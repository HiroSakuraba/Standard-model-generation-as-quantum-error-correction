"""Five-face SU(3) native-Hamiltonian-compatible code search, explicitly limited.

Source: Chen, Mueller, Yao, Phys Rev D 114, 074502 (2026).
Works in the complete 3**5=243 zero-exterior-flux qutrit gauge sector.
Scans unordered *pairs* of the first 18 energy eigenstates at g²=1.
No globally optimal code, certified continuum result or hardware decoder claim.

Reproduce: OPENBLAS_NUM_THREADS=1 python su3_native_eigen_encoding.py
Requires sympy, numpy, scipy and parent repo modules.
"""
from __future__ import annotations
from itertools import combinations
import argparse
import json
from pathlib import Path
from su3_three_plaquette_chain import verify_source_tables
import numpy as np
from scipy.linalg import eigh
from su3_four_five_plaquette import dense_matrices, electric_links, one_site_actions
from su3_channel_recovery import (
    source_data, syndrome_decoder, recover_fidelity, wilson_jump,
    evolved_two_event_channel,
)

N=5
M=18
G2=1.0
TIME=0.6


def problem():
    labels,diag,B = dense_matrices(N)
    pos={q:i for i,q in enumerate(labels)}
    dim=len(labels)
    assert dim==243
    electric=np.asarray(diag,float)
    magnetic=sum(B)
    def hamiltonian(g2):
        # Drop 6N/g²*I: global phase never affects the code/channel.
        return g2*np.diag(electric)-magnetic/g2
    local=[];names=[]
    local.append(np.eye(dim));names.append('I')
    for e in range(3*N+1):
        # E²/(4/3): occupation projector, operator norm 1.
        occ=np.array([float(electric_links(q)[e])*3/4 for q in labels])
        assert set(np.unique(occ)) <= {0.,1.}
        local.append(np.diag(occ));names.append(f'n{e}')
    for z,b in enumerate(B):
        assert np.linalg.norm(b/2,2)<=1+1e-12
        local.append(b/2);names.append(f'B{z}/2')
    wilson={}
    for z in range(N):
        w=np.zeros((dim,dim))
        for col,q in enumerate(labels):
            for dest,a in one_site_actions(q,z,'W').items():
                w[pos[dest],col]=float(a)
        assert np.linalg.norm(w,2)<=1+1e-12
        assert np.array_equal(w+w.T,B[z])
        local += [w,w.T];names += [f'W{z}',f'W{z}^dag']
        wilson[z]=w
    assert len(local)==1+(3*N+1)+N+2*N==32
    return labels,hamiltonian,local,names,wilson


def scan_pairs(H,errors,names):
    eig,V=eigh(H,check_finite=False)
    assert np.max(np.abs(H@V[:,:M]-V[:,:M]*eig[:M]))<1e-11
    U=V[:,:M]
    images=np.stack([a@U for a in errors])
    overlap=np.einsum('adi,bdj->abij',images.conj(),images,optimize=True)
    out=[]
    for i,j in combinations(range(M),2):
        # Every KL pair error overlap K†E_a†E_bK, including E=identity.
        # Center a 2x2 block by subtracting (trace/2) I_2.
        xx=(overlap[:,:,i,i]-overlap[:,:,j,j])/2
        off01=overlap[:,:,i,j];off10=overlap[:,:,j,i]
        defects=np.sqrt(2*np.abs(xx)**2+np.abs(off01)**2+np.abs(off10)**2)
        worst_id=np.unravel_index(np.argmax(defects),defects.shape)
        out.append({'pair':[i,j],'delta_energy':float(eig[j]-eig[i]),
                    'KL_max_centered_Frobenius':float(defects[worst_id]),
                    'worst_pair':[names[worst_id[0]],names[worst_id[1]]]})
    assert len(out)==153 and all(r['KL_max_centered_Frobenius']>1e-6 for r in out)
    best_kl=min(out,key=lambda x:x['KL_max_centered_Frobenius'])
    best_gap=min(out,key=lambda x:x['delta_energy'])
    assert tuple(best_kl['pair'])==(12,15) and tuple(best_gap['pair'])==(15,16)
    assert 0.13<best_kl['KL_max_centered_Frobenius']<0.15
    assert 1e-6<best_gap['delta_energy']<1e-4
    frontier=[]
    for r in out:
        if not any(s is not r and s['delta_energy']<=r['delta_energy']
                   and s['KL_max_centered_Frobenius']<=r['KL_max_centered_Frobenius']
                   and (s['delta_energy']<r['delta_energy']
                        or s['KL_max_centered_Frobenius']<r['KL_max_centered_Frobenius'])
                   for s in out):
            frontier.append(r)
    frontier.sort(key=lambda x:x['delta_energy'])
    for r in out:
        if tuple(r['pair']) == (0,1):first=r
    print('SCAN',json.dumps({'pairs':len(out),'candidates':{'best_KL':best_kl,'smallest_gap':best_gap,'low_energy_01':first}},sort_keys=True))
    return eig,V,{'best_KL':best_kl,'smallest_gap':best_gap,'low_energy_01':first,'candidate_pairs_tested':len(out),'pareto_frontier':frontier,'all_pairs':out}


def projector_decoder_fidelity(kraus_images,K,frame=None):
    """CPTP recovery rho -> K†rho K + tr((I-KK†)rho) I_2/2.

    For eigenpairs optionally rotate the OUTPUT logical qubit by a fixed,
    nominally calibrated virtual Z frame. This is NOT a free hardware gate.
    """
    if frame is None:frame=np.eye(2,dtype=complex)
    norm=sum(float(np.vdot(m,m).real) for m in kraus_images)
    assert abs(norm-2)<2e-9,norm
    coherent=0.;in_code=0.
    for m in kraus_images:
        a=frame@(K.conj().T@m)
        coherent += abs(np.trace(a))**2/4
        in_code += np.linalg.norm(a)**2
    lost=max(0.,norm-in_code)
    fe=coherent+lost/8
    assert -1e-12<=fe<=1+1e-10,(fe,lost)
    return float(fe),float(lost/2)


def unitary_at(H,t,K):
    eig,vec=eigh(H,check_finite=False)
    return (vec*np.exp(-1j*t*eig))@(vec.conj().T@K)


def two_event_images(H,K,W,p,t):
    e,U=eigh(H,check_finite=False)
    S=(U*np.exp(-1j*e*t/3))@U.T
    q=np.sum(W**2,axis=0)
    assert max(q)<=1+1e-11
    kraus=[np.diag(np.sqrt(1-p*q)),np.sqrt(p)*W]
    assert np.linalg.norm(sum(a.T@a for a in kraus)-np.eye(len(q)))<1e-12
    return [S@b@S@a@S@K for a in kraus for b in kraus]


def main():
    verify_source_tables()
    labels,makeH,errors,names,W=problem()
    H=makeH(G2)
    eig,V,scanned=scan_pairs(H,errors,names)
    labels2,pos,K_flux=source_data()
    assert tuple(labels)==tuple(labels2)
    syndrome,_images,_gram=syndrome_decoder(labels2,pos,K_flux)
    results=[]
    for tag,row in [('lowest','low_energy_01'),('best_KL','best_KL'),('closest','smallest_gap')]:
        i,j=scanned[row]['pair']
        K=V[:,[i,j]]
        assert np.linalg.norm(K.conj().T@K-np.eye(2))<1e-13
        frame=np.diag(np.exp(1j*TIME*eig[[i,j]]))
        U_nominal=unitary_at(H,TIME,K)
        f_native,_=projector_decoder_fidelity([U_nominal],K)
        f_tracked,out=projector_decoder_fidelity([U_nominal],K,frame)
        assert abs(f_native-np.cos(scanned[row]['delta_energy']*TIME/2)**2)<2e-11
        assert abs(f_tracked-1)<2e-11
        # Two independent equally likely quasi-static calibration drifts.
        drift={}
        for uncertainty in (0.02,0.05,0.10):
            outputs=[unitary_at(makeH(G2*(1+sign*uncertainty)),TIME,K)/np.sqrt(2)
                     for sign in (-1,1)]
            f,leak=projector_decoder_fidelity(outputs,K,frame)
            drift[str(uncertainty)]={'Fe':f,'out_of_selected_pair':leak}
        jumps={}
        for face in (3,4):
            maps=two_event_images(H,K,W[face],0.10,TIME)
            f,leak=projector_decoder_fidelity(maps,K,frame)
            jumps[str(face)]={'Fe':f,'out_of_selected_pair':leak}
        # Exact local electric occupation jump (one selected link): since n_e^2=n_e,
        # both no-jump and jump Kraus lie in span{I,n_e}.
        n4=errors[1+4]
        electric_local={}
        for t in (0.0,TIME):
            maps=two_event_images(H,K,n4,0.10,t)
            f,leak=projector_decoder_fidelity(maps,K,frame if t else None)
            electric_local[str(t)]={'Fe':f,'out_of_selected_pair':leak}
        result={'name':tag,'pair':[i,j],'delta_energy':scanned[row]['delta_energy'],
                'energy_above_ground':[float(eig[i]-eig[0]),float(eig[j]-eig[0])],
                'KL_max_centered_Frobenius':scanned[row]['KL_max_centered_Frobenius'],
                'no_fault_t06_no_frame_Fe':f_native,'no_fault_t06_virtual_frame_Fe':f_tracked,
                'drift_t06':drift,'two_W_jumps_t06':jumps,
                'two_link4_electric_jumps':electric_local}
        results.append(result)
        print('EIGENPAIR',json.dumps(result,sort_keys=True))
    # Independent regression: previously published *flux* syndrome decoder,
    # physically uncorrected native evolution, NO ideal inverse.
    maps=two_event_images(H,K_flux,W[3],0.,TIME)
    f_flux,raw_flux,out_flux=recover_fidelity(maps,syndrome,K_flux)
    assert abs(f_flux-0.6327487288)<1e-5,(f_flux,raw_flux)
    # Under the same two electric-only jumps at t=0, the original flux code's
    # five-syndrome decoder should still be exact by the prior KL certificate.
    electric_maps=two_event_images(H,K_flux,errors[1+4],0.10,0.0)
    f_electric,raw_electric,out_electric=recover_fidelity(electric_maps,syndrome,K_flux)
    assert abs(f_electric-1)<3e-11
    baseline={'flux_syndrome_F':f_flux,'flux_raw_F':raw_flux,'flux_syndrome_outside':out_flux,
              'flux_t0_two_electric_jumps_recovered':f_electric}
    print('FLUX_BASELINE',json.dumps(baseline,sort_keys=True))
    print('EXACT INPUT ALGEBRA/REGRESSION + NUMERICAL EIGENPAIR CHANNEL COMPARISONS PASS')
    result={'source_citation':'Chen-Mueller-Yao PRD 114 074502 (2026)',
          'dimensions':len(labels),'scan':scanned,'evaluated':results,'baseline':baseline,
          'limitations':['First 18 energy eigenstate pairs only, not all subspaces',
                         'KL metric uses specified 32 dimensionless normalized errors',
                         'Projected decoder and perfect virtual Z frame are abstract CPTP/known-calibration operations',
                         'No modeled preparation/readout, spontaneous decay or continuum glueball physics']}
    print('JSON_RESULTS='+json.dumps(result,sort_keys=True))
    return result

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--output-json',help='Optional reproducible results JSON path')
    args=parser.parse_args()
    result=main()
    if args.output_json:
        path=Path(args.output_json)
        path.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n',encoding='utf-8')
        print('WROTE_JSON',str(path))
