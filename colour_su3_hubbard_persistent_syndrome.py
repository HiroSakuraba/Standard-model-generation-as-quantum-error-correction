"""Persistent occupation leakage in the charged-bus S3 logical syndrome.

Four three-colour active fermion orbitals A0,A1,B0,M with 4 fermions:
choose(12,4)=495 states, INCLUDING ALL intermediate doublon/hole states.
B1 remains a spectator colour qutrit. One group-valued S3 gauge link,
mobile charge location A/B, stationary conjugate reference, and neutral
three-level pointer are also explicit. Charged mediator transport acts on
any mediator occupation by the second-quantized signed-permutation D(g).

Twelve pointer-detuned finite-U Hubbard local pulses replace all 12
ideal pointer-controlled charged-matter swaps of PR53. Four exact link
transports remain; five ideal feedback operations are separately stated.
No occupation sector is reset or projected between gates. The resulting
entanglement fidelity is unconditioned, counting final outside-code
population as loss. This is a microscopic CONTROL module, not full SU3
lattice gauge dynamics or original nine-quark code recovery.
"""
from __future__ import annotations
from functools import lru_cache
from itertools import combinations
import json
from pathlib import Path
import numpy as np
from scipy.linalg import eigh
from scipy.sparse import csr_matrix

from colour_gauge_link_s3 import (
    G, DIM, MATTER, PERMS, REP, compose, inverse, physical_isometry,
    dressed_exchange,
)
from colour_s3_edge_logical_syndrome import code_isometry
from colour_su3_hubbard_virtual_exchange import optimize_pulse

N_FERMIONS=4
N_MODES=12
F_STATES=tuple(sum(1<<k for k in c) for c in combinations(range(N_MODES),N_FERMIONS))
F_INDEX={state:i for i,state in enumerate(F_STATES)}
N_FOCK=len(F_STATES)
assert N_FOCK==495
# orbital 0=A0, 1=A1, 2=B0, 3=M. B1 spectator in own qutrit.
SITE={'A0':0,'A1':1,'B0':2}
SCHEDULE=(("edge",(1,2)),("onsite",(1,2)),("edge",(2,)),("onsite",(2,)))


def fock_index(a0,a1,b0,m):
    return F_INDEX[(1<<a0)|(1<<(3+a1))|(1<<(6+b0))|(1<<(9+m))]


def hopp_term(s,src,dst):
    if not (s>>src)&1 or (s>>dst)&1:return None
    v=(-1)**((s&((1<<src)-1)).bit_count())
    emptied=s^(1<<src)
    v*=(-1)**((emptied&((1<<dst)-1)).bit_count())
    return emptied|(1<<dst),v


@lru_cache(maxsize=12)
def local_hamiltonian_blocks(site,U,delta,t=1.):
    """Exact shared 495D Hamiltonian in sectors of spectator occupations.

    Restricting U to each spectator block (at most 20D) makes repeated
    finite-time propagators efficient without projecting away doublons.
    """
    if site not in SITE or not (U>abs(delta) and t>0):raise ValueError('Hubbard inputs')
    q=SITE[site]
    pair_mask=sum(1<<(3*orb+k) for orb in (q,3) for k in range(3))
    groups={}
    for index,s in enumerate(F_STATES):
        groups.setdefault(s&~pair_mask,[]).append(index)
    output=[]
    for spec,rows in groups.items():
        pos={F_STATES[i]:k for k,i in enumerate(rows)}
        dim=len(rows);H=np.zeros((dim,dim),complex)
        for col,index in enumerate(rows):
            state=F_STATES[index]
            nq=sum((state>>(3*q+a))&1 for a in range(3))
            nm=sum((state>>(9+a))&1 for a in range(3))
            H[col,col]=U*(nq*(nq-1)+nm*(nm-1))/2+delta*(nq-nm)/2
            for a in range(3):
                for src,dst in ((3*q+a,9+a),(9+a,3*q+a)):
                    term=hopp_term(state,src,dst)
                    if term:
                        new,sign=term
                        H[pos[new],col]-=t*sign
        assert np.linalg.norm(H-H.conj().T)<1e-12
        eig,vec=eigh(H)
        output.append((np.array(rows,dtype=int),eig,vec))
    assert sum(len(rows) for rows,_,_ in output)==N_FOCK
    return tuple(output)


@lru_cache(maxsize=40)
def local_fock_pulse(site,U,delta,duration,t=1.):
    """Exact 495x495 sparse full-space unitary, no occupation reset."""
    r=[];c=[];vals=[]
    for indices,energies,V in local_hamiltonian_blocks(site,U,delta,t):
        evolved=(V*np.exp(-1j*duration*energies))@V.conj().T
        nz=np.where(abs(evolved)>1e-14)
        r.extend(indices[nz[0]].tolist())
        c.extend(indices[nz[1]].tolist())
        vals.extend(evolved[nz].tolist())
    out=csr_matrix((vals,(r,c)),shape=(N_FOCK,N_FOCK))
    assert np.linalg.norm((out.getH()@out-csr_matrix(np.eye(N_FOCK))).toarray())<2e-10
    return out


def subgroup_wedge_map(h,orbital_subset):
    """Gauge rotation of occupied colours, as the fermionic exterior power.

    D(h) is a determinant-one SIGNED PERMUTATION; simple exact signed
    Fock-basis permutations handle n=0,1,2,3 and all shared leakage.
    """
    perm=PERMS[h]
    phase=(-1)**sum(perm[i]>perm[j] for i in range(3) for j in range(i+1,3))
    mapped=np.empty(N_FOCK,int); coeff=np.empty(N_FOCK,int)
    for i,s in enumerate(F_STATES):
        rest=s;sign=1
        for orb in orbital_subset:
            cols=[a for a in range(3) if (s>>(3*orb+a))&1]
            for a in cols:rest^=1<<(3*orb+a)
            renamed=[perm[a] for a in cols]
            sign*=phase**len(cols)
            sign*=(-1)**sum(renamed[a]>renamed[b]
                for a in range(len(renamed)) for b in range(a+1,len(renamed)))
            for a in renamed:rest|=1<<(3*orb+a)
        mapped[i]=F_INDEX[rest];coeff[i]=sign
    return mapped,coeff


@lru_cache(maxsize=24)
def _maps(h,subset):
    return subgroup_wedge_map(h,tuple(subset))


def initialize(K):
    """(pointer,link,B1,mediator-location,Fock,reference,logical)."""
    assert K.shape[0]==DIM
    R=K.shape[1]
    z=np.zeros((3,G,3,2,N_FOCK,3,R),complex)
    bits=K.reshape(G,3,3,3,3,R)
    for a0 in range(3):
        for a1 in range(3):
            for b0 in range(3):
                for m in range(3):
                    i=fock_index(a0,a1,b0,m)
                    for b1 in range(3):
                        # 1/sqrt(3) mediator-reference singlet and
                        # 1/sqrt(3) uniform pointer.
                        z[:,:,b1,0,i,m,:]=bits[:,a0,a1,b0,b1,:][None,:,:]/3
    assert abs(np.vdot(z,z).real-np.vdot(K,K).real)<1e-10
    return z


def apply_microscopic_gate(z,site,active,U,detuning,T):
    """Exact full-Fock Hubbard pulse at the mediator's present vertex.

    Other mediator-location sector remains untouched. Pointer levels
    not active see delta_off=0, not a mathematically ideal identity.
    """
    loc=0 if site.startswith('A') else 1
    output=z.copy()
    for p in range(3):
        d=detuning if p in active else 0.
        gate=local_fock_pulse(site,U,d,T)
        block=z[p,:,:,loc,:,:,:]
        # (G,B1,Fock,ref,logical); sparse dot along Fock axis.
        moved=np.moveaxis(block,2,0)
        operated=(gate@moved.reshape(N_FOCK,-1)).reshape(moved.shape)
        output[p,:,:,loc,:,:,:]=np.moveaxis(operated,0,2)
    return output


def link_transport(z):
    """Extend exact PR53 charged-bus transport to ANY mediator occupancy.

    With multiple mediator fermions, D(g) acts through the exterior
    power of its fundamental representation; all leaked sectors kept.
    """
    out=np.empty_like(z)
    for g in range(G):
        for loc in (0,1):
            h=inverse(g) if loc==0 else g
            ix,phase=_maps(h,(3,))
            target=out[:,g,:,1-loc,:,:,:]
            source=z[:,g,:,loc,:,:,:]
            target[:,:,ix,:,:]=phase[None,None,:,None,None]*source
    return out




def gauge_transform(z,h,vertex):
    """Full finite-S3 local gauge action INCLUDING all virtual Fock states.

    Gauges all charged occupied colours through exterior powers of D(h),
    including arbitrarily doubled/empty mediator orbitals. The B1 matter
    qutrit and stationary conjugate reference transform separately.
    """
    if vertex not in ('A','B'):raise ValueError(vertex)
    perm=PERMS[h]
    phase=(-1)**sum(perm[i]>perm[j] for i in range(3) for j in range(i+1,3))
    out=np.zeros_like(z)
    for loc in (0,1):
        selected=((0,1,3) if loc==0 else (0,1)) if vertex=='A' else (
                  (2,3) if loc==1 else (2,))
        ix,sign=_maps(h,selected)
        for g in range(G):
            gd=compose(h,g) if vertex=='A' else compose(g,inverse(h))
            for b1 in range(3):
                bdest=(perm[b1] if vertex=='B' else b1)
                bsign=phase if vertex=='B' else 1
                for ref in range(3):
                    rdest=(perm[ref] if vertex=='A' else ref)
                    rsign=phase if vertex=='A' else 1
                    target=out[:,gd,bdest,loc,:,rdest,:]
                    source=z[:,g,b1,loc,:,ref,:]
                    target[:,ix,:]=sign[None,:,None]*bsign*rsign*source
    return out


def gauss_acceptance(z):
    """Acceptance of the exact full 4-fermion local finite-group Gauss projectors."""
    norm=np.vdot(z,z).real
    if norm==0:raise ValueError('zero input')
    values={}
    for vertex in ('A','B'):
        proj=sum((gauge_transform(z,h,vertex) for h in range(G)))/G
        values[vertex]=float(np.vdot(proj,proj).real/norm)
    return values

def microscopic_schedule():
    ops=[]
    for kind,active in SCHEDULE:
        if kind=='edge':
            ops.extend([('A0',active),('T',None),('B0',active),
                        ('T',None),('A0',active)])
        else:
            ops.extend([('A0',active),('A1',active),('A0',active)])
    assert len(ops)==16
    return ops


def exact_full_channel(K,U,detuning,T,feedback='ideal'):
    """Evaluate entire syndrome, no per-pulse projection/postselection.

    Final code fidelity counts only fully computational four-fermion
    occupation states; all other occupation configurations are
    orthogonal lost/erased states. Feedback uses ideal dressed SWAP
    on code sectors only (optimistic upper bound), and DOES NOT undo
    leaked states: a faithful microscopic feedback is future work.
    """
    if feedback!='ideal':raise ValueError('Only ideal subspace feedback yet')
    results=[]
    source_p=.37
    for E in (np.sqrt(1-source_p)*K,np.sqrt(source_p)*dressed_exchange(K)):
        z=initialize(E)
        trace=[]
        for op,active in microscopic_schedule():
            if op=='T': z=link_transport(z)
            else:z=apply_microscopic_gate(z,op,active,U,detuning,T)
            trace.append(float(occupation_leakage(z)))
        # pointer inverse QFT and classical readout
        z=np.fft.fft(z,axis=0)/np.sqrt(3)
        # For branch label 2, apply ideal logical dressed SWAP to each
        # 486-vector inside singly-occupied full-Fock sectors; zero
        # occupation leakage is NOT artificially restored.
        results.append((z,trace))
    Fe=0.; total=0.;final_leakage=0.;syndrome=np.zeros(3)
    for source_index,(z,trace) in enumerate(results):
        source_weight=(1-source_p if source_index==0 else source_p)
        total+=np.vdot(z,z).real/4
        syndrome+=np.array([np.vdot(q,q).real/4 for q in z])
        # Build code-computational branch for each mediator colour and
        # charge reference, including wrong-charge environment outcomes.
        for label in range(3):
            for loc in (0,):  # ideal recovered charged bus must return to A
                for ref in range(3):
                    for med in range(3):
                        projected=project_matter(z[label,:,: ,loc,:,:,:],med,ref)
                        if label==2:
                            projected=dressed_exchange(projected)
                        Fe+=abs(np.trace(K.conj().T@projected))**2/16
        final_leakage+=source_weight*trace[-1]
    assert abs(total-1)<3e-10
    return {'unconditional_F_e_ideal_feedback':float(Fe),
        'mean_final_occupation_leakage':float(final_leakage),
        'mean_max_transient_occupation_leakage':float((1-source_p)*max(results[0][1])+source_p*max(results[1][1])),
        'pointer_probabilities':syndrome.tolist(),
        'total_probability':float(total),
        'source_swap_probability':source_p,
        'feedback_scope':'Ideal colour-subspace conditional dressed SWAP only; leaked Fock sectors never recovered'}


def project_matter(arr,med,ref):
    """From one pointer/loc branch, project active Fock to occupied code.

    Does not normalize. arr shape (link,B1,Fock,reference,logical).
    """
    R=arr.shape[-1]
    target=np.zeros((G,3,3,3,3,R),complex)
    for a0 in range(3):
        for a1 in range(3):
            for b0 in range(3):
                ix=fock_index(a0,a1,b0,med)
                target[:,a0,a1,b0,:,:]=arr[:,:,ix,ref,:]
    return target.reshape(DIM,R)


def occupation_leakage(z):
    """Fraction outside N_A0=N_A1=N_B0=N_M=1 at any pointer/loc."""
    ok=np.array([all(sum((s>>(3*k+c))&1 for c in range(3))==1
                     for k in range(4)) for s in F_STATES])
    return 1-float(np.vdot(z[:,:,:,:,ok,:,:],z[:,:,:,:,ok,:,:]).real/
                   np.vdot(z,z).real)


def run():
    W=physical_isometry()
    K,F,_=code_isometry(W)
    physical_initial=initialize(K)
    initial_gauss=gauss_acceptance(physical_initial)
    assert min(initial_gauss.values())>1-1e-10
    arms=[]
    for U in (20.,30.,40.):
        opt=optimize_pulse(U)
        T=opt['time'];d=opt['active_detuning']
        metrics=exact_full_channel(K,U,d,T)
        arms.append({'U_over_t':U,'single_pulse_F_e':opt['optimized_metrics']['entanglement_fidelity'],
                     'duration_in_inverse_hopping_t':T,'detuning':d,
                     'complete_syndrome_ideal_feedback_reference':metrics})
    return {'initial_extended_Gauss_acceptance':initial_gauss,
            'active_fermionic_Fock_dimension':N_FOCK,
            'full_register_dimensions':{'pointer':3,'link':G,'spectator_B1':3,
                'mediator_location':2,'four_fermion_active_Fock':N_FOCK,
                'charge_reference':3},
            'all_four_active_orbital_occupations_persist_across_16_operations':True,
            'Hubbard_local_exchange_pulses':12,'gauge_link_transports':4,
            'source_channel':'known dressed swap with p=0.37',
            'arms':arms,'limitations':[
              'Feedback is an ideal code-colour dressed SWAP reference, not a synthesized Hubbard feedback; doublon/hole occupation leakage is not recovered',
              'Two-site Hubbard pulse acts only when mediator location agrees with local matter; link transport is exact wedge representation without physical hopping pulses',
              'Full Gauss covariance of all Fock states must be separately checked before a gauge-invariance claim',
              'No particle loss or dephasing evolution, no ramp or pulse amplitude constraint; four-fermion loss survival is only an assumed factor',
              'SU3 Hubbard hopping is onsite/common-vertex only; finite S3 links are distinct gauge vertices',
              'Not the original nine-quark 21D/65-syndrome recovery code']}


def main():
    out=run()
    Path('colour_su3_hubbard_persistent_syndrome_results.json').write_text(
         json.dumps(out,indent=2,sort_keys=True)+'\n')
    for x in out['arms']:
        print('PERSISTENT_FOCK_SYNDROME',json.dumps({'U_over_t':x['U_over_t'],
             'single_pulse_F_e':x['single_pulse_F_e'],
             'whole_syndrome':x['complete_syndrome_ideal_feedback_reference']},sort_keys=True))

if __name__=='__main__':main()