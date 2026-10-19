"""Hardware-aware LOGICAL qutrit gate compilation for the published SU(3) chain.

An exact conditioned 3x3 local term is decomposed into adjacent-level Givens
rotations, *conditional on neighboring qutrit state(s)*. These are logical
multi-qutrit operations, NOT calibrated native transmon entangling gates.
A Strang product formula provides a real finite-depth approximation to the
five-face Z3 refocusing schedule. No hardware leakage is simulated.
"""
import json
from collections import Counter
from itertools import product
from functools import lru_cache
import numpy as np
from scipy.linalg import eigh
from su3_three_plaquette_chain import PUBLISHED
from su3_native_control_compilation import construct, toggled_matrices, evolution, THETA, entanglement_fidelity, T
from su3_channel_recovery import syndrome_decoder, source_data

N=5
MASK=(0,2,3,4)
CYCLES=8
ORDER=(0,1,2,1,0)
WEIGHTS=(1/6,1/6,1/3,1/6,1/6)


def givens_adjacent(U):
    """Exact QR synthesis: U=G_12† G_01† G_12† D, three adjacent-level rotations.

    The conditional implementation of these Givens rotations requires a
    further multi-qutrit compilation; do not call them native hardware gates.
    """
    Q=np.asarray(U,dtype=complex).copy()
    assert Q.shape==(3,3) and np.linalg.norm(Q.conj().T@Q-np.eye(3))<1e-11
    elim=[]
    for a,b,col in ((1,2,0),(0,1,0),(1,2,1)):
        x,y=Q[a,col],Q[b,col]
        r=np.hypot(abs(x),abs(y))
        if r<1e-13:continue
        g=np.array([[np.conj(x),np.conj(y)],[-y,x]],complex)/r
        G=np.eye(3,dtype=complex);G[np.ix_((a,b),(a,b))]=g
        Q=G@Q
        elim.append(G)
    assert np.max(np.abs(Q-np.diag(np.diag(Q))))<4e-12
    diag=np.diag(np.diag(Q))
    assert np.linalg.norm(diag.conj().T@diag-np.eye(3))<4e-12
    reconstructed=diag
    for G in reversed(elim):reconstructed=G.conj().T@reconstructed
    assert np.linalg.norm(reconstructed-U)<4e-12
    return len(elim)


def local_loop_block(neighbor_left,neighbor_right,phase):
    """Hermitian local magnetic block for fixed neighbor qutrit values."""
    M=np.asarray(PUBLISHED[(neighbor_left,neighbor_right)],dtype=float)
    W=np.zeros((3,3),dtype=complex)
    for a in range(3):W[(a+1)%3,a]=M[(a+1)%3,a]
    B=np.exp(1j*phase)*W+np.exp(-1j*phase)*W.conj().T
    assert np.linalg.norm(B-B.conj().T)<1e-13
    return B


def index_groups(states):
    pos={q:i for i,q in enumerate(states)}
    grouped={}
    for site in range(N):
        blocks=[]
        for other in product(range(3),repeat=N-1):
            vals=[]; k=0
            for i in range(N):
                if i==site:vals.append(0)
                else: vals.append(other[k]);k+=1
            qs=[]
            for a in range(3):
                q=list(vals);q[site]=a;qs.append(pos[tuple(q)])
            l=vals[site-1] if site else 0
            r=vals[site+1] if site<N-1 else 0
            blocks.append((np.array(qs,dtype=int),l,r))
        assert len(blocks)==3**(N-1)
        assert len(set(int(i) for ids,l,r in blocks for i in ids))==len(states)
        grouped[site]=blocks
    return grouped


def crosscheck_all_blocks(states,W,blocks):
    for site in range(N):
        w=W[site]
        for ix,l,r in blocks[site]:
            M=np.asarray(PUBLISHED[(l,r)],dtype=float)
            for a in range(3):
                b=(a+1)%3
                assert abs(w[ix[b],ix[a]]-M[b,a])<1e-13
        # 243 basis transitions per site and 5 sites, including cases 1/9.
    return 3**N * N


def local_unitaries(dt,site,phase,drift=0):
    """Cache 3x3 local unitary per (left,right) sector; coherent drive drift.

    Drift is *deterministic* relative amplitude bias by plaquette, NOT
    stochastic calibrated error, leakage or full transmon dynamics.
    """
    out={}
    for l,r in product(range(3),repeat=2):
        gain=1+drift*((-1)**(site+l+2*r))
        B=local_loop_block(l,r,phase)*gain
        val,vec=eigh(B)
        U=(vec*np.exp(1j*dt*val))@vec.conj().T
        assert givens_adjacent(U)<=3
        out[l,r]=U
    return out


def strang_segments(K,states,d,blocks,mask=MASK,cycles=CYCLES,rounds=1,drift=0):
    """Second-order product formula for each of 5 toggling stages/cycle.

    Physical duration is T; each interval subdivided 'rounds' times.
    Each elementary site unitary is a NEIGHBOR-CONDITIONED U(3), with
    no assumption of a free 3-body native gate.
    """
    y=np.asarray(K,complex).copy()
    total=Counter()
    for _ in range(cycles):
        for toggle,weight in zip(ORDER,WEIGHTS):
            dt=T*weight/(cycles*rounds)
            phase={i:(THETA*toggle if i in mask else 0.) for i in range(N)}
            exp_e=np.exp(-.5j*dt*d)
            schedule=list(range(N))+list(range(N-2,-1,-1))
            cache={}
            for site in range(N):
                cache[site,.5]=(local_unitaries(dt/2,site,phase[site],drift))
                if site==N-1:cache[site,1.]=(local_unitaries(dt,site,phase[site],drift))
            for _ in range(rounds):
                y=exp_e[:,None]*y
                for site in schedule:
                    units=cache[site,1. if site==N-1 else .5]
                    for ix,l,r in blocks[site]:
                        y[ix,:]=units[l,r]@y[ix,:]
                    total['site_unitaries']+=1
                    total['context_blocks']+=3 if site in (0,N-1) else 9
                y=exp_e[:,None]*y
                total['electric_phases']+=2
            total['stage_count']+=1
    assert np.linalg.norm(y.conj().T@y-np.eye(2))<2e-10
    total['givens_max']=3*total['context_blocks']
    total['context_diagonals']=total['context_blocks']
    total['three_qutrit_context_blocks']=sum(1 for _ in ())*0 + (total['site_unitaries']-2*2*cycles*5*rounds)*9  # 2 boundary sites 2 occurrences each
    return y, dict(total)


def main():
    states,K,d,W,H=construct()
    blocks=index_groups(states)
    transitions=crosscheck_all_blocks(states,W,blocks)
    assert transitions==1215
    # Distinct coefficient values / all 99 conditioned oriented transition types.
    coeffs={str(x) for a in PUBLISHED.values() for x in a if x!=0}
    assert len(PUBLISHED)==9
    counts={}
    for i in range(N):
        contexts=3 if i in (0,N-1) else 9
        counts[str(i)]={'neighbor_contexts':contexts,'directed_two_level_terms':3*contexts,
                        'control_arity':1+(i>0)+(i<N-1)}
    assert sum(x['directed_two_level_terms'] for x in counts.values())==99
    for phase in (0,THETA,2*THETA):
        U=local_unitaries(.019,2,phase)
        assert len(U)==9
    decoder,_,gram=syndrome_decoder(*source_data())
    assert gram.rank()==5
    ideal=evolution(K,toggled_matrices(d,W,MASK),T,CYCLES)
    ideal_raw=entanglement_fidelity([ideal],K)
    ideal_dec=entanglement_fidelity([ideal],K,decoder)
    rows=[]
    for rounds in (1,2,4):
        y,cost=strang_segments(K,states,d,blocks,rounds=rounds)
        err=np.linalg.norm(y-ideal,ord='fro')
        row={'trotter_rounds_per_stage':rounds,'matrix_state_error_frobenius':float(err),
             'F_raw':entanglement_fidelity([y],K),
             'F_fixed_decoder':entanglement_fidelity([y],K,decoder),
             **cost,
             'no_event_p_1e-4_per_site_unitary':(1-1e-4)**cost['site_unitaries'],
             'no_event_p_1e-3_per_site_unitary':(1-1e-3)**cost['site_unitaries']}
        for level,count in [('site_unitaries',cost['site_unitaries']),('context_blocks',cost['context_blocks']),('givens_max',cost['givens_max'])]:
            for rate in (1e-5,1e-4,1e-3):
                survival=(1-rate)**count
                # Explicit effective CPTP toy: any failure is a flagged event
                # followed by completely depolarizing replacement of logical qubit.
                # THIS IS NOT an experimentally measured gate infidelity model.
                row[f'{level}_p{rate:g}_event_free']=survival
                row[f'{level}_p{rate:g}_replaced_qubit_F']=survival*row['F_fixed_decoder']+(1-survival)/4
        rows.append(row)
        print('COMPILER',json.dumps(row,sort_keys=True),flush=True)
    assert rows[-1]['matrix_state_error_frobenius']<rows[0]['matrix_state_error_frobenius']
    driftrow=[]
    for fraction in (.01,.05):
        y,cost=strang_segments(K,states,d,blocks,rounds=2,drift=fraction)
        row={'fixed_context_gain_bias':fraction,'F_raw':entanglement_fidelity([y],K),
             'F_fixed_decoder':entanglement_fidelity([y],K,decoder)}
        driftrow.append(row);print('DRIFT',json.dumps(row,sort_keys=True),flush=True)
    data={'source':'Chen Mueller Yao PRD 114 074502 Table VI; abstract qutrit hardware compiler',
          'n_qutrits':N,'dimension':3**N,'mask':list(MASK),'cycles':CYCLES,'time':T,
          'published_conditioned_transition_types':99,'published_contexts_by_site':counts,
          'coefficient_values':sorted(coeffs),'crosschecked_transitions':transitions,
          'syndrome_measurement_outcomes':6,'minimum_pointer_ancilla_qutrits':2,'minimum_readout_bits':3,
          'gate_fault_model':'independent flagged catastrophic error probability p per assumed logical operation; replaced encoded qubit I/2 after any fault, strictly a CPTP illustrative model',
          'ideal_sparse_expm_F_raw':ideal_raw,'ideal_sparse_expm_F_fixed_decoder':ideal_dec,
          'trotter_rows':rows,'context_gain_bias_rows':driftrow,
          'warning':'conditional qutrit unitaries and controlled Givens NOT synthesized into calibrated two-transmon native gates; no physical leakage/T1/T2 simulation'}
    with open('su3_qutrit_gate_compiler_results.json','w') as f:json.dump(data,f,indent=2,sort_keys=True)
    print('ALL CONDITIONAL QUTRIT BLOCK/GIVENS/STRANG VALIDATIONS PASSED')


if __name__=='__main__':main()
