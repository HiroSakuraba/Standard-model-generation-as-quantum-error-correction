"""Exact two-register synthesis of the finite-S3 logical-syndrome circuit.

The high-level PR48 phase-pointer gadget uses two pointer-controlled
dressed edge swaps and two pointer-controlled onsite swaps. Each
controlled dressed swap has an EXACT five-two-register factorization
using register pairs (link,B0), (A0,B0), (pointer,B0).
Each controlled onsite swap uses three two-register factors.

These are GENERAL two-register unitaries: a selective pointer–qutrit
phase, an arbitrary two-qutrit eigenbasis transform Q, and an 18x18
link-matter controlled representation C. They are NOT native
two-level entanglers, two-qubit gates, locality-preserving at every
intermediate step, or calibrated hardware primitives.

The 16-factor syndrome gadget is compared on the entire 486D system
against colour_s3_edge_logical_syndrome.pointer_circuit(). Each
intermediate stage is checked for Gauss-constraint leakage.

Noise sensitivity uses a synthetic, exact CPTP stochastic mixed-unitary
two-register Weyl fault channel after each factor. It includes three
more two-register factors for conditional known-pair swap recovery.
This is a finite-sample model, not leakage, coherent control noise or
physical time/distance/routing analysis.
"""
from __future__ import annotations

from functools import lru_cache
from pathlib import Path
import json
import numpy as np

from colour_gauge_link_s3 import (
    DIM, G, MATTER, REP, physical_isometry, dressed_exchange, projector,
    swap_matter,
)
from colour_s3_edge_logical_syndrome import (
    pointer_circuit, code_isometry, W3,
)

REGISTERS=("pointer","link","A0","A1","B0","B1")
DIMS=(3,6,3,3,3,3)
RANK=4


def exchange_eigenbasis_q():
    """Q maps the 3 antisymmetric two-qutrit SWAP eigenvectors to |a,0>.

    SWAP_AB = Q^dag (I-2 I_A tensor |0><0|_B) Q.
    Symmetric eigenstates map to the six B != 0 basis states.
    """
    eye=np.eye(9,dtype=complex)
    anti=[]
    sym=[]
    for i in range(3):
        sym.append(eye[3*i+i])
    for i,j in ((0,1),(0,2),(1,2)):
        anti.append((eye[3*i+j]-eye[3*j+i])/np.sqrt(2))
        sym.append((eye[3*i+j]+eye[3*j+i])/np.sqrt(2))
    out=[3*a for a in range(3)] + [3*a+b for b in (1,2) for a in range(3)]
    Q=np.zeros((9,9),complex)
    for row,vec in zip(out,anti+sym):
        Q[row,:]=vec.conj()
    assert np.linalg.norm(Q.conj().T@Q-np.eye(9))<1e-13
    diag=np.diag([(-1 if j%3==0 else 1) for j in range(9)])
    swap=swap_matter(0,1).reshape((3,3,3,3))
    swap2=np.zeros((9,9),complex)
    for i in range(3):
        for j in range(3):
            swap2[3*j+i,3*i+j]=1
    assert np.linalg.norm(Q.conj().T@diag@Q-swap2)<1e-13
    return Q


@lru_cache(maxsize=1)
def link_b_matrix():
    """C=sum_g |g><g|_link tensor D(g)_B0, 18x18 pair operation."""
    C=np.zeros((18,18),complex)
    for g in range(G):
        C[g*3:(g+1)*3,g*3:(g+1)*3]=REP[g]
    assert np.linalg.norm(C.conj().T@C-np.eye(18))<1e-13
    return C


def pointer_phase_matrix(active):
    """Selective two-register -1 phase when p in active and target digit=0."""
    active=tuple(active)
    if any(k not in (0,1,2) for k in active):
        raise ValueError(active)
    D=np.eye(9,dtype=complex)
    for p in active:D[3*p+0,3*p+0]=-1
    return D


def pair_gate(state,axes,gate):
    """Apply any two-register U(d_i d_j) to 3x486xr or 486xr, no 3-body gate."""
    x=np.asarray(state,complex)
    has_pointer=x.shape[0]==3 and x.ndim==3
    if has_pointer:
        dims=DIMS
    else:
        assert x.shape[0]==DIM and x.ndim==2
        dims=DIMS[1:]
        axes=(axes[0]-1,axes[1]-1)
        assert min(axes)>=0
    assert 0<=axes[0]<len(dims) and 0<=axes[1]<len(dims) and axes[0]!=axes[1]
    assert gate.shape==(dims[axes[0]]*dims[axes[1]],)*2
    arr=x.reshape(dims+(x.shape[-1],))
    front=np.moveaxis(arr,axes,(0,1))
    s=front.shape
    transformed=(gate@front.reshape(s[0]*s[1],-1)).reshape(s)
    back=np.moveaxis(transformed,(0,1),axes)
    return back.reshape(x.shape)


def controlled_dressed_factors(active):
    """C† * Q† * controlled-sign_B0 * Q * C; time order."""
    Q=exchange_eigenbasis_q()
    C=link_b_matrix()
    return [
        ("C_link_B", (1,4), C),
        ("Q_A_B", (2,4), Q),
        ("CP_pointer_B", (0,4), pointer_phase_matrix(active)),
        ("Qdg_A_B", (2,4), Q.conj().T),
        ("Cdg_link_B", (1,4), C.conj().T),
    ]


def controlled_onsite_factors(active):
    Q=exchange_eigenbasis_q()
    return [
        ("Q_A0_A1", (2,3), Q),
        ("CP_pointer_A1", (0,3), pointer_phase_matrix(active)),
        ("Qdg_A0_A1", (2,3), Q.conj().T),
    ]


def measured_syndrome_factors():
    """Two controlled U=A S, one extra controlled U at pointer value 2."""
    return (controlled_dressed_factors((1,2))
            +controlled_onsite_factors((1,2))
            +controlled_dressed_factors((2,))
            +controlled_onsite_factors((2,)))


def correction_factors():
    """S=C† SWAP_AB C, three exact unconditioned two-register gates."""
    C=link_b_matrix()
    swap=np.zeros((9,9),complex)
    for a in range(3):
        for b in range(3):swap[3*b+a,3*a+b]=1
    return [
        ("C_link_B", (1,4), C),
        ("SWAP_A_B", (2,4), swap),
        ("Cdg_link_B", (1,4), C.conj().T),
    ]


def prepared_pointer(y):
    x=np.asarray(y,complex)
    assert x.shape[0]==DIM and x.ndim==2
    return np.stack((x,x,x),axis=0)/np.sqrt(3)


def phase_readout(y):
    """Single-qutrit inverse Fourier; output label m=omega^m eigenphase."""
    return np.fft.fft(y,axis=0)/np.sqrt(3)


def apply_syndrome(y,noise=0.,rng=None,history=False):
    """Physical-register factor schedule, optional CPTP Weyl fault samples."""
    z=prepared_pointer(y)
    states=[]
    for name,pair,gate in measured_syndrome_factors():
        z=pair_gate(z,pair,gate)
        if noise:
            if rng is None: raise ValueError("rng needed when noise is nonzero")
            if rng.random()<noise:
                z=pair_gate(z,pair,sample_two_register_weyl(pair,rng))
        if history:states.append((name,pair,z.copy()))
    result=phase_readout(z)
    return (result,states) if history else result


def local_weyl(d,a,b):
    basis=np.arange(d)
    X=np.eye(d)[(basis-a)%d]   # cyclic shift up a on column states
    Z=np.diag(np.exp(2j*np.pi*b*basis/d))
    U=X@Z
    assert np.linalg.norm(U.conj().T@U-np.eye(d))<1e-12
    return U


def sample_two_register_weyl(pair,rng):
    d1,d2=DIMS[pair[0]],DIMS[pair[1]]
    # Uniform non-identity product generalized Weyl on this pair.
    n=int(rng.integers(1,(d1*d2)**2))
    x2=n%d2;n//=d2
    z2=n%d2;n//=d2
    x1=n%d1;n//=d1
    z1=n%d1
    return np.kron(local_weyl(d1,x1,z1),local_weyl(d2,x2,z2))


def apply_correction_branch(y,noise=0.,rng=None):
    """Unconditioned edge S on a post-measurement branch (486×r)."""
    z=y
    for _,pair,gate in correction_factors():
        z=pair_gate(z,pair,gate)
        if noise and rng.random()<noise:
            z=pair_gate(z,pair,sample_two_register_weyl(pair,rng))
    return z


def gauss_probability(y):
    """Pair of vertex projector acceptance values on pointer-entangled state."""
    x=np.asarray(y,complex)
    assert x.shape[0]==3
    norm=np.vdot(x,x).real
    return {v:float(sum(np.vdot(projector(z,v),projector(z,v)).real
                         for z in x)/norm)
            for v in ('A','B')}


def circuit_gauge_diagnostics(K):
    outcome,history=apply_syndrome(K,history=True)
    target=pointer_circuit(K)
    assert np.linalg.norm(outcome-target)<1e-10
    arr=[]
    pre=prepared_pointer(K)
    arr.append({'after':'pointer_preparation','Gauss_accept':gauss_probability(pre)})
    for i,(name,pair,z) in enumerate(history):
        arr.append({'after':f'{i+1:02d}_{name}',
                    'registers':[REGISTERS[a] for a in pair],
                    'Gauss_accept':gauss_probability(z)})
    arr.append({'after':'pointer_fourier',
                'Gauss_accept':gauss_probability(outcome)})
    endpoint=arr[-1]['Gauss_accept']
    assert all(abs(p-1)<1e-10 for p in endpoint.values())
    minimum={v:min(r['Gauss_accept'][v] for r in arr) for v in ('A','B')}
    return {'trace':arr,'minimum_intermediate_Gauss_accept':minimum,
            'final_Gauss_accept':endpoint,
            'exact_ref_error_fro':float(np.linalg.norm(outcome-target))}


def channel_fidelity(K,known_error_probability=.37,fault_probability=0.,
                     shots=48,seed=20261019):
    """Sample physical two-register Weyl errors after each compiled factor.

    Full logical entanglement fidelity: sum_{error sources, pointer labels}
    |Tr(K† E_branch)|² / d², then average stochastic gate faults.
    Correction m=2 is decomposed into C,SWAP,C† and also faulted.
    """
    if not (0<=known_error_probability<=1 and 0<=fault_probability<1 and shots>0):
        raise ValueError("probabilities/samples")
    assert K.shape==(DIM,RANK)
    F=dressed_exchange(K)
    rng=np.random.default_rng(seed)
    samples=[]
    for _ in range(shots):
        total=0.
        for E in (np.sqrt(1-known_error_probability)*K,
                  np.sqrt(known_error_probability)*F):
            branches=apply_syndrome(E,noise=fault_probability,rng=rng)
            for label,b in enumerate(branches):
                recovered=(apply_correction_branch(b,fault_probability,rng)
                           if label==2 else b)
                reduced=K.conj().T@recovered
                total+=abs(np.trace(reduced))**2/RANK**2
        assert -1e-8<=total<=1+1e-8,total
        samples.append(float(total))
    return {'source_dressed_swap_probability':known_error_probability,
            'per_two_register_factor_Weyl_fault_probability':fault_probability,
            'shots':shots,'seed':seed,
            'mean_recovered_entanglement_fidelity':float(np.mean(samples)),
            'standard_error':(float(np.std(samples,ddof=1)/np.sqrt(shots))
                              if shots>1 else None),
            'individual_shot_fidelity':samples}


def run():
    W=physical_isometry()
    K,F,U=code_isometry(W)
    rng=np.random.default_rng(31)
    random_x=rng.normal(size=(DIM,2))+1j*rng.normal(size=(DIM,2))
    err=float(np.linalg.norm(apply_syndrome(random_x)-pointer_circuit(random_x)))
    assert err<1e-10,err
    # Check factorized unconditioned S against source dressed exchange.
    correrr=float(np.linalg.norm(apply_correction_branch(random_x)
                                -dressed_exchange(random_x)))
    assert correrr<1e-10
    diagnostics=circuit_gauge_diagnostics(K)
    runs=[channel_fidelity(K,fault_probability=p,shots=48,seed=3100+i)
          for i,p in enumerate((0.,.02,.08))]
    assert abs(runs[0]['mean_recovered_entanglement_fidelity']-1)<1e-10
    return {'system':'finite S3 two-vertex gauge link, 486D, 14D Gauss sector, 4D phase logical code',
            'two_register_pair_ops_syndrome':len(measured_syndrome_factors()),
            'two_register_pair_ops_feedback_if_error':len(correction_factors()),
            'pointer_single_register_operations':2,
            'uses_arbitrary_two_register_18x18_and_9x9_gates':True,
            'native_two_qubit_or_two_qutrit_gate_count':'NOT derived',
            'routed_depth_and_gate_durations':'NOT derived',
            'controlled_swap_basis_Q_unitary_error':float(np.linalg.norm(
                exchange_eigenbasis_q().conj().T@exchange_eigenbasis_q()-np.eye(9))),
            'full_pointer_unitary_ref_error_fro':err,
            'full_feedback_ref_error_fro':correrr,
            'gauge_checks':diagnostics,
            'noise_scenarios':runs,
            'limitations':[
                '18x18 link–matter and 9x9 arbitrary two-qutrit transforms are two-register but NOT native low-cost gates',
                'Intermediate operations are not necessarily gauge covariant; gauge probabilities traced after each gate',
                'No gauge-penalty Hamiltonian or physical time/noise during gates',
                'Weyl pair faults sampled stochastically, not calibrated transmon or SU3 noise',
                'Only known dressed swap in 4D S3 phase code, NOT nine-quark 65-syndrome recovery',
                'No physical coupling graph/routing, transduction or pulse synthesis']}


def main():
    data=run()
    Path('colour_s3_two_register_syndrome_results.json').write_text(
        json.dumps(data,indent=2,sort_keys=True)+'\n')
    print('S3_TWO_REGISTER_SYNDROME',json.dumps({
        k:v for k,v in data.items() if k not in ('noise_scenarios','gauge_checks')
    },sort_keys=True),flush=True)
    for r in data['noise_scenarios']:
        print('S3_TWO_REGISTER_NOISE',json.dumps({k:v for k,v in r.items()
              if k!='individual_shot_fidelity'},sort_keys=True),flush=True)
    print('S3_TWO_REGISTER_INTERMEDIATE_GAUSS',json.dumps(
        data['gauge_checks']['minimum_intermediate_Gauss_accept'],sort_keys=True),flush=True)

if __name__=='__main__':
    main()
