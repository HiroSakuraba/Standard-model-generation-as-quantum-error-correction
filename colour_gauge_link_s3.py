"""Minimal non-Abelian gauge-link exchange and vertex syndrome experiment.

G=S3 embeds in SU(3) as determinant-one signed permutation matrices.
Four matter qutrits (two per site), one six-level oriented link:
3^4*6=486 dimensions and a 14D physical Gauss-invariant subspace.

Dressed across-link colour swap is gauge covariant. Bare cross-site
swap violates independent vertex Gauss constraints. Local gauge checks
are implemented with an explicit six-state coherent pointer circuit:
uniform preparation, controlled vertex gauge action, inverse Fourier,
readout. They cannot diagnose dressed gauge-invariant swap errors.

A four-dimensional gauge-invariant *toy* code corrects a KNOWN dressed
swap by a separate abstract NONLOCAL logical syndrome measurement. This
is NOT the nine-quark 21D code, full SU(3) lattice, or calibrated gate set.
"""
from __future__ import annotations

from itertools import permutations
from pathlib import Path
import json
import numpy as np

PERMS=tuple(permutations(range(3)))
G=6
MATTER=3**4
DIM=G*MATTER
ID=PERMS.index((0,1,2))


def compose(a:int,b:int)->int:
    p,q=PERMS[a],PERMS[b]
    return PERMS.index(tuple(p[q[i]] for i in range(3)))


def inverse(a:int)->int:
    return next(b for b in range(G) if compose(a,b)==ID and compose(b,a)==ID)


def signed_matrix(a:int)->np.ndarray:
    p=PERMS[a]
    parity=(-1)**sum(p[i]>p[j] for i in range(3) for j in range(i+1,3))
    M=np.zeros((3,3),complex)
    for i,j in enumerate(p):M[j,i]=parity
    assert abs(np.linalg.det(M)-1)<1e-12
    return M


REP=tuple(signed_matrix(i) for i in range(G))


def four_rep(a:int,axes:tuple[int,...])->np.ndarray:
    factors=[REP[a] if k in axes else np.eye(3) for k in range(4)]
    return np.kron(np.kron(factors[0],factors[1]),np.kron(factors[2],factors[3]))


def swap_matter(i:int,j:int)->np.ndarray:
    order=np.arange(MATTER).reshape((3,3,3,3))
    axes=list(range(4));axes[i],axes[j]=axes[j],axes[i]
    return np.eye(MATTER,dtype=complex)[np.transpose(order,axes).ravel()]


def transform(y:np.ndarray,h:int,vertex:str)->np.ndarray:
    """Local gauge action on two matter qutrits plus incident link."""
    assert vertex in ('A','B')
    yy=np.asarray(y,complex)
    was_vector=yy.ndim==1
    if was_vector:yy=yy[:,None]
    assert yy.shape[0]==DIM
    out=np.zeros_like(yy)
    local=four_rep(h,(0,1) if vertex=='A' else (2,3))
    for g in range(G):
        new=compose(h,g) if vertex=='A' else compose(g,inverse(h))
        out[new*MATTER:(new+1)*MATTER,:]=local@yy[g*MATTER:(g+1)*MATTER,:]
    return out[:,0] if was_vector else out


def projector(y:np.ndarray,vertex:str)->np.ndarray:
    return sum((transform(y,h,vertex) for h in range(G)))/G


def physical_isometry()->np.ndarray:
    """Build all simultaneous site-Gauss invariants, rank 14."""
    P=sum((four_rep(h,(0,1,2,3)) for h in range(G)))/G
    assert np.linalg.norm(P@P-P)<1e-11
    ev,vec=np.linalg.eigh(P)
    V=vec[:,ev>.5]
    assert V.shape==(MATTER,14)
    W=np.zeros((DIM,14),complex)
    for g in range(G):
        W[g*MATTER:(g+1)*MATTER]=four_rep(inverse(g),(2,3))@V/np.sqrt(G)
    assert np.linalg.norm(W.conj().T@W-np.eye(14))<1e-11
    for v in ('A','B'):
        for h in range(G):assert np.linalg.norm(transform(W,h,v)-W)<1e-10
    return W


def dressed_exchange(y:np.ndarray,dressed:bool=True)->np.ndarray:
    arr=np.asarray(y,complex)
    vector=arr.ndim==1
    if vector:arr=arr[:,None]
    assert arr.shape[0]==DIM
    out=np.zeros_like(arr)
    bare=swap_matter(0,2)
    for g in range(G):
        q=bare@arr[g*MATTER:(g+1)*MATTER]
        if dressed:
            q=four_rep(g,(0,))@four_rep(inverse(g),(2,))@q
        out[g*MATTER:(g+1)*MATTER]=q
    return out[:,0] if vector else out


def dft_ancilla_kraus(y:np.ndarray,vertex:str)->np.ndarray:
    """Coherent six-state pointer: uniform -> controlled Gauss -> inverse DFT.

    Pointer output k=0 is exactly P_Gauss; the other five outcomes
    make the instrument trace preserving. The controlled gauge operation
    is supported on the vertex matter and its incident link, but it is
    not decomposed into a calibrated two-body gate network.
    """
    x=np.asarray(y,complex)
    if x.ndim==1:x=x[:,None]
    assert x.shape[0]==DIM
    controlled=np.stack([transform(x,h,vertex) for h in range(G)],axis=0)/np.sqrt(G)
    return np.fft.ifft(controlled,axis=0)*np.sqrt(G)


def local_gauss_probabilities(y:np.ndarray)->dict:
    x=np.asarray(y,complex)
    denom=np.vdot(x,x).real
    assert denom>0
    return {v:float(np.vdot(projector(x,v),projector(x,v)).real/denom)
            for v in ('A','B')}


def logical_swap_code(W:np.ndarray):
    """Separate gauge-invariant four-dimensional code correcting one known swap."""
    physical_swap=W.conj().T@dressed_exchange(W)
    assert np.linalg.norm(physical_swap-physical_swap.conj().T)<1e-10
    assert np.linalg.norm(physical_swap@physical_swap-np.eye(14))<1e-9
    eigen,V=np.linalg.eigh(physical_swap)
    pos=V[:,eigen>.5];neg=V[:,eigen<-.5]
    k=min(pos.shape[1],neg.shape[1])
    assert k>=1
    code=W@(pos[:,:k]+neg[:,:k])/np.sqrt(2)
    errored=dressed_exchange(code)
    assert np.linalg.norm(code.conj().T@code-np.eye(k))<1e-9
    assert np.linalg.norm(code.conj().T@errored)<1e-9
    Q0=code@code.conj().T
    Q1=errored@errored.conj().T
    assert np.linalg.norm(Q0@code-code)<1e-10
    assert np.linalg.norm(Q1@errored-errored)<1e-10
    if k==7: assert np.linalg.norm(W@W.conj().T-Q0-Q1)<1e-9
    return code,errored,{'swap_positive_multiplicity':pos.shape[1],
                         'swap_negative_multiplicity':neg.shape[1],
                         'logical_code_dimension':k,
                         'logical_syndrome_labels':3,
                         'logical_syndrome_physical_gate_network':'not compiled'}


def study()->dict:
    W=physical_isometry()
    assert W.shape==(DIM,14)
    covariance={}
    for dressed,name in ((True,'dressed'),(False,'bare')):
        largest=0.
        for v in ('A','B'):
            for h in range(G):
                diff=(transform(dressed_exchange(W,dressed),h,v)
                      -dressed_exchange(transform(W,h,v),dressed))
                largest=max(largest,float(np.linalg.norm(diff)))
        covariance[name]=largest
    assert covariance['dressed']<1e-10
    assert covariance['bare']>.1
    dressedW=dressed_exchange(W)
    assert np.linalg.norm(W@W.conj().T@dressedW-dressedW)<1e-10
    ancilla={}
    for v in ('A','B'):
        out=dft_ancilla_kraus(W,v)
        probs=np.array([np.vdot(x,x).real/14 for x in out])
        assert abs(probs.sum()-1)<1e-10
        assert abs(probs[0]-1)<1e-10
        assert np.linalg.norm(out[0]-W)<1e-10
        assert max(probs[1:])<1e-10
        ancilla[v]={'gauss_pass_probability':float(probs[0]),
                    'other_probability':float(probs[1:].sum())}
    bareW=dressed_exchange(W,False)
    bare_rates={v:float(np.linalg.norm(projector(bareW,v))**2/14) for v in ('A','B')}
    dress_rates={v:float(np.linalg.norm(projector(dressedW,v))**2/14) for v in ('A','B')}
    assert min(dress_rates.values())>1-1e-10
    assert min(bare_rates.values())<1-1e-3
    assignment_p=.01
    pointer_accept={'physical_or_dressed':1-assignment_p,
       'bare_swap_A':bare_rates['A']*(1-assignment_p)+(1-bare_rates['A'])*assignment_p/5,
       'bare_swap_B':bare_rates['B']*(1-assignment_p)+(1-bare_rates['B'])*assignment_p/5}
    assert abs(pointer_accept['bare_swap_A']-(.25*.99+.75*.01/5))<1e-10
    code,swapped,stats=logical_swap_code(W)
    k=stats['logical_code_dimension']
    p=.37
    E=(np.sqrt(1-p)*code,np.sqrt(p)*swapped)
    # Ideal NONLOCAL logical syndrome {K, S K, complement}, with swap feedback.
    coherent=(abs(np.trace(code.conj().T@E[0]))**2
              +abs(np.trace(swapped.conj().T@E[1]))**2)/k**2
    assert abs(coherent-1)<1e-10
    return {'gauge_group':'S3 embedded in SU3 via determinant-one signed 3x3 permutations',
        'link_hilbert_dimension':G,'vertices':2,'matter_qutrits':4,
        'full_hilbert_dimension':DIM,'physical_Gauss_singlet_dimension':W.shape[1],
        'gauge_covariance_commutator_norms_on_physical_basis':covariance,
        'single_vertex_gauss_pointer_ancilla_dimension':6,
        'min_qutrits_for_this_six_label_pointer':2,
        'controlled_group_branches_per_vertex_in_this_circuit':6,
        'vertex_gauge_measurement':ancilla,
        'illustrative_symmetric_pointer_misassignment_probability':assignment_p,
        'illustrative_reported_gauss_accept_probabilities':pointer_accept,
        'gauge_pass_after_bare_swap':bare_rates,
        'gauge_pass_after_dressed_swap':dress_rates,
        'dressed_swap_is_gauss_undetectable':True,
        'small_physical_swap_code':stats,
        'known_swap_channel_ideal_nonlocal_decoder_F_e':float(coherent),
        'separation_of_measurements':('Local vertex Gauss circuits detect gauge-breaking faults but cannot detect '
            'gauge-invariant dressed swaps; full logical syndrome remains nonlocal.'),
        'scope':('Exact finite non-Abelian S3 subgroup test only: not nine-colour-qutrit 21D '
            'sign-twist code, continuous SU3 Wilson links, physical quarks, calibrated entangling gates, '
            'or fault-tolerant logical syndrome extraction.')}


def main():
    result=study()
    Path('colour_gauge_link_s3_results.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print('S3_GAUGE_LINK_AND_GAUSS_POINTER',json.dumps(result,sort_keys=True))


if __name__=='__main__':main()
