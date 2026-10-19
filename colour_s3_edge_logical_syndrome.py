"""Edge-supported logical syndrome for finite S3 subgroup gauge simulator.

A=onsite SWAP(A0,A1), S=link-dressed SWAP(A0,B0) on the prior 486D
two-vertex, four-matter-qutrit model. Both commute with local Gauss
transformations; they obey A^2=S^2=(AS)^3=I. On the physical 14D
space U=AS has eigenphase multiplicities [6,4,4] for [1,w,w^2].
The +w eigenphase is a new 4D code K. S flips w <-> w^2.
A qutrit phase-estimation pointer for U resolves {clean,S error}
using four high-level POINTER-CONTROLLED operations: two controlled
onsite A and two controlled dressed-edge S, plus pointer prep/DFT.
No native two-body gate compilation or hardware noise is implied.
This differs from the previous arbitrary-eigenvector-pairing 4D code
and does NOT implement the nine-quark 21D code/65-syndrome decoder.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from colour_gauge_link_s3 import (
    DIM,G,MATTER,transform,physical_isometry,dressed_exchange,
    swap_matter,projector,
)

W3=np.exp(2j*np.pi/3)
SWAP_A=swap_matter(0,1)
DLOG=4

def onsite_swap(y):
    x=np.asarray(y,complex)
    vector=x.ndim==1
    if vector:x=x[:,None]
    z=np.einsum('ij,gjr->gir',SWAP_A,x.reshape(G,MATTER,-1))
    z=z.reshape(DIM,-1)
    return z[:,0] if vector else z

def edge_cycle(y):
    """U=A S, first dressed edge S then onsite swap A."""
    return onsite_swap(dressed_exchange(y))

def projectors(y):
    """Mathematical Fourier spectral projectors of the 3-cycle."""
    z1=edge_cycle(y);z2=edge_cycle(z1)
    return np.stack([(np.asarray(y)+W3**(-m)*z1+W3**(-2*m)*z2)/3
                     for m in range(3)])

def pointer_circuit(y):
    """Qutrit coherent Fourier phase estimation with edge-local supports.

    Pointer prepared (|0>+|1>+|2>)/sqrt3.
    Controlled S then A on pointer labels 1,2; repeat on label 2;
    finally inverse DFT. The controlled-S high-level primitive acts on
    pointer, A0, link and B0, and is not decomposed into 2-body gates.
    """
    x=np.asarray(y,complex)
    vector=x.ndim==1
    if vector:x=x[:,None]
    assert x.shape[0]==DIM
    branch=np.stack([x,x,x])/np.sqrt(3)
    for k in (1,2):
        branch[k]=dressed_exchange(branch[k])
        branch[k]=onsite_swap(branch[k])
    branch[2]=dressed_exchange(branch[2])
    branch[2]=onsite_swap(branch[2])
    z=np.fft.fft(branch,axis=0)/np.sqrt(3)
    return z[:,:,0] if vector else z

def code_isometry(W):
    """Entire w-eigenspace is a 4D gauge-invariant code."""
    v=edge_cycle(W)
    u=W.conj().T@v
    assert np.linalg.norm(edge_cycle(edge_cycle(v))-W)<1e-10
    P=(np.eye(14)+W3**(-1)*u+W3**(-2)*u@u)/3
    assert np.linalg.norm(P-P.conj().T)<1e-10
    ev,V=np.linalg.eigh(P)
    K=W@V[:,ev>.5]
    F=dressed_exchange(K)
    assert K.shape==(DIM,DLOG)
    assert np.linalg.norm(K.conj().T@K-np.eye(DLOG))<1e-10
    assert np.linalg.norm(edge_cycle(K)-W3*K)<1e-10
    assert np.linalg.norm(edge_cycle(F)-W3**2*F)<1e-10
    assert np.linalg.norm(K.conj().T@F)<1e-10
    return K,F,u

def recover_channel(K,p=.37,assignment=0.):
    """Explicit logical Kraus calculation for E={sqrt(1-p)I,sqrt(p)S}.

    m=1 needs no correction; m=2 gets gauge-covariant edge S feedback.
    Misreported pointer outcomes are REPLACED by I/4; this simple
    CPTP readout-error scenario is not a microscopic noisy device.
    """
    if not (0<=p<=1 and 0<=assignment<1):raise ValueError('probabilities')
    error=(np.sqrt(1-p)*K,np.sqrt(p)*dressed_exchange(K))
    fe=0.;sums=np.zeros(3)
    for E in error:
        branches=pointer_circuit(E)
        for m,b in enumerate(branches):
            sums[m]+=np.vdot(b,b).real/DLOG
            z=dressed_exchange(b) if m==2 else b
            logical=K.conj().T@z
            fe+=abs(np.trace(logical))**2/DLOG**2
    assert np.linalg.norm(sums-np.array([0,1-p,p]))<1e-10
    assert abs(fe-1)<1e-10
    return {'dressed_swap_probability':p,'correctly_reported_F_e':float(fe),
            'true_pointer_probabilities_0_1_2':sums.tolist(),
            'assumed_symmetric_assignment_error':assignment,
            'flag_and_replace_F_e':float((1-assignment)*fe+assignment/16)}

def run():
    W=physical_isometry()
    K,F,U=code_isometry(W)
    vals=np.linalg.eigvals(U)
    multiplicities=[int(sum(abs(vals-W3**m)<1e-8)) for m in range(3)]
    assert multiplicities==[6,4,4]
    assert np.linalg.norm(edge_cycle(edge_cycle(edge_cycle(W)))-W)<1e-10
    covariance={}
    for name,op in [('onsite_A_swap',onsite_swap),
                    ('dressed_edge_swap',dressed_exchange),
                    ('edge_cycle_U',edge_cycle)]:
        worst=0.
        for v in ('A','B'):
            for h in range(G):
                worst=max(worst,float(np.linalg.norm(
                    transform(op(W),h,v)-op(transform(W,h,v)))))
        assert worst<1e-10
        covariance[name]=worst
    rng=np.random.default_rng(37)
    x=rng.normal(size=(DIM,2))+1j*rng.normal(size=(DIM,2))
    actual=pointer_circuit(x)
    expect=projectors(x)
    difference=float(np.linalg.norm(actual-expect))
    norm_error=float(abs(np.vdot(actual,actual)-np.vdot(x,x)))
    assert difference<1e-10 and norm_error<1e-10
    for z in (K,F):
        for b in pointer_circuit(z):
            for v in ('A','B'):
                assert np.linalg.norm(projector(b,v)-b)<1e-10
    good=np.array([np.linalg.norm(z)**2/DLOG for z in pointer_circuit(K)])
    bad=np.array([np.linalg.norm(z)**2/DLOG for z in pointer_circuit(F)])
    assert np.linalg.norm(good-[0,1,0])<1e-10
    assert np.linalg.norm(bad-[0,0,1])<1e-10
    rec=recover_channel(K,.37,.01)
    return {'model':'486D S3 subgroup finite gauge-link model',
       'physical_gauss_dimension':14,'logical_code_dimension':DLOG,
       'edge_cycle_order':3,'eigenphases_1_w_w2_multiplicities':multiplicities,
       'gauge_covariance_residuals':covariance,
       'qutrit_pointer_Kraus_projector_residual':difference,
       'pointer_norm_preservation_residual':norm_error,
       'clean_pointer_0_1_2':good.tolist(),
       'dressed_swap_pointer_0_1_2':bad.tolist(),
       'known_swap_recovery':rec,
       'high_level_cost':{
           'pointer_qutrits':1,'uniform_preparations':1,'qutrit_inverse_Fourier':1,
           'pointer_controlled_dressed_edge_swaps':2,
           'pointer_controlled_onsite_swaps':2,
           'conditional_recovery_dressed_edge_swaps_max':1,
           'four_register_pointer_link_edge_primitive_native_compilation':'missing',
           'pointer_movement_or_connectivity_cost':'not estimated'},
       'main_result':'Edge/vertex-supported symmetry-preserving controlled primitives resolve LOGICAL swap syndrome, unlike Gauss-only measurements.',
       'scope':'New four-dimensional phase code, NOT the prior arbitrary-pairing code or 21D nine-quark SU3 code. High-level edge-local gates are not native two-body gates.'}

def main():
    data=run()
    Path('colour_s3_edge_logical_syndrome_results.json').write_text(
        json.dumps(data,indent=2,sort_keys=True)+'\n')
    print('S3_EDGE_LOGICAL_SYNDROME',json.dumps(data,sort_keys=True))

if __name__=='__main__':main()
