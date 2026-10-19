"""Extend finite S3 gauge-admissible noise algebra beyond matter swaps.

Exact *representation* facts on the 14D Gauss sector:
  A=onsite SWAP(A0,A1), B=onsite SWAP(B0,B1), and
  S=link-dressed SWAP(A0,B0) generate the complete 24-element
  matter-permutation group S4. Its common fixed (+1) subspace has
  dimension FOUR, furnishing an exact passive ququart against ALL
  six allowed matter transpositions and all their words.

A finite-group link electric-energy term is provided by the CENTRAL
standard-irrep Peter-Weyl projector P_std of the SIX-state link.
It commutes with independent local vertex Gauss transformations and
is supported on the link alone. H_E=P_std, and the electric phase
fault R=exp(-i pi H_E)=I-2P_std is a gauge-invariant unitary.

Adding H_E to the S4 exchange algebra breaks the four-dimensional
fixed subspace. Numerical commutant dimension on the entire 14D
Gauss Hilbert space is SIX with a clear singular-value gap; because
the noiseless-subsystem logical algebra M4(C) requires at least
dimension 16 in the commutant, no four-dimensional noiseless
subsystem exists for this SPECIFIC broadened error *-algebra.
This last classification depends on DOUBLE-PRECISION matrix rank,
not an exact computer-certified rational proof.

The existing active four-dimensional phase code does NOT satisfy
Knill-Laflamme for both dressed swap and this link electric error,
and its existing ideal syndrome does NOT correct the electric error
perfectly. Report bounded comparisons without asserting a QEC win.

This is the finite S3 toy gauge model, not continuous SU3 or the
original nine-quark 21D/65-syndrome architecture.
"""
from __future__ import annotations

from pathlib import Path
import json

import numpy as np

from colour_gauge_link_s3 import (
    DIM,G,MATTER, physical_isometry,dressed_exchange,
    swap_matter,transform,
)
from colour_s3_edge_logical_syndrome import (
    code_isometry,onsite_swap,pointer_circuit,
)
from colour_s3_strict_gauss_pair_nogo import (
    apply_link,link_central_projectors,
)

DLOG=4
SOURCE_P=.37


def physical_onsite_swap(y,i,j):
    """Unconditional local matter-pair permutation on A or B vertex."""
    if (i,j) not in ((0,1),(2,3)):
        raise ValueError("Expected a single-vertex matter-pair swap")
    x=np.asarray(y,complex)
    vec=x.ndim==1
    if vec:x=x[:,None]
    q=swap_matter(i,j)
    out=np.einsum("ab,gbr->gar",q,x.reshape((G,MATTER,-1)))
    result=out.reshape(x.shape)
    return result[:,0] if vec else result


def compress(W,op):
    return W.conj().T@op(W)


def signed_exchange_generators(W):
    return {
       "onsite_A_01":compress(W,lambda x:physical_onsite_swap(x,0,1)),
       "onsite_B_23":compress(W,lambda x:physical_onsite_swap(x,2,3)),
       "dressed_edge_A0_B0":compress(W,dressed_exchange),
    }


def permutation_group_closure(generators,tolerance=1e-8):
    """Enumerate all 24 group elements, not only short products."""
    found=[np.eye(14,dtype=complex)]
    idx=0
    while idx<len(found):
        x=found[idx]
        for op in generators:
            y=op@x
            if not any(np.linalg.norm(y-z)<tolerance for z in found):
                found.append(y)
                if len(found)>24:
                    raise AssertionError("Expected complete S4, not larger group")
        idx+=1
    assert len(found)==24,len(found)
    return found


def s4_fixed_projector(W):
    gen=signed_exchange_generators(W)
    for name,S in gen.items():
        assert np.linalg.norm(S@S-np.eye(14))<2e-10,name
        assert np.linalg.norm(S-S.conj().T)<2e-10,name
    group=permutation_group_closure(tuple(gen.values()))
    P=sum(group)/len(group)
    assert np.linalg.norm(P@P-P)<2e-10
    assert np.linalg.norm(P-P.conj().T)<2e-10
    rank=int(round(np.trace(P).real))
    assert rank==4,rank
    return P,gen,group


def deterministic_basis(P,dimension):
    basis=[]
    for j in range(P.shape[0]):
        v=P[:,j].astype(complex).copy()
        for _ in range(2):
            for q in basis:
                v-=q*np.vdot(q,v)
        norm=np.linalg.norm(v)
        if norm>1e-9:
            basis.append(v/norm)
        if len(basis)==dimension:break
    if len(basis)!=dimension:
        raise AssertionError("Insufficient rank")
    return np.stack(basis,axis=1)


def electric_projector(W):
    """Link-only S3 standard-irrep electric Hamiltonian in physical 14D."""
    plink=link_central_projectors()["standard"]
    assert plink.shape==(G,G)
    E=compress(W,lambda x:apply_link(x,plink))
    assert np.linalg.norm(E-E.conj().T)<1e-10
    assert np.linalg.norm(E@E-E)<1e-10
    assert abs(np.trace(E).real-9)<1e-10
    return E


def electric_fault(y):
    """Full 486D single-link electric π phase R=I-2 P_standard."""
    plink=link_central_projectors()["standard"]
    return np.asarray(y)-2*apply_link(y,plink)


def commutant_dimension(ops,threshold=1e-8):
    """Numerical linear algebra on M_14(C), explicitly NOT exact proof.

    Every commuting X solves (O^T tensor I - I tensor O) vec(X)=0.
    Either row/column vectorization changes basis but not rank.
    If numerical commutant dimension is below 16, no M4 noiseless
    logical-subsystem factor can be contained in the commutant.
    """
    n=ops[0].shape[0]
    assert n==14 and all(op.shape==(n,n) for op in ops)
    I=np.eye(n)
    constraints=np.concatenate([
        np.kron(op.T,I)-np.kron(I,op)
        for op in ops],axis=0)
    singular=np.linalg.svd(constraints,compute_uv=False)
    rank=int(np.sum(singular>threshold))
    dimension=n*n-rank
    nonzero=singular[singular>threshold]
    zero=singular[singular<=threshold]
    return {
        "commutant_complex_dimension":dimension,
        "numerical_rank":rank,
        "smallest_nonzero_singular_value":float(nonzero[-1]),
        "largest_null_singular_value":float(zero[0]) if len(zero) else None,
        "rank_threshold":threshold,
        "no_exact_rational_rank_certificate":True,
        "no_fitted_four_dim_noiseless_subsystem":bool(dimension<16),
    }


def code_kl_diagnostics(K,E,swap):
    """KL tests for joint unknown channel {I,S,R} on 4D active code."""
    # K and physical operators already live in 14D via K_14.
    R=np.eye(14)-2*E
    ops=[np.eye(14),swap,R]
    compressed={}
    worst=0.
    for i,Oi in enumerate(ops):
        for j,Oj in enumerate(ops):
            c=K.conj().T@Oi.conj().T@Oj@K
            scalar=np.trace(c)/DLOG
            deviation=float(np.linalg.norm(c-scalar*np.eye(DLOG)))
            compressed[f"{i}{j}"]={
                "best_scalar_re":float(scalar.real),
                "best_scalar_im":float(scalar.imag),
                "non_scalar_Fro":deviation}
            worst=max(worst,deviation)
    return {"fault_set":["I","dressed_swap_A0B0","electric_phase_I_minus_2Pstd"],
            "all_Knill_Laflamme_scalar_conditions":bool(worst<1e-10),
            "max_non_scalar_KL_compression_Fro":worst,
            "compressions":compressed}


def entanglement_fidelity_known_error(code14,operator,p=SOURCE_P):
    """No-feedback F_e of error channel {(sqrt(1-p)I),(sqrt(p)U)}."""
    if not 0<=p<=1:raise ValueError(p)
    m=code14.conj().T@operator@code14
    return float((1-p)+p*abs(np.trace(m)/DLOG)**2)


def existing_ideal_phase_syndrome_electric_response(K486,p=SOURCE_P):
    """Original non-microscopic exact PR48 pointer instrument.

    A valid COMPLETE pointer KRAUS decomposition is used; not assuming
    that knowing the Pauli-like dressed error also corrects an electric
    flux-phase fault. Ideal feedback applies S for reported label 2.
    """
    if not 0<=p<=1:raise ValueError(p)
    fe=0.
    pointer=np.zeros(3)
    for E in (np.sqrt(1-p)*K486,np.sqrt(p)*electric_fault(K486)):
        branches=pointer_circuit(E)
        for m,branch in enumerate(branches):
            pointer[m]+=np.vdot(branch,branch).real/DLOG
            corrected=dressed_exchange(branch) if m==2 else branch
            fe+=abs(np.trace(K486.conj().T@corrected))**2/DLOG**2
    assert abs(sum(pointer)-1)<2e-10
    return {"ideal_existing_K_phase_pointer_F_e_for_electric_channel":float(fe),
            "pointer_outcome_probabilities_0_1_2":pointer.tolist(),
            "feedback":"Known dressed swap S only on pointer label 2, not optimized for electric error",
            "no_matched_microscopic_decoder_for_electric_fault":True}


def study():
    W=physical_isometry()
    P,gen,group=s4_fixed_projector(W)
    C14=deterministic_basis(P,DLOG)
    C=W@C14
    for g in group:
        assert np.linalg.norm(g@C14-C14)<2e-10
    for h in range(G):
        for v in ("A","B"):
            assert np.linalg.norm(transform(C,h,v)-C)<2e-10
    E=electric_projector(W)
    R=np.eye(14)-2*E
    for h in range(G):
        for v in ("A","B"):
            test=np.random.default_rng(h+13).normal(size=(DIM,2))
            test=test+1j*np.random.default_rng(h+29).normal(size=test.shape)
            assert np.linalg.norm(
                transform(electric_fault(test),h,v)
                -electric_fault(transform(test,h,v)))<2e-10
    # S4-fixed passive ququart is exactly 4D. E is not scalar there.
    compressed_E=C14.conj().T@E@C14
    eigen=np.linalg.eigvalsh(compressed_E)
    assert np.linalg.norm(eigen-np.array([0.,1/3,2/3,1.]))<2e-10
    leak=float(np.linalg.norm((np.eye(14)-P)@E@C14))
    assert leak>.1
    dims_exchange=commutant_dimension(list(gen.values()))
    dims_electric=commutant_dimension([*gen.values(),E])
    assert dims_exchange["commutant_complex_dimension"]==24
    assert dims_electric["commutant_complex_dimension"]==6
    assert dims_electric["smallest_nonzero_singular_value"]>.1
    assert dims_electric["no_fitted_four_dim_noiseless_subsystem"]
    K486,F,_=code_isometry(W)
    K=W.conj().T@K486
    KL=code_kl_diagnostics(K,E,gen["dressed_edge_A0_B0"])
    assert not KL["all_Knill_Laflamme_scalar_conditions"]
    passive_fe=entanglement_fidelity_known_error(C14,R)
    active_raw=entanglement_fidelity_known_error(K,R)
    ideal=existing_ideal_phase_syndrome_electric_response(K486)
    # Supplementary exact rank observation: P_std has rank nine in
    # 14D, and its complement rank five. A DFS for JUST this one
    # electric operator can hold 4D, but not jointly with all swaps.
    return {
        "model":"finite S3 two-vertex non-Abelian gauge-link 486D / Gauss 14D",
        "gauge_group":"signed-permutation S3 subset SU3",
        "matter_swaps":"onsite A01, onsite B23, dressed A0B0 generate all S4",
        "generated_matter_swap_group_order":len(group),
        "full_swap_group_fixed_4D_subspace_rank":int(round(np.trace(P).real)),
        "permutation_group_full_subspace_is_decoherence_free":True,
        "electric_link_generator":"P_standard = (d/|S3|) sum_h chi_standard(h^-1) L_h; one-link central projector",
        "electric_link_physical_rank":int(round(np.trace(E).real)),
        "electric_link_nonstandard_rank":14-int(round(np.trace(E).real)),
        "electric_phase_fault":"R=e^{-i pi P_standard}=I-2P_standard",
        "electric_commutes_with_both_local_Gauss_groups":True,
        "electric_phase_nontrivial_in_S4_fixed_4D_sector":True,
        "compressed_electric_eigenvalues_on_passive_S4_fixed_4D":eigen.tolist(),
        "electric_leakage_out_of_S4_fixed_sector_Fro":leak,
        "swap_only_commutant":dims_exchange,
        "swap_plus_electric_commutant":dims_electric,
        "no_passive_four_logical_dim_noiseless_subsystem_numerical":True,
        "KL_existing_active_K_three_error_set":KL,
        "electric_channel_probability":SOURCE_P,
        "passive_S4_fixed_no_feedback_F_e_under_electric_channel":passive_fe,
        "active_K_no_feedback_F_e_under_electric_channel":active_raw,
        "existing_active_ideal_phase_syndrome_on_electric_channel":ideal,
        "single_electric_fault_alone_admits_4D_decoherence_free_subspace":True,
        "important_bounds":[
            "The S4 exchange group contains ALL matter transpositions in this four-qutrit model and STILL leaves exactly a 4D passive singlet code",
            "Adding the local gauge-invariant finite-S3 electric LINK irrep projector breaks that passive fixed subspace and reduces the 14D commutant to numerical dimension six",
            "No 4D passive noiseless subsystem under the generated *-algebra follows at the numerical rank resolution; not an exact rational certificate",
            "The existing active 4D K does NOT satisfy exact KL for error set I,S,R: not a constructed active correction for this broadened physical error algebra",
            "The ideal active pointer response to R is not microscopic 495D feedback, and no finite-U/full-noise end-to-end electric fault performance is demonstrated",
            "Physical electric projector is a meaningful finite-group lattice electric Hamiltonian term, but chosen pi pulse strength and fault probability are synthetic",
            "No code state preparation, link energy pulse/charged-transport synthesis, physical T1/Tphi calibration, device-specific leakage or matched hardware costs",
            "This is finite S3 only, not full continuous SU3 lattice, Standard Model quarks, or nine-quark 21D/65-syndrome physical recovery",
        ]}


def main():
    result=study()
    Path("colour_s3_electric_error_algebra_results.json").write_text(
        json.dumps(result,indent=2,sort_keys=True)+"\n")
    print("S4_PASSIVE_FIXED_AND_ELECTRIC_COMMUTANT",json.dumps({
        "full_S4_order":result["generated_matter_swap_group_order"],
        "fixed_rank":result["full_swap_group_fixed_4D_subspace_rank"],
        "electric_rank":result["electric_link_physical_rank"],
        "electric_compression_spectrum":result["compressed_electric_eigenvalues_on_passive_S4_fixed_4D"],
        "algebra":result["swap_plus_electric_commutant"],
        "joint_KL":result["KL_existing_active_K_three_error_set"]["all_Knill_Laflamme_scalar_conditions"],
        "passive_electric_F_e":result["passive_S4_fixed_no_feedback_F_e_under_electric_channel"],
        "ideal_active_electric_F_e":result["existing_active_ideal_phase_syndrome_on_electric_channel"]["ideal_existing_K_phase_pointer_F_e_for_electric_channel"],
    },sort_keys=True))


if __name__=="__main__":
    main()
