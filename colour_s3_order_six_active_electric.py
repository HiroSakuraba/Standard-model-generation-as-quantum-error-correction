"""Construct a four-logical-state ACTIVE code against S3-subgroup exchanges
and one gauge-invariant link-electric phase, with an exact noiseless-
subsystem obstruction and a verified ideal recovery.

Each of the proper order-six matter-swap groups <A,S> and <B,S> has
a SIX-dimensional trivial subrepresentation of the 14D Gauss sector.
The electric projector E restricted to that six-dimensional space has
the *exact* spectral multiset {0,1/3,2/3,2/3,1,1}. No four-dimensional
passive scalar common eigenspace exists (the largest is dimension two).

A 4D code is instead selected from four orthogonal vectors with
C† E C = (2/3) I4, while every g in the selected order-six group
acts as g C=C. For the electric pi-phase R=I-2E:
  C†R C = -(1/3) I4,
  V = (R C + C/3)/(sqrt(8)/3),
  C†V=0, V†V=I4.

A syndrome with projectors C C† and V V†, followed by a partial
isometry V -> C, recovers any mixture of these six permutations and
one electric phase EXACTLY on the four-dimensional code. The error
branches are NOT necessarily themselves orthogonal: the component of
R in the original code has amplitude -1/3. Two outcome probabilities
for a pure R event are 1/9 and 8/9.

Unlike a construction on just one exchange and R, here there is
NO 4D noiseless subsystem for the joint algebra: over Q the subgroup
commutant has 52 integer orbit indicators; commuting with E imposes
an exact rank-41 rational matrix, leaving an 11D commutant (<16).
That is a genuine active-vs-passive distinction within this precise
7-operator error set. It DOES NOT prove a hardware advantage.

Does NOT necessarily correct sequential compound errors gR or Rg
as additional Kraus errors, arbitrary all-24 matter permutations,
continuous SU3 noise, or real quark/gluon field dynamics.

We also clarify the 'virtual particles' question: 495D Hubbard
simulations include high-energy doubly occupied fermion states and
their virtual superexchange; this exact finite-group calculation is
a fixed-particle, fixed-link Hilbert space with no virtual gluon/
quark-antiquark vacuum loops or QCD renormalization.
"""
from __future__ import annotations
from itertools import product
from pathlib import Path
import json
import numpy as np
import sympy as sp

from colour_s3_exact_electric_commutant import (
    DGAUSS, four_colour_orbits, exact_generators,
    physical_orbit_isometry,
)
from colour_s3_exact_active_kl_bound import (
    integer_permutation_rep, generator_actions,
)

D=14
LOGICAL=4
assert DGAUSS==D


def closure_matrices(generators):
    """Full exact integer group from selected involutory generators."""
    group=[sp.eye(D)]
    i=0
    while i<len(group):
        current=group[i]
        for generator in generators:
            next_op=generator*current
            if not any(next_op==g for g in group):
                group.append(next_op)
                assert len(group)<=6
        i+=1
    assert len(group)==6
    return group


def ordered_pair_centralizer_basis(group):
    """Exact commutant of the subgroup action on 14 colour-orbit labels.

    The group acts by permutation matrices, so every X commuting with
    the group is constant on simultaneous subgroup orbits of ordered
    pairs (i,j). Each orbit indicator is an integral independent basis.
    """
    actions=[]
    for M in group:
        p=[]
        for j in range(D):
            nonzero=[i for i in range(D) if M[i,j]!=0]
            assert len(nonzero)==1
            p.append(nonzero[0])
        actions.append(tuple(p))
    unseen=set(product(range(D),repeat=2))
    basis=[]
    while unseen:
        i,j=min(unseen)
        orb={(g[i],g[j]) for g in actions}
        assert orb<=unseen
        unseen-=orb
        X=sp.zeros(D)
        for a,b in orb:X[a,b]=1
        basis.append(X)
    assert len(basis)==52
    assert sum(basis,sp.zeros(D))==sp.ones(D,D)
    return basis


def exact_restricted_algebra_certificate(label,group,E):
    """Exact 196x52 rational E commutator rank 41, kernel 11."""
    basis=ordered_pair_centralizer_basis(group)
    constraint=sp.zeros(D*D,len(basis))
    for col,X in enumerate(basis):
        comm=X*E-E*X
        for i in range(D):
            for j in range(D):
                constraint[i*D+j,col]=comm[i,j]
    rank=constraint.rank()
    assert rank==41
    nullity=len(basis)-rank
    assert nullity==11
    return {
        "group":label,
        "gauge_physical_dimension":D,
        "exchange_group_order":len(group),
        "exchange_commutant_exact_orbit_basis_dimension":len(basis),
        "electric_constraints_shape":[D*D,len(basis)],
        "electric_constraint_exact_rational_rank":rank,
        "joint_error_commutant_exact_dimension":nullity,
        "no_four_logical_noiseless_subsystem":bool(nullity<LOGICAL**2),
        "rank_proven_using":"exact SymPy rational arithmetic, not numerical singular values",
    }


def exact_sector_spectrum(group,E):
    """Restriction to the 6D group-trivial sector: rational charpoly.

    The 14D characteristic polynomial includes 8 automatic zero
    roots of PEP outside the 6D fixed subspace. Inside it, E
    compressed has 0, 1/3, 2/3 twice, 1 twice.
    """
    P=sum(group,sp.zeros(D))/len(group)
    assert P*P==P and sp.trace(P)==6
    Ecomp=P*E*P
    z=sp.Symbol("z")
    char=sp.factor(Ecomp.charpoly(z).as_expr())
    expected=z**9*(3*z-1)*(3*z-2)**2*(z-1)**2/27
    assert sp.expand(char-expected)==0, (char, expected)
    return P, {
        "group_fixed_sector_dimension":6,
        "exact_compressed_electric_eigenvalues_and_multiplicities":{
            "0":1,"1/3":1,"2/3":2,"1":2},
        "characteristic_polynomial_14D":str(char),
    }


def gauge_fixed_normalized_matrix(M,orbit_sizes):
    """Convert rational orbit-SUM matrices into Hermitian/unitary basis."""
    sizes=np.sqrt(np.asarray(orbit_sizes,dtype=float))
    return np.array(M.tolist(),dtype=np.float64)*(
        sizes[:,None]/sizes[None,:])


def construct_active_code(label,group,E,orbit_sizes):
    """One rank-four code with C† E C=(2/3)I within 6D fixed sector.

    Full rank-four compression is obtained by two exact eigenvalue
    2/3 vectors, plus separate combinations of (0,1) and (1/3,1).
    This construction is deterministic up to eigenspace basis/signs;
    verification is at double precision, while multiplicities and
    existence follow from the exact rational characteristic polynomial.
    """
    P,spectrum=exact_sector_spectrum(group,E)
    phys_P=gauge_fixed_normalized_matrix(P,orbit_sizes)
    phys_E=gauge_fixed_normalized_matrix(E,orbit_sizes)
    assert np.linalg.norm(phys_P-phys_P.T)<2e-12
    assert np.linalg.norm(phys_E-phys_E.T)<2e-12
    vals,V=np.linalg.eigh(phys_P)
    fixed=V[:,vals>.5]
    assert fixed.shape==(D,6)
    evalues,Q=np.linalg.eigh(fixed.T@phys_E@fixed)
    target=np.array([0.,1/3,2/3,2/3,1.,1.])
    assert np.linalg.norm(evalues-target)<2e-12,(evalues,target)
    orth=fixed@Q
    e0,e13,e23a,e23b,e1a,e1b=[orth[:,i] for i in range(6)]
    C=np.column_stack((
        e23a, e23b,
        np.sqrt(1/3)*e0+np.sqrt(2/3)*e1a,
        np.sqrt(.5)*e13+np.sqrt(.5)*e1b
    ))
    R=np.eye(D)-2*phys_E
    # R maps code to itself with scalar overlap alpha=-1/3, but also
    # to a distinct orthogonal 4D syndrome space with amplitude sqrt8/3.
    overlap=-1/3
    beta=np.sqrt(8)/3
    Verror=(R@C-overlap*C)/beta
    residuals={
        "C_orthonormal":float(np.linalg.norm(C.T@C-np.eye(4))),
        "electric_compressed_scalar_two_thirds":float(np.linalg.norm(
            C.T@phys_E@C-(2/3)*np.eye(4))),
        "R_compressed_scalar_minus_one_third":float(np.linalg.norm(
            C.T@R@C-overlap*np.eye(4))),
        "electric_syndrome_orthogonality":float(np.linalg.norm(
            C.T@Verror)),
        "electric_syndrome_normalization":float(np.linalg.norm(
            Verror.T@Verror-np.eye(4))),
        "R_component_reconstruction":float(np.linalg.norm(
            R@C-overlap*C-beta*Verror)),
        "invariant_under_every_selected_exchange":max(
            float(np.linalg.norm(gauge_fixed_normalized_matrix(g,orbit_sizes)@C-C))
            for g in group),
    }
    assert max(residuals.values())<5e-12,residuals
    return C,Verror, {
        "group":label,
        **spectrum,
        "logical_dimension":LOGICAL,
        "electric_phase_error_code_overlap_amplitude":"-1/3",
        "electric_phase_error_orthogonal_syndrome_amplitude":"sqrt(8)/3",
        "ideal_R_fault_syndrome_probabilities":{"code":1/9,"electric":8/9},
        "numerical_Knill_Laflamme_residuals":residuals,
        "recovery_type":"ideal two-outcome 4D+4D subspace measurement followed by syndrome-dependent partial isometry; no native circuit",
    }


def verify_five_plus_one_kraus_recovery(C,Verror,group,E,sizes):
    """Independently evaluate a complete ideal CPTP channel on code states.

    For each of six independent group errors and R, measure P_C, P_V,
    recover P_V to C, and trace any orthogonal complement as a flagged
    erasure with a fixed logical replacement. Since the admitted
    channel leaves the C+V syndrome space, flagged probability=0.

    Branch Kraus on encoded states are matrices (C† F† E_i C):
       F0=P_C; F1=C V†
    and give a multiple of I4 for EVERY branch and every listed error.
    """
    P=C@C.T
    Q=Verror@Verror.T
    assert np.linalg.norm(P@Q)<5e-12
    phys_E=gauge_fixed_normalized_matrix(E,sizes)
    R=np.eye(D)-2*phys_E
    output=[]
    for name,Op in (
        [(f"exchange_{i}",gauge_fixed_normalized_matrix(g,sizes))
         for i,g in enumerate(group)]
        +[("electric_phase",R)]
    ):
        k0=C.T@P@Op@C
        k1=Verror.T@Q@Op@C
        recovered=sum(abs(np.trace(k)/LOGICAL)**2 for k in (k0,k1))
        trace_prob=sum(float(np.linalg.norm(k,"fro")**2/LOGICAL)
                       for k in (k0,k1))
        leaked=1-trace_prob
        nonscalar=max(float(np.linalg.norm(k-np.trace(k)/LOGICAL*np.eye(4)))
                      for k in (k0,k1))
        assert abs(recovered-1)<2e-11,(name,recovered)
        assert abs(trace_prob-1)<2e-11,(name,trace_prob)
        assert abs(leaked)<2e-11
        assert nonscalar<3e-11
        output.append({
            "error":name,
            "after_ideal_CPTP_recovery_entanglement_fidelity":float(recovered),
            "decoded_probability":trace_prob,
            "unrecoverable_flagged_probability":float(leaked),
            "max_logical_Kraus_nonscalar_residual":nonscalar,
            "code_syndrome_probability":float(
                np.linalg.norm(k0,"fro")**2/LOGICAL),
            "electric_syndrome_probability":float(
                np.linalg.norm(k1,"fro")**2/LOGICAL),
        })
    return output


def compound_exchange_after_electric_diagnostic(C,Verror,group,E,sizes):
    """EXPLORATORY NEGATIVE CONTROL: a second fault gR, outside the set.

    Compute the SAME decoder's syndrome-captured probability and
    entanglement fidelity for R followed by each independent g.
    Lost syndrome weight is FLAGGED, not renormalized into a false
    perfect score. These gR composite operators were NEVER claimed
    to be in the single-step seven-Kraus correctable error set.
    """
    P=C@C.T
    Q=Verror@Verror.T
    R=np.eye(D)-2*gauge_fixed_normalized_matrix(E,sizes)
    rows=[]
    for i,g in enumerate(group):
        compound=gauge_fixed_normalized_matrix(g,sizes)@R
        k0=C.T@P@compound@C
        k1=Verror.T@Q@compound@C
        captured=sum(np.linalg.norm(k,"fro")**2/LOGICAL for k in (k0,k1))
        fidelity=sum(abs(np.trace(k)/LOGICAL)**2 for k in (k0,k1))
        rows.append({
            "group_element_index":i,
            "two_step_operator":"g R (electric phase then matter exchange)",
            "recovered_F_e_without_postselection":float(fidelity),
            "captured_syndrome_probability":float(captured),
            "flagged_outside_syndrome_probability":float(1-captured),
        })
    assert abs(rows[0]["recovered_F_e_without_postselection"]-1)<2e-11
    assert any(x["recovered_F_e_without_postselection"]<.2 for x in rows)
    assert any(x["flagged_outside_syndrome_probability"]>.5 for x in rows)
    return rows


def common_gauss_physical_crosscheck(classes,sizes,C,Verror,label):
    """Embed constructed logical and syndrome subspaces into 486D gauge W.

    Independent local gauge transformations were checked for W in PR47;
    this crosscheck uses precisely the PR61 canonical orbit W map, not
    a different Hilbert space or different link representation.
    """
    W,sizes_from_embedding=physical_orbit_isometry(classes)
    assert np.array_equal(sizes_from_embedding,np.array(sizes,dtype=float))
    logical=W@C
    syndrome=W@Verror
    assert logical.shape==(486,4) and syndrome.shape==(486,4)
    assert np.linalg.norm(logical.conj().T@syndrome)<2e-11
    assert np.linalg.norm(logical.conj().T@logical-np.eye(4))<2e-11
    from colour_gauge_link_s3 import transform,G
    for vertex in ("A","B"):
        for h in range(G):
            assert np.linalg.norm(transform(logical,h,vertex)-logical)<3e-11
            assert np.linalg.norm(transform(syndrome,h,vertex)-syndrome)<3e-11
    return {
        "physical_code_embedding_shape":[486,4],
        "logical_and_syndrome_both_independent_vertex_Gauss_invariant":True,
        "physical_logical_to_syndrome_overlap_norm":
            float(np.linalg.norm(logical.conj().T@syndrome)),
        "group":label,
    }


def study():
    classes,lookup,ops=exact_generators()
    A,B,S,E=ops
    sizes=[len(x) for x in classes]
    results=[]
    for label,gens in (("AS",(A,S)),("BS",(B,S))):
        group=closure_matrices(gens)
        proof=exact_restricted_algebra_certificate(label,group,E)
        C,V,code=construct_active_code(label,group,E,sizes)
        branches=verify_five_plus_one_kraus_recovery(C,V,group,E,sizes)
        compounds=compound_exchange_after_electric_diagnostic(C,V,group,E,sizes)
        physical=common_gauss_physical_crosscheck(classes,sizes,C,V,label)
        results.append({
            "label":label,
            "exact_no_4D_passive_noiseless_subsystem":proof,
            "active_four_logical_code":code,
            "ideal_recovered_error_branches":branches,
            "outside_error_set_sequential_fault_negative_control":compounds,
            "486D_physical_Gauss_embedding":physical,
            "passive_common_scalar_sector_max_dimension":2,
        })
    return {
        "system":"four colour qutrits and one six-state finite S3 gauge link, physical 14D subspace",
        "operator_error_set":"each of six subgroup exchange unitaries OR one electric central-irrep pi-phase R individually; all 7 Kraus operators are allowed",
        "results":results,
        "active_advantage_against_passive_at_abstract_KL_level":True,
        "full_instrument_native_hardware_compiled":False,
        "virtual_intermediate_physics":{
            "Hubbard_previous_PR55_to_PR58":"explicit high-energy doublon-hole occupations and virtual superexchange retained across active pulses",
            "this_exact_finite_group_model":"fixed four-colour-matter and six-state link; no extra particle number, virtual gluon fields or QCD loops",
            "claim":"Virtual QFT contributions may change microscopic dynamics and fault distributions but do not affect this exact conditional finite-model KL theorem without modifying its Hamiltonian/noise assumptions"
        },
        "limitations":[
            "A mathematical four-state actively correctable code is established against exactly 7 alternatives for each proper six-element exchange subgroup plus R; it FAILS for some sequential gR compound faults and does not cover Rg or the full 24-swap family",
            "The exact 11D error-algebra commutant rules out a four-logical noiseless subsystem, stronger than just absence of a four-dimensional passive common eigenspace",
            "The recovery uses ideal projectors onto 4D code and 4D electric syndrome subspaces and a partial isometry; no native gauge-preserving gate decomposition is supplied",
            "The two codewords' electric branches overlap with amplitude minus one third; treating error outcomes as orthogonal before syndrome projection would be wrong",
            "All 7 errors are abstract finite-S3 gauge-admissible unitaries, with uncalibrated relative physical rates and unknown correlated fault paths",
            "Code preparation, pulse time, Hubbard virtual doublon leakage, gauge-link transport, active fermion T1/Tphi and reference/pointer overhead not included in the new KL advantage result",
            "No demonstration that this active code beats an optimized passive code under a realistic matched microscopic noise distribution",
            "The finite-S3 representation is NOT continuous SU3; no virtual gluon/quark-antiquark loops, actual Standard Model gauge dynamics or original nine-quark 21D 65-syndrome physical recovery",
        ],
    }


def main():
    data=study()
    Path("colour_s3_order_six_active_electric_results.json").write_text(
        json.dumps(data,indent=2,sort_keys=True)+"\n")
    for result in data["results"]:
        print("ORDER_SIX_ACTIVE_KL_CODE",json.dumps({
            "group":result["label"],
            "commutant":result["exact_no_4D_passive_noiseless_subsystem"]["joint_error_commutant_exact_dimension"],
            "largest_passive_common_scalar_sector":result["passive_common_scalar_sector_max_dimension"],
            "logical_dim":result["active_four_logical_code"]["logical_dimension"],
            "max_KL_residual":max(result["active_four_logical_code"][
                "numerical_Knill_Laflamme_residuals"].values()),
            "min_recovered_F_e":min(x[
                "after_ideal_CPTP_recovery_entanglement_fidelity"]
                for x in result["ideal_recovered_error_branches"]),
        },sort_keys=True),flush=True)


if __name__=="__main__":
    main()
