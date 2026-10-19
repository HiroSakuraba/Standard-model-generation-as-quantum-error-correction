"""Finite-S3 LOCAL gauge-covariant control alphabet for PR64's 42D recovery.

Central progress:
- The 14D physical Gauss sector of the 486D four-colour-qutrit /
  six-link model is fully controllable using FOUR physically supported
  gauge-invariant Hermitian Hamiltonians:
  AS: H=(onsite A0-A1 swap, dressed A0-B0 swap S,
         one-link electric P_standard, local B0-B1 equal-colour projector).
  BS: H=(onsite B0-B1 swap, dressed S, P_standard,
         local A0-A1 equal-colour projector).
- Independent finite-field Lie-closure ranks 196 (91 REAL skew plus
  105 IMAGINARY symmetric dimensions), exactly saturating u(14).
  Field 101 and 103 calculations use only RATIONAL generator matrices
  and nested Lie commutators: full controllability has an EXACT
  rational-rank certificate, not a floating-point approximation.
- WITHOUT the opposite-vertex equal-colour projector, the same
  generator Lie algebra has dimension only 21; using the same-vertex
  projector reaches 106, not full control.
- The nonlocal-looking three-register dressed swap exp(-i theta(I-S)/2)
  can be IMPLEMENTED EXACTLY with FIVE gauge-covariant TWO-register
  charged-mediator operations:
        SWAP(m,A0), transport(m), exp(-i theta(I-SWAP(m,B0))/2),
        transport(m), SWAP(m,A0).
  The charged mediator returns in its original singlet with reference.
- Add arbitrary neutral qutrit single-register drives and a tunable
  two-register neutral-pointer/link interaction |1><1|⊗P_standard.
  Standard Lie ideal arguments give full su(42) on 14D Gauss × pointer
  3, hence the PR64 global 42D coherent recovery is REACHABLE IN
  PRINCIPLE from strictly gauge-preserving gates of support <=2
  physical registers PLUS a mobile colour-charged bus.

No finite optimized pulse list for PR64's full recovery is produced.
Lie controllability does NOT imply low gate count, feasible coupling
strengths, acceptable T1/time, or a native device with an exactly
controllable charged mediator. The equal-colour onsite pair projector
is gauge invariant for DISCRETE S3, not general continuous SU3.
"""
from __future__ import annotations

from pathlib import Path
import json
import numpy as np
import sympy as sp

from colour_gauge_link_s3 import (
    DIM,G,transform,physical_isometry,dressed_exchange,
)
from colour_s3_exact_electric_commutant import (
    exact_generators,physical_orbit_isometry,
)
from colour_s3_order_six_active_electric import (
    closure_matrices,ordered_pair_centralizer_basis,
    construct_active_code,gauge_fixed_normalized_matrix,
)
from colour_s3_correlated_recovery_cadence import (
    global_gauge_preserving_instrument,
)
from colour_s3_charged_bus_transport import (
    initialize,local_swap,transport,involution_pulse,
    gauss_acceptance,
)
from colour_s3_strict_gauss_pair_nogo import (
    apply_link,link_central_projectors,
)

SYSTEM_DIM=14
NLOGICAL=4
FIELD_PRIMES=(101,103)
MATTER_DIMS=(3,3,3,3)
FUNDAMENTAL_GATE_SUPPORT_MAX=2


def same_colour_local_projector(classes,site_pair):
    """Exact 14x14 rational/INTEGER local two-matter equality projector.

    Q_v=sum_colour |cc><cc| for the two matters at ONE vertex.
    A simultaneous S3 signed permutation maps equal colour pairs to
    equal colour pairs; Q_v is locally gauge invariant. Crucially,
    Q_v is *not* a continuous SU3-invariant colour observable.
    """
    if site_pair not in ((0,1),(2,3)):
        raise ValueError("Only two matters on one gauge vertex allowed")
    i,j=site_pair
    Q=sp.zeros(SYSTEM_DIM)
    for k,orbit in enumerate(classes):
        outcomes={word[i]==word[j] for word in orbit}
        assert len(outcomes)==1
        Q[k,k]=int(outcomes.pop())
    assert Q*Q==Q and Q==Q.T
    return Q


def same_colour_full_486(y,site_pair):
    """Apply a true two-register vertex-local equal-colour projector."""
    if site_pair not in ((0,1),(2,3)):
        raise ValueError(site_pair)
    x=np.asarray(y,dtype=complex)
    single=x.ndim==1
    if single:x=x[:,None]
    if x.ndim!=2 or x.shape[0]!=DIM:raise ValueError(x.shape)
    shape=(G,3,3,3,3,x.shape[-1])
    colours=np.indices(MATTER_DIMS)
    i,j=site_pair
    mask=(colours[i]==colours[j]).reshape((1,)+MATTER_DIMS+(1,))
    out=(x.reshape(shape)*mask).reshape(x.shape)
    return out[:,0] if single else out


def admissible_controls(classes,ops,label):
    """Specified support <=3, before charged-bus substitution for S."""
    A,B,S,E=ops
    if label=="AS":
        Q=same_colour_local_projector(classes,(2,3))
        H=(A,S,E,Q)
        inventory=(
            ("onsite_swap_A0_A1",2),
            ("dressed_swap_link_A0_B0",3),
            ("electric_link_standard_projector",1),
            ("same_colour_projector_B0_B1",2))
    elif label=="BS":
        Q=same_colour_local_projector(classes,(0,1))
        H=(B,S,E,Q)
        inventory=(
            ("onsite_swap_B0_B1",2),
            ("dressed_swap_link_A0_B0",3),
            ("electric_link_standard_projector",1),
            ("same_colour_projector_A0_A1",2))
    else:raise ValueError(label)
    return H,inventory


def rational_mod_p(M,p):
    """Convert a SymPy rational matrix to the exact prime field F_p."""
    if p<=3:raise ValueError("Prime must avoid thirds")
    return np.array([
        [(int(x.p)*pow(int(x.q),-1,p))%p for x in row]
        for row in M.tolist()],dtype=np.int64)


def finite_field_rank(M,p):
    """Exact rank mod p, bounded 14x14 matrix products/integers."""
    a=np.asarray(M,dtype=np.int64).copy()%p
    nr,nc=a.shape
    row=0
    for col in range(nc):
        inds=np.flatnonzero(a[row:,col])
        if len(inds)==0:continue
        pivot=row+int(inds[0])
        a[[row,pivot]]=a[[pivot,row]]
        a[row,col:]=(a[row,col:]*pow(int(a[row,col]),-1,p))%p
        if row+1<nr:
            multiples=a[row+1:,col].copy()
            a[row+1:,col:]=(a[row+1:,col:]-
                             multiples[:,None]*a[row,col:])%p
        row+=1
        if row==nr:break
    return row


def exact_modular_lie_closure(ops,p):
    """Return exact Lie dimension of {i*H_j} via prime-field witnesses.

    Each H_j is rational and Hermitian in the PHYSICAL normalized orbit
    basis (their unnormalized orbit-sum similarity changes no rank).
    All nested commutators have rational coefficient matrices times
    i^parity. Real-skew and imaginary-symmetric anti-Hermitian blocks
    are independent. We maintain two disjoint pivot bases over F_p
    for the real/imaginary coefficient types, and repeatedly bracket
    each newly independent word with every original iH_j.

    Once dimensions sum to 14^2=196, the exact rational Lie rank is
    AT LEAST 196 (minor nonzero mod p) and at MOST 196 (subalgebra of
    u(14)). Thus it is EXACTLY u(14) in characteristic zero.
    A smaller mod-p rank is only a lower bound, never a no-go.
    """
    matrices=[rational_mod_p(H,p) for H in ops]
    n=matrices[0].shape[0]
    assert n==SYSTEM_DIM and all(M.shape==(n,n) for M in matrices)
    # parity 0 = real commutator word; 1 = i times rational matrix
    pivots=[{},{}]
    words=[]
    def add(M,parity):
        v=(M.reshape(-1)%p).copy()
        for pivot,u in sorted(pivots[parity].items()):
            if v[pivot]:
                v=(v-int(v[pivot])*u)%p
        nz=np.flatnonzero(v)
        if len(nz)==0:return False
        pivot=int(nz[0])
        v=(v*pow(int(v[pivot]),-1,p))%p
        pivots[parity][pivot]=v
        words.append((M.copy(),parity))
        return True
    for M in matrices:add(M,1)
    index=0
    while index<len(words):
        M,parity=words[index]
        index+=1
        for H in matrices:
            # [i^parity M, i H] toggles parity, with an irrelevant sign.
            comm=(M@H-H@M)%p
            add(comm,1-parity)
        if len(words)>n*n:
            raise AssertionError("More than u(n) anti-Hermitian degrees")
    return {
        "prime":p,
        "Lie_rank_mod_prime_lower_bound":len(words),
        "real_skew_rank":len(pivots[0]),
        "imaginary_symmetric_rank":len(pivots[1]),
        "full_u14_rational_rank_certified":len(words)==SYSTEM_DIM**2,
        "rank_method":"nested rational Lie-word commutators, mod-prime row elimination; exact 196 lower bound",
    }


def algebra_commutant_prime(group,E,Q,p):
    """Independent exact module commutant lower-rank crosscheck.

    Exchange subgroup centralizer has 52 integer orbit indicators.
    Adding E and Q gives 196x52 linear constraints twice. Rank 51
    mod p, while identity is an exact rational nullvector => rank<=51.
    Therefore the full commutant dimension is EXACTLY one.
    """
    basis=ordered_pair_centralizer_basis(group)
    assert len(basis)==52
    columns=[]
    for O in (E,Q):
        rows=[sp.Matrix(list(X*O-O*X)) for X in basis]
        coefficients=sp.Matrix.hstack(*rows)
        columns.append(rational_mod_p(coefficients,p))
    rank=finite_field_rank(np.vstack(columns),p)
    assert rank<=51
    return {
        "prime":p,
        "exchange_commutant_basis_dimension":len(basis),
        "electric_plus_pair_equality_constraint_rank":rank,
        "joint_commutant_dimension_upper_bound":len(basis)-rank,
        "identity_is_exact_nullvector":True,
        "exact_scalar_only_commutant":rank==51,
    }


def gauge_covariance_local_projector(pair):
    """Check full nonphysical 486D input, both local gauge actions."""
    rng=np.random.default_rng(6419+sum(pair))
    x=rng.normal(size=(DIM,2))+1j*rng.normal(size=(DIM,2))
    P=same_colour_full_486
    assert np.linalg.norm(P(P(x,pair),pair)-P(x,pair))<2e-12
    worst=0.
    for vertex in ("A","B"):
        for h in range(G):
            diff=transform(P(x,pair),h,vertex)-P(transform(x,h,vertex),pair)
            worst=max(worst,float(np.linalg.norm(diff)))
    assert worst<2e-11
    return worst


def controlled_link_electric_full_486(z,theta):
    """Exact TWO-register (neutral pointer + link) electric phase.

    U=exp(-i theta |1><1|_pointer tensor E_link). Central E_link is a
    projector and invariant under both independent vertex Gauss actions.
    Input shape (3,486,arbitrary logical reference columns).
    """
    x=np.asarray(z,dtype=complex)
    if x.ndim!=3 or x.shape[0]!=3 or x.shape[1]!=DIM:
        raise ValueError(x.shape)
    if not np.isfinite(theta):raise ValueError("finite theta")
    E=link_central_projectors()["standard"]
    out=x.copy()
    out[1]+=(np.exp(-1j*theta)-1)*apply_link(x[1],E)
    return out


def charged_mediator_fractional_dressed_swap(y,theta,history=False):
    """EXACT five TWO-register pulse implementation for arbitrary angle.

    Original full physical unitary is exp(-i theta (I-S_dressed)/2).
    Use W=T*SWAP(m,A0), so S=W† SWAP(m,B0) W and
    U=W† exp(-i theta(I-SWAP(m,B0))/2) W.

    Charged bus and stationary reference must be in a gauge singlet
    with mediator at A at INPUT. All intermediate ops are the already
    PR53-verified gauge-covariant local 2-register primitives. At
    OUTPUT the charged bus returns, without state-dependent residual.
    """
    if not np.isfinite(theta):raise ValueError("finite angle required")
    z=initialize(y,pointer_uniform=False)
    intermediate=[]
    z=local_swap(z,"A0",active=(0,))
    intermediate.append(z.copy() if history else None)
    z=transport(z)
    intermediate.append(z.copy() if history else None)
    z=involution_pulse(z,lambda x:local_swap(x,"B0",active=(0,)),theta)
    intermediate.append(z.copy() if history else None)
    z=transport(z)
    intermediate.append(z.copy() if history else None)
    z=local_swap(z,"A0",active=(0,))
    f=np.exp(-1j*theta)
    expected=((1+f)/2)*y+((1-f)/2)*dressed_exchange(y)
    target=initialize(expected,pointer_uniform=False)
    residual=float(np.linalg.norm(z-target))
    assert residual<3e-10
    out={
        "angle":float(theta),
        "two_register_primitive_count":5,
        "two_register_conjugating_swaps_and_transports":4,
        "one_fractional_two_register_charged_mediator_local_exchange":1,
        "original_three_register_S_fractional_pulse_reconstruction_Fro":
            residual,
        "stationary_reference_singlet_restored":True,
        "charged_mediator_returns_to_A":True,
    }
    if history:
        accept=[gauss_acceptance(v) for v in intermediate+[z]]
        assert all(min(x.values())>1-5e-11 for x in accept)
        out["intermediate_Gauss_acceptance"]=accept
    return out


def pointer_universality_argument(E):
    """Analytic su(42) controllability from system u14 + one entangler.

    Controls: local system u14, neutral-pointer su3, and nontrivial
    H_int=P_pointer_1 tensor E_system, physically pointer + 6-state
    link only. Lie ideals of the simple algebras su14 and su3 show
    [su14,H_int] generates su14 tensor traceless P_1; pointer
    commutators generate su14 tensor su3. Together with the two local
    Lie algebras these span su(42). Trace/identity controls generate
    u(42) up to arbitrary global phase.

    This is a theorem for independently tunable ideal couplings with
    both signs or compact-group recurrence; not an output pulse list.
    """
    assert E*E==E
    rank=int(sp.trace(E))
    assert 0<rank<SYSTEM_DIM
    return {
        "physical_system_dimension":SYSTEM_DIM,
        "neutral_pointer_dimension":3,
        "physical_plus_pointer_dimension":SYSTEM_DIM*3,
        "pointer_single_register_local_su3_controls_assumed":True,
        "pointer_link_entangling_two_register_H":"|1><1|_pointer tensor P_standard(link)",
        "coupled_link_projector_rank_on_physical_system":rank,
        "coupled_link_projector_non_scalar":True,
        "generated_traceless_Lie_algebra_by_ideal_argument":"su(42)",
        "traceless_dimension":(SYSTEM_DIM*3)**2-1,
        "full_unitary_reachable_in_ideal_continuous_control_model":True,
        "explicit_solved_pulse_sequence_for_PR64_global_recovery":False,
        "control_signs_and_unbounded_precision_assumed":True,
    }


def study():
    classes,_,ops=exact_generators()
    A,B,S,E=ops
    W,sizes=physical_orbit_isometry(classes)
    results=[]
    for label,generators in (("AS",(A,S)),("BS",(B,S))):
        H,inventory=admissible_controls(classes,ops,label)
        Q=H[-1]
        group=closure_matrices(generators)
        for h in H:
            # Rational orbit-sum matrices are self-adjoint in the
            # physical orbit-norm, not necessarily Euclidean transpose.
            physical=gauge_fixed_normalized_matrix(h,sizes)
            assert np.linalg.norm(physical-physical.conj().T)<3e-11
        rank_checks=[exact_modular_lie_closure(H,p) for p in FIELD_PRIMES]
        assert all(v["full_u14_rational_rank_certified"] for v in rank_checks)
        no_Q=[exact_modular_lie_closure(H[:-1],p) for p in FIELD_PRIMES]
        assert all(v["Lie_rank_mod_prime_lower_bound"]==21 for v in no_Q)
        pair=[algebra_commutant_prime(group,E,Q,p) for p in FIELD_PRIMES]
        assert all(x["exact_scalar_only_commutant"] for x in pair)
        pair_sites=(2,3) if label=="AS" else (0,1)
        pair_cov=gauge_covariance_local_projector(pair_sites)
        C,V,_=construct_active_code(label,group,E,sizes)
        instrument=global_gauge_preserving_instrument(C,V)
        assert instrument["joint_unitary_dimension"]==42
        results.append({
            "group":label,
            "gauge_invariant_control_Hamiltonian_inventory":[
                {"name":name,"register_support_before_charged_bus":support}
                for name,support in inventory],
            "opposite_vertex_same_colour_projector_pair":list(pair_sites),
            "full_486D_same_colour_projector_Gauss_covariance_Fro":
                pair_cov,
            "independent_exact_finite_field_system_Lie_ranks":rank_checks,
            "without_opposite_vertex_local_projector_lie_rank":no_Q,
            "joint_associative_commutant_exact_scalar_rank_certificate":pair,
            "existing_three_outcome_recovery_target":instrument,
        })
    # Fractional edge exchange pulse uses 5 mediator-pair gates.
    # Columns chosen generically from 486D physical W rather than
    # from only an S=+1 eigen-sector: exercise actual nontrivial S.
    bus=[charged_mediator_fractional_dressed_swap(
        W[:,:3],angle,history=(j==2))
         for j,angle in enumerate((0.,np.pi/5,np.pi/2,np.pi))]
    # Gauge covariance of pointer-link electric coupling off physical
    # sector, including both local independent actions.
    rng=np.random.default_rng(9921)
    z=rng.normal(size=(3,DIM,2))+1j*rng.normal(size=(3,DIM,2))
    theta=.37
    electric_worst=0.
    for vertex in ("A","B"):
        for h in range(G):
            from colour_gauge_link_s3 import transform
            lhs=controlled_link_electric_full_486(
                np.stack([transform(y,h,vertex) for y in z]),theta)
            rhs=np.stack([
                transform(y,h,vertex) for y in
                controlled_link_electric_full_486(z,theta)])
            electric_worst=max(electric_worst,float(np.linalg.norm(lhs-rhs)))
    assert electric_worst<2e-10
    pointer=pointer_universality_argument(E)
    return {
        "system":"two independent local Gauss constraints, finite S3, physical 14D and neutral pointer 3",
        "charged_bus_source":"PR53 existing 2-register charged-mediator transport and vertex-local swaps",
        "subgroup_local_control_results":results,
        "continuous_dressed_swap_five_pair_gate_exact_identity":bus,
        "neutral_pointer_link_electric_pair_covariance_Fro":
            electric_worst,
        "physical_plus_pointer_controllability_argument":pointer,
        "core_scientific_claim":"An explicitly supported strictly gauge-invariant control alphabet reaches all U14 on the physical system and U42 with neutral pointer in the ideal tunable Hamiltonian model; a finite numerical pulse compiler has NOT been built.",
        "claimed_finite_native_compiled_42D_PR64_recovery":False,
        "limitations":[
            "The finite-field nested Lie-word rank 196 is exact for the rational 14D finite-S3 matrices and independent controllable Hamiltonians; it says NOTHING about pulse count, bandwidth or robustness.",
            "A local opposite-vertex two-matter SAME-COLOUR projector is required by this particular minimal four-Hamiltonian universality proof. It respects discrete S3 but NOT arbitrary continuous SU3 gauge transformations.",
            "The 3-register link-dressed exchange continuous pulse is expressed as 5 gauge-covariant TWO-register charged-mediator steps, but requires a fully controllable mobile charged carrier, colour singlet preparation, and tunable onsite mediator exchange.",
            "The pointer-link central-electric coupling is one native ASSUMED gauge-invariant 2-register Hamiltonian, not shown to follow from a specific hardware interaction.",
            "Neutral pointer arbitrary SU3 control and independently tunable +/- pulse times or recurrence are assumed; universal controllability on U42 is by Lie-ideal theorem, not a 42D numerical pulse optimization.",
            "No explicit pulse sequence, depth, wall-clock time, calibration, residual leakage, or T1/T2 fidelity for the COMPLETE three-outcome PR64 recovery instrument has been produced.",
            "The exact charged-bus five-gate continuous S pulse restores the mediator only for the prepared gauge-invariant charged-reference singlet; virtual doublon-hole or QCD virtual gluon fields are NOT included here.",
            "The physical Hilbert space remains finite-S3 four matter qutrits and one group link, NOT full SU3 or the nine-quark 21D/65-syndrome standard-model code.",
        ],
    }


def main():
    result=study()
    Path("colour_s3_local_control_universality_results.json").write_text(
        json.dumps(result,indent=2,sort_keys=True)+"\n")
    for case in result["subgroup_local_control_results"]:
        print("GAUGE_LOCAL_U14_EXACT_CERTIFICATE",json.dumps({
            "group":case["group"],
            "Lie_ranks":[x["Lie_rank_mod_prime_lower_bound"]
                         for x in case["independent_exact_finite_field_system_Lie_ranks"]],
            "Lie_ranks_without_onsite_projector":[x["Lie_rank_mod_prime_lower_bound"]
                    for x in case["without_opposite_vertex_local_projector_lie_rank"]],
            "commutant_ranks":[x["electric_plus_pair_equality_constraint_rank"]
                               for x in case["joint_associative_commutant_exact_scalar_rank_certificate"]],
        },sort_keys=True))
    print("GAUGE_LOCAL_CHARGED_BUS_FRACTIONAL_PULSE",json.dumps([
        {"angle":x["angle"],
         "two_register_primitives":x["two_register_primitive_count"],
         "residual":x["original_three_register_S_fractional_pulse_reconstruction_Fro"]}
         for x in result["continuous_dressed_swap_five_pair_gate_exact_identity"]
    ],sort_keys=True))


if __name__=="__main__":
    main()
