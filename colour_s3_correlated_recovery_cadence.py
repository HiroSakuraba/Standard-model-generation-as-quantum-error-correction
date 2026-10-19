"""Order-sensitive correlated exchange/electric faults and recovery cadence.

Study the *same* two six-element gauge-admissible exchange groups and
four-dimensional code from PR #63, without changing its logical code.
The error R=I-2P_standard acts on the link; g is any of the six group
permutations. The order is consequential:

  Rg C = R C,  because g C=C: existing one-round recovery works.
  gR C generally exits the existing 8D two-syndrome space:
      for four of six g, P_C+P_V captures only 11/24 and
      unpostselected F_e=11/96, not 1.

Inserting an IDEAL syndrome-and-recovery ROUND after the first fault
restores C before the second fault: the two-round channel perfectly
corrects *any such two-step sequence*, in exact noiseless arithmetic.
No claim of fault-tolerance beyond separated single-step events,
or of a native implementation of the recovery.

A coherent, explicitly GAUGE-PRESERVING 14D-plus-neutral-pointer
unitary is given:
  U_synd = (I-P_V) tensor I_2 + P_V tensor X,
  U_feed = I_14 tensor |0><0| + F tensor |1><1|,
  F = C V† + V C† + (I-P_C-P_V).
This proves a finite-dimensional isometric syndrome instrument
and recovery, but these projectors operate on the GLOBAL gauge
physical space and are NOT synthesized from affordable local gates.

Compare independent electric probability p_e and exchange probability
p_g, with f of double events ordered electric-before-exchange (gR),
the rest Rg, and a uniform conditional choice of six g.
The channels are CPTP mixtures of exact unitaries, and ALL branches
are included. For this uniform group model:
  one_round_F_e = 1 - (85/144)*p_e*p_g*f,
  two_round_F_e = 1
before physical loss or extra recovery imperfections.

Four active fermions with hypothetical independent T1 have survival
eta=exp(-4 * duration/T1). One extra equal-duration recovery round
multiplies by another eta. A passive reference restricted to the
six-dimensional group-trivial space can choose the four eigenvectors
maximizing |Tr(C†R C)|, for no-feedback electric F_e=4/9.
This is the OPTIMUM among four-dimensional G-INVARIANT subspaces
under that one electric fault, NOT proof of globally optimal passive
performance under arbitrary finite-T noise.

THIS MODEL HAS NO DEVICE SPECIFIC PULSE TIMES, calibrated rate, native
pointer coupling, microscopic link electric fault or QCD loops.
"""
from __future__ import annotations

from pathlib import Path
import json
import numpy as np

from colour_s3_exact_electric_commutant import exact_generators
from colour_s3_order_six_active_electric import (
    D,LOGICAL,closure_matrices,construct_active_code,
    gauge_fixed_normalized_matrix,
)

N_GROUP=6
FOUR_FERMIONS=4
COMPOUND_GOOD_FRACTION=2/6
BAD_COMPOUND_F_E=11/96
BAD_COMPOUND_CAPTURE=11/24
GROUP_AVERAGE_BAD_F_E=(2+4*BAD_COMPOUND_F_E)/6
COMPOUND_LOGICAL_LOSS=1-GROUP_AVERAGE_BAD_F_E # 85/144
assert abs(COMPOUND_LOGICAL_LOSS-85/144)<1e-15


def branch_kraus(C,V,op):
    """Return 2 corrected logical branch Kraus and flagged residue.

    The global code and syndrome projections give P_C and P_V.
    After partial-isometry feedback V->C, the logical matrices are
    C† Op C and V† Op C. Do not renormalize the discarded complement.
    """
    k0=C.conj().T@op@C
    k1=V.conj().T@op@C
    capture=sum(float(np.linalg.norm(x,"fro")**2/LOGICAL) for x in (k0,k1))
    fe=sum(float(abs(np.trace(x)/LOGICAL)**2) for x in (k0,k1))
    return (k0,k1),{
        "unconditional_F_e":fe,
        "captured_probability":capture,
        "flagged_probability":float(1-capture),
    }


def corrected_two_round(C,V,early,late):
    """Measure/recover between TWO errors; independent pointer per round."""
    ka,_=branch_kraus(C,V,early)
    kb,_=branch_kraus(C,V,late)
    compound=[right@left for left in ka for right in kb]
    fidelity=sum(float(abs(np.trace(x)/LOGICAL)**2) for x in compound)
    captured=sum(float(np.linalg.norm(x,"fro")**2/LOGICAL) for x in compound)
    return {"unconditional_F_e":fidelity,
            "captured_probability":captured,
            "flagged_probability":float(1-captured),
            "decoder_rounds":2}


def global_gauge_preserving_instrument(C,V):
    """Explicit 28x28 syndrome+feedback unitary (NOT a local gate circuit).

    Neutral pointer, initially |0>. U_synd flips it on V; U_feed
    applies F on pointer=1. Both are mathematically unitary. Since C
    and V sit in the simultaneous physical Gauss sector, their
    projectors lifted as W P W† commute with independent local Gauss
    actions. That does not make the gates few-body or fast.
    """
    I=np.eye(D)
    P=C@C.conj().T
    Q=V@V.conj().T
    X=np.array([[0.,1.],[1.,0.]])
    z=np.diag([1.,0.])
    o=np.diag([0.,1.])
    F=C@V.conj().T+V@C.conj().T+I-P-Q
    Us=np.kron(I-Q,np.eye(2))+np.kron(Q,X)
    Uf=np.kron(I,z)+np.kron(F,o)
    U=Uf@Us
    assert np.linalg.norm(F.conj().T@F-I)<5e-11
    assert np.linalg.norm(U.conj().T@U-np.eye(28))<5e-11
    input_code=np.kron(C,np.array([[1.],[0.]]))
    input_error=np.kron(V,np.array([[1.],[0.]]))
    expected_code=input_code
    expected_error=np.kron(C,np.array([[0.],[1.]]))
    residual_code=float(np.linalg.norm(U@input_code-expected_code))
    residual_error=float(np.linalg.norm(U@input_error-expected_error))
    assert max(residual_code,residual_error)<5e-11
    return {
        "pointer_levels":2,
        "logical_code_rank":LOGICAL,
        "orthogonal_syndrome_rank":LOGICAL,
        "joint_unitary_dimension":28,
        "unitary_unverified_as_native_few_body":True,
        "full_logical_code_path_residual":residual_code,
        "full_error_syndrome_path_residual":residual_error,
        "coherent_syndrome_operator":"(I-P_V) tensor I_2 + P_V tensor X",
        "pointer_controlled_feedback_operator":"I_14 tensor |0><0| + F tensor |1><1|",
        "physical_register_count_before_pointer":"four matter qutrits plus one 6-level gauge link",
        "only_global_projector_implementation_shown":True,
    }


def exact_order_diagnostic(C,V,group,E,sizes):
    """Enumerate all six gR and six Rg, one versus two recovery rounds."""
    R=np.eye(D)-2*gauge_fixed_normalized_matrix(E,sizes)
    rows=[]
    for i,g_rational in enumerate(group):
        g=gauge_fixed_normalized_matrix(g_rational,sizes)
        for order,early,late in (
            ("electric_then_exchange_gR",R,g),
            ("exchange_then_electric_Rg",g,R),
        ):
            _,one=branch_kraus(C,V,late@early)
            two=corrected_two_round(C,V,early,late)
            rows.append({
                "group_index":i,"physical_fault_order":order,
                "one_round":one,"two_round":two,
            })
    for order in ("electric_then_exchange_gR","exchange_then_electric_Rg"):
        subset=[v for v in rows if v["physical_fault_order"]==order]
        if order=="electric_then_exchange_gR":
            bad=[x for x in subset if x["one_round"]["unconditional_F_e"]<.2]
            assert len(bad)==4
            for x in bad:
                assert abs(x["one_round"]["unconditional_F_e"]-11/96)<3e-11
                assert abs(x["one_round"]["captured_probability"]-11/24)<3e-11
        else:
            assert all(abs(x["one_round"]["unconditional_F_e"]-1)<3e-11
                       for x in subset)
        assert all(abs(x["two_round"]["unconditional_F_e"]-1)<3e-11
                   for x in subset)
    return rows


def valid_probability(v,name):
    if not np.isfinite(v) or not 0<=v<=1:
        raise ValueError(f"{name} must be finite and within [0,1]")


def matched_channel(C,V,group,E,sizes,p_e=.37,p_g=.4,
                    electric_first_fraction=.5):
    """EXACT CPTP mixture, including no fault and all single/double faults.

    Bernoulli electric event probability p_e; independent Bernoulli
    exchange event probability p_g; each g from full six-element
    subgroup equally likely if exchange occurs. Conditional on BOTH,
    the order is gR with fraction f, Rg otherwise. When only one
    type occurs, its usual single unitary fault is applied.
    """
    for n,x in (("p_e",p_e),("p_g",p_g),
                ("electric_first_fraction",electric_first_fraction)):
        valid_probability(x,n)
    R=np.eye(D)-2*gauge_fixed_normalized_matrix(E,sizes)
    I=np.eye(D)
    events=[("none",I,(1-p_e)*(1-p_g)),
            ("electric_only",R,p_e*(1-p_g))]
    for i,g_rational in enumerate(group):
        g=gauge_fixed_normalized_matrix(g_rational,sizes)
        events.extend([
            (f"exchange_only_{i}",g,(1-p_e)*p_g/N_GROUP),
            (f"gR_{i}",g@R,p_e*p_g*electric_first_fraction/N_GROUP),
            (f"Rg_{i}",R@g,p_e*p_g*(1-electric_first_fraction)/N_GROUP),
        ])
    total=sum(w for _,_,w in events)
    assert abs(total-1)<5e-12
    one=0.;capture=0.;single=0.
    event_rows=[]
    for name,op,w in events:
        _,record=branch_kraus(C,V,op)
        one+=w*record["unconditional_F_e"]
        capture+=w*record["captured_probability"]
        event_rows.append({
            "event":name,"probability":w,
            **record
        })
    expected=1-(85/144)*p_e*p_g*electric_first_fraction
    assert abs(one-expected)<2e-11,(one,expected)
    assert abs(single-0)<1e-12
    return {
        "p_electric":p_e,"p_exchange":p_g,
        "conditional_electric_before_exchange_fraction":electric_first_fraction,
        "independent_Bernoulli_probabilities":True,
        "exchange_draw_includes_identity_among_six":True,
        "single_round_conditional_F_e":float(one),
        "single_round_captured_probability":float(capture),
        "single_round_flagged_probability":float(1-capture),
        "two_rounds_with_ideal_recovery_between_fault_windows_F_e":1.,
        "analytic_single_round_F_e":"1-(85/144)*p_e*p_g*f",
        "enumerated_event_rows":event_rows,
        "all_21_events_probability_sum":float(total),
    }


def passive_group_trivial_best(C,group,E,sizes):
    """EXACT highest no-correction electric F_e within G-FIXED 4D codes.

    On 6D group-fixed space, R=I-2E compressed eigenvalues (sorted)
    are (-1,-1,-1/3,-1/3,1/3,1). Ky Fan trace principle:
    maximum absolute 4D trace is 8/3 (four lowest eigenvalues),
    giving |Tr C† R C/4|²=(2/3)²=4/9.

    This assumes the passive code is a four-state SUBSPACE inside
    G's trivial sector, not an arbitrary noiseless subsystem or
    unconstrained approximate passive implementation.
    """
    sizes_np=np.sqrt(np.asarray(sizes))
    P=sum(group)/N_GROUP
    E_phys=gauge_fixed_normalized_matrix(E,sizes)
    P_phys=gauge_fixed_normalized_matrix(P,sizes)
    eigs,V=np.linalg.eigh(P_phys)
    fixed=V[:,eigs>.5]
    assert fixed.shape==(D,6)
    Rfix=fixed.conj().T@(np.eye(D)-2*E_phys)@fixed
    spec=np.linalg.eigvalsh(Rfix)
    expected=np.array([-1,-1,-1/3,-1/3,1/3,1])
    assert np.linalg.norm(spec-expected)<3e-11
    four_negative=np.sum(spec[:4])
    four_positive=np.sum(spec[-4:])
    best=max(abs(four_negative),abs(four_positive))/LOGICAL
    assert abs(best-2/3)<3e-11
    return {
        "passive_reference_class":"all four-dimensional subspaces of 6D exchange-group FIXED sector",
        "spectrum_R_on_group_fixed_sector":[float(x) for x in spec],
        "max_abs_normalized_trace_electric":float(best),
        "best_conditional_F_e_for_electric_event":float(best**2),
        "matched_no_feedback_channel_F_e":"1-(5/9)*p_e, independent of exchange probability and order",
        "no_claim_globally_optimal_passive_gauge_code":True,
    }


def loss_survival(duration_us,T1_us,n_fermions=FOUR_FERMIONS):
    if not np.isfinite(duration_us) or duration_us<0 or not np.isfinite(T1_us) or T1_us<=0:
        raise ValueError("duration/T1 must be finite with nonnegative duration and positive lifetime")
    return float(np.exp(-n_fermions*duration_us/T1_us))


def lifetime_cadence_comparison(channel,passive,p_e,
                                round_duration_us=1.,T1_us=30.):
    """Fair timing: passive idles for same clock time as each option.

    One decoder recovery occurs after both possible faults, taking
    one assumed duration. Two recoveries bracket the fault windows
    and take two assumed durations. A no-gate passive comparator
    is independently idled for BOTH corresponding full durations.
    No actual physical pointer circuit durations are known; 1us is
    only a hypothetical sensitivity unit, not derived gate timing.
    """
    if not np.isfinite(round_duration_us) or round_duration_us<=0:
        raise ValueError("positive finite ideal recovery round duration")
    valid_probability(p_e,"p_e")
    eta=loss_survival(round_duration_us,T1_us)
    active_one=eta*channel["single_round_conditional_F_e"]
    active_two=eta**2*channel[
        "two_rounds_with_ideal_recovery_between_fault_windows_F_e"]
    passive_cond=1-(5/9)*p_e
    critical=None
    if channel["single_round_conditional_F_e"]<1:
        # Two rounds worth extra time iff eta > F_one.
        threshold=-FOUR_FERMIONS*round_duration_us/np.log(
            channel["single_round_conditional_F_e"])
        critical=threshold
    return {
        "hypothetical_duration_per_recovery_round_us":round_duration_us,
        "hypothetical_active_fermion_T1_us":T1_us,
        "equal_four_fermions_lost_as_flagged_erasures":True,
        "active_one_round_unconditional_F_e":float(active_one),
        "active_two_round_unconditional_F_e":float(active_two),
        "passive_same_single_round_duration_unconditional_F_e":
            float(eta*passive_cond),
        "passive_same_two_round_duration_unconditional_F_e":
            float(eta**2*passive_cond),
        "two_active_rounds_win_over_one_under_this_loss_only_model":
            bool(active_two>active_one),
        "break_even_T1_us_two_rounds_vs_one":critical,
        "not_a_device_calibrated_runtime":True,
        "no_correlated_faults_during_recovery_are_modelled":True,
    }


def study():
    classes,_,(A,B,S,E)=exact_generators()
    sizes=[len(c) for c in classes]
    cases=[]
    for label,gen in (("AS",(A,S)),("BS",(B,S))):
        G=closure_matrices(gen)
        C,V,_=construct_active_code(label,G,E,sizes)
        order=exact_order_diagnostic(C,V,G,E,sizes)
        instr=global_gauge_preserving_instrument(C,V)
        passive=passive_group_trivial_best(C,G,E,sizes)
        runs=[]
        for p_e,p_g,f in (
            (.37,.4,0.),(.37,.4,.5),(.37,.4,1.),
            (.37,1.,1.),(.05,.4,.5),
        ):
            channel=matched_channel(C,V,G,E,sizes,p_e,p_g,f)
            budgets=[lifetime_cadence_comparison(
                channel,passive,p_e,1.,T1) for T1 in (10.,30.,100.)]
            runs.append({"channel":channel,"illustrative_lifetime_scenarios":budgets})
        cases.append({
            "group":label,"four_logical_dimension":LOGICAL,
            "ordered_fault_unitary_audit":order,
            "explicit_GLOBAL_coherent_syndrome_unitary":instr,
            "passive_G_fixed_four_dim_optimum":passive,
            "two_time_window_stochastic_scenarios":runs,
        })
    return {
        "experiment":"correlated fault ordering, ideally interleaved syndrome and four-fermion lifetime tradeoff",
        "logical_code":"PR63 same 4D active finite S3 gauge-admissible code, unchanged",
        "comparisons":cases,
        "ideal_two_round_recovery_compound_faults":True,
        "hardware_native_recovery_compiled":False,
        "limitations":[
            "The seven single-step Kraus errors are individually corrected, but four of six electric-before-exchange gR composites are not corrected by ONE round",
            "Exchange-before-electric Rg is perfectly recovered in one round because every exchange g fixes the logical code",
            "The two-round perfect result assumes a perfect GLOBAL ideal recovery between fault windows; it does not mean arbitrary repeated faults within one round are corrected",
            "The 28D pointer and feedback unitaries are explicitly constructed and gauge-preserving when lifted to the full physical Gauss sector, but they rely on nonlocal 14D projectors and are NOT natively compiled",
            "The passive 4/9 optimum is exact only inside the six-dimensional exchange-invariant subspace, not a global all-14D code optimization",
            "Loss model is uniform four active fermions, flagged erasure; no loss of spectator, pointer, link or reference and no T2, calibration or explicit microscopic post-jump branches",
            "Round duration 1us is illustrative only; physical control depth, synthesis errors and resource costs of the new global syndrome are unknown",
            "Bernoulli electric/exchange rates and temporal-order fractions are synthetic; no QCD-induced rates, correlated baths or actual virtual gluon/sea-quark effects",
            "No native finite S3 link dynamics, continuous SU3, or nine-quark 21D/65-syndrome physical recovery",
        ],
    }


def main():
    data=study()
    Path("colour_s3_correlated_recovery_cadence_results.json").write_text(
        json.dumps(data,indent=2,sort_keys=True)+"\n")
    for case in data["comparisons"]:
        print("CORRELATED_ORDERED_FAULT_CADENCE",json.dumps({
            "group":case["group"],
            "gR_bad_cases":sum(
                r["one_round"]["unconditional_F_e"]<.2
                for r in case["ordered_fault_unitary_audit"]
                if r["physical_fault_order"]=="electric_then_exchange_gR"),
            "Rg_min_F_e":min(
                r["one_round"]["unconditional_F_e"]
                for r in case["ordered_fault_unitary_audit"]
                if r["physical_fault_order"]=="exchange_then_electric_Rg"),
            "two_round_min_F_e":min(
                r["two_round"]["unconditional_F_e"]
                for r in case["ordered_fault_unitary_audit"]),
            "passive_R_F_e":case["passive_G_fixed_four_dim_optimum"][
                "best_conditional_F_e_for_electric_event"],
        },sort_keys=True),flush=True)


if __name__=="__main__":
    main()
