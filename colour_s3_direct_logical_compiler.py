"""Direct seven-fault logical-fidelity pulse optimization for the finite-S3 code.

Research followup to merged PR66: the prior numerical compiler optimized
42D target-isometry overlap and performed WORSE than doing nothing on the
actual seven-fault logical channel. This module instead optimizes the
UNCONDITIONAL (no postselection) logical entanglement fidelity directly,
with optional nonnegative angle/charged-bus effort regularization.

On the PR63 four-logical code C, EACH of six subgroup exchanges g obeys
g C = C, so their uniform Kraus average is identical to clean code
with weight 1-p. The seventh fault is R=I-2 E with weight p. Thus
the exactly equivalent objective needs only eight quantum input states:
X_0=C tensor |0>, X_1=R C tensor |0>, four columns each.

For b=0,1 (recognized pointer outcomes), define
 a_jb = Tr[(C tensor |b>)† U(theta) X_j] (complex scalar).
The trace-decreasing recovered logical subchannel has
 F_e = sum_j w_j sum_b |a_jb|² / 16.
All pointer-2 outcomes AND physical states outside C in pointer0/1
are failures, never postselected or renormalized. The objective can
accept DIFFERENT coherent phases for the two orthogonal pointer
measurement outcomes, as physically required by a measured instrument.

Use analytic reverse-adjoint differentiation through each existing
42D physical/pointer gauge-compatible pulse. A nonnegative, smooth
angle-weighted effort term is optional; final expanded operations and
hypothetical serial duration are always computed separately.

The DO-NOTHING policy is retained as an exact warm-start/fallback
regardless of numerical optimizer success. No unearned correction
claim and no SU3/hardware-native equivalence.
"""
from __future__ import annotations

from pathlib import Path
import json

import numpy as np
from scipy.optimize import minimize

from colour_s3_recovery_pulse_compiler import (
    DIM,POINTER,TOTAL,SpectralCompiler,
    target_three_outcome,native_operation_expansion,
    single_fault_logical_recovery_fidelity,
)
from colour_s3_exact_electric_commutant import exact_generators
from colour_s3_order_six_active_electric import (
    closure_matrices,construct_active_code,gauge_fixed_normalized_matrix,
)


def logical_fault_states(label,p_electric=.37):
    """8 exact code-input columns and two pointer readout projections.

    Six possible g errors reduce EXACTLY to the same gC=C block.
    This reduction is only valid for PR63's restricted single-fault
    family; NOT for compound errors gR, Rg or a general noisy channel.
    """
    if not np.isfinite(p_electric) or not 0<=p_electric<=1:
        raise ValueError("p_electric must be a finite probability")
    classes,_,(A,B,S,E)=exact_generators()
    sizes=[len(c) for c in classes]
    group=closure_matrices((A,S) if label=="AS" else (B,S))
    C,_,_=construct_active_code(label,group,E,sizes)
    R=np.eye(DIM)-2*gauge_fixed_normalized_matrix(E,sizes)
    ptr=np.eye(POINTER)
    input_clean=np.kron(C,ptr[:,[0]])
    input_fault=np.kron(R@C,ptr[:,[0]])
    X=np.column_stack((input_clean,input_fault))
    projectors=np.stack([
        np.kron(C,ptr[:,[0]]),
        np.kron(C,ptr[:,[1]])
    ]) # (2,42,4)
    assert X.shape==(TOTAL,8)
    assert projectors.shape==(2,TOTAL,4)
    assert np.linalg.norm(X.conj().T@X-np.block([
        [np.eye(4),-np.eye(4)/3],
        [-np.eye(4)/3,np.eye(4)]
    ]))<2e-10
    return {
        "inputs":X,
        "readout_code_isometries":projectors,
        "weights":np.array([1-p_electric,p_electric]),
        "electric_probability":p_electric,
        "identity_reference_logical_F_e":float((1-p_electric)+p_electric/9),
        "reduced_fault_blocks":2,
        "individual_Kraus_errors_represented":7,
    }


def effort_weight_for_name(name):
    """Heuristic, NOT calibrated hardware ns or a native gate count."""
    if name.startswith("dressed_swap_"):return 5.
    if name=="pointer_link_electric_entangler":return 2.
    return 1.


class DirectLogicalCompiler:
    def __init__(self,label,cycles,p_electric=.37,
                 effort_lambda=0.,effort_epsilon=.02):
        if not np.isfinite(effort_lambda) or effort_lambda<0:
            raise ValueError("effort_lambda must be finite nonnegative")
        if not np.isfinite(effort_epsilon) or effort_epsilon<=0:
            raise ValueError("effort_epsilon must be positive finite")
        self.base=SpectralCompiler(label,cycles)
        self.label=label
        self.cycles=cycles
        self.channel=logical_fault_states(label,p_electric)
        self.initial=self.channel["inputs"]
        self.readout=self.channel["readout_code_isometries"]
        self.weights=self.channel["weights"]
        self.p_electric=p_electric
        self.effort_lambda=float(effort_lambda)
        self.effort_epsilon=float(effort_epsilon)
        self.names=[self.base.named[i][0] for i in self.base.order]
        self.weighted_effort=np.array([
            effort_weight_for_name(name) for name in self.names
        ],dtype=float)
        self.n_angles=len(self.names)
        self.calls=0

    def _terminal_scores_and_adjoint(self,F):
        """Compute F_e and adjoint for the weighted trace-square objective.

        Input F has columns [U C|0>, U R C|0>]. alpha_jb is trace
        of four-column input block overlapped with encoded |b> output.
        A combined terminal adjoint W_j = w_j sum_b alpha_jb Y_b
        satisfies 2 Re <W, dF> /16 = d(F_e).
        """
        W=np.zeros_like(F)
        total=0.
        amplitudes=np.zeros((2,2),dtype=complex)
        for j in range(2):
            block=F[:,4*j:4*j+4]
            for b in (0,1):
                basis=self.readout[b]
                a=np.vdot(basis,block)
                amplitudes[j,b]=a
                total+=self.weights[j]*(abs(a)**2)/16
                W[:,4*j:4*j+4]+=self.weights[j]*a*basis
        return float(total),W,amplitudes

    def smooth_effort_and_gradient(self,angles):
        """Smooth, zero-at-identity angle effort; transport-aware proxy.

        Cost = Σ w_j [sqrt(theta_j²+eps²)-eps]. Unlike count of
        nonzero native operations, this function has a derivative.
        The weights are heuristics; they are not physical durations.
        """
        angles=np.asarray(angles,float)
        root=np.sqrt(angles*angles+self.effort_epsilon**2)
        quantity=float(np.dot(
            self.weighted_effort,root-self.effort_epsilon))
        gradient=self.weighted_effort*angles/root
        return quantity,gradient

    def objective_and_gradient(self,angles):
        self.calls+=1
        gates=self.base.matrices(angles)
        forward=[self.initial]
        for U in gates:
            forward.append(U@forward[-1])
        fidelity,back,amplitudes=self._terminal_scores_and_adjoint(
            forward[-1])
        grad=np.zeros(self.n_angles)
        for k in range(self.n_angles-1,-1,-1):
            H=self.base.H[k]
            derivative=-1j*(H@forward[k+1])
            grad[k]=-2*np.real(np.vdot(back,derivative))/16
            back=gates[k].conj().T@back
        effort,effort_grad=self.smooth_effort_and_gradient(angles)
        grad+=self.effort_lambda*effort_grad
        score=1-fidelity+self.effort_lambda*effort
        return float(score),grad

    def raw_fidelity(self,angles):
        final=self.base.propagate(angles,input_states=self.initial)
        fe,_,_=self._terminal_scores_and_adjoint(final)
        return fe

    def full_report(self,angles):
        U=self.base.propagate(
            angles,input_states=np.eye(TOTAL,dtype=complex))
        direct=self.raw_fidelity(angles)
        independent=single_fault_logical_recovery_fidelity(
            U,self.label,self.p_electric)
        expected=independent[
            "unconditional_seven_fault_mixture_logical_F_e"]
        if abs(direct-expected)>2e-10:
            raise AssertionError(
                f"two-block gradient objective != seven-Kraus recovery: {direct} vs {expected}")
        ref=self.base.diagnostics(angles)
        effort,_=self.smooth_effort_and_gradient(angles)
        return {
            "actual_unconditional_logical_entanglement_F_e":direct,
            "independent_all_seven_Kraus_recovery":independent,
            "do_nothing_baseline_logical_F_e":self.channel[
                "identity_reference_logical_F_e"],
            "unconditional_logical_improvement_over_identity":float(
                direct-self.channel["identity_reference_logical_F_e"]),
            "isometry_coherent_overlap":ref[
                "three_outcome_input_isometry_coherent_fidelity"],
            "entire_42D_extension_overlap":ref[
                "entire_42D_unitary_fidelity_to_selected_extension"],
            "per_branch":ref["per_branch"],
            "effort_proxy":effort,
            "regularized_objective":float(
                1-direct+self.effort_lambda*effort),
            "no_postselection":True,
        }


def direct_gradient_check(label="AS",cycles=1,seed=812,
                          p_electric=.37,effort_lambda=.0005,
                          epsilon=1e-6):
    c=DirectLogicalCompiler(label,cycles,p_electric,effort_lambda)
    angles=np.random.default_rng(seed).normal(0,.2,c.n_angles)
    val,g=c.objective_and_gradient(angles)
    rows=[]
    for k in (0,1,4,8,12):
        offset=np.zeros(c.n_angles)
        offset[k]=epsilon
        difference=(c.objective_and_gradient(angles+offset)[0]-
                    c.objective_and_gradient(angles-offset)[0])/(2*epsilon)
        rows.append({"index":k,
                     "adjoint_derivative":float(g[k]),
                     "finite_difference_derivative":float(difference),
                     "absolute_error":float(abs(g[k]-difference))})
    return {"group":label,"loss":val,
            "max_error":max(r["absolute_error"] for r in rows),
            "sampled_coordinate_derivatives":rows}


def compile_direct(label,cycles,maxiter=100,seed=66,
                   p_electric=.37,effort_lambda=.0005,
                   warmstart=None):
    """Train *logical fidelity*, never only 42D isometry overlap.

    Identity is always a valid candidate; numerical solver failure
    cannot lower the reported best unpostselected recovered F_e.
    """
    if not isinstance(maxiter,int) or maxiter<1:
        raise ValueError("maxiter must be positive integer")
    c=DirectLogicalCompiler(label,cycles,p_electric,effort_lambda)
    zero=np.zeros(c.n_angles)
    rng=np.random.default_rng(seed+cycles+(label=="BS")*97)
    starts=[("zero_identity",zero),
            ("seeded_small_random",rng.normal(0,.3,c.n_angles)),
            ("seeded_medium_random",rng.normal(0,.65,c.n_angles))]
    if warmstart is not None:
        prior=np.asarray(warmstart,dtype=float)
        if not np.all(np.isfinite(prior)) or not 0<len(prior)<c.n_angles:
            raise ValueError("previous depth must have fewer finite angles")
        starts.append(("previous_depth",np.pad(
            prior,(0,c.n_angles-len(prior)))))
    baseline=c.objective_and_gradient(zero)[0]
    best=(float(baseline),zero.copy(),"no_gate_identity")
    attempts=[]
    for name,start in starts:
        result=minimize(c.objective_and_gradient,start,method="L-BFGS-B",
                        jac=True,bounds=[(-np.pi,np.pi)]*c.n_angles,
                        options={"maxiter":maxiter,"maxls":20,
                                 "gtol":1e-7,"ftol":1e-10})
        val=float(c.objective_and_gradient(result.x)[0])
        attempts.append({
            "initialization":name,
            "final_regularized_loss":val,
            "optimizer_iterations":int(result.nit),
            "optimizer_success":bool(result.success),
            "optimizer_message":str(result.message),
        })
        if val<best[0]-1e-11:
            best=(val,result.x.copy(),name)
    loss,angles,winner=best
    outcome=c.full_report(angles)
    assert outcome["actual_unconditional_logical_entanglement_F_e"]>= (
        outcome["do_nothing_baseline_logical_F_e"]-2e-10)
    assert loss<=baseline+1e-10
    report={
        "group":label,
        "cycles":cycles,
        "controls_per_cycle":13,
        "pulse_count":len(angles),
        "angles_rad":angles.tolist(),
        "pulse_names":c.names,
        "selected_initialization":winner,
        "effort_lambda":effort_lambda,
        "p_electric":p_electric,
        "analytic_objective_evaluations":c.calls,
        "loss":loss,
        "attempts":attempts,
        "logical_fidelity":outcome,
        "identity_policy_retained_even_if_optimizer_fails":True,
    }
    report["native_primitive_cost"]=native_operation_expansion(report)
    duration=report["native_primitive_cost"][
        "hypothetical_serial_duration_ns_including_readout"]
    report["optimistic_four_active_fermion_30us_T1_weighted_logical_F_e_not_physical"]=float(
        outcome["actual_unconditional_logical_entanglement_F_e"]*
        np.exp(-4*duration/(1000*30.)))
    return report


def benchmark(cycles=(1,2,4,8),maxiter=100,labels=("AS","BS"),
              p_electric=.37,effort_lambda=.0005):
    if not cycles or any(type(i)!=int or i<1 for i in cycles):
        raise ValueError("positive depth cycles required")
    if tuple(cycles)!=tuple(sorted(set(cycles))):
        raise ValueError("cycles must be strictly increasing")
    reports=[]
    for label in labels:
        previous=None
        for depth in cycles:
            r=compile_direct(label,depth,maxiter=maxiter,
                p_electric=p_electric,effort_lambda=effort_lambda,
                warmstart=previous)
            reports.append(r)
            previous=r["angles_rad"]
    return {
        "experiment":"seven-single-fault UNCONDITIONAL logical entanglement-fidelity pulse optimization",
        "physical_hilbert_space":14,
        "neutral_pointer_levels":3,
        "reduced_clean_and_electric_input_columns":8,
        "original_single_step_Kraus_error_count":7,
        "equality_of_two_block_and_full_seven_error_channel":True,
        "p_electric":p_electric,
        "effort_lambda":effort_lambda,
        "gradient_check":direct_gradient_check(
            effort_lambda=effort_lambda),
        "depth_cycles":list(cycles),
        "benchmark":reports,
        "claim_limits":[
            "The optimizer's loss uses unconditional corrected logical F_e, not the former PR66 target-unitary overlap; identity is retained as a non-degrading baseline at EVERY tested depth.",
            "The two-block reduction applies because the six selected subgroup exchange unitaries each fix this encoded C; no two-fault gR/Rg or other gauge noise is contained.",
            "The smooth angle/charged-bus effort penalty is heuristic rather than calibrated physical duration; actual hypothetical primitive ledger is reported separately.",
            "Pointer-2 and outside-C outputs are genuine errors/flags, never postselected, and all seven individual Kraus alternatives are checked independently.",
            "Maximizing conditional seven-error F_e may yield an approximate or purely passive policy; the report does NOT presume the optimized circuit has learned genuine electric syndrome extraction.",
            "The default optimizer remains a cyclic 13-control gate template, not a provably minimal circuit or hardware-native implementation.",
            "All controls gauge-preserving only in the signed-permutation finite S3 model; opposite-vertex same-colour projector fails general continuous SU3 symmetry.",
            "Pulse lengths/rates are synthetic and neglect pointer, charged-carrier, reference, and spectator loss, electronic dephasing, compound faults during recovery, and virtual QCD fluctuations.",
            "No physical quantum-memory advantage or Standard-Model nine-quark 21D/65-syndrome code established.",
        ]
    }


def main():
    result=benchmark()
    Path("colour_s3_direct_logical_compiler_results.json").write_text(
        json.dumps(result,indent=2,sort_keys=True)+"\n")
    for r in result["benchmark"]:
        logic=r["logical_fidelity"]
        native=r["native_primitive_cost"]
        print("DIRECT_LOGICAL_FIDELITY_PULSE",json.dumps({
            "group":r["group"],
            "pulses":r["pulse_count"],
            "unconditional_F_e":logic["actual_unconditional_logical_entanglement_F_e"],
            "do_nothing_F_e":logic["do_nothing_baseline_logical_F_e"],
            "gain_over_identity":logic["unconditional_logical_improvement_over_identity"],
            "isometry_overlap":logic["isometry_coherent_overlap"],
            "native_expanded_count":native["two_register_native_or_abstract_primitive_count"],
            "synthetic_serial_time_ns":native["hypothetical_serial_duration_ns_including_readout"],
            "best_initialization":r["selected_initialization"],
        },sort_keys=True),flush=True)


if __name__=="__main__":
    main()
