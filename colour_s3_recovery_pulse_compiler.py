"""Target-specific, reproducible numerical pulse compilation of PR64 recovery.

This is the first FINITE executable pulse schedule attempt (not merely
a Lie-controllability certificate) for the finite-S3 14D Gauss sector
plus a neutral three-level pointer (dimension 42).

The pulse library is precisely PR65's independently controlled
Gauss-compatible Hamiltonians:
  onsite swap, dressed edge swap, link standard-electric projector,
  opposite-vertex equal-colour projector, 8 neutral-pointer su3 drives,
  and pointer/link controlled-electric entangling interaction.

Optimize *one chosen extension* of PR64's three-outcome instrument
on the physically relevant input subspace H_Gauss x |pointer=0>:
clean C -> C/0; electric V -> C/1; six-dimensional complement
-> itself/2. All three branches are scored with a COMMON coherent
global phase, not cherry-picked or postselected. Also independently
report residual against the entire 42x42 unitary extension, which is
NOT optimized in this mode.

Pulse angles are genuine stored reproducible numbers. Analytic
adjoint gradient is independently checked against finite differences.
Schedules are expanded for every nontrivial dressed-edge exchange
using PR65's EXACT five charged-mediator two-register pulse identity.
There is no native compilation of the charged transport or the
pointer-link projector coupling on a specified physical device.

NO success threshold is assumed: report achieved error as-is, even
when poor, across bounded depths. The ideal SU(42) controllability
certificate from PR65 does NOT imply a short sequence exists.
"""
from __future__ import annotations

from pathlib import Path
import json
import numpy as np
from scipy.optimize import minimize

from colour_s3_exact_electric_commutant import exact_generators
from colour_s3_order_six_active_electric import (
    closure_matrices,construct_active_code,
    gauge_fixed_normalized_matrix,
)
from colour_s3_local_control_universality import admissible_controls
from colour_s3_correlated_recovery_cadence import (
    global_gauge_preserving_instrument,
)

DIM=14
POINTER=3
TOTAL=DIM*POINTER
NINPUT=DIM
LABELS=("AS","BS")


def pointer_generators():
    """Explicit eight Hermitian, traceless one-qutrit drives."""
    out=[]
    for i,j in ((0,1),(0,2),(1,2)):
        x=np.zeros((3,3),complex)
        y=np.zeros((3,3),complex)
        x[i,j]=x[j,i]=.5
        y[i,j]=-.5j
        y[j,i]=.5j
        out.extend([(f"pointer_X{i}{j}",x),
                    (f"pointer_Y{i}{j}",y)])
    out.extend([
        ("pointer_Z01",np.diag([.5,-.5,0.])),
        ("pointer_Z12",np.diag([0.,.5,-.5])),
    ])
    assert len(out)==8
    for _,m in out:
        assert np.linalg.norm(m-m.conj().T)<1e-14
        assert abs(np.trace(m))<1e-14
    return out


def controls_for_group(label):
    """13 addressable 42D Hamiltonians, each locally gauge compatible.

    Coordinates are the normalized 14D orbit-sum basis, with order
    SYSTEM tensor POINTER in all matrix Kronecker products.
    """
    classes,_,ops=exact_generators()
    if label not in LABELS:raise ValueError(label)
    sizes=[len(c) for c in classes]
    raw,inventory=admissible_controls(classes,ops,label)
    name_pairs=[]
    for (name,_),matrix in zip(inventory,raw):
        physical=gauge_fixed_normalized_matrix(matrix,sizes)
        name_pairs.append((name,np.kron(physical,np.eye(POINTER))))
    for name,ptr in pointer_generators():
        name_pairs.append((name,np.kron(np.eye(DIM),ptr)))
    E=gauge_fixed_normalized_matrix(ops[-1],sizes)
    pp1=np.diag([0.,1.,0.])
    name_pairs.append(("pointer_link_electric_entangler",
                       np.kron(E,pp1)))
    assert len(name_pairs)==13
    for name,H in name_pairs:
        if np.linalg.norm(H-H.conj().T)>5e-12:
            raise AssertionError(f"{name} non-Hermitian")
    return classes,ops,sizes,name_pairs


def target_three_outcome(label):
    """Construct full PR64 U42 and 14-dimensional pointer-ready target.

    Preserve the 6D flagged-complement basis; an isometry-target
    objective acts on EVERY one of the 14 physical Gauss directions.
    """
    classes,ops,sizes,named=controls_for_group(label)
    A,B,S,E=ops
    grp=closure_matrices((A,S) if label=="AS" else (B,S))
    C,V,_=construct_active_code(label,grp,E,sizes)
    target_reference=global_gauge_preserving_instrument(C,V)
    assert target_reference["pointer_levels"]==3
    I=np.eye(DIM,dtype=complex)
    P=C@C.conj().T
    Q=V@V.conj().T
    T=I-P-Q
    eye=np.eye(3,dtype=complex)
    swap01=eye[[1,0,2]]
    swap02=eye[[2,1,0]]
    F=C@V.conj().T+V@C.conj().T+T
    synd=(np.kron(P,eye)+np.kron(Q,swap01)+np.kron(T,swap02))
    feedback=(np.kron(I,np.diag([1.,0.,0.]))+
              np.kron(F,np.diag([0.,1.,0.]))+
              np.kron(I,np.diag([0.,0.,1.])))
    Uideal=feedback@synd
    assert np.linalg.norm(Uideal.conj().T@Uideal-np.eye(TOTAL))<5e-11
    # tensor order 14D system by 3D pointer; pointer initialized |0>
    initial=np.kron(np.eye(DIM),np.eye(POINTER)[:,[0]])
    goal=Uideal@initial
    assert goal.shape==(TOTAL,NINPUT)
    eigvec,eig=np.linalg.eigh(T)
    complement=eig[:,eigvec>.5]
    assert complement.shape==(DIM,6)
    branches={
        "code":C,
        "electric":V,
        "flag":complement,
    }
    return named,Uideal,initial,goal,branches


def single_fault_logical_recovery_fidelity(U,label,p_electric=.37):
    """Unconditional recovered logical entanglement fidelity, NO postselection.

    For each of the six admitted matter-exchange errors g and electric
    R=I-2E from PR63, initialize code state C and neutral pointer |0>,
    apply the EXACT fault operator, then the ACTUAL compiled U42.
    Measure pointer and count only output in C with pointer 0 OR 1 as
    recovered logical information. Pointer 2 and leakage outside C
    are orthogonal flagged failure, never renormalized away.

    This gives true logical F_e for the trace-decreasing recovered
    subchannel with erasures flagged to orthogonal states. Use a
    separately specified synthetic 7-Kraus mixture: R with p=.37,
    each of six g with (1-p)/6. It does not model double faults,
    hardware leakage, T1/T2 or measurement assignment errors.
    """
    if U.shape!=(TOTAL,TOTAL):
        raise ValueError("Expected physical-system×pointer 42D unitary")
    if not np.isfinite(p_electric) or not 0<=p_electric<=1:
        raise ValueError("p_electric must be in [0,1]")
    classes,ops,sizes,_=controls_for_group(label)
    A,B,S,E=ops
    subgroup=closure_matrices((A,S) if label=="AS" else (B,S))
    C,V,_=construct_active_code(label,subgroup,E,sizes)
    Ephys=gauge_fixed_normalized_matrix(E,sizes)
    faults=[(f"exchange_{i}",gauge_fixed_normalized_matrix(g,sizes),
             (1-p_electric)/6) for i,g in enumerate(subgroup)]
    faults.append(("electric_phase",np.eye(DIM)-2*Ephys,p_electric))
    pointer0=np.eye(POINTER)[:,[0]]
    P=C@C.conj().T
    records=[]
    weighted=0.
    for name,fault,weight in faults:
        evolved=U@np.kron(fault@C,pointer0)
        split=evolved.reshape(DIM,POINTER,4)
        Kraus=[C.conj().T@split[:,i,:] for i in (0,1)]
        recovered=float(sum(abs(np.trace(k)/4)**2 for k in Kraus))
        captured=float(sum(np.linalg.norm(k,"fro")**2/4 for k in Kraus))
        flag=float(np.linalg.norm(split[:,2,:],"fro")**2/4)
        # Unflagged output outside C in pointer 0 or 1 is also a
        # logical failure; preserve it in the unconditional ledger.
        missed=1-captured-flag
        assert -1e-10<=missed<=1+1e-10,(name,missed)
        assert -1e-10<=recovered<=captured+1e-10
        record={
            "error":name,"synthetic_fault_probability":weight,
            "unconditional_recovered_logical_entanglement_fidelity":
                recovered,
            "logical_code_capture_probability":captured,
            "pointer_2_flagged_probability":flag,
            "unflagged_outside_logical_code_probability":float(
                max(0.,missed)),
        }
        weighted+=weight*recovered
        records.append(record)
    assert abs(sum(r["synthetic_fault_probability"] for r in records)-1)<1e-12
    return {
        "fault_model":"R probability p, each of six permitted individual exchange unitaries probability (1-p)/6; NO compound faults",
        "electric_probability":p_electric,
        "unconditional_seven_fault_mixture_logical_F_e":float(weighted),
        "worst_individual_fault_logical_F_e":min(x[
            "unconditional_recovered_logical_entanglement_fidelity"]
            for x in records),
        "individual_faults":records,
        "pointer_2_or_outside_code_events_not_postselected":True,
        "unsimulated_hardware_ancilla_noise":True,
    }


def pulse_pattern(n_controls,cycles):
    if not isinstance(cycles,int) or cycles<1:
        raise ValueError("cycles must be positive integer")
    return tuple(range(n_controls))*cycles


class SpectralCompiler:
    """Exact 42D gate exponentials plus analytic reverse-mode gradient."""

    def __init__(self,label,cycles,full_unitary=False):
        named,Utarget,initial,goal,branches=target_three_outcome(label)
        self.label=label
        self.cycles=cycles
        self.named=named
        self.target_unitary=Utarget
        self.initial=initial if not full_unitary else np.eye(TOTAL)
        self.target=goal if not full_unitary else Utarget
        self.branches=branches
        self.full_unitary=bool(full_unitary)
        self.order=pulse_pattern(len(named),cycles)
        self.dim=self.initial.shape[1]
        self.H=[named[i][1] for i in self.order]
        self.eigen=[]
        for _,H in named:
            evals,v=np.linalg.eigh(H)
            self.eigen.append((evals,v))
        self.eps=np.finfo(float).eps
        self.objective_calls=0

    def matrices(self,angles):
        angles=np.asarray(angles,dtype=float)
        if angles.shape!=(len(self.order),) or not np.all(np.isfinite(angles)):
            raise ValueError("Pulse angles must be finite with proper schedule length")
        out=[]
        for index,angle in zip(self.order,angles):
            eig,v=self.eigen[index]
            U=(v*np.exp(-1j*angle*eig))@v.conj().T
            out.append(U)
        return out

    def propagate(self,angles,input_states=None):
        operators=self.matrices(angles)
        vec=self.initial.copy() if input_states is None else input_states.copy()
        for U in operators:
            vec=U@vec
        return vec

    def objective_and_gradient(self,angles):
        """J=1-|Tr(T† U_angles Phi0)|²/d²; no per-branch rephasing.

        Reverse-mode derivative of U_k=exp(-i a_k H_k) is
        dU_k/da_k = -i H_k U_k, since it commutes with H_k.
        This reduces optimization to cheap 42x14 propagation and
        avoids differentiating the 42x42 matrix exponential numerically.
        """
        self.objective_calls+=1
        gates=self.matrices(angles)
        forwards=[self.initial]
        for U in gates:
            forwards.append(U@forwards[-1])
        amplitude=np.vdot(self.target,forwards[-1])
        fidelity=float(abs(amplitude/self.dim)**2)
        back=self.target.copy()
        grad=np.zeros(len(self.order))
        for k in range(len(self.order)-1,-1,-1):
            derivative=-1j*(self.H[k]@forwards[k+1])
            d_amp=np.vdot(back,derivative)
            grad[k]=-2*np.real(np.conj(amplitude)*d_amp)/(self.dim**2)
            back=gates[k].conj().T@back
        return float(1-fidelity),grad

    def objective(self,angles):
        return self.objective_and_gradient(angles)[0]

    def diagnostics(self,angles):
        U=self.propagate(angles,input_states=np.eye(TOTAL,dtype=complex))
        raw=self.initial
        achieved=U@raw
        amplitude=np.vdot(self.target,achieved)
        isometry_f=float(abs(amplitude/NINPUT)**2)
        entire_f=float(abs(np.trace(
            self.target_unitary.conj().T@U)/TOTAL)**2)
        # Inspect each physically relevant input branch independently.
        ptr0=np.eye(POINTER)[:,[0]]
        cases={}
        goals={"code":0,"electric":1,"flag":2}
        for name,basis in self.branches.items():
            source=np.kron(basis,ptr0)
            out=U@source
            target=self.target_unitary@source
            rank=basis.shape[1]
            ptr=(abs(out.reshape(DIM,POINTER,rank))**2).sum(axis=(0,2))/rank
            target_overlap=float(abs(np.vdot(target,out)/rank)**2)
            cases[name]={
                "input_subspace_rank":rank,
                "intended_pointer_label":goals[name],
                "intended_pointer_probability":float(ptr[goals[name]]),
                "other_pointer_probability":float(1-ptr[goals[name]]),
                "coherent_target_branch_overlap_squared":target_overlap,
                "raw_target_path_frobenius_residual":float(
                    np.linalg.norm(target-out)/np.sqrt(rank)),
                "pointer_outcome_probabilities":ptr.tolist(),
            }
        return {
            "three_outcome_input_isometry_coherent_fidelity":isometry_f,
            "entire_42D_unitary_fidelity_to_selected_extension":entire_f,
            "three_outcome_input_isometry_error":float(1-isometry_f),
            "per_branch":cases,
            "physical_global_phase_irrelevant":True,
            "no_postselection_used":True,
            "synthetic_seven_single_fault_logical_recovery":
                single_fault_logical_recovery_fidelity(U,self.label),
        }


def compile_depth(label,cycles,maxiter=60,seed=20261019,
                  warmstart=None,full_unitary=False):
    """Bounded L-BFGS-B with deterministic random and zero starts.

    Include zero-angle identity and (where supplied) warm-start from
    smaller-depth schedule. Select best FINAL result, including the
    zero baseline; never report a worse optimized value as improvement.
    """
    if not isinstance(maxiter,int) or maxiter<1:
        raise ValueError("maxiter must be a positive integer")
    compiler=SpectralCompiler(label,cycles,full_unitary)
    baseline=np.zeros(len(compiler.order))
    rng=np.random.default_rng(seed+cycles+int(label=="BS")*101)
    starts=[("identity",baseline),
            ("seeded_small_random",rng.normal(
                0,.15,len(compiler.order)))]
    if warmstart is not None:
        warm=np.asarray(warmstart,float)
        if len(warm)>=len(baseline):
            raise ValueError("warmstart length must be smaller than new circuit")
        starts.append(("smaller_depth_warmstart",np.pad(
            warm,(0,len(baseline)-len(warm)))))
    initial_score=float(compiler.objective(baseline))
    best=(initial_score,baseline,"identity_unoptimized")
    attempts=[]
    for name,start in starts:
        result=minimize(
            compiler.objective_and_gradient,start,
            method="L-BFGS-B",jac=True,
            bounds=[(-np.pi,np.pi)]*len(start),
            options={"maxiter":maxiter,"ftol":1e-10,"gtol":1e-7,
                     "maxls":20},
        )
        score=float(compiler.objective(result.x))
        attempts.append({
            "initialization":name,
            "objective":score,
            "iterations":int(result.nit),
            "function_evaluations":int(result.nfev),
            "optimizer_success":bool(result.success),
            "optimizer_message":str(result.message),
        })
        if score<best[0]-1e-12:
            best=(score,result.x.copy(),name)
    err,angles,winner=best
    report=compiler.diagnostics(angles)
    coherent=report["three_outcome_input_isometry_coherent_fidelity"]
    if not full_unitary:
        assert abs(err-(1-coherent))<1e-10
    return {
        "label":label,"cycles":cycles,
        "pulse_count":len(compiler.order),
        "angles_rad":[float(x) for x in angles],
        "pulse_names":[compiler.named[i][0] for i in compiler.order],
        "fidelity_diagnostics":report,
        "zero_angle_identity_baseline_fidelity":float(1-initial_score),
        "improved_over_identity":bool(err<initial_score-1e-8),
        "winning_initialization":winner,
        "attempts":attempts,
        "objective_calls_total":compiler.objective_calls,
        "objective_mode":"full_selected_42D_unitary" if full_unitary
            else "coherent_14D_input_three_outcome_isometry",
        "optimization_converged_to_perfect_target":bool(err<1e-8),
        "explicit_compiled_schedule_exists":True,
        "success_as_accurate_recovery_not_assumed":True,
    }


def native_operation_expansion(report,exchange_mhz=10.,
                               onsite_mhz=10.,electric_mhz=10.,
                               equality_mhz=10.,pointer_mhz=10.,
                               entangler_mhz=10.,transport_ns=100.,
                               readout_ns=300.,skip_below=1e-9):
    """Expand every dressed S into 5 elementary gauge-preserving ops.

    Convert Hamiltonian-angle magnitude to hypothetical execution time
    via t[ns] = |theta|*1000/(2pi*f[MHz]). Exchange dressed SWAP S
    phase is equivalent up to scalar global phase to an involution
    exchange rotation with angle (-2theta), whose coloured-bus
    implementation costs four conjugation primitives plus 1 local
    fractional mediator/matter exchange.

    Finite S3 charged mediator transport itself is still an abstract
    native gate. No reset/transport scheduling overlap is assumed.
    """
    rates=[exchange_mhz,onsite_mhz,electric_mhz,equality_mhz,
           pointer_mhz,entangler_mhz]
    if not all(np.isfinite(v) and v>0 for v in rates):
        raise ValueError("all drive rates must be positive finite")
    if not np.isfinite(transport_ns) or transport_ns<0:
        raise ValueError("transport_ns must be nonnegative finite")
    if not np.isfinite(readout_ns) or readout_ns<0:
        raise ValueError("readout_ns must be nonnegative finite")
    steps=[]
    def add(operation,angle,ns,category):
        steps.append({"name":operation,"rotation_rad":float(angle),
                      "duration_ns":float(ns),"category":category,
                      "register_support_upper_bound":2})
    for name,theta in zip(report["pulse_names"],report["angles_rad"]):
        if abs(theta)<=skip_below:
            continue
        if name.startswith("dressed_swap_"):
            # Full SWAP(m,A0) = exp(-i pi (I-SWAP)/2).
            fixed_ns=1000*np.pi/(2*np.pi*exchange_mhz)
            frac=-2*theta
            frac_ns=1000*abs(frac)/(2*np.pi*exchange_mhz)
            add("charged_mediator_swap_A0",np.pi,fixed_ns,"charged_bus")
            add("charged_transport_A_to_B",0.,transport_ns,"transport")
            add("charged_mediator_fractional_swap_B0",frac,
                frac_ns,"charged_bus")
            add("charged_transport_B_to_A",0.,transport_ns,"transport")
            add("charged_mediator_swap_A0",np.pi,fixed_ns,"charged_bus")
        else:
            if name.startswith("onsite_swap_"):
                speed=onsite_mhz
                # exp(-i theta SWAP) matches exchange-projector angle
                # -2theta up to physically irrelevant global phase.
                physical_angle=-2*theta
            elif name.startswith("electric_link_"):
                speed=electric_mhz;physical_angle=theta
            elif name.startswith("same_colour_"):
                speed=equality_mhz;physical_angle=theta
            elif name.startswith("pointer_link_"):
                speed=entangler_mhz;physical_angle=theta
            elif name.startswith("pointer_"):
                speed=pointer_mhz;physical_angle=theta
            else:raise ValueError("Unknown control name "+name)
            add(name,physical_angle,
                1000*abs(physical_angle)/(2*np.pi*speed),
                "pointer" if name.startswith("pointer_") else "system")
    native_count=len(steps)
    ntransport=sum(s["category"]=="transport" for s in steps)
    nbus=sum(s["category"]=="charged_bus" for s in steps)
    duration=sum(s["duration_ns"] for s in steps)+readout_ns
    return {
        "compiled_continuous_control_pulse_count":report["pulse_count"],
        "nonzero_angle_control_pulses":sum(
            abs(t)>skip_below for t in report["angles_rad"]),
        "two_register_native_or_abstract_primitive_count":native_count,
        "charged_bus_swap_and_fractional_steps":nbus,
        "charged_link_transport_steps":ntransport,
        "hypothetical_serial_duration_ns_including_readout":duration,
        "hypothetical_readout_ns":readout_ns,
        "all_coupling_rates_MHz_assumed":True,
        "all_charged_link_transports_are_abstract_not_device_native":True,
        "assumed_no_parallel_gate_execution":True,
        "expanded_primitive_schedule":steps,
    }


def finite_gradient_error(label="AS",cycles=1,seed=14,epsilon=1e-6):
    compiler=SpectralCompiler(label,cycles)
    angles=np.random.default_rng(seed).normal(0,.2,len(compiler.order))
    value,g=compiler.objective_and_gradient(angles)
    checks=[]
    for index in (0,3,4,9,len(angles)-1):
        changed=np.zeros(len(angles));changed[index]=epsilon
        fd=(compiler.objective(angles+changed)-
            compiler.objective(angles-changed))/(2*epsilon)
        checks.append({"index":index,"analytic":float(g[index]),
                       "finite_difference":float(fd),
                       "absolute_error":float(abs(g[index]-fd))})
    return {"objective":value,"gradient_checks":checks,
            "max_gradient_error":max(x["absolute_error"] for x in checks)}


def study(cycles=(1,2,4,8,16),maxiter=80,labels=LABELS):
    if not cycles or any(not isinstance(x,int) or x<1 for x in cycles):
        raise ValueError("positive cycles required")
    if tuple(sorted(cycles))!=tuple(cycles):
        raise ValueError("cycles must be nondecreasing")
    results=[]
    for label in labels:
        previous=None
        for depth in cycles:
            result=compile_depth(label,depth,maxiter=maxiter,
                                 warmstart=previous)
            expansion=native_operation_expansion(result)
            angles=np.asarray(result["angles_rad"])
            # Common coherent fractional overrotation on EVERY angle:
            # diagnostic only, NOT random independent gate-noise.
            compiler=SpectralCompiler(label,depth)
            errors={}
            for delta in (-.01,.01):
                diags=compiler.diagnostics(angles*(1+delta))
                errors[str(delta)]=diags[
                    "three_outcome_input_isometry_coherent_fidelity"]
            result["uniform_coherent_angle_miscalibration"]=errors
            result["hypothetical_native_two_register_expansion"]=expansion
            result["optimistic_four_active_particle_T1_30us_survival_scaled_isometry_score"]=(
                result["fidelity_diagnostics"][
                    "three_outcome_input_isometry_coherent_fidelity"]*
                np.exp(-4*expansion[
                    "hypothetical_serial_duration_ns_including_readout"]/
                    (1000*30.)))
            result["no_device_performance_claim"]=True
            results.append(result)
            previous=angles
    return {
        "research_target":"PR64 three-outcome neutral-qutrit Gauss-invariant coherent instrument",
        "physical_hilbert_dimension":DIM,
        "pointer_dimension":POINTER,
        "target_unitary_dimension":TOTAL,
        "optimized_input_isometry_columns":NINPUT,
        "system_control_alphabet":"four PR65 gauge-preserving controls, eight neutral-qutrit drives, pointer/link electric entangler",
        "depth_cycles":list(cycles),"iterations_cap_each_start":maxiter,
        "two_physical_gauge_exchange_subgroups":list(labels),
        "gradient_audit":finite_gradient_error(),
        "benchmark":results,
        "limitations":[
            "This is a FINITE reproducible numerical schedule, but optimized objective is only the 14D system × initial-pointer-0 ISOMETRY; report the separate full 42D unitary overlap rather than claiming the whole unitary was compiled.",
            "Input-isometry coherent overlap is NOT the same as logical recovery: a separate unpostselected seven-operator synthetic channel F_e is computed by projecting both recognized pointer branches back into C and treating all pointer-2/outside-code outcomes as orthogonal failures.",
            "The seven-error logical recovery channel uses electric p=.37 and uniform six-element exchange alternatives, not a device-derived noise distribution or multiple correlated electric-exchange faults.",
            "A low benchmark fidelity means the approximate recovery is NOT usable; Lie universal reachability does not guarantee convergence at short pulse depth.",
            "The expansion realizes dressed-S fractionals with five abstract charged-mediator two-register operations; mobile charge transport and local exchange are assumed coherent, tunable primitives, NOT derived device pulses.",
            "The opposite-vertex equal-colour local control is gauge-invariant only under finite signed S3, NOT continuous SU3.",
            "All GHz/MHz rates and gate timing inputs are hypothetical. Charge/reference preparation, spectator/reference/pointer T1, correlated faults during the compiled circuit, reset and readout errors are omitted.",
            "The T1-scaled isometry score is a conditional, optimistic diagnostic, NOT the physical QEC logical entanglement fidelity or evidence of coding advantage.",
            "The optimization uses bounded angles, fixed cyclic gate order, a small number of starts and finite iterations. It is not a globally optimal depth lower bound.",
            "No virtual gluon, sea-quark or nonperturbative QCD field degrees of freedom; this is the four matter-qutrit finite S3 toy model, not the nine-quark 21D/65-syndrome code.",
        ]
    }


def main():
    data=study()
    out=Path("colour_s3_recovery_pulse_compiler_results.json")
    out.write_text(json.dumps(data,indent=2,sort_keys=True)+"\n")
    for item in data["benchmark"]:
        print("RECOVERY_NUMERICAL_PULSE_DEPTH",json.dumps({
            "group":item["label"],"cycles":item["cycles"],
            "control_pulses":item["pulse_count"],
            "fidelity":item["fidelity_diagnostics"][
                "three_outcome_input_isometry_coherent_fidelity"],
            "whole_42D_fidelity":item["fidelity_diagnostics"][
                "entire_42D_unitary_fidelity_to_selected_extension"],
            "synthetic_seven_fault_unconditional_logical_F_e":item["fidelity_diagnostics"]["synthetic_seven_single_fault_logical_recovery"]["unconditional_seven_fault_mixture_logical_F_e"],
            "branch_scores":{k:v["coherent_target_branch_overlap_squared"]
                for k,v in item["fidelity_diagnostics"]["per_branch"].items()},
            "native_expanded_count":item[
                "hypothetical_native_two_register_expansion"][
                    "two_register_native_or_abstract_primitive_count"],
            "serial_assumed_ns":item[
                "hypothetical_native_two_register_expansion"][
                    "hypothetical_serial_duration_ns_including_readout"],
            "ideal_target_reached":item["optimization_converged_to_perfect_target"],
        },sort_keys=True),flush=True)


if __name__=="__main__":
    main()
