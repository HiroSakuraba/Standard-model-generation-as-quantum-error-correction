"""Four-level pointer + finite-S3 effective open-system recovery audit.

CPTP local pointer bosonic loss and Gaussian number dephasing are imported
from the existing four-level transmon probe. A separate gauge-compatible
14D ELECTRIC-SYNDROME phase-dephasing channel and independently flagged
four-active-particle survival are explicit model assumptions. The 14D
Gauss sector is NOT a physical transmon tensor product, so this module
MUST NOT call system erasures 'ququart amplitude damping' or invent
a mapping to an underlying physical CP22 device.

The exact 42D numerical pulse product is embedded in a 14x4=56 space:
pointer |3> is an additional leakage level. Channels are applied between
the two equal-pulse-count halves and after the second half; this is
operator splitting, NOT microscopic gate-time Lindblad evolution.

All Kraus outcomes are summed, NEVER postselected. Same-time passive
and immediate no-operation controls are evaluated with the same fault
inputs and particle loss model. Logical output overlap is an observable
(not an executable free code projector or physical decoder).

This is a bounded feasibility proxy, not device-calibrated QEC or QCD.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path
import json
import numpy as np

from su3_transmon_noise_probe import loss_kraus, dephasing_kraus
from colour_s3_recovery_pulse_compiler import (
    DIM, TOTAL, SpectralCompiler, target_three_outcome,
    native_operation_expansion, controls_for_group,
)
from colour_s3_order_six_active_electric import (
    gauge_fixed_normalized_matrix,
)
from colour_s3_continuation_search import run_continuation

P4 = 4
FULL = DIM * P4
IDENTITY_SYSTEM = np.eye(DIM, dtype=complex)
IDENTITY_POINTER = np.eye(P4, dtype=complex)
COMPUTATIONAL_INDICES = np.array(
    [P4*i+j for i in range(DIM) for j in range(3)], int)


@dataclass(frozen=True)
class NoiseScenario:
    """All times are ns; every numerical value is illustrative."""
    electric_probability: float = .37
    pointer_T1_ns: float = 30000.
    pointer_Tphi_ns: float = 45000.
    electric_Tphi_ns: float = 45000.
    active_particle_T1_ns: float = 30000.
    pointer_leak_angle_per_window: float = .01
    readout_wrong_label_probability: float = .005
    readout_ns: float = 300.
    reset_ns: float = 300.

    def validate(self):
        values = (
            self.electric_probability, self.pointer_T1_ns,
            self.pointer_Tphi_ns, self.electric_Tphi_ns,
            self.active_particle_T1_ns,
            self.pointer_leak_angle_per_window,
            self.readout_wrong_label_probability,
            self.readout_ns, self.reset_ns,
        )
        if not all(np.isfinite(v) for v in values):
            raise ValueError("all parameters must be finite")
        if (not 0 <= self.electric_probability <= 1 or
            not 0 <= self.readout_wrong_label_probability <= 1 or
            any(x <= 0 for x in (
                self.pointer_T1_ns, self.pointer_Tphi_ns,
                self.electric_Tphi_ns, self.active_particle_T1_ns)) or
            any(x < 0 for x in (
                self.pointer_leak_angle_per_window,
                self.readout_ns, self.reset_ns))):
            raise ValueError("probabilities, T1/Tphi, angle or timing invalid")
        return self


def embed_pointer_ququart(u42):
    """Acts as U42 on levels 0,1,2, identity on pointer |3>."""
    u42 = np.asarray(u42, complex)
    if u42.shape != (TOTAL, TOTAL):
        raise ValueError("expected 42x42 control")
    full = np.eye(FULL, dtype=complex)
    full[np.ix_(COMPUTATIONAL_INDICES, COMPUTATIONAL_INDICES)] = u42
    return full


def pointer_leakage_unitary(angle):
    """Toy coherent |2> <-> |3> pointer mixing, not a CP22 transmon pair."""
    if not np.isfinite(angle) or angle < 0:
        raise ValueError("invalid leakage rotation")
    u = IDENTITY_POINTER.copy()
    c = np.cos(angle)
    s = -1j * np.sin(angle)
    u[2,2] = u[3,3] = c
    u[2,3] = u[3,2] = s
    return u


def apply_pointer(state, k):
    """Operate on the pointer while preserving system + logical reference."""
    state = state.reshape(DIM, P4, 4)
    return np.einsum("pq,aqr->apr", k, state,
                     optimize=True).reshape(FULL,4)


def apply_system(state, k):
    state = state.reshape(DIM, P4, 4)
    return np.einsum("ab,bpr->apr", k, state,
                     optimize=True).reshape(FULL,4)


def electric_phase_kraus(dt, tphi, R):
    """Gauge-compatible Z_R random-phase dephasing; NOT site-level Tphi."""
    if not np.isfinite(dt) or dt < 0 or not np.isfinite(tphi) or tphi <= 0:
        raise ValueError("bad electric dephasing time")
    decay = np.exp(-dt/tphi)
    return [np.sqrt((1+decay)/2)*IDENTITY_SYSTEM,
            np.sqrt((1-decay)/2)*R]


def _nonzero_ops(ops):
    return [op for op in ops if np.linalg.norm(op, "fro") > 1e-13]


def channel_branch(states, duration_ns, sc, R):
    """Exact sum of Kraus branches for this coarse time window.

    For dt>0: system Z_R dephasing followed by pointer pure loss followed
    by pointer phase diffusion. Kraus completeness makes this CPTP.
    Branch state norms are intentionally *not* renormalized.
    """
    if not np.isfinite(duration_ns) or duration_ns < 0:
        raise ValueError("invalid duration")
    if duration_ns == 0:
        return list(states)
    sys = _nonzero_ops(electric_phase_kraus(
        duration_ns, sc.electric_Tphi_ns, R))
    ptrloss = _nonzero_ops(loss_kraus(duration_ns, sc.pointer_T1_ns))
    ptrphase = _nonzero_ops(dephasing_kraus(
        duration_ns, sc.pointer_Tphi_ns))
    result = []
    for y in states:
        for s in sys:
            v = apply_system(y, s)
            for l in ptrloss:
                lv = apply_pointer(v, l)
                for ph in ptrphase:
                    z = apply_pointer(lv, ph)
                    if np.linalg.norm(z, "fro") > 1e-13:
                        result.append(z)
    return result


def _fault_state(code, fault):
    return np.kron(fault@code, IDENTITY_POINTER[:, [0]])


def _readout_acceptance(pointer_label, error):
    """Symmetric wrong-label confusion between three reported labels.

    The physical fourth level has no valid report and always erases.
    """
    if pointer_label == 3:
        return 0.
    if pointer_label in (0,1):
        return 1-error/2
    return error


def _kraus_logical_ledger(states, code, error):
    score = 0.
    code_capture = 0.
    leak3 = 0.
    total = 0.
    branch_p = np.zeros(P4)
    for state in states:
        split = state.reshape(DIM,P4,4)
        total += np.linalg.norm(state, "fro")**2/4
        for b in range(P4):
            v = split[:,b,:]
            p = np.linalg.norm(v, "fro")**2/4
            branch_p[b] += p
            if b==3:
                leak3 += p
            a = code.conj().T@v
            accept = _readout_acceptance(b,error)
            score += accept * abs(np.trace(a)/4)**2
            code_capture += accept * np.linalg.norm(a,"fro")**2/4
    if abs(total-1.)>2e-8:
        raise AssertionError(f"Kraus sum not trace preserving: {total}")
    if score > code_capture+2e-9:
        raise AssertionError("fidelity exceeds recognized code probability")
    return {
        "F_e_before_active_particle_erasure":float(score),
        "accepted_code_probability":float(code_capture),
        "fourth_level_pointer_probability":float(leak3),
        "physical_pointer_outcome_probabilities":[float(x) for x in branch_p],
        "total_Kraus_trace":float(total),
        "Kraus_branches":len(states),
    }


def run_instrument(label, pulses, pulse_ns, scenario=NoiseScenario(),
                   nominal_pulse_count=None):
    """Exact two-window CPTP pointer model and an erasure survival ledger.

    'pulses' are the ACTUAL 42D unitary control steps, in chronological
    order, not an idealized synthesized measurement. For an ideal-global
    mathematical upper bound pass [U_ideal] and mark it accordingly.
    """
    sc = scenario.validate()
    if not np.isfinite(pulse_ns) or pulse_ns < 0:
        raise ValueError("pulse_ns must be nonnegative finite")
    if nominal_pulse_count is not None and nominal_pulse_count != len(pulses):
        raise ValueError("pulse count mismatch")
    classes,ops,sizes,_=controls_for_group(label)
    E=gauge_fixed_normalized_matrix(ops[-1],sizes)
    R=IDENTITY_SYSTEM-2*E
    _,_,_,_,branches=target_three_outcome(label)
    code=branches["code"]
    if np.linalg.norm(R.conj().T@R-IDENTITY_SYSTEM)>2e-9:
        raise AssertionError("electric reflection must be unitary")

    mid=len(pulses)//2
    first=[embed_pointer_ququart(u) for u in pulses[:mid]]
    second=[embed_pointer_ququart(u) for u in pulses[mid:]]
    leak=pointer_leakage_unitary(
        sc.pointer_leak_angle_per_window if pulses else 0.)
    fault_rows=[]
    duration=pulse_ns+sc.readout_ns+sc.reset_ns
    for fault_name, fault, weight in (
        ("six_exchange_equivalent_on_code", IDENTITY_SYSTEM,
         1-sc.electric_probability),
        ("electric_reflection", R, sc.electric_probability),
    ):
        input_state=_fault_state(code,fault)
        y=input_state
        for u in first:
            y=u@y
        y=apply_pointer(y,leak)
        # Equal *pulse time* halves, not equal per-gate durations.
        paths=channel_branch([y],pulse_ns/2,sc,R)
        for idx in range(len(paths)):
            for u in second:
                paths[idx]=u@paths[idx]
            paths[idx]=apply_pointer(paths[idx],leak)
        paths=channel_branch(paths,pulse_ns/2+
                             sc.readout_ns+sc.reset_ns,sc,R)
        row=_kraus_logical_ledger(
            paths, code,sc.readout_wrong_label_probability)
        row["fault"]=fault_name
        row["probability"]=float(weight)
        fault_rows.append(row)

    fidelity=sum(r["probability"]*r[
        "F_e_before_active_particle_erasure"] for r in fault_rows)
    ptr3=sum(r["probability"]*r[
        "fourth_level_pointer_probability"] for r in fault_rows)
    # Four independently lost active particles; any event an ORTHOGONAL
    # flagged erasure. No microscopic loss jump is implemented in 14D.
    survival=float(np.exp(-4*duration/sc.active_particle_T1_ns))
    return {
        "group":label,"number_of_physical_42D_unitaries":len(pulses),
        "pulse_time_ns":float(pulse_ns),
        "readout_and_reset_time_ns":float(sc.readout_ns+sc.reset_ns),
        "total_elapsed_ns":float(duration),
        "pointer_hilbert_levels":4,"pointer_leakage_level":3,
        "Gauss_sector_dimension":14,
        "channel":"exact Kraus-summed two-window pointer loss and phase diffusion; 14D electric reflection dephasing",
        "four_active_particle_no_erasure_probability":survival,
        "active_loss_is_flagged_erasure_NOT_microscopic_amplitude_damping":True,
        "unconditional_logical_F_e_before_active_erasure":float(fidelity),
        "unconditional_logical_F_e_including_flagged_erasure":float(fidelity*survival),
        "mean_pointer_fourth_level_population":float(ptr3),
        "fault_details":fault_rows,
        "no_postselection":True,
        "no_device_specific_CP22_realization":True,
        "no_ideal_code_projection_applied_as_a_gate":True,
        "noise_placement":"lumped at mid-sequence and after sequence",
        "scenario":sc.__dict__,
    }


def compare_policies(label, angles, scenario=NoiseScenario()):
    """Noise-corrected pulse schedule vs same-time passive vs zero-time.

    The passive channel experiences the SAME total elapsed time, but
    does not pay an invented readout/reset gate or pointer-label mistake.
    This is a finite-S3 physical-noise *proxy*, not a transmon QEC result.
    """
    sc=scenario.validate()
    compiler=SpectralCompiler(label,len(angles)//13)
    if len(angles)!=len(compiler.order):
        raise ValueError("angles length must be a multiple of 13")
    pulses=compiler.matrices(angles)
    cost=native_operation_expansion({
        "pulse_names":[compiler.named[i][0] for i in compiler.order],
        "angles_rad":list(map(float,angles)), "pulse_count":len(angles),
    },readout_ns=0.)
    pulse_ns=cost["hypothetical_serial_duration_ns_including_readout"]
    candidate=run_instrument(label,pulses,pulse_ns,sc)
    total=candidate["total_elapsed_ns"]
    passive_sc=replace(sc,readout_ns=0.,reset_ns=0.,
                       readout_wrong_label_probability=0.,
                       pointer_leak_angle_per_window=0.)
    matched=run_instrument(label,[],total,passive_sc)
    immediate=run_instrument(
        label,[],0.,replace(passive_sc,active_particle_T1_ns=
                           sc.active_particle_T1_ns))
    return {
        "group":label,"control_angles":len(angles),
        "candidate":candidate,
        "same_elapsed_time_passive":matched,
        "zero_elapsed_time_no_gate_reference":immediate,
        "candidate_minus_same_time_passive_including_erasure":float(
            candidate["unconditional_logical_F_e_including_flagged_erasure"]
            - matched["unconditional_logical_F_e_including_flagged_erasure"]),
        "candidate_beats_matched_passive":bool(
            candidate["unconditional_logical_F_e_including_flagged_erasure"]
            > matched["unconditional_logical_F_e_including_flagged_erasure"]
            + 1e-10),
        "native_expansion_scenario":{
            k:v for k,v in cost.items() if k!="expanded_primitive_schedule"
        },
        "limitations":[
            "Finite-S3 Gauss 14D is NOT five ququarts or continuous SU3.",
            "Pointer four-level pure loss, diffusion and coherent leakage are modeled; mediator and system 4-level microscopic jump operators are NOT.",
            "Only a collective gauge-invariant electric phase noise acts on 14D system; active-fermion loss is independent flagged erasure.",
            "Exact Kraus CPTP noise is lumped into two windows, not during individual gates.",
            "Pulse durations, pointer hardware and readout/reset rates are illustrative.",
            "No physical CP22 gate-to-finite-S3 charged bus correspondence is established.",
            "The logical code overlap is a fidelity diagnostic; no executable decoder projection is credited.",
            "Analog refocusing comparison requires a distinct gauge/device-compatible Hamiltonian mapping.",
        ],
    }


def main():
    # Use actual previously implemented, independently certified PR70
    # pulse products rather than substitute an ideal recovery U42.
    search=run_continuation(depths=(1,2,4,8),labels=("AS","BS"),
                            capture_iter=35,coherent_iter=50)
    rows=[]
    for row in search["results"]:
        if row["cycles"] != 8:
            continue
        for scenario_name,noise in (
            ("illustrative", NoiseScenario()),
            ("no_pointer_leakage", NoiseScenario(
                pointer_leak_angle_per_window=0.)),
        ):
            result=compare_policies(row["group"],row["angles_rad"],noise)
            result["scenario_label"]=scenario_name
            result["source_continuation_depth"]=row["cycles"]
            rows.append(result)
            print("FINITE_S3_NOISE",json.dumps({
                "group":row["group"],"scenario":scenario_name,
                "candidate_F_e":result["candidate"][
                    "unconditional_logical_F_e_including_flagged_erasure"],
                "matched_passive_F_e":result["same_elapsed_time_passive"][
                    "unconditional_logical_F_e_including_flagged_erasure"],
                "candidate_beats_matched":result[
                    "candidate_beats_matched_passive"],
                "pointer_leak":result["candidate"][
                    "mean_pointer_fourth_level_population"],
            },sort_keys=True),flush=True)
    out={
        "research_target":"finite-S3 104-control physical-noise proxy with matched passive",
        "reproduction":"PR70 13->26->52->104 branch-aware continuation",
        "scenarios":rows,
        "hardware_mapping":"not device calibrated; no native CP22 mapping",
        "analog_comparator":"not evaluated without identical physical simulator",
    }
    Path("colour_s3_open_system_results.json").write_text(
        json.dumps(out,indent=2,sort_keys=True)+"\n")
    return out


if __name__=="__main__":
    main()
