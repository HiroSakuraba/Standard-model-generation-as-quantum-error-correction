"""Duration-matched four-logical-state DIRECT storage vs finite-S3 QEC.

This is a CONTROLLED COMPARISON of clearly delimited logical/noise
channels, not a claim of hardware or microscopic fault equivalence.

Encoded memory: 4D finite-S3 gauge phase code with 4 active fermions;
15 finite-U 495D Hubbard pulses, 6 abstract gauge-link transports,
300ns ideal readout, and known dressed swap probability p=.37.
The encoded no-particle-loss fidelity is calculated by the merged
PR57 microscopic-feedback simulator, not estimated from isolated gate
fidelity.

Unencoded memory: SAME 4D LOGICAL Hilbert space stored in TWO direct
colour-qubit orbitals (two occupied fundamental qutrits, colours 0/1).
No encoding, syndrome, correction, or control gates. Idle for precisely
the encoded TOTAL elapsed time. Identical per-active-fermion T1, same
per-colour-mode gamma if relevant.

Two different fair-to-label logical challenges:
  clean bare: NO logical error, best direct-storage noise reference.
  challenged bare: one specified logical X on first colour-qubit with
    probability p=.37, UNCORRECTED. This matches the encoded model's
    ERROR PROBABILITY AND logical dimension, NOT the microscopic
    dressed-exchange quantum channel (which has no bare analogue).
Neither logical comparison is a hardware-equivalent error benchmark.

Exact analytical bare colour-resolved dephasing from L_n=sqrt(2gamma)n:
  c=exp(-2gamma*time_dimensionless) per colour-qubit coherence,
  F_e(no_loss)=((1+c)/2)^2 with no logical X error,
  F_e(no_loss, X-probability=p)=(1-p)*((1+c)/2)^2.
Colour-blind orbital charge has n_total=1 identically on a direct
single-occupation qubit, therefore F_e(no_loss)=1 without logical X.

Uniform fermion erasure:
  eta_bare=exp(-2*time/T1); eta_code=exp(-4*time/T1).
Durations use the same hypothetical t/h, link transport and readout
assumptions from PR57; code/noise parameters are sensitivity scenarios.

Crossover eta_code*F_code > eta_bare*(1-p)*F_bare holds iff
  F_code > (1-p)*F_bare and
  T1 > [2*time]/log(F_code/[(1-p)*F_bare]).
With ZERO bare logical faults and ZERO dephasing, F_bare=1 and
F_code<=1, so the code cannot beat perfect direct idle storage.
"""
from __future__ import annotations

from pathlib import Path
import json
import numpy as np

from colour_gauge_link_s3 import physical_isometry
from colour_s3_edge_logical_syndrome import code_isometry
from colour_su3_hubbard_virtual_exchange import optimize_pulse
from colour_su3_hubbard_loss_dephase_feedback import (
    process_sample, time_unit_ns,
    DEFAULT_T1_US, HOPPING_OVER_H_MHZ,
    LINK_TRANSPORT_NS, POINTER_READOUT_NS,
)

DLOG=4
N_CODE_FERMIONS=4
N_BARE_FERMIONS=2
SOURCE_ERROR_PROBABILITY=.37


def physical_duration_ns(pulse_duration_inverse_t,
                         hopping_mhz=HOPPING_OVER_H_MHZ,
                         transport_ns=LINK_TRANSPORT_NS,
                         readout_ns=POINTER_READOUT_NS):
    """15 Hubbard pulses and six ideal link transports + readout/idle.

    Relative U/t, optimized dimensionless pulse and error/no-loss
    reference remain fixed when t/h is varied. This does NOT model
    a corresponding change of U, bandwidth, gate/control noise, or T1.
    """
    args=(pulse_duration_inverse_t,hopping_mhz,transport_ns,readout_ns)
    if not all(np.isfinite(v) for v in args):
        raise ValueError("All timing parameters must be finite")
    if pulse_duration_inverse_t<=0 or hopping_mhz<=0 or transport_ns<0 or readout_ns<0:
        raise ValueError("Positive pulse and hopping; nonnegative other durations")
    return 15*pulse_duration_inverse_t*time_unit_ns(hopping_mhz) + (
        6*transport_ns + readout_ns)


def bare_colour_dephasing_fidelity(elapsed_inverse_t,gamma=0.,
                                   model="orbital_charge"):
    """Analytic physical two-qutrit direct memory, 4-dimensional subspace.

    Two qutrits carry a logical qubit in colour |0>,|1> each. One
    occupied fundamental per orbital, two independent occupied
    fermions. Under colour-resolved n_{orb,colour} noise each logical
    offdiagonal accumulates exp(-2 gamma t), exactly. Four independent
    colour-mode operators in two colour qutrits; unused third colour
    mode is zero and has no action. Under colour-blind N_orbital the
    occupation is always one, so the channel is identity.
    """
    if model not in ("orbital_charge","colour_resolved"):
        raise ValueError(model)
    if not np.isfinite(gamma) or gamma<0 or not np.isfinite(elapsed_inverse_t) or elapsed_inverse_t<0:
        raise ValueError("Nonnegative finite gamma and elapsed time")
    if model=="orbital_charge":
        return 1.
    coherence=np.exp(-2*gamma*elapsed_inverse_t)
    return float(((1+coherence)/2)**2)


def exactly_flagged_fermion_survival(n_fermions,elapsed_ns,T1_us):
    """Uniform independent c_mode jumps -> exact no-jump probability."""
    if not isinstance(n_fermions,int) or n_fermions<0 or not np.isfinite(elapsed_ns) or elapsed_ns<0 or not np.isfinite(T1_us) or T1_us<=0:
        raise ValueError("Valid fermion count, elapsed time and lifetime needed")
    return float(np.exp(-n_fermions*elapsed_ns/(1000*T1_us)))


def bare_reference(duration_ns, hopping_mhz=HOPPING_OVER_H_MHZ,
                   T1_us=DEFAULT_T1_US, gamma=0.,
                   model="orbital_charge",source_p=SOURCE_ERROR_PROBABILITY):
    """Both an unchallenged direct memory and matched logical-fault analogue.

    The 'challenged' reference uses X_on_first_colour_qubit with p,
    chosen because its trace is zero and therefore the entanglement
    fidelity of its uncorrected branch vanishes under diagonal number
    dephasing. NOT a microscopically identical dressed-SWAP process.
    """
    if not np.isfinite(source_p) or not 0<=source_p<=1:
        raise ValueError("Fault probability must be between 0 and 1")
    elapsed=duration_ns/time_unit_ns(hopping_mhz)
    fe=bare_colour_dephasing_fidelity(elapsed,gamma,model)
    survival=exactly_flagged_fermion_survival(2,duration_ns,T1_us)
    return {
        "equal_logical_dimension":DLOG,
        "physical_colour_qutrits":2,
        "active_fermions":2,
        "storage_duration_ns":float(duration_ns),
        "T1_us":float(T1_us),
        "no_loss_probability":survival,
        "no_fault_before_loss_F_e":fe,
        "no_fault_unconditional_F_e":float(fe*survival),
        "uncorrected_logical_X_probability":float(source_p),
        "matched_probability_X_fault_before_loss_F_e":float((1-source_p)*fe),
        "matched_probability_X_fault_unconditional_F_e":float((1-source_p)*fe*survival),
        "logical_channel_equivalence":"Probability and dimension only, NOT same microscopic gauge-dressed exchange",
        "comparison_assumptions":"Bare two-colour-qubit memory uses same time/T1 and mode-number dephasing rate, no control or readout errors"}


def required_T1_crossover_us(code_no_loss_fidelity,bare_no_loss_fidelity,
                              duration_ns,bare_fault_probability):
    """Exact analytical threshold for 4-active code vs 2-active bare.

    All other parameters and deterministic coherent code errors held
    constant as lifetime changes. >= can tie only asymptotically in
    the infinite-lifetime case; return None if impossible at finite T1.
    """
    a=float(code_no_loss_fidelity)
    b=float(bare_no_loss_fidelity)
    p=float(bare_fault_probability)
    if not all(np.isfinite(x) for x in (a,b,p,duration_ns)) or not (0<=a<=1+1e-9 and 0<=b<=1+1e-9 and 0<=p<=1 and duration_ns>=0):
        raise ValueError("Invalid probabilities or duration")
    target=(1-p)*b
    if target==0:
        return 0. if a>0 else None
    if a<=target or duration_ns==0:
        return None
    return float(2*duration_ns/(1000*np.log(a/target)))


def analyze_anchor(U,K,p_error=SOURCE_ERROR_PROBABILITY):
    """Recompute merged PR57 495D full-circuit deterministic benchmark."""
    opt=optimize_pulse(U)
    pulse=opt["time"]; d=opt["active_detuning"]
    r=process_sample(K,U,d,pulse,shots=1,seed=int(U),gamma=0.,
                     substeps=1,T1_us=DEFAULT_T1_US)
    physical=physical_duration_ns(pulse)
    assert abs(physical-r["nominal_total_duration_ns"])<2e-8
    code_no_loss=r["mean_F_e_conditional_on_no_loss"]
    bare=bare_reference(physical,source_p=p_error)
    code_conditional=code_no_loss*exactly_flagged_fermion_survival(
        N_CODE_FERMIONS,physical,DEFAULT_T1_US)
    assert abs(code_conditional-r["mean_F_e_unconditional_including_flagged_loss"])<1e-10
    crossover=required_T1_crossover_us(
        code_no_loss,1.,physical,p_error)
    return {"U_over_t":float(U),
            "optimized_15state_pulse_inverse_t":pulse,
            "optimized_active_pointer_detuning":d,
            "total_duration_us_at_hypothetical_hopping_10MHz":physical/1000,
            "encoded_microscopic_feedback_no_loss_F_e":code_no_loss,
            "encoded_unconditional_F_e_at_T1_30us":code_conditional,
            "bare_no_fault_unconditional_F_e_at_T1_30us":
                bare["no_fault_unconditional_F_e"],
            "bare_37pct_X_fault_unconditional_F_e_at_T1_30us":
                bare["matched_probability_X_fault_unconditional_F_e"],
            "encoded_vs_bare_no_fault_difference_at_T1_30us":
                code_conditional-bare["no_fault_unconditional_F_e"],
            "encoded_vs_challenged_bare_difference_at_T1_30us":
                code_conditional-bare["matched_probability_X_fault_unconditional_F_e"],
            "break_even_T1_us_vs_37pct_bare_X_fault_colourblind":
                crossover,
            "break_even_T1_us_vs_perfect_bare_idle":None}


def duration_matched_grid(anchors,source_p=SOURCE_ERROR_PROBABILITY):
    """Closed-form sensitivity without rerunning the 495D Fock dynamics."""
    rows=[]
    for a in anchors:
        U=a["U_over_t"];F=a["encoded_microscopic_feedback_no_loss_F_e"]
        pulse=a["optimized_15state_pulse_inverse_t"]
        for hopping in (5.,10.,20.,40.):
            for readout in (100.,300.,1000.):
                time_ns=physical_duration_ns(
                    pulse,hopping_mhz=hopping,readout_ns=readout)
                threshold=required_T1_crossover_us(F,1.,time_ns,source_p)
                for T1 in (30.,50.,75.,100.,300.):
                    bare=bare_reference(time_ns,hopping_mhz=hopping,T1_us=T1,
                                        source_p=source_p)
                    code=F*exactly_flagged_fermion_survival(4,time_ns,T1)
                    rows.append({
                        "U_over_t":U,"hopping_t_over_h_MHz":hopping,
                        "readout_delay_ns":readout,"T1_us":T1,
                        "equal_elapsed_time_us":time_ns/1000,
                        "encoded_unconditional_F_e":float(code),
                        "bare_no_fault_unconditional_F_e":
                            bare["no_fault_unconditional_F_e"],
                        "bare_X_challenged_unconditional_F_e":
                            bare["matched_probability_X_fault_unconditional_F_e"],
                        "encoded_beats_clean_bare":bool(code>bare["no_fault_unconditional_F_e"]),
                        "encoded_beats_challenged_bare":bool(code>bare["matched_probability_X_fault_unconditional_F_e"]),
                        "analytic_break_even_T1_us_vs_challenged_bare":threshold,
                    })
    return rows


def noisy_colour_bare_comparison(K,U=30.,gamma=.0002,shots=6):
    """New fair-in-time bare benchmark for both noisy PR57 environments.

    The PR57 noise ensembles use independent physical Lindblad
    operators; same per-operator gamma is not equal total bath power.
    The bare memory has fewer occupied colour modes. Keep the full
    expected no-loss channel analytic (no stochastic baseline error)
    and encoded estimates explicit Monte Carlo.
    """
    opt=optimize_pulse(U)
    time_ns=physical_duration_ns(opt["time"])
    reports=[]
    for i,model in enumerate(("orbital_charge","colour_resolved")):
        r=process_sample(K,U,opt["active_detuning"],opt["time"],
            model=model,gamma=gamma,shots=shots,seed=126+i,
            T1_us=DEFAULT_T1_US,substeps=2)
        bare=bare_reference(time_ns,model=model,gamma=gamma,
                            T1_us=DEFAULT_T1_US)
        threshold=required_T1_crossover_us(
            r["mean_F_e_conditional_on_no_loss"],
            bare["no_fault_before_loss_F_e"],time_ns,SOURCE_ERROR_PROBABILITY)
        reports.append({
            "model":model,"rate_per_Lindblad_occupation_operator":gamma,
            "encoded_trajectory_count":shots,
            "encoded_mean_conditional_on_no_loss_F_e":
                r["mean_F_e_conditional_on_no_loss"],
            "encoded_conditional_MC_standard_error":
                r["Monte_Carlo_standard_error_conditional"],
            "encoded_unconditional_loss_dephasing_F_e":
                r["mean_F_e_unconditional_including_flagged_loss"],
            "bare_no_fault_unconditional_F_e":
                bare["no_fault_unconditional_F_e"],
            "bare_37pct_X_fault_unconditional_F_e":
                bare["matched_probability_X_fault_unconditional_F_e"],
            "threshold_vs_bare_X_us_from_stochastic_code_estimate":threshold,
            "same_gamma_is_not_equal_total_noise_power":True
        })
    return reports


def run():
    W=physical_isometry()
    K,_,_=code_isometry(W)
    anchors=[analyze_anchor(U,K) for U in (20.,30.,40.)]
    grid=duration_matched_grid(anchors)
    assert len(grid)==3*4*3*5
    assert all(not r["encoded_beats_clean_bare"] for r in grid)
    no_physical_winner=[
        a["U_over_t"] for a in anchors
        if a["encoded_unconditional_F_e_at_T1_30us"]>
            a["bare_37pct_X_fault_unconditional_F_e_at_T1_30us"]]
    assert not no_physical_winner
    noisy=noisy_colour_bare_comparison(K)
    return {
        "title":"Duration-matched two-fermion bare memory vs four-fermion finite-S3 4D logical code",
        "logical_dimension_per_architecture":4,
        "source_known_error_probability":SOURCE_ERROR_PROBABILITY,
        "bare_challenge":"uncorrected X on one of two colour qubits with 37% probability (phenomenological analogue, not a microscopic dressed exchange)",
        "clean_bare_baseline":"No logical fault, no gate/readout, same elapsed total time, T1 and per-mode dephasing",
        "encoded_reference":"PR57 microscopic 495D Hubbard syndrome and correction, including coherent infidelities and occupation leakage",
        "anchors_no_dephasing":anchors,
        "duration_T1_hopping_readout_grid":grid,
        "grid_points":len(grid),
        "encoded_better_than_clean_bare_any_grid":False,
        "encoded_better_than_bare_with_challenge_at_T1_30us":no_physical_winner,
        "analytic_bare_noise_vs_encoded_stochastic_noise_at_U30":noisy,
        "limits":[
            "Direct bare ququart uses two occupied fermionic colour-qutrits with two colours per orbital; coded finite-S3 memory uses four active fermions plus ideal spectators/reference",
            "Bare 37% logical X is probability/dimension-matched but not the same physical noise process as the encoded dressed exchange",
            "Bare clean-no-logical-error memory is the ideal direct storage reference, not a challenged error-correction experiment",
            "Both memories use exactly the same encoded total elapsed duration per grid cell, not the same gate count",
            "Encoded microscopic gate dynamics recalculated only at original 10MHz anchor; the hopping sensitivity holds the dimensionless pulse waveform, coherent F_e, and T1 fixed",
            "Readout-delay sensitivity assumes no extra measurement assignment error, and ignores ancilla preparation/routing time",
            "For bare colour-resolved dephasing, rate gamma is per colour occupation operator; four bare occupied-colour operators vs twelve code-mode operators; equal gamma is not equal total bath power",
            "Loss is exact independently flagged fermion erasure on occupied modes but no post-jump propagation; only active code fermions lose",
            "Only finite S3 four-dimensional phase code; not full SU3 gauge-link hardware or original nine-quark 21D 65-syndrome decoder",
            "No hardware measurements or general QEC advantage/threshold follows from this synthetic comparison",
        ]}


def main():
    result=run()
    Path("colour_su3_duration_matched_bare_results.json").write_text(
        json.dumps(result,indent=2,sort_keys=True)+"\n")
    for a in result["anchors_no_dephasing"]:
        print("BARE_VS_CODE",json.dumps(a,sort_keys=True),flush=True)
    for n in result["analytic_bare_noise_vs_encoded_stochastic_noise_at_U30"]:
        print("MATCHED_COLOUR_NOISE",json.dumps(n,sort_keys=True),flush=True)
    print("COMPARED",result["grid_points"],"matched hardware-scenario grid cells",
          flush=True)


if __name__=="__main__":
    main()
