"""Continuous gauge-covariant three-register dressed swap for finite S3 code.

This is an explicit exact Hamiltonian construction:
    H_S=(I-S_dressed)/2, acting on matter A0, group link, matter B0.
    H_A=(I-SWAP_A0,A1)/2, acting on both matter qutrits at A.
For every time t, exp(-i t H_S) and exp(-i t H_A) commute with the
two vertex gauge groups. At t=pi each equals its corresponding swap.

The PR48 syndrome pointer needs conditional pulses:
|active><active|_pointer tensor H_S (FOUR physical registers),
and the corresponding conditional onsite H_A (THREE registers).
It has two such edge-control and two onsite-control pulse windows.
This is NOT a circuit made only of gauge-preserving three-register
operations including the pointer. No device realizing these multi-body
Hamiltonians is claimed.

Unlike the PR49 16 arbitrary two-register gates, this continuous pulse
path preserves Gauss law at *all* sampled times and analytically at all
times. Compare the two implementations with:
1. exact noiseless syndrome and recovery, and
2. equal TOTAL notional time and equal independent per-register Gaussian
   number-phase diffusion rate (CPTP after ensemble averaging).
The random diagonal link/matter dephasing itself need NOT be gauge
preserving. This is a sensitivity comparison, not a physical cost win.
"""
from __future__ import annotations

from pathlib import Path
import json
import numpy as np

from colour_gauge_link_s3 import (
    DIM, G, MATTER, physical_isometry, dressed_exchange,
    transform, projector,
)
from colour_s3_edge_logical_syndrome import (
    onsite_swap, pointer_circuit, code_isometry,
)
from colour_s3_two_register_syndrome import (
    prepared_pointer, phase_readout, pair_gate, measured_syndrome_factors,
    correction_factors, circuit_gauge_diagnostics,
)

# Register order: pointer (3), link (6), matter A0,A1,B0,B1 (3 each).
PHYSICAL_DIMS=(6,3,3,3,3)
POINTER_DIMS=(3,)+PHYSICAL_DIMS
LOGICAL=4
PULSE_TIME=16.0
FEEDBACK_TIME=3.0


def involution_evolution(y,op,theta):
    """e^{-i theta (I-op)/2}y with op^2=I. True continuous unitary."""
    if not np.isfinite(theta):raise ValueError("finite pulse angle required")
    phase=np.exp(-1j*theta)
    return (1+phase)/2*y+(1-phase)/2*op(y)


def dressed_hamiltonian(y):
    """H_S y=(I-S)y/2 has support on A0-link-B0 (not A1/B1)."""
    return (y-dressed_exchange(y))/2


def onsite_hamiltonian(y):
    """H_A y=(I-A)y/2 supported on vertex A0-A1."""
    return (y-onsite_swap(y))/2


def physical_pulse(y,kind,angle=np.pi):
    if kind=="edge":
        return involution_evolution(y,dressed_exchange,angle)
    if kind=="onsite":
        return involution_evolution(y,onsite_swap,angle)
    raise ValueError(kind)


def controlled_pulse(z,kind,active,angle=np.pi):
    """Gauge-neutral pointer controls physical H, making 4-/3-register H.

    Non-active pointer branches are identity. No global code projectors
    are invoked, and full physical 486D amplitudes are propagated.
    """
    x=np.asarray(z,complex)
    assert x.ndim==3 and x.shape[0]==3 and x.shape[1]==DIM
    active=set(active)
    if not active.issubset({0,1,2}):raise ValueError(active)
    y=x.copy()
    for label in active:
        y[label]=physical_pulse(x[label],kind,angle)
    return y


def control_schedule():
    """Exactly the PR48 pointer branches (S then A, repeated on label 2)."""
    return (
        ("edge",(1,2)),
        ("onsite",(1,2)),
        ("edge",(2,)),
        ("onsite",(2,)),
    )


def gauge_acceptance(z):
    """Gauss acceptance for each vertex including the unmeasured pointer."""
    x=np.asarray(z,complex)
    assert x.ndim==3 and x.shape[0]==3 and x.shape[1]==DIM
    norm=np.vdot(x,x).real
    assert norm>1e-12
    return {v:float(sum(np.vdot(projector(x[p],v),projector(x[p],v)).real
                         for p in range(3))/norm) for v in ("A","B")}


def continuous_gauge_diagnostics(K,substeps=(0.,.25,.5,.75,1.)):
    """Sample every Hamiltonian pulse, including intermediate pulse times."""
    z=prepared_pointer(K)
    rows=[{"step":"pointer_preparation","fraction":0.,
           "Gauss_accept":gauge_acceptance(z)}]
    maximum_cov=0.
    for step,(kind,active) in enumerate(control_schedule(),1):
        for f in substeps:
            z_mid=controlled_pulse(z,kind,active,np.pi*f)
            rows.append({"step":step,"pulse_type":kind,
                         "register_support":(4 if kind=="edge" else 3),
                         "fraction":float(f),"Gauss_accept":gauge_acceptance(z_mid)})
            for v in ("A","B"):
                for h in range(G):
                    # On physical input, invariance throughout the pulse.
                    rotated=np.stack([transform(a,h,v) for a in z],axis=0)
                    lhs=controlled_pulse(rotated,kind,active,np.pi*f)
                    rhs=np.stack([transform(a,h,v) for a in z_mid],axis=0)
                    maximum_cov=max(maximum_cov,float(np.linalg.norm(lhs-rhs)))
        z=controlled_pulse(z,kind,active)
    result=phase_readout(z)
    assert np.linalg.norm(result-pointer_circuit(K))<2e-10
    minima={v:float(min(row["Gauss_accept"][v] for row in rows))
            for v in ("A","B")}
    assert min(minima.values())>1-1e-10
    assert maximum_cov<1e-10
    return {"trace":rows,"minimum_Gauss_acceptance":minima,
            "maximum_Gauss_covariance_residual":maximum_cov,
            "noiseless_reference_residual_fro":
                float(np.linalg.norm(result-pointer_circuit(K)))}


def gaussian_phase_noise(y,dt,rate,rng):
    """Physical-register independent Gaussian number-phase diffusion.

    Applied to ALL pointer, link and matter local registers. Exact random
    unitary => CPTP average. Each register contributes a phase with
    variance 2*dt*rate, so off-diagonal n,m decay as
    exp[-dt*rate*(n-m)^2]. Link basis-number phase is a synthetic,
    generally GAUSS-BREAKING environment, not calibrated physics.
    """
    if dt<0 or rate<0:raise ValueError("nonnegative time/rate")
    if dt==0 or rate==0:return y
    data=np.asarray(y,complex)
    dims=POINTER_DIMS if data.ndim==3 else PHYSICAL_DIMS
    assert data.shape[:-1]==((3,DIM) if data.ndim==3 else (DIM,))
    reshaped=data.reshape(dims+(data.shape[-1],)).copy()
    for axis,d in enumerate(dims):
        angle=rng.normal(0.,np.sqrt(2*dt*rate))
        diag=np.exp(1j*np.arange(d)*angle)
        shape=[1]*(len(dims)+1)
        shape[axis]=d
        reshaped*=diag.reshape(shape)
    return reshaped.reshape(data.shape)


def analog_gadget(y,rate,rng):
    """Four continuous conditional pulses, matching 16 time units total."""
    z=prepared_pointer(y)
    for kind,active in control_schedule():
        z=controlled_pulse(z,kind,active)
        z=gaussian_phase_noise(z,PULSE_TIME/4,rate,rng)
    return phase_readout(z)


def pair_gadget(y,rate,rng):
    """Exact PR49 pair synthesis, same total notional pulse time."""
    z=prepared_pointer(y)
    for _,axes,gate in measured_syndrome_factors():
        z=pair_gate(z,axes,gate)
        z=gaussian_phase_noise(z,PULSE_TIME/16,rate,rng)
    return phase_readout(z)


def analog_feedback(y,rate,rng,correction):
    """One 3-body edge pulse on errored branch; idle otherwise."""
    z=physical_pulse(y,"edge") if correction else y.copy()
    return gaussian_phase_noise(z,FEEDBACK_TIME,rate,rng)


def pair_feedback(y,rate,rng,correction):
    """Three two-register factor gates on correction; equal idle otherwise."""
    z=y.copy()
    for _,axes,gate in correction_factors():
        if correction:z=pair_gate(z,axes,gate)
        z=gaussian_phase_noise(z,FEEDBACK_TIME/3,rate,rng)
    return z


def circuit_fidelity(K,architecture,rate,shots=32,seed=20261019,p=.37):
    """Ensemble of CPTP phase-noise trajectories on encoded maximally mixed state.

    Each trajectory samples independent local random phases after each
    gate, including feedback/idling. Both architectures spend 19
    notional time units. Error channel = {(1-p)I,p S}, measured and
    corrected via actual pointer outcomes with no postselection.
    """
    if architecture not in ("covariant_multibody","pair_decomposed"):
        raise ValueError(architecture)
    if not (0<=p<=1 and rate>=0 and shots>=1):raise ValueError("scenario")
    if architecture=="covariant_multibody":
        gadget,feedback=analog_gadget,analog_feedback
    else:
        gadget,feedback=pair_gadget,pair_feedback
    F=dressed_exchange(K)
    samples=[]
    rng=np.random.default_rng(seed)
    for _ in range(shots):
        fe=0.
        for src in (np.sqrt(1-p)*K,np.sqrt(p)*F):
            outcomes=gadget(src,rate,rng)
            for label,branch in enumerate(outcomes):
                final=feedback(branch,rate,rng,label==2)
                fe+=abs(np.trace(K.conj().T@final))**2/LOGICAL**2
        assert -1e-8<=fe<=1+1e-8,fe
        samples.append(float(fe))
    samples=np.array(samples)
    return {"architecture":architecture,"noise_rate_per_register_per_unit":rate,
            "source_known_dressed_swap_probability":p,"shots":shots,"seed":seed,
            "notional_syndrome_time":PULSE_TIME,
            "notional_feedback_or_idle_time":FEEDBACK_TIME,
            "mean_entanglement_fidelity":float(samples.mean()),
            "monte_carlo_standard_error":
                float(samples.std(ddof=1)/np.sqrt(shots)) if shots>1 else None,
            "individual_trajectory_fidelities":samples.tolist()}


def study(shots=32):
    W=physical_isometry()
    K,F,_=code_isometry(W)
    random_input=np.random.default_rng(4).normal(size=(DIM,2))+1j*np.random.default_rng(5).normal(size=(DIM,2))
    # Full arbitrary 486D coherent input, not just logical basis.
    z=prepared_pointer(random_input)
    for kind,active in control_schedule():
        z=controlled_pulse(z,kind,active)
    arbitrary_error=float(np.linalg.norm(phase_readout(z)-pointer_circuit(random_input)))
    assert arbitrary_error<1e-10
    for kind,op in (("edge",dressed_exchange),("onsite",onsite_swap)):
        x=np.random.default_rng(72).normal(size=(DIM,2))+1j*np.random.default_rng(73).normal(size=(DIM,2))
        assert np.linalg.norm(physical_pulse(x,kind,np.pi)-op(x))<1e-10
        h=(x-op(x))/2
        assert np.linalg.norm(
            (physical_pulse(x,kind,.00001)-x)/.00001+1j*h)<2e-4
        assert np.linalg.norm(physical_pulse(physical_pulse(x,kind,.2),kind,.3)
                  -physical_pulse(x,kind,.5))<1e-10
        # Covariance on arbitrary (not just Gauss-invariant) input.
        for v in ("A","B"):
            for g in range(G):
                assert np.linalg.norm(transform(physical_pulse(x,kind,.37),g,v)
                    -physical_pulse(transform(x,g,v),kind,.37))<1e-10
    gauge3=continuous_gauge_diagnostics(K)
    pair=circuit_gauge_diagnostics(K)
    # Two types of gates; do not equate hardware cost to factor count.
    costs={"covariant_multibody":{
        "conditional_edge_windows":2,"conditional_onsite_windows":2,
        "conditional_edge_support_registers":4,
        "conditional_onsite_support_registers":3,
        "unconditional_edge_feedback_support_registers":3,
        "total_high_level_syndrome_windows":4,
        "physical_multibody_hamiltonian_synthesized":False},
        "pair_decomposed":{
        "general_U9_or_U18_two_register_factors":16,
        "additional_conditional_feedback_pair_factors":3,
        "native_1_or_2_qudit_decomposition":False}}
    records=[]
    for i,rate in enumerate((0.,.002,.01)):
        for arch in ("covariant_multibody","pair_decomposed"):
            records.append(circuit_fidelity(K,arch,rate,shots=shots,seed=100+i*17))
    for x in records[:2]:
        assert abs(x["mean_entanglement_fidelity"]-1)<1e-10
    return {"model":"finite S3 two-vertex four-qutrit link gauge simulator, new 4D phase code",
            "physical_dimension":DIM,"physical_Gauss_rank":14,
            "logical_dimension":4,"source_error":"one known gauge-invariant dressed-swap pair",
            "continuous_gauge_covariant_path":gauge3,
            "prior_pair_path_minimum_Gauss_accept":pair["minimum_intermediate_Gauss_accept"],
            "exact_full_system_reference_residual":arbitrary_error,
            "resource_comparison":costs,
            "matched_noise_comparison":records,
            "assumed_notional_syndrome_time":PULSE_TIME,
            "assumed_notional_feedback_time":FEEDBACK_TIME,
            "main_scientific_result":"H=(I-S_dressed)/2 is a gauge-covariant three-register Hamiltonian throughout; pointer control is FOUR-register and not yet synthesized.",
            "limitations":[
                "No actual four-register neutral-pointer control construction, physical coupling graph, pulse generator, or calibrated native gate times",
                "Pair factor and multibody gates are assigned equal total notional time by assumption, not device measurement",
                "Only Gaussian computational-number phase diffusion, not T1, non-Markovian gauge bath, coherent overrotations, or link leakage",
                "Local link and matter dephasing faults generally violate gauge constraints even though the ideal pulse does not",
                "Four-dimensional finite S3 phase code, NOT nine-quark 21-dimensional code or 65-syndrome instrument",
                "No architecture win inferred from a toy noise comparison"]}


def main():
    out=study()
    Path("colour_s3_covariant_multibody_control_results.json").write_text(
        json.dumps(out,indent=2,sort_keys=True)+"\n")
    print("S3_COVARIANT_MULTIBODY",json.dumps({
        "exact_reference_residual":out["exact_full_system_reference_residual"],
        "new_path_minimum_Gauss_accept":
            out["continuous_gauge_covariant_path"]["minimum_Gauss_acceptance"],
        "old_pair_path_minimum_Gauss_accept":
            out["prior_pair_path_minimum_Gauss_accept"],
        "resources":out["resource_comparison"]},sort_keys=True),flush=True)
    for x in out["matched_noise_comparison"]:
        print("S3_MATCHED_PHASE_NOISE",json.dumps({
            k:v for k,v in x.items() if k!="individual_trajectory_fidelities"
        },sort_keys=True),flush=True)

if __name__=="__main__":
    main()
