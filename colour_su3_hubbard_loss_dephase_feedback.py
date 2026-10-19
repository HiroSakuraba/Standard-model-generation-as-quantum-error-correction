"""Microscopic charged-bus feedback and CPTP dephasing/fermion-loss audit.

Builds directly on the MERGED persistent four-fermion 495D Fock model
(PR56). No reset/projection of doublons between syndrome *or feedback*
pulses. Five feedback operations: Hubbard(A0), link transport,
Hubbard(B0), link transport, Hubbard(A0), with detuning active only on
measured syndrome pointer label 2. Inactive feedback labels undergo
the actual off-detuned Hubbard pulses, NOT magically identity.

Dephasing is a random-unitary unraveling of independent, Hermitian
number-operator Lindblad channels: L_j=sqrt(2 gamma) N_j.
We sample Gaussian phases N(0,2 gamma dt) at pulse subdivisions and
transport/measurement intervals. Two distinct environments:
  "orbital_charge": N_j=sum_colour n_orbital, locally SU3-invariant.
  "colour_resolved": n_orbital,colour, generally GAUSS-BREAKING.

Uniform single-particle loss has L_c=sqrt(kappa) c_c for all 12 active
modes. On the fixed N=4 fermion sector, sum L†L=4*kappa*I.
The exact NO-JUMP propagator is scalar exp(-2*kappa*T), and the
probability of at least one jump is 1-exp(-4*kappa*T).
All jumped trajectories occupy fewer than 4 particles, so are
represented by an ORTHOGONAL FLAGGED ERASURE SECTOR (zero logical
fidelity). This is an exact CPTP loss-vs-no-loss instrument under
uniform, independent loss; it does NOT evolve the 3-,2-,1-fermion
post-jump dynamics or any particle repumping.

Ideal readout/projective pointer and exact assumed finite S3 link
transport remain. No calibrated device, no continuous SU3 gauge field
or original 21D nine-quark 65-syndrome decoder.
"""
from __future__ import annotations

from pathlib import Path
import json
import numpy as np

from colour_gauge_link_s3 import physical_isometry, dressed_exchange
from colour_s3_edge_logical_syndrome import code_isometry
from colour_su3_hubbard_virtual_exchange import optimize_pulse
from colour_su3_hubbard_persistent_syndrome import (
    F_STATES, N_FOCK, SITE, N_MODES, G, DIM,
    initialize, local_fock_pulse, apply_microscopic_gate,
    link_transport, microscopic_schedule, project_matter,
    occupation_leakage, gauss_acceptance, exact_full_channel
)

DLOG=4
LINK_TRANSPORT_NS=100.
POINTER_READOUT_NS=300.
HOPPING_OVER_H_MHZ=10.
DEFAULT_T1_US=30.
OCC=np.array([[(state>>k)&1 for k in range(N_MODES)]
              for state in F_STATES],dtype=np.float64)
ORBITAL_OCC=OCC.reshape(N_FOCK,4,3).sum(axis=2)
assert all(ORBITAL_OCC.sum(axis=1)==4)
assert np.allclose(OCC.sum(axis=1),4.)


def time_unit_ns(hopping_over_h_MHz=HOPPING_OVER_H_MHZ):
    """t=1 is an angular hopping amplitude: t/h=f MHz."""
    if not np.isfinite(hopping_over_h_MHz) or hopping_over_h_MHz<=0:
        raise ValueError("positive finite hopping")
    return 1000./(2*np.pi*hopping_over_h_MHz)


def flagged_uniform_loss(total_ns,T1_us=DEFAULT_T1_US):
    """EXACT uniform independent fermion loss/erasure on the N=4 sector.

    Kraus no-jump operator has norm sqrt(eta), eta=e^(-4*T/T1).
    Remaining Kraus outcomes map into mutually orthogonal flagged
    missing-particle sectors, and have no overlap with the encoded code.
    The stored scalar eta is therefore the exact branch probability,
    not a Monte Carlo approximation or postselected memory fidelity.
    """
    if not np.isfinite(total_ns) or total_ns<0 or not np.isfinite(T1_us) or T1_us<=0:
        raise ValueError("nonnegative time and positive T1")
    eta=float(np.exp(-4*total_ns/(1000*T1_us)))
    return {"no_jump_probability":eta,
            "flagged_particle_loss_probability":float(-np.expm1(
                -4*total_ns/(1000*T1_us))),
            "lost_particle_branch_logical_fidelity":0.,
            "loss_jumps":["sqrt(1/T1) c_mode for each of the 12 active colour modes"],
            "total_active_fermions":4,
            "scope":"Exact uniform independent loss as flagged erasure; no evolution of fewer-than-four-particle jump states"}


def number_phase_kick(z,dt,gamma,rng,model):
    """One Gaussian Lindblad-number-dephasing unraveling step.

    For L_j=sqrt(2 gamma) N_j, a Gaussian random phase with variance
    2 gamma dt gives exp[-gamma dt (n_j-n'_j)^2] after averaging.
    Every realization is unitary, so averaging gives a CPTP channel.
    """
    if model not in ("orbital_charge","colour_resolved"):
        raise ValueError("unknown noise model")
    if not np.isfinite(gamma) or gamma<0 or not np.isfinite(dt) or dt<0:
        raise ValueError("nonnegative finite dt and gamma")
    if gamma==0 or dt==0:
        return z.copy()
    occupations=ORBITAL_OCC if model=="orbital_charge" else OCC
    angles=rng.normal(0.,np.sqrt(2*gamma*dt),size=occupations.shape[1])
    phasor=np.exp(-1j*occupations@angles)
    return z*phasor[None,None,None,None,:,None,None]


def apply_split_microscopic_gate(z,site,active,U,detuning,T,gamma,rng,
                                model,substeps):
    """Same 495D state throughout; phase environment during every pulse."""
    if not isinstance(substeps,int) or substeps<1:
        raise ValueError("positive integer substeps")
    out=z
    for _ in range(substeps):
        out=apply_microscopic_gate(out,site,active,U,detuning,T/substeps)
        out=number_phase_kick(out,T/substeps,gamma,rng,model)
    return out


def full_no_jump_branch(K,U,detuning,T,source_p=.37,
                        gamma=0.,model="orbital_charge",substeps=2,
                        rng=None,transport_ns=LINK_TRANSPORT_NS,
                        readout_ns=POINTER_READOUT_NS):
    """All 12 Hubbard syndrome pulses AND 3 Hubbard feedback pulses.

    Pointer branch 2 is conditionally corrected; other branches evolve
    under the real off-detuning pulse rather than a fictitious identity.
    Dephasing acts on persistent Fock occupations even between pulses.
    """
    if rng is None:rng=np.random.default_rng(20261019)
    if not 0<=source_p<=1:raise ValueError("source_p")
    unit=time_unit_ns()
    link_time=transport_ns/unit
    read_time=readout_ns/unit
    fe=0.;norm=0.;final_occup=0.;max_occup=0.
    pointer=np.zeros(3)
    final_gauss={v:0. for v in ("A","B")}
    for weight,E in (
            (1-source_p,K),(source_p,dressed_exchange(K))):
        z=initialize(np.sqrt(weight)*E)
        max_source=0.
        for name,active in microscopic_schedule():
            if name=="T":
                z=link_transport(z)
                z=number_phase_kick(z,link_time,gamma,rng,model)
            else:
                z=apply_split_microscopic_gate(
                    z,name,active,U,detuning,T,gamma,rng,model,substeps)
            max_source=max(max_source,occupation_leakage(z))
        z=np.fft.fft(z,axis=0)/np.sqrt(3)
        z=number_phase_kick(z,read_time,gamma,rng,model)
        # Actual five-step charge-bus feedback, including the real
        # off-pointer Hubbard pulses on labels 0/1.
        for name in ("A0","T","B0","T","A0"):
            if name=="T":
                z=link_transport(z)
                z=number_phase_kick(z,link_time,gamma,rng,model)
            else:
                z=apply_split_microscopic_gate(
                    z,name,(2,),U,detuning,T,gamma,rng,model,substeps)
            max_source=max(max_source,occupation_leakage(z))
        norm+=np.vdot(z,z).real/DLOG
        pointer+=np.array([np.vdot(x,x).real/DLOG for x in z])
        final_occup+=weight*occupation_leakage(z)
        max_occup+=weight*max_source
        if gamma:
            gau=gauss_acceptance(z)
            for v in ("A","B"):
                final_gauss[v]+=weight*gau[v]
        # Logical process fidelity, sum over all pointer and charge
        # environment outcomes, without normalizing escaped population.
        for p in range(3):
            for loc in (0,): # source code and charge reference at A
                for ref in range(3):
                    for med in range(3):
                        block=project_matter(z[p,:,:,loc,:,:,:],med,ref)
                        fe+=abs(np.trace(K.conj().T@block))**2/DLOG**2
    assert abs(norm-1)<2e-9,(norm,source_p)
    assert 0<=fe<=1+1e-9,fe
    assert abs(sum(pointer)-1)<2e-9
    if not gamma:final_gauss={"A":None,"B":None}
    return {"F_e_conditional_on_no_particle_loss":float(fe),
            "final_occupation_leakage_conditional":float(final_occup),
            "max_stage_occupation_leakage_conditional":float(max_occup),
            "pointer_probability_0_1_2":pointer.tolist(),
            "full_no_jump_trace":float(norm),
            "final_extended_Gauss_acceptance":final_gauss}


def process_sample(K,U,detuning,T,model="orbital_charge",gamma=0.,
                   shots=4,seed=12345,T1_us=DEFAULT_T1_US,substeps=2,
                   loss_enabled=True):
    """Average NO-LOSS stochastic trajectories and exact CPTP erasure branch."""
    if shots<1:raise ValueError("positive shots")
    if model not in ("orbital_charge","colour_resolved"):
        raise ValueError(model)
    duration_ns=(15*T*time_unit_ns()+6*LINK_TRANSPORT_NS
                 +POINTER_READOUT_NS)
    loss=flagged_uniform_loss(duration_ns,T1_us) if loss_enabled else {
        "no_jump_probability":1.,"flagged_particle_loss_probability":0.,
        "lost_particle_branch_logical_fidelity":0.,
        "scope":"loss disabled"}
    rng=np.random.default_rng(seed)
    trajectories=[]
    for _ in range(shots):
        trajectories.append(full_no_jump_branch(
            K,U,detuning,T,gamma=gamma,model=model,
            substeps=substeps,rng=rng))
    fidelity=np.array([v["F_e_conditional_on_no_particle_loss"]
                       for v in trajectories])
    loss_eta=loss["no_jump_probability"]
    gauss=[v["final_extended_Gauss_acceptance"] for v in trajectories]
    gauss_av={v:float(np.mean([x[v] for x in gauss]))
              if gamma else None for v in ("A","B")}
    return {"U_over_t":float(U),"fermionic_dephasing_model":model,
            "per_operator_Lindblad_dephasing_rate_per_inverse_t":gamma,
            "dephasing_L":"sqrt(2 gamma) times N_orbital or n_orbital_colour",
            "source_known_dressed_swap_probability":.37,
            "n_trajectories":shots,"seed":seed,"pulse_substeps":substeps,
            "total_microscopic_Hubbard_pulses":15,
            "total_ideal_link_transports":6,
            "nominal_total_duration_ns":duration_ns,
            "measurement_delay_ns_assumed":POINTER_READOUT_NS,
            "loss_channel":loss,
            "mean_F_e_conditional_on_no_loss":float(fidelity.mean()),
            "mean_F_e_unconditional_including_flagged_loss":
                float(loss_eta*fidelity.mean()),
            "Monte_Carlo_standard_error_conditional":float(fidelity.std(ddof=1)/np.sqrt(shots))
                    if shots>1 else None,
            "mean_final_occupation_leakage_conditional":
                float(np.mean([v["final_occupation_leakage_conditional"]
                               for v in trajectories])),
            "mean_final_Gauss_acceptance_conditional":gauss_av,
            "individual_F_e_conditional":fidelity.tolist()}


def run():
    W=physical_isometry()
    K,_,_=code_isometry(W)
    output=[]
    # Coherent feedback with NO environmental noise at three interaction strengths
    for U in (20.,30.,40.):
        opt=optimize_pulse(U)
        output.append(process_sample(K,U,opt["active_detuning"],
            opt["time"],gamma=0.,shots=1,seed=int(U),substeps=1))
    # Same physical U/t=30 pulse schedule, different dephasing operators.
    opt=optimize_pulse(30.)
    for i,model in enumerate(("orbital_charge","colour_resolved")):
        output.append(process_sample(K,30.,opt["active_detuning"],
            opt["time"],gamma=.0002,model=model,shots=6,
            seed=126+i,substeps=2))
    # Benchmark against old optimistic PR56 feedback, at U/t=30.
    old=exact_full_channel(K,30.,opt["active_detuning"],opt["time"])
    return {"model":"Finite S3, 4D phase code, persistent 495D four-fermion Hubbard dynamics",
        "added_over_PR56":[
            "3 microscopic pointer-conditioned 495D Hubbard feedback pulses and 2 link transports (5 operations)",
            "CPTP occupation-number dephasing during pulses, transport and readout",
            "exact uniform active-fermion particle-loss versus no-loss flagged-erasure instrument"],
        "PR56_ideal_feedback_reference_U30":old,
        "loss_model":{"operators":"sqrt(kappa)c_{orbital_colour}, kappa=1/T1; equal for all 12 modes",
                     "N_particles":4,"joint_no_jump_amplitude":
                        "exp[-(4 kappa total_duration)/2]",
                     "eventual_fermion_loss_is_flagged_erasure":True,
                     "explicit_N3_post_jump_dynamics":False},
        "runs":output,
        "limitations":[
            "Local 15D Hubbard superexchange is embedded in persistent 495D four-particle register; calibrated fermionic device coupling is not available",
            "Two-body microscopic local exchanges only; S3 gauge-link charge transport remains exact high-level unitary, not native Hamiltonian",
            "Loss jump trajectories are flagged erasures with zero decoded fidelity, not explicitly evolved 3-,2-,1-fermion states; no active loss correction or repumping",
            "Gaussian stochastic dephasing is a trajectory estimate for a specified Lindblad channel; finite trajectories and Trotter steps introduce sampling/time-step uncertainty",
            "Only the four ACTIVE fermions suffer loss, while spectator B1/reference and link-pointer hardware loss are unmodeled",
            "Pointer Fourier/readout assumed exact; its duration contributes idle dephasing and particle loss but no assignment-error model",
            "No true physical local-SU3 link, continuous gauge group, original nine-quark 21D / 65-syndrome readout, native pulses, fault threshold or unencoded matched benchmark",
            "Hardware units t/h=10MHz, link transport 100ns, readout 300ns, and T1=30us are hypothetical"]}


def main():
    result=run()
    Path("colour_su3_hubbard_loss_dephase_feedback_results.json").write_text(
        json.dumps(result,indent=2,sort_keys=True)+"\n")
    for r in result["runs"]:
        print("FULL_HUBBARD_FEEDBACK_NOISE",json.dumps({
            k:v for k,v in r.items() if k!="individual_F_e_conditional"
        },sort_keys=True),flush=True)


if __name__=="__main__":main()
