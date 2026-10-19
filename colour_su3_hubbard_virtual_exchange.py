"""Microscopic SU(3)-symmetric two-orbital Hubbard control of colour exchange.

Two spinless fermions with three colours on TWO ORBITALS AT ONE GAUGE
VERTEX: 6 fermionic modes, choose(6,2)=15 physical Fock states.
9 states have one particle/orbital; 6 have an empty and doubly occupied
orbital. The 15x15 Hamiltonian exactly propagates virtual doublons.

H_p = -t sum_colour (c_L† c_R + h.c.)
      + (U/2) sum_site n_site(n_site-1)
      + (delta_p/2)(n_L-n_R).

The neutral qutrit pointer shifts the SITE NUMBER difference: delta_0=0,
delta_1=delta_2=delta. This is a two-register pointer+orbital-occupation
coupling, not a native pointer-controlled colour exchange.

In the single-occupation subspace the second-order effective H is
-J(delta) (I-SWAP)/2 with
J(delta)=2t²[1/(U+delta)+1/(U-delta)].
Choose a pulse with J(0) T=2pi and J(delta) T=3pi: off-pointer ~I,
on-pointer ~SWAP. Exact Fock dynamics reveal finite-U errors and doublon
leakage. SU3 colour invariance holds for both orbitals sharing ONE
gauge vertex, not for independent gauge vertices without a link.

No full charged-bus circuit, no calibrated hardware, no actual fermionic
gauge-link transport or physical fault-tolerant correction is claimed.
"""
from __future__ import annotations

from functools import lru_cache
from itertools import combinations
from pathlib import Path
import json

import numpy as np
from scipy.linalg import eigh
from scipy.optimize import minimize

from colour_singlet_code import gellmann

T_COL=3
ORBITALS=2
MODES=6
TWO_FERMION_DIM=15
COMPUTATIONAL=9

FOCK=tuple(sum(1<<k for k in pair) for pair in combinations(range(MODES),2))
INDEX={s:i for i,s in enumerate(FOCK)}
W=np.zeros((TWO_FERMION_DIM,COMPUTATIONAL),complex)
for a in range(3):
    for b in range(3):
        W[INDEX[(1<<a)|(1<<(b+3))],3*a+b]=1

SWAP=np.zeros((9,9),complex)
for a in range(3):
    for b in range(3):
        SWAP[3*b+a,3*a+b]=1
Q_ANTISYM=(np.eye(9)-SWAP)/2


def creation_annihilation(dst:int,src:int,state:int):
    """Matrix sign of c†_dst c_src in ascending-mode Fock order."""
    if not (state>>src)&1 or (state>>dst)&1:
        return None
    sign=(-1)**((state&((1<<src)-1)).bit_count())
    emptied=state^(1<<src)
    sign*=(-1)**((emptied&((1<<dst)-1)).bit_count())
    return emptied|(1<<dst),sign


def one_body(dst,src):
    if dst==src:
        return np.diag([float((state>>src)&1) for state in FOCK]).astype(complex)
    result=np.zeros((15,15),complex)
    for j,s in enumerate(FOCK):
        term=creation_annihilation(dst,src,s)
        if term is not None:
            state,sign=term
            result[INDEX[state],j]=sign
    return result


def _building_blocks():
    hopping=np.zeros((15,15),complex)
    interaction=np.zeros((15,15),complex)
    number_imbalance=np.zeros((15,15),complex)
    for a in range(3):
        hopping-=one_body(a,a+3)+one_body(a+3,a)
    for i,s in enumerate(FOCK):
        left=sum((s>>k)&1 for k in range(3))
        right=2-left
        interaction[i,i]=int(left==2 or right==2)
        number_imbalance[i,i]=(left-right)/2
    assert np.linalg.norm(hopping-hopping.conj().T)<1e-13
    assert int(np.trace(interaction).real)==6
    assert np.linalg.norm(W.conj().T@W-np.eye(9))<1e-12
    return hopping,interaction,number_imbalance


HOP,DOUBLON,IMBALANCE=_building_blocks()


def hamiltonian(U,delta,t=1.):
    if not all(np.isfinite(x) for x in (U,delta,t)) or U<=0 or t<=0 or abs(delta)>=U:
        raise ValueError("require finite U>0, t>0, abs(delta)<U")
    return t*HOP+U*DOUBLON+delta*IMBALANCE


def su3_total_generator(T):
    return sum((T[a,b]*(one_body(a,b)+one_body(a+3,b+3))
                for a in range(3) for b in range(3)),
               np.zeros((15,15),complex))


def invariant_residual(U=20.,delta=0.,t=1.):
    H=hamiltonian(U,delta,t)
    return max(float(np.linalg.norm(H@su3_total_generator(T)
                                  -su3_total_generator(T)@H))
               for T in gellmann())


def effective_exchange(U,delta,t=1.):
    """Exact second-order Schrieffer-Wolff block and analytic SU3 value."""
    H=hamiltonian(U,delta,t)
    high=np.flatnonzero(np.diag(DOUBLON).real>.5)
    B=H[np.ix_(high,range(15))]@W
    E=np.diag(H)[high]
    Heff=-B.conj().T@np.diag(1/E)@B
    J=2*t*t*(1/(U+delta)+1/(U-delta))
    error=float(np.linalg.norm(Heff+J*Q_ANTISYM))
    assert error<1e-10
    return J,Heff,error


def projected_evolution(U,delta,time,t=1.):
    """Physical 15D unitary; return 9x9 computational projection, leakage."""
    if not np.isfinite(time) or time<0:raise ValueError("duration")
    e,V=eigh(hamiltonian(U,delta,t))
    projected=W.conj().T@((V*np.exp(-1j*e*time))@V.conj().T)@W
    leak=float(1-np.linalg.norm(projected,"fro")**2/9)
    return projected,leak


def gate_metrics(U,time,delta,t=1.):
    """27D coherent pointer-controlled isometry entanglement fidelity.

    Target = diag(I,SWAP,SWAP); no projection/postselection correction
    for escaped physical doublon states.
    """
    if time<=0 or abs(delta)>=U:return {"entanglement_fidelity":0.}
    idle,l0=projected_evolution(U,0.,time,t)
    active,l1=projected_evolution(U,delta,time,t)
    overlap=np.trace(idle)+2*np.trace(SWAP.conj().T@active)
    Fe=float(abs(overlap/27)**2)
    leak=float((l0+2*l1)/3)
    return {"entanglement_fidelity":Fe,
            "average_final_doublon_leakage":leak,
            "off_pointer_leakage":l0,
            "on_pointer_leakage":l1,
            "off_I_error_fro_per_sqrt9":float(np.linalg.norm(idle-np.eye(9))/3),
            "on_swap_error_fro_per_sqrt9":float(np.linalg.norm(active-SWAP)/3)}


def perturbative_schedule(U,t=1.):
    """Relative phase 2pi off, 3pi on. δ/U = sqrt(1/3)."""
    d=U/np.sqrt(3)
    J0=effective_exchange(U,0,t)[0]
    J1=effective_exchange(U,d,t)[0]
    T=2*np.pi/J0
    assert abs(J1/J0-1.5)<1e-12
    return {"time_in_inverse_t_angular_units":T,
            "on_pointer_detuning":d,
            "off_exchange_rate":J0,
            "on_exchange_rate":J1,
            "on_over_off_exchange_ratio":J1/J0}


def optimize_pulse(U,t=1.):
    """Bounded two-variable fit against EXACT 15D Fock evolution.

    Fixed 1-pulse rectangular hopping profile, δ_off=0. Optimized
    only T and δ_on within bounds near perturbative 2pi/3pi proposal.
    """
    proposal=perturbative_schedule(U,t)
    T0=proposal["time_in_inverse_t_angular_units"]
    d0=proposal["on_pointer_detuning"]
    def loss(par):
        T,d=par
        return 1-gate_metrics(U,T,d,t)["entanglement_fidelity"]
    r=minimize(loss,[T0,d0],method="L-BFGS-B",
        bounds=[(.8*T0,1.2*T0),(.2*U,.88*U)],
        options={"maxiter":140,"ftol":1e-12})
    T,d=map(float,r.x)
    return {"time":T,"active_detuning":d,
            "optimized_success":bool(r.success),
            "objective_function":"1 - exact leakage-inclusive controlled-SWAP entanglement fidelity",
            "baseline_metrics":gate_metrics(U,T0,d0,t),
            "optimized_metrics":gate_metrics(U,T,d,t),
            "fractional_time_vs_perturbative":T/T0}


def max_transient_leakage(U,delta,time,t=1.,steps=41):
    """Doublon occupancy at sampled intermediate points, averaged over 9 inputs."""
    values=[projected_evolution(U,delta,float(f)*time,t)[1]
            for f in np.linspace(0,1,steps)]
    return {"number_time_samples":steps,
            "max_transient_doublon_occupation":float(max(values)),
            "mean_transient_doublon_occupation":float(np.mean(values)),
            "final_doublon_occupation":float(values[-1])}


def illustrative_loss_budget(Fe,time,t_over_h_MHz=10.,T1_us=30.):
    """Erasure/loss-only multiplicative survival scenario for two particles.

    Angular model parameter t = 2pi * (t/h in MHz) * .001 rad/ns.
    Any fermion loss is treated as flagged erasure, never corrected;
    phenomenological independent loss, NOT a Hubbard Lindblad solve.
    """
    if t_over_h_MHz<=0 or T1_us<=0:raise ValueError("positive physical assumptions")
    duration_ns=1000*time/(2*np.pi*t_over_h_MHz)
    survival=float(np.exp(-2*duration_ns/(1000*T1_us)))
    return {"hypothetical_t_over_h_MHz":t_over_h_MHz,
            "hypothetical_two_fermion_T1_us":T1_us,
            "pulse_duration_ns":float(duration_ns),
            "two_particle_no_loss_probability":survival,
            "leakage_and_loss_survival_weighted_F_e":float(Fe*survival),
            "scope":"Flagged independent particle erasure estimate; no physical dissipative evolution, Tphi or encoded recovery"}


def run():
    invariance=[invariant_residual(20.,d) for d in (0.,4.,10.)]
    assert max(invariance)<1e-10
    arms=[]
    for U in (20.,30.,40.):
        p=perturbative_schedule(U)
        optimum=optimize_pulse(U)
        T=optimum["time"]; d=optimum["active_detuning"]
        leakage={
            "pointer_off":max_transient_leakage(U,0.,T),
            "pointer_on":max_transient_leakage(U,d,T)}
        arms.append({"U_over_t":U, "analytic":p, "optimization":optimum,
            "transient_doublon":leakage,
            "loss_scenario":illustrative_loss_budget(
                optimum["optimized_metrics"]["entanglement_fidelity"],T)})
    return {"basis":"two identical fermions, 6 single-particle orbitals, exactly choose(6,2)=15",
        "one_per_orbital_colour_subspace_dimension":9,
        "virtual_doublon_hole_dimension":6,
        "neutral_pointer_dimension":3,
        "SU3_global_same_vertex_commutator_residuals":invariance,
        "second_order_effective_exchange_identity":"H_eff = -2t^2*(1/(U+delta)+1/(U-delta))*(I-SWAP)/2",
        "pointer_onsite_detuning":"delta_p/2 (n_L - n_R), delta_off=0, delta_on for pointer=1,2",
        "benchmark_arms":arms,
        "limitations":[
            "Both fermionic orbitals are COLOCATED under ONE vertex SU3 Gauss action: colour-blind hopping is NOT invariant under two independent vertex gauges without a link",
            "No dynamical gauge link, no mobile-charged-bus end-to-end physical syndrome, no nine-quark 21D/65-syndrome recovery",
            "Second-order superexchange is perturbative; the 15-state pulse fidelities include all doublon states and are NOT obtained from H_eff",
            "Each pointer condition changes orbital detuning with spin-blind density coupling; Hubbard interaction itself is two-body in fermionic operators",
            "No ramp-on/off, pointer dispersive dressing calibration, Tphi/Lindblad pulse-time channel, particle loss dynamics, transport or reset",
            "The hypothetical T1 erasure survival estimate is not a simulated physical noise channel",
            "No calibrated fermion orbital hopping/interaction rates or actual quantum-memory advantage"]}


def main():
    out=run()
    Path("colour_su3_hubbard_virtual_exchange_results.json").write_text(
        json.dumps(out,indent=2,sort_keys=True)+"\n")
    for x in out["benchmark_arms"]:
        print("MICROSCOPIC_SU3_HUBBARD",json.dumps({
            "U_over_t":x["U_over_t"],
            "perturbative_gate_F_e":x["optimization"]["baseline_metrics"]["entanglement_fidelity"],
            "optimized_gate_F_e":x["optimization"]["optimized_metrics"]["entanglement_fidelity"],
            "optimized_final_leakage":x["optimization"]["optimized_metrics"]["average_final_doublon_leakage"],
            "on_max_transient_doublon":x["transient_doublon"]["pointer_on"]["max_transient_doublon_occupation"],
            "hypothetical_duration_ns":x["loss_scenario"]["pulse_duration_ns"],
            "loss_weighted_F_e":x["loss_scenario"]["leakage_and_loss_survival_weighted_F_e"],
        },sort_keys=True))
    print("SU3_HUBBARD_MICROSCOPIC_TWO_ORBITAL_PROBE_PASSED")


if __name__=="__main__":
    main()
