"""Device-parameterized two-transmon CP22 pulse and physical noise probe.

A four-level Duffing-plus-exchange pair with a smooth *detuning* excursion
is propagated exactly (16x16). A diagonal product-Z frame is extracted and
only physically local phases are removed. The pulse is NOT a calibrated or
guaranteed selective CP22; the reported residual error is the result.

Amplitude damping is the exact four-level bosonic pure-loss CPTP channel.
Dephasing is an exact Gaussian number-phase Schur CPTP channel. They act
AFTER the modeled coherent pulse (lumped-channel approximation, not full
pulse-time Lindblad). Nothing here is calibrated to a particular device.

Primary context: Phys Rev A 108 022614; PRX Quantum 4 030327;
Phys Rev Lett 134 050601; npj Quantum Information 10 2024 00892-z.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass, replace, asdict
from math import comb
import json
from pathlib import Path
import numpy as np
from scipy.linalg import expm

D = 4
PAIR_DIM = D * D
COMP = np.array([D*i+j for i in range(3) for j in range(3)], int)
PROJECT = np.eye(PAIR_DIM, dtype=complex)[:, COMP]
NUMBER = np.diag(np.arange(D))
ANNIHILATOR = np.zeros((D,D), complex)
for n in range(1,D):
    ANNIHILATOR[n-1,n] = np.sqrt(n)
EYE = np.eye(D, dtype=complex)
EYE_PAIR = np.eye(PAIR_DIM, dtype=complex)
TOTAL_EXCITATION = np.kron(NUMBER,EYE)+np.kron(EYE,NUMBER)


@dataclass(frozen=True)
class PairParameters:
    """All entries illustrative; frequencies in MHz, time in ns."""
    delta_idle_mhz: float = 650.
    delta_pulse_mhz: float = 290.
    anharmonicity_a_mhz: float = -230.
    anharmonicity_b_mhz: float = -260.
    exchange_mhz: float = 9.
    duration_ns: float = 120.
    steps: int = 48
    t1_ns: float = 30000.
    tphi_ns: float = 45000.


def duffing_hamiltonian(p: PairParameters, detuning_mhz: float) -> np.ndarray:
    """Rotating frame of transmon A: -i 2*pi*1e-3 H_MHz t_ns."""
    n = NUMBER
    anh = n @ (n-EYE)/2
    a = ANNIHILATOR
    h = (p.anharmonicity_a_mhz*np.kron(anh,EYE)
         + np.kron(EYE,detuning_mhz*n+p.anharmonicity_b_mhz*anh)
         + p.exchange_mhz*(np.kron(a.conj().T,a)+np.kron(a,a.conj().T)))
    assert np.linalg.norm(h-h.conj().T) < 1e-12
    assert np.linalg.norm(h@TOTAL_EXCITATION-TOTAL_EXCITATION@h)<1e-10
    return h


def smooth_pulse(p: PairParameters) -> np.ndarray:
    """A detuning pulse followed by ideal removal of identical idle evolution.

    The idle inverse is a REFERENCE FRAME diagnostic, not an executable free
    operation in a compiled quantum circuit.
    """
    if p.steps < 2 or p.duration_ns <= 0:
        raise ValueError("Positive pulse time and at least 2 time steps required")
    h_idle = duffing_hamiltonian(p,p.delta_idle_mhz)
    dt = p.duration_ns/p.steps
    u = EYE_PAIR.copy()
    for step in range(p.steps):
        t = (step+.5)/p.steps
        delta = p.delta_idle_mhz+(p.delta_pulse_mhz-p.delta_idle_mhz)*np.sin(np.pi*t)**2
        h = duffing_hamiltonian(p,delta)
        u = expm(-2j*np.pi*.001*dt*h) @ u
    idle = expm(-2j*np.pi*.001*p.duration_ns*h_idle)
    u = idle.conj().T@u
    assert np.linalg.norm(u.conj().T@u-EYE_PAIR)<1e-10
    assert np.linalg.norm(u@TOTAL_EXCITATION-TOTAL_EXCITATION@u)<1e-9
    return u


def local_phase_frame(u: np.ndarray):
    """Remove row/column phases inferable from |i0> and |0j> amplitudes.

    Returns the compensated 16D matrix and the remaining correlated |22>
    phase. Does NOT use a nonlocal correction or artificially remove leakage.
    """
    angle = np.angle(np.diag(u)).reshape((D,D))
    first = angle[:,0]-angle[0,0]
    second = angle[0,:]
    correction = np.exp(-1j*(first[:,None]+second[None,:])).reshape(-1)
    compensated = correction[:,None]*u
    realized = float(np.angle(compensated[2*D+2,2*D+2]))
    return compensated, realized


def ideal_cp22(angle: float) -> np.ndarray:
    u = EYE_PAIR.copy()
    u[2*D+2,2*D+2] = np.exp(1j*angle)
    return u


def loss_kraus(duration_ns: float, t1_ns: float) -> list[np.ndarray]:
    """Bosonic pure loss: |n> -> |n-l> at eta=exp(-t/T1), n<=3."""
    if duration_ns < 0 or t1_ns <= 0:
        raise ValueError("Invalid loss time")
    eta = float(np.exp(-duration_ns/t1_ns))
    k = []
    for jumps in range(D):
        op = np.zeros((D,D), complex)
        for n in range(jumps,D):
            op[n-jumps,n] = (np.sqrt(comb(n,jumps))
                * (1-eta)**(jumps/2)*eta**((n-jumps)/2))
        k.append(op)
    assert np.linalg.norm(sum((a.conj().T@a for a in k),np.zeros((D,D),complex))-EYE)<1e-11
    return k


def dephasing_kraus(duration_ns: float, tphi_ns: float) -> list[np.ndarray]:
    """Exact Gaussian phase diffusion: rho_nm -> exp[-t(n-m)^2/Tphi] rho_nm."""
    if duration_ns < 0 or tphi_ns <= 0:
        raise ValueError("Invalid dephasing time")
    indices = np.arange(D)
    schur = np.exp(-duration_ns/tphi_ns*(indices[:,None]-indices[None,:])**2)
    eig,vec = np.linalg.eigh(schur)
    if np.min(eig)<-1e-10:
        raise AssertionError("Dephasing Schur kernel is not positive")
    kraus = [np.diag(np.sqrt(max(0.,eig[i]))*vec[:,i]) for i in range(D)]
    assert np.linalg.norm(sum((k.conj().T@k for k in kraus),np.zeros((D,D),complex))-EYE)<1e-11
    return kraus


def single_site_noise(duration_ns: float, t1_ns: float, tphi_ns: float):
    return [phase@loss for phase in dephasing_kraus(duration_ns,tphi_ns)
            for loss in loss_kraus(duration_ns,t1_ns)]


def two_site_noise(duration_ns: float, t1_ns: float, tphi_ns: float):
    single = single_site_noise(duration_ns,t1_ns,tphi_ns)
    return [np.kron(a,b) for a in single for b in single]


def channel_metrics(unitary: np.ndarray, target_angle: float,
                    noise: list[np.ndarray] | None=None):
    """Entanglement fidelity including leakage (no postselection).

    The target is a 16x9 ideal isometry. Both coherence and leakage from
    computational levels are included in the unconditioned process overlap.
    """
    if noise is None:
        noise=[EYE_PAIR]
    target = ideal_cp22(target_angle)@PROJECT
    coherent = 0.
    inspace = 0.
    norm = 0.
    for k in noise:
        acted=k@unitary@PROJECT
        overlap=target.conj().T@acted
        coherent+=abs(np.trace(overlap))**2
        inspace+=float(np.vdot(PROJECT.conj().T@acted,PROJECT.conj().T@acted).real)
        norm+=float(np.vdot(acted,acted).real)
    assert abs(norm-9)<2e-8, norm
    fe=float(coherent/81)
    leakage=float(max(0.,1-inspace/9))
    assert -1e-10<=fe<=1+1e-10
    return {'unconditional_entanglement_fidelity':fe,
            'final_fourth_level_leakage_probability':leakage,
            'logical_survival_probability':float(inspace/9)}


def pulse_candidate(p: PairParameters, target=np.pi, include_noise=False):
    unit=smooth_pulse(p)
    calibrated,realized=local_phase_frame(unit)
    ideal=channel_metrics(calibrated,target)
    at_realized=channel_metrics(calibrated,realized)
    comp=PROJECT.conj().T@calibrated@PROJECT
    diagonal=np.angle(np.diag(comp))
    spectator=float(max(abs(np.angle(np.exp(1j*diagonal[i])))
                        for i in range(9) if i != 8))
    result={'duration_ns':p.duration_ns,'delta_pulse_mhz':p.delta_pulse_mhz,
            'realized_cp22_angle_rad':realized,
            'target_cp22_angle_rad':float(target),
            'coherent_target_F_e':ideal['unconditional_entanglement_fidelity'],
            'coherent_realized_angle_F_e':at_realized['unconditional_entanglement_fidelity'],
            'coherent_leakage_probability':ideal['final_fourth_level_leakage_probability'],
            'max_other_computational_diagonal_phase_rad':spectator,
            'physical_time_model':'smooth detuning pulse plus nonphysical idle-frame subtraction'}
    if include_noise:
        noisy=channel_metrics(calibrated,target,two_site_noise(p.duration_ns,p.t1_ns,p.tphi_ns))
        result.update({'post_pulse_T1_Tphi_target_F_e':noisy['unconditional_entanglement_fidelity'],
                       'post_pulse_T1_Tphi_final_leakage_probability':noisy['final_fourth_level_leakage_probability'],
                       't1_ns_scenario':p.t1_ns,'tphi_ns_scenario':p.tphi_ns,
                       'noise_composition':'exact CPTP local channels LUMPED AFTER pulse'})
    return result


def study(p: PairParameters, durations=(40.,80.,160.), detunings=(260.,290.,320.)):
    grid=[]
    for duration in durations:
        for detuning in detunings:
            trial=replace(p,duration_ns=duration,delta_pulse_mhz=detuning)
            grid.append(pulse_candidate(trial))
    ranked=sorted(grid,key=lambda r:r['coherent_target_F_e'],reverse=True)
    best=ranked[0]
    best_p=replace(p,duration_ns=best['duration_ns'],
                   delta_pulse_mhz=best['delta_pulse_mhz'])
    return {'model':'four-level two-Duffing-transmon detuning excursion plus exchange',
            'illustrative_parameters':asdict(p),
            'target_cp22_angle_rad':float(np.pi),
            'candidate_count':len(grid),'candidates':grid,
            'best_from_limited_grid':pulse_candidate(best_p,include_noise=True),
            'important_limitation':('This grid does not synthesize arbitrary-angle '
               'selective CP22. A high-fidelity target requires calibrated control '
               'and physical pulse optimization; even a promising grid point is '
               'not a characterized native two-transmon gate.')}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',default='su3_transmon_noise_probe_results.json')
    args=parser.parse_args()
    data=study(PairParameters())
    Path(args.output).write_text(json.dumps(data,indent=2,sort_keys=True)+'\n')
    print('DEVICE_PULSE_PROBE',json.dumps(data['best_from_limited_grid'],sort_keys=True),flush=True)
    print('PHYSICAL_CPTP_TWO_TRANSMON_PROBE_PASS',flush=True)
    return data


if __name__=='__main__':
    main()
