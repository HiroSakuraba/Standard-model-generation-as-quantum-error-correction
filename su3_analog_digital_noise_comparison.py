"""SU3 full-circuit digital vs conditional analog cubic-clock comparison.

Analog arm: actual five-plaquette native H remains on during finite cubic
electric pulses, using existing repository implementation. Digital arm:
expanded five-qutrit CP22 circuit from PR39. A MHz conversion assumes
H/h = f_MHz * H_dimensionless; there is no measured device calibration.
There is NO full 4^5 T1/Tphi simulation or compiled noisy decoder here.
"""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path
import numpy as np
from scipy.sparse.linalg import expm_multiply
from su3_cubic_casimir_echo import construct, entanglement_fidelity
from su3_native_control_compilation import (
    lab_cubic_pulses, toggled_matrices, evolution, THETA, T
)
from su3_channel_recovery import syndrome_decoder, source_data
from su3_qutrit_gate_compiler import MASK, CYCLES
from su3_full_refocusing_cp22 import propagate, toy_risk


def dimensional_mapping(scale_mhz=0.01, tau=.002):
    """Hypothetical conversion: H_phys/h [MHz] = scale_mhz * H_dimless."""
    if scale_mhz<=0 or tau<=0:raise ValueError('positive scale and pulse width')
    factor_ns=1/(2*np.pi*.001*scale_mhz)
    n_pulses=4*CYCLES
    return {'scale_mhz_assumed':scale_mhz,
            'pulse_width_dimensionless':tau,
            'pulse_width_ns':float(tau*factor_ns),
            'total_analog_elapsed_us':float((T+n_pulses*tau)*factor_ns/1000),
            'native_no_control_elapsed_us':float(T*factor_ns/1000),
            'clock_drive_strength_mhz_per_addressed_qutrit':float(THETA*scale_mhz/tau),
            'control_pulses_global':n_pulses,
            'addressed_face_clock_operations':n_pulses*len(MASK),
            'warning':'Assumed spectral mapping, NOT a calibrated gauge simulator.'}


def syndrome_scenarios(readout_ns=800.,reset_ns=300.,p_readout=.001,
                       assumed_extraction_gates=(0,24,128)):
    """Resource assumptions; actual six-outcome instrument is NOT synthesized."""
    assert readout_ns >= 0 and reset_ns>=0 and 0<=p_readout<1
    n_outcomes=6
    n_pointers=2
    assert 3**n_pointers >= n_outcomes
    assert 2**2 < n_outcomes <= 2**3
    rows=[]
    for gates in assumed_extraction_gates:
        rows.append({
            'assumed_extraction_entanglers':gates,
            'two_pointer_readouts':n_pointers,
            'minimum_classical_bits':3,
            'minimum_distinct_syndrome_outcomes':n_outcomes,
            'optimistic_pointer_cycle_us':float((readout_ns+reset_ns)/1000),
            'pointer_readout_event_free':float((1-p_readout)**n_pointers),
            'not_implemented':'This is a scenario, not a syndrome-extraction circuit.'})
    return rows


def analog_flagged_replacement(f_noiseless, p_local_pulse, readout_p=.001):
    """NOT the physical four-level decoherence model: logical CPTP toy only."""
    if not 0<=p_local_pulse<1 or not 0<=readout_p<1:
        raise ValueError('probabilities')
    count=4*CYCLES*len(MASK)
    survival=float((1-p_local_pulse)**count*(1-readout_p)**2)
    return {'clock_action_fault_probability':p_local_pulse,
            'clock_operations':count,'pointer_fault_probability_each':readout_p,
            'event_free_probability':survival,
            'logical_replacement_F_e':float(survival*f_noiseless+(1-survival)*.25),
            'model':'CPTP logical replacement risk only, NOT physical noise'}


def run(tau=.002,scale_mhz=.01):
    states,K,d,W,H=construct()
    decoder,_,gram=syndrome_decoder(*source_data())
    assert gram.rank()==5 and len(MASK)==4
    uncontrolled=expm_multiply(-1j*T*H,K,traceA=0)
    uncontrolled_f=entanglement_fidelity([uncontrolled],K,decoder)
    ideal_analog=evolution(K,toggled_matrices(d,W,MASK),T,CYCLES)
    finite_analog=lab_cubic_pulses(states,H,K,MASK,CYCLES,tau)
    analog_f=entanglement_fidelity([finite_analog],K,decoder)
    ideal_analog_f=entanglement_fidelity([ideal_analog],K,decoder)
    digital=[]
    for approximate in (False,True):
        y,counts=propagate(K,d,1,approximate=approximate)
        f=entanglement_fidelity([y],K,decoder)
        digital.append({'variant':'approximate_32_CP22' if approximate else 'exact_40_CP22',
                        'recovered_F_e_no_gate_noise':f,**counts,
                        'assumed_serial_gate_noise_example':
                          toy_risk(f,uncontrolled_f,counts,1e-5,1e-5,
                                   readout_us=.8,p_readout=.001)})
    return {'model':'separate hypothetical digital and analog realizations of same truncated five-plaquette model',
            'gauge_model_state_dimension':243,'control_mask':list(MASK),
            'native_uncontrolled_noiseless_recovered_F_e':uncontrolled_f,
            'ideal_instantaneous_analog_recovered_F_e':ideal_analog_f,
            'finite_pulse_analog_recovered_F_e':analog_f,
            'finite_pulse_analog_raw_F_e':entanglement_fidelity([finite_analog],K),
            'analog_time_conversion':dimensional_mapping(scale_mhz,tau),
            'digital':digital,
            'syndrome_extraction':syndrome_scenarios(),
            'analog_flagged_replacement_scenarios':[
                analog_flagged_replacement(analog_f,p) for p in (0,1e-4,1e-3,1e-2)],
            'limitations':[
                'Uncalibrated cubic-Casimir control and native SU3 Wilson Hamiltonian',
                'Analog elapsed time based on arbitrary MHz energy scale',
                'Digital CP22/local operation times are hypothetical serial schedules',
                'Device pulse probe is only on a two-transmon 4x4 pair; it is NOT applied to full five-plaquette run',
                'No full-chain physical amplitude damping, dephasing, leakage or CPTP noisy syndrome extraction',
                'Analog and digital processor assumptions may be mutually incompatible',
                'Uncontrolled storage benchmark is noiseless, not a same-device noisy comparison']}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--tau',type=float,default=.002)
    parser.add_argument('--scale-mhz',type=float,default=.01)
    parser.add_argument('--output',default='su3_analog_digital_noise_comparison_results.json')
    args=parser.parse_args()
    result=run(args.tau,args.scale_mhz)
    Path(args.output).write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print('ANALOG_DIGITAL_REFERENCE',json.dumps({
        'uncontrolled':result['native_uncontrolled_noiseless_recovered_F_e'],
        'ideal_analog':result['ideal_instantaneous_analog_recovered_F_e'],
        'finite_analog':result['finite_pulse_analog_recovered_F_e'],
        'exact_digital':result['digital'][0]['recovered_F_e_no_gate_noise'],
        'approx_digital':result['digital'][1]['recovered_F_e_no_gate_noise'],
        'analog_time':result['analog_time_conversion']},sort_keys=True),flush=True)
    print('SEPARATE_DIGITAL_AND_ANALOG_BUDGETS_PASS',flush=True)
    return result


if __name__=='__main__':
    main()
