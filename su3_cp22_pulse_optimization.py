"""Bounded synthetic optimization of the PR40 four-level Duffing CP22 pulse.

Optimizes (duration, pulse detuning) in a smooth sin² detuning excursion.
The existing 16x16 exact unitary and locally-removable phase decomposition
are reused. Objective includes unconditional process entanglement fidelity,
so leakage/spectator phases are penalized automatically. A scan at target
pi/2 demonstrates the difference between target angles; it does NOT
certify a single arbitrary-angle native gate family.

Differential evolution is deterministic with seed; small budget. An
independent finer step discretization is required before trusting results.
Neither the pulse shape nor the parameters represent a calibrated chip.
"""
from __future__ import annotations
import argparse
import json
from dataclasses import asdict,replace
from pathlib import Path
import numpy as np
from scipy.optimize import differential_evolution
from su3_transmon_noise_probe import PairParameters,pulse_candidate

SEED=20261019


def optimize_angle(angle, iterations=7, steps=24, seed=SEED):
    base=PairParameters()
    def evaluate(x):
        duration, detuning=map(float,x)
        pars=replace(base,duration_ns=duration,
                     delta_pulse_mhz=detuning,steps=steps)
        return 1-pulse_candidate(pars,target=angle)['coherent_target_F_e']
    bounds=[(32.,200.),(210.,420.)]
    result=differential_evolution(evaluate,bounds,
         seed=seed,popsize=5,maxiter=iterations,polish=False,
         tol=1e-3,updating='immediate',workers=1)
    old=replace(base,duration_ns=80.,delta_pulse_mhz=290.,steps=steps)
    baseline=pulse_candidate(old,target=angle)
    candidate=replace(base,duration_ns=float(result.x[0]),
                      delta_pulse_mhz=float(result.x[1]),steps=steps)
    low=pulse_candidate(candidate,target=angle,include_noise=True)
    high=pulse_candidate(replace(candidate,steps=2*steps),target=angle,include_noise=True)
    baseline_high=pulse_candidate(replace(old,steps=2*steps),target=angle,include_noise=True)
    # Never call an optimizer output an improvement if the incumbent wins
    # under the SAME fine time-discretization.
    improved=high['coherent_target_F_e']>baseline_high['coherent_target_F_e']
    chosen_low=low if improved else baseline
    chosen_high=high if improved else baseline_high
    assert chosen_high['coherent_target_F_e']>=baseline_high['coherent_target_F_e']-1e-12
    return {
      'target_angle_radians':float(angle),
      'search_bounds_duration_ns':[32.,200.],
      'search_bounds_detuning_mhz':[210.,420.],
      'candidate_baseline':baseline,
      'coarse_step_optimization':low,
      'optimizer_finer_step_evaluation':high,
      'baseline_finer_step_evaluation':baseline_high,
      'optimizer_improved_over_baseline':bool(improved),
      'finer_step_evaluation':chosen_high,
      'coarse_minus_fine_fidelity':float(chosen_low['coherent_target_F_e']-chosen_high['coherent_target_F_e']),
      'candidate_evaluations':int(result.nfev),
      'search_seed':seed,
      'iterations':iterations,
      'coarse_time_steps':steps,
      'certifies_arbitrary_angle_gate_family':False,
      'limitation':'Small deterministic numerical search, exact 16D unitary for a synthetic pair; no calibrated waveform or leakage suppression hardware.'}


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--iterations',type=int,default=7)
    ap.add_argument('--steps',type=int,default=24)
    ap.add_argument('--output',default='su3_cp22_pulse_optimized_results.json')
    args=ap.parse_args()
    if args.steps<8 or args.iterations<1:raise ValueError('small step count/iterations')
    rows=[optimize_angle(angle,args.iterations,args.steps,SEED+i)
          for i,angle in enumerate((np.pi/2,np.pi))]
    out={'model':'four-level exchange-coupled Duffing transmons; single smooth detuning control, synthetic parameters',
         'optimizer':'bounded scipy differential evolution; explicit fine-discretization crosscheck',
         'base_parameters':asdict(PairParameters()),
         'rows':rows,
         'warning':'Even a substantial simulated local improvement cannot establish any calibrated arbitrary-angle CP22 or guarantee that 10,080 gates form a functioning memory.'}
    Path(args.output).write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    for row in rows:print('PULSE_OPT',json.dumps({
       'target_angle':row['target_angle_radians'],
       'baseline_F':row['candidate_baseline']['coherent_target_F_e'],
       'optimized_F':row['finer_step_evaluation']['coherent_target_F_e'],
       'duration_ns':row['finer_step_evaluation']['duration_ns'],
       'detuning_mhz':row['finer_step_evaluation']['delta_pulse_mhz'],
       'leakage':row['finer_step_evaluation']['coherent_leakage_probability'],
       'evaluations':row['candidate_evaluations']}),flush=True)
    return out


if __name__=='__main__':main()
