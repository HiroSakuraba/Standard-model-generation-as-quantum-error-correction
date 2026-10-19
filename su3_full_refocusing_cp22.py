"""Five-plaquette SU(3) Z3 refocusing using EXPLICIT 40/32-CP22 schedules.

Reproduces the Strang splitting in su3_qutrit_gate_compiler.py, replacing
neighbor-conditioned U(3) blocks with actual gate lists on a five-qutrit
statevector. Device-native CP22 and local qutrit gates remain ASSUMED, not
calibrated. Logical replacement-risk estimates are CPTP toy scenarios, not
T1/T2 or four-level physical-channel propagation.

Run (from repository root):
    OPENBLAS_NUM_THREADS=1 python su3_full_refocusing_cp22.py --rounds 1
"""
from __future__ import annotations
import argparse
import json
from collections import Counter
from functools import lru_cache
from pathlib import Path
import numpy as np
from scipy.sparse.linalg import expm_multiply
from su3_cubic_casimir_echo import construct, entanglement_fidelity
from su3_native_control_compilation import toggled_matrices, evolution, THETA, T
from su3_channel_recovery import syndrome_decoder, source_data
from su3_qutrit_gate_compiler import (
    MASK, CYCLES, ORDER, WEIGHTS, index_groups, local_unitaries,
    strang_segments,
)
from su3_multiplexed_commutator_compiler import (
    compile_mux, emit_local, controlled_u3,
)
from su3_selective_cp_compiler import target_block

N = 5
DIM = 3 ** N
PATH = (0, 1, 2, 3, 4, 3, 2, 1, 0)
SKIPPED = ((1, 2),)


def compile_boundary(site: int, dt: float, phase: float):
    """One-neighbor multiplexor, using 4 CP22s (rather than 40)."""
    if site == 0:
        letter = 'R'
        blocks = [target_block(0, r, dt, phase)[1] for r in range(3)]
    elif site == N-1:
        letter = 'L'
        blocks = [target_block(l, 0, dt, phase)[1] for l in range(3)]
    else:
        raise ValueError('Boundary compilation applies only to sites 0 and 4')
    baseline = blocks[0]
    sequence = []
    emit_local(sequence, 'T', baseline)
    for state in (1, 2):
        controlled_u3(sequence, letter, state, blocks[state] @ baseline.conj().T)
    assert sum(kind == 'CP22' for kind, _ in sequence) == 4
    return sequence


@lru_cache(maxsize=128)
def compiled_site(site: int, dt: float, phase: float, approximate: bool):
    if site in (0, N-1):
        return tuple(compile_boundary(site, dt, phase))
    seq, _, _ = compile_mux(dt, phase, skipped=SKIPPED if approximate else ())
    assert sum(tag == 'CP22' for tag, _ in seq) == (32 if approximate else 40)
    return tuple(seq)


def site_mapping(site):
    if site == 0:
        return {'T': 0, 'R': 1}
    if site == N-1:
        return {'L': N-2, 'T': N-1}
    return {'L': site-1, 'T': site, 'R': site+1}


def apply_gate_schedule(y, site, gates, cp_phase_fraction=0.):
    """Apply actual local U(3)/nearest-neighbor CP22 gates in time order."""
    mapping = site_mapping(site)
    z = np.asarray(y, complex).reshape((3,) * N + (y.shape[-1],)).copy()
    for tag, payload in gates:
        if tag == 'LOCAL':
            label, matrix = payload
            axis = mapping[label]
            z = np.moveaxis(np.tensordot(matrix, z, axes=(1, axis)), 0, axis)
        elif tag == 'CP22':
            controller, angle = payload
            axes = [mapping[controller], mapping['T']]
            assert abs(axes[0]-axes[1]) == 1
            select = [slice(None)] * (N+1)
            for axis in axes:
                select[axis] = 2
            z[tuple(select)] *= np.exp(1j * angle * (1 + cp_phase_fraction))
        else:
            raise ValueError(tag)
    return z.reshape(DIM, y.shape[-1])


def exact_local_action(y, states, site, dt, phase):
    """Independent conditioned 3x3 direct block application (not compiled)."""
    groups = index_groups(states)
    units = local_unitaries(dt, site, phase)
    out = y.copy()
    for ix, left, right in groups[site]:
        out[ix, :] = units[left, right] @ y[ix, :]
    return out


def local_gate_regression(states):
    rng = np.random.default_rng(20261019)
    y = rng.standard_normal((DIM,2)) + 1j*rng.standard_normal((DIM,2))
    y /= np.linalg.norm(y, axis=0)[None,:]
    tests = []
    for site, phase in ((0, 0.), (1, THETA), (2, 2*THETA), (3, 0.), (4, THETA)):
        dt = .013
        ideal = exact_local_action(y, states, site, dt, phase)
        got = apply_gate_schedule(y, site, compiled_site(site, dt, phase, False))
        error = float(np.linalg.norm(got-ideal, ord='fro'))
        assert error < 2e-10, (site, error)
        tests.append({'site':site, 'phase':float(phase), 'fro_error':error,
                      'cp22':sum(tag=='CP22' for tag,_ in compiled_site(site,dt,phase,False))})
    return tests


def propagate(K, d, rounds, approximate=False, cp_phase_fraction=0.):
    """40 Z3 stages, 9 site operations per round/stage with explicit CP22s."""
    assert rounds >= 1
    out = np.array(K, dtype=complex, copy=True)
    counter = Counter()
    per_site = Counter()
    for _ in range(CYCLES):
        for toggle, weight in zip(ORDER, WEIGHTS):
            dt = T * weight / (CYCLES * rounds)
            electric = np.exp(-.5j * dt * d)
            for _ in range(rounds):
                out *= electric[:,None]
                for site in PATH:
                    phase = THETA * toggle if site in MASK else 0.
                    step = dt if site == N-1 else dt/2
                    gates = compiled_site(site, step, phase, approximate)
                    out = apply_gate_schedule(out, site, gates, cp_phase_fraction)
                    cnt = Counter(tag for tag,_ in gates)
                    counter['cp22'] += cnt['CP22']
                    counter['local_qutrit'] += cnt['LOCAL']
                    counter['magnetic_site_blocks'] += 1
                    per_site[site] += cnt['CP22']
                out *= electric[:,None]
                counter['electric_half_phases'] += 2
            counter['stages'] += 1
    assert abs(np.linalg.norm(out[:,0]) - 1) < 1e-9
    assert np.linalg.norm(out.conj().T@out - np.eye(2)) < 2e-9
    assert counter['stages'] == CYCLES * len(ORDER)
    assert counter['magnetic_site_blocks'] == CYCLES * len(ORDER) * len(PATH) * rounds
    assert counter['cp22'] == (9600 if not approximate else 7680)*rounds + 480*rounds
    return out, {**dict(counter), 'cp22_by_plaquette':{str(k):per_site[k] for k in range(N)}}


def threshold_error(fidelity, baseline, n_cp, floor_survival=1.):
    """Break-even for independently flagged gate faults and I/2 replacement."""
    if fidelity <= baseline:
        return None
    required = (baseline - .25) / (fidelity - .25)
    if floor_survival <= required:
        return None
    return float(-np.expm1(np.log(required/floor_survival)/n_cp))


def toy_risk(fidelity, baseline, cost, p_cp, p_local,
             entangler_ns=580., local_ns=30.,
             readout_us=0., p_readout=.0, t1_us=None, dephase_jump_us=None):
    """Logical CPTP replacement channel, NOT physical qutrit noise simulation.

    Independent flagged event replaces the output qubit by I/2 (F_e=1/4).
    Exponential hazards are conservative illustrative replacement-event rates.
    """
    assert 0 <= p_cp < 1 and 0 <= p_local < 1 and 0 <= p_readout < 1
    gate_us = (cost['cp22']*entangler_ns + cost['local_qutrit']*local_ns)/1000.
    elapsed_us = gate_us + readout_us
    log_survival = (cost['cp22']*np.log1p(-p_cp)
                    + cost['local_qutrit']*np.log1p(-p_local)
                    + 2*np.log1p(-p_readout))
    if t1_us is not None:
        assert t1_us > 0
        log_survival -= N*elapsed_us/t1_us
    if dephase_jump_us is not None:
        assert dephase_jump_us > 0
        log_survival -= N*elapsed_us/dephase_jump_us
    survival = float(np.exp(log_survival))
    recovered = float(survival * fidelity + (1-survival)*.25)
    return {'p_cp22_flagged':p_cp, 'p_local_flagged':p_local,
            'p_each_pointer_readout_flagged':p_readout,
            'readout_us_assumed':readout_us, 't1_replacement_hazard_us':t1_us,
            'dephasing_replacement_hazard_us':dephase_jump_us,
            'gate_time_us_serial':gate_us, 'total_time_us_serial_with_readout':elapsed_us,
            'no_flag_survival':survival, 'logical_replacement_F_e':recovered,
            'beats_uncontrolled_reference':recovered > baseline}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--rounds', type=int, default=1)
    parser.add_argument('--output', default='su3_full_refocusing_cp22_results.json')
    parser.add_argument('--skip-regression', action='store_true')
    args=parser.parse_args()
    states,K,d,W,H=construct()
    decoder,_,gram=syndrome_decoder(*source_data())
    assert gram.rank() == 5
    tests=[] if args.skip_regression else local_gate_regression(states)
    ideal = evolution(K, toggled_matrices(d,W,MASK), T,CYCLES)
    ideal_f = entanglement_fidelity([ideal],K,decoder)
    uncontrolled=expm_multiply(-1j*T*H,K,traceA=0)
    baseline=entanglement_fidelity([uncontrolled],K,decoder)
    rows=[]
    target_product, abstract_counts = strang_segments(K,states,d,index_groups(states),rounds=args.rounds)
    for approx in (False, True):
        y, counts=propagate(K,d,args.rounds,approximate=approx)
        f=entanglement_fidelity([y],K,decoder)
        refdiff=float(np.linalg.norm(y-target_product,ord='fro'))
        if not approx:
            assert refdiff < 3e-8, refdiff
        ncp=counts['cp22']
        fid_baseline=threshold_error(f,baseline,ncp)
        scenario=[toy_risk(f,baseline,counts,p_cp,p_local,readout_us=2.,p_readout=.001)
                  for p_cp,p_local in ((0.,0.),(1e-5,1e-5),(1e-4,1e-4),(1e-3,1e-3))]
        scenario.append(toy_risk(f,baseline,counts,0.,0.,readout_us=2.,
                                 p_readout=.001,t1_us=1e5,dephase_jump_us=1e5))
        row={'compiler':'approximate_32_CP22' if approx else 'exact_40_CP22',
             'noiseless_recovered_F_e':f,
             'raw_F_e':entanglement_fidelity([y],K),
             'state_error_vs_ideal_stage_evolution_fro':float(np.linalg.norm(y-ideal,ord='fro')),
             'state_error_vs_direct_Strang_fro':refdiff,
             'cp22_flagged_break_even_probability_per_gate_other_ops_perfect':fid_baseline,
             **counts, 'toy_replacement_scenarios':scenario}
        rows.append(row)
        print('COMPILER',json.dumps({k:v for k,v in row.items() if k!='toy_replacement_scenarios'},sort_keys=True),flush=True)
        for entry in scenario:print('TOY_NOISE',json.dumps(entry,sort_keys=True),flush=True)
        if approx:
            row['worst_case_encoded_input_vector_difference_vs_exact'] = float(
                np.linalg.svd(y-target_exact, compute_uv=False)[0])
        else:
            target_exact=y
    data={'protocol':'five-face SU3 refocusing; 8 cycles, 5 symmetric stages per cycle',
          'source':'Chen-Mueller-Yao 2026 Appendix B Table VI, via existing repository',
          'n_qutrits':N,'logical_dimension':2,'rounds_per_stage':args.rounds,
          'mask':list(MASK),'total_dimensionless_time':T,
          'uncontrolled_recovered_F_e':baseline,
          'ideal_exact_stages_recovered_F_e':ideal_f,
          'abstract_logical_Strang_count':abstract_counts,
          'local_243D_gate_regression':tests,'rows':rows,
          'minimum_syndrome_pointer_qutrits':2, 'minimum_classical_readout_bits':3,
          'uncompiled_costs':['gauge-preserving encoding','clock controls','native CP22 pulse synthesis',
                              'syndrome extraction','ancilla reset','logical recovery'],
          'noise_limitations':('CPTP flagged logical replacement channel only; not physical qutrit T1/T2, '
                              'not leakage to a fourth level, not calibrated gate or readout error. '
                              'Two qutrit readouts and serial durations are assumptions. '
                              'Real magnetic coupling/refocusing pulse engineering is unsolved.'),
          'warning':('No claim of a physical SU3 quantum memory or gluonic error correction. '
                     'Gate angle-selectivity and local rotations are hypothetical.')}
    Path(args.output).write_text(json.dumps(data,indent=2,sort_keys=True)+'\n')
    print('FULL FIVE-PLAQUETTE CP22 CIRCUIT VERIFIED',json.dumps({
        'base_recovered':baseline,'ideal_evolution_recovered':ideal_f,
        'exact_gate_count':rows[0]['cp22'],'approx_gate_count':rows[1]['cp22']}),flush=True)
    return data


if __name__ == '__main__':
    main()
