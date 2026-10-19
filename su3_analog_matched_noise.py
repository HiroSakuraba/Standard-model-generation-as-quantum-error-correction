"""Matched-time analog SU(3) refocusing vs free evolution with noisy syndrome pointers.

The exact 243D truncated published five-face H is propagated throughout 8
cycles of finite-duration cubic-Casimir pulses. Uncontrolled storage is
propagated for the SAME 0.6+32*tau dimensionless duration with identical
segment boundaries and sampled *gauge-register* three-level loss and phase
diffusion. Physical clock controls and native SU(3) Wilson interactions
remain hypothetical; this is NOT a superconducting transmon simulation.

Six mutually orthogonal syndrome sectors are coherently tagged by TWO
qutrit-pointer basis states through the mathematically exact nonlocal
isometry sum_r P_r tensor |pointer_r>. Pointer amplitude damping and
assignment error give a stochastic readout confusion matrix. Only correctly
reported syndrome outcomes are decoded; all other outcomes are replaced by
I_2/2. The nonlocal sector-controlled operation and decoder are NOT
compiled to native local gates.

State 243x2 purifies the encoded maximally mixed qubit. Each conditional
Kraus sample is properly weighted and normalized. Number dephasing uses
Gaussian stochastic phases, and 3-level pure-loss CPTP Kraus operators
are sampled after EACH native free/pulse segment (not pulse-time Lindblad).
"""
from __future__ import annotations
import argparse
import json
from dataclasses import dataclass, asdict
from math import comb
from pathlib import Path
import numpy as np
from scipy.linalg import eigh
from scipy.sparse.linalg import expm_multiply

from su3_cubic_casimir_echo import construct, entanglement_fidelity
from su3_channel_recovery import source_data, syndrome_decoder
from su3_native_control_compilation import THETA,T,WEIGHTS
from su3_qutrit_gate_compiler import MASK, CYCLES

N=5
LOCAL=3
DIM=LOCAL**N
POINTERS=((0,0),(0,1),(0,2),(1,0),(1,1),(1,2))
POINTER_INDEX=np.array([3*a+b for a,b in POINTERS],dtype=int)
LABELS=np.arange(DIM).reshape((LOCAL,)*N)
NUMBER=np.arange(LOCAL)


@dataclass(frozen=True)
class Scenario:
    """Hypothetical physical register parameters; not hardware-calibrated."""
    scale_mhz:float=.01
    tau:float=.002
    t1_ns:float=30000.
    tphi_ns:float=45000.
    pointer_t1_ns:float=30000.
    pointer_read_ns:float=800.
    pointer_assignment_error:float=.001
    pointer_reset_ns:float=300.
    extraction_ns:float=0.
    clock_angle_fractional_bias:float=0.
    shots:int=48
    seed:int=20261019


def time_ns(dimless,scale_mhz):
    if dimless<0 or scale_mhz<=0:raise ValueError('invalid conversion')
    return dimless/(2*np.pi*.001*scale_mhz)


def number_clock(states,subset=MASK):
    return np.array([sum(0 if q[i]==0 else (1 if q[i]==1 else -1)
                       for i in subset) for q in states],float)


def exponentiate_hermitian(H,t):
    eig,vec=eigh(H)
    return (vec*np.exp(-1j*eig*t))@vec.conj().T


def propagation_schedule(states,H,scenario):
    """72 time-matched segments. Native H ON throughout analog sequence.

    Analog free segments include the physical H. Four pulse windows per
    cycle use H +/- THETA*(1+bias)/tau times a *diagonal cubic clock*.
    Uncontrolled branch evolves under the SAME H for pulse-window time.
    """
    if scenario.tau<=0:raise ValueError('positive finite pulse width needed')
    dense=H.toarray()
    clock=np.diag(number_clock(states))
    options={}
    schedule=[]
    for cycle in range(CYCLES):
        for k,weight in enumerate(WEIGHTS):
            dt=T*weight/CYCLES
            tag=('free',round(dt,15))
            if tag not in options:
                options[tag]=exponentiate_hermitian(dense,dt)
            schedule.append((options[tag],options[tag],time_ns(dt,scenario.scale_mhz)))
            if k<4:
                sign=1 if k<2 else -1
                pulse=('pulse',sign)
                if pulse not in options:
                    control=(THETA/scenario.tau)*(1+scenario.clock_angle_fractional_bias)
                    options[pulse]=exponentiate_hermitian(dense+sign*control*clock,scenario.tau)
                idlekey=('idle',scenario.tau)
                if idlekey not in options:
                    options[idlekey]=exponentiate_hermitian(dense,scenario.tau)
                schedule.append((options[pulse],options[idlekey],
                                 time_ns(scenario.tau,scenario.scale_mhz)))
    assert len(schedule)==CYCLES*9
    assert abs(sum(x[2] for x in schedule)
               -time_ns(T+4*CYCLES*scenario.tau,scenario.scale_mhz))<1e-7
    return schedule


def loss_kraus(dt,t1):
    if dt<0 or t1<=0: raise ValueError('invalid damping parameters')
    eta=float(np.exp(-dt/t1))
    out=[]
    for jump in range(LOCAL):
        mat=np.zeros((LOCAL,LOCAL),complex)
        for n in range(jump,LOCAL):
            mat[n-jump,n]=np.sqrt(comb(n,jump))*(1-eta)**(jump/2)*eta**((n-jump)/2)
        out.append(mat)
    assert np.linalg.norm(sum((k.conj().T@k for k in out),np.zeros((3,3),complex))-np.eye(3))<2e-11
    return out


def apply_site(mat,site,gate):
    data=mat.reshape((LOCAL,)*N+(2,))
    result=np.tensordot(gate,data,axes=(1,site))
    return np.moveaxis(result,0,site).reshape(DIM,2)


def noise_segment(y,duration_ns,t1,tphi,rng):
    """Single time segment: CPTP product pure loss plus Gaussian phase diffusion.

    The sampled two-column (system,reference) state remains normalized in
    squared Frobenius norm to 2. This is not a conditional pure logical state.
    """
    if duration_ns==0:return y,0
    if t1<=0 or tphi<=0:raise ValueError('invalid coherence times')
    ops=loss_kraus(duration_ns,t1)
    jumps=0
    for site in range(N):
        branches=[apply_site(y,site,k) for k in ops]
        weights=np.array([max(0.,np.vdot(m,m).real/2) for m in branches])
        assert abs(weights.sum()-1)<1e-8,(site,weights)
        weights/=weights.sum()
        choice=int(rng.choice(LOCAL,p=weights))
        y=branches[choice]/np.sqrt(weights[choice])
        jumps+=choice
        phase=rng.normal(0,np.sqrt(2*duration_ns/tphi))
        y=apply_site(y,site,np.diag(np.exp(1j*NUMBER*phase)))
    assert abs(np.vdot(y,y).real-2)<1e-7
    return y,jumps


def syndrome_projectors(syndromes):
    """Nonlocal six-sector PVM: 5 rank-2 sectors and rank-233 complement."""
    ps=[v@v.conj().T for v in syndromes]
    complement=np.eye(DIM)-sum(ps)
    ps.append(complement)
    for i,p in enumerate(ps):
        assert np.linalg.norm(p@p-p)<2e-10
        assert np.linalg.norm(p-p.conj().T)<1e-10
        if i<5:assert abs(np.trace(p)-2)<2e-10
    assert np.linalg.norm(sum(ps)-np.eye(DIM))<2e-10
    assert abs(np.trace(ps[-1]).real-(DIM-10))<2e-10
    return ps


def tag_syndrome_pointers(y,projectors):
    """Exact PVM Stinespring isometry into two pointer qutrits, 9 labels.

    Assumes a *nonlocal controlled interaction*. The isometry is
    mathematically explicit but not a compiled hardware gate network.
    """
    if len(projectors)!=6 or y.shape[0]!=DIM:raise ValueError('shape')
    tagged=np.zeros((9,DIM,y.shape[1]),complex)
    for r,p in enumerate(projectors):
        tagged[POINTER_INDEX[r]]=p@y
    assert abs(np.vdot(tagged,tagged).real-np.vdot(y,y).real)<2e-9
    assert np.linalg.norm(tagged.sum(axis=0)-y)<2e-9
    return tagged


def pointer_single_transition(read_ns,t1_ns,assignment):
    """P(output basis digit | tagged input digit).

    Pointer qutrit decays with exact 3-level pure-loss channel; we then
    make a symmetric readout assignment error. Phase noise commutes with
    population readout in the tagged computational basis and is irrelevant.
    """
    if not 0<=assignment<1:raise ValueError('pointer assignment probability')
    assign=np.full((LOCAL,LOCAL),assignment/(LOCAL-1))
    np.fill_diagonal(assign,1-assignment)
    transition=np.zeros((LOCAL,LOCAL))
    for tagged in range(LOCAL):
        for k in loss_kraus(read_ns,t1_ns):
            prob=np.abs(k[:,tagged])**2
            transition[:,tagged]+=assign@prob
    assert np.max(abs(transition.sum(axis=0)-1))<1e-12
    assert np.min(transition)>=0
    return transition


def pointer_confusion(scenario):
    """Full 9 observed two-pointer outcomes by 6 tagged syndromes."""
    t=pointer_single_transition(scenario.pointer_read_ns,
                                scenario.pointer_t1_ns,
                                scenario.pointer_assignment_error)
    out=np.zeros((9,6))
    for idx,(a,b) in enumerate(POINTERS):
        out[:,idx]=np.kron(t[:,a],t[:,b])
    assert np.max(abs(out.sum(axis=0)-1))<1e-12
    return out


def syndrome_recovered_fidelity(y,syndromes,confusion):
    """Exact CPTP instrument and recovery fidelity for each encoded trajectory.

    The outcome r is tagged in the two pointers through P_r; a corrupted
    pointer label is never interpreted as the original r. Correct labels
    apply V_r^dag, while bad labels and complement replace with I/2.
    """
    assert y.shape==(DIM,2)
    assert abs(np.vdot(y,y).real-2)<2e-7
    good=0.;supported=0.
    for r,v in enumerate(syndromes):
        a=v.conj().T@y
        p=float(confusion[POINTER_INDEX[r],r])
        good+=p*abs(np.trace(a))**2/4
        supported+=p*np.vdot(a,a).real
    fidelity=float(good+(2-supported)/8)
    assert -1e-8<=fidelity<=1+1e-8,fidelity
    return fidelity, float(supported/2)


def run_one(y,segments,arm,scenario,syndromes,confusion,rng):
    jumps=0
    for active,idle,dt_ns in segments:
        y=(active if arm=='analog' else idle)@y
        y,loss=noise_segment(y,dt_ns,scenario.t1_ns,scenario.tphi_ns,rng)
        jumps+=loss
    # Both paths experience identical syndrome extraction, read/reset idle
    # time. Projective measurement itself is an ideal, nonlocal 6-sector PVM.
    extra=scenario.extraction_ns+scenario.pointer_read_ns+scenario.pointer_reset_ns
    y,loss=noise_segment(y,extra,scenario.t1_ns,scenario.tphi_ns,rng)
    jumps+=loss
    f,weight=syndrome_recovered_fidelity(y,syndromes,confusion)
    return f,weight,jumps


def study(scenario=Scenario()):
    if scenario.shots<1 or scenario.extraction_ns<0 or scenario.pointer_reset_ns<0:
        raise ValueError('invalid sampling parameters')
    states,K,d,W,H=construct()
    syndrome,_,gram=syndrome_decoder(*source_data())
    assert gram.rank()==5
    ps=syndrome_projectors(syndrome)
    tagged=tag_syndrome_pointers(K,ps)
    assert abs(np.vdot(tagged,tagged).real-2)<1e-9
    confusion=pointer_confusion(scenario)
    segments=propagation_schedule(states,H,scenario)
    # The explicit isometry K -> sum_r P_r K ⊗ |pointer_r> is
    # norm-preserving by the rank-six projector resolution of identity.
    rng=np.random.default_rng(scenario.seed)
    shots={'analog':[],'uncontrolled':[]}
    weights={'analog':[],'uncontrolled':[]}
    jumps={'analog':[],'uncontrolled':[]}
    for shot in range(scenario.shots):
        for arm in ('analog','uncontrolled'):
            f,w,loss=run_one(K.copy(),segments,arm,scenario,syndrome,confusion,rng)
            shots[arm].append(f);weights[arm].append(w);jumps[arm].append(loss)
    # Independently crosscheck the noiseless references against published
    # lab_cubic_pulses for analog, and expm(-i H total) for control.
    from su3_native_control_compilation import lab_cubic_pulses
    analog_ref=lab_cubic_pulses(states,H,K,MASK,CYCLES,scenario.tau,
                               phase_error={i:scenario.clock_angle_fractional_bias for i in MASK})
    a=K.copy();u=K.copy()
    for active,idle,_ in segments:
        a=active@a;u=idle@u
    assert np.linalg.norm(a-analog_ref,'fro')<2e-9
    no_control=expm_multiply(-1j*(T+4*CYCLES*scenario.tau)*H,K,traceA=0)
    assert np.linalg.norm(u-no_control,'fro')<2e-9
    ideal={'analog':entanglement_fidelity([a],K,syndrome),
           'uncontrolled':entanglement_fidelity([u],K,syndrome)}
    rows=[]
    for arm in ('analog','uncontrolled'):
        samples=np.asarray(shots[arm])
        rows.append({'arm':arm,
                     'noiseless_ideal_decoder_F_e':ideal[arm],
                     'noisy_pointer_channel_samples':shots[arm],
                     'noisy_F_e_mean':float(samples.mean()),
                     'noisy_F_e_standard_error':float(samples.std(ddof=1)/np.sqrt(len(samples))) if len(samples)>1 else None,
                     'mean_correct_syndrome_weight':float(np.mean(weights[arm])),
                     'mean_sampled_local_decay_quanta':float(np.mean(jumps[arm]))})
    meta={'scenario':asdict(scenario),
          'duration_ns_each_arm':sum(x[2] for x in segments)+scenario.extraction_ns+scenario.pointer_read_ns+scenario.pointer_reset_ns,
          'control_pulse_windows':4*CYCLES,'addressed_clock_sites':len(MASK),
          'pointer_labels':[list(z) for z in POINTERS],
          'pointer_self_report_probs':[float(confusion[POINTER_INDEX[r],r]) for r in range(6)],
          'syndrome_sector_ranks':[2]*5+[DIM-10],
          'pointer_register_qutrits':2,'classical_readout_bits_min':3,
          'random_seed':scenario.seed,'noise_events_per_trial':(len(segments)+1)*N,
          'reference_difference_fro_max':float(max(np.linalg.norm(a-analog_ref,'fro'),np.linalg.norm(u-no_control,'fro'))),
          'caveat':('Uncalibrated native SU3 analog Hamiltonian and selective cubic-Casimir clocks; '
                    'gauge-fixed qutrit register amplitude loss/phase diffusion is illustrative, '
                    'not physical gauge-field noise; local-noise Kraus after each physical segment, '
                    'not pulse-time Lindblad; six-outcome syndrome-pointer tag is nonlocal '
                    'and not implemented by local gates. Equal elapsed time and readout budget only '
                    'within this analog-versus-no-control model, not versus digital hardware.') }
    return {'meta':meta,'results':rows}


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--shots',type=int,default=48)
    p.add_argument('--scale-mhz',type=float,default=.01)
    p.add_argument('--t1-ns',type=float,default=30000.)
    p.add_argument('--tphi-ns',type=float,default=45000.)
    p.add_argument('--pointer-t1-ns',type=float,default=30000.)
    p.add_argument('--assignment',type=float,default=.001)
    p.add_argument('--bias',type=float,default=0.)
    p.add_argument('--extraction-ns',type=float,default=0.)
    p.add_argument('--seed',type=int,default=20261019)
    p.add_argument('--output',default='su3_analog_matched_noise_results.json')
    args=p.parse_args()
    settings=Scenario(shots=args.shots,scale_mhz=args.scale_mhz,
         t1_ns=args.t1_ns,tphi_ns=args.tphi_ns,pointer_t1_ns=args.pointer_t1_ns,
         pointer_assignment_error=args.assignment,
         clock_angle_fractional_bias=args.bias,
         extraction_ns=args.extraction_ns,seed=args.seed)
    output=study(settings)
    Path(args.output).write_text(json.dumps(output,indent=2,sort_keys=True)+'\n')
    for row in output['results']:
        print('MATCHED_NOISE',json.dumps({k:v for k,v in row.items()
             if k!='noisy_pointer_channel_samples'},sort_keys=True),flush=True)
    print('ANALOG_MATCHED_NOISE_COMPLETE',json.dumps({
        'shots':args.shots,'duration_ns_each_arm':output['meta']['duration_ns_each_arm'],
        'clock_pulse_windows':output['meta']['control_pulse_windows']}),flush=True)


if __name__=='__main__':main()
