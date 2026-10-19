"""Full five-QUQUART (4^5=1024) stochastic physical-channel pilot for SU(3).

The existing 5-qutrit Strang compiler's complete 40/32-CP22 gate list is
applied to a 1024x2 encoded-entanglement purification. Every CP22 can mix
|22> <-> |31> using an explicitly unitary two-level Hamiltonian. After
each symmetric refocusing stage, independent four-level exact bosonic
pure-loss Kraus events and Gaussian phase-diffusion random unitaries are
sampled for all FIVE physical ququarts.

This is a finite-sample quantum-trajectory unraveling of explicitly CPTP
channels, NOT a 1024x1024 density-matrix evolution, calibrated transmon
pulse, full pulse-time Lindblad solver, or implemented syndrome hardware.
Syndrome confusion is a six-outcome instrument followed by replacement on
misreported sectors, not merely multiplying recovery fidelity by readout.

Full circuit CP durations and local gate durations are serial *scenarios*.
The three-level Hamiltonian's action on leaked states is arbitrarily
extended as identity for the electric step; clock-pulse gates, preparation,
and coherent noisy syndrome extraction remain uncompiled.
"""
from __future__ import annotations
import argparse
import json
from collections import Counter
from functools import lru_cache
from itertools import product
from pathlib import Path
from math import comb
import numpy as np

from su3_cubic_casimir_echo import construct
from su3_channel_recovery import syndrome_decoder, source_data
from su3_qutrit_gate_compiler import MASK, CYCLES, ORDER, WEIGHTS
from su3_native_control_compilation import THETA, T
from su3_full_refocusing_cp22 import compiled_site, PATH, N

LOCAL_DIM = 4
SIZE = LOCAL_DIM ** N
COMPUTATIONAL = np.array([
    sum(v * LOCAL_DIM**(N-1-i) for i,v in enumerate(q))
    for q in product(range(3), repeat=N)
], dtype=int)
assert len(COMPUTATIONAL)==243 and len(set(COMPUTATIONAL))==243
IN_CODE = np.zeros(SIZE, bool)
IN_CODE[COMPUTATIONAL] = True
QNUM = np.arange(LOCAL_DIM)


def initial_entangled_map(K):
    """Matrix sqrt(2)|Phi> stored as two physical state columns."""
    mat = np.zeros((SIZE,2), complex)
    mat[COMPUTATIONAL,:] = K
    assert abs(np.vdot(mat,mat).real-2)<1e-12
    return mat


def site_axes(site):
    if site==0: return {'T':0,'R':1}
    if site==N-1: return {'L':N-2,'T':N-1}
    return {'L':site-1,'T':site,'R':site+1}


def apply_one(mat, axis, op):
    """Single-ququart 4x4 gate/Kraus action, preserving the reference index."""
    reshaped = mat.reshape((LOCAL_DIM,)*N+(2,))
    out = np.tensordot(op, reshaped, axes=(1,axis))
    return np.moveaxis(out,0,axis).reshape(SIZE,2)


@lru_cache(maxsize=8192)
def cp_mixing(angle, leakage, detuning=0.):
    """Exact exp(i [[theta,eps],[eps,detune]]) in |22>,|31> pair sector.

    Identical geometric sector convention on either nearest-neighbor edge.
    When eps=detuning=0, this is a perfect extended CP22. Leakage amplitude
    and detuning are radians PER gate, not a measured transmon drive.
    """
    angle,leakage,detuning=float(angle),float(leakage),float(detuning)
    mean=(angle+detuning)/2
    z=(angle-detuning)/2
    norm=np.hypot(z,leakage)
    factor=(np.sin(norm)/norm if norm>1e-15 else 1.)
    part=np.array([[z,leakage],[leakage,-z]],complex)
    return np.exp(1j*mean)*(np.cos(norm)*np.eye(2)+1j*factor*part)


def apply_cp22_mixing(mat, pair, angle, leakage=0., detuning=0.):
    a,b=sorted(pair)
    assert b-a==1
    u=cp_mixing(float(angle),float(leakage),float(detuning))
    z=mat.reshape((LOCAL_DIM,)*N+(2,)).copy()
    pa=[slice(None)]*(N+1)
    pb=[slice(None)]*(N+1)
    pa[a]=pa[b]=2
    pb[a],pb[b]=3,1
    first=z[tuple(pa)].copy()
    second=z[tuple(pb)].copy()
    z[tuple(pa)]=u[0,0]*first+u[0,1]*second
    z[tuple(pb)]=u[1,0]*first+u[1,1]*second
    return z.reshape(SIZE,2)


def apply_compiled_site(mat, site, dt, phase, approximate, leakage):
    axes=site_axes(site)
    schedule=compiled_site(site,float(dt),float(phase),bool(approximate))
    counts=Counter()
    result=mat
    for tag,payload in schedule:
        if tag=='LOCAL':
            label,gate=payload
            enlarged=np.eye(LOCAL_DIM,dtype=complex)
            enlarged[:3,:3]=gate
            result=apply_one(result,axes[label],enlarged)
            counts['local']+=1
        elif tag=='CP22':
            label,angle=payload
            result=apply_cp22_mixing(result,(axes[label],axes['T']),angle,leakage)
            counts['cp22']+=1
        else:raise ValueError(tag)
    return result,counts


def stage_iter(cycles=CYCLES):
    for _ in range(cycles):
        for stage,weight in zip(ORDER,WEIGHTS):
            yield stage,weight


def loss_kraus(dt,t1):
    """Exact 4-level bosonic pure-loss Kraus family, including multiple losses."""
    if dt<0 or t1<=0: raise ValueError("invalid decay parameters")
    eta=np.exp(-dt/t1)
    ops=[]
    for l in range(LOCAL_DIM):
        k=np.zeros((LOCAL_DIM,LOCAL_DIM),complex)
        for n in range(l,LOCAL_DIM):
            k[n-l,n]=np.sqrt(comb(n,l)*(1-eta)**l*eta**(n-l))
        ops.append(k)
    assert np.linalg.norm(sum((k.conj().T@k for k in ops),np.zeros((4,4),complex))-np.eye(4))<5e-11
    return ops


def sample_noise(mat,dt,t1,tphi,rng):
    """Unravel product CPTP loss + diffusion channel, preserving reference.

    Loss: one of four exact Kraus outcomes sampled on the joint (system,
    reference) state. Phase: random unitary with variance 2 dt/Tphi,
    reproducing exp[-dt*(n-m)^2/Tphi] after averaging trajectories.
    """
    if dt==0:
        return mat,0
    ops=loss_kraus(float(dt),float(t1))
    result=mat
    jumps=0
    for site in range(N):
        branches=[apply_one(result,site,k) for k in ops]
        probs=np.array([max(0.,np.vdot(b,b).real/2) for b in branches])
        if abs(sum(probs)-1)>1e-7: raise AssertionError((site,probs))
        probs/=sum(probs)
        choice=int(rng.choice(LOCAL_DIM,p=probs))
        result=branches[choice]/np.sqrt(probs[choice])
        jumps+=int(choice)
        phi=rng.normal(0,np.sqrt(2*dt/tphi))
        local=np.diag(np.exp(-1j*QNUM*phi))
        result=apply_one(result,site,local)
    assert abs(np.vdot(result,result).real-2)<1e-7
    return result,jumps


def readout_fidelity(mat,syndromes,pointer_error=0.,backaction_error=0.):
    """CPTP six-sector projective instrument + classical readout confusion.

    Each of 5 correctable syndrome sectors is distinguished perfectly before
    two independent noisy pointer outputs. A misclassified outcome prepares
    I/2, as does the syndrome complement. Hence this evaluates a physically
    meaningful CPTP *abstract instrument*, but does NOT compile measurement.
    """
    if not 0<=pointer_error<1 or not 0<=backaction_error<1:\n        raise ValueError('readout/backaction probabilities')
    assert abs(np.vdot(mat,mat).real-2)<2e-7
    correct=(1-pointer_error)**2*(1-backaction_error)
    comp=mat[COMPUTATIONAL,:]
    coherent=0.
    supported=0.
    for v in syndromes:
        a=v.conj().T@comp
        coherent+=abs(np.trace(a))**2*correct/4
        supported+=np.vdot(a,a).real*correct
    value=float(coherent+(2-supported)/8)
    leakage=float(np.vdot(mat[~IN_CODE,:],mat[~IN_CODE,:]).real/2)
    assert -1e-8<=value<=1+1e-8
    return value,leakage


def simulate(approximate,leakage=0.01,t1_ns=30000.,tphi_ns=45000.,
             cp_ns=580.,local_ns=30.,shots=4,seed=19,cycles=CYCLES,
             pointer_error=.001,backaction_error=.0,\n             readout_ns=800.,reset_ns=300.,extraction_entanglers=0):
    """Stochastic full-five-ququart digital run; no silent extrapolation.

    40 symmetric stages by default, total 1024-dimensional Hilbert space.
    Each stage has 9 compiled plaquette-site blocks. Duration/decoherence
    accrued serially; physical CP22 and local-qutrit pulse parameters
    remain uncalibrated.
    """
    if shots<1 or cycles<1 or t1_ns<=0 or tphi_ns<=0 or cp_ns<0 or local_ns<0\n       or readout_ns<0 or reset_ns<0 or extraction_entanglers<0:
        raise ValueError("invalid scenario")
    states,K,d,W,H=construct()
    V,_,gram=syndrome_decoder(*source_data())
    assert gram.rank()==5
    phase_diag=np.ones(SIZE, complex)
    qd=np.asarray(d)
    rng=np.random.default_rng(seed)
    metrics=[]
    summary=None
    for trajectory in range(shots):
        y=initial_entangled_map(K)
        total=Counter()
        for stage,weight in stage_iter(cycles):
            dt=T*weight/(CYCLES)
            phase_diag[:]=1
            phase_diag[COMPUTATIONAL]=np.exp(-.5j*dt*qd)
            y=phase_diag[:,None]*y
            counts=Counter()
            for site in PATH:
                phase=THETA*stage if site in MASK else 0.
                small_dt=dt if site==N-1 else dt/2
                y,c=apply_compiled_site(y,site,small_dt,phase,approximate,leakage)
                counts.update(c)
            y=phase_diag[:,None]*y
            duration=counts['cp22']*cp_ns+counts['local']*local_ns
            total.update(counts)
            total['time_ns']+=duration
            y,jumps=sample_noise(y,duration,t1_ns,tphi_ns,rng)
            total['loss_quanta']+=jumps
        # Time spent extracting and measuring a syndrome exposes ALL data\n        # ququarts to idle noise. Entanglers are a COST assumption only.\n        instrument_ns=readout_ns+reset_ns+extraction_entanglers*cp_ns\n        y,instrument_jumps=sample_noise(y,instrument_ns,t1_ns,tphi_ns,rng)\n        total['loss_quanta']+=instrument_jumps\n        total['time_ns']+=instrument_ns\n        f,leak=readout_fidelity(y,V,pointer_error,backaction_error)
        metrics.append({'trajectory':trajectory,'F_e_sample':f,
                        'population_in_fourth_level_end':leak,
                        'loss_quanta_sampled':total['loss_quanta']})
        if summary is None:summary=dict(total)
        assert total['cp22']==((8160 if approximate else 10080)*cycles//CYCLES)
    fs=np.array([x['F_e_sample'] for x in metrics])
    ls=np.array([x['population_in_fourth_level_end'] for x in metrics])
    return {'variant':'approximate_32_cp22' if approximate else 'exact_40_cp22',
            'hilbert_dimension':SIZE,'full_cycles':cycles,
            'trajectory_count':shots,'random_seed':seed,'per_CP22_mixing_radians':leakage,
            'CP22_leakage_pair':'|22> <-> |31> on addressed physical neighbor edge',
            'assumed_T1_ns':t1_ns,'assumed_Tphi_ns':tphi_ns,
            'per_CP22_duration_ns':cp_ns,'per_local_duration_ns':local_ns,
            'readout_pointer_error_each':pointer_error,\n            'syndrome_backaction_replacement_error':backaction_error,\n            'assumed_extraction_entanglers_NOT_IMPLEMENTED':extraction_entanglers,\n            'two_pointer_readout_duration_ns':readout_ns,\n            'ancilla_reset_duration_ns':reset_ns,\n            'minimum_pointer_qutrits':2,'minimum_syndrome_bits':3,
            'serial_resources':summary,
            'sample_mean_recovered_entanglement_fidelity':float(fs.mean()),
            'sample_standard_error':float(fs.std(ddof=1)/np.sqrt(shots)) if shots>1 else None,
            'sample_mean_final_leakage_probability':float(ls.mean()),
            'trajectories':metrics,
            'warning':'Exploratory finite-trajectory CPTP unraveling. Noise is lumped after every stage; no pulse-time Lindblad, gauge-preserving leaked-sector dynamics, calibrated CP22, or noisy extraction circuit.'}


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--shots',type=int,default=4)
    p.add_argument('--leakage',type=float,default=.01)
    p.add_argument('--t1-ns',type=float,default=30000.)
    p.add_argument('--tphi-ns',type=float,default=45000.)
    p.add_argument('--cp-ns',type=float,default=580.)
    p.add_argument('--local-ns',type=float,default=30.)
    p.add_argument('--pointer-error',type=float,default=.001)
    p.add_argument('--backaction-error',type=float,default=.0)\n    p.add_argument('--extraction-entanglers',type=int,default=0)\n    p.add_argument('--readout-ns',type=float,default=800.)\n    p.add_argument('--reset-ns',type=float,default=300.)\n    p.add_argument('--seed',type=int,default=19)
    p.add_argument('--output',default='su3_five_ququart_trajectories_results.json')
    args=p.parse_args()
    rows=[simulate(approximate=x,leakage=args.leakage,t1_ns=args.t1_ns,
                   tphi_ns=args.tphi_ns,cp_ns=args.cp_ns,local_ns=args.local_ns,
                   pointer_error=args.pointer_error,backaction_error=args.backaction_error,\n                   readout_ns=args.readout_ns,reset_ns=args.reset_ns,\n                   extraction_entanglers=args.extraction_entanglers,shots=args.shots,
                   seed=args.seed+i*1000) for i,x in enumerate((False,True))]
    out={'study':'five-ququart encoded reference trajectories; exact vs approx digital CP22',
         'noise_model':'stage-lumped bosonic loss + stochastic Gaussian dephasing; explicit coherent CP22 leakage',
         'rows':rows,'limitations':[
             'No physical gate calibrations, no true pulse-time dissipator.',
             'First-order stroboscopic channel splitting and finite sample size.',
             'Arbitrary extension of electric phase on states containing |3>.',
             'Measurement is an abstract six-sector instrument with independent mislabeling, not native syndrome hardware.',
             'No device-matched uncontrolled/analog physical-noise trajectory comparison.']}
    Path(args.output).write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    for row in rows:
        print('FIVE_QUQUART',json.dumps({k:v for k,v in row.items() if k!='trajectories'},sort_keys=True),flush=True)
    return out

if __name__=='__main__':main()
