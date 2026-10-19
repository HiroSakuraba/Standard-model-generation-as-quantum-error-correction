"""Exact nine-context SU(3) magnetic multiplexor compiled into nearest-neighbor CP22.

Uses factorized two-neighbor multiplexor and U(3) commutator identity.
Assumes arbitrary-angle |22> controlled-phase on L-T or T-R and ideal local
qutrit U(3) (decomposed/verified to adjacent-level SU(2) and diagonal phases).
The gate set is an ASSUMPTION, not a demonstrated transmon calibration.
"""
from __future__ import annotations
import argparse
import json
from collections import Counter
from pathlib import Path
import numpy as np
from scipy.linalg import expm, schur
from su3_selective_cp_compiler import target_block, embedded_target, site_gate, fixed_cp22_gate, adjacent_swaps_to_two, simulate_fixed_22, compile_full_interior, PUBLISHED

EYE3=np.eye(3,dtype=complex)
EYE27=np.eye(27,dtype=complex)


def unitary_error(u):
    return float(np.linalg.norm(u.conj().T@u-np.eye(u.shape[0]),ord=2))


def local_qr(u):
    """Decompose an arbitrary U(3) into <=3 adjacent SU2 mixers + one diagonal.
    The diagonal is a qutrit virtual-Z frame (three phases, one unobservable global).
    """
    q=u.copy();elim=[]
    for a,b,col in ((1,2,0),(0,1,0),(1,2,1)):
        x,y=q[a,col],q[b,col]
        norm=np.hypot(abs(x),abs(y))
        if norm<1e-13:continue
        g=np.array([[np.conj(x),np.conj(y)],[-y,x]],complex)/norm
        big=EYE3.copy();big[np.ix_((a,b),(a,b))]=g
        q=big@q;elim.append(big)
    assert np.linalg.norm(q-np.diag(np.diag(q)),ord='fro')<3e-11
    diag=np.diag(np.diag(q))
    assert unitary_error(diag)<3e-11
    schedule=[diag]+[m.conj().T for m in reversed(elim)]
    z=EYE3.copy()
    for gate in schedule:z=gate@z
    assert np.linalg.norm(z-u,ord=2)<3e-11
    return schedule,len(elim),1


def emit_local(out,site,u,split=True):
    if split:
        schedule,_,_=local_qr(u)
        for mat in schedule:out.append(('LOCAL',(site,mat)))
    else:out.append(('LOCAL',(site,u)))


def emit_cp22(out,controller,ctrl_value,target_level,angle):
    """Route a computational control/target level to |22> via adjacent swaps."""
    for g in adjacent_swaps_to_two(ctrl_value):out.append(('LOCAL',(controller,g)))
    for g in adjacent_swaps_to_two(target_level):out.append(('LOCAL',('T',g)))
    out.append(('CP22',(controller,float(angle))))
    for g in reversed(adjacent_swaps_to_two(target_level)):out.append(('LOCAL',('T',g)))
    for g in reversed(adjacent_swaps_to_two(ctrl_value)):out.append(('LOCAL',(controller,g)))


def controlled_u3(out,controller,value,u,split=True):
    """Exact C_value(U) with 2 selective CP22 and local qutrit phases.

    Schur U=V diag(lambda) Vdag; use global eigenphase lambda_0 as
    a **controller-only** phase gate, then two relative target eigenphases.
    This does not introduce an uncounted entangling controlled global phase.
    """
    assert unitary_error(u)<3e-10
    t,v=schur(u,output='complex')
    assert np.linalg.norm(t-np.diag(np.diag(t)))<1e-9
    lam=np.diag(t)
    assert np.max(np.abs(abs(lam)-1))<3e-10
    emit_local(out,'T',v.conj().T,split)
    control_diag=np.eye(3,dtype=complex)
    control_diag[value,value]=lam[0]
    emit_local(out,controller,control_diag,split=False)
    for k in (1,2):
        emit_cp22(out,controller,value,k,np.angle(lam[k]/lam[0]))
    emit_local(out,'T',v,split)


def su3_commutator(u):
    """Construct A,B with A B A^dag B^dag = U, det(U)=1.

    In U's eigenbasis, A cyclically permutes the three eigenvectors;
    B=diag(1,1/lambda_1,1/(lambda_1 lambda_2)).
    """
    t,v=schur(u,output='complex')
    lam=np.diag(t)
    assert np.linalg.norm(t-np.diag(lam))<1e-9
    assert abs(np.prod(lam)-1)<1e-9, np.prod(lam)
    shift=np.roll(EYE3,1,axis=0)
    b=np.diag([1.,1/lam[1],1/(lam[1]*lam[2])])
    a=v@shift@v.conj().T
    b=v@b@v.conj().T
    error=np.linalg.norm(a@b@a.conj().T@b.conj().T-u,ord=2)
    assert error<3e-9,error
    return a,b


def double_controlled_su3(out,left,right,u,split=True):
    a,b=su3_commutator(u)
    # Chronological sequence gives [C_L(A), C_R(B)] on intersection sector.
    for ctrl,val,op in [('R',right,b.conj().T),('L',left,a.conj().T),
                        ('R',right,b),('L',left,a)]:
        controlled_u3(out,ctrl,val,op,split)


def wanted(dt,phase):
    blocks={(l,r):target_block(l,r,dt,phase)[1] for l in range(3) for r in range(3)}
    full=EYE27.copy()
    for (l,r),u in blocks.items():full+=embedded_target(u,l,r)-EYE27
    return blocks,full


def compile_mux(dt=.19,phase=2*np.pi/3,split=True,skipped=()):
    blocks,reference=wanted(dt,phase)
    base=blocks[0,0]
    left={l:blocks[l,0]@base.conj().T for l in (1,2)}
    right={r:blocks[0,r]@base.conj().T for r in (1,2)}
    out=[]
    emit_local(out,'T',base,split)
    for l in (1,2):controlled_u3(out,'L',l,left[l],split)
    for r in (1,2):controlled_u3(out,'R',r,right[r],split)
    residuals={}
    for l in (1,2):
        for r in (1,2):
            predicted=right[r]@left[l]@base
            delta=blocks[l,r]@predicted.conj().T
            residuals[l,r]=float(np.linalg.norm(delta-EYE3,ord=2))
            if (l,r) not in skipped:double_controlled_su3(out,l,r,delta,split)
    count=Counter(tag for tag,_ in out)
    assert count['CP22']==40-8*len(skipped)
    # 4 single-control U3 plus 4 doubly-controlled residual U3;
    # each double uses 4 single controls. => 20 controlled SU3 x 2 CP22.
    return out,reference,residuals


def calculate_unitary(schedule):
    result=EYE27.copy()
    for tag,arg in schedule:
        if tag=='LOCAL':
            site,g=arg
            result=site_gate(g,site)@result
        else:result=fixed_cp22_gate(*arg)@result
    return result


def verify(schedule,target,tolerance=2e-9):
    result=calculate_unitary(schedule)
    e=np.linalg.norm(result-target,ord=2)
    assert e<tolerance,e
    return float(e)


def local_counts(schedule):
    cnt=Counter(tag for tag,_ in schedule)
    su2=0;diagonal=0;swaps=0
    for tag,arg in schedule:
        if tag!='LOCAL':continue
        _,m=arg
        if np.linalg.norm(m-np.diag(np.diag(m)))<1e-11:diagonal+=1
        else:
            dif=m-np.eye(3)
            active=tuple(i for i in range(3) if np.linalg.norm(dif[i,:])>1e-10 or np.linalg.norm(dif[:,i])>1e-10)
            assert active in ((0,1),(1,2)),active
            su2+=1
            if np.linalg.norm(m@m-np.eye(3))<1e-11 and np.max(abs(m-np.abs(m)))<1e-10:swaps+=1
    return {'cp22':cnt['CP22'],'single_qutrit_ops':cnt['LOCAL'],
            'adjacent_su2_rotations_or_swaps':su2,'diagonal_qutrit_phase_frames':diagonal,
            'adjacent_level_swaps':swaps}


def crosscheck_reference():
    """Second, original repository's independent published Table VI transcription."""
    from su3_three_plaquette_chain import PUBLISHED as independently_transcribed
    for pair,mat in PUBLISHED.items():
        assert np.linalg.norm(np.asarray(mat,float)-np.asarray(independently_transcribed[pair],float))<1e-13
    return 9


def leakage_audit(schedule,target):
    """Reuse PR36's same deliberately phenomenological four-level model.

    It does NOT include T1/T2 decoherence during the schedule, calibrated
    waveform bandwidth, leakage on single-qutrit operations or readout.
    """
    old=[('LOCAL_AT',v) if t=='LOCAL' else (t,v) for t,v in schedule]
    out=[]
    for leak,xtalk,bias in ((0.,0.,0.),(.01,.002,0.),(.03,.006,.02)):
        row=simulate_fixed_22(old,None,leak=leak,crosstalk=xtalk,
                              angle_fraction=bias,target_27=target)
        row.update({'leak_amplitude':leak,'spectator_crosstalk_phase':xtalk,
                    'fractional_entangler_phase_error':bias})
        assert row['counts']['fixed_CP22_entanglers']==40
        if leak==0:assert row['unconditional_entanglement_fidelity']>1-1e-11
        out.append(row)
    return out


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',default='su3_multiplexed_commutator_results.json');args=p.parse_args()
    assert crosscheck_reference()==9
    rows=[]
    baseline=compile_full_interior(run_fault_test=False)
    assert baseline['full_CP22_count']==360 and baseline['full_single_qutrit_ops']==1728
    for dt,phase in ((.19,2*np.pi/3),(.08,0.),(.3,.23),(.047,-1.1)):
        seq,ideal,residual=compile_mux(dt,phase)
        error=verify(seq,ideal)
        counts=local_counts(seq)
        row={'dt':dt,'magnetic_phase':float(phase),'operator_error':error,
             'residual_norms':{f'{l},{r}':v for (l,r),v in residual.items()},
             **counts}
        assert counts['cp22']==40
        if not rows:
            row['same_four_level_model']=leakage_audit(seq,ideal)
        cheap,_,_=compile_mux(dt,phase,skipped=((1,2),))
        approx=calculate_unitary(cheap)
        row['approx_32_cp22']={
          'operator_norm_error':float(np.linalg.norm(approx-ideal,ord=2)),
          'process_entanglement_fidelity':float(abs(np.trace(ideal.conj().T@approx))**2/27**2)}
        rows.append(row)
        print('MUX',json.dumps(row,sort_keys=True))
    # Transparent approximate 10x-or-better budget: omit one or more
    # double-control residual blocks, report the actual unitary error.
    approximations=[]
    for skipped in (((1,2),),((1,2),(2,1)),((1,1),(1,2),(2,1),(2,2))):
        sequence,target,_=compile_mux(skipped=skipped)
        comp=calculate_unitary(sequence)
        op_err=float(np.linalg.norm(comp-target,ord=2))
        process=float(abs(np.trace(target.conj().T@comp))**2/27**2)
        ncp=sum(tag=='CP22' for tag,arg in sequence)
        assert ncp==40-8*len(skipped)
        assert 0<op_err<.3 and process<1
        row={'skipped_double_control_contexts':[list(v) for v in skipped],
             'cp22':ncp,'two_qutrit_reduction_factor':360/ncp,
             'operator_norm_error':op_err,'process_entanglement_fidelity':process,
             'single_qutrit_ops':local_counts(sequence)['single_qutrit_ops']}
        approximations.append(row)
        print('APPROXIMATE_OPTION',json.dumps(row,sort_keys=True))
    data={'source':'Chen Mueller Yao PRD 114 074502, Appendix B Table VI; previous 360-CP22 baseline',
          'assumed_connectivity':'line L-T-R with arbitrary-angle fixed |22> CP and adjacent SU2 plus diagonal qutrit phases',
          'baseline_cp22':360,'baseline_local_ops':1728,'optimized_cp22':rows[0]['cp22'],
          'optimized_local_ops':rows[0]['single_qutrit_ops'],
          'illustrative_entangler_duration_ns':580,'illustrative_local_duration_ns':30,
          'baseline_serial_duration_ns':360*580+1728*30,
          'optimized_serial_duration_ns':40*580+rows[0]['single_qutrit_ops']*30,
          'reduction_factor':360/rows[0]['cp22'],
          'one_control_SU3':4,'double_control_SU3':4,
          'independent_published_blocks_checked':9,
          'double_control_method':'SU3 cyclic-eigenbasis commutator, four single-control U3 operations, two CP22 per single control',
          'test_rows':rows,'approximate_options':approximations,
          'warning':'Not calibrated transmon gates; arbitrary-angle fixed CP22 remains an assumed primitive. Local rotations and swaps explicitly counted. No noise or pulse physics assumed for correctness.'}
    Path(args.output).write_text(json.dumps(data,indent=2,sort_keys=True)+'\n')
    print('COMPLETE 9-CONTEXT EXACT MULTIPLEXOR SU3 COMPILATION PASS',json.dumps({'cp22':rows[0]['cp22'],'reduction_factor':data['reduction_factor']}))

if __name__=='__main__':main()
