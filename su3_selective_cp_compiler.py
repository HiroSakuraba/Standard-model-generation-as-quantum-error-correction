"""Constructive three-qutrit SU(3) magnetic-block compilation to two-body gates.

Assumed gate set (not necessarily realized as a universal primitive on any
specific chip): arbitrary-angle, level-selective two-qutrit CP on L-T or T-R,
plus one-qutrit adjacent-level SU(2) rotations / diagonal phase operations.
No nearest-neighbor three-body gates and no direct L-R coupler are used.
Includes a three-transmon 4-level leakage/crosstalk toy model, not a calibrated
pulse simulation. Requires numpy and scipy.
"""
import argparse
import json
from collections import Counter
from pathlib import Path
import numpy as np
from scipy.linalg import expm

SQRT3 = 1/np.sqrt(3)
THIRD = 1/3
PUBLISHED = {
    (0,0):[[0,1,1],[1,0,1],[1,1,0]],
    (0,1):[[0,THIRD,SQRT3],[THIRD,0,SQRT3],[SQRT3,SQRT3,0]],
    (0,2):[[0,SQRT3,THIRD],[SQRT3,0,SQRT3],[THIRD,SQRT3,0]],
    (1,1):[[0,1/9,THIRD],[1/9,0,THIRD],[THIRD,THIRD,0]],
    (1,2):[[0,SQRT3/3,SQRT3/3],[SQRT3/3,0,THIRD],[SQRT3/3,THIRD,0]],
    (2,2):[[0,THIRD,1/9],[THIRD,0,THIRD],[1/9,THIRD,0]],
}
PUBLISHED.update({(r,l):v for (l,r),v in list(PUBLISHED.items()) if l!=r})


def target_block(left=1,right=2,dt=.19,phase=2*np.pi/3):
    m=np.asarray(PUBLISHED[left,right],complex)
    w=np.zeros((3,3),complex)
    for a in range(3): w[(a+1)%3,a]=m[(a+1)%3,a]
    b=np.exp(1j*phase)*w + np.exp(-1j*phase)*w.conj().T
    assert np.linalg.norm(b-b.conj().T)<1e-13 and abs(np.trace(b))<1e-13
    u=expm(1j*dt*b)
    assert abs(np.linalg.det(u)-1)<1e-12
    return b,u


def embed_rotation(two,levels):
    x=np.eye(3,dtype=complex)
    x[np.ix_(levels,levels)]=two
    return x


def decompose_su3(u):
    """U = G0^dag G1^dag G2^dag Z01 Z12; 5 SU(2) factors max.

    Each G only rotates adjacent 1-2 or 0-1 levels. Z rotations are also
    adjacent. This works because exp(i*traceless B) has det exactly 1.
    """
    q=u.copy(); eliminated=[]
    for a,b,col in ((1,2,0),(0,1,0),(1,2,1)):
        x,y=q[a,col],q[b,col]
        rad=np.hypot(abs(x),abs(y))
        if rad<1e-14:continue
        small=np.array([[np.conj(x),np.conj(y)],[-y,x]],complex)/rad
        g=embed_rotation(small,(a,b));q=g@q
        assert abs(np.linalg.det(small)-1)<2e-12
        eliminated.append(((a,b),small))
    assert np.max(np.abs(q-np.diag(np.diag(q))))<1e-12
    d=np.diag(q)
    alpha=np.angle(d)
    theta01=alpha[0]
    theta12=-alpha[2]
    z01=np.diag([np.exp(1j*theta01),np.exp(-1j*theta01)])
    z12=np.diag([np.exp(1j*theta12),np.exp(-1j*theta12)])
    assert np.linalg.norm(embed_rotation(z01,(0,1)) @ embed_rotation(z12,(1,2))-q)<3e-12
    chronological=[((1,2),z12),((0,1),z01)] + [((levels),g.conj().T) for levels,g in reversed(eliminated)]
    check=np.eye(3,dtype=complex)
    for levels,g in chronological: check=embed_rotation(g,levels)@check
    assert np.linalg.norm(check-u)<3e-12
    return chronological


def su2_commutator(rot):
    """Construct A,B in SU(2) with A B A^dag B^dag = rot."""
    assert abs(np.linalg.det(rot)-1)<1e-11
    vals,vec=np.linalg.eig(rot)
    # Choose orthonormal eigenvectors (normal matrix; degeneracy excluded here).
    vec[:,1]-=vec[:,0]*np.vdot(vec[:,0],vec[:,1]);vec[:,1]/=np.linalg.norm(vec[:,1])
    assert np.linalg.norm(vec.conj().T@vec-np.eye(2))<1e-11
    theta=np.angle(vals[0])
    b=vec@np.diag([np.exp(-.5j*theta),np.exp(.5j*theta)])@vec.conj().T
    a=vec@np.array([[0,1],[-1,0]],complex)@vec.conj().T
    assert np.linalg.norm(a@b@a.conj().T@b.conj().T-rot)<4e-11
    return a,b


def compile_controlled_su2(controller,state,levels,r):
    """One two-qutrit conditional SU2 using two selective CP entanglers.

    V diag(e^{i alpha},e^{i beta}) V^dag = r; single qutrit V cancels
    outside selected control sector. Time order: V^dag, CP_a, CP_b, V.
    """
    lam,vec=np.linalg.eig(r)
    vec[:,1]-=vec[:,0]*np.vdot(vec[:,0],vec[:,1]);vec[:,1]/=np.linalg.norm(vec[:,1])
    assert np.linalg.norm(vec.conj().T@vec-np.eye(2))<2e-11
    local=embed_rotation(vec,levels)
    return [('LOCAL', local.conj().T),('CP',(controller,state,levels[0],float(np.angle(lam[0])))),
            ('CP',(controller,state,levels[1],float(np.angle(lam[1])))),('LOCAL',local)]


def compile_doubly_controlled_su2(left_value,right_value,levels,rot):
    """Use commutator: [C_L(A),C_R(B)] implements C_L C_R(rot).

    Chronological order C_R(B^dag), C_L(A^dag), C_R(B), C_L(A).
    No ancilla or three-body primitive used.
    """
    a,b=su2_commutator(rot)
    schedule=[]
    for c,s,r in [('R',right_value,b.conj().T),('L',left_value,a.conj().T),
                  ('R',right_value,b),('L',left_value,a)]:
        schedule.extend(compile_controlled_su2(c,s,levels,r))
    return schedule


def compile_block(left=1,right=2,dt=.19,phase=2*np.pi/3):
    _,u=target_block(left,right,dt,phase)
    schedule=[]
    for levels,r in decompose_su3(u):
        schedule.extend(compile_doubly_controlled_su2(left,right,levels,r))
    return u,schedule


def cp_dense(controller,state,level,phi,dim=3):
    diag=np.ones(dim**3,dtype=complex)
    for l in range(dim):
        for t in range(dim):
            for r in range(dim):
                if t==level and ((controller=='L' and l==state) or (controller=='R' and r==state)):
                    diag[(l*dim+t)*dim+r]=np.exp(1j*phi)
    return np.diag(diag)


def local_dense(v,dim=3):
    k=np.eye(dim,dtype=complex)
    k[:3,:3]=v
    return np.kron(np.eye(dim),np.kron(k,np.eye(dim)))


def embedded_target(u,left=1,right=2):
    result=np.eye(27,dtype=complex)
    for i in range(3):
        for j in range(3):
            result[(left*3+i)*3+right,(left*3+j)*3+right]=u[i,j]
    return result


def ideal_compile_test(u,schedule,left,right):
    result=np.eye(27,dtype=complex)
    for typ,val in schedule:
        gate=local_dense(val) if typ=='LOCAL' else cp_dense(*val)
        result=gate@result
    expected=embedded_target(u,left,right)
    err=np.linalg.norm(result-expected,ord=2)
    assert err<2e-10,err
    return err


def two_ququart_leaky_cp(controller,state,level,phi,leak=.0,detune=1.4):
    """A phenomenological 4x4 x 4x4 unitary for a level-selective CP.

    The addressed pair uses states (ctrl,tgt) for L and (tgt,ctrl) for R.
    Leakage mixes |ctrl=2,tgt=1> with |ctrl=3,tgt=0> (modeled |21> <-> |30>),
    with a detuned level 3; unitary exponent, not a non-TP loss map.
    No claim that this represents a measured driven-transmon Hamiltonian.
    """
    h=np.zeros((16,16),complex)
    a=(state*4+level) if controller=='L' else (level*4+state)
    h[a,a]=phi
    if controller=='L':i,j=2*4+1,3*4+0
    else:i,j=1*4+2,0*4+3
    h[i,j]=h[j,i]=leak
    h[j,j]+=detune
    # Off-resonant detuning on |3> is idle-picture subtracted when leak=0.
    g=expm(1j*h)
    if abs(leak)<1e-15:
        # unwanted detuned |3> phase is outside computational basis and harmless
        assert abs(g[a,a]-np.exp(1j*phi))<1e-12
    if controller=='L':return np.kron(g,np.eye(4))
    return np.kron(np.eye(4),g)


def simulate_ququart(schedule,u,leak=0,crosstalk=0,angle_fraction=0):
    """Propagate 27 logical columns embedded into the full 4^3=64 space.

    Each native two-ququart CP pulse receives a coherent |21>--|30>
    leakage term and a separate unwanted spectator-pair phase after it.
    The latter is an explicit sequential crosstalk model, NOT driven dynamics.
    """
    computational=np.array([(l*4+t)*4+r for l in range(3) for t in range(3) for r in range(3)])
    k=np.zeros((64,27),complex);k[computational,np.arange(27)]=1
    counter=Counter()
    for typ,val in schedule:
        if typ=='LOCAL':
            gate=local_dense(val,4)
            counter['single_qutrit_U2_mixers']+=1
        else:
            ctrl,state,level,phi=val
            gate=two_ququart_leaky_cp(ctrl,state,level,phi*(1+angle_fraction),leak)
            counter['selective_two_ququart_CP']+=1
        k=gate@k
        if typ=='CP' and crosstalk:
            # A parasitic diagonal CP on the *other* neighbor edge.
            # Spectator couples to target level 2 while it is in level 1.
            diag=np.ones(64,complex)
            for l in range(4):
                for t in range(4):
                    for r in range(4):
                        bad=(t==2 and r==1) if val[0]=='L' else (l==1 and t==2)
                        if bad:diag[(l*4+t)*4+r]=np.exp(1j*crosstalk)
            k*=diag[:,None]
    target=embedded_target(u)
    logical=k[computational,:]
    leak_prob=float(1-np.linalg.norm(logical,ord='fro')**2/27)
    process_f=float(abs(np.trace(target.conj().T@logical))**2/27**2)
    conditioned_f=process_f/(1-leak_prob) if leak_prob<1 else 0
    return {'leakage_probability':leak_prob,'unconditional_entanglement_fidelity':process_f,
            'no_leakage_conditioned_entanglement_fidelity':conditioned_f,
            'native_counts':dict(counter),
            'operator_error_logical_fro':float(np.linalg.norm(logical-target,ord='fro'))}



def adjacent_swaps_to_two(level):
    """Chronological swaps taking a chosen logical qutrit level to 2."""
    swap01=np.array([[0,1,0],[1,0,0],[0,0,1]],complex)
    swap12=np.array([[1,0,0],[0,0,1],[0,1,0]],complex)
    if level==2:return []
    if level==1:return [swap12]
    if level==0:return [swap01,swap12]
    raise ValueError(level)


def compile_fixed_22(schedule):
    """Replace arbitrary-level CP_{s,t}(phi) by fixed CP_22(phi).

    Per pair: rotate requested control and target levels into |2>, apply
    a nearest-neighbor CP22, undo both permutations. Swaps only 01/12.
    Still assumes *arbitrary angle* for CP22 and addressable single levels.
    """
    out=[]
    for kind,val in schedule:
        if kind=='LOCAL':
            out.append(('LOCAL_AT',('T',val)))
        else:
            ctrl,state,target,phi=val
            for gate in adjacent_swaps_to_two(state):out.append(('LOCAL_AT',(ctrl,gate)))
            for gate in adjacent_swaps_to_two(target):out.append(('LOCAL_AT',('T',gate)))
            out.append(('CP22',(ctrl,phi)))
            for gate in reversed(adjacent_swaps_to_two(target)):out.append(('LOCAL_AT',('T',gate)))
            for gate in reversed(adjacent_swaps_to_two(state)):out.append(('LOCAL_AT',(ctrl,gate)))
    return out


def site_gate(g,site,dim=3):
    x=np.eye(dim,dtype=complex);x[:3,:3]=g
    factors=[x if letter==site else np.eye(dim) for letter in 'LTR']
    return np.kron(np.kron(factors[0],factors[1]),factors[2])


def fixed_cp22_gate(controller,phi,dim=3):
    return cp_dense(controller,2,2,phi,dim)


def ideal_fixed_22_test(schedule,u,left=1,right=2):
    # Target all 27 computational vectors, not only |1,t,2>.
    k=np.eye(27,dtype=complex)
    for typ,val in schedule:
        if typ=='LOCAL_AT':gate=site_gate(val[1],val[0]);k=gate@k
        else:gate=fixed_cp22_gate(*val);k=gate@k
    error=np.linalg.norm(k-embedded_target(u,left,right),ord=2)
    assert error<1e-10,error
    return error


def two_ququart_fixed_22(controller,phi,leak=0,detune=1.4):
    """4-level selective |22> CP with coherent |22> <-> |31> leakage.

    Excitation-conserving 2+2 = 3+1 two-body mixing is an illustrative
    transmon-like avoided-crossing channel, not derived from a device.
    """
    h=np.zeros((16,16),complex)
    h[2*4+2,2*4+2]=phi
    i,j=2*4+2,3*4+1
    h[i,j]=h[j,i]=leak
    h[j,j]=detune
    gate=expm(1j*h)
    return np.kron(gate,np.eye(4)) if controller=='L' else np.kron(np.eye(4),gate)


def simulate_fixed_22(schedule,u,leak=0,crosstalk=0,angle_fraction=0,target_27=None):
    comp=np.array([(l*4+t)*4+r for l in range(3) for t in range(3) for r in range(3)])
    k=np.zeros((64,27),complex);k[comp,np.arange(27)]=1
    counts=Counter()
    for typ,val in schedule:
        if typ=='LOCAL_AT':
            site,g=val
            k=site_gate(g,site,4)@k
            counts['single_qutrit_local_or_swap']+=1
            # Swap count is computed structurally by compile_fixed_22,
            # not guessed from whether a general mixer happens to be involutory.
        else:
            ctrl,phi=val
            k=two_ququart_fixed_22(ctrl,phi*(1+angle_fraction),leak)@k
            counts['fixed_CP22_entanglers']+=1
            if crosstalk:
                diag=np.ones(64,complex)
                for l in range(4):
                    for t in range(4):
                        for r in range(4):
                            bad=t==2 and r==1 if ctrl=='L' else l==1 and t==2
                            if bad:diag[(l*4+t)*4+r]=np.exp(1j*crosstalk)
                k*=diag[:,None]
    logical=k[comp,:]
    target=embedded_target(u) if target_27 is None else target_27
    leakage=float(max(0,1-np.linalg.norm(logical,'fro')**2/27))
    fe=float(abs(np.trace(target.conj().T@logical))**2/27**2)
    return {'leakage_probability':leakage,'unconditional_entanglement_fidelity':fe,
            'conditional_no_leakage_entanglement_fidelity':fe/(1-leakage) if leakage<1 else 0,
            'operator_error_logical_fro':float(np.linalg.norm(logical-target,'fro')),
            'counts':dict(counts)}


def compile_full_interior(dt=.19,phase=2*np.pi/3,run_fault_test=True):
    """All nine neighbor-conditioned blocks of one three-qutrit magnetic U.

    The context pieces commute logically because their neighbor projectors
    are disjoint; all gates remain nearest-neighbor 2-local + 1-local.
    """
    full=np.zeros((27,27),complex)
    big=[];context_counts={}
    for l in range(3):
        for r in range(3):
            u,schedule=compile_block(l,r,dt,phase)
            full+=embedded_target(u,l,r)-np.eye(27)
            built=compile_fixed_22(schedule)
            c=Counter(x[0] for x in built)
            swap_count=sum(2*(len(adjacent_swaps_to_two(state))+len(adjacent_swaps_to_two(target)))
                           for typ, payload in schedule if typ=='CP' for ctrl,state,target,phi in [payload])
            assert c['LOCAL_AT']==len([x for x in schedule if x[0]=='LOCAL'])+swap_count
            context_counts[f'{l},{r}']={'cp22':c['CP22'],'local':c['LOCAL_AT'],'adjacent_swaps':swap_count}
            big.extend(built)
    # full=I+sum_context(U_context-I), exactly one context active per basis state
    full+=np.eye(27)
    bm=np.zeros((27,27),complex)
    for l in range(3):
        for r in range(3):
            b,_=target_block(l,r,dt,phase)
            for a in range(3):
                for c in range(3):bm[(l*3+a)*3+r,(l*3+c)*3+r]=b[a,c]
    target_full=expm(1j*dt*bm)
    assert np.linalg.norm(full-target_full)<3e-12
    k=np.eye(27,dtype=complex)
    for typ,val in big:
        k=(site_gate(val[1],val[0]) if typ=='LOCAL_AT' else fixed_cp22_gate(*val))@k
    error=np.linalg.norm(k-target_full,ord=2)
    assert error<2e-10,error
    cost=Counter(x[0] for x in big)
    fault=simulate_fixed_22(big,None,leak=.01,crosstalk=.002,target_27=target_full) if run_fault_test else None
    return {'full_27D_operator_error':float(error),'context_gate_counts':context_counts,
            'full_CP22_count':cost['CP22'],'full_single_qutrit_ops':cost['LOCAL_AT'],
            'full_swap_count':sum(z['adjacent_swaps'] for z in context_counts.values()),
            'total_contexts':len(context_counts),'published_blocks':len(PUBLISHED),
            'full_4_level_toy_error':fault}

def repo_table_crosscheck():
    """Independent previous publication transcription, when present."""
    try:
        from su3_three_plaquette_chain import PUBLISHED as source_table
    except ImportError:
        print('NOTE: previous publication reference file not importable; internal 9-block tests remain')
        return False
    assert set(source_table)==set(PUBLISHED)
    for context, block in PUBLISHED.items():
        assert np.linalg.norm(np.asarray(source_table[context],float)-np.asarray(block,float))<1e-13
    return True


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',default='su3_selective_cp_compiler_results.json');args=parser.parse_args()
    parent_benchmark_verified=repo_table_crosscheck()
    # Independent regression across different blocks, couplings, and phases.
    for ll,rr,tt,pp in ((0,0,.08,0),(1,1,.3,.23),(2,2,.2,-1.1),(0,2,.047,2.094395)):
        uu,ss=compile_block(ll,rr,tt,pp)
        ideal_compile_test(uu,ss,ll,rr)
        ideal_fixed_22_test(compile_fixed_22(ss),uu,ll,rr)
    left,right=1,2;dt=.19;phase=2*np.pi/3
    b,u=target_block(left,right,dt,phase)
    factors=decompose_su3(u)
    sched_u,sequence=compile_block(left,right,dt,phase)
    assert np.linalg.norm(sched_u-u)<1e-12
    exact_error=ideal_compile_test(u,sequence,left,right)
    count=Counter(x[0] for x in sequence)
    fixed_schedule=compile_fixed_22(sequence)
    fixed_error=ideal_fixed_22_test(fixed_schedule,u)
    fixed_count=Counter(x[0] for x in fixed_schedule)
    assert fixed_count["CP22"]==count["CP"]
    total_swaps=sum(2*(len(adjacent_swaps_to_two(state))+len(adjacent_swaps_to_two(target)))
                    for typ,payload in sequence if typ=='CP' for ctrl,state,target,phi in [payload])
    assert fixed_count['LOCAL_AT']==count['LOCAL']+total_swaps
    assert count['CP']==8*len(factors) and count['LOCAL']==8*len(factors)
    # All gates use only L-T or T-R connectivity and act on <=2 physical sites.
    assert all(x[1][0] in ('L','R') for x in sequence if x[0]=='CP')
    checks=[simulate_ququart(sequence,u,leak=x,crosstalk=y,angle_fraction=z)
            | {'leak_amplitude':x,'spectator_crosstalk_phase':y,'relative_CP_angle_error':z}
            for x,y,z in [(0,0,0),(.01,0,0),(.03,0,0),(.01,.002,0),(.01,0,.01),(.03,.006,.02)]]
    assert abs(checks[0]['unconditional_entanglement_fidelity']-1)<1e-10
    assert abs(checks[0]['leakage_probability'])<1e-10
    assert checks[-1]['leakage_probability']>0
    fixed_checks=[simulate_fixed_22(fixed_schedule,u,leak=x,crosstalk=y,angle_fraction=z)
                  | {'leak_amplitude':x,'spectator_crosstalk_phase':y,'relative_CP_angle_error':z}
                  for x,y,z in [(0,0,0),(.01,0,0),(.03,0,0),(.01,.002,0),(.01,0,.01),(.03,.006,.02)]]
    assert abs(fixed_checks[0]['unconditional_entanglement_fidelity']-1)<1e-10
    assert fixed_checks[-1]['leakage_probability']>0
    full_audit=compile_full_interior(dt=dt,phase=phase)
    # Conservative time-accounting scenario, NOT a calibrated arbitrary-angle selective pulse.
    times={'selective_CP_ns':580.,'one_qutrit_SU2_mixer_ns':30.}
    nominal_duration_ns=fixed_count['CP22']*times['selective_CP_ns']+fixed_count['LOCAL_AT']*times['one_qutrit_SU2_mixer_ns']
    # Fixed CP22 routing via only local adjacent-level swaps is explicitly counted.
    # Each arbitrary adjacent-level local SU2 is still counted as a composite mixer;
    # its virtual-Z/drive calibration and waveform duration are NOT compiled.
    result={'source':'Chen, Mueller, Yao PRD 114 074502 Table VI; experimental qutrit CP references are not the assumed universal selective CP',
            'target_context':[left,right],'dt':dt,'magnetic_phase_rad':phase,
            'local_block_B':np.round(b,12).tolist() if False else [[(float(z.real),float(z.imag)) for z in row] for row in b],
            'previous_published_transcription_crosschecked':parent_benchmark_verified,
            'number_su2_factors':len(factors),'native_assumed_selective_CP':count['CP'],'native_single_qutrit_mixers':count['LOCAL'],
            'max_nonadjacent_single_qutrit_mixers':0,'fixed_CP22_native_gate_counts':dict(fixed_count),
            'fixed_CP22_27D_error':float(fixed_error),'fixed_CP22_adjacent_level_swap_count':total_swaps,
            'full_interior_plaquette':full_audit,'exact_27D_operator_error':float(exact_error),
            'four_level_states':64,'computational_states':27,
            'illustrative_unverified_gate_times_ns':times,'illustrative_total_nominal_duration_ns':nominal_duration_ns,
            'noise_model':'Each CP uses a 4-level pair Hamiltonian with resonant-style |2,1> <-> |3,0> coupling (dimensionless), separate spectator conditional phase, and fractional intended CP-phase bias. No stochastic relaxation, waveform, crosstalk calibration or qutrit coherence simulation.',
            'test_rows':checks,'fixed_CP22_four_level_rows':fixed_checks,'warning':'Assumed arbitrary-angle two-qutrit CP22 and uncounted internal arbitrary-SU2 waveform decomposition. CP22 duration 580ns is only an illustrative choice anchored to a different demonstrated CZ gate; not a physical gate-time lower bound or calibration.'}
    Path(args.output).write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
    print('PARENT PUBLISHED TABLE VI CROSSCHECK:',parent_benchmark_verified)
    print('SU3 3-QUTRIT UNITARY MATRIX MATCH:',exact_error)
    print('Gates:',dict(count),'CP22 compiled:',dict(fixed_count),'SU2 factors:',len(factors),'scenario_gate_time_ns:',nominal_duration_ns)
    print('ALL NINE NEIGHBOR BLOCKS FULL UNITARY:',json.dumps(full_audit,sort_keys=True))
    for row in checks:print('TEST_SELECTIVE',json.dumps(row,sort_keys=True))
    for row in fixed_checks:print('TEST_FIXED_CP22',json.dumps(row,sort_keys=True))
    print('ALL CONSTRUCTIVE 3-QUTRIT COMPILER/4-LEVEL UNITARY TESTS PASSED')

if __name__=='__main__':main()
