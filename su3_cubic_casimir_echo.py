"""Three-phase SU(3) electric cubic-Casimir refocusing on 5-face minimal chain.

Sources: Chen, Mueller, Yao, Phys Rev D 114, 074502 (2026) for H.
This work's clock pulses are an *assumed additional control*, not from paper.
Only modeled gauge-invariant open-chain qutrits; no physical hardware gate.
Reproduce: OPENBLAS_NUM_THREADS=1 python su3_cubic_casimir_echo.py
"""
import json
import numpy as np
import sympy as sp
from scipy.sparse import csr_matrix, diags
from scipy.sparse.linalg import expm_multiply
from su3_four_five_plaquette import basis_states, electric_links, one_site_actions, dense_matrices
from su3_channel_recovery import source_data, syndrome_decoder, recover_fidelity

N=5
DIM=3**N
G2=1.0
TIME=0.6
ROOT=np.exp(2j*np.pi/3)


def construct():
    states=tuple(basis_states(N))
    pos={q:i for i,q in enumerate(states)}
    K=np.zeros((DIM,2),dtype=complex)
    for j,q in enumerate(((0,0,0,0,1),(0,0,0,0,2))):K[pos[q],j]=1
    d=np.array([float(sum(electric_links(q))) for q in states])
    W=[]
    for site in range(N):
        rows=[];cols=[];data=[]
        for j,q in enumerate(states):
            for p,a in one_site_actions(q,site,'W').items():
                rows.append(pos[p]);cols.append(j);data.append(float(a))
        W.append(csr_matrix((data,(rows,cols)),shape=(DIM,DIM),dtype=complex))
    H=diags(d,format='csr',dtype=complex)-sum((w+w.getH() for w in W),csr_matrix((DIM,DIM),dtype=complex))
    # Additive source magnetic 6N I/g² is global phase; omit.
    # Native parameter g²=1 and independent validation vs published matrices.
    _, old_D, old_B=dense_matrices(N)
    assert np.max(np.abs(d-old_D))<1e-12
    assert np.linalg.norm(H.toarray()-(np.diag(old_D)-sum(old_B)),ord=np.inf)<1e-12
    assert np.linalg.norm(H.toarray()-H.toarray().conj().T,ord=np.inf)<1e-12
    return states,K,d,W,H


def exact_clock_test(states):
    # SU(3) C3(p,q)=(p-q)(2p+q+3)(p+2q+3)/18
    def cubic(p,q):return sp.Rational((p-q)*(2*p+q+3)*(p+2*q+3),18)
    assert [cubic(0,0),cubic(1,0),cubic(0,1)]==[0,sp.Rational(10,9),-sp.Rational(10,9)]
    # The pulse is Q=exp[(2pi i/3)*(9/10)sum_i C3(top_link_i)].
    # Every directed W_i raises the corresponding q_i by +1 mod 3.
    crossings=0
    for q in states:
        for site in range(N):
            for new,a in one_site_actions(q,site,'W').items():
                if a:
                    assert (sum(new)-sum(q))%3==1
                    crossings+=1
    assert crossings>500
    # sum_{k=0,1,2} Q^k H Q^{-k}/3 = D exactly,
    # because the directed loop picks up phase omega^k; its adjoint omega^-k.
    z=(-1+sp.I*sp.sqrt(3))/2
    assert sp.simplify(1+z+z*z)==0
    # Electric Casimirs only distinguish occupied / unoccupied states,
    # so a quadratic-electric pulse cannot change the relative phase of |+> and |->.
    q_plus=(0,0,0,0,1);q_minus=(0,0,0,0,2)
    assert electric_links(q_plus)==electric_links(q_minus)
    # Yet the last B4 has a nonzero X matrix element in this pair.
    assert one_site_actions(q_plus,4,'B')[q_minus]==1
    assert one_site_actions(q_minus,4,'B')[q_plus]==1
    print('EXACT: C3(3)=+10/9, C3(bar3)=-10/9; all',crossings,'directed W_i transitions have Z3 clock charge 1')
    print('EXACT: 3-phase average cancels every directed fundamental magnetic W_i and adjoint; quadratic electric cannot refocus logical X')


def cycle_hamiltonians(d,W,states,amp_error=0.0):
    s=np.array([sum(0 if x==0 else (1 if x==1 else -1) for x in q) for q in states],dtype=float)
    # Clock eigenphase: +2pi/3 * sum (9/10 C3) per top link.
    phase=(2*np.pi/3)*(1+amp_error)*s
    H0=diags(d,format='csr',dtype=complex)-sum((w+w.getH() for w in W),csr_matrix((DIM,DIM),dtype=complex))
    coo=H0.tocoo()
    out=[]
    for k in range(3):
        phase_factor=np.exp(1j*k*(phase[coo.row]-phase[coo.col]))
        hk=csr_matrix((coo.data*phase_factor,(coo.row,coo.col)),shape=(DIM,DIM))
        assert np.max(np.abs((hk-hk.getH()).data),initial=0)<1e-12
        out.append(hk)
    return out


def evolution(y,hams,time,cycles,scheme='symmetric'):
    assert cycles>=1
    if scheme=='symmetric':order=(0,1,2,1,0);weights=(1/6,1/6,1/3,1/6,1/6)
    elif scheme=='first_order':order=(0,1,2);weights=(1/3,1/3,1/3)
    else:raise ValueError(scheme)
    out=np.array(y,dtype=complex,copy=True)
    dt=time/cycles
    for _ in range(cycles):
        for k,a in zip(order,weights):
            out=expm_multiply((-1j*dt*a)*hams[k],out,traceA=0)
    return out


def finite_pulse_sequence(y,H,states,time,cycles,tau):
    """Real-time lab-frame sequence: native H remains ON during 4 pulses/cycle.

    Each signed cubic-Casimir pulse has duration tau, magnitude theta/tau
    on each of the five top links; net top-link clock is exactly identity
    in the zero-pulse-width limit. Elapsed duration is time+4*cycles*tau.
    """
    s=np.array([sum(0 if q==0 else (1 if q==1 else -1) for q in face) for face in states],float)
    clock=diags(s,format='csr')
    theta=2*np.pi/3
    out=np.array(y,dtype=complex,copy=True)
    dt=time/cycles
    for _ in range(cycles):
        for k,w in enumerate((1/6,1/6,1/3,1/6,1/6)):
            out=expm_multiply(-1j*(w*dt)*H,out,traceA=0)
            if k<4:
                # Physical toggle: Q^dag, Q^dag, Q, Q, net identity.
                sign=(1 if k<2 else -1)
                pulse=H+(sign*theta/tau)*clock
                out=expm_multiply(-1j*tau*pulse,out,traceA=0)
    return out


def binary_refocusing(y,H,states,time,cycles):
    """Negative control: quadratic-electric phase (-1)^{occupied face links}."""
    occ=np.array([sum(q!=0 for q in face) for face in states])
    phase=np.pi*occ
    coo=H.tocoo()
    H1=csr_matrix((coo.data*np.exp(1j*(phase[coo.row]-phase[coo.col])),
                    (coo.row,coo.col)),shape=H.shape)
    out=np.array(y,dtype=complex,copy=True)
    for _ in range(cycles):
        for mat,w in zip((H,H1,H),(1/4,1/2,1/4)):
            out=expm_multiply(-1j*(time/cycles*w)*mat,out,traceA=0)
    return out


def entanglement_fidelity(maps,K,V=None):
    if V is None:
        fidelity=sum(abs(np.trace(K.conj().T@m))**2/4 for m in maps)
        return float(fidelity)
    return recover_fidelity(maps,V,K)[0]


def jump_operator(states,p,site,directed='W'):
    pos={q:i for i,q in enumerate(states)}
    rows=[];cols=[];data=[]
    for src,q in enumerate(states):
        for dst,a in one_site_actions(q,site,directed).items():
            rows.append(pos[dst]);cols.append(src);data.append(float(a))
    W=csr_matrix((data,(rows,cols)),shape=(DIM,DIM))
    effect=np.asarray(W.power(2).sum(axis=0)).flatten()
    assert effect.max()<=1+1e-12
    K0=diags(np.sqrt(1-p*effect),format='csr')
    K1=np.sqrt(p)*W
    completeness=K0.getH()@K0+K1.getH()@K1
    assert np.max(np.abs((completeness-diags(np.ones(DIM),format='csr')).data),initial=0)<2e-12
    return K0,K1


def run():
    states,K,d,W,H=construct()
    exact_clock_test(states)
    V,_,g=syndrome_decoder(*source_data())
    assert g.rank()==5
    controlled=cycle_hamiltonians(d,W,states)
    Dmat=diags(d,format='csr')
    averaged=sum(controlled,csr_matrix((DIM,DIM),dtype=complex))/3
    assert np.max(np.abs((averaged-Dmat).data),initial=0)<2e-13
    print('NUMERICAL twirled Hamiltonian: max residual',np.max(np.abs((averaged-Dmat).data),initial=0))
    rows=[]
    nat=expm_multiply((-1j*TIME)*H,K,traceA=0)
    before_raw=entanglement_fidelity([nat],K)
    before_recovered=entanglement_fidelity([nat],K,V)
    print('NATIVE no-control t=.6: F_raw',before_raw,'F_fixed_recovery',before_recovered)
    for m in (1,2,4,8,16,32):
        seq=evolution(K,controlled,TIME,m)
        fraw=entanglement_fidelity([seq],K)
        frec=entanglement_fidelity([seq],K,V)
        leakage=1-np.real(np.trace((K.conj().T@seq).conj().T@(K.conj().T@seq)))/2
        row={'cycles':m,'global_clock_pulses':4*m,'t':TIME,'clock_calibration_error':0.0,'F_raw':fraw,'F_fixed_recovery':frec,'code_leakage':float(leakage)}
        rows.append(row)
        print('SYMMETRIC',json.dumps(row,sort_keys=True))
        assert fraw>before_raw
    binary=binary_refocusing(K,H,states,TIME,8)
    binary_row={'control':'quadratic_occupation_Z2','cycles':8,'F_raw':entanglement_fidelity([binary],K),'F_fixed_recovery':entanglement_fidelity([binary],K,V)}
    print('NEGATIVE_CONTROL_BINARY',json.dumps(binary_row,sort_keys=True))
    rows.append(binary_row)
    assert binary_row['F_fixed_recovery']<0.8
    ideal8=evolution(K,controlled,TIME,8)
    ultra_short=finite_pulse_sequence(K,H,states,TIME,8,1e-6)
    # Source magnetic constant and unrefocused electric control-phase time
    # give a physically irrelevant global phase as pulse duration tends to zero.
    global_phase=np.angle(np.vdot(ideal8[:,0],ultra_short[:,0]))
    assert np.linalg.norm(ultra_short*np.exp(-1j*global_phase)-ideal8)<4e-5
    # Nonzero-duration cubic electric pulses with full published H left on.
    for tau in (0.0005,0.002,0.01):
        real=finite_pulse_sequence(K,H,states,TIME,8,tau)
        row={'pulse_duration':tau,'control_amplitude_per_top_link':float(2*np.pi/(3*tau)),
             'cycles':8,'pulse_count':32,'elapsed_time':TIME+32*tau,
             'F_raw':entanglement_fidelity([real],K),
             'F_fixed_recovery':entanglement_fidelity([real],K,V)}
        print('FINITE_DURATION_PULSE',json.dumps(row,sort_keys=True))
        rows.append(row)
        assert row["F_raw"]>0.98
    # Fitted clock control errors; pulse gate acts instantaneously but has a coherent angle error.
    for error in (0.01,0.05,0.10):
        f=evolution(K,cycle_hamiltonians(d,W,states,error),TIME,8)
        row={'clock_calibration_error':error,'cycles':8,'F_raw':entanglement_fidelity([f],K),'F_fixed_recovery':entanglement_fidelity([f],K,V)}
        print('PULSE_CALIBRATION',json.dumps(row,sort_keys=True))
        rows.append(row)
    # Show physical CPTP faults are NOT erased by coherent Hamiltonian refocusing.
    for site in (3,4):
        jump=jump_operator(states,0.1,site,'W')
        first=evolution(K,controlled,TIME/2,8)
        maps=[evolution(op@first,controlled,TIME/2,8) for op in jump]
        row={'jump_site':site,'jump_probability':0.10,'cycles_per_half':8,'t':TIME,'F_raw':entanglement_fidelity(maps,K),'F_fixed_recovery':entanglement_fidelity(maps,K,V)}
        print('UNCONDITIONED_DIRECTED_JUMP',json.dumps(row,sort_keys=True))
        rows.append(row)
        if site==4:assert row['F_fixed_recovery']<0.95
    assert before_recovered<0.7 and max(r.get('F_raw',0) for r in rows)>0.99
    print('ALL CUBIC-CASIMIR Z3 DECOUPLING TESTS PASS')
    result={'dimension':DIM,'g2':G2,'t':TIME,'native_F_raw':before_raw,'native_F_fixed_recovery':before_recovered,'control_is_instantaneous':True,'rows':rows}
    with open('su3_cubic_casimir_echo_results.json','w') as f:json.dump(result,f,indent=2,sort_keys=True)
    return result

if __name__=='__main__':run()
