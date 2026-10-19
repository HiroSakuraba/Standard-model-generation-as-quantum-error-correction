"""Explicit recovery and CPTP physical-channel stress for a five-face SU(3) toy code.

Builds the syndrome isometries from the exactly KL-certified 21-element
error span in su3_four_five_plaquette.py. Tests the recovery against (a)
a code-conditioned CPTP demonstration channel supported in that span,
(b) non-conditioned, gauge-invariant directed-Wilson jump channels, and
(c) time evolution under the published five-plaquette Hamiltonian before,
between and after two faults.

Requires numpy scipy sympy, existing research source files. No hardware
recovery circuit or glueball mass claim. The optional inverse-Hamiltonian
operation is an explicitly idealized *benchmark*, not free recovery.
"""
import json
import numpy as np
from scipy.linalg import eigh
import sympy as sp
from su3_four_five_plaquette import (
    basis_states, dense_matrices, one_site_actions, applied_error, electric_links,
    kl_check,
)

N=5
DIM=3**N
CODE=((0,0,0,0,1),(0,0,0,0,2))
ERRORS=[('I',None)]+[('E',i) for i in range(3*N+1)]+[('B',j) for j in range(N-1)]


def source_data():
    labels=basis_states(N)
    pos={s:i for i,s in enumerate(labels)}
    K=np.zeros((DIM,2))
    for i,q in enumerate(CODE):K[pos[q],i]=1
    return labels,pos,K


def matrix_error(labels,pos,kind,site):
    mat=np.zeros((DIM,DIM))
    for col,q in enumerate(labels):
        for state,a in applied_error(q,kind,site).items():
            mat[pos[state],col]=float(a)
    return mat


def syndrome_decoder(labels,pos,K):
    """KL diagonalization yields 5 mutually orthogonal two-column isometries."""
    passed,exact_gram=kl_check(CODE,ERRORS)
    assert passed and exact_gram.rank()==5 and len(ERRORS)==21
    images=[matrix_error(labels,pos,*e)@K for e in ERRORS]
    left=np.column_stack([x[:,0] for x in images])
    right=np.column_stack([x[:,1] for x in images])
    assert np.max(np.abs(left.T@left-right.T@right))<1e-13
    assert np.max(np.abs(left.T@right))<1e-13
    gram=left.T@left
    eig,rotation=eigh(gram)
    select=np.flatnonzero(eig>1e-9)
    assert len(select)==5, eig
    V=[]
    for k in select:
        a=left@rotation[:,k]/np.sqrt(eig[k]);b=right@rotation[:,k]/np.sqrt(eig[k])
        V.append(np.column_stack((a,b)))
    stack=np.column_stack(V)
    assert np.linalg.norm(stack.T@stack-np.eye(2*len(V)))<2e-12
    assert np.linalg.norm(stack@stack.T@stack@stack.T-stack@stack.T)<2e-12
    return V,images,exact_gram


def recover_fidelity(output_kraus,V,K):
    """Logical recovered entanglement fidelity for a qubit CPTP channel.

    Kraus arguments are DIM x 2 matrices A_j=K_j K_code. Explicit CPTP
    recovery is R(X)=sum_r V_r^dag X V_r + Tr[(I-S)X] I_2/2, S=sum V_rV_r^dag.
    Entanglement fidelity uses all Kraus operators of the recovery including
    the positive complement replacement; the complement is NOT discarded.
    """
    total_norm=sum(float(np.vdot(m,m).real) for m in output_kraus)
    assert abs(total_norm-2)<2e-10, total_norm
    coherent=0.0;inside=0.0
    for m in output_kraus:
        for v in V:
            a=v.T.conj()@m
            coherent+=abs(np.trace(a))**2/4
            inside+=np.linalg.norm(a)**2
    outside=max(0., total_norm-inside)
    fidelity=coherent+outside/8
    raw=sum(abs(np.trace(K.T.conj()@m))**2 for m in output_kraus)/4
    assert -1e-10<=fidelity<=1+1e-10
    return float(fidelity),float(raw),float(outside/2)


def certified_cptp_channel(images,exact_gram,K):
    """EXPLICITLY code-conditioned channel realizing only certified errors.

    E_a P are Kraus building blocks, scaled by the exact code-restricted
    norm c_aa; K_perp=I-P gives a trace-preserving extension on full space.
    This nonlocal code projection is NOT a physical local-noise mechanism.
    """
    c=np.array([float(exact_gram[i,i]) for i in range(len(ERRORS))])
    active=np.flatnonzero(c>1e-14)
    assert len(active)>=5
    probs=np.zeros(len(ERRORS))
    probs[0]=.55
    others=[i for i in active if i!=0]
    probs[others]=.45/len(others)
    assert abs(probs.sum()-1)<1e-13
    output=[np.sqrt(probs[i]/c[i])*images[i] for i in active]
    effects=sum(m.T@m for m in output)
    assert np.linalg.norm(effects-np.eye(2))<1e-12
    # Global extension: K_perp=I-KK^dag; other K_j=(E_j K)K^dag.
    # Sum K_j^dag K_j = KK^dag + (I-KK^dag)=I.
    return output,len(active)


def recovered_matrix_unit(i,j,maps,V):
    """Independently verify the full logical superoperator, not fidelity alone."""
    output=np.zeros((2,2),complex)
    for m in maps:
        sub=[v.T.conj()@m for v in V]
        for a in sub:
            output+=np.outer(a[:,i],np.conj(a[:,j]))
        complement=np.vdot(m[:,j],m[:,i])-sum(np.vdot(a[:,j],a[:,i]) for a in sub)
        output+=np.eye(2)*complement/2
    return output


def hermitian_sign_kick(hermitian,theta):
    """Two equally likely coherent gauge-invariant ±theta B_i rotations."""
    e,v=eigh(hermitian)
    result=[]
    for sign in (1,-1):
        U=(v*np.exp(1j*sign*theta*e))@v.T
        assert np.linalg.norm(U.conj().T@U-np.eye(DIM))<1e-10
        result.append(U/np.sqrt(2))
    assert np.linalg.norm(sum(x.conj().T@x for x in result)-np.eye(DIM))<1e-10
    return result


def wilson_jump(labels,pos,p,site,directed='W'):
    """Local 2-Kraus CPTP Wilson jump, K1=sqrt(p)W_i, K0=sqrt(I-pW_i^dag W_i).

    Source-truncated W_i is NOT unitary, nor is this dissipative jump a
    measured rate from the published Hamiltonian. Values 0<=p<=1 work
    because its source transition amplitudes have magnitude <=1.
    """
    W=matrix_error(labels,pos,directed,site)
    diagonal=np.sum(W**2,axis=0)
    assert np.linalg.norm(W.T@W-np.diag(diagonal))<1e-13
    assert np.max(diagonal)<=1+1e-13
    assert 0<=p<=1
    K0=np.diag(np.sqrt(1-p*diagonal));K1=np.sqrt(p)*W
    assert np.linalg.norm(K0.T@K0+K1.T@K1-np.eye(DIM))<2e-12
    return (K0,K1)


def evolved_two_event_channel(K,noise,step,nominal_compensation):
    """Three equal ideal-evolution segments; two local quantum jumps between them."""
    # Evolved physical circuit U_step E_b U_step E_a U_step.
    U=step
    ideal=U@U@U
    undo=ideal.T.conj() if nominal_compensation else np.eye(DIM)
    results=[]
    for a in noise:
        middle=U@a@U@K
        for b in noise:
            results.append(undo@U@b@middle)
    return results


def run():
    labels,pos,K=source_data()
    V,images,gram=syndrome_decoder(labels,pos,K)
    allowed,nactive=certified_cptp_channel(images,gram,K)
    F,raw,out=recover_fidelity(allowed,V,K)
    assert abs(F-1)<5e-12
    for i in range(2):
        for j in range(2):
            expected=np.zeros((2,2),complex);expected[i,j]=1
            assert np.linalg.norm(recovered_matrix_unit(i,j,allowed,V)-expected)<2e-12
    print('EXACT-LINEAR-ERROR-SPAN RECOVERY:',len(ERRORS),'operators, Gram rank',len(V),'active channel branches',nactive,'F_e',F)

    basis,D,B=dense_matrices(N)
    assert basis==labels
    # Physical square-chain H at g²=1; additive 30 I can be ignored in evolution.
    H=np.diag(D)-sum(B)
    eigen,U=eigh(H)
    rows=[]
    for time in (0.,.15,.6):
        timestep= np.exp(-1j*eigen*(time/3))
        step=(U*timestep)@U.T
        if time==0: assert np.linalg.norm(step-np.eye(DIM))<3e-12
        # Baseline native Hamiltonian is NOT a quantum memory: compare with
        # ideal Hamiltonian inverse, which requires global calibrated control.
        for fault in [('W',3),('Wdag',3),('W',4)]:
            for p in (0.,.02,.1):
                noise=wilson_jump(labels,pos,p,fault[1],fault[0])
                kraus=evolved_two_event_channel(K,noise,step,True)
                fidelity,raw,outs=recover_fidelity(kraus,V,K)
                row={'t':time,'jump':fault[0]+str(fault[1]),'p_per_event':p,'F_recovered':fidelity,'F_raw_after_ideal_inverse':raw,'syndrome_outside_probability':outs}
                rows.append(row)
                print('CHANNEL',json.dumps(row,sort_keys=True))
                if p==0:assert abs(fidelity-1)<3e-10
        if time in (0.,.6):
            for theta in (.1,.3):
                kicks=hermitian_sign_kick(B[3],theta)
                result=evolved_two_event_channel(K,kicks,step,True)
                f,raw,out=recover_fidelity(result,V,K)
                print('HERMITIAN_KICK',json.dumps({'t':time,'theta':theta,'plaquette':3,'F_recovered':f,'F_raw':raw,'outside':out},sort_keys=True))
            jump=wilson_jump(labels,pos,0.,3,'W')
            uncorrected=evolved_two_event_channel(K,jump,step,False)
            f,raw,out=recover_fidelity(uncorrected,V,K)
            print('NATIVE_EVOLUTION_WITHOUT_INVERSE',json.dumps({'t':time,'F_recovered':f,'F_raw':raw,'outside':out},sort_keys=True))
    print('ALL SU3 FIVE-PLAQUETTE EXPLICIT CPTP CHANNEL/RECOVERY CHECKS PASSED')
    print('JSON_RESULTS='+json.dumps({'dimension':DIM,'exact_error_count':len(ERRORS),'syndromes':len(V),'code_conditioned_F':F,'rows':rows},sort_keys=True))

if __name__=='__main__':run()
