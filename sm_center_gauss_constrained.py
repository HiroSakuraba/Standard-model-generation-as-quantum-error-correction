"""Finite-group Gauss-constrained four-link model for the SM centre echo.

A closed oriented four-link loop carries A-flux in Z_(6/q), q=|Gamma|,
with Gauss div(A)=0, plus B-flux in Z2 and matter occupations m_v with
Gauss div(B)=m_v mod 2. Matter hops flip the connecting B-link, not A.
The Z6 centre character acts on A-flux. This is an explicitly constrained
FINITE-GROUP analogue, not a full SU3xSU2xU1 theory, nor real SM fermions.

For Gamma=Z6 only A=0 is allowed. All centre pulses become the identity,
and centre-charged magnetic loops are absent (negative control).
"""
from __future__ import annotations
from itertools import product
import json
from pathlib import Path
import numpy as np
from scipy.linalg import expm

ORDERS=(1,2,3,6)
V=4

def valid(a,b,m,d):
    return (all((a[i]-a[(i-1)%V])%d==0 for i in range(V))
            and all((b[i]-b[(i-1)%V]-m[i])%2==0 for i in range(V)))

def basis(q):
    if q not in ORDERS:raise ValueError(q)
    d=6//q
    result=[]
    for v in range(d):
        a=(v,)*V
        for m in product((0,1),repeat=V):
            if sum(m)%2:continue
            for first in (0,1):
                b=[first]
                for i in range(1,V): b.append((b[-1]+m[i])%2)
                st=(a,tuple(b),tuple(m))
                assert valid(*st,d)
                result.append(st)
    assert len(result)==16*d
    return result

def operators(q):
    states=basis(q); pos={s:i for i,s in enumerate(states)}
    d=6//q;n=len(states)
    E=np.zeros((n,n),complex);M=np.zeros_like(E);G=np.zeros_like(E)
    phase=np.zeros(n)
    hop_count=0
    for col,(a,b,m) in enumerate(states):
        flux=a[0]
        phase[col]=q*flux
        E[col,col]=.18*min(flux,d-flux)**2+.09*sum(b)+.17*sum(m)+.037*flux*sum(m)
        if d>1:
            nxt=(tuple((x+1)%d for x in a),b,m)
            assert valid(*nxt,d)
            G[pos[nxt],col]-=.63
        for v in range(V):
            j=(v+1)%V
            if m[v]!=m[j]:
                mm=list(m);mm[v],mm[j]=mm[j],mm[v]
                bb=list(b);bb[v]^=1
                nxt=(a,tuple(bb),tuple(mm))
                assert valid(*nxt,d) and nxt in pos
                M[pos[nxt],col]-=.31
                hop_count+=1
    if d>1:G+=G.conj().T
    P=np.diag(np.exp(2j*np.pi*phase/6))
    for H in (E,M,G): assert np.linalg.norm(H-H.conj().T)<1e-10
    assert np.linalg.norm(np.linalg.matrix_power(P,6)-np.eye(n))<1e-10
    assert np.linalg.norm(P.conj().T@M@P-M)<1e-10
    return E,M,G,P,{'dimension':n,'flux_group_order':d,'matter_hopping_transitions':hop_count}

def experiment(q,time=.9,cycles=(4,16,64)):
    E,M,G,P,stats=operators(q);H=E+M+G
    n=len(E);frames=[np.linalg.matrix_power(P,k) for k in range(6)]
    twirl=sum((F.conj().T@H@F for F in frames))/6
    mean_res=float(np.linalg.norm(twirl-E-M))
    assert mean_res<1e-10,mean_res
    target=expm(-1j*time*(E+M if q!=6 else H))
    errors=[]
    for ncy in cycles:
        U=expm(-1j*H*time/(6*ncy))
        cyc=np.eye(n, dtype=complex)
        for F in frames:cyc=F.conj().T@U@F@cyc
        errors.append(float(np.linalg.norm(np.linalg.matrix_power(cyc,ncy)-target,'fro')/np.sqrt(n)))
    if q==6:
        assert np.linalg.norm(P-np.eye(n))<1e-12
        assert np.linalg.norm(G)<1e-12
        assert max(errors)<1e-9
    else:
        assert np.linalg.norm(P-np.eye(n))>1e-3
        assert np.linalg.norm(G)>1e-3
        assert errors[-1]<errors[0]
    return {'gamma_order':q,**stats,'allowed_centre_classes_e6':list(range(0,6,q)),
            'number_of_distinct_pulse_phases':6//q,
            'pulse_identity':q==6,'magnetic_loop_fro':float(np.linalg.norm(G)),
            'matter_pulse_commutator_fro':float(np.linalg.norm(P@M-M@P)),
            'averaged_magnetic_cancellation_residual_fro':mean_res,
            'cycles':list(cycles),'finite_echo_errors_per_sqrt_dimension':errors}

def main():
    output={'model':'Gauge-constrained Z_(6/q) central loop and independent Z2 matter-transport sector on four-vertex ring',
            'arms':[experiment(q) for q in ORDERS],
            'limitations':['Not full SU(3)xSU(2)xU(1) or actual quark/lepton dynamics',
                'Gauge constraints explicitly imposed and checked on each transition',
                'No native pulse Hamiltonian, noise, or measurable global-form discrimination',
                'Gamma=Z6 excludes charged Wilson shifts; apparent echo is trivially the identity']}
    Path('sm_center_gauss_constrained_results.json').write_text(json.dumps(output,indent=2)+'\n')
    for row in output['arms']:print('QUOTIENT_NEGATIVE_CONTROL',json.dumps(row,sort_keys=True))
    print('GAUSS_CONSTRAINED_CENTRE_ECHO_CHECKS_PASSED')

if __name__=='__main__':main()
