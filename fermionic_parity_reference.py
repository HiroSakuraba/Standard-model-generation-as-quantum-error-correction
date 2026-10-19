"""Parity-reference resolution for the SO(10)/Reed-Muller 15/16-mode codes.

The punctured 15-mode [[15,1,3]] binary Reed-Muller code has basis
|0_L> (even data parity), |1_L> (odd data parity).
Parity superselection forbids coherent use of their isolated superposition.

Adding a reference occupation mode makes TWO even-total-parity codewords
|0_L,0_ref>, |1_L,1_ref>. This 2D *relational code* has one logical qubit:
15-mode rank-14 checks plus ONE total parity check => 16-15=1.
Logical Z is reference number parity; logical X flips all 16 qubit-mode
occupation bits, commutes with total parity, and is extremely nonlocal.

Adding the unpunctured RM(1,4) extra X stabilizer fixes the + logical
state, leaving 16 independent checks => k=0. Thus its unique state is
NOT one freely encoded qubit. A reference mode plus parity constraint
and the full unpunctured code are DIFFERENT constructions.

Uses existing spinor_reed_muller.kets() as a source of punctured code
basis states, and checks logical operations independently on 2^16 vectors.
All of this is a formal fermionic-mode register, not the dynamics of the SM.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from spinor_reed_muller import kets,rank2,rm,span2,spinor_bits

MODES=16
DATA_MODES=15
DIM=2**MODES

def construct_even_code():
    n,checks,z0,z1=kets()
    assert n==DATA_MODES and len(z0)==2**15
    zero=np.zeros_like(z0)
    v0=np.concatenate((z0,zero))
    v1=np.concatenate((zero,z1))
    K=np.stack((v0,v1),axis=1)
    assert K.shape==(DIM,2)
    assert np.linalg.norm(K.conj().T@K-np.eye(2))<1e-11
    return checks,z0,z1,K


def bit_parity(n):
    """0/1 parity of all mode occupations (MSB = reference mode)."""
    out=np.zeros(n,dtype=np.int8)
    for j in range(MODES):
        out^=((np.arange(n,dtype=np.uint32)>>j)&1).astype(np.int8)
    return out


def logical_dressed_X(vector):
    """Occupation complement on ALL 16 modes, an even-parity nonlocal gate."""
    return vector[np.arange(DIM,dtype=np.uint32)^(DIM-1)]


def total_parity(vector):
    return vector*np.where(bit_parity(DIM)==0,1.,-1.)


def reference_Z(vector):
    return vector*np.concatenate((np.ones(2**DATA_MODES),-np.ones(2**DATA_MODES)))


def rank_certificate(checks):
    """Binary CSS rank, with physical total parity constraint.

    15-data X checks = 4 independent rows. Z checks = 10 independent
    rows. Total parity on all 16 = an independent 11th Z check.
    Extra X_ALL from unpunctured RM(1,4) removes the relational logical.
    """
    x=[list(row)+[0] for row in checks]
    assert rank2(x)==4
    L=np.asarray(checks,dtype=int)
    z=[]
    # Z checks of punctured [[15,1,3]]: orthogonal to X and EVEN parity.
    for item in span2([[int(i==j) for i in range(DATA_MODES)] for j in range(DATA_MODES)]):
        if not np.any(L@np.array(item)%2) and sum(item)%2==0:
            z.append(list(item)+[0])
    assert rank2(z)==10
    parity=[1]*MODES
    assert rank2(z+[parity])==11
    assert all(sum(int(a)&int(b) for a,b in zip(xx,parity))%2==0 for xx in x)
    all_X=[1]*MODES
    assert rank2(x+[all_X])==5
    assert all(sum(int(a)&int(b) for a,b in zip(zz,all_X))%2==0 for zz in z+[parity])
    return {'modes':16,'punctured_15_X_check_rank':4,
            'punctured_15_Z_check_rank':10,
            'extra_total_parity_check_rank':1,
            'relational_16_mode_total_stabilizer_rank':15,
            'relational_logical_qubits':1,
            'full_unpunctured_16_mode_rank':16,
            'full_unpunctured_logical_qubits':0}


def audit():
    checks,z0,z1,K=construct_even_code()
    p_data=np.array([bin(i).count('1')%2 for i in range(2**DATA_MODES)])
    assert np.linalg.norm(z0[p_data==1])<1e-12
    assert np.linalg.norm(z1[p_data==0])<1e-12
    for j in range(2):
        assert np.linalg.norm(total_parity(K[:,j])-K[:,j])<1e-12
        assert np.linalg.norm(logical_dressed_X(K[:,j])-K[:,1-j])<1e-12
        assert np.linalg.norm(reference_Z(K[:,j]) - ((-1)**j)*K[:,j])<1e-12
    alpha=np.sqrt(.3);beta=np.sqrt(.7)
    psi=alpha*K[:,0]+beta*K[:,1]
    assert abs(np.vdot(psi,psi)-1)<1e-11
    assert np.linalg.norm(total_parity(psi)-psi)<1e-11
    Z=float(np.vdot(psi,reference_Z(psi)).real)
    X=float(np.vdot(psi,logical_dressed_X(psi)).real)
    assert abs(Z-(alpha**2-beta**2))<1e-11
    assert abs(X-2*alpha*beta)<1e-11
    # The isolated 15-mode parity SSR twirl removes off-diagonal terms
    # between |0_L> and |1_L>; its expectation of parity-odd X is zero.
    naive=alpha*z0+beta*z1
    parity_even=np.where(p_data==0,naive,0.)
    parity_odd=np.where(p_data==1,naive,0.)
    naive_twirl_X=(np.vdot(parity_even,parity_even[::-1])+
                    np.vdot(parity_odd,parity_odd[::-1]))
    assert abs(naive_twirl_X)<1e-12
    # The *full unpunctured* extra X_ALL projects to the unique + state.
    plus=(K[:,0]+K[:,1])/np.sqrt(2)
    minus=(K[:,0]-K[:,1])/np.sqrt(2)
    assert np.linalg.norm(logical_dressed_X(plus)-plus)<1e-11
    assert np.linalg.norm(logical_dressed_X(minus)+minus)<1e-11
    # Both relational logical basis states have total occupations 0,8,16.
    weights=np.array([int(x).bit_count() for x in range(DIM)])
    assert all(set(weights[abs(K[:,j])>1e-10])<={0,8,16} for j in range(2))
    phases=np.exp(1j*np.pi*weights/4)
    assert np.linalg.norm(phases[:,None]*K-K)<1e-11
    report={'citation':'Mudassar, Schuckert and Gottesman, PRX Quantum 7, 020315 (2026), DOI 10.1103/bglx-ssjg; comparison only',
            **rank_certificate(checks),
            'relational_logical_Z_expectation':Z,
            'relational_logical_X_expectation':X,
            'isolated_15_mode_parity_twirled_X_expectation':float(np.real(naive_twirl_X)),
            'parity_fixed_logical_X_is_nonlocal':True,
            'full_unpunctured_unique_state':'(|0L,0ref>+|1L,1ref>)/sqrt(2)',
            'major_limitation':'One 16-mode reference parity subspace is a formal logical qubit; the full unpunctured RM stabilizer state is k=0. No gate-locality or physical SM dynamics established.'}
    return report


def main():
    data=audit()
    Path('fermionic_parity_reference_results.json').write_text(json.dumps(data,indent=2)+'\n')
    print('FERMIONIC_PARITY_REFERENCE',json.dumps(data,sort_keys=True))
    print('EVEN_PARITY_REFERENCE_AND_K0_STATE_DISTINGUISHED')


if __name__=='__main__':main()
