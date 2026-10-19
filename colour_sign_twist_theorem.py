"""Independent 42-dimensional Young-seminormal theorem for nine-colour-qutrit code.

No import from the expensive 3^9 colour_singlet_code implementation.
For SU(3) nine-qutrit colour singlets the S9 multiplicity module is the
Specht module of self-conjugate partition (3,3,3), dimension 42.

On standard Young tableaux T, the adjacent transposition matrices use
Young's orthogonal form. The signed tableau-transpose map J sends |T>
to parity(rowword(T)) |T^transpose>. It is real orthogonal, skew and
anticommutes with EVERY transposition. Thus the +i eigenspace K has
dimension 21 and P_K Swap_ij P_K = 0.

Together with SU(3) invariance and the fact that every two-qutrit
SU(3)-invariant operator lies in span{I,Swap}, this proves exact
Knill-Laflamme (KL) for ALL single-qutrit errors, and two-qutrit erasure.
Proof is valid independent of floating-point tests.

The explicit recovery below corrects the channel {sqrt(1-p)I, sqrt(p)S12}
by measuring the sign-twist sector and applying S12 on the -i sector.
It is NOT a recovery for arbitrary mixtures of unknown pair exchanges:
products Sij Skl have non-scalar KL compressions. A general one-qutrit
recovery exists from KL polar decomposition, but no physical gate-level
measurement is synthesized here.
"""
from __future__ import annotations
from itertools import combinations
import json
from pathlib import Path
import numpy as np

SHAPE=(3,3,3)
N=9
DIM=42


def tableaux():
    """Recursively add 1..9 to Young diagrams without breaking row order."""
    out=[]
    def rec(rows,n):
        if n==N+1:
            out.append(tuple(tuple(r) for r in rows))
            return
        for i in range(3):
            j=len(rows[i])
            if j<3 and (i==0 or len(rows[i-1])>j):
                next_rows=[list(r) for r in rows]
                next_rows[i].append(n)
                rec(next_rows,n+1)
    rec([[],[],[]],1)
    assert len(out)==DIM
    return tuple(out)


def adjacent_generators(ts):
    """S_i |T> = 1/d |T> + sqrt(1-d^-2)|s_i T>, d=c(i+1)-c(i)."""
    index={t:i for i,t in enumerate(ts)}
    generators=[]
    for k in range(1,N):
        M=np.zeros((DIM,DIM))
        for j,t in enumerate(ts):
            pos={v:(row,col) for row,line in enumerate(t) for col,v in enumerate(line)}
            a,b=pos[k],pos[k+1]
            d=(b[1]-b[0])-(a[1]-a[0])
            assert d!=0
            M[j,j]=1/d
            changed=tuple(tuple(k+1 if v==k else k if v==k+1 else v
                                for v in line) for line in t)
            if changed in index:
                M[index[changed],j]=np.sqrt(1-1/d**2)
            else:
                assert abs(d)==1
        assert np.linalg.norm(M@M-np.eye(DIM))<2e-13
        generators.append(M)
    return generators


def signed_tableau_transpose(ts):
    """Explicit sign-twist intertwiner; no numerical null-space search."""
    index={t:i for i,t in enumerate(ts)}
    J=np.zeros((DIM,DIM))
    for j,t in enumerate(ts):
        word=sum((list(r) for r in t),[])
        inversions=sum(word[i]>word[k] for i in range(N) for k in range(i+1,N))
        dual=tuple(tuple(t[k][i] for k in range(3)) for i in range(3))
        J[index[dual],j]=(-1)**inversions
    return J


def swaps(generators):
    """All 36 (ij), using S_(i,j) = s_i...s_(j-1)...s_i."""
    out={}
    for i,j in combinations(range(1,N+1),2):
        M=np.eye(DIM)
        # matrix application order: walk i to j and back.
        for k in list(range(i,j))+list(range(j-2,i-1,-1)):
            M=generators[k-1]@M
        assert np.linalg.norm(M@M-np.eye(DIM))<2e-12
        out[(i,j)]=M
    return out


def code_isometry(J):
    """Explicit +i sector projector, rank=21; numerical ONB for diagnostics."""
    assert np.linalg.norm(J@J+np.eye(DIM))<2e-12
    projector=(np.eye(DIM)-1j*J)/2
    assert np.linalg.norm(projector@projector-projector)<2e-12
    eig,vec=np.linalg.eigh(projector)
    K=vec[:,eig>.5]
    assert K.shape==(DIM,21)
    assert np.linalg.norm(J@K-1j*K)<2e-12
    return K


def exchange_survival(t):
    """H=(S12-I/3)/2. Exact code survival = cos^2(t/2)."""
    return float(np.cos(t/2)**2)


def exchange_channel_recovery(K,J,S,p):
    """Exact isometric flagged recovery for {I, S12}, tested as a quantum channel.

    Measure P_+=(I-iJ)/2 and P_-=(I+iJ)/2. Conditional S12 maps the
    -i sector back to +i. Then decode via K^dag. On the rest of the
    Hilbert space replace by maximally mixed logical state.
    """
    assert 0<=p<=1
    Pp=(np.eye(DIM)-1j*J)/2
    Pm=(np.eye(DIM)+1j*J)/2
    E0=np.sqrt(1-p)*K
    E1=np.sqrt(p)*S@K
    # Kraus recovery R_+=Kdag P_+, R_-=Kdag S P_-.
    R0=K.conj().T@Pp
    R1=K.conj().T@S@Pm
    r0,r1=K.conj().T@K, K.conj().T@S@S@K
    assert np.linalg.norm(r0-np.eye(21))<2e-12
    assert np.linalg.norm(r1-np.eye(21))<2e-12
    coherent=0
    for R in (R0,R1):
        for E in (E0,E1):
            # entanglement fidelity of maximally mixed encoded qudit
            coherent+=abs(np.trace(R@E))**2/(21**2)
    # complement sector not reached by this particular channel; general
    # CPTP extension replaces uncorrected remainder with maximally mixed.
    combined=R0.conj().T@R0+R1.conj().T@R1
    assert np.min(np.linalg.eigvalsh(np.eye(DIM)-combined))>-2e-12
    assert abs(coherent-1)<2e-12
    return {'p':p,'entanglement_fidelity_after_explicit_recovery':float(coherent),
            'recovery_kraus_count_before_complement':2,
            'complement_rank':int(round(np.trace(np.eye(DIM)-combined).real))}


def run():
    ts=tableaux()
    S=adjacent_generators(ts)
    J=signed_tableau_transpose(ts)
    assert np.linalg.norm(J.T@J-np.eye(DIM))==0
    assert np.linalg.norm(J+J.T)==0
    assert np.linalg.norm(J@J+np.eye(DIM))==0
    worst_braid=max(np.linalg.norm(S[i]@S[i+1]@S[i]-S[i+1]@S[i]@S[i+1])
                    for i in range(N-2))
    assert worst_braid<2e-12
    worst_sign=max(np.linalg.norm(J@s+s@J) for s in S)
    assert worst_sign<2e-12
    all_swaps=swaps(S)
    K=code_isometry(J)
    worst_compression=max(np.linalg.norm(K.conj().T@M@K) for M in all_swaps.values())
    assert worst_compression<2e-12
    # Three-cycle is a weight-3 operator and is not scalar on code.
    cycle=all_swaps[(1,2)]@all_swaps[(2,3)]
    reduced=K.conj().T@cycle@K
    nontrivial=np.linalg.norm(reduced-np.trace(reduced)/21*np.eye(21))
    assert nontrivial>1
    # Exact controlled exchange moves the code into orthogonal J sector.
    t=np.pi/2
    U=np.cos(t/2)*np.eye(DIM)-1j*np.sin(t/2)*all_swaps[(1,2)]
    survival=np.linalg.norm(K.conj().T@U@K)**2/21
    assert abs(survival-exchange_survival(t))<1e-12
    # Exact KL for error family {I,S12}; general 1-site KL proof via
    # two-site commutant span{I,Sij} is documented separately.
    rec=exchange_channel_recovery(K,J,all_swaps[(1,2)],0.37)
    report={
       'singlet_specht_shape':SHAPE,
       'singlet_dimension':DIM,'code_dimension':K.shape[1],
       'number_of_transpositions':len(all_swaps),
       'sign_twist_J_squared_minus_identity':True,
       'adjacent_braid_max_fro_residual':float(worst_braid),
       'exchange_anticommutation_max_fro_residual':float(worst_sign),
       'all_exchange_KLK_max_fro_residual':float(worst_compression),
       'three_cycle_non_scalar_fro_norm':float(nontrivial),
       'single_exchange_time_pi_over_2_survival':float(survival),
       'one_pair_channel':rec,
       'claims':{
          'KL_all_single_qutrit_errors':'proved by singlet SU3 Schur lemma plus two-qutrit commutant and PswapP=0',
          'maximality_21':'only among zero-exchange-compression subspaces, NOT globally proved',
          'survival_exchange':'cos^2(t/2); compression of Hamiltonian is not dynamical immunity',
          'recovery':'explicit {I,S12} recovery; general 1-local KL polar recovery exists but not gate-compiled'}
    }
    return report


def main():
    r=run()
    Path('colour_sign_twist_theorem_results.json').write_text(json.dumps(r,indent=2)+'\n')
    print('INDEPENDENT_SIGN_TWIST_VERIFIED',json.dumps(r,sort_keys=True))


if __name__=='__main__':
    main()
