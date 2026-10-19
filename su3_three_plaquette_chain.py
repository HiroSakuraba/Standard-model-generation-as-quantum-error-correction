"""Exact published minimally truncated SU(3) square-chain Hamiltonian.

Chen, Mueller & Yao, Physical Review D 114, 074502 (2026),
DOI 10.1103/1qt7-3ywl, Eqs. (3),(5),(10),(11),
Appendix B Table VI (journal); independent vertex Table I cross-check.

Zero external flux, 1 + 3 + anti-3 irreps on links, open chain.
For N=3, the gauge-invariant face-qutrit Hilbert space has 3**3=27
states. This is not a continuum glueball spectrum or new QEC code.
"""
from itertools import product
import json
import sympy as sp
import numpy as np

SQ = sp.sqrt(3) / 3
TH = sp.Rational(1,3)
BASIC = sp.ones(3) - sp.eye(3)

# Table I: independent *vertex* recoupling factors (ordered neighbors a,b).
# Source: published PRD PDF p. 4, Table I.
VERTEX = {
 (0,0): BASIC,
 (0,1): sp.Matrix([[0,1,1],[TH,0,SQ],[SQ,1,0]]),
 (0,2): sp.Matrix([[0,1,1],[SQ,0,1],[TH,SQ,0]]),
 (1,0): sp.Matrix([[0,TH,SQ],[1,0,1],[1,SQ,0]]),
 (1,1): BASIC,
 (1,2): sp.Matrix([[0,SQ,1],[1,0,1],[SQ,TH,0]]),
 (2,0): sp.Matrix([[0,SQ,TH],[1,0,SQ],[1,1,0]]),
 (2,1): sp.Matrix([[0,1,SQ],[SQ,0,TH],[1,1,0]]),
 (2,2): BASIC,
}

# Independently transcribed *whole plaquette* matrices from published
# PRD Appendix B Table VI, rows and columns labeled neighbors 0,1,2.
PUBLISHED = {
 (0,0): BASIC,
 (0,1): sp.Matrix([[0,TH,SQ],[TH,0,SQ],[SQ,SQ,0]]),
 (0,2): sp.Matrix([[0,SQ,TH],[SQ,0,SQ],[TH,SQ,0]]),
 (1,1): sp.Matrix([[0,sp.Rational(1,9),TH],[sp.Rational(1,9),0,TH],[TH,TH,0]]),
 (1,2): sp.Matrix([[0,SQ/3,SQ/3],[SQ/3,0,TH],[SQ/3,TH,0]]),
 (2,2): sp.Matrix([[0,TH,sp.Rational(1,9)],[TH,0,TH],[sp.Rational(1,9),TH,0]]),
}
PUBLISHED.update({(b,a):m for (a,b),m in list(PUBLISHED.items()) if a!=b})


def verify_source_tables():
    assert len(VERTEX)==len(PUBLISHED)==9
    for l,r in product(range(3), repeat=2):
        # Eq. (8) around 4 CCW corners with PB=PT=0.
        predicted=VERTEX[(0,r)].multiply_elementwise(VERTEX[(r,0)])
        predicted=predicted.multiply_elementwise(VERTEX[(0,l)])
        predicted=predicted.multiply_elementwise(VERTEX[(l,0)])
        assert predicted == PUBLISHED[(l,r)], (l,r,predicted)
        assert predicted == predicted.T and all(predicted[i,i]==0 for i in range(3))
    assert PUBLISHED[(1,1)][0,1]==sp.Rational(1,9)
    assert PUBLISHED[(1,2)][0,1]==sp.sqrt(3)/9
    assert PUBLISHED[(2,2)][0,2]==sp.Rational(1,9)
    print('SOURCE AUDIT: all 9 published Table VI 3x3 blocks (81 entries) reproduced by four exact Table I vertex factors')


def states(n):
    return tuple(product(range(3),repeat=n))


def flux_reps(q):
    """Orient top/bottom and vertical links; 2n horizontal+n+1 vertical."""
    top=tuple(q)
    bottom=tuple((-v)%3 for v in q)
    vertical=((-q[0])%3,)+tuple((q[i]-q[i+1])%3 for i in range(len(q)-1))+(q[-1],)
    return top+bottom+vertical


def electric_casimir(q):
    links=flux_reps(q)
    assert len(links)==3*len(q)+1
    return sp.Rational(4,3)*sum(v!=0 for v in links)


def operators(n):
    basis=states(n); ix={v:i for i,v in enumerate(basis)}
    dim=len(basis)
    D=sp.diag(*(electric_casimir(q) for q in basis))
    W=[]
    for site in range(n):
        w=sp.zeros(dim)
        for i,q in enumerate(basis):
            left=q[site-1] if site else 0
            right=q[site+1] if site<n-1 else 0
            dest=list(q);dest[site]=(q[site]+1)%3
            w[ix[tuple(dest)],i]=PUBLISHED[(left,right)][(q[site]+1)%3,q[site]]
        W.append(w)
    B=[w+w.T for w in W]
    assert all(x==x.T for x in B)
    if n >= 2:
        assert all(w*w.T != sp.eye(dim) for w in W)  # nontrivial recoupling makes W nonunitary
    return basis,D,W,B


def charge_conjugation(basis):
    ix={q:i for i,q in enumerate(basis)};C=sp.zeros(len(basis))
    for j,q in enumerate(basis):C[ix[tuple((-v)%3 for v in q)],j]=1
    return C


def reflection(basis):
    ix={q:i for i,q in enumerate(basis)};R=sp.zeros(len(basis))
    for j,q in enumerate(basis):R[ix[q[::-1]],j]=1
    return R


def exact_validation():
    verify_source_tables()
    for n in (1,2,3):
        basis,D,W,B=operators(n)
        assert len(basis)==3**n
        # Each active face excitation has uniquely specified link irrep.
        link_patterns=[flux_reps(q) for q in basis]
        assert len(set(link_patterns))==len(basis)
        # Discrete charge-conservation at all trivalent (possibly trivial) vertices
        # follows from each edge label being a difference of neighboring faces.
        assert all(len(z)==3*n+1 and set(z)<=set(range(3)) for z in link_patterns)
        C=charge_conjugation(basis)
        assert C*C==sp.eye(3**n)
        assert all(C*a==a*C for a in [D]+B)
        if n==3:
            R=reflection(basis)
            assert R*R==sp.eye(27) and R*C==C*R
            assert R*D==D*R and R*B[0]==B[2]*R and R*B[1]==B[1]*R
            # 27 states resolve exactly into C-even (14) and C-odd (13).
            assert (sp.eye(27)+C).rank()==14 and (sp.eye(27)-C).rank()==13
            sectors={(c,r):sp.trace((sp.eye(27)+c*C)*(sp.eye(27)+r*R))/4
                     for c,r in product((1,-1),repeat=2)}
            assert sum(sectors.values())==27
            print('THREE-PLAQUETTE C/REFLECTION SYMMETRY SECTORS', sectors)
        print(f'States n={n}: {len(basis)}; links={3*n+1}; directed Wilson transitions='+str([sum(x!=0 for x in w) for w in W]))

    # Reproduce independently benchmarked N=2 zero-external-flux operators.
    from su3_published_square_chain_benchmark import (
        PLAQUETTE, NEW_LABELS, published_electric, published_matrix)
    basis2,D2,W2,B2=operators(2)
    perm=[basis2.index(PLAQUETTE[name]) for name in NEW_LABELS]
    assert D2.extract(perm,perm)==published_electric()
    for face,expected in [(0,'L'),(1,'R')]:
        assert W2[face].extract(perm,perm)==published_matrix(expected,oriented=True)
        assert B2[face].extract(perm,perm)==published_matrix(expected)
    print('N=2 PARENT BENCHMARK: 81/81 magnetic matrix entries for each directed loop and every electric energy PASS')

    basis3,D,W,B=operators(3)
    # Freeze rightmost q2=0. Project every left/middle operator back to N=2.
    embedding=[basis3.index(q+(0,)) for q in basis2]
    assert D.extract(embedding,embedding)==D2
    assert all(W[i].extract(embedding,embedding)==W2[i] for i in (0,1))
    assert W[2].extract(embedding,embedding)==sp.zeros(9)
    # Interior i=1 sees all pairs (left,right); all 9 Table VI blocks are used.
    assert set((q[0],q[2]) for q in basis3)==set(product(range(3),repeat=2))
    for left,right in product(range(3),repeat=2):
        sites=[basis3.index((left,m,right)) for m in range(3)]
        assert B[1].extract(sites,sites)==PUBLISHED[(left,right)]
    print('MIDDLE FACE: all 9 conditioned published Table VI blocks embedded exactly in 27-state Hamiltonian')

    # Encoded pair e.g. (0,+1,0),(0,-1,0) is electric-indistinguishable
    # even for ALL individual link Casimirs, but magnetic Bmiddle is X.
    pair=[basis3.index((0,1,0)),basis3.index((0,2,0))]
    assert D.extract(pair,pair)==sp.Rational(16,3)*sp.eye(2)
    links=[flux_reps(basis3[i]) for i in pair]
    assert all((links[0][j]!=0)==(links[1][j]!=0) for j in range(10))
    assert B[1].extract(pair,pair)==sp.Matrix([[0,1],[1,0]])
    # A single magnetic Wilson is admissible and violates KL with identity.
    # This is a candidate-specific obstruction, not active-QEC no-go.
    print('QEC: local electric-only pair KL PASS; middle fundamental Wilson compresses to logical X -> KL FAIL')
    return basis3,D,W,B,C,R


def spectral_audit(basis,D,B,C,R):
    # H = g² D +(6N I - sum B_i)/g² (g² is physical parameter).
    n=3;m=len(basis);N=np.array(D,dtype=float)
    S=np.array(sum(B,sp.zeros(m)),dtype=float)
    Cn=np.array(C,dtype=float);Rn=np.array(R,dtype=float)
    occupation=np.array([[int(v!=0) for v in q] for q in basis],dtype=float)
    _,D2,_,B2=operators(2)
    N2=np.array(D2,dtype=float)
    S2=np.array(sum(B2,sp.zeros(9)),dtype=float)
    summary=[]
    for g2 in (0.5,1.0,1.5,2.0):
        H=g2*N+(6*n*np.eye(m)-S)/g2
        ev,U=np.linalg.eigh(H)
        assert np.linalg.norm(H@U-U*ev,ord=np.inf)<2e-12
        psi=U[:,0]
        c=float(psi@Cn@psi);r=float(psi@Rn@psi)
        assert abs(abs(c)-1)<1e-10 and abs(abs(r)-1)<1e-10
        H2=g2*N2+(12*np.eye(9)-S2)/g2
        ev2=np.linalg.eigvalsh(H2)
        gap2=float(ev2[1]-ev2[0])
        occ=occupation.T@(psi**2)
        occ_lr=float(np.dot((psi**2),occupation[:,0]*occupation[:,2])-occ[0]*occ[2])
        ent=np.linalg.svd(psi.reshape((3,9)),compute_uv=False)
        vnent=float(-sum(p*np.log(p) for p in ent**2 if p>1e-15))
        row={'g2':g2,'E0':float(ev[0]),'gap':float(ev[1]-ev[0]),
             'gap_two_plaquettes':gap2,'gap_three_minus_two':float(ev[1]-ev[0])-gap2,
             'C_ground':int(round(c)),'R_ground':int(round(r)),
             'occ_left':float(occ[0]),'occ_center':float(occ[1]),
             'outer_occupation_connected':occ_lr,'left_cut_entropy_nats':vnent}
        summary.append(row)
        print('SPECTRUM', json.dumps(row,sort_keys=True))
    return summary


def main():
    basis,D,W,B,C,R=exact_validation()
    rows=spectral_audit(basis,D,B,C,R)
    print('JSON_RESULTS='+json.dumps({'dimension':27,'rows':rows},sort_keys=True))
    print('ALL EXACT PUBLISHED THREE-PLAQUETTE SU3 CHECKS PASSED')


if __name__=='__main__':
    main()
