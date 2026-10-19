"""SU(3) published qutrit square chain: four/five-plaquette spectra and KL errors.

Source normalization: Chen, Mueller, Yao, Phys. Rev. D 114 074502 (2026),
Eqs. (3),(5),(10),(11), Appendix B Table VI. Open, zero-flux boundary.
Exact SymPy validation / Knill-Laflamme; NumPy floating spectra.
NOT a glueball mass, continuum limit, or physical recovery circuit.
"""
from itertools import product
import json
import numpy as np
import sympy as sp

from su3_three_plaquette_chain import (
    PUBLISHED, VERTEX, verify_source_tables, operators as prior_operators,
    charge_conjugation as prior_conjugation, reflection as prior_reflection,
)

CAS = sp.Rational(4,3)
IDENTITY = sp.Integer(1)

def basis_states(n):
    assert n >= 1
    return tuple(product(range(3), repeat=n))

def link_fluxes(q):
    # Top n, bottom n, vertical n+1. Values modulo 3: 0=1, 1=3, 2=bar3.
    return tuple(q)+tuple((-x)%3 for x in q)+((-q[0])%3,) + tuple((a-b)%3 for a,b in zip(q,q[1:])) + (q[-1],)

def electric_links(q):
    return tuple(CAS if x else sp.Integer(0) for x in link_fluxes(q))

def one_site_actions(q,site,kind="B"):
    """Return sparse exact output (new face state, matrix element)."""
    left = q[site-1] if site else 0
    right = q[site+1] if site+1 < len(q) else 0
    m = PUBLISHED[(left,right)]
    old = q[site]
    if kind == "W":
        changes = ((old+1)%3,)
    elif kind == "Wdag":
        changes = ((old-1)%3,)
    elif kind == "B":
        changes = ((old+1)%3,(old-1)%3)
    else:
        raise ValueError(kind)
    out = {}
    for new in changes:
        p = q[:site]+(new,)+q[site+1:]
        a = sp.simplify(m[new,old])
        if a:
            out[p] = out.get(p,sp.Integer(0))+a
    return out

def dense_matrices(n):
    labels=basis_states(n)
    pos={q:i for i,q in enumerate(labels)}
    diag=np.array([float(sum(electric_links(q))) for q in labels])
    hparts=[]
    for site in range(n):
        block=np.zeros((len(labels),len(labels)))
        for col,q in enumerate(labels):
            for q2,a in one_site_actions(q,site).items():
                block[pos[q2],col]=float(a)
        assert np.array_equal(block,block.T)
        hparts.append(block)
    return labels, diag, hparts

def conjugate(q):return tuple((-x)%3 for x in q)

def transformed(q):return q[::-1]

def symmetry_matrix(labels,transform):
    ix={q:i for i,q in enumerate(labels)}
    out=np.zeros((len(labels),len(labels)))
    for j,q in enumerate(labels):out[ix[transform(q)],j]=1
    return out

def validate_exact():
    verify_source_tables()
    for n in (2,3,4,5):
        labels=basis_states(n)
        assert len(labels)==3**n
        assert len(set(link_fluxes(q) for q in labels))==len(labels)
        assert all(len(link_fluxes(q))==3*n+1 for q in labels)
        # C conjugation, reflection. All magnetic coefficients exactly real,
        # invariant under global C, reflected left/right by chain reversal.
        for q in labels:
            assert electric_links(conjugate(q)) == electric_links(q)
            for site in range(n):
                x=one_site_actions(q,site)
                cx=one_site_actions(conjugate(q),site)
                rx=one_site_actions(transformed(q),n-1-site)
                assert {conjugate(z):v for z,v in x.items()}==cx
                assert {transformed(z):v for z,v in x.items()}==rx
        if n <= 3:
            p,D,W,B=prior_operators(n)
            assert p==labels
            for q in labels:
                i=p.index(q)
                assert D[i,i]==sum(electric_links(q))
                for site in range(n):
                    for dest,a in one_site_actions(q,site,"W").items():
                        assert W[site][p.index(dest),i]==a
                    for dest,a in one_site_actions(q,site).items():
                        assert B[site][p.index(dest),i]==a
        print(f'EXACT n={n}: {len(labels)} unique singlets, {3*n+1} links; all local Wilson entries and C/reflection actions PASS')


def applied_error(q,kind,site=None):
    if kind == 'I':return {q:IDENTITY}
    if kind == 'E':return {q:electric_links(q)[site]}
    return one_site_actions(q,site,kind)

def dot_sparse(a,b):
    return sp.simplify(sum(sp.conjugate(v)*b.get(q,0) for q,v in a.items()))

def kl_check(codewords,errors):
    # Each error is a tuple (kind,site), with site None for I.
    imgs=[[applied_error(q,*e) for q in codewords] for e in errors]
    # P E_a^dagger E_b P must be scalar as a k x k matrix.
    grams={}
    for i,a in enumerate(imgs):
        for j,b in enumerate(imgs):
            m=sp.Matrix(len(codewords),len(codewords),lambda p,q:dot_sparse(a[p],b[q]))
            if m != m[0,0]*sp.eye(len(codewords)):
                return False,(errors[i],errors[j],m)
            grams[i,j]=m[0,0]
    return True,sp.Matrix(len(errors),len(errors),lambda i,j:grams[i,j])

def kl_tests():
    n=5
    q0=(0,0,0,0,1)
    q1=(0,0,0,0,2)
    code=[q0,q1]
    # The following is spatially restricted *single-event* noise; it is
    # neither uniform error immunity nor correction of the generated algebra.
    electric=[('E',e) for e in range(3*n+1)]
    limited=[('I',None)]+electric+[('B',j) for j in range(n-1)]
    passed,mat=kl_check(code,limited)
    assert passed and mat.shape==(3*n+1+n,)*2
    assert mat==mat.T.conjugate() and mat.rank()==5
    print(f'KL PASS: n=5 separated orientation code, {len(limited)} elementary errors (identity + {3*n+1} local link Casimirs + B0..B3), {len(limited)**2} exact overlaps, Gram rank 5')
    with_all=limited+[('B',4)]
    passed,issue=kl_check(code,with_all)
    assert not passed and issue[2]==sp.Matrix([[0,1],[1,0]]),issue
    print('KL FAIL: adding local B4 acts as logical X, despite exact gauge invariance')
    passed,issue=kl_check(code,[('I',None),('W',3)])
    assert not passed and issue == (('W',3),('W',3),sp.diag(sp.Rational(1,9),sp.Rational(1,3)))
    print('KL FAIL: W3 by itself (one plaquette away), P W3^dag W3 P = diag(1/9,1/3)')
    # Mixed offsite adjacent noise terms can discriminate code via the boundary.
    # Important: no claim that time evolution generated by physical H protects code.
    return {'codewords':[list(q) for q in code],'correctable_error_count':len(limited),
            'KL_matrix_tests':len(limited)**2,'gram_rank':5,'B4_witness':[[0,1],[1,0]],
            'W3_failure':str(issue)}

def parity_sectors(labels,H,n):
    C=symmetry_matrix(labels,conjugate)
    R=symmetry_matrix(labels,transformed)
    assert np.array_equal(C@C,np.eye(len(labels)))
    assert np.array_equal(R@R,np.eye(len(labels)))
    assert np.allclose(H@C,C@H,atol=1e-12)
    assert np.allclose(H@R,R@H,atol=1e-12)
    dims={(c,r):round(np.trace((np.eye(len(labels))+c*C)@(np.eye(len(labels))+r*R))/4)
          for c,r in product((1,-1),repeat=2)}
    assert sum(dims.values())==len(labels)
    assert sum(d for (c,r),d in dims.items() if c==1)==(3**n+1)//2
    assert sum(d for (c,r),d in dims.items() if c==-1)==(3**n-1)//2
    return C,R,dims

def spectral_study():
    results=[]
    for n in (2,3,4,5):
        labels,D,B=dense_matrices(n)
        additive=6*n
        C,R=symmetry_matrix(labels,conjugate),symmetry_matrix(labels,transformed)
        assert np.allclose(C@C,np.eye(3**n)) and np.allclose(R@R,np.eye(3**n))
        for g2 in (0.5,1.,1.5,2.):
            H=g2*np.diag(D)+(additive*np.eye(len(labels))-sum(B))/g2
            cn,rn,dims=parity_sectors(labels,H,n)
            eig,V=np.linalg.eigh(H)
            assert np.max(np.abs(H@V-V*eig)) < 2e-11
            ground=V[:,0]
            c=float(ground@cn@ground)
            r=float(ground@rn@ground)
            assert abs(c-1)<1e-9 and abs(r-1)<1e-9
            excited=V[:,1]
            excited_c=float(excited@cn@excited)
            excited_r=float(excited@rn@excited)
            assert abs(abs(excited_c)-1)<1e-9 and abs(abs(excited_r)-1)<1e-9
            weights=(ground*ground)
            left=np.array([int(q[0]!=0) for q in labels])
            right=np.array([int(q[-1]!=0) for q in labels])
            conn=float(weights@(left*right)-(weights@left)*(weights@right))
            # The ground state is nondegenerate numerically at tested couplings.
            assert eig[1]-eig[0] > 1e-8
            row={'plaquettes':n,'dim':len(labels),'g2':g2,'ground_energy':float(eig[0]),
                 'first_gap':float(eig[1]-eig[0]),'second_gap':float(eig[2]-eig[0]),
                 'ground_C':round(c),'ground_R':round(r),
                 'first_excited_C':round(excited_c),'first_excited_R':round(excited_r),
                 'end_occupation_cov':conn}
            results.append(row)
            print('SPECTRUM',json.dumps(row,sort_keys=True))
        print('SYMMETRY C/R sector dims',n,{str(k):v for k,v in dims.items()})
    return results

def main():
    validate_exact()
    kl=kl_tests()
    spectra=spectral_study()
    print('JSON_RESULTS='+json.dumps({'spectra':spectra,'bounded_error_KL':kl},sort_keys=True))
    print('ALL SU3 N=4/5 PUBLISHED-CHAIN AND BOUNDED-KL CHECKS PASSED')

if __name__=='__main__':main()
