"""Exact finite SU(3) character-basis gauge-invariant error-algebra audit.

Single-loop class-function sector, NOT a lattice glueball spectrum. Uses
Peter-Weyl characters and a hard cutoff in SU(3) highest weights p+q<=L.
Fundamental multiplication via the SU(3) Pieri rule; adjoint multiplication
from full-character fusion (including intermediate reps beyond cutoff).

Requires numpy and sympy. Run: python su3_single_plaquette_error_algebra.py
"""
from __future__ import annotations

import argparse
import json
import numpy as np
import sympy as sp

# Use an intentionally interpretable cutoff-two ordering.
LOW_REPS = ((0, 0), (1, 0), (0, 1), (1, 1), (2, 0), (0, 2))


def reps_at_cutoff(level):
    reps = [(p, q) for k in range(level + 1) for p in range(k, -1, -1) for q in (k-p,)]
    if level == 2:
        return LOW_REPS
    return tuple(reps)


def fundamental_fusion(label):
    p,q=label
    out=[(p+1,q)]
    if p >= 1: out.append((p-1,q+1))
    if q >= 1: out.append((p,q-1))
    return tuple(out)


def casimir(label):
    p,q=label
    return sp.Rational(p*p+q*q+p*q+3*p+3*q,3)


def matrices(reps):
    index={x:i for i,x in enumerate(reps)}
    n=len(reps)
    w=sp.zeros(n)
    for j,rep in enumerate(reps):
        for target in fundamental_fusion(rep):
            if target in index:
                w[index[target],j]+=1
    d=sp.diag(*(casimir(rep) for rep in reps))
    # chi_8 = chi_3 chi_bar3 - chi_1. Perform fusion BEFORE truncation.
    # Intermediate reps are the union of bar3*touched characters, thus
    # may lie outside the retained space (boundary-return paths matter).
    a8=sp.zeros(n)
    for j,rep in enumerate(reps):
        p,q=rep
        intermediates=[(p,q+1)]
        if q>=1: intermediates.append((p+1,q-1))
        if p>=1: intermediates.append((p-1,q))
        for mid in intermediates:
            for dest in fundamental_fusion(mid):
                if dest in index:
                    a8[index[dest],j]+=1
        a8[j,j]-=1
    assert a8 == a8.T
    assert (w+w.T)==(w+w.T).T
    if (1,1) in index:
        assert a8[index[(1,1)],index[(0,0)]]==1
    return d,w,a8


def commutant_dimension(ms):
    n=ms[0].rows
    eye=sp.eye(n)
    constraints=sp.Matrix.vstack(*(sp.kronecker_product(eye,m)-sp.kronecker_product(m.T,eye) for m in ms))
    return n*n - sp.polys.matrices.DomainMatrix.from_Matrix(constraints).convert_to(sp.QQ).rank()


def algebra_basis(gens):
    n=gens[0].rows
    mats=[sp.eye(n)]
    vectors=[sp.Matrix(n*n,1,list(mats[0]))]
    while True:
        old=len(mats)
        for x in tuple(mats):
            for g in gens:
                new=x*g
                col=sp.Matrix(n*n,1,list(new))
                if sp.Matrix.hstack(*vectors,col).rank()>len(vectors):
                    mats.append(new)
                    vectors.append(col)
        if len(mats)==old:
            return mats


def scalar_on_code(mat, indexes):
    sub=mat.extract(indexes,indexes)
    return sub == sub[0,0]*sp.eye(len(indexes))


def low_level_exact():
    reps=LOW_REPS
    d,w,a8=matrices(reps)
    assert d==sp.diag(0,sp.Rational(4,3),sp.Rational(4,3),3,sp.Rational(10,3),sp.Rational(10,3))
    assert all((p+2*q)%3 == (r+2*s)%3 for j,(p,q) in enumerate(reps) for i,(r,s) in enumerate(reps) if a8[i,j])
    assert a8[0,3]==1 and a8[3,3]==2  # 8 x 8 contains 1 and two 8s
    assert a8[1,1]==a8[2,2]==1  # 8 x 3 contains 3; conjugate pair agrees
    assert w[2,1]==1, 'fundamental loop flips 3 to anti-3'
    # Distinct triality sectors and conjugation-paired copies of the same
    # 2x2 'gauge' module under adjoint-only operators.
    t0=[0,3]
    t1=[1,5]
    t2=[2,4]
    assert d.extract(t1,t1)==d.extract(t2,t2)
    assert a8.extract(t1,t1)==a8.extract(t2,t2)
    assert a8.extract(t1,t2)==sp.zeros(2)
    assert a8.extract(t2,t1)==sp.zeros(2)
    assert a8.extract(t0,t1)==sp.zeros(2)
    dims={
        'Casimir_only':commutant_dimension([d]),
        'Casimir_and_adjoint_loop':commutant_dimension([d,a8]),
        'Casimir_and_C_even_fundamental_loop':commutant_dimension([d,w+w.T]),
        'Casimir_and_fundamental_and_conjugate_loops':commutant_dimension([d,w,w.T]),
    }
    alg_adj=algebra_basis([d,a8])
    alg_real=algebra_basis([d,w+w.T])
    alg_full=algebra_basis([d,w,w.T])
    assert len(alg_adj)==8, len(alg_adj)
    assert len(alg_real)==20, len(alg_real)
    assert len(alg_full)==36, len(alg_full)
    assert dims=={'Casimir_only':10,'Casimir_and_adjoint_loop':5,'Casimir_and_C_even_fundamental_loop':2,'Casimir_and_fundamental_and_conjugate_loops':1},dims
    # All products B^* C in the full 8D noise algebra compress to scalar,
    # proving full-algebra KL for encoded |3>,|bar3> under D,A8 errors.
    code=(1,2)
    for a in alg_adj:
        for b in alg_adj:
            assert scalar_on_code(a.T*b,code)
    assert not scalar_on_code(w,code)
    # Since I is in the Kraus linear span, w being nonscalar already
    # violates KL for an enlarged error set {I,w}.
    assert w.extract(code,code)==sp.Matrix([[0,0],[1,0]])
    assert (w+w.T).extract(code,code)==sp.Matrix([[0,1],[1,0]])
    assert d.extract((0,3),(0,3))==sp.diag(0,3)
    summary={
        'reps':[list(r) for r in reps],
        'casimir':[str(d[i,i]) for i in range(6)],
        'adjoint_algebra_dim':len(alg_adj),
        'real_fundamental_algebra_dim':len(alg_real),
        'full_algebra_dim':len(alg_full),
        'commutant_dims':dims,
        'conditional_code':'|3> and |anti-3> exactly satisfy KL against entire alg(C2, chi8)',
        'fundamental_breaks':'P chi3 P = [[0,0],[1,0]]',
    }
    return summary


def spectrum_truncation(level, mag=1.0):
    reps=reps_at_cutoff(level)
    d,w,a8=matrices(reps)
    h=np.array(d,float)-mag*np.array(w+w.T,float)/2
    vals, vecs = np.linalg.eigh(h)
    result={'cutoff':level,'dim':len(reps),'lowest_energies':[round(float(v),9) for v in vals[:min(4,len(vals))]],'first_gap':round(float(vals[1]-vals[0]),9)}
    if level == 2:
        q=vecs[:, :2]
        pd=q.T @ np.array(d,float) @ q
        scalar_error=np.linalg.norm(pd-np.trace(pd)/2*np.eye(2),ord=2)
        result['low_two_eigenstates_Casimir_KL_traceless_norm']=round(float(scalar_error),9)
        assert scalar_error > 0.8
    return result


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--max-cutoff',type=int,default=4)
    args=parser.parse_args()
    s=low_level_exact()
    print('Six-character exact SU(3) Casimir, fusion and error-algebra tests: PASS')
    print('Adjoint-only algebra dimension 8, commutant 5; KL-encoded qubit PASS')
    print('Fundamental gauge-invariant Wilson loop grows algebra to M6; KL broken: PASS')
    # A small cutoff check: energy gaps are NOT glueball masses.
    for level in range(1,args.max_cutoff+1):
        print('CUTOFF:',spectrum_truncation(level))
    # Independent full-matrix algebras persist at the next two cutoffs.
    for level in (3,4):
        r=reps_at_cutoff(level)
        d,w,_=matrices(r)
        dim=commutant_dimension([d,w,w.T])
        assert dim == 1, (level,dim)
        print(f'CUTOFF {level}: full fundamental/electric algebra M_{len(r)} verified by exact commutant')
    print('JSON_RESULTS='+json.dumps(s,sort_keys=True))


if __name__=='__main__':
    main()
