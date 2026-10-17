#!/usr/bin/env python3
"""Exact finite-field checks of three SO(10) generations and a flavor-mixing obstruction.

Pure Python 3.9+. Formal code coordinates represent fermion species, not
physical particle-register Hilbert spaces. No external libraries.

Run: python3 three_generation_flavor.py
"""
from __future__ import annotations
from collections import Counter
from fractions import Fraction
from itertools import combinations, product
from matter_qec_study import (weights, Q0, rank_cols, triorthogonal_integer,
    columns, z_logical_supports, rows_full, dot)

# The first two family labels reproduce the original (+,-) construction;
# the third supplies an independent charge-label direction.
FAMILIES = ((1, 0), (-1, 0), (0, 1))
B3 = ((1,0,0,0,0),(0,1,0,0,1),(0,0,1,1,0))
B5 = ((1,0,0,0,0),(0,1,1,1,1))


def charge_rows(B, labels=FAMILIES):
    q = weights()
    common = [[w[j] for lab in labels for w in q] for j in range(5)]
    differential = [[lab[k]*dot(b,w) for lab in labels for w in q]
                    for k in range(2) for b in B]
    return common+differential


def puncture(rows, n_families=3):
    cs=columns(rows)
    removed_indices=set(16*g for g in range(n_families))
    removed=[cs[i] for i in sorted(removed_indices)]
    physical=[cs[i] for i in range(len(cs)) if i not in removed_indices]
    return removed, physical


def counting(rows,p,n_families=3, max_weight=4):
    removed,physical=puncture(rows,n_families)
    r=rank_cols(columns(rows),p)
    k=rank_cols(removed,p)
    assert r==len(rows), 'chosen rows not independent'
    assert rank_cols(physical,p)==r, 'puncture loses row rank'
    assert k==n_families, 'punctured singlets fail to encode independently'
    good,bad=triorthogonal_integer(rows)
    assert good, ('cubic anomaly obstruction', bad)
    result={'prime':p,'n':len(physical),'k':k,'charge_rank':r,
            'X_stabilizers':r-k,'Z_stabilizers':len(physical)-r}
    for w in range(1,max_weight+1):
        supp=z_logical_supports(removed,physical,p,w)
        partitions=Counter(tuple(sum(i//15==g for i in c) for g in range(n_families))
                           for c in supp)
        result['weight_'+str(w)+'_logical_supersets']=len(supp)
        result['weight_'+str(w)+'_partition']={str(a):v for a,v in sorted(partitions.items())}
        if w==3:
            assert len(supp)==25*n_families
            assert all(max(parts)==3 for parts in partitions)
        if w==4:
            for g,h in combinations(range(n_families),2):
                kk=tuple(2 if i in (g,h) else 0 for i in range(n_families))
                assert partitions[kk] in (88,396)
                result[f'irreducible_mixed_2plus2_{g}{h}']=partitions[kk]
    assert result['weight_1_logical_supersets']==0
    assert result['weight_2_logical_supersets']==0
    return result


def x_distance(rows,p):
    """Exact optimized X-distance for the three-family rows.

    Each row combination is (a.q + s1*beta.Bq + s2*gamma.Bq).
    Minimize pair (first two families, opposite beta) independently of third
    family (gamma) for fixed common vector a. This searches *every* row
    combination without p^(5+2m) iteration.
    """
    q=weights()
    B=(B3 if p==3 else B5)
    # Integer row generators reduce modulo p; no need to enumerate duplicate
    # coefficients since these B have independent rows modulo the chosen p.
    images=[]
    for coeff in product(range(p), repeat=len(B)):
        w=tuple(sum(coeff[j]*B[j][i] for j in range(len(B)))%p for i in range(5))
        images.append(w)
    assert len(set(images))==len(images)
    def evaluate(v):
        vals=[sum(v[j]*w[j] for j in range(5))%p for w in q]
        return sum(z!=0 for z in vals[1:]), vals[0]!=0
    stats={a:evaluate(a) for a in product(range(p),repeat=5)}
    ans=100000
    witness=None
    for a in product(range(p),repeat=5):
        w3base, l3base=stats[a]
        best_pair_any=(100000,None)
        best_pair_log=(100000,None)
        best_third_any=(100000,None)
        best_third_log=(100000,None)
        for b in images:
            plus=tuple((a[i]+b[i])%p for i in range(5))
            minus=tuple((a[i]-b[i])%p for i in range(5))
            wp,lp=stats[plus];wm,lm=stats[minus]
            if wp+wm<best_pair_any[0]:best_pair_any=(wp+wm,b)
            if (lp or lm) and wp+wm<best_pair_log[0]:best_pair_log=(wp+wm,b)
            if wp<best_third_any[0]:best_third_any=(wp,b)
            if lp and wp<best_third_log[0]:best_third_log=(wp,b)
        # Any puncture nonzero in either the first pair or third generation.
        for pair,third in ((best_pair_any,best_third_log),(best_pair_log,best_third_any)):
            if pair[1] is not None and third[1] is not None and pair[0]+third[0]<ans:
                ans=pair[0]+third[0]
                witness=(a,pair[1],third[1])
    assert witness is not None
    return ans,witness


def cabibbo_check():
    """Exact rational two-flavor rotation check; not a physical rate.

    An orthogonal rotation with cos=4/5, sin=3/5 is a clean toy example.
    For D=diag(d1,d2), the offdiagonal of V.T D V has magnitude
    abs(d1-d2)*abs(cos*sin). It vanishes for nonzero mixing only if d1=d2.
    """
    c,s=Fraction(4,5),Fraction(3,5)
    V=((c,s),(-s,c))
    d1,d2=Fraction(1),Fraction(-1)
    dprime=[[sum(V[k][i]*[d1,d2][k]*V[k][j] for k in range(2))
             for j in range(2)] for i in range(2)]
    assert dprime[0][1] == (d1-d2)*c*s
    assert dprime[1][0] == (d1-d2)*c*s
    assert dprime[0][1] == Fraction(24,25)
    # Algebraic observation: for sin(theta)*cos(theta) != 0 and d1 != d2,
    # a flavor diagonal generator is not diagonal after rotation.
    return {'illustrative_rotation': 'cos=4/5, sin=3/5',
            'charge_eigenvalues':(int(d1),int(d2)),
            'offdiagonal_magnitude':str(abs(dprime[0][1])),
            'meaning':'broken exact generation charge; not a qudit error rate'}


def main():
    print('FLAVOR MIXING:',cabibbo_check())
    for p,B in ((3,B3),(5,B5)):
        rows=charge_rows(B)
        data=counting(rows,p)
        dx,witness=x_distance(rows,p)
        data['X_distance']=dx
        data['X_witness']={'a':witness[0], 'beta_Bq':witness[1], 'gamma_Bq':witness[2]}
        assert dx>=3
        print('THREE GENERATIONS:',data)
    print('ALL THREE-GENERATION CHECKS PASSED')

if __name__=='__main__':main()
