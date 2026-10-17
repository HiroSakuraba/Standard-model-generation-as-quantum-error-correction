#!/usr/bin/env python3
"""Two Standard Model generations as triorthogonal quantum CSS codes.

A standalone, standard-library-only exact arithmetic check of constructions
in the accompanying research note.  Python 3.9+.

The result is a mathematical code constructed from chosen charge lattices.
It does not assert that Standard Model fermions physically store qudits.
"""
from __future__ import annotations
from itertools import combinations, product
from collections import Counter

Q0 = (1, 1, 1, 1, 1)


def weights():
    """Return punctured SO(10) spinor weights in SU(5) 10 + 5bar order."""
    ten = []
    for i,j in combinations(range(5),2):
        q = list(Q0)
        q[i]=q[j]=-1
        ten.append(tuple(q))
    fivebar = []
    for i in range(5):
        q = [-1]*5
        q[i]=1
        fivebar.append(tuple(q))
    return [Q0]+ten+fivebar


def rank_cols(cols,p):
    if not cols:return 0
    a=[[col[r]%p for col in cols] for r in range(len(cols[0]))]
    n=len(a); m=len(cols); r=0
    for j in range(m):
        pivot=next((k for k in range(r,n) if a[k][j]),None)
        if pivot is None:continue
        a[r],a[pivot]=a[pivot],a[r]
        inv=pow(a[r][j],-1,p)
        a[r]=[(x*inv)%p for x in a[r]]
        for k in range(n):
            if k != r and a[k][j]:
                v=a[k][j]
                a[k]=[(x-v*y)%p for x,y in zip(a[k],a[r])]
        r+=1
        if r==n:break
    return r


def dot(a,b):return sum(x*y for x,y in zip(a,b))


def rows_full(B=(),two_generations=True,flavor=False):
    q=weights()
    common=[[u[i] for u in q] for i in range(5)]
    if not two_generations:return common
    common=[r+r for r in common]
    diff=[[dot(b,u) for u in q]+[-dot(b,u) for u in q] for b in B]
    if flavor:diff.append([1]*16+[-1]*16)
    return common+diff


def triorthogonal_integer(rows):
    for i in range(len(rows)):
        for j in range(i,len(rows)):
            for k in range(j,len(rows)):
                if sum(x*y*z for x,y,z in zip(rows[i],rows[j],rows[k])) !=0:
                    return False,(i,j,k)
    return True,None


def columns(rows):return [tuple(row[j] for row in rows) for j in range(len(rows[0]))]


def punctured(rows,p,two_generations=True):
    allcols=columns(rows)
    omitted=(0,16) if two_generations else (0,)
    removed=[allcols[i] for i in omitted]
    phys=[allcols[j] for j in range(len(allcols)) if j not in omitted]
    r=rank_cols(allcols,p)
    k=rank_cols(removed,p)
    assert rank_cols(phys,p)==r,'rank loss after puncture'
    return removed,phys,r,k


def z_logical_supports(removed,phys,p,w=3):
    pr=rank_cols(removed,p)
    res=[]
    for J in combinations(range(len(phys)),w):
        cs=[phys[j] for j in J]
        if rank_cols(removed+cs,p) < pr+rank_cols(cs,p):
            res.append(J)
    return res


def x_distance(rows,p,two_generations=True):
    removed,physical,r,k=punctured(rows,p,two_generations)
    minwt=len(physical)+1; num=0
    for coeff in product(range(p),repeat=len(rows)):
        if all(dot(coeff,u)%p==0 for u in removed):continue
        w=sum(dot(coeff,u)%p!=0 for u in physical)
        if w<minwt:minwt=w;num=0
        if w==minwt:num+=1
    return minwt,num


def dimensions(rows,p,two=True):
    removed,physical,r,k=punctured(rows,p,two)
    return {'n':len(physical),'k':k,'charge_rank':r,
            'X_stabilizer_rank':r-k,'Z_stabilizer_rank':len(physical)-r}


def rref_subspaces(m,p):
    """Exactly one RREF basis for every m-plane in F_p^5."""
    for pivots in combinations(range(5),m):
        free=[(i,j) for i in range(m) for j in range(pivots[i]+1,5)
              if j not in pivots]
        for params in product(range(p),repeat=len(free)):
            B=[[0]*5 for _ in range(m)]
            for i,j in enumerate(pivots):B[i][j]=1
            for (i,j),v in zip(free,params):B[i][j]=v
            yield B


def bad_cross_generation_pairs(B,p):
    """Indices of matched species giving weight-2 logical Z errors.

    For p odd, a 2-column relation with two generations must involve
    identical SO(10) weights, with opposite coefficients. Such a pair
    is logically nontrivial iff B*q is a nonzero scalar times B*q0.
    """
    b0=tuple(dot(b,Q0)%p for b in B)
    assert any(b0)
    failure=[]
    for i,q in enumerate(weights()[1:]):
        u=tuple(dot(b,q)%p for b in B)
        if any(u) and rank_cols([b0,u],p)==1:failure.append(i)
    return failure


def basic_check(name,rows,p,two=True,need_x=False):
    if p != 2:
        valid,wit=triorthogonal_integer(rows)
        assert valid,('cubic anomaly',name,wit)
    data=dimensions(rows,p,two)
    rem,phys,_,_=punctured(rows,p,two)
    d=next(w for w in (1,2,3) if z_logical_supports(rem,phys,p,w))
    data['Z_distance']=d
    if d==3:
        third=z_logical_supports(rem,phys,p,3)
        data['weight_3_logical_supports']=len(third)
        if two:
            data['weight_3_partition']={str(a):n for a,n in sorted(Counter(
                sum(i<15 for i in support) for support in third).items())}
    if need_x:
        data['X_distance'],data['min_X_count']=x_distance(rows,p,two)
    assert d >= 2 or name=='invalid'
    print(name, 'p=',p, data)
    return data


def main():
    q=weights()
    assert len(q)==16 and len(set(q))==16
    assert all(sum(x==-1 for x in u)%2==0 for u in q)
    assert all(sum(u[i] for u in q)==0 for i in range(5))
    assert all(sum(u[i]*u[j]*u[k] for u in q)==0
               for i in range(5) for j in range(5) for k in range(5))

    # Paper's one-generation code; the muon/charm/strange family is identical
    # to electron/up/down under gauge charges, so only labels change.
    for p in (3,5):
        d=basic_check('second generation alone',rows_full(two_generations=False),p,False)
        assert (d['n'],d['k'],d['Z_distance'],d['weight_3_logical_supports'])==(15,1,3,25)

    # Shared gauge charges do not resolve both punctured neutrinos.
    d=basic_check('shared gauge only',rows_full(),5)
    assert (d['n'],d['k'],d['Z_distance'])==(30,1,3)

    # Anomaly-free generation-number difference resolves both, but gives
    # undetectable Z errors between matched generation species.
    d=basic_check('shared gauge + family number',rows_full(flavor=True),5)
    assert (d['n'],d['k'],d['Z_distance'])==(30,2,2)

    # New family-sensitive differential Cartan constructions.
    B3=[[1,0,0,0,0],[0,1,0,0,1],[0,0,1,1,0]]
    B5=[[1,0,0,0,0],[0,1,1,1,1]]
    for p,B in ((3,B3),(5,B5),(7,B5)):
        assert bad_cross_generation_pairs(B,p)==[]
        d=basic_check('coupled differential Cartan',rows_full(B=B),p,need_x=(p in (3,5)))
        assert (d['n'],d['k'],d['Z_distance'],d['weight_3_logical_supports'])==(30,2,3,50)
        assert d['weight_3_partition']=={'0':25,'3':25}

    # The full separated symmetry is the direct product of two copies.
    eye=[[int(i==j) for j in range(5)] for i in range(5)]
    d=basic_check('direct sum',rows_full(B=eye),3)
    assert (d['n'],d['k'],d['Z_distance'],d['weight_3_logical_supports'])==(30,2,3,50)

    # Finite-field exhaustiveness for rank-1 and rank-2 differential
    # Cartan subspaces in the qutrit family: only subspaces that resolve
    # the two removed neutrino coordinates are included.
    for m,expected in ((1,81),(2,1080)):
        tested=0
        for B in rref_subspaces(m,3):
            if not any(dot(b,Q0)%3 for b in B):continue
            tested+=1
            assert bad_cross_generation_pairs(B,3),('counterexample',B)
        assert tested==expected,(m,tested)
        print(f'EXHAUSTIVE p=3 m={m}: {tested} admissible row-spaces, all have distance-2 logical pairs')

    # Binary case uses the paper's bit construction, not the odd-prime
    # reduction of +/-1 charges.
    bitrows=[[1]*16]+[[(1-u[i])//2 for u in q] for i in range(4)]
    assert all(sum((sum(c*row[j] for c,row in zip(coeff,bitrows))%2)
                   for j in range(16))%8==0 for coeff in product((0,1),repeat=5))
    d=basic_check('second generation qubit Reed-Muller',bitrows,2,False,True)
    assert (d['n'],d['k'],d['Z_distance'],d['X_distance'])==(15,1,3,7)
    print('ALL CHECKS PASSED')

if __name__=='__main__':main()
