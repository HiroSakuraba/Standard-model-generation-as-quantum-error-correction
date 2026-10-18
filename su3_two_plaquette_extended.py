"""SU(3) two-holonomy Wilson basis enlargement; exact SU(3) Haar moments.

Gauge-fixed two plaquette holonomies U,V with simultaneous conjugation.
Basis includes seven original loop states, four simultaneous-loop traces,
and higher face characters (6,bar6,8). It is not a complete link-irrep
truncation or physical glueball spectrum. Exact color fusion and Wilson
matrix elements; numerical comparisons separate from exact algebra.
"""
from functools import lru_cache
from itertools import permutations, product
from collections import defaultdict
from math import factorial, prod
import sympy as sp
import numpy as np

N=3


def inv(p):
    x=[0]*len(p)
    for i,v in enumerate(p):x[v]=i
    return tuple(x)

def cycle_type(p):
    seen=set();out=[]
    for i in range(len(p)):
        if i not in seen:
            j=i;n=0
            while j not in seen:
                seen.add(j); n+=1;j=p[j]
            out.append(n)
    return tuple(sorted(out,reverse=True))

def partitions(n,maxpart=None):
    if n==0: yield ();return
    for i in range(min(n,maxpart or n),1-1,-1):
        for tail in partitions(n-i,i):yield (i,)+tail

def hook_data(lam):
    hooks=[];content=1
    for i,row in enumerate(lam):
        for j in range(row):
            hook=row-j+sum(1 for below in lam[i+1:] if below>j)
            hooks.append(hook)
            content*=N+j-i
    dim=sp.Rational(content, prod(hooks))
    f=factorial(sum(lam))//prod(hooks)
    return dim,f

def rim_char(lam,cyc):
    # Murnaghan-Nakayama via subpartition strips, tiny k <= 4.
    if not cyc: return int(not lam)
    target=cyc[0]; old=sum(lam)
    out=0
    for mu in partitions(old-target):
        if len(mu)>len(lam) or any(mu[i]>lam[i] for i in range(len(mu))):continue
        cells={(i,j) for i,row in enumerate(lam) for j in range(row)}-{
            (i,j) for i,row in enumerate(mu) for j in range(row)}
        if len(cells)!=target:continue
        # no 2x2 block, connected in grid
        if any((i+1,j) in cells and (i,j+1) in cells and (i+1,j+1) in cells for i,j in cells):continue
        found={next(iter(cells))}
        while True:
            extended=found|{q for i,j in found for q in ((i-1,j),(i+1,j),(i,j-1),(i,j+1)) if q in cells}
            if extended==found:break
            found=extended
        if found!=cells:continue
        height=len(set(i for i,j in cells))-1
        out+=(-1)**height*rim_char(mu,cyc[1:])
    return out

@lru_cache(None)
def wg_map(k):
    if not k:return {():sp.Integer(1)}
    vals={}
    for p in permutations(range(k)):
        typ=cycle_type(p)
        if typ not in vals:
            z=sp.Rational(0)
            for lam in partitions(k):
                if len(lam)>N:continue
                dim,f=hook_data(lam)
                z+=sp.Rational(f*f, factorial(k)**2)/dim*rim_char(lam,typ)
            vals[typ]=sp.factor(z)
    return vals


def signed_permutations(k):
    for x in permutations(range(k)):
        parity=sum(x[i]>x[j] for i in range(k) for j in range(i+1,k))
        yield x,(-1)**parity


# Trace words are built with symbols shared along their cyclic index lines.
# Each word is (sequence of group symbols, conjugate 0/1).
# Polynomial represented as list of terms (coefficient, tuple(words)).

OPS=('vac','L+','L-','R+','R-','O+','O-',
     'L+R+','L+R-','L-R+','L-R-',
     'L6','L6b','R6','R6b','L8','R8','X+','X-')

def polynomial(label):
    if label=='vac':return [(1,())]
    if label=='L+':return [(1,(('U',False),))]
    if label=='L-':return [(1,(('U',True),))]
    if label=='R+':return [(1,(('V',False),))]
    if label=='R-':return [(1,(('V',True),))]
    if label=='O+':return [(1,(('UV',False),))]
    if label=='O-':return [(1,(('UV',True),))]
    if label=='X+':return [(1,(('Uv',False),))]
    if label=='X-':return [(1,(('uV',False),))]
    if len(label)==4 and label[0]=='L' and label[2]=='R':
        return [(1,(('U',label[1]=='-'),('V',label[3]=='-')))]
    group="U" if label[0]=="L" else "V"
    if label.endswith('6') or label.endswith('6b'):
        star=label.endswith('6b')
        return [(1,((group,star),(group,star))),(-1,((group,not star),))]
    if label.endswith('8'):
        return [(1,((group,False),(group,True))),(-1,())]
    raise KeyError(label)

def monomial(words):
    factors=[];i=0
    for word,star in words:
        # Each word W_0...W_m cyclic trace.
        indices=[('v',i+j) for j in range(len(word))];i+=len(word)
        for j,group in enumerate(word):
            factors.append((group.upper(),star ^ group.islower(),indices[j],indices[(j+1)%len(word)]))
    return factors

@lru_cache(maxsize=None)
def moment_terms(group_factors):
    # For one group, return a linear combination of Kronecker constraints.
    plus=[(r,c) for star,r,c in group_factors if not star]
    minus=[(r,c) for star,r,c in group_factors if star]
    p,q=len(plus),len(minus)
    if (p-q)%3:return ()
    if p==q==0:return ((sp.Integer(1),()),)
    copies=abs(p-q)//3
    assert copies<=1,(p,q)
    if p>q:
        minus.extend([(('det',i),i) for i in range(3)])
        twist=True
    elif q>p:
        plus.extend([(('det',i),i) for i in range(3)])
        twist=True
    else: twist=False
    k=len(plus)
    assert k==len(minus) and k<=4,(p,q,k)
    ps=list(permutations(range(k)))
    wg=wg_map(k)
    terms=defaultdict(lambda:sp.Integer(0))
    for dperm,dsgn in (list(signed_permutations(3)) if twist else [((),1)]):
        def convert(pair):
            a,b=pair
            if isinstance(a,tuple) and a[:1]==('det',):
                return (('fix',a[1]),('fix',dperm[a[1]]))
            return (a,b)
        aa=list(map(convert,plus));bb=list(map(convert,minus))
        for sigma in ps:
            s_inv=inv(sigma)
            left=[(aa[i][0],bb[sigma[i]][0]) for i in range(k)]
            for tau in ps:
                pi=tuple(tau[s_inv[i]] for i in range(k))
                coeff=dsgn*wg[cycle_type(pi)]
                if not coeff:continue
                edges=left+[(aa[i][1],bb[tau[i]][1]) for i in range(k)]
                # constraints may be redundant; canonical equivalence handled on union.
                terms[tuple(edges)]+=coeff
    return tuple((c,e) for e,c in terms.items() if c)

def partition_count(nodes, equations):
    """Sum free SU(3) color indices under exact Weingarten delta constraints."""
    root={x:x for x in nodes}
    def find(x):
        while root[x]!=x:
            root[x]=root[root[x]];x=root[x]
        return x
    for a,b in equations:
        ra,rb=find(a),find(b)
        if ra!=rb:root[ra]=rb
    groups=defaultdict(set)
    for node in nodes:
        r=find(node)
        if isinstance(node,tuple) and node[0]=='fix':groups[r].add(node[1])
        else:groups[r]
    if any(len(values)>1 for values in groups.values()):return 0
    return 3**sum(len(values)==0 for values in groups.values())


@lru_cache(maxsize=None)
def integrate_words(words):
    factors=monomial(words)
    u=tuple(x[1:] for x in factors if x[0]=='U')
    v=tuple(x[1:] for x in factors if x[0]=='V')
    left=moment_terms(u);right=moment_terms(v)
    if not left or not right:return sp.Integer(0)
    nodes={x for _,_,r,c in factors for x in (r,c)}
    result=sp.Rational(0)
    for cu,eu in left:
        for cv,ev in right:
            edges=eu+ev
            allnodes=nodes|{z for a,b in edges for z in (a,b)}
            val=partition_count(allnodes,edges)
            if val:result+=cu*cv*val
    return sp.factor(result)

def integrand(ops):
    ops=tuple(ops)
    terms=[(sp.Integer(1),())]
    for label,dagger in ops:
        upcoming=[]
        for coefficient,words in polynomial(label):
            new=tuple((word,not star) for word,star in words) if dagger else words
            for a,b in terms:upcoming.append((a*coefficient,b+new))
        terms=upcoming
    accum=sp.Rational(0)
    for coeff,words in terms:
        accum+=coeff*integrate_words(words)
    return sp.factor(accum)

def audit(selected=OPS):
    import time
    t0=time.time()
    m=len(selected)
    G=sp.Matrix(m,m,lambda i,j:integrand(((selected[i],True),(selected[j],False))))
    print('GRAM',m,'rank',G.rank(),'time',round(time.time()-t0,1),flush=True)
    assert G==G.T and G.det()!=0
    print('G off-diagonal',[(selected[i],selected[j],G[i,j]) for i in range(m) for j in range(i+1,m) if G[i,j]],flush=True)
    WL=sp.Matrix(m,m,lambda i,j:integrand(((selected[i],True),('L+',False),(selected[j],False))))
    print('WLEFT',m,'time',round(time.time()-t0,1),flush=True)
    WR=sp.Matrix(m,m,lambda i,j:integrand(((selected[i],True),('R+',False),(selected[j],False))))
    print('WRIGHT',m,'time',round(time.time()-t0,1),flush=True)
    assert all(x.is_Rational for x in G) and all(x.is_Rational for x in WL)
    print('WL[:,vac]',[(selected[i],WL[i,0]) for i in range(m) if WL[i,0]])
    assert WL[1,0]==1 and WR[3,0]==1
    # generalized projected operator in nonorthonormal character basis
    WU=G.inv()*WL; WV=G.inv()*WR
    assert G*WU==WL and G*WV==WR
    print('PROJECTED COMMUTATOR RANK', (WU*WV-WV*WU).rank(),flush=True)
    return G,WU,WV

def casimir(label):
    if label.endswith(('6','6b')):return sp.Rational(10,3)
    if label.endswith('8'):return sp.Integer(3)
    if label=='vac':return sp.Integer(0)
    return sp.Rational(4,3)


def electric_cov(G,labels):
    n=len(labels);idx={lab:i for i,lab in enumerate(labels)}
    def reps(lab):
        if lab=='vac':return (0,0)
        if lab[0]=='L' and not lab.startswith('L+') and not lab.startswith('L-'):
            return (casimir(lab),0)
        if lab[0]=='R' and not lab.startswith('R+') and not lab.startswith('R-'):
            return (0,casimir(lab))
        if lab in ('L+','L-'):return (sp.Rational(4,3),0)
        if lab in ('R+','R-'):return (0,sp.Rational(4,3))
        return (sp.Rational(4,3),sp.Rational(4,3))
    vals=[reps(x) for x in labels]
    U=G*sp.diag(*(a for a,b in vals));V=G*sp.diag(*(b for a,b in vals))
    assert U==U.T and V==V.T
    E=sp.zeros(n)
    for lab in labels:
        i=idx[lab]
        if lab.startswith('O') or lab=='vac':value=0
        elif lab in ('L+R+','L-R-','L+R-','L-R+','X+','X-'):value=None
        else:value=sum(vals[i])
        if value is not None:E[i,i]=value
    # Physical shared-edge irreps at opposite-direction strands:
    #    3 x bar3 = 1 + 8, Casimir 0 or 3;
    # orthogonal outer trace is the shared-link singlet.
    for lab in ('L+R+','L-R-'):
        if lab in idx:E[idx[lab],idx[lab]]=sp.Rational(8,3)
    # Physical shared-edge irreps with same-directed strands:
    #   3 x 3 = bar3 + 6; E=7/3 I + swap on the two strands.
    # In the nonorthogonal basis A=trU trV*, B=tr(UV*),
    #   Gram = [[1,1/3],[1/3,1]], so covariant E=G(7/3 I+swap).
    for a,b in (('L+R-','X+'),('L-R+','X-')):
        if a not in idx or b not in idx:
            if a in idx:E[idx[a],idx[a]]=sp.Rational(8,3)
            continue
        i,j=idx[a],idx[b]
        E[i,i]=E[j,j]=sp.Rational(8,3)
        E[i,j]=E[j,i]=sp.Rational(16,9)
    assert E==E.T
    assert all(x.is_Rational for x in E)
    # Exact seven-link electric Casimirs commute after gauge reduction.
    Uop,Vop,Eop=(G.inv()*x for x in (U,V,E))
    assert Uop*Vop==Vop*Uop
    assert Uop*Eop==Eop*Uop and Vop*Eop==Eop*Vop
    # Check positive semidefinite electric as a G-selfadjoint operator.
    e=G.inv()*E
    assert sp.factor(e.charpoly().as_expr()) is not None
    eigen=e.eigenvals()
    assert all(z>=0 for z in eigen)
    for a,b in (('O+','L+R+'),('O-','L-R-')):
        if a not in idx or b not in idx:continue
        i,j=idx[a],idx[b]
        S=G.extract([i,j],[i,j]);M=E.extract([i,j],[i,j])
        assert sorted((S.inv()*M).eigenvals())==[0,3]
    for a,b in (('L+R-','X+'),('L-R+','X-')):
        if a not in idx or b not in idx:continue
        i,j=idx[a],idx[b]
        S=G.extract([i,j],[i,j]);M=E.extract([i,j],[i,j])
        assert sorted((S.inv()*M).eigenvals())==[sp.Rational(4,3),sp.Rational(10,3)]
    return 3*U+3*V+E, U,V,E


def commutant_dimension(gens):
    eye=sp.eye(gens[0].rows);n=eye.rows
    eq=sp.Matrix.vstack(*[sp.kronecker_product(eye,a)-sp.kronecker_product(a.T,eye) for a in gens])
    return n*n-sp.polys.matrices.DomainMatrix.from_Matrix(eq).convert_to(sp.QQ).rank()


def capacity_audit(labels=OPS):
    import json
    from scipy.linalg import eigh
    G,WL,WR=audit(labels)
    n=len(labels)
    electric,U,V,E=electric_cov(G,labels)
    print('TOTAL E COV DIAG',[(label,electric[i,i]) for i,label in enumerate(labels)])
    # Convert covariant matrices to linear operators in nonorthogonal basis.
    Eop=G.inv()*electric
    Bleft=WL+G.inv()*WL.T*G
    Bright=WR+G.inv()*WR.T*G
    C=sp.zeros(n)
    for j,label in enumerate(labels):
        if len(label)==4 and label[0]=='L' and label[2]=='R':
            conj=label[0]+('+' if label[1]=='-' else '-')+'R'+('+' if label[3]=='-' else '-')
        elif label.endswith('+'):conj=label[:-1]+'-'
        elif label.endswith('-'):conj=label[:-1]+'+'
        elif label.endswith('6b'):conj=label[:-1]
        elif label.endswith('6'):conj=label+'b'
        else:conj=label
        assert conj in labels,(label,conj)
        C[labels.index(conj),j]=1
    assert C*C==sp.eye(n)
    assert all(C*a==a*C for a in (Eop,Bleft,Bright))
    print('CHARGE CONJUGATION',((sp.eye(n)+C)/2).rank(),((sp.eye(n)-C)/2).rank())
    for tag,gens in [('hermitian',[Eop,Bleft,Bright]),('full',[Eop,WL,WR,G.inv()*WL.T*G,G.inv()*WR.T*G]),('hermitian-plus-shared',[Eop,G.inv()*E,Bleft,Bright])]:
        dimension=commutant_dimension(gens)
        assert dimension == (1 if tag=='full' else 2),(tag,dimension)
        print('COMMUTANT',tag,dimension,flush=True)
    ix=[labels.index('L+'),labels.index('L-')]
    for gen in [U,V,E]:
        assert G.extract(ix,ix)==sp.eye(2)
        val=gen.extract(ix,ix)
        assert val[0,0]==val[1,1] and val[0,1]==val[1,0]==0
    assert G.extract(ix,range(n))*Bleft.extract(range(n),ix) == sp.Matrix([[0,1],[1,0]])
    assert G.extract(ix,range(n))*Bright.extract(range(n),ix) == sp.zeros(2)
    print('ELECTRIC-ONLY code exact PASS; Hermitian left Wilson logical X exact FAIL')
    # Generalized Hermitian eigenproblem uses positive exact G.
    Hcov=electric-sp.Rational(1,2)*(G*Bleft+G*Bright)
    A=np.asarray(Hcov,float);B=np.asarray(G,float)
    vals,_=eigh(A,B)
    oldlabels=labels[:7]
    oldids=[labels.index(x) for x in oldlabels]
    vals7,_=eigh(A[np.ix_(oldids,oldids)],B[np.ix_(oldids,oldids)])
    print('EIGENVALUES first 6:',vals[:6]);print('FULL FIRST GAP',vals[1]-vals[0])
    print('SEVEN-STATE FIRST GAP',vals7[1]-vals7[0]);print('FULL GROUND ENERGY',vals[0],'SEVEN GROUND',vals7[0])
    assert vals[0]<=vals7[0]+1e-10
    # Seven-state covariant matrices reproduce their prior exact Wilson projectors.
    if labels[:7] == OPS[:7]:
        expected=sp.zeros(7)
        for i,j,value in ((1,0,1),(2,1,1),(0,2,1),(5,3,sp.Rational(1,3)),(4,6,sp.Rational(1,3))):
            expected[i,j]=value
        assert (G*WL)[:7,:7]==expected  # independent previous 7D exact transitions
    return G,WL,WR,electric

def microscopic_checks():
    # 3x3 fundamental SU(3) Lie-algebra generators, trace T^a T^b=delta/2.
    gs=[]
    for i,j in ((0,1),(0,2),(1,2)):
        A=sp.zeros(3);B=sp.zeros(3)
        A[i,j]=A[j,i]=sp.Rational(1,2)
        B[i,j]=-sp.I/2;B[j,i]=sp.I/2
        gs.extend([A,B])
    gs.extend([sp.diag(1,-1,0)/2,sp.diag(1,1,-2)/(2*sp.sqrt(3))])
    assert len(gs)==8
    for i,j,k,l in product(range(3),repeat=4):
        rhs=(sp.KroneckerDelta(i,l)*sp.KroneckerDelta(j,k)
             -sp.KroneckerDelta(i,j)*sp.KroneckerDelta(k,l)/3)/2
        assert sp.simplify(sum(T[i,j]*T[k,l] for T in gs)-rhs)==0
    assert integrand((('O+',True),('O+',False)))==1
    assert integrand((('O+',True),('L+R+',False)))==sp.Rational(1,3)
    assert integrand((('X+',True),('L+R-',False)))==sp.Rational(1,3)
    assert integrand((('L6',True),('L6',False)))==1
    assert integrand((('L8',True),('L8',False)))==1
    assert integrand((('L+',False),)*3)==1  # SU3 epsilon determinant moment
    print('MICROSCOPIC FIERZ + SU3 HAAR/determinant CHECKS PASSED')


def main():
    microscopic_checks()
    for n in (7,11,17,19):
        print('\n====== TRUNCATION',n,'======',flush=True)
        G,WL,WR,D=capacity_audit(OPS[:n])
    print('ALL EXACT EXTENDED TWO-PLAQUETTE SU3 TESTS PASSED')


if __name__=='__main__':
    main()
