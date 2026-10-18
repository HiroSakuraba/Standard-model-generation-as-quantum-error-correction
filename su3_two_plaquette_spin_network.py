"""Two adjacent SU(3) square plaquettes: exact six-vertex Gauss certificates.

Seven-link pure-gauge graph; Peter-Weyl link representation labels
1, 3, anti-3. Retains vacuum and oriented simple flux loops L/R/outer.
Gauss checked by direct local SU(3) generator actions on the COLOR INDICES.
Projected Wilson matrix elements are Haar-exact in this subspace (including
shared-link perimeter amplitude 1/3), not invented vacuum-only hopping.
This is a seven-dimensional gauge-invariant variational TRUNCATION, not the
complete minimal-truncated link theory or a glueball mass calculation.
"""
from collections import defaultdict
from itertools import product
import json
import sympy as s
import numpy as np

# Graph: 0--1--2 on top, 3--4--5 bottom. Edge 5 is shared.
EDGES=((0,1),(1,2),(3,4),(4,5),(0,3),(1,4),(2,5))
CYCLES={
 'L':((0,1),(5,1),(2,-1),(4,-1)),
 'R':((1,1),(6,1),(3,-1),(5,-1)),
 'O':((0,1),(1,1),(6,1),(3,-1),(2,-1),(4,-1)),
}
LABELS=('vac','L+','L-','R+','R-','O+','O-')

def states_for_loop(cycle):
    """Sparse full seven-link wavefunction, UNNORMALIZED exact integer coefficients."""
    if cycle is None: return {tuple([(0,0,0)]*7):1}
    edges=list(cycle)
    vertex_order=[]
    for ei,_ in edges:
        for v in EDGES[ei]:
            if v not in vertex_order: vertex_order.append(v)
    out={}
    for colors in product(range(3),repeat=len(vertex_order)):
        by_vertex=dict(zip(vertex_order,colors))
        links=[(0,0,0)]*7
        for ei,orientation in edges:
            v,w=EDGES[ei]
            links[ei]=(orientation,by_vertex[v],by_vertex[w])
        key=tuple(links)
        assert key not in out
        out[key]=1
    return out

def wavefunctions():
    out={'vac':states_for_loop(None)}
    for k,cycle in CYCLES.items():
        for sign in (1,-1):
            out[k+('+' if sign==1 else '-')]=states_for_loop([(edge,d*sign) for edge,d in cycle])
    return out

def local_gauss_action(vec,vertex,kind):
    """Chevalley raising/lowering and two independent Cartan generators.
    Anti-fundamental action -T^transpose is correct at right link ends.
    """
    pairs={'E01':(0,1),'E10':(1,0),'E12':(1,2),'E21':(2,1)}
    diag={'H01':(1,-1,0),'H12':(0,1,-1)}
    out=defaultdict(int)
    for config,amp in vec.items():
        for ei,(a,b) in enumerate(EDGES):
            if vertex not in (a,b):continue
            rep,left,right=config[ei]
            if rep==0:continue
            source=vertex==a
            fundamental=(source and rep==1) or ((not source) and rep==-1)
            idx=1 if source else 2
            color=config[ei][idx]
            if kind in diag:
                val=diag[kind][color]*(1 if fundamental else -1)
                if val: out[config]+=amp*val
            else:
                r,c=pairs[kind]
                # anti repr -T^transpose
                from_index=c if fundamental else r
                to_index=r if fundamental else c
                if color!=from_index:continue
                moved=list(config)
                entry=list(moved[ei]);entry[idx]=to_index;moved[ei]=tuple(entry)
                out[tuple(moved)]+=amp*(1 if fundamental else -1)
    return {k:v for k,v in out.items() if v}

def exact_gauss_test():
    vectors=wavefunctions()
    norms={label:len(psi) for label,psi in vectors.items()}
    assert norms == dict(vac=1,**{'L+':81,'L-':81,'R+':81,'R-':81,'O+':729,'O-':729})
    assert len(set(frozenset(v.keys()) for v in vectors.values()))==7
    checks=0
    for name,psi in vectors.items():
        for v in range(6):
            for kind in ('E01','E10','E12','E21','H01','H12'):
                residual=local_gauss_action(psi,v,kind)
                assert not residual,(name,v,kind,len(residual))
                checks+=1
    print('Gauss constraints: 7 basis vectors x 6 vertices x 6 exact Lie generators =',checks,'zero residuals')
    return norms

# Independent Haar-polynomial integration in two fundamental SU(3) holonomies.
# In maximal tree gauge, L holonomy is U and R holonomy V, while the exterior
# contour is UV. The 7 basis states are corresponding normalized characters.
# Each term represents product of (group, conjugation, row, col) matrix entries.

def character(name):
    if name=='vac': return [((),1)]
    group,sign=name[:-1],name[-1]
    star=(sign=='-')
    if group=='L': return [((('U',star,i,i),),1) for i in range(3)]
    if group=='R': return [((('V',star,i,i),),1) for i in range(3)]
    if group=='O': return [((('U',star,i,j),('V',star,j,i)),1) for i in range(3) for j in range(3)]
    raise ValueError(name)

def dagger(monomial):return tuple((g,not star,i,j) for g,star,i,j in monomial)

def epsilon(x):
    if len(set(x))<3:return 0
    seq=list(x)
    return (-1)**sum(seq[i]>seq[j] for i in range(3) for j in range(i+1,3))

def haar_one(factors):
    plus=[(i,j) for g,star,i,j in factors if not star]
    minus=[(i,j) for g,star,i,j in factors if star]
    p,q=len(plus),len(minus)
    if p==q==0: return s.Integer(1)
    if p==q==1:
        return s.Rational(int(plus[0]==minus[0]),3)
    if p==3 and q==0:
        return s.Rational(epsilon(tuple(a for a,b in plus))*epsilon(tuple(b for a,b in plus)),6)
    if p==0 and q==3:
        return s.Rational(epsilon(tuple(a for a,b in minus))*epsilon(tuple(b for a,b in minus)),6)
    # Only degrees <=3 occur in this projected Wilson calculation.
    assert max(p,q)<=3,(p,q)
    return s.Integer(0)

def integrate(factors):
    return haar_one([f for f in factors if f[0]=='U'])*haar_one([f for f in factors if f[0]=='V'])

def w_matrix(which):
    op=character(which)
    mat=s.zeros(7)
    for i,bra in enumerate(LABELS):
        for j,ket in enumerate(LABELS):
            score=s.Integer(0)
            for b,bc in character(bra):
                for o,oc in op:
                    for k,kc in character(ket):
                        score+=bc*oc*kc*integrate(dagger(b)+o+k)
            mat[i,j]=score
    return mat

def commutant_dim(matrices):
    n=matrices[0].rows
    equations=s.Matrix.vstack(*[s.kronecker_product(s.eye(n),a)-s.kronecker_product(a.T,s.eye(n)) for a in matrices])
    rank=s.polys.matrices.DomainMatrix.from_Matrix(equations).convert_to(s.QQ).rank()
    return n*n-rank

def algebra_dim(generators):
    n=generators[0].rows
    queue=[s.eye(n)];basis=[];vec=s.zeros(n*n,0)
    while queue:
        mat=queue.pop(0)
        col=s.Matrix(n*n,1,list(mat))
        if s.Matrix.hstack(vec,col).rank()==vec.cols:continue
        vec=vec.row_join(col);basis.append(mat)
        if len(basis)==n*n:break
        queue.extend([g*mat for g in generators])
    return len(basis)

def test():
    norms=exact_gauss_test()
    # Independent two-holonomy Haar integration confirms orthonormal spin networks.
    gram=s.Matrix(7,7,lambda i,j:sum(
        a*b*integrate(dagger(left)+right)
        for left,a in character(LABELS[i])
        for right,b in character(LABELS[j])))
    assert gram==s.eye(7)
    Wl,Wr=w_matrix('L+'),w_matrix('R+')
    assert Wl.T==w_matrix('L-') and Wr.T==w_matrix('R-')
    assert Wl[1,0]==1 and Wr[3,0]==1
    assert Wl[5,3]==s.Rational(1,3) and Wr[5,1]==s.Rational(1,3)
    assert Wl[4,6]==s.Rational(1,3) and Wr[2,6]==s.Rational(1,3)
    # Two independent plaquette insertions create an outer perimeter flux.
    assert (Wr*Wl)[5,0] == (Wl*Wr)[5,0] == s.Rational(1,3)
    # Full Wilson multiplication operators commute. Their hard projected
    # 7D representatives need not commute: truncation, not SU(3), causes it.
    assert Wl*Wr != Wr*Wl
    print('W_left (exact projected fundamental Wilson multiplication):');print(Wl)
    print('W_right:');print(Wr)
    # Exact local electric link occupation operators; E²=(4/3)N.
    N=[]
    for ei in range(7):
        vals=[int(label!='vac' and ei in dict(CYCLES[label[0]]).keys()) for label in LABELS]
        N.append(s.diag(*vals))
    D=sum((s.Rational(4,3)*m for m in N),s.zeros(7))
    assert D==s.diag(0,*[s.Rational(16,3)]*4,*[s.Integer(8)]*2)
    assert N[5][1,1]==N[5][3,3]==1 and N[5][5,5]==0
    # Different local positions distinguish the two face loops in a gauge invariant way.
    assert N[0][1,1]==1 and N[0][3,3]==0
    Bleft=Wl+Wl.T;Bright=Wr+Wr.T
    print('Electric E =',list(D.diagonal()))
    assert commutant_dim([D,Bleft,Bright])==2
    assert algebra_dim([D,Bleft,Bright])==25
    charge=s.zeros(7);charge[0,0]=1
    for i,j in [(1,2),(3,4),(5,6)]:charge[i,j]=charge[j,i]=1
    assert charge*charge==s.eye(7) and all(charge*a==a*charge for a in (D,Bleft,Bright))
    assert (s.eye(7)+charge).rank()==4 and (s.eye(7)-charge).rank()==3
    # 4+3 parity sectors and 25=4²+3² certify M4 ⊕ M3, with
    # NO noiseless logical multiplicity under this FULL generated algebra.
    assert commutant_dim([D,Bleft,Bright]+N)==2
    print('Hermitian-electric/Wilson algebra: M4 + M3, dimension=25, commutant=2; no noiseless logical qubit')
    print('commutant(E,Bleft,Bright)',commutant_dim([D,Bleft,Bright]))
    print('commutant(E,Bleft,Bright,Nedges)',commutant_dim([D,Bleft,Bright]+N))
    assert commutant_dim([D,Wl,Wr,Wl.T,Wr.T])==1
    assert algebra_dim([D,Wl,Wr,Wl.T,Wr.T])==49
    print('commutant(E,Wl,Wr,adjoints)',commutant_dim([D,Wl,Wr,Wl.T,Wr.T]))
    # This pair is exactly decoherence-free against ANY products of
    # diagonal local electric-Casimir measurements, including the shared link.
    # The standard gauge-invariant real plaquette magnetic term destroys it.
    K=s.eye(7)[:,[1,2]]
    for A in N+[D]:
        Q=K.T*A*K
        assert Q==Q[0,0]*s.eye(2)
    assert K.T*Bleft*K==s.Matrix([[0,1],[1,0]])
    print('Electric-only oriented left-loop qubit: PASS; local magnetic Wilson loop acts as logical X: FAIL')
    print('algebra(E,Wl,Wr,adjoints)',algebra_dim([D,Wl,Wr,Wl.T,Wr.T]))
    M=np.array(D-s.Rational(1,2)*(Bleft+Bright),float)
    vals=np.linalg.eigvalsh(M)
    eigs,V=np.linalg.eigh(M)
    Q=V[:,:2]
    projected_e=Q.T@np.array(D,float)@Q
    electric_violation=np.linalg.norm(projected_e-np.trace(projected_e)/2*np.eye(2),ord=2)
    assert electric_violation > 2.56
    print('First two levels fail electric-only Knill-Laflamme; traceless compression norm',electric_violation)
    print('Spectrum =',vals,'gap=',vals[1]-vals[0])
    print('EXACT TWO-PLAQUETTE GAUSS AND WILSON TESTS PASSED')
    return Wl,Wr,D,N,vals

if __name__=='__main__': test()
