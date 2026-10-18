"""Complete SU(3) two-plaquette Gauss-singlet sector with link irreps 1,3,bar3.

Exact enumeration of EVERY per-link representation labeling and local
intertwiner multiplicity. Wilson multipliers are projected from the exact
SU(3) Haar measure in two holonomies. This is a specific 2-plaquette
Peter-Weyl truncation, not a continuum glueball or a calibrated device.

Requires sympy, numpy and su3_two_plaquette_spin_network.py (Gauss helper).
"""
from collections import Counter, defaultdict
from itertools import product, permutations
from functools import lru_cache
import json
import sympy as sp
import numpy as np
from su3_two_plaquette_spin_network import (EDGES, CYCLES, LABELS as OLD_LABELS,
    wavefunctions, local_gauss_action, w_matrix, exact_gauss_test)

NEW_LABELS = OLD_LABELS + ('B+','B-')


def vertex_irreps(labels, v):
    return tuple((i, labels[i] if v==a else -labels[i])
                 for i,(a,b) in enumerate(EDGES) if v in (a,b))


def invariant_multiplicity(irreps):
    used=[rep for _,rep in irreps if rep]
    if not used:return 1
    if len(used)==2:return int(used[0]==-used[1])
    if len(used)==3:return int(abs(sum(used))==3)
    return 0


def enumerate_basis():
    valid=[];weights=Counter();bad=Counter()
    for labels in product((-1,0,1),repeat=len(EDGES)):
        local=[invariant_multiplicity(vertex_irreps(labels,v)) for v in range(6)]
        multiplicity=np.prod(local)
        if multiplicity:valid.append(labels);weights[sum(x!=0 for x in labels)]+=int(multiplicity)
        else:bad[sum(x!=0 for x in labels)]+=1
    assert sum(weights.values())==9
    assert len(valid)==9
    assert sum(weights.values())+sum(bad.values())==3**len(EDGES)==2187
    return valid,weights


def loop_pattern(name):
    if name=='vac':return (0,)*7
    path=CYCLES[name[0]];sgn=1 if name[1]=='+' else -1
    out=[0]*7
    for edge,orient in path:out[edge]=sgn*orient
    return tuple(out)


def boundary_network(label_assignment):
    """Exact nonzero link-matrix coefficients from delta/epsilon intertwiners.

    Construct independent color assignments to vertices, then link-end
    representations, not an inferred trace-polynomial Gauss test.
    """
    choices=[]
    for v in range(6):
        active=vertex_irreps(label_assignment,v)
        active=tuple((e,ir) for e,ir in active if ir)
        if not active:choices.append([((),1)]);continue
        if len(active)==2:
            assert active[0][1]==-active[1][1]
            choices.append([((c,c),1) for c in range(3)])
        elif len(active)==3:
            assert abs(sum(ir for _,ir in active))==3
            def parity(perm):
                return (-1)**sum(perm[i]>perm[j] for i in range(3) for j in range(i+1,3))
            choices.append([(p,parity(p)) for p in permutations(range(3))])
        else:raise AssertionError((v,active))
    out={}
    for chosen in product(*choices):
        data={}; amp=1
        for v,((colors,weight)) in enumerate(chosen):
            amp*=weight
            edges=[e for e,r in vertex_irreps(label_assignment,v) if r]
            for ei,color in zip(edges,colors):data[ei,v]=color
        link=[]
        for i,(a,b) in enumerate(EDGES):
            rep=label_assignment[i]
            link.append((rep,data[i,a],data[i,b]) if rep else (0,0,0))
        key=tuple(link)
        assert key not in out
        out[key]=amp
    return out


def microscopic_certificate(valid):
    old={loop_pattern(name):name for name in OLD_LABELS}
    new=[x for x in valid if x not in old]
    assert len(new)==2 and all(all(t!=0 for t in x) for x in new)
    # Choose B+ as the configuration in which the first perimeter edge
    # carries a fundamental color orientation, B- as its conjugate.
    positive=next(x for x in new if x[0]==1)
    negative=next(x for x in new if x[0]==-1)
    assert tuple(-x for x in positive)==negative
    mapped={**{name:loop_pattern(name) for name in OLD_LABELS}, 'B+':positive,'B-':negative}
    assert set(mapped.values())==set(valid)
    oldvectors=wavefunctions()
    exact_gauss_test()  # independently recheck the previous 252 constraints
    for lab in OLD_LABELS:
        assert set(boundary_network(mapped[lab]))==set(oldvectors[lab])
        vals=boundary_network(mapped[lab]);reference=oldvectors[lab]
        assert any(vals[k]!=reference[k] for k in vals) is False or all(vals[k]==-reference[k] for k in vals)
    gauss_tests=0; norms={}
    for lab in ('B+','B-'):
        state=boundary_network(mapped[lab]);norms[lab]=sum(abs(v)**2 for v in state.values())
        assert norms[lab]==2916
        for v in range(6):
            for gen in ('E01','E10','E12','E21','H01','H12'):
                residual=local_gauss_action(state,v,gen)
                assert not residual,(lab,v,gen,len(residual))
                gauss_tests+=1
    assert gauss_tests==72
    # Existing repo checks all 7 old states separately (252 constraints).
    return mapped,norms,gauss_tests


# Analytic SU(3)xSU(3) Haar moments needed for 9-basis elementary Wilson
# matrix entries: at most 3 fundamental matrix factors on either group.
def epsilon(xs):
    if len(set(xs))<3:return 0
    return (-1)**sum(xs[i]>xs[j] for i in range(3) for j in range(i+1,3))

@lru_cache(maxsize=None)
def haar_group(args):
    plus=[(r,c) for starred,r,c in args if not starred]
    minus=[(r,c) for starred,r,c in args if starred]
    p,q=len(plus),len(minus)
    if p==q==0:return sp.Integer(1)
    if p==q==1:
        return sp.Rational(int(plus[0]==minus[0]),3)
    if (p,q)==(3,0):
        return sp.Rational(epsilon([r for r,c in plus])*epsilon([c for r,c in plus]),6)
    if (p,q)==(0,3):
        return sp.Rational(epsilon([r for r,c in minus])*epsilon([c for r,c in minus]),6)
    assert p+q<=3,('Unexpected Haar degree',p,q)
    return sp.Integer(0)

@lru_cache(maxsize=None)
def integrate_monomial(factors):
    u=tuple((star,r,c) for group,star,r,c in factors if group=='U')
    v=tuple((star,r,c) for group,star,r,c in factors if group=='V')
    return haar_group(u)*haar_group(v)


def trace_parts(which):
    if which=='vac':return [(sp.Integer(1),())]
    if which in ('L+','L-','R+','R-'):
        group='U' if which.startswith('L') else 'V'
        dagger=which.endswith('-')
        return [(sp.Integer(1),((group,dagger,i,i),)) for i in range(3)]
    if which in ('O+','O-'):
        dagger=which.endswith('-')
        # tr(U V): sum U_ij V_ji; (tr(UV))* = sum conj U_ij conj V_ji.
        return [(sp.Integer(1), (('U',dagger,i,j),('V',dagger,j,i)))
                for i in range(3) for j in range(3)]
    assert which in ('B+','B-')
    dagger=which=='B-'
    parts=[]
    for i in range(3):
        for j in range(3):
            # trU tr Vdagger - tr(U Vdagger) and its conjugate.
            parts.append((sp.Integer(1), (('U',dagger,i,i),('V',not dagger,j,j))))
            parts.append((-sp.Integer(1), (('U',dagger,i,j),('V',not dagger,i,j))))
    return parts


def conjugate_monomial(mono):
    return tuple((group,not dagger,i,j) for group,dagger,i,j in mono)


def inner(op=None):
    """Exact covariant matrix <basis_i| (optional trU/trV) |basis_j>."""
    states=[trace_parts(lab) for lab in NEW_LABELS]
    middle=[(sp.Integer(1),())] if op is None else trace_parts(op)
    cov=sp.zeros(9)
    for i,bra in enumerate(states):
        for j,ket in enumerate(states):
            ans=sp.Rational(0)
            for a,b in bra:
                bc=conjugate_monomial(b)
                for o,m in middle:
                    for z,k in ket:
                        ans+=a*o*z*integrate_monomial(bc+m+k)
            cov[i,j]=ans
    return cov


def commutant_dim(mats):
    n=mats[0].rows;I=sp.eye(n)
    eq=sp.Matrix.vstack(*[sp.kronecker_product(I,m)-sp.kronecker_product(m.T,I) for m in mats])
    return n*n-sp.polys.matrices.DomainMatrix.from_Matrix(eq).convert_to(sp.QQ).rank()


def closure_dim(gens):
    n=gens[0].rows
    basis=[];flat=sp.zeros(n*n,0);queue=[sp.eye(n)]
    while queue:
        a=queue.pop(0);v=sp.Matrix(n*n,1,list(a))
        if flat.row_join(v).rank()==flat.cols:continue
        basis.append(a);flat=flat.row_join(v)
        if len(basis)==n*n:break
        queue.extend(g*a for g in gens)
    return len(basis)


def run():
    basis,distribution=enumerate_basis()
    mapped,norms,gauss_tests=microscopic_certificate(basis)
    print('Complete per-link irrep cutoff 1+3+anti3: 3^7=2187 rep assignments, Gauss singlets =',len(basis))
    print('Gauge-invariant states by occupied links:',dict(sorted(distribution.items())))
    print('New epsilon-junction configurations:',mapped['B+'],mapped['B-'])
    print('Both baryon color networks: 2916 integer amplitudes each; 252 old + 72 new = 324 exact Gauss checks PASS')
    G=inner()
    assert G==sp.diag(1,1,1,1,1,1,1,sp.Rational(4,3),sp.Rational(4,3)),G
    assert G.det()==sp.Rational(16,9)
    print('Trace-polynomial Gram: exact diagonal entries',list(G.diagonal()))
    Wl_cov=inner('L+');Wr_cov=inner('R+')
    Wl=G.inv()*Wl_cov;Wr=G.inv()*Wr_cov
    old_wl=w_matrix('L+');old_wr=w_matrix('R+')
    assert Wl_cov[:7,:7]==old_wl and Wr_cov[:7,:7]==old_wr
    assert Wl_cov.T==inner('L-') and Wr_cov.T==inner('R-')
    assert Wl_cov.T*G.inv()*G==Wl_cov.T
    print('Exact Haar fundamental Wilson matrices reproduce all 7-state entries')
    # Change to the actual normalized basis: each new raw baryon has norm^2=4/3.
    # The omitted block has two nonzero entries of magnitude 1/sqrt(3);
    # its induced 2-norm is exactly 1/sqrt(3), not a tiny truncation tail.
    S=sp.diag(*[sp.Integer(1)]*7,sp.sqrt(3)/2,sp.sqrt(3)/2)
    magnetic_w_l=S.inv()*Wl*S
    magnetic_w_r=S.inv()*Wr*S
    assert sp.simplify(magnetic_w_l[7:, :7]*magnetic_w_l[7:, :7].T-sp.eye(2)/3)==sp.zeros(2)
    assert sp.simplify(magnetic_w_r[7:, :7]*magnetic_w_r[7:, :7].T-sp.eye(2)/3)==sp.zeros(2)
    assert Wl*Wr != Wr*Wl # compression artifact, not full Wilson noncommutation
    print('Wilson discarded-block operator norm from 7-state host to junction states: 1/sqrt(3) exactly')
    print('Additional magnetic transitions from original states to baryon channels:')
    print('W_L cov[7:,:7] =',Wl_cov[7:,:7]);print('W_R cov[7:,:7] =',Wr_cov[7:,:7])
    NL=[]
    for edge in range(7):
        NL.append(sp.diag(*[int(mapped[lab][edge]!=0) for lab in NEW_LABELS]))
    D=sum((sp.Rational(4,3)*N for N in NL),sp.zeros(9))
    assert D==sp.diag(0,*[sp.Rational(16,3)]*4,*[sp.Integer(8)]*2,*[sp.Rational(28,3)]*2)
    assert all(g*D==D*g for g in NL)
    B_l=Wl+G.inv()*Wl_cov.T
    B_r=Wr+G.inv()*Wr_cov.T
    assert G*B_l == (G*B_l).T and G*B_r==(G*B_r).T
    # Charge conjugation maps all oriented flux and baryon pairs. It
    # commutes with Hermitian Wilson/electric operators but not all
    # complex fundamental Wilson multipliers.
    C=sp.eye(9)
    for i,j in ((1,2),(3,4),(5,6),(7,8)):
        C[i,i]=C[j,j]=0;C[i,j]=C[j,i]=1
    assert C*C==sp.eye(9)
    assert ((sp.eye(9)+C)/2).rank()==5
    assert ((sp.eye(9)-C)/2).rank()==4
    assert all(C*M==M*C for M in (D,B_l,B_r))
    print('Hermitian error algebra preserves charge-conjugation sectors: dimensions 5+4')
    print('Electric Casimir diag:',list(D.diagonal()))
    code=sp.eye(9)[:,[1,2]]
    assert all(code.T*N*code == N[1,1]*sp.eye(2) for N in NL)
    assert code.T*G*B_l*code==sp.Matrix([[0,1],[1,0]])
    print('Electric-only exact KL PASS, fundamental Wilson logical X remains FAIL')
    for tag,gens in [('Hermitian',[D,B_l,B_r]),('Hermitian+link-Casimirs',[D,B_l,B_r]+NL),('full complex',[D,Wl,Wr,G.inv()*Wl_cov.T,G.inv()*Wr_cov.T])]:
        d=commutant_dim(gens);print('Exact commutant',tag,d)
    # Numeric generalized Hermitian eigenproblems at several magnetic strengths.
    sq=np.diag([float(sp.sqrt(x)) for x in G.diagonal()]);sqinv=np.linalg.inv(sq)
    d=np.array(D,float);b=np.array(B_l+B_r,float)
    data=[]
    for coupling in (0.25,0.5,1,2,4):
        H=np.array(sq@ (d-coupling/2*b) @ sqinv,float)
        vals=np.linalg.eigvalsh((H+H.T)/2)
        h7=np.array((d-coupling/2*b)[:7,:7],float)
        v7=np.linalg.eigvalsh((h7+h7.T)/2)
        row=dict(g=float(coupling),nine_ground=float(vals[0]),nine_gap=float(vals[1]-vals[0]),seven_ground=float(v7[0]),seven_gap=float(v7[1]-v7[0]))
        data.append(row)
        print('Spectrum:',row)
    output=dict(dim=len(basis),state_labels=NEW_LABELS,configs={x:list(y) for x,y in mapped.items()},distribution=dict(distribution),gauss_new=gauss_tests,
                magnetic_left=[[str(Wl[i,j]) for j in range(9)] for i in range(9)],
                magnetic_right=[[str(Wr[i,j]) for j in range(9)] for i in range(9)],spectrum=data)
    print('JSON_RESULTS='+json.dumps(output,sort_keys=True))
    print('ALL COMPLETE MINIMAL-LINK SU(3) GAUSS AND WILSON TESTS PASSED')
    return output

if __name__=='__main__':run()
