"""Nine-quark colour code: physically allowed exchange algebra and syndrome rank.

Representation-theoretic statements checked on independent 42D Young module:
* swaps on a connected coupling graph generate S9.
* its (3,3,3) singlet Specht representation is irreducible, hence Burnside's
  theorem gives a FULL M42(C) complex associative error algebra.
* every individual exchange anticommutes with J, hence is detectable, but
  two distinct exchanges (when unknown) fail simultaneous Knill-Laflamme.
* two exchange couplings can produce a non-scalar second-order logical
  return amplitude even though P_K H P_K is scalar.

For all 9 independently addressable fundamental-colour qutrits, the
single-site error span is {I, T_i^a}, 1+9*8=73 labeled errors. On the
21D code, the 73x73 KL Gram matrix is explicitly known from SU(3)
invariance and P_K Swap_ij P_K=0:
  C_00=1, C_0,ia=0,
  C_ia,jb = delta_ab * [1/6 (i=j), -1/48 (i!=j)].
It has rank 65: 1 singlet syndrome and 8 adjoint components times 8
independent relative-site modes. Thus an orthogonal syndrome-label
implementation requires 65 labels, >=4 qutrit pointers (3**4=81).
This is an INFORMATIONAL projective-syndrome requirement, NOT a general
physical gate-count lower bound or a proof that colour-octet labels can
be locally measured in a gauge invariant experiment.

Physical locality: a bare colour exchange S_ij commutes with the
DIAGONAL SU(3) action of two colocated quarks, not with INDEPENDENT
local SU(3) rotations at two separated sites. The separated-site
interaction requires link transport/dressing and dynamical gauge links.
"""
from __future__ import annotations

from itertools import combinations
from math import ceil, log
from pathlib import Path
import json
import numpy as np
from scipy.linalg import expm
from colour_sign_twist_theorem import (
    tableaux, adjacent_generators, signed_tableau_transpose,
    code_isometry, swaps,
)
from colour_singlet_code import gellmann

N=9
COLORS=8
DIM=42
LOGICAL=21

def allowed_graph_components(edges,vertices=N):
    """Connected transposition edges generate product of S_(component size)."""
    neighbors={i:set() for i in range(1,vertices+1)}
    for i,j in edges:
        if i==j or i not in neighbors or j not in neighbors:
            raise ValueError((i,j))
        neighbors[i].add(j)
        neighbors[j].add(i)
    seen=set()
    comps=[]
    for i in neighbors:
        if i not in seen:
            todo=[i]
            group=set()
            while todo:
                k=todo.pop()
                if k in group:continue
                group.add(k);seen.add(k);todo.extend(neighbors[k]-group)
            comps.append(tuple(sorted(group)))
    return sorted(comps,key=lambda x:x[0])

def exact_single_colour_gram():
    """73x73 real rational KL Gram (returned as float), proven analytically."""
    C=np.zeros((1+N*COLORS,1+N*COLORS),float)
    C[0,0]=1.
    for a in range(COLORS):
        for i in range(N):
            for j in range(N):
                C[1+i*COLORS+a,1+j*COLORS+a]=1/6 if i==j else -1/48
    eigen=np.linalg.eigvalsh(C)
    assert np.min(eigen)>-1e-12
    assert np.count_nonzero(eigen>1e-10)==65
    assert abs(eigen[-1]-1)<1e-12
    assert np.count_nonzero(np.isclose(eigen,3/16,atol=1e-12))==64
    assert np.count_nonzero(abs(eigen)<1e-12)==8
    # Eight missing site-mode syndromes = exact Gauss/global-singlet relations.
    for a in range(COLORS):
        v=np.zeros(1+N*COLORS)
        for i in range(N):v[1+i*COLORS+a]=1
        assert np.linalg.norm(C@v)<1e-12
    return C

def orthogonal_modes(C):
    """Construct 64 orthogonal relative-site generator modes without full 3^9."""
    modes=[]
    # Helmert orthonormal contrasts for nine quark sites, sum_i f_ui = 0.
    for k in range(1,N):
        f=np.zeros(N)
        f[:k]=1/np.sqrt(k*(k+1))
        f[k]=-k/np.sqrt(k*(k+1))
        assert abs(sum(f))<1e-12
        modes.append(f)
    hel=np.stack(modes)
    assert np.linalg.norm(hel@hel.T-np.eye(N-1))<1e-12
    # For each a, T_a(f) has KL norm sqrt(3/16). Normalize to unit.
    vectors=[]
    for a in range(COLORS):
        for f in modes:
            vec=np.zeros(1+N*COLORS)
            for i in range(N):vec[1+i*COLORS+a]=f[i]*np.sqrt(16/3)
            vectors.append(vec)
    B=np.array(vectors).T
    assert B.shape==(73,64)
    assert np.linalg.norm(B.T@C@B-np.eye(64))<2e-12
    assert np.linalg.norm(B.T@C[:,0])<1e-12
    return hel,B

def physical_gauge_commutator():
    """Exchange is SU3-diagonal invariant, not invariant at independent sites."""
    S=np.array([[1 if a//3==b%3 and a%3==b//3 else 0
                 for b in range(9)] for a in range(9)],complex)
    I=np.eye(3)
    gens=gellmann()
    independent=[]
    global_rot=[]
    for T in gens:
        left=np.kron(T,I)
        right=np.kron(I,T)
        independent.append(float(np.linalg.norm(S@left-left@S)))
        global_rot.append(float(np.linalg.norm(S@(left+right)-(left+right)@S)))
    assert min(independent)>1e-2
    assert max(global_rot)<1e-12
    return {'min_independent_vertex_gauss_commutator_fro':min(independent),
            'max_common_vertex_colour_commutator_fro':max(global_rot),
            'interpretation':'Bare swap is allowed for co-located colour fundamentals under one vertex gauge action; separated-site swap requires gauge-link dressing.'}

def exchange_obstructions():
    S=adjacent_generators(tableaux())
    J=signed_tableau_transpose(tableaux())
    K=code_isometry(J)
    sw=swaps(S)
    graph=[(i,i+1) for i in range(1,9)]
    connected=allowed_graph_components(graph)
    assert connected==[tuple(range(1,10))]
    # All single swaps individually have zero code compression.
    worst_individual=max(np.linalg.norm(K.conj().T@A@K) for A in sw.values())
    assert worst_individual<2e-12
    # But the error set {I,S12,S23} cannot be jointly corrected:
    # compressed E1^dag E2 contains logical information.
    xy=K.conj().T@(sw[(1,2)]@sw[(2,3)])@K
    xy_non_scalar=float(np.linalg.norm(xy-np.trace(xy)/LOGICAL*np.eye(LOGICAL)))
    assert xy_non_scalar>1
    # A Hamiltonian with two noncommuting native pair exchanges has
    # scalar first-order compression, but non-scalar order t^2 return.
    H=.5*(sw[(1,2)]+sw[(2,3)])
    reduced=K.conj().T@H@K
    assert np.linalg.norm(reduced)<2e-12
    second=K.conj().T@(H@H)@K
    second_non_scalar=float(np.linalg.norm(second-np.trace(second)/LOGICAL*np.eye(LOGICAL)))
    assert second_non_scalar>.1
    time=.5
    returned=K.conj().T@expm(-1j*time*H)@K
    return_non_scalar=float(np.linalg.norm(returned-np.trace(returned)/LOGICAL*np.eye(LOGICAL)))
    assert return_non_scalar>1e-3
    # Formal theorem: connected graph edges generate S9. The Specht
    # representation (3,3,3) is irreducible over C; Burnside implies
    # the complex associative algebra generated by these swaps is M42.
    # Do NOT confuse this with the real Lie algebra of controls.
    return {'connected_graph_components':[list(c) for c in connected],
            'connected_exchange_group':'S9',
            'complex_associative_exchange_algebra_dimension':DIM*DIM,
            'proof':'connected edge transpositions generate S9; irreducible Specht (3,3,3) plus Burnside',
            'maximum_single_exchange_KLK_fro':float(worst_individual),
            'two_unknown_swap_error_KL_non_scalar_fro':xy_non_scalar,
            'second_order_native_exchange_return_non_scalar_fro':second_non_scalar,
            'finite_time_half_native_exchange_return_non_scalar_fro':return_non_scalar,
            'passively_protected_from_all_allowed_exchanges':False,
            'note':'Associative algebra M42 is NOT a claim that all U(42) pulse schedules are efficiently physically controllable.'}

def study():
    gram=exact_single_colour_gram()
    hel,syndrome_basis=orthogonal_modes(gram)
    gauge=physical_gauge_commutator()
    exch=exchange_obstructions()
    rank=65
    pointer_qutrits=next(n for n in range(1,10) if 3**n>=rank)
    pointer_qubits=next(n for n in range(1,12) if 2**n>=rank)
    assert pointer_qutrits==4 and pointer_qubits==7
    return {'project':'nine-qutrit global colour singlet sign-twist, operational error-algebra cost',
            'gauge_locality':gauge,'exchange_error_algebra':exch,
            'single_site_error_gram':{
                'number_of_labeled_generators_including_identity':73,
                'rank_of_KL_gram':rank,
                'gram_eigenvalue_one_multiplicity':1,
                'gram_eigenvalue_three_sixteenths_multiplicity':64,
                'gram_zero_eigenvalue_multiplicity':8,
                'site_mode_rank_per_su3_generator':8,
                'orthogonal_relative_site_modes':int(hel.shape[0]),
                'orthonormalized_error_basis_residual_fro':float(
                    np.linalg.norm(syndrome_basis.T@gram@syndrome_basis-np.eye(64))),
                'code_dimension':LOGICAL,
                'orthogonal_error_images_dimension':rank*LOGICAL,
                'physical_nine_qutrit_hilbert_dimension':3**N,
                'minimal_orthogonal_syndrome_pointer_qutrits':pointer_qutrits,
                'minimal_orthogonal_syndrome_pointer_qubits':pointer_qubits,
                'warning':'Pointer dimension is for resolving the full abstract correctable generator-error span. It is not a universal lower bound on gates, measurements, or physically gauge-invariant readout.'},
            'decision':'No passive gauge-only protection under connected co-located colour exchange graph; active gauge-covariant syndrome extraction or restrictions on error generators are necessary.',
            'limitations':[
                'Full 3^9 spatial-colour channel recovery is NOT synthesized; Gram plus KL provides abstract syndrome isometries.',
                'Error span of single-site T_i^a is not itself gauge invariant; syndrome environment must carry adjoint charge to preserve covariance.',
                'Only simultaneous colour SU3 at a single vertex makes bare swap gauge invariant; distant vertices need parallel-transport links.',
                'No operation counts for nonlocal J projector or 65-label measurement; four qutrit pointers are an information dimension bound only.',
                'Logical second-order exchange return does not by itself imply a measured environmental fault rate or decoder performance.']}
def main():
    out=study()
    Path('colour_allowed_error_algebra_budget_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print('COLOUR_ALLOWED_ERROR_ALGEBRA',json.dumps(out,sort_keys=True))

if __name__=='__main__':main()
