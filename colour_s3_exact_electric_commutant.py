"""Exact S3 electric+S4 matter-swap commutant and a restricted active code.

Replaces the FLOATING-POINT-only 6D commutant finding in merged PR60
with a rational certificate on the 14 physical Gauss states.

Four signed-permutation colour fundamentals have (sign(g))**4=1.
The globally colour-invariant, gauge-fixed matter space therefore has
one basis orbit-sum for each equality pattern of four labels under
simultaneous S3 colour relabelling: EXACTLY 14 patterns.

All 24 physical matter permutations act as integer 14x14 permutation
matrices on those patterns. Their commutant has one integer basis
matrix for each diagonal S4-orbit on ordered pairs of patterns:
EXACTLY 24 indicator matrices (disjoint support).

The standard-irrep link electric projector reduces in gauge-fixed
coordinates to the RATIONAL 14x14 matrix
  E = (1/3) sum_{h in S3} chi_standard(h) D_B(h)^tensor2.
This is a link-only gauge-invariant electric Hamiltonian, not a
microscopic native implementation. Requiring a linear combination of
the 24 S4-centralizer indicator matrices to commute with E gives a
196x24 RATIONAL matrix of EXACT rank 18. Thus the joint complex
commutant has exactly dimension 24-18 = SIX. Four logical noiseless
subsystem dimensions would require an M4(C) commutant (dimension 16),
which is impossible in this precise finite S3 error algebra.

Independent full 784x196 commutator matrices have rank 190 modulo
101 and 65521; in characteristic zero these ranks give a matching
lower bound 190, agreeing with the rational 24-18 derivation.

A further restricted example: for errors {I, S_A0B0, R_electric}, one
can construct a 4D KL-perfect phase code with S acting trivially and
R mapping the code to an orthogonal syndrome space. HOWEVER another
4D gauge-invariant DFS exists in the five-dimensional joint
eigenspace S=+1, E=+1, so active recovery is unnecessary against this
restricted set. No four-dimensional active code for the FULL S4+E
error set is constructed or ruled out by the commutant alone.

No real SU3 gauge-link device, nine-quark 21D code or 65 syndrome.
"""
from __future__ import annotations

from itertools import permutations, product
from pathlib import Path
import json

import numpy as np
import sympy as sp

from colour_gauge_link_s3 import (
    DIM, G, MATTER, REP, four_rep, inverse, physical_isometry,
    dressed_exchange,
)
from colour_s3_strict_gauss_pair_nogo import (
    apply_link, link_central_projectors,
)

NCOLOUR=3
NSITES=4
DGAUSS=14
DSYS=4


def four_colour_orbits():
    """Colour-orbit sums B_j are INTEGER bases of the Gauss-fixed space.

    The global 4-fold signed permutation is an unsigned simultaneous
    relabelling because the sign occurs four times.
    """
    relabels=tuple(permutations(range(NCOLOUR)))
    words=tuple(product(range(NCOLOUR),repeat=NSITES))
    orbit=lambda x:tuple(sorted({
        tuple(h[c] for c in x) for h in relabels}))
    classes=sorted({orbit(x) for x in words})
    assert len(classes)==DGAUSS
    lookup={x:i for i,items in enumerate(classes) for x in items}
    assert len(lookup)==len(words)==81
    return classes,lookup,relabels


def exact_generators():
    """Integer S4 positional action and rational gauge-fixed electric E."""
    classes,lookup,relabels=four_colour_orbits()
    def permute_sites(t,pi):
        return tuple(t[pi[i]] for i in range(NSITES))
    def transposition(i,j):
        pi=list(range(NSITES))
        pi[i],pi[j]=pi[j],pi[i]
        M=sp.zeros(DGAUSS)
        for k,cl in enumerate(classes):
            dest={lookup[permute_sites(t,pi)] for t in cl}
            assert len(dest)==1
            M[dest.pop(),k]=1
        assert M*M==sp.eye(DGAUSS)
        return M
    A=transposition(0,1)
    B=transposition(2,3)
    S=transposition(0,2)
    # E is single-link standard-irrep central projector. In link=g=e
    # gauge it is the S3 character average acting on the two B colours.
    # sign(h)^2=+1 for two B fundamental matter qutrits.
    E=sp.zeros(DGAUSS)
    for k,cl in enumerate(classes):
        for t in cl:
            for h in relabels:
                chi=2 if h==(0,1,2) else -1 if h in (
                    (1,2,0),(2,0,1)) else 0
                if not chi:continue
                target=(t[0],t[1],h[t[2]],h[t[3]])
                i=lookup[target]
                E[i,k]+=sp.Rational(chi,3*len(classes[i]))
    assert E*E==E
    assert sp.trace(E)==9
    return classes,lookup,(A,B,S,E)


def positional_s4(classes,lookup):
    """All 24 permutations of FOUR matter site positions."""
    group=[]
    for pi in permutations(range(NSITES)):
        action=tuple(lookup[tuple(cl[0][pi[j]]
                                  for j in range(NSITES))]
                     for cl in classes)
        assert sorted(action)==list(range(DGAUSS))
        group.append(action)
    assert len(set(group))==24
    return group


def invariant_orbital_matrices(group):
    """EXACT basis of centralizer of the S4 14-point permutation action.

    A matrix commutes with all site permutations iff its entries are
    constant along simultaneous S4-orbits on ordered index pairs.
    Each indicator matrix is integral and supports a DISJOINT orbit,
    so all returned matrices are linearly independent automatically.
    """
    occupied=set()
    basis=[]
    for i in range(DGAUSS):
        for j in range(DGAUSS):
            if (i,j) in occupied:continue
            orbit={(g[i],g[j]) for g in group}
            assert not (orbit&occupied)
            occupied.update(orbit)
            X=sp.zeros(DGAUSS)
            for a,b in orbit:X[a,b]=1
            basis.append(X)
    assert len(occupied)==DGAUSS**2
    assert len(basis)==24
    return basis


def electric_centralizer_constraint(basis,E):
    """196x24 EXACT rational linear map X -> XE-EX."""
    D=sp.zeros(DGAUSS**2,len(basis))
    for col,X in enumerate(basis):
        comm=X*E-E*X
        for i in range(DGAUSS):
            for j in range(DGAUSS):
                D[i*DGAUSS+j,col]=comm[i,j]
    return D


def modular_rank(M,p):
    """Exact prime-field row elimination; no floating-point operations."""
    m=np.array(M,dtype=np.int64)%p
    row=0
    for col in range(m.shape[1]):
        nonzero=np.flatnonzero(m[row:,col])
        if len(nonzero)==0:continue
        k=row+int(nonzero[0])
        m[[row,k]]=m[[k,row]]
        inv=pow(int(m[row,col]),-1,p)
        m[row,col:]=m[row,col:]*inv%p
        factors=m[row+1:,col].copy()
        m[row+1:,col:]=(m[row+1:,col:]-
                        factors[:,None]*m[row,col:])%p
        row+=1
        if row==m.shape[0]:break
    return row


def full_commutator_ranks(ops):
    """Independent 784x196 certificate from rational operators.

    Multiply all operators by 3 so integer modular reduction is exact.
    Over Q, any nonzero minor mod p implies a nonzero integer minor,
    hence the rational rank is >= the mod-p rank.
    """
    I=np.eye(DGAUSS,dtype=np.int64)
    integer_ops=[]
    for O in ops:
        scaled=3*O
        assert all(v.q==1 for v in scaled)
        integer_ops.append(np.array(scaled,dtype=np.int64))
    stacked=np.concatenate([np.kron(O.T,I)-np.kron(I,O)
                            for O in integer_ops])
    assert stacked.shape==(784,196)
    return {str(p):modular_rank(stacked,p) for p in (101,65521)}


def physical_orbit_isometry(classes):
    """Independent 486D test; map orbit sums to actual two-vertex Gauss W."""
    V=np.zeros((MATTER,DGAUSS),complex)
    sizes=np.array([len(cl) for cl in classes],dtype=float)
    for j,cl in enumerate(classes):
        for t in cl:
            V[np.ravel_multi_index(t,(3,3,3,3)),j]=1/np.sqrt(len(cl))
    W=np.zeros((DIM,DGAUSS),complex)
    for g in range(G):
        W[g*MATTER:(g+1)*MATTER]=four_rep(inverse(g),(2,3))@V/np.sqrt(G)
    assert np.linalg.norm(W.conj().T@W-np.eye(DGAUSS))<2e-12
    return W,sizes


def crosscheck_existing_physical_model(classes,ops):
    """Both electrically physical models agree, not just their ranks."""
    W,n=physical_orbit_isometry(classes)
    A,B,S,E=ops
    diag=np.sqrt(n)
    def unit(unnormalized):
        return np.array(unnormalized.tolist(),dtype=float)*(
            diag[:,None]/diag[None,:])
    from colour_s3_electric_error_algebra import (
        physical_onsite_swap,
    )
    full=(
        W.conj().T@physical_onsite_swap(W,0,1),
        W.conj().T@physical_onsite_swap(W,2,3),
        W.conj().T@dressed_exchange(W),
        W.conj().T@apply_link(W,link_central_projectors()["standard"]),
    )
    errors=[float(np.linalg.norm(ref-unit(x)))
            for ref,x in zip(full,ops)]
    assert max(errors)<2e-11
    return {"compressed_orbit_vs_486D_residuals":errors,
            "max_residual":max(errors),
            "orbit_basis_normalizations":n.astype(int).tolist()}


def restricted_swap_electric_code(ops,sizes):
    """Restricted exact-KL CONSTRUCTION; also finds a larger passive DFS.

    {I,S,R=I-2E} is a genuine proper subset of the full exchange+E
    algebra. S has +1 dimension 10. On this 10D subspace, compressed
    E has eigenvalues 0,0,1/2,1/2,1/2,1,1,1,1,1.

    Two 0/1 eigenvector superpositions and two 1/2 eigenvectors give
    a four-dimensional C satisfying S C=C, C†R C=0, R C ⟂ C.
    Syndrome projection onto C versus RC, then R feedback, perfectly
    recovers {I,R} in an IDEAL NON-NATIVE global measurement.

    BUT the simultaneous S=+1,E=+1 sector has dimension FIVE;
    therefore an equally large passive code also avoids the entire
    restricted {I,S,R} set. No active advantage is established.
    """
    _,_,S,E=ops
    n=np.sqrt(np.array(sizes,dtype=float))
    Sm=np.array(S.tolist(),dtype=float)*(n[:,None]/n[None,:])
    Em=np.array(E.tolist(),dtype=float)*(n[:,None]/n[None,:])
    assert np.linalg.norm(Sm-Sm.T)<1e-12
    assert np.linalg.norm(Em-Em.T)<1e-12
    vals,vecs=np.linalg.eigh(Sm)
    V=vecs[:,vals>.5]
    assert V.shape==(14,10)
    evals,Q=np.linalg.eigh(V.T@Em@V)
    low=V@Q[:,np.abs(evals)<1e-9]
    middle=V@Q[:,np.abs(evals-.5)<1e-9]
    high=V@Q[:,np.abs(evals-1)<1e-9]
    assert (low.shape[1],middle.shape[1],high.shape[1])==(2,3,5)
    C=np.column_stack(((low[:,0]+high[:,0])/np.sqrt(2),
                       (low[:,1]+high[:,1])/np.sqrt(2),
                       middle[:,0],middle[:,1]))
    R=np.eye(DGAUSS)-2*Em
    residuals={
        "code_orthonormality":float(np.linalg.norm(C.T@C-np.eye(4))),
        "swap_acts_as_identity":float(np.linalg.norm(Sm@C-C)),
        "electric_fault_orthogonal_syndrome":float(np.linalg.norm(C.T@R@C)),
        "electric_fault_image_orthonormal":float(np.linalg.norm(
            (R@C).T@(R@C)-np.eye(4))),
        "passive_4D_joint_eigenspace_error":float(np.linalg.norm(
            Sm@high[:,:4]-high[:,:4])+np.linalg.norm(
            Em@high[:,:4]-high[:,:4])),
    }
    assert max(residuals.values())<2e-11
    return {
        "restricted_error_set":["I","single dressed swap S","electric phase R=I-2E"],
        "S_plus_space_dimension":10,
        "electric_compressed_on_S_plus_eigenspace":{
            "zero":2,"one_half":3,"one":5},
        "KL_correctable_active_code_dimension":4,
        "active_syndrome_space_rank":8,
        "exact_KL_residuals_numeric_construction":residuals,
        "passive_joint_S_plus_E_one_space_dimension":5,
        "passive_4D_code_also_exists":True,
        "ideal_syndrome_is_global_not_hardware_compiled":True,
    }


def exact_study():
    classes,lookup,ops=exact_generators()
    A,B,S,E=ops
    group=positional_s4(classes,lookup)
    basis=invariant_orbital_matrices(group)
    D=electric_centralizer_constraint(basis,E)
    assert D.shape==(196,24)
    assert all((3*x).q==1 for x in D)
    # Ranks are computed over the exact rational field by SymPy.
    rank=D.rank()
    nullity=len(D.nullspace())
    assert rank==18 and nullity==6
    for v in D.nullspace():
        assert D*v==sp.zeros(196,1)
    allr=full_commutator_ranks(ops)
    assert all(x==190 for x in allr.values())
    assert DGAUSS**2-190==nullity
    swaps=[A,B,S]
    for X in basis:
        assert all(M*X==X*M for M in swaps)
    P=sp.zeros(DGAUSS)
    for pi in group:
        Q=sp.zeros(DGAUSS)
        for j,i in enumerate(pi):Q[i,j]=1
        P+=Q/24
    assert P*P==P
    assert sp.trace(P)==4
    assert (P*E-E*P)!=sp.zeros(DGAUSS)
    I=sp.eye(DGAUSS)
    joint_e_one_nullity=14-sp.Matrix.vstack(S-I,E-I).rank()
    assert joint_e_one_nullity==5
    comparison=crosscheck_existing_physical_model(classes,ops)
    restricted=restricted_swap_electric_code(ops,[len(x) for x in classes])
    return {
        "exact_arithmetic_engine":"SymPy rational and integer matrices; no floating-point ranks in final certificate",
        "physical_Gauss_space_dimension":14,
        "gauge_fixed_global_S3_colour_orbit_count":len(classes),
        "colour_orbit_cardinalities":[len(x) for x in classes],
        "signed_fourfold_action_is_unsigned_colour_permutation":True,
        "matter_S4_order":len(group),
        "S4_fixed_common_subspace_rank_exact":int(sp.trace(P)),
        "integer_S4_centralizer_orbital_basis_dimension":len(basis),
        "electric_projector_exact_rank":E.rank(),
        "electric_projector_idempotent_rational":E*E==E,
        "reduced_electric_commutator_matrix_shape":[D.rows,D.cols],
        "reduced_commutator_exact_rational_rank":rank,
        "joint_commutant_exact_complex_dimension":nullity,
        "independent_full_modular_commutator_ranks":allr,
        "numerical_full_486D_crosscheck":comparison,
        "exact_no_four_dim_noiseless_subsystem_under_S4_plus_E":True,
        "reason":"A four-dimensional noiseless subsystem would require a full 16D matrix algebra M4(C) in the exact 6D commutant.",
        "restricted_I_S_R_active_and_passive_code":restricted,
        "claim_limits":[
            "This EXACT dimension-six commutant certificate replaces PR60's numerical-only claim for this precise finite S3 error algebra",
            "The S4 orbit-commutant uses colour orbit SUMS (not normalized amplitudes); rational E is similarity-conjugate to physical Hermitian projector, confirmed independently",
            "No 4D NOISELESS SUBSYSTEM for *full algebra* does not rule out a 4D actively correctable code for a smaller error SET",
            "The constructed restricted active 4D KL code for {I,S,R} has an IDEAL 8D syndrome, but a passive 4D DFS also exists in that restricted setting",
            "No claim that a 4D active code for full {I,A,B,S,R} exists or is impossible; a separate KL/packing analysis is needed",
            "No native Fourier pointer, link-electric syndrome, physical link error rates, finite-U recovery or overall device performance",
            "This is a finite signed S3 subgroup with one link, not full continuous SU3 or nine-quark 21D 65-syndrome correction"
        ]}


def main():
    r=exact_study()
    Path("colour_s3_exact_electric_commutant_results.json").write_text(
        json.dumps(r,indent=2,sort_keys=True)+"\n")
    print("EXACT_ELECTRIC_COMMUTANT",json.dumps({
        "orbit_count":r["gauge_fixed_global_S3_colour_orbit_count"],
        "S4_commutant":r["integer_S4_centralizer_orbital_basis_dimension"],
        "reduced_rational_rank":r["reduced_commutator_exact_rational_rank"],
        "joint_complex_commutant_dimension":r["joint_commutant_exact_complex_dimension"],
        "modular_ranks":r["independent_full_modular_commutator_ranks"],
        "restricted_code":r["restricted_I_S_R_active_and_passive_code"]
    },sort_keys=True),flush=True)


if __name__=="__main__":
    main()
