"""Strictly Gauss-preserving two-register logical-syndrome obstruction.

Follows colour_s3_two_register_syndrome.py (PR #49).
For neutral ancillas, any one/two-register gate commuting with BOTH local
S3 Gauss actions preserves the irreducible Peter-Weyl label r of the
single six-level link. Why:
  - without the link, the link label is untouched;
  - with link + matter/neutral pointer at one vertex, invariance under
    the OTHER vertex's independent action implies commutation with
    left/right-regular central irrep projectors on the link.

In the Gauss-invariant 14D physical space the link-irrep dimensions are
4 (trivial), 1 (sign), and 9 (standard). The dressed exchange S mixes
these sectors. The healthy four-dimensional phase-code K and errored SK
have NONZERO blockwise overlaps K^dag P_r SK. Thus no instrument formed
from strictly invariant <=2-register gates + neutral pointers can
perfectly and non-destructively report whether S occurred.

Positive test: a *single* gauge-invariant two-register pointer–link
unitary can read out r into one qutrit pointer. This is not the logical
error syndrome: K and SK have the SAME r probabilities.

Independent algebra check: all Gauss-twirled physical-pair-supported
operators compress to a 29D linear span, whose associative closure on
the 14D invariant space is dimension 98 = 4^2 + 1^2 + 9^2.
These are numerical rank checks accompanying the algebraic sector
preservation theorem, not exact arithmetic certificates.

NOT a no-go for three-register gauge-invariant primitives, gauge-charged
ancillas, temporary gauge violation, nonlocal link operations, or
continuous-SU3 physical circuits.
"""
from __future__ import annotations

from itertools import combinations
from pathlib import Path
import json

import numpy as np

from colour_gauge_link_s3 import (
    DIM, G, MATTER, PERMS, REP, compose, inverse, transform, projector,
    physical_isometry, dressed_exchange,
)
from colour_s3_edge_logical_syndrome import code_isometry

DIMS = (6, 3, 3, 3, 3)     # link, A0, A1, B0, B1
IRREPS = ("trivial", "sign", "standard")
IRREP_DIMS = (1, 1, 2)
LOGICAL = 4


def parity(g):
    p = PERMS[g]
    return (-1)**sum(p[i] > p[j] for i in range(3) for j in range(i+1, 3))


def link_central_projectors():
    """Peter-Weyl sector P_r = (d_r/|G|) sum_h chi_r(h^-1) L_h.

    S3 irreps have characters 1, sign, standard. Signed permutation D
    decomposes as sign + standard, giving chi_standard=Tr(D)-sign.
    """
    trivial = np.ones(G)
    sign = np.array([parity(h) for h in range(G)])
    standard = np.array([np.trace(REP[h]).real - sign[h]
                         for h in range(G)])
    chars = (trivial, sign, standard)
    out = {}
    for name, d, chi in zip(IRREPS, IRREP_DIMS, chars):
        link_p = np.zeros((G,G),complex)
        for h in range(G):
            weight = d * chi[inverse(h)] / G
            for g in range(G):
                link_p[compose(h,g),g] += weight
        assert np.linalg.norm(link_p@link_p-link_p) < 1e-12
        assert np.linalg.norm(link_p-link_p.conj().T) < 1e-12
        assert abs(np.trace(link_p).real-d*d) < 1e-12
        out[name] = link_p
    assert np.linalg.norm(sum(out.values())-np.eye(G)) < 1e-12
    return out


def apply_link(y, op):
    """Apply a 6x6 link matrix on a 486xR input; no dense 486x486."""
    x = np.asarray(y,complex)
    vector=x.ndim==1
    if vector:x=x[:,None]
    assert x.shape[0]==DIM and op.shape==(G,G)
    new=np.einsum("hg,gmr->hmr",op,x.reshape((G,MATTER,-1)))
    new=new.reshape(x.shape)
    return new[:,0] if vector else new


def pointer_link_irrep_gate():
    """One 18x18 strictly gauge-invariant pair unitary:
       V = sum_r X_pointer^r tensor P_r(link), r=0,1,2.
    Every constituent commutes with both vertex S3 actions.
    """
    p=link_central_projectors()
    out=np.zeros((3*G,3*G),complex)
    for r,name in enumerate(IRREPS):
        shift=np.roll(np.eye(3,dtype=complex),r,axis=0)
        out+=np.kron(shift,p[name])
    assert np.linalg.norm(out.conj().T@out-np.eye(18)) < 1e-12
    return out


def irrep_pointer_kraus(y):
    """Ancilla |0> then a SINGLE 18x18 pointer-link pair gate."""
    x=np.asarray(y,complex)
    vector=x.ndim==1
    if vector:x=x[:,None]
    assert x.shape[0]==DIM
    z=np.zeros((3,G,MATTER,x.shape[1]),complex)
    z[0]=x.reshape(G,MATTER,-1)
    gate=pointer_link_irrep_gate()
    output=(gate@z.reshape((18,-1))).reshape((3,DIM,x.shape[1]))
    if vector:return output[:,:,0]
    return output


def compressed_pair_matrix_units(W,i,j):
    """Physical compression of all E_ab supported on pair (i,j).

    W invariant => W^dag E W = W^dag (Gauss-twirl(E)) W.
    Gauss twirling a two-register operator leaves its support unchanged,
    so every matrix in this span has a strictly gauge-invariant
    pair-supported representative. No arbitrary noninvariant gate is
    used in the conceptual algebra generation.
    """
    assert i!=j and 0<=i<5 and 0<=j<5
    shape=DIMS+(W.shape[1],)
    a=np.moveaxis(W.reshape(shape),(i,j),(0,1))
    d=DIMS[i]*DIMS[j]
    a=a.reshape(d,-1,W.shape[1])
    return np.einsum("ark,brl->abkl",a.conj(),a).reshape((-1,14*14))


def generated_pair_algebra(W):
    """Numerical associative algebra rank; cost small (14x14 matrices)."""
    pair_stacks=[compressed_pair_matrix_units(W,i,j)
                 for i,j in combinations(range(5),2)]
    allpairs=np.concatenate(pair_stacks)
    _,sv,Vh=np.linalg.svd(allpairs,full_matrices=False)
    first=int(np.sum(sv>1e-9))
    assert first==29,first
    generators=Vh[:first].reshape(first,14,14)
    current=generators.copy()
    layers=[first]
    for _ in range(2):
        products=np.einsum("aij,bjk->abik",current,generators).reshape(-1,196)
        stack=np.concatenate((current.reshape(-1,196),products),axis=0)
        _,sing,Vh2=np.linalg.svd(stack,full_matrices=False)
        rank=int(np.sum(sing>1e-8))
        current=Vh2[:rank].reshape(rank,14,14)
        layers.append(rank)
    assert layers==[29,98,98],layers
    return {"pair_linear_span_dimension":first,
            "associative_closure_ranks":layers,
            "associative_algebra_dimension":layers[-1],
            "sector_block_matrix_dimension_sum":
                sum(n*n for n in (4,1,9))}


def main_result(include_algebra=True):
    W=physical_isometry()
    K,F,_=code_isometry(W)
    p=link_central_projectors()
    physical_ranks={}
    k_weights={}
    f_weights={}
    block_overlap={}
    swap_sector_commutators={}
    for name,P in p.items():
        # These link projectors commute with both *full* independent
        # Gauss actions and are physical, not arbitrary gauge projectors.
        QW=apply_link(W,P)
        for v in ("A","B"):
            for h in range(G):
                assert np.linalg.norm(transform(QW,h,v)-QW)<2e-10
        physical_ranks[name]=float(np.trace(W.conj().T@QW).real)
        k_weights[name]=float(np.linalg.norm(apply_link(K,P))**2/LOGICAL)
        f_weights[name]=float(np.linalg.norm(apply_link(F,P))**2/LOGICAL)
        block_overlap[name]=float(np.linalg.norm(K.conj().T@apply_link(F,P)))
        swap_sector_commutators[name]=float(np.linalg.norm(
            apply_link(dressed_exchange(W),P)
            -dressed_exchange(apply_link(W,P))))
    assert [round(physical_ranks[x]) for x in IRREPS]==[4,1,9]
    assert max(abs(k_weights[k]-f_weights[k]) for k in IRREPS)<1e-12
    assert all(block_overlap[k]>.4 for k in IRREPS)
    assert all(swap_sector_commutators[k]>.5 for k in IRREPS)

    rng=np.random.default_rng(1927)
    random_y=rng.normal(size=(DIM,2))+1j*rng.normal(size=(DIM,2))
    branches=irrep_pointer_kraus(random_y)
    assert abs(np.vdot(branches,branches)-np.vdot(random_y,random_y))<1e-10
    want=np.stack([apply_link(random_y,p[n]) for n in IRREPS])
    assert np.linalg.norm(branches-want)<1e-10
    pointer_gauge_errors={}
    for v in ("A","B"):
        worst=0.
        for h in range(G):
            for branch in branches:
                # For arbitrary input the *covariance* of the pointer
                # instrument matters; physical inputs remain physical.
                test=irrep_pointer_kraus(transform(random_y,h,v))
                ref=np.stack([transform(z,h,v) for z in branches])
                worst=max(worst,float(np.linalg.norm(test-ref)))
        pointer_gauge_errors[v]=worst
        assert worst<1e-10

    healthy=irrep_pointer_kraus(K)
    faulted=irrep_pointer_kraus(F)
    good_probs=[float(np.linalg.norm(z)**2/LOGICAL) for z in healthy]
    bad_probs=[float(np.linalg.norm(z)**2/LOGICAL) for z in faulted]
    assert np.linalg.norm(np.array(good_probs)-np.array(bad_probs))<1e-12

    result={
        "system":"two-vertex finite S3 subgroup gauge model, one 6-state link",
        "physical_gauss_dimension":14,
        "logical_phase_code_dimension":4,
        "link_irreps":list(IRREPS),
        "link_irrep_dimensions":[1,1,2],
        "physical_sector_ranks":physical_ranks,
        "healthy_link_sector_probabilities":k_weights,
        "dressed_swap_link_sector_probabilities":f_weights,
        "healthy_vs_dressed_syndrome_block_overlap_fro":block_overlap,
        "dressed_swap_link_sector_commutator_fro":swap_sector_commutators,
        "one_pointer_qutrit_irrep_readout":{
            "two_register_unitary_dimension":18,
            "one_arbitrary_pointer_link_pair_unitary":1,
            "healthy_pointer_probabilities":good_probs,
            "dressed_swap_pointer_probabilities":bad_probs,
            "measurement_kraus_completeness_residual":
                float(abs(np.vdot(branches,branches)-np.vdot(random_y,random_y))),
            "gauge_covariance_residuals":pointer_gauge_errors,
            "resolves_known_dressed_swap":False,
            "native_gate_compilation":"not provided"},
        "no_go":{
            "conditions":[
                "all intermediate gates commute with BOTH vertex Gauss representations",
                "at most two physical registers per elementary gate, pointer included",
                "ancilla pointers transform trivially under the gauge group",
                "no nonlocal joint link operation or gauge-charged mediator"],
            "exact_reason":"Each such gate conserves the link Peter-Weyl irrep; nonzero K†P_r S K for each sector excludes a deterministic nondemolition code-vs-error syndrome.",
            "does_not_rule_out":[
                "gauge-invariant operations on three or more registers",
                "gauge-charged ancillas spanning the edge",
                "temporary Gauss violation as in prior PR49 decomposition",
                "other encodings or approximate/partial measurements"]},
        "scope":"Finite S3 subgroup and new 4D eigenphase code, not full SU3 or original 21D nine-quark code."}
    if include_algebra:
        result["pair_operator_algebra"]=generated_pair_algebra(W)
    return result


def main():
    result=main_result()
    Path("colour_s3_strict_gauss_pair_nogo_results.json").write_text(
        json.dumps(result,indent=2,sort_keys=True)+"\n")
    print("STRICT_GAUSS_TWO_REGISTER_RESULT",json.dumps({
        "physical_sector_ranks":result["physical_sector_ranks"],
        "block_overlap_fro":result["healthy_vs_dressed_syndrome_block_overlap_fro"],
        "algebra":result["pair_operator_algebra"],
        "healthy_link_sector_probabilities":result["healthy_link_sector_probabilities"],
        "errored_link_sector_probabilities":result["dressed_swap_link_sector_probabilities"],
        "pointer_covariance":result["one_pointer_qutrit_irrep_readout"]["gauge_covariance_residuals"]
    },sort_keys=True),flush=True)


if __name__=="__main__":
    main()
