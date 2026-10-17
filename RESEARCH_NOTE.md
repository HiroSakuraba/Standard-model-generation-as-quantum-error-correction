# Second-generation matter as a quantum error-correcting code

**Independent exploratory calculation, 17 October 2026.**  Starting point: Shiroman Prakash and K. L. Ramya, *One generation of the Standard Model as a quantum error-correcting code*, arXiv:2610.21879v1 (15 October 2026), especially Theorem 1, Table III, the paragraph "Other representations," and Supplement Part C. Original authors' code: https://github.com/spdei/standard-model-code.

## Executive finding

The isolated second generation (charm, strange, muon and muon neutrino) gives exactly the **same** `[[15,1,3]]_p` triorthogonal code as the first, for all primes (with the paper's separate binary construction at p=2). The masses do not change SO(10) weights or gauge anomalies. A nontrivial *two-generation*, *two-logical-qudit* construction is possible when new, **generation-sensitive**, anomaly-free charge rows are included:

- For `p=3`, an eight-row charge lattice, puncturing both conjugate right-handed neutrinos, yields **`[[30,2,3]]_3`** with six independent X stabilizers. Its X distance is 7.
- For every odd prime `p >= 5`, a seven-row charge lattice yields **`[[30,2,3]]_p`** with five independent X stabilizers. At `p=5` its X distance is 9.
- These constructions are *not merely the block-diagonal direct sum as charge matrices*: they retain all five common SO(10) Cartan charge rows but use only 3 or 2 extra differential rows, respectively, rather than five. No claim of new code inequivalence under arbitrary Clifford transformations has been made.
- By contrast, adding only a uniform family-difference U(1) charge gives `[[30,2,2]]_p` for odd `p`: matched fermions in the two families cause weight-two logical Z errors. Using only the identical, flavor-blind SO(10) rows gives `[[30,1,3]]_p` (for odd primes), **not** two logical qudits.

All statements above can be reproduced by running `python3 matter_qec_study.py`. The script checks exact integer cubic anomaly cancellation, finite-field ranks, Z-logical supports up to weight three, and exhaustive X distances for the specified p=3/5 coupled codes.

## Precise construction

Each SO(10) generation supplies 16 left-handed Weyl fields labelled by the even-parity sign weights

`Q = {q in {+1,-1}^5 : the number of negative signs is even}`.

The conjugate right-handed neutrino has weight `q0 = (1,1,1,1,1)`. Let `q` run over all 16 weights in each family, with generation index `s = +1` (electron family) or `s = -1` (muon family). Define an integral charge column

`h(q,s) = (q_1,...,q_5, s (B q)_1, ..., s (B q)_m)`.

Here `B` is an `m x 5` integer matrix. The full charge matrix has `32` columns and `5+m` rows. Puncture the two coordinates `(q0,+1)` and `(q0,-1)` associated with `nu_e^c` and `nu_mu^c`. The remaining 30 fields are the formal physical qudits. The rank-two projection onto those punctures supplies two logical qudits.

The cubic anomaly vanishes **over the integers** for every such `B`, since the full contribution from each of the two 16-weight families vanishes separately:

`sum_{s=+/-} sum_{q in Q} (a.q + s b.(Bq))^3 = 0`, for all integer coefficient vectors `a,b`.

The logical phase from a transversal cubic phase gate for p>=5 is

`omega^{-(a^3+b^3)}` on a logical basis state `|a,b>`.

The analogous qutrit 9th-root gate is exact by the same integer-lift argument given in the paper. This is a **separable** logical phase, not an entangling two-logical-qudit gate. Binary p=2 is not included in this odd-prime construction: the sign-charge lattice collapses there and requires the paper's different Reed–Muller construction.

### How to retain distance three

Let `b0 = B q0` mod p. The two punctured columns have common components `q0` and opposite differential components `+b0` and `-b0`. They are independent when `b0 != 0` and p is odd.

Within the 15 other SO(10) weights, `q0` lies in the span of no one or two distinct non-singlet weights over any odd prime, as proved in the paper. Thus any potentially nontrivial weight-two Z logical in the two-family construction must use *the same weight in both families*, with opposite coefficients. Such a pair has charge `(0, 2 Bq)`. It is logical exactly when `Bq` is a **nonzero scalar multiple** of `Bq0` modulo p.

Consequently, distance >=3 is guaranteed if every non-singlet weight has `Bq=0` or `Bq` linearly independent from `Bq0`. A three-field logical is always present inside each family because the original one-generation 15-qudit code has distance three; hence the distance is exactly three when this condition holds. This establishes Z distance three; the separately computed X distances establish full CSS distance three in the tested cases. For all odd p, the usual bound for the triorthogonal construction, `d_Z <= d_X`, also applies.

**Qutrit choice (p=3)**

```
B_3 = [ 1  0  0  0  0 ]
      [ 0  1  0  0  1 ]
      [ 0  0  1  1  0 ]
```

Here `B_3 q0 = (1,2,2)` modulo 3. If `B_3 q` is a nonzero multiple of `B_3 q0`, its first component forces the multiplier to be `q1 = +/-1`, and the other two components force `q2 = q5 = q1` and `q3 = q4 = q1`. Thus `q` is all plus or all minus. All plus is the punctured singlet; all minus has odd parity and is absent. The code has `[[30,2,3]]_3`, six independent X stabilizers, 22 independent Z stabilizers, and X distance seven.

**Prime dimensions p>=5**

```
B_>=5 = [ 1  0  0  0  0 ]
        [ 0  1  1  1  1 ]
```

Here `B q0 = (1,4)`. If `B q` is proportional to `(1,4)`, then `q2+q3+q4+q5 = 4 q1` modulo p. Write `k` for the number of minus signs among coordinates 2,...,5, so the left side is `4-2k` where `0 <= k <= 4`. If `q1=+1`, congruence requires `p | 2k`, and thus k=0. If `q1=-1`, it requires `p | 2(4-k)`, and thus k=4. Either yields the all-plus singlet or the excluded odd-parity all-minus vector. Thus the condition is satisfied for *all* primes p>=5. This produces `[[30,2,3]]_p`, five independent X stabilizers and 23 independent Z stabilizers. Its X distance at p=5 is nine.

The above proofs also exclude weight-one nontrivial errors: a non-singlet sign weight is never proportional to the removed all-plus weight for odd p.

### Minimal differential-row count in this specific architecture

With one differential Cartan row, there are always non-singlet weights with nonzero projection, and all such projections are proportional, giving Z-distance two. For p=3 and two differential rows, a complete enumeration of the `1080` two-dimensional row spaces in F_3^5 that resolve the two punctures shows every space retains at least one weight-two logical pair. The three-row qutrit choice is therefore minimal **within the architecture consisting of all five common SO(10) rows plus differential Cartan rows**. For p>=5, one row fails and the explicit two-row construction succeeds, making two minimal within the same architecture. This is not a universal minimality statement for arbitrary quantum error-correcting codes or arbitrary anomaly-free charge assignments.

## Exact computation summary

| Construction | prime | n | logical k | X stabilizers | Z stabilizers | Z dist. | X dist. | weight-3 Z-logical supports |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| One muon/charm/strange family | 3 | 15 | 1 | 4 | 10 | 3 | >=3 | 25 |
| Two families, common Cartan only | 5 | 30 | 1 | 4 | 25 | 3 | >=3 | 200 |
| Two families, common Cartan + uniform family U(1) | 5 | 30 | 2 | 4 | 24 | 2 | -- | -- |
| Two families, `B_3` | 3 | 30 | 2 | 6 | 22 | 3 | **7** | **50** |
| Two families, `B_>=5` | 5 | 30 | 2 | 5 | 23 | 3 | **9** | **50** |
| Two families, `B_>=5` | 7 | 30 | 2 | 5 | 23 | 3 | >=3 | **50** |
| Direct sum of two one-family codes | 3 | 30 | 2 | 8 | 20 | 3 | >=3 | **50** |

For each reduced-rank, distance-three two-family code, the 50 weight-three logical supports split into **25 exclusively within each family**; no mixed-family three-field logical supports remain. Thus there is **no demonstrated improvement in the leading cubic term of magic-state distillation error** over two independent 15-qudit codes. In the qutrit setting with the same symmetric independent Z-error model as the paper, the leading infidelity per output remains `6.25 epsilon^3`. The acceptance rates and higher-order error terms could differ; these have not been evaluated here.

As a separate qubit cross-check, the second generation gives the known 15-qubit Reed–Muller construction at p=2, with Z distance 3 and X distance 7. The binary construction has 35 weight-three logical supports, unlike the odd-prime count of 25; this is consistent with the paper's `35 epsilon^3` qubit result.

## Physical interpretation and limitations

1. **This is an algebraic correspondence, not evidence for particles performing error correction.** A fermion species labels a *formal code coordinate*. To turn this into hardware would require a physical Hilbert-space mapping, controllable stabilizer measurements, encoding/decoding, and noise characterization, none of which follows from anomaly cancellation.
2. **Flavor charges are extra assumptions.** The common five charges originate from SO(10) Cartan directions. The generation-sensitive rows `s Bq` are mathematical flavor-dependent U(1) charge assignments, not gauge symmetries of the observed Standard Model. If exact, they constrain intergeneration Yukawa couplings and are generally incompatible with unrestricted observed CKM/PMNS flavor mixing without symmetry breaking or other new dynamics.
3. **The extra charge rows are anomaly free, but this does not define a realistic theory.** Abelian anomaly cancellation is a necessary consistency check, not a phenomenological model. Full model building would need symmetry breaking, gauge bosons, neutrino masses, and constraints from flavor-changing interactions.
4. **Novelty is not established.** The p=3 and p>=5 reduced-row constructions were derived and computationally checked here, but have not undergone a literature survey for equivalence to existing generalized triorthogonal codes or independent peer review. The authors' supplemental discussion already warns that simply duplicating identical generations fails to produce useful distance; our construction adds distinguishing charges rather than merely duplicating the exact same charge column.

## Next tests worth doing

- Compute full weight enumerators and **exact acceptance/error polynomials** of the coupled codes, especially whether fewer X-stabilizer checks help or hurt distillation when realistic noisy measurements are included.
- Classify these `[[30,2,3]]_p` codes under monomial and local-Clifford equivalence against the direct sum, and evaluate possible correlated logical error channels at weight four and beyond.
- Model flavor breaking and CKM/PMNS rotations explicitly. Ask whether approximate flavor constraints can support an approximate code; with ordinary Standard Model mass mixing, the new flavor-dependent charges are not exact symmetries.
- Search alternative integral charge assignments satisfying cubic *and mixed* anomaly constraints that might offer distance four or a transversal **entangling** logical phase without simply multiplying generations.

## Reproduction

Run `python3 matter_qec_study.py` from this directory. No packages or network access are needed. `test_results.txt` records a completed run. Code tests exact integer triorthogonality in the odd-prime lattices, not a numerical approximation. Full Z distance at p>=5 follows from the analytic argument plus rank tests for p=5 and p=7; it is not based on brute-force enumeration for infinitely many primes. The p=3 two-row minimality test enumerates every eligible rank-two row space over F3 in RREF form.
