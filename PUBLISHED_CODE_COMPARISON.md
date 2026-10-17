# Explicit published-code comparison: Prakash–Saha (2025)

**Research status, 17 October 2026:** reproducible exclusion of one published odd-prime triorthogonal *family*. This is a comparison with external generator matrices, not a general novelty claim.

## The reference and its precise construction

Shiroman Prakash and Tanay Saha, [*Low Overhead Qutrit Magic State Distillation*, Quantum **9**, 1768 (2025)](https://doi.org/10.22331/q-2025-06-12-1768), [arXiv:2403.06228](https://arxiv.org/abs/2403.06228), Eqs. (8)–(11), construct an unpunctured qutrit triorthogonal space of length 9m. Its 3m generator rows are a repeated row (0,1,2,0,1,2,...) and 3m−1 block rows, each carrying (1,1,1) in one three-site block and (2,2,2) in the final block.

Puncture the first coordinate of each of k distinct early blocks. This produces their explicit family

    [[9m-k, k, 2]]_3,  1 <= k <= 3m-2.

**Independent audit:** the new script constructs their matrix from the published equations, and verifies that at (m,k)=(2,4) the resulting **complete 6×14 matrix** is *entry-for-entry identical* to Eq. (11), with the published logical/stabilizer row separation. This is stronger than merely checking its reported parameters.

## Reproduced examples

| Published m | k | Physical n | Encoded k | Distance |
|---:|---:|---:|---:|---:|
| 2 | 4 | 14 | 4 | 2 |
| 3 | 7 | 20 | 7 | 2 |
| 4 | 6 | 30 | 6 | 2 |

For all three examples, the script:

1. Checks every cubic triorthogonality identity of the unpunctured matrix over F3, including repeated row indices.
2. Computes full matrix rank and X-stabilizer rank exactly over F3.
3. Confirms that all X-stabilizer rows are orthogonal to all generator rows.
4. Exhaustively excludes weight-one logical Z and finds a weight-two logical Z witness.
5. Excludes weight-one logical X, proving full quantum-code distance exactly two.

### The same-length comparison

The published choice m=4,k=6 and our coupled two-generation qutrit construction both use 30 physical qutrits, and each has **six independent X-stabilizer checks**. Despite those coincidences, the codes differ:

| Invariant | Prakash–Saha m=4,k=6 | Coupled matter B3 |
|---|---:|---:|
| p | 3 | 3 |
| physical n | 30 | 30 |
| independent X checks | 6 | 6 |
| logical k | **6** | **2** |
| distance | **2** | **3** |
| minimum nonzero X-stabilizer weight | 6 | 10 |

Their *complete* X-stabilizer weight enumerators (weight: count) are

    Prakash–Saha: 0:1, 6:30, 9:40, 12:90, 15:60, 18:22, 24:486
    Coupled B3:    0:1, 10:4, 11:4, 12:8, 16:52, 18:128,
                   19:48, 20:200, 21:40, 22:148, 23:48, 24:48

Each sums to 3^6=729. The different X-stabilizer enumerators exclude any CSS-X-preserving coordinate permutation/rescaling between these presentations; they are **not** invariants of unrestricted local-Clifford or arbitrary local-unitary equivalence. The different (n,k,d) already excludes *all* unitary code equivalence between this particular pair.

### Exclusion of the published family

The Prakash–Saha family has **d=2 for every allowed (m,k)**, while our B3 two- and three-generation examples have **d=3**. Equivalence of codes preserves their distance under local Clifford transformations and coordinate permutations; in this case a stronger exclusion follows from the absence of any matching (n,k):

- For [[30,2,3]]_3, a reference family match would require 9m=n+k=32, impossible.
- For [[45,3,3]]_3, a match would require 9m=n+k=48, impossible.

Thus *no member of this explicitly defined published family* has the same (n,k) as either matter construction. These statements do **not** exclude subcodes, shortenings, nonstandard puncturings, or other generalized Reed–Muller families outside the specified construction.

## Reproduction

    python -m pip install numpy
    python published_qutrit_family.py

Dependencies: Python standard library for the published reference generator, plus repository-native \`charge_lattice_tools.py\` (which imports NumPy) for the matter-code comparison. The repository's GitHub Actions workflow runs this script.

## What this does and does not establish

**Established:** a faithful generator-matrix reconstruction of a published qutrit family, exactly verified code parameters, an equal-length but inequivalent published example, and a family-wide (n,k) exclusion.

**Still open:** direct equivalence comparison against a previously published **indecomposable [[30,2,3]]_3, [[45,3,3]]_3, [[30,2,3]]_5, or [[45,3,3]]_5** generator matrix, should one be identified. The broader odd-prime punctured Reed–Muller families in [Saha–Prakash, *Physical Review A* **113**, 042405 (2026)](https://doi.org/10.1103/tzkz-9p2f) have not been exhaustively searched for matching parameters.

### Same-looking parameters in the binary literature

Shi, Lu, Kim and Solé, [*Triorthogonal codes and self-dual codes* (2024)](https://arxiv.org/abs/2408.09685), construct a binary triorthogonal example denoted **[[30,2,3]]**. This is a **qubit** code (p=2), not a qutrit code (p=3), and the authors use the third slot as a Z-distance in their triorthogonal-code tables. It is not a candidate for a p=3 local-unitary equivalence comparison. Identical numerical brackets without the subscript prime and a consistent distance convention are insufficient for code identification.

These are formal stabilizer-code comparisons; they do not provide a physical flavor-symmetry mechanism or a Standard Model quantum recovery map.
