# E8 and the generation codes: the chirality obstruction in code language

**Status, 18 October 2026:** short exploratory comparison, not peer reviewed. Verified by `python3 e8_comparison.py` (needs numpy; runs in about a second; ends with `ALL E8 COMPARISON CHECKS PASSED`).

**Question.** Does Garrett Lisi's E8 unification proposal ([*An Exceptionally Simple Theory of Everything*](https://arxiv.org/abs/0711.0770), 2007) bear on the generation codes in this repository?

**Short answer.** It doesn't help, but it mirrors our results closely. E8 meets the same two walls these codes meet, **chirality** and **family labels**, and each wall has an exact code-theoretic counterpart.

## Shared starting point

One generation here is the 16 SO(10) spinor weights: the sign vectors in {±1}⁵ with an even number of minus signs.

Under the standard branching E8 ⊃ SO(16) ⊃ SO(10) × SO(6), we have 248 = 120 + 128, and the 128 half-integer roots (scaled by 2) are the even-parity sign vectors in {±1}⁸. Splitting the 8 coordinates as 5 + 3 gives

```
128 = (16, 4) + (16̄, 4̄),
```

that is, four copies of one generation, labelled by SU(4) ≅ SO(6) weights, plus four mirror copies. The script checks that the second block is exactly the set of negatives of the first.

## The external result

Distler and Garibaldi, [*There is no "Theory of Everything" inside E8*](https://arxiv.org/abs/0905.2658) (Commun. Math. Phys., [doi:10.1007/s00220-010-1006-y](https://doi.org/10.1007/s00220-010-1006-y)), prove two things. Under their stated assumptions (Lorentz group SL(2,ℂ), compact gauge group commuting with it, restricted spins), embeddings of gravity plus the Standard Model into any real form of E8 give only self-conjugate, non-chiral fermion representations. And three generations cannot fit as distinct weight vectors: they would need at least 180 dimensions, but only 112 or 128 are available.

## Wall 1: chirality ⟺ anomaly cancellation ⟺ distance

| Charge lattice (8 Cartan rows of E8) | Triorthogonal? | Code after puncturing the four ν^c columns |
|---|---|---|
| Chiral block (16, 4) alone | **No.** The cross cubic sum over the three SU(4) Cartan rows is 64, which fails over ℤ and mod 3, 5, 7. | Not a valid triorthogonal code |
| Full non-chiral 128 | Yes, over ℤ | **[[124, 4, 1]]_p** for p = 3, 5, 7 |

The repository's construction works because a single *chiral* generation cancels its cubic anomalies non-trivially, and that cancellation is exactly triorthogonality.

The chiral E8 block carries an uncancelled SU(4)³-type anomaly. Restoring cancellation the only way E8 allows, by adding the mirror block, makes cancellation trivial: every column has its negative present. A punctured ν^c column then reappears as its mirror among the physical columns. That mirror is a weight-one logical operator, so the code has **distance 1**.

The Distler–Garibaldi non-chirality obstruction thus has a direct counterpart here: there is no anomaly-free *and* protecting lattice inside this E8 decomposition.

## Wall 2: additive family labels give only the family-badge code

Ignore the anomaly and use E8's extra SU(4) Cartan directions as family labels on the chiral block. Every species in family g then receives the same label w_g, the column for species q being (q, w_g). This is exactly the "uniform family badge" construction (attempt B in the README). Twin pairs (q, w_g) − (q, w_h) coincide with differences of punctured columns, giving **[[60, 4, 2]]_p** at p = 3, 5, 7.

The distance-3 "smart badge" constructions need charges of the form f_g · (b·q), which *multiply* a family label by an SO(10) weight. The Cartan subalgebra of any Lie group containing SO(10) × G_family as a direct product acts *additively*, as (x, y)·(q, w) = x·q + y·w, so it can never produce such products.

Family-dependent SO(10)-Cartan charges in general need something like one SO(10) per family (rank 15 for three families). That exceeds E8's rank of 8.

## Scope

- **The decomposition used.** The checks use the standard textbook branching E8 ⊃ SO(10) × SO(6). Lisi's actual particle assignments differ and also place gravity inside E8. Distler–Garibaldi covers his framework under their assumptions; the computations here cover this decomposition. Other E8-inspired constructions outside those assumptions are not excluded.
- **Integer roots not used.** The 112 integer roots (inside SO(16)'s adjoint) contain no copy of the 16 and play no role here.
- **What the comparison shows.** It concerns charge lattices and their codes. It makes no claim about the physics of E8 models beyond what the cited theorem states.

## Takeaway

The family-coupled codes in this repository and the E8 programme fail for the same structural reason. A family structure that is **both chiral and sharply labelled** is not available, whether you look for it in one large symmetry group (E8: mirrors and additive labels) or in the low-energy Standard Model (flavour mixing; see `OPEN_STEPS_NOTE.md` and `EXACT_SYMMETRY_AND_CLASSIFICATION.md`).
