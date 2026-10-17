# Exact comparison: full-evaluation punctured Reed–Muller codes (2026)

**Status:** source-checked family-specific exclusion, plus a separate characterization of the unpunctured charge lattice. Not a general proof that the coupled quantum codes are new. This follows \`PUBLISHED_CODE_COMPARISON.md\`.

## Primary source and what it actually assumes

Tanay Saha and Shiroman Prakash, [*Sublogarithmic distillation in all prime dimensions using punctured Reed–Muller codes*](https://doi.org/10.1103/tzkz-9p2f), *Physical Review A* **113**, 042405 (2026); [arXiv:2510.10852](https://arxiv.org/pdf/2510.10852).

- **Section 2, Definition 1 (PDF p. 5):** the *unpunctured* "classical triorthogonal space" is both **cubic-orthogonal** and **pairwise self-orthogonal** over Fp.
- **Section 2, Eq. (1) (PDF p. 5):** puncturing \`t\` independent coordinates of a classical \`N\`-position matrix gives a CSS code of length \`N-t\` encoding \`t\` logical qudits, **provided the retained generator matrix has full rank**, as explicitly assumed there.
- **Section 3.1 (PDF p. 6):** a full generalized Reed–Muller evaluation matrix on Fp^m has **N = p^m** columns.
- **Section 3.2 (PDF p. 7):** the paper's Reed–Muller codes result from puncturing \`t\` evaluation points. The combinatorial puncture locations may vary, but the starting number of columns stays \`p^m\`.

Do not confuse the first paper's **unpunctured** evaluation code with a full Reed–Muller code of length 15 or with an arbitrary shortened/affine-subvariety construction.

## Proposition A — whole-family parameter exclusion, under the published rank assumption

Every code arising from the **full-evaluation, independent-\`t\`-puncture construction as specified in the paper** obeys

\`\`\`
n + k = p^m
\`\`\`

for an integer \`m ≥ 1\`, because \`n = p^m - t\` and \`k = t\`.

| Matter code | n+k | Required p^m | Result |
|---|---:|---|---|
| \`[[30,2,3]]_3\` | 32 | a power of 3 | impossible |
| \`[[45,3,3]]_3\` | 48 | a power of 3 | impossible |
| \`[[30,2,3]]_5\` | 32 | a power of 5 | impossible |
| \`[[45,3,3]]_5\` | 48 | a power of 5 | impossible |

This exclusion is **independent of Reed–Muller degree, puncture coordinates, or distance**, so it covers *all puncture locations and all admissible orders within the paper's full-rank framework*, including its computationally searched puncture schemes.

**Limits:** It does **not** exclude arbitrary generalized Reed–Muller *subcodes*, full-rank *failures* after puncturing, extra shortenings, non-full affine-variety evaluation sets, or a subsequently modified code with different \`n+k\`. Codes obtained by concatenation or by other constructions are not covered. In particular, a rank-deficient puncturing might encode fewer logical qudits than the number removed, invalidating the equation \`k=t\`; it is outside this proposition.

The executable check \`python3 reed_muller_family_exclusion.py\` verifies the entire table and positive controls: \`[[26,1,...]]_3\`, \`[[25,2,...]]_3\`, \`[[24,1,...]]_5\`, \`[[23,2,...]]_5\` all pass the necessary \`n+k=p^m\` test (which does not by itself assert they attain any particular distance).

## Proposition B — cubic anomaly cancellation is not the same as the paper's *two-condition* definition

The Standard Model representation here is built from the **even-parity sign set**

\`\`\`
S = {q in {−1,+1}^5 : product(q_i)=+1}, |S|=16.
\`\`\`

Its integer moments satisfy, for coordinate indices i,j,k,

\`\`\`
sum_{q in S} q_i = 0
sum_{q in S} q_i q_j = 16 δ_ij
sum_{q in S} q_i q_j q_k = 0 .
\`\`\`

Every multi-family **Cartan-linear** row evaluates a linear function of \`q\` on each family. Consequently the complete unpunctured charge matrix has vanishing *integer cubic* products of arbitrary rows, including repeated indices, and thus vanishing cubic anomalies.

However, in both of the B3/B5 family-coupled lattices the first differential row evaluates as \`(q_1, −q_1)\` on the two generations (or \`(q_1,−q_1,0)\` for three generations), with all 16 spinor weights in each block. Its quadratic self-product is

\`\`\`
sum_j (G_{row,j})^2 = 16+16 = 32.
\`\`\`

**For every odd prime p, 32 ≠ 0 (mod p)**. Thus the supplied unpunctured multi-family charge span does **not** satisfy the *quadratic* self-orthogonality demanded by Saha–Prakash Definition 1, even though its cubic anomaly sums vanish exactly.

This has a precise interpretation:

- The present matrices are **cubic-triorthogonal** charge lattices with a different quadratic inner-product structure.
- The Saha–Prakash unpunctured Reed–Muller spaces obey **both** the cubic and quadratic conditions.
- This difference is a property of these particular unpunctured matrix representations. It is **not an invariant of the resulting punctured quantum code under arbitrary local Clifford/unitary transformations**, and must not be cited as an independent no-equivalence theorem.

### Exact evaluation-code viewpoint

The 16 SO(10) sign columns can also be described as evaluations of five linear coordinate functions on **16 selected points of Fp^5**, not as evaluations on *all* p^m affine points. Equivalently, four independent binary signs select the fifth by parity, \`q_5 = q_1 q_2 q_3 q_4\`. In character-language, each coordinate row is a character of the 4-dimensional Boolean group over Fp (p odd).

The family-coupled code uses *linear transforms of this restricted evaluation design in each block*, and punctures the singlet q=(1,1,1,1,1) from each block. This is a useful classification route: **codes on character/evaluation subsets**, not necessarily full generalized Reed–Muller codes over Fp^m. It does not imply the constructed code is equivalent to a previously published affine-variety code.

## Reproduction

\`\`\`bash
python -m pip install numpy
python reed_muller_family_exclusion.py
\`\`\`

The script checks all source-scoped parameter obstructions, 16-point parity moments, full integer cubic constraints, and the quadratic witness in p=3/5, 2/3-family examples. The GitHub Actions smoke workflow includes it.

## Best remaining comparison

The results now exclude (i) the 2025 Prakash–Saha \`[[9m-k,k,2]]_3\` construction, and (ii) the **2026 Saha–Prakash full-evaluation, full-rank punctured RM family**.

For a true novelty claim, we still need an **explicit published code of matching \`(p,n,k)\`** or an exhaustive published classification for the relevant odd-prime lengths. The next useful target is an *affine-variety evaluation code* or an indecomposable odd-prime triorthogonal construction, accompanied by its actual check matrices and local-Clifford/monomial equivalence criteria.
