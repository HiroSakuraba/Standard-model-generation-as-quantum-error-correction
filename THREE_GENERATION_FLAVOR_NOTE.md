# Three generations, mixed logical errors, and the obstruction from flavor mixing

**Research extension, October 17, 2026.** Authorial status: independent exploratory calculation, **not peer reviewed**. Continues the calculations in `RESEARCH_NOTE.md`. Underlying paper: Shiroman Prakash and K. L. Ramya, *One generation of the Standard Model as a quantum error-correcting code*, [arXiv:2610.21879v1](https://arxiv.org/abs/2610.21879), posted October 15, 2026. The authors derive a formal `[[15,1,3]]_p` code from an SO(10) generation with a conjugate right-handed neutrino. We extend that *mathematical construction*; we do not claim ordinary matter actually corrects quantum errors. These are formal qudits indexed by fermion species, not direct physical subsystems.

## Results: a compact statement

Two three-generation triorthogonal **CSS** codes, with 45 physical qudits and three logical qudits, are constructed from 48 SO(10) spinor weights by deleting one conjugate right-handed neutrino from each generation:

- `[[45,3,3]]_3`: 11 independent charge rows, 8 independent X stabilizers, 34 Z stabilizers, X distance **7** and Z distance **3**.
- `[[45,3,3]]_5`: 9 independent charge rows, 6 independent X stabilizers, 36 Z stabilizers, X distance **9** and Z distance **3**.
- The distance-3 construction extends algebraically to every **odd prime**: use the first B matrix for `p=3` and the second for `p>=5`. X distances are computed exactly for p=3 and 5, not for every p.
- For either explicitly checked prime there are exactly **75** three-coordinate logical-Z supports: **25 per generation**, none crossing generations.
- For four-coordinate supports split across two generations as 2+2, there are **88 per generation pair at p=3** and **396 per pair at p=5**. Since no weight-three logical support fits inside such a 2+2 support, these are genuine *minimum-support-size-four* cross-generation logical errors. The corresponding 2+2 supports are absent for the block-diagonal direct-sum construction.
- The logical transversal cubic phase is a **tensor product** of one-qudit phase gates, not an entangling logical gate.

These observations **do not** show a better quantum error correction code than known constructions; fewer X stabilizers mean more Z stabilizers, not less total stabilizer rank, and no circuit-level overhead, threshold, or decoder study has been performed.

## Algebraic construction

Let the 16 weights in each SO(10) spinor be

```
Q = {q in {+1,-1}^5 : the number of negative entries is even}
q0 = (1,1,1,1,1)  # conjugate right-handed neutrino
```

Give the three generations the **integer** family vectors

```
f_e   = (+1,  0)
f_mu  = (-1,  0)
f_tau = ( 0, +1)
```

They are affinely independent modulo **every odd prime**, because the determinant of their two differences is ±2.

Choose one of the differential charge matrices

```
B3 (for p=3) = [1 0 0 0 0]      (3 x 5)
               [0 1 0 0 1]
               [0 0 1 1 0]

B5 (for p>=5)= [1 0 0 0 0]      (2 x 5)
               [0 1 1 1 1]
```

For a generation label `g` with family vector `f_g`, assign its SO(10) weight `q` the column

```
h(q,g) = ( q , f_g[0] * Bq , f_g[1] * Bq )^T
```

The complete `48`-column integer charge matrix has `5+2m` rows. Deleting the three `q0` columns gives 45 physical coordinates, while the independent deleted columns give three logical degrees of freedom.

### Cubic anomaly cancellation over the integers

Every row combination restricts to a linear functional `a_g · q` on each generation's 16 SO(10) weights. The underlying 16-weight SO(10) set obeys

`sum_(q in Q) (a_g · q)^3 = 0`

**over integers for arbitrary a_g**. The full 48-column lattice therefore has identically vanishing mixed cubic anomalies for any B and any family labels. The script tests all distinct triples of generator rows, not merely diagonal anomalies, and reports PASS for both B choices.

The integer lift also gives exactly separable transversal cubic phases after puncturing:

`U_logical = U^(-1) on each of the 3 logical qudits`,

with the appropriate ninth-root-of-unity qutrit gate. Nothing here implies a transversal *entangling* logical gate.

### Why the Z distance remains three

For an odd prime, no SO(10) non-singlet weight `q != q0` is a scalar multiple of `q0`, and no two distinct non-singlet weights span `q0` (as used in the original paper's Theorem 1 proof). So a one- or two-column logical error, if any, would have to be a matched pair of identical non-singlet weights in two distinct generations, with coefficients summing to zero in the common five charge coordinates.

Its remaining differential component is

`(f_g - f_h) tensor Bq`.

It lies in the span of the removed singlet columns, with zero common component, exactly when `Bq` is proportional to `Bq0` modulo p (with nonzero proportionality). Our B matrices exclude this for **every** non-singlet q:

- For `B3`, the first projected coordinate fixes the proportionality and the remaining two imply that all signs are equal; the all-negative weight has forbidden odd parity.
- For `B5`, if `k` is the number of negative entries among q2...q5, then proportionality requires `4 - 2k = 4q1 (mod p)`. For p>=5, `0<=k<=4` forces either all positive or all negative (the latter forbidden).

Thus `d_Z >= 3`; the familiar three-fermion relation wholly within a single generation gives `d_Z <= 3`. The accompanying code explicitly enumerates all 45 and all 990 one- and two-coordinate supports, then all 14,190 three-coordinate supports, over F3 and F5. Full CSS distances for p=3,5 follow from separately computed `d_X=7,9`. For arbitrary odd p the standard triorthogonal property `d_Z <= d_X` gives full distance three.

The result should be understood as an **algebraic family of codes**, not an assertion that the extra charge symmetries occur in nature.

## Weight-four errors: the coupling is real but appears later

A simple direct sum of three distance-three single-generation codes cannot have a Z logical error whose support contains only **two** coordinates in generation A and **two** in generation B: each block's individual Z distance is three. The reduced-charge three-generation codes instead admit such mixed four-coordinate logical errors:

| Exact finite-field enumeration | `p=3` | `p=5` |
|---|---:|---:|
| Three-coordinate logical supports, total | 75 | 75 |
| Three-coordinate supports spanning generations | 0 | 0 |
| Mixed 2+2 four-coordinate logical supports, **per family pair** | 88 | 396 |
| Mixed 2+2 four-coordinate logical supports, all 3 pairs | 264 | 1188 |
| All four-coordinate supports *containing some logical*, including 3+1 supersets | 3579 | 4488 |

**Definition caution:** The last row enumerates four-element *sets* that support at least one nontrivial logical vector, which may be supported on just three of those elements. The 2+2 counts are different: they cannot hide a three-element logical, because no weight-three logical spans generations. The enumerator counts supports, not distinct Pauli error coefficient patterns or actual probabilities.

The direct-sum comparison is a reason the coupled construction is structurally nontrivial at the level of physical coordinate supports. It is **not** a proof of novelty under arbitrary quantum-code equivalences.

## Exact flavor-mixing obstruction

Quark and lepton mass terms mix family labels. The observed Cabibbo–Kobayashi–Maskawa (**CKM**) quark matrix and Pontecorvo–Maki–Nakagawa–Sakata (**PMNS**) neutrino mixing matrix make flavor-sensitive diagonal charges physically problematic. The Particle Data Group's 2026 reviews document these effects: https://pdg.lbl.gov/2026/reviews/contents_sports.html, entries 12 and 14. The quark sector requires nontrivial flavor mixing; the lepton sector does as well if neutrino masses are incorporated.

For a two-flavor illustration, take a diagonal family charge `D = diag(d1,d2)` and a rotation

```
       [ cos(theta)   sin(theta) ]
V(theta) =
       [-sin(theta)   cos(theta) ]
```

In the rotated basis,

```
(V.T D V)_12 = (d1-d2)*cos(theta)*sin(theta).
```

If mixing is nonzero and nontrivial (`sin(theta)*cos(theta) != 0`), the family charge stays diagonal in both bases **only when `d1=d2`**. Hence any proposed *exact* conserved U(1) that distinguishes the two generations is incompatible with a generic flavor-mixing mass interaction, unless extra charged fields/spurions change the assumed symmetry transformation.

For the B matrices above, every SO(10) weight has at least one nonzero projected B coordinate, because B's first row is q1, which is always ±1. Thus two generations have distinct charge vectors for that field species under at least one extra charge row. The extra family-sensitive symmetries **cannot** all remain exact in the presence of nontrivial mixing connecting those generations. Further, arbitrary family-dependent SO(10) Cartan rows are not automatically generators of a gauge group compatible with the observed single SU(3) x SU(2) x U(1): many do not commute with the nonabelian gauge action on a family multiplet. Vanishing torus cubic sums is necessary but not sufficient to embed these hypothetical charges into a consistent interacting gauge theory.

An exact rational arithmetic toy test included in the script takes `cos=4/5`, `sin=3/5`, `D=diag(1,-1)` and produces `|(V.T D V)_12|=24/25`. This is **not** a realistic CKM/PMNS angle and **not** a quantum-code failure probability. It is simply an exact check of the matrix identity.

With only truly *generation-blind* charge rows, the three removed neutrino columns coincide, so their rank falls from three to one. For this charge-lattice construction, the three separate logical symbols collapse to **one** logical degree of freedom. This is an algebraic obstruction to treating flavor-sensitive rows as physically exact conserved symmetries, not a no-go theorem for alternative quantum codes or symmetry-protected constructions.

**Important distinction:** The original SO(10) charge lattice is a representation-theoretic structure. Some of its generators are broken in the low-energy theory. Calling any of its generator rows a *physical stabilizer measurement* would require an explicit microscopic Hilbert space and dynamics, neither provided by the original paper nor by this extension.

## Next scientific tests

1. **Work in a concrete flavor Lagrangian.** Write actual Yukawa matrices and Higgs/spurion transformation rules, then identify the exact residual symmetries and approximate ones. Use charged-lepton and neutrino basis conventions consistently. A Yukawa-covariant *formal* code may exist, but merely transforming species labels does not supply a physical QEC channel or recoverability theorem.
2. **Classify correlated logical channels**, not just their supports. Compute exact Pauli-coefficient enumerators and distillation acceptance/output polynomials, including correlated two-generation four-coordinate events.
3. **Look for a physical encoding map** or an operational reinterpretation in quantum simulation, where charge and gauge constraints can actually be measured. Without this, protection remains combinatorial, not dynamical.
4. **Investigate higher-distance extensions**, avoiding naïve identical family duplication. Derive an integer anomaly-free chiral charge lattice with higher distance and study whether nonseparable logical cubic phases are possible using a different puncturing structure; the present construction is separable by definition.

## Reproduce

Run from this directory:

```
python3 matter_qec_study.py
python3 three_generation_flavor.py
```

Both programs use only the Python standard library. `three_generation_results.txt` records the exact run. This extension tests numerical field dimensions `p=3,5` and proves the odd-prime Z-distance analytically within the stated construction; it does not claim exhaustive code-classification novelty or an experimentally realized error-protection mechanism.
