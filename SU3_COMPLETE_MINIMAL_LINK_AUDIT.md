# Complete SU(3) two-plaquette minimal-link Gauss sector

**18 October 2026 — independently executable exact finite-lattice certificate.** The companion program is [`su3_complete_minimal_links.py`](su3_complete_minimal_links.py), which imports the preceding [`su3_two_plaquette_spin_network.py`](su3_two_plaquette_spin_network.py). Python 3.9+; SymPy and NumPy. This is **not** a physical glueball spectrum, Yang–Mills mass-gap proof, or demonstrated quantum error-correcting device.

## Research decision: stop growing unrelated character-trial bases

The preceding nested **7/11/17/19-state Wilson trace** calculation retained some irreducible representations `6`, `bar6` and `8` on individual links. It did *not* systematically retain every gauge-invariant state with a prescribed link cutoff. Thus its convergence between 17 and 19 states cannot be interpreted as convergence of the **minimal link truncation** used in Chen, Müller and Yao (2026). We instead enumerate the complete gauge-singlet Hilbert space of **exactly the same seven-link, six-vertex open rectangle as in the earlier seven-state test**, imposing **per-link Peter–Weyl irreducible representations `1 ⊕ 3 ⊕ bar3`**. No dynamical matter is included.

The graph is two adjacent square faces, vertices `0--1--2` on the top and `3--4--5` below, with links (in code index order):

```
(0,1), (1,2), (3,4), (4,5), (0,3), (1,4), (2,5)
```

All six boundary/corner vertices have degree two except shared-face vertices 1 and 4, which have degree three. Endpoints of an oriented `3` link carry `3` at the source and `bar3` at the target; an oriented `bar3` link carries the opposite representations.

## Exact completeness certificate: nine physical states

Enumerate all `3^7 = 2187` possible assignments of the link representation labels `-1,0,+1` (`bar3,1,3`). At a **degree-two** vertex, a singlet exists only for no flux or two conjugate representations, and the intertwiner multiplicity is one. At a **degree-three** vertex, a singlet exists only with no flux; a conjugate pair and one trivial representation; or three fundamentals (`epsilon_{ijk}`), or three antifundamentals (`epsilon^{ijk}`). In each allowed case the singlet intertwiner space has dimension **one**; no higher-multiplicity vertex has been omitted on this particular graph.

The executable directly performs this representation-theory test for **all 2,187** assignments and obtains **nine** admissible gauge-singlet networks:

| Nontrivial occupied links | Number of states | Interpretation |
|---:|---:|---|
| 0 | 1 | Vacuum |
| 4 | 4 | Two orientations each on left/right square |
| 6 | 2 | Two orientations on the outer perimeter |
| 7 | 2 | Simultaneous-flux `epsilon` junction networks at both degree-three vertices |
| **Total** | **9** | Complete *within the chosen link-irrep cutoff* |

The last two configurations, in the fixed seven-edge orientation, are `(+,-,-,+,-,-,-)` and its color conjugate `(-,+,+,-,+,+,+)`. They are distinct sectors of the representation-label Hilbert space, not higher-`6` or `8` link excitations. Their three-strand `epsilon` junctions are *pure-gauge invariant tensors*; **they do not represent dynamical matter baryons**.

### Independent microscopic Gauss checks

The script constructs each new junction vector by contracting the explicit SU(3) invariant `delta` at every degree-two vertex and `epsilon` at the two degree-three vertices. Each junction wavefunction has **2,916 nonzero exact integer color-index components**, with squared unnormalized norm 2,916. It acts with four complex root generators and two Cartan generators separately at **each of six vertices**, using fundamental `T` and antifundamental `-T^T` link-end actions. Both new networks have **72** zero-residual generator–state checks. The existing seven networks also pass their **252** independent microscopic checks, for **324 exact checks** total. The root and Cartan generators span complexified `su(3)`, so all eight real generators vanish as well. This is a *direct*, not assumed, local Gauss certificate.

## A complete and correctly normalized nine-state basis

After a maximal-tree gauge choice, the two independent fundamental face holonomies are `U,V ∈ SU(3)`. Write `u=tr U` and `v=tr V`. The seven original normalized basis functions are

```
1, tr U, tr U†, tr V, tr V†, tr(UV), tr(UV)†.
```

The additional **unnormalized** network representatives are

\[
B_+=\operatorname{tr}U\operatorname{tr}V^\dagger-
\operatorname{tr}(UV^\dagger),\quad
B_-=\operatorname{tr}U^\dagger\operatorname{tr}V-
\operatorname{tr}(U^\dagger V).
\]

These are the antisymmetric `bar3` and `3` shared-link fusion channels of the simultaneously excited plaquettes. Their **squared norms are exactly `4/3`**, and they are orthogonal to each other and all seven original states. Thus the raw Gram matrix is `diag(1,1,1,1,1,1,1,4/3,4/3)`; normalized junction states are `(sqrt3/2) B_±`. The code independently obtains this Gram matrix through exact Haar integration over SU(3)×SU(3), rather than falsely assuming that every polynomial is normalized.

Every new state has seven occupied fundamental/antifundamental links, so its electric energy in the declared units is **`7*(4/3)=28/3`**. The complete electric diagonal is

\[
D=\operatorname{diag}(0,16/3,16/3,16/3,16/3,8,8,28/3,28/3).
\]

The code also constructs each of the seven individual **link** Casimirs. Unlike the 19-state character basis, this minimal-link model has no `6` or `8` link Casimir eigenvalues.

## Exact fundamental Wilson matrices — why seven states are incomplete

The fundamental face multiplication operators are `W_L=tr U` and `W_R=tr V`, compressed into the complete nine-dimensional cutoff. Every one of the 81 entries of each **covariant** matrix `S_ij=<f_i|W|f_j>` is evaluated by **exact SU(3) Haar moments**, including the `epsilon epsilon/6` determinant moment. For the unnormalized basis the actual coefficient operator is `G^{-1}S`, and the Hilbert-space adjoint is `G^{-1}S†`, **not** ordinary transpose of `G^{-1}S`.

Their original seven-by-seven principal submatrices reproduce **every** exact matrix element of `su3_two_plaquette_spin_network.py`. The newly exposed covariant matrix elements are

\[
\langle B_+|W_L|R_-\rangle=\frac23,\quad
\langle B_-|W_L|O_+\rangle=\frac23,\\
\langle B_+|W_R|O_+\rangle=\frac23,\quad
\langle B_-|W_R|L_-\rangle=\frac23.
\]

With the **normalized** junction wavefunctions each transition coefficient becomes `1/sqrt(3)`. The full `Q W P` block from the old seven states `P` to the omitted two-state complement `Q` has **operator norm exactly `1/sqrt(3)`** for either elementary fundamental Wilson insertion. That is direct, nonzero truncation leakage: seven states are *not* closed under the very same local magnetic interaction. (The full unprojected multiplication operators commute; any failure of commutation of hard-projected Wilson matrices is a cutoff artifact.)

## Error-correction verdict and complete-algebra structure

The code `|0_L>=tr(U)`, `|1_L>=tr(U†)` still corrects arbitrary errors in the **complete algebra of the seven link-electric Casimirs**, because those operators act with identical eigenvalues on the two states. However the left fundamental **Hermitian** Wilson magnetic error `B_L=W_L+W_L†` still compresses onto the code as exact logical `X`: 

\[
P B_L P=\begin{pmatrix}0&1\\1&0\end{pmatrix}.
\]

Thus the physically permitted magnetic operator defeats the proposed *electric-only* code without ever violating Gauss's law. The code distinguishes correcting the **complete generated error algebra** from correcting a small, physically bounded *linear error set*.

For exact rational nine-by-nine error generators:

| Error algebra generated by | Exact commutant dimension | Complete-algebra interpretation |
|---|---:|---|
| `D`, real Hermitian left/right magnetic plaquettes | **2** | Two inequivalent charge-conjugation sectors of dimensions **5 and 4**, hence algebra `M5⊕M4`; **no multiplicity qubit** |
| Above plus **all seven** link electric Casimirs | **2** | Same noiseless-subsystem conclusion |
| `D`, both complex fundamental Wilson multipliers and their adjoints | **1** | Full `M9`; **no noiseless subsystem** |

Charge conjugation is an explicit central involution for the real/Hermitian set, with ranks **5 and 4**, checked exactly. These claims apply only to **arbitrary errors from the full specified generated algebra** on the nine-state host. They are not a theorem excluding all finite-weight active quantum codes, noiseless algebras on larger lattices, or all possible physical SU(3) encodings.

## Finite-volume spectral sensitivity across magnetic couplings

Illustrative Hamiltonian `H(g)=D−g(B_L+B_R)/2`. Double-precision diagonalization uses the normalized Gram metric. This is the **same coefficient convention** across the seven- and nine-state calculations but **not a calibrated identification** with the Hamiltonian parameters used by Chen–Müller–Yao.

| Dimensionless magnetic coefficient `g` | 7-state ground | 9-state ground | 7-state gap | 9-state gap |
|---:|---:|---:|---:|---:|
| 0.25 | -0.0119761 | -0.0119761 | 5.2203068 | 5.2203094 |
| 0.5 | -0.0487297 | -0.0487734 | 5.1320631 | 5.1321067 |
| 1 | -0.1989836 | -0.1997103 | 5.0323169 | 5.0330436 |
| 2 | -0.7853342 | -0.7957794 | 5.1186676 | 5.1291127 |
| 4 | -2.6924164 | -2.7846110 | 6.0257498 | 6.1179443 |

The nine-state ground energy is no higher than the seven-state one, as expected from the variational principle. The **first excitation gap has no monotone variational guarantee**. These small-lattice energy gaps have arbitrary units, no continuum/thermodynamic scaling, no physical `J^PC` glueball classification, and **are not glueball masses**.

## Relation to primary literature; remaining engineering work

Chen, Müller and Yao, *Minimally truncated SU(3) lattice gauge theory and string tension*, *Physical Review D* **114**, 074502, published **5 October 2026**, https://doi.org/10.1103/1qt7-3ywl, use link irreps `1,3,bar3` and construct a minimally truncated square-plaquette-chain Hamiltonian. Our **complete link representation inventory for two adjacent squares** now uses the same declared local irrep *set*, but **we have not established numerical equivalence to their published plaquette-transition amplitudes, boundary conventions or coupling normalization**. Do **not** advertise this as a replication of their full matrix or spectrum. Their construction is the correct next source for an independent coefficient-by-coefficient comparison and extension to more plaquettes.

Research priorities after this milestone:

1. Reproduce their precise two-plaquette Hamiltonian using its original boundary, orientation, representation, and magnetic-operator normalization conventions. Compare every nonzero transition against our nine-state exact Haar matrix, and reconcile any discrepancy **before** doing larger volumes.
2. Quantify projected Wilson truncation effects relative to a larger **systematic per-link** irrep cutoff such as `1,3,bar3,6,bar6,8`, including all possible intertwiner multiplicities. Do not substitute a selected 19-polynomial variational basis for this complete cutoff.
3. Build a three-plaquette and larger gauge-invariant Hilbert space using transfer matrices/dynamic programming over link irreps, then compute finite-size scaling and gauge-invariant connected-loop correlation functions.
4. Only with significantly larger geometry and spin/parity/charge-conjugation sectors attempt glueball interpolators, or physically calibrated noise models and recovery operations.

Other primary source: Kogut and Susskind, *Hamiltonian formulation of Wilson's lattice gauge theories*, *Physical Review D* **11**, 395 (1975), https://doi.org/10.1103/PhysRevD.11.395. Mathematical correction standard: Knill and Laflamme, *Theory of quantum error-correcting codes*, *Physical Review A* **55**, 900 (1997), https://doi.org/10.1103/PhysRevA.55.900.

## Reproduce

```sh
python -m pip install sympy numpy
OPENBLAS_NUM_THREADS=1 python su3_two_plaquette_spin_network.py
OPENBLAS_NUM_THREADS=1 python su3_complete_minimal_links.py
python -m py_compile su3_complete_minimal_links.py
```

Exact assertions: exhaustive representation sectors and intertwiner count, direct 324 SU(3) Gauss generator residuals on the explicit sparse link-color states, exact Haar Gram, seven-link electric Casimirs, fundamental Wilson matrices, normalized omitted-state coupling, charge-conjugation ranks, commutant dimensions and the local magnetic logical-X obstruction. Only eigenvalues/gaps are double-precision numerical and not interval-certified. The companion `.txt` and `.json` files are reproduced by the script's output.