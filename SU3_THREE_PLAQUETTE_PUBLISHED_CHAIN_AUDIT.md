# Three SU(3) square plaquettes: published Hamiltonian and exact error algebra

**18 October 2026 — reproducible finite-lattice research calculation.** Companion programs: [`su3_three_plaquette_chain.py`](su3_three_plaquette_chain.py) and [`su3_three_plaquette_error_algebra.py`](su3_three_plaquette_error_algebra.py). Exact tests use SymPy 1.14 and a finite-field modular rank witness; numerical spectra use NumPy double precision. These are **not** measurements of glueball masses, a continuum Yang–Mills mass-gap proof, a hardware quantum memory, or a new Standard Model code.

## 1. Primary source and continuity with the two-plaquette benchmark

Vincent Chen, Berndt Müller and Xiaojun Yao, *Minimally truncated SU(3) lattice gauge theory and string tension*, **Physical Review D 114**, 074502 (published 5 October 2026), DOI [10.1103/1qt7-3ywl](https://doi.org/10.1103/1qt7-3ywl), is the source for the minimal `1⊕3⊕bar3` link cutoff, the open square-plaquette chain, the electric/magnetic normalization, and the recoupling matrices. See Eqs. (3), (5), (8), (10), (11), **Table I** (published PDF p. 4) and **Table VI** (published PDF p. 14). The comparable table in the January arXiv v1 preprint was numbered Table IV; the source references in this report follow the published version.

PR #27 established the full **nine-dimensional** Gauss-singlet space of *two* plaquettes and PR #28, now merged into PR #27's branch, checked its full Hamiltonian against the source. PR #27 remains an open parent pull request; this extension is stacked on that branch. Our script explicitly imports its independent two-face benchmark and compares all directed Wilson/electric entries before constructing the new model.

## 2. The complete 27-state model at the fixed minimal-link cutoff

The open chain has **three square plaquettes**, four upper and four lower vertices, **ten links** (three top, three bottom, four vertical), with top/bottom exterior flux and left/right boundary flux fixed to zero. The face-qutrit labels `q=(q0,q1,q2)` range independently over `Z3^3`, giving **3³ = 27 gauge-invariant states** under the published minimal-link, flux-free boundary prescription.

The oriented link representations are fixed by the face differences:

\[
\ell(q)=(q_0,q_1,q_2,-q_0,-q_1,-q_2, -q_0,q_0-q_1,q_1-q_2,q_2)\pmod3.
\]

`0` means trivial representation; `1` and `2` mean fundamental `3` and antifundamental `bar3`. The link labels are not independent once the Gauss constraints and zero exterior flux are fixed: every local trivalent flux balance follows algebraically from the difference rule. **Unlike the previous two-face microscopic spin-network verifier, this three-face program does not explicitly enumerate all eight SU(3) Lie-generator actions on every link-color configuration.** It verifies the complete published face-to-link prescription, its one-to-one mapping and the already source-benchmarked two-face restrictions.

The gauge-invariant electric Hamiltonian is the sum of quadratic Casimirs of these ten links,

\[
D(q)=\frac43\left[2\sum_{i=0}^2[q_i\ne0]+[q_0\ne0]+[q_2\ne0]+[q_0\ne q_1]+[q_1\ne q_2]\right].
\]

Here `[condition]` is one if true and zero otherwise. It follows from the stated link cutoff and reproduces the parent two-plaquette Casimir exactly on `q2=0`.

## 3. Exact source-to-implementation magnetic operator verification

For each dynamical plaquette `i`, its immediate left/right neighbor labels are `q_(i−1)` and `q_(i+1)`; at the two ends the missing neighbor is fixed to `0`. The **source's complete published Table VI** provides a 3×3 Hermitian magnetic matrix `M_(left,right)` for **every one of the nine pairs** of neighbor labels.

The source verifier independently transcribes both published **Table I**, which contains nine directed trivalent-vertex recoupling matrices, and the **Table VI** complete-plaquette blocks. Using the four vertices around a square face, it forms the elementwise product

\[
M_{(l,r)}=V_{(0,r)}\odot V_{(r,0)}\odot V_{(0,l)}\odot V_{(l,0)}.
\]

**All nine 3×3 matrices (81 entries) agree exactly** with the separately transcribed Table VI. The matrices exercise the nontrivial interior blocks including `M_(1,1)[0,1]=1/9`, `M_(1,2)[0,1]=1/(3sqrt(3))`, and `M_(2,2)[0,2]=1/9`. These were **not used** by the original boundary-zero two-face model, which needs only `(0,n)`.

Let `W_i` shift the face label `q_i → q_i+1 mod 3` by the exact matrix element in the appropriate `M_(left,right)`. The real magnetic term is `B_i=W_i+W_i†`. The complete source-normalized 27×27 Hamiltonian, including the constant, is

\[
\boxed{H(g^2)=g^2D+\frac{18I_{27}-B_0-B_1-B_2}{g^2}.}
\]

**Independent regression:** Freeze `q2=0` and project the three-face operators to the remaining nine states. Their electric matrix and both original directed Wilson matrices agree *entry by entry* with the independently benchmarked two-face source Hamiltonian; the third-face Wilson compresses to zero. Every one of the nine interior neighbor blocks also appears as a literal 3×3 submatrix of `B_1`.

The magnetic operators in the published minimal link truncation are **not unitaries** in general. They are Hermitian or adjoint pairs as explicitly constructed; do not call a finite-dimensional projected `W_i` a unitary fundamental Wilson holonomy.

## 4. Exact algebra result, no floating-point rank assumption

Let `C` globally conjugate all three `Z3` face variables, `q_i→−q_i`; it is an exact involution commuting with `D,B_0,B_1,B_2`. Its eigenvalue sectors have dimensions **14 (+)** and **13 (−)**. Left/right reflection `R` also commutes with the full Hamiltonian **with equal plaquette couplings**. The finer `C,R` symmetry-sector dimensions are **10, 4, 8, 5** for `(+, +), (+, −), (−, +), (−, −)` respectively. Reflection exchanges individual `B_0` and `B_2` and is **not** a commutant member of the complete individually generated local-error algebra.

For an exact rank certificate, every commutator equation `[X,M]=0` is a homogeneous linear system with coefficients in `Q(sqrt(3))`, with denominators invertible modulo 97. Since `10²≡3 mod97`, the specialization `sqrt(3)→10` is a ring homomorphism to `F_97`. A nonzero minor modulo 97 is a nonzero minor in characteristic zero, so **modular rank lower-bounds the characteristic-zero rank**, giving an upper bound on the exact commutant dimension. The verifier performs sparse Gaussian elimination, with no numerical tolerances, obtaining:

| Specified complete error *-algebra | Verified modular commutator rank | Exact commutant dimension | Algebra structure |
|---|---:|---:|---|
| Generated by `D,B_0,B_1,B_2` | **727** of 729 | **2** (explicit independent `I,C` attain bound) | `M_14(C) ⊕ M_13(C)` (dimension **365**) |
| Generated by `D,W_i,W_i†` for all 3 faces | **728** of 729 | **1** (explicit identity attains bound) | `M_27(C)` (dimension **729**) |

Because the algebras are finite-dimensional *-algebras, exact commutant dimension and the explicit charge-conjugation sectors give the stated full block decomposition by the double-commutant theorem. The **complete specified algebra** has no noiseless multiplicity subsystem of logical dimension 2 or more. This is a **strong error-model assumption** allowing arbitrary products: it does **not** rule out active correction for a restricted single-fault error span, other encodings on bigger lattices, or a quantum error-correcting code in actual SU(3) quantum chromodynamics.

### A concrete code-specific counterexample

Consider `|0_L> = |0,1,0>` and `|1_L> = |0,2,0>`. Every individual electric-link Casimir measures the same occupation pattern on these states, hence **the full commuting electric-error algebra obeys the exact Knill–Laflamme condition** for this code.

But the middle **physically gauge-invariant fundamental magnetic** error compresses to `P B_1 P = [[0,1],[1,0]]`—logical `X`. If both the identity and this magnetic error are permitted, the Knill–Laflamme scalarity test fails. This does not invalidate the restricted electric-only code; it demonstrates why local gauge symmetry alone fails to protect it.

## 5. Numerical finite-chain spectrum and new spatial diagnostics

Diagonalizing the published Hamiltonian with NumPy double precision gives the following values, all in **dimensionless source units**, with an additive shift `6N/g²` for `N` plaquettes:

| Source parameter `g²` | First gap, 2 plaquettes | First gap, 3 plaquettes | 3-plaquette ground energy |
|---:|---:|---:|---:|
| 0.5 | 4.618270 | **4.493288** | 30.223535 |
| 1.0 | 5.129113 | **5.108866** | 16.823851 |
| 1.5 | 7.569155 | **7.567657** | 11.647555 |
| 2.0 | 10.264213 | **10.264018** | 8.853848 |

The ground state is **charge-conjugation even and reflection even** at all four sampled couplings. The script also reports a gauge-invariant connected covariance of *face occupation* on the leftmost and rightmost faces (0.00032248 at `g²=1`) and the von Neumann entropy for the left-face versus remaining-faces bipartition (0.00733597 nats at `g²=1`). These characterize correlations **inside this finite gauge-fixed qutrit model**. They are not physically calibrated glueball correlators and the entropy depends on the chosen tensor-product face variables.

The near equality of two- and three-plaquette gaps at larger `g²` is a limited strong-electric-coupling observation, **not** finite-volume extrapolation, an interval proof or a physical mass measurement. Serious glueball spectroscopy requires increasing volumes and link-representation cutoffs, classifying spin/parity/charge-conjugation interpolating operators, and separately controlling lattice spacing, discretization and confinement dynamics.

## 6. Next research steps

1. Continue to four and five faces with **the same published operator convention**, taking care to keep the model's open square-chain boundary. Measure finite-volume trends over several couplings rather than extrapolate from two points.
2. Reproduce the published paper's multi-plaquette numerical spectra (where published values/parameters permit direct comparison), and compare the fermion/static-charge sectors **only after** adding corresponding boundary flux and Gauss sources.
3. Specify fault location and timing, e.g. local electric errors, Hermitian magnetic errors, and correlated adjacent-plaquette faults. Test **finite error spans** and explicit recovery rather than inferring active-code impossibility from the full-algebra commutant.
4. Move to larger two-dimensional spatial geometry and consistent `J^PC` glueball-like interpolators before any actual glueball mass calculation. Do not assert a mass-gap theorem from a finite spectral gap.

## Reproduce

From the repository folder with parent PR #27 / merged PR #28 files available:

```bash
python -m pip install 'sympy==1.14.0' numpy
OPENBLAS_NUM_THREADS=1 python su3_three_plaquette_chain.py
OPENBLAS_NUM_THREADS=1 python su3_three_plaquette_error_algebra.py
python -m py_compile su3_three_plaquette_chain.py su3_three_plaquette_error_algebra.py
```

The chained benchmark imports `su3_published_square_chain_benchmark.py`, which in turn uses the prior `su3_complete_minimal_links.py` and `su3_two_plaquette_spin_network.py`. Exact symbolic assertions and both finite-field ranks passed in the test environment. NumPy energies use double precision without interval enclosure. The two scripts and separately recorded outputs form the reproducible research package.