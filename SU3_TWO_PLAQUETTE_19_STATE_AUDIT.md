# SU(3) adjacent-plaquette expansion: simultaneous flux and shared-link fusion

**Status:** 18 October 2026, exact finite-dimensional symbolic calculation with an illustrative numerical spectrum. Executable: [`su3_two_plaquette_extended.py`](su3_two_plaquette_extended.py); recorded full console output: [`su3_two_plaquette_extended_results.txt`](su3_two_plaquette_extended_results.txt). **Not** a glueball mass, a continuum mass-gap proof, a quantum memory, or a reproduction of a published full link-truncation Hamiltonian.

## 1. Why this experiment matters

The existing seven-state SU(3) model of two neighboring square plaquettes retains only the vacuum and single oriented loops around the left, right, and exterior contours. It checks local Gauss law on all six vertices, but discards the possibility of simultaneous plaquette flux, higher representations, and multiple fusion channels along the shared link. This extension asks whether the restricted-code success and full-error-algebra obstruction survive a **nested, larger gauge-invariant variational basis**.

The physical graph still contains seven links and six vertices (see `SU3_TWO_PLAQUETTE_GAUSS_ERROR_AUDIT.md`). Gauge-fix a maximal tree so the independent oriented face holonomies are `U,V ∈ SU(3)` and both transform by simultaneous conjugation under the residual gauge rotation. Every retained state is a gauge-invariant trace polynomial in the closed based holonomies. In particular, **all six original local Gauss constraints are respected by construction for these Wilson functions**: at each internal vertex, adjacent open-link color transformations cancel, while transformations at the common basepoint conjugate whole words and disappear under the trace. This is an algebraic gauge-invariance argument, **not** a new explicit sparse 19×6×6 generator-action enumeration. The earlier seven-state calculation supplies that independent enumeration for its original states.

## 2. Exactly specified nested basis

Write `χ3(U)=tr(U)`, `χbar3(U)=tr(U†)`, `χ6(U)=(tr U)^2−tr(U†)`, and `χ8(U)=tr U tr(U†)−1`. The latter two identities follow from `3⊗3=6⊕bar3` and `3⊗bar3=1⊕8`.

| Retained dimension | States added in this step |
|---:|---|
| **7** | `1, tr U, tr U†, tr V, tr V†, tr(UV), tr(UV)†` (previous model) |
| **11** | Four simultaneous two-plaquette products `tr U^(±) tr V^(±)` |
| **17** | `χ6(U), χbar6(U), χ6(V), χbar6(V), χ8(U), χ8(V)` |
| **19** | Mixed winding traces `tr(UV†), tr(U†V)` |

The basis is **not orthogonal**. For example,

\[
\langle \mathrm{tr}(UV),\mathrm{tr}U\,\mathrm{tr}V\rangle=\frac13,
\qquad
\langle \mathrm{tr}(UV^\dagger),\mathrm{tr}U\,\mathrm{tr}V^\dagger\rangle=\frac13.
\]

Those are the four and **only** off-diagonal overlap pairs across the 19 retained trace polynomials; their partners are the charge conjugates. The complete Gram matrix has exact rational entries and **rank 19**. Calculations use the Gram matrix explicitly instead of incorrectly treating trace polynomials as orthogonal basis states.

These states resolve the representation channels **on the common link**:

- For `tr U tr V`, the shared link has oppositely oriented color strands, `3⊗bar3=1⊕8`. The normalized exterior contour `tr(UV)` carries the singlet component, and the orthogonal linear combination carries the adjoint `8`.
- For `tr U tr V†` and `tr(UV†)`, the strands are aligned, `3⊗3=bar3⊕6`. Their symmetric and antisymmetric combinations resolve the `6` and `bar3` sectors, respectively.
- The extra individual-face characters involve link labels `6`, `bar6`, or `8` around an entire plaquette. Thus this goes **beyond** a per-link `1,3,bar3` truncation. It does not systematically include *all* admissible intertwiners of a link-truncated theory.

## 3. Exact SU(3) magnetic matrices

The elementary Wilson multiplication operators are `W_L ψ=tr(U)ψ` and `W_R ψ=tr(V)ψ`, followed by orthogonal projection into the chosen trace-polynomial space. Each covariant matrix entry `S^L_{ij}=⟨f_i|tr(U)|f_j⟩` is calculated by **exact Haar integration over SU(3)×SU(3)**. The script implements a small finite-degree SU(3) Weingarten integral, including the determinant-twisted sectors when the number of fundamental and antifundamental factors differs by three. It uses rational numbers; no Monte Carlo or float threshold enters these matrix elements.

For nonorthogonal Gram matrix `G`, the actual projected linear operator is `G⁻¹S^L`, not `S^L`. The script checks Hermitian adjoints with this metric. Its original seven-state Wilson submatrix **exactly** reproduces the five nonzero left-loop transitions in the previous file, including the shared-link amplitude `1/3`.

An important truncation warning persists: **unprojected Wilson multiplication operators commute**, whereas their projected matrices need not. The projected commutator's rank in the stages 7/11/17/19 is respectively **2/0/6/6**. This nonmonotonic behavior is caused by changing which intermediate representation channels remain available. It is not a new non-Abelian commutator between physical Wilson multiplication operators.

## 4. Shared-link electric energy from SU(3) fusion, not guessed coefficients

In a consistent loop orientation, an infinitesimal gauge-electric variation of the common link acts on `U→Ue^{itT}` and `V→e^{-itT}V`. The shared-edge Casimir is therefore

\[
E_s^2=(R_U-L_V)^2.
\]

The six **unshared** links contribute three independent copies of the quadratic Casimir on each face holonomy, so the projected seven-link electric Hamiltonian is

\[
\boxed{D=3C_U+3C_V+E_s^2}.
\]

On a character `χ_(p,q)` of one holonomy, `C=(p²+q²+pq+3p+3q)/3` (hence `C3=4/3`, `C6=10/3`, `C8=3`). On the two-strand shared link, the SU(3) Fierz identity

\[
\sum_{a=1}^{8}T^a_{ij}T^a_{kl}
=\tfrac12(\delta_{il}\delta_{jk}-\tfrac13\delta_{ij}\delta_{kl})
\]

fixes the remaining recoupling coefficients; the verifier confirms **all 81 component identities** with exact symbolic Gell-Mann generators.

For same-oriented face flux, set `A=tr U tr V`, `B=tr(UV)`. Exact shared-electric action gives

\[
E_s^2B=0,\qquad E_s^2A=3A-B.
\]

Thus the Gram-orthogonal sectors have shared-link Casimir `0` (singlet) and `3` (adjoint). For mixed orientations, `A=tr U tr V†` and `B=tr(UV†)` obey

\[
E_s^2A=\tfrac73 A+B,\qquad
E_s^2B=A+\tfrac73 B.
\]

The symmetric/antisymmetric combinations have shared-link eigenvalues **10/3 (`6`)** and **4/3 (`bar3`)**. The exact Gram-adjusted matrix diagonalizations verify both decompositions and nonnegative local-electric energy. All three independently specified local electric operators (`C_U`, `C_V`, shared `E_s²`) commute in this retained representation. No arbitrary energy values were assigned to the added flux configurations.

## 5. Quantum information: robust *conclusion*, not robust *memory*

Let the logical pair be the two opposite orientations of one left-plaquette fundamental loop, `|0_L⟩=tr(U)` and `|1_L⟩=tr(U†)`. Both are orthonormal and eigenvectors with **identical eigenvalues** for all seven local electric-link Casimirs. Therefore their span corrects every operator in the **entire algebra generated by those electric Casimirs** (an exact restricted Knill–Laflamme certificate).

However, the standard physically gauge-invariant Hermitian left Wilson operator `B_L=W_L+W_L†` compresses to

\[
P B_L P=\begin{pmatrix}0&1\\1&0\end{pmatrix}=X_L.
\]

So this candidate code fails Knill–Laflamme when `I` and the magnetic plaquette disturbance are both allowed. The result persists **in all four nested basis sizes**.

A stronger, basis-wide statement comes from exact **rational** commutant calculations:

| Stage | Algebra from `D,B_L,B_R` | Complex commutant | Algebra from `D,W_L,W_R,W_L†,W_R†` | Commutant |
|---:|---|---:|---|---:|
| 7 | `M4⊕M3` | 2 | `M7` | 1 |
| 11 | `M6⊕M5` | 2 | `M11` | 1 |
| 17 | `M10⊕M7` | 2 | `M17` | 1 |
| 19 | `M11⊕M8` | 2 | `M19` | 1 |

Charge-conjugation parity is exact for the real/Hermitian Wilson terms, with two **inequivalent, irreducible** matrix blocks. Therefore **no multiplicity-two noiseless subsystem** is protected against arbitrary operations in that *complete generated algebra*. The complex fundamental plaquette multipliers are also individually gauge invariant and generate the full matrix algebra in every tested finite truncation. These conclusions are **not** no-go theorems for all active error correction against bounded fault sets, other larger encodings, or continuum SU(3).

## 6. Finite-volume spectral comparison (illustrative only)

Use the same dimensionless normalization for every nested basis:

\[
H=D-\tfrac12(B_L+B_R).
\]

Diagonalize its **generalized Hermitian eigenproblem** `H_cov v=E G v`, rather than treating the nonorthogonal functions as orthonormal. The exact matrix elements are rational, but the eigenvalues are obtained with double-precision SciPy diagonalization:

| Retained states | Ground energy | First excitation gap |
|---:|---:|---:|
| 7 | −0.1989835531 | 5.0323168864 |
| 11 | −0.2026695952 | 5.0307804479 |
| 17 | −0.2051695463 | 4.9318717963 |
| 19 | −0.2051883138 | 4.9318706386 |

The monotone decrease in the ground energy is the expected variational result for nested trial subspaces. The near agreement of the 17- and 19-state **gaps** is only an observation for one arbitrary electric/magnetic ratio; eigenvalue-gap convergence is not ensured by the variational principle. No lattice-spacing or spatial-volume extrapolation, `J^PC` glueball spectroscopy, or physical calibration is included. These values are **not glueball masses** and say nothing rigorous about a continuum Yang–Mills mass gap.

## 7. Reproduce and evaluate

```bash
python -m pip install sympy==1.14.0 numpy scipy
OPENBLAS_NUM_THREADS=1 python su3_two_plaquette_extended.py
python -m py_compile su3_two_plaquette_extended.py
```

Exact assertions: SU(3) generator Fierz completeness, determinant-charge Haar moments, Gram ranks, shared-link Casimir eigenvalues, positive electric operator, mutual electric commutation, magnetic matrix identities, charge-conjugation symmetry and parity dimensions, commutants, and exact KL pass/failure. Every exact assertion passed locally using SymPy 1.14. Spectral diagonalizations use ordinary double precision and are **not interval-certified**.

**Evaluation criteria before claiming additional physics:**

1. Compare this **selected-loop** 19-state projection with the complete gauge-singlet Hilbert space for a chosen finite set of per-link irreps, especially all vertex intertwiner multiplicities; report missed states and Wilson-operator leakage.
2. Build the actual minimally truncated square-chain Hamiltonian from Chen, Müller and Yao on the **same geometry and normalization**, then quantify the difference. That paper retains per-link `1,3,bar3`, unlike this variational basis, which includes some `6`/`8` states but omits many others.
3. Measure how the spectral and algebraic conclusions vary with magnetic-to-electric coupling ratio, boundary conditions, and additional states. Do not extrapolate one toy gap.
4. Only then study gauge-invariant connected-loop correlation functions on **larger lattices**, lattice parity, charge conjugation and rotation-sector projection; a two-plaquette box cannot support serious physical glueball spin identification.
5. For a hardware QEC claim, specify actual local noise Kraus operators, correlated magnetic faults, syndrome extraction and recovery; the no-noiseless-subsystem result for a complete algebra is not automatically an active-code no-go.

## Primary references

- Kogut and Susskind, *Hamiltonian formulation of Wilson's lattice gauge theories*, Phys. Rev. D **11**, 395 (1975): https://doi.org/10.1103/PhysRevD.11.395.
- Chen, Müller and Yao, *Minimally truncated SU(3) lattice gauge theory and string tension*, Phys. Rev. D **114**, 074502 (published 5 October 2026): https://doi.org/10.1103/1qt7-3ywl; arXiv:2601.10065. This source presents a **different, complete minimally link-truncated** operator construction for small lattices, **not the present selective character basis**.
- Mathur and Sreeraj, *Canonical transformations and loop formulation of SU(N) lattice gauge theories*, Phys. Rev. D **92**, 125018 (2015): https://doi.org/10.1103/PhysRevD.92.125018.
- Knill and Laflamme, *Theory of quantum error-correcting codes*, Phys. Rev. A **55**, 900 (1997): https://doi.org/10.1103/PhysRevA.55.900.
