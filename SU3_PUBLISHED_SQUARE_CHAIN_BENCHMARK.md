# Published two-plaquette SU(3) Hamiltonian: exact matrix-by-matrix benchmark

**Status: 18 October 2026; exact nine-dimensional source-to-code comparison.** Executable: [`su3_published_square_chain_benchmark.py`](su3_published_square_chain_benchmark.py). Parent PR #27 supplies the complete two-face microscopic SU(3) link model. **This is a finite-lattice operator benchmark, not a physical glueball mass, a continuum theorem, or a new quantum-error-correcting code.**

## Primary-source comparison

Vincent Chen, Berndt Müller, and Xiaojun Yao, *Minimally truncated SU(3) lattice gauge theory and string tension*, **Physical Review D 114**, 074502 (published 5 October 2026), [DOI 10.1103/1qt7-3ywl](https://doi.org/10.1103/1qt7-3ywl), [publisher's PDF](https://journals.aps.org/prd/pdf/10.1103/1qt7-3ywl). The comparison uses the **published version**: Section IIIA–B, equations (3)–(5), (10)–(11), and Appendix B **Table VI** (PDF page 14). In the earlier [arXiv:2601.10065v1](https://arxiv.org/pdf/2601.10065v1), the equivalent Appendix B table is numbered **Table IV**. The changes in table numbering are intentional: do not cite the preprint number as the final journal number.

The published square-chain formulation retains link representations `1`, `3`, `anti-3` and uses plaquette electric-flux qutrit numbers `0,1,2` (where `2 ≡ −1 mod 3`). Each link has the representation difference of its adjacent plaquette fluxes; exterior plaquettes and the top/bottom background have flux zero here. The previously developed **two adjacent square faces, seven links, six vertices** correspond to this *zero-exterior-flux two-face case*. Boundary vertices of degree two may be regarded as trivalent vertices with an externally attached trivial representation. **This equivalence depends on the fixed external singlet boundary condition**; it has not been tested for flux sectors, other boundary conditions, or larger lattices.

## Exact state mapping

Write the pair of dynamical face labels as `(qL,qR)` in `Z3²`. For the previous repo's state order:

| Spin-network label | Paper's `(qL,qR)` | Electric sum of link Casimirs |
|---|---|---|
| `vac` | `(0,0)` | `0` |
| `L+`, `L-` | `(1,0)`, `(2,0)` | `16/3` |
| `R+`, `R-` | `(0,1)`, `(0,2)` | `16/3` |
| `O+`, `O-` | `(1,1)`, `(2,2)` | `8` |
| `B+`, `B-` | `(1,2)`, `(2,1)` | `28/3` |

Oriented seven-link labels (edge indices matching the repo's `EDGES` ordering) are

\[
\ell(q_L,q_R)=(q_L,q_R,-q_L,-q_R,-q_L,q_L-q_R,q_R)\pmod3.
\]

All nine images are distinct and reproduce **all nine** link-representation configurations in the parent verifier's exhaustive `3^7=2187` irrep-assignment enumeration. In particular the two `B±` junction states are not invented excess states: they are the valid opposite-flux pairs `(1,2)` and `(2,1)`. This is a bijection, not a heuristic dimension match.

## Independent, published magnetic transitions

For two faces the left plaquette has a frozen-left neighbor `0` and active-right neighbor `n=qR`; the right plaquette has `n=qL` and frozen-right neighbor `0`. By the reflection symmetry stated in the published appendix, only the **three** boundary-zero blocks of Table VI are required:

\[
M_0=\begin{pmatrix}0&1&1\\1&0&1\\1&1&0\end{pmatrix},\quad
M_1=\begin{pmatrix}0&\frac13&\frac1{\sqrt3}\\\frac13&0&\frac1{\sqrt3}\\\frac1{\sqrt3}&\frac1{\sqrt3}&0\end{pmatrix},
\]

\[
M_2=\begin{pmatrix}0&\frac1{\sqrt3}&\frac13\\\frac1{\sqrt3}&0&\frac1{\sqrt3}\\\frac13&\frac1{\sqrt3}&0\end{pmatrix}.
\]

These matrices are **hardcoded transcriptions of the published table**; they are not fitted to our previous matrices. The directed fundamental plaquette operator shifts `q→q+1 (mod 3)` as in the paper's equation (4); its amplitude is the corresponding element of `M_n`. Summing the directed shift and its adjoint reproduces `M_n` on each fixed-neighbor sector.

The independent benchmark constructs both **9×9 directed Wilson operators** using only these published data and the qutrit face labels. It then constructs the previous link-spin-network operators from **independently evaluated SU(3) Haar integrals**, correctly accounting for the raw B± states' squared norm `4/3` by using normalized `|B±⟩ = sqrt(3)/2 |B±raw⟩`. Both operators agree **entry by entry**:

\[
W_L^{\rm publication}=W_L^{\rm link},\qquad
W_R^{\rm publication}=W_R^{\rm link}.
\]

All `81+81` directed entries agree exactly in SymPy; likewise the two Hermitian `W+W†` matrices. The only nonzero directed amplitudes are `1`, `1/3`, and `1/sqrt(3)`. This comparison is not merely spectral similarity and is not sensitive to numerical floating-point thresholds.

## Electric term and full coupling convention

All seven links have quadratic electric Casimir `0` for `1` or `4/3` for `3/anti-3`. The paper's electric operator (per-link Casimir sum) is

\[
D=\mathrm{diag}(0,16/3,16/3,16/3,16/3,8,8,28/3,28/3).
\]

The benchmark constructs this diagonal *independently* from the published qutrit link-difference formula, then compares it exactly with all nine link configurations of the prior microscopic enumeration.

The published square-lattice Hamiltonian, its Eq. (3) with `n_l=2` and `Σ=1` and **two plaquettes**, is

\[
H_{\mathrm{PRD}}(g^2)=g^2D+\frac1{g^2}\Big(12I_9-B_L-B_R\Big),
\quad B_\alpha=W_\alpha+W_\alpha^\dagger.
\]

The term `12 I/g²` is an energy-independent constant, which the paper explicitly drops in subsequent calculations. Prior repo scans instead used `H_prior(c)=D-(c/2)(B_L+B_R)`. These are **exactly related**, provided

\[
c=\frac{2}{(g^2)^2},\qquad
H_{\mathrm{PRD}}(g^2)=g^2H_{\mathrm{prior}}(c)+\frac{12}{g^2}I_9.
\]

The verifier checks all 81 matrix entries of this identity for `g² ∈ {1/2,1,3/2,2}` using symbolic rational arithmetic. It also diagonalizes the source-normalized 9×9 Hamiltonian numerically at `g²=1` and `2` as a **secondary** check, not an exact eigenvalue certificate.

| Source `g²` | First energy | First gap |
|---:|---:|---:|
| 1 | 11.2042206179 | 5.1291127154 |
| 2 | 5.9024532821 | 10.2642133846 |

All energies are dimensionless and depend on the specified finite box. The additive `12/g²` shift changes ground energies but not gaps; multiplying by `g²` changes both gaps and nonconstant energies. Do not compare earlier `c`-scan numbers to paper `g²` by setting their numerical parameters equal.

## Scientific conclusion and next steps

**We now have an exact published-operator benchmark for the fully enumerated nine-state two-plaquette model in the zero-external-flux sector.** The underlying Gauss-singlet completeness and Haar/fusion calculations are independent of the benchmark's published transition table. This materially strengthens confidence in our finite-link implementation. The old seven-state projection remains physically incomplete (it omits nonzero Wilson transitions), and the 19-state selected higher-irrep basis uses a different truncation and is not part of this equivalence.

The previously demonstrated conditional **electric-only code**, and failure under the physically allowed **fundamental magnetic Wilson loop**, apply to this benchmarked Hamiltonian at the stated cutoff. This is still an operator-algebra result, not an implemented error-correction recovery or a new Standard Model code.

Next prioritize **three and more square plaquettes** in the same published `Z3` face-qutrit formulation, with the published Table VI for all relevant `(PL,PR)` sectors and fixed external flux boundary, or build the complete link-intertwiner version independently and compare again. Add electric/magnetic coupling sweeps and spatial-correlation diagnostics before making *any* glueball or continuum mass-gap statements. The published benchmark does not establish equivalence for arbitrary boundary flux, static quarks, mixed representation cutoffs, or physical fault rates.

## Reproduce

```bash
python -m pip install sympy numpy
OPENBLAS_NUM_THREADS=1 python su3_published_square_chain_benchmark.py
python -m py_compile su3_published_square_chain_benchmark.py
```

The verifier imports the parent PR #27's `su3_complete_minimal_links.py` and its `su3_two_plaquette_spin_network.py` dependency. **Merge PR #27 before separately merging a dependent follow-up PR**, or review this as a stacked PR. All deterministic operator-identity checks use exact symbolic arithmetic; spectra use NumPy double precision.
