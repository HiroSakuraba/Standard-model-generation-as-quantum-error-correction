# Four- and five-plaquette SU(3) square chain: spectra and bounded-error KL tests

**Research status — 18 October 2026.** Companion exact and numerical verifier: [`su3_four_five_plaquette.py`](su3_four_five_plaquette.py). Run output: [`su3_four_five_plaquette_results.txt`](su3_four_five_plaquette_results.txt). These results concern a *finite minimally link-truncated pure SU(3) gauge simulator*, not observed glueballs, the continuum mass gap, full Standard Model anomalies, or a physical fault-tolerant memory.

## Source and test scope

We extend the published Hamiltonian of Vincent Chen, Berndt Müller and Xiaojun Yao, *Minimally truncated SU(3) lattice gauge theory and string tension*, **Physical Review D 114**, 074502 (5 October 2026), DOI [10.1103/1qt7-3ywl](https://doi.org/10.1103/1qt7-3ywl), square-chain equations (3), (5), (10), (11), and Appendix B **Table VI** (journal version). The repository's [`su3_three_plaquette_chain.py`](su3_three_plaquette_chain.py) contains *independent* symbolic transcriptions of the nine 3×3 magnetic blocks and the trivalent-vertex factors that generate them; its parent tests also reproduce the published **two-plaquette 9×9** Hamiltonian element by element. The present program **imports those published blocks**, checks their independent vertex decomposition, and uses precisely the same open, zero-exterior-flux convention.

Each of N faces carries a qutrit flux label \(q_i\in\{0,1,2\}\simeq \mathbb Z_3\); 1 and 2 indicate the fundamental and antifundamental SU(3) link representations, respectively. The complete gauge-singlet face model has \(3^N\) states and \(3N+1\) spatial links for an open strip. The links' representation labels are determined by top \(q_i\), bottom \(-q_i\), and vertical differences \((-q_0,q_0-q_1,\ldots,q_{N-2}-q_{N-1},q_{N-1})\), all modulo 3. Each nontrivial link has quadratic electric Casimir \(4/3\). Thus the Hamiltonian is

\[
H_N(g^2)=g^2 D_N+\frac{6N I-\sum_{i=0}^{N-1}(W_i+W_i^\dagger)}{g^2},
\quad D_N(q)=\tfrac43\#\{\text{nontrivial link representations}\}.
\]

Each directed Wilson operator shifts one plaquette \(q_i\mapsto q_i+1\) with the exact coefficient in the published \(M_{(q_{i-1},q_{i+1})}\) Table VI block, using exterior neighbors zero. These projected Wilson operators are not unitary in general. All magnetic matrix elements are real in the chosen basis; the Hermitian \(B_i=W_i+W_i^\dagger\) is a standard gauge-invariant local interaction.

## Verification

The exact script tests, for N=2,3,4,5:

- all \(3^N\) unique qutrit-to-link assignments and the \(3N+1\) link count;
- every one-site magnetic transition under global charge conjugation \(C:q_i\mapsto -q_i\) and chain reflection \(R:q_i\mapsto q_{N-1-i}\), using symbolic fractions and \(\sqrt3\), and all electric energies under conjugation;
- complete **entrywise** agreement of its new N=2 and N=3 operators with the existing exact parent Hamiltonians;
- an independently assembled numerical Hamiltonian satisfying \([H,C]=[H,R]=0\), orthogonal eigenvector residuals better than \(2\times10^{-11}\), and the simultaneous C/R sector sizes.

At N=4, the four (C,R) sector dimensions \((+,+),(+,-),(-,+),(-,-)\) are **25,16,20,20**, totaling 81. At N=5 they are **70,52,65,56**, totaling 243. These sector counts follow exactly from the finite permutation involutions; their numerical use only verifies the assembled matrices. The computed ground state has \(C=R=+1\) at all sampled couplings; the first excited state has \(C=+1\) and R parity \(-1\) for even N and \(+1\) for odd N among N=2,...,5 at these parameters.

## First excitation gaps (double precision, dimensionless source units)

| Source \(g^2\) | N=2 (9 states) | N=3 (27) | N=4 (81) | N=5 (243) |
|---:|---:|---:|---:|---:|
| 0.5 | 4.618270 | 4.493288 | 4.225976 | 4.132531 |
| 1.0 | 5.129113 | 5.108866 | 5.098354 | **5.091694** |
| 1.5 | 7.569155 | 7.567657 | 7.566912 | 7.566458 |
| 2.0 | 10.264213 | 10.264018 | 10.263920 | 10.263861 |

At \(g^2=1\), the connected endpoint-occupation covariance \(\langle n_0n_{N-1}\rangle-\langle n_0\rangle\langle n_{N-1}\rangle\) is approximately \(-3.2442\times10^{-3}\) (N=2), \(+3.2248\times10^{-4}\) (N=3), \(-2.9694\times10^{-5}\) (N=4), and \(+2.9612\times10^{-6}\) (N=5). The sign alternation and decrease in magnitude are observations for one geometry, cutoff, and coupling; they do **not** establish a continuum correlation length or universal mass gap. A finite open chain has neither the requisite 3D spatial rotation spectrum for glueball spin \(J\) nor systematic lattice-spacing and representation-cutoff extrapolations. Numerical eigenvalues are double precision, not interval-certified.

## A mathematically exact **bounded-error** two-dimensional code at N=5

Consider normalized orthogonal face basis states

\[
|0_L\rangle=|0,0,0,0,1\rangle,\qquad
|1_L\rangle=|0,0,0,0,2\rangle.
\]

These differ only by reversing the fundamental versus antifundamental flux on the **rightmost** plaquette. Their 16 link *electric Casimirs* all have identical eigenvalues: each only tests the link's representation magnitude, not its orientation.

Define the **single-event linear noise space** (not its algebraic closure)

\[
\mathcal E=\operatorname{span}\{I, E_0^2,\ldots,E_{15}^2, B_0,B_1,B_2,B_3\},
\quad B_i=W_i+W_i^\dagger.
\]

The script forms each \(E_a|0_L\rangle\), \(E_a|1_L\rangle\) as an **exact sparse qutrit vector**, computes **all 441** overlaps, and verifies the full Knill–Laflamme scalarity identity

\[
\langle p_L|E_a^\dagger E_b|q_L\rangle=\lambda_{ab}\delta_{pq}
\]

for every one of the **21 listed elementary operators**. The error Gram matrix is Hermitian and has rank **5**: the nominal 16 electric operators are highly degenerate on the code, so 21 is not the number of independent nontrivial syndrome outcomes. The certificate ensures mathematical existence of a recovery for arbitrary noise channels with Kraus operators in this particular span. **It is not a distance-five code, a compiled recovery circuit, or protection against unrestricted native Hamiltonian evolution.** All five magnetic terms occur in the physical source Hamiltonian, and the last one is deliberately excluded from \(\mathcal E\).

Two **gauge-invariant counterexamples** delimit the success:

1. **Direct logical flip at the encoded face.** \(P B_4 P=X_L\), so permitting \(B_4\) alongside identity fails KL at once. The native Hamiltonian includes \(B_4\).
2. **Directional error one face closer.** Although the Hermitian sum \(B_3\) passes KL, the *directed* fundamental Wilson operator \(W_3\) is a valid gauge-invariant non-Hermitian error operator (e.g. as a candidate jump/Kraus operator, with appropriate channel normalization). Its error weight already distinguishes the two codewords:

\[
P W_3^\dagger W_3 P=\operatorname{diag}(1/9,1/3),
\]

so even the two-operator set \(\{I,W_3\}\) fails KL. The difference in orientation-sensitive recoupling amplitudes is physically meaningful within the published truncated operator model. This shows why grouping directed magnetic faults into one Hermitian sum is a **substantive noise assumption**, not a harmless change of notation.

This conditional KL success is essentially spatial separation, not intrinsic protection from SU(3) Gauss constraints. No explicit unitary or dissipative decoder, ancilla construction, syndrome measurement, recovery time, hardware noise distribution, fault-tolerance threshold or code lifetime is supplied. We did **not** compute the complete operator-algebra commutant at N=4 or N=5; previous N=3 complete-algebra obstructions should not be silently extrapolated as verified here.

## Reproduction

Run from this repo branch after the three-plaquette source files have been merged or checked out:

```bash
python -m pip install sympy==1.14.0 numpy
OPENBLAS_NUM_THREADS=1 python su3_four_five_plaquette.py
python -m py_compile su3_four_five_plaquette.py
```

The new script imports only the existing `su3_three_plaquette_chain.py` source (which contains the independent published Table VI and Table I transcriptions). The original 3-plaquette script itself imports the two-plaquette benchmark *inside* its optional test path; it is not required for the new N=4,5 run. All symbolic assertions and numerical checks passed locally on 18 October 2026 with SymPy 1.14 and NumPy 2.3.5. The recorded console output includes all sampled energies, the exact KL witnesses, and machine-readable JSON.

## Next research milestone

A physically substantive next step is to replace the custom 21-generator fault span with a stated **probabilistic local jump/Kraus model**, incorporating both Wilson orientations, correlated adjacent-plaquette faults, and time evolution by the *native* Hamiltonian during storage/recovery. Compare optimal approximate recovery with a constructible gauge-invariant decoder. Separately, enlarge the spatial lattice beyond an open one-plaquette-wide strip and study genuine glueball interpolators and discrete spatial symmetry representations. Do not equate a finite-chain energy gap with a glueball mass or a Yang–Mills mass-gap theorem.
