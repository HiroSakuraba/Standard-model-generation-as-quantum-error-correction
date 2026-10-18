# SU(2) native-vertex exchange audit: Gauss protection versus logical protection

**18 October 2026 · finite four-color-spin study · exact symbolic + numerical crosschecks**  
Executable: [`su2_native_vertex_controls.py`](su2_native_vertex_controls.py). This is an extension to the gauge-error-algebra work, not a Standard Model quantum code, a fabricated simulator, or a validated experiment.

## 1. Why change the graph?

The prior eight-spin star graph was a **noise model chosen to preserve a correctable multiplicity**; it was never mapped to the required controls of a physical gauge simulator. In contrast, Halimeh, Homeier, Bohrdt and Grusdt, *Spin Exchange-Enabled Quantum Simulator for Large-Scale Non-Abelian Gauge Theories*, **PRX Quantum 5**, 030358 (2024), DOI [10.1103/PRXQuantum.5.030358](https://doi.org/10.1103/PRXQuantum.5.030358), propose strong **equal intravertex antiferromagnetic Heisenberg interactions** to enforce the SU(2) color-singlet Gauss-law condition energetically. Their proposed ultracold-molecule and ultracold-atom realizations are physical proposals, not demonstrations of the independent high-bandwidth coupler settings tested below.

A honeycomb vertex contains a matter mode and up to three adjacent rishon modes. This study conditions on **one matter spin half and three occupied local rishon color spins**, with one rishon per link conserved. That gives one four-spin vertex. The full dynamical matter/link Hilbert space, other vertex orientations, and additional density/Hubbard interactions of the cited proposal are **not simulated**. Even the four-spin sector need not be dynamically isolated in the device.

## 2. Exact local gauge and logical Hilbert spaces

For four co-located spin halves, let \(\mathbf G=\sum_{i=0}^{3}\mathbf S_i\), \(h_{ij}=\mathbf S_i\cdot\mathbf S_j=(\mathrm{SWAP}_{ij}-I/2)/2\), and \(P\) project onto the two-dimensional total-spin-zero subspace. The full 16-dimensional spin Hilbert space decomposes into two spin-0 irreps, three spin-1 irreps (9 states), and one spin-2 irrep (5 states). The verifier constructs the entire 16×16 exchange and Gauss matrices with exact rational entries and explicit normalized singlets, and tests \([h_{ij},G^z]=[h_{ij},G^+]=0\) on every pair.

The complete six-edge (\(K_4\)) uniform antiferromagnetic Hamiltonian is

\[
H_{\rm G}=J\sum_{i<j}h_{ij}=\frac J2(\mathbf G^2-3I),\quad J>0.
\]

The singlets have energy \(-3J/2\), the spin-1 states \(-J/2\), and the spin-2 states \(3J/2\). Hence the exact Gauss excitation gap is **\(J\)**.

The **three-edge star alone**, \(J(h_{01}+h_{02}+h_{03})\), is *not* an equivalent Gauss-penalty Hamiltonian: it has spin-one states at \(-5J/4\), below the singlets at \(-3J/4\). Thus importing the eight-spin star result directly as the **native gauge-protection** mechanism is invalid. A separate protected code under a star-shaped **error** algebra is still mathematically valid, but the star alone does not enforce the source's four-spin Gauss-law ground sector.

For the full 16-state four-spin space, the identity and the six distinct pair exchanges have exact linear rank **7**. Consequently, implementing an operator proportional to \(\mathbf G^2\) through *only these pair Heisenberg terms and a scalar offset* requires equal pair coefficients. This is stronger than simply balancing logical energies inside the singlet sector.

## 3. Equal Gauss protection does not remove gauge-invariant logical errors

Choose an orthonormal singlet basis where the exchange SWAP(0,1) projects to \(-Z_L\). The physical pair exchange projects exactly to

\[
P h_{01}P = -\frac14 I_L-\frac12 Z_L
=\begin{pmatrix}-3/4&0\\0&1/4\end{pmatrix}.
\]

A coupling error \(\delta h_{01}\) causes an **exact logical splitting \(|\delta|\)**, independent of the ideal protection strength \(J\), since \(H_{\rm G}\) acts as an identity multiple within the singlets and commutes with all the exchanges. Moreover, \(\delta h_{01}\) commutes with all Gauss generators, so **ideal Gauss-charge monitoring sees no violation**.

For arbitrary six Heisenberg couplings, define opposite-edge sums

\[
A=J_{01}+J_{23},\quad B=J_{02}+J_{13},\quad C=J_{03}+J_{12}.
\]

The exact logical Hamiltonian is

\[
H_L=-\frac{A+B+C}{4}I_L+\frac{-2A+B+C}{4}Z_L+\frac{\sqrt3(C-B)}4 X_L.
\]

The singlet energies are degenerate if and only if **\(A=B=C\)**. Thus *correlated* calibration errors satisfying those two linear conditions do not disturb the logical energy splitting even when individual pairs differ. But such a balanced distortion need not realize the exact full-space Gauss penalty. These statements replicate/extend the prior exact six-coupling exchange calculation, not claim an original representation-theoretic theorem.

## 4. A specified stochastic error channel with no Gauss syndrome

Model a single unknown run-to-run coupling error taking either sign \(\pm\delta\) with probability one half, constant within a trial. The stored logical qubit experiences

\[
\mathcal N_{\delta,t}(\rho)=\frac12\sum_{s=\pm1} e^{-it(H_G+s\delta h_{01})}\rho\,e^{it(H_G+s\delta h_{01})}.
\]

The code remains entirely gauge singlet, and the exact **uncorrected entanglement fidelity** to the identity is

\[
F_e^{\rm bare}=\frac{1+\cos(\delta t)}2=\cos^2(\delta t/2).
\]

This follows from the logical phase difference \(\delta t\), without perturbation theory, and is **independent of \(J\)**. Unknown signs make recovery by tracking one deterministic rotation impossible without further information, though we do **not** claim to have optimized over every possible quantum recovery map.

## 5. An exchange-only refocusing control, under additional assumptions

The four-spin singlet has an exact gauge-invariant logical \(X\) generator,

\[
\widetilde X=\frac{2}{\sqrt3}(h_{03}-h_{02}),\qquad
P\widetilde X P=X_L,
\]

and the 16×16 symbolic verifier confirms that the operator **preserves** the two singlets, not just that its compression looks like \(X\). A controlled pulse \(H_{\rm pulse}=J_b\sum h_{ij}+\Omega\widetilde X\) therefore induces \(\Omega X_L\) up to a global phase, while commuting with local SU(2) Gauss generators. Choosing \(J_b>2\Omega/\sqrt3\) makes **every individual pair coupling positive**:

\[
J_{02}=J_b-2\Omega/\sqrt3,\quad J_{03}=J_b+2\Omega/\sqrt3,
\]

with the other four bonds at \(J_b\). The calculation uses \(J_b=12, \Omega=5\) in dimensionless units (couplings 6.2265, 17.7735, and 12 for the remaining four). A nominal logical \(X\) pulse has duration \(\tau=\pi/(2\Omega)=\pi/10\approx0.314159\). The two-pulse echo is \(U_pU_fU_pU_f\), where each free interval lasts \(t/2\) and both finite pulses feel the same unknown error \(\pm\delta h_{01}\). Ideal instantaneous \(X\) pulses would cancel quasi-static \(Z_L\) errors exactly. **Positive exchange strength alone does not establish hardware feasibility:** the required *independent and fast* positive pair-coupler modulations, on/off switching, uniform background calibration, and preservation of link/rishon constraints have not been demonstrated in the cited simulator.

### Finite-pulse benchmark

All entries are computed with NumPy double-precision diagonalization of the **exact two-dimensional singlet block**, not an approximate Schrieffer–Wolff gate. The verifier independently confirms the 16-dimensional full-spin pulse action equals the projected two-dimensional action on the code to numerical precision. \(t=3\) is free storage time; echo pulses add \(2\tau\approx0.6283\) for total protocol time 3.6283. In these tests no extra idle-bath noise accrues by waiting the additional pulse time.

| \(\delta\) | Bare \(F_e\), t=3 | Finite echo, \(\pm\delta h_{01}\) static | Finite echo, \(\pm\delta h_{01}\) independently flipped by half | Finite echo, \(\pm\delta h_{02}\) static |
|---:|---:|---:|---:|---:|
| 0.00 | 1.000000000 | 1.000000000 | 1.000000000 | 1.000000000 |
| 0.02 | 0.999100270 | 0.999999996 | 0.999542138 | 0.999013000 |
| 0.05 | 0.994385539 | 0.999999828 | 0.997142896 | 0.993842794 |
| 0.10 | 0.977668245 | 0.999997260 | 0.988636136 | 0.975535431 |
| 0.20 | 0.912667807 | 0.999956473 | 0.955565359 | 0.904729496 |

When the error flips independently between the two halves, echo cancellation is incomplete. The \(h_{02}\) fault contains an \(X_L\) component parallel to the control axis, so the same echo does not cancel it and is **slightly worse** for this selected case. The table is not a benchmark of arbitrary noise or an experimentally fitted process.

## 6. What "gauge-compatible recovery" would require

The pulse above is a **unitary open-loop mitigation control**, not recovery from an arbitrary unknown logical error. In the two-dimensional singlet space, local SU(2)-invariant exchange operators realize the full matrix algebra \(M_2\). One unknown exchange error is already a logical observable, and no recovery on this unencoded two-dimensional qubit can correct \(\{I,\mathrm{SWAP}_{01}\}\) while preserving both arbitrary logical basis states. A Gauss-law syndrome measurement always gives the same outcome. Measuring the pair SWAP eigenvalue instead reveals \(Z_L\) (the protected information) and can destroy an unknown superposition; it is **not** an ordinary stabilizer syndrome.

The existing [five-block SU(2) recovery](FIVE_BLOCK_SU2_RECOVERY.md) provides **mathematical redundancy**, but its multi-module stabilizer readout still needs a source-compatible elementary interaction scheme. The existing [two-body mediated gate study](TWO_BODY_SU2_MEDIATED_GATE.md) demonstrates an engineered perturbative route with leakage; it does not validate the real ultracold-molecule platform's ancilla, pulse, measurement or full error budget.

## 7. Research decision

**Supported:** exact local four-spin model and exact diagnosis of why Casimir Gauss protection does not correct intra-singlet exchange noise; first-order sensitivity and correlated balance criterion; an engineered positive-Heisenberg exchange control that ideally refocuses one *specified* slow error; full-spin versus logical-block pulse agreement.

**Not supported:** physically implemented independent six-coupler controls; dynamical matter/rishon constraints outside the selected orientation; stochastic noise spectrum measured in any device; realistic pulse switching, T1/T2, gate readout; multi-module recovery circuit; Standard Model quantum error correction or hardware benefit.

**Next implementable comparison:** import an explicit site-resolved proposal's microscopic Hamiltonian including occupation projectors and per-link tunneling; check Gauss-constraint protection and leakage at finite coupling and finite drive frequency, while auditing whether the six independent exchange controls used by the echo are actually available. For operational correction, compile and test one full five-module stabilizer measurement from elementary gates with explicit noisy ancilla and readout rather than treating logical Pauli checks as free.

## Reproduce

```bash
python -m pip install sympy==1.14.0 numpy
OPENBLAS_NUM_THREADS=1 python su2_native_vertex_controls.py
```

Outputs all exact assertion passes, a finite-pulse logical-fidelity table, and a JSON results line. Calculations are small (16-dimensional exact symbolic, 2D numerical and 16D numerical crosschecks); no interval verification claim is made.

**Primary source:** Jad C. Halimeh, Lukas Homeier, Annabelle Bohrdt, Fabian Grusdt, *Spin Exchange-Enabled Quantum Simulator for Large-Scale Non-Abelian Gauge Theories*, **PRX Quantum 5**, 030358 (2024). [Journal](https://doi.org/10.1103/PRXQuantum.5.030358) · [arXiv:2305.06373](https://arxiv.org/abs/2305.06373).
