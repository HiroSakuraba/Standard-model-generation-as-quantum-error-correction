# Smooth exchange pulses for SU(2) singlet gates: an eight-spin numerical optimization

**Research status — 18 October 2026.** Reproducible numerical pulse shaping, *not* an implemented quantum processor, experimentally validated fault-tolerant gate, or a Standard Model anomaly-derived quantum code. The sole new control family in this milestone is a synchronized smooth envelope of two specified on-site Heisenberg exchanges. Source: [\`su2_shaped_pulse.py\`](su2_shaped_pulse.py). The original square-pulse interaction and exact perturbative analysis are in [\`TWO_BODY_SU2_MEDIATED_GATE.md\`](TWO_BODY_SU2_MEDIATED_GATE.md).

## 1. What was optimized

The microscopic model contains two four-spin SU(2) singlet registers, one playing the *encoded ancilla* and the other the *encoded data qubit*. All eight constituent spins are assumed co-located at a single gauge vertex. Each register's spin-zero subspace has dimension two. Their product defines a four-dimensional **logical** subspace of the entire \(2^8=256\)-dimensional Hilbert space.

Only **pairwise SU(2)-invariant Heisenberg exchange** plus a block-spin Casimir penalty is used:

\[
H(t)=\Delta(\mathbf S_A^2+\mathbf S_D^2)
+J_1(t)\mathbf S_{A_0}\cdot\mathbf S_{D_0}
+J_2(t)\mathbf S_{A_1}\cdot\mathbf S_{D_1}.
\]

The first study used abrupt, constant \(J_1=J_2=J\). Its second-order effective Hamiltonian produces a \(Z_A Z_D\) phase with nominal \(\pi/4\) duration \(t=6\pi\Delta/J^2\). The abrupt finite pulse suffers appreciable **leakage out of the product-of-two-singlets code subspace**, even though the total SU(2) Gauss charge is conserved.

We instead fix \(\hbar=\Delta=1\), choose \(J_{\max}=1\), and apply

\[
\boxed{J_1(t)=J_2(t)=J_{\max}\sin^2(\pi t/T),\quad 0\le t\le T.}
\]

This is a **positive, simultaneous, finite-bandwidth-in-principle envelope** that turns both couplings on and off continuously. The implementation assumes its analog waveform is perfectly available; it neither simulates a limited device bandwidth nor provides a digital hardware compilation.

We perform a bounded one-dimensional numerical search over \(30\le T\Delta\le46\) to minimize the code-compressed operator-norm deviation from a specified **phase-calibrated controlled-Z (CZ) representative**, keeping \(J_{\max}\) fixed. This is **local optimization in this declared family**, not a global optimal-control theorem, and the resulting gate is only *locally equivalent* to CZ after individual logical phase adjustments.

## 2. Exact invariant-sector reduction, numerical propagation and metric

The full eight-spin Hilbert space has dimension 256. Since the physical Hamiltonian is a global SU(2) scalar, any initial global spin singlet stays in the \(S=0\) invariant sector. The program constructs that sector instead of assuming its basis:

1. Construct all 28 microscopic pair-exchange matrices in the 256-dimensional computational spin basis, and the full \(S_{\mathrm{tot}}^2\) Casimir.
2. Restrict to total \(S_z=0\) (70 states), diagonalize \(S_{\mathrm{tot}}^2\), and verify that exactly **14 states** belong to \(S=0\), in agreement with the fourth Catalan number.
3. Build the \(256\times14\) isometric embedding \(V_S\) of these singlets and verify separately that **each of the three time-independent Hamiltonian terms preserves this subspace**, to norm less than \(10^{-10}\).
4. Embed the **four explicit logical valence-bond codewords**, verify their orthonormality and that they lie entirely inside the 14-dimensional SU(2) singlet sector.
5. Integrate the *actual time-ordered Hamiltonian* by a midpoint piecewise-constant product of exact \(14\times14\) matrix exponentials; **do not** integrate the second-order Schrieffer–Wolff approximation.
6. Independently perform a 256-dimensional sparse-matrix exponential propagation of the same nontrivial pulse at 24 midpoint slices. The projected full-versus-reduced unitary difference has Frobenius norm approximately **\(3.1\times10^{-13}\)**.

The last comparison validates the symmetry reduction against independent full-spin dynamics. Fine-mesh temporal convergence is checked separately at 256, 512 and 1,024 bins; this is necessary because coarse time slicing can produce accidental false leakage minima.

Let \(Q\) be the four-dimensional code isometry and \(U(T)\) the full time-ordered unitary. The exact code restriction is \(W=Q^\dagger U(T)Q\), which need not be unitary if leakage occurs.

**Leakage:** The worst-case probability of leaving this four-dimensional product-of-singlets code at the end of the pulse is

\[
L_{\max}=1-\lambda_{\min}(W^\dagger W).
\]

**Logical gate deviation:** For a diagonal logical gate, the phase combination \(\varphi_{00}-\varphi_{01}-\varphi_{10}+\varphi_{11}\) is invariant under global and independent logical-\(Z\) phase adjustments. A CZ-class target has this phase equal to \(\pi\) modulo \(2\pi\). The program constructs a *specified* CZ-class diagonal target \(D\) by distributing the computed entangling-phase mismatch equally over the four diagonal phases and measures \(\lVert W-D\rVert_2\). This is a reproducible **candidate phase-calibrated operator-norm deviation**, not a proven globally nearest target, a physical average gate infidelity, a stochastic Pauli error probability, or a demonstrated physical phase-calibration procedure.

The code separately reports the off-diagonal norm, which is at numerical roundoff for this particular symmetric Hamiltonian.

## 3. Results with independent convergence checks

The bounded search with 512 midpoint bins returns

\[
\boxed{J_{\max}/\Delta=1,\quad T\Delta\approx38.234506218.}
\]

| Pulse protocol | \(T\Delta\) | Phase-calibrated operator deviation | Worst-case code-space leakage probability |
|---|---:|---:|---:|
| Previous rectangular pulse, \(J/\Delta=0.4\), nominal perturbative time | 117.809725 | 0.124337 | 0.019131 |
| Rectangular pulse, same **peak** \(J/\Delta=1\), nominal perturbative time | 18.849556 | 0.297560 | 0.222458 |
| Smooth \(\sin^2\), \(J_{\max}/\Delta=1\), retimed as above (1,024-bin verification) | **38.234506** | **0.00000290** | **0.0000000977** |

The comparison to the older \(J/\Delta=0.4\) pulse includes a substantial increase in permitted peak coupling. The same-peak \(J/\Delta=1\) square pulse is included to make that distinction explicit; it is evaluated at its nominal perturbative pulse time, **not independently optimized for the best square-pulse duration**.

An unrounded worst-case leakage probability \(L=9.77\times10^{-8}\) corresponds to \(9.77\times10^{-6}\%\), not a 9.77%-level result.

### Numerical time-step convergence

| Midpoint time slices | Phase-calibrated operator deviation | Worst code leakage probability |
|---:|---:|---:|
| 256 | 0.00001167 | \(9.508\times10^{-8}\) |
| 512 (optimization grid) | \(4.86\times10^{-8}\) | \(9.715\times10^{-8}\) |
| **1,024** (independent verification) | **0.00000290** | **\(9.769\times10^{-8}\)** |

The objective is deliberately evaluated on a *different, finer discretization* after optimization. Its finite difference from the 512-bin optimization grid must be reported: a near-zero score **at the optimized grid** is not a continuum-time certificate. For the chosen pulse, the off-diagonal logical mixing remains below \(10^{-9}\), and the full-versus-reduced 256-state propagation agrees to better than \(10^{-9}\) for a separate 24-slice test. Time-discretization and hardware waveform errors beyond this check have not been bounded rigorously.

**Speed-versus-control comparison:** The shaped \(J_{\max}=1\) gate has duration about **38.23**, versus **117.81** for the older lower-peak \(J=0.4\) square pulse—a factor of \(\sim3.1\). That is *not* a same-coupling speedup. At the same \(J_{\max}=1\), the shaped pulse takes longer than the nominal square pulse but is vastly more accurate on the declared gate/leakage metrics.

## 4. Calibration robustness: still poor

A small nominal leakage score does **not** imply a robust gate. Retain the optimized waveform and duration, but alter the amplitude of **only the first physical cross exchange** \(J_1(t)\rightarrow(1+\epsilon)J_1(t)\). This fault is gauge-invariant at every instant.

| Coupling amplitude error \(\epsilon\) | Phase-calibrated operator deviation | Worst code leakage |
|---:|---:|---:|
| \(-10\%\) | 0.086337 | \(1.01\times10^{-7}\) |
| \(-5\%\) | 0.043368 | \(9.97\times10^{-8}\) |
| \(+5\%\) | 0.043734 | \(9.41\times10^{-8}\) |
| \(+10\%\) | 0.087801 | \(9.06\times10^{-8}\) |

The main problem is a **wrong accumulated logical entangling phase**, not leakage. For example, the \(+10\%\) perturbation changes the logical nonlocal phase by about \(0.351\) radians. Merely improving the pulse's adiabatic return to the code space cannot remove that error. No stochastic gate error rate or hardware noise model is inferred from these deterministic perturbations.

This is a concrete **release blocker** for a practical error-correcting machine: the nominal calibrated pulse is excellent in this model, but no robustness to physically allowed exchange-strength fluctuations has been established. A follow-on composite-pulse/echo design must explicitly suppress the derivative of the nonlocal phase with respect to such perturbations, *not just optimize leakage*.

## 5. What changed relative to prior work

- The source physics is the same specified two-body Hamiltonian from PR #15; no new matter representation is quietly introduced.
- The pulse is a real time-dependent pair-exchange schedule, simulated beyond perturbation theory.
- The exact global SU(2) singlet sector enables efficient but *verified* 14-dimensional evolution; the full 256-dimensional comparison checks the reduction.
- The new gate should be understood as a calibrated **CZ-class entangler**, not a literal native controlled-\(Z\) acting perfectly on the full 256-dimensional physical space.
- The amplitude-coupling fault test continues to show substantial coherent logical errors that cannot be caught by Gauss-law measurement alone.

This construction uses established ideas of exchange-only quantum control and smooth pulse shaping. Relevant prior work includes [Heinz et al., *Fast Quantum Gates for Exchange-Only Qubits Using Simultaneous Exchange Pulses* (PRX Quantum 2025)](https://doi.org/10.1103/njq3-fcdd), and [*Singlet-Only Always-On Gapless Exchange Qubits with Baseband Control* (PRL 2025)](https://doi.org/10.1103/j1hn-k5m8). These papers concern related physical architectures; **neither validates the particular SU(2) gauge-rishon mediator or numerical waveform constructed here**.

## 6. Reproduction and next decision

\`\`\`bash
python -m pip install numpy sympy==1.14.0 scipy==1.17.0
OPENBLAS_NUM_THREADS=1 python su2_shaped_pulse.py
\`\`\`

The GitHub Actions research workflow installs the numerical dependency and runs the independent 256-dimensional crosscheck.

**Next high-value step:** Develop an exchange-based **composite phase gate** that cancels first-order sensitivity to calibration in either physical cross coupling without increasing leakage or duration excessively. It may require noncommuting logical rotations and thus additional controllable SU(2)-scalar exchanges. Compare full-spin pulse errors under \(\epsilon=\pm1\%,\pm5\%,\pm10\%\), not only best-case nominal fidelity. The second milestone is to expose this *physical pulsed interaction* to the five-module flagged extraction simulator with explicit leakage states and to test a fault during each elementary interaction.

Do not promote this result to a full Standard Model quantum-error-correction mechanism: no physical Standard Model encoding, executable gauge-compatible multi-module stabilizer measurement, physical noise model, or device-level fault tolerance has been demonstrated.
