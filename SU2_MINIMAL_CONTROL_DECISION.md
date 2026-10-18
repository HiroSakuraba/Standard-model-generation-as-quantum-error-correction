# Minimal-control decision for the SU(2) singlet entangling gate

**18 October 2026 — verified numerical engineering comparison, not a physics discovery about the full Standard Model.** This report compares a **single smooth exchange pulse without additional physical logical \(X/Y\) rotations** against the longer six-envelope SK1-inspired composite sequence. It includes calibration faults in *all control families required by each design*, rather than only the two cross-register exchanges. Reproduce with [\`su2_minimal_pulse_audit.py\`](su2_minimal_pulse_audit.py).

## Executive decision

**For the declared worst-case deterministic calibration tests, prefer the original single smooth pulse as the less complex reference design.** A six-pulse composite indeed suppresses coherent errors in its **two main cross-register exchange couplers**, but is **more** sensitive to even small systematic errors in the many additional **intra-register exchange rotations** required to implement its correction axes. The single pulse has no such additional local rotations to miscalibrate. The composite is also approximately nine times longer in cross-exchange evolution alone.

**This is a conditional engineering decision, not a claim that the single pulse is intrinsically error-resistant or is a physical CZ device.** Its own cross-coupler drift remains a serious logical error. A device with extremely accurate intra-register exchange control but imprecise cross couplers could still favor the composite. The tests count a finite set of coherent fault scenarios; they do **not** infer statistical likelihoods, error thresholds or a hardware noise budget.

## Physical model and direct comparison

The physical system is unchanged from the two preceding studies: **eight spin-\(\tfrac12\) SU(2) fundamentals co-located at one vertex**, with an engineered block-spin Casimir penalty. Each set of four spins has a two-dimensional local color-singlet subspace, and their product defines the **four-dimensional logical code**. Gauge-invariant *pairwise* Heisenberg interactions mediate the two-qubit gate:

\[
H(t)=\Delta(\mathbf S_A^2+\mathbf S_D^2)
+J_1(t)\,\mathbf S_{A_0}\cdot\mathbf S_{D_0}
+J_2(t)\,\mathbf S_{A_1}\cdot\mathbf S_{D_1}.
\]

The microscopic eight-spin Hilbert space has dimension 256, and its exact total-spin-zero invariant subspace has dimension 14. Previous scripts verified this symmetry reduction against direct evolution in all 256 spin states.

**Minimal design:** A single, simultaneous waveform \(J_1(t)=J_2(t)=\Delta\sin^2(\pi t/T)\), \(T\Delta=38.234506218\), calibrated to a controlled-Z-equivalent nonlocal phase. Unlike SK1, it requires **zero additional exchange-generated logical \(X/Y\) rotations during the entangling pulse**. A fixed global or local logical \(Z\)-frame adjustment is allowed when declaring its nominal equivalence to CZ. Whether such phase tracking can be implemented reliably during *physical five-block syndrome extraction* has not been established, so this is not a free verified hardware operation.

**Composite design:** Three echoed entangling blocks (six smooth cross-exchange envelopes) plus many **extra local SU(2)-scalar exchange rotations**, at total cross-exchange duration \(T\Delta=344.208400\), excluding all local-rotation durations. Each pulse is propagated through the **full time-dependent 14D Hamiltonian**; none is replaced by a two-qubit perturbative effective-gate approximation.

Both nominal gates have tiny model leakage: approximately \(9.72\times10^{-8}\) for the minimal pulse versus \(1.73\times10^{-6}\) for the composite. The shorter pulse has lower *computed* leakage, but that alone does not imply a device can execute it fault tolerantly.

### Error metric and calibration rules

For any pulse \(U\), let \(Q\) embed the four physical codewords into the 14D singlet sector. Its compressed action is \(W=Q^\dagger UQ\), which can be nonunitary due to code-space leakage.

For each **fixed physical pulse design**, compare a faulted \(W_\epsilon\) with its **own fixed nominal** \(W_0\), using only one physically irrelevant global-phase alignment:

\[
\delta_\epsilon
=\left\|e^{-i\arg\operatorname{Tr}(W_0^\dagger W_\epsilon)}
 W_\epsilon-W_0\right\|_2.
\]

This is a reproducible **spectral-norm deviation**; the chosen trace-based phase alignment is not certified to give the mathematically closest phase. **No fault-dependent recalibration of local logical phases is permitted.** The initial nominal CZ-class comparison is stated separately from this perturbation metric.

Fault scenarios are coherent **multiplicative amplitude errors** \(\pm1\%\), \(\pm5\%\), \(\pm10\%\), applied independently to (a) the first cross coupler, (b) the second, or (c) both together with the same sign. For the composite only, also apply a common coherent angle miscalibration to **all its required intra-register exchange rotations**. This is a deliberately specified finite scenario set, not exhaustive over arbitrary errors, waveforms, drifts or correlations.

## Reproducible robustness results

For each uncertainty level, the table reports the **largest deviation found among the sampled faults** in that column. It is *not* an estimated failure probability or physical average fidelity.

| Allowed amplitude uncertainty | Single smooth: worst cross-coupler fault | Composite: worst cross-coupler fault | Composite: worst extra-local-rotation fault |
|---:|---:|---:|---:|
| \(\pm1\%\) | **0.034573** | 0.001218 | **0.101316** |
| \(\pm5\%\) | **0.177005** | 0.031567 | **0.517882** |
| \(\pm10\%\) | **0.364117** | 0.127834 | **1.034373** |

Under *cross-coupler errors alone*, the composite is the clear winner at all three levels. Once its **extra local control operations** are admitted as possible calibration faults with comparable fractional uncertainty, its **worst sampled logical deviation exceeds the simple pulse's** at every tested level. Its approximately 9× greater cross-exchange duration also amplifies exposure to any time-dependent decoherence (not simulated here).

For clarity, the short gate's \(\pm1\%\) error in **one** cross coupler is already approximately **0.0172**, while a same-sign error in **both** is approximately **0.0346**. The short gate therefore does **not** somehow become intrinsically calibration-immune; it wins this comparison by *avoiding an additional much more fragile family of controls*.

The composite-local-rotation test is one common coherent systematic error applied to all local rotation angles, rather than a single isolated faulty local exchange. Faults on individual local rotations, finite pulse widths, always-on gap Hamiltonian during rotations, ancilla preparation/readout and stochastic idle errors are not included in this comparison. The effect of those omissions could change the design decision.

##  A precise weak-coupling obstruction to a single-envelope calibration sweet spot

The earlier exact second-order SU(2) mediator calculation gives

\[
J_{ZZ}^{(2)}=-\frac{J_1J_2}{24\Delta}.
\]

For a fixed exchange-waveform shape and duration, write the corresponding **leading-order nonlocal phase** as \(\theta=C\,a_1a_2\), where \(a_i\) multiply the physical coupler envelopes and \(C\) is a fixed nonzero integral. Independent multiplicative calibration errors give

\[
\theta(\epsilon_1,\epsilon_2)
=\theta(0,0)(1+\epsilon_1)(1+\epsilon_2).
\]

Consequently

\[
\left.\frac{\partial\theta}{\partial\epsilon_1}\right|_0
=\left.\frac{\partial\theta}{\partial\epsilon_2}\right|_0
=\theta(0,0)\ne0
\]

for any nontrivial target entangling angle at this order. The new executable verifies this symbolic identity with exact SymPy algebra.

**What this theorem says:** Simply resculpting **one fixed weak-coupling envelope** cannot remove the leading response to each amplitude scaling while retaining a nonzero entangling phase *within this second-order effective model*.

**What it does not say:** It does not rule out a fully nonperturbative sweet spot, sign-changing couplings, noncommuting multi-pulse sequences, alternative interactions, measurement-based recovery, or all single-pulse control designs. In particular the preceding \(J/\Delta=1\) gate is outside a rigorously established small-\(J/\Delta\) asymptotic error bound; the practical conclusion for that pulse is based on its actual numerical simulation, not extrapolation from this derivative identity.

## Small bounded nonperturbative waveform search

To test whether a *modest* change of a positive one-pulse waveform solves the calibration problem at fixed peak bound, define

\[
J_{1,2}(t)=
\Delta\,\frac{\sin^2(\pi t/T)\,[1+a\cos(2\pi t/T)]}
                  {1+|a|},
\quad a\in\{-0.2,0,+0.2\}.
\]

This shape is nonnegative, uses **no extra local exchange rotations**, and keeps the physical coupler amplitudes bounded by \(\Delta\). At each fixed \(a\) optimize only \(T\Delta\in[30,75]\) for a *nominal CZ-class* phase-calibrated gate. Evaluate nominal accuracy and worst sampled \(\pm5\%\) coupler deviations independently at **512 midpoint steps** (after a 256-step duration calibration).

| Waveform \(a\) | Optimized \(T\Delta\) | Nominal phase-calibrated operator deviation | Nominal code leakage | Worst sampled \(\pm5\%\) deviation | Valid nominal CZ under threshold 0.001? |
|---:|---:|---:|---:|---:|---|
| \(-0.2\) | 42.782325 | \(1.36\times10^{-5}\) | \(1.70\times10^{-8}\) | **0.176940** | Yes |
| \(0\) | 38.234506 | \(4.86\times10^{-8}\) | \(9.72\times10^{-8}\) | **0.177007** | Yes |
| \(+0.2\) | 74.999995 *(search boundary)* | **0.034566** | \(5.18\times10^{-9}\) | 0.167586 | **No** |

The \(a=-0.2\) waveform suppresses nominal leakage a little further but takes longer and does **not** materially improve calibration tolerance. The \(a=+0.2\) waveform ends at the declared **search boundary** without meeting the required nominal CZ accuracy; its marginally smaller fault response **does not qualify as an improved gate**. This is a small sampled search, not an exhaustive optimization over physical waveforms or all pulse durations.

## The decision and next engineering check

The evidence supports **simplifying rather than extending SK1** as the near-term research baseline:

1. Keep the **single, calibrated \(\sin^2\) exchange pulse** as the reference gate for subsequent system-level studies: one analog envelope family, no mandatory corrective logical X/Y exchange pulses, shorter duration, and lower modeled nominal code leakage.
2. Use the composite sequence only as a **conditional comparison** for experiments in which local exchange rotation accuracy is substantially better than cross-coupler accuracy. The present project has no physical evidence such an asymmetric noise advantage exists.
3. Stop claiming that nominal leakage suppression or a complicated composite automatically provides fault-tolerant quantum memory. The **single gate is still very sensitive to coherent amplitude drift**, and none of the present studies supplies physical fault-tolerant stabilizer measurements.
4. Before replacing ideal controlled-Paulis in the five-block decoder, model **the elementary exchange pulses and their gauge-preserving noise**, possible leakage, and the ancilla's preparation and readout.

**Next useful experiment:** retain the simple waveform but introduce a concrete, physically permitted *independent calibration monitor* or a jointly optimized bias/phase control that can reduce effective amplitude uncertainty **without adding fragile noncommuting local rotations**. Then compare the *cost of that calibration mechanism* against the extra composite-control burden. Without new hardware capabilities or extra observables, the weak-coupling derivative obstruction limits what simple fixed envelopes can guarantee.

## Reproduction

\`\`\`bash
python -m pip install numpy==2.3.3 sympy==1.14.0 scipy==1.17.0
OPENBLAS_NUM_THREADS=1 python su2_minimal_pulse_audit.py
\`\`\`

This is an ordinary GitHub Actions test in the existing research workflow. It imports the previously merged full-spin SU(2) mediator and both pulse-control verifiers; earlier test results are retained unchanged.

**Related primary literature:** Brown, Harrow and Chuang, [*Arbitrarily accurate composite pulse sequences*](https://doi.org/10.1103/PhysRevA.70.052318), *Physical Review A* **70**, 052318 (2004), provides established noncommuting-axis correction methods. The present numerical gate calculations are about a separately engineered matter/rishon SU(2) truncation. They do **not** establish natural quantum error correction in the full Standard Model.
