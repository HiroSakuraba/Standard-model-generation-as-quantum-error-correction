# SU(3) five-plaquette memory: explicit recovery versus CPTP Wilson faults

**18 October 2026 — reproducible finite-lattice simulation, not experimental noise calibration or a glueball mass.**

Executable: [`su3_channel_recovery.py`](su3_channel_recovery.py); structured output: [`su3_channel_recovery_results.json`](su3_channel_recovery_results.json). Builds on [`su3_four_five_plaquette.py`](su3_four_five_plaquette.py) and the Chen–Müller–Yao square-plaquette model (Phys. Rev. D **114**, 074502, DOI [10.1103/1qt7-3ywl](https://doi.org/10.1103/1qt7-3ywl)). The original finite-error correction criterion is Knill–Laflamme (Phys. Rev. A **55**, 900, DOI [10.1103/PhysRevA.55.900](https://doi.org/10.1103/PhysRevA.55.900)).

## 1. Purpose and claim boundaries

We previously verified all 441 Knill–Laflamme scalar overlap equations for the two five-plaquette qutrit states

\[
 |0_L\rangle=|0,0,0,0,1\rangle,\qquad
 |1_L\rangle=|0,0,0,0,2\rangle
\]

against the 21-element **linear single-event** error space

\[
 \mathcal E=\operatorname{span}\{I,E_0^2,\ldots,E_{15}^2,B_0,B_1,B_2,B_3\},
 \quad B_i=W_i+W_i^\dagger.
\]

The rank of the exact Knill–Laflamme Gram matrix is **5**, not 21. The present work constructs the associated **explicit completely positive trace-preserving (CPTP) decoder**, rather than stopping at recovery existence; then it stress-tests that fixed decoder under two fully specified channels that go *beyond* this error span. All state spaces remain within the **243-dimensional gauge-invariant minimal-link qutrit truncation**. No implementation on a physical lattice or gauge-compatible experimental recovery circuit is constructed.

## 2. Concrete recovery and its exact target

Let \(K:\mathbb C^2\to\mathbb C^{243}\) embed the logical basis. Let \(A_a=E_aK\) be the 243×2 error images, and \(G_{ab} I_2=K^\dagger E_a^\dagger E_b K\). Exact SymPy checks from the previous implementation certify that \(G\) is Hermitian positive semidefinite of rank five. Diagonalizing its **numerical representation** produces five orthogonal syndrome isometries \(V_r:\mathbb C^2\to\mathbb C^{243}\), satisfying \(V_r^\dagger V_s=\delta_{rs}I_2\). Set \(S=\sum_rV_rV_r^\dagger\). The explicit recovery is

\[
 \boxed{\mathcal R(X)=\sum_{r=1}^{5}V_r^\dagger XV_r+
  \operatorname{tr}[(I-S)X]\,\frac{I_2}{2}.}
\]

This maps the 243-state system to a logical qubit. Both terms are CP, and the trace is preserved since the syndrome projectors plus their orthogonal complement sum to identity. The program verifies the ten syndrome vectors are orthonormal and checks **all four logical matrix units** after a demonstration channel supported on the error span. The reported entanglement fidelity is \(F_e=1\) up to double-precision rounding.

**Important conditioning:** To define that exact recovery test as a globally trace-preserving channel, the demonstration Kraus operators on the code are proportional to \(E_aK\), with the extra Kraus operator \(I-KK^\dagger\) completing the map outside the code. It explicitly **uses the nonlocal code projector** and is *not* evidence that a physical local noise source has precisely that action. The exact mathematical result is that the recovery reverses any channel supported by the declared error span *on initial code inputs*.

## 3. Physically motivated but uncalibrated CPTP directed-Wilson jumps

The projected fundamental plaquette operator \(W_i\) is not unitary. To construct a genuine channel for an **unconditioned local Wilson jump**, use two Kraus operators:

\[
 K_1=\sqrt{p}\,W_i,\qquad
 K_0=\sqrt{I-pW_i^\dagger W_i},\qquad 0\le p\le1.
\]

The squared transition amplitudes in this published minimal-link Hamiltonian ensure \(0\le W_i^\dagger W_i\le I\). The executable checks \(K_0^\dagger K_0+K_1^\dagger K_1=I_{243}\) directly on the full model space, including states outside the code. The model parameter \(p\) is **not** a laboratory fault probability or a measured Lindblad rate. The no-jump branch depends on \(W_i^\dagger W_i\), as it must for trace preservation. It cannot be discarded when computing fidelity.

Apply **two independent instances** of the same noise channel separated by *three equal evolution intervals* of the source Hamiltonian at \(g^2=1\), using the dimensionless storage duration \(t\):

\[
 \mathcal N_{t,p} =\operatorname{Ad}_{U(t)}^\dagger\circ
 \operatorname{Ad}_{U(t/3)}\circ\mathcal J_p\circ
 \operatorname{Ad}_{U(t/3)}\circ\mathcal J_p\circ
 \operatorname{Ad}_{U(t/3)},
\]

where the leading \(\operatorname{Ad}_{U(t)}^\dagger\) is **perfect inversion of the known entire Hamiltonian evolution**. This highly favorable compensation is included only to isolate fault sensitivity; an actual device would need to realize it and incur its own errors. The script also reports the un-compensated evolution and does **not** call ideal inversion a free practical control.

All branches of the two-event channel are retained, including no-jump/jump orderings and two-jump trajectories. The same fixed recovery \(\mathcal R\) is applied in every case. We record both **raw** logical entanglement fidelity after the ideal inversion but before recovery and **recovered** entanglement fidelity. The failure-space weight is included in recovery via maximally mixed replacement, never silently postselected away.

### Results: two Wilson jump opportunities, p = 0.10 per event

| Dimensionless t | Error at face | Raw \(F_e\) after ideal inverse | Recovered \(F_e\) |
|---:|---|---:|---:|
| 0 | \(W_3\) | 0.956049 | **0.982037** |
| 0 | \(W_3^\dagger\) | 0.956049 | **0.982037** |
| 0 | \(W_4\) | 0.810000 | **0.833750** |
| 0.15 | \(W_3\) | 0.955501 | **0.981602** |
| 0.15 | \(W_4\) | 0.810412 | **0.834104** |
| 0.6 | \(W_3\) | 0.951999 | **0.977861** |
| 0.6 | \(W_4\) | 0.823288 | **0.842620** |

At \(p=0\), every ideally inverse-compensated channel is the identity and recovers with \(F_e=1\). At \(p=.02\), zero-duration two-event recovered fidelity is **0.996393** for \(W_3\) and **0.965350** for \(W_4\). The violation at \(W_3\) is expected: even a **single** directed loop error has the exact code-space syndrome obstruction \(K^\dagger W_3^\dagger W_3 K=\operatorname{diag}(1/9,1/3)\). By contrast, a Hermitian combination \(B_3\) belongs to the prior *single-event correctable span*. **A Hermitian sum passing Knill–Laflamme is not equivalent to its directed components passing separately.**

## 4. Finite coherent Hermitian kicks are not single-event errors

For a second valid, gauge-invariant CPTP stress model, use equal-probability random unitaries

\[
 K_\pm=2^{-1/2}\exp(\pm i\theta B_3),
\]

at both fault opportunities, with evolution and ideal inverse as above. At a single event's *first order* in \(\theta\), its nontrivial term is proportional to the exactly correctable \(B_3\). At finite angle, powers \(B_3^2,B_3^3,\ldots\) are introduced; with intervening Hamiltonian dynamics, time-dressed errors add further operators.

| t | Kick angle \(\theta\) | Raw \(F_e\) | Recovered \(F_e\) |
|---:|---:|---:|---:|
| 0 | 0.1 | 0.991173 | **0.999994** |
| 0 | 0.3 | 0.924931 | **0.999571** |
| 0.6 | 0.1 | 0.989092 | **0.998051** |
| 0.6 | 0.3 | 0.908284 | **0.983784** |

This quantitatively illustrates near-exact correction for sufficiently weak coherent kicks and its finite-amplitude/time limitations. The tested duration and angle are arbitrary dimensionless scales; **there is no inferred experimental memory lifetime or quantum fault-tolerance threshold**.

## 5. Native Hamiltonian itself is not protected storage

Without the favorable full-Hamiltonian inverse, the *zero-fault* native evolution at \(t=.6\) has **raw** embedded-code entanglement fidelity \(\approx0.254693\). Applying the same recovery improves it to only \(\approx0.632749\), even with **no** jump faults. This is the fundamental practical obstacle: the native magnetic term at the **encoded last plaquette** is allowed and moves information outside the error subspace our five-syndrome decoder corrects. A known unitary is theoretically invertible, but the missing issue is implementing reliable, economical **gauge-invariant** controls and recovery under real noise.

## 6. Next engineering decision

The constrained result is a useful **mathematical code and explicit decoder for a narrow, spatially selected channel**; it is not a passive gauge-protected memory. The more realistic stress tests fail at finite strength or when the native Hamiltonian is left active. The next implementation should compare syndrome readout and noisy control costs; consider alternative code locations with bounds on realistic directed-loop jumps and finite-duration recovery, and an implementable protection mechanism against the encoded plaquette's native \(B_4\). If no credible native way to suppress \(B_4\) is available, the present two-orientation encoding should **not** be prioritized for hardware implementation.

## Reproduction

```bash
python -m pip install sympy==1.14.0 numpy scipy
OPENBLAS_NUM_THREADS=1 python su3_channel_recovery.py
python -m py_compile su3_channel_recovery.py
```

The tests use exact SymPy for the 441 KL relations and Gram rank; orthogonal recovery is formed in double precision; the 243-state Hamiltonian is diagonalized in double precision. Asserted CPTP completeness and logical matrix-unit recovery are numerical checks, not formally certified interval inequalities. Inputs are the published minimal-link qutrit matrices already independently benchmarked in this repo. In contrast to an experimentally calibrated device, the "fault probabilities" and angles here are **chosen test parameters**.