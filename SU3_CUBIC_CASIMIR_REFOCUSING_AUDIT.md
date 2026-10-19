# Gauge-invariant SU(3) cubic-electric refocusing of the five-plaquette memory

**October 18, 2026 — finite-model control experiment and exact operator identity.** Companion: `su3_cubic_casimir_echo.py`. The Kogut–Susskind five-plaquette minimal-link Hamiltonian and its magnetic coefficients come from Chen, Müller and Yao (*Physical Review D* **114**, 074502, 2026, DOI 10.1103/1qt7-3ywl). **The cubic-electric pulse capability is our hypothetical additional control, not an established capability of their hardware proposal or of real QCD.** This is not a glueball mass, fault-tolerance threshold, or new Standard Model quantum code.

## Why this approach

Earlier experiments proved that a specified 21-operator, bounded magnetic/electric error span can be decoded on the end-flux logical pair, but the same code is not stored by the native five-face Hamiltonian: with dimensionless time `t=.6`, embedded-code entanglement fidelity drops to `0.25469270`, and the fixed five-syndrome decoder raises it only to `0.63274935`. Our new question is whether a **local gauge-invariant electric pulse family** can suppress that coherent native evolution without using a perfect inverse of the complete 243×243 Hamiltonian.

## Exact SU(3) cubic-Casimir control

Use SU(3) Dynkin labels `(p,q)` and generators normalized so the quadratic electric Casimir of the fundamental is `4/3`. The cubic Casimir has eigenvalues

\[
 C_3(p,q)=\frac{(p-q)(2p+q+3)(p+2q+3)}{18},\quad
 C_3(\mathbf{1})=0,\quad C_3(\mathbf 3)=+10/9,\quad C_3(\bar{\mathbf3})=-10/9.
\]

On each **top link** of the published square-chain gauge-fixed basis, representation `q_i∈{0,1,2}` is respectively `1,3,bar3` (view `2` as `−1`). Define the gauge-invariant link-phase operator

\[
 Q_i=\exp\left[\frac{2\pi i}{3}\frac9{10}C_3(E_{\mathrm{top},i})\right],\qquad Q=\prod_{i=0}^{4} Q_i.
\]

The cubic Casimir of a link is invariant under the independent SU(3) transformations acting on either endpoint: it is constructed from local electric generators and acts as a representation scalar. Thus `Q` preserves every local Gauss-law sector *kinematically*. This statement does **not** establish that an experimental platform can efficiently pulse the cubic operator or can address five links independently.

For the directed elementary Wilson-loop operator `W_i`, which changes the face label `q_i→q_i+1 mod 3`, we have exactly

\[
 Q W_i Q^\dagger=\omega W_i,\qquad
 Q W_i^\dagger Q^\dagger=\omega^{-1}W_i^\dagger,
 \quad \omega=e^{2\pi i/3}.
\]

The executable verifies this charge for **all 1,215 nonzero directed Wilson transition entries** of the five-plaquette model. With the source-normalized native Hamiltonian at `g²=1`, `H=D−Σ_i(W_i+W_i†)` up to the physically irrelevant additive `30I`, the exact algebraic clock twirl is

\[
 \frac13\sum_{k=0}^{2}Q^k H Q^{-k}=D.
\]

Therefore, in the **infinitely rapid control limit**, the engineered effective Hamiltonian contains only electric energy. The pair `|00001>`,`|00002>` have identical electric energy, so they acquire the same dynamical phase and retain their logical superposition. This is dynamical decoupling (Hamiltonian engineering), **not autonomous gauge protection or quantum error correction**. It also removes the desired native magnetic glue-field dynamics, so the same device is no longer performing faithful Yang–Mills evolution during memory storage.

A two-phase `(-1)^occupation` control from quadratic electric information does not distinguish `3` from `bar3`. The within-encoded-pair Hermitian `B_4` bit-flip survives that control; the numerical binary negative control confirms the deficit.

## Finite-cycle test, no physical fault channel

Symmetric toggling cycle uses frames `k=(0,1,2,1,0)` and normalized free-evolution durations `(1/6,1/6,1/3,1/6,1/6)`. Between segments are four `Q`/`Q†` pulses which return the lab-frame toggle to identity at the end of each cycle. Matrix exponentials are propagated numerically on the **complete 243-dimensional** finite model with SciPy sparse `expm_multiply`. Total dimensionless free-evolution time is `t=.6`. The fixed recovery is the previously committed explicit five-syndrome completely positive trace-preserving (CPTP) map, with its positive complement branch. `F_raw` is entanglement fidelity relative to the original two-state code **without** a recovery step; `F_fixed_recovery` includes that specific fixed decoder, not an optimized one.

| Symmetric cycles | Global Q pulses | F_raw | F_fixed_recovery | Original subspace leakage |
|---:|---:|---:|---:|---:|
| No control | 0 | 0.254693 | 0.632749 | Not separately reported |
| 1 | 4 | 0.901759 | 0.980559 | 0.098022 |
| 2 | 8 | 0.995500 | 0.999216 | 0.004493 |
| 4 | 16 | 0.999741 | 0.999956 | 0.000259 |
| 8 | 32 | 0.999984 | 0.999997 | 0.00001583 |
| 16 | 64 | 0.999999014 | 0.999999833 | 0.000000984 |
| 32 | 128 | 0.999999938 | 0.999999990 | 0.0000000614 |

The exact group average is not equal to a finite pulse sequence because toggled Hamiltonians do not all commute. The finite-cycle results above include their commutators. No extrapolation to a physical dephasing time is justified.

## Controls are not free: calibration and finite durations

Assume every global ideal clock pulse has a coherent fractional angle calibration error `η` and keep 8 cycles / 32 ideal instantaneous pulses:

| η | F_raw | F_fixed_recovery |
|---:|---:|---:|
| 0 | 0.999984 | 0.999997 |
| 0.01 | 0.999637 | 0.999680 |
| 0.05 | 0.990481 | 0.992129 |
| 0.10 | 0.961211 | 0.968781 |

For the next test, model every clock pulse as **a finite interval in which the original native Hamiltonian remains switched ON** and a signed control `±ΩΣ_i(9/10)C3(top_i)` is applied. Set `Ω=2π/(3τpulse)` per top link, with 8 cycles, **32 pulse intervals** and free-evolution time `t=.6`:

| Pulse width τ | Ω per top link (relative model units) | Total elapsed time | F_raw | F_fixed_recovery |
|---:|---:|---:|---:|---:|
| .0005 | 4,188.79 | .616 | .999825 | .999925 |
| .002 | 1,047.20 | .664 | .998036 | .998879 |
| .01 | 209.44 | .920 | .984995 | .989380 |

These controls are **very strong compared with the original O(1) magnetic interaction**. No experimental strength, bandwidth, crosstalk, finite rise-time or anharmonicity constraint has been established. The simulation verifies the lab-frame finite-pulse construction approaches the ideal toggling-frame evolution as pulse duration tends to zero (up to a global phase); it does not demonstrate hardware feasibility. Pulses during idle intervals and decoder steps were not modeled.

## Decoherence not repaired by Hamiltonian refocusing

An unconditioned gauge-invariant directed-Wilson jump channel at a chosen plaquette is specified by

\[
 K_1=\sqrt{p}\,W_j,\qquad K_0=\sqrt{I-pW_j^\dagger W_j},\quad p=.1.
\]

The program checks trace preservation over **all 243 states**, inserts one jump opportunity halfway between two 8-cycle echo blocks, then measures logical entanglement fidelity after the fixed decoder. This is **one** event opportunity, rather than the *two* used in the earlier channel experiment; do not compare their numbers as if identical conditions.

| Jump location | F_raw | F_fixed_recovery |
|---|---:|---:|
| Neighboring plaquette 3 | 0.977745 | 0.990612 |
| Encoded plaquette 4 | 0.899999 | 0.912500 |

Thus our coherent refocusing succeeds against the native Hamiltonian yet remains exposed to quantum jumps at the encoded plaquette. A fixed decoder and fast control are not equivalent to a fault-tolerant memory. The simpler **8-cycle binary quadratic-electric** negative control gives `F_raw=0.680898`, `F_fixed_recovery=0.680908` even with no injected jumps.

## Precise scope and next engineering decision

**Established:** exact clock-charge identity and cancellation of the complete *declared* minimal-link magnetic Hamiltonian under an ideal group twirl; finite-cycle unitary simulation; explicit phase calibration and finite-duration control stress; full-state CPTP test for one directed Wilson jump; benchmark with the earlier fixed recovery.

**Not established:** cubic-Casimir pulses as native gates, sufficiently high available control bandwidth (required strength reaches thousands of uncalibrated simulation units in the fast-pulse test), suppression of real high-frequency noise, recovery measurements or their locality, fault tolerance, error thresholds, robustness to higher-link irreps, hardware device fidelity, continuum glueball masses or any new Standard Model generation QEC code.

**Next priority:** build a device-native cubic-electric gate synthesis or conclude the cubic control is prohibitive; include pulse errors at each individual link and the real time/energy cost. Separate a *memory mode* which suppresses gauge dynamics from a *physics simulation mode* which retains it. Do not treat both as the same operating Hamiltonian.

## Reproduction

```bash
python -m pip install 'sympy==1.14.0' numpy scipy
OPENBLAS_NUM_THREADS=1 python su3_cubic_casimir_echo.py
python -m py_compile su3_cubic_casimir_echo.py
```

The script requires existing repo files `su3_four_five_plaquette.py`, `su3_channel_recovery.py`, `su3_three_plaquette_chain.py`, `su3_published_square_chain_benchmark.py`, `su3_complete_minimal_links.py`, `su3_two_plaquette_spin_network.py`. All reference data are from this project except the referenced published Kogut–Susskind Hamiltonian.

**Primary references:** Chen, Müller, Yao, *Minimally truncated SU(3) lattice gauge theory and string tension*, *Phys. Rev. D* **114**, 074502 (2026), DOI https://doi.org/10.1103/1qt7-3ywl; Knill and Laflamme, *Theory of quantum error-correcting codes*, *Phys. Rev. A* **55**, 900 (1997), DOI https://doi.org/10.1103/PhysRevA.55.900.