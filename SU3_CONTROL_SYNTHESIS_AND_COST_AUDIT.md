# SU(3) square-chain clock refocusing: conditional native compilation and control costs

**Status: 18 October 2026 — reproducible finite-lattice simulation and finite-model operator identities.**

Companions: `su3_native_control_compilation.py`, `su3_control_mask_scan.py`, `su3_native_control_compilation_results.json`, `su3_control_mask_scan_results.json` and logs. This extends `SU3_CUBIC_CASIMIR_REFOCUSING_AUDIT.md` without altering its data. This is **not** a compiled physical SU(3) link gate, a realized hardware device, a quantum error-correction threshold, or a glueball mass calculation.

## 1. Model and what changed

The model is the **published Chen–Müller–Yao** pure SU(3) square plaquette chain, five open faces, zero external flux and the minimal per-link `1⊕3⊕bar3` cutoff. The gauge-invariant Hilbert space is `3^5=243`. At `g²=1`, after omitting a physically immaterial scalar, `H=D−Σ_i(W_i+W_i†)`; `D` is the sum of 16 link quadratic-electric Casimirs and the directed fundamental plaquette `W_i` shifts face-qutrit flux `q_i→q_i+1 mod 3` with the published source-table amplitude. The logical states are `|00001>` and `|00002>`. The fixed **nonlocal**, five-syndrome trace-preserving recovery is inherited unchanged from `su3_channel_recovery.py`.

We ask whether we must act on all five representation-sensitive top-link sectors, and whether a faster alternative to pulsing a strong cubic SU(3) electric control exists on suitable quantum-simulator architectures.

## 2. Exact algebra: necessary representation resolution and two conditional implementations

The quadratic Casimir has `C2(1)=0`, `C2(3)=C2(bar3)=4/3`. The cubic Casimir has `C3(1)=0`, `C3(3)=+10/9`, `C3(bar3)=−10/9`. Consequently, **on the truncated representation space**

`Q_i = exp[(2π i/3)(9/10) C3(top-link i)] = diag(1,ω,ω²)` with `ω=exp(2π i/3)`.

Let `n_i` denote the face qutrit's representation orientation, `0,+1,−1`. The equivalent **encoded** diagonal operator is `Z12_i=|1><1|−|2><2| = diag(0,1,−1)`. On the gauge-fixed square-chain basis, **the operator identity** `Z12_i=(9/10)C3(top-link i)` holds. This is an equality *within this cutoff and qutrit encoding*, not a claim that a single microscopic fundamental SU(3) color qubit implements the cubic operator.

For each directed Wilson-loop term `W_j`, `Q_i W_j Q_i† = ω W_j` if `i=j`, and otherwise `W_j` is unchanged. The modulo-3 wrap `2→0` gives the same phase because `ω^(−2)=ω`. Hence the three-phase twirl of `H` for any selected control subset `S` is **exactly**

`(H + Q_S H Q_S† + Q_S² H Q_S†²)/3 = D − Σ_(j∉S)(W_j+W_j†)`.

Both implementations below give these same toggled Hamiltonians:

1. **Physical representation-sensitive electric pulses:** Generate ±`(2π/3) Σ_(i∈S) Z12_i` using a real, gauge-invariant cubic-electric control. This is not part of the **native published Hamiltonian**. A symmetric refocusing cycle uses four physical phase pulses, and finite-width pulses keep native `H` on simultaneously.
2. **Programmable magnetic-drive phases:** *If and only if the simulator can physically program both Wilson quadratures*, use `H_k=D−Σ_{j∈S}(ω^k W_j+ω^(−k) W_j†)−Σ_{j∉S}(W_j+W_j†)`. Thus phase updates synthesize the toggling frames **without high-amplitude electric pulses**. The individual `i(W_j−W_j†)` Hermitian quadrature is gauge invariant but odd under charge conjugation. A logical/virtual qutrit `Z12` frame update is not sufficient when magnetic interactions are fixed and always on: the **physical** drive phase must change. No cited SU(3) platform is shown to provide all required magnetic quadratures and selective timing.

**Lower-order control limit:** Both codewords have *identical values for all 16 local quadratic Casimirs*, yet `P B_4 P=X_L`. Any diagonal phase pulse synthesized **solely from functions of these commuting quadratic Casimirs** acts as the same phase on both codewords and cannot alter the `B_4` logical-X matrix element in a toggling-frame average. This is a **specific direct-matrix-element obstruction**, not a universal no-go for every use of quadratic controls plus auxiliary interactions.

## 3. All 32 addressed-link subsets, eight symmetric cycles

Every one of the `2^5=32` choices was evaluated at the **same** dimensionless `g²=1`, nominal evolution time `t=0.6`, eight palindromic `0,1,2,1,0` cycles, and the **same preexisting decoder**. Five toggled Hamiltonian segments per cycle are weighted `(1/6,1/6,1/3,1/6,1/6)`; each cycle has four phase transitions. Reported numerical fidelity uses the full 243D evolution, with no ideal full-Hamiltonian inverse.

| Number of controlled top-link sectors | Best subset (zero-based) | Recovered entanglement fidelity | Undecoded fidelity |
| ---: | --- | ---: | ---: |
| 0 | none | 0.6327493502 | 0.2546927035 |
| 1 | `4` | 0.8795137884 | 0.3834190853 |
| 2 | `0,4` | 0.9489613842 | 0.4911774032 |
| 3 | `0,1,4` | 0.9874114825 | 0.6324512439 |
| 4 | `0,2,3,4` | **0.9999971825** | 0.7012590414 |
| 5 | `0,1,2,3,4` | **0.9999973034** | 0.9999841483 |

In this **finite exhaustive search** no 3-or-fewer control subset reaches recovered fidelity 0.99, and four sectors suffice for `>0.99999`. This is **not** a proof of a minimum controller count when pulse timings, amplitudes, encoding, decoder, or circuit design can be optimized. The four-versus-five result reduces **addressed sites** and the total number of link-phase operations 20% (128 versus 160 over eight cycles), but it still needs **32 simultaneous transition events**, and relies critically on nonlocal recovery of the remaining native magnetic interaction.

## 4. Timing, amplitude, calibration, and a physically distinct drive-phase option

The detailed timing benchmarks use the **last four controlled sites `S={1,2,3,4}`**, which achieve nearly the same fidelity (`0.9999964604`) as the best four-site subset. This choice keeps an intuitive leftmost uncontrolled face. In the ideal sequence, their undecoded fidelity is only `0.7012983704`: the very high **decoded** number must not be presented as passive protection.

Physical **finite-width cubic-electric pulses** use 32 pulses of width `τ` while the source magnetic interaction remains active, giving elapsed time `0.6+32τ`. The necessary magnitude of the additional control on each addressed top link is `Ω=(2π/3)/τ` in *dimensionless energy units*. The integrated absolute control action is `32×4×(2π/3)≈268.08` link-radians for `S={1,2,3,4}`, versus `335.10` for five links.

| Pulse implementation (`S=1,2,3,4`) | Elapsed time | Required extra control or phase-slew | Recovered fidelity |
| --- | ---: | ---: | ---: |
| Instantaneous ideal pulses | 0.600 | unbounded ideal bandwidth | 0.99999646 |
| Electric control, `τ=0.002` | 0.664 | `Ω≈1047.20` per addressed link | **0.99909224** |
| Magnetic drive-phase *ramp*, `τ=0.002` | 0.664 | `|dφ/dt|≈1047.20` (drive-phase rate) | **0.99857980** |
| Magnetic drive-phase *ramp*, `τ=0.01` | 0.920 | `|dφ/dt|≈209.44` | **0.97860565** |

The electric-pulse bound and the magnetic-phase-slew bound have **different meanings**, and neither is a calibrated physical hertz/volt quantity. The drive-ramp calculation approximates every transition by sixteen midpoint-held Hermitian Hamiltonians while leaving magnetic interactions active. At τ=.002 the same computation gave recovered fidelity 0.99854884 with 4 substeps, 0.99857367 with 8, and 0.99857980 with 16: a numerical integration convergence check, **not** a rigorous error bound. Do not infer reliable analog electronics performance.

**Different error mechanisms, separate tests:**

| Eight-cycle systematic fault (`S=1,2,3,4`) | Recovered fidelity |
| --- | ---: |
| Physical cubic-electric clock angle miscalibrated by 5% | 0.99264175 |
| Physically programmable magnetic Wilson drive phase miscalibrated by 5% | 0.99721201 |
| Magnetic drive phase crosstalk 5% onto unaddressed plaquette | 0.99895518 |

These are not interchangeable: a miscalibrated electric clock acts differently on the wrapped `2→0` face-qutrit transition than a homogeneous error in the phase of a complex directed Wilson drive. The simulation implements them as different Hamiltonians.

**Noise not eliminated:** A single unconditional Kraus-complete directed Wilson jump of nominal `p=0.10` on the encoded end plaquette at the midpoint produces recovered fidelity `0.91249689` (4 controlled sectors), compared with `0.91249764` (5 controlled sectors). The control suppresses coherent magnetic evolution; it cannot correct that fault. The neighboring jump on plaquette 3 gives `0.98842489` (4 controls) or `0.99060914` (5 controls).

## 5. Engineering judgment / what is and is not native

- **Link-irrep hardware:** Explicitly generating `C3(top-link)` is a *higher-order*, representation-sensitive electric coupling. It cannot be produced by a polynomial of the single-link `C2` alone, which is blind to `3↔bar3`. A hardware realization would need new interactions or a gauge-compatible encoded representation sensor. We do not synthesize such a microscopic circuit here.
- **Gauge-fixed face-qutrit hardware:** Within the restricted face-qutrit formulation, the required unitary is an **ordinary diagonal one-qutrit `Z12` phase gate**. Superconducting qutrit experiments have implemented virtual `Z` operations in two-level subspaces, but their existence does **not** prove the published multi-plaquette magnetic couplings, addressability, or phase programmability are available. Virtual frames alone cannot change a fixed always-on Hamiltonian.
- **Driven magnetic-coupler hardware:** If `W+W†` and `i(W−W†)` can actually be driven with independent quadratures and phase-switch timing, the exactly equivalent Hamiltonian-engineered protocol can remove the need for a physical C3 pulse. This is a **conditional control compilation**, not an achieved device implementation. Driving a C-odd quadrature breaks charge-conjugation symmetry *during pulses* while respecting local gauge invariance.
- **Cost limit:** The 4-site strategy lowers site-address and aggregate phase actions by 20%, but not transition count or maximum slew/strength. At τ=.002, the 32 transition events add 10.67% elapsed time; at τ=.01 they add 53.3%. No noise on phase controllers, leakage into representations beyond `1,3,bar3`, cross-plaquette correlated faults, or local physical syndrome extraction is modeled.

## 6. Reproduce / exact and approximate claim boundary

```bash
python -m pip install sympy numpy scipy
OPENBLAS_NUM_THREADS=1 python su3_native_control_compilation.py
OPENBLAS_NUM_THREADS=1 python su3_control_mask_scan.py
python -m py_compile su3_native_control_compilation.py su3_control_mask_scan.py
```

Exact algebraic assertions verify the SU(3) cubic eigenvalues, qutrit `Z12` identification, gauge-invariant Wilson qutrit shifts, the necessary quadratic-only occupation degeneracy, and cancellation of each selected `W_i` in the group average. All 243D full-Hamiltonian propagation results use double-precision SciPy matrix exponential actions, **not** rigorous interval enclosures. The 32-subset search is exhaustive only over *this fixed protocol* and fixed decoder.

## Primary sources

1. V. Chen, B. Müller and X. Yao, *Minimally truncated SU(3) lattice gauge theory and string tension*, **Physical Review D 114**, 074502 (2026), https://doi.org/10.1103/1qt7-3ywl. The Hamiltonian source; **does not claim the new cubic control**.
2. Y. Cho et al., *Direct pulse-level compilation of arbitrary quantum logic gates on superconducting qutrits*, **Physical Review Applied 22**, 034066 (2024), https://doi.org/10.1103/PhysRevApplied.22.034066. Evidence of general qutrit gate control on specific hardware, **not an SU(3) link simulator**.
3.  *Universal Qudit Gate Synthesis for Transmons*, **PRX Quantum 4**, 030327 (2023), https://doi.org/10.1103/PRXQuantum.4.030327. Shows virtual qudit phase frames; may not refocus fixed undriven interactions.
4. E. Knill and R. Laflamme, *Theory of quantum error-correcting codes*, **Physical Review A 55**, 900 (1997), https://doi.org/10.1103/PhysRevA.55.900. Defines the stipulated error-subspace recovery criterion.
