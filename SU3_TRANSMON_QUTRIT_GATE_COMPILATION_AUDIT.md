# SU(3) minimal-link memory: conditional transmon-qutrit gate compiler and fault-cost stress

**18 October 2026 · Finite simulator, exact block construction and floating-point time propagation.** Reproducer: [`su3_qutrit_gate_compiler.py`](su3_qutrit_gate_compiler.py); JSON and console results alongside. This does **not** establish that an SU(3) magnetic plaquette term has been built in superconducting hardware. No fabricated device parameters, gate calibrations, decoherence times or glueball masses are inferred.

## 1. Actual reference architecture

Consider **five physical superconducting transmon qutrits** in a nearest-neighbor chain (levels `0,1,2`), abstractly encoding the five gauge-reduced plaquette-flux qutrits. This differs from seven/ten/sixteen *physical SU(3) gauge links*; we are compiling the **gauge-reduced simulator**, not showing that independently controllable trivalent Gauss-law link hardware exists. The zero-exterior-flux `1⊕3⊕bar3` square-chain Hamiltonian is sourced to Chen, Müller and Yao, *Phys. Rev. D* **114**, 074502 (5 October 2026), DOI **10.1103/1qt7-3ywl**, Table VI and Eq. (3).

Experimentally established building blocks include single-transmon-qutrit two-level rotations and virtual phase advances, and some **two-qutrit** entangling gates. For example: `Universal Qudit Gate Synthesis for Transmons`, *PRX Quantum* **4**, 030327 (2023), DOI **10.1103/PRXQuantum.4.030327**; `High-fidelity qutrit entangling gates for superconducting circuits`, *Nature Communications* **13**, 7479 (2022), DOI **10.1038/s41467-022-34851-z**; `Extending the computational reach of a superconducting qutrit processor`, *npj Quantum Information* (2024), DOI **10.1038/s41534-024-00892-z**. These sources **do not demonstrate** independently programmable three-qutrit conditional SU(3) plaquette interactions.

### The missing hardware primitive

For plaquette `i`, let `l,r∈{0,1,2}` be the neighboring flux labels, zero outside the strip, and `a∈{0,1,2}` the local label. The published *directed* magnetic operator is

\[
W_i=\sum_{l,r,a}M_{l,r}[(a+1)\bmod3,a]\,P_l^{(i-1)}\otimes|a+1\rangle\langle a|_i\otimes P_r^{(i+1)},
\]

where `M_{l,r}` is the source's exact 3×3 matrix (Table VI), and boundary projectors are omitted. The Hermitian driven interaction is `B_i(φ)=e^{iφ}W_i+e^{-iφ}W_i†`. Our Z3 refocusing uses physical toggling phases `0,2π/3,4π/3` on a selected set of these **conditional** interactions.

There are **99 source-controlled directed transition types** in five plaquettes: 9 each for the two endpoints (3 neighbor contexts × 3 qutrit transitions), and 27 each for the three interior sites (9 neighbor contexts × 3 transitions). The program cross-checks all `5×243=1,215` corresponding full-Hilbert-space nonzero matrix elements against the already validated published Hamiltonian.

One interior term is **three-qutrit conditioned** (`P_l⊗G_i⊗P_r`), not a free one-qutrit drive. Virtual qutrit Z/frame updates change reference phases for subsequent drives; **they cannot change an always-on physical Hamiltonian by bookkeeping alone**. Complex `W_i` quadrature control would need to be implemented and calibrated for every participating neighbor condition, possibly through substantially longer native entangling-gate decompositions. This is the chief open engineering issue.

## 2. Exact local single-qutrit compilation, and where it stops

The script computes each local conditional unitary `U_{i,l,r}(dt,φ)=exp(+i dt B_i(φ)|_{l,r})` by 3×3 Hermitian eigendecomposition. For every such block it then performs an independent constructive unitary QR elimination using only adjacent levels `(1,2)`, `(0,1)`, `(1,2)`. This verifies the identity

\[
U=G_{12}^{\dagger}G_{01}^{\dagger}G_{12}^{\dagger}D,
\]

with up to three exact **Givens two-level rotations** and a diagonal qutrit phase `D`. It verifies reconstruction error <4×10^-12. Therefore a direct microwave **0↔2 dipole drive** is *not mathematically required* for the isolated three-level unitary. **But each Givens and phase is still conditioned on one or two neighboring qutrits.** The work does **not** synthesize these conditional primitives into calibrated two-transmon gates or provide their entangling-gate count.

One complete **second-order symmetric Trotter step** consists of electric half-steps surrounding the plaquette sequence `B0/2,B1/2,B2/2,B3/2,B4,B3/2,B2/2,B1/2,B0/2`. Time-ordering errors from distinct neighboring `B_i` are therefore explicit, not ignored. The 8-cycle five-stage symmetric echo has **40 stages**. At one Trotter step per stage the gate description requires:

| Logical resource | Count |
|---|---:|
| Full 5-qutrit conditioned site unitaries | **360** |
| Neighbor-context 3×3 controlled blocks | **2,520** |
| Of those, two-neighbor / three-qutrit contexts | **1,800** |
| Adjacent-level **controlled** Givens (upper bound, 3 per context) | **7,560** |
| Controlled diagonal qutrit phases | **2,520** |
| Global electric diagonal half-step layers | **80** |

A **controlled context block** is a logical operation in an abstract instruction set, **not a demonstrated native gate**, and the 7,560 count is an *upper bound on the 2-level factors at that layer*, **not** a proven optimal number of physical two-qutrit entangling gates. Physical decomposition could be much more expensive. Some diagonal terms and repeated half-steps could be optimized together; no hardware optimizer is claimed.

## 3. Finite digital time-step accuracy (noiseless)

Use `N=5`, `3^5=243` gauge-reduced dimensions, full published `g²=1` Hamiltonian without its global constant, the previously best four-site mask `(0,2,3,4)`, `t=0.6`, **eight** symmetric Z3 cycles and the prior explicit five-syndrome CPTP recovery. The reference is exact sparse **matrix exponential propagation** of the *piecewise constant conditioned Hamiltonians* (within double-precision numerical accuracy). For the digital compiler the only approximation is the finite Trotter product formula and ordinary floating-point arithmetic:

| Trotter steps per phase segment | ∥ψ_comp−ψ_reference∥_F | Recovered entanglement fidelity | Conditioned site unitaries | Max controlled Givens |
|---:|---:|---:|---:|---:|
| 1 | 1.08165×10^-3 | **0.999996859** | 360 | 7,560 |
| 2 | 2.70352×10^-4 | **0.999997109** | 720 | 15,120 |
| 4 | 6.75843×10^-5 | **0.999997164** | 1,440 | 30,240 |

A factor-of-two increase in digital subdivisions reduces the state error by about a factor of four, consistent with second-order splitting. It simultaneously **doubles gate opportunities for faults**, so increasing digital depth indefinitely is not a good noise-aware strategy.

Raw (no decoder) storage fidelity remains approximately **0.701** for this selective four-site protocol: the high recovered number specifically depends on the **previously assumed ideal nonlocal decoder**.

## 4. Explicit, deliberately pessimistic gate-fault/erasure channel

Because there are no measured failure rates or hardware-native gate decompositions, do not pass off a made-up fidelity number as a prediction. Instead the code includes an **explicit hypothetical CPTP effective logical channel**:

- Every counted operation independently has probability `p` of a *flagged catastrophic failure*.
- If **no** event occurs (probability `s=(1-p)^G`), apply the computed noiseless compiled channel and fixed decoder.
- If **any** event occurs, replace the logical output by maximally mixed state `I_2/2`. This depolarizing replacement is CPTP; its entanglement fidelity is exactly `1/4`.
- The resulting *toy* fidelity is `s F_noiseless + (1-s)/4`.

For the **one-step** compiler:

| Hypothetical independent fault location | Number of opportunities G | p = 10^-4: event-free probability | p = 10^-4: toy recovered fidelity |
|---|---:|---:|---:|
| Each conditioned site unitary | 360 | 0.964639 | ~0.97348 |
| Each conditioned neighbor-context block | 2,520 | ~0.7772 | ~0.8329 |
| Each controlled Givens upper-bound primitive | 7,560 | ~0.4695 | **0.60214** |

These are **three different hypothetical fault granularities**, not simultaneous faults that should be multiplied together. They do not model actual transmon leakage, amplitude damping, or a hardware-available two-qutrit gate. In particular there is **no evidence** each context or Givens would be a primitive whose physical error rate is 10^-4. The table identifies which native-gate performance would need to be measured to make the apparently high ideal fidelity meaningful. A finer two-step compilation has more sites (720), and under the `p=10^-4 per max-Givens` toy would return **0.41534** rather than 0.60214: more accurate Hamiltonian compilation can perform *worse* under gate faults.

**Coherent context drift:** a separate deterministic local test changes each context's magnetic coupling by alternating `±1%` or `±5%` (a reproducible pattern, not random transmon calibration data). At two Trotter steps per stage its recovered fidelities are **0.999997065** and **0.999996888**; the corresponding raw fidelities change because the decoder corrects some of this coherent disturbance in the selected code. Do not extrapolate these particular patterns to arbitrary crosstalk or correlated noise. Actual **T1/T2, fourth-level leakage, readout and real cross-drive crosstalk** remain unmodeled; treating the toy depolarizing branch as a validated leakage simulator would be incorrect.

## 5. Syndrome extraction is also nonlocal

The existing five-plaquette decoder reconstructs **five orthogonal two-dimensional syndrome subspaces** plus their orthogonal complement. Implementing its projective syndrome extraction with *distinguishable pointer outcomes* thus requires at least **six measurement results**, which in turn require at least **three classical readout bits** (or **two qutrit pointer ancillas**, since 3²≥6). These are information-theoretic **ancilla/readout lower bounds**, not a circuit synthesis or proof of enough ancilla gates. The syndrome projectors live in the 243-state multi-qutrit host and no local circuit, fault propagation, readout fidelity or schedule is provided. Thus even achieving the compiled unitary with high fidelity would **not** suffice to claim a working quantum memory.

## 6. Decision and missing checks

**Supported:** Exact source-local Hamiltonian block decomposition and qutrit adjacent-rotation factorization; successful finite-depth digital unitary compilation in the 243-state space; reproducible resources and an explicitly defined catastrophic-fault sensitivity model.

**Not supported:** Implementation using a demonstrated superconducting device gate set; native entangler count, T1/T2/leakage channel calibration, 3-body conditional coupling fidelity, digital gauge-symmetry preservation outside the gauge-reduced code, physical SU(3) link control, syndrome extraction and recovery circuit, code threshold, continuum mass-gap estimate or glueball measurement.

The next sensible work is **one 3-qutrit interior plaquette block** compiled into a named two-qutrit native gate (e.g. controlled-phase plus 0↔1/1↔2 rotations), verified in the full qutrit Hilbert space with a modeled fourth-level spectator, crosstalk, and its actual gate count. Only after a complete, benchmarked controlled-block compilation should large-chain simulations consume that gate-cost model.

## Reproduce

```bash
python -m pip install sympy numpy scipy
OPENBLAS_NUM_THREADS=1 python su3_qutrit_gate_compiler.py
python -m py_compile su3_qutrit_gate_compiler.py
```

The program imports the previously merged research-chain files (included in the standalone ZIP). Its `su3_qutrit_gate_compiler_results.json` records all row-level fidelity, operator counts, exact-source regression count and hypothetical fault assumptions. No interval arithmetic or calibrated hardware parameters are used.