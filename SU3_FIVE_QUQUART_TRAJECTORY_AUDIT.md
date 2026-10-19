# SU(3) full five-ququart noise trajectories and pulse optimization

Date: 19 October 2026. Continues merged PRs #39 (explicit five-qutrit circuit) and #40 (four-level two-transmon pulse, pair CPTP channels and conditional analog comparison).

## What is newly computed

1. **Bounded transmon pulse search** over an explicit 16×16 Duffing/exchange Hamiltonian with smooth detuning control; two target angles (π/2 and π), independent finer-step numerical checks, and fidelity including fourth-level leakage. This is a **synthetic two-device model**, not measured hardware, a guarantee of arbitrary-angle selective CP22, or a demonstration of transmon-native QEC. The candidate parameter budget and search seed are recorded.
2. **Full five-ququart Hilbert space** (4⁵ = 1024 basis states), with both logical codeword columns retained as a purification of a maximally mixed encoded qubit. It propagates all 40 stages of the eight-cycle five-site digital Strang schedule. Each CP22 is an *explicit* nearest-neighbor gate with coherent |22⟩–|31⟩ mixing in a two-dimensional leakage subblock, and each single-qutrit operation is extended as identity on |3⟩. The full 10,080-gate exact and 8,160-gate approximate layouts are simulated individually.
3. **Stochastic completely positive, trace-preserving (CPTP) physical noise:** after each of 40 refocusing stages, four-level bosonic loss is sampled from its complete four-Kraus set for *each of five devices*; independent random number-phase unitary diffusion samples a Gaussian channel reproducing ρ_nm → exp[−dt(n−m)²/Tφ]ρ_nm. The purified state is reweighted by the correctly sampled joint Kraus outcome probability. This is a **finite-sample channel unraveling**, not a dense 1024² density-matrix computation or a continuous-time Lindblad integration during pulses.
4. **Six-outcome noisy syndrome recovery instrument:** five orthogonal two-dimensional corrected syndrome sectors and a complementary sector are measured ideally as an *abstract nonlocal projective instrument*. Two independent pointer misclassification events and an independent *abstract* measurement-backaction replacement fault cause randomization to I/2 if recovery fails; otherwise the original correction map is used. The physical data ququarts also accumulate loss/dephasing during the assumed measurement/reset time. No ancilla coupling network or fault-tolerant extraction circuit is actually compiled.

## Most important scientific distinctions

- The **pulse optimizer** and the **full-register leakage model** are **two separately parameterized experiments**. The fitted synthetic two-transmon pulse is not substituted for arbitrary-angle CP22 gates in the full circuit. In particular, a mixing strength quoted in radians/gate is not derived from a characterized transmon coupling or from optimized pulse process tomography.
- The full-register physical noise is **stage-lumped**. Loss and phase noise are CPTP for each stage, but commuting them past the thousands of internal gates is an approximation. The script does not simulate continuous dissipative evolution during each gate.
- The applied electric diagonal acts with the previously published SU(3) phases **only** on computational qutrit basis vectors; its extension to states with at least one |3⟩ is identity. This is an arbitrary leaked-sector Hamiltonian extension. It is not a gauge-preserving four-level Hamiltonian.
- The 4⁵ system is a **gauge-fixed simulator register**. Adding an extra transmon level does not make |3⟩ a physical SU(3) representation. Leaked populations are counted as errors, not gluonic states.
- The *time budget* is serial: illustrative 580 ns per CP22 and 30 ns per local operation. No parallel-depth optimization, device-matched noise baseline, control-pulse synthesis, preparation or gauge-restoring corrections have been added. Long calculated serial durations compared with illustrative T1 values are not calibrated experimental forecasts.
- **Syndrome extraction** is represented as a CPTP abstract instrument with classical errors, replacement backaction, and idle-time noise. The requirement of ≥2 qutrit pointers or ≥3 classical bits follows from six outcomes. It is not a compiled local extraction circuit. Assumed additional entanglers increase time only and must never be presented as implemented gates.
- Entanglement fidelity estimates from **few trajectories carry large sampling uncertainty**. Standard error calculated across trajectories describes Monte Carlo noise only and does not cover model error or confidence in unknown device parameters.

## Implementation and interpretation checks

- A 1024×2 encoded map has squared Frobenius norm 2 and embeds the exact two original codewords.
- The |22⟩–|31⟩ two-state unitary is checked for unitarity; zero mixing gives exactly the ideal CP22 on the computational sector.
- The exact loss Kraus family sums to identity on four levels; trajectory sampling preserves norm and the reference coherence.
- Embedded site circuits reproduce PR #39's three-qutrit action when leakage is disabled.
- The readout confusion/backaction model is CPTP by construction. With perfect pointer/backaction and identity storage it agrees with the existing ideal decoder.
- Regression tests run first. The complete numerical study runs second and publishes JSON containing all trajectories, settings and sample standard errors. **A passing workflow does not mean that physical hardware advantage is demonstrated.**

## Commands

```bash
python -m pip install numpy scipy sympy
OPENBLAS_NUM_THREADS=1 python -m unittest -v test_su3_five_ququart_trajectories
OPENBLAS_NUM_THREADS=1 python su3_cp22_pulse_optimization.py --iterations 7 --steps 24
OPENBLAS_NUM_THREADS=1 python su3_five_ququart_trajectories.py --shots 4 --leakage .01 \
  --t1-ns 30000 --tphi-ns 45000 --pointer-error .001 \
  --backaction-error .001 --extraction-entanglers 24
```

These are **synthetic scenarios**, not recommended physical parameter values.

## Initial numerical validation and sensitivity design

An initial complete-trajectory run on this PR executed eight regression tests
and all eight cycles on both circuit variants. In the **illustrative serial**
configuration (580 ns per CP22, 30 ns per local gate, 30 µs T1 and 45 µs Tφ),
mean recovered entanglement fidelity was effectively 0.25 for both variants.
This is the maximally mixed-qubit recovery floor, consistent with decay over
approximately 6–8 milliseconds of serial gate execution. It **does not**
establish a model-independent quantum memory no-go result.

We specifically distinguish **transient** fourth-level occupation before
and after each stage's loss/dephasing from *final* leakage. Relaxation may
empty the fourth level by the final readout while irreversibly erasing the
encoded information. The latter cannot be inferred from a small final
leakage number.

The workflow now adds a separate **deliberately speculative optimistic arm**
with 25 ns CP22, 1 ns local gates, 1 ms T1 and Tφ, and lower leakage
coupling. These are stress-test assumptions, **not** plausible measured
parameters supplied by the cited papers. Both arms run exact 40-CP22 and
approximate 32-CP22 layouts. The standard serial arm uses four sampled
trajectories per compiler because it tests a strong-decay regime. The
speculative high-coherence arm uses 24 trajectories per compiler. Even
24 samples may be insufficient to rank circuits reliably; report sample
standard error and never infer a precise physical threshold from these
toy parameter scenarios.

The pulse optimizer retains the *independently verified finer-time-step
baseline* whenever the bounded random search does not improve fidelity.
An apparently different optimum must not be called a gain unless it
outperforms the same benchmark under equal time discretization.

## Research decisions

A useful physical claim requires (i) calibrated native pulse or robust optimized multi-control waveform across *all* required CP22 angles, including local frames and spectators, (ii) pulse-time Lindblad quantum dynamics or a benchmarked lower-cost trajectory method per gate, (iii) a gauge-compatible treatment of leaked sectors, (iv) an actually synthesized and noisy six-outcome syndrome extraction network, and (v) a **same-device, same-time, equally noise-budgeted uncontrolled and analog baseline**. A low full-circuit fidelity under present scenarios is a meaningful cost warning, not a general impossibility proof.

Primary literature context: [Universal Qudit Gate Synthesis for Transmons](https://doi.org/10.1103/PRXQuantum.4.030327) (2023, numerical compilation); [Exploring ququart computation on a transmon using optimal control](https://doi.org/10.1103/PhysRevA.108.062609) (2023, experiment); [Leakage in restless quantum gate calibration](https://doi.org/10.1103/PhysRevA.108.022614) (2023, multi-level relaxation models); [Qudit Dynamical Decoupling on a Superconducting Quantum Processor](https://doi.org/10.1103/PhysRevLett.134.050601) (2025, experimental decoupling). These do **not** calibrate our synthetic SU(3) construction.
