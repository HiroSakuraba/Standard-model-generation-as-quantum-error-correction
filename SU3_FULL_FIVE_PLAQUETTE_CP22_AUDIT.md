# SU(3) five-plaquette CP22 refocusing: digital-circuit audit

**19 October 2026.** Source branch: main, after the integration in PR #38.
**Status:** An executable proposed experiment, **not** proof of a working memory,
calibrated CP22 pulse, physical leakage rate, or physical Yang–Mills spectroscopy.

## What the executable actually does

\`su3_full_refocusing_cp22.py\` takes the published minimal-link SU(3) source
Hamiltonian and five-face encoded codewords from this repository. It reproduces
the *existing* symmetric Strang/Trotter decomposition for the eight-cycle
four-addressed-site Z3 refocusing protocol. Each interior three-qutrit
conditioned magnetic evolution is replaced by the **expanded** circuit from
\`su3_multiplexed_commutator_compiler.py\`. A separately synthesized boundary
multiplexor uses one controller with four CP22 gates. No free conditional
three-qutrit operator remains inside the magnetic part.

The electric diagonal half steps are still ideal diagonal operations.
So are the switching of clock/phase controls, native CP22 pulse synthesis,
code preparation, syndrome extraction, ancilla reset and logical recovery.
They are **not** included in the following gate/time numbers.

The 32-entangler approximation is tested through sequential circuit action,
not inferred from its one-interaction process fidelity.

## Structural counts at one Trotter subdivision per stage

There are 8 cycles × 5 symmetric toggling stages = **40 stages**.
Each stage uses site order \`0,1,2,3,4,3,2,1,0\`.
Thus the 3 interior plaquettes (1,2,3) each appear twice per stage,
for 240 interior magnetic blocks total. The two boundary plaquettes
appear 3 times per stage, for 120 boundary blocks total.

| Variant | Interior CP22 | Boundary CP22 | All CP22 |
|---|---:|---:|---:|
| Exact 40 per interior | 240 × 40 = 9,600 | 120 × 4 = 480 | **10,080** |
| Approximate 32 per interior | 240 × 32 = 7,680 | 120 × 4 = 480 | **8,160** |

Each additional Strang subdivision multiplies these *unoptimized* counts
by the subdivision count. Consecutive local gate/frame cancellations and
physical parallel scheduling are not optimized by this module.

At an **illustrative, uncalibrated** 580 ns per CP22, the CP22 portions
alone occupy 5.8464 ms (exact) or 4.7328 ms (approximate) in a serial
schedule. This is neither an experimental estimate nor a hardware lower
bound: the assumed arbitrary-angle CP22 primitive is not a demonstrated
calibrated gate and parallelization has not been analyzed.

## Regression criteria

1. Five independent 243-dimensional local block comparisons, including
   both boundaries, must match direct source-block exponentials.
2. Full exact compiled state evolution must reproduce
   \`su3_qutrit_gate_compiler.strang_segments\` on the **two encoded
   basis states** to Frobenius norm < 3e-8.
3. Unitary/isometry checks and the full 10,080/8,160 gate counts must pass.
4. Report direct ideal-toggled-Hamiltonian recovery fidelity,
   uncompensated native-Hamiltonian recovery, exact full circuit
   recovery, approximate full circuit recovery, coherent state-vector
   differences and a worst-case encoded-input deviation (singular value
   of the 243×2 difference).
5. Run the GitHub Actions workflow and record real JSON outputs before
   asserting any new numerical fidelity results. The benchmark runs
   numerical double precision; it is not interval certified.

## Risk model, scope and physical limitations

A separate *logical-level* replacement-risk scenario assumes independent
flagged faults. If **no** flagged event occurs, the logical output is the
calculated noiseless output. Otherwise it is replaced by the completely
mixed qubit \`I/2\`. The resulting channel is completely positive and
trace preserving (CPTP) by construction. Its entanglement fidelity is

\`F_e = s F_noiseless + (1-s)/4\`,

where \`s\` is the event-free survival probability. This is a **pessimistic
phenomenological bound/scenario**, not a native physical CP22, T1/T2,
fourth-level-leakage or noisy-decoder simulation. It also cannot be used to
deduce a hardware threshold without specifying a physical device.

Using an ideal controlled fidelity of about 0.999997, the report's
uncontrolled fidelity about 0.632749, other gates perfect and 10,080
independently flagged CP22s, this model needs error probability below
roughly **6.67×10^-5** per CP22 to outperform the uncontrolled baseline.
This is a *conditional calculation*, not an observed threshold.

The script also supplies *illustrative* gate times, local-gate fault
probabilities, pointer-readout fault probabilities and exponential
replacement hazards (labeled by assumed microsecond lifetimes). Such
hazards must **not** be represented as a simulation of physical qutrit
amplitude damping or dephasing.

**The digital-versus-analog distinction matters.** The CP22 circuit
costs concern a digital simulator implementing magnetic interactions.
An analog model with genuine native magnetic evolution and physically
addressable cubic-Casimir control may have a completely different cost.
Therefore these counts do **not** establish a mandatory physical cost of
the underlying gauge theory or of every refocusing realization.

## Still required before a claim of physical advantage

- Specify and optimize a device Hamiltonian implementing arbitrary-angle
  selective CP22, including gate durations, pulse bandwidth and topology.
- Full five-site \`|3>\` leakage/crosstalk trajectories and trace-preserving
  physical qutrit damping/dephasing channels.
- Account for clock pulses, electric phases, parallel gate depth and drift.
- Implement or explicitly lower-bound *native noisy* syndrome extraction,
  ancilla preparation/reset and recovery.
- Compare to an equally hardware-budgeted uncontrolled baseline, not
  solely ideal uncontrolled evolution.
- Vary Trotter subdivision, coherent approximation residual selection,
  and gauge-breaking errors.

### Reproduce

\`\`\`bash
python -m pip install numpy scipy sympy
OPENBLAS_NUM_THREADS=1 python -m unittest -v test_su3_full_refocusing_cp22
OPENBLAS_NUM_THREADS=1 python su3_full_refocusing_cp22.py --rounds 1 \
  --output su3_full_refocusing_cp22_results.json
\`\`\`

The workflow \`.github/workflows/su3_full_refocusing_cp22.yml\` runs both
commands and uploads the JSON. Tests must complete before the PR is
described as validated.
