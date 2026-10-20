# Finite-S3 recovery: longer-pulse continuation and sector-capture warmup

October 19, 2026. Follow-up to merged PRs #68–69.

## Research question

Can an actual product of gauge-compatible, charged-mediator-aware *abstract*
Hamiltonian pulses approach the **complete three-outcome quantum recovery
instrument**, rather than just improving one logical fidelity figure under
an unrealistically free final code projector?

PR #69 tried 13 and 26 tunable controls. Six dedicated unit tests passed, but
none of its four experiments met the independent three-outcome acceptance
criterion (0.99 per branch). The weakest coherent branch transfer reached
approximately 0.269 (AS, 26 controls) and 0.298 (BS, 26 controls), while the
projected seven-error logical fidelity remained below the 0.671111 no-gate
baseline. These are **limited numerical search outcomes**, not impossibility
theorems.

## New numerical approach

**Stage A: capture a complete physical syndrome sector.** For each of the
three input branches `b`, compute

`q_b = || T_b† U(theta) X_b ||_F² / rank_b`.

This counts probability ending **in both the correct pointer state and the
correct system sector**. A coherent target overlap from PR #69 can be nearly
zero even if some useful amplitude has reached the right subspace. The
sector-capture surrogate can have a different optimization gradient and
helps find candidate paths out of that plateau.

However, `q_b` ignores unwanted rotations *within* the code and flagged
subspaces, so a high `q_b` is not evidence of protected quantum information.
The new optimizer explicitly treats it only as a **warmup objective**.

**Stage B: protect the actual target map.** From each Stage-A endpoint,
optimize PR #69's *coherent branch transfer*

`f_b = |Tr(T_b† U(theta) X_b)/rank_b|²`

with the balanced penalty for weak branches. This rejects in-sector
scrambling even when the pointer label is right. Both objectives have analytic
reverse-adjoint gradients, with Stage A independently checked against
centered finite differences.

**Search continuation:** Test fixed cyclic schedules of 13, 26, 52 and
104 abstract Hamiltonian rotations for AS and BS. Use deterministic seeded
starts, the strict zero-gate control, and a copy of the preceding shallower
sequence extended by zero-angle pulses. Appending zeros preserves the
*entire* physical unitary, not just a summary score.

Select the candidate with the highest **worst branch coherent transfer**,
breaking ties with the rank-weighted average. This explicitly protects the
weakest extraction outcome against being sacrificed to improve another.
The previous-depth candidate is retained at deeper searches, so the
worst-branch metric cannot worsen across these depths. This is a property
of the *selection method*, not proof of algorithmic convergence.

## Independent verification and physical limitations

Every *selected* pulse product is sent to the separate PR #68
`certify_three_outcome` implementation, testing all 14 physical input
directions in the 42-dimensional system-pointer space. It checks intended
pointer probabilities, branch-coherent transfer, independent seven-Kraus
single-fault scores and ordered `gR/Rg` counterexamples. The new code
also asserts that the independently measured sector probability matches
Stage A's diagnostic and upper-bounds coherent target transfer.

Acceptance remains **0.99 for each of the three branches**. A better loss,
better sector capture, more pulses, or a passing CI run *does not itself
mean* a decoder has been demonstrated. The searched pulses form real finite
matrix products but their device-native transport, pointer interaction,
calibration and simultaneous scheduling are not established. The original
unconditional projected seven-error `F_e` still assumes an ideal final
code subspace selection; we report it as an **optimistic diagnostic** and
never let it substitute for complete-instrument certification.

Operation count, serial duration and four-active-fermion no-jump lifetime
figures remain **hypothetical**. This code does not yet model a specified
CP22 device, fourth-level leakage, amplitude damping, dephasing, reset,
measurement assignment errors, mediator loss or feedback latency. Neither
does it model continuum SU(3), QCD glueballs or virtual QCD excitations.

## Execute

```bash
python -m pip install numpy scipy sympy
OPENBLAS_NUM_THREADS=1 python -m unittest -v test_colour_s3_continuation_search
OPENBLAS_NUM_THREADS=1 python colour_s3_continuation_search.py
```

The benchmark writes `colour_s3_continuation_results.json` containing all
selected pulse angles, independent syndrome outcomes, cost assumptions,
candidate metrics and optimizer convergence states. CI uploads it as an
artifact. The benchmark is intentionally bounded in CPU work. Failure to
reach 0.99 at 104 controls means only **not achieved by this search**.

## Decision after this experiment

If full coherent transfer remains low, investigate the accessible control
alphabet, reparameterization, pulse ordering and local gauge/ancilla coupling
before spending large-scale computation on longer blindly repeated patterns.
If higher-fidelity candidates emerge, next build a complete device-specific
open-system instrument and count readout, feedback and reset in a
duration-matched comparison against both passive memory and analog
refocusing. **No memory advantage is claimed in this PR.**

## First completed GitHub Actions benchmark (19 October 2026, Eastern)

The first complete dedicated [CI run](https://github.com/HiroSakuraba/Standard-model-generation-as-quantum-error-correction/actions/runs/39351069597) passed **7/7 new tests** and archived all eight selected schedules/independent instrument certificates. The following are **actual bounded numerical outcomes**, not upper bounds or hardware fidelity estimates. Values rounded to three decimals.

| Cyclic controls | AS minimum coherent branch | BS minimum coherent branch | AS projected seven-fault F_e* | BS projected seven-fault F_e* |
|---:|---:|---:|---:|---:|
| 13 | 0.002 | 0.002 | 0.630 | 0.639 |
| 26 | 0.284 | 0.272 | 0.306 | 0.295 |
| 52 | 0.458 | 0.473 | 0.479 | 0.478 |
| 104 | **0.578** | **0.556** | **0.579** | **0.555** |

`*F_e` here is an **optimistic diagnostic** that assumes an ideal final logical-code projection. No selected circuit passed the independent full-instrument 0.99 threshold. The zero-gate projected reference under this toy seven-fault channel is **0.671111**.

At 104 controls the expanded abstract primitive counts were **136 (AS)** and **136 (BS)**; assumed serial durations including a nominal readout were ~**3.992 µs (AS)** and ~**4.002 µs (BS)**. These numbers do **not** include a fully calibrated pointer, charged mediator, feedback/reset circuit and must not be interpreted as physical gate counts on a specific device.

**Interpretation:** the two-stage search makes a measurable difference to the *worst branch of full instrument transfer* as depth grows, but improving that metric has not produced higher recovered logical fidelity than doing nothing. This experiment rules out **these eight selected numerical schedules** as demonstrated full-instrument recoveries; it does not rule out better pulse orderings, devices, control libraries, or deeper schedules.
