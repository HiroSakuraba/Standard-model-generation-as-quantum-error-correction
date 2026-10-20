# Fourteen-direction branch-aware finite-S3 pulse optimizer

19 October 2026 — next step after merged PR #68. This is a **pulse-search experiment, not a demonstrated physical quantum memory**.

## Why this experiment?

PR #67 optimized *unconditional seven-single-fault logical entanglement fidelity* (`F_e`) after assuming a final ideal code-subspace projection. It found numerical gains over doing nothing, but the full three-outcome syndrome extraction could still be badly implemented. In particular, the six-dimensional orthogonal failure subspace was not an optimizer objective. PR #68 added an independent 14-direction certificate requiring code, electric-error, and failure outcomes to be correctly recognized.

This PR actually **optimizes a score on all 14 input directions** instead of evaluating that certificate only after the fact. No code projection is applied to the states during the optimization.

## Physical/encoded Hilbert spaces and target

Within the **finite S3** toy system's gauge-physical 14-state basis, initialize a neutral three-level pointer to `|0>`. The ideal *reference* operation from PR #64 sends:

| Subspace | Rank | Target physical and pointer output |
|---|---:|---|
| Code `C` | 4 | Code `C`, pointer 0 |
| Electric syndrome `V` | 4 | Code `C`, pointer 1 |
| Orthogonal complement `T` | 6 | Unchanged flagged complement, pointer 2 |

The actual operator `U(theta)` is **only** a product of PR #65's 13 gauge-compatible 42D Hamiltonian controls, repeated for a bounded number of cycles. The ideal global target unitary is never inserted into the compiled gate product.

## Objective

For each branch `b`, let `X_b` be its orthonormal input columns, `T_b=U_ideal X_b` and `n_b` its rank. Measure the **branch-coherent transfer fidelity**

`f_b = |Tr(T_b† U(theta) X_b)/n_b|²`.

This is not a ratio conditioned on observing pointer `b`; all lost or wrong-pointer amplitudes lower the score. By Cauchy–Schwarz, `f_b` cannot exceed the raw probability of the intended pointer label. Each pointer result is measured separately, so the three branches are allowed independent overall phases.

The search minimizes

`L = Σ_b (n_b/14)(1-f_b) + (η/3)Σ_b max(0,τ-f_b)² + λ Σ_k q_k(sqrt(θ_k²+ε²)-ε)`.

The first term weights all 14 physical directions fairly. The second penalizes failing to transfer a weak branch: `η=2`, `τ=0.90` are merely **search heuristics**. The third is the existing smooth charged-bus-aware cost *proxy* from PR #67; it does not measure real pulse duration.

Analytic reverse differentiation through the full 42D pulse sequence is verified against finite differences. The all-zero gate policy, two seeded random starts and a previous-depth warm start are considered, while retaining the zero-gate candidate rather than reporting a worse result as an improvement.

## Independent acceptance and negative controls

**Success is NOT defined by an optimizer's decreasing loss.** Every selected unitary is independently run through PR #68's `certify_three_outcome` on all 14 input directions with **0.99** intended-pointer and branch-coherent transfer thresholds. It also checks all seven single-error Kraus alternatives, and all 12 ordered `gR` and `Rg` double faults, including the known four `gR` failures of the *ideal one-round* decoder.

The loss minimized here is **not** the previous projected seven-error `F_e`. That `F_e` remains a separate, explicitly optimistic diagnostic with an ideal final subspace projection, so **it cannot establish native end-to-end quantum correction**. Branch-transfer certification would show that a finite idealized gauge-preserving Hamiltonian sequence *mathematically* resembles the full three-outcome syndrome-and-feedback instrument, not that a pulse sequence has been physically realized on a device.

The numerical report expands charged-mediator abstract operations and supplies assumed serial durations. An optional four-active-fermion `T1` no-jump sensitivity check is retained, but no hardware-level `T1/T2`, phase damping, fourth-level leakage, CP22 calibration, pointer misassignment, charged transport noise, meter reset or control-latency simulation is performed. Crucially, the reported duration therefore **does not include a complete device-specific physical syndrome extraction cost**.

## Run

```bash
python -m pip install numpy scipy sympy
OPENBLAS_NUM_THREADS=1 python -m unittest -v test_colour_s3_branch_aware_optimizer
OPENBLAS_NUM_THREADS=1 python colour_s3_branch_aware_optimizer.py
```

The second command runs bounded 13- and 26-control searches for AS and BS, writes `colour_s3_branch_aware_results.json` including angles and exact branch diagnostics, and prints a concise status per run. No *claim* of reaching the 0.99 certificate threshold is made in advance.

## Interpret outcomes

- If no tested depth passes the certificate, report **not achieved at these search depths**, not a Lie-algebra impossibility theorem.
- If a candidate passes, it is a concrete 42D **toy-model sequence** with a full mathematical instrument audit. Device timing, control decomposition, physical noise, readout, feedback and erasure costs are still necessary.
- If only the optimistic projected seven-error `F_e` improves while electric or flag fidelity is poor, **do not** claim full syndrome extraction.
- This research is **not** a demonstration of protection against continuous SU(3) gauge noise, Standard Model dynamics or virtual QCD fields.

Next scientific step: independently implement pointer readout/feedback/reset and device-specific noise channels (amplitude damping, dephasing, fourth-level leakage and correlated CP22 faults) with a duration-matched passive and identity comparison; compare against analog echo/refocusing.
