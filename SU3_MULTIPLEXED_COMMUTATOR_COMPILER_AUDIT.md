# SU(3) interior magnetic interaction: shared-control multiplexor and commutator compilation

**18 October 2026. Status: numerically verified exact logical 27-dimensional gate synthesis, and a separate error-bounded approximate option.** This file documents a follow-up to the previous independent nine-context |22> two-qutrit compiler (PR #36). No claims of calibrated experimental gates, continuum SU(3), glueball masses, or an operational Standard Model quantum memory follow from this calculation.

## Problem and primary input

The published source is Vincent Chen, Berndt Müller and Xiaojun Yao, *Minimally truncated SU(3) lattice gauge theory and string tension*, **Physical Review D 114, 074502 (published 5 October 2026)**, https://doi.org/10.1103/1qt7-3ywl. Its square-chain minimal link cutoff consists of the representations `1`, `3`, and `anti-3`. We use the independently transcribed **journal Appendix B Table VI** from `su3_three_plaquette_chain.py`, cross-checked against the previous 360-entangler compiler's own transcription.

For an *interior* target qutrit T with adjacent neighbors L and R, each of the nine `(l,r) ∈ {0,1,2}²` contexts has a published magnetic block `M_(l,r)`. With directed fundamental Wilson shift `W` and an externally chosen phase `phi`, the 3×3 Hamiltonian block is

`B_(l,r)(phi) = exp(i phi) W_(l,r) + exp(-i phi) W_(l,r)†`.

We compile the full block-controlled unitary `U = direct_sum_(l,r) exp(i dt B_(l,r))` on the 27-dimensional register `(L,T,R)`. The overall sign on `dt` is conventional. The independently published matrices are the *input*, not fitted to the compiled output.

## 1. Exact factorization across the nine contexts

Let `U_lr = exp(i dt B_(l,r))`. Define one **unconditional** target operation `U_00`, two left-neighbor corrections `L_l = U_l0 U_00†` for `l=1,2`, and two right-neighbor corrections `R_r = U_0r U_00†` for `r=1,2`.

The remaining two-neighbor corrections occur only when **both** neighbors are nonzero:

`C_lr = U_lr (R_r L_l U_00)†`, for `(l,r) ∈ {(1,1),(1,2),(2,1),(2,2)}`.

The chronological circuit applies `U_00`, then the controlled `L_l`, then the controlled `R_r`, then the controlled `C_lr`; its action on any of the nine sectors is **exactly `U_lr`**. This avoids duplicating the entire target-level decomposition nine times. Every factor is unitary and has determinant one (SU(3)) within numerical tolerance.

## 2. Constructive SU(3) commutator and two-neighbor control

For each residual `C`, Schur-decompose `C = V diag(lambda_0,lambda_1,lambda_2) V†`, with `lambda_0 lambda_1 lambda_2 = 1`. Define the cyclic three-level shift `S|j⟩=|j+1 mod 3⟩` and

`A = V S V†`, `B = V diag(1, lambda_1^-1, (lambda_1 lambda_2)^-1) V†`.

Then **exactly** `A B A† B† = C`, because conjugation by `S` cyclically permutes the three diagonal entries. The circuit

`controlled_R(B†), controlled_L(A†), controlled_R(B), controlled_L(A)`

in *chronological* order applies `C` only when **both** prescribed neighbor values occur; on all other sectors the paired operations cancel. It requires no ancilla qutrit and no direct L–R coupling.

A singly neighbor-controlled arbitrary U(3) target unitary needs **two** entangling controlled phase gates, plus local qutrit operations: diagonalize its three eigenphases; the common eigenphase is a *controller-only* diagonal phase gate, while the other two relative eigenphases are selective CP interactions. Every selective CP is routed onto a fixed `|22⟩` CP primitive by adjacent-level 0↔1/1↔2 swaps before and after. An arbitrary local U(3) is itself decomposed into at most **three adjacent-level SU(2) mixers plus one diagonal phase frame**, and the verifier reconstructs each local U(3) before using it.

**Entangler count, full nine-context interaction:**

- Four singly controlled `L_l/R_r` corrections × two CP = **8 CP22**;
- Four doubly controlled residuals × four singly controlled `A/B` factors × two CP = **32 CP22**;
- Unconditional target operation: no entanglers;
- **Total: 40 CP22**, down from **360 CP22** in PR #36 (exactly **9.0× fewer**).

There are **264 separately counted local qutrit operations**, broken into 203 adjacent-level SU(2) operations (including 80 explicit nearest-level swaps) and 61 diagonal qutrit phase frames. PR #36's nine-independent-context construction used 1,728 local operations; comparison assumes the same local-operation accounting conventions.

## 3. Verification, not just an algebraic gate count

The executable `su3_multiplexed_commutator_compiler.py` constructs the **entire sequential 27×27 matrix** by multiplying every compiled CP22 and every local operation. It compares that matrix in operator norm against an independent direct blockwise exponential. All four sampled `dt,phi` pairs pass:

| `dt` | `phi` (rad) | 27D unitary operator-norm difference |
|---:|---:|---:|
| 0.19 | 2π/3 | 7.26×10⁻¹⁵ |
| 0.08 | 0 | 6.65×10⁻¹⁵ |
| 0.30 | 0.23 | 1.37×10⁻¹⁴ |
| 0.047 | −1.1 | 1.43×10⁻¹⁴ |

These checks are **double-precision numerical tests**, not interval-certified equalities. The algebraic identity is constructive, but device pulse response has not been characterized.

## 4. A 32-entangler *approximation* clears the tenfold target at a cost

The exact factorization allows omission of one or more **residual** doubly controlled corrections. For the reference `dt=0.19`, `phi=2π/3`, omitting `(l,r)=(1,2)` gives 32 CP22 (**11.25×** fewer than 360), with:

- Single-interaction 27D process entanglement fidelity against the exact target: **0.9994898873**;
- Worst-case (operator-norm) unitary difference: **0.0945825483**;
- Local operations: **212**.

Omitting both `(1,2)` and `(2,1)` uses 24 CP22, process fidelity 0.9989799048 and the same worst-case norm difference because they affect mutually orthogonal neighbor sectors. Omitting all four doubly controlled residuals uses only eight CP22, process fidelity **0.9975871187**, but a worst-case norm error **0.1059583722**. None of these approximate circuits are exact published-Hamiltonian implementations.

The 32-entangler option was repeated across all four test points. Its worst operator-norm difference among them was **0.1463** (at `dt=0.30`, `phi=0.23`), with process fidelity **0.998715**. These single-interaction values **cannot simply be multiplied or extrapolated to a long refocusing sequence**; systematic coherent error can accumulate across repeated blocks.

## 5. Matched phenomenological fourth-level noise model

We reused the previous compiler's **identical phenomenological** 64-dimensional three-ququart propagation: each CP22 is a 16D unitary with coherent `|22⟩↔|31⟩` leakage mixing and a separate spectator conditional phase. Local operations are perfect and leave fourth levels untouched. In this intentionally simplified model, at leakage-amplitude parameter `0.01` and spectator phase `0.002` per CP:

| Implementation | CP22 | Leakage probability | Unconditional process entanglement fidelity |
|---|---:|---:|---:|
| PR #36 nine separate blocks | 360 | 0.00097949 | 0.99102580 |
| This 40-CP multiplexor | **40** | **0.00051271** | **0.99938436** |

For a more severe test `(leak=0.03, spectator phase=0.006, CP-angle bias=2%)`, the new method has leakage **0.004569** and process fidelity **0.994377**. All fault parameters are invented simulation parameters; no transmon measurements are implied. The numerical model omits T1/T2 decay, pulse spectral bandwidth, dressing, residual link gauge-constraint leakage, and readout/noisy recovery. In particular, a shorter simulated circuit is not proof of a practical gate.

## 6. Timing illustration, not a physical benchmark

If the same **illustrative** serial durations used in PR #36 (580 ns per arbitrary-angle CP22 and 30 ns per one-qutrit operation) are imposed mechanically:

- Earlier circuit: `360×580 + 1728×30 = 260,640 ns` (260.64 μs).
- New exact circuit: `40×580 + 264×30 = 31,120 ns` (31.12 μs).

This is an **8.38× reduction in the simplistic serialized-duration model**, without justified CP-angle timing or T1/T2. The 580 ns number was taken from a *different* experimental qutrit entangler, not demonstrated for this arbitrarily selective CP22. No real hardware time or fidelity estimate is established.

## 7. Next validation priorities

1. Optimize the exact 40-CP circuit further **without sacrificing arbitrary-angle correctness**; do not call 32 CP exact. Consider shared-control cancellation, control-label Gray-code traversal, and numerical two-qutrit circuit synthesis with independently verified operator norms.
2. Calibrate or at least specify a concrete **fixed-angle/echoed native two-qutrit gate** with realistic pulse schedules, instead of an assumed arbitrary-angle CP22.
3. Model T1/T2 leakage, dressed fourth-level dynamics, crosstalk, state preparation and measurement on the full gate sequence. Report implementation-level error and time, not only ideal process fidelity.
4. Retest the optimized blocks **inside the full five-plaquette refocusing sequence**. A single-interaction gate process fidelity does not establish long-time memory fidelity.
5. Before a quantum-memory claim, compile the previously idealized **five-syndrome recovery** to physical ancilla measurements and characterize logical fault propagation.

## Reproduce

```bash
python -m pip install numpy scipy sympy
OPENBLAS_NUM_THREADS=1 python su3_multiplexed_commutator_compiler.py \
  --output su3_multiplexed_commutator_results.json
python -m py_compile su3_multiplexed_commutator_compiler.py
```

Requires the previous repository files `su3_selective_cp_compiler.py` and `su3_three_plaquette_chain.py` in the same directory. The script verifies their published transition table agreement, constructs all nine contexts, checks full compiled matrix products, and performs a matched fourth-level toy-noise propagation. JSON results and the recorded console output are committed alongside the script. Earlier research implementations are **not changed**.
