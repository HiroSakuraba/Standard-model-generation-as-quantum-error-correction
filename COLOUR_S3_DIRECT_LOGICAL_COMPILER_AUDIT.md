# Direct logical-fidelity recovery optimization — finite S3 experiment

19 October 2026. Follows merged [PR #66](https://github.com/HiroSakuraba/Standard-model-generation-as-quantum-error-correction/pull/66). That experiment produced actual 13–208-pulse Hamiltonian schedules and found that maximizing resemblance to the ideal 42D syndrome target did **not** maximize the quantum memory's corrected logical fidelity: at 208 controls, the AS and BS circuits scored ~0.557 and ~0.563 on the synthetic seven-error logical channel, **below the identity/do-nothing channel's ~0.6711**. The new work directly optimizes the corrected information rather than cosmetic target overlap.

## 1. Precisely which error channel?

PRs #62–64 established ideal four-logical-state active recovery within the Gauss-physical 14D space of four matter qutrits and one finite signed-permutation-S3 gauge link. There are six allowed matter-permutation faults from either proper order-six group AS or BS, and a link-electric phase \(R=I-2P_{\rm standard}\).

For a synthetic single-step stochastic channel, take electric error probability \(p_e\) and each of six exchange errors with probability \((1-p_e)/6\). **Every exchange g fixes all four encoded columns C**, so \(gC=C\) exactly. This reduces the numerical fidelity optimization *without approximation* to just two input blocks:

\[
X_{\rm clean}=C\otimes|0\rangle_p,\qquad
X_{\rm electric}=RC\otimes|0\rangle_p.
\]

Thus the optimizer carries **eight** 42-dimensional state vectors instead of a full 42×42 unitary. The reduction is **not valid for two sequential exchange/electric faults** of type \(gR\) or \(Rg\), or an arbitrary gauge-compatible noise channel.

## 2. Objective: unconditional entanglement fidelity after the actual circuit

Let \(U(\theta)\) be the product of tunable, independently gauge-compatible control pulses from PR #65 (same 13 controls per cycle as PR #66). Pointer outcomes 0 and 1 denote recognized recovery branches; pointer outcome 2 is the unrecoverable flag.

For \(j\in\{\rm clean,electric\}\) and recognized pointer outcome \(b\in\{0,1\}\), define

\[
a_{jb}=\operatorname{Tr}\left[
(C\otimes|b\rangle)^\dagger U(\theta)X_j
\right].
\]

The **true unconditional logical entanglement fidelity** of the recovery subchannel (erasure/flag outcomes orthogonal to the logical output) is

\[
F_e(\theta)=\sum_{j,b}\frac{w_j}{16}|a_{jb}|^2,\qquad
(w_{\rm clean},w_{\rm electric})=(1-p_e,p_e).
\]

This expression follows from the trace-of-Kraus definition and treats the two pointer measurement outcomes **incoherently**, with no falsely imposed common phase between physically distinct classical outcomes. The objective gives zero credit for pointer 2, and zero credit to physical output outside the logical code C despite a recognized pointer label. Such events are *not postselected away or renormalized*.

The implementation **independently recomputes this score using all seven original physical Kraus operators** and compares it to the two-block reduced objective for every candidate reported. Both calculations must agree to numerical precision.

At \(p_e=0.37\), a circuit doing **nothing** has

\[
F_{e,\mathrm{identity}}=(1-p_e)+p_e/9
=0.671111\ldots,
\]

because \(C^\dagger R C=-I_4/3\). By contrast, the unsynthesized ideal global PR #64 syndrome/feedback has \(F_e=1\). Both baselines are asserted by dedicated regression tests.

## 3. Differentiable pulse compiler and a non-degradation fallback

All controls remain the same finite-S3 gauge-preserving Hamiltonians: local onsite exchange, charged-bus-decomposable dressed exchange, link-electric projector, opposite-vertex local two-matter equal-colour projector, eight neutral-qutrit pointer drives, and the assumed pointer+link electric entangler.

The optimizer derives an **analytic reverse/adjoint gradient of the actual channel fidelity**, not of 42D unitary-target overlap. Denote \(F=U(\theta)[X_{\rm clean},X_{\rm electric}]\). At the final time, each weighted complex trace amplitude contributes a terminal adjoint term \(w_j a_{jb}(C\otimes|b\rangle)\). Back-propagating those terms through the variable-angle unitaries yields exact first-order derivatives with the same cost as a small collection of state propagations. These gradients are checked against centered finite differences, including the optional cost penalty.

The regularized objective is

\[
\mathcal L(\theta)
=1-F_e(\theta)+\lambda\sum_k q_k\bigl(
\sqrt{\theta_k^2+\varepsilon^2}-\varepsilon
\bigr).
\]

Here \(q_k=5\) for dressed cross-link rotations (each maps to five charged-mediator two-register operations), \(q_k=2\) for the pointer/link entangler, and \(q_k=1\) for the others. These are **heuristic complexity weights**, *not* measured pulse duration, energy, or experimental fault rates. The JSON also records true abstract operation counts and a separately parameterized **hypothetical** serial gate-duration estimate via the PR #66 expansion.

The compiler uses bounded L-BFGS-B, deterministic small/medium random starts, the exact all-zero initialization, and warm starts from shallower pulse schedules. **Identity is always explicitly included as a selectable, zero-control policy**; no failed numerical optimization can be reported as outperforming it. It also incurs no fictitious readout duration if selected.

A selected control sequence can be better than doing nothing while still not approximating the *specific ideal 42D projector unitary*. Direct logical optimization only enforces good code recovery for the **specified restricted single-fault channel**; the diagnostic output separately records complete 14-column isometry overlap, full 42D unitary overlap, and all three original syndrome-branch metrics. Improving F_e alone does not establish proper three-outcome flagged classification, arbitrary sequential-fault correction, or a universal syndrome extractor.

## 4. Validation and experiment

Twelve dedicated tests cover the clean/electric input-state reduction, no-postselection single-step channel, ideal and do-nothing baselines, independent seven-Kraus numerical agreement, analytic-gradient checks at 13 and 26 controls with/without effort regularization, finite control angles, deterministic repeatability, identity/no-readout fallback, input validation and depth warm starts.

The benchmark enumerates both AS and BS at **13, 26, 52, 104** tunable controls, with an explicit maximum optimizer-iteration budget. For each candidate it archives its full floating-point angles and names, initialized-start ledger, unconditional seven-error logical F_e and gain over do-nothing, expanded charged-mediator operations and transports, tentative timing and loss-sensitivity figures, and all three ideal-target diagnostic overlaps.

**The benchmark must not be called a success until its verified F_e exceeds 0.671111 at a useful depth and survives independently simulated device noise and fair resource accounting.** If no improvement is found, the result is an honest numerical negative for the tested cyclic template/search budget, not a Knill–Laflamme theorem or a lower bound on all circuits.

The previously derived four-dimensional code and its ideal recovery remain exactly correct against the **seven individual faults**. Neither this compiler nor its timing assumptions model nonperturbative QCD, virtual gluon or sea-quark loops, the original nine-quark 21D/65-syndrome code, realistic electronic charge transport, pointer/reference/spectator loss, syndrome assignment errors, or correlated faults occurring during recovery. The opposite-vertex equal-colour projector is only invariant under the finite S3 subgroup, not full continuous SU3.

### Reproduce

\`\`\`bash
python -m pip install numpy scipy sympy
OPENBLAS_NUM_THREADS=1 python -m unittest -v test_colour_s3_direct_logical_compiler
OPENBLAS_NUM_THREADS=1 python colour_s3_direct_logical_compiler.py
\`\`\`

A separate GitHub Actions workflow preserves actual schedules and no-postselection channel evidence as a JSON artifact. Only the green scientific and repository-wide workflows authorize the PR merge.
