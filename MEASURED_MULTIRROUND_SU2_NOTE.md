# Two measured flagged rounds with adaptive verification: an exact one-fault containment result

**Research status, 17 October 2026:** This is an exact binary-symplectic analysis of a complete *four-check* extraction schedule on the **five-qubit outer code**, whose five qubits are themselves embedded into the fixed-endpoint SU(2) singlet modules from this repository. The decoding algorithm uses the **actually recorded ancilla measurements**, not an externally supplied perfect syndrome. It does **not** implement the SU(2) controlled-exchange gates on native hardware or establish complete fault tolerance. Run [\`measured_multiround_su2.py\`](measured_multiround_su2.py).

The corresponding preceding investigation [\`FLAGGED_SU2_STABILIZER_NOTE.md\`](FLAGGED_SU2_STABILIZER_NOTE.md) verified individual faulted *gadgets* while assuming a perfectly measured *later* syndrome. This release replaces that oracle with an executable adaptive measurement step.

## Circuit and feedback algorithm

The outer code is the established \([[5,1,3]]_2\) five-qubit code. Its four check operators on the five **SU(2) singlet modules** are

\`\`\`text
g0 = XZZXI
g1 = IXZZX
g2 = XIXZZ
g3 = ZXIXZ
\`\`\`

For each check, prepare a fresh syndrome ancilla \(|+\rangle\) and a fresh flag ancilla \(|0\rangle\), perform four ideal controlled-Pauli data couplings, insert \(\mathrm{CX}_{m\to f}\) after data coupling 1 and after data coupling 3, measure the syndrome ancilla in **X** and the flag in **Z**, then discard/reset both ancillas. Each individual check gadget contains six two-qubit gates.

One *full flagged round* consists of **all four** check gadgets in order. The actual classical records are the four syndrome measurement bits \((b_0,\ldots,b_3)\) and four flag bits \((f_0,\ldots,f_3)\). A complete two-round baseline contains eight fresh-ancilla check gadgets and **48 abstract controlled two-qubit gates**, plus eight syndrome and eight flag readouts. These are not native SU(2) hardware gate counts.

The feedback protocol is:

1. Run two full flagged rounds and store all **actually measured bits** as \((b^{(1)},f^{(1)})\) and \((b^{(2)},f^{(2)})\).
2. If **any** flag is raised, or \(b^{(1)}\ne b^{(2)}\), run a **third complete flagged round** and record its actually measured four-bit syndrome \(b^{(3)}\).
3. For a flagged event in check \(g_i\), use the saved flag location together with the actual third-round syndrome to select the precomputed *flag-specific* Pauli recovery. The flag lookup is from the independently verified single-gadget fault enumeration in PR #11.
4. For an unflagged mismatch, use the third-round measured syndrome with the standard one-module lookup. If neither flag nor mismatch occurs, use the **second-round measured syndrome** with the same lookup.
5. Correct the data and terminate, or retain any one-module error that appeared too late to be measured consistently; such an error remains correctable by a subsequent round.

The conditional logic is executed by the script's \`recover()\` function. **It does not call \`syndrome_data()\`**, inspect hidden data-error labels, or use a noiseless syndrome oracle. The verifier separately calls \`syndrome_data()\` **after decoding**, solely to classify the residual Pauli by its ideal logical equivalence class.

## Exactly declared fault and input model

- The initial encoded data are a perfect outer-code codeword. **Exactly one** circuit fault, or none, occurs **anywhere in the first two full measurement rounds**.
- After **each** of the \(8\times6=48\) two-qubit gates, allow **all 15** nonidentity two-qubit Pauli errors on the participating pair.
- Before each individual gadget, allow an \(X\), \(Y\), or \(Z\) preparation fault on either fresh ancilla.
- At each gadget's end, allow one classical X-basis syndrome-readout bit flip, or one Z-basis flag-readout bit flip.
- This gives **\(6\times15+6+2=98\)** fault labels per gadget, **\(8\times98=784\)** mutually exclusive single-fault cases for the full baseline. Fault labels are enumerated, **not weighted by noise rates**.
- With no circuit fault, additionally test all **15 nontrivial single-module Pauli errors** already on the input data.
- The optional \`--stress\` diagnostic combines each preexisting single-module error with each of the 784 circuit faults. Its result is **not** included in the certified clean-input one-fault claim.

The Clifford propagation is exact. Each data Pauli is transported through every later controlled-Pauli gate and acquires the appropriate ancilla measurement bit-flip signs. All ancillas are freshly reset between gadgets and rounds; faults are never propagated through a falsely reused ancilla.

**Why the third round may be assumed fault-free under the enumerated model:** When the input is clean and no fault has occurred, the first two rounds give matching all-zero syndromes and no flags. Therefore every invocation of the third round is evidence that the **one** allowed fault already occurred in the first two rounds. Given the global one-fault restriction, the third round contains no additional fault. That inference is a property of this **limited model**, not a realistic claim of perfect repeat measurements under continuous noise. The code actually executes the third round's measurement gates and takes their observed syndrome bits.

## Exact results on the clean-input single-fault model

The initial GitHub Actions run established:

| Experiment | Result |
|---|---:|
| Full two-round single-fault cases | **784** |
| Completely corrected (residual stabilizer) | **744** |
| Left one *correctable* single-module Pauli-equivalence class | **40** |
| Left any uncorrectable logical error | **0** |
| Conditional third round triggered | **653** |
| Accepted after just two rounds | **131** |
| Third-round-triggered cases corrected exactly | **653 / 653** |
| Preexisting one-module Pauli, **no** new fault | **15 / 15 fully corrected** |

The residual criterion is not that the remaining *physical error weight* is literally one; it is that the residual is **equivalent modulo a code stabilizer** to an error of weight at most one. That is precisely what can be removed by an ideal future syndrome-decoding round.

Thus the protocol proves, in the stated clean-input / at-most-one-fault model:

\[
\text{one circuit fault}\Longrightarrow
\text{residual error equivalent to at most one module fault}.
\]

It does **not** prove that every one-fault execution leaves a zero-error output: 40 cases leave a correctable residual, and the result is therefore **one-fault containment** rather than unconditional complete correction.

### Why recorded bits matter

A single readout fault may produce a false syndrome while leaving data intact. A single ancilla Pauli fault may both corrupt a recorded syndrome and propagate through the remaining data couplings. The protocol contrasts two actual rounds, then performs a third if their records disagree. This directly addresses the original oracle-syndrome weakness of PR #11: **the decoder receives measured classical information, not a hidden true syndrome**.

The price is extra measurement depth; under repeated real noise this depth produces further faults, which are not covered by the one-fault case analysis.

## Combined prior-error and circuit-fault stress test: **FAIL**

The optional \`--stress\` path exhaustively tests **15 possible single-module input Pauli errors, each together with each of 784 single-circuit-fault labels**: 11,760 cases. This is a *strictly harder model* than the passing clean-input one-fault test. The first GitHub Actions run of this diagnostic returned:

| Outcome of the current adaptive decoder | Exact cases |
|---|---:|
| Fully corrected | 4,160 |
| One correctable residual | 600 |
| **Uncorrectable logical failure** | **4,760** |
| **Recovery lookup undefined** | **2,240** |
| Total | **11,760** |

These figures count uniformly enumerated *fault labels*, not physical probabilities, and the last two rows are **not passed tests**. GitHub Actions succeeds because it checks that the diagnostic is **executed and reports** its shortcomings, while requiring full success only for the explicitly supported clean-input one-fault theorem.

An example begins with a preexisting \(X\) error on module 0 and a further two-body error after data gate 2 of the first extraction gadget. The first two rounds disagree, correctly triggering a third measurement, but the standard single-error decoder leaves a logical fault; another paired fault raises a flag whose (flag, later-syndrome) lookup key is **undefined**. Repeating a noiseless syndrome by itself cannot convert a two-error pattern into a correctable one-module pattern.

The difference between **one fault from a clean input** and **one incoming data error plus one fresh fault** is therefore decisive. The latter must not be silently inferred from the former. The full-output recovery of all 11,760 cases would require a substantially stronger decoder/schedule and may be obstructed by measurement-record ambiguities; the subsequent exact transcript-collision analysis is intended to establish which failures are fundamentally unresolvable without extra measurements.

## Fault-tolerance limits and next decisions

The construction is **not** an end-to-end fault-tolerant error-correction implementation because:

1. The certified property assumes **no preexisting data error when the new circuit fault occurs**. Handling one prior error **plus** one new fault requires the separate stress analysis and perhaps additional decoding logic.
2. The third round, although physically simulated gate-by-gate in the noise-free branch, is **not itself exposed to an independent fault** under the single-fault-global-budget assumption; nor are there failures in the conditional feed-forward or recovery Pauli gates.
3. The model omits idle noise, multi-location faults, correlated two-module noise, leakage out of the two-dimensional inner singlet sectors, and imperfections of ancilla reset.
4. The controlled-\(X_\mathrm{in}\) and controlled-\(Z_\mathrm{in}\) operations are ideal encoded Clifford gates, **not yet decomposed into gauge-invariant, locality-respecting matter/rishon interactions**. The earlier SU(2) singlet analysis demonstrates only unconditionally applied single-module X/Z *at the algebraic level*, not their ancilla-controlled versions.
5. It neither predicts a memory lifetime nor claims the actual Standard Model realizes the code. The outer \([[5,1,3]]\) code and the flag idea are published prior art.

**Next milestone:** (a) analyze the combined *one preexisting data error plus one new circuit fault* diagnostic and patch observable recovery ambiguities; (b) include recovery-gate faults and repeated-measurement faults with a time-dependent decoder; (c) compile ancilla-controlled SU(2) singlet operations into native color-scalar interactions and re-enumerate physical gate faults. Only then consider quantitative noise thresholds.

## Reproduce

\`\`\`bash
python measured_multiround_su2.py
python measured_multiround_su2.py --stress  # additional combined-error diagnostic
\`\`\`

Runs with the Python 3.11 standard library plus the two existing local files [\`flagged_su2_measurement.py\`](flagged_su2_measurement.py) and [\`five_block_su2_recovery.py\`](five_block_su2_recovery.py), with no numerical approximations. GitHub Actions executes the full diagnostic version.

**Primary literature:** Rui Chao and Ben W. Reichardt, [*Quantum Error Correction with Only Two Extra Qubits*](https://doi.org/10.1103/PhysRevLett.121.050502), *Physical Review Letters* **121**, 050502 (2018); Christopher Chamberland and Michael E. Beverland, [*Flag fault-tolerant error correction with arbitrary distance codes*](https://doi.org/10.22331/q-2018-02-08-53), *Quantum* **2**, 53 (2018). Neither paper validates a physical implementation of our SU(2) matter/rishon module.
