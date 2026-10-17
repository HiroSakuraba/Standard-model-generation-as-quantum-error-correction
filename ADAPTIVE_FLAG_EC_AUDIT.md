# Source-informed adaptive flag recovery for five SU(2) singlet modules

**17 October 2026 — scope:** reproducible **abstract Clifford-level** adaptive recovery experiments on the known \([[5,1,3]]_2\) five-qubit code, with each qubit represented by a two-dimensional SU(2)-gauge-singlet matter/rishon module from this repository. This is neither a new code nor an implementation of the Standard Model. See [\`paper_style_flag_ec.py\`](paper_style_flag_ec.py).

## 1. Primary-source check and an important scientific correction

The primary reference is Rui Chao and Ben W. Reichardt, [*Quantum Error Correction with Only Two Extra Qubits*, *Phys. Rev. Lett.* **121**, 050502 (2018)](https://doi.org/10.1103/PhysRevLett.121.050502), especially **Figure 1(c)** and the numbered error-correction procedure on pages 1–2 of the [accepted manuscript](https://link.aps.org/accepted/10.1103/PhysRevLett.121.050502). Their published **five-qubit fault-tolerant adaptive control flow** is:

1. Measure each stabilizer using a flagged gadget in a specified order.
2. At the *first* observed flag **or** nontrivial syndrome bit, **stop** the remaining flagged checks.
3. Extract **all four** stabilizer syndromes again using *unflagged* gadgets.
4. If the flag fired, correct from the allowable flagged-error set for that particular check; otherwise correct using a standard weight-one syndrome lookup.
5. If all initial checks are trivial, **apply no correction**.

This source provides a constructive full fault-tolerance argument **for the circuits and procedures actually studied in that paper**.

**What is implemented here:** exactly that **branching/measurement/recovery strategy**, but reusing our *independently developed*, optimized controlled-\(P\) flag gadget from [\`flagged_su2_measurement.py\`](flagged_su2_measurement.py), in which a syndrome ancilla is prepared in \(|+\rangle\), measured in X, and two flag CNOTs are inserted after data interactions 1 and 3. The conditional lookup is taken from our previous **exact simulator enumeration**. This is **not** a gate-for-gate transcription or independent reproduction of Chao–Reichardt Figure 1(c); no paper circuit is silently substituted by assumption.

**Correction to interpretation of PR #12.** The earlier experiment found **344 pairs of identical classical records + true final stabilizer syndromes with differing logical Pauli actions** when a **preexisting input error and a *new* circuit fault were both present**. That establishes limits on exact preservation of the original logical state for that *two-event* model and observation schedule. But a standard distance-three **one-fault** criterion requires preservation of the input logical information when the *total* input-error weight and the new fault count satisfy \(r+s\le1\)—**not** when \(r=s=1\). The cited literature has a separate output-containment condition when the EC unit itself suffers one fault. Consequently these 344 ambiguity groups **do not alone disprove standard one-fault fault tolerance**. They remain genuine stress-test data and should not be erased.

Standard two-part extended-rectangle EC criteria are also stated in the error-correction section of [*The small stellated dodecahedron code and friends* (2018)](https://pmc.ncbi.nlm.nih.gov/articles/PMC5990658/), section 5.1. This program explicitly verifies two important \(r+s\le1\) cases and reports \(r+s=2\) separately. It does **not** claim a fully general fault-tolerance theorem for arbitrary channels, arbitrary incoming states, noisy recovery gates or SU(2) hardware.

## 2. The tested schedule and its exact fault alphabet

Outer stabilizers, written in the effective two-dimensional color-singlet basis of each of the five modules:

\`\`\`text
XZZXI
IXZZX
XIXZZ
ZXIXZ
\`\`\`

Each *flagged* check contains four ideal encoded controlled-\(X\)/controlled-\(Z\) data operations plus two syndrome-ancilla-to-flag CNOTs: **six two-qubit gates**. Each *unflagged* check has four ideal encoded controlled-data operations. Fresh ancillas are prepared before each gadget and measured afterward. The simulator tracks **actual ancilla Pauli-sign flips** and *does not* inspect the hidden error's true final syndrome during correction.

At most one fault is injected in the initial flagged series for the **clean-input one-fault** theorem:

- 15 nonidentity two-qubit Pauli errors injected after *each* two-qubit gate,
- all 3 single-qubit Pauli preparation faults on each of the two ancillas,
- one classical syndrome and one classical flag readout flip.

This gives \(6\times15 + 6 + 2 = 98\) distinct Pauli *fault labels* per flagged stabilizer. Since a clean input never triggers before the intentionally injected fault, all four locations are reachable, giving \(4\times98=\mathbf{392}\) cases.

The circuit includes **24 abstract two-body gates on an untriggered run**, versus **48** for our earlier two-complete-flagged-round baseline. If a flag/syndrome triggers after the \(j\)-th flagged check (1-based), the entire execution contains \(6j+16\) abstract two-body gates, between 22 and 40. The 2018 paper's own detailed gate configuration and gate counts may differ from our separately designed gadget.

The follow-up unflagged measurements are actually executed gate-by-gate and their **recorded outcomes** determine correction. Under a clean initial codeword with **one** fault restricted to the initial flagged stage, a detected trigger means that the sole allowed fault has already happened; hence the follow-up is assumed fault-free by the global one-fault budget. If an input error is already present, a new fault can occur during follow-up, and that is explicitly included in the separate stress diagnostic rather than misrepresented as correctable by the clean-input theorem.

## 3. Exact test results — what passed

From the completed GitHub Actions run of \`paper_style_flag_ec.py --stress\`:

| Experiment | Cases |
|---|---:|
| Incoming weight-one Pauli on one module, with fault-free extraction | **15 / 15 fully corrected** |
| Clean input + one injected fault in initial flagged gadgets | **392** |
| Of those: fully corrected | **352** |
| Of those: left at most a weight-one-equivalent residual | **40** |
| Of those: nontrivial logical error | **0** |
| Of those: flagged/nontrivial-syndrome triggers | **328** |
| Of those: no trigger | **64** |

All **328 triggered cases** are fully corrected via the *actually measured* unflagged follow-up syndrome. The **40 residual cases** are *correctable one-module Pauli equivalence classes*, appropriate for a subsequent syndrome cycle; they are **not** claimed to be zero-error outputs.

The verifier calls the true-syndrome oracle **only after the protocol has returned**, to classify residual output. The runtime uses observed ancilla bits, observed flags, and a precomputed flag-specific syndrome-to-correction lookup. All Pauli propagation and comparisons are exact binary symplectic arithmetic.

These are counts of **possible fault labels**, not probabilistic failure rates or a physical noise threshold.

## 4. Stronger two-event stress diagnostic — not a one-fault release criterion

The \`--stress\` option adds **one incoming single-module Pauli** *and* **one separate new fault**. It enumerates only **physically reached gadget locations**, including unflagged follow-up gates that become reachable because the incoming error triggers the protocol. Consequently its case count is different from simply multiplying 15 by 392:

| Outcome relative to the **original input logical information** | Reached test cases |
|---|---:|
| Fully corrected | 1,216 |
| One correctable residual | 1,142 |
| Original logical state altered | 3,558 |
| Flag-specific decoder lookup undefined | 472 |
| **Total reached cases** | **6,388** |

These are **not** 6,388 required single-fault passes. In every row the **total event budget \(r+s=2\)** exceeds what a distance-three one-fault correction condition guarantees. The script prints explicit failure examples rather than suppressing them. In particular, some observed flag/syndrome combinations were absent from the original *single-fault* lookup, so the decoder is undefined under two simultaneous independent events.

The distinction is essential: a state with a different logical value after an **uncorrectable two-event pattern** is not proof that the circuit violates the *one-fault* definition. This is why the negative PR #12 results should not be marketed as a theorem that no proper five-qubit FT correction exists.

## 5. What is still missing

- Reproduce the **literal published Figure 1(c) gate ordering** and its error table, not just the branching strategy.
- Check the exact extended-rectangle output-containment condition and more general initial states; the present implementation exhausts a declared Pauli basis under restricted initial errors, not full noisy recovery.
- Include independent faults in the conditional follow-up at \(r+s=1\) if any other valid fault-free trigger mechanism exists; model noisy correction operations and one fault during ancilla reset.
- Implement fault-tolerant **controlled** logical-\(X\)/logical-\(Z\) operations on SU(2) matter/rishon blocks; our earlier projects only provide *uncontrolled* on-site gauge-invariant exchange matrices. There is still **no local physical hardware compilation**.
- Extend the simulator to within-module leakage and correlated faults induced by actual native gates, not only idealized encoded-module Pauli labels.

The new result is a verified **adaptive Clifford-level** improvement over PR #12's original measurement schedule. It is **not** a demonstration that Standard Model particles perform quantum recovery.

## 6. Reproduce

\`\`\`bash
python paper_style_flag_ec.py
python paper_style_flag_ec.py --stress
\`\`\`

Python 3.11 standard library plus the repository's previously merged \`flagged_su2_measurement.py\` and \`measured_multiround_su2.py\`. The GitHub Actions research workflow runs the stress version with the existing physics/code verifiers.
