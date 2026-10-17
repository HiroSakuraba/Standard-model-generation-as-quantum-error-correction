# Fault propagation in outer stabilizer measurements on SU(2) singlet modules

**Status — 17 October 2026.** Exact single-fault Clifford *gadget* enumeration and a limited, explicitly conditional correction lookup. Not a full fault-tolerant recovery protocol. Not a physical SU(2) hardware implementation. Use [\`flagged_su2_measurement.py\`](flagged_su2_measurement.py) to reproduce everything below.

## 1. Why this step was needed

The [five-module singlet construction](FIVE_BLOCK_SU2_RECOVERY.md) proved that the established \([[5,1,3]]_2\) code, applied to five **independent SU(2)-gauge-singlet inner qubits**, corrects any arbitrary **single-module** error with perfect syndrome measurements.

That algebraic theorem says nothing about faults occurring *during* the four multi-module stabilizer measurements. A faulty syndrome ancilla can propagate its error into **two or more modules**, exceeding the code's distance-three single-error capability. To evaluate that problem instead of assuming it away, this step analyzes explicit **controlled-Pauli extraction gadgets** for the code's four outer generators:

\`\`\`text
XZZXI    IXZZX    XIXZZ    ZXIXZ
\`\`\`

The operators \(X,Z\) here refer to **entire gauge-invariant singlet modules**, as established from actual on-site SU(2)-invariant exchanges in \`five_block_su2_recovery.py\`. The ancilla couplings to those operators are *ideal abstract Clifford gates*. A physically compiled controlled-exchange implementation remains to be designed.

## 2. Circuits and ancilla preparation

The syndrome qubit \(m\) starts in \(|+\rangle\) and is measured in the **X basis**. For each of the four nonidentity positions of a stabilizer \(g=P_1P_2P_3P_4\), perform an ideal controlled-\(P_j\) with \(m\) as the computational-basis control and data module \(j\) as target. \(P_j=X\) is represented by CNOT; \(P_j=Z\) by controlled-\(Z\) (CZ).

Because different factors act on different modules, their controlled operations commute. In the *fault-free* circuit, the ancilla-X measurement measures \(g\) nondestructively.

**Bare circuit:**

\`\`\`text
m: |+> ── C(P1) ── C(P2) ── C(P3) ── C(P4) ── X measurement
\`\`\`

**Two-coupling flagged circuit:** prepare flag qubit \(f\) in \(|0\rangle\), insert a CNOT \(\mathrm{CX}_{m\to f}\) after the **first** and after the **third** data coupling, and measure \(f\) in the Z basis:

\`\`\`text
m: |+> ── C(P1) ──●── C(P2) ── C(P3) ──●── C(P4) ── X measure
                   │                      │
f: |0> ─────────────X──────────────────────X─────────────── Z measure
\`\`\`

Without faults the two flag CNOTs cancel *coherently*, and the flag measurement is \(0\), preserving the intended stabilizer measurement.

For comparison, an **early-flag placement** puts the first flag CNOT *before* \(C(P_1)\) and the second after \(C(P_3)\). This also catches dangerous hook errors, but does not in general provide enough information to select a unique correction.

## 3. Exact enumerated fault model

Each gate acts on precisely two **abstract qubits**—the syndrome/flag ancillas or an encoded SU(2) data module. Fault injection includes:

- **Every nonidentity two-qubit Pauli after every two-qubit gate** (15 choices: \(\{I,X,Y,Z\}^{\otimes2}\setminus\{II\}\)).
- All \(X,Y,Z\) preparation Pauli errors on each present ancilla.
- A single independent classical readout bit flip for the X-basis syndrome readout; and for the Z-basis flag readout if present.

At most **one listed fault occurs per execution**. This is exhaustive for this *abstract circuit* fault set. All conjugation of a Pauli error through subsequent CNOT/CZ gates is tracked exactly in the binary symplectic representation; Pauli global phases are ignored because they cannot change syndromes or Knill–Laflamme classes.

There are **64 fault cases** for the four-gate bare gadget (60 gate faults + three syndrome-preparation faults + one readout fault), and **98 cases** for the flagged gadget (90 gate faults + six preparation faults + two readout faults).

The measured syndrome from the potentially faulty gadget is **not** trusted. For the purpose of *classifying* errors, the script additionally assumes a **subsequent fault-free measurement of all four stabilizers**, yielding the true four-bit data syndrome. This trusted later syndrome is *an oracle assumption*, not a circuit implementation or a result of this test.

## 4. Results for every outer-code stabilizer

| Stabilizer measured | Bare: faults leading to logical failure under naive one-site decoding | Poor flag: ambiguous recovery groups | Improved flag: dangerous cases that raise flag | Improved flag: ambiguous recovery groups |
|---|---:|---:|---:|---:|
| \`XZZXI\` | 12 | 1 | **12/12** | **0** |
| \`IXZZX\` | 12 | 1 | **12/12** | **0** |
| \`XIXZZ\` | 12 | 1 | **12/12** | **0** |
| \`ZXIXZ\` | 12 | 1 | **12/12** | **0** |

These are **numbers of enumerated fault locations and Pauli labels**, **not statistical failure rates**. To make them into rates one must specify physical gate, preparation, measurement and idling noise rates and their correlations.

In the improved gadget, **no fault with a trivial flag produces a residual logical error under the ordinary one-module decoder**. For flagged faults, the pair

\[
(\hbox{flag outcome},\ \hbox{trusted subsequent four-bit syndrome})
\]

identifies a unique data-error class up to a code stabilizer among all enumerated single faults. The script computes a recovery representative for each occupied pair of classical labels and verifies that its product with *every* corresponding propagated data error has (i) zero stabilizer syndrome and (ii) zero logical \(X\) and \(Z\) action. There are **21** observed (flag, syndrome) correction keys per measured stabilizer.

By contrast, the early-flag gadget has one group with the **same flag and the same subsequent true syndrome but different logical residuals**. Those errors cannot both be corrected by a single recovery choice. **Flag detection alone is not enough: placement and correction information matter.**

The original bare circuit also has ambiguous syndrome groups and produces 12 unflagged, naively miscorrrected single-fault cases per check.

### What is—and is not—established

**Established conditional circuit property:** Given an initially fault-free five-module codeword, no more than one of the listed abstract gate/ancilla faults **inside one extraction gadget**, and a later perfectly measured full data syndrome, the optimized flag placement provides enough information to choose a data correction up to stabilizer. This is a reproducible one-gadget, single-fault result.

**Not established:** A complete **one-fault-tolerant error-correction round**, which must handle at least (a) preexisting data errors, (b) faults in the other stabilizer gadgets, (c) errors or flags in the *repeated* syndrome extraction, (d) any flag/postselection scheduling and decoder validation, and (e) correlated errors or leakage introduced while implementing the required ancilla-to-module controlled operators.

The model tracks single Pauli faults at the **encoded-module** level, not arbitrary physical faults within the underlying six-spin matter/rishon Hilbert space. A physically elementary gate fault might spread across an entire module or multiple modules. Moreover, local on-site spin exchanges implement **logical X/Z generators in the gauge-singlet block**, but this does not automatically supply fault-tolerant controlled versions of these logical rotations or ancilla readout in a gauge-theory simulator.

## 5. Literature and relation to known flag correction

The purpose and structure of the extra flag qubit follow the established **flag fault-tolerance** paradigm. In particular:

- Rui Chao and Ben W. Reichardt, [*Quantum Error Correction with Only Two Extra Qubits*](https://doi.org/10.1103/PhysRevLett.121.050502), *Phys. Rev. Lett.* **121**, 050502 (2018), developed two-ancilla fault-tolerant procedures explicitly for the \([[5,1,3]]\) code. This repository does **not** claim to reproduce or improve their full experimental protocol.
- Christopher Chamberland and Michael E. Beverland, [*Flag fault-tolerant error correction with arbitrary distance codes*](https://doi.org/10.22331/q-2018-02-08-53), *Quantum* **2**, 53 (2018), gives the broader framework for identifying dangerous fault propagation through flag circuits.
- Chao and Reichardt, [*Flag Fault-Tolerant Error Correction for any Stabilizer Code*](https://doi.org/10.1103/PRXQuantum.1.010302), *PRX Quantum* **1**, 010302 (2020), discusses the difference between a fault-detecting **flag gadget** and an entire fault-tolerant syndrome-extraction protocol.

Our contribution is a **reproducible exact comparison of specific ideal Clifford gadgets** in the context of the SU(2) inner modules previously constructed in this repository. The outer five-qubit code, the flag paradigm and its general theory are all established prior work.

## 6. Next milestone

The next release criterion should be an **end-to-end single-fault correction round**, with every syndrome gadget and repeat included in the same circuit, a decoder using only actually measured bits (no perfect-syndrome oracle), and exhaustive one-fault propagation through the *whole round*. That is a larger experiment than the current conditional classification.

At the physical layer, a separate task is to replace ideal encoded controlled-Pauli interactions with **explicit SU(2)-gauge-invariant, locality-respecting native gates**, then enumerate faults in those elementary gates and test for multi-module propagation. This would address the central missing operational ingredient before any claim about real particles.

## Reproduce

\`\`\`bash
python flagged_su2_measurement.py
\`\`\`

Python 3.11+ standard library only. The existing GitHub Actions research workflow runs it alongside the full charge-code and SU(2) test collection.
