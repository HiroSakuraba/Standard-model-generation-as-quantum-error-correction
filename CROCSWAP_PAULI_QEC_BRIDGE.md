# CrocSwap integer-mult-bounds: a genuine six-qubit code hiding inside a Pauli measurement partition

**18 October 2026 — independent research extension of a pinned third-party mathematical fixture.** This is a **real abstract stabilizer quantum error-correcting code**, rather than an analogy between integer multiplication and quantum errors. It does **not** prove a new quantum code family, lower the cost of the existing five-module construction, or show that Standard Model anomaly cancellation generates this code.

The executable independent test is [\`crocswap_pauli_qec_bridge.py\`](crocswap_pauli_qec_bridge.py).

## The third-party source and what it actually proves

CrocSwap's large repository [\`integer-mult-bounds\`](https://github.com/CrocSwap/integer-mult-bounds), maintained by Douglas Colkitt and contributors, investigates **conditional bounds on deterministic integer multiplication** using bit and complex transformation networks. Its selected multiplication result retains upstream hypotheses. Neither the multiplication exponent nor the project's use of coding-theory vocabulary is evidence for natural physical quantum error correction.

A more directly relevant mathematical subproject appears at pinned commit **\`20ecf2891d80b5acf18d980b2164fb44ea7faa90\`**, under the [Atlas applications package](https://github.com/CrocSwap/integer-mult-bounds/tree/20ecf2891d80b5acf18d980b2164fb44ea7faa90/research/atlas-frontier-applications/applications/adapters):

- [\`pauli_recovery.py\`](https://github.com/CrocSwap/integer-mult-bounds/blob/20ecf2891d80b5acf18d980b2164fb44ea7faa90/research/atlas-frontier-applications/applications/adapters/pauli_recovery.py) reconstructs a missing **maximal commuting Pauli measurement class** from the other classes.
- [The six-qubit measurement fixture](https://github.com/CrocSwap/integer-mult-bounds/blob/20ecf2891d80b5acf18d980b2164fb44ea7faa90/research/atlas-frontier-applications/applications/adapters/six-qubit-recovery/FULL-MEASUREMENT.json) provides 65 commuting measurement classes on six qubits.
- [Its independent checker](https://github.com/CrocSwap/integer-mult-bounds/blob/20ecf2891d80b5acf18d980b2164fb44ea7faa90/research/atlas-frontier-applications/applications/adapters/verify_pauli_recovery.py) validates the class partition, literal Pauli phases and mutually unbiased rank-one projectors. [The recorded receipt](https://github.com/CrocSwap/integer-mult-bounds/blob/20ecf2891d80b5acf18d980b2164fb44ea7faa90/research/atlas-frontier-applications/applications/adapters/INDEPENDENT-PAULI-RECEIPT.json) reports six qubits, 65 classes, 4,095 nonidentity Pauli observables, and 4,160 rank-one measurement projectors.

**Important distinction:** the source's “recovery” means *recovering an erased measurement class from a complete finite-geometric partition*. It does not implement a noisy physical quantum channel, a logical encoding of an unknown state, syndrome-based correction, or Knill–Laflamme recovery. Its full six-generator classes have 64 distinct joint sign outcomes, each rank one in a 64-dimensional Hilbert space; *each individual joint eigenspace encodes zero qubits*.

The following code construction and stabilizer-distance census are **our derived analysis of the pinned source data**, not results claimed by CrocSwap.

## 1. From source Pauli classes to explicit codes

Ignoring physically irrelevant center phases, represent a six-qubit Pauli operator by binary \((x,z)\in\mathbb F_2^{12}\). It commutes with another \((x',z')\) exactly when the *binary symplectic product*

\[
[x,z;x',z']=x\cdot z'+z\cdot x'=0\quad(\bmod\,2).
\]

Each CrocSwap measurement class is a **maximal isotropic subspace** \(L\) of dimension six under this product: all its elements commute. The source produces **65 disjoint maximal classes** whose union (apart from the shared identity) exhausts all 4,095 nonidentity Pauli labels.

A maximal class's six independent Hermitian +1 checks leave a **one-dimensional** code space, \([[6,0,\dots]]\). The key transformation is to choose a **five-dimensional hyperplane** \(S\subset L\), rather than using all six checks. The stabilizer projector is

\[
P_S=\frac{1}{2^5}\prod_{j=1}^{5}(I+\widehat g_j).
\]

Because the checks commute and are independent, \(\operatorname{rank}P_S=2^{6-5}=2\): it encodes **one logical qubit**. Every Pauli label can be screened for whether it (a) anticommutes with any check, (b) belongs to the stabilizer, or (c) commutes with the stabilizer but acts nontrivially on the code. Its distance is therefore

\[
\boxed{d(S)=\min\{\mathrm{wt}(E):
E\in S^\perp\setminus S\}.}
\]

Unlike simply discarding an arbitrarily listed generator, one must enumerate *all 63 hyperplanes* of a six-dimensional binary vector space. For a nonzero six-bit coefficient functional \(h\), a basis of \(\ker(h)\) can be constructed without a numerical rank heuristic by choosing the first pivot \(p\) with \(h_p=1\) and taking \(g_j+h_j g_p\) for \(j\ne p\).

The independent script reconstructs the source six-qubit fixture from its primitive binary-field polynomial \(x^6+x+1\), field trace pairing, and pinned symplectic transvections; checks the first, erased 23rd, and last maximal classes against the source; and verifies the **entire 65-class, 4,095-label disjoint partition** before looking for codes.

## 2. Complete 4,095-subspace census (exact)

There are 65 maximal commuting classes and 63 possible hyperplanes per class: **4,095 distinct candidate \([[6,1,d]]\) stabilizer subspaces** in this specified fixture.

| Distance | Candidate codes | Capability |
|---|---:|---|
| \(d=1\) | **1,539** | Does not detect every single-qubit Pauli error |
| \(d=2\) | **2,539** | Detects all weight-one Pauli errors, but does not correct an arbitrary unknown single-qubit error |
| \(d=3\) | **17** | Corrects arbitrary errors confined to one qubit |
| **Total** | **4,095** | Exactly exhaustive for the specified finite source fixture |

The classification uses **integer-only binary linear algebra**, no floating point or statistical sampling. It scans all 4,095 nonidentity six-qubit Pauli labels for each candidate, in increasing Pauli weight. The 17 distance-three subspaces are distinct subspaces of *this* symplectic partition. **They are not asserted to be 17 inequivalent quantum-code families or newly discovered stabilizer codes.**

The straightforward act of deleting one of the six supplied generators performs badly: for all \(65\times6=390\) five-generator subspaces obtained this way, the distance is **one**. This makes the distinction between a source measurement basis and an error-correcting stabilizer code operationally important. The correctly chosen *linear combinations* of checks are the entire improvement.

## 3. One explicit, degenerate \([[6,1,3]]\) certificate

Take the **first** source measurement class and the coefficient hyperplane \(h=51\). The resulting five independent commuting Hermitian Pauli checks are:

| Integer symplectic label | Pauli word, qubits 0–5 |
|---:|---|
| 3094 | \`IXXIYZ\` |
| 584 | \`ZIIYII\` |
| 283 | \`XXZXXI\` |
| 2265 | \`YZIXXZ\` |
| 2164 | \`ZIXIXY\` |

A Pauli label \(v\) is the usual **Hermitian** operator \(i^{\mathrm{wt}(x\wedge z)}X^xZ^z\), so a \`Y\` means the Hermitian matrix \(Y=iXZ\), not the real operator \(XZ\) without its phase. The five independent generators commute, leaving a rank-two code projector and **one encoded qubit**.

An explicit anticommuting logical operator pair is

\[
\overline X=\texttt{ZXXIII},\qquad
\overline Z=\texttt{IZXIXI}.
\]

Each has physical weight three, commutes with every stabilizer, and anticommutes with the other. Exhausting all lower-weight Paulis proves that **no weight-one or weight-two operator acts as a nontrivial logical operator**, so the distance is *exactly three*.

For the full 19-element error basis

\[
\mathcal E=\{I\}\cup
\{X_i,Y_i,Z_i:i=0,\ldots,5\},
\]

the script checks **all 361 ordered pairs** \(E,F\). In every case, either \(E^\dagger F\) anticommutes with at least one stabilizer, giving

\[
P_SE^\dagger FP_S=0,
\]

or it equals a stabilizer *up to a Pauli phase*, giving

\[
P_SE^\dagger FP_S=c_{EF}P_S.
\]

These are exactly the **Knill–Laflamme quantum error-correction conditions** for an arbitrary one-qubit quantum channel. Syndrome-based recovery is explicit: measure the five stabilizers, choose one representative Pauli for each observed syndrome, and apply its inverse. The source model provides abstract Pauli observables, **not** a physical nondemolition measurement circuit.

This code is **degenerate**: two distinct single-qubit errors, \`ZIIIII\` and \`IIIYII\`, have the same syndrome because their product is the weight-two stabilizer \`ZIIYII\`. A correct recovery must allow that collision; a false requirement that every error have a different syndrome would incorrectly reject the code. The 19 errors occupy **18** syndrome values.

## 4. Relation to the Standard Model project

The new six-qubit code can be *formally placed* on **six independent copies of our two-dimensional SU(2) matter/rishon singlet module**, interpreting each \(X_i,Y_i,Z_i\) as operators acting on the \(i\)th gauge-invariant inner qubit. Our previous exact on-site spin-exchange constructions give a possible abstract representation of each logical \(X,Z\), and \(Y\) can be synthesized algebraically from them.

Under the **same idealized assumptions** as the earlier five-module recovery construction, the six-module code corrects an error confined to one encoded singlet module. This is a legitimate mathematical stabilizer concatenation, not a natural physical encoding supplied by anomaly cancellation.

The comparison is:

| Property | Previously tested five-module code | This source-derived six-module code |
|---|---|---|
| Outer code | Published \([[5,1,3]]\) perfect code | Explicit source-derived \([[6,1,3]]\) degenerate code |
| SU(2) singlet inner modules needed | **5** | **6** |
| Independent outer stabilizer checks | **4** | **5** |
| General single-module Pauli recovery | Yes (exactly checked) | **Yes (361 exact KL pair conditions)** |
| Measurement hardware / extraction fidelity | Still idealized | **Still idealized** |
| Direct link to family anomaly cancellation | **None established** | **None established** |
| Claim of better overhead | No | **No; this uses more modules/checks** |

The extra syndrome capacity (five check bits, 32 values, only 18 occupied by the identity and first-order error classes) *may* prove useful in designing flagged extraction to distinguish propagated correlated faults, but that is a research **hypothesis**, not established by this code. Some checks have substantial weight, and the necessary nonlocal SU(2)-compatible measurements are unimplemented. No new quantum memory lifetime, physical threshold or device advantage is inferred.

## 5. Other items inspected in CrocSwap and their real scope

Three more constructions are worth remembering, but they are **not themselves operational quantum codes**:

- [Fano network-coding audit](https://github.com/CrocSwap/integer-mult-bounds/blob/20ecf2891d80b5acf18d980b2164fb44ea7faa90/docs/research/fano-completion-audit.md): a genuine **classical** characteristic-two network code transmits three source bits by XOR through a constrained graph. The source proves a coding-versus-routing property for that graph and shows limitations of certain reversible completions. It is not a stabilizer recovery gadget. Fano geometry also underlies familiar classical Hamming and quantum Steane-code constructions, but that is **outside the source's stated result**.
- [Exact Clifford–Pauli anchor comparison](https://github.com/CrocSwap/integer-mult-bounds/blob/20ecf2891d80b5acf18d980b2164fb44ea7faa90/research/clifford-pauli-anchor/PROOF.md): for validated Clifford gate words, **signed conjugation of all Pauli generators plus one exact reference column** determines equality of the whole unitary (including global phase). The repo implements a substantially faster circuit verifier. This is directly relevant to **checking our five/six-module Clifford extraction circuits** without simulating every input vector, but it certifies circuit equality, **not error recovery**.
- [Coded-carry audit](https://github.com/CrocSwap/integer-mult-bounds/blob/20ecf2891d80b5acf18d980b2164fb44ea7faa90/docs/research/coded-carry-audit.md): branches integer carries into masked affine translations and restores overwritten gap bits exactly. Here “coding” means reversible computational restructuring, not redundancy capable of detecting/correcting noise. Its authors report a scoped cost obstruction rather than a new multiplication improvement.
- [Binary/rational rank-product core](https://github.com/CrocSwap/integer-mult-bounds/blob/20ecf2891d80b5acf18d980b2164fb44ea7faa90/docs/research/rank-product-core.md): low-rank incidence over \(\mathbb F_2\) combined with separate rational fitting constraints helps compress multiplication networks. Those are useful finite-algebraic conditions, but **different fields cannot silently share one rank argument**, and no physical noise model or quantum recovery follows from that result.

All these projects emphasize a useful methodological lesson: **an algebraic identity, a reversible computation, an erased measurement reconstruction, and a quantum recovery code are four different mathematical claims**. The new six-qubit construction is firmly in the last category at the *abstract stabilizer code* level.

## 6. Reproduction, source binding and next decision

\`\`\`bash
python crocswap_pauli_qec_bridge.py
\`\`\`

Python standard library only. The script has no runtime network dependency, imports no code from CrocSwap, and reproduces its pinned binary-field spread algorithm independently with explicit expected source generators, class cardinalities and a full Pauli partition check. It searches **all 4,095 candidate subspaces**, verifies the selected code and checks every single-qubit Pauli error pair.

**Provenance:** CrocSwap's source, fixture, finite-field construction and measured-class erasure recovery belong to the credited Atlas contributors and are Apache-2.0. The hyperplane-distance search, example stabilizer extraction, and explicit KL verification here are a separate derived analysis. This is not a novelty or authorship claim over earlier stabilizer-code literature.

**Next milestone:** Test the six-module code in our existing **flagged stabilizer-measurement simulator**: synthesize its five check measurements from physically allowed SU(2) singlet operations, enumerate a single fault in a check ancilla and propagated errors on multiple modules, and compare the **actual recovery behavior and hardware overhead** with the established five-module code under the **same noise assumptions**. Do not conclude that additional syndromes help until that test succeeds.

Source-independent scientific grounding: the stabilizer formalism and Knill–Laflamme criterion are standard, e.g. [Knill and Laflamme, *Theory of quantum error-correcting codes*, PRA 55, 900 (1997)](https://doi.org/10.1103/PhysRevA.55.900), and [Gottesman, *Stabilizer Codes and Quantum Error Correction*, 1997](https://arxiv.org/abs/quant-ph/9705052).
