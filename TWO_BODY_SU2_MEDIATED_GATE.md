# Two-body SU(2) exchange mediation: a physical gate-set test and its leakage cost

**Status (18 October 2026): exact representation-theory and second-order perturbation results, plus full finite eight-spin numerical evolution.** This is a deliberately engineered and truncated non-Abelian matter/rishon **model**. It neither reproduces a specific device Hamiltonian nor produces an exact native controlled gate, nor does it connect the original Standard Model anomaly-charge code to physical particles.

Executable verification: [\`two_body_su2_mediator.py\`](two_body_su2_mediator.py); the previous ideal ancilla-conditioned Hamiltonians are in [\`SU2_CONTROLLED_EXCHANGE_IMPLEMENTATION.md\`](SU2_CONTROLLED_EXCHANGE_IMPLEMENTATION.md).

## 1. A rigorous obstruction: why an ordinary neutral ancilla is insufficient

Take a **color-neutral ancilla** \(\mathcal H_a=\mathbb C^2\) and a **single SU(2) fundamental** color spin \(\mathcal H_c=\mathbb C^2\). Local SU(2) acts as \(I_a\otimes U_c\). For any two-body Hermitian coupling \(H\) to be invariant under *every* gauge rotation,

\[
[H,I_a\otimes U_c]=0\qquad \forall U_c\in SU(2).
\]

By irreducibility of the fundamental representation and Schur's lemma, the entire invariant operator algebra is

\[
\operatorname{End}(\mathcal H_a)\otimes I_c.
\]

It contains **no entangling interaction** between the neutral ancilla and one color fundamental. This is a *precise restriction on this specified two-particle, two-body gate set*, not a universal ban on gauge-invariant controlled operations. Gauge-scalar couplings involving **two** colored spins, or additional mediators/higher-body interactions, remain possible.

The verifier independently writes all 16 basis matrices on \(\mathbb C^2\otimes\mathbb C^2\), imposes the three SU(2) commutator equations, and confirms **rank 12** and **invariant dimension 4**, exactly.

This explains why the explicit neutral-ancilla Hamiltonian \(|1\rangle\langle1|_a\otimes \mathbf S_i\cdot\mathbf S_j\) proposed in PR #14 was inherently *three-body*: the neutral ancilla needs to control a **color-scalar operator involving at least two** fundamental spins.

## 2. An alternate elementary interaction set using an encoded colored ancilla

Instead of an external neutral physical qubit, add **four more SU(2) fundamental spins** at the same gauge vertex as the previous four-doublet rich matter/rishon vertex. Their total-spin-zero subspace is two-dimensional, and so encodes one **collectively neutral** ancilla logical qubit.

- Ancilla register \(A\): four co-located spin-\(1/2\) fundamentals, \((A_0,A_1,A_2,A_3)\).
- Data register \(D\): the existing four co-located matter/rishon doublets, \((D_0,D_1,D_2,D_3)\).
- Physical local Hilbert space: **eight SU(2) spins = \(2^8=256\) dimensions**, acted on by one local SU(2) color transformation.
- Each four-spin singlet register has logical dimension two; the desired code block \(P\) therefore has dimension four.

A simple *native interaction proposal* is only **pairwise Heisenberg exchange**,

\[
h_{ij}=\mathbf S_i\cdot\mathbf S_j
=\tfrac12(\mathrm{SWAP}_{ij}-\tfrac12 I),
\]

plus a block-spin Casimir energy penalty, also expressible entirely using two-body exchange.

The concrete Hamiltonian is

\[
\boxed{
H=\Delta(\mathbf S_A^2+\mathbf S_D^2)
  +J_1\,\mathbf S_{A_0}\cdot\mathbf S_{D_0}
  +J_2\,\mathbf S_{A_1}\cdot\mathbf S_{D_1}.
}
\]

The \(\Delta\) term includes six intra-block exchange pairs for **each** four-spin register, and the two terms proportional to \(J_1,J_2\) link one individual spin from each register. All constituent terms are **SU(2)-invariant two-spin operators at the same gauge vertex**, with no direct three-body interaction assumed. The local penalty is **additional engineered physics**, not simply the original Gauss-law Casimir: all states under these interactions can obey the *combined* vertex Gauss law even when individual four-spin registers leave their singlet spaces.

**Physical proviso:** the extra four color fundamentals are additional gauge-charged constituents. They must be supplied by an actual matter/rishon simulator and combined into an SU(2) singlet. Calling the encoded logical qubit a "neutral ancilla" does not mean its constituent spins transform trivially. Its preparation, placement and independent control have not been demonstrated, and this construction changes the original model's particle inventory.

## 3. Exact second-order entangling Hamiltonian

Let \(P\) project onto the tensor product of the two four-spin singlet subspaces. The two cross exchanges together are \(V=J_1h_{A_0D_0}+J_2h_{A_1D_1}\).

Because a single spin operator takes a spin-zero state into spin-one, **\(PVP=0\)**. Each cross exchange takes both sub-registers from total spin zero to total spin one; the penalty costs \(2\Delta\) in each register, or \(4\Delta\) total.

The exact leading, nonzero Schrieffer–Wolff effective Hamiltonian is therefore

\[
H_{\mathrm{eff}}^{(2)}
=-\frac1{4\Delta}PV^2P.
\]

The verifier constructs the complete four-dimensional logical basis as explicit **eight-spin valence-bond wavefunctions** and computes the exact matrix elements of these operators, with no numerical tolerance.

In each four-spin logical basis the internal pair exchange obeys

\[
P\mathbf S_0\cdot\mathbf S_1P
=-\tfrac14 I-\tfrac12 Z.
\]

Consequently, for arbitrary real couplings, the exact *second-order coefficient* of the entangling term is

\[
\boxed{J_{ZZ}=-\frac{J_1J_2}{24\Delta}.}
\]

For equal \(J_1=J_2=J\),

\[
\boxed{
H_{\mathrm{eff}}^{(2)}
=-\frac{J^2}{4\Delta}\left[
 \frac5{12}I+
 \frac1{12}(Z_A+Z_D)+
 \frac16 Z_AZ_D\right].
}
\]

This entangles the two singlet qubits in the **effective model**. The identity gives a global phase, and the single-\(Z\) components give local phases. If these local phases can be controlled or compensated, evolving under the \(ZZ\) term through magnitude \(\pi/4\) gives a gate locally equivalent to CZ; an additional logical Hadamard produces CNOT.

The **nominal duration** is

\[
t_{\pi/4}=\frac{6\pi\Delta}{|J_1J_2|}
\quad(\hbar=1).
\]

This is an **approximate perturbative gate**, not a finite-pulse equality. Higher-order terms and leakage invalidate a claim of exact CZ.

## 4. Full eight-spin evolution: quantified gate error and leakage

To determine whether the second-order model is predictive in a non-asymptotic regime, the script constructs and diagonalizes the **actual 256×256 Hamiltonian** for \(\Delta=1\), \(J_1=J_2=J\), with all 12 intra-block and both cross-block Heisenberg couplings.

It computes the exact (to numerical eigensolver precision) finite-time evolution \(e^{-iHt_{\pi/4}}\), compresses that unitary to the four-dimensional singlet-pair code space, and compares it against \(e^{-iH_{\mathrm{eff}}^{(2)}t_{\pi/4}}\).

| \(J/\Delta\) | \(t_{\pi/4}\Delta\) | Operator-norm difference from second-order prediction | Maximum leakage probability from the logical space |
|---:|---:|---:|---:|
| 0.4 | 117.810 | 0.362404 | 0.019131 |
| 0.2 | 471.239 | 0.176107 | 0.002408 |
| 0.1 | 1884.956 | 0.088395 | 0.000533 |

The operator-norm comparison is between **raw projected evolution matrices**, with no calibration of global or logical local phases. It is **not** a reported physical gate infidelity. The leakage quantity is \(1-\lambda_{\min}(P U^\dagger P U P)\) on the four-dimensional code, the **worst-case probability to leave that code** at the end of the nominal pulse. It is not a stochastic Pauli error rate or complete leakage-correction benchmark.

**Interpretation:** Lowering \(J/\Delta\) improves the perturbative approximation and suppresses finite-pulse leakage, but makes the gate slow: the desired phase accumulates only at order \(J^2/\Delta\). None of the tested finite pulses is shown to be exactly equivalent to a CZ after any specified calibration. The effect of hardware dephasing, pulse noise and finite ancilla lifetime is not modeled. We should not claim a useful error-correction clock rate or threshold from these results.

## 5. Gauge-invariant faults and physical error propagation

The simulator checks both cross exchanges commute exactly with the **full eight-spin Gauss-law generators** on all \(2^8\) basis vectors (2,048 exact commutator comparisons). Every intra-register Heisenberg exchange is likewise a color scalar.

The full finite-time unitary **remains inside the total SU(2)-singlet sector** when initialized there, but generally leaks **out of the chosen tensor product of two four-spin singlet registers**. This is an especially important distinction: leakage from the *logical code* is not necessarily a violation of the *gauge constraint*. Monitoring only the total Gauss charge would miss it.

Two different physical error classes remain to be modeled:

1. **Gauge-invariant coupling faults:** \(\delta J_1,\delta J_2,\delta\Delta\), timing errors and unintended pair exchanges. These change the entangling strength and can drive leakage **without producing Gauss violations**. They do not automatically fit a single logical Pauli fault.
2. **Gauge-violating imperfections:** noise acting directly on one colored constituent, which a gauge-simulator implementation may suffer and which can violate local Gauss law.

The script's full-pulse deviation and leakage form a **baseline for native-operation-level error modeling**, but it does not yet include stochastic faults at each individual elementary coupling and does not establish the flagged extraction circuit's physical fault tolerance.

## 6. Scientific value and research decision

The present result closes a narrowly defined gap: **two-body, SU(2)-invariant interactions can mediate an entangling gate between two gauge-invariant encoded singlet qubits**, albeit perturbatively and with leakage. This demonstrates a pathway around the strict neutral-ancilla/single-color-spin two-body obstruction.

What it does **not** establish is a native hardware gate set available in a realistic matter/rishon device, an exact Clifford gate under finite pulses, low-overhead cross-module stabilizer measurement, physical recovery under gate faults, or anything like a naturally occurring Standard Model error-correcting code.

The next two experiments are now unambiguous:

- **Pulse synthesis:** optimize controllable intra-block exchanges and the two cross-block couplings for a near-ideal controlled-phase gate *and explicitly penalize both leakage and long duration*. Then test that pulse sequence against coherent \(\delta J\) and stochastic SU(2)-invariant coupling fluctuations.
- **System-level circuit:** substitute these physical pulse-derived gates into the flagged five-module measurement circuit, and test whether one faulty pulse can spread errors across both the ancilla block and more than one data module. Model the larger encoded ancilla's own preparation and readout, not as ideal two-level operations.

Do not proceed to Standard Model anomaly-code claims without a compatible **physical Hilbert-space encoding, conserved/measurable checks and demonstrated recovery**.

## Reproduce

\`\`\`bash
python -m pip install numpy sympy==1.14.0
OPENBLAS_NUM_THREADS=1 python two_body_su2_mediator.py
\`\`\`

The accompanying GitHub Actions step runs this deterministic exact-plus-diagonalization check; the notebook/code never allocates any \(1024^5\) state for the previous five-module system.

## Prior literature and novelty boundary

This construction is an illustrative combination of established four-spin singlet encodings, exchange-only interactions and second-order effective Hamiltonian methods. It is **not** claimed to be a novel universal quantum gate family.

- Bacon, Kempe, Lidar and Whaley, [*Universal Fault-Tolerant Quantum Computation on Decoherence-Free Subspaces*](https://doi.org/10.1103/PhysRevLett.85.1758), *Phys. Rev. Lett.* **85**, 1758 (2000).
- Weinstein and Hellberg, [*Energetic suppression of decoherence in exchange-only quantum computation*](https://doi.org/10.1103/PhysRevA.72.022319), *Phys. Rev. A* **72**, 022319 (2005), already studied encoded qubits, gap protection and exchange-only controlled phase gates.
- [*Singlet-Only Always-On Gapless Exchange Qubits with Baseband Control*](https://doi.org/10.1103/j1hn-k5m8), *Phys. Rev. Lett.* (2025), already proposes four-electron spin-singlet qubits and two-qubit exchange gates. The present construction is a model-specific exact audit and numerical baseline, not a priority claim over this literature.
