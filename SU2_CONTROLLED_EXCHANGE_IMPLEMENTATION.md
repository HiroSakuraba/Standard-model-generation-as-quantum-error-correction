# Gauge-invariant ancilla-controlled exchange gates for the SU(2) singlet modules

**Status, 18 October 2026.** Exact constructive Hamiltonians, microscopic gauge-commutator checks and encoded-Clifford verification, plus an independent primary-source audit of the original five-qubit flag error catalog. **This is not a native gate compilation, not a complete fault-tolerant circuit, and not evidence for quantum error correction by the Standard Model.**

Executable proofs: [\`su2_controlled_exchange.py\`](su2_controlled_exchange.py) and [\`published_flag_figure2_catalog.py\`](published_flag_figure2_catalog.py).

## 1. The missing controlled operation

Previously, each physical inner module was a pair of color-singlet states in a **four-link, two-matter-site, single-rishon-per-link SU(2) truncation**. Five independent modules with the published five-qubit code provide exact register-level correction against one-module noise. The outer flagged measurement circuit nevertheless assumed access to ideal ancilla-controlled \(X_{\rm in}\) and \(Z_{\rm in}\) gates. Proving that *uncontrolled* local exchanges represent Pauli operators did not establish their *controlled* counterparts.

Here the ancilla is an **external SU(2)-color-neutral two-level system**, assumed to be brought to the four-doublet vertex. We construct explicit **gauge-invariant Hamiltonian operators** whose evolution produces the required controlled logical gates **on the singlet subspace**. This proves *kinematic reachability if the specified ancilla-conditioned interactions can be synthesized*. It does **not** prove such a coupling is a native two-body Hamiltonian of a particular gauge simulator.

## 2. Three physical color-scalar exchange operators

Retain the reference endpoint pattern \(o=(0,1,1,1)\) from [\`FOUR_LINK_SU2_DFS.md\`](FOUR_LINK_SU2_DFS.md). Matter doublets \(m_1\) and rishons \(r_1,r_2,r_3\), labeled as spins \(1,3,4,5\), are co-located at vertex 1. The normalized states \(|0_L\rangle, |1_L\rangle\) are the two inequivalent SU(2) total-spin-zero intertwiners.

Let \(S_{1j}\) swap the color spin of \(m_1\) and rishon spin \(j\). For co-located spins, it commutes with all local SU(2) Gauss-law generators because exchange is invariant under simultaneous rotations. Within the singlet pair,

\[
\begin{aligned}
P S_{13}P &= - Z_L,\\
P S_{14}P &= (Z_L-\sqrt3 X_L)/2,\\
P S_{15}P &= (Z_L+\sqrt3 X_L)/2.
\end{aligned}
\]

Consequently the **physical Hermitian** combinations

\[
\boxed{\widetilde Z=-S_{13},\qquad
\widetilde X=\frac{S_{15}-S_{14}}{\sqrt3}}
\]

act as **exact logical Pauli \(Z_L\) and \(X_L\)**. The script verifies these as **equalities of full six-spin vectors** acting on *both* logical basis states, not merely equality of \(2\times2\) projections. This ensures the exchanges do not leak a codeword out of the chosen two-dimensional singlet block under the ideal operation.

Using the spin-half identity \(S_{ij}=2\mathbf S_i\!\cdot\!\mathbf S_j+\frac12I\),

\[
\widetilde X=
\frac{2}{\sqrt3}\left(
\mathbf S_1\!\cdot\!\mathbf S_5-
\mathbf S_1\!\cdot\!\mathbf S_4\right),\qquad
\widetilde Z= -2\mathbf S_1\!\cdot\!\mathbf S_3-\tfrac12 I.
\]

**Crucial:** No *single* co-located permutation \(S_{13}\), \(S_{14}\), or \(S_{15}\) equals logical \(X_L\) up to sign in this basis. \(X_L\) requires a **balanced difference of two exchanges**, or a separate pulse-synthesis scheme. The difference is Hermitian, but is not claimed to square to the identity outside the protected block.

## 3. Gauge invariance on *all* orientation sectors, not just selected states

A naked swap of spins at **different vertices** would be unphysical: it would not commute with independently acting SU(2) gauge transformations at those vertices. The model includes 16 endpoint orientations, so we must not use the reference-sector swaps blindly.

Define \(N_{j,R}\) as the gauge-invariant occupation projector recording whether the corresponding rishon of link \(j\) is at vertex 1, co-located with matter spin 1. Define the **endpoint-conditioned** physical spin exchange

\[
\mathsf T_{1j}=N_{j,R}S_{1j}+(I-N_{j,R})I,\qquad j=3,4,5.
\]

On a sector in which the two spins are co-located, this acts as the usual SU(2) exchange. On a sector where the rishon is at the opposite vertex, it acts as **identity**, never swapping colors across gauge sites. It conserves every rishon-end occupation.

This construction needs only **local endpoint occupancy at vertex 1**, rather than a projector onto the complete four-link pattern. The resulting direct-sum operator is Hermitian and unitary on the complete truncated kinematic Hilbert space.

The executable checks **all 16 orientations times 64 spin configurations = 1,024 kinematic basis states** and all three exchanges against both local generators \(2G_v^z,G_v^+\). This is **12,288 exact commutator checks**, all passed in CI. The lowering-generator and Hermiticity relations complete the Lie-algebra commutation argument. Linear combinations

\[
\widetilde X_{\rm all}=
(\mathsf T_{15}-\mathsf T_{14})/\sqrt3,\qquad
\widetilde Z_{\rm all}=-\mathsf T_{13}
\]

are therefore gauge invariant in *every* occupation sector. On the target singlet pair they reduce exactly to logical \(X_L,Z_L\). The ancilla factor below is a gauge singlet, so conditioning these operators on ancilla state preserves gauge invariance.

**Limits:** The model's fixed matter occupations and one-rishon-per-link truncation are stipulated. This commutator proof does not supply full Yang–Mills dynamics or determine which physical interaction strengths an experiment can realize.

## 4. A constructive ancilla-controlled logical Hamiltonian

Let the ancilla projector be \(Q_a=|1\rangle\langle1|_a\). For \(P\in\{X,Z\}\) use the just-defined **physical Hermitian SU(2)-scalar** \(\widetilde P_{\rm all}\). Set \(\hbar=1\) and consider

\[
\boxed{
K_P=\frac{\pi}{2}\;Q_a\otimes
       (I-\widetilde P_{\rm all}),\qquad
U_P=e^{-iK_P}.
}
\]

**All-sector property.** \(K_P\) is Hermitian and commutes with every local Gauss generator across the entire 1,024-dimensional truncated module, including when the rishon ends are elsewhere. Thus \(U_P\) is unitary and gauge-invariant globally, not only on the selected code. We do **not** claim \(U_X\) is a controlled-Pauli outside the chosen two-dimensional encoding; \(\widetilde X_{\rm all}\) can have additional eigenvalues there.

**On the encoded singlets:** \(P_L^2=I\). Therefore the encoded Hamiltonian block satisfies \(K_P^2=\pi K_P\), and the exponential is exactly

\[
e^{-iK_P}=I-\frac2\pi K_P
\quad\Rightarrow\quad
\boxed{U_P=
|0\rangle\langle0|_a\otimes I_L+
|1\rangle\langle1|_a\otimes P_L.}
\]

This is the required **controlled-\(X_L\)** or **controlled-\(Z_L\)**. No numerical matrix exponential is used: the script proves the exact minimal-polynomial identity and unitary \(4\times4\) target matrix symbolically.

The Hamiltonian terms are explicit. For example,

\[
K_X=\frac{\pi}{2} Q_a\otimes I-
\frac{\pi}{2\sqrt3} Q_a\otimes
(\mathsf T_{15}-\mathsf T_{14}).
\]

These are **ancilla-controlled exchange couplings**. If the flag/syndrome ancilla can coherently modulate the relative strength of two on-site SU(2)-invariant spin exchanges, the gate follows. But the terms \(Q_a \otimes \mathbf S_i\cdot\mathbf S_j\) are **effective three-body interactions**, not an automatically available native two-body exchange gate. The program does not synthesize them into hardware gates or model microscopic errors during synthesis.

## 5. Consistency with the existing fault-propagation simulator

The earlier \`flagged_su2_measurement.py\` abstracted the ancilla-controlled inner operations as ordinary CNOT and CZ gates. To test that this abstraction is correct **on the encoded Hilbert space**, the new verifier:

1. Constructs both nontrivial controlled gate matrices exactly from \(K_X\) and \(K_Z\).
2. Verifies both act correctly on the *full singlet codeword vectors* and that the ancilla-0 block is identity.
3. Compares the conjugation of **all 16 two-qubit Pauli words** (control ancilla and target singlet qubit) against the simulator's binary-symplectic \`CX\` and \`CZ\` transformations, up to physically irrelevant Pauli signs.
4. Checks a genuine Bell-state entangling witness \(CX(|+\rangle_a|0_L\rangle)=(|0_a0_L\rangle+|1_a1_L\rangle)/\sqrt2\).

All 16 signatures for **each** gate match. This gives an **exact algebraic bridge from SU(2) matter/rishon operators to the ideal encoded Clifford simulator**. It is **not** a fault-propagation equivalence between hypothetical elementary physical ancilla-spin gates and one abstract logical controlled-P gate: that missing compilation may introduce correlated physical faults.

## 6. Primary-source error catalog: the published five-qubit flag circuit

A separate script, [\`published_flag_figure2_catalog.py\`](published_flag_figure2_catalog.py), is grounded in the visually inspected **arXiv:1705.02329v1, PDF page 2, Figure 2(c,d)** for Chao and Reichardt, [*Quantum Error Correction with Only Two Extra Qubits*](https://doi.org/10.1103/PhysRevLett.121.050502), *Physical Review Letters* **121**, 050502 (2018).

**Figure correction:** The relevant flagged \([[5,1,3]]\) circuit is **Figure 2(c)** in that original arXiv manuscript, and its fault catalog is **Figure 2(d)**. The previous [\`ADAPTIVE_FLAG_EC_AUDIT.md\`](ADAPTIVE_FLAG_EC_AUDIT.md) referred to “Figure 1(c)” and should not be read as an accurate figure identifier for that specific version.

The original paper lists four two-Pauli fault labels after circuit point **b** and four after **c**; the eight rows yield **seven** distinct nontrivial data errors:

| Original location | Fault labels | Resulting data errors |
|---|---|---|
| b | IZ, XZ, YZ, ZZ | IIZXI, IXZXI, IYZXI, IZZXI |
| c | IZ, XZ, YZ, ZZ | IIIXI, IIXXI, IIYXI, IIZXI |

The program independently checks that the seven distinct data Pauli errors have **seven distinct nonzero four-bit stabilizer syndromes** for the stated five-qubit code. This reproduces the specific correctability fact used by the published flag construction.

**Strict separation:** This catalog is a **verified transcription of one source figure and a syndrome-check regression**. It is **not** a gate-by-gate simulation of Figure 2(c), and our ancilla-controlled \(|+\rangle/X\)-readout gadget uses a different ancilla-basis convention from the source's \(|0\rangle/Z\)-readout syndrome ancilla. The flag ancilla conventions differ accordingly. Even a basis-equivalent ideal Clifford gadget cannot automatically be assumed to have the same physical fault locations or native gate noise.

## 7. Achieved, not achieved, and next step

| Research criterion | Status |
|---|---|
| SU(2)-singlet logical inner \(X,Z\) from co-located exchanges | **Exact verified** |
| Gauge-invariant extension across all truncated orientation sectors | **Exact verified** |
| Explicit Hermitian ancilla-controlled Hamiltonian for logical CX/CZ | **Exact constructed** |
| Controlled logical gates agree with ideal Clifford error propagation | **32 exact Pauli conjugations passed** |
| Published Figure 2(d) seven flagged errors have distinct syndromes | **Exact verified** |
| Native two-body hardware decomposition of conditioned exchange | **Not implemented** |
| Gate-for-gate reproduction of published Figure 2(c) | **Not implemented** |
| Physical noise/threshold analysis and leakage-safe recovery | **Not implemented** |
| Standard Model SU(3) color, Yukawa dynamics, anomaly-derived code realized | **Not established** |

The next milestone should specify an *actual elementary interaction set* (e.g. ancilla-controlled exchange mediated by a neutral bus, permitted interactions and geometry), synthesize \(U_X,U_Z\) using that set, and enumerate physical faults **at the elementary-gate level**, not at the ideal controlled-P logical level. In parallel, transcribe Figure 2(c)'s native gate directions and measurement-basis conjugations exactly, then compare its fault catalog under the **same** noise assumptions. Until both are complete, the result is an explicit **gauge-invariant Hamiltonian proposal** with algebraic certificates, not a physical fault-tolerant implementation.

## Reproduction

\`\`\`bash
python -m pip install sympy==1.14.0
python su2_controlled_exchange.py
python published_flag_figure2_catalog.py
\`\`\`

GitHub Actions runs these tests alongside the prior charge-code, non-Abelian singlet, recovery and flag experiments.

**Primary sources:** Chao and Reichardt, [arXiv:1705.02329v1 (see Fig. 2)](https://arxiv.org/pdf/1705.02329), published as [Phys. Rev. Lett. **121**, 050502 (2018)](https://doi.org/10.1103/PhysRevLett.121.050502); Bacon, Kempe, Lidar and Whaley, [Phys. Rev. Lett. **85**, 1758 (2000)](https://doi.org/10.1103/PhysRevLett.85.1758), on encoded operations within SU(2)-symmetric spin-singlet subspaces.
