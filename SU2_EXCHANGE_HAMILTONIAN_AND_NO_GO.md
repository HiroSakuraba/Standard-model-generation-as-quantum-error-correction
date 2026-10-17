# SU(2) matter encoding: a gapped symmetric Hamiltonian and an exact noise-algebra obstruction

**Status, 17 October 2026:** source-informed, mathematically reproducible research for the *specific* two-matter-site / four-parallel-link / one-rishon-per-link SU(2) truncation. This is **not the full Standard Model**, not the original anomaly-derived charge code, and not a newly invented four-spin decoherence-free encoding. The main executable check is [\`su2_exchange_algebra.py\`](su2_exchange_algebra.py); preceding constructions are in [\`FOUR_LINK_SU2_DFS.md\`](FOUR_LINK_SU2_DFS.md).

## Executive result

We can engineer a **gauge-invariant Hamiltonian with a two-dimensional, spectrally gapped ground space** corresponding to the exact SU(2) matter/rishon singlet qubit. Uniform collective spin exchange does not distinguish the logical states. In fact, equal individual couplings are not necessary: **three precise pair-sum equalities** are sufficient and necessary within a six-coupling Heisenberg exchange model.

But the encoded states remain **unprotected against generic local SU(2)-invariant exchange perturbations**. Two on-site exchanges act as noncommuting qubit observables and generate the full matrix algebra on the logical space. Including the occupancy operators, they generate a full **eight-block error algebra** on the 16-dimensional gauge-singlet sector. No encoding of a two-dimensional logical Hilbert space corrects *every* error in that **complete specified error algebra**.

The distinction is between **protected against a restricted noise class** and **protected by gauge symmetry itself**. The result is neither a universal SU(2) no-go theorem nor a physical lifetime estimate.

## 1. Physical sector and exact local exchange matrices

Each of the four links holds one fundamental rishon at one endpoint. There is one fixed matter doublet per vertex. The 16 endpoint configurations each have \(2^6=64\) spin assignments: the full truncated kinematic Hilbert space has 1,024 dimensions.

Exactly **eight endpoint configurations** are color-singlet-admissible, each with **two independent** singlet channels. The full gauge-invariant space is therefore \(\mathcal H_{\rm phys}=\bigoplus_{o=1}^{8}\mathbb C^2\), dimension 16.

Within the chosen endpoint configuration \((0,1,1,1)\), vertex 0 has two spin-half degrees of freedom and vertex 1 has four. Let the four at vertex 1 be ordered \((1,3,4,5)\): matter spin 1 followed by rishons 1,2,3. The exact orthonormal basis is

\[
|0_L\rangle = |s(0,2)\rangle|s(1,3)\rangle|s(4,5)\rangle,\qquad
|1_L\rangle = (2|s(0,2)\rangle|s(1,4)\rangle|s(3,5)\rangle-|0_L\rangle)/\sqrt3.
\]

The script verifies **each** of the eight sectors independently, constructs normalized valence-bond basis states, checks both vertices' non-Abelian Gauss conditions exactly, and projects all local spin SWAPs onto the chosen logical basis.

Let \(X_L\) and \(Z_L\) be ordinary Pauli matrices with \(Z_L=\operatorname{diag}(1,-1)\). The six pair exchanges \(\mathrm{SWAP}_{ij}\) project to:

| Pair | Exact encoded matrix |
|---|---|
| \((1,3),(4,5)\) | \(-Z_L\) |
| \((1,4),(3,5)\) | \((Z_L-\sqrt3 X_L)/2\) |
| \((1,5),(3,4)\) | \((Z_L+\sqrt3 X_L)/2\) |

They are involutions, and their total sum is **zero**. Already \(\mathrm{SWAP}_{13}\) and \(\mathrm{SWAP}_{14}\) fail to commute and, together with the identity and their product, span the **four complex dimensions of \(M_2(\mathbb C)\)**. These are SU(2)-invariant, *co-located* spin exchanges; no gauge charge is violated.

The same facts hold for **all eight** allowed endpoint configurations, with their four-spin vertex relabeled accordingly. The code checks the rank-four span in each sector.

## 2. Arbitrary local Heisenberg couplings: exact splitting condition

For local spin halves, \(\mathbf S_i\cdot\mathbf S_j=(\mathrm{SWAP}_{ij}-I/2)/2\). Therefore the six-coupling intravertex Hamiltonian

\[
H_{\rm ex}=\sum_{i<j\in\{1,3,4,5\}}J_{ij}\,\mathbf S_i\cdot\mathbf S_j
\]

has an exact **two-by-two** representation. Define the three sums

\[
A=J_{13}+J_{45},\qquad B=J_{14}+J_{35},\qquad C=J_{15}+J_{34}.
\]

Then

\[
\boxed{
H_L=-\frac{A+B+C}{4}I_L+
\frac{-2A+B+C}{4}Z_L+
\frac{\sqrt3(C-B)}4 X_L.
}
\]

The logical energy gap squared is

\[
\boxed{(E_+-E_-)^2=
\frac{(-2A+B+C)^2+3(C-B)^2}{4}.}
\]

For real couplings, the encoded energies are equal **if and only if \(A=B=C\)**. The equality can hold even with unequal individual couplings; the new program verifies an explicit nonuniform example. This is a **balance condition**, not a generic symmetry guarantee.

For completely uniform exchange \(J_{ij}=J\), the projected Hamiltonian is

\[
H_L=-\frac32J\,I_L.
\]

Uniform exchange only measures the total spin Casimir and hence cannot distinguish two SU(2) spin-zero states. This familiar mechanism is directly related to the established *four-spin decoherence-free subspace* literature (Bacon, Kempe, Lidar and Whaley 2000/2001). The implementation here is an independent **matter/rishon sector bookkeeping calculation**, not a claim to a new code.

A perturbation \(\delta J_{13}=\delta\) and all other deviations zero produces \(\operatorname{diag}(-3\delta/4,\delta/4)\), so the **first-order logical splitting is \(|\delta|\)**. Increasing a penalty on non-singlet or wrong-endpoint configurations cannot remove this in-code perturbation.

## 3. Positive example: an exactly gapped but fine-tuned gauge-invariant Hamiltonian

Let \(N_{e,L}\) be the left-link-end occupation and take the target pattern \(n^*=(1,0,0,0)\). Let \(G_v^2\) be the local SU(2) total-spin Casimir. Consider the explicitly **engineered** kinematic Hamiltonian

\[
H_0=\Delta\sum_{e=0}^3(N_{e,L}-n_e^*)^2+
\lambda\sum_{v=0}^1 G_v^2+
J\sum_{v=0}^1\sum_{i<j\ \mathrm{at}\ v}\mathbf S_i\cdot\mathbf S_j.
\]

When a rishon is absent at a vertex, its spin does not enter that vertex's exchange term. This is a well-defined, color-scalar operator on every endpoint occupation sector. Take \(\Delta>0\) and \(g:=\lambda+J/2>0\).

In a sector with \(n_v\) local color doublets,

\[
\sum_{i<j\ \mathrm{at}\ v}\mathbf S_i\cdot\mathbf S_j=
\frac12\left(G_v^2-\frac34n_v\right).
\]

Since \(n_0+n_1=6\) in every sector,

\[
H_0=-\frac{9J}{4}I+
\Delta\,d_{\rm Hamming}(o,n^*)+
g(G_0^2+G_1^2).
\]

This is **gauge invariant**, diagonal in endpoint configurations, and positive above its constant baseline. It has **exactly two** ground states, the two target-sector SU(2) singlets, both with energy \(-9J/4\).

The complete excitation gap is

\[
\boxed{\Gamma=\min\{2g,\ \Delta+\tfrac32g,\ 2\Delta\}>0.}
\]

The three minima are attained respectively by (i) a local spin-1 excitation with the correct endpoints, (ii) one changed endpoint with local half-integer spins, and (iii) another color-singlet sector two endpoint flips away. The script enumerates the exact minimum-spin Casimirs of all 16 endpoint configurations and verifies this formula and sample parameter choices. No \(1024\times1024\) numerical diagonalization or uncontrolled approximation is required.

**Crucial limit:** This is a *spectral gap*, not an error-correcting distance. It protects against leaving the ground manifold and yields **no suppression at all** of a gauge-invariant local exchange perturbation that has a nonscalar two-by-two restriction inside the ground manifold. Energetic isolation does not eliminate first-order logical dephasing. The construction is fine-tuned and does not demonstrate a self-correcting memory or a natural Standard Model Hamiltonian.

## 4. Stronger obstruction for a declared **complete local error algebra**

The local SU(2)-invariant **endpoint occupation operators** \(N_{e,L}\), together with on-site **matter–rishon and rishon–rishon color exchanges whenever both spins are present at a vertex**, form a naturally specified error algebra.

Within the **gauge-singlet space**:

1. The four commuting number operators resolve all eight admissible endpoint configurations. Their products yield the eight orthogonal projectors \(P_o\).
2. Inside each sector \(P_o\mathcal H_{\rm phys}\cong\mathbb C^2\), two local SU(2)-invariant spin swaps are noncommuting and span, with their product and identity, all of \(M_2(\mathbb C)\).
3. Multiplying these exchanges by endpoint-occupation projectors, which commute with the Gauss generators, yields the full block-diagonal algebra

\[
\boxed{\mathcal A=\bigoplus_{o=1}^8 M_2(\mathbb C),\qquad
\dim_{\mathbb C}\mathcal A=8\times4=32.}
\]

**No-two-dimensional-code proposition.** Suppose errors may contain **arbitrary operators in \(\mathcal A\)** (equivalently an error basis whose products span this algebra), and suppose a two-dimensional code \(V:\mathbb C^2\to\mathcal H_{\rm phys}\) satisfies Knill–Laflamme for that entire error algebra. Write \(V_o=P_o V\). For any rank-one matrix \(|w\rangle\langle w|\) supported in block \(o\), KL would require

\[
V_o^\dagger|w\rangle\langle w|V_o=c_{o,w}I_2.
\]

The left side has rank at most one, while \(I_2\) has rank two. Consequently \(c_{o,w}=0\), and because this must hold for every \(w\), \(V_o=0\). Doing so in **all eight blocks** contradicts \(V^\dagger V=I_2\). Therefore **no 2D subspace of the 16D physical singlet sector can correct arbitrary noise from this complete specified algebra**.

The program independently checks the full rank-four algebra generated by two projected swaps in each of all eight sectors and confirms the 32-dimensional block structure; the rank-one KL argument is a proof for **any** putative encoded qubit, not just the fixed-endpoint code.

**Very important scope:** An error *algebra* is much stronger than a narrow physical error set such as one endpoint dephasing channel or one exchange perturbation. This theorem **does not** exclude subspace codes under smaller channels, higher-dimensional redundant physical systems, extra stabilizers, autonomous correction, active measurement, or arbitrary SU(2) quantum-memory designs. It does not assign a likelihood or rate to the individual processes in \(\mathcal A\).

## 5. Research implications

This provides a constructive positive result **and** a bounded obstruction:

- A gauge-invariant Hamiltonian with a genuine **two-fold degenerate, gapped ground subspace** exists in our tractable model.
- Collective or exactly pair-balanced exchange is scalar on the encoding, so the model is a decoherence-free subspace for this restricted class.
- Generic **local SU(2)-invariant pair exchange** generates logical \(X\)- and \(Z\)-type operators, acting entirely within the gapped ground subspace.
- For the complete local occupancy-and-exchange noise algebra, **no protected encoded qubit exists anywhere in the fixed four-link gauge-singlet sector**.
- No feature of the anomaly-derived Standard Model charge matrices supplied this encoding or its Hamiltonian. Their operational implementation remains an open, different problem.

**Next useful milestone:** either identify a *smaller physically justified noise model* that preserves the pair-sum equality to measured precision, or enlarge the Hilbert space with independent redundant color-singlet blocks and test the full KL equations under explicitly local exchange errors. Those are separate, falsifiable directions; neither has been done here.

## Reproduce and source comparisons

\`\`\`bash
python -m pip install sympy==1.14.0
python su2_exchange_algebra.py
\`\`\`

All results use exact representation-theory identities, exact rational/square-root matrix entries, and finite matrix spans. GitHub Actions runs the new test with previous anomaly-code and SU(2) tests.

**Relevant primary literature:**

- Bacon, Kempe, Lidar and Whaley, [*Universal Fault-Tolerant Quantum Computation on Decoherence-Free Subspaces*](https://doi.org/10.1103/PhysRevLett.85.1758), *Physical Review Letters* **85**, 1758 (2000): an already-known four-spin singlet DFS and its exchange interactions.
- Kempe, Bacon, Lidar and Whaley, [*Theory of decoherence-free fault-tolerant universal quantum computation*](https://doi.org/10.1103/PhysRevA.63.042307), *Physical Review A* **63**, 042307 (2001): exchange-generated logical operations and collective symmetry protection.
- [*Spin Exchange-Enabled Quantum Simulator for Large-Scale Non-Abelian Gauge Theories*](https://doi.org/10.1103/PRXQuantum.5.030358), *PRX Quantum* **5**, 030358 (2024): a relevant rishon/matter SU(N) gauge-simulator proposal explicitly using spin exchange and local singlet constraints; it does **not** imply a noise-protected Standard Model memory.
- Knill and Laflamme, [*Theory of quantum error-correcting codes*](https://doi.org/10.1103/PhysRevA.55.900), *Physical Review A* **55**, 900 (1997): the exact algebraic quantum-recovery criterion.
