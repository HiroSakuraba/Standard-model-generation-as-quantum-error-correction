# A finite non-Abelian gauge theory with matter: an exact Gauss-law recovery counterexample

**Research status, 17 October 2026.** This is an exactly solvable, minimal **SU(2)-kinematic quantum-link/rishon truncation with matter**, not a lattice discretization of the full Standard Model. The result is an explicit *negative* Knill–Laflamme (KL) test: genuine non-Abelian gauge invariance can admit a nontrivial logical degree of freedom **and** an undetectable, physically gauge-invariant local error. It does **not** claim that no quantum error correction is possible for SU(2), or that the anomaly-derived code cannot be engineered on qudit hardware.

Executable proof: [\`su2_matter_singlet.py\`](su2_matter_singlet.py). It uses exact rational singlet amplitudes and exact SymPy ranks, without stochastic simulation.

## 1. Why this problem differs from anomaly cancellation

The prior [operational Abelian benchmark](OPERATIONAL_GAUGE_QEC.md) separated genuine physical gauge constraints from the **abstract stabilizer operators** defined by the Standard Model charge matrix. The outstanding question is whether the same obstruction survives when Gauss law is non-Abelian, color singlets are genuine representation-theoretic objects, and **matter** is present.

It does, in a small model. The point is not to construct an impressive code from two sites; it is to exhibit **the exact noncorrectable error** and **the encoded basis** in a model where SU(2) gauge transformations act explicitly.

## 2. Explicit microscopic Hilbert space and Gauss generators

Take a periodic two-vertex, two-link geometry (two parallel oriented edges joining vertices 0 and 1). Fix:

- One distinguishable matter color doublet (spin-1/2, no empty/doubly occupied matter sectors) at **each** vertex.
- One color doublet called a *rishon* per link. For each link the rishon may be at the vertex-0 endpoint or the vertex-1 endpoint, but not both. This is the fixed-single-rishon-per-link kinematic truncation.
- A color-up/down basis for the two matter doublets and two rishons.

There are four link-end assignments \((o_0,o_1)\in\{0,1\}^2\) and \(2^4=16\) color states per assignment: \(\dim\mathcal H_{\rm kin}=64\).

For every vertex \(v\), the local gauge generators in this model are the ordinary SU(2) spin generators

\[
G_v^a=S^a_{m_v}+\sum_{e:o_e=v}S^a_{r_e},
\qquad a=x,y,z.
\]

A physical **color singlet** satisfies \(G_v^a|\psi\rangle=0\) at both vertices. The script constructs the exact integer matrices for \(2G_v^z\) and \(G_v^+\) and computes their *joint kernel*, with no assumption that a simple parity test is sufficient. For a unitary SU(2) representation, a zero-weight state annihilated by the raising operator lies in a spin-zero irrep.

The results are:

| Link 0 end | Link 1 end | Color-singlet dimension |
|---|---|---:|
| 0 | 0 | 0 |
| 0 | 1 | **1** |
| 1 | 0 | **1** |
| 1 | 1 | 0 |

Why? If both rishons end at the same vertex, the *other* vertex contains only one matter doublet, which cannot form a singlet. When exactly one rishon ends at each vertex, each local matter–rishon pair has one singlet channel. The total gauge-invariant space therefore has **dimension two**.

## 3. Two exact orthonormal gauge-invariant states

For any matter spin \(m\) and rishon spin \(r\), define the SU(2) singlet

\[
|s(m,r)\rangle=
\frac{|\uparrow_m\downarrow_r\rangle-|\downarrow_m\uparrow_r\rangle}{\sqrt2}.
\]

The two physical basis states are

\[
|0_L\rangle=|o_0=0,o_1=1\rangle\otimes
  |s(m_0,r_0)\rangle|s(m_1,r_1)\rangle,
\]

\[
|1_L\rangle=|o_0=1,o_1=0\rangle\otimes
  |s(m_0,r_1)\rangle|s(m_1,r_0)\rangle.
\]

These live in distinct link-end sectors. Each product expands into four computational-basis states with coefficients \(\pm 1/2\). The program verifies, exactly, both normalization and orthogonality, **every local \(G_v^z\) constraint**, and every non-Abelian raising constraint \(G_v^+\). Combined with the exact joint-kernel dimension calculation, this proves they span the *entire* gauge-singlet sector of the stipulated 64-state model.

The encoded qubit exists as a **kinematic logical degree of freedom**. That alone is not evidence of quantum error correction.

## 4. A *local, gauge-invariant error* acts as logical Z

Define \(N_{0,L}\) to be the occupation number of the **left endpoint of link 0**. It is a local operator on a link-end degree of freedom, diagonal in endpoint occupancy, and **commutes with all local SU(2) generators** because it is a color scalar.

In the orthonormal gauge-singlet basis above, its exact matrix is

\[
P N_{0,L} P=
\begin{pmatrix}1&0\\0&0\end{pmatrix}.
\]

The corresponding local SU(2) electric Casimir is \(C^2_{0,L}=\tfrac34N_{0,L}\); it has restricted matrix \(\mathrm{diag}(3/4,0)\). Thus a gauge-invariant fluctuation in the local electric energy **distinguishes the two would-be logical states** without violating Gauss law.

Even more explicitly, the local **unitary** phase error

\[
U=\exp(i\pi N_{0,L})=(-1)^{N_{0,L}}
\]

is fully gauge invariant and preserves the exact physical singlet subspace, but acts there as

\[
PUP=\begin{pmatrix}-1&0\\0&1\end{pmatrix}.
\]

For the error set \(\{I,U\}\), the Knill–Laflamme condition would require \(PI^\dagger U P=cP\) for some scalar \(c\). This **fails exactly**, since the matrix is not proportional to the identity. Gauss-law measurements report the same zero-charge syndrome for \(I\) and \(U\).

This is a **physically meaningful class of errors in the defined quantum-link/rishon model**, not an arbitrary non-gauge-invariant Pauli fault. Its presence does not depend on choosing the basis of encoded logical states: non-scalar operator restriction is basis invariant on the full two-dimensional sector.

## 5. Interpretation boundaries

**What this proves:** In this explicit finite SU(2)-kinematic model with matter, the *entire Gauss-singlet space* cannot correct even the particular local gauge-invariant dephasing error \(\{I,(-1)^{N_{0,L}}\}\). Non-Abelian color-singlet constraints alone do not make that entire two-dimensional sector a quantum error-correcting code against that error set.

**What it does not prove:**

- A no-go theorem against all SU(2) matter encodings, code subspaces, additional magnetic/plaquette checks, or engineered Hamiltonians.
- A full quantum-link time evolution, fermionic hopping dynamics, or Kogut–Susskind SU(2) Hamiltonian. The truncation fixes one matter doublet per vertex and one rishon per link; it intentionally omits vacuum and double matter occupancy, other link representations, and an independent link-determinant term.
- That \(\mathrm{SU(2)}\) is color in the actual Standard Model (the color gauge group is \(\mathrm{SU(3)}_c\)); this is a tractable non-Abelian precursor, not a reconstruction of physical SM matter.
- That the two logical basis states are protected from generic errors, or that SU(2) symmetry forbids an additional independent syndrome observable.

The explicit microscopic single-rishon model should also **not** be confused with the **pure-gauge** error-correcting constructions in Yao (2026), which use a different truncated Hilbert space and additional encoded checks.

## 6. Reproduce and extend

\`\`\`bash
python -m pip install sympy==1.14.0
python su2_matter_singlet.py
\`\`\`

The script verifies all four exact singlet-sector dimensions, normalized singlet states, gauge generators, local SU(2)-invariant occupancy and Casimir, and an explicit failed KL matrix identity.

**Next meaningful experiment:** enlarge to 3–4 links and vary matter occupancy; construct the gauge-invariant local *error algebra* (occupancy, Casimir, and gauge-covariant hopping combinations) and compute its matrix restrictions on candidate code subspaces. Search for a proper subspace with \(\dim C\ge 2\) satisfying the full KL matrix equalities, and report the physical check operations and Hamiltonian needed to preserve it. This will test *whether additional redundancy can work*, rather than mistake Gauss invariance for error correction.

## References and relation to existing results

- Xiaojun Yao, [*Quantum error correction codes for truncated SU(2) lattice gauge theories*](https://doi.org/10.1103/5yg7-6mph), *Physical Review D* **113**, 114512 (2026). **Pure** truncated SU(2) gauge theory; not the matter/rishon model constructed here.
- [*Spin Exchange-Enabled Quantum Simulator for Large-Scale Non-Abelian Gauge Theories*](https://doi.org/10.1103/PRXQuantum.5.030358), *PRX Quantum* **5** (2024), describes rishon matter-site constructions and local singlet constraints. It informs the representation choice but does not provide this specific two-link KL calculation.
- [*Digital Quantum Simulation of a (1+1)D SU(2) Lattice Gauge Theory with Ion Qudits*](https://doi.org/10.1103/PRXQuantum.5.040309), *PRX Quantum* **5** (2024), constructs gauge-invariant dressed-site bases with matter and rishons.
- E. Knill and R. Laflamme, [*Theory of quantum error-correcting codes*](https://doi.org/10.1103/PhysRevA.55.900), *Physical Review A* **55**, 900 (1997), gives the general operator-correction criterion.

**Decision:** Do not expand the anomaly-derived code's physical claims on the basis of this result. Continue only with explicit gauge-invariant encodings, specified errors, and checkable recovery criteria.
