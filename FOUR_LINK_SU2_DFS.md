# A gauge-invariant SU(2) logical qubit: what it protects and what breaks it

**Status — exact constructive result for a finite, deliberately truncated matter/rishon model.** This is not a Standard Model implementation or a discovery of a new universal quantum code. The executable reference is [\`four_link_su2_subspace.py\`](four_link_su2_subspace.py), requiring SymPy 1.14. It builds on the exact two-site Gauss-law calculation in [\`SU2_MATTER_GAUSS_KL.md\`](SU2_MATTER_GAUSS_KL.md).

## Main result

A two-site **SU(2)** quantum-link *kinematic* model, with one color-doublet matter excitation at each vertex and **four** distinguishable links carrying one spin-half rishon each, has a **16-dimensional gauge-singlet sector**. Eight allowed assignments of link endpoints each supply two independent color-singlet states.

**Positive result:** Fix the endpoint assignment \((0,1,1,1)\). Its entire **two-dimensional** singlet space is a *decoherence-free subspace* for any error operator that acts **only as a diagonal function of the four link-end occupancies**. For this noise algebra, Knill–Laflamme is satisfied exactly; the recovery is the identity and no syndrome needs to be measured.

**Negative result:** An SU(2)-invariant **local color-spin exchange** at the vertex containing four color doublets acts as a *nontrivial logical phase flip* on that very same code. Thus the code **does not** correct the physically allowed error set \(\{I,\mathrm{SWAP}\}\), despite obeying Gauss's law throughout.

The result is narrower but more instructive than simply saying non-Abelian gauge symmetry fails: one can obtain a genuine, precisely delimited passive code *within* the physical singlet space, but its noise selectivity is essential.

## 1. Exact Hilbert space and multiplicity

The model has two vertices \(v=0,1\), four parallel links, two spin-\(1/2\) matter degrees of freedom, and four spin-\(1/2\) rishons. Each rishon occupies exactly one endpoint of its link. No matter-empty/doubly-occupied sectors, link representation changes, or dynamical Hamiltonian are included.

There are \(2^4\) link-end configurations and \(2^6\) color states per configuration: \(\dim\mathcal H_\mathrm{kin}=1,024\). For a configuration with \(l\) left-end rishons, the two vertices contain \(1+l\) and \(5-l\) fundamental SU(2) doublets. Local Gauss law projects each site's tensor product onto its total-spin-zero subspace.

A collection of \(2r\) fundamental spin-half degrees of freedom has singlet multiplicity \(C_r\), the \(r\)-th Catalan number; an odd collection has none. In the four-link model, a physical configuration requires \(l\in\{1,3\}\): there are four choices for each value. The vertex with two doublets has one singlet; the vertex with four doublets has **two** linearly independent singlets. Therefore

\[
\dim\mathcal H_\mathrm{phys}=4(1\times2)+4(2\times1)=16.
\]

The verifier computes every configuration's exact SU(2) spin-zero multiplicity, checks the total 16, and confirms the earlier four-link result \((8,8)\) split by first-link orientation.

## 2. Explicit codewords in a fixed endpoint sector

Fix the link-end assignment \(o=(0,1,1,1)\), where rishon \(r_0\) is at vertex 0 and \(r_1,r_2,r_3\) are at vertex 1. Let \(m_v\) denote the matter doublet on vertex \(v\). The normalized SU(2) singlet of two doublets is

\[
|s(a,b)\rangle =
\frac{|\uparrow_a\downarrow_b\rangle-|\downarrow_a\uparrow_b\rangle}{\sqrt2}.
\]

Form two normalized states, in the *same* endpoint sector:

\[
\begin{aligned}
|A\rangle &= |s(m_0,r_0)\rangle|s(m_1,r_1)\rangle|s(r_2,r_3)\rangle,\\
|B\rangle &= |s(m_0,r_0)\rangle|s(m_1,r_2)\rangle|s(r_1,r_3)\rangle.
\end{aligned}
\]

Their exact overlap is \(\langle A|B\rangle=\tfrac12\), so they are linearly independent but not orthogonal. An orthonormal **logical basis** is

\[
\boxed{|0_L\rangle=|A\rangle,\qquad
|1_L\rangle=(2|B\rangle-|A\rangle)/\sqrt3.}
\]

Both are exact SU(2) singlets at **both** vertices. The program expands each into the six-doublet computational basis, verifies orthonormality, and applies all local \(G_v^z\) and \(G_v^+\) operators exactly (with no numeric tolerances).

## 3. A precisely protected class of physically gauge-invariant noise

Let \(N_{e,L}\) denote the left-end occupation number of the rishon on link \(e\). Each \(N_{e,L}\) is a **color-scalar, gauge-invariant, locally measurable** number operator in this model.

Since both logical states share the *same* endpoint pattern, on the code projector \(P\) all occupancies are scalars:

\[
P N_{e,L}P = \begin{cases}P,&e=0\\0,&e=1,2,3.\end{cases}
\]

Indeed every **diagonal function** \(f(N_{0,L},N_{1,L},N_{2,L},N_{3,L})\), including finite combinations and phases of link-end electric energies, acts as \(f(1,0,0,0)\,I_L\) on the logical qubit. Thus for any two Kraus operators \(E_\alpha,E_\beta\) in this algebra,

\[
P E_\alpha^\dagger E_\beta P = c_{\alpha\beta}P.
\]

This is an exact **decoherence-free subspace (DFS)**, a restricted form of quantum error protection. The logical qubit does not dephase under arbitrary fluctuations coupling *only* to these occupations. No active correction is required.

This is **not** protection against all SU(2)-invariant operators. It also assumes a separate preparation/confinement mechanism keeps states inside the fixed-pattern singlet sector; the model does not supply a preserving Hamiltonian.

## 4. A local, gauge-invariant exchange breaks the code

On vertex 1 there are **four** doublets \(m_1,r_1,r_2,r_3\). Exchange the *color-spin labels* of \(m_1\) and \(r_1\) without moving rishons or changing link occupancy:

\[
E_{\mathrm{swap}}=\mathrm{SWAP}_{m_1,r_1}.
\]

Because a simultaneous SU(2) rotation acts identically on these two doublets, swapping them **commutes with all Gauss generators**. The script checks this commutation exactly against every six-spin computational-basis vector.

In the logical basis above,

\[
P E_{\mathrm{swap}}P =
\begin{pmatrix}-1&0\\0&+1\end{pmatrix}= -Z_L
\quad
\text{(up to the chosen sign convention).}
\]

The associated physical local spin-exchange term is

\[
\mathbf S_{m_1}\cdot\mathbf S_{r_1}
=\frac12\mathrm{SWAP}_{m_1,r_1}-\frac14 I,
\]

with exact logical restriction

\[
P(\mathbf S_{m_1}\cdot\mathbf S_{r_1})P =
\begin{pmatrix}-3/4&0\\0&+1/4\end{pmatrix}.
\]

These **distinct eigenvalues** produce logical dephasing while leaving every SU(2) Gauss-law syndrome unchanged. For the error set \(\{I,E_{\mathrm{swap}}\}\), Knill–Laflamme would require \(PE_{\mathrm{swap}}P=cP\), which fails.

This is an exact *counterexample to broader protection claims*, not an argument that physical SU(2) implementations necessarily experience this noise at any particular rate.

## 5. Key distinction from earlier results

- The two-link gauge-singlet space had dimension two and was **not protected even against the local link-occupation phase error**: both endpoint occupations occurred in the physical basis.
- In the four-link model, one can choose a **proper two-dimensional subspace of the 16-dimensional gauge-invariant space** in a fixed endpoint sector. That subspace *is* protected against all occupation-diagonal dephasing.
- It nevertheless fails against an SU(2)-invariant *spin-exchange error*, despite both errors respecting the same gauge group.
- Neither the two-link nor four-link encoding has been shown to realize the original anomaly-derived charge-code matrices, the observed Standard Model flavor mixing, or a complete circuit-level recovery implementation.

This illustrates why the error algebra must be specified **before** describing a gauge-invariant structure as useful quantum error correction. The non-Abelian symmetry identifies singlet channels, but it does not automatically suppress physical interactions *among* those channels.

## 6. Reproduce and subsequent experiment

\`\`\`bash
python -m pip install sympy==1.14.0
python four_link_su2_subspace.py
\`\`\`

All comparisons use exact SymPy expressions and small integer representations; no stochastic or floating-point approximation is needed. The GitHub Actions workflow runs this check alongside the prior matter, charge-code, and lattice-gauge tests.

**Next research target:** incorporate a specific SU(2)-gauge-invariant Hamiltonian on these four links, compute its matrix on the 16-dimensional singlet space, and identify whether any symmetry sector or energetic separation can suppress the dangerous exchange relative to occupation-diagonal noise. Only then would this DFS have an argued physical lifetime. Beyond that, build the same analysis with dynamical matter and more realistic link representations. Those tests remain unperformed.

### References

- E. Knill and R. Laflamme, [*Theory of quantum error-correcting codes*](https://doi.org/10.1103/PhysRevA.55.900), *Physical Review A* **55** (1997), 900.
- [The preceding exact SU(2) singlet/rishon result](SU2_MATTER_GAUSS_KL.md), including its scope and source references.
