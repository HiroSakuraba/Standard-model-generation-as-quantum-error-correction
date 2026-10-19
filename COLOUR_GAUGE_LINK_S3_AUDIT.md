# Two-vertex gauge-link-dressed exchange and local Gauss pointer: finite non-Abelian audit

19 October 2026. Follows the merged nine-quark colour algebra PR #46.

## Research question

Can a local gauge-field link make an exchange interaction between separately gauged colour registers gauge invariant? Can one **actually construct** an ancilla measurement of the local Gauss constraints, and will that measurement reveal errors that preserve gauge symmetry?

**Yes, yes, and no, respectively — in the finite-group construction below.** The third answer is the important barrier for attempting to recover the original nine-quark quantum memory.

## Precise minimal system

Take the six-element non-Abelian permutation group \(S_3\), embedded faithfully into \(SU(3)\) via \(D(g)=\operatorname{sgn}(g)P(g)\), where \(P(g)\) permutes the three standard colour axes. Every \(D(g)\) is a \(3\times3\) special unitary matrix. This **does not** reproduce continuous SU(3) gauge symmetry: the representation restricted to \(S_3\) is reducible, and finite-group charge sectors differ from physical colour multiplets.

Place two fundamental colour qutrits at vertex A, two at vertex B, and one oriented group-element link \(|g\rangle\), \(g\in S_3\), between them. The total Hilbert space has dimension \(3^4\cdot6=486\). Independent local transformations are

\[
R_A(h): |a_0,a_1,g,b_0,b_1\rangle
 \mapsto D(h)^{\otimes2}|a_0,a_1\rangle\otimes|hg\rangle\otimes|b_0,b_1\rangle
\]
\[
R_B(h): |a_0,a_1,g,b_0,b_1\rangle
 \mapsto |a_0,a_1\rangle\otimes|gh^{-1}\rangle\otimes D(h)^{\otimes2}|b_0,b_1\rangle.
\]

The matter and link factors are ordered in code as \(|g\rangle\otimes|a_0a_1b_0b_1\rangle\); the formulas above are illustrative orderings. The vertex transformations commute with one another. Both Gauss projectors are exact finite-group averages: \(P_v=|S_3|^{-1}\sum_h R_v(h)\).

The intersection of the two invariant subspaces has dimension **14**, constructed with an explicit isometry W: \(\mathbb C^{14}\to\mathbb C^{486}\). This is not the original 21-dimensional nine-quark code.

## Dressed edge exchange, bare exchange negative control

Define the edge-local but three-register (matter A0, link, matter B0) operation

\[
S_{\text{dress}} |g\rangle\otimes|a,b\rangle
   = |g\rangle\otimes (D(g)\otimes D(g)^\dagger)\operatorname{SWAP}|a,b\rangle,
\]

where the other local matter qutrits are untouched. It is unitary and its square is identity, and **commutes with both independent vertex gauge actions**. It can be described at a high level as a cross-edge SWAP plus two link-controlled matter rotations. Those conditional operations are *not* decomposed into a universal/native entangling-gate set and might be expensive.

A bare SWAP between A0 and B0 is not invariant under independent local gauge transformations. In the actual calculation, dressed-SWAP covariance residual on all 14 physical basis states is exactly zero, while the bare-SWAP covariance residual is about **5.196** (Frobenius norm, not a probability).

A perfect Gauss projector accepts a uniformly mixed state on the physical 14D sector with probability 1. After a **bare** cross-vertex SWAP, either single-vertex Gauss projector accepts with probability **0.25**. After the gauge-dressed SWAP, either accepts with probability **1**. Thus, a valid gauge-invariant exchange is **invisible to Gauss-only checks**. Gauss checks alone are not an error-correction syndrome for the logical quantum memory.

## Explicit locally supported vertex syndrome circuit

For each vertex, introduce a six-level coherent group-label pointer. The circuit is:

1. Prepare a uniform pointer \(|u\rangle=6^{-1/2}\sum_{h=0}^{5}|h\rangle\).
2. Apply the **one multiplexed local controlled-gauge operation** \(U_v=\sum_h |h\rangle\langle h|\otimes R_v(h)\). Every branch acts only on the two matter qutrits at that vertex and its incident link. The multiplexor has six branches; this is **not six proven native gates**.
3. Apply the inverse normalized six-point discrete Fourier transform to the pointer and measure its six output basis labels.

The \(k=0\) Kraus operator is exactly \(P_v\), whereas \(k=1,\ldots,5\) together contain the other gauge-charge outcomes. The script computes all output branches, verifies \(\sum_k M_k^\dagger M_k=I\) via norm conservation on arbitrary states, and verifies acceptance of every exact physical basis state. The Fourier transform is a *label-space* transform, not an irreducible Fourier transform on the non-Abelian group.

To encode all **six labels** of this implementation needs six pointer levels (e.g. at least two qutrits); if one merely wants a binary pass/fail Gauss test, this is **not** a lower bound of two qutrits. There are no hardware/native two-qudit-gate counts for controlled vertex transformations, pointer preparation, Fourier operations, or readout.

As a clearly synthetic illustration, a uniform 1% **six-way classical pointer assignment error** yields 99% reported acceptance on a clean or dressed-swapped gauge state, and about 24.9% acceptance for a bare-swapped state under a single-vertex check. That is not a calibrated readout model.

## The remaining logical-syndrome obstruction

On the **14-dimensional** gauge-invariant sector, the dressed swap has eigenvalue +1 with multiplicity 10 and -1 with multiplicity 4. Pairing four opposite-sign eigenvectors yields an artificial **4-dimensional** gauge-invariant code K such that K† S_dress K=0.

For the **specific known-pair exchange channel** with Kraus operators \(\sqrt{1-p}I\) and \(\sqrt p\,S_{\text{dress}}\), the mathematical three-outcome instrument \(\{P_K,\,P_{SK},\,I-P_K-P_{SK}\}\), conditional S_dress feedback, and decoding recover the logical four-level state perfectly: ideal entanglement fidelity **1**.

**This three-outcome logical instrument is nonlocal and not compiled.** Its exact success is a mathematical baseline, not hardware protection. The implemented local vertex Gauss-pointer circuit gives the *same trivial gauge syndrome* for both K and S_dress K. It cannot replace the nonlocal logical syndrome.

This toy code and recovered four-level logical space must not be identified with the 21-dimensional sign-twist nine-quark code or its **65** independent one-quark error syndromes. The 65-label recovery remains unimplemented.

## External source context

The finite group-link basis, its left/right gauge transformations, and gauge-invariant matter-sector projector are standard in lattice gauge theory. See *Counting gauge-invariant states with matter fields and finite gauge groups*, European Physical Journal C (2025), https://doi.org/10.1140/epjc/s10052-025-15050-8; and *Quantum simulation of lattice gauge theories in more than one space dimension—requirements, challenges and methods* (2022), https://pmc.ncbi.nlm.nih.gov/articles/PMC8886423/. A separate source, *Quantum circuits for SU(3) lattice gauge theory*, Physical Review D (2025), https://doi.org/10.1103/k8f6-yft8, makes clear that *actual* continuous SU(3) Gauss constraints and link truncation require additional representation theory. None verifies our device-dependent costs.

## Test and next decision

    OPENBLAS_NUM_THREADS=1 python -m unittest -v test_colour_gauge_link_s3
    OPENBLAS_NUM_THREADS=1 python colour_gauge_link_s3.py

Six regression tests check the finite SU(3)-subgroup representation and multiplication, independent vertex actions, 14D physical Gauss space, gauge-covariant dressed SWAP, Kraus completeness of the **explicit** local pointer circuit, and the known-error logical counterexample.

**Next useful step:** either (a) compile the local vertex controlled-group multiplexor to a particular hardware gate set and count its depth/error; or (b) search for a *locally implemented gauge-covariant measurement of the dressed-exchange logical syndrome*, rather than just Gauss-charge detection. Only after (b) could the nine-quark 65-syndrome physical recovery question be addressed in earnest. Simulating the *full SU(3)* Wilson link is a separate, substantially larger step.
