# Exact pair-register compilation of the finite-S3 logical syndrome

19 October 2026. Builds on merged PRs #47 (two-vertex non-Abelian gauge link) and #48 (one-qutrit logical phase measurement).

## Main result

The PR #48 syndrome circuit used **two controlled link-dressed swaps** and **two controlled on-site matter swaps**, where each controlled edge operation acted on four registers (pointer, matter A0, group link, matter B0). This batch derives an **exact circuit using only two-register unitary gates**, checks the entire resulting pointer map on the 486-dimensional matter/link model, and records Gauss-law violations **after every constituent gate**.

### Two-register identity for an ordinary qutrit controlled swap

On two qutrits A and B, the SWAP is +1 on the six-dimensional symmetric sector and −1 on the three-dimensional antisymmetric sector. Choose any exact two-register unitary Q_AB taking the antisymmetric sector into the three computational states with **B=0**, and the symmetric sector into the six basis states with B≠0. Then

\[
\mathrm{SWAP}_{AB}=Q_{AB}^{\dagger}[I-2(I_A\otimes|0\rangle\langle0|_B)]Q_{AB}.
\]

Consequently, SWAP_AB conditioned on pointer values in a specified set M⊂{0,1,2} is the product of exactly **three general two-register operations** in time order:

1. Q_AB on the A–B pair.
2. A selective controlled phase −1 on the pointer–B pair only when pointer∈M and B=0.
3. Q_AB† on A–B.

This is a constructive *exact* factorization, not an approximate numerical optimum. It uses an arbitrary two-qutrit basis transform Q_AB and an assumed selective two-qutrit phase operation; **neither is yet compiled to a restricted device-native gate basis**.

### Dressed edge exchange

Let the oriented group link hold g∈S3 and D(g) be its signed-permutation SU(3) subgroup representation. Set

\[
C_{LB}=\sum_g|g\rangle\langle g|_{\rm link}\otimes D(g)_{B0},
\]

a two-register 18×18 unitary on link and B0. The gauge-covariant dressed exchange from PR #47 obeys

\[
S_{\rm dress}=C_{LB}^\dagger\,\mathrm{SWAP}_{A0,B0}\,C_{LB}.
\]

A *pointer-controlled* dressed exchange therefore decomposes into **five** general two-register operations, in time order: C_LB, Q_A0B0, selective pointer–B0 phase, Q†_A0B0, C†_LB.

An on-site pointer-controlled swap needs only the three factors involving Q, pointer phase, Q†. Thus the entire logical-syndrome measurement's four controlled operations become

\[
2\times5+2\times3=\boxed{16\ \text{two-register factors}},
\]

plus one single-register pointer uniform preparation and one single-register Fourier transform. A conditional un-controlled dressed-exchange recovery requires **three** two-register factors (C, ordinary SWAP, C†).

**This is a count of general U(9)/U(18) pair-register factors, not native CP, CNOT, SUM, two-qubit, or two-qutrit entanglers.** Link–matter interactions involve a 6-level system. Compiling the 18×18 controlled representation, the Q transform and selective pointer phase into a realistic gate basis, including routing, remains open.

## Intermediate Gauss-law exposure

Both constituent operations Q_AB and C_LB generally fail to commute individually with vertex Gauss operators even though their complete product **does**. The program checks each step on the exact 3-pointer×486-system amplitudes and reports the probability that vertex A and B Gauss projectors accept the state. It verifies that both acceptance probabilities return to 1 at the ends of all four controlled high-level operations and at final pointer readout.

This provides a quantitative measure of *transient* gauge-constraint violation under this algebraic gate factorization. It does not simulate a continuous-time gauge-penalty Hamiltonian or error introduced by gauge-breaking virtual states during a physical control pulse. **Intermediate Gauss failure can matter physically even if the ideal final product is gauge invariant.**

## Fault-injection experiment

The implementation applies the exact two-register operations directly to pointer+link+matter state tensors and compares the result to the original high-level PR #48 circuit on an arbitrary 486×2 reference input. Separately it Monte Carlo samples a **completely positive, trace-preserving** mixed-unitary error channel after each two-register factor. With probability p, an independently uniformly sampled nonidentity generalized two-register Weyl product is applied; otherwise no fault occurs. This is a synthetic noisy-pair scenario, not calibrated physical decoherence or a gauge-invariant noise channel.

The physical gauge-link logical code remains four-dimensional. The source known-dressed-swap error channel {√(1−q)I,√q S} is propagated with the actual compiled phase pointer, measured into three branches, and corrected conditionally with a further three two-register gates. The resulting **unconditional logical entanglement fidelity** is estimated over explicit stochastic sampled fault trajectories. The mean and Monte Carlo standard error are reported. The ideal zero-noise fidelity must be 1.

The compiled circuit and noise propagation do not establish a hardware threshold. No noisy comparable uncontrolled-storage protocol, coherent calibration model, gauge penalty, idle time, routing, or physical control of arbitrary U(18) is included.

## Strict evidence limits

- An **exact two-register decomposition** has been derived, not a native two-body physical gate synthesis. General U(18) link-matter operations may themselves cost many elementary physical entanglers; quoting "16 native gates" would be false.
- The finite S3 subgroup representation with four matter qutrits and one group link is not a continuous SU3 lattice, nine-quark 21-dimensional sign-twist code or its 65-syndrome decoder.
- Logical syndrome measurement here corrects **one known dressed-swap error**, not arbitrary unknown exchanges.
- The controlled operation and the final recovery preserve gauge symmetry **as total unitaries**, not necessarily at each intermediate two-register factor. The diagnostic reports intermediate violation explicitly.
- Only after choosing a *native* gate set and a gauge-penalty/control Hamiltonian can one quantify physical gate times, pulse-time gauge leakage, noise, or architecture competitiveness.

## Reproduce

Run from the repository root:

    python -m pip install numpy scipy
    OPENBLAS_NUM_THREADS=1 python -m unittest -v test_colour_s3_two_register_syndrome
    OPENBLAS_NUM_THREADS=1 python colour_s3_two_register_syndrome.py

GitHub Actions publishes the machine-readable gate counts, intermediate Gauss probabilities, exact unitary comparison, and per-shot stochastic fidelity estimates.

**Next:** synthesize/link the U(18) and U(9) factors into a specific device-native gate alphabet with explicit graph constraints, and test whether only using gauge-covariant native operators is possible or whether temporary gauge violation must be controlled by an energy penalty or echo.
