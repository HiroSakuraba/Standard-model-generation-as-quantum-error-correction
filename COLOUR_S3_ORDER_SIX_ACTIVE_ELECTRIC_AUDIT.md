# Active four-state quantum code for proper order-six gauge-compatible exchange groups and electric-link fault

19 October 2026. This follows merged PR #62, which **exactly excluded all four-dimensional active codes against the complete 24-swap Kraus family plus the electric fault**. That no-go did *not* apply to proper six-element matter permutation subgroups. This experiment resolves the two formerly open subgroup cases, **AS** and **BS**, and identifies a genuine mathematical active-versus-passive distinction in the existing finite-S3 two-vertex gauge model.

## 1. The exact, bounded physical error model

Four colour fundamentals and a six-level finite group S3 gauge link form a 486D Hilbert space, with a **14D** simultaneous local-Gauss-invariant sector. Use the exact integer/rational 14D gauge-fixed colour-orbit-sum representation from PR #61, independently embedded back into 486D with normalized orbit amplitudes.

There are three gauge-invariant matter exchanges:

- A: onsite SWAP(A0,A1);
- B: onsite SWAP(B0,B1);
- S: gauge-link-dressed SWAP(A0,B0).

Either proper subgroup \(\langle A,S\rangle\) or \(\langle B,S\rangle\) has **six** elements. The electrical link phase is \(R=e^{-i\pi E}=I-2E\), where \(E=P_{\mathrm{standard}}\) is the single-link central standard-irrep projector that commutes with both Gauss actions.

The precise allowed **single-step error operators** are the six elements of **one** selected subgroup, **plus R as a seventh alternative**. A stochastic channel with nonzero weights on those alternatives is correctable by the constructed code, regardless of their relative weights.

**The error set does not include arbitrary sequential composite faults \(gR\) or \(Rg\), the full 24-element matter permutation group, or all physically gauge-compatible quantum operations.** Those additions would require distinct Knill–Laflamme analysis; PR #62 rules out a four-dimensional code when all 24 permutations and R are individually in the error span.

## 2. Passive storage fails, including passive noiseless subsystems

The representation of each order-six exchange subgroup on the 14D physical space contains a **six-dimensional trivial sector** (every subgroup exchange acts as the identity), and a two-dimensional irreducible representation with multiplicity four. The sign representation is absent.

For this subgroup plus R, PR #62 exhaustively checked the ranks of all common scalar eigenspaces: the largest has **dimension two**. Therefore there is no passive four-dimensional decoherence-free *subspace*.

A stronger criterion rules out *noiseless subsystems* as well. The group acts by 14×14 integer permutations on the gauge-fixed orbit patterns. Its commutant has an **exact 52-element integer basis**, one indicator matrix for each simultaneous group orbit on ordered pairs of physical basis indices. Requiring a general linear combination to commute with the **rational electric projector E** yields an exact **196×52 matrix of rank 41**. Consequently, for **each** of AS and BS,

\[
\dim_{\mathbb C} \{ \langle A,S\rangle,E\}' =
\dim_{\mathbb C}\{ \langle B,S\rangle,E\}'=52-41=\mathbf{11}.
\]

A four-dimensional noiseless logical *subsystem* would require a complete \(M_4(\mathbb C)\) algebra of operators inside this commutant, which itself would need dimension at least **16**. Hence neither joint error algebra admits such a passive subsystem. This conclusion is from **exact rational arithmetic**, not floating-point singular-value rank.

The passive obstruction does *not* preclude active error correction, because the latter relies on conditions for the linear **single-step error span**, rather than requiring the code to be invariant under the entire generated algebra.

## 3. Construct a four-logical-state code that corrects the entire seven-error set

Write \(P_G=(1/6)\sum_{g\in G}g\) for the **six-dimensional exchange-invariant sector** of either order-six group. It is a rational projector. The compression \(P_G E P_G\) has exact eigenvalue multiset on that sector:

\[
\underbrace{0}_{1},
\quad\underbrace{1/3}_{1},
\quad\underbrace{2/3,2/3}_{2},
\quad\underbrace{1,1}_{2}.
\]

The full 14D characteristic polynomial includes eight additional zero roots from outside the fixed sector; it is verified exactly by SymPy.

Select an orthonormal basis \(e_0,e_{1/3},e_{2/3}^{(1)},e_{2/3}^{(2)},e_1^{(1)},e_1^{(2)}\) of the compressed E. Define the four logical columns

\[
\begin{aligned}
 c_0&=e_{2/3}^{(1)},\\
 c_1&=e_{2/3}^{(2)},\\
 c_2&=\sqrt{1/3}\,e_0+\sqrt{2/3}\,e_1^{(1)},\\
 c_3&=\sqrt{1/2}\,e_{1/3}+\sqrt{1/2}\,e_1^{(2)}.
\end{aligned}
\]

For \(C=[c_0,c_1,c_2,c_3]\), these columns have mutually disjoint spectral support except where the eigenvalue is already degenerate. Consequently

\[
C^\dagger C=I_4,\qquad
gC=C\quad (\forall g\in G),\qquad
C^\dagger E C=\frac23 I_4.
\]

The electric phase error \(R=I-2E\) then satisfies

\[
C^\dagger R C=-\frac13 I_4.
\]

All Knill–Laflamme matrices \(C^\dagger F_i^\dagger F_j C\) for \(F_i,F_j\in G\cup\{R\}\) are scalar multiples of identity: subgroup permutations act trivially on the code, and \(R\) has scalar compression \(-I/3\). This is **perfect active correctability of four logical dimensions** for the entire specified seven-operator fault family.

The exact group and E spectral identities establish the existence theorem without numerical search. Numerical double precision is used only to choose an orthonormal realization of the codewords and verify their residuals in the original complex Hilbert space.

## 4. A real logical syndrome and recovery, not an orthogonality shortcut

The original error and code are **not fully orthogonal**: \(C^\dagger R C=-I/3\). Define the orthogonal error-syndrome isometry

\[
V=\frac{RC+C/3}{\sqrt8/3}.
\]

Then \(C^\dagger V=0\), \(V^\dagger V=I_4\), and

\[
RC=-\frac13 C+\frac{\sqrt8}{3}V.
\]

An ideal nondemolition projective measurement of \(P_C=CC^\dagger\) and \(P_V=VV^\dagger\), with a third unused/complementary outcome, followed by a partial isometry \(V\to C\) on the second branch, exactly restores the encoded four-state information. The error R gives code-syndrome probability **1/9** and orthogonal-syndrome probability **8/9**. Subgroup permutation errors always produce the original code branch.

The simulator checks all **seven possible error Kraus operators** separately, retaining both syndrome branches and accounting for probabilities rather than postselecting. Each branch's decoded logical Kraus operator is a scalar multiple of \(I_4\), and the **unconditional recovered entanglement fidelity is 1** in the ideal model.

This is an **abstract/global recovery specification**, not a compiled finite-S3 or physical SU3 syndrome circuit. The earlier PRs #48–57 construct hardware-inspired syndrome implementations for *different code spaces*; one cannot transfer those gate counts or fidelity values to this new code without a fresh compilation.

### Explicit two-step fault negative control

The program also tests **electric phase R followed by each possible matter permutation g**, giving composite \(gR\), which is **outside** the original seven-alternative *single-fault* set. The same ideal recovery now fails on some such sequences. In the evaluated group representation, four nontrivial composites retain only about **0.4583** of their probability inside the intended code-plus-electric-syndrome projectors; the remaining **0.5417** is flagged rather than being renormalized away. Their unconditional recovered entanglement fidelity is approximately **0.1146**, not 1. Other subgroup elements happen to commute suitably and preserve recovery, but that does **not** imply correction of all repeated errors.

This makes the boundary operational: **perfect one-fault correction does not imply fault tolerance under two sequential errors**. A microscopic error process that naturally produces correlated exchange-electric sequences would require a redesigned decoder or another code. Virtual high-energy processes and omitted field fluctuations could be relevant to the frequency and structure of those compound errors in hardware.

## 5. How virtual particles fit—and what remains absent

The entire research program currently has two different dynamical resolutions:

1. **Explicit virtual superexchange (PRs #55–58):** the Fermi-Hubbard model contains high-energy doublon–hole occupations and actual fermionic hopping. Those intermediate states mediate effective colour exchange. Its 495D full-Fock extension **retains their coherent occupation and residual leakage throughout multiple operations**. The term *virtual excitation* is standard in the Hubbard effective-Hamiltonian literature; it does not mean the simulator creates new Standard Model particles. As a relevant source, see *Schrieffer-Wolff transformations for experiments: Dynamically suppressing virtual doublon-hole excitations in a Fermi-Hubbard simulator*, Physical Review A **106**, 012428 (2022), DOI https://doi.org/10.1103/PhysRevA.106.012428.
2. **Exact symmetry/encoding test here:** a finite S3 gauge Hilbert space with fixed matter-particle count and a six-level link. There is **no dynamical QCD gluon field, no quark–antiquark pair creation, no virtual gluon or vacuum-polarization Feynman loops, and no renormalization-group running**. Its seven allowed errors are explicitly specified gauge-compatible operators, not derived rates from the QCD vacuum.

In quantum field theory, internal Feynman diagram lines represent *virtual contributions*, not ordinarily detectable free particles. These loop and exchange effects matter for the **correct effective Hamiltonian and natural error spectrum** of a real particle-physics realization; absent effects might generate new allowed noise that defeats our code. But they do **not invalidate** an exact theorem about a precisely specified finite error family. See CERN's QCD review of quark/gluon loops: https://cds.cern.ch/record/724415/files/p173.pdf.

A future full-Standard-Model program would need physical quark/gluon degrees of freedom, continuum/regulated lattice dynamics, treatment of sea quarks and vacuum polarization, and nonperturbative gauge dynamics. Whether that changes the error-correctability of the **original nine-quark 21D/65-syndrome encoding** remains unknown.

## 6. Scientific result and next release criterion

This PR supplies a mathematically **positive active-vs-passive result in the finite toy model**: a four-dimensional active code is perfect against a seven-alternative gauge-compatible fault set for which **no equally large noiseless subspace or subsystem exists**.

This is **not** a demonstrated advantage for physical quantum memory. The conditioned projector/partial-isometry recovery is global and currently unsynthesized, the charged mediator and link-electric terms are not native calibrated device Hamiltonians, and the seven error rates are synthetic and restricted to individual faults. We haven't analyzed repeated/compound errors, spectator/reference loss, state preparation, measurement errors, T1/T2, or optimized passive codes under a matched noisy device.

**Next engineering milestone:** synthesize the new \(P_C/P_V\) syndrome and \(V\to C\) feedback with local Gauss-preserving operations, preferably by exploiting a charged-mediator architecture; then include the full Hubbard microscopic error, link-electric noise strength, and matched passive controls. **Next fundamental milestone:** determine whether the code survives sequential exchange-plus-electric errors and whether an enlarged continuous SU3 theory retains any analogous syndrome.

### Reproduce

\`\`\`bash
python -m pip install numpy scipy sympy
OPENBLAS_NUM_THREADS=1 python -m unittest -v test_colour_s3_order_six_active_electric
OPENBLAS_NUM_THREADS=1 python colour_s3_order_six_active_electric.py
\`\`\`

The dedicated GitHub workflow archives exact group/algebra ranks, exact compressed-electric characteristic polynomials, the two explicit four-logical isometries, all seven ideal recovered Kraus branches and 486D Gauss compatibility, along with scientific limitations.
