# Exact active four-logical-state obstruction for all S4 matter faults plus electric-link error

19 October 2026. Follow-up to [PR #61](https://github.com/HiroSakuraba/Standard-model-generation-as-quantum-error-correction/pull/61), which established by rational arithmetic that the commutant of the entire S4 matter-exchange algebra plus one standard-irrep link-electric projector has **complex dimension six**. That ruled out a four-dimensional *noiseless subsystem*, but not active quantum error correction.

**This experiment proves a stronger, exact, and carefully scoped statement: no four-dimensional code satisfies the Knill–Laflamme (KL) conditions when the error span includes every one of the 24 S4 matter-permutation operators and a gauge-invariant electric-link phase error.** It also classifies passive four-state subspaces for eight smaller swap-generator subgroups combined with the electric error.

## 1. The precise fault set, versus its three generators

The two-vertex finite-S3 gauge model has four matter qutrits, one six-state group link, a 486D total Hilbert space and a 14D local-Gauss-invariant subspace.

The three physical gauge-covariant operators A=onsite SWAP(A0,A1), B=onsite SWAP(B0,B1) and S=dressed SWAP(A0,B0) generate an order-24 group isomorphic to the symmetric group S4. They account for all pair permutations of the four matter qutrits, including gauge-dressed cross-link permutations.

**The no-go error set is every one of those 24 elements as separately possible Kraus error operators,** plus the link-electric phase fault \(R=I-2E\), where \(E=P_{\rm standard}\) is the single-link central projector of PR #60. A stochastic probability distribution assigning nonzero support to all 25 unitary errors defines a physically legitimate but presently *synthetic* quantum channel. The no-go concerns this **25-dimensional operator-error span**, not a sparse channel containing only \(\{I,A,B,S,R\}\). Distinct generators may generate the same algebra as Hamiltonians without generating the same span of single-step Kraus errors, and exact KL conditions depend on the latter.

There are no physical calibration rates or microscopic jump operators making all 24 events equally likely. We do **not** claim that realistic hardware produces this whole fault set; it is a completeness/possibility bound within the finite-group model.

## 2. Exact S4 representation decomposition

In the 14D Gauss sector, use the exact integer colour-orbit representation of PR #61: four signed-permutation fundamental factors cancel the common sign, so the global S3 gauge-fixed invariant space has exactly **14 integer orbit sums** of four colour labels.

All 24 position permutations act as 14×14 integer matrices. The five irreducible characters of S4 are integer-valued. Construct the exact rational central idempotents

\[
P_\lambda=\frac{d_\lambda}{24}\sum_{g\in S_4}
\chi_\lambda(g^{-1})\,U(g).
\]

Using the ordinary S4 conjugacy classes \(1,\ (12),\ (12)(34),\ (123),\ (1234)\), SymPy evaluates every \(P_\lambda\) with exact fractions, verifies all idempotence/orthogonality identities and yields:

| S4 irreducible type | Irrep dimension | Multiplicity | Total isotypic dimension |
|---|---:|---:|---:|
| Trivial | 1 | **4** | 4 |
| Sign | 1 | 0 | 0 |
| Standard [3,1] | 3 | **2** | 6 |
| Sign-twisted standard [2,1,1] | 3 | 0 | 0 |
| Two-dimensional [2,2] | 2 | **2** | 4 |
| **Total** | | | **14** |

Thus the physical Hilbert space decomposes exactly as \( (\mathbb C^1\otimes\mathbb C^4) \oplus (\mathbb C^3\otimes\mathbb C^2) \oplus (\mathbb C^2\otimes\mathbb C^2)\).

## 3. Exact active Knill–Laflamme impossibility for the full error span

Suppose \(K:\mathbb C^4\to\mathcal H_{\rm physical}\) encodes four logical states and perfectly corrects each of the 24 permutations. KL applied to \(I\) and every other group error requires

\[
K^\dagger U(g)K=c_g I_4,\qquad g\in S_4.
\]

By linearity, the same scalar-compression condition holds for **every operator in the complex group-algebra image** \(\operatorname{span}\{U(g)\}\), including all irreducible-block matrix units \( |a\rangle\langle b|\otimes I_{m_\lambda}\).

For each irrep \(\lambda\), resolve the code's isometry into maps \(V_{\lambda,a}:\mathbb C^4\to\mathbb C^{m_\lambda}\). The matrix-unit KL relations imply

\[
V_{\lambda,a}^\dagger V_{\lambda,b}
=c_{\lambda,ab}I_4.
\]

If any component is nonzero, at least one \(c_{\lambda,aa}>0\), so \(V_{\lambda,a}\) is injective and has rank **four**. Hence the multiplicity space must satisfy \(m_\lambda\ge4\).

**Only the trivial S4 representation has multiplicity four.** Every 4D code correcting the full matter-permutation error set must therefore occupy the unique four-dimensional *trivial S4 sector*—the passive all-swap invariant memory already constructed in PR #60.

But the electric projector on that unique 4D sector has exactly

\[
\operatorname{spec}(P_{\rm triv} E P_{\rm triv}|_{\operatorname{ran}P_{\rm triv}})
= \left\{0,\frac13,\frac23,1\right\}.
\]

The exact rational characteristic polynomial on the full 14D space is
\(x^{11}(x-1)(3x-1)(3x-2)/9\); the ten additional zero eigenvalues arise from the complement of the four-dimensional trivial sector. As a result \(E\), and \(R=I-2E\), **do not compress to a scalar on that full four-dimensional space**. The KL condition for the error pair \(I,R\) fails.

Therefore

\[
\boxed{\text{No perfect active 4D code corrects all 24 matter permutations plus this link-electric phase fault.}}
\]

This proof uses exact integer/rational character and polynomial identities plus the elementary injectivity argument above. It is stronger than the absence of a noiseless subsystem and has no floating-point rank threshold. It is **not** a no-go for an error set with only three individual swap generators plus R, nor for approximate error correction, changing the code dimension, adding independent physical resources, or full SU(3).

## 4. Exhaustive exact classification of smaller exchange subgroups plus R

For every subset of the three swap generators A, B, S, take the **whole generated subgroup** as a possible exchange error family, and add the electric-phase R. A passive decoherence-free *subspace* is simultaneously scalar under every group element and R. Each selected involutory generator therefore has scalar eigenvalue ±1, and E has eigenvalue 0 or 1.

The code exhaustively enumerates these sign assignments and solves their common eigenspace ranks using **exact rational nullities** of stacked 14×14 constraints. The maxima are:

| Subgroup generator choice | Group order | Largest passive common eigenspace | Four-state passive code? |
|---|---:|---:|---|
| None | 1 | 9 | Yes |
| A alone | 2 | 6 | Yes |
| B alone | 2 | 6 | Yes |
| S alone | 2 | 5 | Yes |
| A, B (commuting disjoint swaps) | 4 | 4 | Yes |
| A, S (order-six subgroup) | 6 | **2** | **No** |
| B, S (order-six subgroup) | 6 | **2** | **No** |
| A, B, S (full S4) | 24 | **1** | **No** |

This is an exhaustive passive **common-eigenvector subspace** census for these particular error families, not a classification of all noiseless subsystems, KL-correctable codes, or every gauge-compatible interaction. Some proper order-six subgroups lack passive four-state protection against E, so they are natural targets to test **active** KL conditions. Their actively correctable four-state capacity remains **open** in this PR. Do not infer that an active code exists merely because passive common eigenstates do not.

For the full S4 case the result is stronger: the separate theorem above excludes even perfect active four-state correction when all 24 swaps are possible Kraus errors.

## 5. What this means for the research program

We can now draw an exact boundary around the small gauge model:

- **All swaps alone:** four logical dimensions enjoy complete passive protection.
- **All swaps plus a link-electric fault:** neither passive four-state storage nor any perfect active four-state code works for the *full 25-operator error span*.
- **Selected proper subgroups plus electric:** the first mathematically interesting unresolved cases are the order-six exchange subgroups generated by A,S or B,S, for which the largest passive common eigenspace has only two dimensions.
- **Changing error assumptions can change the conclusion.** The presence of a symmetry alone has not yielded a physically advantageous code; actual device noise and resource costs must determine whether a restricted subgroup error family is meaningful.

The sensible next experiment is to check the **exact KL conditions and maximal logical dimension** for the two order-six subgroup-plus-electric fault sets, preferably with a constructive encoding/recovery when one exists or an exact code-dimension argument when it does not. Then introduce a *single device-motivated stochastic channel* with calibrated rates, preparation, transport and readout costs before comparing physical performance.

**Still unresolved:** continuous SU(3), real Standard-Model matter, electric-link hardware synthesis, nine-quark 21-dimensional/65-syndrome physical recovery, or any measured gauge-QEC advantage.

## Reproduce

\`\`\`bash
python -m pip install numpy scipy sympy
OPENBLAS_NUM_THREADS=1 python -m unittest -v test_colour_s3_exact_active_kl_bound
OPENBLAS_NUM_THREADS=1 python colour_s3_exact_active_kl_bound.py
\`\`\`

A dedicated CI job archives exact S4 character projector multiplicities, electric polynomial, active KL impossibility certificate, all eight passive subgroup tables, and explicit scientific limitations in JSON.
