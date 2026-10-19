# Exact finite-S3 swap-plus-electric commutant certificate and restricted active-code construction

19 October 2026. Builds on merged [PR #60](https://github.com/HiroSakuraba/Standard-model-generation-as-quantum-error-correction/pull/60), where a numerically six-dimensional commutant under all matter swaps plus a link-electric projector was reported, with an explicitly stated lack of exact arithmetic certification.

## 1. Correctly reduce the physical gauge Hilbert space using integer orbit sums

The model has four signed-permutation-colour matter qutrits (two at each vertex), one six-element group-valued S3 link, full dimension 486 and a 14-dimensional simultaneous Gauss-invariant space. Gauge fixing the link to identity identifies the physical states with global *diagonal* S3 invariants in \((\mathbf3)^{\otimes4}\).

The signed-permutation representation contributes \(\operatorname{sgn}(g)\) on every fundamental. With **four** tensor factors the total sign is \(\operatorname{sgn}(g)^4=1\). Consequently the global invariant subspace has an exact integer basis: sums over orbits of length-four colour words under simultaneous permutation of the three colours. There are **14** orbits, one of size three and thirteen of size six. These correspond to partitions of four labelled sites into at most three indistinguishable colour classes.

This construction does not use floating point to select a degenerate 14D eigenspace. It also establishes a concrete, independent correspondence to the original 486D physical gauge embedding.

## 2. All matter-exchange words form S4, and its commutant has an exact 24-element basis

On the orbit-sum basis, every elementary matter-site exchange permutes the *positions* of a length-four colour word. Thus the on-site exchanges A0↔A1 and B0↔B1, together with the gauge-dressed A0↔B0 exchange, act as **integer permutation matrices**. They generate the entire 24-element S4 positional permutation group.

For a permutation representation on 14 basis labels, a 14×14 operator commutes with every S4 matrix if and only if its entries are constant on each diagonal S4-orbit of **ordered pairs** of those 14 labels. There are **24 such orbits**. Their characteristic/indicator matrices have disjoint support and are therefore an **exact linearly independent basis** of the S4 commutant, over either rationals or complexes.

This provides a theoretical basis for PR #60's swap-only commutant dimension 24, without relying on double-precision singular values.

## 3. Link electric energy has rational entries in the orbit-sum basis

The central standard-irrep Peter–Weyl projector on the six-level link is

\[
E=P_{\mathrm{standard}}=\frac{2}{6}
\sum_{g\in S_3}\chi_{\mathrm{standard}}(g)\,L_g.
\]

At the identity-link gauge slice, \(L_g\) is equivalent to the inverse group action on the two B-vertex colour qutrits. The two signed-representation factors cancel, and the standard character takes values **2 on identity, 0 on transpositions, −1 on 3-cycles**. Thus the gauge-fixed operator \(E\) is rational with denominators dividing three. It is an **exact idempotent of rank nine**, agreeing with the compressed six-link numerical operator in PR #60.

In the **unnormalized** orbit-sum basis this rational E need not be symmetric under the usual Euclidean inner product, because the orbit sums have unequal sizes. In the normalized orbit basis, a diagonal similarity transformation with entries \(\sqrt{|{\cal O}_i|}\) restores the physical Hermitian representation. The program verifies this normalized model against the original 486D gauge-link operator for A, B, S, and E, with numerical Frobenius residuals below \(2\times10^{-11}\). The exact commutant dimension is invariant under similarity.

## 4. An exact rational, independently cross-checked commutant certificate

Write a general operator commuting with S4 as \(X=\sum_{j=1}^{24}x_j C_j\), where \(C_j\) are the exact integer ordered-pair-orbit indicator matrices.

The extra condition \(XE=EX\) is a **196×24 matrix over the rational field**, because every entry of E is rational. Its **exact SymPy rational rank is 18**, with an exactly verified **six-dimensional rational nullspace**. Therefore:

\[
\boxed{\dim_{\mathbb C}\{A,B,S,E\}'=24-18=6.}
\]

As an independent integer check, clear the denominators in A, B, S, E and form the **784×196** full commutator linear system. Its rank is **190 modulo 101** and again **190 modulo 65521**, giving a rigorous lower bound of 190 on the characteristic-zero rank. This agrees with the rational argument, which gives characteristic-zero rank \(196-6=190\).

Unlike PR #60, the **specific six-dimensional commutant finding is now exact**, conditional only on the code implementing the explicitly stated finite-group model. No interval arithmetic or floating-point threshold is used in the proof itself.

A noiseless subsystem containing a four-dimensional logical factor needs a full \(M_4(\mathbb C)\) operator algebra in the commutant, of complex dimension at least \(4^2=16\). An exact six-dimensional commutant cannot contain it. Consequently:

**No four-dimensional noiseless subsystem exists under the *-algebra generated by all four matter permutations and this specific link-electric projector in the finite S3 model.**

This statement is **not** a no-go theorem against four-dimensional *active* quantum error-correcting codes for restricted error sets; Knill–Laflamme conditions must be assessed separately.

## 5. A constructive but nonadvantageous restricted active code

As a controlled positive example, consider only the proper restricted error set

\[
\{I,\ S_{A0,B0},\ R\},\qquad R=I-2E.
\]

The dressed swap S has a **ten-dimensional +1 eigenspace**. Compress E to it. The characteristic eigenvalues are

\[
0,0,\quad \tfrac12,\tfrac12,\tfrac12,\quad 1,1,1,1,1.
\]

Construct four orthonormal logical vectors: two equal superpositions of mutually orthogonal zero- and one-eigenvectors, and two independent one-half-eigenvectors. The resulting code K satisfies

\[
SK=K,\qquad K^\dagger R K=0,\qquad
K^\dagger K=I_4,\qquad (RK)^\dagger RK=I_4.
\]

Thus the binary syndrome subspaces K and RK are orthogonal. An **ideal, global**, nondemolition measurement followed by R feedback exactly corrects a known-probability R fault and is also immune to S (which acts as identity on K). The numerical construction confirms the KL residuals at approximately machine precision, but it is *not* a physical recovery circuit.

**However, this is not evidence for active correction being superior.** The intersection of S=+1 and E=+1 eigenspaces already has **dimension five**, supporting a four-dimensional *passive* code on which S=+I and R=−I. Hence the restricted channel does not need active syndrome extraction at all. In particular, constructing *some* active code is easier than demonstrating that active correction achieves protection absent from all passive choices.

The comparison faithfully retains the distinction between a logical code, a mathematically valid syndrome instrument, and a physically synthesized circuit.

## 6. Unresolved questions

The exact commutant obstruction applies to the **full algebra** generated by S4 matter permutations and the single-link standard-irrep electric projector. It does not prove that the entire physically permitted operator algebra is the same or that no active four-dimensional code can correct selected microscopic error operators.

The next meaningful tasks are to identify an **error set outside the small passive intersection**, seek an optimized active code and explicit recovery satisfying every KL compression for that set, and benchmark that recovered channel with *the same microscopic control and physical fault model* as a properly optimized passive memory. It may be more fruitful to prove a four-dimensional packing/linear-algebra obstruction than to search for ever larger circuits without a code.

No native gauge-electric pulse, no link-transport Hamiltonian, no real device T1/Tphi, no physical code preparation, no 21D nine-quark/65-syndrome implementation, and **no measured gauge-QEC advantage** are provided.

### Reproduction

\`\`\`bash
python -m pip install numpy scipy sympy
OPENBLAS_NUM_THREADS=1 python -m unittest -v test_colour_s3_exact_electric_commutant
OPENBLAS_NUM_THREADS=1 python colour_s3_exact_electric_commutant.py
\`\`\`

The exact rational matrix, nullspace and independent modular ranks are reconstructed from the group definitions at runtime. The machine-readable report includes every check and every claim boundary.
