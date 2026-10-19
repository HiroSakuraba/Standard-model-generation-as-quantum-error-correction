# The fermionic 16: one generation as a Reed–Muller code

**Status, 19 October 2026:** exploratory note, not peer reviewed. Verified by `python3 spinor_reed_muller.py` (numpy; about a second; ends with `ALL FERMIONIC-GENERATION CHECKS PASSED`).

**Question.** One generation has 16 Weyl fermions, equivalently 32 Majoranas. The Standard Model with a ℤ4 symmetry is anomaly-free only when the number of fermions is a multiple of 16 (García-Etxebarria and Montero, [arXiv:1808.00009](https://arxiv.org/abs/1808.00009)). Does 16 matter for the code? And what happens when the code is built from fermions rather than qudits?

**Answer.**
- The generation's binary shadow is the 16-point Reed–Muller code, and puncturing ν^c gives the 15-qubit code with a transversal T gate.
- SO(10) is the smallest SO(2r) whose generation gives a non-Clifford transversal gate.
- Built from fermions, the 15-mode code's logical qubit is fermion parity, so it is forbidden by superselection. Adding ν^c repairs it.
- The mod-16 anomaly itself does not appear in these codes. A mod-8 count does.

## The binary shadow

Over odd primes, the spinor charges ±1 give the [[15,1,3]]_p code of the main paper. Mod 2, both +1 and −1 are 1, and the lattice collapses. So write each weight as bits, b = (1 − s)/2.

For SO(2r), the 2^{r−1} spinor weights become the even-weight subspace of 𝔽₂^r, and the affine functions on that subspace form the Reed–Muller code RM(1, r − 1). For SO(10), the weight enumerator is 1 + 30z⁸ + z¹⁶, and ν^c (all signs +) is the origin.

Puncturing the origin gives the quantum Reed–Muller codes. Transversal R_k = diag(1, e^{2πi/2^k}) acts logically, as R_k†, exactly up to k = r − 2:

| Group | Weights | Code | X / Z checks | Transversal gate |
|---|---|---|---|---|
| SO(6) | 4 | [[3,1,1]] | 2 / 0 | Z |
| SO(8) | 8 | [[7,1,3]] (Steane) | 3 / 3 | S (Clifford) |
| **SO(10)** | **16** | **[[15,1,3]]** | **4 / 10** | **T (non-Clifford)** |
| SO(12) | 32 | [[31,1,3]] | 5 / 25 | √T |

SO(10) is the first rung with a non-Clifford gate. Its code has the same n, k and check counts as the paper's [[15,1,3]]_p. Prakash and Ramya call their code the qudit analogue of the 15-qubit Reed–Muller code; this makes the match exact at p = 2.

For SO(10), the script also checks the 2¹⁵-dimensional state vectors directly: T^{⊗15} is logical T†, and the Knill–Laflamme conditions hold for all 46 Pauli errors of weight at most one.

## Fermions instead of qubits

Give each of the 16 Weyl fields a complex mode with Majoranas γ_j and γ′_j. CSS checks become Majorana operators:

- an X check on a support S is the product of the γ_j for j in S;
- a Z check on a support T is the parity of the modes in T.

Two Majorana products on supports A and B commute iff |A||B| − |A ∩ B| is even. The script validates this rule with explicit Jordan–Wigner matrices on every support of 6 Majoranas, then checks that all checks commute.

**Without ν^c (15 modes).**
- The X checks have rank 4, each a product of 8 Majoranas; the parity checks have rank 10.
- The logical X is a product of 15 Majoranas, so it is fermion-odd.
- The logical Z is the fermion parity (−1)^F itself: |0_L⟩ has N ∈ {0, 8} and |1_L⟩ has N ∈ {7, 15}.

A superposition of logical states mixes fermion parities, and superselection forbids it. The spin-ℤ4 generator i^N, whose square is (−1)^F, acts as logical S†.

**With ν^c (16 modes).** The unpunctured RM(1,4) code has X checks of rank 5 and parity checks of rank 11, so k = 0 and the code is a single state:

  (|0_L⟩|0⟩_ν + |1_L⟩|1⟩_ν) / √2.

Every branch has N ∈ {0, 8, 16}, so the state is fermion-even. It is invariant under i^N and under e^{iπN/4}, but not under e^{iπN/8}. The right-handed neutrino is the reference mode that makes the logical qubit physical. This is the fermionic version of the superselection-with-reference-frame point in the paper: ν^c is the frame.

## A mod-8 count, not mod 16

Take time reversal with γ → γ. Quadratic terms iγγ are odd and quartic terms are even. A gapped symmetric stabilizer state on 2n Majoranas then needs commuting real products of 4, 8, … Majoranas generating a maximal group, which is a doubly-even self-dual binary code. Exhaustive search finds none for 2n = 2, 4, 6 and one at 2n = 8. Gleason's theorem says such codes exist exactly when 2n ≡ 0 mod 8. This is the Fidkowski–Kitaev ℤ8 ([arXiv:0904.2197](https://arxiv.org/abs/0904.2197)) in code language.

One generation, with 32 Majoranas, passes: RM(2,5) is doubly-even and self-dual, and the script checks both. Fifteen Weyl fields give 30 Majoranas and fail. But the condition is mod 8 in Majoranas, which is mod 4 in Weyl fields. Four or eight Weyl fields would pass too.

## What 16 does and does not do

| Claim | Status |
|---|---|
| 16 = 2⁴ weights → RM(1,4) → transversal T after puncturing ν^c | proved by computation |
| ν^c restores fermion-parity superselection | proved by computation |
| the 16-mode code state is invariant under ℤ4 and ℤ8 phases | proved; it fails for ℤ16 |
| the Dai–Freed ℤ16 anomaly appears as a code property | **not found**; these zero-dimensional codes see only mod 8 |

The match between the 16 of the spinor and the 16 of the anomaly is not explained by these codes.

## Scope

The fermionic codes are formal. The modes are the 16 fields of one generation at one point, with no dynamics. Time reversal here is a toy antiunitary symmetry, not the Standard Model's CPT.
