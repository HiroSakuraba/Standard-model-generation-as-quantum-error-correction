# Confinement as a code: nine quarks protect themselves

**Status, 19 October 2026:** exploratory note, not peer reviewed. Verified by `python3 colour_singlet_code.py` (numpy and scipy; about a minute; ends with `ALL COLOUR-SINGLET CODE CHECKS PASSED`).

**Question.** Confinement says physical states are colour singlets. Treat each quark's colour as a qutrit and the singlet subspace as a code. What does it protect against?

**Answer.**
- Every colour-singlet space detects any error on a single quark.
- None of them corrects one, because quark exchange acts non-trivially inside the singlet space.
- The exception is nine quarks. Their singlet space contains a 21-dimensional code that corrects any single-quark colour error, or any two lost quarks. It is fixed by the Pauli principle for nine quarks filling three colours and three flavours.

## Singlet spaces

The colour-singlet subspace of n quarks has the dimension of the number of standard Young tableaux of shape (n/3, n/3, n/3):

| Quarks | Singlet dimension | Example |
|---|---|---|
| 3 | 1 | one baryon |
| 6 | 5 | "hidden colour" in the deuteron channel |
| 9 | 42 | three baryons |

The script builds each space from products of ε tensors and checks that it is killed by all eight total colour generators.

## Detection: distance at least 2, for every n

Any operator E on one quark decomposes into a multiple of the identity plus an octet. An octet acting on a singlet gives a state in the octet, orthogonal to every singlet. So

  ⟨s| E_i |s′⟩ = (tr E / 3) ⟨s|s′⟩

for all singlets s, s′. Every single-quark colour error is detected, and one erased quark can be restored. The script checks this with random complex 3 × 3 errors on every site for n = 6 and 9.

## The obstruction: exchange

Correcting single errors needs the Knill–Laflamme conditions for products E_i† F_j on two different quarks. Projected onto singlets, the only invariant two-quark operators are the identity and the exchange P_ij, so the condition becomes

  Π P_ij Π ∝ Π   for every pair i < j.

By Schur–Weyl duality, the singlet space is an irreducible module of the permutation group S_n. P_ij has eigenvalues ±1 on it, so the full singlet space never corrects one error. Its distance is exactly 2. This is the SU(3) version of the obstruction met for SU(2) singlets in [`SU2_EXCHANGE_HAMILTONIAN_AND_NO_GO.md`](SU2_EXCHANGE_HAMILTONIAN_AND_NO_GO.md).

For six quarks the exchanges span a 9-dimensional space of traceless operators on the 5-dimensional singlet space. A numerical search found no 2-dimensional subcode on which they all compress to scalars: across 20 random starts, the smallest residual was 0.621.

## Nine quarks: the sign twist

The singlet shape for nine quarks is the 3 × 3 square (3,3,3), and the square is the only self-conjugate shape of the form (m, m, m). When a shape equals its transpose, tensoring the S_n module with the sign representation returns the same module. So there is an operator J with

  J P_ij = −P_ij J   for every exchange.

The script finds J in two independent ways:

1. **Linear algebra.** The equations J s_k + s_k J = 0 for the eight adjacent exchanges have a one-dimensional solution space.
2. **Pauli principle.** Put nine quarks in the nine states colour ⊗ flavour (3 colours, 3 flavours, one spatial orbital). The only fully antisymmetric state is the filled shell, and it is a colour singlet and a flavour singlet. Its amplitudes in the singlet ⊗ singlet basis define a matrix J_mn. Antisymmetry under swapping two quarks, colour and flavour together, is exactly J P_ij = −P_ij J.

The two agree up to normalization. J is real, orthogonal and antisymmetric, so J² = −1.

**The code.** Let K be the +i eigenspace of J. It has dimension 21, half the singlet space. For u, v in K,

  ⟨u|P_ij|v⟩ = −⟨u|J P_ij J^{−1}|v⟩ = −⟨u|P_ij|v⟩ = 0.

Every exchange has zero matrix elements on K, so all 36 Knill–Laflamme conditions hold with constant 0.

**Full check in 3⁹ dimensions.** The script embeds K in the 19,683-dimensional space of nine qutrits. For each site and each pair of sites, it computes the reduced overlap tensor, which covers all 9 one-site and all 81 two-site operator directions at once. Every one compresses to a scalar on K, with maximum deviation 1.1 × 10⁻¹⁵. So K corrects any single-quark colour error, and any two erased quarks.

**Distance exactly 3.** A 3-cycle of quarks commutes with J and acts non-trivially on K, so some three-quark error is undetectable.

**Colour exchange is invisible.** One-gluon exchange between quarks i and j is T_i·T_j = ½(P_ij − ⅓). On K, every interaction Σ c_ij T_i·T_j is a constant, for any couplings c_ij. The script checks this with random c_ij.

**Size.** A search for a 22-dimensional correcting subcode failed (residual 0.864), so 21 appears to be the largest. In qutrit units this is [[9, log₃ 21 ≈ 2.77, 3]]₃. The quantum Singleton bound allows 5 qutrits, so the code is not efficient. What it offers is symmetry: it is invariant under every global colour rotation and also corrects local errors.

## Why nine, and why colour

For SU(N), N² particles form the self-conjugate square:

| Group | Particles | Singlet dimension | Half |
|---|---|---|---|
| SU(2) | 4 | 2 | 1 |
| SU(3) | 9 | 42 | 21 |
| SU(4) | 16 | 24,024 | 12,012 |

For SU(2) the sign-twist code is one-dimensional and stores nothing. SU(3) is the smallest gauge group where the construction protects information. For three colours, nine quarks is the only quark number where it works.

## Scope

- **Errors are model errors.** A single-quark colour rotation is not gauge invariant; here it stands for an arbitrary error on a qutrit that simulates one quark's colour. Spin, flavour and position are not part of the code.
- **Literature.** Codes inside single SU(3) irreps exist ([arXiv:2312.00162](https://arxiv.org/abs/2312.00162)). A recent framework ([arXiv:2604.05871](https://arxiv.org/abs/2604.05871)) gets Knill–Laflamme conditions from SU(d) decoupling groups. Neither builds this sign-twist code. A full survey of SU(d)-invariant codes has not been done.
- **Maximality.** "21 is the largest" is a numerical statement, not a proof.
