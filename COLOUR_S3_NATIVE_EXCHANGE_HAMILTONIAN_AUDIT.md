# From an abstract charged-bus Fredkin to a specified SU(3)-isotropic exchange Hamiltonian

19 October 2026. Follow-up to merged PR #53. This models the required **pointer-controlled local matter–mediator exchange** as a concrete gauge-invariant Hamiltonian and checks its control requirements. It does not claim an experimental realization.

## 1. Exact generator and quantum control

For two co-located SU(3) fundamental qutrits, the mobile charged mediator \(m\) and matter qutrit \(q\), let \(T^a=\lambda^a/2\) be the eight Gell-Mann generators, normalized \(\operatorname{Tr}(T^aT^b)=\delta_{ab}/2\). The exchange identity is

\[
\sum_{a=1}^8 T_m^a\otimes T_q^a =
\frac12\left(\mathrm{SWAP}_{mq}-\frac13 I_9\right).
\]

The positive-semidefinite antisymmetric-sector projector

\[
H_{\mathrm{ex}} = \frac13 I_9 - \sum_a T_m^a\otimes T_q^a
                = \frac12 (I_9-\mathrm{SWAP}_{mq})
\]

is a **fully SU(3)-invariant two-qutrit Hamiltonian**. Consequently \(e^{-i\pi H_{\mathrm{ex}}}=\mathrm{SWAP}_{mq}\). This is a rigorous algebraic identity, not an optimization result. Its sign and scalar offset are essential for obtaining precisely SWAP, rather than an uncontrolled relative phase between pointer branches.

For a neutral qutrit pointer and selected active labels \(A\subset\{0,1,2\}\), the desired coherent control is

\[
H_{\mathrm{ctrl}}=
\left(\sum_{p\in A}|p\rangle\langle p|\right)\otimes H_{\mathrm{ex}},
\quad
e^{-i\pi H_{\mathrm{ctrl}}}=
\text{pointer-controlled-SWAP}_{mq}.
\]

The operator has **three-register support (pointer, mediator colour, target matter colour)** and is strictly invariant under the shared SU(3) colour gauge action at every pulse time. It therefore reproduces the high-level three-register local exchanges from PR #53 without assuming a black-box Fredkin matrix. The *location*-dependent charge transport on the link remains the separately assumed gauge-covariant **two-register** unitary of PR #53.

## 2. Why strictly invariant two-register native Hamiltonians cannot synthesize this control under full SU(3)

This obstruction is more general and more precise than an unsuccessful pulse search. The fundamental \(3\) is an **irreducible representation** of continuous SU(3). By Schur's lemma, an SU(3)-invariant operator on a gauge-neutral pointer \(p\) and one fundamental colour qutrit \(q\) has the form \(A_p\otimes I_q\); its complex commutant dimension is **9** for a three-state pointer. For two fundamental colour qutrits \(m,q\), the invariant algebra is \(\mathrm{span}\{I,\mathrm{SWAP}\}\), dimension **2**.

If the only available operations are SU(3)-invariant **one- and two-register** Hamiltonians on the three existing registers \((p,m,q)\), the pointer never couples to the coloured pair: the whole allowed Lie algebra factorizes between pointer-only operators and SU(3)-invariant pair exchange. Commutators cannot generate \(P_p\otimes\mathrm{SWAP}_{mq}\), so no gate count suffices.

The new code verifies the exact decomposition and independently computes both 9×9 commutant dimensions by numerical null spaces. The proof is Schur/Schur–Weyl, **not** double-precision singular-value calculations. The no-go assumes no extra degrees of freedom, gauge-neutral pointer, a *full* SU(3) irreducible fundamental, and every native Hamiltonian commuting with the gauge action at each instant.

**Careful exception:** The simulator's actual finite signed-permutation \(S_3\) subgroup acts **reducibly** on a three-dimensional colour qutrit (sign ⊕ standard). Consequently this exact *SU(3)* two-register no-go does **not automatically apply** to the finite S3 charged-bus model, where extra pair couplings could exist. The extended gauge-link restrictions of PR #50 and charged-mediator escape of PR #53 have different assumptions. Additional charged mediators, encoded gauge multiplicity pointers, gauge-breaking virtual pulses, or an explicitly engineered three-register exchange controller likewise fall outside the obstruction.

## 3. Full-system pulse and anisotropy tests

The code directly exponentiates the 9×9 local colour Hamiltonian (no logical projector inserted) and embeds its action in the **full PR #53 register**: neutral qutrit pointer, six-level gauge link, four matter qutrits, six-level mobile charged mediator, stationary conjugate-charge reference, and four logical-reference columns. It checks exact equality to the previous controlled local SWAP for arbitrary off-Gauss inputs, strict covariance at arbitrary pulse angles and every local SU3 generator, the full 16-operation logical pointer circuit, and exact *unconditional* recovery of the known dressed-swap channel after feedback.

A second, explicitly synthetic sensitivity scenario adds a colour-anisotropic exchange perturbation

\[
H_{\mathrm{local}}=H_{\mathrm{ex}}+\epsilon\, T_m^3\otimes T_q^3
\]

during each selected local exchange pulse and during the feedback. Here \(T^3=\mathrm{diag}(1,-1,0)/2\). This coupling **breaks full SU(3) and the selected finite S3 gauge invariance**; the code measures its commutator residuals with all eight infinitesimal SU3 generators and all six finite S3 matrices. It then traces all pointer, mobile-mediator location/colour, and reference outcomes to estimate unconditional recovered entanglement fidelity at \(\epsilon=0,0.01,0.05\).

The perturbation is coherent colour-axis miscalibration, **not a measured physical stochastic error rate**. The zero-anisotropy result must agree with PR #53's exact recovered fidelity one; the nonzero results assess sensitivity but are not a physical error threshold.

## 4. Illustrative timing, not a demonstrated interaction

A scenario mapping \(H_{\mathrm{ctrl}}/h=J_{\mathrm{eff}}\, P_p\otimes H_{\mathrm{ex}}\), with hypothetical effective three-register coupling \(J_{\mathrm{eff}}=5\) MHz, gives a conditional local SWAP pulse of \(1/(2J)=100\) ns. Assuming equally hypothetical two-register charged-link transport operations of 100 ns, the 12 controlled onsite swaps plus four transports require **1.6 μs** of serial syndrome time; the five-operation feedback adds **0.5 μs**, before pointer preparation, Fourier, readout, reference-charge preparation, movement, and reset.

**These numbers are explicitly scenario inputs, not measurements, engineering feasibility, hardware lower bounds, or physical predictions.** The interaction \(P_p\otimes \sum_a T_m^aT_q^a\) is itself a three-register conditional exchange. No native two-body device implementation or effective coupling \(J_{\mathrm{eff}}\) has been derived or calibrated.

## 5. Decision

The new work identifies a precise native-interaction specification: the hardware must realize a **quantum-conditional, colour-isotropic exchange amplitude** without selecting a preferred colour direction. For a full SU3 fundamental, a neutral pointer cannot provide such control using only strictly invariant two-register couplings on those three registers. A device must therefore supply a genuine effective multibody term, additional gauge-charged encoded ancillary degrees of freedom, or a carefully characterized temporary gauge-violating gadget.

The **next legitimate hardware milestone** is deriving \(J_{\mathrm{eff}}\), fidelity and cost from a particular microscopic Hamiltonian (e.g. a specified microscopic Hubbard-type model with extra charged ancillary modes), not assuming a free three-body control. Then include full physical noise, finite transport link errors, charge reference preparation and comparison with an unencoded storage benchmark. This is still finite-S3 **four-dimensional** logical code research, **not** the nine-quark SU3 21D/65-syndrome memory.

## Reproduction

    python -m pip install numpy scipy
    OPENBLAS_NUM_THREADS=1 python -m unittest -v test_colour_s3_native_exchange_hamiltonian
    OPENBLAS_NUM_THREADS=1 python colour_s3_native_exchange_hamiltonian.py

The dedicated CI workflow archives an audit JSON of algebraic identities, numerical commutants, exact circuit residuals, Gauss acceptance, anisotropy/logical-fidelity samples, and all assumptions.
