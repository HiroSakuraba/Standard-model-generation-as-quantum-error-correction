# Correlated exchange–electric faults: timing, ideal recovery, and the physical-cost boundary

19 October 2026. Follow-up to merged [PR #63](https://github.com/HiroSakuraba/Standard-model-generation-as-quantum-error-correction/pull/63).

## Question

The finite non-Abelian S3 subgroup code in PR #63 protects four logical quantum states against every *individual* matter permutation in a six-element group, plus the link-electric phase R=I−2P_standard. The decoder has two four-dimensional syndrome branches and restores every admitted single error perfectly. But if R occurs and then a matter permutation g occurs before readout, four of six gR composites escape the intended syndrome space. Is that an unavoidable failure, or can the recovery timing prevent it? Does the extra recovery round cost more in particle loss than it saves in logical errors?

This experiment answers that question **inside the same 14D Gauss-constrained toy Hilbert space** for both six-element error groups <A,S> and <B,S>. The original 486D embedding, exact four-dimensional code, and electric-phase error are unchanged. This is a mathematical fault-order and resource sensitivity study, **not** a device-level quantum-memory implementation.

## 1. Fault order changes correctability

For either order-six subgroup G and PR #63's four-dimensional code C,

\[
gC=C\quad(\forall g\in G),\qquad
RC=-\tfrac13 C+\tfrac{\sqrt8}{3}V,\quad C^\dagger V=0.
\]

Here V is the four-dimensional orthogonal electric-syndrome isometry. The original recovery projects onto C and V, maps V→C, and flags any remaining probability outside those eight dimensions, **without postselection or renormalization**.

A composite error **gR** means the electric phase happened first and the exchange second; **Rg** means the exchange happened first. Thus \(RgC=RC\) remains perfectly recoverable in one round. But \(gRC\) need not lie in the eight-dimensional C+V syndrome space.

For both groups, the exact six-element enumeration (excluding none of the group elements) gives:

| Two-fault operator | Group elements | One-round unconditional F_e | Syndrome capture | Flagged outside syndrome |
|---|---:|---:|---:|---:|
| gR, group elements compatible with syndrome | 2 of 6 | 1 | 1 | 0 |
| gR, remaining elements | 4 of 6 | **11/96 ≈ 0.114583** | **11/24 ≈ 0.458333** | **13/24 ≈ 0.541667** |
| Rg, all group elements | 6 of 6 | **1** | 1 | 0 |

These are from explicitly represented **14×14 gauge-invariant operators** and the same original decoder. They are not generic claims about the order of quantum errors, and the gR failure is *not* part of PR #63's seven-single-fault theorem.

## 2. Why a second recovery round works ideally

Apply the original perfect recovery once after the first fault, then again after the second. The first admitted fault maps each syndrome branch to a **scalar multiple of I₄** on the logical space, so the logical state is restored independently of the classical syndrome. The second allowed fault then also acts on a freshly encoded state and is corrected.

The simulator composes all four branchwise logical Kraus products **without postselection** and checks for every gR and Rg pair that:

\[
F_e^{\text{two ideal rounds}}=1,\qquad
P_{\text{captured}}=1
\]

to numerical tolerance, for both subgroups.

This is a **time-separated single-fault correction protocol**. It assumes no additional errors occur during recovery, and it does not make the one-round code correct an arbitrary product of error operators. Whether a sufficiently fast physical syndrome/readout can occur between real microscopic events is unresolved.

## 3. Explicit gauge-preserving three-outcome coherent recovery—but globally supported

Let \(P_C=CC^\dagger\) and \(P_V=VV^\dagger\) on the 14D physical Gauss Hilbert space, and \(P_{\mathrm{rest}}=I_{14}-P_C-P_V\), of rank six. Use a gauge-neutral **three-level pointer** starting in \(|0\rangle\), not just a two-level pointer. This distinction matters because the compound errors can populate the six-dimensional orthogonal complement.

Define pointer swaps \(X_{01}\) between \(|0\rangle,|1\rangle\) and \(X_{02}\) between \(|0\rangle,|2\rangle\). The coherent syndrome operation is

\[
U_{\rm synd}=P_C\otimes I_3+
 P_V\otimes X_{01}+P_{\mathrm{rest}}\otimes X_{02}.
\]

Define \(F=CV^\dagger+VC^\dagger+P_{\mathrm{rest}}\) and apply feedback controlled on pointer 1:

\[
U_{\rm feed}=I_{14}\otimes|0\rangle\langle0|+
F\otimes|1\rangle\langle1|+
I_{14}\otimes|2\rangle\langle2|.
\]

Both are unitary on the entire **42-dimensional** physical-plus-pointer space. The tests verify three complete input routes: \(C|0\rangle\to C|0\rangle\); \(V|0\rangle\to C|1\rangle\); and any state in \(P_{\rm rest}\) goes to the same flagged physical state with pointer \(|2\rangle\), with no postselection. Because C and V lie inside the original full 486D simultaneous Gauss sector, their projector lifts \(W P W^\dagger\) commute with both vertex gauge actions.

**This is not a local pulse implementation.** Projecting onto the four-dimensional code, four-dimensional error syndrome and six-dimensional complement generally requires operators involving all four matter qutrits and the gauge link. Unitarity and gauge invariance do *not* imply few-body physical implementability. The number, duration, strength and failure probability of elementary charged-carrier operations required to realize this 42D unitary are unknown. In particular the illustrative 1 μs round below is not an estimate derived from this construction.

## 4. A normalized 20-event CPTP fault model with an order parameter

Use independently specified Bernoulli fault probabilities:

- Electric-phase event with probability \(p_e\).
- Matter permutation event with probability \(p_g\). Conditional on occurrence, select uniformly among **all six** subgroup elements, including the identity. This convention deliberately matches the exact six-element group average.
- Conditional on both events, fraction \(f\) has the harmful order gR (electric first), while fraction \(1-f\) has Rg.

There are **20 listed weighted unitary alternatives**: no error, electric only, six exchange-only, six gR, and six Rg. The probability sum equals 1 exactly in the model. Because each alternative is unitary, this mixture is completely positive and trace preserving. The simulator explicitly weights all alternatives and retains all flagged syndrome leakage.

Only four of six gR elements are harmful; their F_e is 11/96. The uniform group mean F_e after a gR is \(59/144\). Therefore, before particle loss,

\[
\boxed{F_e^{(1\ \rm round)}=1-\frac{85}{144}p_e p_g f}
\]

and the **two-round ideal protocol** has \(F_e^{(2\ \rm rounds)}=1\). These are exact rational consequences of the admitted channel and tested code, not a universal quantum error-correction threshold.

### Example: \(p_e=0.37\), \(p_g=0.40\)

| Harmful-order fraction f | One-round conditional F_e | Two-round ideal conditional F_e |
|---:|---:|---:|
| 0 | 1.000000 | 1.000000 |
| 0.5 | 0.956319 | 1.000000 |
| 1 | 0.912639 | 1.000000 |

The result is strongly sensitive to *time ordering*, not merely to total fault count. The input probabilities are illustrative and **are not measured rates** of gluons, nuclear matter, or any actual quantum device.

## 5. A mathematically grounded passive reference

For context, a four-state passive memory that is invariant under every g in a six-element subgroup can be selected inside its six-dimensional trivial sector. On that space, the electric phase R has exact real eigenvalues

\[
-1,-1,-\tfrac13,-\tfrac13,\tfrac13,1.
\]

The Ky Fan extremal trace principle says that the maximal absolute trace of R restricted to any four-dimensional subspace of that six-dimensional exchange-fixed sector is **8/3** (the four most negative eigenvalues). The best normalized logical trace therefore has magnitude \(2/3\), giving exact no-feedback entanglement fidelity

\[
F_{e,\mathrm{passive\ fixed}}(R)=\frac49.
\]

For any stochastic channel of the form in Section 4, exchanges have no effect on this G-fixed passive code, regardless of temporal order. Its optimum within that restricted class is \(F_e=1-\frac59 p_e\), before particle loss.

**Important:** this is an exact best result only among **four-dimensional subspaces inside the six-dimensional G-trivial sector**. It is *not* an optimization over arbitrary 14D passive codes, approximate noiseless memories, different physical inventories, or different noise strengths. The exact nonexistence of a passive 4D noiseless *subsystem* from PR #63 remains a separate error-algebra theorem.

## 6. Include lifetime and equal-clock-time overhead

In an explicitly hypothetical scenario, assume one idealized recovery round takes \(\tau=1\ \mu\text{s}\); two rounds take \(2\tau\). The code and the passive reference each carry the **same four active fermions**, with independent particle-loss lifetime \(T_1\), and loss of any one is a flagged erasure. Thus an exact no-jump factor is \(\eta=e^{-4\tau/T_1}\).

The single-round unconditional fidelity is \(\eta F_e^{(1)}\), and the two-round unconditional fidelity is \(\eta^2 F_e^{(2)}=\eta^2\). The passive G-fixed reference can be idled for \(\tau\) or \(2\tau\) to match whichever active scheme is evaluated, with the same \(\eta\) factors. No T2, ancilla loss, link error during recovery, preparation/reset/readout errors or physical pulse calibration is included.

**Two rounds beat one after loss iff**

\[
\eta>1-\frac{85}{144}p_e p_g f.
\]

For the illustrative \(p_e=.37,p_g=.40\) and one-microsecond rounds:

| f | One round after loss at T₁=30 μs | Two rounds after loss | Required T₁ for two rounds to win |
|---:|---:|---:|---:|
| 0 | 0.875173 | 0.765928 | No finite crossover |
| 0.5 | 0.836945 | 0.765928 | Above ~89.56 μs |
| 1 | 0.798717 | 0.765928 | Above ~43.76 μs |

This is a **resource threshold inside a synthetic noise model**. Real syndrome extraction might take considerably longer than one microsecond and introduce noise of its own; then the crossover shifts or disappears.

## 7. What this does and does not resolve

**Established mathematically:** with two distinct time windows and an ideal recovery in between, the new 4D finite-S3 code repairs all tested two-fault gR/Rg sequences; one recovery at the end fails for particular fault orders. The exact order-dependent Kraus mixture, exact group-fixed passive bound, and four-active-fermion lifetime crossover are independently testable.

**Not established physically:** an affordable native implementation, virtual-QCD corrections, microscopic correlated fault rates and waiting-time distributions, independent loss of spectator/reference/pointer/link, recovery during noisy fault events, or a direct nine-quark (21D/65-syndrome) implementation. The original PR #55–58 Hubbard experiments include virtual high-energy doublon–hole intermediates, but their pulse optimization and physical T1 estimates cannot be transferred automatically to this newly constructed four-state code.

**Next priority:** compile the 28D gauge-preserving global pointer instrument into specified local/charged-mediator primitives; determine its actual pulse duration and leakage; then use that measured or simulated duration, not 1 μs, to evaluate the fault-order crossover. A microscopic finite-S3 link-electric noise mechanism and repeated stochastic jump trajectories would determine whether such protective timing is realistic.

### Reproduce

\`\`\`bash
python -m pip install numpy scipy sympy
OPENBLAS_NUM_THREADS=1 python -m unittest -v test_colour_s3_correlated_recovery_cadence
OPENBLAS_NUM_THREADS=1 python colour_s3_correlated_recovery_cadence.py
\`\`\`

The dedicated GitHub workflow archives both groups, all 12 ordered two-fault cases per group, the explicit 42D three-outcome coherent syndrome unitary tests, 20-event weighted channels, exact passive bound and no-loss/lifetime comparisons in JSON.
