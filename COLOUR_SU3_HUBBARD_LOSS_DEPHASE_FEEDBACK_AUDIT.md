# Microscopic feedback, actual number-dephasing channels, and flagged fermion loss

19 October 2026. This is the next experiment after merged PR #55 (the exact 15-state two-fermion Hubbard gate) and PR #56 (persistent leakage through the full 495-state charged-mediator syndrome).

## Question

PR #56 applied twelve microscopic Hubbard colour-exchange pulses and four link transports in the complete charged-bus syndrome, retaining all occupation leakage. However, syndrome-dependent **correction** was still an ideal dressed colour SWAP confined to the computational sector. It also applied an external survival factor instead of a dynamical environmental-noise channel.

This PR makes two changes: **microscopic Hubbard correction** using the *same* 495D persistent Fock state, and an explicitly defined **completely positive, trace-preserving (CPTP)** channel for colour-mode or orbital-number dephasing and fermion loss.

The result remains a finite-group, gauge-inspired four-dimensional **toy code**, not a hardware proof or the original nine-quark 21D code.

## 1. Actual microscopic conditional feedback

After the twelve Hubbard exchanges and four ideal link transports of PR #56, the coherent qutrit Fourier measurement produces labels 0, 1 and 2. On label 2, correct by the charge-bus sequence

    Hubbard(A0) -> transport -> Hubbard(B0) -> transport -> Hubbard(A0).

The three local Hubbard operations use **exactly the same finite-U, pointer-dependent detuning and pulse duration** as the syndrome. Importantly, pointer labels 0 and 1 still undergo the corresponding **off-detuned Hubbard evolution**; they are NOT assigned an ideal identity gate. Both link transports act on all pointer labels, cancel ideally on the inactive branches, and retain every doubled/empty orbital configuration. The program never resets a leaked occupation sector between syndrome or correction operations.

This makes **15 microscopic Hubbard pulses and six link transports** across the complete syndrome and feedback, instead of twelve and four with an ideal decoder.

The logical entanglement fidelity is computed by summing the complete (unpostselected) code overlap after tracing every syndrome-pointer label, mediator location/colour, and stationary reference. Only physical amplitudes with exactly one particle per active orbital can overlap the input code, but those amplitudes are **never renormalized** to hide leakage.

## 2. A genuine stochastic dephasing channel at intermediate times

During Hubbard pulses, each pulse is divided into \(n\) finite-time intervals (default two). The exact 495D occupation-space Hubbard propagator is applied for each interval, followed by independent Gaussian phase rotations

\[
U_{\theta}=\exp\left[-i\sum_j\theta_j N_j\right],\qquad
\theta_j\sim {\cal N}(0,2\gamma \Delta t).
\]

After averaging over these random-unitary trajectories, the channel obeys the Lindblad equation with Hermitian jump operators \(L_j=\sqrt{2\gamma}\,N_j\). This is CPTP by construction and is **not** merely an end-of-circuit fidelity penalty. The phase environment is also sampled during idealized charge-link transports and a hypothetical pointer measurement delay, using the same dimensionless rates.

Two noise couplings are compared:

- **Orbital-charge dephasing:** \(N_j=n_{j,1}+n_{j,2}+n_{j,3}\), for each of the four active orbitals. Because it counts total occupation independent of colour, each noise realization is invariant under the extended finite-S3 gauge action, and under common-vertex SU3 colour rotations.
- **Colour-resolved dephasing:** \(N_j=n_{j,a}\) separately for all twelve colour modes. A realization generally selects privileged colour axes and can violate the vertex Gauss constraints, even if the *ideal Hubbard Hamiltonian* is gauge invariant.

The experiment records final extended-Gauss acceptance for the resulting noisy states. Random-unitary quantum trajectories estimate conditional-on-no-loss fidelity, with a reported Monte Carlo standard error. Finite subdivision and sampling uncertainty remain; there is no claim that a few trajectories constitute a high-precision Lindblad integration.

## 3. A loss channel, not an unexplained scalar adjustment

The microscopic active register has **exactly four identical fermions** in twelve colour modes. Suppose every active mode has the same one-particle loss rate \(\kappa=1/T_1\), with annihilation jump \(L_c=\sqrt\kappa\,c_c\). On the fixed-four-particle sector,

\[
\sum_c L_c^\dagger L_c=\kappa\sum_c n_c=4\kappa I.
\]

Hence the no-jump trajectory has a **state-independent, exactly known** survival probability

\[
\eta=\exp[-4\kappa T_{\mathrm{total}}],
\]

even when the four fermions become redistributed among doublon/hole configurations. Any actual particle loss changes the total number below four. Since the Hubbard exchanges and charge transports conserve particle number, no jumped state can return to the original fixed-four-fermion code without an extra particle-injection operation.

We therefore use a **CPTP flagged erasure instrument**: with probability \(\eta\), propagate the complete four-fermion channel, and with probability \(1-\eta\), output an orthogonal flagged missing-particle sector with **zero encoded fidelity**. The unconditional logical fidelity is \(\eta\) times the (not postselected) full four-particle dephasing/feedback fidelity. This branch decomposition is exactly justified by the uniform-loss Lindblad jump operators and fixed particle number; it is not a claim that the 495D model also propagates all explicit three-, two- and one-particle post-jump states.

**Scope:** unequal colour-dependent losses, loss-induced repumping, recovery after erasure, photon-assisted virtual-state decay, and microscopic post-jump dynamics remain absent. As in earlier experiments, the four active fermions are the only modes assigned loss. The spectator colour B1, reference qutrit, pointer, and finite group link have ideal survival; this is an optimistic assumption.

## 4. Fair comparison and reproducibility

The nominal matched-time scenario uses the same hypothetical values as PRs #55–56: \(t/h=10\) MHz; ideal charged-link transports of 100 ns each; a 300 ns pointer measurement/readout delay; and \(T_1=30\,\mu\mathrm{s}\) for each of four active fermions. **Total time includes feedback**, so the calculation covers fifteen Hubbard pulses, six link transports, and the readout delay. This is substantially longer than the PR #56 syndrome-only loss factor.

The report includes:

- U/t = 20, 30, 40 with **microscopic correction and zero dephasing**, compared against the corresponding PR #56 **ideal-feedback** no-noise reference.
- At U/t = 30, a small seeded noisy trajectory set for **orbital-charge** versus **colour-resolved** Lindblad-number dephasing, both at the same per-Lindblad-operator rate.
- The full no-jump logical fidelity, exact loss-erasure probability, **unconditional** loss-and-dephasing recovered fidelity, final occupation leakage, gauge acceptance, pulse count and hypothetical total physical duration.

The numbers from this synthetic model must never be called *measured hardware fidelity*. Unequal channel counts (four total-orbital operators versus twelve mode operators) also mean equal per-operator \(\gamma\) is **not equal total bath power**; the report makes no ranking of environments on that basis.

## 5. Remaining physical and mathematical limitations

- The local 495D Hubbard exchanges conserve colour gauge symmetry at a shared vertex. Link transport is still an **assumed** exact finite-S3 gauge-covariant register operation, not microscopic gauge-link hopping.
- The pointer-detuning actuator, hopping strengths, readout delay, link duration and active-particle lifetime are hypothetical.
- The stochastic dephasing environment is diagonal in occupation number. No amplitude noise, coherent hopping calibration drift, virtual-state-specific relaxation, coloured loss, measurement assignment error, or non-Markovian environment is modeled.
- The erasure-flagged loss channel is fully trace preserving but has **no explicit N=3 post-loss propagation**. All such events are unrecoverable for this code in the absence of repumping.
- End-to-end ideal syndrome readout remains a high-level Fourier transformation and ideal measurement.
- No duration-matched unencoded **reference memory** is implemented here, so there is no demonstration of physical advantage, threshold, or practically useful gauge-protected quantum error correction.
- The code corrects one *known* dressed-swap error in a 4D **finite-S3 toy** encoding. It is neither the full SU3 non-Abelian lattice with physical Standard Model fermions nor the original 21D nine-quark/65-syndrome problem.

## Verification

\`\`\`bash
python -m pip install numpy scipy
OPENBLAS_NUM_THREADS=1 python -m unittest -v test_colour_su3_hubbard_loss_dephase_feedback
OPENBLAS_NUM_THREADS=1 python colour_su3_hubbard_loss_dephase_feedback.py
\`\`\`

Eight dedicated regressions check occupation conservation, genuine stochastic unitary dephasing, gauge-invariant orbital noise versus colour-breaking mode noise, the exact loss-channel probability decomposition, actual Hubbard corrective gates and pointer branch trace conservation, pulse subdivision consistency, and seeded reproducibility.

**Next meaningful evaluation:** a full microscopic Lindblad solver with explicit post-jump sectors and a matched-duration **unencoded storage baseline** (including the same charged-reference and pointer preparation if required). Only then should physical viability and overhead of this finite-S3 scheme be judged.
