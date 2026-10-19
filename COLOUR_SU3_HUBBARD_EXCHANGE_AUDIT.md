# Microscopic virtual exchange: exact SU(3)-symmetric two-orbital Hubbard test

19 October 2026. Follow-up to PR #54, which specified but did **not** synthesize the Hamiltonian \(P_{\rm pointer}\otimes(I-\mathrm{SWAP})/2\) in a charged mediator/matter pair. This experiment supplies an underlying *two-fermion microscopic Hamiltonian* with occupation-dependent pointer control, and computes the finite-coupling error **without substituting the effective exchange for full dynamics**.

## Experimental Hilbert space and gauge boundary

There are **six fermion modes**: orbital L has colours \(\alpha=1,2,3\) and orbital R has the same three colours. With exactly two fermions, the Fock sector has \(\binom62=15\) states. Nine computational states have one fermion per orbital, and **six virtual doubly occupied / empty orbital states** are propagated explicitly.

**Both orbitals here must be interpreted as two co-located subsystems sharing ONE SU(3) gauge transformation at a common vertex.** The colour-blind fermionic hopping commutes with that common SU(3), as verified for all eight colour generators on the full 15D Fock space. It does *not* commute with independently specified SU(3) gauge transformations at distinct spatial vertices in the absence of a gauge link. This test is therefore a **microscopic onsite exchange module for PR #53's charged mediator**, not a full gauge-link-transport simulation or an independently gauged two-vertex Hamiltonian.

The fermionic signs in hopping and Fock occupancy are treated exactly; virtual doublons are neither projected out during pulses nor deemed 'corrected' when they remain at the end.

## Microscopic Hamiltonian, including ordinary pointer–density control

For two fermions in this six-mode register,

\[
H_p=-t\sum_{\alpha=1}^{3}
  (c^\dagger_{L\alpha}c_{R\alpha}+c^\dagger_{R\alpha}c_{L\alpha})
  +\frac U2\sum_{s=L,R} n_s(n_s-1)
  +\frac{\delta_p}{2}(n_L-n_R).
\]

The neutral **three-level pointer** changes the orbital *detuning* \(\delta_p\): pointer value 0 gives \(\delta_0=0\) and values 1/2 give \(\delta_1=\delta_2=\delta\). This is a colour-blind pointer-to-density energy coupling, rather than a pointer-controlled colour-exchange or Fredkin Hamiltonian. The onsite \(U\) is the usual two-fermion density interaction. The proposal consequently has a genuinely more microscopic **two-body interaction structure** at the fermion-mode/pointer–orbital level, **but a coherent pointer-conditioned detuning actuator and its strength/range are assumed, not experimentally established**.

The pointer is not used to pick a preferred colour axis, so the Hamiltonian commutes with common-vertex SU(3) for all three pointer blocks.

## Superexchange prediction and full finite-U result

Eliminate the six high-energy double-occupancy states to second order only to derive an *analytical reference*:

\[
H_{\rm eff}^{(2)}
 =-J(\delta)\frac{I-\mathrm{SWAP}}2,\qquad
J(\delta)=2t^2\left[\frac{1}{U+\delta}+\frac{1}{U-\delta}\right]
 =\frac{4t^2U}{U^2-\delta^2}.
\]

The program verifies the 9×9 matrix identity exactly to floating-point precision by independently calculating the second-order block contraction. The genuine gate and leakage results use the **exact exponential of the 15×15 Hubbard Hamiltonian**, not the 9×9 effective approximation.

With \(\delta_{\rm off}=0\), and \(\delta_{\rm on}=U/\sqrt3\), the effective exchange rates satisfy \(J_{\rm on}/J_{\rm off}=3/2\). A proposed single rectangular pulse of duration \(T=2\pi/J_{\rm off}\) gives an **off** antisymmetric phase \(2\pi\) (identity) and an **on** phase \(3\pi\) (a swap). The symmetric colour sector is dark. This is a coherent pointer-controlled SWAP proposal based on *detuning*, not a native 3-body controlled-exchange assumption.

The pulse is then optimized over precisely **two** variables (duration and on-pointer detuning) within a documented bounded window near the analytic starting point. There is no unrestricted waveform discovery or hidden tuning of independent pointer phases. The objective is the **leakage-inclusive** entanglement fidelity of the full 27-dimensional pointer × 9-colour target isometry \(\operatorname{diag}(I_9,\mathrm{SWAP}_9,\mathrm{SWAP}_9)\).

Crucially, fidelity is computed from the projection of **actual full 15D evolution** back to singly occupied subspace, without renormalization:

\[
F_e = \left|\frac{\mathrm{Tr}(M_{\rm off})+
                    2\,\mathrm{Tr}(\mathrm{SWAP}^\dagger M_{\rm on})}{27}
       \right|^2.
\]

The probability remaining in all six doublon/hole sectors is accounted for explicitly as final occupation leakage, and the maximum intermediate virtual occupation is sampled at 41 pulse times. Even when the final leakage is small, transient double occupancy may be material for decoherence and breakdown of a truncated computational-register model.

The program tests interaction ratios \(U/t=20,30,40\), reports both untuned and optimized fidelity, and *does not* call either one a hardware result.

## Illustrative time and particle loss

For context only, if the microscopic hopping is assumed to have \(t/h=10\) MHz, the dimensionless pulse duration \(T\) maps to \(1000T/(2\pi\times10)\) nanoseconds. This is an **assumption**, not a calibrated hopping amplitude or gate duration. The report also includes a labelled **independent-particle loss/erasure illustration** assuming \(T_1=30\,\mu s\) for *each of the two fermions*, with no-loss probability \(\exp[-2t_{\rm phys}/T_1]\). All lost-particle events are counted as erasures (never corrected), giving a simple survival-weighted fidelity figure.

That illustration is **not** a solution of the Lindblad master equation or a calibrated decoherence model. It does not include dephasing, virtual doublon decay, a physical gauge penalty, readout/feedback, spectator couplings, transport, or noise during pointer preparation. It also does not include the **twelve local quantum-controlled swaps plus link transports** used by the entire PR #53 syndrome, and does not assume multiplying one-gate fidelity twelve times is an accurate process calculation.

## Interpretation and follow-up

The real achievement would be an accurate **two-body microscopic route to one local controlled exchange**, with a physically plausible detuning actuator, while maintaining the common-vertex SU(3) symmetry. The 15D test can reveal its finite-U colour-operation errors and double-occupancy leakage. It is **not** equivalent to a full device-level quantum memory.

The next milestones should be (a) coherent physical pointer–site-detuning Hamiltonian, realistic finite ramps and native control bandwidth; (b) evolve the enlarged charge/occupation space through the **entire** charged-mediator syndrome, allowing doublon leakage to persist between pulses; (c) introduce consistent microscopic particle loss, pure dephasing and link/mediator decoherence, followed by a matched unencoded reference memory; (d) only then consider transferring any construction to continuous SU(3) link Hilbert spaces and the original nine-quark 21D/65-syndrome recovery problem.

Related primary sources on SU(N) Hubbard Mott and strong-coupling exchange include **Phys. Rev. A 104, 043316 (2021)**, DOI https://doi.org/10.1103/PhysRevA.104.043316, and the SU(N) Fermi-Hubbard Hamiltonian in *Synthetic gauge fields in synthetic dimensions: interactions and chiral edge modes*, *New Journal of Physics* **18**, 035010 (2016), DOI https://doi.org/10.1088/1367-2630/18/3/035010. Neither supports the specific device timing or practical memory performance assumed here.

## Reproduce

    python -m pip install numpy scipy
    OPENBLAS_NUM_THREADS=1 python -m unittest -v test_colour_su3_hubbard_virtual_exchange
    OPENBLAS_NUM_THREADS=1 python colour_su3_hubbard_virtual_exchange.py

GitHub Actions archives every scan point, leakage rate, exact optimised fidelity, and all claim boundaries as JSON.
