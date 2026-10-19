# Persistent virtual doublons across the complete charged-mediator syndrome

19 October 2026. Follow-up to PR #55 (microscopic controlled exchange) and PR #53 (charged mediator transported over an S3 gauge link). The key new test is **whole-circuit accumulation**, not another isolated controlled gate.

## Question and Hilbert space

PR #55 calculated one pointer-detuned two-orbital SU(3) Hubbard controlled swap in the exact **15D two-fermion** Fock sector. Applying that 15D map separately to each gate would discard the crucial history: virtual doublon/hole populations that survive a pulse can affect later operations.

We now retain **four active three-colour fermionic orbitals**—matter A0, matter A1, matter B0, and the mobile mediator M—with **four identical fermions in twelve modes**. The persistent active-sector Hilbert dimension is

\[\binom{12}{4}=495.\]

The full model also retains a spectator B1 colour qutrit, a six-state finite S3 gauge link, a two-position mediator-location register, a stationary reference colour qutrit paired with the charged mediator, and the three-level syndrome pointer. No Fock sectors are projected, reinitialized, normalized or discarded **between gates**. All occupation configurations—doublon/hole, triply occupied and even fully redistributed configurations allowed by four-fermion statistics—remain available.

The finite S3 link is *not* a continuous SU(3) gauge link. **Within one site/vertex**, each microscopic two-orbital hopping Hamiltonian is common-SU(3) symmetric. The S3 gauge transformations at the two vertices act on all occupied colours through the corresponding fermionic exterior-power representations, including leaked occupations. Charged-mediator link transport rotates *every occupied mediator colour* with the link and moves its location. This generalized transport is an exact **assumed gauge-covariant unitary**, not a microscopic hopping implementation.

## Persistent local dynamics

For each of the twelve pointer-controlled exchanges required by PR #53, evolve the **same shared four-fermion state** with

\[H_p=-t\sum_{a=1}^3(c^\dagger_{M,a}c_{q,a}+c^\dagger_{q,a}c_{M,a})+\tfrac U2[n_M(n_M-1)+n_q(n_q-1)]+\tfrac{\delta_p}{2}(n_q-n_M).\]

Here q is A0, A1 or B0, *only if the charged mediator is at that vertex*. The neutral pointer chooses detuning \(\delta_p=0\) for inactive labels and \(\delta_p=\delta\) for active labels. At the mathematical Hubbard level the pointer couples to site-number imbalance, not to an assumed Fredkin Hamiltonian. The four link transports from PR #53 are applied as exact unitary transformations of the **entire four-fermion Fock space**, including virtual occupations of M. The charged reference singlet is kept throughout.

The active Fock Hamiltonian is 495×495, but spectator orbital occupation numbers are conserved by a pair pulse. The code diagonalizes independent blocks of size **at most 20** and assembles the exact sparse propagator, rather than truncating to 9D colour code space. Fermionic anticommutation signs are computed from the actual twelve-mode ordering. The resulting computation is modest enough to test all sixteen stages at full state-vector resolution.

## Error channel and metrics

The input is the **four-dimensional finite S3 gauge-invariant eigenphase code** from PR #48. The tested source error is a **known dressed-swap** channel with probabilities 0.63 (no swap) and 0.37 (swap). Perform the complete neutral-pointer three-cycle syndrome: **twelve Hubbard pulses plus four mediator-link transports**, followed by the qutrit Fourier measurement. Measure final syndrome probabilities without postselecting.

The readout analysis computes **unconditional entanglement fidelity** by tracing all pointer, mobile-mediator location/colour, stationary reference and leaked Fock-sector environments. Only amplitude returned to one fermion per active orbital can overlap the original code. Final occupation leakage is reported separately, as is the *largest intermediate* leakage. There is **no per-gate success conditioning**.

**Correction boundary:** syndrome-2 feedback is currently an **ideal dressed colour swap on the computational occupation sector**. No microscopic Hubbard implementation of correction is included, and leaked sectors are **not recovered**. This is a useful *ideal-feedback reference* for the syndrome, **not** a general upper bound on all possible optimized decoders or a full device recovery channel.

Test the same optimized single-pulse parameters from PR #55 for U/t=20, 30, 40. The whole-circuit fidelity must be independently measured; it is invalid to extrapolate directly from one gate's fidelity. Even at a moderate local interaction ratio, small transient doublons can accumulate, generate state-dependent phases, fail to return at the end, or distort the three-cycle pointer outcomes. The Gauss symmetry can stay exact at the same time: **gauge invariance is not equivalent to correctability**.

## Separate illustrative loss-time cost

As in PR #55, take *illustrative only* t/h=10 MHz and assume each transport takes 100 ns; then syndrome time is 12 times the microscopic optimized gate duration plus 4 transport durations, **excluding** the reference preparation, pointer Fourier, measurement, and feedback. With independent active-fermion T1=30 microseconds and four active fermions, quote the *hypothetical no-loss factor*

\[S=\exp(-4T_\mathrm{syndrome}/T_1)\]

and its product with the ideal-feedback reference fidelity. These are labelled sensitivity figures, **not** a Lindblad calculation or experimental data. They exclude spectator and reference loss, dephasing, doubled-occupation-specific decay, noise during link transport, and all syndrome readout and feedback delays. No architecture win can be inferred from them.

## Tests, limitations and next experiment

Dedicated tests cover the complete 495D four-fermion basis, unitary pair propagators, reproduction of PR #55's 15D two-orbital projected unitary in the single-occupation spectator sector, exterior-power S3 transport including occupation leakage, exact extended local Gauss invariance, and the non-postselected complete logical syndrome output. In particular, the leakage fraction must **not reset when transported**. CI archives all computed Fock leakage and fidelity data.

Main limitations: finite S3 link rather than full SU3 lattice gauge theory; only the four-dimensional toy code, **not** the nine-quark 21-dimensional 65-syndrome memory; idealized gauge-link transport, charge-reference singlet, feedback and pointer Fourier; rectangular Hubbard pulses with optimizations imported from the isolated-gate benchmark; no microscopic pointer actuator calibration, finite ramps, particle loss/dephasing quantum trajectories, hardware topology, or verified advantage over a matched unencoded memory.

The next experimentally relevant program step is replacing the **ideal feedback, link-transport unitary, and stochastic survival factor** with the same microscopic occupation and noise model used for syndrome extraction, then measuring recovery against a duration-matched unencoded storage baseline. In parallel, consider a continuous SU3 link formulation. The current numerical experiment is an explicit intermediate result, not the end of physical realization.

### Reproduce

```bash
python -m pip install numpy scipy
OPENBLAS_NUM_THREADS=1 python -m unittest -v test_colour_su3_hubbard_persistent_syndrome
OPENBLAS_NUM_THREADS=1 python colour_su3_hubbard_persistent_syndrome.py
```