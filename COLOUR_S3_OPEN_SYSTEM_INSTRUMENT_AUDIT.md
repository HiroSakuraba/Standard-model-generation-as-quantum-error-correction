# Finite-S3 recovery under four-level pointer noise: matched-time audit

**19 October 2026 — following merged PR #70. Status: executable
phenomenological open-system experiment, not a calibrated physical
decoder.**

## Why this experiment is different

PR #70 produced explicit 104-control products for both finite-S3 codes.
The weakest coherent branch reached ~0.578 (AS) and ~0.556 (BS), well
below the 0.99 mathematical three-outcome acceptance threshold. Its
projected seven-fault logical fidelity remained below the zero-control
0.671111 baseline, *before* physical decoherence.

This PR asks a narrower and falsifiable question: **how do those
unpostselected pulse products compare to passive storage when an
independently noisy four-level pointer, gauge-compatible electric
dephasing, measurement overhead and four-particle erasure are charged
over the same elapsed time?** It does not treat CP22 as an available
finite-S3 gate or confuse the Gauss-reduced 14D basis with a physical
14-level transmon.

## Actual operations simulated

The 14D gauge-physical system and its three-level neutral pointer
(original 42D coherent controls) are embedded into a 14×4 = 56
dimensional Hilbert space, where level `|3>` is an explicit **additional
pointer leakage state**. The existing compiled 42D Hamiltonian
unitaries act on pointer levels 0,1,2; they are extended as identity
on pointer 3. This extension is a declared toy model, *not* an actual
device Hamiltonian.

The compiler's pulse product is split into two successive halves.
After the first half and after the second half, we apply the complete
four-level bosonic amplitude-damping Kraus channel already verified in
`su3_transmon_noise_probe.py`, followed by its Gaussian number-phase
Schur-Kraus dephasing channel. No outcomes are sampled or thrown away:
**the program sums the entire Kraus tree**, keeping the logical
reference index and individual recovery outcomes.

The noisy windows have times `t_pulse/2` and
`t_pulse/2 + t_readout + t_reset`. Pointer leakage arises from an
independent coherent `|2> <-> |3>` rotation at each window boundary.
This is **not** a physically calibrated CP22 `|22> <-> |31>`
transition. The associated fourth-level population is measured after
noise and included in the unpostselected fidelity.

The 14D system has no tensor-factorization into four-level devices.
Instead we model **a collective gauge-compatible electric-phase
dephasing channel** with two Kraus operators

`K_0 = sqrt((1+exp(-t/Tphi_E))/2) I`,
`K_1 = sqrt((1-exp(-t/Tphi_E))/2) R`,

where `R=I-2P_E` is the exact 14D electric-reflection fault.
This is a fully trace-preserving mixed-unitary channel because
`R†R=I`; its rate is a phenomenological coupling, **not**
microscopic per-quark/local-transmon dephasing.

Four independent active fermions are also assigned a **flagged
particle-loss probability** `1-exp(-4t/T1)`. A loss is an orthogonal
eraser and contributes zero logical entanglement fidelity, so the
reported scalar becomes `F_e_full = exp(-4t/T1) F_e_surviving`.
This is *not* a bosonic ladder amplitude-damping channel on the
fermions, nor a simulated post-loss gauge state. Pointer damping is
the actual four-level amplitude-damping channel; active-particle
loss is explicitly a separate erasure approximation.

## Fair negative controls

Every input includes the same previously defined synthetic
single-fault mixture: electric reflection with probability 0.37;
six distinct physical exchange errors jointly carry weight 0.63.
Because all six exchanges fix the input code exactly, two input blocks
represent the complete **restricted one-fault** channel.
Compound errors, repeated syndrome rounds and hardware cross-talk
remain outside this experiment.

The same unrecovered logical entanglement fidelity is computed for:

1. The actual PR #70 compiled pulse sequence (including modeled
   pointer losses, confusion and readout/reset times).
2. **Same-elapsed-time passive storage:** no syndrome operations,
   no invented pointer readout or resets, but **the same total elapsed
   time** under the same collective electric dephasing, pointer
   relaxation assumption and active-particle erasure.
3. The **zero-elapsed-time no-gate reference** at the same initial
   fault distribution (0.671111 before noise).

For the candidate, a symmetric three-label readout confusion matrix
with probability `p` of being assigned an incorrect label counts
whether a pointer outcome is accepted as code/electric; fourth-level
leakage is always rejected. Code-projector overlaps appear **only as
measured diagnostics of the logical output**. There is **no physical
code projection or corrective gate inserted for free**.

We label the result a *matched-time noise proxy*, not a
device-matched experimental advantage or a truly compiled measurement
network.

## Device-specific CP22 and analog refocusing

This repository already includes independent 16D Duffing/exchange
pulse studies with fourth-level transmon noise
(`su3_transmon_noise_probe.py`) and a five-ququart (4^5 = 1024)
quantum-trajectory study
(`su3_five_ququart_trajectories.py`) for the separate
**five-plaquette digital SU(3)** construction. They are not silently
reused as charged-bus controls in this *finite-S3* model.

To honestly insert a device-specific CP22 into this finite-S3
recovery, we would first need a **physical register encoding** for
its 14 Gauss states, reference, charged mediator and gauge link,
and a decomposition of each of the 13 42D control Hamiltonians into
native drive channels. Neither encoding/decomposition exists yet.
Importing a calibrated gate error probability without that mapping
would be physically invalid. Likewise an analog comparison requires
a common gauge-simulator Hamiltonian, control amplitudes and storage
horizon; it is not manufactured by multiplying fidelity by a chosen
exponential.

## Usage

```bash
python -m pip install numpy scipy sympy
OPENBLAS_NUM_THREADS=1 python -m unittest -v test_colour_s3_open_system_instrument
OPENBLAS_NUM_THREADS=1 python colour_s3_open_system_instrument.py
```

The benchmark exactly reproduces the 13→26→52→104-control
continuation for AS and BS, then evaluates the selected 104-control
products under two illustrative noise scenarios (with and without
coherent pointer level-3 mixing). The JSON artifact includes
unconditional Kraus-summed logical scores, per-fault and
per-pointer probability ledgers, population in level 3, all
assumed times, active-particle erasure survival, passive and
zero-time controls.

The dedicated regression tests cover pointer-four extension,
coherent leakage unitarity, trace-preserving noise,
noiseless perfect-global target, the 0.671111 no-gate baseline,
real pulse products, readout assignment and **exact equal-time
comparison**, and parameter rejection.

## Interpretation boundaries

- Even a coherent global ideal instrument is an **unsynthesized
  mathematical upper bound**; only the actual numerical pulse products
  count for circuit-level claims.
- The two-window open-system simulation is not Lindblad evolution
  during each gate. Actual 104-control schedules have strongly
  nonuniform durations; equal-window approximations may matter.
- A CPTP amplitude damping channel on the pointer is **not** a
  device-specific CP22 implementation or a source of SU(3) particles.
- Recorded `F_e` is unconditional over all measured pointer/Kraus
  outcomes and separately over flagged loss; it does not postselect
  the observed code population to boost fidelity.
- The candidate has **not** passed PR #68's 0.99 full-instrument
  certificate, and a favorable noise-proxy result alone would not
  supersede that failure.
- No conclusion about real QCD, a continuum SU(3) simulator, glueball
  dynamics, physical transmon errors inside a gauge simulator, or a
  reliable error-correcting quantum memory follows from this proxy.

**Next decision:** only if an effective-noise region is favorable
should we invest in full physical-control synthesis. Otherwise prioritize
a native gauge-register/charged-bus encoding and a same-device
digital-versus-analog hardware model before more blind pulse-depth growth.

## First completed numerical results (dedicated CI)

The [initial workflow](https://github.com/HiroSakuraba/Standard-model-generation-as-quantum-error-correction/actions/runs/39352785779)
passed **10/10 new regression tests** and completed both 104-control
pulse schedules under the assumed `T1=30 μs` active-particle loss
and pointer decay, `Tphi=45 μs` pointer/electric dephasing,
300 ns readout, 300 ns reset and 0.005 readout-label error.
The reported values are the **unconditional logical entanglement
fidelity including separately flagged four-particle erasure**, NOT
postselected fidelity. They are reproducible within the specified
phenomenological model.

| Schedule | Noisy compiled recovery | Same-time passive storage | Difference |
|---|---:|---:|---:|
| AS, illustrative coherent pointer level-3 mixing | 0.304206 | 0.376978 | **−0.072772** |
| BS, illustrative coherent pointer level-3 mixing | 0.311722 | 0.372601 | **−0.060878** |
| AS, pointer level-3 mixing disabled | 0.304254 | 0.376978 | **−0.072724** |
| BS, pointer level-3 mixing disabled | 0.311751 | 0.372601 | **−0.060849** |

For these settings, coherent leakage into the extra pointer
level is tiny (approximately **3.36×10⁻⁵ AS** and **1.53×10⁻⁵ BS**
of the final physical pointer population). Its suppression changes
`F_e` only by about **0.000048** and **0.000029** respectively.
Those small numbers should **not** be generalized to genuine
CP22 ququart leakage in another device: no hardware mapping exists.

**Decision:** neither compiled schedule beats matched passive storage
in this explicit noise proxy, even with pointer-level mixing disabled.
The dominant observed issue is that the compiled recovery already
fails to reproduce the ideal three-outcome instrument. This is a
negative result for **these two fixed schedules and assumptions**;
it is not evidence that symmetry-protected recovery is impossible
with a different encoding, control alphabet or hardware.
