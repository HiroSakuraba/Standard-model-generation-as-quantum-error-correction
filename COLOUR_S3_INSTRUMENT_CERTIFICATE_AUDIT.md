# Finite-S3 three-outcome instrument certification (research audit)

19 October 2026. Follows merged PR #67. **This is a stricter acceptance test, not a compiled physical quantum decoder.**

## Problem

The four-logical-state active code corrects seven *individual* admitted faults
under an ideal syndrome-and-recovery operation. PR #67 optimized the numerical
unconditional entanglement fidelity of the restricted seven-fault channel.
It reported approximately 0.808 (AS) and 0.823 (BS) at 104 controls before
physical loss. Critically, the calculated recovered logical fidelity
**still credits an ideal final code-subspace projection**. It does not establish
that the circuit extracts all syndromes and performs conditional feedback.

In addition, any circuit should handle a six-dimensional *flagged complement*
without quietly discarding it. A high logical score alone is insufficient.

## New test: all 14 initial Gauss-physical states

The function `certify_three_outcome(U, label)` accepts an **explicit 42×42
unitary**, acting on 14 Gauss-invariant physical system states times a neutral
three-level pointer initialized to 0. It audits three mutually orthogonal
input subspaces:

| Input | Rank | Required output after syndrome and feedback |
|---|---:|---|
| Clean code `C` | 4 | `C` with pointer 0 |
| Electric syndrome `V` | 4 | `C` with pointer 1 |
| Orthogonal flag `T` | 6 | Same flagged space with pointer 2 |

For **every** branch it reports the full pointer-probability vector,
probability of jointly occupying the proper pointer and physical output
sector, overlap with the *actual target output mapping*, and raw Frobenius
distance. The certificate requires both intended pointer probability
and branch coherent transfer fidelity to exceed 0.99 by default. That
threshold is an explicit research acceptance setting, not a hardware
error-correction threshold.

The ideal global unitary from PR #64 passes the mathematical certificate
as an upper bound; the literal no-gate identity fails the electric and flag
requirements. Neither finding claims new hardware implementation.

All seven original single-fault Kraus alternatives are independently
evaluated. The ledger **never postselects** pointer 2 or falsely allocates
unrecognized output outside `C` to recovered information.

## Essential correlated-fault counterexamples

The verifier also independently applies all 12 `gR` and `Rg` ordered
two-fault strings, for each order-six group (`AS`, `BS`).
Here `gR` means **electric first, then exchange** when operators act on
the ket from right to left. With the mathematically ideal one-round
instrument, `Rg` succeeds for all six; `gR` fails on four of the six
(`F_e=11/96`, pointer-2 flag probability `13/24`). This verifies the
known PR #64 limitation and prevents an incorrect generic two-fault claim.
Testing the ideal global instrument is an **algebraic regression**, not a
demonstration of fault-tolerant repeated extraction.

## Hardware/resource boundary

The 14D Gauss-preserving model does not contain native pulse durations,
complete dynamical leakage, charged-mediator transport noise, charge/reference
singlet preparation, pointer assignment noise, measurement latency,
conditional feedback time or reset. Those must be accounted for in a
subsequent local decoder implementation and against duration-matched passive
and no-gate controls. **The diagnostics use code projectors as observables
only**, not as cost-free physical gates; no physical code-projection or
correction is inserted by this new module.

The pass flag therefore has the deliberately narrow name
`mathematical_three_outcome_instrument_pass`. It is **not** a hardware
readiness verdict. A numerical 42D product of permitted Hamiltonians can
be checked with the same API, and may fail even if its single-fault
entanglement fidelity rises.

## Reproduce

```bash
python -m pip install numpy scipy sympy
OPENBLAS_NUM_THREADS=1 python -m unittest -v test_colour_s3_instrument_certificate
OPENBLAS_NUM_THREADS=1 python colour_s3_instrument_certificate.py
```

The executable prints two controls: **no gates** and **uncompiled
global ideal**; no optimizer is run and no misleading device parameters
are silently assumed.

## Next step

Use the certificate in a *branch-constrained optimizer*, retaining all
14 input directions alongside the seven-Kraus logical score, and demand
an independently recorded charged-mediator/pointer/feedback/measurement
resource ledger. Only a finite sequence that passes the complete
instrument checks, survives device-specific noise and beats matched controls
would support operational quantum-memory advantage.

**Scope:** finite signed-permutation `S3` link with four colour qutrits;
not continuous SU(3) QCD, glueball spectroscopy, or a Standard Model memory.
