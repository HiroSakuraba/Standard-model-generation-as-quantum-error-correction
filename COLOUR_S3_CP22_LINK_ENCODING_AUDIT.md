# Finite-S3 physical link encoding: exact selective CP22 mappings and the Gauss-path obstruction

October 19–20, 2026. Follows PR #71's matched-time physical-noise result. **This is an exact gate-algebra witness, not a calibrated transmon decoder.**

## The first native-gate bridge that can actually be checked

The S3 gauge link has six states, and its regular representation is naturally structured as a parity bit and the three elements of a cyclic subgroup. Choose the group-element ordering

`[e,r,r²,s,sr,sr²]`

with `r³=e` and `s²=e`. This is a computational relabeling of the six link states, allowing a **parity qubit × cycle qutrit** hardware register.

Fourier-transform *within each three-element cycle*. In this precisely specified encoding, the **central standard-representation electric projector** of the existing 14D Gauss model satisfies the *exact matrix identity*

`U P_standard U† = I_parity ⊗ diag(0,1,1)`.

The standard representation has rank four in the 6D link Hilbert space, the other two states being the trivial and sign irreducible sectors. The script compares this decomposition with the independently implemented Peter–Weyl character projector from PR #50 and asserts machine-precision agreement.

It follows that the **neutral-pointer × link standard-electric entangler**,

`exp(-i θ |1><1|_pointer ⊗ P_standard)`,

has an exact finite decomposition in the encoded link basis:

`CP(pointer=1,cycle=1,θ) × CP(pointer=1,cycle=2,θ)`,

independent of the parity bit. Each conditioned phase is a native-*style* selective `|22>` interaction, conjugated by local qutrit level swaps to address the selected levels. This gives **two ideal arbitrary-angle CP22 primitives** and local basis permutations. It is a **specific constructive mapping of one of the 13 abstract Hamiltonian controls**, not the full recovery or a measured implementation.

The unconditioned link-electric rotation uses **two local cycle-qutrit phases and zero entangling CP22s**, after this encoding. The same-colour equality projector on a single matter vertex uses **three selective CP22s** on `|00>, |11>, |22>`. An onsite matter–matter qutrit SWAP has an explicit but deliberately **nonoptimal upper bound of 12 selective CP22s**: three modular SUM transformations, each built from four selective phases and local qutrit Fourier transforms, followed by a local modular negation. The gate identities are checked on every basis vector and against the finite signed-S3 local action.

| Model-level Hamiltonian/control | Exact idealized CP22 count | Other required resources |
|---|---:|---|
| Pointer-conditioned central link electric `|1><1| ⊗ P_standard` | **2** | Parity×cycle encoding, per-gate qutrit level swaps |
| Unconditioned standard-electric link rotation | **0** | Two local cycle-qutrit phase shifts |
| Vertex-local equal-colour projector `Q` | **3** | Pair-local qutrit level swaps |
| Onsite matter qutrit SWAP | **12** (constructive upper bound) | Six qutrit Fourier transforms, one modular negation, local permutations |
| Charged mediator transport and preparation | **uncompiled** | Six-state link, charged matter, reference singlet, gauge-covariant movement |
| Complete pointer measurement / feedforward / reset | **uncompiled** | Control timing, assignment, resets, hardware-native encoder |

All CP22s here are **ideal arbitrary-angle two-qutrit conditional phases**. The separate four-level Duffing pulse from earlier five-plaquette SU(3) work does **not** establish that these arbitrary phases are accurately achievable on this device geometry. Mapping the link to a parity qubit and cycle qutrit is not a free physical preparation/routing operation if the real device starts in the six group-element basis.

## Key negative discovery: intermediate gates break the Gauss constraint

A complete projector `P_standard` commutes with **both** independent vertex gauge transformations. However, addressing only the Fourier-cycle mode 1 or 2 does **not** commute with all right-vertex S3 transformations. The script explicitly constructs both six-state link regular actions, checks all local commutators and evolves all 14 gauge-invariant code-basis states in the full 486D system.

Consequently:

- A one-mode selective CP22 phase can cause intermediate Gauss-sector leakage.
- Completing the second equal-phase mode restores the exact central gauge-invariant rotation.
- The complete **algebraic gate identity is not sufficient to call the sequence strictly Gauss-preserving at every pulse**. A device with a strong gauge-energy penalty would require a separately analyzed perturbative treatment, or a native simultaneous two-mode interaction commuting with the gauge generators at all intermediate times.

This distinction is particularly important because previous Lie-algebra results established reachability using *gauge-preserving abstract controls*. We have **not** shown that each of those controls has a strictly gauge-preserving realization in the selective-CP22 native alphabet.

A similar issue occurs for finite-S3 same-colour equality: the complete three-phase product is invariant, but a single colour-specific CP22 factor is not. It also occurs for the 12-CP22 onsite swap construction. **Do not use the exact final gate decomposition as a substitute for a continuous-time gauge-covariant pulse certificate.**

## How to reproduce

```bash
python -m pip install numpy scipy sympy
OPENBLAS_NUM_THREADS=1 python -m unittest -v test_colour_s3_cp22_link_mapping
OPENBLAS_NUM_THREADS=1 python colour_s3_cp22_link_mapping.py
```

The output is a machine-readable JSON with exact residuals at θ=π/7, π/2, π, two-vertex Gauss checks, CP22 operation-count upper bounds and a complete list of unsupported hardware steps. The dedicated CI workflow saves the JSON artifact.

## Next explicit engineering milestones

1. **Gauge-preserving elementary gates:** find or rule out a strictly invariant native route to the two-mode electric interaction; model the penalty-dependent intermediate Gauss excitations if strict invariance is relaxed.
2. **Charged mediator:** encode link×mobile charge into qubit/qutrit/ququart registers and decompose its two-register transport, with charge/reference singlet preparation, into addressable CP22 and single-device operations. Quantify routing.
3. **Syndrome pointer:** explicitly compile controlled electric interactions, pointer preparation, readout, feedback and reset into a repeatable *end-to-end* three-outcome instrument. No global free code projection.
4. **Noise and comparison:** use the same hardware native pulse times, amplitude damping/dephasing, leakage and gauge penalty for passive, digital and analog echo approaches.

**Scientific status:** an exact native-style decomposition for specific finite-S3 Hamiltonians and an equally exact intermediate-symmetry warning. **Not** a complete SU(3) quantum error-correcting memory, physical device calibration, or a demonstration of advantage over passive storage.
