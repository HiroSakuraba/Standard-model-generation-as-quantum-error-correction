# Mobile gauge-charged mediator: local link transport and logical syndrome

19 October 2026. Builds on PRs #47–52 of the **finite non-Abelian S3 subgroup** gauge-link simulator. The model contains four matter qutrits at two independently gauged vertices and a six-element oriented gauge link. This is **not continuous SU(3)** and **not** the original nine-quark 21D/65-syndrome recovery architecture.

## Precise engineering question

PR #50 established that a neutral syndrome pointer combined with strictly gauge-invariant one-/two-register primitives cannot measure the desired dressed-swap syndrome: the link representation label is conserved. PR #51 constructed an ideal gauge-invariant 3-register dressed-swap Hamiltonian, but its neutral-pointer-controlled version required a **four-register elementary Hamiltonian**. PR #52 tested alternative mediator/penalty assumptions without actually implementing link transport.

**This experiment changes the register inventory:** introduce a *gauge-charged* mobile carrier M, with location A/B and a three-dimensional colour state (six carrier basis states), together with a stationary conjugate-charge qutrit R at vertex A. The initial local singlet \((|00\rangle+|11\rangle+|22\rangle)/\sqrt3\) between M and R is gauge invariant because the chosen signed-permutation S3 subgroup matrices are real orthogonal. Preparing an isolated pure charged fundamental without a balancing reference would not be locally Gauss invariant.

The neutral syndrome pointer remains a qutrit. Both independent vertex Gauss actions are explicitly extended to act on the charged mediator *at its current location*, plus its stationary conjugate reference.

## An exact TWO-register charge-transport operator

For the oriented group link \(|g\rangle\), the transport T acts on **link+mobile mediator only**:

\[
T\,|g,A,v\rangle=|g,B,D(g^{-1})v\rangle,\qquad
T\,|g,B,w\rangle=|g,A,D(g)w\rangle.
\]

The matter qudits and reference are spectators. T is an involutory unitary. If \(g\mapsto h_A g h_B^{-1}\) and the carrier transforms as \(D(h_A)\) at A or \(D(h_B)\) at B, the link rotation ensures T commutes exactly with **both** vertex Gauss transformations. The tests check both covariance and T²=I on arbitrary states.

A local swap between the mediator colour and an individual matter qudit at the carrier's current vertex acts on **mediator+one matter qudit** (two registers) and commutes with the vertex gauge action; the mediator's location selects its active site. Controlling the local swap on the neutral qutrit pointer increases its support to **three registers**. No gate with four-register simultaneous support is used in this high-level decomposition.

## Reversible exchange bus and pointer readout

Starting with the mediator at A, the exact sequence

1. locally swap M with matter A0;
2. transport M from A to B;
3. locally swap M with matter B0;
4. transport M from B to A;
5. locally swap M with matter A0 again

implements the **link-dressed SWAP** \(S_{\rm dressed}\) between A0 and B0 while returning M to vertex A and its original colour, still entangled only with the stationary reference R. This is an exact identity for arbitrary initial matter amplitudes. The link itself is unchanged. The nonlocal three-register dressed exchange has become **five individually Gauss-covariant two-register operations**, using a new mobile charged degree of freedom.

The coherent qutrit logical-syndrome circuit of PR #48 measures the eigenphase of \(U=A_{\rm onsite}S_{\rm dressed}\). For each pointer-conditioned edge exchange, condition the three *local matter–mediator swaps* on pointer outcomes; perform the two link transports unconditionally. In inactive pointer branches, two successive transports undo each other, so they cause no logical action. To condition the on-site A0–A1 SWAP, use three pointer-conditioned mediator-to-matter exchanges entirely at A.

For two edge-S windows and two onsite-A windows, the circuit comprises:

| Primitive | Number | Register support |
|---|---:|---:|
| Unconditional link–charged-mediator transports | 4 | 2 |
| Neutral-pointer-controlled local charged-mediator/matter swaps | 12 | 3 |
| **Total coherent syndrome operations** | **16** | **at most 3** |
| Classical-outcome-conditioned dressed-swap feedback | 5 | 2 |
| Pointer preparation/Fourier | 2 | 1 |

These are exact *register-locality counts*, **not** calibrated native pulse counts. In particular, the three-register neutral-pointer-controlled local swap remains **an assumed elementary controlled operation**. It has not been decomposed into a universal *gauge-invariant two-register* alphabet, so this is not a contradiction of PR #50's restriction. A mobile carrier with charged reference and location coherence also requires its own preparation, transport, coherence times and readout/reset. The result is a constructive **max-register-support reduction from 4 to 3**, not a demonstrated native-gate breakthrough.

## Checks on the complete gauge-constrained system

The numerical implementation uses **all 486** original gauge-system amplitudes plus a neutral pointer, two-position/three-colour charged mediator, conjugate reference and all logical reference columns. It compares the full pointer instrument with PR #48 on all four code basis columns; it does *not* synthesize or apply the global code projectors.

The instrument is traced over the pointer and the **entire mobile/reference environment** when computing logical entanglement fidelity. The tests explicitly verify:

- singlet preparation satisfies both independent Gauss constraints;
- each link transport and local swap commutes with both extended vertex actions even on arbitrary off-Gauss inputs;
- the mediator and reference return to a disentangled singlet for the ideal circuit;
- every one of 16 intermediate gates preserves Gauss law;
- clean and known dressed-swapped states produce the correct distinct pointer labels, while leaving arbitrary code amplitudes undisturbed;
- the ideal source channel \(\{\sqrt{1-p}I,\sqrt p S_{\rm dressed}\}\) is recovered with **unconditional entanglement fidelity 1** (with all ancilla outcomes included), for p=0,0.37,1.

An additional **deterministic coherent pulse-angle error** sensitivity scan applies the same fractional overrotation to every mediator swap/transport and the classically commanded feedback. Overrotations of gauge-commuting involutions still preserve Gauss symmetry, but can degrade recovered logical fidelity and fail to return the mobile reference cleanly. This is not T1/T2, leakage or measured hardware calibration.

## Limitations and release criteria

No waveform/native physical realization of the location-coherent mobile charged carrier, selective three-register pointer-controlled exchange, SU(3)-compatible link transport, or reference-singlet preparation is established. Link/mediator movement, routing, speed, noise, cross-talk and reset are unmeasured. A true dynamical gauge-field transport may have different allowed charges and interactions. This does **not** implement the nine-quark code's 65-syndrome readout, nor prove practical quantum-memory advantage over ordinary codes.

**Most important next hardware test:** synthesize the three-register pointer–charged-mediator–matter controlled SWAP from a *specified* local gauge-respecting native Hamiltonian and compare its noise and duration with the four-register ideal pulse and the gauge-breaking U(9)/U(18) pair decomposition. Do not simply assume all 16 high-level operations are cheap.

Reproduce:

    python -m pip install numpy scipy
    OPENBLAS_NUM_THREADS=1 python -m unittest -v test_colour_s3_charged_bus_transport
    OPENBLAS_NUM_THREADS=1 python colour_s3_charged_bus_transport.py

The dedicated workflow archives the complete Gauss trace, exact residuals and unconditional coherent-overrotation sensitivity results.
