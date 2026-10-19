# Local-support gauge-covariant logical syndrome via a three-cycle

19 October 2026. Continues merged PR #47: finite S3 subgroup of SU3, two locally gauged vertices, four colour-qutrit matter registers, and one six-state edge link. **This is a new four-dimensional eigenphase code in that finite-group testbed; it does not implement the nine-quark 21D/65-syndrome SU(3) memory.**

## Why this is different from the Gauss-check circuit

The previous work established that a **gauge-invariant dressed exchange** is accepted by the Gauss projectors with probability 1 whether it occurred or not. A measurement of local Gauss charge cannot diagnose it. This batch constructs a different, logical syndrome from two gauge-invariant operations with bounded spatial support:

- \(A\): SWAP the two qutrits **within vertex A**. It commutes with independent vertex gauge actions.
- \(S\): SWAP matter qutrit A0 with B0, **dressed by their connecting finite-group link**, as constructed and verified in PR #47. It commutes with both vertex Gauss operators, but must interact with the link as well as matter on each side.

Both satisfy \(A^2=S^2=I\). For the modeled representation they obey the braid relation \(ASA=SAS\), so \(U=AS\) satisfies \(U^3=I\). Each application of U consists of S followed by A; its factors have vertex or edge support. In the **14-dimensional physical Gauss sector**, U has eigenphases \(1,\omega,\omega^2\) with multiplicities **6,4,4**, where \(\omega=e^{2\pi i/3}\).

Choose \(K=\ker(U-\omega I)\), a **four-dimensional** gauge-invariant logical code. Since \(SUS=U^{-1}\), the error S sends K into the **orthogonal** four-dimensional \(\omega^2\) eigenspace, without changing its internal encoded amplitudes. This choice is *not* the arbitrary balanced ±S-eigenvector code from PR #47.

## An explicit pointer circuit, not an assumed projector

Use **one qutrit pointer**. Prepare the uniform superposition \((|0\rangle+|1\rangle+|2\rangle)/\sqrt3\). Apply U conditioned on pointer labels 1 and 2, then U once more conditioned only on pointer label 2. Apply the inverse three-point discrete Fourier transform and measure the pointer. This exactly realizes the three spectral Kraus operators

\[
M_m=\frac13\left[I+\omega^{-m}U+\omega^{-2m}U^2\right],\quad m=0,1,2.
\]

In particular a state in K produces **pointer label 1 with certainty**; the same state after S produces **label 2 with certainty**; the pointer does not reveal which of the four logical basis states was encoded.

The conditioned U sequence expands into **four high-level controlled operations**: **two pointer-controlled link-dressed swaps S, two pointer-controlled on-site swaps A**, in addition to one pointer preparation, one three-point inverse Fourier, and at most one conditional corrective S. The program executes these controlled branches directly on all **486 system amplitudes**, preserves the norm of arbitrary mixed-reference inputs, and verifies exact gauge covariance at both vertices. No global logical \(P_K\) or \(P_{SK}\) is applied by fiat.

**Important architectural limitation:** A *pointer-controlled dressed edge swap* acts jointly on the pointer, matter qutrit A0, the gauge link and matter qutrit B0. We have **not** broken this four-register primitive into one- and two-register gates, assigned pulse strengths, counted routing steps, or modeled resulting noise. Calling this "native local implementation" would be incorrect. It is a **constructive circuit using gauge-covariant high-level operations with edge-local support**, rather than a compiled native gate sequence.

## Exact recovery of one known-pair exchange

Consider the CPTP error channel \(\mathcal E(\rho)=(1-p)\rho+pS\rho S\). Measure the three-level pointer. On syndrome 1, leave the system alone. On syndrome 2, apply one dressed S; this maps \(SK\to K\). Syndrome 0 is the complement and is unused for this noise channel. The state is returned to K and decoded, with entanglement fidelity exactly **1** for any \(p\in[0,1]\). The measured states retain arbitrary four-dimensional logical superpositions.

An illustrative symmetric 1% classical pointer-report error, under a **flag-and-replace** policy for misreports, yields \((1-0.01)\cdot1+0.01/16=0.990625\) entanglement fidelity. This is a **specified logical CPTP scenario**, not a calibrated qutrit readout error, error propagation through conditional swaps, or pulse-level decoherence.

## What this result does not show

- The finite S3 representation is **not** the continuous local SU(3) Gauss group, and physical matter spin, statistics and link dynamics are absent.
- The four-dimensional phase code is **not** the original **21-dimensional nine-quark** sign-twist code, whose arbitrary single-colour-qutrit error span requires **65 independent syndromes**.
- A single known-pair dressed swap is correctable by this instrument. It does **not** prove correctability of unknown exchange locations or general gauge-invariant errors. Connected exchange algebra can still destroy passive protection, as established in PR #46.
- The pointer's spatial movement/coupler layout and decomposition of its controlled four-register edge swap into native primitives are open. Thus the **resource count is in high-level gauge-covariant gates**, not device-calibrated two-body gates or runtime.
- A Gauss-syndrome measurement alone still cannot detect the dressed swap. The new measurement extracts **an eigenphase of a gauge-invariant loop of operations**, not a Gauss-charge violation.

## Reproduce and next engineering check

\`\`\`bash
OPENBLAS_NUM_THREADS=1 python -m unittest -v test_colour_s3_edge_logical_syndrome
OPENBLAS_NUM_THREADS=1 python colour_s3_edge_logical_syndrome.py
\`\`\`

The dedicated GitHub workflow runs both commands and archives the report. The next substantive step is to **compile the controlled edge swap into physically specified adjacent-register two-body gates** while preserving its gauge covariance at each intermediate operation, or explicitly quantify the gauge-breaking introduced by individual compilation gates and their correction. Only after that should physical T1/T2, link leakage, fault probability and end-to-end recovery be compared with a matched unencoded baseline.

This limited, constructive success is compatible with the earlier passive-protection no-go: **active logical syndrome measurement adds a resource beyond gauge symmetry**.
