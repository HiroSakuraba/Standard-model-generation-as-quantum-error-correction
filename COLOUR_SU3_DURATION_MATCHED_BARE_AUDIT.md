# Does the microscopic finite-S3 code beat an unencoded four-dimensional memory?

**19 October 2026.** Follow-up to merged PR #57. This is a **duration-matched direct-storage baseline and parameter-sensitivity analysis**, not a demonstration that the two physical implementations are equivalent or that a quantum-error-correcting advantage has been established.

## 1. What exactly is compared?

The **encoded candidate** is the established four-dimensional phase code \(K\) in the finite \(S_3\) gauge-link model, with the 495-dimensional four-fermion persistent Hubbard simulator from PRs #55–57. It applies fifteen microscopic local Hubbard pulses (including actual feedback), six ideal charge-link transports and a nominal 300 ns syndrome readout interval. Its environmental no-jump entanglement fidelity \(F_{\rm code}\) is recomputed directly by the PR #57 model for \(U/t=20,30,40\). Its four active fermions are subject to uniform loss. Any jump is counted as an orthogonal flagged erasure.

The **unencoded control** stores the *same four-dimensional logical Hilbert space* in **two fundamental colour qutrits**, each using colours \(|0\rangle\) and \(|1\rangle\) as a two-level memory. Exactly **two fermions** are occupied. The direct memory is held idle for **precisely the same total elapsed duration** as the encoded candidate. It receives the same per-fermion loss lifetime and, for environmental comparisons, the same *per-number-operator* dephasing rate. There is no assumed encoding, correcting or syndrome readout operation in direct storage.

We report **two different unencoded logical channels**, which must not be confused:

1. **Clean direct storage:** no imposed logical error. This measures the cost of attempting to correct a fault that is not present. It is a deliberately strong unencoded reference, not an error-correction challenge.
2. **Matched 37% *logical-fault probability*:** with probability 0.37 apply \(X\) to the first of the two stored colour qubits and leave it uncorrected. The \(X\) branch is traceless in the logical four-dimensional space. The encoded candidate instead receives the known-pair **dressed-exchange** channel with probability 0.37, as in PRs #48–57. **These are not the same microscopic physical noise process.** The matching is **logical dimension and source-event probability only**, so any crossover is a controlled *phenomenological comparison*, not evidence of a physically superior gauge device.

No third-party external device benchmarks or measured T1 values are presented; all physical-scale parameters are hypothetical.

## 2. Exact two-fermion loss and dephasing channel

Let \(T\) be equal elapsed storage time, \(T_1\) the lifetime of each active fermion, and \(F_{\rm code}\) the encoded entanglement fidelity computed from the full microscopic circuit **conditioned on no particle loss**. Independent uniform annihilation jumps give exact no-loss weights:

\[
\eta_{\rm code}=e^{-4T/T_1},\qquad
\eta_{\rm bare}=e^{-2T/T_1}.
\]

These are properly **unconditional** performance factors when missing particles are flagged and unrecoverable. The encoded model omits particle loss on the extra spectator B1/reference registers, so its quoted loss performance remains optimistic.

For **colour-blind orbital-charge dephasing**, the direct memory has fixed total occupation \(n_{\rm orbital}=1\), so the noise operator is a scalar on its code space, giving \(F^{(0)}_{\rm bare}=1\) when no logical fault occurs. The encoded model must still *propagate its actual Hubbard dynamics* for the same colour-blind bath; its \(F_{\rm code}\) is not replaced by 1.

For **colour-resolved occupation dephasing**, \(L_{j,a}=\sqrt{2\gamma}n_{j,a}\). Two occupied colour modes per direct-memory qubit are alternative logical basis states. Each qubit's off-diagonal coherence decreases by \(c=e^{-2\gamma t}\), where \(t=T/[1000/(2\pi(t/h)_{\rm MHz})]\) is dimensionless hopping time. The exact bare two-qubit channel has no-fault entanglement fidelity

\[
F^{(0)}_{\rm bare}=\left(\frac{1+c}{2}\right)^2.
\]

With the known-probability logical-\(X\) challenge but no error correction, every \(X\) error Kraus term remains traceless under diagonal dephasing, so the direct memory has

\[
F^{(X)}_{\rm bare}=(1-p)F^{(0)}_{\rm bare},\qquad p=0.37.
\]

The program independently verifies this analytic expression using the explicit \(4\times4\) Kraus matrices of two independent qubit dephasing channels and a first-qubit \(X\) fault. There is **no stochastic sampling uncertainty on the bare baseline**.

## 3. Analytic threshold: precisely when can extra correction overhead pay?

For the baseline colour-blind/no-dephasing anchor, the **unconditional** quantities are

\[
F^{\rm uncond}_{\rm code}=F_{\rm code}e^{-4T/T_1},
\quad
F^{\rm uncond}_{\rm bare,X}= (1-p)e^{-2T/T_1}.
\]

Thus the encoded candidate exceeds the uncorrected challenged bare memory **iff**

\[
F_{\rm code}>1-p
\quad\text{and}\quad
T_1>\frac{2T}{\ln[F_{\rm code}/(1-p)]}.
\]

If \(F_{\rm code}\le 1-p\), then there is **no finite-lifetime crossover**. The program also evaluates the general analytic threshold with colour-resolved \(F_{\rm bare}^{(0)}\), replacing \(1-p\) with \((1-p)F_{\rm bare}^{(0)}\).

By contrast, the **perfect idle direct memory** has no logical fault and has \(F_{\rm bare}^{(0)}=1\). With no dephasing and a two-versus-four active-fermion loss penalty, an encoded circuit with \(F_{\rm code}\le 1\) cannot beat it at any finite \(T_1\). This is an exact restriction of the chosen comparison, **not a no-go theorem for all quantum error-correcting codes**.

### Numerical anchors (from PR #57, recomputed on main)

At hypothetical \(t/h=10\) MHz, each link transport 100 ns, readout 300 ns and active-particle \(T_1=30\,\mu{\rm s}\):

| \(U/t\) | Total time (μs) | Code \(F_e\) before loss | Code after loss | Bare clean after loss | Bare with 37% uncorrected \(X\) | Break-even \(T_1\) (μs) |
|---:|---:|---:|---:|---:|---:|---:|
| 20 | 8.419 | 0.7727 | 0.2515 | 0.5705 | 0.3594 | 82.46 |
| 30 | 12.160 | 0.9119 | 0.1802 | 0.4446 | 0.2801 | 65.76 |
| 40 | 15.907 | 0.9714 | 0.1165 | 0.3463 | 0.2182 | 73.48 |

**None of the tested encoded circuits wins even the matched-37% challenged direct-storage comparison at \(T_1=30\,\mu s\)**, and none beats clean direct storage at any finite lifetime under this ideal orbital-charge bath.

This reverses the misleading impression that improving one local Hubbard controlled-SWAP or the loss-free recovered entanglement fidelity necessarily improves system performance. Device lifetime and total circuit time dominate the stated scenario.

## 4. Sensitivity surfaces and what is actually held constant

The code calculates a **180-point closed-form sensitivity grid**:

- \(U/t \in \{20,30,40\}\)
- hypothetical hopping \(t/h\in\{5,10,20,40\}\) MHz;
- readout delay \(\in\{100,300,1000\}\) ns;
- active-fermion \(T_1\in\{30,50,75,100,300\}\) μs.

Each cell uses the **same encoded end-to-end elapsed time** for the bare memory and the coded run. The optimized *dimensionless* Hubbard pulse and its **coherent** \(F_{\rm code}\) are held constant while converting frequency to physical time. Changing hopping without changing U, decoherence rates, actuator bandwidth, or coupling strengths on real hardware is **not assumed physically available**; this is a sensitivity analysis.

Separately, at \(U/t=30\), the script reruns the encoded **full microscopic feedback** with the two PR #57 stochastic Lindblad-number baths (\(\gamma=0.0002\) per operator) and evaluates the bare storage channel analytically at that same total elapsed time:

- four active **orbital total-number** operators for coded dynamics versus two scalar orbital operators for direct storage;
- twelve active **colour-resolved mode-number** operators for the code versus four nontrivial logical-basis colour occupations in direct storage.

**Equal \(\gamma\) per operator is not matched total environmental noise power.** Therefore these results only demonstrate sensitivity to a specified bath structure; they do not constitute a device-equivalent fairness result. Encoded estimates use **six stochastic trajectories** and report Monte Carlo standard error; the analytic bare values do not.

## 5. Conclusions and release criteria

1. The current finite-S3, four-dimensional Hubbard-syndrome construction is **not a demonstrably advantageous quantum memory** at the nominal hypothetical particle lifetimes and frequencies.
2. The comparison is explicit and falsifiable: only when lifetime, hopping and readout parameters cross a **known analytic threshold** can this circuit beat an uncorrected, probability-matched bare logical-\(X\) challenge under the simplifying assumptions.
3. The complete code still assumes **ideal finite-S3 charge-link transport**, stationary reference preparation and survival, perfect pointer Fourier/readout and no loss on spectators. Physicalizing these likely makes the encoded cost worse.
4. The source dressed swap is not a bare logical-\(X\) physical channel. **No claim of a physically equivalent error environment or genuine encoding advantage** is made.
5. Neither the original **21-dimensional nine-quark** Standard-Model-generation code nor its **65 logical-syndrome** instrument is implemented.

### Reproduce

\`\`\`bash
python -m pip install numpy scipy
OPENBLAS_NUM_THREADS=1 python -m unittest -v test_colour_su3_duration_matched_bare
OPENBLAS_NUM_THREADS=1 python colour_su3_duration_matched_bare.py
\`\`\`

The dedicated GitHub workflow archives exact bare Kraus and lifetime thresholds, 180 scenario cells, seed-controlled encoded dephasing estimates, and every scientific limitation as JSON.

**Next experiment:** to establish a physically meaningful code-versus-unencoded claim, simulate the **same microscopic fault/environment** on both real architectures, adding native gauge-link transport, encoded auxiliary/reference losses and the actual encoder/decoder timing. Alternatively study **which additional physically allowed errors** the code corrects that a bare colour-ququart cannot without an incomparable phenomenological fault substitution.
