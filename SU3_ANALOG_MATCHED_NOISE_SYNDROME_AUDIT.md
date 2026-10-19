# SU(3) analog refocusing: equal-time noise and two-qutrit syndrome-pointer model

19 October 2026. Continues merged PRs #39–#41.

## Research question

Does finite-width cubic-Casimir analog refocusing beat *uncontrolled native evolution on the same hypothetical five-qutrit register*, at **the same elapsed time and the same sampled local noise**? Does an explicit two-pointer, six-outcome measurement instrument retain the predicted recovery advantage?

Our reference is the source-truncated 243-state SU(3) Hamiltonian and prior finite-width pulse routine in su3_native_control_compilation.py. We propagate **both arms through the same 72 segment durations**, including 32 finite pulse windows. During windows, the controlled arm evolves under H plus the signed cubic-clock generator, while the uncontrolled arm evolves under H alone. This fixes the earlier unequal-time comparator, which evolved the baseline for only T=0.6 rather than T+32*tau.

After each of the 72 segments, **both** arms receive identical parameterized channel families: three-level bosonic pure-loss with eta=exp(-time/T1), and qutrit Gaussian phase diffusion with variance 2*time/Tphi. We sample one Kraus outcome per site on a 243x2 encoded purification; averaging is a completely positive trace-preserving (CPTP) channel. These are *gauge-fixed register-level* noise assumptions, **not** calibrated physical gauge-field noise. Dissipation is segment-lumped, not driven gate-time Lindblad dynamics.

## Explicit six-sector syndrome model

The earlier exact Knill–Laflamme analysis provides five mutually orthogonal two-column syndrome isometries V_r. Define P_r=V_r V_r^dagger, each rank 2. Their sum has rank 10; P_5=I-sum P_r is the rank-233 complement.

Two pointer qutrits provide nine orthogonal basis labels. Tag the six sectors with labels (0,0), (0,1), (0,2), (1,0), (1,1), (1,2). The **mathematically exact** controlled isometry maps input |psi>|00> to sum_r P_r|psi>|pointer_r>. Orthogonal projectors and a resolution of identity make this norm preserving. The script constructs and tests this explicitly.

**The sector-controlled isometry is nonlocal and NOT compiled to local native hardware gates.** Two qutrits suffice as a dimension/information bound, but that is not a proof of low-cost extraction.

Pointer errors are **not** just an overall fidelity multiplier. Each qutrit undergoes an exact three-level pure-loss channel over the assumed readout time and then symmetric classical assignment error, yielding a normalized 9x6 confusion matrix. Correctly reported sectors are corrected by the original V_r^dagger; wrong sectors and complement are replaced by I_2/2. This defines a CPTP logical instrument (no postselection). Both arms also incur the same data-register loss/dephasing during assumed extraction, readout, and reset time. Ancilla preparation, nonlocal extraction dynamics, and additional measurement backaction remain unmodeled.

## Hardware and time conversion assumptions

Converting the dimensionless source Hamiltonian to clock time assumes H_phys/h=f_MHz * H_dimensionless. The default f=0.01 MHz, tau=0.002, register T1=30us, Tphi=45us, pointer T1=30us and 800ns readout are **scenarios, not device measurements**. The native Wilson interactions and individually addressable cubic-Casimir clock strengths have not been physically demonstrated here.

The analog and no-control arms are **same-duration/same-noise within this hypothetical simulator**. Comparing either to the digital transmon circuit of PR #41 is **not a controlled hardware comparison**.

## Verification

- Exact noiseless controlled sequence matches the existing finite-width lab_cubic_pulses calculation for both logical inputs.
- Free reference matches exp(-iH*(T+32*tau)), not the shorter prior T-only baseline.
- The 5 rank-2 syndrome projectors plus rank-233 complement sum to identity; pointer-tag isometry preserves arbitrary state norms.
- Pure-loss Kraus operators satisfy completeness; sampled channel trajectories preserve reference-state norm.
- Pointer assignment and relaxation confusion columns sum to one. Perfect channel agrees with previously certified decoder.
- GitHub Actions executes all tests, then produces trajectory JSON with actual sample values, sampling standard errors, timing and assumptions.

Reproduce from the repo root by installing numpy scipy sympy, then run:

    OPENBLAS_NUM_THREADS=1 python -m unittest -v test_su3_analog_matched_noise
    OPENBLAS_NUM_THREADS=1 python su3_analog_matched_noise.py --shots 48

**Do not claim a hardware quantum memory on these results.** If noisy analog refocusing beats the matched baseline, it is evidence for this synthetic control/noise scenario, not an engineered gate implementation. A negative result is likewise scenario-dependent.

Next: derive local controlled-projector extraction hardware and a same-device pulse/control model including noise during control, followed by broader parameter scans and tighter statistical error bounds.
