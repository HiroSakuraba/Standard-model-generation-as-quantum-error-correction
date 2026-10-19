# SU(3) research: transmon pulse, physical pair channels, and analog alternative

19 October 2026. Follows successfully merged five-plaquette circuit PR #39.

## Scope and evidence levels

This is an **experiment-design/physics feasibility probe**, not a verified
superconducting implementation or a claim of quantum error correction
from QCD. The source gauge simulation has 243 three-level gauge-invariant
states. The physical ququart pair model has 16 states and **does not**
replace the full five-plaquette quantum channel.

| Investigation | Actually simulated | Not supplied |
|---|---|---|
| Digital five-plaquette | Explicit ideal 10,080/8,160-CP22 sequences and ideal recovery, reused from PR #39 | Full 4^5 physical channels; noisy recovery; calibrated CP22 pulses |
| Physical two-transmon pulse | 4 levels per transmon, Duffing anharmonicities, exchange, smooth frequency-detuning excursion, 16x16 exact unitary | Fabricated chip parameters; optimized selective arbitrary-angle CP22; spectator and drive imperfections |
| Damping/dephasing | Explicit trace-preserving four-level local Kraus channels, tensor-product two-device channel, including final leakage and unconditional process entanglement fidelity | Lindblad dissipation **during** pulse; cross-correlated noise; gate-by-gate 4^5 trajectories |
| Analog SU(3) refocusing | Full source five-plaquette H remains on during 32 finite-width cubic electric-clock pulses; same ideal decoder | Proof that four addressed cubic-Casimir clocks, or tunable Wilson quadratures, exist on the chosen device |
| Syndrome extraction | Six distinguishable outcomes imply 2 pointer qutrits minimum (9 states) or 3 bits; parameterized readout/reset scenarios | A local physical implementation of the 5 two-dimensional syndrome subspaces plus complement, noisy extraction fidelity |

### Device-parameterized pulse model

In a rotating frame of transmon A, we simulate a model

H/h = (alpha_A/2) n_A(n_A-1) + delta(t) n_B +
      (alpha_B/2) n_B(n_B-1) + g(a_A† a_B + a_A a_B†),

with four levels retained in each subsystem; delta(t) smoothly traverses
an avoided crossing with a sin² detuning excursion. Time evolution uses
exp[-i 2π H_MHz (time_ns/1000)] over many steps. A corresponding idle
unitary is **mathematically divided out as a diagnostic**, not as a
free, executable hardware operation. Only **local diagonal** qutrit
frames, estimated from the |i0> and |0j> phases, are removed.

The script measures, for the 9 computational input states:

1. Residual cross-level spectator phases.
2. Correlated conditional |22> phase.
3. Unconditional process entanglement fidelity to the **target CP22(pi)**,
   including unwanted coherent mixing and fourth-level leakage.
4. Same fidelity when the target is the *actually produced* conditional
   angle, to separate angle-calibration error from other pulse errors.
5. Final computational-subspace leakage after the **post-pulse** local
   amplitude damping/dephasing Kraus channels.

The parameter values are **illustrative scenarios**:
idle detuning 650 MHz, driven detuning candidate 260/290/320 MHz,
anharmonicities -230/-260 MHz, exchange 9 MHz, candidate pulse
durations 40/80/160 ns, with T1=30 microseconds and Tphi=45 microseconds.
None is taken from a characterized device owned by this project.

A high score at one point would not establish a universal
**arbitrary-angle**, high-fidelity selective CP22, which the 40-CP22
compiler assumes. A low score falsifies only these toy pulses, not
all plausible controls.

### Exact local CPTP noise construction

Pure loss uses four bosonic Kraus operators with eta=exp(-dt/T1):

K_l |n> = sqrt(C(n,l)) (1-eta)^(l/2) eta^((n-l)/2) |n-l>.

This has sum_l K_l† K_l = I_4. Phase diffusion uses the positive Schur
kernel G_nm=exp[-dt(n-m)^2/Tphi], factored into diagonal Kraus matrices
through eigen-decomposition. The composed local and product two-site maps
remain CPTP. These channels act **after** the pulse (operator-splitting
limitation); they do not jointly evolve during the driven interaction.

The unconditional process entanglement fidelity is computed as the
Kraus trace sum over the 9 input computational states; it is *not*
a postselected fidelity that hides leaked population.

### Analog/digital common logical target, different hardware assumptions

The five-plaquette analog test uses the already audited
su3_native_control_compilation.lab_cubic_pulses() method. It leaves the
published truncated native Hamiltonian on during 32 clock pulse windows
and decodes with the same idealized syndrome recovery used in the digital
refocusing work. Four addressed sites correspond to **128 single-qutrit
clock actions** (assuming they can be applied in parallel within each of
the 32 windows), instead of 10,080 exact or 8,160 approximate two-qutrit
conditional phase entanglers in the digital construction.

This is **not** an equivalent physically realized machine: the analog
arm assumes direct native SU(3) Wilson interactions and switchable,
gauge-compatible cubic-Casimir clocks. Such controls have not been
demonstrated here. For comparison of clock duration with gate time, we
parameterize H_physical/h = f_MHz H_dimensionless. Example choice
f=0.01 MHz and dimensionless pulse width tau=0.002 is *arbitrary* and
implies clock strength about 10.47 MHz and pulse width about 31.83 ns.
Changing f changes the physical comparison; no continuum/QCD physical
energy conversion is being claimed.

The analog logical replacement noise and prior digital flagged
replacement noise are a **separate illustrative CPTP scenario**.
They are not derived from the two-transmon channel and must not be
mistaken for a like-for-like same-device experiment.

### Reproducibility and decision criteria

Run:

\`\`\`bash
python -m pip install numpy scipy sympy
OPENBLAS_NUM_THREADS=1 python -m unittest -v test_su3_transmon_noise_probe
OPENBLAS_NUM_THREADS=1 python su3_transmon_noise_probe.py
OPENBLAS_NUM_THREADS=1 python su3_analog_digital_noise_comparison.py
\`\`\`

The GitHub Actions workflow runs all three and archives both JSON outputs.
**Do not quote derived numerical results until this workflow actually
finishes and the JSON is inspected.**

The next go/no-go criteria should be:

- Can a physically implementable Hamiltonian produce selective CP22
  at arbitrary angles with very small *unconditional* error, over a
  usable angle range, after realistic waveform constraints?
- Under pulse-time Lindblad dynamics and **repeated full digital sequence**
  noise, can any real gate/time/decoherence regime beat a same-hardware
  uncontrolled baseline after noisy syndrome recovery?
- Can analog cubic-clock and gauge-interaction controls be implemented
  natively with bounded control strength, finite pulses and gauge leakage?
- Is the extraction of six outcomes with negligible logical backaction
  achievable at a cost below the protected-memory advantage?
- Does the smaller 32-CP22 circuit remain preferable after **coherent**
  model mismatch and *measured* gate imperfections?

No pass on this first pair-channel probe implies a viable five-plaquette
memory; no fail disproves hypothetical analog gauge protection.

## Primary background (examples, not direct calibration evidence)

- *Universal Qudit Gate Synthesis for Transmons*, PRX Quantum 4,
  030327 (2023), DOI https://doi.org/10.1103/PRXQuantum.4.030327
- *Leakage in restless quantum gate calibration*, Physical Review A 108,
  022614 (2023), DOI https://doi.org/10.1103/PhysRevA.108.022614
- *Qudit Dynamical Decoupling on a Superconducting Quantum Processor*,
  Physical Review Letters 134, 050601 (2025),
  DOI https://doi.org/10.1103/PhysRevLett.134.050601
- *Extending the computational reach of a superconducting qutrit
  processor*, npj Quantum Information (2024),
  https://www.nature.com/articles/s41534-024-00892-z
