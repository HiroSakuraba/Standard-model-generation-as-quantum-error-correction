# First numerical pulse compiler for the finite-S3 three-outcome quantum recovery

**19 October 2026 — research experiment, not a practical device claim.**

This report follows [PR #65](https://github.com/HiroSakuraba/Standard-model-generation-as-quantum-error-correction/pull/65). That PR proved that a specifically supported set of finite-\(S_3\), gauge-invariant Hamiltonians generates the complete physical \(\mathfrak u(14)\) control algebra and, with a neutral qutrit pointer, reaches \(\mathfrak{su}(42)\) in principle. The logical next step is a **finite, numerical pulse-sequence compiler** for the three-outcome recovery from [PR #64](https://github.com/HiroSakuraba/Standard-model-generation-as-quantum-error-correction/pull/64).

## 1. What is actually compiled?

The original model has four matter-colour qutrits, an oriented finite-group \(S_3\) link, and **14 states satisfying the two independent Gauss laws**. A neutral three-level pointer extends this physical space to **42 dimensions**.

PR #63 supplies the four-dimensional logical code \(C\) and the four-dimensional correctable electric-syndrome space \(V\) for each of the two order-six matter-exchange groups \(AS\) and \(BS\). The remaining orthogonal **six-dimensional subspace** must be flagged, not silently normalized out.

The ideal three-outcome instrument from PR #64 acts on pointer initialized to \(|0\rangle\):

\[
C\otimes|0\rangle\mapsto C\otimes|0\rangle,\quad
V\otimes|0\rangle\mapsto C\otimes|1\rangle,\quad
P_{\rm rest}\otimes|0\rangle\mapsto P_{\rm rest}\otimes|2\rangle.
\]

The target is an **isometry on all fourteen pointer-ready physical states**; the code also constructs the complete, specified \(42\times42\) unitary extension \(U_{\rm target}\). **The default optimizer only fits those fourteen relevant input columns.** Its report also evaluates the achieved **full 42D unitary fidelity** independently. Fitting the initial-pointer-0 subspace must *not* be called an exact compilation of the entire global unitary.

## 2. Thirteen independently tunable gauge-compatible controls

Each cycle contains exactly the following **13 elementary control rotations**:

- Four physical-system Hamiltonians proven gauge compatible in PR #65: local on-site matter exchange \(A\) (or \(B\)); dressed cross-link matter exchange \(S\); single-link electric \(E=P_{\rm standard}\); and the **opposite-vertex** local two-matter equal-colour projector \(Q_B\) (or \(Q_A\)).
- Eight neutral-pointer single-register Hermitian \(SU(3)\) generators (three X-like, three Y-like and two independent Z-like controls).
- One assumed tunable **two-register pointer+link** interaction \(|1\rangle\langle1|_{\rm pointer}\otimes P_{\rm standard}(\mathrm{link})\).

All gates are represented on the normalized 14D Gauss slice × neutral qutrit pointer. For a control \(H_j\) and angle \(\theta_j\),

\[
U_j(\theta_j)=e^{-i\theta_j H_j},\qquad
U(\boldsymbol\theta)=U_N(\theta_N)\cdots U_1(\theta_1).
\]

Each Hamiltonian is diagonalized once to evaluate variable-angle pulses without approximating a matrix exponential on each optimizer call. The compiler uses **bounded L-BFGS-B**, deterministic random starts, the identity start and warm-starting larger circuits from smaller results. It archives the full reproducible angle list, optimizer termination status and all three branch metrics.

### Honest objective and gradient

Let \(\Phi_0\) be the 42×14 isometry initializing the pointer to \(|0\rangle\), and \(Y=U_{\rm target}\Phi_0\). The default objective is the *coherent phase-sensitive-together* target overlap, with one common irrelevant global phase:

\[
J(\boldsymbol\theta)=
1-\left|
\frac{\operatorname{Tr}(Y^\dagger U(\boldsymbol\theta)\Phi_0)}{14}
\right|^2.
\]

The code differentiates every rotation analytically,

\[
\frac{\partial U_j}{\partial\theta_j}=-iH_jU_j,
\]

and computes the full gradient using forward states and backward-adjoint states. Five randomly chosen coordinate derivatives are independently checked by centered finite differences. An optional full-42D unitary objective exists for diagnostic comparisons, but the depth sweep defaults to the more physically relevant 14-input-column problem.

**This is a meaningful coherent target-isometry score**, but it is not automatically the *unconditional quantum-memory entanglement fidelity* under a particular stochastic electric/exchange noise mixture. The tool separately evaluates code, electric and six-dimensional flagged input branches, their intended pointer-label probabilities, off-label probability, and target-overlap fidelity without postselection.

### Separate operational logical-channel fidelity, not just target overlap

To test whether the ACTUAL compiled pulse sequence corrects anything, apply to the same four-dimensional code C each of the **seven admitted physical error unitaries** from PR #63: all six individual exchange permutations \(g\) of one chosen subgroup, plus the electric phase \(R=I-2E\). Initialize the pointer to \(|0\rangle\), apply the real compiled 42D pulse product, and read the three pointer outcomes.

For recognized pointer outcomes 0 and 1, project the *physical output* back onto the intended code C. Each output defines an unnormalized logical Kraus matrix \(K_j=C^\dagger \langle j| U_{\rm compiled}(F_i C\otimes|0\rangle)\). The **unconditional recovered logical entanglement fidelity** for fault \(F_i\) is

\[
F_e(F_i)=\sum_{j=0,1}\frac{|\operatorname{Tr} K_j|^2}{4^2}.
\]

Probability appearing in **pointer outcome 2** or in physical states **outside C even with pointer 0 or 1** is counted as loss of the logical state. No renormalization or postselection is allowed. The output archive reports each fault separately, all code-capture and flag probabilities, the worst individual-fault fidelity, and a synthetic **normalized seven-operator stochastic channel** with electric-phase probability \(p_e=0.37\) and uniform probability \((1-p_e)/6\) over the six subgroup permutations.

This channel is a controlled mathematical test, NOT the device-correlated Bernoulli fault mixture from PR #64 and NOT a measured physical distribution. The ideal global target must obtain **F_e=1 for all seven faults**; the do-nothing identity circuit obtains exactly \(F_e=(1-p_e)+p_e/9\), or about **0.671111** at \(p_e=0.37\), because the electric fault has scalar overlap \(-I/3\). Both analytical reference values are asserted in dedicated tests.

The new **operational channel F_e** is more relevant to correctness of the approximate recovery than coherent isometry overlap, but even it excludes pulse noise, virtual Hubbard leakage, charge/reference loss, pointer assignment error and any faults during recovery. It is not a physical quantum-memory demonstration.

## 3. Finite depth sweep: report failures rather than assuming success

The default experiment fits both \(AS\) and \(BS\) recovery targets at **1, 2, 4, 8 and 16 cycles**, or respectively **13, 26, 52, 104 and 208** tunable rotations, with a bounded iteration budget. It tries identity, reproducible small random angles, and a warm start from the preceding depth when available. It always keeps the best verified result, including the unoptimized identity schedule as a reference.

This intentionally does **not** assume even the 208-pulse depth reaches a satisfactory recovery. The research output prints and records: coherent input-isometry fidelity, entire-42D unitary fidelity, per-branch target overlaps, intended-versus-wrong pointer probabilities, actual angles and optimization termination. A result not sufficiently close to the ideal target is a **failed or incomplete compilation at that depth**, regardless of the exact \(\mathfrak{su}(42)\) controllability theorem from PR #65.

The choices of cycle ordering, angle limits and optimization iterations are **search heuristics**, not lower bounds on any best possible circuit.

## 4. Expand every dressed exchange to charged-mediator pair operations

The dressed cross-link \(S\) initially couples a gauge link and two matter qutrits. PR #65 established that its fractional pulse can be performed by exactly five **abstract two-register gauge-covariant charged-mediator steps**:

\[
\mathrm{SWAP}(m,A0),\quad T(m,\mathrm{link}),\quad
e^{-i\phi(I-\mathrm{SWAP}(m,B0))/2},\quad
T(m,\mathrm{link}),\quad\mathrm{SWAP}(m,A0).
\]

For a compiler control \(e^{-i\theta S}\), take \(\phi=-2\theta\); this differs from the dressed-exchange pulse only by a common global phase and gives the correct relative logical phases. The simulator expands each nonzero dressed-\(S\) rotation to these **five** primitive entries in the JSON archive, with two transport operations per dressed edge control. Other angles are recorded as one assumed tunable one- or two-register gate each.

The resulting estimate counts **actual control rotations**, **expanded one/two-register operations**, **charged-link transports** and **serial nominal elapsed time**. Rates for all controls and the readout are independently configurable; the default values (10 MHz, 100 ns/transport, 300 ns/readout) are **hypothetical sensitivity assumptions only**. The two-register charged transport itself, the fractional mediated exchange interaction, the opposite-vertex \(S_3\) pair-colour projector and the pointer+link projector are **not native hardware gates demonstrated on a device**.

As an explicit optimistic diagnostic, the report also multiplies the achieved isometry target score by a uniform four-active-fermion no-jump factor \(e^{-4T/T_1}\) with an assumed \(T_1=30\) μs. **This must not be called the physical fidelity**: it omits spectator/reference loss, the extra charged mediator during transports, pulse dephasing, stochastic correlated faults, compilation leakage, pointer readout errors and repeated jump dynamics.

## 5. Quality and correctness checks

The regression suite checks, for **both** selected gauge-exchange subgroups where relevant:

1. The thirteen 42D control Hamiltonians are Hermitian.
2. The exact target is 42D unitary, with clean/electric/flag ranks 4/4/6.
3. Every variable-angle control pulse is unitary.
4. The analytic adjoint gradient matches finite-difference derivatives, in both 14-input and full-42D objective modes.
5. The identity circuit does **not** spuriously report perfect recovery or a valid flagged outcome.
6. The bounded optimizer reports reproducible, finite angles with no regression relative to the identity baseline.
7. The 5-step bus expansion is explicitly counted and its dressed-edge phase angle is \(-2\theta\).
8. Nonphysical pulse rates or time parameters are rejected.
9. Reproducible seeded optimization and smaller-depth warm starts work.
10. Actual seven-error recovered logical entanglement fidelity is evaluated without postselection, with ideal-target F_e=1 and identity-circuit F_e=(1-p_e)+p_e/9 regression baselines.
11. The result retains all physical limitations in the JSON audit rather than asserting a hardware-native recovery.

The exact finite-\(S_3\) mathematical controllability certificate remains in PR #65; this PR addresses **numerical constructive efficiency**.

## 6. What counts as a scientific win?

The short-term win is an **honestly reproducible finite circuit**, even if initially poor, accompanied by an explicit residual-versus-depth comparison and strict separation of intended inputs, flagged leakage and full-unitary performance.

The longer-term objective is a circuit whose **unconditional noisy logical recovery** exceeds optimized equal-hardware passive storage when both are subjected to the same device-derived fault process, with realistic native transport, pointer coupling, charged-reference survival, intermediate doubly occupied virtual Hubbard states and recovery cadence.

**No result here establishes actual full \(SU(3)\) gauge-code recovery or the original nine-quark, 21-dimensional/65-syndrome implementation.** The opposite-vertex equal-colour projector is a legal *discrete* \(S_3\) Hamiltonian but is not invariant under arbitrary continuous \(SU(3)\). The virtual gluon, sea-quark and QCD vacuum sectors still remain absent.

### Reproduce

\`\`\`bash
python -m pip install numpy scipy sympy
OPENBLAS_NUM_THREADS=1 python -m unittest -v test_colour_s3_recovery_pulse_compiler
OPENBLAS_NUM_THREADS=1 python colour_s3_recovery_pulse_compiler.py
\`\`\`

The dedicated workflow archives full per-depth schedules and angle lists, independent analytic-gradient checks, every pointer outcome, an unexpanded-versus-expanded gate cost ledger, and explicitly hypothetical noise-duration estimates.
