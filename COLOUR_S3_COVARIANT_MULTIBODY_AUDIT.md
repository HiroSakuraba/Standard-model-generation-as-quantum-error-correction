# Continuous gauge-covariant multibody controls versus gauge-breaking pair compilation

19 October 2026. Follows PR #50's **strict two-register, neutral-pointer obstruction** and PR #49's exact 16-factor (but temporarily gauge-breaking) compilation. This is a small finite \(\mathrm S_3\subset SU(3)\) testbed, **not** a construction of the nine-quark 21D code.

## 1. A genuinely gauge-invariant three-register Hamiltonian

On the two-vertex system from PR #47, let \(S\) be the **link-dressed exchange** of matter at A0 and B0, with the oriented six-state gauge link between them. It has support on exactly these **three physical registers**. We have already verified \(S^2=I\) and \([S,G_A(h)]=[S,G_B(h)]=0\) for every vertex transformation \(h\in S_3\).

Define the Hermitian three-register projector Hamiltonian

\[
H_S=\tfrac12(I-S), \qquad
U_S(t)=e^{-it H_S}=\tfrac{1+e^{-it}}2 I+\tfrac{1-e^{-it}}2 S.
\]

For every \(t\), \(U_S(t)\) commutes with **both** independent vertex gauge transformations. At \(t=\pi\), \(U_S(\pi)=S\). This *continuous* path never leaves the Gauss-invariant subspace when starting in that subspace. Similarly, the on-site swap \(A\) of qutrits A0 and A1 gives \(H_A=(I-A)/2\), a two-register gauge-invariant interaction, with \(U_A(\pi)=A\).

**Important limitation:** this is a mathematically defined Hamiltonian, not an engineered physical three-body coupling, native two-qutrit pulse, or a demonstrated analog control. Its coefficients, energy scale and availability on a real gauge simulator have not been established.

## 2. Exact logical-syndrome pulse sequence and register cost

PR #48's measurement uses a qutrit pointer to Fourier-measure the three-cycle \(U=AS\). This requires **two pointer-conditioned dressed exchanges** and **two pointer-conditioned on-site swaps**, then a pointer Fourier transform.

Conditioning \(H_S\) on a **gauge-neutral pointer** makes its control generator \(\Pi_p\otimes H_S\), with support on:

- pointer register,
- matter A0,
- oriented group link,
- matter B0.

That is a **four-register** interaction. Conditioning \(H_A\) touches three registers (pointer, A0, A1). Thus the ideal logical-syndrome circuit has **two four-register controlled pulse windows plus two three-register controlled pulse windows**, followed by single-pointer Fourier/readout. Conditional correction after a recognized dressed swap uses **one unconditioned three-register** \(U_S(\pi)\).

This explicitly does **not** evade PR #50's no-go using a gauge-preserving *two*-register native gate alphabet. A gauge-invariant three-register matter–link–matter operation resolves the **link-irrep mixing** obstruction; controlling its execution by a neutral pointer increases the necessary support by one register. No method of synthesizing the resulting four-body generator out of **gauge-preserving lower-body gates** is claimed.

## 3. Gauss-law path check

The code directly propagates each pulse in the full 486-dimensional system plus a three-level pointer; it does not start by applying logical projectors. For every pulse it checks fractions \(t/\pi=0,1/4,1/2,3/4,1\), on every encoded basis vector, against both independent six-element local gauge actions. It also verifies covariance for arbitrary 486D inputs. The continuous analytic commutator identity ensures preservation for all pulse times, not just the sampled fractions.

The earlier PR #49 two-register gate path returns to the physical sector at the boundaries of each high-level operation, but its intermediate snapshots have Gauss acceptance near zero. The present multibody Hamiltonian path keeps acceptance \(\simeq1\) at every tested fraction. This is a **path property**, not proof of lower gate cost, greater noise tolerance, or a feasible control amplitude.

## 4. Matched-noise comparison

The prior 16-factor pair circuit and new four-window multibody circuit represent different *hardware alphabets*. A count of **four** multibody controls versus **16** general \(U(9)/U(18)\) pair factors is **not a cost ranking**.

For a controlled but explicitly artificial diagnostic, we assume **exactly the same 16 notional units of syndrome time** in both cases and **exactly the same 3 additional units** of feedback/idle time. The four multibody windows last four units each; the sixteen pair factors last one unit each. The three-body feedback lasts three units; pair recovery uses three one-unit factors. The pointer preparation/Fourier, readout and idle routing are treated as ideal common overhead, not free on actual hardware.

At the end of each unitary segment, apply independent Gaussian random **computational-number phase rotations** to the pointer, six-level link, and four matter qutrit registers, with per-register variance \(2t\gamma\). The ensemble defines a **completely positive, trace-preserving (CPTP)** phase-diffusion channel. The physical dephasing channel is the *same* in both architectures, with three sensitivity rates \(\gamma=0,0.002,0.01\) per notional time unit. It is applied to the **entire system**, not merely a logical register. Each trajectory computes unconditional entanglement fidelity for the known dressed-swap channel \(\{\sqrt{1-p}I,\sqrt p\,S\}\), with \(p=0.37\), and feedback for reported pointer outcome 2. The pure unitary results must agree with the exact PR #48 instrument and yield fidelity 1.

**Interpretation boundaries:** the sample means reflect a toy *same-time number-dephasing* comparison, **not** a matched physical device. The dephasing faults on link/matter generally **do not commute with Gauss law**, even though ideal multibody control does. There is no T1, spontaneous decay, leakage, pulse waveform, time-dependent Lindblad solver, physical gauge penalty, coherent pulse overrotation, crosstalk, calibrated four-body Hamiltonian, or dynamical gauge enforcement. Monte Carlo standard errors quantify only trajectory sampling uncertainty. A more favorable fidelity does not establish an architecture or hardware advantage.

## 5. Test and next decision

The new suite checks: projector Hamiltonians and unitary time composition, covariance on arbitrary states for all local group elements, Gauss acceptance at all intermediate pulse fractions, exactly matched pointer circuit for arbitrary inputs, 3-versus-4-register support of elementary controlled controls, CPTP stochastic phase-noise norm preservation, perfect noise-free known-swap recovery, deterministic seeded nonzero-noise comparisons, and explicit scientific scope/cost reporting.

\`\`\`bash
python -m pip install numpy scipy
OPENBLAS_NUM_THREADS=1 python -m unittest -v test_colour_s3_covariant_multibody_control
OPENBLAS_NUM_THREADS=1 python colour_s3_covariant_multibody_control.py
\`\`\`

Next work should determine whether **a neutral-pointer-controlled gauge-invariant four-register Hamiltonian can be physically generated** with useful strength and noise characteristics. Competing approaches are gauge-charged mediators, perturbative three-body gadgets with explicit intermediate gauge penalties, and deliberately gauge-breaking pair gates whose leakage is controlled. Any final cost comparison needs **the same implementable native gate set, routing, time and channel parameters**. The original nine-quark 65-syndrome recovery problem remains unsolved.
