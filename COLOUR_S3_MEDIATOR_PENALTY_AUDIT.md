# Feasibility comparison: charged mediator, detuned gadget, and gauge penalty

19 October 2026. Continues merged PRs #49–51 of the finite non-Abelian S3-subgroup quantum-memory simulator. All reported results are synthetic and deliberately separate a **physically specified Hamiltonian** from a merely **mathematically exact Hamiltonian**.

## Problem statement

PR #51 generated a Gauss-covariant, three-register dressed exchange with H_S=(I-S_dressed)/2. A *neutral-pointer-controlled* syndrome gate still requires a four-register interaction. PR #50 proved that any strictly Gauss-preserving sequence of neutral-pointer two-register gates cannot implement it, regardless of gate count.

This PR tests three proposed escapes from that bottleneck. **No route yet provides a calibrated two-body/native device realization.** The results can therefore guide architecture selection without falsely claiming a solved hardware problem.

## 1. Gauge-charged mediator: allowed kinematics, incomplete circuit

Introduce a new three-level charge carrier M transforming as the signed-permutation SU(3)-subgroup fundamental under the local gauge group at vertex A. The simulator extends the exact local gauge action to the auxiliary mediator and constructs SWAP(M,A0), a **two-register** local exchange. Its covariance under both independently acting vertex gauge groups is tested for all six group elements on arbitrary complex states. This is an example of a local **gauge-covariant operation involving a charged mediator**, rather than assuming a costless neutral pointer.

It does **not** produce the logical three-cycle eigenphase syndrome, transport the charged mediator across the gauge link, or synthesize its controlled readout. Moving M to vertex B without violating Gauss law would require an explicit link-coupled charge-transfer operation and proper preparation/disposal of the gauge charge. The original PR #50 no-go does not automatically apply once charged mediator degrees of freedom are introduced, but it is not bypassed merely by adding one.

## 2. Detuned two-level mediator gadget: exact dynamics and resource tradeoff

Let Q=(I-S_dressed)/2 be the gauge-invariant **three-register** edge projector and P be a projector onto the selected qutrit pointer labels. Add a two-level neutral mediator m with gap Delta. Consider the exact Hamiltonian

    H = Delta |1><1|_m + g X_m (P+Q) + (g^2/Delta)(P+Q).

The last term is a counterterm. Since P and Q are commuting projectors, the system decomposes into three scalar sectors L=P+Q with eigenvalues 0,1,2, and each sector has an **exact 2×2 Rabi Hamiltonian**. The script evaluates its exact ground/excited mediator amplitudes rather than using an uncontrolled truncated perturbation series.

For Delta >> g, eliminating the excited mediator gives the low-energy effective Hamiltonian

    H_eff ≈ -(g^2/Delta)[(P+Q)^2-(P+Q)]
          = -(2g^2/Delta) P Q.

A pulse of duration t=pi*Delta/(2g^2) thus tends toward the **controlled dressed swap** because exp(i pi P Q)=I−2P Q. Reset and trace out the mediator after each of the four PR #48 syndrome windows. Enumerate all 2^4=16 mediator environment histories **without postselection** and compute unconditional recovered entanglement fidelity for the exact known-pair dressed-swap error channel. Record mediator excitation, coherent errors, and the pulse time versus detuning.

**Critical unresolved hardware issue:** the g X_m Q term acts jointly on the **mediator and the three-register edge** (FOUR registers). The construction trades a pointer-controlled four-register Hamiltonian for a mediator-controlled four-register Hamiltonian plus mediator–pointer interactions. It is *not* a reduction to two-body interactions. Real perturbative gadget reductions to bounded-strength two-body couplings, e.g. Bravyi, DiVincenzo, Loss and Terhal, Phys. Rev. Lett. **101**, 070503 (2008), DOI https://doi.org/10.1103/PhysRevLett.101.070503, require further ancillas, interaction graphs and explicit accuracy/strength analysis which are **not** provided here.

Increasing detuning can improve the weak-gadget approximation but increases syndrome time proportionally, even when ignoring T1/T2 and reset errors. No physical energy-scale calibration is assumed.

## 3. Gauge-penalty-assisted pair pulse: an exact local-constraint experiment

Take **one** of PR #49's gauge-breaking 9×9 two-qutrit pair factors, Q_A0,B0. Choose the Hermitian generator H_pair=i log Q_pair, and add two independent, vertex-local Gauss penalties:

    H(lambda) = H_pair + lambda[(I−P_A)+(I−P_B)].

P_A and P_B are the exact S3-averaged Gauss projectors at vertices A and B. Their individual support includes a link plus two vertex matter qutrits (THREE registers); the two terms do not use a nonlocal *product* projector as the physical penalty.

The script propagates the **complete 486-dimensional register** for a fixed pulse time with sparse-matrix exponential action and measures Gauss return at each vertex, joint physical return probability and deviation from the intended unpenalized factor for lambda=0,2,8,32. The same pulse time means large lambda may reject gauge-violating transitions at the cost of distorting the desired gate. **This is not a complete 16-factor penalized syndrome and gives no logical recovery fidelity for the penalty route.** It also does not account for building vertex-local Gauss penalties from native two-body interactions.

Related primary literature: Ball and Cohen, *Zeno effect suppression of gauge drift in quantum simulations*, Phys. Rev. A **110**, 022417 (2024), DOI https://doi.org/10.1103/PhysRevA.110.022417; *Protecting gauge symmetries in the dynamics of SU(3) lattice gauge theories*, Communications Physics (2025), https://www.nature.com/articles/s42005-025-02230-x. These provide context, not a calibration or validation of our S3 pulse.

## Decision criteria and open problems

| Route | Gauge invariance | Missing engineering | Key tradeoff |
|---|---|---|---|
| Charged mediator local swap | Exact for demonstrated vertex-A exchange | Link-mediated charge transfer and logical readout | Gauge charge, link connectivity, ancilla preparation |
| Detuned two-level gadget | Exact Hamiltonian commutes with both Gauss groups | Four-register mediator-edge interaction still assumed | Larger detuning reduces gadget errors but lengthens storage |
| Gauss-penalty pulse | Vertex penalties gauge invariant but added pair drive violates gauge symmetry | Native construction of penalties and whole compiled syndrome | Suppressed off-sector dynamics versus wrong target gate |

None should be promoted to a usable quantum memory on the basis of a toy-simulator result. The **next hard step** is either synthesizing a *genuinely lower-body* gauge-covariant charged-mediator transport and logical readout, or computing a **whole-circuit pulse-time** benchmark with specified implementable native gates and penalties, including environment noise.

## Reproduce

    python -m pip install numpy scipy
    OPENBLAS_NUM_THREADS=1 python -m unittest -v test_colour_s3_mediator_penalty_feasibility
    OPENBLAS_NUM_THREADS=1 python colour_s3_mediator_penalty_feasibility.py

CI archives a JSON report containing all actual numerical results and explicit limitations. Numerical eigenvalue/exponential residuals verify this implementation, not a computer-assisted physical feasibility theorem.
