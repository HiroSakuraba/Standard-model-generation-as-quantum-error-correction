# Compiling one interior SU(3) square-chain magnetic interaction to nearest-neighbor qutrit gates

**18 October 2026 · Status: constructive logical compilation and phenomenological 4-level error simulation.**

Runnable code: [`su3_selective_cp_compiler.py`](su3_selective_cp_compiler.py). Machine output: [`results.json`](results.json), [`output.txt`](output.txt). Reference-table dependency: [`su3_three_plaquette_chain.py`](su3_three_plaquette_chain.py).

## Question and decisive result

The prior SU(3) square-chain studies treat an interior plaquette's neighboring qutrits as *controls*, but did not implement its neighbor-conditioned three-qutrit magnetic action in a two-qutrit gate set. Can one compile that particular three-qutrit unitary exactly using only nearest-neighbor two-qutrit entanglers, local operations, and no clean ancilla?

**Yes, within a specified idealized gate set.** The implementation explicitly compiles the nine neighbor-conditioned SU(3) blocks into adjacent-level SU(2) operations, synthesizes doubly controlled rotations from commutators of two-body controlled gates, converts each two-body controlled rotation into level-selective conditional phases, then compiles each level-selective phase into one fixed-level `CP_22(phi)` using adjacent-level qutrit permutations. All entanglers couple only left–center or center–right. There is no native three-body gate, no direct left–right link, and no ancillary qutrit in this circuit.

An independent full 27×27 unitary simulation gives an operator-norm reconstruction error approximately **9.76×10^-15** for one full interior-plaquette step at a dimensionless duration **0.19**, magnetic phase **2π/3**. This number is floating-point verification of a constructive analytic circuit identity; it is **not a numerical tolerance certified with interval arithmetic**.

The **unoptimized** exact circuit uses **360 fixed `CP_22` two-qutrit gates**, **1,728 local qutrit operations** (of which **1,368 are adjacent-level swaps**), for all nine neighbor sectors. A *single selected neighbor-conditioned block* uses **40 `CP_22` gates** and **152 local operations** (112 swaps, 40 eigenbasis transformations). These counts are a constructive upper bound **for this decomposition**, not a lower bound on the optimal circuit or a calibrated gate schedule. Local eigenbasis transformations are adjacent-level `U(2)` mixers and may need additional physical drive/virtual-Z primitives; the code treats each as one *abstract local operation*.

## Source SU(3) Hamiltonian and independent source continuity

Vincent Chen, Berndt Müller, and Xiaojun Yao, *Minimally truncated SU(3) lattice gauge theory and string tension*, **Physical Review D 114**, 074502 (2026), DOI [10.1103/1qt7-3ywl](https://doi.org/10.1103/1qt7-3ywl), published Appendix B **Table VI**. In the zero-exterior-flux square-chain, let `(l,r)` be the qutrit labels of the plaquettes neighboring the interior face. Its published conditional magnetic block is a real symmetric three-by-three matrix `M_(l,r)`, with entries from `1,1/3,1/9,1/sqrt(3),1/(3sqrt(3))` and zero diagonal. Write

\[
W_{(l,r)}=\sum_{a=0}^{2} [M_{(l,r)}]_{a+1\bmod3,a}\,|a+1\rangle\langle a|,
\quad B_{(l,r)}(\varphi)=e^{i\varphi}W_{(l,r)}+e^{-i\varphi}W_{(l,r)}^\dagger.
\]

The desired **three-qutrit** unitary on left, center, right is

\[
U_{\mathrm{face}}(\tau,\varphi)=\sum_{l,r=0}^{2}|l\rangle\langle l|_L\otimes
e^{i\tau B_{(l,r)}(\varphi)}\otimes|r\rangle\langle r|_R.
\]

This is the magnetic term *alone* at one interior face, not a full time step of the five-face electric-plus-magnetic Hamiltonian. The script checks its table entries against the **previous independent published transcription** in `su3_three_plaquette_chain.py`, as well as performing a separate direct 27×27 exponentiation of the full block-diagonal generator.

## Exact compilation identities

First, `tr B_(l,r)=0`, so every conditional `U_(l,r)=exp(iτ B_(l,r))` lies in `SU(3)`. Adjacent-level Givens elimination decomposes it into at most **three** rotations on `(1,2),(0,1),(1,2)` plus **two** diagonal `SU(2)` phase rotations on adjacent subspaces `(0,1)` and `(1,2)`. There is **no conditional overall U(1) phase** and hence no hidden direct left–right conditional-phase gate.

For any target `R∈SU(2)`, diagonalize it to choose `B∈SU(2)` satisfying `B²=R†`. Choose `A∈SU(2)` exchanging B's two eigenvectors, so `A B A†=B†`. Then

\[
A B A^\dagger B^\dagger=R.
\]

The chronological gate sequence `C_R(B†), C_L(A†), C_R(B), C_L(A)` implements **exactly** the double-control operation `C_{L=l,R=r}(R)`, since gates on a control sector where the control is not selected cancel. Each `C_L` or `C_R` acts on only the control and middle qutrit. For a two-level `S=V diag(e^{iα},e^{iβ}) V†`, the **two-body** controlled-S is `V†`, `CP_(control=l,target=a)(α)`, `CP_(control=l,target=b)(β)`, `V` in time order; outside the selected control state the local `V` factors cancel. Thus one double-control SU(2) uses **eight selective two-qutrit phases** and **eight local basis operations**.

Finally each arbitrary-level phase `CP_(control=s,target=t)(φ)` becomes one fixed `CP_22(φ)` by swapping logical levels `s→2` and `t→2` using only adjacent `(0,1)` or `(1,2)` swaps, applying `CP_22`, and undoing them. No nonadjacent drive is needed by this routing. This last stage assumes a **continuously tunable angle** for the fixed conditional phase primitive, which has **not been established for the specific proposed hardware**. The code verifies both the level-selective and fixed-level versions **for the whole 27-state space**, not just for their intended control sector.

## Explicit four-level transmon stress model

Three physical transmons are represented by **four energy levels apiece**, yielding a 64-dimensional Hilbert space; the computational subspace has 3³=27 states. Each fixed `CP_22` includes a **phenomenological** coherent mixing of the computational `|22>` pair with leakage `|31>` on that coupled pair, plus a detuning on the fourth-level partner. An independent spectator-edge conditional phase is applied after each entangler to represent systematic crosstalk. A fractional coherent bias multiplies the intended conditional-phase angle. All simulated operations are full unitaries on the 64-dimensional space, preserving total norm; leakage is measured as probability outside the 27-state computational sector.

Representative results for the complete nine-context interaction:

| Illustrative error parameters | Mean computational leakage | Unconditional entanglement fidelity |
| --- | ---: | ---: |
| None | zero within precision | 1.000000 |
| Coherent leakage amplitude 0.01 **and** off-target phase 0.002 per entangler | **0.0009795** | **0.991026** |

For the selected `(1,2)` block, leakage amplitudes 0.01 and 0.03 alone give leakage probabilities approximately **0.000106** and **0.000948**. A combined `(0.03, 0.006, 2% phase bias)` test has unconditioned process fidelity approximately **0.997493** for that selected block.

These parameters are **dimensionless inventions**. The leakage matrix element and detuning are not derived from a physical transmon's anharmonicity, coupler, pulse shape, bandwidth, flux noise, relaxation, or measured crosstalk. The apparently high simulated fidelities do not include T1/T2 decoherence or stochastic infidelity of hundreds of native gates. They must **not** be presented as hardware predictions.

## Hardware feasibility and meaningful costs

We chose this gate set because experimental superconducting-qutrit work has demonstrated **two-qutrit entangling cross-Kerr interactions**, not because the exact fixed `CP_22(φ)` primitive is known to be available with arbitrary angle and high fidelity. Relevant primary results:

- Luo *et al.*, *Experimental Realization of Two Qutrits Gate with Tunable Coupling in Superconducting Circuits*, **Phys. Rev. Lett. 130**, 030603 (2023), DOI https://doi.org/10.1103/PhysRevLett.130.030603. The reported two-qutrit gate fidelity is **89.3%**; its specific implemented operation should not be confused with the unlimited arbitrary-angle primitive assumed here.
- *High-fidelity qutrit entangling gates for superconducting circuits*, **Nature Communications** (2022), DOI https://doi.org/10.1038/s41467-022-34851-z. Reports a specific two-qutrit CZ with **580 ns** and related CZ† with **783 ns** gate durations, not arbitrary-parameter fixed-`|22>` conditional phases.
- *Exploring ququart computation on a transmon using optimal control*, **Phys. Rev. A 108**, 062609 (2023), DOI https://doi.org/10.1103/PhysRevA.108.062609. Experimentally demonstrates four-level operation, without validating this three-device leakage model.

As a deliberately explicit **timing scenario**, suppose that each ideal arbitrary-angle `CP_22` consumes 580 ns and each abstract local mixer/swap consumes 30 ns. Then **one selected `(1,2)` block** requires `40×580 +152×30=27,760 ns =27.76 µs`. The **full interior unitary** takes `360×580+1728×30 =260,640 ns =260.64 µs` in a completely serial schedule. This is *not* a proven duration, minimum or calibration. Many local operations are only symbolic U(2) rotations, and 580 ns refers to a different demonstrated gate. The long modeled time suggests that decoherence and calibration costs could overwhelm any benefit; this study does **not** claim a feasible quantum memory.

There may be far shorter compilations using simultaneous controls, direct analog plaquette interactions, multi-qutrit optimal control, alternative native gates, or controlled unitary factoring. No optimality theorem is offered.

## Scientific decision and release criteria

The **conceptual compilation obstacle is solved** for the declared abstract two-body gate set, but the **physical hardware feasibility remains unresolved**. The next sensible work is not more circuit repetition: optimize one interior interaction's two-body entangling count, restrict to measured Hamiltonian terms and available coupler connectivity, include calibrated T1/T2 and leakage/readout statistics, and implement a complete five-face refocusing cycle to see whether the QEC advantage survives realistic time and errors.

A physically supportable claim requires (a) a measured calibrated variable-angle `CP_22` implementation or a verified decomposition into actual native gates, (b) verified single-qutrit adjacent-transition pulse synthesis, (c) calibrated leakage/crosstalk and timing model, and (d) compiled gauge-invariant syndrome recovery. A native 3-level face clock still does not imply a native SU(3) gauge **link** control. This work proves no physical glueball mass, continuum Yang–Mills gap, Standard Model generation code, or fault-tolerance threshold.

## Reproduce

```bash
python -m pip install numpy scipy sympy
OPENBLAS_NUM_THREADS=1 python su3_selective_cp_compiler.py --output results.json
python -m py_compile su3_selective_cp_compiler.py
```

`su3_three_plaquette_chain.py` is included so the code can compare the published magnetic blocks against the earlier independent source transcription. The program includes assertions for each SU(2) commutator, Givens reconstruction, full 27-state 3-body unitary, fixed-`|22>` conversion, and 4-level unitary/no-leakage baseline. All tested parameters and results are saved in the companion machine-readable JSON.
