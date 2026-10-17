# Exact recovery of one local exchange fault using five SU(2) singlet modules

**17 October 2026 — research status:** a reproducible *positive* Knill–Laflamme result for a carefully specified modular SU(2) matter/rishon model. It uses the **existing, published five-qubit perfect code** as an outer error-correcting code. It is **not** a new code from Standard Model anomaly cancellation, an implementation of the full Standard Model, or a physical syndrome-extraction circuit.

## Executive statement

Five **independent copies** of the four-link matter/rishon SU(2) module developed in [FOUR_LINK_SU2_DFS.md](FOUR_LINK_SU2_DFS.md) each provide one logical qubit within a fixed link-end configuration and a two-dimensional color-singlet subspace. An individual module is **not** protected against generic local SU(2)-invariant color exchange; two distinct swaps generate all of \(M_2(\mathbb C)\) on that singlet pair.

However, a standard **\([[5,1,3]]_2\)** outer code across five such modules protects one global logical qubit against **arbitrary errors acting on one unknown module**, given coherent preparation, gauge-compatible syndrome measurement, and recovery operations. This includes errors built from noncommuting local color exchanges and endpoint occupation fluctuations. The result is exact, not a perturbative or mean-noise-rate estimate.

The gap from earlier work is clear: one module is a *decoherence-free subspace only for a narrow occupation noise algebra*; five independent modules plus active cross-module checks provide **actual correction of a complete single-module logical error algebra**.

## 1. Physical inner modules and redundancy

One module is the four-parallel-link SU(2) kinematic truncation, with:

- Two sites, one spin-half **matter color doublet at each site**.
- Four distinguishable links, one spin-half **rishon** per link located at one endpoint.
- A fixed end assignment \((0,1,1,1)\), giving a two-dimensional product of local SU(2) color-singlet sectors.
- Explicit orthonormal inner vectors \(|0\rangle_\mathrm{in}\) and \(|1\rangle_\mathrm{in}\), constructed exactly in [\`four_link_su2_subspace.py\`](four_link_su2_subspace.py) using valence-bond singlets.

Five independent modules amount, in this kinematic model, to **10 matter doublets and 20 rishon doublets, associated with 20 four-link connections**. The fixed-endpoint, gauge-invariant inner Hilbert space is \((\mathbb C^2)^{\otimes5}\), dimension **32**, inside the much larger physical singlet sector of all five modules. The gauge laws of distinct modules commute.

The earlier [SU(2) exchange-algebra analysis](SU2_EXCHANGE_HAMILTONIAN_AND_NO_GO.md) shows how co-located, SU(2)-invariant swaps act on the inner qubit. At the rich four-spin vertex, with spin labels \(1,3,4,5\) and Pauli matrices \(X_\mathrm{in}\), \(Z_\mathrm{in}\),

\[
\begin{aligned}
S_{13}&=-Z_\mathrm{in},\\
S_{14}&=(Z_\mathrm{in}-\sqrt3 X_\mathrm{in})/2,\\
S_{15}&=(Z_\mathrm{in}+\sqrt3 X_\mathrm{in})/2.
\end{aligned}
\]

In particular,

\[
Z_\mathrm{in}=-S_{13},\qquad
X_\mathrm{in}=(S_{15}-S_{14})/\sqrt3.
\]

Every \(S_{ij}\) swaps two **co-located** fundamental color doublets and therefore commutes with the local non-Abelian Gauss generators. The verifier imports the *actual projected matrices* from the SU(2) module and checks these identities exactly; it does not simply assume Pauli operators acting on fictional matter bits.

## 2. The known five-qubit code

We apply the **1996 five-qubit perfect code** of Laflamme, Miquel, Paz and Zurek, [*Perfect Quantum Error Correcting Code*, Physical Review Letters **77**, 198](https://doi.org/10.1103/PhysRevLett.77.198), [arXiv:quant-ph/9602019](https://arxiv.org/abs/quant-ph/9602019).

Here each position of that established code is an **entire gauge-invariant SU(2) inner module**.

Its four commuting, independent stabilizer generators are

\[
\begin{array}{ccccc}
X&Z&Z&X&I\\
I&X&Z&Z&X\\
X&I&X&Z&Z\\
Z&X&I&X&Z.
\end{array}
\]

These \(X,Z\) are *inner singlet-qubit operators*, not bare spin-half Pauli operators on color-charged matter. A representative logical \(X\) is \(XXXXX\), and a logical \(Z\) is \(ZZZZZ\).

Each generator is a **gauge-invariant observable when restricted to the chosen inner-module Hilbert spaces**, because it consists of products of operators that commute with local SU(2) gauge transformations. But its nonlocal, multi-module measurement and compatible fault-tolerant gates have **not** been supplied. The stabilizer-code theorem is an algebraic existence statement conditional on these operations.

### Concrete encoded states

The pure-Python exact verifier generates codewords using the commuting positive stabilizer projectors:

\[
|\overline0\rangle=\tfrac14
(I+g_1)(I+g_2)(I+g_3)(I+g_4)|00000\rangle.
\]

In this conventional representation each resulting codeword has **16 nonzero signed amplitudes**, exactly \(\pm1/4\); the companion state is \(|\overline1\rangle=XXXXX|\overline0\rangle\). The script verifies that the two vectors are orthonormal and fixed by each \(g_i\).

Because the module basis is made of color singlets, **every** superposition of these five-module codewords remains a gauge singlet on every matter site. No gauge-charge superselection sector is mixed.

## 3. Exact one-module Knill–Laflamme proof

Consider the following **16** error basis elements:

\[
\mathcal E=\{I\}\cup
\{X_j,\ Z_j,\ X_j Z_j: j=0,\dots,4\}.
\]

The \(XZ\) matrix is the real-orthogonal representative of the usual \(Y\) Pauli error up to a global factor \(i\); this phase does not change quantum recovery or syndromes.

Every error has a distinct **four-bit syndrome** under the four outer checks. The script checks all 16 values and explicitly constructs an inverse-error correction lookup table.

Writing \(|u\rangle\) for either encoded state, the verifier establishes **all 1,024 exact scalar identities**

\[
\boxed{\langle u|E^\dagger F|v\rangle
=\delta_{EF}\delta_{uv}}
\quad
(E,F\in\mathcal E,\ u,v\in\{0,1\}).
\]

There are \(16^2\cdot2^2=1,024\) such matrix elements. Every amplitude is an exact rational number and every comparison is exact equality—no floating point, Monte Carlo or unimplemented interval tools.

This is the Knill–Laflamme condition \(PE^\dagger FP=c_{EF}P\) in its nondegenerate perfect-code form. All **32 states** \(E|\overline u\rangle\) are orthonormal, filling the five-module 32-dimensional logical register. For every single-module channel with Kraus operators expanded in \(\mathcal E\), measuring the four syndrome bits and applying the corresponding inverse Pauli recovers the encoded state.

**Important noise scope:** The error operators need not be Pauli *in nature*: the operator basis spans all \(2\times2\) matrices on any **one** inner singlet qubit, so arbitrary *coherent* single-module exchange perturbations are covered exactly. The scheme does **not** correct simultaneous independent errors on two modules, correlated intermodule exchanges, or faulty syndrome extraction without further analysis.

### A stronger mathematical extension: module leakage

Let \(V:\mathbb C^2\to\mathcal H_\mathrm{module}\) embed an inner singlet qubit into its full physical module Hilbert space, and let \(W=V^{\otimes5}W_{\mathrm{outer}}\).

For two arbitrary one-module physical errors \(E_i,F_j\), project \(W^\dagger E_i^\dagger F_j W\) into the outer five-qubit register:

- If \(i=j\), the result is the projection \(V^\dagger E_i^\dagger F_iV\), an *arbitrary single-qubit operator*. The distance-three outer code restricts it to a scalar.
- If \(i\ne j\), it is the tensor product of \(V^\dagger E_i^\dagger V\) on module \(i\) and \(V^\dagger F_j V\) on module \(j\), an operator with support at most two blocks. The distance-three outer code detects its nonidentity part.

Thus the mathematical outer code actually satisfies Knill–Laflamme for **arbitrary errors acting on at most one physical module**, even ones that leak outside the chosen two-dimensional singlet block, as long as errors are single-module operators and *ideal global recovery* is allowed. The provided executable checks cover the complete Pauli basis **inside** the 32-dimensional inner-register space. The leakage extension is an algebraic corollary, not a simulated SU(2) leakage/recovery circuit. If only gauge-invariant operations are allowed, their detailed implementation remains to be established.

## 4. Why five modules?

With each SU(2) singlet block used strictly as a **two-level subsystem**, correcting arbitrary unknown-location single-block errors requires an encoded \([[n,1,3]]_2\) code. The quantum Singleton bound states

\[
n-k\ge 2(d-1).
\]

For \(k=1,d=3\), this requires \(n\ge5\). The standard five-qubit code achieves the bound and is the **smallest number of two-level modules for this full one-module-noise requirement**.

This minimum claim does **not** exclude using fewer physical modules while exploiting their larger, 16-dimensional gauge-singlet sectors, different error constraints, subsystem encodings or approximate correction. The five-module construction is optimal only within the declared two-level-module architecture.

One can use three blocks for **pure one-axis phase errors**, but that biased-noise repetition construction fails for the *other* noncommuting exchanges already identified inside our SU(2) inner blocks. The five-module outer code closes this particular gap.

## 5. What has, and has not, been achieved

**Mathematical achievement:**

- Explicit gauge-singlet inner states and physically meaningful on-site SU(2)-invariant exchange operators.
- A full gauge-invariant multi-module encoded qubit **at the level of an abstract projective subspace**.
- A complete **exact Knill–Laflamme recovery certificate** for any one failed module within the fixed-orientation logical register.
- A general single-module-leakage correction theorem by composition with the standard five-qubit code.

**Still missing for a physical implementation:**

1. A **gauge-compatible protocol to prepare the highly entangled five-module codewords** from local matter/rishon states.
2. Realizable, non-demolition measurements of the four **cross-module stabilizers** without destroying the internal singlets or spreading one error across multiple blocks.
3. A gauge-invariant, controllable recovery circuit with a realistic noisy measurement/gate model.
4. A native matter/rishon Hamiltonian that suppresses higher-weight and correlated intermodule errors during storage and correction.
5. A route from this deliberately truncated SU(2) model to the actual Standard Model \(\mathrm{SU}(3)_c\times\mathrm{SU}(2)_L\times\mathrm U(1)_Y\) interactions and flavor mixing.
6. Any connection between these **outer five-qubit stabilizers** and the cubic anomaly-cancellation-derived charge matrices studied in the rest of this repository.

**Novelty boundary:** The recovery mechanism is the **known five-qubit code concatenated with an SU(2) singlet embedding**. It is a legitimate, explicit existence construction for the stipulated mathematical model, but not evidence of a naturally emergent code in real matter.

## 6. Reproduce

\`\`\`bash
python -m pip install sympy==1.14.0
python five_block_su2_recovery.py
\`\`\`

The new script is registered in the GitHub Actions research workflow and runs after the earlier SU(2) results. It performs all exact basis/recovery checks without creating a \(1024^5\)-dimensional matrix; the very large physical Hilbert space enters only through a mathematically explicit five-fold **isometric embedding** of five 2D singlet blocks.

## 7. Next scientifically useful milestone

The immediate challenge is no longer the Knill–Laflamme condition—it **passes** in this modular model—but the **physicality of syndrome extraction and error locality**. The next step should design a cross-module stabilizer-measurement gadget expressed in the color-singlet exchange basis, analyze which required interactions couple distant SU(2) matter sites, and quantify whether a *single* faulty extraction operation can propagate into two modules. If it can, the [[5,1,3]] outer code does not guarantee fault tolerance; extra verification or flags will be necessary.

Primary references: [Laflamme, Miquel, Paz, Zurek (1996)](https://doi.org/10.1103/PhysRevLett.77.198), [Knill and Laflamme (1997)](https://doi.org/10.1103/PhysRevA.55.900). The SU(2) matter/rishon physical model assumptions and exchange projection are stated in [FOUR_LINK_SU2_DFS.md](FOUR_LINK_SU2_DFS.md) and [SU2_EXCHANGE_HAMILTONIAN_AND_NO_GO.md](SU2_EXCHANGE_HAMILTONIAN_AND_NO_GO.md).
