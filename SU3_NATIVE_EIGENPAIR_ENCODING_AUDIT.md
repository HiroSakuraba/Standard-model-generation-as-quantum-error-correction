# SU(3) five-plaquette storage without inverse evolution: eigenstate-pair trade-off

**18 October 2026 — reproducible finite-model numerical audit, not a new quantum memory.**

Primary source: V. Chen, B. Müller, X. Yao, *Minimally truncated SU(3) lattice gauge theory and string tension*, *Physical Review D* **114**, 074502 (5 October 2026), DOI [10.1103/1qt7-3ywl](https://doi.org/10.1103/1qt7-3ywl). The Hamiltonian and all nine magnetic transition blocks have already been matched exactly against the source in this repository. Here we study possible passive encodings inside the same **five-plaquette, 243-state, gauge-singlet qutrit chain**, without an assumed inverse of the entire Hamiltonian.

Executable: [`su3_native_eigen_encoding.py`](su3_native_eigen_encoding.py). Full deterministic 153-pair scan, including a ten-point Pareto frontier and all modeled channels: [`su3_native_eigen_encoding_results.json`](su3_native_eigen_encoding_results.json). Console output: [`su3_native_eigen_encoding_output.txt`](su3_native_eigen_encoding_output.txt).

## The scientific question

The previous code stores one qubit in the two oppositely oriented flux states at the right end of the chain. An explicit decoder corrects a restricted linear error span exactly. But the **native Hamiltonian itself** evolves that code: without artificial full-Hamiltonian inverse, its recovered entanglement fidelity fell to about **0.63275** at dimensionless time `t=0.6`, even with no injected fault. This motivates states that are compatible with the actual free dynamics.

A pair of eigenstates of a fixed Hamiltonian automatically stays in its two-dimensional span. Its remaining native effect is a *relative logical phase*. If the energies coincide, the pair is exactly stationary up to global phase. If they differ, the relative phase accumulates. We therefore examine whether the eigenstates with small energy splitting are also indistinguishable to the *local gauge-invariant errors* that a simulator could experience.

## Model and declared error operators

The published chain has five dynamical face qutrits, open zero-exterior-flux boundaries, and 16 links with representations restricted to `1`, `3`, `anti-3`. Define the source Hamiltonian, dropping only an irrelevant scalar energy shift:

\[
H(g^2)=g^2 D-\frac{1}{g^2}\sum_{i=0}^{4}B_i,
\quad B_i=W_i+W_i^\dagger,
\quad D=\sum_{e=0}^{15}E_e^2.
\]

At nominal `g²=1`, diagonalize the full real symmetric `243×243` Hamiltonian. Let `|v_i>` be its ordered eigenvectors, `E_i` their energies, and `K_ij=[v_i,v_j]` the logical isometry for an unordered pair `0≤i<j<18`. This is a **restricted 153-candidate survey of energy-eigenstate pairs**, not an optimization over arbitrary rank-two projectors or all eigenstates. The low-temperature feasibility of using excited eigenstates is not established.

The **declared finite, normalized single-event error span** has 32 operators:

\[
\mathcal E=\operatorname{span}\{I,\ n_e=\tfrac34 E_e^2\ (e=0,\ldots,15),\ B_i/2,\ W_i,\ W_i^\dagger\ (i=0,\ldots,4)\}.
\]

Here `n_e` is the link-occupation projector (`0` or `1`). `B_i/2` and all directed `W_i` have operator norms at most one. Operator scaling is explicitly declared because any approximate Knill–Laflamme (KL) metric depends on the chosen noise normalization; it is **not** an experimentally fitted probability or a unique physical distance. Every candidate is tested against **all 32×32 = 1,024 pairs** in the *linear error span*, rather than just the 32 individual compressions. For each overlap, define the traceless `2×2` defect:

\[
A_{ab}=K^\dagger E_a^\dagger E_b K,\qquad
\delta_{ab}(K)=\|A_{ab}-\tfrac12\mathrm{tr}(A_{ab})I_2\|_F.
\]

The reported score is `max_ab δ_ab` (Frobenius norm). **Exact** KL against the full declared span requires this score to be zero; a small **nonzero** number is neither a recovery certificate nor a rigorous upper bound on actual infidelity. The source-matrix lookup and eigenpair scores are reproducible numerically with double precision, not exact algebraic claims.

## Tested results at `g²=1`

| Pair (zero-indexed energy rank) | Energy splitting `ΔE` | Worst normalized KL defect | Largest defect comes from |
|---|---:|---:|---|
| Lowest-energy `0,1` | 5.091694 | 0.547301 | `W2`, `W2†` |
| Minimum splitting of first 18: `15,16` | 0.00001176 | 0.665892 | `n3`, `n4` |
| Minimum KL defect of first 18: `12,15` | 0.555643 | **0.139363** | `I`, `W4†` |

There are **ten nondominated pairs** in the two-dimensional (splitting, KL-defect) trade-off among these 153 candidates; the full list is in the JSON. **No pair among the tested 153 has zero KL defect**. The nearly degenerate pair is a poor error-correcting code under the tested physical single-event span. The candidate with lower KL defect pays for it with larger splitting and higher excitation energy: its two energies are about `8.5709` and `9.1265` above the ground state in these source units.

## Channel-level tests: actual trace-preserving quantum maps

Unlike the code's earlier *globally conditioned* demonstration errors, these models act unconditionally on **all 243 basis states**. Probabilities, time and couplings are deliberately chosen stress-test parameters, **not measured decoherence rates**.

### Native Hamiltonian without inverse

At `t=0.6`, native undisturbed evolution of an eigenpair only rotates its logical phase. Under the explicit **trace-preserving projection-and-replacement decoder**

\[
\mathcal R_K(\rho)=K^\dagger\rho K+\mathrm{tr}[(I-KK^\dagger)\rho]I_2/2,
\]

no inverse of the full Hamiltonian is required. Without phase compensation, eigenpair entanglement fidelity is exactly the analytic finite-matrix value `cos²(ΔE t/2)`. If the nominal energy difference is *known perfectly*, a subsequent **virtual logical Z-frame correction** `diag(exp(+it E_i),exp(+it E_j))` can ideally restore the relative phase. This frame must eventually be accounted for in physical readout/gates and is **not** a demonstrated free recovery control. For a perfectly known, time-independent Hamiltonian, it is not an error-correction achievement to undo its deterministic logical phase.

| Encoding | `t=0.6` no inverse, no frame correction | Ideal virtual logical phase frame |
|---|---:|---:|
| Lowest energy pair `0,1` | 0.001873 | 1.000000 |
| Nearly degenerate pair `15,16` | **1.000000** | 1.000000 |
| Lowest KL defect pair `12,15` | **0.972470** | 1.000000 |

The original flux code, with its **different and stronger five-syndrome CPTP decoder**, reaches `F_e = 0.632749` under the same native evolution without inversion or frame tracking. These are **different decoders and reference frames**; the comparison is designed to expose the native-Hamiltonian incompatibility of the old code, **not** to establish a fair hardware performance ranking.

### Quasi-static uncertainty in the native couplings

Two equally likely Hamiltonians `H(g²(1+ε))` and `H(g²(1−ε))` act for the full `t=0.6`. Decode with `R_K` and the **fixed nominal** logical frame, with no fault-dependent phase adjustment. At nominal `g²=1`:

| Encoding | `ε=0.02` | `ε=0.05` | `ε=0.10` |
|---|---:|---:|---:|
| `0,1` | 0.998373 | 0.989886 | 0.960313 |
| `15,16` | 0.999283 | 0.995522 | 0.982167 |
| `12,15` | 0.999144 | 0.994660 | 0.978757 |

An eigenpair of the **nominal** Hamiltonian need not remain exactly invariant when the coupling changes, so the projection-and-replacement complement is essential to preserve trace even under leakage. These are calibration-uncertainty channels, not environmental field measurements.

### Directed Wilson jumps and electric faults

For each selected plaquette, use the same globally valid two-Kraus channel from the preceding experiment:

\[
K_1=\sqrt p\,W_i,\qquad K_0=\sqrt{I-pW_i^\dagger W_i}.
\]

For a local electric occupation `n_e`, use the analogous pair `(sqrt(p)n_e, sqrt(I-pn_e))`. Each is **CPTP on all 243 states**, not a postselected jump. Apply two fault opportunities separated by **three equal segments** of the native Hamiltonian at total `t=0.6`, and decode with `R_K` and the fixed nominal virtual logical frame. The shot probability parameter `p=0.10` is defined **per fault opportunity**. No full-Hamiltonian inversion is used.

| Encoding | Two `W3` jumps, `p=.10` | Two `W4` jumps, `p=.10` | Two link-4 electric faults, `p=.10`, `t=0` |
|---|---:|---:|---:|
| Lowest `0,1` | 0.897348 | 0.886667 | 0.979978 |
| Nearly degenerate `15,16` | 0.916734 | 0.913949 | 0.952485 |
| Lowest KL defect `12,15` | **0.919992** | **0.917517** | 0.962541 |

The old flux code with the actual five-syndrome decoder **exactly corrects** the same two link-4 electric faults at `t=0`, yielding `F_e=1` (to rounding), as guaranteed by its exact KL proof. Therefore the alternative eigenpairs are **not uniformly better**: passive compatibility trades against correction of real local errors.

## What has been accomplished and what remains open

**Accomplished:** a reproducible native-compatible energy-eigenpair search, including every local electric projector and both directions of all five magnetic plaquette Wilson operators; an explicit, globally CPTP code-projection decoder; finite-time channels with no inverse Hamiltonian; and a finite-list Pareto trade-off between splittings and errors. The script crosschecks the original five-face syndrome-decoder benchmark (approximately `0.632749` at `t=.6`) independently.

**Not established:** any exact active-code correction for the eigenpairs, an optimal decoder, a physically executable energy-eigenstate encoder/readout or logical phase tracking, immunity to arbitrary local errors, a stable thermally accessible ground-space memory, full Standard Model implementation, glueball spectroscopy, or a continuum mass-gap proof. Importantly, the number **0.139363** is *not* an infidelity and should never be reported as an error rate.

**Recommended decision:** Do not promote any scanned eigenpair to a protected memory. The near-degenerate candidate passes the Hamiltonian compatibility test but fails the chosen local-error scalarity test. The most promising low-KL pair has an appreciable logical energy splitting and is highly excited. The next hard engineering question is whether a **bounded local control or measurement scheme** can simultaneously preserve a two-dimensional subspace and suppress *both* directed Wilson error orientations under realistic calibration uncertainty. Consider small joint optimization of pulse/control and code only after declaring what operations the quantum-link hardware actually permits.

## Reproduction

```bash
python -m pip install sympy==1.14.0 numpy scipy
OPENBLAS_NUM_THREADS=1 python su3_native_eigen_encoding.py --output-json su3_native_eigen_encoding_results.json
python -m py_compile su3_native_eigen_encoding.py
```

The new script uses existing repository modules `su3_three_plaquette_chain.py`, `su3_four_five_plaquette.py`, and `su3_channel_recovery.py` (which itself imports the validated five-face model). Numeric spectrum and channel results are double precision, not Arb/SageMath interval enclosures. Primary cited equations use the fully published APS version from **5 October 2026**, not an older preprint table numbering.