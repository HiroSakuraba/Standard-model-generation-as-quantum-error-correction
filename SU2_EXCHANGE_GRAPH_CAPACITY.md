# Exact SU(2) exchange-graph error capacity (four and eight spins)

**Status, 18 October 2026.** Exact finite-model tests, not an experimentally calibrated lifetime or a Standard Model anomaly-derived code. Companion verifiers: `su2_bridge_error_geometry.py` and `su2_graph_qutrit_witness.py` (require `su2_error_algebra_capacity.py`). SymPy 1.14; no floating-point acceptance criteria.

## Physical setup and two distinct error demands

Eight co-located spin-half SU(2) fundamentals at **one common gauge vertex**; registers A=(0,1,2,3) and D=(4,5,6,7). Their independent four-spin singlets form a four-dimensional product-singlet code space \(P\) inside the *fourteen-dimensional* global singlet sector \(H_0\). Every pair Heisenberg exchange \(h_{ab}=\mathbf S_a\cdot\mathbf S_b=(\mathrm{SWAP}_{ab}-I/2)/2\) is SU(2)-invariant at this vertex. The graph of permitted couplers is a **stipulated noise restriction**, not a result of SU(2) symmetry.

For a chosen graph G of cross-register edges compare (i) the **single-application** linear error space \(\mathcal E_G=\operatorname{span}\{I,h_{ab}:(a,b)\in G\}\), for which Knill–Laflamme tests products \(E^\dagger F\); and (ii) the **full generated algebra** \(\mathcal A_G=*\text{-alg}(\mathcal E_G)\), which allows arbitrary products/sequences. Passing (i) does not establish (ii). Protection against (ii) requires a suitable multiplicity/noiseless subsystem or another exact recovery code.

## Exact first-order Knill–Laflamme tests inside the 4D product-singlet host

Set \(S^A_{ac}=P_A(\mathbf S_a\cdot\mathbf S_c)P_A\) and similarly \(S^D_{bd}\), using 0-based indices within each four-spin register. SU(2) rotational invariance gives

\[
P h_{ab}P=0,\quad
P h_{ab}h_{cd}P=\frac13 S^A_{ac}\otimes S^D_{bd},\quad
S_{aa}=\frac34I_2.
\]

The script confirms this with five independent exact contractions in the full 256-dimensional eight-spin computational basis.

| Allowed cross-register exchange locations | Maximum correctable code dimension **inside P** | Method |
|---|---:|---|
| Single edge (0,4) | **4** | \(PhP=0\), \(Ph^2P=3P/16\) |
| Two disjoint edges (0,4),(1,5) | **2** | Encoded \(\{|01\rangle,|10\rangle\}\), plus Hermitian interlacing optimality |
| Four-edge star (0,4),(0,5),(0,6),(0,7) | **2** | Encode on A, freeze D to one singlet; exact KL and interlacing |
| Three-edge L (0,4),(0,5),(1,4) | **1** | Four rank-one spectral projectors lie in the KL operator system |
| Four-edge square containing L | **1** | Inherits the L obstruction |
| All 16 bridges (complete K4,4) | **1** | Inherits L; KL operator-system rank is 9 |

A dimension-one code stores no nontrivial quantum information. The upper bounds in this table apply **only to P**, not to every encoding in \(H_0\).

For the disjoint pair,

\[
P h_{04}h_{15}P=\operatorname{diag}(3/16,-1/16,-1/16,1/48).
\]

The middle two eigenstates give an exact 2D correctable code for \(\{I,h_{04},h_{15}\}\). The second and third ordered eigenvalues, \(-1/16\) and \(1/48\), differ, forbidding a rank-three scalar compression by Cauchy interlacing. Its error Gram matrix is

\[
\begin{pmatrix}1&0&0\\0&3/16&-1/16\\0&-1/16&3/16\end{pmatrix}.
\]

**Star, stronger conclusion:** All star exchange operators, and arbitrary products of them, act only on physical spin A0 and the four D spins. For the complete four-spin singlet encoding on A, \(P_A\sigma_{A0}P_A=0\) for Pauli \(\sigma=X,Y,Z\), and \(P_A I_{A0}P_A=I_A\). Hold D in one known singlet. Then the resulting 2D code corrects **every** operator supported on the five-spin subsystem A0+D, not just single faults: its compression to the code is scalar. This is the standard fixed-support erasure/noiseless multiplicity principle, not a demonstration of passive dephasing immunity or a constructed physical recovery gadget. A rank-three code is excluded by the nondegenerate spin-pair observable on D.

**L graph no-go within P:** The four matrices \(I\), \(S^A_{01}\otimes I\), \(I\otimes S^D_{01}\), \(S^A_{01}\otimes S^D_{01}\) all belong to the complex linear span of \(P E^\dagger F P\). Both 2×2 spin-pair operators have eigenvalues \(-3/4,1/4\). Their span contains all four orthogonal rank-one product spectral projectors \(R_j\) with \(\sum_jR_j=P\). If K encodes k≥2 dimensions in P and satisfies KL, \(K^\dagger R_j K=c_j I_k\). Since each LHS has rank at most one, all \(c_j=0\); summing forces \(K^\dagger K=0\), impossible. This applies to arbitrary complex 2D subspaces of P. One extra nonstar bridge can turn the star into this obstructed pattern.

## Full fourteen-dimensional singlet space: graph-algebra multiplicities

Exchange swaps on a connected component of m spins generate \(S_m\), the symmetric permutation group. Spin/permutation representation theory produces matrix-algebra blocks acting on the component multiplicities. Spins not connected by allowed errors supply *protected multiplicity factors*.

| Allowed connected exchange support | Full generated algebra on the 14D singlet | Commutant dimension | Maximum protected code dimension for complete algebra |
|---|---|---:|---:|
| L graph: four interacting spins, four untouched | \(M_2\otimes I_2\;\oplus\;M_3\otimes I_3\;\oplus\;\mathbb C\) | **14** | **3 (qutrit)** |
| Star graph: five interacting spins, three untouched | \(M_5\otimes I_2\;\oplus\;M_4\) | **5** | **2 (qubit)** |
| Connected exchange graph on all eight spins | \(M_{14}\) | **1** | **1 (none)** |

For the L graph, the four interacting spins have spin \(j=0,1,2\) representation multiplicities \(2,3,1\), matched by the four spectator spins. The global-spin-zero sectors therefore have dimensions \(4,9,1\). The 9D sector factors into a **three-dimensional error factor and a three-dimensional protected factor**. For the star, the 5-spin component couples to three spectator spins in spin-1/2 and spin-3/2 sectors of dimensions 10 and 4; the protected factor is two-dimensional.

The code verifies these dimensions with exact 14×14 rational permutation representations, Casimir eigenvalue ranks and commutant nullspaces. The maximum correctable dimensions follow from the matrix blocks' multiplicities: applying the KL rank bound to rank-one projectors inside each full-matrix block gives \(k\le \max m_\alpha\). Connected all-eight exchange gives \(M_{14}\), excluding any nontrivial code correcting **the full algebra**. That is not a no-go for sparse single-fault error spaces.

## Explicit three-level quantum code in the full singlet sector

The companion verifier `su2_graph_qutrit_witness.py` supplies three **explicit 256-component, exact wavefunctions** for the L-graph algebra. Let C=(0,1,4,5) be the four-spin connected component and R=(2,3,6,7) its four unaffected spectator spins. In each group, form three independent spin-1 multiplicity copies from one singlet pair and one triplet pair; orthogonalize exactly. Fix connected-copy index \(c=0\), vary spectator-copy index \(r=0,1,2\), and couple the two spin-1 irreducible representations into total spin zero:

\[
|\overline r\rangle=\frac1{\sqrt3}\big(
|1,+1;c\rangle_C|1,-1;r\rangle_R-
|1,0;c\rangle_C|1,0;r\rangle_R+
|1,-1;c\rangle_C|1,+1;r\rangle_R
\big).
\]

The script verifies all three states are orthonormal, have \(G^z=0\), and are annihilated by total \(G^+\). It then tests **all 24 permutations of C**; each has a scalar compression to the qutrit basis. Since these permutations form a group, the full complex algebra they generate satisfies Knill–Laflamme for arbitrary error Kraus operators in that algebra. This is a constructive error-correction certificate, although no physical preparation or recovery operation is implemented.

This qutrit lies in the larger 14D singlet sector, **not inside the previously excluded 4D product-singlet host P**. There is no contradiction: the available noise multiplicities change when the encoding space changes.

## Interpretation, limits, next test

**A gauge-invariant Hilbert space is not enough. The error algebra's representation and interaction graph decide which multiplicity sectors can protect information.** Restricted graphs can leave a correctable subsystem; unrestricted exchanges on a connected set of all eight spins generate every operator on the global singlet. This does **not** mean a real gauge theory naturally enforces those restrictions.

The eight spins here share *one* SU(2) transformation; a spatial gauge theory with multiple vertices/links may forbid naked cross-vertex exchange or require gauge dressing. Noise on spectators, extra exchange couplers, correlated errors, measurement faults, and their rates are not modeled. No lifetime, circuit-level threshold, self-correcting memory, Standard Model particle implementation or new anomaly-code relation is claimed. CrocSwap's arithmetic frame gauges are not physical SU(2) gauge transformations.

**Next:** constrain exchange graphs using the couplers of a proposed matter/rishon quantum-link device. Add small symmetry-preserving **off-graph** couplings, calculate how their amplitudes destroy the KL conditions and protected multiplicities, and compare exact full-spin unitary leakage and logical infidelity against suitable recovery maps. Do not infer quantitative noise suppression from algebra alone.

## Reproduce

```bash
python -m pip install sympy==1.14.0
python su2_error_algebra_capacity.py
python su2_bridge_error_geometry.py
python su2_graph_qutrit_witness.py
```

All assertions passed locally on 18 October 2026 with SymPy 1.14.0; no sampling or floating-point tolerances.

**Primary background:** Knill and Laflamme, *Theory of quantum error-correcting codes*, Phys. Rev. A **55**, 900 (1997), https://doi.org/10.1103/PhysRevA.55.900; Knill, Laflamme and Viola, *Theory of quantum error correction for general noise*, Phys. Rev. Lett. **84**, 2525 (2000), https://doi.org/10.1103/PhysRevLett.84.2525; Bacon, Kempe, Lidar and Whaley, *Universal fault-tolerant quantum computation on decoherence-free subspaces*, Phys. Rev. Lett. **85**, 1758 (2000), https://doi.org/10.1103/PhysRevLett.85.1758. Also see the existing repository's `SU2_ALLOWED_ERROR_ALGEBRA.md`, `SU2_EXCHANGE_HAMILTONIAN_AND_NO_GO.md`, and `TWO_BODY_SU2_MEDIATED_GATE.md`.
