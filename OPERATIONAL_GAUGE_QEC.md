# From anomaly-derived matrices to operational quantum recovery: an exact gauge-theory benchmark

**Status:** exact computational checks plus standard stabilizer-code arguments for a toy, finite **Z_p lattice gauge theory**; not a claim that real Standard Model particles implement any of these checks. All runs are deterministic; nothing is peer-reviewed here. Source scripts: \`operational_gauge_qec.py\` (no dependencies) and \`abstract_stabilizer_recovery.py\` (NumPy, the existing charge-lattice toolkit).

## The actual question

Anomaly cancellation is a **consistency condition on quantum field-theory charge assignments**. Cubic triorthogonality turns some of those assignments into generator matrices of valid *abstract* quantum stabilizer codes. Neither fact supplies the **physical ingredients** of error correction:

1. An encoded subspace inside a specified physical Hilbert space, with a clear meaning for logical superpositions.
2. A declared physically motivated noise algebra or list of Kraus error operators.
3. Measurable observables whose syndromes distinguish those errors without revealing the logical state.
4. A recovery channel (a syndrome-dependent corrective operation) that satisfies the **Knill–Laflamme conditions**.
5. A dynamics and symmetry analysis checking that these operations make sense for the theory, not simply for qudit hardware.

The last item is the principal missing step for the full Standard Model. Physical gauge charges impose local constraints and superselection restrictions; treating charge labels as freely manipulated qudit degrees of freedom without addressing those restrictions is not an operational model.

## A general obstruction: charge-only syndromes

**Lemma.** Let \(P\) project onto a proposed code space stabilized by a commuting collection of gauge/charge measurements \(\{G_a\}\). Suppose a physical unitary error \(E\) commutes with every \(G_a\) and preserves the code space, but its restriction \(PEP\) is **not** a scalar multiple of \(P\). Then the error set \(\{I,E\}\) **cannot** be corrected on that code space, since the Knill–Laflamme condition \(PI^\dagger EP=cP\) fails. Measuring the \(G_a\) also produces the same syndrome before and after \(E\).

The assumptions matter. This does **not** rule out adding independent checks that distinguish an \(E\) which preserves the gauge-physical subspace but moves states outside the *smaller* quantum code, as the magnetic-check toric example below demonstrates. It also does not say that every gauge-invariant error has nontrivial logical action.

## Result A. Gauge invariance by itself does not protect a quantum memory

Take an elementary closed ring of three p-level links, with electric-flux computational basis \`|e_0,e_1,e_2>\`. No matter or background charge is included. Impose local Gauss laws

\`\`\`
G_i = Z_i (Z_(i-1))^-1 = +1, for i in Z_3,
Z|j> = omega^j|j>, omega=exp(2*pi*i/p).
\`\`\`

The **gauge-invariant physical Hilbert space** is exactly

\`\`\`
C_Gauss = span{ |j,j,j> : j in F_p }.
\`\`\`

It therefore has p dimensions and could *formally* encode one logical p-level system.

- An isolated link shift \`X_i^a\`, \`a != 0\`, changes the electric flux at its two endpoints. Each of the \`3(p−1)\` single-link shift errors has its own nonzero Gauss-law syndrome. Measuring \`G_i\` and applying \`X_i^(-a)\` repairs this **specific gauge-violating error set**.
- But a local electric phase \`Z_i^b\`, \`b != 0\`, **commutes with all the Gauss generators** and acts on the p-dimensional physical subspace as
  \`Z_i^b |j,j,j> = omega^(b*j)|j,j,j>\`.
  This is a **nontrivial logical Z**, not a scalar. Thus \`P Z_i^b P\` is not proportional to \`P\`; the Knill–Laflamme condition for the error pair \`{I,Z_i^b}\` fails.

**Conclusion:** Gauss-law conservation detects certain failures to stay in the physical sector. It does **not**, by itself, detect or correct all local, gauge-invariant errors *within* that sector. In a fundamental theory enforcing gauge invariance exactly, the gauge-violating \`X_i^a\` errors may not be allowed physical operations; they are relevant instead to an imperfect **simulator** whose gates can leak out of its target gauge sector.

This is a concrete obstruction, not a universal no-go against using gauge theories to encode quantum information.

## Result B. Adding a genuinely independent gauge-invariant check makes recovery possible

Use the standard **Z_p toric/quantum-double code** on a periodic \`3×3\` square lattice with **18 directed link qudits**.

Let \`Z_star(v)\` be the product of link-Z operators with \`+1\` outgoing and \`−1\` incoming exponents (the **Gauss-law/electric** constraint). Let \`X_face(f)\` be the oriented product of four link-X operators around face f (a gauge-invariant **magnetic** shift). These commute since the divergence of a closed loop is zero.

There are nine vertices and nine faces, but one dependence of each type:

\`\`\`
rank(Gauss Z-star) = 8
rank(magnetic X-face) = 8
n = 18,  k = n - 8 - 8 = 2.
\`\`\`

**An explicit physical basis:** let \`C_face\` be the eight-dimensional space of oriented face boundaries in \`F_p^18\`, and let \`gamma_h, gamma_v\` be unit-flux configurations around the two independent noncontractible cycles. Then two logical qudits are encoded by the normalized states

\`\`\`
|a,b>_L = p^-4 sum_(c in C_face) |a*gamma_h + b*gamma_v + c>,
a,b in F_p.
\`\`\`

All basis states in these superpositions have zero divergence and hence obey the Gauss constraints. The uniform sum is invariant under adding a face boundary and therefore obeys all magnetic checks. The \`p^2\` distinct homology cosets are orthogonal. Noncontractible X and Z loops serve as logical Weyl operators; their symplectic intersection pairing is verified in the script.

**Exact recovery test:** For each edge i and each nonidentity single-edge Weyl error \`E_i(a,b)=X_i^a Z_i^b\`, the script computes the full **18-component** syndrome (nine Gauss and nine face-check values; each nine-component group has one linear dependence and rank eight). For p=3 and p=5:

| p | Encoded code | Number of identity + single-edge Pauli syndromes | Collisions |
|---:|---|---:|---:|
| 3 | \`[[18,2,3]]_3\` | 145 | 0 |
| 5 | \`[[18,2,3]]_5\` | 433 | 0 |

Identity has zero syndrome; every other listed error has its own *nonzero* syndrome. Consequently each \`E_i^\dagger E_j\` for distinct errors has a nontrivial syndrome, so its code-space projector vanishes. For identical errors \`E_i^\dagger E_i=I\`. Thus the necessary-and-sufficient **Knill–Laflamme criterion** holds for all one-site Pauli operators. Since they span the space of arbitrary one-site error operators, the code corrects **arbitrary noise affecting one unknown link**. Measure the complete syndrome and apply the inverse of the uniquely identified one-site Pauli error.

This is a **standard, previously known quantum code**, not an invention or a prediction of the Standard Model. The extra plaquette checks are essential: a single \`Z_i\` commutes with **every Gauss law** but triggers the magnetic face measurements.

The code distance is exactly 3: the syndrome uniqueness excludes logical errors of weight 1 or 2; the explicitly constructed noncontractible length-three logical loop witnesses show distance at most 3.

## Result C. Our charge-matrix codes genuinely pass the *abstract* recovery condition

The separate \`abstract_stabilizer_recovery.py\` uses **the actual repository's charge matrices**. It constructs all the CSS X-check generators and the Z-check nullspace, verifies their commutation, and enumerates every identity and single-site error \`X_i^a Z_i^b\`.

For p=3 and p=5, it checks:

| Abstract code | Independent X/Z syndrome classes for identity and one-site Weyl errors |
|---|---:|
| \`[[15,1,3]]_3\` | 121 |
| \`[[30,2,3]]_3\` | 241 |
| \`[[45,3,3]]_3\` | 361 |
| \`[[15,1,3]]_5\` | 361 |
| \`[[30,2,3]]_5\` | 721 |
| \`[[45,3,3]]_5\` | 1081 |

Every syndrome is nonzero except identity and **distinct**. The resulting lookup gives an exact inverse-Pauli recovery rule for arbitrary single-qudit errors in the *formal p-level register* model. This is more operational than checking only Z-distance: it explicitly verifies both X- and Z-error syndromes and arbitrary mixed single-site Weyl errors.

However, this **still does not** place these formal qudits in the Hilbert space of particles or show that the Standard Model Hamiltonian preserves their code space. The charge checks \`S_X\` in the mathematical stabilizer construction are **not automatically conserved gauge observables or locally measurable particle operations**. That distinction is non-negotiable.

## How this updates the broader program

We now have a clean division of three claims:

**Established mathematical claim:** odd-prime, two- and three-family charge matrices define formal \`[[30,2,3]]_p\` and \`[[45,3,3]]_p\` qudit CSS codes (at the primes used here), with explicit one-site Weyl syndrome and Knill–Laflamme recovery in the register model.

**Established operational example in a different theory:** a **Z_p lattice-gauge simulator**, with independently imposed magnetic constraints, has an explicit physical encoded subspace, nontrivial logical degrees of freedom, a declared error model, measurable check algebra, and a complete recovery rule. Gauge constraints alone fail for gauge-invariant phase errors.

**Not established for the full Standard Model:** a charge-lattice encoding inside the physical \`SU(3)_c × SU(2)_L × U(1)_Y\` Hilbert space, exact check operations consistent with CKM/PMNS mixing and Yukawa couplings, physically allowed noise operators, and a recovery map.

The next *research* task should choose a manageable intermediate theory, e.g. a truncated **SU(2)** lattice gauge model with matter, and test its **physical gauge-invariant error algebra** and implementable recovery measurements. Recent independent work addresses the right category of problem: Rajput, Roggero and Wiebe ([npj Quantum Information 9, 41 (2023)](https://doi.org/10.1038/s41534-023-00706-8)) use Gauss-law redundancy together with phase-flip correction; Yao ([Physical Review D 113, 114512 (2026)](https://doi.org/10.1103/5yg7-6mph)) develops codes for truncated SU(2) lattice gauge theories. Neither paper establishes that our SO(10)-charge code is physically realized.

## Reproduction

\`\`\`bash
python3 operational_gauge_qec.py
python3 -m pip install numpy
python3 abstract_stabilizer_recovery.py
\`\`\`

Both scripts use exact modular arithmetic and symbolic stabilizer logic; neither attempts an exponentially large \`p^18\` statevector allocation, Monte Carlo trials, or approximate interval verification. The GitHub Actions research workflow runs them with the existing deterministic checks.

**References:** Emanuel Knill and Raymond Laflamme, [*Theory of quantum error-correcting codes*, Physical Review A **55**, 900 (1997)](https://doi.org/10.1103/PhysRevA.55.900); Rajput, Roggero and Wiebe (2023), and Yao (2026), linked above.
