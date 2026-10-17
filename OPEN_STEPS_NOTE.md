# The five open steps: flavour physics, logical channels, equivalence, operational setting, and new lattices

**Research extension, 17 October 2026.** Independent exploratory calculation, **not peer reviewed**. This note continues `RESEARCH_NOTE.md` and `THREE_GENERATION_FLAVOR_NOTE.md`, and works through the five "next steps" listed in the README. Everything below is reproduced by `python3 open_steps.py`. It needs Python 3.9+ and numpy, takes about 4 minutes (`--quick` skips the p = 13 search), and asserts every number quoted here. `open_steps_results.txt` records a complete run. The shared toolkit is `charge_lattice_tools.py`.

Source paper: Shiroman Prakash and K. L. Ramya, *One generation of the Standard Model as a quantum error-correcting code*, [arXiv:2610.21879](https://arxiv.org/abs/2610.21879).

## Executive summary

| Step | Question | Answer |
|---|---|---|
| 1 | Do realistic flavour Lagrangians leave any of these charges intact? | **No exact survivor.** The B3/B5 badge rows already fail to commute with SU(3)c × SU(2)L. The only gauge-compatible family charges are per-family **(B−L)** and **T₃R**. They give new `[[30,2,3]]_p` and `[[45,3,3]]_p` codes for every p ≥ 5, but an exhaustive Yukawa scan finds no realistic texture that preserves any family-differential charge (1HDM or 2HDM). Distance 3 needs charges on lepton doublets, which neutrino mixing breaks at O(1). The best approximate structure (right-handed T₃R per family, broken at 0.4–6 %) only gives distance 2. |
| 2 | What do the logical error channels look like? | Exact enumerators for all codes. Single-output error is the same as for direct sums at leading order. However, the coupled codes have **correlated two-output logical errors at order ε⁴**, while direct sums only get them at **ε⁶**. These come entirely from the 2+2 cross-family events. |
| 3 | Are the coupled codes equivalent to direct sums? | **No, not even under arbitrary local unitaries with qudit permutations.** Every coupled stabilizer group is indecomposable: it has one matroid component, while a direct sum always has two. This covers monomial and local-Clifford equivalence. |
| 4 | Where could the charge constraints be measured? | In a conjugate basis the X checks are **Z_p charge measurements** and the logical value is "the charge of the missing right-handed neutrinos." That is natural in engineered qudit registers and gauge-theory simulators, but impossible with real gauge charges because of superselection. Measuring the coupled codes costs about 18 % more total check weight than the direct sum. The correlated ε⁴ error is a concrete experimental signature. |
| 5 | Distance 4, or a transversal entangling gate? | **No distance 4** in any lattice searched. A proof rules it out for every family-wise Cartan-linear lattice, and exhaustive searches cover integer-anomaly-free nonlinear rows. Mod-p anomaly-free nonlinear rows exist for p = 7 and 13 and cut the weight-3 weak spots from 50 to **8**, lowering the leading distillation error 6.25× at p = 7. **No transversal diagonal gate** (cubic or quadratic, any per-qudit powers) acts as an entangling logical gate in any code tested. |

---

## Step 2 — Exact logical channels and distillation polynomials

**Model.** This is the source paper's symmetric model. Each input qudit independently suffers `Z^j` (j = 1…p−1), each with probability ε/(p−1). Distillation post-selects on a trivial X-stabilizer syndrome. An accepted error `e` has `G_P e = G_R c` and applies the logical Pauli `Z^c` (c ∈ F_p^k).

**Method (exact).** For every `s`,

```
Σ_{e : G_P e = s} x^{wt e} = p^{-r} Σ_a ω^{-a·s} (1+(p−1)x)^{n−w(a)} (1−x)^{w(a)},   w(a) = wt(a G_P).
```

Put `s = G_R c`, so that `a·s = b·c` with `b = a G_R`. Grouping all `a ∈ F_p^r` by `(w(a), b)` gives integer histograms. The character sum collapses because `a → u·a` (u ≠ 0) permutes the classes `b·c = t ≠ 0`. This yields the exact weight distribution `A_w(c)` of every logical class, for every code. Direct sums are obtained by exact convolution, which was checked against direct computation. The toolkit reproduces the known qubit 15-to-1 value `35ε³` and the paper's qutrit `6.25ε³`. All classes agree with an independent brute-force enumeration of **every** error of weight ≤ 4.

**Leading terms** (the probability that an accepted output block carries a logical error, summed over outputs):

| Code | single-output logical error | correlated (≥ 2 outputs) | per-output infidelity at ε = 10⁻³ | P(correlated \| accept) at ε = 10⁻³ |
|---|---|---|---|---|
| Coupled B3, 2 families, p = 3 | 25/2 ε³ − 1271/4 ε⁴ … | **5 ε⁴** − 83 ε⁵ … | 6.2838×10⁻⁹ | 5.1×10⁻¹² |
| Direct sum, 2 × [[15,1,3]]₃ | 25/2 ε³ − 1295/4 ε⁴ … | 625/16 ε⁶ … | 6.2757×10⁻⁹ | 3.9×10⁻¹⁷ |
| Coupled B5, 2 families, p = 5 | 25/8 ε³ − 649/8 ε⁴ … | **9/2 ε⁴** − 3027/32 ε⁵ … | 1.5734×10⁻⁹ | 4.5×10⁻¹² |
| Direct sum, 2 × [[15,1,3]]₅ | 25/8 ε³ − 1325/16 ε⁴ … | 625/256 ε⁶ … | 1.5680×10⁻⁹ | 2.5×10⁻¹⁸ |
| Coupled B3, 3 families, p = 3 | 75/4 ε³ … | **15 ε⁴** … | 6.2920×10⁻⁹ | 1.5×10⁻¹¹ |
| Coupled B5, 3 families, p = 5 | 75/16 ε³ … | **27/2 ε⁴** … | 1.5788×10⁻⁹ | 1.4×10⁻¹¹ |
| Direct sums, 3 families | as above ×3/2 | ∝ ε⁶ | 6.2757×10⁻⁹ / 1.5680×10⁻⁹ | 1.2×10⁻¹⁶ / 7.4×10⁻¹⁸ |

**Where the correlated term comes from.** In the coupled two-family codes, the weight-4 logical error patterns split as follows:

| | 2+2 across families, correlated | 2+2 across families, single output | 4+0 within one family |
|---|---:|---:|---:|
| p = 3 (B3) | **80** | 96 | 220 |
| p = 5 (B5) | **1152** | 432 | 400 |

The 4+0 patterns are exactly those of the separate one-family codes. The 2+2 patterns are new, and **most of them hit both logical outputs at once**. This is the channel-level content of the "four-tile cross-family events" from the three-generation note.

**Interpretation.**
- Through order ε³ the coupled codes are indistinguishable from direct sums: same acceptance, same single-output error. Fewer X checks never helps the leading term.
- At ε = 10⁻³ the per-output infidelity is only 0.13 % (p = 3) to 0.69 % (3 families, p = 5) worse than separate codes.
- The *correlation* structure is qualitatively different. Correlated output errors are 10⁵–10⁶ times more likely than for separate codes, and appear two orders earlier in ε.
- This matters for protocols that consume both outputs together or feed them into a second distillation round as a block. It is invisible to protocols that consume the outputs independently.
- **Net:** no distillation advantage, and a small, structurally interesting disadvantage.

## Step 3 — Equivalence to direct sums

**Invariant.** For a stabilizer code with group S and a bipartition A|B of the qudits, the code space factorizes as `C_A ⊗ C_B` **iff** `S = S_A × S_B`, where `S_A` is the subgroup supported on A.

*Proof sketch.* Expand the projector `Π = |S|⁻¹ Σ_{s∈S} s` in the Pauli basis. A product projector `Π_A ⊗ Π_B` has Pauli support `T_A × T_B`, and since 1 ∈ T_B every `(t_A, 1)` lies in S.

This factorization property is preserved by any local unitary combined with any qudit permutation. Hence it is preserved by monomial and local-Clifford equivalence. A direct sum of codes always factorizes.

**Computation.** For a CSS group, `S_A = S_{X,A} × S_{Z,A}`. Linear-code splittings are exactly unions of **matroid connected components**, and `S_Z = ker G_P` has the same components as its dual `rowspace(G_P)`. So the finest factorization of S is the join of the component partitions of `S_X` and `rowspace(G_P)`. The components come from fundamental circuits of an RREF basis.

| Code | components of S | dim S − (dim S_A + dim S_B) across the family cut |
|---|---|---:|
| Direct sum (p = 3, 5) | 15 + 15 | 0 |
| Coupled B3, 2 families, p = 3 | **30** | 4 |
| Coupled B5, 2 families, p = 5 and 7 | **30** | 6 |
| Coupled B3 / B5, 3 families | **45** | 4 / 6 |
| Gauge-compatible (B−L, T₃R) lattice, 2 / 3 families, p = 5 | **30 / 45** | 6 |

**Conclusion.** No coupled code is equivalent to a direct sum, or to *any* tensor product of smaller codes, under local unitaries and qudit permutations. That rules out the monomial and local-Clifford cases asked about.

Still open: whether the coupled codes are equivalent to some *other* known indecomposable qudit code, for example from the qudit Reed–Muller family of [Campbell, Anwar & Browne](https://arxiv.org/abs/1205.3104). That needs a literature comparison, not just this invariant.

## Step 1 — A concrete flavour Lagrangian

### 1a. The smart-badge rows already conflict with the gauge group

A family-dependent charge is a symmetry of the QCD and electroweak kinetic terms only if, within each family, it is constant on every SU(3)c × SU(2)L multiplet. Gluons and W bosons are shared by all families, so a family-specific colour or isospin Cartan rotation cannot be undone.

None of the B3 or B5 rows passes this test: every one already splits the quark doublet Q. Solving the constraint mod 3, 5, 7 and 101 gives the complete gauge-compatible space per family:

```
S3 = (1,1,1,0,0)·q = 3(B−L),     S2 = (0,0,0,1,1)·q = −4 T₃R      (equivalently B−L and hypercharge Y = S3/6 − S2/4)
```

| multiplet | Q | u^c | d^c | L | e^c | ν^c |
|---|---|---|---|---|---|---|
| (S3, S2) | (1,0) | (−1,2) | (−1,−2) | (−3,0) | (3,−2) | (3,2) |

### 1b. A gauge-compatible version of the smart badges exists for p ≥ 5

Use rows `f_g·S3` and `f_g·S2` with the family labels from before: (±1) for two families, and (1,0), (−1,0), (0,1) for three.

| Lattice | p = 3 | p = 5, 7, 11 | X distance (p = 5, 7) |
|---|---|---|---|
| (B−L)_g + T₃R_g, 2 families | [[30,2,**2**]] | **[[30,2,3]]** | 7 |
| (B−L)_g + T₃R_g, 3 families | [[45,3,**2**]] | **[[45,3,3]]** (p = 5, 7) | 7 |
| (B−L)_g only | [[30,1,3]] | [[30,2,2]] | 14 |
| T₃R_g only | [[30,2,2]] | [[30,2,2]] | 7 |
| Y_g only ("tri-hypercharge") | — | [[30,1,3]] / [[45,1,3]] | 14 / 21 |

*Why distance 3 for all p ≥ 5.* By the original argument, a weight-2 logical must be a matched twin pair whose differential vector (S3, S2) is a nonzero multiple of ν^c's (3, 2). Check each multiplet:
- Q: (1,0) would need 1 ≡ 0.
- u^c: (−1,2) would need λ = 1 and 4 ≡ 0.
- d^c: λ = −1 and 2 ≡ 0.
- L: λ = 0 and 3 ≡ 0.
- e^c: λ = −1 and 6 ≡ 0.

All are impossible for p ≥ 5. For p = 3 the e^c twins fail, which is the single weight-2 support found.

Two consequences:
- **Tri-hypercharge alone gives k = 1**, because the right-handed neutrino has zero hypercharge. The family-gauged hypercharge models of [Fernández Navarro & King](https://arxiv.org/abs/2305.07690) therefore cannot by themselves separate the logical qudits. Family-dependent B−L (or T₃R) is required.
- Both rows are needed for distance 3: B−L_g alone has lepton-doublet twins, and T₃R_g alone has right-handed twins.

### 1c. Exhaustive Yukawa scan: no exact survivor

Take every family-dependent charge assignment `θ_g = a_g S3 + b_g S2` mod p, together with every Higgs charge. That's one Higgs doublet, and separately two doublets with independent up and down charges. Each Yukawa entry is allowed iff its total charge vanishes mod p. Call a texture **realistic** if:
- Y_u, Y_d and Y_e can all be full rank (a perfect matching in the allowed pattern, then confirmed numerically with generic entries);
- the CKM matrix has no zero entries for generic entries;
- some neutrino mass mechanism (Weinberg operator, Dirac, or type-I seesaw) gives at least two massive neutrinos and a PMNS matrix with no zero entries.

| Families | p | Higgs | realistic assignments | with any family-differential charge |
|---|---|---|---:|---:|
| 2 | 3, 5, 7 | 1HDM | 9, 25, 49 | **0** |
| 3 | 3, 5, 7 | 1HDM | 9, 25, 49 | **0** |
| 2 / 3 | 5 / 3, 5 | 2HDM | 25 / 9, 25 | **0** |

The survivors are exactly the p² family-universal assignments. For those all punctures coincide and **k = 1**.

This is the expected result: the only exact abelian symmetries left by generic Yukawas are family-universal. Here it is checked exhaustively within this class.

### 1d. How badly broken? (approximate symmetries)

Because exact symmetry is impossible, the question becomes how small the breaking can be. A family charge acting on **doublets** with distinct values on two blocks of families forces the mixing matrix to be block-diagonal. Its breaking is therefore at least the largest unavoidable off-block mixing element, minimised over block assignments:

| charge separates | quark doublets (CKM) | lepton doublets (PMNS) |
|---|---:|---:|
| family 3 from (1,2) | 0.041 | 0.66 |
| family 1 from (2,3) | 0.23 | 0.55 |
| all three | 0.23 | 0.66 |

A family charge acting only on **right-handed** fields (T₃R_g) is unaffected by CKM or PMNS mixing. With a single Higgs it allows one Yukawa column per sector, so the breaking is set by mass ratios: m_c/m_t ≈ 0.004, m_s/m_b ≈ 0.019 and m_μ/m_τ ≈ 0.059 at M_Z. This is a standard Froggatt–Nielsen-type approximate symmetry.

(Inputs: CKM global fit from the [PDG review](https://pdg.web.cern.ch/2020/reviews/rpp2020-rev-ckm-matrix.pdf); PMNS from [NuFIT 6.0](https://arxiv.org/abs/2410.05380), normal ordering; running masses at M_Z from [Huang & Zhou](https://arxiv.org/abs/2009.04851). These are illustrative central values, and the conclusions don't depend on their precision.)

**Verdict on an "approximate code".**
- The best-protected structure, T₃R_g with breaking of a few percent, gives only **distance 2**. It can detect a single error but not correct one.
- Every distance-3 option needs B−L_g, which charges the lepton doublets (S3 = −3). Neutrino mixing breaks that at **O(1)**.
- So no distance-3 code is protected by a small symmetry-breaking parameter in the observed flavour sector.
- Majorana neutrino masses would also break B−L itself, even family-universally.

Here "approximate" means an approximate selection rule of the Lagrangian. It does not define a quantum noise channel on the code.

## Step 4 — An operational setting

**Charge-register dictionary.** A codeword is

```
|c⟩ ∝ Σ_{a : a·G_R = c} | a·G_P ⟩ .
```

Each computational basis state records, for every physical species j, the charge `a·h_j` (mod p) of that species under a "probe" charge combination a. The logical label c is the charge carried by the punctured right-handed neutrinos.

- **X stabilizers** shift the probe a by combinations under which every ν^c is neutral.
- **Z stabilizers** `Z(z)`, `z ∈ ker G_P`, correspond to charge-neutral multisets of species.
- **Z logicals** are multisets carrying exactly the charge of c right-handed neutrinos.
- The transversal cubic phase gives `ω^{Σ_j (a·h_j)³} = ω^{−Σ_i c_i³}`. This is *precisely* where anomaly cancellation `Σ_all (a·h)³ = 0` enters.

**Conjugate basis.** After a Fourier transform on every qudit, `X(aG_P)` becomes `exp(2πi a·Q(n)/p)` with `Q(n) = Σ_j n_j h_j`. In that basis:
- the X checks are **measurements of conserved Z_p charges** of an occupation pattern n;
- the code space is the set of charge sectors with `Q(n) ∈ span(G_R)`;
- the logical value is "how many right-handed-neutrino units of charge the other species carry, mod p."

**Where this can be measured, and where it can't.**
- **Real gauge charges:** no. Logical basis states sit in different total-charge sectors, and superselection forbids coherent superpositions of exact gauge charges. So the code cannot live in physical particle charges, even before considering flavour mixing.
- **Engineered qudit registers:** yes. The codes are ordinary CSS codes. Trapped-ion processors have demonstrated qudits up to d = 7, with qutrit Clifford errors of about 2×10⁻³ ([Ringbauer et al.](https://arxiv.org/abs/2109.06903)). The 15-, 30- and 45-qudit codes need larger registers than that demonstration, plus ancillas.
- **Gauge-theory quantum simulators:** these already measure and enforce Gauss-law charge constraints ([Yang et al., 71-site simulator](https://arxiv.org/abs/2003.08945); [Halimeh & Hauke review](https://arxiv.org/abs/2204.13709)). The Z_p-charge form of the X checks fits that toolbox naturally, as global, non-gauged conserved quantities.

**Measurement cost** (a minimum-weight check basis, found by a greedy matroid basis over all stabilizers; p = 5):

| Code | X-check weights | Z checks | total check weight |
|---|---|---|---:|
| [[15,1,3]]₅ | 8, 8, 8, 8 | 10 × weight 4 | 72 |
| Direct sum, 2 families | 8 × weight 8 | 20 × weight 4 | 144 |
| Coupled B5 / gauge-compatible (B−L, T₃R) | 15, 15, 16, 16, 16 | 23 × weight 4 | **170** |

The coupled codes have fewer X checks but much heavier ones. Their total circuit cost is about 18 % higher.

**A concrete experimental signature.** Inject independent random Z errors at rate ε, post-select, and record the *joint* logical outcome of both outputs. The coupled code predicts correlated logical errors at `(9/2)ε⁴` (p = 5), where the direct sum predicts `(625/256)ε⁶`. This distinguishes the two codes without any physics interpretation.

## Step 5 — Distance 4 and entangling gates

### 5a. No-go for every Cartan-linear lattice

**Proposition.** Suppose that within each family g, every charge row is a linear function of the SO(10) weight, so the column is `M_g q`. Suppose also that each punctured column `M_g q0` is nonzero and the punctures are independent. Then `d_Z ≤ 3`.

*Proof.* The source paper's Theorem 1 gives `q0 = αq_a + βq_b + γq_c` over F_p for three non-singlet weights. Then `M_g q0 = αM_g q_a + βM_g q_b + γM_g q_c`, so the weight-3 error (α, β, γ) on those species of family g has syndrome-free image `G_R e_g ≠ 0`. That is a nontrivial logical.

This covers every lattice in the series so far: B3, B5, (B−L, T₃R), direct sums, and any family-dependent Cartan charges.

**A single family is rigid.** The only rows compatible with the five SO(10) Cartan rows, even at the linear level of triorthogonality, are the Cartan rows themselves (checked for p = 3, 5, 7). So any improvement must come from **nonlinear rows that couple families**.

### 5b. Nonlinear, family-antisymmetric rows

Rows of the form `(y, −y)` are triorthogonal with the common Cartan rows iff y is *isotropic*: `Σ_q q_a y_q y'_q = 0` for a = 1…5 and every pair y, y' in the set.

- **Exactly anomaly-free (over the integers), entries in {−1, 0, 1}.** There are 221 zero-sum supports, 132,928 isotropic vectors, and 1,913,966 isotropic pairs (reduced by the S₅ symmetry). **No** pair protects all twin pairs. Of the 25 integer ternary rows compatible with the B5 rows, **none** is nonlinear.
- **Mod-p anomaly-free** (Z_p triorthogonality, which is all the code and its transversal gate need for p ≥ 5), extending B5. Results by prime:

| p | nonlinear isotropic rows | orbits | resulting codes |
|---|---:|---:|---|
| 5 | 0 | 0 | — |
| 7 | 588 | 17 | all [[30,2,3]]₇ with **8** weight-3 supports |
| 11 | 0 | 0 | — |
| 13 | 4056 | 57 | all [[30,2,3]]₁₃ with **8** weight-3 supports |

Nonlinear solutions appeared only for p ≡ 1 (mod 3). That is an observation, not a theorem.

**Best new code.** At p = 7, one extra row (for example `y = (4,6,6,6,6,1,0,5,5,0,1,0,4,4,4,4)` in the species order of `matter_qec_study.py`) gives `[[30,2,3]]₇` with X distance 9.
- It has only 8 weight-3 logical supports (4 per family, all through d^c₁), against 50 for B5.
- The exact leading distillation error is `P(any logical) = (2/9)ε³`, against `(25/18)ε³` for the B5 code at p = 7: **6.25× lower**.
- Caveats: the row is triorthogonal mod 7 but **not over the integers**, so it isn't a consistent continuous U(1) charge. It is also not constant on SM multiplets. It's a coding-theory object, not a physics one.

**Distance 4** was not found in any search. No compatible *second* nonlinear row exists at p = 7.

### 5c. Transversal entangling gates: none

Consider all transversal diagonal gates `⊗_j D^{w_j}` with `D = ω^{x³}` (cubic) or `ω^{x²}` (quadratic) and arbitrary per-qudit powers w_j. The gate preserves the code iff `T_w(x,y,s) = Σ_j w_j (x·h_j)(y·h_j)(s·h_j) = 0` for every s in the X-stabilizer coefficient space (for the quadratic case, the analogous bilinear condition). This holds for p ≥ 5, where cubic polynomial functions are faithful. The condition is *linear* in w, so the whole gate space can be computed exactly. The induced logical gate is the descended cubic or quadratic form on c-space, and it is entangling iff that form has mixed monomials.

| Code | valid cubic power vectors | logical rank | **entangling rank** | quadratic: valid / logical / entangling |
|---|---:|---:|---:|---|
| Direct sum, p = 5 | 2 | 2 | **0** | 10 / 2 / 0 |
| Coupled B5, 2 families, p = 5 and 7 | 2 | 2 | **0** | 12 / 2 / 0 |
| Coupled B5, 3 families, p = 5 | 3 | 3 | **0** | 19 / 3 / 0 |
| (B−L, T₃R), 2 / 3 families, p = 5 | 3 / 5 | 2 / 3 | **0** | 13 / 2 / 0, 21 / 3 / 0 |
| Nonlinear-row code, p = 7 | 1 | 1 | **0** | 8 / 2 / 0 |
| Family-number code, p = 5 | 6 | 2 | **0** | 16 / 2 / 0 |

No code supports a transversal diagonal entangling logical gate of degree ≤ 3, at the Clifford or the non-Clifford level. The cross-family coupling shows up in *errors* (Step 2), never in *gates*.

Not covered: non-diagonal transversal gates, transversal-plus-permutation ("fold") gates, and the qutrit ninth-root gates at p = 3.

## What changed in the overall picture

1. The coupling between families is real at the level of code equivalence (Step 3) and of correlated errors (Step 2). It is absent at the level of logical gates (Step 5c), and it never helps the leading distillation term.
2. The physics route closes more tightly than before:
   - the original smart badges violate the gauge structure itself;
   - the gauge-compatible replacement, family-dependent (B−L, T₃R), works as a code for p ≥ 5;
   - but realistic Yukawas break it, and the part that distance 3 needs is broken at O(1) by neutrino mixing.
3. Distance 4 is excluded for the whole Cartan-linear class. Nonlinear Z_p charges can thin out the weak spots (50 → 8) but did not remove them.

## Reproduce

```
pip install numpy
python3 open_steps.py            # ~4 min, ends with ALL OPEN-STEP CHECKS PASSED
python3 open_steps.py --quick    # skips the p = 13 search
```

The scripts check every claim at the stated primes. Statements "for all p ≥ 5" rest on the short proofs given above, together with explicit checks at p = 5, 7 and 11.
