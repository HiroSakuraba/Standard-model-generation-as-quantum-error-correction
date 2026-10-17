# Standard Model generations as quantum error correction

**The rule that keeps the Standard Model of particle physics mathematically consistent is the same rule a special quantum error-correcting code needs. A recent paper showed this for one family of particles. This project asks what happens with two families, then all three, and runs into a wall that nature itself puts up: flavour mixing. A follow-up works through the open questions that raised.**

📖 **Plain-language guide with an interactive "sneak an error past the checks" game:** [hirosakuraba.github.io/Standard-model-generation-as-quantum-error-correction](https://hirosakuraba.github.io/Standard-model-generation-as-quantum-error-correction/)

> Independent exploratory study by Benjamin John Schulz, October 2026. Not peer-reviewed.
> Builds on Shiroman Prakash & K. L. Ramya, *One generation of the Standard Model as a quantum error-correcting code*, [arXiv:2610.21879](https://arxiv.org/abs/2610.21879) ([their code](https://github.com/spdei/standard-model-code)).

---

## The background, in plain terms

**Generations.** Ordinary matter is built from up and down quarks, the electron, and the electron neutrino. Count every variety, with each quark in three "colours" and separate left- and right-handed versions, and you get 16 pieces: one *generation*. Nature repeats the set twice as heavier copies: charm, strange, muon and muon neutrino, then top, bottom, tau and tau neutrino. The copies feel identical forces; only their masses differ.

**The surprise.** Each particle carries charges, numbers saying how it responds to the forces. For the theory to be consistent, the cubes of those charges must add up to exactly zero (*anomaly cancellation*). Prakash and Ramya noticed that if you write the charges as a table, one column per particle, this rule is *exactly* the condition (*triorthogonality*) for the table to define a quantum error-correcting code with a built-in "magic" operation. Pull out the right-handed neutrino, and the other 15 particles become the 15 slots of the well-known code used in magic-state distillation, a key ingredient for fully capable quantum computers.

**Error-correcting codes in one paragraph.** Spread a secret across many lockers and add checksums. A code's *distance* is the fewest lockers a tamperer must touch to change the secret without tripping any check. Distance 3 means any single fault can be found and fixed. A *qudit* is a quantum digit with *p* settings, where *p* is a prime (qubit: p = 2, qutrit: p = 3, and so on).

## Part 1: two generations

**The second family on its own** gives *exactly the same* 15-slot code. The code reads charges, not masses, and the muon family's charges match the electron family's column for column. Mass is invisible at this level.

**Both families together** gives 30 slots, with both right-handed neutrinos pulled out, and the hope is two protected digits. The obstacle is that every particle has an identical-looking twin in the other family. Three ways of handling it:

| Approach | Result | Problem |
|---|---|---|
| **A.** Shared charges only | 1 protected digit | The two hidden spots look identical, so they can only hold one secret |
| **B.** Add a "family badge" (+1 for family 1, −1 for family 2) | 2 digits, **distance 2** | Nudging a particle one way and its twin the opposite way slips past every check |
| **C.** **"Smart badges"**: family-flipped copies of some charges, combined carefully | **2 digits, distance 3** | None of these twin-pair blind spots. Needs 3 extra rows for p = 3, and only 2 for every prime p ≥ 5 |

**Minimality.** One extra row never works. For p = 3 the script checks all 1,080 ways of adding two rows; every one leaves a blind spot. So three rows is minimal *within this recipe* (keep the five shared charges and add family-flipped rows). It is not a claim about every possible code.

## Part 2: three generations

Give each family a two-number label instead of a single ±1: **(+1, 0)**, **(−1, 0)** and **(0, +1)**. These point in genuinely different directions for every odd prime, so all three hidden neutrino spots can be told apart. Then apply the smart-badge rows along each direction.

- **Result:** 45 slots, **3 protected digits, distance 3**. That's `[[45,3,3]]_3` with 11 charge rows and 8 checks, and `[[45,3,3]]_5` with 9 rows and 6 checks. Distance 3 holds for every odd prime, and the bit-flip-type distance is 7 (p = 3) and 9 (p = 5). Anomaly cancellation holds exactly over the integers.
- **The three-tile weak spots stay inside single families:** 25 per family, 75 in all, none crossing.
- **But the families really are coupled, at four tiles.** There are silent errors using just **two tiles in each of two families**: 88 per family pair at p = 3 and 396 at p = 5. Three *separate* codes can't have these. So the joined code is structurally different from stacking copies, but the coupling sits one tile beyond the distance-3 guarantee.

## Part 3: the flavour-mixing catch

The particles with a definite **mass** aren't the ones with a definite **family label**. Each is a blend (this is measured: the CKM matrix for quarks and the PMNS matrix for neutrinos). A blend can't carry a fixed family badge.

Concretely, rotate two families by a mixing angle θ. A diagonal badge `diag(d₁, d₂)` picks up off-diagonal entries `(d₁ − d₂)·sin θ·cos θ`. These vanish only if θ = 0 (no mixing) or d₁ = d₂ (the badge can't tell the families apart). So **any badge that distinguishes the families can't be an exact symmetry** while there is real mixing, unless new physics protects it. Keep only the family-blind charges, and the three hidden spots collapse back to **one** protected digit.

The structure the multi-family code needs is precisely the structure flavour mixing erases.

### Key numbers

| Construction | p | Slots | Digits | Distance | Checks |
|---|---:|---:|---:|---:|---:|
| 2nd family alone | 3, 5 | 15 | 1 | 3 | 4 |
| Two families, shared charges | 5 | 30 | 1 | 3 | 4 |
| Two families + family badge | 5 | 30 | 2 | **2** | 4 |
| **Two families + smart badges (3 rows)** | 3 | 30 | 2 | **3** | 6 |
| **Two families + smart badges (2 rows)** | ≥ 5 | 30 | 2 | **3** | 5 |
| Two separate copies (baseline) | 3 | 30 | 2 | 3 | 8 |
| **Three families + smart badges** | 3 | 45 | 3 | **3** | 8 |
| **Three families + smart badges** | 5 | 45 | 3 | **3** | 6 |
| Three families, family-blind charges only | odd | 45 | 1 | – | – |

## Part 4: the follow-up on the open steps

The five open steps from earlier versions are now worked through in [`OPEN_STEPS_NOTE.md`](OPEN_STEPS_NOTE.md), with every number checked by [`open_steps.py`](open_steps.py). In plain terms:

1. **Real flavour physics.** The original smart badges have a deeper problem than mixing: they treat the three quark colours differently, so the strong force itself breaks them. The only badges the forces allow are family-by-family versions of **B−L** and **right-handed isospin**. Those *do* give working codes (`[[30,2,3]]` and `[[45,3,3]]` for every prime p ≥ 5). But an exhaustive scan of possible Yukawa couplings finds that none of these badges survives exactly in a realistic Standard Model. The part needed for distance 3 is broken at order one by neutrino mixing. The best approximately preserved structure gives only distance 2.
2. **Error channels, exactly.** Exact formulas for every logical error class. To leading order the joined codes behave like separate copies. But they produce **correlated errors that hit both outputs at once** at order ε⁴, which separate copies don't do until ε⁶. These come from the 2+2 cross-family weak spots.
3. **Genuinely different codes.** The joined codes cannot be turned into direct sums even by arbitrary local unitaries and qudit relabelling: the projector cannot factorize into smaller code blocks, as proved using the connected stabilizer-support components. This does not establish novelty among other indecomposable codes. Their stabilizer groups don't split.
4. **Where it could be measured.** In a Fourier-rotated basis, the checks are measurements of conserved charges mod p, and the protected digit is "the charge of the missing right-handed neutrinos." That fits engineered qudit hardware and gauge-theory simulators. It can't work with real gauge charges, because nature forbids superposing different charge sectors. The joined codes cost about 18 % more to measure. The ε⁴ correlated errors are a concrete experimental signature.
5. **Distance 4 and entangling gates.** Distance 4 is impossible for every lattice built this way (there's a short proof) and wasn't found in wider searches. Exotic mod-7 and mod-13 charges cut the weakest error patterns from 50 to 8, which lowers the leading distillation error 6.25× at p = 7, but they aren't consistent physical charges. No tested code has a diagonal transversal entangling logical gate within the degree-2 and degree-3 families searched.

## What this does *not* show

- **Particles aren't computing anything.** Each particle is a label for a column in a table. Real hardware would need a physical system, measurements and noise modelling, and none of that follows from the mathematics.
- **The smart badges aren't real forces.** They are extra charges added by hand. Part 3 shows flavour mixing rules them out as exact symmetries, and Part 4 shows the original B3/B5 rows even conflict with the strong and weak forces. The force-compatible replacement (family B−L and right-handed isospin) is also broken by realistic Yukawa couplings.
- **No advantage over separate codes.** The leading magic-state-distillation error rate equals that of independent 15-slot codes. The joined codes add correlated two-output errors at order ε⁴ and cost about 18 % more to measure. No decoder, threshold or circuit-level study has been done.
- **The magic operation acts on each digit separately.** It doesn't entangle them.
- **Novelty and review.** The joined codes are proven inequivalent to direct sums, but haven't been compared against other known qudit codes. None of this work has been peer-reviewed.

## Files

| File | What it is |
|---|---|
| [`index.html`](index.html) | Plain-language explainer: interactive error game (2 or 3 families) and a flavour-mixing slider (the GitHub Pages site) |
| [`RESEARCH_NOTE.md`](RESEARCH_NOTE.md) | Technical note, part 1: two-generation construction, distance-3 proofs, minimality, limitations |
| [`THREE_GENERATION_FLAVOR_NOTE.md`](THREE_GENERATION_FLAVOR_NOTE.md) | Technical note, part 2: three-generation codes, weight-four cross-family errors, the flavour-mixing obstruction |
| [`matter_qec_study.py`](matter_qec_study.py) | Exact-arithmetic verification for the one- and two-generation results |
| [`three_generation_flavor.py`](three_generation_flavor.py) | Exact-arithmetic verification for the three-generation and mixing results (imports `matter_qec_study.py`) |
| [`test_results.txt`](test_results.txt), [`three_generation_results.txt`](three_generation_results.txt) | Recorded outputs of completed runs |
| [`interactive_charge_map.html`](interactive_charge_map.html) | Technical explorer: two-generation charge matrices with per-species pair diagnostics |
| [`interactive_flavor_map.html`](interactive_flavor_map.html) | Technical explorer: three-generation charges, flavour-rotation demo, support counts |
| [`OPEN_STEPS_NOTE.md`](OPEN_STEPS_NOTE.md) | Technical note, part 3: flavour Lagrangian, exact logical channels, equivalence, operational setting, distance-4 / entangling-gate searches |
| [`open_steps.py`](open_steps.py), [`charge_lattice_tools.py`](charge_lattice_tools.py) | Verification script for part 3 (asserts every reported number) and its exact mod-p toolkit (needs numpy) |
| [`open_steps_results.txt`](open_steps_results.txt) | Recorded output of a completed run |
| [`OPERATIONAL_GAUGE_QEC.md`](OPERATIONAL_GAUGE_QEC.md) | Exact Knill–Laflamme recovery versus gauge constraints: what works for a simulator and what remains unproved for the Standard Model |
| [`operational_gauge_qec.py`](operational_gauge_qec.py) | Pure-Python gauge-ring Gauss-only failure and 3×3 lattice-gauge toric-code recovery verification |
| [`abstract_stabilizer_recovery.py`](abstract_stabilizer_recovery.py) | Full single-site Weyl syndrome and recovery verification for the original charge CSS codes (requires NumPy) |
| [`SU2_MATTER_GAUSS_KL.md`](SU2_MATTER_GAUSS_KL.md) | Explicit SU(2) matter+rishon color-singlet basis and failed Knill–Laflamme test for local gauge-invariant noise |
| [`su2_matter_singlet.py`](su2_matter_singlet.py) | Exact 64-state, two-vertex SU(2) Gauss law kernels, singlet wavefunctions, and logical-Z counterexample (requires SymPy) |

## Reproduce

To reproduce the new deterministic and comparison checks:

```bash
python3 exact_yukawa.py
python3 code_fingerprints.py            # requires numpy
python3 code_fingerprints.py --extended # adds three-generation codes
python3 local_purity.py                 # exact four-site reduced-projector purities; requires numpy
python3 published_qutrit_family.py     # reproduces external 2025 paper matrices
python3 reed_muller_family_exclusion.py # source-scoped 2026 family no-go and exact charge moments
python3 walsh_character_lattice.py       # all-prime rigidity and masked direct-sum theorem
python3 operational_gauge_qec.py         # Gauss-only failure vs genuine lattice-gauge recovery
python3 abstract_stabilizer_recovery.py  # all single-site Pauli errors on charge-code qudits
python3 su2_matter_singlet.py           # SU(2) with matter: Gauss-invariant error fails KL (requires sympy)
```

Existing checks:

Requires Python 3.9+. The first two scripts need nothing else; the third needs numpy. Run from the repository folder:

```bash
python3 matter_qec_study.py          # a few seconds  → ALL CHECKS PASSED
python3 three_generation_flavor.py   # about 20 seconds → ALL THREE-GENERATION CHECKS PASSED
python3 open_steps.py                # about 4 minutes → ALL OPEN-STEP CHECKS PASSED  (--quick skips the p = 13 search)
```

Between them, the scripts check cubic anomaly cancellation in exact integer arithmetic, compute finite-field ranks, logical-qudit counts and distances, run the exhaustive 1,080-case p = 3 minimality search, and enumerate every logical support of size 1–4 in the three-generation codes over F₃ and F₅. They also find X distances exactly and check the two-flavour mixing identity with exact fractions (cos θ = 4/5, sin θ = 3/5 gives an off-diagonal of 24/25; this is a toy angle, not a physical rate).

## Open next steps

The original five steps are addressed in Part 4 and [`OPEN_STEPS_NOTE.md`](OPEN_STEPS_NOTE.md). What's left:

- **New:** [exact Yukawa theorem and classification fingerprints](EXACT_SYMMETRY_AND_CLASSIFICATION.md), with [deterministic single-Higgs verification](exact_yukawa.py) and [code-fingerprint comparison](code_fingerprints.py). The latter distinguishes the joined codes from direct sums but does not establish novelty among all qudit codes.
- **New:** [four-qudit reduced-state purity invariant](LOCAL_PURITY_NOTE.md), computed by [local_purity.py](local_purity.py). For the two-family coupled codes there are 184 (p=3) or 249 (p=5) four-site purity hotspots, versus 150 for independent copies. The three-family counts are 327/522, versus 225. These are local-unitary invariants, not logical error rates; matching existing indecomposable qudit codes remains open.
- **Published matrix comparison:** [Prakash–Saha 2025 qutrit code family, exact reconstruction and comparison](PUBLISHED_CODE_COMPARISON.md), verified with [published_qutrit_family.py](published_qutrit_family.py). Their m=4,k=6 code has length 30 and six X checks, but encodes six logical qutrits at distance 2; the coupled matter construction encodes two at distance 3. The entire specified [[9m−k,k,2]]₃ family is excluded by parameters, **not** all odd-prime Reed–Muller/triorthogonal codes.
- **Further source-checked exclusion:** [Saha–Prakash 2026 full-evaluation punctured Reed–Muller comparison](FULL_REED_MULLER_COMPARISON.md), verified by [reed_muller_family_exclusion.py](reed_muller_family_exclusion.py). Under that paper's full-rank puncture construction, n+k must be p^m; ours have n+k=32 or 48, neither an odd-prime power at p=3 or 5. The original charge lattice is cubic-orthogonal but not pairwise self-orthogonal, unlike the unpunctured construction in that paper. This **does not exclude arbitrary Reed–Muller variants or all known qudit codes**.
- **Walsh-character theorems:** [exact single-family rigidity and the disguised direct-sum construction](WALSH_CHARACTER_THEOREMS.md), independently checked by [walsh_character_lattice.py](walsh_character_lattice.py). No additional single-family row can preserve all cubic Cartan constraints over any odd prime. All compatible five-row, pure-Walsh two-family extensions are exactly two original 15-qudit codes after a coordinate permutation—even when 44 of their 50 minimum-weight logical supports span the *original* family labels. Cross-family support counts alone are not proofs of genuine coupling; the earlier B3/B5 indecomposability proof remains essential.
- **Operational interpretation:** [physical gauge-state and Knill–Laflamme benchmark](OPERATIONAL_GAUGE_QEC.md) shows that the *formal* generation-based codes correct arbitrary single-qudit register errors, but this is not physical Standard Model protection. A three-link gauge ring fails for local gauge-invariant phase noise; an engineered 3×3 toric code succeeds only with additional magnetic checks. A physically justified SM-Hilbert-space encoding and measurable SM-compatible syndromes remain unconstructed.
- **Next non-Abelian matter result:** [two-vertex SU(2) matter singlets and exact Gauss-only KL obstruction](SU2_MATTER_GAUSS_KL.md), verified in [su2_matter_singlet.py](su2_matter_singlet.py). The full color-singlet sector of a deliberately truncated 64-state rishon/matter model has dimension two, but a link-local SU(2)-invariant electric Casimir distinguishes its logical states and produces an undetectable phase error. This is *not* a Standard Model encoding or a general SU(2) no-go theorem.
- Compare the joined codes against **same-(n,k,p)** explicit check matrices of other published odd-prime codes, especially shortened/affine-variety constructions. That broader novelty question is open.
- Explain why nonlinear charge rows appear only for p ≡ 1 (mod 3), and whether the floor of 8 weight-3 weak spots can be broken.
- Extend the gate search to non-diagonal and permutation-assisted ("fold") transversal gates, and to the p = 3 ninth-root gates.
- Run a circuit-level simulation, with noisy checks of weight 15–16, to see whether the ε⁴ correlated errors matter in a full distillation factory.

## The whole story in one chain

Anomaly cancellation = triorthogonality → one generation is a 15-slot magic code → the second generation alone is the identical code (charges, not masses) → two together hit the twins problem → a family badge gives 2 digits but distance 2 → smart badges restore distance 3 → three families with independent labels give 45 slots, 3 digits, distance 3, coupled only at four tiles → but flavour mixing blends the families, so no family-distinguishing badge can be exact, and the digits collapse back to one → the follow-up shows the only force-compatible badges (family B−L and right-handed isospin) still need charges that neutrino mixing breaks at order one, that the joined codes are inequivalent to independent copies whose coupling shows up as correlated ε⁴ errors rather than new gates, and that distance 4 is out of reach for this whole class.

Two levels of description each refuse to supply a stable family label. At the level of symmetry (charges), the generations are perfect twins. At the level of mass, they are blended. A multi-family code needs exactly what neither level provides.
