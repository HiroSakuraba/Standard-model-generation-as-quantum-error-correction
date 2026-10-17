# Standard Model generations as quantum error correction

**The rule that keeps the Standard Model of particle physics mathematically consistent is the same rule a special quantum error-correcting code needs. A recent paper showed this for one family of particles. This project asks what happens with two families, then all three, and runs into a wall that nature itself puts up: flavour mixing.**

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

## What this does *not* show

- **Particles aren't computing anything.** Each particle is a label for a column in a table. Real hardware would need a physical system, measurements and noise modelling, and none of that follows from the mathematics.
- **The smart badges aren't real forces.** They are extra charges added by hand, and Part 3 shows real flavour mixing rules them out as exact symmetries. Many family-dependent charge rows also wouldn't fit consistently inside the Standard Model's gauge structure.
- **No proven advantage over separate codes.** The weakest (three-tile) errors stay inside single families, so the leading magic-state-distillation error rate equals that of independent 15-slot codes. The joined codes use fewer checks and have extra four-tile cross-family weak spots. Whether that trade helps is open. No decoder, threshold or circuit study has been done.
- **The magic operation acts on each digit separately.** It doesn't entangle them.
- **Novelty and review.** No literature search has been done for equivalent known codes, and the work hasn't been peer-reviewed.

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

## Reproduce

Requires Python 3.9+ and nothing else. Run both from the repository folder:

```bash
python3 matter_qec_study.py          # a few seconds  → ALL CHECKS PASSED
python3 three_generation_flavor.py   # about 20 seconds → ALL THREE-GENERATION CHECKS PASSED
```

Between them, the scripts check cubic anomaly cancellation in exact integer arithmetic, compute finite-field ranks, logical-qudit counts and distances, run the exhaustive 1,080-case p = 3 minimality search, and enumerate every logical support of size 1–4 in the three-generation codes over F₃ and F₅. They also find X distances exactly and check the two-flavour mixing identity with exact fractions (cos θ = 4/5, sin θ = 3/5 gives an off-diagonal of 24/25; this is a toy angle, not a physical rate).

## Open next steps

- Work in a concrete flavour Lagrangian (actual Yukawa matrices and spurions) to see which approximate symmetries survive, and whether an approximate code makes sense.
- Classify the correlated logical channels, not just supports: exact Pauli-coefficient enumerators and distillation acceptance/output polynomials, including the four-tile cross-family events.
- Classify the `[[30,2,3]]_p` and `[[45,3,3]]_p` codes against direct sums under monomial and local-Clifford equivalence.
- Look for an operational setting (e.g. quantum simulation) where these charge constraints could actually be measured.
- Search other anomaly-free charge lattices for distance 4 or a transversal *entangling* logical gate.

## The whole story in one chain

Anomaly cancellation = triorthogonality → one generation is a 15-slot magic code → the second generation alone is the identical code (charges, not masses) → two together hit the twins problem → a family badge gives 2 digits but distance 2 → smart badges restore distance 3 → three families with independent labels give 45 slots, 3 digits, distance 3, coupled only at four tiles → but flavour mixing blends the families, so no family-distinguishing badge can be exact, and the digits collapse back to one.

Two levels of description each refuse to supply a stable family label. At the level of symmetry (charges), the generations are perfect twins. At the level of mass, they are blended. A multi-family code needs exactly what neither level provides.
