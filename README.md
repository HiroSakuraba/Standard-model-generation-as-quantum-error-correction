# Standard Model generations as quantum error correction

**The rule that keeps the Standard Model of particle physics mathematically consistent is the same rule a special quantum error-correcting code needs. A recent paper showed this for one family of particles. This project asks what happens with two.**

📖 **Plain-language guide with an interactive "sneak an error past the checks" game:** [hirosakuraba.github.io/Standard-model-generation-as-quantum-error-correction](https://hirosakuraba.github.io/Standard-model-generation-as-quantum-error-correction/)

> Independent exploratory study by Benjamin John Schulz, October 2026. Not peer-reviewed.
> Builds on Shiroman Prakash & K. L. Ramya, *One generation of the Standard Model as a quantum error-correcting code*, [arXiv:2610.21879](https://arxiv.org/abs/2610.21879) ([their code](https://github.com/spdei/standard-model-code)).

---

## The background, in plain terms

**Generations.** Ordinary matter is built from up and down quarks, the electron, and the electron neutrino. Count every variety, with each quark in three "colours" and separate left- and right-handed versions, and you get 16 pieces: one *generation*. Nature then repeats the set twice as heavier copies. The second generation is charm, strange, the muon and the muon neutrino. The copies feel identical forces; only their masses differ.

**The surprise.** Each particle carries charges, numbers saying how it responds to the forces. For the theory to be consistent, the cubes of those charges must add up to exactly zero (*anomaly cancellation*). Prakash and Ramya noticed that if you write the charges as a table, one column per particle, this rule is *exactly* the condition (*triorthogonality*) for the table to define a quantum error-correcting code with a built-in "magic" operation. Pull out the right-handed neutrino, and the other 15 particles become the 15 slots of the well-known code used in magic-state distillation, a key ingredient for fully capable quantum computers.

**Error-correcting codes in one paragraph.** Spread a secret across many lockers and add checksums. A code's *distance* is the fewest lockers a tamperer must touch to change the secret without tripping any check. Distance 3 means any single fault can be found and fixed. A *qudit* is a quantum digit with *p* settings, where *p* is a prime (qubit: p = 2, qutrit: p = 3, and so on).

## What this project does

**1. The second family on its own** gives *exactly the same* 15-slot code. The code reads charges, not masses, and the muon family's charges match the electron family's column for column. Mass is invisible at this level.

**2. Both families together** is where it gets interesting. That gives 30 slots, with both right-handed neutrinos pulled out, and the hope is two protected digits. The obstacle is that every particle has an identical-looking twin in the other family. Three ways of handling it:

| Approach | Result | Problem |
|---|---|---|
| **A.** Shared charges only | 1 protected digit | The two hidden spots look identical, so they can only hold one secret |
| **B.** Add a "family badge" (+1 for family 1, −1 for family 2) | 2 digits, **distance 2** | Nudging a particle one way and its twin the opposite way slips past every check |
| **C.** **"Smart badges"** (this study): family-flipped copies of some charges, combined carefully | **2 digits, distance 3** | None of these twin-pair blind spots. Needs 3 extra rows for p = 3, and only 2 for every prime p ≥ 5 |

**3. Minimality.** One extra row never works. For p = 3 the script checks all 1,080 ways of adding two rows; every one leaves a blind spot. So three rows is minimal *within this recipe* (keep the five shared charges and add family-flipped rows). It is not a claim about every possible code.

### Key numbers

| Construction | p | Slots | Digits | Distance | Checks |
|---|---:|---:|---:|---:|---:|
| 2nd family alone | 3, 5 | 15 | 1 | 3 | 4 |
| Two families, shared charges | 5 | 30 | 1 | 3 | 4 |
| Two families + family badge | 5 | 30 | 2 | **2** | 4 |
| **Two families + smart badges (3 rows)** | 3 | 30 | 2 | **3** | 6 |
| **Two families + smart badges (2 rows)** | ≥ 5 | 30 | 2 | **3** | 5 |
| Two separate copies (baseline) | 3 | 30 | 2 | 3 | 8 |

In code notation, the smart-badge codes are `[[30,2,3]]_p`. Their distance against the other error type (X distance) is 7 at p = 3 and 9 at p = 5.

## What this does *not* show

- **Particles aren't computing anything.** Each particle is a label for a column in a table. Real hardware would need a physical system, measurements and noise modelling, and none of that follows from the mathematics.
- **The smart badges aren't real forces.** They are extra charges added by hand. Real quarks and neutrinos mix between families, which would break them unless new physics protects them.
- **No proven advantage over two separate codes.** The weakest error patterns stay entirely inside one family (25 + 25, none mixed). So the leading magic-state-distillation error rate equals that of two independent 15-slot codes. The combined code uses fewer checks (5–6 vs 8); whether that helps is open.
- **The magic operation acts on each digit separately.** It doesn't entangle them.
- **Novelty and review.** No literature search has been done for equivalent known codes, and the work hasn't been peer-reviewed.

## Files

| File | What it is |
|---|---|
| [`index.html`](index.html) | Plain-language explainer with the interactive error game (the GitHub Pages site) |
| [`interactive_charge_map.html`](interactive_charge_map.html) | Technical explorer: full charge matrices for each construction, with per-species pair diagnostics |
| [`RESEARCH_NOTE.md`](RESEARCH_NOTE.md) | The full technical note: precise construction, proofs of distance 3, minimality argument, limitations, next steps |
| [`matter_qec_study.py`](matter_qec_study.py) | Standalone exact-arithmetic verification script (standard library only) |
| [`test_results.txt`](test_results.txt) | Output of a completed run |

## Reproduce

Requires Python 3.9+ and nothing else. Runs in a few seconds.

```bash
python3 matter_qec_study.py
```

The script checks cubic anomaly cancellation in exact integer arithmetic, computes finite-field ranks, logical-qudit counts and Z distances, finds X distances exhaustively for the p = 3 and p = 5 coupled codes, runs the exhaustive 1,080-case p = 3 minimality search, and cross-checks the binary (p = 2) Reed–Muller case. It ends with `ALL CHECKS PASSED`.

## Open next steps

- Exact weight enumerators and acceptance/error polynomials for distillation with the coupled codes.
- Classify the `[[30,2,3]]_p` codes against the direct sum under monomial and local-Clifford equivalence.
- Model flavour breaking (CKM/PMNS mixing) explicitly: can an approximate symmetry support an approximate code?
- Search other anomaly-free charge assignments (cubic *and* mixed) for distance 4 or a transversal *entangling* logical gate.

## The whole story in one chain

Anomaly cancellation = triorthogonality → one generation is a 15-slot magic code → the second generation alone is the identical code (charges, not masses) → two together hit the twins problem → a family badge gives 2 digits but distance 2 → smart badges (3 rows at p = 3, 2 rows at p ≥ 5) restore distance 3 → but there is no distillation gain over two separate codes yet, and the badges aren't observed physics.

The code lives entirely at the level of symmetry, where only charges exist. At that level the generations are perfect twins, so getting a genuinely two-family code means adding structure by hand that nature, as far as we know, doesn't supply.
