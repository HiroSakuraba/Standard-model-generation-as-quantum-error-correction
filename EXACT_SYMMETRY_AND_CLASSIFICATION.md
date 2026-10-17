# Exact Yukawa symmetry theorem and code-comparison fingerprints

**Status (17 October 2026):** exploratory follow-up, not peer-reviewed. These charge-lattice CSS codes are mathematical objects. No physical encoding, recovery map, or dynamical protection is asserted.

## 1. Replace the stochastic one-Higgs screen with an exact argument

The original `open_steps.py` enumerates discrete family charges and allowed Yukawa textures, then uses random complex matrices and singular-value decompositions to test generic mixing. That was useful exploratory evidence but is not a deterministic proof.

**Proposition (exact unbroken, family-diagonal symmetry; single-Higgs model).** Let `G >= 2` quark families have an exactly conserved, unbroken family-diagonal Abelian symmetry, with charges in the gauge-compatible span of `S3 = 3(B-L)` and `S2 = -4 T3R`. Write the generation parameters as `(a_i,b_i)` over the rationals, reals, or a field of odd prime characteristic. Suppose the up- and down-quark mass matrices have full rank and their left mass matrices admit **no nontrivial common family-block decomposition** (generic non-block-diagonal CKM mixing). Then `a_i` and `b_i` are independent of generation.

**Proof.** The left quark doublet carries charge `a_i`; `u^c_j` and `d^c_j` carry `-a_j+2b_j` and `-a_j-2b_j`. An up-type Yukawa entry requires

```
a_i - a_j + 2 b_j + eta_up = 0
```

where `eta_up` is generation-independent. Down-type entries obey the corresponding rule with `-2b_j` and the appropriate Higgs charge. Invariance implies `H_u=Y_u Y_u†` and `H_d=Y_d Y_d†` both commute with `diag(a_i)`. If there are distinct `a_i`, they share a nontrivial block decomposition into eigenspaces of that charge, contrary to the mixing hypothesis. Therefore every `a_i=a`. Now each column `j` of `Y_u` is either fully allowed or identically zero depending on `2 b_j + eta_up`. Full rank forbids a zero column. Since 2 is invertible, all `b_j` are equal. QED.

**Scope.** This is about exact *unbroken* family charge selection rules. If multiple Higgs doublets have different charges and both acquire vacuum expectation values, the candidate family symmetry is broken; it must not be described as a surviving exact low-energy conservation law. Spurions, extra particles and non-Abelian flavor groups require their own analysis. This proposition is not a no-go theorem for every two-Higgs or beyond-Standard-Model theory.

### Deterministic verification

Run `python3 exact_yukawa.py` (standard-library only). It enumerates all single-Higgs charge assignments for G=2 at p=3,5,7 and G=3 at p=3,5.

- For each Yukawa support, **full rank is possible exactly when its bipartite graph has a perfect matching**.
- If the union of left-row overlap graphs from up and down Yukawas is disconnected, then both `Y Y†` matrices are block-diagonal; connectivity is a **necessary**, not sufficient, condition for generic CKM mixing.
- All surviving assignments are family-universal: exactly `p²` in each enumerated case. No numerical random matrices or tolerance thresholds are used.

The existing `open_steps.py` 2HDM screen remains a finite, restricted-model computational search involving numerical trials, not a general theorem.

## 2. Comparison fingerprints: distinguish direct sums, not all known codes

Run `python3 code_fingerprints.py` (requires NumPy and `charge_lattice_tools.py`). Use `--extended` to include the three-generation codes. It reports:

1. qudit dimension p, code length n, logical dimension k, X and Z stabilizer ranks;
2. the **complete X-stabilizer Hamming-weight enumerator**;
3. connected components of the Pauli-stabilizer support decomposition using the repo's existing `split_components` routine;
4. weight-three Z-logical support partitions by family;
5. genuinely mixed 2+2 weight-four Z-logical support counts.

The X-stabilizer enumerator is invariant under coordinate permutations and local nonzero rescalings **preserving the CSS X structure**. It is not invariant under unrestricted local Clifford transformations that exchange X and Z.

The Pauli-stabilizer support-component sizes distinguish the tested coupled models (one 30-coordinate component, or one 45-coordinate component) from the tested direct sums (two 15-coordinate components). The existing `OPEN_STEPS_NOTE.md` supplies the stronger **arbitrary-local-unitary** argument: the stabilizer-code projector factorizes across a bipartition precisely when its Pauli support (the stabilizer group) factorizes there. A local tensor-product unitary preserves projector factorization, so an indecomposable joined code cannot become a direct sum under local unitaries and a coordinate permutation. This does **not** prove inequivalence to every other indecomposable qudit code.

The weight-four mixed 2+2 Z-logical support counts are 88 per family pair at p=3 and 396 at p=5 for the reported coupled examples, against zero for direct sums. This is an explicit useful comparison, not a claim of complete classification.

### Relevant published / preprint comparison families

- Prakash and Ramya, [*One generation of the Standard Model as a quantum error-correcting code*](https://arxiv.org/abs/2610.21879) (2026): source construction.
- Nezami and Haah, [*Classification of small triorthogonal codes*](https://doi.org/10.1103/PhysRevA.106.012437) (2022): **qubit** codes, `n+k <= 38`.
- Wills, Jain and Singh, [*Classification of Generalised Triorthogonal Codes through Length 54*](https://arxiv.org/abs/2609.30860) (2026): enlarged **qubit T-state** classification; its length coverage does *not* constitute a p=3 or p=5 qudit catalogue.
- Prakash, [*Cyclic triorthogonal codes in prime dimension*](https://arxiv.org/abs/2610.08012) (2026): odd-prime **cyclic** subclass.
- Saha and Prakash, [*Sublogarithmic distillation in all prime dimensions using punctured Reed-Muller codes*](https://doi.org/10.1103/tzkz-9p2f) (2026): odd-prime comparison family; equivalence against specific check matrices remains untested.

**Open:** retrieve external qudit code generator/check matrices, normalize conventions for non-self-orthogonal triorthogonality, compare invariants, and then attempt canonical equivalence. Until that work is completed, describe the family-coupled codes as **non-direct-sum constructions whose novelty relative to other qudit codes remains unverified**.

## Reproduce

```bash
python3 exact_yukawa.py          # standard library; deterministic
python3 -m pip install numpy
python3 code_fingerprints.py     # two-family codes
python3 code_fingerprints.py --extended  # includes three-family codes
```

These tests supplement rather than replace the existing `matter_qec_study.py`, `three_generation_flavor.py` and `open_steps.py` checks.
