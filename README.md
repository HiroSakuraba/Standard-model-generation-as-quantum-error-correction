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
| [`FOUR_LINK_SU2_DFS.md`](FOUR_LINK_SU2_DFS.md) | Constructive two-dimensional SU(2) matter/rishon decoherence-free subspace against occupation noise, and exact failure for local color exchange |
| [`four_link_su2_subspace.py`](four_link_su2_subspace.py) | Exact four-link singlet states, Gauss checks, noise algebra and projected local SU(2) spin-exchange matrices (requires SymPy) |
| [`SU2_EXCHANGE_HAMILTONIAN_AND_NO_GO.md`](SU2_EXCHANGE_HAMILTONIAN_AND_NO_GO.md) | Exact gapped two-dimensional SU(2) ground-space Hamiltonian, three exchange-coupling balance conditions, and complete local-algebra Knill–Laflamme obstruction |
| [`su2_exchange_algebra.py`](su2_exchange_algebra.py) | Verifies all eight spin-singlet sectors, on-site swap matrices, Hamiltonian gap, and full endpoint/exchange operator algebra (requires SymPy) |
| [`FIVE_BLOCK_SU2_RECOVERY.md`](FIVE_BLOCK_SU2_RECOVERY.md) | Five independent SU(2) singlet modules plus the known [[5,1,3]] outer code; exact recovery proof, physical limitations and two-module fault |
| [`five_block_su2_recovery.py`](five_block_su2_recovery.py) | 1,024 exact Knill–Laflamme overlaps, 16 error syndromes, local SU(2) exchange-to-Pauli map, correlated-noise counterexample |
| [`FLAGGED_SU2_STABILIZER_NOTE.md`](FLAGGED_SU2_STABILIZER_NOTE.md) | Exhaustive single-fault Clifford extraction-gadget comparison, flag placement and recovery ambiguities; not a full fault-tolerant measurement protocol |
| [`flagged_su2_measurement.py`](flagged_su2_measurement.py) | Pure-Python exact propagation of 64 bare / 98 flagged fault cases for each outer check and conditional recovery validation |
| [`MEASURED_MULTIRROUND_SU2_NOTE.md`](MEASURED_MULTIRROUND_SU2_NOTE.md) | Two full measured flagged rounds plus adaptive third; 784 single-fault containment cases and an exact failure/observation-collision result for combined input+fault noise |
| [`measured_multiround_su2.py`](measured_multiround_su2.py) | Standalone exact Pauli circuit/feedback simulator with only observed syndrome bits, optional 11,760-case combined-fault stress diagnostic |
| [`ADAPTIVE_FLAG_EC_AUDIT.md`](ADAPTIVE_FLAG_EC_AUDIT.md) | Source-checked Chao–Reichardt *adaptive decision tree* with explicit non-identical circuit and correction of one-fault release criteria |
| [`paper_style_flag_ec.py`](paper_style_flag_ec.py) | Observed-only adaptive flagged-until-trigger then full unflagged extraction; 392 clean-input fault cases, 15 incoming-only errors, additional two-event stress |
| [`SU2_CONTROLLED_EXCHANGE_IMPLEMENTATION.md`](SU2_CONTROLLED_EXCHANGE_IMPLEMENTATION.md) | Gauge-invariant ancilla-controlled exchange Hamiltonians on four-link SU(2) modules, exact Clifford bridge, physical limitations, and primary-source Figure 2 correction |
| [`su2_controlled_exchange.py`](su2_controlled_exchange.py) | 12,288 full-kinematic SU(2) Gauss commutator checks and exact ancilla-controlled X/Z Clifford verification on the singlet subspace |
| [`published_flag_figure2_catalog.py`](published_flag_figure2_catalog.py) | Independently checks the seven distinct nontrivial syndromes in Chao–Reichardt arXiv:1705.02329v1 Figure 2(d) |
| [`TWO_BODY_SU2_MEDIATED_GATE.md`](TWO_BODY_SU2_MEDIATED_GATE.md) | Schur-lemma limitation for neutral ancilla–single fundamental coupling, explicit two-body exchange-mediated ZZ, full eight-spin leakage and miscalibration tests |
| [`two_body_su2_mediator.py`](two_body_su2_mediator.py) | Exact SU(2) invariants and second-order Hamiltonian, 256-state pulse simulation, coherent single-exchange calibration faults |
| [`SU2_SMOOTH_PULSE_OPTIMIZATION.md`](SU2_SMOOTH_PULSE_OPTIMIZATION.md) | Pulse-shaping experiment, 256-state vs 14-singlet consistency, convergence, coherent gate deviations and exchange calibration sensitivity |
| [`su2_shaped_pulse.py`](su2_shaped_pulse.py) | Deterministic bounded smooth-pulse optimization, 1024-bin verification, SU(2) invariant-sector reduction and full-spin independent crosscheck |
| [`SU2_SK1_COMPOSITE_PULSE.md`](SU2_SK1_COMPOSITE_PULSE.md) | SK1-inspired SU(2)-invariant exchange-composite pulse, all calibrated error cases, full-spin audit and ninefold duration cost |
| [`su2_sk1_composite_pulse.py`](su2_sk1_composite_pulse.py) | Full 8-spin modeled noncommuting SU(2) singlet control, independent coupler faults, 256/14D crosscheck, pulse convergence |
| [`SU2_MINIMAL_CONTROL_DECISION.md`](SU2_MINIMAL_CONTROL_DECISION.md) | Worst-case simple-versus-composite control budget, modest waveform search, analytic weak-coupling sensitivity limitation |
| [`su2_minimal_pulse_audit.py`](su2_minimal_pulse_audit.py) | Exact 14D microscopic Hamiltonian simulation, both cross-coupler and extra-local control drift, fixed nominal comparisons |
| [`SU2_ALLOWED_ERROR_ALGEBRA.md`](SU2_ALLOWED_ERROR_ALGEBRA.md) | Exact four/eight-spin gauge-invariant exchange algebra and restricted-error Knill–Laflamme tests |
| [`su2_error_algebra_capacity.py`](su2_error_algebra_capacity.py) | Exact singlet algebra verifier (SymPy) |
| [`SU2_EXCHANGE_GRAPH_CAPACITY.md`](SU2_EXCHANGE_GRAPH_CAPACITY.md) | Exchange-graph geometry and constructive 14D-singlet qutrit code |
| [`su2_bridge_error_geometry.py`](su2_bridge_error_geometry.py), [`su2_graph_qutrit_witness.py`](su2_graph_qutrit_witness.py) | Exact bridge-error tests and explicit gauge-invariant qutrit witness |
| [`SU2_OFFGRAPH_ROBUSTNESS.md`](SU2_OFFGRAPH_ROBUSTNESS.md) | Weak parasitic exchange: exact error-algebra obstruction vs finite-time, recoverable logical fidelity |
| [`su2_offgraph_robustness.py`](su2_offgraph_robustness.py) | Full 256D spin evolution, independent 14D crosscheck and explicit CPTP erasure recovery |
| [`SU2_NATIVE_VERTEX_CONTROL_AUDIT.md`](SU2_NATIVE_VERTEX_CONTROL_AUDIT.md) | Published exchange-enabled quantum-link vertex comparison: six-pair Gauss enforcement, coupling mismatch and selective finite-pulse echo |
| [`su2_native_vertex_controls.py`](su2_native_vertex_controls.py) | Exact 16-state SU(2) Casimir/exchange tests, numerical finite-pulse channel and independently checked full-spin evolution |
| [`SU3_PLAQUETTE_GAUGE_ERROR_ALGEBRA.md`](SU3_PLAQUETTE_GAUGE_ERROR_ALGEBRA.md) | SU(3) character-truncated one-loop error-algebra test: adjoint-only conditional qubit and fundamental Wilson-loop obstruction |
| [`su3_single_plaquette_error_algebra.py`](su3_single_plaquette_error_algebra.py), [`su3_single_plaquette_results.txt`](su3_single_plaquette_results.txt) | Exact 6/10/15-state SU(3) commutants, conditional Knill–Laflamme overlaps and illustrative eigenvalue gaps |
| [`SU3_TWO_PLAQUETTE_GAUSS_ERROR_AUDIT.md`](SU3_TWO_PLAQUETTE_GAUSS_ERROR_AUDIT.md) | Exact two-face SU(3) link-color spin networks: local Gauss-law tests, shared-link Wilson amplitudes and KL/noiseless-subsystem limits |
| [`su3_two_plaquette_spin_network.py`](su3_two_plaquette_spin_network.py), [`su3_two_plaquette_results.txt`](su3_two_plaquette_results.txt) | Reproduce 252 exact SU(3) local Gauss checks, Haar Wilson transitions, matrix-algebra tests and toy finite-volume spectrum |
| [`SU3_TWO_PLAQUETTE_19_STATE_AUDIT.md`](SU3_TWO_PLAQUETTE_19_STATE_AUDIT.md) | Nested 7→11→17→19 SU(3) loop basis; exact shared-link 1/8 and bar3/6 fusion, electric/Wilson error algebras and truncation study |
| [`su3_two_plaquette_extended.py`](su3_two_plaquette_extended.py), [`su3_two_plaquette_extended_results.txt`](su3_two_plaquette_extended_results.txt) | Exact SU(3) determinant-twisted Haar integrals, Fierz checks, rational commutants and illustrative generalized spectra |
| [`SU3_COMPLETE_MINIMAL_LINK_AUDIT.md`](SU3_COMPLETE_MINIMAL_LINK_AUDIT.md) | **Complete** SU(3) two-square minimal-link `1⊕3⊕bar3` Gauss singlets: nine total; epsilon junctions; exact Wilson leakage and exchange-error obstruction |
| [`su3_complete_minimal_links.py`](su3_complete_minimal_links.py), [`su3_complete_minimal_links_results.txt`](su3_complete_minimal_links_results.txt) | Exhaustive `3^7` irrep enumeration; 324 local Gauss checks; 9-state Haar matrices and five-coupling finite-model spectrum |
| [`SU3_PUBLISHED_SQUARE_CHAIN_BENCHMARK.md`](SU3_PUBLISHED_SQUARE_CHAIN_BENCHMARK.md) | **Exact published comparison:** Chen–Müller–Yao zero-flux two-plaquette SU(3) Hamiltonian (PRD Appendix B Table VI, Eqs. 3–5, 10–11) matches the complete nine-state model |
| [`su3_published_square_chain_benchmark.py`](su3_published_square_chain_benchmark.py), [`su3_published_square_chain_benchmark_results.txt`](su3_published_square_chain_benchmark_results.txt) | Independent qutrit-to-link bijection, all 162 directed Wilson entries, electric term and complete source-normalized 9×9 Hamiltonian identity |
| [`SU3_THREE_PLAQUETTE_PUBLISHED_CHAIN_AUDIT.md`](SU3_THREE_PLAQUETTE_PUBLISHED_CHAIN_AUDIT.md) | **Three-square SU(3)** complete minimal qutrit chain (27 states, 10 links), all nine exact published magnetic blocks, electric/magnetic error-algebra classification and finite-volume comparison |
| [`su3_three_plaquette_chain.py`](su3_three_plaquette_chain.py), [`su3_three_plaquette_results.txt`](su3_three_plaquette_results.txt) | Exact PRD Tables I/VI crosscheck, 27D Hamiltonian, 9D parent regression and dimensionless finite-size spectra/correlations |
| [`su3_three_plaquette_error_algebra.py`](su3_three_plaquette_error_algebra.py), [`su3_three_plaquette_error_algebra_results.txt`](su3_three_plaquette_error_algebra_results.txt) | Exact finite-field commutant-rank witness: Hermitian M14⊕M13, complex Wilson M27; not an active-code no-go |
| [`SU3_FOUR_FIVE_PLAQUETTE_BOUNDED_ERRORS.md`](SU3_FOUR_FIVE_PLAQUETTE_BOUNDED_ERRORS.md) | Four/five-face SU(3) finite-size gaps and charge-conjugation/reflection sectors; exact 21-element local-error-span code and Wilson counterexamples |
| [`su3_four_five_plaquette.py`](su3_four_five_plaquette.py), [`su3_four_five_plaquette_results.txt`](su3_four_five_plaquette_results.txt) | Exact source-table validation, 81/243-state Hamiltonian diagonalization and 441 exact Knill–Laflamme overlaps |

## Reproduce

To reproduce the new deterministic and comparison checks:

```bash
python3 exact_yukawa.py
python3 code_fingerprints.py            # requires numpy
python3 code_fingerprints.py --extended # adds three-generation codes
python3 local_purity.py                 # exact four-site reduced-projector purities; requires numpy
python3 published_qutrit_family.py     # reproduces external 2025 paper matrices
python3 reed_muller_family_exclusion.py # source-scoped 2026 family no-go and exact charge moments
python3 e8_comparison.py               # E8 chirality obstruction as a code obstruction; requires numpy
python3 walsh_character_lattice.py       # all-prime rigidity and masked direct-sum theorem
python3 operational_gauge_qec.py         # Gauss-only failure vs genuine lattice-gauge recovery
python3 abstract_stabilizer_recovery.py  # all single-site Pauli errors on charge-code qudits
python3 su2_matter_singlet.py           # SU(2) with matter: Gauss-invariant error fails KL (requires sympy)
python3 four_link_su2_subspace.py       # fixed-sector SU(2) DFS, failure under local exchange
python3 su2_exchange_algebra.py        # engineered SU(2) ground gap, balance conditions, full local-noise algebra
python3 five_block_su2_recovery.py     # exact one-module recovery via known five-qubit perfect code
python3 flagged_su2_measurement.py    # flag-ancilla scheduling and single-gadget fault enumeration
python3 measured_multiround_su2.py --stress  # 2 flagged rounds, adaptive 3rd; stress with prior data errors
python3 paper_style_flag_ec.py --stress    # source-inspired adaptive first-trigger recovery; t=1 vs r+s=2
python3 su2_controlled_exchange.py       # explicit SU(2)-invariant ancilla-conditioned gate Hamiltonians
python3 published_flag_figure2_catalog.py # paper Figure 2(d) original flagged error syndromes
OPENBLAS_NUM_THREADS=1 python3 two_body_su2_mediator.py # two-body SU(2) mediator, leakage and coherent gate faults
OPENBLAS_NUM_THREADS=1 python3 su2_shaped_pulse.py # optimized SU(2) smooth waveform, leakage and coupler robustness
OPENBLAS_NUM_THREADS=1 python3 su2_sk1_composite_pulse.py # SK1-style noncommuting gauge-invariant robustness sequence
OPENBLAS_NUM_THREADS=1 python3 su2_minimal_pulse_audit.py # decision on simpler pulse with all modeled calibration faults
python3 su2_error_algebra_capacity.py       # exact symmetry-error algebra (requires SymPy)
python3 su2_bridge_error_geometry.py        # exact graph-limited Knill-Laflamme conditions (requires SymPy)
python3 su2_graph_qutrit_witness.py         # explicit 3-state gauge-singlet code (requires SymPy)
OPENBLAS_NUM_THREADS=1 python3 su2_offgraph_robustness.py # 256D Hamiltonian + CPTP recovery (SymPy, NumPy, SciPy)
OPENBLAS_NUM_THREADS=1 python3 su2_native_vertex_controls.py # source-motivated four-spin vertex, coupling/echo audit (SymPy, NumPy)
OPENBLAS_NUM_THREADS=1 python3 su3_single_plaquette_error_algebra.py # SU3 character fusion/error-algebra (SymPy, NumPy); not a glueball mass
OPENBLAS_NUM_THREADS=1 python3 su3_two_plaquette_spin_network.py # two adjacent SU(3) plaquettes, explicit local Gauss constraints (SymPy, NumPy)
OPENBLAS_NUM_THREADS=1 python3 su3_two_plaquette_extended.py # exact 19-state SU3 flux-fusion audit (SymPy, NumPy, SciPy)
OPENBLAS_NUM_THREADS=1 python3 su3_complete_minimal_links.py # full 9-state minimal-link SU3 Gauss sector (SymPy, NumPy)
OPENBLAS_NUM_THREADS=1 python3 su3_published_square_chain_benchmark.py # exact source-to-code match to PRD 114 074502 Table VI (SymPy, NumPy)
OPENBLAS_NUM_THREADS=1 python3 su3_three_plaquette_chain.py # exact 27-state source-benchmarked SU3 chain (SymPy, NumPy)
OPENBLAS_NUM_THREADS=1 python3 su3_three_plaquette_error_algebra.py # F97 exact commutant certificate (SymPy)
OPENBLAS_NUM_THREADS=1 python3 su3_four_five_plaquette.py # four/five-face spectra and restricted local-error KL (SymPy, NumPy)
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
- **Minimal-control decision with complete modeled calibration families:** [`SU2_MINIMAL_CONTROL_DECISION.md`](SU2_MINIMAL_CONTROL_DECISION.md), verified by [`su2_minimal_pulse_audit.py`](su2_minimal_pulse_audit.py). Reuse the prior short single `sin²` SU(2) two-body gate without the composite's many extra local X/Y rotations. Under ±1% coherent uncertainty in either/both cross-register couplers, it has a largest sampled deviation of ~`0.0346`; the long SK1 composite's worst sampled control error is ~`0.1013` when its requisite local exchange rotations may also be wrong by ±1%. At ±5% these become ~`0.177` vs ~`0.518`; SK1 **still wins on cross-couplers alone**, but not across the full tested control set. The simple gate has shorter active duration (~38.2 vs ~344.2 plus extra rotations) and lower nominal leakage. A bounded three-shape pulse search fails to improve drift materially, while exact second-order analysis explains why merely reshaping one weak-coupling pulse cannot cancel its amplitude derivative. **Neither gate is intrinsically noise-proof, hardware ready, or evidence of the actual Standard Model's physical QEC.**
- **SK1-inspired composite exchange control:** [`SU2_SK1_COMPOSITE_PULSE.md`](SU2_SK1_COMPOSITE_PULSE.md), verified by [`su2_sk1_composite_pulse.py`](su2_sk1_composite_pulse.py). Noncommuting intra-block SU(2)-invariant exchange rotations plus six smoothly shaped cross-register exchange pulses cancel the leading coherent calibration response: a ±1% fault on either physical cross coupler changes a fixed-calibration logical gate by ~`0.00030` rather than ~`0.0087` for a single echoed primitive; at ±10%, ~`0.03` versus ~`0.087`. Independent physical 256D and symmetry-reduced 14D composite propagation agree to ~`1.1e-11`. **Major limitations:** ~`344/Delta` cross-exchange duration (~9x the prior short gate) *plus* uncounted intra-block rotations; ~`1.7e-6` nominal code leakage (larger than PR #16). **A separate exact numerical stress test shows merely ±1% common miscalibration of those extra local rotations produces ~`0.10` logical deviation (±5% produces ~`0.48–0.52`).** Thus this is NOT a robust full control sequence: local rotation calibration, physical pulse timing, bandwidth, and integrated stabilizer extraction remain unresolved. This is deterministic Hamiltonian simulation, not a viable device-level threshold or a Standard Model code.

The original five steps are addressed in Part 4 and [`OPEN_STEPS_NOTE.md`](OPEN_STEPS_NOTE.md). What's left:

- **New:** [exact Yukawa theorem and classification fingerprints](EXACT_SYMMETRY_AND_CLASSIFICATION.md), with [deterministic single-Higgs verification](exact_yukawa.py) and [code-fingerprint comparison](code_fingerprints.py). The latter distinguishes the joined codes from direct sums but does not establish novelty among all qudit codes.
- **New:** [four-qudit reduced-state purity invariant](LOCAL_PURITY_NOTE.md), computed by [local_purity.py](local_purity.py). For the two-family coupled codes there are 184 (p=3) or 249 (p=5) four-site purity hotspots, versus 150 for independent copies. The three-family counts are 327/522, versus 225. These are local-unitary invariants, not logical error rates; matching existing indecomposable qudit codes remains open.
- **Published matrix comparison:** [Prakash–Saha 2025 qutrit code family, exact reconstruction and comparison](PUBLISHED_CODE_COMPARISON.md), verified with [published_qutrit_family.py](published_qutrit_family.py). Their m=4,k=6 code has length 30 and six X checks, but encodes six logical qutrits at distance 2; the coupled matter construction encodes two at distance 3. The entire specified [[9m−k,k,2]]₃ family is excluded by parameters, **not** all odd-prime Reed–Muller/triorthogonal codes.
- **Further source-checked exclusion:** [Saha–Prakash 2026 full-evaluation punctured Reed–Muller comparison](FULL_REED_MULLER_COMPARISON.md), verified by [reed_muller_family_exclusion.py](reed_muller_family_exclusion.py). Under that paper's full-rank puncture construction, n+k must be p^m; ours have n+k=32 or 48, neither an odd-prime power at p=3 or 5. The original charge lattice is cubic-orthogonal but not pairwise self-orthogonal, unlike the unpunctured construction in that paper. This **does not exclude arbitrary Reed–Muller variants or all known qudit codes**.
- **E8 comparison:** [Lisi's E8 proposal and the Distler–Garibaldi chirality theorem in code language](E8_COMPARISON.md), verified by [e8_comparison.py](e8_comparison.py). Inside E8 ⊃ SO(10)×SO(6), the chiral four-generation block (16,4) fails the anomaly/triorthogonality test. Adding the mirror block restores it but gives distance 1, and E8's additive family labels only reproduce the distance-2 family-badge code.
- **Walsh-character theorems:** [exact single-family rigidity and the disguised direct-sum construction](WALSH_CHARACTER_THEOREMS.md), independently checked by [walsh_character_lattice.py](walsh_character_lattice.py). No additional single-family row can preserve all cubic Cartan constraints over any odd prime. All compatible five-row, pure-Walsh two-family extensions are exactly two original 15-qudit codes after a coordinate permutation—even when 44 of their 50 minimum-weight logical supports span the *original* family labels. Cross-family support counts alone are not proofs of genuine coupling; the earlier B3/B5 indecomposability proof remains essential.
- **Operational interpretation:** [physical gauge-state and Knill–Laflamme benchmark](OPERATIONAL_GAUGE_QEC.md) shows that the *formal* generation-based codes correct arbitrary single-qudit register errors, but this is not physical Standard Model protection. A three-link gauge ring fails for local gauge-invariant phase noise; an engineered 3×3 toric code succeeds only with additional magnetic checks. A physically justified SM-Hilbert-space encoding and measurable SM-compatible syndromes remain unconstructed.
- **Next non-Abelian matter result:** [two-vertex SU(2) matter singlets and exact Gauss-only KL obstruction](SU2_MATTER_GAUSS_KL.md), verified in [su2_matter_singlet.py](su2_matter_singlet.py). The full color-singlet sector of a deliberately truncated 64-state rishon/matter model has dimension two, but a link-local SU(2)-invariant electric Casimir distinguishes its logical states and produces an undetectable phase error. This is *not* a Standard Model encoding or a general SU(2) no-go theorem.
- **Restricted positive physical-code result:** [four-link SU(2) decoherence-free qubit](FOUR_LINK_SU2_DFS.md), verified by [four_link_su2_subspace.py](four_link_su2_subspace.py). In the four-link model a **proper 2D subspace** of the 16D color-singlet Hilbert space exactly protects a logical qubit from all diagonal link-end occupancy noise. Nevertheless a gauge-invariant matter–rishon color exchange acts as a logical Z phase error. This makes the protected error model explicit but leaves Hamiltonian-based protection and the full SM interpretation unproved.
- **Engineered Hamiltonian and no-go result:** [exact SU(2) exchange Hamiltonian and full-noise algebra](SU2_EXCHANGE_HAMILTONIAN_AND_NO_GO.md), verified by [su2_exchange_algebra.py](su2_exchange_algebra.py). An engineered SU(2)-invariant penalty plus collective exchange Hamiltonian makes the encoded qubit an exactly degenerate, gapped ground sector. The six local pair couplings preserve degeneracy iff three opposite-pair sums match; generic allowed exchanges split it at first order despite the gap. Complete local occupancy + conditional exchange errors generate eight full 2×2 operator blocks, precluding any 2D code correcting *that entire error algebra*. This is not a universal gauge-code no-go.
- **First exact active-recovery SU(2) example:** [five singlet modules with the established [[5,1,3]] outer quantum code](FIVE_BLOCK_SU2_RECOVERY.md), checked by [five_block_su2_recovery.py](five_block_su2_recovery.py). The verifier constructs actual local SU(2)-invariant exchange representations for module X/Z, two orthonormal five-module codewords, 16 distinct error syndromes and all 1,024 Knill–Laflamme overlaps. A one-module fault is correctable **assuming ideal cross-module syndrome measurement and recovery**. An explicit two-module fault is miscorrected. This is *concatenation of a published code*, not a new Standard Model code or a physical circuit.
- **Measured-check gadget diagnosis:** [flagged outer stabilizer extraction fault enumeration](FLAGGED_SU2_STABILIZER_NOTE.md), checked by [flagged_su2_measurement.py](flagged_su2_measurement.py). Each of four bare controlled-Pauli gadgets has 12 single-gate/ancilla-fault cases that become logical errors under naive decoding. Adding two flag-ancilla CNOTs after data couplings 1 and 3 flags all dangerous cases and yields an unambiguous recovery lookup **assuming a later perfect four-bit data syndrome**; a different flag placement has ambiguous correction classes. This is a *conditional one-gadget* result, **not** a complete one-fault-tolerant extraction round or local SU(2) hardware compilation.
- **Actual measured multi-round decoding:** [two flagged full rounds plus a measured adaptive third](MEASURED_MULTIRROUND_SU2_NOTE.md), verified by [measured_multiround_su2.py](measured_multiround_su2.py). From clean input, all 784 exhaustive one-fault cases leave no logical error (744 fully corrected, 40 at most one correctable residual). The **stronger** 11,760-case test of one preexisting module error plus one new circuit fault produced 4,760 logical alterations and 2,240 undefined lookup keys, with 344 identical-transcript-and-true-syndrome ambiguity groups. **Important correction:** this stress test has (r+s=2), exceeding the standard distance-three (t=1) *original-logical-preservation* condition, so these results do **not** disprove standard one-fault fault tolerance. It still shows those two-event histories cannot all be decoded from the same transcript to the original logical state.
- **Source-checked adaptive schedule:** [Chao–Reichardt branching strategy and exact one-fault criteria audit](ADAPTIVE_FLAG_EC_AUDIT.md), with executable [paper_style_flag_ec.py](paper_style_flag_ec.py). Reuses this repo's **own** flagged controlled-P gadget (NOT a literal transcription of the published Figure 2(c) of arXiv:1705.02329v1); at the first measured flag/nonzero syndrome it stops flagged checks and remeasures all four unflagged. For 392 exact clean-input single-circuit-fault cases: 352 corrected, 40 at most a correctable module residual, zero uncorrectable logical errors. All 15 incoming single-module Pauli errors with no circuit fault are corrected. This requires idealized encoded controls and recovery; physical SU(2) gadget compilation remains open.
- **Explicit SU(2)-gauge-invariant controlled singlet gates:** [Hamiltonian-level construction and primary-source audit](SU2_CONTROLLED_EXCHANGE_IMPLEMENTATION.md). [`su2_controlled_exchange.py`](su2_controlled_exchange.py) produces encoded ancilla-controlled X and Z from endpoint-co-located local color exchanges, verifies **12,288** exact full-kinematic-space Gauss commutators, and checks **32** Clifford Pauli conjugation identities against the original ideal gate simulator. [`published_flag_figure2_catalog.py`](published_flag_figure2_catalog.py) checks all seven distinct original Chao–Reichardt Figure 2(d) flagged syndromes. **Not solved:** elementary hardware synthesis of ancilla-conditioned two-spin exchange, physical gate faults, or literal gate-for-gate replication of Figure 2(c). The charge-matrix Standard Model program still lacks a physical encoding.
- **Two-body exchange-mediated SU(2) entangler:** [exact mathematical construction and full finite-pulse baseline](TWO_BODY_SU2_MEDIATED_GATE.md), verified by [`two_body_su2_mediator.py`](two_body_su2_mediator.py). Schur's lemma rules out a nontrivial gauge-invariant 2-body interaction of a neutral qubit with only one SU(2) fundamental. Re-encoding the ancilla into four additional color spins and using only same-vertex pairwise Heisenberg exchange plus block-spin penalties creates a **second-order logical ZZ interaction** `-J1 J2/(24 Delta)`. Full 256-state evolution measures nonzero coherent gate deviations and leakage; +/-10% gauge-invariant elementary exchange faults produce significant logical changes. **This remains a perturbative proposal, not an exact native gate, a tested device, or fault-tolerant extraction**, and it adds nontrivial charged-ancilla preparation overhead.
- **Numerical smooth-pulse improvement:** [`SU2_SMOOTH_PULSE_OPTIMIZATION.md`](SU2_SMOOTH_PULSE_OPTIMIZATION.md) with executable [`su2_shaped_pulse.py`](su2_shaped_pulse.py). Using the *same* specified SU(2)-scalar eight-spin Hamiltonian, a synchronized `Jpeak sin²(pi t/T)` exchange pulse with `Jpeak/Delta=1`, locally optimized duration `T Delta=38.234506`, and ideal individual logical-Z phase compensation gives ~`2.90e-6` operator-norm deviation and `9.77e-8` maximum product-singlet-code leakage (1024-bin verification). The 14-dimensional total SU(2) singlet reduction is cross-checked against full 256-dimensional evolution. **Crucial negative result:** a +10% error on one physical coupler still gives ~`0.0878` coherent logical-gate deviation. This is numerical pulse-control evidence, not a hardware-ready, calibration-robust CZ or actual Standard Model code.
- Compare the joined codes against **same-(n,k,p)** explicit check matrices of other published odd-prime codes, especially shortened/affine-variety constructions. That broader novelty question is open.
- Explain why nonlinear charge rows appear only for p ≡ 1 (mod 3), and whether the floor of 8 weight-3 weak spots can be broken.
- Extend the gate search to non-diagonal and permutation-assisted ("fold") transversal gates, and to the p = 3 ninth-root gates.
- Run a circuit-level simulation, with noisy checks of weight 15–16, to see whether the ε⁴ correlated errors matter in a full distillation factory.

## The whole story in one chain

Anomaly cancellation = triorthogonality → one generation is a 15-slot magic code → the second generation alone is the identical code (charges, not masses) → two together hit the twins problem → a family badge gives 2 digits but distance 2 → smart badges restore distance 3 → three families with independent labels give 45 slots, 3 digits, distance 3, coupled only at four tiles → but flavour mixing blends the families, so no family-distinguishing badge can be exact, and the digits collapse back to one → the follow-up shows the only force-compatible badges (family B−L and right-handed isospin) still need charges that neutrino mixing breaks at order one, that the joined codes are inequivalent to independent copies whose coupling shows up as correlated ε⁴ errors rather than new gates, and that distance 4 is out of reach for this whole class.

Two levels of description each refuse to supply a stable family label. At the level of symmetry (charges), the generations are perfect twins. At the level of mass, they are blended. A multi-family code needs exactly what neither level provides.
