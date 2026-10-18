# SU(3) two-plaquette spin networks: local Gauss law and error-algebra test

**18 October 2026. Exact finite-state algebra calculations plus illustrative double-precision diagonalization.** Executable: [`su3_two_plaquette_spin_network.py`](su3_two_plaquette_spin_network.py); recorded output: [`su3_two_plaquette_results.txt`](su3_two_plaquette_results.txt). Requires Python 3.9+, SymPy and NumPy. **Not a glueball mass, Yang–Mills mass-gap proof, quantum memory implementation, or reproduction of a published hardware simulator.**

## Why the construction differs from the earlier single-loop experiment

The earlier `su3_single_plaquette_error_algebra.py` modeled one plaquette holonomy directly as a finite collection of SU(3) characters. It could not test independently acting gauge transformations at different vertices or interactions between spatially neighboring plaquettes. Here we start with a graph having **six vertices, seven links and two adjacent square faces**, linked by one shared edge:

```text
        0 -------- 1 -------- 2
        |  LEFT    |  RIGHT   |
        |          |          |
        3 -------- 4 -------- 5
```

Oriented links indexed 0–6 are `(0,1), (1,2), (3,4), (4,5), (0,3), (1,4), (2,5)`; link 5 is shared by both faces. The left face follows edges `(0,+),(5,+),(2,-),(4,-)`; the right follows `(1,+),(6,+),(3,-),(5,-)`; their outer perimeter follows `(0,+),(1,+),(6,+),(3,-),(2,-),(4,-)`.

Every link carries a **restricted SU(3) Peter–Weyl representation**: trivial `1`, fundamental `3` (left/source index transforms as `3`, right/target as `anti-3`), or antifundamental `anti-3` (opposite endpoint representations). A colored link in `3` or `anti-3` has two color indices from 0,1,2; the kinematic link space of these irreps has dimension `1+9+9=19`. We **do not** tensor together and diagonalize all `19^7` states. We retain just seven explicitly written gauge-singlet spin-network functions, and apply a hard projection to this seven-dimensional variational subspace.

## Exactly gauge-invariant states and verified Gauss operators

Basis order: `vac`, `L+`, `L-`, `R+`, `R-`, `O+`, `O-`: empty vacuum and clockwise/counterclockwise oriented fundamental color flux around left, right and outer contours. Each loop with length `ell` has an unnormalized wavefunction obtained by summing over all `3^ell` assignments of endpoint colors subject to the invariant `3 ⊗ anti-3 → 1` Kronecker contraction at each loop vertex. The script builds these **explicit physical link-color coefficients** as sparse dictionaries: 1, 81, 81, 81, 81, 729, 729 nonzero amplitudes, respectively. Different flux patterns are orthogonal; normalizing by `3^(-ell/2)` gives unit states. An independent Haar-character Gram test gives exactly `I_7`.

At each of the six vertices the SU(3) charge generator is the sum of link-end generators: `T` on a fundamental endpoint, `-T^transpose` on an antifundamental endpoint. Acting on **every sparse color configuration** (rather than just asserting a Wilson-loop trace is invariant), the verifier checks the four root generators `E01,E10,E12,E21` and two Cartan generators `H01,H12`. Because these generate the full complexified Lie algebra `sl3`, they suffice to check the eight real `su3` generators. All **7 × 6 × 6 = 252** complete generator–state actions vanish with exact integer coefficients. This is a **direct local Gauss test**, not an argument that a single global color charge happens to vanish.

## Exact plaquette multiplication includes a shared-link interaction

Gauge fix a maximal tree for analytic matrix elements. The two independent face holonomies are `U,V ∈ SU(3)`, and the retained gauge-invariant functions are

`1, χ3(U), χbar3(U), χ3(V), χbar3(V), χ3(UV), χbar3(UV)`.

Define **fundamental Wilson-loop multiplication** by `χ3(U)` or `χ3(V)` and project to these seven states. The code explicitly integrates *every matrix element* over `SU(3) × SU(3)` using exact Haar moments (`∫U_ij U*_{kl}=δ_ik δ_jl/3`, `∫UUU=εε/6`), rather than setting magnetic transitions by hand.

For example, the exact left multiplication operator `W_L` has only these nonzero actions within the retained basis:

```
W_L |vac> = |L+>
W_L |L+> = |L->
W_L |L-> = |vac>
W_L |R+> = (1/3)|O+>
W_L |O-> = (1/3)|R->
```

The right operator `W_R` is the analogous left-right swap of the labels. These include the **shared-link fusion amplitude 1/3** between an individual face loop and the outer perimeter. In particular, two face insertions generate outer flux with amplitude `⟨O+|W_R W_L|vac⟩ = ⟨O+|W_L W_R|vac⟩ = 1/3`. This is a spatial two-plaquette effect unavailable to an isolated single-loop basis. SU(3) three-fundamental `ε` contraction also allows `W_L|L+>=|L->`, unlike simple U(1) loop raising.

**Projection warning:** Full unprojected Wilson multiplication operators commute as functions of the link holonomies. The **hard projected** `7 × 7` matrices `P W_L P` and `P W_R P` do not commute on all retained states, because discarded intermediate flux sectors matter. This is a truncation artifact, **not** evidence of intrinsically noncommuting physical magnetic Wilson loops. The calculation computes *projected single-insertion operators*, and their products must not be confused with `P W_L W_R P` in general.

## Local electric energy and quantum error correction

Local gauge-invariant electric Casimir `E_e²` is `4/3` on an occupied fundamental or antifundamental link and zero on vacuum. On our seven-state spin-network basis, the total electric operator is exactly

`D = diag(0, 16/3, 16/3, 16/3, 16/3, 8, 8)`.

Each local link Casimir is diagonal on these flux states. The shared link has occupation 1 for `L±` and `R±`, but 0 for `O±`. A left-only edge separates left from right plaquette excitations. Thus **the physical geometry is resolvable by local electric observables**.

For the particular two-dimensional code `|0_L>=|L+>` and `|1_L>=|L->`, *every local electric Casimir and arbitrary products of them* have identical action on the two orientations. Their entire diagonal electric-error algebra obeys Knill–Laflamme on this code; in this restricted noise model it is a decoherence-free logical pair. Yet the Hermitian **physically gauge-invariant** magnetic plaquette operator `B_L = W_L+W_L†` compresses to

```
             | 0  1 |
P B_L P =    | 1  0 | = logical X.
```

Adding this valid local magnetic error immediately destroys that proposed code's electric-only protection. The local Gauss generators still annihilate both codewords and all states evolved by the projected magnetic dynamics.

For the **complete generated operator algebras on this seven-state host**:

| Allowed generators | Complex algebra | Commutant | Noiseless logical qubit |
|---|---|---:|---|
| `D, B_L, B_R` | `M4(C) ⊕ M3(C)` (dimension **25**) | **2** | No |
| `D, B_L, B_R` plus all seven individual link Casimirs | same | **2** | No |
| `D, W_L, W_R, W_L†, W_R†` | full `M7(C)` (dimension **49**) | **1** | No |

The hermitian-generator algebra splits according to charge conjugation, with exact even/odd sector dimensions 4 and 3. The two sectors are **inequivalent full matrix blocks**, not equivalent noisy representations with a protected multiplicity. Exact rational commutator ranks and matrix-span closure certify these statements; no floating rank tolerance. **This does not rule out active correction of a bounded error set or a different code in a larger gauge-theory Hilbert space.**

## Illustrative two-face spectrum, not a glueball mass

For an arbitrarily normalized illustrative model `H = D − (B_L+B_R)/2`, seven-level numerical diagonalization produces

```
E ≈ [-0.19898355, 4.83333333, 5.01416210,
      5.80798878, 5.83333333, 8.01815479, 8.02534456]
first gap ≈ 5.03231689 (dimensionless)
```

The first two low-energy eigenstates **do not** correct electric-energy errors: the traceless spectral norm of the electric operator compressed to their span is approximately **2.56498012**. This is a check against the tempting claim that an energetic excitation gap automatically creates an error-correcting logical pair.

The energies are finite-volume, severe spin-network-truncation eigenvalues in arbitrary units. There is no continuum/lattice-spacing extrapolation, no finite-volume scaling, no spatial spin/parity glueball classification, no calibrated dynamics, and no full Kogut–Susskind link truncation. **None of these energies is a predicted glueball mass or a Yang–Mills mass-gap proof.**

## How this extends the research program

The operational finding is now present in a graph with independently verified local SU(3) Gauss transformations: **gauge invariance constrains what states and interactions are physical; it does not by itself prevent gauge-invariant local Wilson loops from acting as logical errors.** In this model both a positive electric-only restricted code and a physically allowed magnetic counterexample are certified.

The next work should **enlarge the spin-network basis**, starting with two-loop occupancies and shared-edge fusions into `8` and `6` representations. This is required to lessen projection artifacts and recover the real Kogut–Susskind local plaquette algebra. Then compare the full minimally link-truncated two-plaquette Hamiltonian with the spin-network approximation at increasing cutoffs, introduce physical noise and recovery operations, and only afterwards compute glueball interpolator correlation functions on larger lattices with proper lattice rotation, parity, and charge conjugation quantum numbers. Separate error-set Knill–Laflamme tests from noiseless-subsystem claims at every step.

## Reproduce

```bash
python -m pip install sympy numpy
OPENBLAS_NUM_THREADS=1 python su3_two_plaquette_spin_network.py
python -m py_compile su3_two_plaquette_spin_network.py
```

The exact certificates are the Gauss residuals, Haar matrix elements/Gram, local electric eigenvalues, commutant and algebra dimensions, and the gauge-invariant logical-X counterexample. The spectrum and low-eigenstate error diagnostic are NumPy double-precision calculations; no claimed interval certification. Reference result is saved in `su3_two_plaquette_results.txt`.

**Primary sources:** Kogut and Susskind, *Hamiltonian formulation of Wilson's lattice gauge theories*, Phys. Rev. D **11**, 395 (1975), https://doi.org/10.1103/PhysRevD.11.395; Chen, Müller and Yao, *Minimally truncated SU(3) lattice gauge theory and string tension*, Phys. Rev. D **114**, 074502 (2026), https://doi.org/10.1103/1qt7-3ywl, arXiv:2601.10065. Their actual minimally link-truncated square-plaquette-chain Hamiltonian is **not identical** to this seven-spin-network projection. Background on gauge-invariant plaquette variables: Anciavarella and Bauer, Phys. Rev. Lett. **133**, 111901 (2024), https://doi.org/10.1103/PhysRevLett.133.111901.
