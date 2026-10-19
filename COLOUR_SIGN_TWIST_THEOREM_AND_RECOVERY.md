# Independent nine-quark sign-twist theorem, dynamics, and recovery

**19 October 2026.** Independent follow-up to COLOUR_SINGLET_CODE.md. This is an exact representation-theoretic argument supplemented by deterministic double-precision checks; not a certified computer-assisted proof of maximal code dimension in the most general case.

## Theorem: a 21-dimensional code correcting arbitrary single-colour-qutrit errors

On nine distinguishable colour-qutrit tensor factors, the SU(3)-singlet sector carries the irreducible S₉ Specht module with Young diagram λ=(3,3,3). Its hook-length dimension is **42**.

Construct the standard Young-tableaux basis of shape (3,3,3). Each adjacent swap s_i=(i,i+1) acts by the *Young orthogonal form*. For a tableau T, let c_i(T)=column(i)-row(i) and d=c_(i+1)-c_i. Then

s_i |T> = (1/d)|T> + sqrt(1-1/d²)|s_iT>,

with second term absent if the exchanged filling is not standard. The matrices obey s_i²=I, distant commutation and braid relations. This construction is independent of the expensive 3⁹ epsilon-tensor embedding from colour_singlet_code.py.

Define J on tableaux as follows: transpose the 3×3 tableau, then multiply by the sign of the permutation given by its row-major word. Since λ is self-conjugate, this map stays in the same module. Direct combinatorial verification gives JᵀJ=I, Jᵀ=−J, J²=−I, and Js_i=−s_iJ for all eight adjacent swaps, hence for **all 36 transpositions**.

Let K=ker(J−iI), complex dimension 21. Every transposition maps K into its orthogonal −i partner; therefore for all u,v∈K, ⟨u|s_ij|v⟩=0.

### From exchange compression to full Knill–Laflamme

The complete global SU(3)-singlet projector kills every nontrivial SU(3) tensor component of a local error. By two-factor Schur–Weyl duality, the invariant part of any operator supported on sites i,j is a linear combination of the identity and exchange s_ij. Its compression to K is therefore scalar. This holds for **every** operator on two sites, including products E_i†F_j of arbitrary single-site errors.

Thus the Knill–Laflamme conditions hold exactly for the family of all single-site colour-qutrit faults at any of nine positions, with **distance at least three**. A three-cycle has nonscalar action on K, so distance is **exactly three**. The earlier 19,683-state independent reduced-tensor numerical check remains a valuable separate crosscheck.

### What maximality means

For a fixed transposition s_ij acting on the 42-dimensional singlet module, J conjugates its ±1 eigenspaces, hence the signature is (21,21). If P_C s_ij P_C=0, the putative C is totally isotropic for this Hermitian form, so dim C≤21. Therefore **21 is maximal among subspaces with zero compression of every exchange**. This does **not** prove there are no larger codes with *nonzero scalar* exchange compressions. The old randomized search for dimension 22 was not a proof of global optimality.

## Constructive CPTP recovery and its limits

For a known exchange fault on one pair, let S=s_12, P_+=(I−iJ)/2, P_-=(I+iJ)/2. The completely positive channel with Kraus operators sqrt(1−p) I and sqrt(p) S has an exact recovery:

- Measure the two orthogonal J sectors.
- Apply S conditioned on the −i sector.
- Decode the +i code sector by the original isometry K†.

Its logical Kraus operators are K†P_+ and K†SP_-; the complementary space is mapped to maximally mixed logical output if necessary. Their effects sum to identity in the 42D singlet module, so this is a CPTP recovery, not a postselected fidelity. The simulation verifies entanglement fidelity 1 for all tested p∈[0,1].

**Crucial distinction:** this does **not** correct a set of *unknown pair exchanges*. Products s_ij s_kl can act nonscalarly on the code. The general arbitrary single-site colour error family *is* correctable by the KL theorem and has an explicit abstract recovery through diagonalizing its Gram matrix C_(a,b)=K†E_a†E_bK=c_(a,b) I and constructing orthogonal syndrome isometries from eigenvectors. This polar/Kraus recovery is not a low-depth hardware circuit, and only the one-known-exchange instance is numerically compiled here.

## Exchange Hamiltonians do not act as scalars dynamically

For SU(3) colour matrices T_i, the one-gluon-exchange colour operator is

T_i·T_j = (s_ij−I/3)/2.

Its first-order compression is indeed constant: P_K(T_i·T_j)P_K=−P_K/6. But s_ij maps the code out of itself. Evolving the *full* singlet module under H=(s_12−I/3)/2 sends a code vector to

e^(it/6) [cos(t/2)|ψ>−i sin(t/2)s_12|ψ>],

so the **uncorrected probability of staying in K is cos²(t/2)**. At t=π, it is zero. Confusing P_K H P_K with exact invariance under e^(−itH) would produce a false physical claim. If S is a known calibrated error mechanism and the J measurement is available, the same syndrome recovery reverses this coherent event, but that is active correction with a potentially nonlocal measurement.

## Physics boundaries

The nine colour qutrits are **distinguishable tensor factors** (for example, a simulator with addressable sites). A physical nine-quark baryon-like system must incorporate spatial wavefunctions, spin, flavour, fermionic antisymmetry, Gauss constraints and allowed interactions. The filled-shell argument is a symmetry motivation, not a demonstration of an isolated physical memory. Single-quark colour rotations are not gauge-invariant laboratory noise operations by themselves. There is no device-level syndrome extraction, threshold, or confinement-derived passive protection claim.

## Reproduce

Run from repository root with numpy and scipy installed:

    OPENBLAS_NUM_THREADS=1 python -m unittest -v test_colour_sign_twist_theorem
    OPENBLAS_NUM_THREADS=1 python colour_sign_twist_theorem.py

The new standalone tableau verifier uses only a 42×42 representation and is independent of colour_singlet_code.py's much larger 3⁹ embedding. The latter is retained for the separate full-tensor numerical validation.
