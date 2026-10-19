# Gauss-constrained finite-group test of the Standard Model centre echo

19 October 2026. Complements SM_CENTER_ECHO.md and GLOBAL_FORM_CODES.md. It is a **structural negative-control experiment**, not a physical Standard Model lattice regulator.

## Gauss constraints included explicitly

Take four vertices connected by four oriented periodic-cycle links. The link register carries an A-centre electric flux a_l in Z_d with d=6/|Gamma|, alongside a second B gauge flux b_l in Z2. At each vertex v, impose the full physical-subspace equations:

  a_v - a_(v-1) = 0 mod d,
  b_v - b_(v-1) = m_v mod 2,

where m_v∈{0,1} represents charged hard-core matter in the auxiliary B sector. The parity sum m_v is even. The enumerated Gauss-invariant basis has dimension 16d: 96,48,32,16 as Gamma goes through 1,Z2,Z3,Z6.

A B-charged matter hop from v to v+1 swaps a particle and a hole and flips the **specific connecting B link**. Every resulting transition is explicitly checked to satisfy Gauss at every vertex. Magnetic A Wilson dynamics move a conserved closed loop a↦a+1, leaving Gauss at every vertex. A diagonal electric Hamiltonian completes the toy.

This deliberately factorizes the *combined-centre* gauge flux A from a different matter gauge sector B. Such matter is centre neutral by construction but still gauge charged under B. **This is not a reconstruction of actual Standard Model SU(3), SU(2), U(1) matter representations, their non-Abelian fusion rules, spinors, or local Yukawa couplings.** It is considerably more constrained than the earlier free link register, but still a finite-group analogue.

## What happens to the six-pulse echo

Represent the allowed global-form electric centre classes by e6=q a mod6, q=|Gamma|, d=6/q. The proposed Z6 generator is P|a,b,m>=exp(2πi q a/6)|a,b,m>. The six toggling frames P^k leave H_matter and H_electric invariant and cancel an allowed nontrivial closed A magnetic loop for d>1.

For full Gamma=Z6 the only admitted centre electric class is e6=0, so P=I on the entire physical Hilbert space. No nontrivial A loop-shift operator remains (the projected, centre-charged Wilson loop is forbidden). The apparent echo has **no effect**, not a miraculous decoupling. This is the essential global-form negative control.

Finite echo error is compared with the evolved matter/electric target using 4,16,64 cycles. Exact average-Hamiltonian cancellation is a matrix identity within each explicitly constrained finite model. Convergence at nonzero cycle length is approximate and can be impacted by noncommuting electric and magnetic terms.

**Crucial limitation:** a real SM closed fundamental Wilson loop in a permitted representation can exist or not depending on the quotient and line spectrum; this toy models only one class shift. Neither the imposed d values nor absence of a loop in q=6 should be taken as a complete characterization of non-Abelian gauge Wilson dynamics. Tong (2017), *Line Operators in the Standard Model*, DOI https://doi.org/10.1007/JHEP07(2017)104 is the primary source for quotient-dependent line spectra.

## Next physical check

Construct gauge-covariant finite representation fusion operators for actual SU(3) triplets, SU(2) doublets, U(1) hypercharge and dynamical matter on a minimal lattice; restrict to genuine reps under each quotient **before** twirling. Include native pulse strengths, time, Gauss-violating errors and a feasible syndrome readout. Only then interpret the echo as a candidate analogue control architecture.

Reproduce:

    OPENBLAS_NUM_THREADS=1 python -m unittest -v test_sm_center_gauss_constrained
    OPENBLAS_NUM_THREADS=1 python sm_center_gauss_constrained.py
