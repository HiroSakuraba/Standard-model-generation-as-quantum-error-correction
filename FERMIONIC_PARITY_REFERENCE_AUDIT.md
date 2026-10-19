# Fermionic generation: a reference qubit is not the unpunctured stabilizer state

**19 October 2026.** Crosscheck of FERMIONIC_GENERATION_CODE.md and spinor_reed_muller.py with a new explicit 2^16-state implementation.

## Two constructions with different numbers of logical qubits

**A. Physical even-parity relational encoding (k=1).** Start with the punctured Reed–Muller [[15,1,3]] binary code. Its two basis states |0_L> and |1_L> have respectively even and odd data-mode fermion parity. A superposition α|0_L>+β|1_L> without a reference violates fermion-parity superselection.

Introduce one reference occupation mode ν and encode

  |0_R> = |0_L>|0_ν>,   |1_R> = |1_L>|1_ν>.

**Both states have even total fermion parity.** Their arbitrary superposition α|0_R>+β|1_R> remains entirely in that allowed sector. The original 15-mode stabilizer group has 14 independent checks (four X-type, ten Z-type); adding **one total-parity Z check** gives 15 independent stabilizers on 16 modes, hence **one encoded logical qubit**. A logical Z is reference mode parity, and a parity-preserving (though nonlocal) logical X complements all 16 occupation bits. The new program checks their Pauli action and superselection consistency explicitly on 65,536 amplitudes.

This extra reference mode does **not** make individual operators local or experimentally available. A logical X implemented as 16 simultaneous occupation flips is parity-even but can be extremely nonlocal in fermionic hardware; additional fault tolerance and measurement circuits would be needed.

**B. Full unpunctured RM(1,4) stabilizer state (k=0).** The extra X check with support on *all 16 modes* is the dressed logical X of A. Imposing that as an additional +1 stabilizer fixes the particular state

  |+_R> = (|0_L,0_ν> + |1_L,1_ν>)/sqrt(2).

The independent checks become 16 on 16 modes (five X and eleven Z). The code space is one-dimensional, **no arbitrary logical qubit remains**. Therefore the original statement that the unique 16-mode state 'makes the qubit physical' conflates a prepared reference state with a two-dimensional relational *encoding*.

## What is independently verified

- 15-mode logical basis states have opposite parity, and an isolated parity twirl destroys their coherent logical X expectation.
- The 16-mode relational basis states each have **even total parity**, and the even logical X interchanges them; logical Z is reference parity.
- An arbitrary encoded state α|0_R>+β|1_R> keeps its coherence under total-parity superselection; a local observer of data only does not access that coherence without the reference.
- The extra all-occupation X stabilizer isolates |+_R>, reducing k=1 to k=0.
- Both relational basis states have total fermion numbers in {0,8,16}. The diagonal operation exp(iπN/4) is identity on the **whole** relational code; the isolated data 15-mode transversal phase had different logical action. Thus the reference construction changes the gate structure.
- The previously quoted mod-8 Type-II binary stabilizer criterion is a restriction on that stabilizer family. It is **not** by itself the full Fidkowski–Kitaev interacting-phase classification, nor the Dai–Freed Z16 anomaly.

## Independent literature context

Mudassar, Schuckert and Gottesman, *Fault-Tolerant Operations in Majorana-Based Quantum Codes: Gates, Measurements, and High-Rate Constructions*, PRX Quantum **7**, 020315 (2026), DOI https://doi.org/10.1103/bglx-ssjg; arXiv:2508.09928. Its parity reference-frame gadgets and even/odd Majorana-code distinction are directly relevant. **Our 16-mode relational construction is a specific formal example, not a reproduction of their full circuit methods or a claim of novelty over their work.**

The SO(10) spinor-weight and 15-qubit Reed–Muller matching remains a separate algebraic correspondence. The model's fermions are formal modes, not Standard Model particles undergoing Yukawa dynamics. This result does not imply a usable physical error-correcting qubit without an extraction/recovery circuit.

### Reproduce

    python -m unittest -v test_fermionic_parity_reference
    python fermionic_parity_reference.py

Both commands use the 15-mode source vectors from spinor_reed_muller.kets(), but the parity and rank assertions over 16 modes are independent.
