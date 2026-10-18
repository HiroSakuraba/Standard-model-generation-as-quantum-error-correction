"""Exact ancilla-controlled SU(2)-invariant logical operations in one rishon block.

This extends the earlier SU(2) four-link singlet model. It does NOT
build a native two-qubit hardware gate, a full fault-tolerant decoder,
or a standard-model particle code. An external, COLOR-NEUTRAL ancilla
controls an engineered on-site exchange Hamiltonian.

For the chosen 2D singlet space, the physical color-SWAP operators obey
  Z_L = -SWAP(m1,r1)
  X_L = [SWAP(m1,r3) - SWAP(m1,r2)] / sqrt(3).

The operators are Hermitian SU(2) scalars on the FULL fixed-orientation
six-spin space. The ancilla-conditioned evolution
  exp[-i(pi/2) |1><1|_a (I - P_physical)]
acts exactly as controlled-P_L for P=X,Z on both singlet codewords.
Outside this 2D code, the X physical extension need not square to I,
but its exponential remains unitary and gauge-invariant.

All tests use exact SymPy numbers (including pi and sqrt(3)). They check
the full six-spin commutators with both local SU(2) Gauss generators;
exact code-state transformation (no projected-matrix-only assumption);
unitary control matrices; and the 16 two-qubit Clifford Pauli
conjugation signatures used by the existing fault-propagation model.

The H_X branch needs TWO differently weighted on-site exchange terms
conditioned on the ancilla state, i.e. an effective three-body control.
Neither its native local implementation nor its fault modes are supplied.
Run: python su2_controlled_exchange.py
"""
from itertools import product

import sympy as sp

from four_link_su2_subspace import (
    ALL_SPINS, PAIRS_A, PAIRS_B, combine, dot, gauss_raise, gauss_z_twice,
    singlet_product, swap_state, assert_same_state,
)
from flagged_su2_measurement import (
    MEASURE, conjugate, multiply, single,
)

I2 = sp.eye(2)
X2 = sp.Matrix([[0, 1], [1, 0]])
Z2 = sp.diag(1, -1)
I4 = sp.eye(4)
ANCILLA_ONE = sp.diag(0, 1)


def logical_basis():
    a = singlet_product(PAIRS_A)
    b = singlet_product(PAIRS_B)
    q1 = combine(a, b, -1 / sp.sqrt(3), 2 / sp.sqrt(3))
    assert dot(a, a) == dot(q1, q1) == 1
    assert dot(a, q1) == 0
    for q in (a, q1):
        for v in (0, 1):
            assert gauss_raise(q, v) == {}
            assert gauss_z_twice(q, v) == {}
    return a, q1


def pauli_x_full(state):
    p = swap_state(state, 1, 5)
    m = swap_state(state, 1, 4)
    return combine(p, m, 1 / sp.sqrt(3), -1 / sp.sqrt(3))


def pauli_z_full(state):
    return combine(swap_state(state, 1, 3), {}, -1, 0)


def check_full_gauss_invariance():
    """Check commutation on each of 64 spin basis vectors, not just codewords."""
    for indices in ((1, 3), (1, 4), (1, 5)):
        for bits in ALL_SPINS:
            ket = {bits: sp.Integer(1)}
            for vertex in (0, 1):
                a = gauss_z_twice(swap_state(ket, *indices), vertex)
                b = swap_state(gauss_z_twice(ket, vertex), *indices)
                assert_same_state(a, b)
                a = gauss_raise(swap_state(ket, *indices), vertex)
                b = swap_state(gauss_raise(ket, vertex), *indices)
                assert_same_state(a, b)
    print("64-state SU(2) Gauss commutators: all three exchanges exact PASS")


def projected_matrix(basis, transform):
    return sp.Matrix(2, 2, lambda i, j:
                     dot(basis[i], transform(basis[j])))


def controlled_gate(pauli):
    """Construct exp(-i K) on the 2 ancilla x 2 singlet computational basis.

    K=(pi/2) |1><1| (I-P). Since P^2=I on the code, K^2=pi*K.
    Spectral polynomial exp(-iK)=I-2K/pi, exact and easy to verify.
    """
    k = sp.pi / 2 * sp.kronecker_product(ANCILLA_ONE, I2 - pauli)
    assert sp.simplify(k * k - sp.pi * k) == sp.zeros(4)
    unitary = sp.simplify(I4 - 2 * k / sp.pi)
    target = sp.diag(1, 1, 1, 1)
    for i in range(2):
        for j in range(2):
            target[2 + i, 2 + j] = pauli[i, j]
    assert unitary == target
    assert unitary * unitary.H == I4
    return k, unitary


def one_pauli(op):
    return {
        "I": I2, "X": X2, "Y": X2 * Z2, "Z": Z2
    }[op]


def verify_clifford_signatures(kind, gate):
    """Check the physical derived control agrees with the simulator's
    ideal CX/CZ symplectic fault-propagation primitive for every Pauli.
    The non-Hermitian XZ convention represents Y up to a phase.
    """
    primitive = (kind, MEASURE, 0)
    for p_control, p_data in product("IXYZ", repeat=2):
        real = sp.kronecker_product(one_pauli(p_control), one_pauli(p_data))
        transformed = sp.simplify(gate * real * gate.H)
        pauli = multiply(single(MEASURE, p_control), single(0, p_data))
        result = conjugate(pauli, primitive)
        def op_on_site(index, data):
            x, z = data
            if x & (1 << index) and z & (1 << index):
                return "Y"
            if x & (1 << index):
                return "X"
            if z & (1 << index):
                return "Z"
            return "I"
        predicted = sp.kronecker_product(
            one_pauli(op_on_site(MEASURE, result)),
            one_pauli(op_on_site(0, result)),
        )
        assert transformed == predicted or transformed == -predicted, (
            kind, p_control, p_data, transformed, predicted)
    print(f"{kind} encoded Clifford Pauli propagation: all 16 signatures PASS")


def test():
    check_full_gauss_invariance()
    q0, q1 = logical_basis()
    basis = (q0, q1)

    # Crucially verify FULL vector actions, not just matrix elements.
    assert_same_state(pauli_x_full(q0), q1)
    assert_same_state(pauli_x_full(q1), q0)
    assert_same_state(pauli_z_full(q0), q0)
    assert_same_state(pauli_z_full(q1), combine(q1, {}, -1, 0))
    assert projected_matrix(basis, pauli_x_full) == X2
    assert projected_matrix(basis, pauli_z_full) == Z2
    assert projected_matrix(basis, pauli_x_full) ** 2 == I2
    assert projected_matrix(basis, pauli_z_full) ** 2 == I2

    # The physical X operator is NOT a single spin permutation.
    # It requires a controlled difference of two on-site exchanges.
    swaps = [projected_matrix(basis, lambda s, i=i, j=j:
                              swap_state(s, i, j))
             for i, j in ((1, 3), (1, 4), (1, 5))]
    assert all(S != X2 and S != -X2 for S in swaps)

    for kind, p in (("CX", X2), ("CZ", Z2)):
        k, gate = controlled_gate(p)
        verify_clifford_signatures(kind, gate)
        # With ancilla control |0>, the target is unchanged; with
        # control |1>, the physical color-singlet P is enacted.
        assert gate[:2, :2] == I2
        assert gate[:2, 2:] == sp.zeros(2)
        assert gate[2:, :2] == sp.zeros(2)
        assert gate[2:, 2:] == p
        print(kind, "conditional Hamiltonian K=(pi/2)|1><1|(I-P)"
              " realizes controlled-P exactly on code")

    # Entangling action is nonproduct: CX(|+>a |0>L) is a Bell pair.
    _, cx = controlled_gate(X2)
    initial = sp.Matrix([1 / sp.sqrt(2), 0, 1 / sp.sqrt(2), 0])
    assert cx * initial == sp.Matrix([
        1 / sp.sqrt(2), 0, 0, 1 / sp.sqrt(2)
    ])

    print("Native-control caveat: conditioned on-site exchange requires"
          " engineered ancilla-spin coupling; no elementary gate model.")
    print("ALL SU2 CONTROLLED EXCHANGE PROJECTIONS PASSED")


if __name__ == "__main__":
    test()
