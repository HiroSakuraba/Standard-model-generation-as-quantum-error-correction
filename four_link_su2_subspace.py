"""A two-dimensional SU(2) singlet decoherence-free subspace with matter.

Four parallel quantum links each carry one spin-1/2 rishon, localized at
one endpoint, and there is one spin-1/2 matter doublet at each vertex.
This is the same deliberately finite kinematic model as
su2_matter_singlet.py, NOT the full SU(2) Hamiltonian or the Standard Model.

For each allowed endpoint pattern there are two SU(2) singlet channels.
A fixed-pattern 2D code is a decoherence-free subspace for ALL diagonal
functions of rishon endpoint occupation. However a local, gauge-invariant
SU(2) color-exchange SWAP on the four-doublet vertex acts as logical Z,
hence the code FAILS Knill--Laflamme for {I, SWAP}.

All matrix elements, singlet states, and KL certificates are exact SymPy
expressions (no floating-point / Monte Carlo / statevector simulation).
Run: python3 four_link_su2_subspace.py
"""
from itertools import product
from math import comb

import sympy as sp

from su2_matter_singlet import (
    parallel_link_singlet_dimensions, singlet_multiplicity
)

# Color-spin order: matter vertex 0, matter vertex 1,
# rishon link 0, link 1, link 2, link 3.
NSPIN = 6
ALL_SPINS = tuple(product((0, 1), repeat=NSPIN))
ORIENTATION = (0, 1, 1, 1)  # one rishon at v0 and three at v1
PAIRS_A = ((0, 2), (1, 3), (4, 5))
PAIRS_B = ((0, 2), (1, 4), (3, 5))


def multiplicity_by_pattern():
    """Representation-theoretic exact dimensions for all 16 patterns."""
    dims = {}
    for orientation in product((0, 1), repeat=4):
        left = orientation.count(0)
        right = orientation.count(1)
        dims[orientation] = (
            singlet_multiplicity(1 + left) *
            singlet_multiplicity(1 + right)
        )
    assert sum(dims.values()) == 16
    assert sum(value > 0 for value in dims.values()) == 8
    assert all((value == 2) == (orientation.count(0) % 2 == 1)
               for orientation, value in dims.items())
    assert parallel_link_singlet_dimensions(4) == (8, 8)
    return dims


def singlet_product(pairing):
    """Sparse normalized valence-bond product, amplitudes +/- 1/sqrt(8)."""
    assert len(pairing) == 3
    assert sorted(i for pair in pairing for i in pair) == list(range(NSPIN))
    result = {}
    for bits in product((0, 1), repeat=3):
        spins = [None] * NSPIN
        parity = 0
        for (a, b), bit in zip(pairing, bits):
            spins[a], spins[b] = bit, 1 - bit
            parity += bit
        result[tuple(spins)] = (-1 if parity % 2 else +1) / sp.sqrt(8)
    assert len(result) == 8
    return result


def combine(left, right, a, b):
    """Exact linear combination of sparse spin states."""
    return {
        key: sp.simplify(a * left.get(key, 0) + b * right.get(key, 0))
        for key in left.keys() | right.keys()
        if sp.simplify(a * left.get(key, 0) + b * right.get(key, 0)) != 0
    }


def dot(left, right):
    return sp.simplify(sum(
        sp.conjugate(amp) * right.get(bits, 0)
        for bits, amp in left.items()
    ))


def swap_state(state, i, j):
    out = {}
    for bits, amp in state.items():
        moved = list(bits)
        moved[i], moved[j] = moved[j], moved[i]
        key = tuple(moved)
        out[key] = sp.simplify(out.get(key, 0) + amp)
    return {key: amp for key, amp in out.items() if amp != 0}


def local_spins(vertex):
    """Spin labels transforming in the local SU(2) irrep."""
    return (vertex,) + tuple(2 + i for i, end in enumerate(ORIENTATION)
                             if end == vertex)


def gauss_raise(state, vertex):
    """Apply G_v^+ exactly, with SU(2) spin-half raising on each factor."""
    out = {}
    for bits, amp in state.items():
        for j in local_spins(vertex):
            if bits[j] == 1:
                shifted = list(bits)
                shifted[j] = 0
                shifted = tuple(shifted)
                out[shifted] = sp.simplify(out.get(shifted, 0) + amp)
    return {key: v for key, v in out.items() if v != 0}


def gauss_z_twice(state, vertex):
    out = {}
    for bits, amp in state.items():
        eigenvalue = sum(1 - 2 * bits[i] for i in local_spins(vertex))
        if eigenvalue:
            out[bits] = sp.simplify(eigenvalue * amp)
    return {key: v for key, v in out.items() if v != 0}


def projected(basis, transform):
    return sp.Matrix(len(basis), len(basis), lambda i, j:
                     dot(basis[i], transform(basis[j])))


def assert_same_state(a, b):
    assert a.keys() == b.keys()
    assert all(sp.simplify(a[k] - b[k]) == 0 for k in a)


def validate_swap_is_gauge_invariant():
    """The on-site SWAP of two v1 doublets commutes with all G_v."""
    for vertex in (0, 1):
        indices = local_spins(vertex)
        assert ((1 in indices) == (3 in indices))
        for bits in ALL_SPINS:
            ket = {bits: sp.Integer(1)}
            left_z = gauss_z_twice(swap_state(ket, 1, 3), vertex)
            right_z = swap_state(gauss_z_twice(ket, vertex), 1, 3)
            assert left_z == right_z
            left_p = gauss_raise(swap_state(ket, 1, 3), vertex)
            right_p = swap_state(gauss_raise(ket, vertex), 1, 3)
            assert left_p == right_p


def test():
    dims = multiplicity_by_pattern()
    print("4-link SU(2) singlet dimensions:", dict(dims))

    A = singlet_product(PAIRS_A)
    B = singlet_product(PAIRS_B)
    assert dot(A, A) == dot(B, B) == 1
    assert dot(A, B) == sp.Rational(1, 2)

    # Orthonormal basis of the two inequivalent spin-zero intertwiner
    # channels for this single endpoint configuration.
    # |0L>=A ; |1L> = (2 B - A)/sqrt(3).
    Q0 = A
    Q1 = combine(A, B, -1 / sp.sqrt(3), 2 / sp.sqrt(3))
    basis = (Q0, Q1)
    for i, state in enumerate(basis):
        for j, other in enumerate(basis):
            assert dot(state, other) == int(i == j)
        for vertex in (0, 1):
            assert gauss_z_twice(state, vertex) == {}
            assert gauss_raise(state, vertex) == {}

    # Every link-end occupation N_e is constant across the two
    # states. Products and arbitrary diagonal functions of all the
    # N_e therefore also act as scalar on the code.
    # This proves a *passive* decoherence-free subspace for a precise
    # noise class, not a universal recovery scheme.
    identity = sp.eye(2)
    for i, end in enumerate(ORIENTATION):
        n_left = int(end == 0)
        N = projected(basis, lambda state, val=n_left: {
            k: val * amp for k, amp in state.items()
        })
        assert N == n_left * identity
        phase = projected(basis, lambda state, val=n_left: {
            k: ((-1) ** val) * amp for k, amp in state.items()
        })
        assert phase == ((-1) ** n_left) * identity
    # All diagonal functions of 4 occupations are scalar because each
    # member has exactly the same fixed link-end occupation pattern.

    # The microscopic SU(2)-invariant SWAP(m1, rishon1) at vertex 1
    # conserves every endpoint pattern and all Gauss constraints.
    validate_swap_is_gauge_invariant()
    swap = projected(basis, lambda state: swap_state(state, 1, 3))
    assert swap == sp.diag(-1, 1)
    # Spin exchange S_m1 · S_r1 = (SWAP - 1/2 I)/2 is a local,
    # gauge-invariant two-body operator on vertex 1.
    exchange = swap / 2 - identity / 4
    assert exchange == sp.diag(sp.Rational(-3, 4), sp.Rational(1, 4))
    # SWAP is unitary on the physical space and acts as logical Z:
    # P I†SWAP P is NOT proportional to P, so KL fails for {I,SWAP}.
    assert swap * swap == identity
    assert swap != swap[0, 0] * identity
    # SWAP preserves both logical states, not merely their span.
    assert_same_state(swap_state(Q0, 1, 3), combine(Q0, {}, -1, 0))
    assert_same_state(swap_state(Q1, 1, 3), Q1)

    print("Two orthonormal fixed-orientation SU(2) singlets: PASS")
    print("All four link-end number / phase checks PEP=cP: PASS")
    print("Gauge-invariant on-site SWAP projected:", swap)
    print("Local SU(2) exchange projected:", exchange)
    print("Knill-Laflamme for {I,SWAP}: FAIL (logical Z)")
    print("ALL FOUR-LINK SU2 SUBSPACE CHECKS PASSED")


if __name__ == "__main__":
    test()
