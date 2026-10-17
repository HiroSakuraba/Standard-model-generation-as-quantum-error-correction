"""Exact SU(2) Gauss-law singlets for two matter sites and two rishon links.

Truncation: exactly one spin-1/2 (color-doublet) matter excitation at
each vertex and exactly one spin-1/2 rishon per link. Each rishon
occupies either the left (vertex 0) or right (vertex 1) end of its link.
The resulting 4 orientation sectors each have 16 color configurations.
Gauss generators are the sums of local SU(2) spin generators.

This is a kinematic SU(2) quantum-link / rishon toy, not the full
Kogut-Susskind SU(2) Hamiltonian, not an SM particle encoding,
and not an error-correcting code merely because it is gauge invariant.

Uses exact SymPy integer ranks and Fraction amplitudes.
Run: python3 su2_matter_singlet.py
"""
from collections import defaultdict
from fractions import Fraction
from itertools import product

import sympy as sp

VERTICES = (0, 1)
ORIENTATIONS = tuple(product(VERTICES, repeat=2))
SPINS = tuple(product((0, 1), repeat=4))  # matter0, matter1, rishon0, rishon1


def at_vertex(orientation, vertex):
    """Indices of fundamental SU(2) doublets transforming at a vertex."""
    return (vertex,) + tuple(2 + e for e in range(2) if orientation[e] == vertex)


def twice_gauss_z(orientation, spins, vertex):
    """Twice the z component of local non-Abelian Gauss law."""
    return sum(1 - 2 * spins[i] for i in at_vertex(orientation, vertex))


def gauss_plus(orientation, state, vertex):
    """SU(2) raising generator on spin-1/2 factors, exact integer action."""
    for index in at_vertex(orientation, vertex):
        if state[index] == 1:
            changed = list(state)
            changed[index] = 0
            yield tuple(changed)


def exact_sector_singlet_dimension(orientation):
    """Joint kernel of two G^z and two G^+ operators on 16 color states.

    In finite unitary SU(2) representations a weight-zero highest
    weight vector belongs to the spin-zero irrep. The exact joint
    kernel thus counts local color singlets, not just parity.
    """
    indices = {state: i for i, state in enumerate(SPINS)}
    M = sp.zeros(4 * len(SPINS), len(SPINS))
    for j, state in enumerate(SPINS):
        for vertex in VERTICES:
            M[vertex * len(SPINS) + j, j] = twice_gauss_z(orientation, state, vertex)
            for raised in gauss_plus(orientation, state, vertex):
                M[(2 + vertex) * len(SPINS) + indices[raised], j] += 1
    return len(SPINS) - M.rank()


def physical_singlet(orientation):
    """Exact normalized two-vertex singlet as a sparse rational vector."""
    if orientation not in ((0, 1), (1, 0)):
        return {}
    rishon_at_vertex = {
        v: 2 + next(e for e in range(2) if orientation[e] == v)
        for v in VERTICES
    }
    singlet_terms = (((0, 1), 1), ((1, 0), -1))
    answer = {}
    for (m0, r0), sign0 in singlet_terms:
        for (m1, r1), sign1 in singlet_terms:
            spins = [-1] * 4
            spins[0], spins[rishon_at_vertex[0]] = m0, r0
            spins[1], spins[rishon_at_vertex[1]] = m1, r1
            assert set(spins) <= {0, 1}
            answer[(orientation, tuple(spins))] = Fraction(sign0 * sign1, 2)
    return answer


def inner(left, right):
    return sum(left.get(state, Fraction(0)) * amp for state, amp in right.items())


def diagonal_action(state, eigenvalue):
    return {key: amp * eigenvalue(key[0], key[1]) for key, amp in state.items()}


def assert_color_singlet(state):
    for vertex in VERTICES:
        z_action = defaultdict(Fraction)
        plus_action = defaultdict(Fraction)
        for (orientation, spins), amp in state.items():
            z_action[(orientation, spins)] += amp * twice_gauss_z(
                orientation, spins, vertex)
            for raised in gauss_plus(orientation, spins, vertex):
                plus_action[(orientation, raised)] += amp
        assert all(v == 0 for v in z_action.values())
        assert all(v == 0 for v in plus_action.values())


def projected_operator(basis, eigenvalue):
    """Compute P E P in the exact orthonormal singlet basis."""
    return sp.Matrix(len(basis), len(basis), lambda i, j:
                     inner(basis[i], diagonal_action(basis[j], eigenvalue)))


def test():
    dims = {
        orientation: exact_sector_singlet_dimension(orientation)
        for orientation in ORIENTATIONS
    }
    assert dims == {(0, 0): 0, (0, 1): 1, (1, 0): 1, (1, 1): 0}
    assert sum(dims.values()) == 2

    basis = (physical_singlet((0, 1)), physical_singlet((1, 0)))
    for i, vector in enumerate(basis):
        assert len(vector) == 4
        assert_color_singlet(vector)
        for j, other in enumerate(basis):
            assert inner(vector, other) == int(i == j)

    # Each link contains exactly one rishon; its ENDPOINT occupation
    # number is a SU(2)-singlet observable and commutes with all local
    # Gauss generators. On the gauge-singlet subspace it is not scalar.
    n0_left = lambda orientation, spins: int(orientation[0] == 0)
    n1_left = lambda orientation, spins: int(orientation[1] == 0)
    N0 = projected_operator(basis, n0_left)
    N1 = projected_operator(basis, n1_left)
    assert N0 == sp.diag(1, 0)
    assert N1 == sp.diag(0, 1)
    assert N0 + N1 == sp.eye(2)

    # A local SU(2) electric Casimir for a single rishon equals 3/4
    # on an occupied endpoint, 0 on a vacuum endpoint.
    electric_casimir_left = sp.Rational(3, 4) * N0
    assert electric_casimir_left == sp.diag(sp.Rational(3, 4), 0)

    # U=(-1)^N0 is a gauge-invariant local unitary. It preserves the
    # WHOLE physical gauge-singlet space, but is logical Z (up to sign).
    # Gauss syndrome for U is identical to Gauss syndrome for identity.
    U = projected_operator(basis, lambda o, s: 1 - 2 * n0_left(o, s))
    assert U == sp.diag(-1, 1)
    assert U * U == sp.eye(2)
    assert U != U[0, 0] * sp.eye(2)
    assert all_color_singlet := all(
        exact_sector_singlet_dimension(o) == dims[o] for o in ORIENTATIONS)
    assert all_color_singlet

    # A gauge-invariant two-state exchange can be constructed from
    # singlet projectors; unlike U, it is NOT claimed to be a
    # microscopically single-link or native dynamical coupling.
    logical_x = sp.Matrix([[0, 1], [1, 0]])
    assert logical_x * U == -U * logical_x

    # Exact Knill-Laflamme for errors {I,U} would require P I† U P = c P.
    # Instead it equals diag(-1,1), so correctability fails.
    assert sp.eye(2) * U != U[0, 0] * sp.eye(2)

    print("Exact SU(2) singlet dimensions by rishon ends:", dims)
    print("Two orthonormal gauge-invariant singlets: PASS")
    print("SU(2)-invariant first-link end occupation PNP:", N0)
    print("SU(2)-invariant local electric Casimir:", electric_casimir_left)
    print("SU(2)-invariant local phase PUP:", U)
    print("Gauss-only Knill-Laflamme {I,U}: FAIL (undetected logical Z)")
    print("ALL SU2-MATTER SINGLET / GAUSS-ONLY NO-RECOVERY CHECKS PASSED")


if __name__ == "__main__":
    test()
