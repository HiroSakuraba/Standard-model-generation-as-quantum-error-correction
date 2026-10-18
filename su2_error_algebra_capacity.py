"""Exact allowed-error algebra checks for four/eight SU(2) spin-halves.

Source model: two four-spin singlet qubits at one SU(2) gauge vertex,
embedded in the fourteen-dimensional global eight-spin singlet sector.

This is an operator-algebra test, NOT a distance/threshold/lifetime claim.
Require: SymPy >= 1.12. Run: python su2_error_algebra_capacity.py
"""
from itertools import product
import sympy as sp

I2 = sp.eye(2)
I4 = sp.eye(4)
Z2 = sp.diag(1, -1)


def noncrossing_matchings(indices):
    """Catalan basis: planar spin singlet pairings of an ordered register."""
    if not indices:
        yield ()
        return
    i = indices[0]
    for jpos in range(1, len(indices), 2):
        j = indices[jpos]
        for left in noncrossing_matchings(indices[1:jpos]):
            for right in noncrossing_matchings(indices[jpos + 1:]):
                yield ((i, j),) + left + right


def valence_state(n, matching):
    """Unnormalized: coefficients +-1; norm squared 2^(n/2)."""
    out = {}
    for bs in product((0, 1), repeat=n // 2):
        state = [0] * n
        for (i, j), bit in zip(matching, bs):
            state[i], state[j] = bit, 1 - bit
        out[tuple(state)] = (-1) ** sum(bs)
    return out


def dot(a, b):
    return sp.expand(sum(sp.conjugate(x) * b.get(k, 0) for k, x in a.items()))


def swap_state(state, i, j):
    out = {}
    for bits, coeff in state.items():
        moved = list(bits)
        moved[i], moved[j] = moved[j], moved[i]
        moved = tuple(moved)
        out[moved] = out.get(moved, 0) + coeff
    return out


def scale_sum(*terms):
    out = {}
    for scalar, state in terms:
        for bits, coeff in state.items():
            out[bits] = out.get(bits, 0) + scalar * coeff
    return {k: sp.simplify(v) for k, v in out.items() if sp.simplify(v) != 0}


def h(state, i, j):
    return scale_sum((sp.Rational(1, 2), swap_state(state, i, j)),
                     (-sp.Rational(1, 4), state))


def four_spin_code():
    a = {k: sp.Rational(v, 2) for k, v in valence_state(4, ((0, 1), (2, 3))).items()}
    b = {k: sp.Rational(v, 2) for k, v in valence_state(4, ((0, 2), (1, 3))).items()}
    assert dot(a, a) == dot(b, b) == 1
    assert dot(a, b) == sp.Rational(1, 2)
    q = scale_sum((-1 / sp.sqrt(3), a), (2 / sp.sqrt(3), b))
    assert dot(q, q) == 1 and dot(a, q) == 0
    return (a, q)


def tensor(a, b):
    return {x + y: u * v for x, u in a.items() for y, v in b.items()}


def projected(basis, transform):
    return sp.Matrix(len(basis), len(basis), lambda i, j:
                     sp.simplify(dot(basis[i], transform(basis[j]))))


def commutant_dimension(matrices):
    """Exact rational linear equations [X,M]=0 in column vec convention."""
    n = matrices[0].rows
    eye = sp.eye(n)
    equations = sp.Matrix.vstack(*[
        sp.kronecker_product(eye, m) - sp.kronecker_product(m.T, eye)
        for m in matrices
    ])
    # Rational row reduction via DomainMatrix, not a floating rank tolerance.
    rank = (equations.rank() if n <= 2 else
            sp.polys.matrices.DomainMatrix.from_Matrix(equations).convert_to(sp.QQ).rank())
    return n * n - rank


def singlet_exchange_representation(n):
    pairings = list(noncrossing_matchings(tuple(range(n))))
    assert len(pairings) == {4: 2, 8: 14}[n]
    states = [valence_state(n, pairing) for pairing in pairings]
    gram = sp.Matrix(len(states), len(states), lambda i, j: dot(states[i], states[j]))
    assert gram.det() != 0
    inv = gram.inv()
    def swap(i, j):
        overlaps = sp.Matrix(len(states), len(states), lambda a, b:
                             dot(states[a], swap_state(states[b], i, j)))
        rep = inv * overlaps
        assert rep * rep == sp.eye(len(states))
        return rep
    return swap


def four_spin_algebra_test():
    pair = four_spin_code()
    swap01 = projected(pair, lambda ket: swap_state(ket, 0, 1))
    swap02 = projected(pair, lambda ket: swap_state(ket, 0, 2))
    assert swap01 == -Z2
    assert swap02 == (Z2 - sp.sqrt(3) * sp.Matrix([[0, 1], [1, 0]])) / 2
    assert commutant_dimension([swap01]) == 2
    assert commutant_dimension([swap01, swap02]) == 1
    # Every two-spin swap commutes with collective SU(2), nevertheless
    # the two noncommuting swaps give the FULL logical M2 error algebra.
    print("FOUR SPINS: one swap commutant=2; two distinct swaps commutant=1; full M2")


def eight_spin_KL_test():
    pair = four_spin_code()
    code = [tensor(a, b) for a in pair for b in pair]
    assert projected(code, lambda ket: ket) == I4
    h04 = lambda ket: h(ket, 0, 4)
    h15 = lambda ket: h(ket, 1, 5)
    assert projected(code, h04) == sp.zeros(4)
    assert projected(code, h15) == sp.zeros(4)
    assert projected(code, lambda ket: h04(h04(ket))) == sp.Rational(3, 16) * I4
    assert projected(code, lambda ket: h15(h15(ket))) == sp.Rational(3, 16) * I4
    # Exact KL passes for E={I, h04}; its two-by-two Gram is scalar.
    # But the two *distinct* elementary cross interactions fail KL:
    cross = projected(code, lambda ket: h04(h15(ket)))
    local_pair = -I2 / 4 - Z2 / 2
    assert sp.simplify(cross - sp.kronecker_product(local_pair, local_pair) / 3) == sp.zeros(4)
    assert cross != cross[0, 0] * I4
    print("EIGHT SPINS: {I,h04} exact KL PASS; {I,h04,h15} KL FAIL")
    print("P h04 h15 P =", cross)


def eight_spin_algebra_test():
    swap = singlet_exchange_representation(8)
    within = [swap(i, i + 1) for i in (0, 1, 2, 4, 5, 6)]
    assert commutant_dimension(within) == 3
    # Two spin-zero four-spin registers occupy dim=2*2=4; other
    # global singlets have spin blocks j_A=j_D=1 (3*3=9), 2 (1*1=1).
    ca = sum((swap(i, j) for i in range(4) for j in range(i + 1, 4)), sp.zeros(14))
    ranks = [14 - (ca - eig * sp.eye(14)).rank() for eig in (0, 2, 6)]
    assert ranks == [4, 9, 1], ranks
    assert (ca * (ca - 2 * sp.eye(14)) * (ca - 6 * sp.eye(14))).is_zero_matrix
    # One *additional* inter-register swap makes the exchange graph connected.
    # Neighbor exchanges then generate S8, acting irreducibly on 8-spin J=0.
    # By Schur and Burnside, its complex group algebra is ALL M14(C).
    assert commutant_dimension(within + [swap(3, 4)]) == 1
    print("EIGHT-SPIN GLOBAL SINGLET: dim=14; jA=jD sectors=4+9+1")
    print("Intrablock swap algebra: M4 + M9 + C, dimension=98; commutant=3")
    print("Intra + one bridge swap: full M14, dimension=196; commutant=1")


def main():
    four_spin_algebra_test()
    eight_spin_KL_test()
    eight_spin_algebra_test()
    print("ALL EXACT SU(2) ALLOWED-ERROR ALGEBRA CHECKS PASSED")


if __name__ == "__main__":
    main()
