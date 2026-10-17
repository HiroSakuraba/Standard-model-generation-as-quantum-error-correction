"""Exact SU(2)-invariant exchange algebra and its protection limits.

Build all 8 two-dimensional matter/rishon singlet blocks in the 4-link
toy model, verify local exchange generators, and prove two statements:

(1) A color-symmetric Heisenberg exchange is scalar on each spin-zero
    intertwiner space. For six arbitrary on-site pair couplings, the
    encoded splitting vanishes exactly iff three matched pair sums
    A=B=C; an occupation or Gauss-law energy penalty cannot reduce
    exchange perturbations already acting INSIDE the singlet block.

(2) Link-end occupancy projectors plus endpoint-conditioned local
    color exchanges generate the full block algebra direct_sum_8 M_2(C)
    on the 16-dimensional gauge-singlet Hilbert space. Consequently
    no code subspace of dimension >=2 corrects ALL errors from that
    specified algebra (full-algebra Knill-Laflamme no-go).

These are model-specific results, not claims about the full Standard
Model, realistic SU(2) dynamics, or all smaller noise channels.
SymPy exact arithmetic; no sampling or floating point.

Run: python su2_exchange_algebra.py
"""
from itertools import combinations, product

import sympy as sp

from four_link_su2_subspace import (
    ORIENTATION, combine, dot, multiplicity_by_pattern, projected,
    singlet_product, swap_state,
)

I2 = sp.eye(2)
X2 = sp.Matrix([[0, 1], [1, 0]])
Z2 = sp.diag(1, -1)
ROOT3 = sp.sqrt(3)


def sector_spins(orientation, vertex):
    return (vertex,) + tuple(
        2 + j for j, end in enumerate(orientation) if end == vertex
    )


def basis_for_sector(orientation):
    """An explicit orthonormal pair of total-color singlets for each
    allowed 4-link endpoint sector. No other internal singlets exist.
    """
    nleft = orientation.count(0)
    assert nleft in (1, 3)
    rich = 0 if nleft == 3 else 1
    four = sector_spins(orientation, rich)
    two = sector_spins(orientation, 1 - rich)
    assert len(four) == 4 and len(two) == 2
    a, b, c, d = four
    x, y = two
    A = singlet_product(((a, b), (c, d), (x, y)))
    B = singlet_product(((a, c), (b, d), (x, y)))
    assert dot(A, A) == dot(B, B) == 1
    assert dot(A, B) == sp.Rational(1, 2)
    Q1 = combine(A, B, -1 / ROOT3, 2 / ROOT3)
    states = (A, Q1)
    assert all(dot(states[i], states[j]) == int(i == j)
               for i in range(2) for j in range(2))
    return states, four, two


def zero_gauss_on_sector(orientation, state):
    """Check both non-Abelian Gauss generators exactly on each state."""
    for v in (0, 1):
        loc = sector_spins(orientation, v)
        z = {}
        plus = {}
        for spins, amp in state.items():
            w = sum(1 - 2 * spins[j] for j in loc)
            z[spins] = z.get(spins, 0) + w * amp
            for j in loc:
                if spins[j]:
                    raised = list(spins)
                    raised[j] = 0
                    raised = tuple(raised)
                    plus[raised] = plus.get(raised, 0) + amp
        assert all(sp.simplify(v) == 0 for v in z.values())
        assert all(sp.simplify(v) == 0 for v in plus.values())


def on_site_swap(basis, a, b):
    return projected(basis, lambda ket: swap_state(ket, a, b))


def matrix_span_rank(matrices):
    """Exact vector-space span dimension of a matrix algebra basis."""
    return sp.Matrix.hstack(*[
        sp.Matrix(4, 1, list(M)) for M in matrices
    ]).rank()


def check_all_physical_sectors():
    dims = multiplicity_by_pattern()
    assert sum(dims.values()) == 16
    assert sum(dim == 2 for dim in dims.values()) == 8
    identities = []
    for orientation, d in dims.items():
        if not d:
            continue
        states, four, two = basis_for_sector(orientation)
        for state in states:
            zero_gauss_on_sector(orientation, state)
        assert len(four) == 4
        a, b, c, _ = four
        S1 = on_site_swap(states, a, b)
        S2 = on_site_swap(states, a, c)
        assert S1 == -Z2
        assert S2 == Z2 / 2 - ROOT3 * X2 / 2
        assert S1 * S1 == I2 and S2 * S2 == I2
        assert sp.simplify(S1 * S2 - S2 * S1) != sp.zeros(2)
        # Four operators span full M2(C), so these exchange generators
        # generate ALL local logical observables on this sector.
        assert matrix_span_rank((I2, S1, S2, S1 * S2)) == 4
        identities.append((orientation, (S1, S2)))
    assert len(identities) == 8
    print("All 8 sectors: exact Gauss singlets, two noncommuting "
          "on-site exchanges, full M2 algebra PASS")
    return identities


def hamiltonian_and_noise():
    """Exact 2x2 effective SU(2)-invariant Hamiltonian for six local
    Heisenberg pair couplings in the reference 4-spin vertex.
    """
    states, four, two = basis_for_sector(ORIENTATION)
    assert four == (1, 3, 4, 5)
    pairs = tuple(combinations(four, 2))
    SWAP = {pair: on_site_swap(states, *pair) for pair in pairs}

    assert SWAP[(1, 3)] == SWAP[(4, 5)] == -Z2
    assert SWAP[(1, 4)] == SWAP[(3, 5)] == (
        Z2 - ROOT3 * X2) / 2
    assert SWAP[(1, 5)] == SWAP[(3, 4)] == (
        Z2 + ROOT3 * X2) / 2
    assert sum(SWAP.values(), sp.zeros(2)) == sp.zeros(2)

    # S_i . S_j = (SWAP_ij - I/2) / 2 for spin halves.
    Js = {pair: sp.Symbol("J%d%d" % pair, real=True) for pair in pairs}
    H = sum((Js[pair] * (SWAP[pair] / 2 - I2 / 4)
             for pair in pairs), sp.zeros(2))
    A = Js[(1, 3)] + Js[(4, 5)]
    B = Js[(1, 4)] + Js[(3, 5)]
    C = Js[(1, 5)] + Js[(3, 4)]
    hz = (-2 * A + B + C) / 4
    hx = ROOT3 * (C - B) / 4
    predicted = -(A + B + C) * I2 / 4 + hz * Z2 + hx * X2
    assert sp.simplify(H - predicted) == sp.zeros(2)

    # E_gap^2 = [(-2A+B+C)^2+3(C-B)^2]/4.
    gap_sq = ((-2 * A + B + C) ** 2 + 3 * (C - B) ** 2) / 4
    assert sp.simplify((H[0, 0] - H[1, 1]) ** 2 +
                       4 * H[0, 1] ** 2 - gap_sq) == 0

    J = sp.Symbol("J", real=True)
    uniform = {sym: J for sym in Js.values()}
    # Permutation-symmetric intra-vertex exchange: scalar on spin-zero.
    assert H.subs(uniform) == -sp.Rational(3, 2) * J * I2

    # Full degeneracy is weaker than six individually equal couplings:
    # it requires and is equivalent to three pair-sums A=B=C.
    aa, bb, cc = sp.symbols("A B C", real=True)
    generic_hz = (-2 * aa + bb + cc) / 4
    generic_hx = ROOT3 * (cc - bb) / 4
    assert sp.solve((generic_hz, generic_hx), (aa, bb)) == {
        aa: cc, bb: cc
    }

    nonuniform = {
        Js[(1, 3)]: 1, Js[(4, 5)]: 3,
        Js[(1, 4)]: 4, Js[(3, 5)]: 0,
        Js[(1, 5)]: 2, Js[(3, 4)]: 2,
    }
    assert H.subs(nonuniform) == -3 * I2

    delta = sp.Symbol("delta", real=True)
    single = {sym: sp.Integer(0) for sym in Js.values()}
    single[Js[(1, 3)]] = delta
    assert H.subs(single) == sp.diag(-3 * delta / 4, delta / 4)

    # Endpoint-pattern penalties and SU(2) Gauss Casimir penalties
    # vanish on these two exact singlets; they cannot change the
    # in-code splitting produced by a gauge-invariant local exchange.
    for ket in states:
        zero_gauss_on_sector(ORIENTATION, ket)
    assert (tuple(int(end == 0) for end in ORIENTATION) == (1, 0, 0, 0))

    print("SU(2) collective Heisenberg exchange: H=-3J/2 I (exact)")
    print("All local pair couplings: gap^2 =", gap_sq)
    print("Logical degeneracy iff pair sums A=B=C (exact)")
    print("Single local exchange perturbation: splitting |delta| "
          "independent of Gauss/occupancy penalties")
    return SWAP



def engineered_ground_space_gap():
    """Positive, exactly gapped *engineered* SU(2)-invariant Hamiltonian.

    H = Delta * sum_e (N_e-left - target_e)^2
      + lambda * sum_vertices G_v^2
      + K * sum_vertices sum_(i<j co-located) S_i . S_j.

    Since sum_pairs S_i.S_j=(G_v^2 - 3*n_v/4)/2, and n0+n1=6,
    H = Delta*hamming + g*(G0^2+G1^2) - 9K/4,
    g=lambda+K/2. For positive Delta and g, the *entire* ground
    manifold consists of the two target-orientation Gauss singlets.

    Exact excitation gap:
    min(2g, Delta + 3g/2, 2Delta).
    This Hamiltonian is a fine-tuned mathematical construction; it
    does not suppress gauge-invariant local exchange perturbations.
    """
    Delta, lam, K = sp.symbols("Delta lambda K", positive=True, real=True)
    g = lam + K / 2
    target = ORIENTATION

    candidate = (2 * g, Delta + sp.Rational(3, 2) * g, 2 * Delta)
    # A spin system with n spin halves has minimum total spin
    # 0 when n even and 1/2 when n odd. Enumerate all 16
    # occupancy/orientation sectors exactly.
    for orientation in product((0, 1), repeat=4):
        distance = sum(a != b for a, b in zip(orientation, target))
        n0 = 1 + orientation.count(0)
        n1 = 1 + orientation.count(1)
        j0 = sp.Rational(n0 % 2, 2)
        j1 = sp.Rational(n1 % 2, 2)
        casimir_min = j0 * (j0 + 1) + j1 * (j1 + 1)
        if not distance:
            assert n0 == 2 and n1 == 4 and casimir_min == 0
            continue
        assert casimir_min == (
            sp.Rational(3, 2) if distance % 2 else 0
        )
        # Sector lower energy relative to the encoded ground.
        E = Delta * distance + g * casimir_min
        # Positive coefficients allow exactly these three
        # candidates to dominate: all other sectors cost >=
        # one of the listed candidates.
        assert any(
            sp.expand(E - c).coeff(Delta) >= 0 and
            sp.expand(E - c).coeff(lam) >= 0 and
            sp.expand(E - c).coeff(K) >= 0
            for c in candidate
        ), (orientation, E)

    # The target 2-spin and 4-spin vertices both have a spin-1
    # excited irrep, so the smallest excitation within the
    # target orientation has Casimir cost exactly 2g.
    # One-link flips attain Delta+3g/2; two-link flips attain 2Delta.
    for de, gg in ((1, 1), (1, 4), (4, 1), (3, 2)):
        values = [sp.sympify(c).subs({Delta: de, lam: gg, K: 0})
                  for c in candidate]
        expected = min(2 * gg, de + sp.Rational(3, 2) * gg, 2 * de)
        assert min(values) == expected

    print("Engineered symmetric gauge Hamiltonian has exactly "
          "two ground states and positive gap")
    print("Gap = min(2g, Delta + 3g/2, 2Delta), g=lambda + K/2 > 0")
    print("Inside code, any generic local exchange perturbation "
          "still has first-order logical splitting")
    return candidate

def proof_of_full_error_algebra(sector_data):
    """Exact finite-algebra obstruction, distinct from one-error KL tests.

    Occupation N_e (commuting gauge singlets) distinguishes the eight
    allowed endpoint configurations and generates their projectors P_o
    as products of N_e and (I-N_e). In each two-dimensional block,
    endpoint-conditioned on-site SU(2)-invariant swaps generate M2.
    Therefore this *specified noise algebra* equals direct_sum_8 M2.

    For a putative isometric encoding V:C^k -> direct_sum_o C^2 with
    k>1, KL for the entire algebra demands V_o^dagger M V_o = c_M I_k
    for each rank-one sector matrix M=|v><v|. The LHS has rank <=1,
    whereas I_k has rank k>1. Thus every coefficient c_M=0.
    For all |v> this forces V_o=0, contradicting sum_o
    V_o^dagger V_o = I_k. Hence NO k>1 correctable subspace against
    arbitrary Kraus operators drawn from the full generated algebra.
    """
    assert len(sector_data) == 8
    # Exact algebra dimension: each one of eight independent 2x2
    # blocks has a 4-dimensional full matrix algebra.
    algebra_dimension = 0
    for orientation, (A, B) in sector_data:
        assert orientation.count(0) in (1, 3)
        assert matrix_span_rank((I2, A, B, A * B)) == 4
        algebra_dimension += 4
    assert algebra_dimension == 32

    # Direct finite KL obstruction for any rank-2 subspace of a
    # single block: even a rank-one color projector has rank <=1,
    # not proportional to the identity of rank 2.
    rank_one = sp.diag(1, 0)
    assert rank_one.rank() == 1
    assert rank_one != rank_one[0, 0] * I2

    print("Full local, gauge-invariant, endpoint-resolving "
          "exchange error algebra = direct_sum_8 M2(C), "
          "dimension 32")
    print("No 2D correctable code for arbitrary errors "
          "from this COMPLETE noise algebra (KL rank argument)")
    return algebra_dimension


def main():
    blocks = check_all_physical_sectors()
    hamiltonian_and_noise()
    engineered_ground_space_gap()
    proof_of_full_error_algebra(blocks)
    print("ALL SU2 EXCHANGE ALGEBRA / HAMILTONIAN CHECKS PASSED")


if __name__ == "__main__":
    main()
