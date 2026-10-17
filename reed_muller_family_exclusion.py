"""Source-scoped comparison with full-evaluation punctured Reed–Muller triorthogonal codes.

Saha & Prakash, arXiv:2510.10852, Secs 2–3, start from ALL p**m
evaluation positions, then puncture t independent positions. In the
nondegenerate, full-rank-after-puncture setup of Eq. (1), n+t=p**m and
the quantum code encodes k=t. So n+k MUST equal p**m.

Checks also identify an important distinction between their
quadratic+cubic triorthogonality and the matter lattice's cubic
anomaly-cancellation condition.

No claim excludes arbitrary RM subcodes, shortenings, affine-variety
evaluation codes, rank-deficient puncturing, or other constructions.
"""
from itertools import combinations, combinations_with_replacement

from charge_lattice_tools import B3, B5, Q0, W16, coupled


SOURCE = "https://arxiv.org/abs/2510.10852"


def full_evaluation_parameter_match(p, n, k):
    """Can (n,k) obey n+k=p**m for some integer m>=1?"""
    assert p >= 2 and n >= 1 and k >= 1
    total, power, m = n + k, p, 1
    while power < total:
        power *= p
        m += 1
    return m if power == total else None


def parity_moments():
    """SO(10) even-sign 16-point evaluation set over any odd F_p.

    Each one-, two- or three-character monomial other than a repeated
    square has zero total; quadratic diagonal moments are 16.
    """
    assert len(W16) == len(set(W16)) == 16
    assert W16[0] == Q0
    assert all(len(q) == 5 and set(q) <= {-1, 1}
               and q[0] * q[1] * q[2] * q[3] * q[4] == 1
               for q in W16)
    assert all(sum(q[i] for q in W16) == 0 for i in range(5))
    for i in range(5):
        for j in range(5):
            second = sum(q[i] * q[j] for q in W16)
            assert second == (16 if i == j else 0)
    for i, j, k in combinations_with_replacement(range(5), 3):
        assert sum(q[i] * q[j] * q[k] for q in W16) == 0
    print("16-point parity-evaluation moments: cubic = 0; quadratic = 16 I_5")


def matter_charge_witness(B, families, p):
    """Check cubic anomaly and a nonzero quadratic inner product.

    The FIRST differential charge row on the first family axis is
    (q_1, -q_1[, 0]) and has squared norm 2*16=32, nonzero modulo
    every odd prime. This is a witness that the *given unpunctured
    charge span* is NOT a quadratic-self-orthogonal RM triorthogonal
    space, although it is cubic-orthogonal over the integers.
    """
    code = coupled(B, families, p)
    G = code.Gint
    full_length = G.shape[1]
    assert full_length == families * 16
    assert code.n == 15 * families
    assert code.k == families
    assert len(code.punct) == families

    # First extra row: families ( +1, -1 ) or (+1,-1,0).
    first_extra = [int(x) for x in G[5]]
    expected = [q[0] for q in W16] + [-q[0] for q in W16]
    if families == 3:
        expected += [0] * 16
    assert first_extra == expected
    self_product = sum(x * x for x in first_extra)
    assert self_product == 32 and self_product % p != 0

    # Exact integer cubic identities, repeated row indices included.
    for i, j, k in combinations_with_replacement(range(len(G)), 3):
        assert sum(int(a) * int(b) * int(c)
                   for a, b, c in zip(G[i], G[j], G[k])) == 0, (
                       families, p, i, j, k)

    assert code.n + code.k == full_length
    assert full_evaluation_parameter_match(p, code.n, code.k) is None
    print(f"[[{code.n},{code.k},3]]_{p}: full 16*{families}={full_length}, "
          f"not a power of {p}; cubic identities exact, "
          f"quadratic witness 32 mod {p} = {32 % p}")


def test():
    # Positive controls taken directly from the full-evaluation
    # construction, showing the obstruction is parameter-specific.
    assert full_evaluation_parameter_match(3, 26, 1) == 3
    assert full_evaluation_parameter_match(3, 25, 2) == 3
    assert full_evaluation_parameter_match(5, 24, 1) == 2
    assert full_evaluation_parameter_match(5, 23, 2) == 2
    assert full_evaluation_parameter_match(3, 30, 2) is None
    assert full_evaluation_parameter_match(3, 45, 3) is None
    assert full_evaluation_parameter_match(5, 30, 2) is None
    assert full_evaluation_parameter_match(5, 45, 3) is None
    parity_moments()
    for p, B in ((3, B3), (5, B5)):
        for families in (2, 3):
            matter_charge_witness(B, families, p)
    print("ALL REED-MULLER FAMILY-EXCLUSION CHECKS PASSED")


if __name__ == "__main__":
    test()
