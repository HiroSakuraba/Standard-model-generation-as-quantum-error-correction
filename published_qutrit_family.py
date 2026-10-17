"""Exact comparison against Prakash-Saha 2025, equations (8)-(11).

Independent reference construction; requires only Python 3 stdlib plus
the repository's charge_lattice_tools.py for the final cross-comparison.
Does not claim a published same-(n,k,p) indecomposable comparator.
"""
from itertools import product, combinations
from collections import Counter

PAPER = 'https://doi.org/10.22331/q-2025-06-12-1768'

# Eq. (11): top 4 logical rows H1, bottom 2 stabilizer rows H0.
PAPER_EQ11 = tuple(tuple(int(c) for c in line) for line in (
    '11000000000222',
    '00110000000222',
    '00001100000222',
    '00000011000222',
    '12121212012012',
    '00000000111222',
))


def rank(rows, p):
    if not rows:
        return 0
    A = [list(map(lambda v: int(v) % p, r)) for r in rows]
    n = len(A[0])
    d = 0
    for col in range(n):
        piv = next((i for i in range(d, len(A)) if A[i][col]), None)
        if piv is None:
            continue
        A[d], A[piv] = A[piv], A[d]
        inv = pow(A[d][col], -1, p)
        A[d] = [(x * inv) % p for x in A[d]]
        for i in range(len(A)):
            if i != d and A[i][col]:
                factor = A[i][col]
                A[i] = [(x - factor * y) % p for x, y in zip(A[i], A[d])]
        d += 1
        if d == len(A):
            break
    return d


def construction(m, k):
    """Prakash-Saha 2025, Eq. (8)-(9), first-coordinate punctures."""
    assert isinstance(m, int) and isinstance(k, int) and m >= 1 and 1 <= k <= 3 * m - 2
    n = 9 * m
    w = tuple(j % 3 for j in range(n))
    v = []
    for a in range(3 * m - 1):
        x = [0] * n
        x[3 * a:3 * a + 3] = [1, 1, 1]
        x[-3:] = [2, 2, 2]
        v.append(tuple(x))
    rows = v + [w]
    removed = set(range(0, 3 * k, 3))
    physical = [j for j in range(n) if j not in removed]
    H1 = tuple(tuple(row[j] for j in physical) for row in rows[:k])
    H0 = tuple(tuple(row[j] for j in physical) for row in [w] + rows[k:-1])
    return rows, H1, H0


def cubic_valid(full):
    """Exact all-index ternary triorthogonality, repetitions included."""
    for i in range(len(full)):
        for j in range(i, len(full)):
            for k in range(j, len(full)):
                if sum(a * b * c for a, b, c in zip(full[i], full[j], full[k])) % 3:
                    return False, (i, j, k)
    return True, None


def z_distance_up_to_two(H1, H0):
    H = H1 + H0
    n = len(H[0])

    def syndrome(rows, positions, coefficients):
        return tuple(sum(a * row[j] for j, a in zip(positions, coefficients)) % 3
                     for row in rows)

    for weight in (1, 2):
        for positions in combinations(range(n), weight):
            for coeff in product((1, 2), repeat=weight):
                if not any(syndrome(H0, positions, coeff)) and any(syndrome(H, positions, coeff)):
                    return weight, (positions, coeff)
    return None, None


def has_single_qudit_x_logical(H1, H0):
    H = H1 + H0
    r = rank(H, 3)
    for j in range(len(H[0])):
        ei = tuple(int(i == j) for i in range(len(H[0])))
        if rank(H + (ei,), 3) == r and rank(H0 + (ei,), 3) > rank(H0, 3):
            return True, j
    return False, None


def x_weight_enumerator(rows, p=3):
    """Complete X row-span enumerator, invariant under CSS-preserving monomial maps."""
    counts = Counter()
    for coeffs in product(range(p), repeat=len(rows)):
        weight = sum(sum(a * row[j] for a, row in zip(coeffs, rows)) % p != 0
                     for j in range(len(rows[0])))
        counts[weight] += 1
    assert sum(counts.values()) == p ** len(rows)
    return dict(sorted(counts.items()))


def compare_30_site_codes():
    from charge_lattice_tools import B3, coupled

    target = coupled(B3, 2, 3)
    _, H1, H0 = construction(4, 6)
    assert target.n == 30 and target.k == 2
    assert rank(H1 + H0, 3) - rank(H0, 3) == 6
    assert rank(target.SX.tolist(), 3) == 6
    external = x_weight_enumerator(H0)
    internal = x_weight_enumerator(target.SX.tolist())
    assert external == {0: 1, 6: 30, 9: 40, 12: 90, 15: 60, 18: 22, 24: 486}, external
    assert internal == {0: 1, 10: 4, 11: 4, 12: 8, 16: 52, 18: 128, 19: 48,
                        20: 200, 21: 40, 22: 148, 23: 48, 24: 48}, internal
    assert external != internal
    print('Same 30 physical qutrits, six X checks, different X row-span enumerators:')
    print('  Prakash-Saha m=4,k=6:', external)
    print('  Coupled matter code:', internal)


def test():
    assert len(PAPER_EQ11) == 6 and all(len(r) == 14 for r in PAPER_EQ11)
    for m, k in ((2, 4), (4, 6), (3, 7)):
        full, H1, H0 = construction(m, k)
        assert cubic_valid(full) == (True, None)
        assert len(full) == 3 * m and len(full[0]) == 9 * m
        assert len(H1) == k and len(H0) == 3 * m - k
        assert rank(H1 + H0, 3) == 3 * m
        assert rank(H0, 3) == 3 * m - k
        assert all(sum(a * b for a, b in zip(x, y)) % 3 == 0
                   for x in H0 for y in (H1 + H0))
        if (m, k) == (2, 4):
            assert H1 + H0 == PAPER_EQ11, (H1 + H0, PAPER_EQ11)
        dz, witness = z_distance_up_to_two(H1, H0)
        assert dz == 2, (m, k, dz)
        has_x1, xwit = has_single_qudit_x_logical(H1, H0)
        assert not has_x1, (m, k, xwit)
        n = 9 * m - k
        print(f'Prakash-Saha 2025 m={m}, k={k}: [[{n},{k},2]]_3; '
              f'Z-distance=2; cubic condition verified; Z witness={witness}')
    # For the published family, n+k=9m must be divisible by 9.
    assert (30 + 2) % 9 != 0 and (45 + 3) % 9 != 0
    print('No same (n,k) in [[9m-k,k,2]]_3 family for (30,2) or (45,3).')
    print('Published distance-2 family cannot be equivalent to matter distance-3 codes.')
    compare_30_site_codes()
    print('ALL PUBLISHED QUTRIT REFERENCE CHECKS PASSED')


if __name__ == '__main__':
    test()
