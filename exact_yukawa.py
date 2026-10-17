"""Deterministic support-graph test for exact unbroken family charges.

No floating-point matrices, SVD, or random choices. No claim about spontaneously
broken multi-Higgs/flavon theories or approximate conservation.
"""
from itertools import product, permutations


def allowed(a, b, higgs_up, p, sector):
    """Multiplet charges: Q=a_i, u^c=-a_j+2b_j, d^c=-a_j-2b_j.
    Each nonzero Yukawa entry must have total charge zero mod p.
    """
    sign = +1 if sector == "up" else -1
    eta = higgs_up if sector == "up" else -higgs_up
    return tuple(tuple((a[i] - a[j] + sign * 2 * b[j] + eta) % p == 0
                       for j in range(len(a))) for i in range(len(a)))


def full_rank_possible(M):
    """Exact structural full-rank criterion: a perfect matching exists."""
    n = len(M)
    return any(all(M[i][perm[i]] for i in range(n))
               for perm in permutations(range(n)))


def combined_row_graph_connected(U, D):
    """Necessary condition for irreducible generic CKM mixing.

    An edge joins two left-handed families when an up- or down-type texture
    has some shared right-handed column. This potentially gives a nonzero
    off-diagonal term in Y Y^dagger. A disconnected graph forces an exact
    common block decomposition of both Y Y^dagger.

    Connectivity is NOT sufficient for measured CKM parameters.
    """
    n = len(U)
    seen = {0}
    changed = True
    while changed:
        changed = False
        for i in tuple(seen):
            for j in range(n):
                if j in seen:
                    continue
                if any((U[i][c] and U[j][c]) or (D[i][c] and D[j][c])
                       for c in range(n)):
                    seen.add(j)
                    changed = True
    return len(seen) == n


def exhaustive_1hdm(generations, p):
    """Enumerate F_p charges and one-Higgs assignments with full-rank
    up/down textures and connected combined row-overlap graph.
    The quark sector alone gives the unbroken-symmetry obstruction.
    """
    survivors = []
    for a in product(range(p), repeat=generations):
        for b in product(range(p), repeat=generations):
            for eta in range(p):
                U, D = allowed(a, b, eta, p, 'up'), allowed(a, b, eta, p, 'down')
                if (full_rank_possible(U) and full_rank_possible(D)
                        and combined_row_graph_connected(U, D)):
                    survivors.append((a, b, eta))
    return survivors


def run():
    for generations, p in ((2, 3), (2, 5), (2, 7), (3, 3), (3, 5)):
        survivors = exhaustive_1hdm(generations, p)
        assert len(survivors) == p * p, (generations, p, len(survivors))
        assert all(len(set(a)) == len(set(b)) == 1
                   and (2 * b[0] + eta) % p == 0
                   for a, b, eta in survivors)
        print(f'G={generations}, p={p}: {len(survivors)} solutions; '
              f'all family-universal (exact deterministic support test)')

    up = ((True, False), (False, True))
    assert full_rank_possible(up) and not combined_row_graph_connected(up, up)
    triangular = ((True, True), (False, True))
    assert full_rank_possible(triangular)
    assert combined_row_graph_connected(triangular, up)
    print('ALL EXACT SUPPORT TESTS PASSED')


if __name__ == '__main__':
    run()
