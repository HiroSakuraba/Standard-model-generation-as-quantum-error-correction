"""Exact four-site reduced-state purity fingerprints of coupled qudit CSS codes.

The state is the *normalized projector* onto the code space, not a chosen
encoded pure state. For every four-coordinate subset A,

  Tr(rho_A^2) = p^(dim S_A - 4),

where S_A consists of Pauli stabilizers supported entirely on A.

All family Cartan-linear charge rows restrict on each family to a linear form
in five SO(10) signs. For p=3,5 a nonzero such form has weight >=7 on
the fifteen unpunctured species (verified exhaustively). Consequently no
X stabilizer can be supported on four sites; dim S_A is exactly the nullity
of the physical charge submatrix G_P|A. These purities are invariant under
arbitrary local unitaries and qudit permutations.

Usage: python3 local_purity.py
Requires numpy via charge_lattice_tools.py. No random sampling.
"""
from collections import Counter
from itertools import combinations, product
from math import comb

from charge_lattice_tools import B3, B5, W16, coupled, direct_sum


def small_rank(columns, p):
    """Rank of a handful of column vectors over the prime field F_p."""
    pivots = {}
    for column in columns:
        v = [int(x) % p for x in column]
        for i in sorted(pivots):
            if v[i]:
                scale = v[i]
                v = [(a - scale * b) % p for a, b in zip(v, pivots[i])]
        pivot = next((i for i, x in enumerate(v) if x), None)
        if pivot is None:
            continue
        inv = pow(v[pivot], -1, p)
        reduced = [(x * inv) % p for x in v]
        for j, row in pivots.items():
            if row[pivot]:
                factor = row[pivot]
                pivots[j] = [(x - factor * y) % p for x, y in zip(row, reduced)]
        pivots[pivot] = reduced
    return len(pivots)


def minimum_nonzero_one_family_cartan_weight(p):
    """Enumerate the five-dimensional Cartan row space on 15 species.

    A nonzero form on a single family has no identically zero punctured
    evaluation (rank = 5). Its minimum weight bounds the support of any
    nonzero multi-family Cartan-linear X stabilizer below.
    """
    phys = W16[1:]
    best = 15
    for coefficients in product(range(p), repeat=5):
        if not any(coefficients):
            continue
        weight = sum(
            sum(a * b for a, b in zip(coefficients, q)) % p != 0
            for q in phys
        )
        best = min(best, weight)
    return best


def four_site_profile(code):
    """Counts of 4-subsets by dim S_A; family partition is diagnostic only."""
    n, p = code.n, code.p
    col = [tuple(int(x) for x in code.GP[:, i]) for i in range(n)]
    families = [int(x) for x in code.families]
    nfam = max(families) + 1
    dimensions = Counter()
    parts = Counter()
    for indices in combinations(range(n), 4):
        nullity = 4 - small_rank([col[i] for i in indices], p)
        dimensions[nullity] += 1
        if nullity:
            fam = tuple(sum(families[i] == g for i in indices)
                        for g in range(nfam))
            parts[fam] += 1
    assert sum(dimensions.values()) == comb(n, 4)
    assert set(dimensions).issubset({0, 1}), dimensions
    return dict(dimensions), dict(sorted(parts.items()))


def run():
    for p, B, mixed in ((3, B3, 34), (5, B5, 99)):
        minimum = minimum_nonzero_one_family_cartan_weight(p)
        assert minimum == 7, (p, minimum)
        print(f'F_{p}: every nonzero family-local Cartan form has weight >= {minimum}')
        for nfam in (2, 3):
            reference = direct_sum(nfam, p)
            joined = coupled(B, nfam, p)
            direct_counts, direct_parts = four_site_profile(reference)
            coupled_counts, coupled_parts = four_site_profile(joined)

            expected_direct = 75 * nfam
            expected_joined = expected_direct + mixed * comb(nfam, 2)
            assert direct_counts.get(1, 0) == expected_direct, direct_counts
            assert coupled_counts.get(1, 0) == expected_joined, coupled_counts
            assert len([x for x in direct_parts if len([i for i in x if i]) > 1]) == 0
            for a, b in combinations(range(nfam), 2):
                key = tuple(2 if f in (a, b) else 0 for f in range(nfam))
                assert coupled_parts.get(key, 0) == mixed, (p, nfam, key, coupled_parts)
                assert direct_parts.get(key, 0) == 0
            print(f'p={p} families={nfam} n={15 * nfam}: '
                  f'4-site purity hotspots direct={expected_direct}, '
                  f'coupled={expected_joined}; mixed 2+2 per pair={mixed}')
            # For these models, all hotspot subsets have dim(S_A)=1.
            # Hence Tr(rho_A^2) = p^-3 rather than p^-4.
    print('ALL LOCAL-PURITY CHECKS PASSED')


if __name__ == '__main__':
    run()
