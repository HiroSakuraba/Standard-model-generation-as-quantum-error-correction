"""Walsh-character rigidity and a deceptive two-family nonlinear lattice.

Pure 16-point character enumeration plus exact finite-field verification
against the repo's ChargeCode. No physical-family symmetry is inferred.

Run: python3 walsh_character_lattice.py (numpy required for ChargeCode).
"""
from collections import Counter
from itertools import combinations, combinations_with_replacement

import numpy as np

from charge_lattice_tools import W16, build, direct_sum, split_components

# The fifth SO(10) sign is the product of the first four:
# q5 = q1*q2*q3*q4 on all 16 even-parity sign vectors.
GAUGE_MASKS = (1, 2, 4, 8, 15)
GAUGE_SET = frozenset(GAUGE_MASKS)


def character(q, mask):
    assert 0 <= mask < 16
    value = 1
    for i in range(4):
        if mask & (1 << i):
            value *= q[i]
    return value


def binary_rank(values):
    pivots = [0] * 4
    rank = 0
    for v in values:
        x = v
        for j in range(3, -1, -1):
            if x & (1 << j):
                if pivots[j]:
                    x ^= pivots[j]
                else:
                    pivots[j] = x
                    rank += 1
                    break
    return rank


def compatible(masks):
    """Cubic orthogonality of the five shared q_i and antisymmetric
    Walsh rows (chi_a,-chi_a). The only nontrivial cubic condition:
    chi_a*chi_b != q_i for every pair a!=b and every i.
    """
    return len(set(masks)) == len(masks) and all(
        (a ^ b) not in GAUGE_SET for a, b in combinations(masks, 2)
    )


def matched_twin_faults(masks):
    """Number of matched non-singlet twins producing weight-2 logical Z.

    For two-family columns (q, +-chi_a(q)), a same-species pair with
    coefficients +1,-1 is logical exactly when all chi_a(q) agree.
    Their relative equations have binary rank r, leaving 2^(4-r)
    sign assignments; subtract the punctured all-plus singlet.
    """
    if not masks:
        return 15
    r = binary_rank([mask ^ masks[0] for mask in masks[1:]])
    return (1 << (4 - r)) - 1


def classify_admissible_character_sets():
    valid = []
    counts = Counter()
    for size in range(1, 7):
        for masks in combinations(range(16), size):
            if compatible(masks):
                counts[size] += 1
                valid.append(masks)
    assert [counts[s] for s in range(1, 7)] == [16, 80, 160, 120, 16, 0]
    maximals = [frozenset(S) for S in valid if len(S) == 5]
    translates = {frozenset(a ^ t for a in GAUGE_MASKS) for t in range(16)}
    assert set(maximals) == translates
    assert len(translates) == 16
    assert all(matched_twin_faults(S) == 0 for S in maximals)
    assert all(matched_twin_faults(S) >= 1 for S in valid if len(S) < 5)
    print("Compatible mask sets by size 1..6:", [counts[s] for s in range(1, 7)])
    print("Every 5-mask solution is a translate t XOR GAUGE_MASKS; "
          "all smaller solutions leave matched-twin weight-2 faults.")


def single_family_rigidity():
    """The characters chi_a span all functions on the 16-point sign set.
    The 11 products q_i*q_j have masks 0 plus ten distinct pair XORs,
    so their annihilator is exactly span of the five q_i, for odd p.
    """
    assert len(W16) == 16
    assert len(set(W16)) == 16
    assert all(q[4] == q[0] * q[1] * q[2] * q[3] for q in W16)
    for a in range(16):
        for b in range(16):
            val = sum(character(q, a) * character(q, b) for q in W16)
            assert val == (16 if a == b else 0)
    products = {a ^ b for a, b in combinations_with_replacement(GAUGE_MASKS, 2)}
    assert len(products) == 11
    assert GAUGE_SET.isdisjoint(products)
    assert products | GAUGE_SET == set(range(16))
    for p in (3, 5, 7, 11):
        assert 16 % p != 0
    print("One-generation rigidity: 11 pair-product characters + "
          "5 Cartan characters form a Walsh basis; no nonlinear "
          "single-family extension for any odd prime.")


def candidate(t, p):
    """Shared Cartan rows; differential rows chi_(t XOR gauge mask)."""
    masks = tuple(t ^ a for a in GAUGE_MASKS)
    def family(sign):
        def f(q):
            return list(q) + [sign * character(q, m) for m in masks]
        return f
    return build([family(1), family(-1)], p), masks


def prove_explicit_block_permutation(t, p):
    """A row transformation sends (q, +-chi_t(q)*q) to (q,0) or
    (0,q). Relabel by that effective block to recover the direct sum.
    This is a constructive equivalence, not a numerical heuristic.
    """
    code, masks = candidate(t, p)
    reference = direct_sum(2, p)
    assert code.n == 30 and code.k == 2 and code.r == 10
    assert len(code.punct) == 2
    inv2 = pow(2, -1, p)
    transformed = np.vstack(((code.GP[:5] + code.GP[5:]) * inv2 % p,
                             (code.GP[:5] - code.GP[5:]) * inv2 % p))
    permutation = []
    for j in range(30):
        family, species = divmod(j, 15)
        q = W16[species + 1]
        switched = character(q, t) == -1
        permuted = (family ^ switched) * 15 + species
        permutation.append(permuted)
        assert np.array_equal(transformed[:, j], reference.GP[:, permuted])
    assert set(permutation) == set(range(30))
    assert sorted(map(len, split_components(code))) == [15, 15]
    return code, masks


def logical_triplet_partitions(code):
    """For this *already proved direct-sum* code, all 50 weight-3
    logical supports are images of the 25 one-family supports under
    the explicit permutation. No cross-output error occurs at weight 3.
    """
    out = Counter()
    for support in code.z_logical_supports(3):
        out[sum(j < 15 for j in support)] += 1
    assert sum(out.values()) == 50
    return out


def integer_triorthogonality(code):
    G = code.Gint
    for i, j, k in combinations_with_replacement(range(G.shape[0]), 3):
        assert sum(int(x) * int(y) * int(z)
                   for x, y, z in zip(G[i], G[j], G[k])) == 0


def test():
    single_family_rigidity()
    classify_admissible_character_sets()

    example = (0, 3, 5, 9, 14)
    assert tuple(1 ^ a for a in GAUGE_MASKS) == example
    for p in (3, 5):
        for t in range(16):
            code, masks = prove_explicit_block_permutation(t, p)
            assert compatible(masks)
            assert matched_twin_faults(masks) == 0
            if t in (0, 1, 3):
                integer_triorthogonality(code)
                parts = logical_triplet_partitions(code)
                if t == 0:
                    assert dict(parts) == {0: 25, 3: 25}
                if t == 1:
                    assert dict(parts) == {0: 3, 1: 22, 2: 22, 3: 3}
                if t == 3:
                    assert dict(parts) == {0: 6, 1: 19, 2: 19, 3: 6}
                print(f"p={p}, mask t={t}: [[30,2,3]]_{p}; "
                      f"3-site supports by original family count {dict(sorted(parts.items()))}; "
                      "exactly two 15-site stabilizer components")
    print("ALL WALSH-CHARACTER LATTICE CHECKS PASSED")


if __name__ == "__main__":
    test()
