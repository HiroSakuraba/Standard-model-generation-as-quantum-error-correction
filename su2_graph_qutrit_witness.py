"""Explicit exact 3-level quantum code against local S4 exchange algebra.

Eight SU(2) spin-halves at one gauge vertex. Four interact through arbitrary
pair exchanges; the four spectator spins are untouched. This implements the
3-dimensional noiseless multiplicity from j=1 triplet channels. Does NOT
claim physical realization, preparation, or gauge-compatible recovery.

Run: python su2_graph_qutrit_witness.py
"""
from itertools import permutations
import sympy as sp

from su2_error_algebra_capacity import dot, projected, scale_sum

C = (0, 1, 4, 5)
R = (2, 3, 6, 7)
PAIRS = ((0, 1), (0, 2), (0, 3))


def local_triplet_with_singlet(singlet_pair, mag):
    """One 4-spin total-j=1 irrep copy, indexed by spin projection m."""
    i, j = singlet_pair
    a, b = tuple(k for k in range(4) if k not in singlet_pair)
    singlet = {(0, 1): sp.sqrt(2) / 2, (1, 0): -sp.sqrt(2) / 2}
    triplet = {
        +1: {(0, 0): sp.Integer(1)},
        0: {(0, 1): sp.sqrt(2) / 2, (1, 0): sp.sqrt(2) / 2},
        -1: {(1, 1): sp.Integer(1)},
    }[mag]
    out = {}
    for (si, sj), x in singlet.items():
        for (ta, tb), y in triplet.items():
            bits = [0] * 4
            bits[i], bits[j], bits[a], bits[b] = si, sj, ta, tb
            out[tuple(bits)] = x*y
    return out


def orthonormal_triplet_irreps():
    """Find same exact Gram-Schmidt coefficients for m=+1,0,-1."""
    originals = {m: [local_triplet_with_singlet(pair, m)
                     for pair in PAIRS] for m in (+1, 0, -1)}
    gram = sp.Matrix(3, 3, lambda i, j: dot(originals[1][i], originals[1][j]))
    assert gram.det() != 0
    for m in (0, -1):
        assert sp.Matrix(3, 3, lambda i, j: dot(originals[m][i], originals[m][j])) == gram
    coeffs = []
    for a in range(3):
        candidate = sp.eye(3)[:, a]
        for old in coeffs:
            candidate = sp.simplify(candidate - (old.T * gram * candidate)[0] * old)
        norm = sp.sqrt((candidate.T * gram * candidate)[0])
        coeffs.append(sp.simplify(candidate / norm))
    assert sp.simplify(sp.Matrix.hstack(*coeffs).T * gram * sp.Matrix.hstack(*coeffs) - sp.eye(3)) == sp.zeros(3)
    result = {}
    for a in range(3):
        for m in (+1, 0, -1):
            result[(a,m)] = scale_sum(*[(coeffs[a][j], originals[m][j]) for j in range(3)])
    return result


def combine_at_sites(left, right):
    out = {}
    for l, x in left.items():
        for r, y in right.items():
            bits = [0] * 8
            for i, site in enumerate(C):
                bits[site] = l[i]
            for i, site in enumerate(R):
                bits[site] = r[i]
            key = tuple(bits)
            out[key] = out.get(key, 0) + x * y
    return out


def gauge_plus(state):
    out = {}
    for bs, amp in state.items():
        for i, val in enumerate(bs):
            if val:
                t = list(bs)
                t[i] = 0
                t = tuple(t)
                out[t] = out.get(t, 0) + amp
    return {k: sp.simplify(v) for k, v in out.items() if sp.simplify(v) != 0}


def permute_C(state, order):
    out = {}
    for bits, coeff in state.items():
        modified = list(bits)
        for target, source in enumerate(order):
            modified[C[target]] = bits[C[source]]
        key = tuple(modified)
        out[key] = out.get(key, 0) + coeff
    return out


def main():
    triplet = orthonormal_triplet_irreps()
    fixed_component_multiplicity = 0
    code = []
    for spectator_multiplicity in range(3):
        vectors = []
        for magnetic, sign in ((+1, 1), (0, -1), (-1, 1)):
            vectors.append((sign / sp.sqrt(3), combine_at_sites(
                triplet[(fixed_component_multiplicity, magnetic)],
                triplet[(spectator_multiplicity, -magnetic)])))
        psi = scale_sum(*vectors)
        assert gauge_plus(psi) == {}
        assert all(sum(1 if bit == 0 else -1 for bit in k) == 0 for k in psi)
        code.append(psi)
    assert projected(code, lambda psi: psi) == sp.eye(3)
    for order in permutations(range(4)):
        p = projected(code, lambda psi, order=order: permute_C(psi, order))
        assert sp.simplify(p - p[0, 0] * sp.eye(3)) == sp.zeros(3), (order, p)
    # Products of group operators are group operators; 24 scalar
    # compressions imply exact KL for the entire complex group algebra.
    print("EXACT 3D QUTRIT: 24 S4 permutations of the connected four spins all compress to scalars")
    print("All 3 codewords are orthonormal global SU(2) singlets; all Gauss raising and z tests PASS")
    print("Arbitrary error algebra on connected four-spin component is correctable (group closure proof)")
    print("NO claim of encoding against spectator errors or physical recovery hardware")


if __name__ == "__main__":
    main()
