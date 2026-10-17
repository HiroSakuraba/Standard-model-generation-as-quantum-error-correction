"""Reproducible, equivalence-scoped fingerprints for family-coupled CSS codes.

Usage: python3 code_fingerprints.py [--extended]
Dependencies: numpy and the repository's charge_lattice_tools.py.

The fingerprints distinguish SELECTED reference codes and support the existing
projector-factorization argument against direct sums under local unitaries;
they do not establish novelty among all known indecomposable codes.
"""
from __future__ import annotations

import json
import sys
from itertools import combinations

import numpy as np

from charge_lattice_tools import (
    B3, B5, coupled, direct_sum, rref, rank, split_components, vectors_range,
)


def stabilizer_x_enumerator(code):
    """Complete X-stabilizer Hamming-weight enumerator.

    Invariant under coordinate permutations and nonzero local rescalings that
    preserve CSS X. Not invariant under general local Clifford X/Z rotations.
    """
    basis, _ = rref(code.SX, code.p)
    dimension = len(basis)
    counts = np.zeros(code.n + 1, dtype=np.int64)
    for start in range(0, code.p ** dimension, 65536):
        coeffs = vectors_range(dimension, code.p, start, start + 65536)
        words = (coeffs @ basis) % code.p
        counts += np.bincount(np.count_nonzero(words, axis=1),
                              minlength=code.n + 1)
    assert int(counts.sum()) == code.p ** dimension
    assert int(counts[0]) == 1
    return {str(i): int(c) for i, c in enumerate(counts) if c}


def family_partition(code, weight):
    counts = {}
    for support in code.z_logical_supports(weight):
        parts = tuple(sum(code.families[j] == fam for j in support)
                      for fam in range(max(code.families) + 1))
        key = '-'.join(str(x) for x in parts)
        counts[key] = counts.get(key, 0) + 1
    return dict(sorted(counts.items()))


def mixed_two_plus_two(code):
    """Check only 2+2 cross-family support candidates (no 3+1 supersets)."""
    nf = max(code.families) + 1
    by_family = [[j for j, f in enumerate(code.families) if f == i]
                 for i in range(nf)]
    r_removed = rank(code.GR, code.p)
    total = 0
    for f1, f2 in combinations(range(nf), 2):
        for ab in combinations(by_family[f1], 2):
            for cd in combinations(by_family[f2], 2):
                cols = code.GP[:, list(ab + cd)]
                if (rank(np.column_stack((code.GR, cols)), code.p)
                        < r_removed + rank(cols, code.p)):
                    total += 1
    return total


def fingerprint(name, code):
    return {
        'name': name,
        'p': code.p,
        'n': code.n,
        'k': code.k,
        'x_stabilizers': rank(code.SX, code.p),
        'z_stabilizers': code.n - code.r,
        'stabilizer_components': sorted(map(len, split_components(code))),
        'x_stabilizer_weight_enumerator': stabilizer_x_enumerator(code),
        'logical_z_supports_weight_3_by_family': family_partition(code, 3),
        'mixed_2_plus_2_logical_z_supports': mixed_two_plus_two(code),
    }


def main():
    models = [
        ('direct-2-qutrit', direct_sum(2, 3)),
        ('coupled-2-qutrit', coupled(B3, 2, 3)),
        ('direct-2-p5', direct_sum(2, 5)),
        ('coupled-2-p5', coupled(B5, 2, 5)),
    ]
    if '--extended' in sys.argv:
        models.extend([
            ('coupled-3-qutrit', coupled(B3, 3, 3)),
            ('coupled-3-p5', coupled(B5, 3, 5)),
        ])
    results = {name: fingerprint(name, code) for name, code in models}
    for p, tag in ((3, 'qutrit'), (5, 'p5')):
        direct = results[f'direct-2-{tag}']
        joined = results[f'coupled-2-{tag}']
        assert direct['stabilizer_components'] == [15, 15]
        assert joined['stabilizer_components'] == [30]
        assert direct['mixed_2_plus_2_logical_z_supports'] == 0
        assert joined['mixed_2_plus_2_logical_z_supports'] == (88 if p == 3 else 396)
        assert direct['k'] == joined['k'] == 2
        assert (direct['x_stabilizer_weight_enumerator']
                != joined['x_stabilizer_weight_enumerator'])
    if '--extended' in sys.argv:
        assert results['coupled-3-qutrit']['stabilizer_components'] == [45]
        assert results['coupled-3-p5']['stabilizer_components'] == [45]
        assert results['coupled-3-qutrit']['mixed_2_plus_2_logical_z_supports'] == 264
        assert results['coupled-3-p5']['mixed_2_plus_2_logical_z_supports'] == 1188
    print(json.dumps(results, indent=2, sort_keys=True))
    print('ALL FINGERPRINT CHECKS PASSED')


if __name__ == '__main__':
    main()
