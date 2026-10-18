#!/usr/bin/env python3
"""E8 comparison: chirality obstruction <-> code obstruction (see E8_COMPARISON.md).

Standard branching E8 > SO(16) > SO(10) x SO(6), with 248 = 120 + 128 and
128 = (16,4) + (16bar,4bar).  The 128 half-integer roots (scaled by 2) are the
sign vectors in {+1,-1}^8 with an even number of minus signs; the first five
coordinates are the SO(10) weight, the last three the SO(6) = SU(4) weight.

Checks (exact integer arithmetic, plus mod-p ranks from charge_lattice_tools):
 1. the chiral block (16,4) is NOT triorthogonal (anomalous) over Z, mod 3, 5, 7;
 2. the full non-chiral 128 IS triorthogonal, but puncturing the four nu^c columns
    leaves their mirror partners as weight-one logicals: distance 1;
 3. ignoring the anomaly, using the SU(4) weights as additive family labels on the
    chiral block reproduces the 'uniform family badge': [[60,4,2]]_p.
Requires numpy.  Ends with ALL E8 COMPARISON CHECKS PASSED.
"""
import itertools
import numpy as np
from charge_lattice_tools import W16, Q0, ChargeCode, build

ROOTS = [v for v in itertools.product((1, -1), repeat=8) if v.count(-1) % 2 == 0]
SPINOR10 = set(W16)
SU4 = [w for w in itertools.product((1, -1), repeat=3) if w.count(-1) % 2 == 0]


def first_cubic_failure(cols, mod=None):
    G = np.array(cols, dtype=np.int64).T
    for i in range(8):
        for j in range(i, 8):
            for k in range(j, 8):
                s = int((G[i] * G[j] * G[k]).sum())
                if (s if mod is None else s % mod) != 0:
                    return (i, j, k, s)
    return None


def main():
    assert len(ROOTS) == 128 and len(SU4) == 4
    chiral = [v for v in ROOTS if tuple(v[:5]) in SPINOR10]
    mirror = [v for v in ROOTS if tuple(v[:5]) not in SPINOR10]
    assert len(chiral) == len(mirror) == 64
    assert set(tuple(-x for x in v) for v in chiral) == set(mirror)   # mirror = negatives
    print('128 = (16,4) + (16bar,4bar); the second block is exactly the negatives of the first')

    # 1. anomaly of the chiral block
    fail = first_cubic_failure(chiral)
    assert fail == (5, 6, 7, 64)
    for p in (3, 5, 7):
        assert first_cubic_failure(chiral, p) is not None
    print('chiral (16,4): cubic sum over the three SU(4) Cartan rows = 64 != 0 '
          '(fails over Z and mod 3, 5, 7) -> not a valid triorthogonal charge lattice')

    # 2. full non-chiral block
    assert first_cubic_failure(ROOTS) is None
    for p in (3, 5, 7):
        punct = [ROOTS.index(tuple(Q0) + w) for w in SU4]
        c = ChargeCode(np.array(ROOTS).T, p, punct)
        assert (c.n, c.k) == (124, 4) and c.z_distance(2) == 1
        print(f'full 128, p={p}: triorthogonal over Z, but [[{c.n},{c.k},1]]_{p} '
              '(each punctured nu^c has its mirror among the physical columns)')

    # 3. additive SU(4) family labels on the chiral block (anomaly ignored)
    for p in (3, 5, 7):
        c = build([(lambda w: (lambda q: list(q) + list(w)))(w) for w in SU4], p)
        assert (c.n, c.k) == (60, 4) and c.z_distance(3) == 2
        print(f'chiral block with additive SU(4) family labels, p={p}: [[60,4,2]]_{p} '
              '(the uniform family-badge twins problem)')
    print('ALL E8 COMPARISON CHECKS PASSED')


if __name__ == '__main__':
    main()
