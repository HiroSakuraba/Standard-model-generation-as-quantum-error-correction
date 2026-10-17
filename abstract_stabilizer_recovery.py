"""Exact Knill--Laflamme certificates for the EXISTING formal charge codes.

This tests the stabilizer quantum codes built by charge_lattice_tools.py.
It does NOT establish that any of their X/Z checks are measurable conserved
observables of the Standard Model. It does NOT assert a physical noise model.

For an (n,k,d>=3) nondegenerate CSS code, distinct syndromes for all
single-site Weyl errors E=X_i^a Z_i^b and identity imply
P E^\dagger F P = delta_{EF} P (up to unimportant phases), and thus
exact correction of arbitrary single-site noise by linearity.
"""
from itertools import product

import numpy as np

from charge_lattice_tools import B3, B5, coupled, nullspace, rank, single


def recovery_certificate(code):
    p, n = code.p, code.n
    sx = code.SX % p
    sz = nullspace(code.GP, p)
    assert rank(sx, p) == code.r - code.k
    assert rank(sz, p) == code.n - code.r
    assert (sx @ sz.T % p == 0).all(), "CSS checks do not commute"
    assert rank(np.vstack((sx, sz)), p) <= n  # independently scoped to each type

    # Record each syndrome as (Z-check responses to X, X-check
    # responses to Z), using matrices generated from the actual
    # charge code. An error basis has p^2 Weyl operators per site.
    zero = (0,) * (len(sz) + len(sx))
    lookup = {zero: ("identity",)}
    for j in range(n):
        for a, b in product(range(p), repeat=2):
            if not (a or b):
                continue
            syndrome = tuple(int(x) for x in np.concatenate((
                (a * sz[:, j]) % p, (b * sx[:, j]) % p)))
            assert any(syndrome), (p, n, j, a, b, "undetected single-site error")
            assert syndrome not in lookup, (
                p, n, j, a, b, "KL collision", lookup.get(syndrome))
            lookup[syndrome] = (j, a, b)
    assert len(lookup) == 1 + n * (p * p - 1)

    # Given a measured syndrome, apply the inverse single-site Weyl.
    # A minus sign modulo p in each exponent corrects the detected
    # Pauli error up to its irrelevant global commutation phase.
    for syndrome, error in lookup.items():
        if error == ("identity",):
            assert syndrome == zero
            continue
        j, a, b = error
        recovery = (j, (-a) % p, (-b) % p)
        assert (a + recovery[1]) % p == (b + recovery[2]) % p == 0

    params = code.params()
    print(f"{code.name}: [[{n},{code.k},>=3]]_{p} "
          f"(X checks={params['X_stab']}, Z checks={params['Z_stab']}); "
          f"{len(lookup)} distinct identity/one-site Weyl syndromes; "
          "exact KL/recovery lookup PASS")
    return len(lookup)


def test():
    for p, B in ((3, B3), (5, B5)):
        recovery_certificate(single(p, name=f"single family p={p}"))
        for nfam in (2, 3):
            c = coupled(B, nfam, p, name=f"B{p}-coupled {nfam} families")
            recovery_certificate(c)
    print("ALL ABSTRACT CHARGE-CODE RECOVERY CHECKS PASSED")


if __name__ == "__main__":
    test()
