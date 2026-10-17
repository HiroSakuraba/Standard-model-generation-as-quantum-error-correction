"""Exact outer recovery for five independent SU(2) matter/rishon singlet modules.

Each module is the two-dimensional fixed-endpoint SU(2) singlet sector
from four_link_su2_subspace.py. Projected local color-spin exchanges
generate its full logical M2 algebra: one module alone is NOT protected
against generic local exchange.

Concatenate five such modules with the *previously known* [[5,1,3]]_2
perfect stabilizer code (Laflamme et al., Phys Rev Lett 77, 198, 1996).
This is a standard outer code, NOT a new Standard Model or SU(2) code.

Certificate: exactly construct two 32-component normalized codewords
from four commuting stabilizers; verify all 16 identity/single-module
Pauli syndromes distinct and all 16**2*2**2 Knill--Laflamme overlaps
against rational amplitudes. A decoder table corrects an arbitrary
one-module exchange channel by linearity in its Weyl/Pauli expansion.

The mathematical result extends to arbitrary errors acting on one
physical module, even if they leak out of the internal 2D singlet
sector, by the standard isometric-encoding KL argument (documented
separately). Physical nonlocal stabilizer measurements, gauge-compatible
control, and a native Hamiltonian have NOT been constructed.

Requires SymPy for the local SU(2) exchange-matrix provenance.
Run: python five_block_su2_recovery.py
"""
from collections import defaultdict
from fractions import Fraction
from itertools import product

import sympy as sp

from four_link_su2_subspace import ORIENTATION
from su2_exchange_algebra import (
    I2, X2, Z2, basis_for_sector, on_site_swap,
)

GENERATORS = ("XZZXI", "IXZZX", "XIXZZ", "ZXIXZ")
LOGICAL_X = "XXXXX"
LOGICAL_Z = "ZZZZZ"
NMODULES = 5


def masks(word):
    """Pauli word represented as X^x Z^z; Y = i XZ.

    The operator XZ is a convenient *real orthogonal* Pauli basis
    element, differing from Hermitian Y only by a global i phase.
    """
    assert len(word) == NMODULES
    x = z = 0
    for i, char in enumerate(word):
        assert char in "IXYZ"
        if char in "XY":
            x |= 1 << i
        if char in "ZY":
            z |= 1 << i
    return x, z


def commutation(a, b):
    """Symplectic parity of Pauli commutation."""
    x1, z1 = masks(a)
    x2, z2 = masks(b)
    return ((x1 & z2).bit_count() + (z1 & x2).bit_count()) % 2


def apply(word, vector):
    x, z = masks(word)
    out = defaultdict(Fraction)
    for state, amplitude in vector.items():
        phase = -1 if (state & z).bit_count() % 2 else 1
        out[state ^ x] += phase * amplitude
    return {state: amp for state, amp in out.items() if amp}


def add(a, b):
    out = defaultdict(Fraction, a)
    for index, amplitude in b.items():
        out[index] += amplitude
    return {index: amp for index, amp in out.items() if amp}


def inner(a, b):
    return sum((amp * b.get(state, 0) for state, amp in a.items()),
               Fraction(0))


def logical_codewords():
    """Project computational |00000> into +1 common eigenspace."""
    ket0 = {0: Fraction(1)}
    for generator in GENERATORS:
        ket0 = add(ket0, apply(generator, ket0))
    # 16 different signed computational states: normalized by 1/4.
    assert len(ket0) == 16
    ket0 = {state: amp / 4 for state, amp in ket0.items()}
    ket1 = apply(LOGICAL_X, ket0)
    assert inner(ket0, ket0) == inner(ket1, ket1) == 1
    assert inner(ket0, ket1) == 0
    assert apply(LOGICAL_Z, ket0) == ket0
    assert apply(LOGICAL_Z, ket1) == {
        state: -amp for state, amp in ket1.items()
    }
    for psi in (ket0, ket1):
        for generator in GENERATORS:
            assert apply(generator, psi) == psi
    return ket0, ket1


def one_module_error_basis():
    """16 error operators: identity and three per module.

    Local XY error is represented by XZ, not iXZ. This changes no
    syndrome and is a complete operator basis for one qubit.
    """
    errors = {"identity": "I" * NMODULES}
    for i in range(NMODULES):
        for op in "XYZ":
            chars = ["I"] * NMODULES
            chars[i] = op
            errors[f"{op}{i}"] = "".join(chars)
    assert len(errors) == 16
    return errors


def syndrome(error):
    return tuple(commutation(error, gen) for gen in GENERATORS)


def exact_knill_laflamme():
    """Check ALL 1,024 matrix elements <u|E†F|v>.

    If the 16 error spaces are mutually orthogonal, recovery applies
    to arbitrary linear combinations of errors on an unknown block,
    including all coherent local exchanges and dephasing operators.
    """
    assert all(not commutation(a, b) for a in GENERATORS for b in GENERATORS)
    basis = logical_codewords()
    errors = one_module_error_basis()
    syndromes = {}
    for name, pauli in errors.items():
        sig = syndrome(pauli)
        assert sig not in syndromes, (name, sig, syndromes.get(sig))
        syndromes[sig] = name
    assert len(syndromes) == 16
    assert syndrome("IIIII") == (0, 0, 0, 0)
    assert all(any(syndrome(op) ) for name, op in errors.items() if name != "identity")
    assert all(not commutation(LOGICAL_X, gen) for gen in GENERATORS)
    assert all(not commutation(LOGICAL_Z, gen) for gen in GENERATORS)
    assert commutation(LOGICAL_X, LOGICAL_Z) == 1

    error_images = {
        name: tuple(apply(word, v) for v in basis)
        for name, word in errors.items()
    }
    tested = 0
    for name_a, va in error_images.items():
        for name_b, vb in error_images.items():
            for u, v in product(range(2), repeat=2):
                overlap = inner(va[u], vb[v])
                expected = Fraction(int(name_a == name_b and u == v))
                assert overlap == expected, (
                    name_a, name_b, u, v, overlap, expected)
                tested += 1
    assert tested == 16 * 16 * 2 * 2
    # All 32 vectors E |u_L> constitute an orthonormal basis of the
    # entire 32-dim five-module logical Hilbert space. This is the
    # saturated nondegenerate quantum Hamming bound.
    print("Five-module [[5,1,3]]_2: 16 unique 4-bit syndromes.")
    print(f"Exact Pauli-basis KL comparisons: {tested} / {tested} PASS")
    print("Error recovery: read syndrome, identify single-module "
          "Pauli, apply its inverse (global phase irrelevant).")
    return syndromes


def source_local_gauge_exchange_matrices():
    """Verify outer X/Z are built from actual on-site SU(2) exchanges.

    S13=SWAP(m1,r1)=-Z, S14=(Z-sqrt(3) X)/2,
    S15=(Z+sqrt(3) X)/2 in the chosen 4-spin singlet basis.
    Each swap is between co-located local doublets, commutes with
    Gauss law and preserves the two-color-singlet subspace.
    """
    states, rich, poor = basis_for_sector(ORIENTATION)
    assert rich == (1, 3, 4, 5)
    assert poor == (0, 2)
    s13 = on_site_swap(states, 1, 3)
    s14 = on_site_swap(states, 1, 4)
    s15 = on_site_swap(states, 1, 5)
    z = -s13
    x = sp.simplify((s15 - s14) / sp.sqrt(3))
    assert z == Z2 and x == X2
    assert x * z == -z * x
    assert sp.simplify(x * x - I2) == sp.zeros(2)
    assert sp.simplify(z * z - I2) == sp.zeros(2)
    print("On-site gauge-invariant swaps synthesize exact "
          "encoded X and Z on each SU(2) singlet module.")
    print("Outer XZZXI etc. are gauge-invariant as products, "
          "but cross-module syndrome measurement is not implemented.")
    return x, z


def test():
    source_local_gauge_exchange_matrices()
    table = exact_knill_laflamme()
    print("Syndrome table:", dict(sorted(table.items())))
    print("ALL FIVE-BLOCK SU2 RECOVERY CHECKS PASSED")


if __name__ == "__main__":
    test()
