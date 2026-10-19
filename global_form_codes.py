#!/usr/bin/env python3
"""The global form of the Standard Model gauge group as a choice of stabilizer code
(see GLOBAL_FORM_CODES.md).

The gauge group is (SU(3) x SU(2) x U(1)_Y)/Gamma with Gamma in {1, Z2, Z3, Z6} (Tong,
arXiv:1705.01853).  Line operators carry electric labels (z_e2, z_e3, q) and magnetic labels
(z_m2, z_m3, g) with Tong's Dirac condition (his eq. 3.1)
        3 z_e2 z_m2 + 2 z_e3 z_m3 - 6 q g  in 6Z.
Checks:
 1. Modulo screening by SM matter, the electric label reduces to a class (a2, a3) in Z2 x Z3 and
    the magnetic label to (b2, b3) = (z_m2, z_m3); Tong's pairing becomes 3 a2 b2 + 2 a3 b3 (mod 6).
 2. A maximal mutually local set of lines is a Lagrangian subgroup of Z6 x Z6 (order 6).  There
    are exactly 12; sorting by the magnetic image gives Gamma = 1, Z2, Z3, Z6 with 1, 2, 3, 6 each.
 3. Theta-angle shifts (Witten effect, a -> a + b) permute them; the orbits are exactly the four
    Gamma, of sizes 1, 2, 3, 6.
 4. Qudit dictionary: line (a, b) -> Weyl operator Z^a X^b on a qubit (x) qutrit = one 6-level
    qudit.  Mutual locality = commutation.  Each global form is a stabilizer basis (6 states);
    theta shifts are the Clifford phase gates S (qubit) and diag(1,1,omega) (qutrit); the
    electric-magnetic swap is the Fourier gate and exchanges Gamma = 1 with Gamma = Z6.
 5. On a spatial 3-torus the flux sectors of each global form are L^3: 216 sectors in every
    case, and the overlap between two forms is |L cap L'|^3.
 6. Link to SM_CENTER_ECHO.md: in the theory with quotient Gamma = Z_q the centre echo acts
    with 6/q distinct phases on the allowed links; for Gamma = Z6 it is the identity.
Requires numpy.  Ends with ALL GLOBAL-FORM CODE CHECKS PASSED.
"""
from fractions import Fraction as F
import itertools
import numpy as np

SM = {  # (z_e2 = SU(2) N-ality, z_e3 = triality, q = 6Y)
    'Q': (1, 1, 1), 'u^c': (0, 2, -4), 'd^c': (0, 2, 2), 'L': (1, 0, -3),
    'e^c': (0, 0, 6), 'nu^c': (0, 0, 0), 'H': (1, 0, 3),
}


def tong_pairing(e, m):
    """Tong eq. (3.1) with 6g given as an integer G = 6g."""
    ze2, ze3, q = e
    zm2, zm3, G = m
    return (3 * ze2 * zm2 + 2 * ze3 * zm3 - q * G) % 6


def magnetic_labels():
    """Magnetic labels local with all SM matter: (z_m2, z_m3, 6g) with 6g = 3 z_m2 + 2 z_m3 (mod 6)."""
    out = []
    for zm2 in range(2):
        for zm3 in range(3):
            for G in range(-12, 13):
                if (G - 3 * zm2 - 2 * zm3) % 6 == 0:
                    out.append((zm2, zm3, G))
    return out


def electric_class(e):
    ze2, ze3, q = e
    return ((ze2 - q) % 2, (ze3 - q) % 3)


def check_reduction():
    mags = magnetic_labels()
    for f, e in SM.items():
        assert electric_class(e) == (0, 0), f
        assert all(tong_pairing(e, m) == 0 for m in mags), f
    # matter is local with *only* these magnetic labels (among a wider set)
    wider = [(a, b, G) for a in range(2) for b in range(3) for G in range(-12, 13)]
    local = [m for m in wider if all(tong_pairing(e, m) == 0 for e in SM.values())]
    assert sorted(local) == sorted(mags)
    for ze2 in range(2):
        for ze3 in range(3):
            for q in range(-12, 13):
                a2, a3 = electric_class((ze2, ze3, q))
                for m in mags:
                    assert tong_pairing((ze2, ze3, q), m) == (3 * a2 * m[0] + 2 * a3 * m[1]) % 6
    print('SM matter has electric class (0,0); magnetic lines local with all matter obey 6g = 3z_m2 + 2z_m3 mod 6')
    print("Tong's pairing reduces to <(a,b),(a',b')> = 3(a2 b2' - b2 a2') + 2(a3 b3' - b3 a3') mod 6 on classes")


# dyonic class x = (a2, b2, a3, b3)
ELEMENTS = [(a2, b2, a3, b3) for a2 in range(2) for b2 in range(2) for a3 in range(3) for b3 in range(3)]


def add(x, y):
    return ((x[0] + y[0]) % 2, (x[1] + y[1]) % 2, (x[2] + y[2]) % 3, (x[3] + y[3]) % 3)


def form(x, y):
    return (3 * (x[0] * y[1] - x[1] * y[0]) + 2 * (x[2] * y[3] - x[3] * y[2])) % 6


def span(gens):
    S = {(0, 0, 0, 0)}
    while True:
        new = {add(s, g) for s in S for g in gens} | S
        if new == S:
            return frozenset(S)
        S = new


def lagrangians():
    subs = set()
    for r in range(1, 3):
        for gens in itertools.combinations(ELEMENTS, r):
            S = span(gens)
            if len(S) == 6 and all(form(x, y) == 0 for x in S for y in S):
                subs.add(S)
    return sorted(subs, key=lambda S: sorted(S))


def gamma_order(L):
    return len({(x[1], x[3]) for x in L})


def theta2(x):
    return ((x[0] + x[1]) % 2, x[1], x[2], x[3])


def theta3(x):
    return (x[0], x[1], (x[2] + x[3]) % 3, x[3])


def check_lagrangians():
    Ls = lagrangians()
    assert len(Ls) == 12
    by_q = {}
    for L in Ls:
        q = gamma_order(L)
        by_q.setdefault(q, []).append(L)
        assert len({(x[0], x[2]) for x in L if x[1] == 0 and x[3] == 0}) == 6 // q
    assert {q: len(v) for q, v in by_q.items()} == {1: 1, 2: 2, 3: 3, 6: 6}
    print('Lagrangian subgroups of Z6 x Z6 (maximal mutually local line sets): 12')
    print('  by magnetic image |Gamma|: q=1: 1, q=2: 2, q=3: 3, q=6: 6;  Wilson lines allowed: 6/q classes')
    # theta orbits
    idx = {L: i for i, L in enumerate(Ls)}
    orbits = []
    seen = set()
    for L in Ls:
        if L in seen:
            continue
        orb = {L}
        frontier = [L]
        while frontier:
            M = frontier.pop()
            for T in (theta2, theta3):
                N = frozenset(T(x) for x in M)
                assert N in idx
                if N not in orb:
                    orb.add(N); frontier.append(N)
        seen |= orb
        orbits.append(orb)
        assert len({gamma_order(M) for M in orb}) == 1
    sizes = sorted(len(o) for o in orbits)
    assert sizes == [1, 2, 3, 6]
    for o in orbits:
        assert len(o) == gamma_order(next(iter(o)))
    print('  theta shifts (Witten effect a -> a + b) have orbits of sizes 1, 2, 3, 6: one orbit per Gamma')
    return Ls


# ------------------------------------------------------------------- qudit dictionary
def weyl(x):
    a2, b2, a3, b3 = x
    Z2 = np.diag([1, -1]).astype(complex); X2 = np.array([[0, 1], [1, 0]], dtype=complex)
    w = np.exp(2j * np.pi / 3)
    Z3 = np.diag([1, w, w * w]); X3 = np.roll(np.eye(3), 1, axis=0).astype(complex)
    mp = np.linalg.matrix_power
    return np.kron(mp(Z2, a2) @ mp(X2, b2), mp(Z3, a3) @ mp(X3, b3))


def same_ray_sets(U, V, tol=1e-9):
    """Columns of U and V span the same set of rays (up to order and phase)."""
    O = np.abs(U.conj().T @ V)
    return np.allclose(np.sort(O, axis=1)[:, -1], 1, atol=tol)


def basis_of(L):
    gens = [x for x in L if len(span([x])) == 6]
    W = weyl(gens[0])
    # make it normal-diagonalisable with distinct eigenvalues: W^6 is a scalar
    W6 = np.linalg.matrix_power(W, 6)
    assert np.allclose(W6, W6[0, 0] * np.eye(6))
    w, V = np.linalg.eig(W)
    assert min(abs(w[i] - w[j]) for i in range(6) for j in range(i + 1, 6)) > 1e-6
    V, _ = np.linalg.qr(V)
    # V diagonalises every element of L
    for x in L:
        D = V.conj().T @ weyl(x) @ V
        assert np.allclose(D, np.diag(np.diag(D)), atol=1e-9)
    return V


def check_qudit(Ls):
    for x in ELEMENTS:
        for y in ELEMENTS:
            comm = np.allclose(weyl(x) @ weyl(y), weyl(y) @ weyl(x))
            assert comm == (form(x, y) == 0)
    print('qudit dictionary: line (a2,b2,a3,b3) -> Z^a2 X^b2 (x) Z3^a3 X3^b3; mutual locality = commutation (36 x 36 pairs)')
    bases = [basis_of(L) for L in Ls]
    print('  each of the 12 global forms is a stabilizer basis of one 6-level qudit (6 joint eigenstates, all 1-dim)')
    S = np.kron(np.diag([1, 1j]), np.eye(3))
    P3 = np.kron(np.eye(2), np.diag([1, 1, np.exp(2j * np.pi / 3)]))
    idx = {L: i for i, L in enumerate(Ls)}
    for L, B in zip(Ls, bases):
        for gate, T in ((S, theta2), (P3, theta3)):
            target = bases[idx[frozenset(T(x) for x in L)]]
            assert same_ray_sets(gate @ B, target)
    print('  theta shifts are Clifford gates: S on the qubit, diag(1,1,omega) on the qutrit')
    H2 = np.array([[1, 1], [1, -1]]) / np.sqrt(2)
    F3 = np.array([[np.exp(2j * np.pi * j * k / 3) for k in range(3)] for j in range(3)]) / np.sqrt(3)
    Fg = np.kron(H2, F3)
    simply = [L for L in Ls if gamma_order(L) == 1][0]
    z6_0 = frozenset((0, b2, 0, b3) for b2 in range(2) for b3 in range(3))
    assert z6_0 in idx
    assert same_ray_sets(Fg @ bases[idx[simply]], bases[idx[z6_0]])
    print('  the electric-magnetic Fourier gate maps the Gamma=1 basis to the Gamma=Z6 (theta=0) basis')
    return bases


def gname(q):
    return '1' if q == 1 else f'Z{q}'


def check_torus(Ls):
    simply = [L for L in Ls if gamma_order(L) == 1][0]
    names = {}
    for L in Ls:
        names.setdefault(gamma_order(L), L)
    print('spatial 3-torus: every global form has 6^3 = 216 flux sectors; sectors shared with Gamma=1:')
    for q in (1, 2, 3, 6):
        L = names[q]
        shared = len(L & simply) ** 3
        assert shared == (6 // q) ** 3
        print(f'  Gamma = {gname(q)}: {len(L) ** 3} sectors, {shared} shared with SU(3)xSU(2)xU(1)')


def check_echo_link(Ls):
    """Distinct Z6 echo phases on the links allowed in the Gamma = Z_q theory."""
    print('centre echo on the link classes each global form allows:')
    for q in (1, 2, 3, 6):
        L = [L for L in Ls if gamma_order(L) == q][0]
        allowed = [(x[0], x[2]) for x in L if x[1] == 0 and x[3] == 0]
        # Z6 echo generator acts on electric class (a2, a3) by exp(2 pi i (a2/2 + a3/3)) (up to sign of a3)
        phases = {(F(a2, 2) - F(a3, 3)) % 1 for a2, a3 in allowed}
        assert len(phases) == 6 // q
        print(f'  Gamma = {gname(q)}: allowed link classes {sorted(allowed)}; Z6 echo takes {len(phases)} distinct value(s)')
    print('  -> the centre echo of SM_CENTER_ECHO.md is a working control only if Gamma != Z6')


def main():
    check_reduction()
    Ls = check_lagrangians()
    check_qudit(Ls)
    check_torus(Ls)
    check_echo_link(Ls)
    print('ALL GLOBAL-FORM CODE CHECKS PASSED')


if __name__ == '__main__':
    main()
