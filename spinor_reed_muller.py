#!/usr/bin/env python3
"""Fermionic 16: the binary shadow of one generation and its Majorana version
(see FERMIONIC_GENERATION_CODE.md).

One generation = the 16 spinor weights of SO(10), sign vectors s in {+1,-1}^5 with an even
number of minus signs.  Over odd primes the charges +-1 give the [[15,1,3]]_p code of the main
paper.  Mod 2, +-1 are both 1 and the lattice collapses, so write each weight as bits
b = (1 - s)/2 instead.  Checks:

 1. For SO(2r), r = 3..6, the bit vectors of the 2^(r-1) spinor weights are the even-weight
    subspace of F2^r, and the affine functions on them form RM(1, r-1).  For SO(10):
    weight enumerator 1 + 30 z^8 + z^16, and nu^c is the origin b = 0.
 2. Puncturing nu^c gives the quantum Reed-Muller codes [[3,1,1]], [[7,1,3]], [[15,1,3]],
    [[31,1,3]] (SO(6), SO(8), SO(10), SO(12)), the same X/Z check counts as the paper's code.
 3. Transversal R_k = diag(1, e^{2 pi i/2^k}) acts logically up to k = r - 2: S for SO(8),
    T for SO(10), sqrt(T) for SO(12).  SO(10) is the smallest SO(2r) whose generation code has
    a transversal non-Clifford gate.  State-vector check for SO(10): T^{x15} = logical T^dag,
    and Knill-Laflamme for all weight-1 Pauli errors on the 2^15-dimensional space.
 4. Fermionic version (16 complex modes = 32 Majoranas).  X-type checks are products of 8
    Majoranas, Z-type checks are mode parities; all commute.  Without nu^c, the logical X is a
    product of 15 Majoranas (fermion-odd) and the logical Z is the fermion parity (-1)^F:
    superselection forbids the logical qubit.  The spin-Z4 generator i^N acts as logical S^dag.
 5. With nu^c as the 16th mode the code has k = 0: a unique state, every branch with N in
    {0, 8, 16}.  It is (|0_L>|0>_nu + |1_L>|1>_nu)/sqrt2, fermion-even and invariant under
    i^N and e^{i pi N/4}, but not e^{i pi N/8}.
 6. Fidkowski-Kitaev-type count: a stabilizer state of 2n Majoranas with only real quartic
    (time-reversal-even) stabilizers needs a doubly-even self-dual code, which exists iff
    2n = 0 mod 8.  One generation (32 Majoranas) passes with RM(2,5); fifteen Weyl fields
    (30 Majoranas) fail.  This is a mod-8 statement, not the mod-16 Dai-Freed anomaly.
Requires numpy.  Ends with ALL FERMIONIC-GENERATION CHECKS PASSED.
"""
import itertools
import numpy as np
from charge_lattice_tools import W16, Q0, single


# ------------------------------------------------------------------- F2 helpers
def rank2(M):
    M = [int(''.join(map(str, row)), 2) for row in M]
    r = 0
    rows = M[:]
    while rows:
        piv = max(rows)
        if piv == 0:
            break
        rows.remove(piv)
        r += 1
        top = piv.bit_length() - 1
        rows = [x ^ piv if (x >> top) & 1 else x for x in rows]
    return r


def span2(rows):
    rows = [np.array(r, dtype=np.int64) % 2 for r in rows]
    out = {tuple([0] * len(rows[0]))}
    for r in rows:
        out |= {tuple((np.array(v) + r) % 2) for v in out}
    return out


def spinor_bits(r):
    pts = [s for s in itertools.product((1, -1), repeat=r) if s.count(-1) % 2 == 0]
    bits = [tuple((1 - x) // 2 for x in s) for s in pts]
    return pts, bits


def rm1(m):
    """RM(1,m) evaluation rows on the points of F2^m."""
    pts = list(itertools.product((0, 1), repeat=m))
    return [[1] * len(pts)] + [[p[i] for p in pts] for i in range(m)], pts


def check_shadow():
    for r in range(3, 7):
        pts, bits = spinor_bits(r)
        assert len(bits) == 2 ** (r - 1) and all(sum(b) % 2 == 0 for b in bits)
        G = [[1] * len(bits)] + [[b[i] for b in bits] for i in range(r)]
        assert rank2(G) == r                                  # 1 + (r - 1) independent rows
        # identify with RM(1, r-1) via dropping the last bit (a linear iso E -> F2^{r-1})
        R, rpts = rm1(r - 1)
        order = [rpts.index(b[:-1]) for b in bits]
        Rp = [[row[j] for j in order] for row in R]
        assert span2(G) == span2(Rp)
        if r == 5:
            assert set(pts) == set(W16) and bits[pts.index(tuple(Q0))] == (0,) * 5
            wts = sorted(sum(c) for c in span2(G))
            enum = {w: wts.count(w) for w in set(wts)}
            assert enum == {0: 1, 8: 30, 16: 1}
            print('SO(10): bits of the 16 spinor weights carry RM(1,4); weight enumerator 1 + 30 z^8 + z^16;'
                  ' nu^c is the origin')
    print('SO(2r), r = 3..6: affine functions on the 2^(r-1) spinor weights = RM(1, r-1)')


def punctured_code(r):
    pts, bits = spinor_bits(r)
    keep = [i for i, b in enumerate(bits) if any(b)]
    lin = [[bits[i][j] for i in keep] for j in range(r)]   # X checks (a_0 = 0)
    n = len(keep)
    return n, lin


def distances(n, lin):
    X = span2(lin)
    ones = np.ones(n, dtype=np.int64)
    dx = min(sum((np.array(v) + ones) % 2) for v in X)
    L = np.array(lin)
    dz = None
    for w in range(1, n + 1):
        for supp in itertools.combinations(range(n), w):
            z = np.zeros(n, dtype=np.int64); z[list(supp)] = 1
            if not (L @ z % 2).any() and w % 2 == 1:
                dz = w
                break
        if dz:
            break
    return dx, dz


def transversal_level(n, lin):
    X = span2(lin)
    coset0 = [sum(v) for v in X]
    coset1 = [n - sum(v) for v in X]       # complement = v + all-ones
    best = 0
    for k in range(1, 8):
        mod = 2 ** k
        if len({w % mod for w in coset0}) == 1 and len({w % mod for w in coset1}) == 1:
            best = k
            phase = (coset1[0] - coset0[0]) % mod
    return best, phase


def check_quantum_rm():
    expect = {3: (3, 1, 1, 1), 4: (7, 3, 3, 2), 5: (15, 7, 3, 3), 6: (31, 15, 3, 4)}
    names = {3: 'SO(6)', 4: 'SO(8)', 5: 'SO(10)', 6: 'SO(12)'}
    gates = {1: 'Z (Pauli)', 2: 'S (Clifford)', 3: 'T (non-Clifford)', 4: 'sqrt T (non-Clifford)'}
    for r in range(3, 7):
        n, lin = punctured_code(r)
        rx = rank2(lin)
        rz = n - rx - 1
        dx, dz = distances(n, lin)
        k, phase = transversal_level(n, lin)
        assert (n, dx, min(dx, dz), k) == expect[r], (r, n, dx, dz, k)
        assert phase == 2 ** k - 1                               # logical R_k^dagger
        print(f'{names[r]:6s}: [[{n},1,{min(dx, dz)}]]  X checks {rx}, Z checks {rz}, X distance {dx}; '
              f'transversal R_{k} = {gates[k]} acts as logical R_{k}^dag')
    c = single(3).params()
    assert (c['n'], c['k'], c['X_stab'], c['Z_stab']) == (15, 1, 4, 10)
    print("  same n, k and check counts as the paper's [[15,1,3]]_p for odd p")


def kets(r=5):
    n, lin = punctured_code(r)
    X = span2(lin)
    pw = 2 ** np.arange(n - 1, -1, -1)
    z0 = np.zeros(2 ** n); z1 = np.zeros(2 ** n)
    for v in X:
        v = np.array(v)
        z0[int(v @ pw)] = 1
        z1[int(((v + 1) % 2) @ pw)] = 1
    return n, lin, z0 / np.linalg.norm(z0), z1 / np.linalg.norm(z1)


def check_statevector():
    n, lin, z0, z1 = kets()
    N = 2 ** n
    weights = np.array([bin(i).count('1') for i in range(N)])
    T = np.exp(1j * np.pi / 4 * weights)
    assert np.allclose(T * z0, z0) and np.allclose(T * z1, np.exp(-1j * np.pi / 4) * z1)
    idx = np.arange(N)

    def pauli(v, site, kind):
        bit = 1 << (n - 1 - site)
        if kind == 'X':
            return v[idx ^ bit]
        if kind == 'Z':
            return v * np.where(idx & bit, -1, 1)
        return 1j * pauli(pauli(v, site, 'Z'), site, 'X')
    ops = [(None, None)] + [(s, k) for s in range(n) for k in 'XYZ']
    code = [z0, z1]
    applied = [[v if s is None else pauli(v, s, k) for v in code] for s, k in ops]
    worst = 0.0
    for a in range(len(ops)):
        for b in range(a, len(ops)):
            G = np.array([[np.vdot(applied[a][i], applied[b][j]) for j in range(2)] for i in range(2)])
            worst = max(worst, abs(G[0, 1]), abs(G[1, 0]), abs(G[0, 0] - G[1, 1]))
    assert worst < 1e-12
    print(f'state vectors (2^15): T^x15 = logical T^dag; KL holds for all {len(ops)} weight<=1 Paulis'
          f' ({len(ops) * (len(ops) + 1) // 2} pairs)')
    return n, lin, z0, z1, weights


# ------------------------------------------------------------------- Majorana layer
def majorana_commute(A, B):
    """Products of Majoranas on supports A, B (sets of indices) commute iff |A||B| - |A&B| is even."""
    return (len(A) * len(B) - len(A & B)) % 2 == 0


def check_majorana_rule():
    """Validate the commutation rule with explicit Jordan-Wigner matrices on 3 modes."""
    n = 3
    X = np.array([[0, 1], [1, 0]]); Y = np.array([[0, -1j], [1j, 0]]); Z = np.diag([1, -1])

    def jw(j, P):
        mats = [Z] * j + [P] + [np.eye(2)] * (n - j - 1)
        out = mats[0]
        for m in mats[1:]:
            out = np.kron(out, m)
        return out
    gam = [jw(j, P) for j in range(n) for P in (X, Y)]
    for A in itertools.chain.from_iterable(itertools.combinations(range(6), w) for w in range(1, 7)):
        for B in itertools.chain.from_iterable(itertools.combinations(range(6), w) for w in range(1, 4)):
            MA = np.eye(8); MB = np.eye(8)
            for i in A:
                MA = MA @ gam[i]
            for i in B:
                MB = MB @ gam[i]
            assert np.allclose(MA @ MB, MB @ MA) == majorana_commute(set(A), set(B))
    print('Majorana commutation rule |A||B| - |A cap B| even, validated on all supports of 6 Majoranas')


def check_fermionic(n, lin, z0, z1, weights):
    # mode j has Majoranas 2j (gamma) and 2j+1 (gamma'); X-type check on S -> {2j : j in S};
    # Z-type (parity) check on T -> {2j, 2j+1 : j in T}
    Xs = [frozenset(2 * j for j in range(n) if row[j]) for row in lin]
    L = np.array(lin)
    zbasis = []
    for v in itertools.product((0, 1), repeat=n):
        v = np.array(v)
        if not (L @ v % 2).any() and v.sum() % 2 == 0:
            zbasis.append(v)
    Zs = [frozenset(x for j in range(n) if z[j] for x in (2 * j, 2 * j + 1)) for z in zbasis]
    checks = Xs + Zs
    for A, B in itertools.combinations(checks, 2):
        assert majorana_commute(A, B)
    assert all(len(A) % 2 == 0 for A in checks)
    logical_X = frozenset(2 * j for j in range(n))
    assert len(logical_X) == 15 and all(majorana_commute(logical_X, A) for A in checks)
    zr = rank2([list(z) for z in zbasis])
    assert rank2(lin) == 4 and zr == 10
    print(f'15-mode fermionic code: X-type checks of rank 4 (8 Majoranas each), parity checks of rank {zr};'
          f' 5 X generators and all 2^10 parity checks commute pairwise;')
    print('  logical X is a product of 15 Majoranas: fermion-odd')
    parity = np.where(weights % 2, -1, 1)
    assert np.allclose(parity * z0, z0) and np.allclose(parity * z1, -z1)
    Z4 = 1j ** weights
    assert np.allclose(Z4 * z0, z0) and np.allclose(Z4 * z1, -1j * z1)
    print('  logical Z = (-1)^F, so a logical superposition mixes fermion parities (forbidden);'
          ' spin-Z4 i^N = logical S^dag')
    # add nu^c as mode 16: state (|0_L>|0> + |1_L>|1>)/sqrt2
    psi = np.concatenate([z0, z1]) / np.sqrt(2)     # nu^c is the leading (most significant) mode
    N16 = np.concatenate([weights, weights + 1])
    supp = np.abs(psi) > 1e-12
    assert set(N16[supp]) == {0, 8, 16}
    for phase, invariant in ((1j, True), (np.exp(1j * np.pi / 4), True), (np.exp(1j * np.pi / 8), False)):
        assert np.allclose(phase ** N16 * psi, psi) == invariant
    # it is the unique state of the unpunctured RM(1,4) stabilizer code on 16 modes
    pts, bits = spinor_bits(5)
    G = [[1] * 16] + [[b[i] for b in bits] for i in range(5)]
    words = span2(G)
    order = [bits.index((0,) * 5)] + [i for i, b in enumerate(bits) if any(b)]
    pw = 2 ** np.arange(15, -1, -1)
    ref = np.zeros(2 ** 16)
    for w in words:
        ref[int(np.array([w[i] for i in order]) @ pw)] = 1
    ref /= np.linalg.norm(ref)
    assert np.isclose(abs(np.vdot(ref, psi)), 1)
    assert rank2(G) == 5          # X checks RM(1,4); Z checks = its dual RM(2,4), rank 11; 5 + 11 = 16 -> k = 0
    print('16 modes (nu^c added): the unpunctured RM(1,4) code has k = 0; its unique state is')
    print('  (|0_L>|0>_nu + |1_L>|1>_nu)/sqrt2, with N in {0, 8, 16}: fermion-even, invariant under i^N and')
    print('  e^{i pi N/4}, not under e^{i pi N/8}.  nu^c is the reference mode that makes the qubit physical.')


# ------------------------------------------------------------------- Fidkowski-Kitaev-type count
def max_doubly_even_self_orthogonal(length):
    words = [w for w in itertools.product((0, 1), repeat=length) if sum(w) % 4 == 0 and any(w)]
    words = [int(''.join(map(str, w)), 2) for w in words]
    best = 0

    def ok(a, b):
        return bin(a & b).count('1') % 2 == 0

    def grow(basis, code, cands):
        nonlocal best
        best = max(best, len(basis))
        if best == length // 2:
            return
        for i, w in enumerate(cands):
            if w in code:
                continue
            if all(ok(w, c) for c in basis):
                new = code | {c ^ w for c in code}
                if all(bin(c).count('1') % 4 == 0 for c in new):
                    grow(basis + [w], new, cands[i + 1:])
                    if best == length // 2:
                        return
    grow([], {0}, words)
    return best


def rm(rdeg, m):
    pts = list(itertools.product((0, 1), repeat=m))
    rows = []
    for d in range(rdeg + 1):
        for S in itertools.combinations(range(m), d):
            rows.append([int(all(p[i] for i in S)) for p in pts])
    return rows


def check_fk():
    print('real quartic Majorana stabilizer states (doubly-even self-dual codes) on 2n Majoranas:')
    for L in (2, 4, 6, 8):
        d = max_doubly_even_self_orthogonal(L)
        exists = d == L // 2
        assert exists == (L % 8 == 0), L
        print(f'  2n = {L:2d}: max doubly-even self-orthogonal dimension {d} of {L // 2} -> '
              f'{"exists" if exists else "none"}')
    print('  (all lengths: Type II codes exist iff the length is 0 mod 8, Gleason / MacWilliams-Sloane ch. 19)')
    G = rm(2, 5)
    assert rank2(G) == 16
    for a, b in itertools.combinations_with_replacement(G, 2):
        assert sum(x & y for x, y in zip(a, b)) % 2 == 0
    assert all(sum(row) % 4 == 0 for row in G)
    print('  one generation, 32 Majoranas: RM(2,5) is doubly-even self-dual -> a T-even gapped stabilizer state exists')
    print('  fifteen Weyl fields, 30 Majoranas: 30 != 0 mod 8 -> none.  (mod 8 in Majoranas = mod 4 in Weyl fields)')


def main():
    check_shadow()
    check_quantum_rm()
    n, lin, z0, z1, weights = check_statevector()
    check_majorana_rule()
    check_fermionic(n, lin, z0, z1, weights)
    check_fk()
    print('ALL FERMIONIC-GENERATION CHECKS PASSED')


if __name__ == '__main__':
    main()
