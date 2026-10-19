#!/usr/bin/env python3
"""Confinement as a code: colour-singlet subspaces of n quarks (see COLOUR_SINGLET_CODE.md).

Each quark's colour is a qutrit.  The colour-singlet subspace of n quarks is the code.

 1. Singlet spaces: n = 3, 6, 9 give dimensions 1, 5, 42 (standard Young tableaux of
    shape (n/3, n/3, n/3)); built from products of epsilon tensors, checked SU(3)-invariant.
 2. Detection: every single-quark operator E acts on singlets as (tr E / 3) * identity, so the
    singlet space detects every single-quark colour error (distance >= 2).
 3. Obstruction: on the full singlet space, quark exchange P_ij is not a scalar (the space is an
    irreducible S_n module), so it does not correct single errors: distance exactly 2.
 4. n = 6: numerical search finds no 2-dimensional subcode correcting single-quark errors.
 5. n = 9: the shape (3,3,3) is self-conjugate, so there is an operator J with
    J P_ij = -P_ij J for every exchange.  J is built two ways: as the one-dimensional solution
    of the anticommutation equations, and as the Pauli filled shell of 9 quarks with 3 colours
    and 3 flavours.  J^2 = -1 and its +i eigenspace K (dimension 21) has <u|P_ij|v> = 0 for all
    u, v in K and all 36 exchanges.
 6. Full Knill-Laflamme check of K in the 3^9 = 19683-dimensional space for every one- and
    two-site operator (all 81 two-site operator directions at once, via reduced overlap tensors):
    K corrects every single-quark error and every pair of erased quarks.  Distance exactly 3.
 7. Every colour-exchange interaction sum_ij c_ij T_i.T_j compresses to a constant on K.
Requires numpy and scipy.  Ends with ALL COLOUR-SINGLET CODE CHECKS PASSED.
"""
import itertools
import math
import numpy as np
from scipy.optimize import minimize

OMEGA = np.exp(2j * np.pi / 3)


def epsilon():
    E = np.zeros((3, 3, 3))
    for p in itertools.permutations(range(3)):
        E[p] = round(np.linalg.det(np.eye(3)[list(p)]))
    return E


def triple_partitions(sites):
    if not sites:
        yield []
        return
    a = sites[0]
    for b, c in itertools.combinations(sites[1:], 2):
        rest = [s for s in sites if s not in (a, b, c)]
        for r in triple_partitions(rest):
            yield [(a, b, c)] + r


def singlet_basis(n):
    """Orthonormal real basis (3^n x D) of the colour-singlet space."""
    E = epsilon()
    vecs = []
    for P in triple_partitions(list(range(n))):
        T = E
        for _ in P[1:]:
            T = np.tensordot(T, E, axes=0)
        order = [s for t in P for s in t]
        vecs.append(np.transpose(T, np.argsort(order)).reshape(-1))
    A = np.array(vecs).T
    U, s, _ = np.linalg.svd(A, full_matrices=False)
    return U[:, : int((s > 1e-9).sum())]


def gellmann():
    L = np.zeros((8, 3, 3), dtype=complex)
    k = 0
    for i, j in itertools.combinations(range(3), 2):
        L[k, i, j] = L[k, j, i] = 1; k += 1
        L[k, i, j] = -1j; L[k, j, i] = 1j; k += 1
    L[6] = np.diag([1, -1, 0]); L[7] = np.diag([1, 1, -2]) / np.sqrt(3)
    return L / 2


def apply_site(vecs, n, site, M):
    """Apply 3x3 matrix M on one site to columns of vecs (3^n x m)."""
    m = vecs.shape[1]
    T = vecs.T.reshape((m,) + (3,) * n)
    T = np.moveaxis(np.tensordot(M, T, axes=([1], [site + 1])), 0, site + 1)
    return T.reshape(m, -1).T


def swap_matrix(B, n, i, j):
    D = B.shape[1]
    T = B.T.reshape((D,) + (3,) * n)
    ax = list(range(n + 1)); ax[i + 1], ax[j + 1] = ax[j + 1], ax[i + 1]
    return B.T @ np.transpose(T, ax).reshape(D, -1).T


def syt_count(shape):
    n = sum(shape)
    hooks = 1
    conj = [sum(1 for r in shape if r > c) for c in range(shape[0])]
    for i, r in enumerate(shape):
        for j in range(r):
            hooks *= (r - j - 1) + (conj[j] - i - 1) + 1
    return math.factorial(n) // hooks


def check_spaces():
    bases = {}
    for n, D in ((3, 1), (6, 5), (9, 42)):
        B = singlet_basis(n)
        assert B.shape[1] == D == syt_count((n // 3,) * 3)
        bases[n] = B
    L = gellmann()
    for n in (6, 9):
        B = bases[n]
        for a in range(8):
            tot = sum(apply_site(B, n, s, L[a]) for s in range(n))
            assert np.abs(tot).max() < 1e-10
    print('colour-singlet dimensions: 3 quarks -> 1, 6 -> 5, 9 -> 42 (all SU(3)-invariant)')
    return bases


def check_detection(bases):
    rng = np.random.default_rng(3)
    for n in (6, 9):
        B = bases[n]
        for site in range(n):
            M = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3))
            C = B.conj().T @ apply_site(B, n, site, M)
            assert np.abs(C - np.trace(M) / 3 * np.eye(B.shape[1])).max() < 1e-10
        P = swap_matrix(B, n, 0, 1)
        evals = np.linalg.eigvalsh(P)
        assert evals.min() < -0.5 and evals.max() > 0.5
        print(f'n={n}: every single-quark operator acts as (tr E/3) 1 on singlets -> detected; '
              f'P_01 has eigenvalues {sorted(set(int(round(e)) for e in evals))} -> distance exactly 2')


def traceless_span(Ms):
    D = Ms[0].shape[0]
    V = np.array([(M - np.trace(M) / D * np.eye(D)).reshape(-1) for M in Ms])
    _, s, Vt = np.linalg.svd(V, full_matrices=False)
    A = Vt[: int((s > 1e-9).sum())].reshape(-1, D, D)
    return (A + A.transpose(0, 2, 1)) / 2


def kl_residual_search(A, k, trials, seed):
    """Minimise sum ||X^dag A X - scalar||^2 + ||X^dag X - 1||^2 over complex D x k matrices X."""
    D = A.shape[1]

    def fg(x):
        X = (x[: D * k] + 1j * x[D * k:]).reshape(D, k)
        Pn = X.conj().T @ X - np.eye(k)
        AX = np.einsum('rij,jb->rib', A, X)
        C = np.einsum('ia,rib->rab', X.conj(), AX)
        C = C - (np.einsum('raa->r', C) / k)[:, None, None] * np.eye(k)
        f = np.sum(np.abs(C) ** 2) + np.sum(np.abs(Pn) ** 2)
        g = 2 * np.einsum('rib,rba->ia', AX, C) + 2 * X @ Pn
        return f, np.concatenate([2 * g.real.ravel(), 2 * g.imag.ravel()])
    rng = np.random.default_rng(seed)
    best = np.inf
    for _ in range(trials):
        res = minimize(fg, rng.normal(size=2 * D * k) / np.sqrt(D), jac=True, method='L-BFGS-B',
                       options={'maxiter': 20000, 'gtol': 1e-15, 'ftol': 1e-20, 'maxcor': 50})
        best = min(best, res.fun)
    return best


def check_six(bases):
    B = bases[6]
    Ms = [swap_matrix(B, 6, i, j) for i, j in itertools.combinations(range(6), 2)]
    A = traceless_span(Ms)
    best = kl_residual_search(A, 2, 20, 0)
    assert best > 0.5
    print(f'n=6: exchange span has dimension {len(A)}; best residual for a 2-dim correcting subcode '
          f'over 20 random starts = {best:.3f} (none found)')


def j_from_anticommutation(B):
    """Solve J s_k + s_k J = 0 for the 8 adjacent exchanges; return the unique solution."""
    D = B.shape[1]
    I = np.eye(D)
    rows = []
    for k in range(8):
        M = swap_matrix(B, 9, k, k + 1)
        rows.append(np.kron(M, I) + np.kron(I, M.T))     # vec(M J + J M), row-major
    _, s, Vt = np.linalg.svd(np.vstack(rows))
    null = int((s < 1e-9).sum())
    return null, Vt[-1].reshape(D, D)


def j_from_filled_shell(B):
    """J_mn = <b_m (colour) x b_n (flavour) | 9 quarks filling 3 colours x 3 flavours>."""
    perms = np.array(list(itertools.permutations(range(9))), dtype=np.int64)
    inv = 0
    for i in range(9):
        for j in range(i + 1, 9):
            inv = inv + (perms[:, i] > perms[:, j])
    sign = np.where(inv % 2 == 0, 1.0, -1.0)
    pw = 3 ** np.arange(8, -1, -1)
    cidx = (perms // 3) @ pw
    fidx = (perms % 3) @ pw
    J = (B[cidx, :] * sign[:, None]).T @ B[fidx, :]
    return J / math.sqrt(math.factorial(9))


def check_nine(bases):
    B = bases[9]
    D = B.shape[1]
    null, J1 = j_from_anticommutation(B)
    assert null == 1
    J2 = j_from_filled_shell(B)
    c = np.vdot(J1.ravel(), J2.ravel()) / np.vdot(J1.ravel(), J1.ravel())
    assert np.abs(J2 - c * J1).max() < 1e-10 and abs(c) > 1e-6
    J = J2 / np.sqrt(np.abs(np.trace(J2.T @ J2)) / D)
    assert np.abs(J.T @ J - np.eye(D)).max() < 1e-10            # orthogonal
    assert np.abs(J @ J + np.eye(D)).max() < 1e-10              # J^2 = -1
    assert np.abs(J + J.T).max() < 1e-10                        # so J is antisymmetric
    print('n=9: anticommutation equations have a 1-dim solution space; it equals the filled-shell')
    print('     operator of 9 quarks in 3 colours x 3 flavours; J is orthogonal, antisymmetric, J^2 = -1')
    w, V = np.linalg.eig(J)
    K = V[:, np.abs(w - 1j) < 1e-8]
    K, _ = np.linalg.qr(K)
    assert K.shape[1] == D // 2 == 21
    Ms = [swap_matrix(B, 9, i, j) for i, j in itertools.combinations(range(9), 2)]
    for M in Ms:
        assert np.abs(J @ M + M @ J).max() < 1e-10
        assert np.abs(K.conj().T @ M @ K).max() < 1e-10
    print(f'     +i eigenspace K has dimension {K.shape[1]}; all 36 exchanges have zero matrix elements on K')
    # the subcode is maximal among those our search finds
    A = traceless_span(Ms)
    best22 = kl_residual_search(A, 22, 3, 1)
    assert best22 > 1e-3
    print(f'     search for a 22-dim correcting subcode: best residual {best22:.3f} (none found)')
    # three-site invariant operator that K does not detect: a 3-cycle
    T = B.T.reshape((D,) + (3,) * 9)
    cyc = np.transpose(T, [0, 2, 3, 1] + list(range(4, 10))).reshape(D, -1).T
    C3 = K.conj().T @ (B.T @ cyc) @ K
    assert np.abs(C3 - np.trace(C3) / 21 * np.eye(21)).max() > 0.1
    print('     a 3-cycle of quarks acts non-trivially on K: distance exactly 3')
    return B @ K


def check_full_kl(code):
    """KL in the full 3^9 space for every one- and two-site operator."""
    n = 9
    m = code.shape[1]
    T = code.T.reshape((m,) + (3,) * n)
    worst = 0.0
    for i in range(n):
        X = np.moveaxis(T, i + 1, 1).reshape(m, 3, -1)
        G = np.einsum('mxr,nyr->mnxy', X.conj(), X)
        for x in range(3):
            for y in range(3):
                g = G[:, :, x, y]
                worst = max(worst, np.abs(g - np.trace(g) / m * np.eye(m)).max())
    for i, j in itertools.combinations(range(n), 2):
        X = np.moveaxis(T, (i + 1, j + 1), (1, 2)).reshape(m, 9, -1)
        G = np.einsum('mxr,nyr->mnxy', X.conj(), X)
        dev = G - np.einsum('mmxy->xy', G)[None, None] / m * np.eye(m)[:, :, None, None]
        worst = max(worst, np.abs(dev).max())
    assert worst < 1e-10
    print(f'full KL in C^(3^9): all 1- and 2-site operators (81 directions per pair, 36 pairs) '
          f'compress to scalars on the 21-dim code; max deviation {worst:.1e}')
    print('     -> corrects any single-quark colour error, or any two erased quarks')
    # colour-exchange interaction
    rng = np.random.default_rng(11)
    L = gellmann()
    v = code[:, :4]
    H = np.zeros_like(v)
    for i, j in itertools.combinations(range(n), 2):
        cij = rng.normal()
        for a in range(8):
            H += cij * apply_site(apply_site(v, n, j, L[a]), n, i, L[a])
    C = v.conj().T @ H
    assert np.abs(C - np.trace(C) / 4 * np.eye(4)).max() < 1e-10
    print('     random colour-exchange Hamiltonian sum c_ij T_i.T_j is a constant on the code')


def check_uniqueness():
    print('self-conjugate singlet shapes (m^N) with m = N: n = N^2 particles of SU(N)')
    for N in (2, 3, 4):
        D = syt_count((N,) * N)
        print(f'  SU({N}): {N * N} particles, singlet dimension {D}, half = {D // 2}')
    assert syt_count((2, 2)) == 2 and syt_count((3, 3, 3)) == 42
    for m in (1, 2, 4, 5):
        shape = (m,) * 3
        conj = tuple(sum(1 for r in shape if r > c) for c in range(m))
        assert conj != shape
    print('  for SU(3), n = 9 is the only quark number whose singlet shape (n/3)^3 is self-conjugate')


def main():
    bases = check_spaces()
    check_detection(bases)
    check_six(bases)
    code = check_nine(bases)
    check_full_kl(code)
    check_uniqueness()
    print('ALL COLOUR-SINGLET CODE CHECKS PASSED')


if __name__ == '__main__':
    main()
