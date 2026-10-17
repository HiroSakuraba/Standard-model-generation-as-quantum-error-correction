#!/usr/bin/env python3
"""Follow-up computations for the five open steps (see OPEN_STEPS_NOTE.md).

    python3 open_steps.py            # full run (a few minutes)
    python3 open_steps.py --quick    # skips the p = 13 nonlinear-row search

Requires Python 3.9+ and numpy.  All code-theoretic arithmetic is exact (integers mod p,
Fractions); numpy floating point is used only for the illustrative Yukawa-texture rank and
mixing tests in Step 1 (random generic entries) and for the measured-data breaking table.
Every reported number is asserted, so a clean run ends with ALL OPEN-STEP CHECKS PASSED.
"""
from __future__ import annotations
import sys, time, itertools
import numpy as np
from fractions import Fraction
from charge_lattice_tools import *

QUICK = '--quick' in sys.argv
T0 = time.time()
BSM = ((1, 1, 1, 0, 0), (0, 0, 0, 1, 1))          # family B-L (S3) and right-handed isospin (S2) rows
IDX = {q: i for i, q in enumerate(W16)}
NAME_IDX = {nm: i for i, nm in enumerate(NAMES1)}
WM = np.array(W16).T


def say(*a):
    print(*a, f'[{time.time() - T0:.0f}s]' if a and str(a[0]).startswith('==') else '')
    sys.stdout.flush()


def fam_code(nf, p, extra):
    """common SO(10) rows + per-family extra rows extra(q, g) (list)."""
    return build([(lambda g: (lambda q: list(q) + extra(q, g)))(g) for g in range(nf)], p)


# =====================================================================================  STEP 2
def step2():
    say('== STEP 2: exact logical channels under i.i.d. Z noise (distillation post-selection)')
    nz = lambda cc: sum(1 for x in cc if x)
    base = {p: class_weight_distributions(single(p)) for p in (3, 5)}

    def series(D, n, p, pred, order=7):
        A = [sum(D[cc][w] for cc in D if pred(cc)) for w in range(n + 1)]
        return A, to_eps(A, p, n, order)
    rows = {}
    for nm, p, nf, c in (('coupled B3 2fam', 3, 2, coupled(B3, 2, 3)), ('direct 2fam', 3, 2, None),
                         ('coupled B3 3fam', 3, 3, coupled(B3, 3, 3)), ('direct 3fam', 3, 3, None),
                         ('coupled B5 2fam', 5, 2, coupled(B5, 2, 5)), ('direct 2fam', 5, 2, None),
                         ('coupled B5 3fam', 5, 3, coupled(B5, 3, 5)), ('direct 3fam', 5, 3, None)):
        n = 15 * nf
        if c is None:
            D = direct_sum_distributions(base[p], 15, nf)
            D = {cc: (v + [0] * (n + 1))[:n + 1] for cc, v in D.items()}
        else:
            D = class_weight_distributions(c)
        acc, Sacc = series(D, n, p, lambda cc: True)
        one, Sone = series(D, n, p, lambda cc: nz(cc) == 1)
        mul, Smul = series(D, n, p, lambda cc: nz(cc) >= 2)
        e = Fraction(1, 1000)
        pa = exact_value(acc, p, n, e)
        per_out = sum(exact_value([sum(D[cc][w] for cc in D if cc[i]) for w in range(n + 1)], p, n, e)
                      for i in range(nf)) / pa / nf
        corr = exact_value(mul, p, n, e) / pa
        rows[(nm, p)] = (one[:7], mul[:7], Sone[3], Smul[4], Smul[6], float(per_out), float(corr), Sacc[:5])
        say(f'  p={p} {nm:16s} single-logical A_w {one[:7]}  correlated A_w {mul[:7]}')
        say(f'           P(single) = {Sone[3]} e^3 + {Sone[4]} e^4 ...   P(correlated) = {Smul[4]} e^4 + {Smul[5]} e^5 + {Smul[6]} e^6 ...')
        say(f'           at eps=1e-3: per-output infidelity {float(per_out):.6e}, P(correlated | accept) {float(corr):.3e}')
    # assertions: leading terms
    assert rows[('coupled B3 2fam', 3)][3] == 5 and rows[('coupled B5 2fam', 5)][3] == Fraction(9, 2)
    assert rows[('coupled B3 3fam', 3)][3] == 15 and rows[('coupled B5 3fam', 5)][3] == Fraction(27, 2)
    for nf in (2, 3):
        for p in (3, 5):
            assert rows[('direct %dfam' % nf, p)][3] == 0 and rows[('direct %dfam' % nf, p)][4] != 0
    assert rows[('coupled B3 2fam', 3)][2] == rows[('direct 2fam', 3)][2] == Fraction(25, 2)
    assert rows[('coupled B5 2fam', 5)][2] == rows[('direct 2fam', 5)][2] == Fraction(25, 8)
    # independent brute-force check of all weight <= 4 errors, with family split
    for nm, c, exp in (('coupled B3 2fam p3', coupled(B3, 2, 3), {('correlated', (2, 2)): 80, ('single', (2, 2)): 96, ('single', (4, 0)): 220}),
                       ('coupled B5 2fam p5', coupled(B5, 2, 5), {('correlated', (2, 2)): 1152, ('single', (2, 2)): 432, ('single', (4, 0)): 400})):
        D = class_weight_distributions(c)
        cnt, split = brute_low_weight(c, 4)
        for w in range(1, 5):
            for cc in D:
                assert cnt.get((w, cc), 0) == D[cc][w]
        agg = {}
        for (w, cc, fs), v in split.items():
            if w == 4:
                key = ('correlated' if nz(cc) >= 2 else 'single', fs)
                agg[key] = agg.get(key, 0) + v
        assert agg == exp, (nm, agg)
        say(f'  brute force (all errors of weight <= 4) agrees for {nm}; weight-4 logical patterns: {sorted(agg.items())}')
    return rows


# =====================================================================================  STEP 3
def step3():
    say('== STEP 3: equivalence to direct sums (stabilizer-group decomposability)')
    out = {}
    for nm, c in (('single p3', single(3)), ('direct 2fam p3', direct_sum(2, 3)), ('direct 2fam p5', direct_sum(2, 5)),
                  ('coupled B3 2fam p3', coupled(B3, 2, 3)), ('coupled B5 2fam p5', coupled(B5, 2, 5)),
                  ('coupled B5 2fam p7', coupled(B5, 2, 7)), ('coupled B3 3fam p3', coupled(B3, 3, 3)),
                  ('coupled B5 3fam p5', coupled(B5, 3, 5)), ('SM (B-L,T3R) 2fam p5', coupled(BSM, 2, 5)),
                  ('SM (B-L,T3R) 3fam p5', coupled(BSM, 3, 5))):
        comps = [len(x) for x in split_components(c)]
        out[nm] = comps
        sd = split_dims(c, list(range(15)))
        say(f'  {nm:22s} components {comps}   family cut: dim S = {sd["dimS"]}, dim S_A + dim S_B = {sd["dimSA"] + sd["dimSB"]}')
    assert out['direct 2fam p3'] == [15, 15] and out['single p3'] == [15]
    for nm in out:
        if nm.startswith(('coupled', 'SM')):
            assert len(out[nm]) == 1
    return out


def split_dims(code, A):
    p = code.p; B = [j for j in range(code.n) if j not in A]
    dSX = rank(code.SX, p); dSZ = code.n - code.r
    sxA = dSX - rank(code.SX[:, B], p); sxB = dSX - rank(code.SX[:, A], p)
    szA = len(A) - rank(code.GP[:, A], p); szB = len(B) - rank(code.GP[:, B], p)
    return dict(dimS=dSX + dSZ, dimSA=sxA + szA, dimSB=sxB + szB)


# =====================================================================================  STEP 1
MULT = {'Q': ['u1', 'u2', 'u3', 'd1', 'd2', 'd3'], 'uc': ['uc1', 'uc2', 'uc3'], 'dc': ['dc1', 'dc2', 'dc3'],
        'L': ['nu', 'e'], 'ec': ['ec'], 'nuc': ['nuc']}
S32 = {'Q': (1, 0), 'uc': (-1, 2), 'dc': (-1, -2), 'L': (-3, 0), 'ec': (3, -2), 'nuc': (3, 2)}


def sm_commuting(c):
    return all(len({dot(c, W16[NAME_IDX[x]]) for x in mem}) == 1 for mem in MULT.values())


def yukawa_scan(G, p, twoHDM, rng):
    def th(a, b, F, g): return (a[g] * S32[F][0] + b[g] * S32[F][1]) % p

    def mask(a, b, F1, F2, eta):
        return np.array([[(th(a, b, F1, g) + th(a, b, F2, h) + eta) % p == 0 for h in range(G)] for g in range(G)])

    def matching(M): return any(all(M[i, s[i]] for i in range(G)) for s in itertools.permutations(range(G)))
    def rnd(M): return np.where(M, rng.normal(size=M.shape) + 1j * rng.normal(size=M.shape), 0)
    def generic(U1, U2): return np.abs(U1.conj().T @ U2).min() > 1e-7
    surv = []
    for a in itertools.product(range(p), repeat=G):
        for b in itertools.product(range(p), repeat=G):
            for eu in range(p):
                for ed in (range(p) if twoHDM else [(-eu) % p]):
                    Mu, Md, Me = mask(a, b, 'Q', 'uc', eu), mask(a, b, 'Q', 'dc', ed), mask(a, b, 'L', 'ec', ed)
                    if not (matching(Mu) and matching(Md) and matching(Me)):
                        continue
                    Uu = np.linalg.svd(rnd(Mu))[0]; Ud = np.linalg.svd(rnd(Md))[0]; Ue = np.linalg.svd(rnd(Me))[0]
                    if not generic(Uu, Ud):
                        continue
                    opts = []
                    K = rnd(mask(a, b, 'L', 'L', 2 * eu)); opts.append(K + K.T)
                    Yn = rnd(mask(a, b, 'L', 'nuc', eu)); opts.append(Yn)
                    MR = rnd(mask(a, b, 'nuc', 'nuc', 0)); MR = MR + MR.T
                    if np.linalg.matrix_rank(MR) == G:
                        opts.append(Yn @ np.linalg.inv(MR) @ Yn.T)
                    for m in opts:
                        U, s, _ = np.linalg.svd(m)
                        if (s > 1e-9).sum() >= 2 and generic(Ue, U):
                            surv.append((a, b, eu, ed)); break
    nontriv = [s for s in surv if len(set(s[0])) > 1 or len(set(s[1])) > 1]
    return len(surv), len(nontriv)


def max_offblock(M, blocks):
    best = None
    lab = {}
    for bi, blk in enumerate(blocks):
        for i in blk:
            lab[i] = bi
    for perm in itertools.permutations(range(3)):
        v = max(M[i, j] for i in range(3) for j in range(3) if lab[i] != lab[perm[j]])
        best = v if best is None else min(best, v)
    return best


def step1():
    say('== STEP 1: concrete flavour Lagrangian')
    for nm, B in (('B3', B3), ('B5', B5)):
        for b in B:
            assert not sm_commuting(b)
    sp = {p: nullspace(np.array([[W16[NAME_IDX[x]][i] - W16[NAME_IDX[mem[0]]][i] for i in range(5)]
                                 for mem in MULT.values() for x in mem[1:]]), p).tolist() for p in (3, 5, 7, 101)}
    for p, v in sp.items():
        assert v == [[1, 1, 1, 0, 0], [0, 0, 0, 1, 1]]
    say('  none of the B3/B5 badge rows commutes with SU(3)c x SU(2)L; the commuting family Cartan charges are')
    say('  exactly span{S3=(1,1,1,0,0), S2=(0,0,0,1,1)} = span{3(B-L), -4 T3R} per family (checked mod 3,5,7,101)')
    # codes from the SM-compatible lattice
    res = {}
    for nm, B in (('(B-L)_g+T3R_g', BSM), ('(B-L)_g only', BSM[:1]), ('T3R_g only', BSM[1:]), ('Y_g only (tri-hypercharge)', ((2, 2, 2, -3, -3),))):
        for nf, ps in ((2, (3, 5, 7, 11)), (3, (3, 5, 7))):
            for p in ps:
                if nm.startswith('Y_g') and p in (3, 11):
                    continue
                c = coupled(B, nf, p)
                d = c.z_distance(3)
                xd = c.x_distance() if (c.r <= 9 and p <= 7) else None
                res[(nm, nf, p)] = (c.k, d, xd)
                say(f'  {nm:26s} fam={nf} p={p:2d}  [[{c.n},{c.k},{d}]]  X-distance {xd}')
    for p in (5, 7, 11):
        assert res[('(B-L)_g+T3R_g', 2, p)][:2] == (2, 3)
    for p in (5, 7):
        assert res[('(B-L)_g+T3R_g', 3, p)][:2] == (3, 3)
        assert res[('(B-L)_g+T3R_g', 2, p)][2] == 7 and res[('(B-L)_g+T3R_g', 3, p)][2] == 7
        assert res[('Y_g only (tri-hypercharge)', 2, p)][0] == 1
        assert res[('T3R_g only', 2, p)][1] == 2 and res[('(B-L)_g only', 2, p)][1] == 2
    assert res[('(B-L)_g+T3R_g', 2, 3)][1] == 2
    # Yukawa scan
    rng = np.random.default_rng(1)
    for G, p, two in ((2, 3, False), (2, 5, False), (2, 7, False), (3, 3, False), (3, 5, False), (3, 7, False),
                      (2, 5, True), (3, 3, True), (3, 5, True)):
        ns, nt = yukawa_scan(G, p, two, rng)
        say(f'  Yukawa scan G={G} p={p} {"2HDM" if two else "1HDM"}: realistic assignments {ns}, with family-differential charges {nt}')
        assert nt == 0 and ns == p * p
    # measured breaking sizes
    V = np.array([[0.97401, 0.22650, 0.00361], [0.22636, 0.97320, 0.04053], [0.00854, 0.03978, 0.999172]])
    s12, s23, s13 = np.sqrt([0.307, 0.561, 0.02195]); dl = np.radians(177)
    c12, c23, c13 = np.sqrt(1 - np.array([s12, s23, s13]) ** 2)
    U = np.abs(np.array([[c12 * c13, s12 * c13, s13 * np.exp(-1j * dl)],
                         [-s12 * c23 - c12 * s23 * s13 * np.exp(1j * dl), c12 * c23 - s12 * s23 * s13 * np.exp(1j * dl), s23 * c13],
                         [s12 * s23 - c12 * c23 * s13 * np.exp(1j * dl), -c12 * s23 - s12 * c23 * s13 * np.exp(1j * dl), c23 * c13]]))
    for nm, blocks in (('3 vs (1,2)', [[0, 1], [2]]), ('1 vs (2,3)', [[0], [1, 2]]), ('all three distinct', [[0], [1], [2]])):
        say(f'  doublet charge separating {nm:18s}: smallest achievable largest off-block mixing  CKM {max_offblock(V, blocks):.4f}   PMNS {max_offblock(U, blocks):.3f}')
    m = dict(u=1.23e-3, c=0.620, t=168.26, d=2.67e-3, s=53.16e-3, b=2.839, e=0.48307e-3, mu=0.101766, tau=1.72856)
    say(f'  right-handed (T3R_g) charges: forbidden light columns, size m2/m3 = up {m["c"]/m["t"]:.4f}, down {m["s"]/m["b"]:.4f}, lepton {m["mu"]/m["tau"]:.4f}')
    return res


# =====================================================================================  STEP 5
def step5():
    say('== STEP 5: distance 4 and entangling transversal gates')
    # (a) single family is rigid
    for p in (3, 5, 7):
        P = [(WM[a] * WM[b]) % p for a in range(5) for b in range(a, 5)]
        assert len(nullspace(np.array(P), p)) == 5
    say('  single family: the only rows compatible with the five SO(10) Cartan rows are the Cartan rows (p = 3, 5, 7)')
    # (b) integer (exactly anomaly-free) ternary antisymmetric rows (y,-y)
    supp = [[j for j in range(16) if m >> j & 1] for m in range(1, 1 << 16)]
    supp = [s for s in supp if not WM[:, s].sum(1).any()]
    c1, c0 = [], []
    for S in supp:
        if 0 in S:
            rest = [j for j in S if j]
            for sg in itertools.product((1, -1), repeat=len(rest)):
                y = np.zeros(16, np.int64); y[0] = 1; y[rest] = sg; c1.append(y)
        else:
            for sg in itertools.product((1, -1), repeat=len(S) - 1):
                y = np.zeros(16, np.int64); y[S[0]] = 1; y[S[1:]] = sg; c0.append(y)
    c1 = np.array(c1); Y2 = np.vstack([c1, np.array(c0)])
    perms = [np.array([IDX[tuple(q[s[i]] for i in range(5))] for q in W16]) for s in itertools.permutations(range(5))]
    pw = 3 ** np.arange(16)
    K = np.full(len(c1), np.iinfo(np.int64).max)
    for P in perms:
        Yp = np.empty_like(c1); Yp[:, P] = c1; K = np.minimum(K, (Yp + 1) @ pw)
    reps = c1[np.unique(K, return_index=True)[1]]
    npairs = 0; passing = 0
    for y1 in reps:
        cand = Y2[~((WM * y1) @ Y2.T).any(0)]; npairs += len(cand)
        for y2 in cand:
            t = y2[0]
            if not any((y1[q], y2[q]) in ((1, t), (-1, -t)) for q in range(1, 16)):
                passing += 1
    say(f'  integer ternary isotropic rows: {len(supp)} zero-sum supports, {len(c1) + len(c0)} vectors, {npairs} isotropic pairs; pairs protecting all twins: {passing}')
    assert passing == 0
    Bl = [np.array([dot(b, q) for q in W16]) for b in B5]
    okB = ~((np.vstack([WM * l for l in Bl]) @ Y2.T).any(0))
    nonlin_int = [y for y in Y2[okB] if rank(np.vstack([WM, y]), 101) == 6]
    say(f'  integer ternary rows compatible with B5: {int(okB.sum())}, of which nonlinear: {len(nonlin_int)}')
    assert len(nonlin_int) == 0
    # (c) mod-p nonlinear rows extending B5
    summary = {}
    best = None
    for p in ((5, 7, 11) if QUICK else (5, 7, 11, 13)):
        L = [l % p for l in Bl]
        V = nullspace(np.vstack([(WM * l) % p for l in L]), p)
        Cs = rref(WM % p, p)[0]; piv = [int(np.nonzero(r)[0][0]) for r in Cs]
        iso = []
        for st in range(0, p ** len(V), 400000):
            Ys = (vectors_range(len(V), p, st, st + 400000) @ V) % p
            iso.append(Ys[~((WM @ ((Ys * Ys) % p).T) % p).any(0)])
        iso = np.vstack(iso); Z = iso.copy()
        for r, pc in zip(Cs, piv):
            Z = (Z - np.outer(Z[:, pc], r)) % p
        nl = iso[Z.any(1)]
        seen = {}
        for y in nl:
            key = None
            for lam in range(1, p):
                z = (lam * y) % p
                for P in perms:
                    zp = np.empty_like(z); zp[P] = z; k = tuple(zp.tolist())
                    key = k if key is None or k < key else key
            seen.setdefault(key, y)
        counts = {}
        for y in seen.values():
            yl = [int(v) for v in y]
            c = fam_code(2, p, lambda q, g, yl=yl: [(1 - 2 * g) * dot(b, q) for b in B5] + [(1 - 2 * g) * yl[IDX[q]]])
            if c.k != 2:
                continue
            d = c.z_distance(3); n3 = len(c.z_logical_supports(3)) if d == 3 else 0
            counts[(d, n3)] = counts.get((d, n3), 0) + 1
            if p == 7 and best is None:
                best = (yl, c)
        if p == 7:   # no compatible second nonlinear row (mutually isotropic, independent of Cartan + first)
            second = 0
            for y1 in seen.values():
                ok = ~((WM @ ((y1 * nl) % p).T) % p).any(0)
                second += sum(1 for y2 in nl[ok] if rank(np.vstack([Cs, y1, y2]), p) == 7)
            assert second == 0
            say(f'  p= 7: compatible second nonlinear rows: {second}')
        summary[p] = (len(nl), len(seen), counts)
        say(f'  p={p:2d}: nonlinear isotropic rows extending B5: {len(nl)} ({len(seen)} orbits); (dZ, #weight-3 supports) -> orbits {counts}')
    assert summary[5][0] == 0 and summary[11][0] == 0
    assert summary[7][:2] == (588, 17) and summary[7][2] == {(3, 8): 17}
    if not QUICK:
        assert summary[13][:2] == (4056, 57) and summary[13][2] == {(3, 8): 57}
    yl, c = best
    G = c.Gint
    assert all(int((G[i] * G[j] * G[k]).sum()) % 7 == 0 for i in range(len(G)) for j in range(i, len(G)) for k in range(j, len(G)))
    assert not all(int((G[i] * G[j] * G[k]).sum()) == 0 for i in range(len(G)) for j in range(i, len(G)) for k in range(j, len(G)))
    D = class_weight_distributions(c)
    anyl = [sum(D[cc][w] for cc in D if any(cc)) for w in range(31)]
    Dref = class_weight_distributions(coupled(B5, 2, 7))
    anyr = [sum(Dref[cc][w] for cc in Dref if any(cc)) for w in range(31)]
    s_new, s_ref = to_eps(anyl, 7, 30, 4), to_eps(anyr, 7, 30, 4)
    xd = c.x_distance()
    say(f'  example p=7 row y = {yl}: [[30,2,3]]_7, X-distance {xd}, triorthogonal mod 7 but not over Z;')
    say(f'      P(any logical) = {s_new[3]} e^3 + ...   vs coupled B5 at p=7: {s_ref[3]} e^3   (ratio {s_ref[3] / s_new[3]})')
    assert s_new[3] == Fraction(2, 9) and s_ref[3] == Fraction(25, 18) and xd == 9
    # (d) transversal diagonal gates
    def nl7(g): return lambda q: list(q) + [(1 - 2 * g) * dot(b, q) for b in B5] + [(1 - 2 * g) * yl[IDX[q]]]
    for nm, cc in (('direct 2fam p5', direct_sum(2, 5)), ('coupled B5 2fam p5', coupled(B5, 2, 5)), ('coupled B5 2fam p7', coupled(B5, 2, 7)),
                   ('coupled B5 3fam p5', coupled(B5, 3, 5)), ('SM (B-L,T3R) 2fam p5', coupled(BSM, 2, 5)),
                   ('SM (B-L,T3R) 3fam p5', coupled(BSM, 3, 5)), ('nonlinear-row 2fam p7', build([nl7(0), nl7(1)], 7)),
                   ('family-number 2fam p5', build([lambda q: list(q) + [1], lambda q: list(q) + [-1]], 5))):
        g3, g2 = gate_space(cc, 3), gate_space(cc, 2)
        say(f'  {nm:24s} cubic gates: {g3["dim_valid"]} valid power vectors, logical rank {g3["rank_logical"]}, entangling rank {g3["rank_mixed"]};'
            f' quadratic: {g2["dim_valid"]}, {g2["rank_logical"]}, {g2["rank_mixed"]}')
        assert g3['rank_mixed'] == 0 and g2['rank_mixed'] == 0
    return summary


# =====================================================================================  STEP 4
def step4():
    say('== STEP 4: operational cost of measuring the charge constraints')
    def xgens(c):
        p = c.p; S = (all_vectors(len(c.K), p) @ c.SX) % p
        w = (S != 0).sum(1); basis = []
        for i in np.argsort(w):
            if w[i] and rank(np.array(basis + [S[i]]), p) > len(basis):
                basis.append(S[i])
            if len(basis) == len(c.K):
                break
        return sorted(int((b != 0).sum()) for b in basis)

    def zgens(c, wmax=4):
        p = c.p; target = c.n - c.r; basis = []; ws = []
        for w in range(2, wmax + 1):
            for J in itertools.combinations(range(c.n), w):
                for v in nullspace(c.GP[:, J], p):
                    if (v != 0).sum() == w:
                        e = np.zeros(c.n, np.int64); e[list(J)] = v
                        if rank(np.array(basis + [e]), p) > len(basis):
                            basis.append(e); ws.append(w)
                if len(basis) == target:
                    return ws
        return ws
    out = {}
    for nm, c in (('single p5', single(5)), ('direct 2fam p5', direct_sum(2, 5)), ('coupled B5 2fam p5', coupled(B5, 2, 5)),
                  ('SM (B-L,T3R) 2fam p5', coupled(BSM, 2, 5))):
        xs, zs = xgens(c), zgens(c)
        assert len(zs) == c.n - c.r and set(zs) == {4}
        out[nm] = (xs, len(zs), sum(xs) + sum(zs))
        say(f'  {nm:22s} X-check weights {xs}; {len(zs)} Z-checks of weight 4; total check weight {sum(xs) + sum(zs)}')
    assert out['coupled B5 2fam p5'][0] == [15, 15, 16, 16, 16] and out['direct 2fam p5'][0] == [8] * 8
    return out


if __name__ == '__main__':
    step2(); step3(); step1(); step5(); step4()
    print('ALL OPEN-STEP CHECKS PASSED', f'({time.time() - T0:.0f}s)')
