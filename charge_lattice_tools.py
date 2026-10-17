"""Exact mod-p toolkit for charge-lattice CSS codes built from Standard Model generations.

Used by open_steps.py.  Requires numpy; all arithmetic is exact integer / Fraction arithmetic.

Contents
  * mod-p linear algebra (rref, rank, nullspace)
  * ChargeCode: CSS code from an integer charge matrix punctured on the right-handed-neutrino columns
  * builders: single generation, coupled (smart-badge) lattices, direct sums, arbitrary family maps
  * exact logical-channel enumerators under i.i.d. Z noise (Fourier / MacWilliams), direct-sum convolution
  * stabilizer-group decomposability (matroid components)
  * transversal diagonal-gate classification (cubic / quadratic per-qudit phase gates)
"""
from __future__ import annotations
import numpy as np
from fractions import Fraction
from math import comb
import itertools
from itertools import combinations, product

Q0 = (1, 1, 1, 1, 1)


def weights():
    """16 SO(10) spinor weights, q0 first, then SU(5) 10, then 5bar (same order as matter_qec_study.py)."""
    ten = []
    for i, j in combinations(range(5), 2):
        q = [1] * 5; q[i] = q[j] = -1; ten.append(tuple(q))
    fb = []
    for i in range(5):
        q = [-1] * 5; q[i] = 1; fb.append(tuple(q))
    return [Q0] + ten + fb


W16 = weights()
NAMES1 = ['nuc', 'uc3', 'uc2', 'd1', 'u1', 'uc1', 'd2', 'u2', 'd3', 'u3', 'ec', 'dc1', 'dc2', 'dc3', 'nu', 'e']


# ---------------------------------------------------------------- linear algebra mod p
def rref(M, p):
    A = np.array(M, dtype=np.int64) % p
    if A.ndim == 1:
        A = A[None, :]
    rows, cols = A.shape
    piv = []
    r = 0
    for c in range(cols):
        if r == rows:
            break
        nz = np.nonzero(A[r:, c])[0]
        if len(nz) == 0:
            continue
        i = r + nz[0]
        A[[r, i]] = A[[i, r]]
        A[r] = (A[r] * pow(int(A[r, c]), -1, p)) % p
        f = A[:, c].copy(); f[r] = 0
        A = (A - np.outer(f, A[r])) % p
        piv.append(c); r += 1
    return A[:r], piv


def rank(M, p):
    M = np.asarray(M)
    if M.size == 0:
        return 0
    return len(rref(M, p)[1])


def nullspace(M, p):
    """Basis (rows) of {x : M x = 0 mod p}."""
    M = np.asarray(M, dtype=np.int64) % p
    n = M.shape[1]
    R, piv = rref(M, p) if M.shape[0] else (np.zeros((0, n), np.int64), [])
    free = [c for c in range(n) if c not in piv]
    basis = []
    for f in free:
        x = np.zeros(n, np.int64); x[f] = 1
        for i, c in enumerate(piv):
            x[c] = (-R[i, f]) % p
        basis.append(x)
    return np.array(basis, dtype=np.int64).reshape(len(basis), n)


def all_vectors(dim, p, chunk=None):
    """All vectors of F_p^dim as an int64 array (dim small)."""
    if dim == 0:
        return np.zeros((1, 0), np.int64)
    idx = np.arange(p ** dim, dtype=np.int64)
    out = np.empty((p ** dim, dim), np.int64)
    for i in range(dim):
        out[:, dim - 1 - i] = (idx // p ** i) % p
    return out


def vectors_range(dim, p, start, stop):
    """Vectors number start..stop-1 of F_p^dim (same ordering as all_vectors), without materialising all."""
    idx = np.arange(start, min(stop, p ** dim), dtype=np.int64)
    out = np.empty((len(idx), dim), np.int64)
    for i in range(dim):
        out[:, dim - 1 - i] = (idx // p ** i) % p
    return out


# ---------------------------------------------------------------- codes
class ChargeCode:
    """CSS code from an integer charge matrix G (rows = charges, columns = species),
    punctured on the columns in `punct` (the logical coordinates)."""

    def __init__(self, G, p, punct, labels=None, families=None, name=''):
        self.Gint = np.array(G, dtype=np.int64)
        self.p = p
        self.name = name
        self.punct = list(punct)
        N = self.Gint.shape[1]
        self.phys = [j for j in range(N) if j not in self.punct]
        G = self.Gint % p
        # drop dependent rows mod p
        R, piv = rref(G.T, p)  # column-space pivots of G^T = independent rows of G
        self.G = G[piv]
        self.r = len(piv)
        self.GR = self.G[:, self.punct]
        self.GP = self.G[:, self.phys]
        self.n = len(self.phys)
        self.k = rank(self.GR.T, p) if self.punct else 0
        assert rank(self.GP, p) == self.r, 'rank loss after puncture'
        self.labels = labels
        self.families = families  # family index of each physical qudit
        # X-stabilizer coefficient space K = {a : a G_R = 0}
        self.K = nullspace(self.GR.T, p) if self.punct else np.eye(self.r, dtype=np.int64)
        self.SX = (self.K @ self.GP) % p  # generators of X stabilizers (rows)

    def params(self):
        return dict(n=self.n, k=self.k, charge_rank=self.r, X_stab=self.r - self.k, Z_stab=self.n - self.r)

    # Z-type error e (length n) : syndrome-free iff G_P e in span(G_R)
    def classify(self, e):
        s = (self.GP @ np.asarray(e)) % self.p
        if not s.any():
            return 'stabilizer'
        if rank(np.column_stack([self.GR, s]), self.p) == rank(self.GR, self.p):
            return 'logical'
        return 'detected'

    def logical_class(self, e):
        """Return c with G_P e = G_R c (or None if detected)."""
        s = (self.GP @ np.asarray(e)) % self.p
        A = np.column_stack([self.GR, s])
        R, piv = rref(A, self.p)
        if self.k in piv:
            return None
        c = np.zeros(self.k, np.int64)
        for i, col in enumerate(piv):
            c[col] = R[i, -1]
        return tuple(int(x) for x in c)

    def z_logical_supports(self, w):
        p = self.p
        rR = rank(self.GR, p)
        out = []
        for J in combinations(range(self.n), w):
            cs = self.GP[:, J]
            if rank(np.column_stack([self.GR, cs]), p) < rR + rank(cs, p):
                out.append(J)
        return out

    def z_distance(self, wmax=4):
        for w in range(1, wmax + 1):
            if self.z_logical_supports(w):
                return w
        return None

    def x_distance(self):
        """Exact: min weight of a G_P over a with a G_R != 0 (enumerates F_p^r)."""
        p = self.p
        best = None
        for st in range(0, p ** self.r, 200000):
            a = vectors_range(self.r, p, st, st + 200000)
            b = (a @ self.GR) % p
            w = ((a @ self.GP) % p != 0).sum(1)
            mask = b.any(1)
            if mask.any():
                m = int(w[mask].min())
                best = m if best is None else min(best, m)
        return best


def build(maps, p, name='', with_labels=True):
    """maps: list over families of functions q -> integer column (list).
    Columns ordered family-major, 16 per family, puncture q0 of each family."""
    cols, labels, fams = [], [], []
    for g, f in enumerate(maps):
        for i, q in enumerate(W16):
            cols.append(list(f(q)))
            labels.append(f'{NAMES1[i]}[{g + 1}]')
            fams.append(g)
    G = np.array(cols, dtype=np.int64).T
    punct = [16 * g for g in range(len(maps))]
    phys_f = [fams[j] for j in range(len(fams)) if j not in punct]
    phys_l = [labels[j] for j in range(len(labels)) if j not in punct]
    return ChargeCode(G, p, punct, labels=phys_l, families=phys_f, name=name)


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


B3 = ((1, 0, 0, 0, 0), (0, 1, 0, 0, 1), (0, 0, 1, 1, 0))
B5 = ((1, 0, 0, 0, 0), (0, 1, 1, 1, 1))
FAM3 = ((1, 0), (-1, 0), (0, 1))


def coupled(B, nfam, p, labels_vec=None, name=''):
    """Common SO(10) rows + differential rows f_g[axis] * (b.q). For nfam=2 uses s=+1,-1."""
    if labels_vec is None:
        labels_vec = ((1,), (-1,)) if nfam == 2 else FAM3
    naxes = len(labels_vec[0])

    def mk(g):
        f = labels_vec[g]
        return lambda q: list(q) + [f[ax] * dot(b, q) for ax in range(naxes) for b in B]
    return build([mk(g) for g in range(nfam)], p, name=name)


def direct_sum(nfam, p, name=''):
    def mk(g):
        return lambda q: [list(q)[i] if gg == g else 0 for gg in range(nfam) for i in range(5)]
    return build([mk(g) for g in range(nfam)], p, name=name)


def single(p, name=''):
    return build([lambda q: list(q)], p, name=name)


# ---------------------------------------------------------------- exact logical-channel enumerators
def histogram(code):
    p, r, k, n = code.p, code.r, code.k, code.n
    H = np.zeros((n + 1, p ** k), dtype=np.int64)
    pw = p ** np.arange(k, dtype=np.int64)
    for st in range(0, p ** r, 250000):
        a = vectors_range(r, p, st, st + 250000)
        w = ((a @ code.GP) % p != 0).sum(1)
        bi = ((a @ code.GR) % p) @ pw
        np.add.at(H, (w, bi), 1)
    return H


def class_weight_distributions(code, H=None):
    """Return {c: [A_0..A_n]} exact counts of Z-errors of each weight in logical class c
    (c=0 includes the trivial class = Z stabilizers)."""
    p, r, k, n = code.p, code.r, code.k, code.n
    if H is None:
        H = histogram(code)
    B = all_vectors(k, p)  # b vectors indexed like bi
    # polynomial in x for each w(a):  (1+(p-1)x)^{n-w} (1-x)^w   (integer coefficients)
    def poly(w):
        c = [0] * (n + 1)
        for i in range(n - w + 1):
            for j in range(w + 1):
                c[i + j] += comb(n - w, i) * (p - 1) ** i * comb(w, j) * (-1) ** j
        return c
    polys = [poly(w) for w in range(n + 1)]
    # M_b(x) = sum_w H[w,b] poly_w(x)
    M = {}
    for bi in range(p ** k):
        acc = [0] * (n + 1)
        for w in range(n + 1):
            h = int(H[w, bi])
            if h:
                pw_ = polys[w]
                for i in range(n + 1):
                    acc[i] += h * pw_[i]
        M[bi] = acc
    out = {}
    for c in all_vectors(k, p):
        c = tuple(int(x) for x in c)
        # sum_b w^{-b.c} M_b ; group by t = b.c
        At = {t: [0] * (n + 1) for t in range(p)}
        for bi, b in enumerate(B):
            t = int(np.dot(b, c)) % p
            for i in range(n + 1):
                At[t][i] += M[bi][i]
        if any(c):
            for t in range(2, p):
                assert At[t] == At[1]
            tot = [At[0][i] - At[1][i] for i in range(n + 1)]
        else:
            tot = [sum(At[t][i] for t in range(p)) for i in range(n + 1)]
        assert all(x % p ** r == 0 for x in tot)
        out[c] = [x // p ** r for x in tot]
    return out


def to_eps(dist, p, n, order):
    """Probability series sum_w A_w (eps/(p-1))^w (1-eps)^(n-w) up to eps^order (exact Fractions)."""
    s = [Fraction(0)] * (order + 1)
    for w, A in enumerate(dist):
        if not A or w > order:
            continue
        base = Fraction(A, (p - 1) ** w)
        for m in range(0, order - w + 1):
            s[w + m] += base * comb(n - w, m) * (-1) ** m
    return s


def exact_value(dist, p, n, eps):
    eps = Fraction(eps)
    return sum(Fraction(A) * (eps / (p - 1)) ** w * (1 - eps) ** (n - w) for w, A in enumerate(dist) if A)


def brute_low_weight(code, wmax):
    """Independent check: enumerate all Z errors of weight <= wmax; return {(w,c): count} and
    family-split counts for logical errors."""
    p, n = code.p, code.n
    K = code.K % p
    # left inverse of G_R on its column span
    counts = {}
    split = {}
    for w in range(1, wmax + 1):
        vals = np.array(list(product(range(1, p), repeat=w)), dtype=np.int64)  # V x w
        for J in combinations(range(n), w):
            s = (code.GP[:, J] @ vals.T) % p  # r x V
            undet = ~((K @ s) % p).any(0) if len(K) else np.ones(s.shape[1], bool)
            if not undet.any():
                continue
            for v in np.nonzero(undet)[0]:
                e = np.zeros(n, np.int64); e[list(J)] = vals[v]
                c = code.logical_class(e)
                key = (w, c)
                counts[key] = counts.get(key, 0) + 1
                if any(c) and code.families is not None:
                    fs = tuple(sorted([sum(code.families[j] == g for j in J) for g in set(code.families)], reverse=True))
                    split[(w, c, fs)] = split.get((w, c, fs), 0) + 1
    return counts, split


def direct_sum_distributions(D1, n1, copies):
    """Class weight distributions of a direct sum of identical codes (exact convolution)."""
    cur = {(): [1]}
    for _ in range(copies):
        nxt = {}
        for a, da in cur.items():
            for b, db in D1.items():
                conv = [0] * (len(da) + len(db) - 1)
                for i, x in enumerate(da):
                    if x:
                        for j, y in enumerate(db):
                            if y:
                                conv[i + j] += x * y
                nxt[a + b] = conv
        cur = nxt
    return cur


# ---------------------------------------------------------------- decomposability
class UF:
    def __init__(s,n): s.p=list(range(n))
    def f(s,x):
        while s.p[x]!=x: s.p[x]=s.p[s.p[x]]; x=s.p[x]
        return x
    def u(s,a,b): s.p[s.f(a)]=s.f(b)

def matroid_components(M,p,uf=None):
    """Connected components of the column matroid of M (rows span the code)."""
    M=np.asarray(M)%p; n=M.shape[1]
    uf=uf or UF(n)
    if M.shape[0]==0: return uf
    R,piv=rref(M,p)
    for j in range(n):
        if j in piv: continue
        for i,c in enumerate(piv):
            if R[i,j]%p: uf.u(j,c)
    return uf

def split_components(code):
    """Finest bipartition-compatible decomposition of S = S_X x S_Z.
    S_Z = ker G_P has the same matroid components as its dual, rowspace(G_P)."""
    uf=UF(code.n)
    matroid_components(code.SX,code.p,uf)
    matroid_components(code.GP,code.p,uf)
    comps={}
    for j in range(code.n): comps.setdefault(uf.f(j),[]).append(j)
    return sorted(comps.values(),key=len,reverse=True)


# ---------------------------------------------------------------- transversal diagonal gates
def gate_space(code, degree):
    p=code.p; H=code.GP; r=code.r; K=code.K % p
    eqs=[]
    if degree==3:
        for i in range(r):
            for l in range(i,r):
                for s in K:
                    eqs.append((H[i]*H[l]*((s@H)%p))%p)
    else:
        for i in range(r):
            for s in K:
                eqs.append((H[i]*((s@H)%p))%p)
    Wsp=nullspace(np.array(eqs),p) if eqs else np.eye(code.n,dtype=np.int64)
    # preimages a^(i) with a G_R = e_i
    A=[]
    for i in range(code.k):
        # solve a G_R = e_i  ->  G_R^T a^T = e_i
        M=np.column_stack([code.GR.T%p, np.eye(code.k,dtype=np.int64)[i]])
        R,piv=rref(M,p); a=np.zeros(r,np.int64)
        for row,c in zip(R,piv): a[c]=row[-1]
        assert ((a@code.GR)%p==np.eye(code.k,dtype=np.int64)[i]).all()
        A.append(a)
    X=[(a@H)%p for a in A]  # values x_j on logical basis preimages
    mixed=[]; forms=[]
    idxs=list(itertools.combinations_with_replacement(range(code.k),degree))
    for w in Wsp:
        coef={}
        for ix in idxs:
            prod=np.ones(code.n,np.int64)
            for i in ix: prod=(prod*X[i])%p
            coef[ix]=int((w*prod).sum()%p)
        forms.append(coef)
    # rank of the map w -> descended form, and of its mixed part
    F=np.array([[f[ix] for ix in idxs] for f in forms]) if forms else np.zeros((0,len(idxs)),np.int64)
    mixed_cols=[j for j,ix in enumerate(idxs) if len(set(ix))>1]
    rk=rank(F,p) if len(F) else 0
    rk_mixed=rank(F[:,mixed_cols],p) if len(F) and mixed_cols else 0
    return dict(dim_valid=len(Wsp), rank_logical=rk, rank_mixed=rk_mixed, monomials=len(idxs))

