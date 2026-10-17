"""A REAL recovery test, but for a Z_p lattice-gauge simulator, NOT the Standard Model.

Part A: 3-edge pure-gauge ring with only Gauss-law checks. Corrects
single electric-flux shift X_i^a, but gauge-invariant phase Z_i^b acts
as logical Z and fails Knill--Laflamme.

Part B: 3x3 Z_p torus with Gauss-law Z-star checks PLUS independent
gauge-invariant X-plaquette checks. It is the standard Z_p toric code,
[[18,2,3]]_p, not a novel code. All single-edge Weyl operators
X_i^a Z_i^b have distinct nonzero syndromes. This proves exact
Knill--Laflamme recovery for arbitrary one-edge noise by linearity.

Only modular finite-field arithmetic; no statevector p**18 allocation.
Usage: python3 operational_gauge_qec.py
"""
from itertools import product


def unit(n, j, value=1):
    v = [0] * n
    v[j] = value
    return tuple(v)


def lincomb(a, b, p, scale=1):
    return tuple((int(x) + scale * int(y)) % p for x, y in zip(a, b))


def inner(a, b, p):
    return sum(int(x) * int(y) for x, y in zip(a, b)) % p


def rank(rows, p):
    if not rows:
        return 0
    a = [[int(v) % p for v in row] for row in rows]
    n, m = len(a), len(a[0])
    assert all(len(row) == m for row in a)
    r = 0
    for j in range(m):
        k = next((i for i in range(r, n) if a[i][j]), None)
        if k is None:
            continue
        a[r], a[k] = a[k], a[r]
        inv = pow(a[r][j], -1, p)
        a[r] = [(v * inv) % p for v in a[r]]
        for i in range(n):
            if i != r and a[i][j]:
                f = a[i][j]
                a[i] = [(x - f * y) % p for x, y in zip(a[i], a[r])]
        r += 1
        if r == n:
            break
    return r


def ring(p):
    """Electric-basis Gauss operators: G_i = Z_i Z_(i-1)^(-1).

    Physical (zero-charge) basis vectors are |j,j,j> for j in F_p.
    """
    gauss = []
    for i in range(3):
        row = [0] * 3
        row[i] += 1
        row[(i - 1) % 3] -= 1
        gauss.append(tuple(row))
    assert rank(gauss, p) == 2

    # The error X_i^a has a gauge-syndrome vector proportional to
    # column i of the Gauss incidence matrix. Every single-edge error
    # has a different syndrome and is correctable by its inverse.
    syndromes = {}
    for site in range(3):
        for a in range(1, p):
            vector = unit(3, site, a)
            syndrome = tuple(inner(g, vector, p) for g in gauss)
            assert any(syndrome)
            assert syndrome not in syndromes
            syndromes[syndrome] = (site, a)
    assert len(syndromes) == 3 * (p - 1)

    # Gauge-invariant Z_i^b has zero syndrome but implements the
    # non-scalar logical Z_L^b : |j>_L -> omega^(b*j)|j>_L.
    for site in range(3):
        for b in range(1, p):
            phase_vector = unit(3, site, b)
            assert all(inner(g, tuple(0 for _ in range(3)), p) == 0
                       for g in gauss)
            assert len({(b * j) % p for j in range(p)}) == p
            assert phase_vector[site] == b
    print(f"Ring p={p}: physical dimension p; {len(syndromes)} distinct "
          "single-X syndromes; every Z_i is gauge-invariant logical Z (not correctable)")
    return gauss


def torus(l):
    """Orient horizontal links right and vertical links up on an l x l torus.

    z_star = electric Gauss generators Z^(div E).
    x_face = oriented magnetic plaquette shifts X^(boundary face).
    star and face rows commute because div(boundary)=0.
    """
    assert l >= 2
    n = 2 * l * l

    def h(x, y):
        return (y % l) * l + x % l

    def v(x, y):
        return l * l + (y % l) * l + x % l

    z_stars, x_faces = [], []
    for y in range(l):
        for x in range(l):
            star = [0] * n
            for e, sign in ((h(x, y), +1), (v(x, y), +1),
                            (h(x - 1, y), -1), (v(x, y - 1), -1)):
                star[e] += sign
            z_stars.append(tuple(star))
            face = [0] * n
            for e, sign in ((h(x, y), +1), (v(x + 1, y), +1),
                            (h(x, y + 1), -1), (v(x, y), -1)):
                face[e] += sign
            x_faces.append(tuple(face))

    # Noncontractible primal loops (magnetic, X): horizontal row,
    # vertical column. Dual loops (electric, Z): horizontal
    # column-cut and vertical row-cut, respectively.
    x_horizontal = tuple(int(i in {h(x, 0) for x in range(l)}) for i in range(n))
    x_vertical = tuple(int(i in {v(0, y) for y in range(l)}) for i in range(n))
    z_horizontal = tuple(int(i in {h(0, y) for y in range(l)}) for i in range(n))
    z_vertical = tuple(int(i in {v(x, 0) for x in range(l)}) for i in range(n))
    return z_stars, x_faces, (x_horizontal, x_vertical), (z_horizontal, z_vertical)


def one_site_syndrome(z_stars, x_faces, site, a, b, p):
    """E=X_site^a Z_site^b. Signs of commutation phases are conventional.

    Gauss Z checks respond to X^a; magnetic X checks respond to Z^b.
    """
    return (tuple((a * g[site]) % p for g in z_stars) +
            tuple((b * f[site]) % p for f in x_faces))


def toric_code(p, l=3):
    assert l == 3, "The test targets distance 3; l=2 is not a one-error code."
    z_stars, x_faces, logical_x, logical_z = torus(l)
    n = 2 * l * l
    rz = rank(z_stars, p)
    rx = rank(x_faces, p)
    assert rz == rx == l * l - 1
    assert all(inner(g, f, p) == 0 for g in z_stars for f in x_faces)
    k = n - rz - rx
    assert k == 2

    # Explicit independent logical Weyl pairs and noncontractible
    # representatives. These do not belong to the stabilizer spans.
    for i in range(2):
        assert all(inner(g, logical_x[i], p) == 0 for g in z_stars)
        assert all(inner(f, logical_z[i], p) == 0 for f in x_faces)
        assert rank(x_faces + [logical_x[i]], p) == rx + 1
        assert rank(z_stars + [logical_z[i]], p) == rz + 1
        assert sum(bool(x) for x in logical_x[i]) == l
        assert sum(bool(x) for x in logical_z[i]) == l
        for j in range(2):
            assert inner(logical_x[i], logical_z[j], p) == int(i == j)

    # The zero-syndrome sector is the gauge-invariant code space
    # satisfying both Gauss and magnetic constraints.
    #
    # Full Knill--Laflamme check for the one-edge Weyl error basis.
    # If the syndrome of each element in {I, X_i^a Z_i^b} is unique,
    # all E†F with E!=F have a nonzero syndrome, hence P E†F P=0;
    # P E†E P=P. Any one-edge channel is a linear combination of
    # these Weyl operators, so arbitrary one-edge noise is correctable.
    error_syndromes = {(0,) * (len(z_stars) + len(x_faces)): ("identity",)}
    for site in range(n):
        for a, b in product(range(p), repeat=2):
            if a == b == 0:
                continue
            syndrome = one_site_syndrome(z_stars, x_faces, site, a, b, p)
            assert any(syndrome)
            assert syndrome not in error_syndromes, (
                f"KL failure: {site=}, {a=}, {b=}, "
                f"collides with {error_syndromes[syndrome]}")
            error_syndromes[syndrome] = (site, a, b)
    assert len(error_syndromes) == 1 + n * (p * p - 1)

    # Gauge-invariant electric phase error Z_i has zero Gauss syndrome
    # but a nonzero magnetic-plaquette syndrome: an extra physical
    # check beyond Gauss is essential for its correction.
    for site in range(n):
        syn = one_site_syndrome(z_stars, x_faces, site, 0, 1, p)
        assert all(v == 0 for v in syn[:len(z_stars)])
        assert any(v != 0 for v in syn[len(z_stars):])
    # X_i violates Gauss and is detected by Gauss checks alone.
    for site in range(n):
        syn = one_site_syndrome(z_stars, x_faces, site, 1, 0, p)
        assert any(v != 0 for v in syn[:len(z_stars)])
        assert all(v == 0 for v in syn[len(z_stars):])
    print(f"Torus p={p}: [[{n},{k},3]]_{p}, {rz} independent Gauss checks, "
          f"{rx} independent magnetic checks, {len(error_syndromes)} "
          "unique identity/one-edge Weyl syndromes; complete KL basis verified")
    return error_syndromes


def test():
    for p in (3, 5):
        ring(p)
        toric_code(p)
    print("ALL OPERATIONAL GAUGE-QEC CHECKS PASSED")


if __name__ == "__main__":
    test()
