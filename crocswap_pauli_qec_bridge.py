"""Extract genuine six-qubit stabilizer codes from CrocSwap's Pauli spread.

Source mathematical fixture (Apache-2.0, CrocSwap/integer-mult-bounds):
commit 20ecf2891d80b5acf18d980b2164fb44ea7faa90
research/atlas-frontier-applications/applications/adapters/pauli_recovery.py
research/atlas-frontier-applications/applications/adapters/
six-qubit-recovery/FULL-MEASUREMENT.json

The source recovers a missing *maximal commuting Pauli measurement class*.
It does NOT claim that this is quantum error correction of an unknown qubit.
Here we independently rebuild its six-qubit symplectic spread and construct
rank-five stabilizer subspaces inside its 65 rank-six measurement classes.

For EVERY hyperplane (65*(2**6-1)=4095), determine the stabilizer code's
distance by exact binary symplectic algebra. A proper rank-five subspace
encodes one logical qubit; a full rank-six class encodes ZERO qubits.

Also verify 19x19 Knill-Laflamme single-qubit Pauli error-pair conditions,
degenerate syndrome recovery, and a pair of anticommuting weight-three
logical Paulis for one fully specified [[6,1,3]] witness.

No numerical tolerances, external dependencies or physical Standard Model
claim. One source-linked code found does not imply 17 inequivalent codes
or physical gauge-compatible measurements/recovery.

Run: python crocswap_pauli_qec_bridge.py
"""
from collections import Counter, defaultdict

N = 6
MASK = (1 << N) - 1
FIELD_MODULUS = 67  # x^6+x+1, exactly as the source fixture
SOURCE_GAUGE_TRANSVECTIONS = (
    3746, 866, 813, 1804, 1971, 43, 3359, 3374, 166,
    123, 2499, 796, 2152, 3391, 1731, 3636, 1621,
)
SOURCE_MISSING_CLASS = (2076, 1069, 564, 314, 133, 83)
SOURCE_FIRST_CLASS = (2125, 1115, 584, 283, 148, 57)
SOURCE_LAST_CLASS = (2102, 1087, 540, 316, 177, 82)
CHOSEN_CLASS = 0
CHOSEN_HYPERPLANE = 51
CHOSEN_CHECKS = (3094, 584, 283, 2265, 2164)


def symplectic(a, b):
    """One iff their Hermitian Pauli representatives anticommute."""
    return (((a & MASK) & (b >> N)).bit_count()
            + ((b & MASK) & (a >> N)).bit_count()) % 2


def weight(p):
    return ((p & MASK) | (p >> N)).bit_count()


def pauli_label(p):
    return "".join(
        "Y" if (p >> i) & 1 and (p >> (N+i)) & 1
        else "X" if (p >> i) & 1
        else "Z" if (p >> (N+i)) & 1
        else "I" for i in range(N)
    )


def gf_multiply(a, b):
    """Polynomial basis for F_(2^6), matching CrocSwap's fixture."""
    result = 0
    while b:
        if b & 1:
            result ^= a
        b >>= 1
        a <<= 1
        if a >> N:
            a ^= FIELD_MODULUS
    return result


def gf_trace(a):
    total = 0
    for _ in range(N):
        total ^= a
        a = gf_multiply(a, a)
    assert total in (0, 1)
    return total


def reduced_binary_basis(vectors):
    """Independent ascending-source-generator -> descending pivot RREF."""
    pivot = {}
    for v in vectors:
        x = v
        for k in sorted(pivot, reverse=True):
            if (x >> k) & 1:
                x ^= pivot[k]
        if x:
            k = x.bit_length()-1
            for j in list(pivot):
                if (pivot[j] >> k) & 1:
                    pivot[j] ^= x
            pivot[k] = x
    return tuple(pivot[k] for k in sorted(pivot, reverse=True))


def span(g):
    out = {0}
    for v in g:
        out |= {u ^ v for u in out}
    return out


def reconstruct_source_spread():
    """Rebuild exact 65-class fixture without importing CrocSwap code.

    The construction is the field-trace symplectic spread, followed
    by the source's fixed sequence of symplectic transvections.
    """
    gram = tuple(tuple(gf_trace(gf_multiply(1 << i, 1 << j))
                       for j in range(N)) for i in range(N))
    assert gram == (
        (0, 0, 0, 0, 0, 1),
        (0, 0, 0, 0, 1, 0),
        (0, 0, 0, 1, 0, 0),
        (0, 0, 1, 0, 0, 0),
        (0, 1, 0, 0, 0, 0),
        (1, 0, 0, 0, 0, 1),
    )

    def dual(x):
        return sum(
            (sum(gram[i][j] * ((x >> j) & 1)
                 for j in range(N)) & 1) << i
            for i in range(N)
        )

    classes = []
    for slope in range(1 << N):
        classes.append([(1 << i) | (
            dual(gf_multiply(slope, 1 << i)) << N
        ) for i in range(N)])
    classes.append([1 << (N+i) for i in range(N)])

    result = []
    for g in classes:
        transformed = []
        for x in g:
            for t in SOURCE_GAUGE_TRANSVECTIONS:
                if symplectic(x, t):
                    x ^= t
            transformed.append(x)
        result.append(reduced_binary_basis(transformed))

    assert len(result) == 65
    assert result[0] == SOURCE_FIRST_CLASS
    assert result[22] == SOURCE_MISSING_CLASS
    assert result[-1] == SOURCE_LAST_CLASS
    occupied = set()
    for generators in result:
        assert len(generators) == N
        assert all(not symplectic(a, b)
                   for a in generators for b in generators)
        labels = span(generators)-{0}
        assert len(labels) == 63
        assert not (occupied & labels)
        occupied.update(labels)
    assert occupied == set(range(1, 1 << (2*N)))
    print("CrocSwap source construction replayed: 65 maximal commuting "
          "classes, 4,095 distinct nonidentity Pauli labels")
    return result


def code_checks(class_generators, functional):
    """All 63 distinct rank-five subspaces inside one rank-six class.

    A binary hyperplane is kernel of nonzero linear functional on
    the six generator coordinates, NOT just deleting one generator.
    """
    assert 1 <= functional < 1 << N
    pivot = next(i for i in range(N) if (functional >> i) & 1)
    checks = tuple(
        class_generators[i] ^ (
            class_generators[pivot] if (functional >> i) & 1 else 0
        ) for i in range(N) if i != pivot
    )
    assert len(span(checks)) == 32
    assert all(not symplectic(a, b) for a in checks for b in checks)
    return checks


PAULI_ORDER = sorted(range(1, 1 << (2*N)), key=weight)


def distance(checks):
    """Smallest weight of Pauli commuting with all checks but not
    in their stabilizer span. Exact symplectic stabilizer distance.
    """
    s = span(checks)
    for p in PAULI_ORDER:
        if p not in s and all(not symplectic(p, g) for g in checks):
            return weight(p), p
    raise AssertionError("no logical Pauli despite rank-five stabilizer")


def single_qubit_error_set():
    errors = [0]
    for i in range(N):
        for x, z in ((1, 0), (0, 1), (1, 1)):
            errors.append((x << i) | (z << (N+i)))
    assert len(set(errors)) == 19
    return errors


def exact_one_error_recovery(checks):
    """Exact Knill-Laflamme for the full span of single-qubit errors.

    P E†F P can be nonzero ONLY when the difference label e^f
    commutes with the checks. In that case it must belong to the
    stabilizer group (up to Pauli phase), hence is scalar on P.
    The same-syndrome errors are degenerate but correctable.
    """
    code_stabilizers = span(checks)
    errors = single_qubit_error_set()
    syndromes = defaultdict(list)
    for e in errors:
        signature = tuple(symplectic(e, s) for s in checks)
        syndromes[signature].append(e)

    pairs = 0
    degeneracies = []
    for e in errors:
        for f in errors:
            difference = e ^ f
            same_syndrome = all(
                symplectic(difference, s) == 0 for s in checks
            )
            if same_syndrome:
                assert difference in code_stabilizers, (
                    "Failed KL: undetectable, nonstabilizer E†F",
                    pauli_label(e), pauli_label(f)
                )
                if e != f:
                    degeneracies.append((pauli_label(e), pauli_label(f)))
            pairs += 1
    assert pairs == 19**2

    # Syndrome-based recovery: choose one representative per group.
    # The correction acts trivially on the logical qubit up to a
    # stabilizer and a physically irrelevant Pauli phase.
    for signature, candidates in syndromes.items():
        representative = candidates[0]
        assert all((representative ^ e) in code_stabilizers
                   for e in candidates)
    assert len(syndromes) <= 19
    print(f"KL error-pair audit: {pairs} exact relations; "
          f"{len(syndromes)} occupied syndromes; "
          f"{len(degeneracies)//2} degenerate unordered error pairs")
    return syndromes, degeneracies


def test():
    spread = reconstruct_source_spread()

    # Each full class has six independent checks on six qubits:
    # 64/64 = 1-dimensional eigenspaces, encoding no logical qubit.
    assert all(len(span(g)) == 64 for g in spread)
    census = Counter()
    distance3 = []
    for class_index, generators in enumerate(spread):
        for hyperplane in range(1, 1 << N):
            checks = code_checks(generators, hyperplane)
            d, logical = distance(checks)
            census[d] += 1
            if d >= 3:
                distance3.append((class_index, hyperplane, checks, logical))
    assert sum(census.values()) == 65*63 == 4095
    assert census == Counter({1: 1539, 2: 2539, 3: 17})
    assert len(distance3) == 17
    print("All 4,095 rank-five stabilizer subspaces:"
          " [[6,1,1]] 1,539; [[6,1,2]] 2,539; [[6,1,3]] 17")
    print("17 refers to distinct subspaces in this fixture,"
          " NOT inequivalent code families.")

    # Explicit source-derived, distance-three code.
    selected = code_checks(spread[CHOSEN_CLASS], CHOSEN_HYPERPLANE)
    assert selected == CHOSEN_CHECKS
    assert distance(selected) == (3, 70)
    assert len(span(selected)) == 32
    assert all(not symplectic(a, b) for a in selected for b in selected)
    # Every generator is Hermitian by canonical Pauli phase choice;
    # five independent commuting +1 checks project onto a 2D space:
    # rank(P)=64/2^5=2 from Pauli trace orthogonality.
    syndrome_table, degeneracy = exact_one_error_recovery(selected)
    assert len(syndrome_table) == 18
    assert len(degeneracy) == 2
    logical_x, logical_z = 70, 148
    assert weight(logical_x) == weight(logical_z) == 3
    assert symplectic(logical_x, logical_z) == 1
    assert all(not symplectic(v, g)
               for v in (logical_x, logical_z) for g in selected)
    assert logical_x not in span(selected)
    assert logical_z not in span(selected)
    print("Representative [[6,1,3]] checks:",
          [pauli_label(v) for v in selected])
    print("Logical anticommuting operators:",
          pauli_label(logical_x), pauli_label(logical_z))
    print("Verified full single-site Pauli-noise recovery algebra.")
    print("No physical gauge implementation, anomaly-code relation,"
          " or improvement over the published five-qubit code claimed.")
    print("ALL CROCSWAP PAULI-QEC BRIDGE CHECKS PASSED")


if __name__ == "__main__":
    test()
