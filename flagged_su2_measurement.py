"""Exact single-fault propagation through one five-block SU(2) stabilizer check.

The data 'qubits' are the two-dimensional color-singlet modules established
by four_link_su2_subspace.py; five_block_su2_recovery.py proves their outer
[[5,1,3]] code and identifies module X/Z from SU(2)-invariant exchanges.

Ideal gadget for each weight-four stabilizer:
  1. Prepare syndrome ancilla |+> and flag |0>.
  2. Controlled-P from syndrome to four distinct data modules.
  3. For a flagged gadget, CNOT syndrome -> flag after data gates 1 and 3.
  4. Measure syndrome in X, flag in Z.

This file exactly tracks all Pauli faults inserted AFTER each 2-qubit gate
(15 nonidentity Paulis on the gate's two supports), plus single ancilla
preparation Pauli faults and readout bit flips. Classical syndrome lookup
uses a subsequently PERFECT four-check syndrome, not the unreliable current
gadget's reported syndrome. This is NOT an end-to-end fault-tolerant circuit
or an implementation of controlled SU(2) exchange on actual hardware.

All calculations are mod-2 symplectic and exact. No dependencies.
Run: python flagged_su2_measurement.py
"""
from collections import defaultdict
from itertools import product

CHECKS = ("XZZXI", "IXZZX", "XIXZZ", "ZXIXZ")
LOGICAL_X, LOGICAL_Z = "XXXXX", "ZZZZZ"
DATA_COUNT, MEASURE, FLAG, N = 5, 5, 6, 7


def from_word(word):
    x = z = 0
    for j, ch in enumerate(word):
        assert ch in "IXYZ"
        if ch in "XY":
            x |= 1 << j
        if ch in "ZY":
            z |= 1 << j
    return x, z


def label(pauli, n=DATA_COUNT):
    x, z = pauli
    return "".join(
        "Y" if x & z & (1 << i) else
        "X" if x & (1 << i) else
        "Z" if z & (1 << i) else "I"
        for i in range(n)
    )


def single(i, op):
    return from_word("I" * i + op + "I" * (N - i - 1))


def multiply(p, q):
    # Global Pauli phases are irrelevant to syndrome and KL classes.
    return p[0] ^ q[0], p[1] ^ q[1]


def commute_bit(a, b):
    x1, z1 = a
    x2, z2 = b
    return ((x1 & z2).bit_count() + (z1 & x2).bit_count()) & 1


def syndrome_data(data):
    return tuple(commute_bit(data, from_word(g)) for g in CHECKS)


def standard_decoder():
    result = {(0, 0, 0, 0): (0, 0)}
    for i in range(DATA_COUNT):
        for op in "XYZ":
            error = single(i, op)
            syndrome = syndrome_data(error)
            assert syndrome not in result
            result[syndrome] = error
    assert len(result) == 16
    return result


def syndrome_circuit(check, flagged, flag_after=(1, 3), order=None):
    """Controlled-P gates implemented abstractly by CX/CZ on inner qubits.

    Ancilla acts as computational-basis control, initially |+>.
    Two CNOTs to |0> flag uncompute in the no-fault circuit.
    """
    positions = [i for i, c in enumerate(check) if c != "I"]
    assert len(positions) == 4
    if order is not None:
        assert sorted(order) == sorted(positions)
        positions = list(order)
    gates = []
    if flagged and 0 in flag_after:
        gates.append(("CX", MEASURE, FLAG))
    for j, i in enumerate(positions, 1):
        assert check[i] in "XZ"
        gates.append(("CX" if check[i] == "X" else "CZ", MEASURE, i))
        if flagged and j in flag_after:
            gates.append(("CX", MEASURE, FLAG))
    assert not flagged or len(flag_after) == 2
    return gates


def conjugate(pauli, gate):
    kind, a, b = gate
    x, z = pauli
    A, B = 1 << a, 1 << b
    if kind == "CX":
        if x & A:
            x ^= B
        if z & B:
            z ^= A
    elif kind == "CZ":
        if x & A:
            z ^= B
        if x & B:
            z ^= A
    else:
        raise ValueError(kind)
    return x, z


def propagate(pauli, gates, start):
    for gate in gates[start:]:
        pauli = conjugate(pauli, gate)
    return pauli


def fault_injections(gates, flagged):
    # A Pauli fault just after each elementary two-body gate.
    for i, (_, a, b) in enumerate(gates):
        for op_a, op_b in product("IXYZ", repeat=2):
            if op_a == op_b == "I":
                continue
            yield (f"after-{i}:{op_a}@{a}-{op_b}@{b}",
                   multiply(single(a, op_a), single(b, op_b)), i + 1)

    # Ancilla preparation faults, including unwanted phases.
    for j in ([MEASURE, FLAG] if flagged else [MEASURE]):
        for op in "XYZ":
            yield (f"prep-{op}@{j}", single(j, op), 0)

    # A classical readout bit flip is represented by an anticommuting
    # Pauli immediately before the intended X or Z basis measurement.
    yield ("readout-X-syndrome-flip", single(MEASURE, "Z"), len(gates))
    if flagged:
        yield ("readout-Z-flag-flip", single(FLAG, "X"), len(gates))


def class_after_simple_decoder(error, decoder):
    corrected = multiply(error, decoder[syndrome_data(error)])
    assert syndrome_data(corrected) == (0, 0, 0, 0)
    return (commute_bit(corrected, from_word(LOGICAL_Z)),
            commute_bit(corrected, from_word(LOGICAL_X)))


def simulate(check, flagged, flag_after=(1, 3), order=None):
    gates = syndrome_circuit(check, flagged, flag_after, order)
    decoder = standard_decoder()
    rows = []
    for position, fault, next_gate in fault_injections(gates, flagged):
        full = propagate(fault, gates, next_gate)
        mask = (1 << DATA_COUNT) - 1
        data = full[0] & mask, full[1] & mask
        rows.append({
            "position": position,
            "data": label(data),
            "syndrome": syndrome_data(data),
            "flag": (full[0] >> FLAG) & 1,
            "syndrome_readout_flip": (full[1] >> MEASURE) & 1,
            "simple_decoder_logical": class_after_simple_decoder(data, decoder),
        })
    return rows


def recovery_classes(rows):
    """Test whether a flag + PERFECT LATER syndrome can select a correction.

    All faults with identical (flag, data syndrome) must have the same
    logical class after applying the standard single-qubit decoder.
    A representative error corrects the common class up to a stabilizer.
    """
    decoder = standard_decoder()
    groups = defaultdict(list)
    for row in rows:
        groups[(row["flag"], row["syndrome"])].append(row)
    ambiguity = {}
    corrections = {}
    for key, cases in groups.items():
        classes = {case["simple_decoder_logical"] for case in cases}
        if len(classes) > 1:
            ambiguity[key] = tuple(sorted(classes))
            continue
        if key[0] == 0 and classes == {(0, 0)}:
            correction = decoder[key[1]]
        else:
            correction = from_word(cases[0]["data"])
        for case in cases:
            residual = multiply(correction, from_word(case["data"]))
            assert syndrome_data(residual) == (0, 0, 0, 0)
            assert commute_bit(residual, from_word(LOGICAL_X)) == 0
            assert commute_bit(residual, from_word(LOGICAL_Z)) == 0
        corrections[key] = label(correction)
    return corrections, ambiguity


def verify():
    for g in CHECKS:
        plain = simulate(g, False)
        early_flag = simulate(g, True, flag_after=(0, 3))
        correctly_placed = simulate(g, True, flag_after=(1, 3))
        p_corrections, p_ambiguity = recovery_classes(plain)
        e_corrections, e_ambiguity = recovery_classes(early_flag)
        f_corrections, f_ambiguity = recovery_classes(correctly_placed)
        del p_corrections, e_corrections

        # Each 2-qubit gate has 15 possible nontrivial Pauli faults.
        # Bare: 4 data controls, 3 prep faults, 1 readout fault =64.
        # Flag: 4 data + 2 flag controls, 6 prep and 2 readout =98.
        assert len(plain) == 64
        assert len(early_flag) == len(correctly_placed) == 98
        dangerous = lambda rows: [
            r for r in rows if r["simple_decoder_logical"] != (0, 0)
        ]
        assert len(dangerous(plain)) == 12
        assert len(dangerous(correctly_placed)) == 12
        assert len([r for r in dangerous(correctly_placed)
                    if r["flag"] == 0]) == 0

        # Correct early flag detects the bad hook errors but does NOT
        # tell the decoder how to correct one collision group. Moving
        # first flag coupling AFTER data-1 resolves this ambiguity.
        assert p_ambiguity
        assert e_ambiguity and len(e_ambiguity) == 1
        assert not f_ambiguity
        assert len(f_corrections) == 21
        assert any(r["position"] == "readout-X-syndrome-flip" and
                   r["syndrome_readout_flip"] and not r["flag"] and
                   r["data"] == "IIIII" for r in correctly_placed)

        print(f"{g}: bare 64 single-fault cases, "
              "12 unflagged naive-decoder logical failures; "
              "poor placement 1 ambiguous correction group; "
              "optimized 98 cases, all 12 dangerous faults flagged, "
              "0 correction ambiguities given a trusted later syndrome")

    print("ALL FLAGGED SU2 STABILIZER GADGET CHECKS PASSED")


if __name__ == "__main__":
    verify()
