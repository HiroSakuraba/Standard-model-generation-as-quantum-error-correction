"""Paper-inspired adaptive flag EC and an exact distance-three criteria audit.

Implements the DECISION TREE of Chao & Reichardt, PRL 121, 050502
(2018): measure successive stabilizers with one flagged gadget each;
at first observed nontrivial flag or syndrome, measure all four
syndromes with *unflagged* gadgets, then decode using the flag-specific
error class (if flagged) or ordinary weight-one lookup.

IMPORTANT provenance: this reuses OUR independently analyzed
controlled-Pauli flag gadget from flagged_su2_measurement.py, NOT a
gate-for-gate transcription of Fig. 2(c) of arXiv:1705.02329v1 from Chao--Reichardt.
The encoded controls remain ideal; no SU(2) physical implementation.

It enumerates the same exact Pauli fault alphabet as the previous
research, injecting faults only at physically reached gate locations.
The runtime decoder sees measured ancilla bits/flags ONLY; ideal
syndromes are used solely in *post-run testing*. All results require
perfect recovery operations and at most one specified new fault.

A key scientific correction: a distance-three one-fault EC is NOT
obliged to preserve the original logical state when both one incoming
error AND one additional gadget fault occur (r+s=2>1). Such tests are
reported separately and not counted as failures of the one-fault
criterion. The source's extended-rectangle criteria are distinguished
in the accompanying note.

Run: python paper_style_flag_ec.py [--stress]
"""
from collections import Counter
from itertools import product

from flagged_su2_measurement import (
    CHECKS, DATA_COUNT, FLAG, LOGICAL_X, LOGICAL_Z, MEASURE,
    commute_bit, conjugate, fault_injections, from_word, label,
    multiply, single, standard_decoder, syndrome_circuit, syndrome_data,
)
from measured_multiround_su2 import (
    MASK, NO_ERROR, classify_residual, flag_recovery_maps,
)


def extract(data, check_index, flagged, injected=None):
    """Exact ancilla-Pauli propagation through an observable check gadget.

    Fresh m ancilla is initialized |+>, measured X; fresh f is |0>,
    measured Z if flagged. 3-tuple injection is a single fault AFTER a
    selected gate, at preparation, or at readout.
    """
    gates = syndrome_circuit(CHECKS[check_index], flagged,
                            flag_after=(1, 3))
    result = data
    for step in range(len(gates) + 1):
        if injected is not None and injected[2] == step:
            result = multiply(result, injected[1])
        if step < len(gates):
            result = conjugate(result, gates[step])
    observed_syndrome_bit = (result[1] >> MEASURE) & 1
    observed_flag = ((result[0] >> FLAG) & 1) if flagged else 0
    remaining = (result[0] & MASK, result[1] & MASK)
    return remaining, observed_syndrome_bit, observed_flag


def run(input_error=NO_ERROR, fault=None, maps=None):
    """Paper-style adaptive schedule, with NO oracle data syndrome.

    fault is (phase, stabilizer_index, injection) where phase is
    'flagged' or 'unflagged'. Only realized injection locations
    are marked 'used'. An injected fault on an unreachable branch
    is not counted by the test. The decoder uses only returned
    measurement outcomes and prerecorded offline lookup tables.
    """
    if maps is None:
        maps = flag_recovery_maps()
    decoder = standard_decoder()
    data = input_error
    trace = {"flagged": [], "unflagged": [], "used": False, "trigger": None}

    for idx, check in enumerate(CHECKS):
        inject = (fault[2] if fault is not None and
                  fault[0] == "flagged" and fault[1] == idx else None)
        if inject is not None:
            trace["used"] = True
        data, bit, flagged = extract(data, idx, True, inject)
        trace["flagged"].append((idx, bit, flagged))

        if not (bit or flagged):
            continue

        trace["trigger"] = (idx, "flag" if flagged else "syndrome")
        readings = []
        # The unflagged full round is actually executed, NOT an oracle.
        for j in range(len(CHECKS)):
            followup_injection = (
                fault[2] if fault is not None and
                fault[0] == "unflagged" and fault[1] == j else None
            )
            if followup_injection is not None:
                trace["used"] = True
            data, b, unused_flag = extract(
                data, j, False, followup_injection)
            assert unused_flag == 0
            readings.append(b)
        observed_syndrome = tuple(readings)
        trace["unflagged"] = list(readings)
        trace["followup_syndrome"] = observed_syndrome

        if flagged:
            correction = maps[idx].get((1, observed_syndrome))
            if correction is None:
                trace["decoder_undefined"] = (idx, observed_syndrome)
                return data, trace
            trace["decoder"] = "flag-specific"
        else:
            correction = decoder[observed_syndrome]
            trace["decoder"] = "standard"
        data = multiply(data, correction)
        trace["correction"] = label(correction)
        return data, trace

    # No nontrivial syndrome/flag: leave the current data as-is.
    trace["decoder"] = "no-correction"
    return data, trace


def all_single_faults(phase, index):
    gates = syndrome_circuit(CHECKS[index], phase == "flagged",
                            flag_after=(1, 3))
    faults = list(fault_injections(gates, phase == "flagged"))
    assert len(faults) == (98 if phase == "flagged" else 64)
    return ((phase, index, fault) for fault in faults)


def verify_one_fault_criterion():
    """Exhaust clean input plus one new gadget fault, and correct r=1,s=0."""
    maps = flag_recovery_maps()
    data, log = run(maps=maps)
    assert data == NO_ERROR and log["trigger"] is None
    assert log["unflagged"] == []

    # The 15 incoming one-module errors and no new faults MUST be
    # recovered without changing the encoded logical state.
    incoming = Counter()
    for i in range(DATA_COUNT):
        for op in "XYZ":
            data, log = run(input_error=single(i, op), maps=maps)
            assert log["trigger"] is not None
            outcome = classify_residual(data)
            assert outcome == "fully_corrected", (i, op, outcome, log)
            incoming[outcome] += 1
    assert incoming == {"fully_corrected": 15}

    cases = Counter()
    witnesses = []
    for index in range(len(CHECKS)):
        for fault in all_single_faults("flagged", index):
            data, trace = run(fault=fault, maps=maps)
            assert trace["used"], ("unreached fault location", fault)
            assert "decoder_undefined" not in trace, (fault, trace)
            outcome = classify_residual(data)
            cases[outcome] += 1
            cases["trigger" if trace["trigger"] else "no_trigger"] += 1
            if (trace["trigger"] and outcome != "fully_corrected") or (
                    outcome == "logical_failure"):
                witnesses.append((index, fault[2][0], outcome, trace))

    assert sum(cases[s] for s in
               ("fully_corrected", "one_error_residual",
                "logical_failure")) == 4 * 98
    assert not witnesses, witnesses[:6]
    assert cases["logical_failure"] == 0
    print("Chao-Reichardt-inspired adaptive decision schedule; "
          "independently designed controlled-Pauli flag gadget.")
    print("15 incoming module errors, no new faults: 15 fully corrected")
    print("392 clean-input single-fault cases:", dict(cases))
    print("All triggered cases fully recovered; untriggered cases "
          "have at most one correctable module residual.")
    print("The decoder used recorded flag/syndrome bits, never oracle syndromes.")
    return cases


def stress_two_independent_faults():
    """Diagnostic: r=1 incoming + s=1 circuit fault, r+s=2>t=1.

    This is NOT a required preserve-original-logical criterion for
    distance-three one-fault EC. The problem is deliberately counted
    separately, including single faults on conditional followup gates
    that become reachable when the incoming error triggers a check.
    """
    maps = flag_recovery_maps()
    counts = Counter()
    examples = []
    for i, op in product(range(DATA_COUNT), "XYZ"):
        incoming = single(i, op)
        for phase in ("flagged", "unflagged"):
            for j in range(len(CHECKS)):
                for fault in all_single_faults(phase, j):
                    data, transcript = run(
                        input_error=incoming, fault=fault, maps=maps)
                    if not transcript["used"]:
                        continue
                    if "decoder_undefined" in transcript:
                        outcome = "decoder_undefined"
                    else:
                        outcome = classify_residual(data)
                    counts[outcome] += 1
                    if outcome in ("logical_failure",
                                   "decoder_undefined") and len(examples) < 3:
                        examples.append(
                            (i, op, phase, j, fault[2][0],
                             outcome, transcript))
    assert sum(counts.values()) > 0
    print("NON-CRITERION diagnostic (incoming data error plus gadget fault):",
          dict(counts))
    print("First such case histories:", examples)
    print("r+s=2 exceeds distance-three one-fault preservation requirement.")
    return counts


def main():
    verify_one_fault_criterion()
    if "--stress" in __import__("sys").argv:
        stress_two_independent_faults()
    print("ADAPTIVE FLAG SCHEDULE CHECKS COMPLETED")


if __name__ == "__main__":
    main()
