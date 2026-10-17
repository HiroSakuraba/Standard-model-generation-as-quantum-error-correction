"""Two observed flagged rounds plus adaptive third: exact single-fault proof.

Unlike flagged_su2_measurement.py, this recovery algorithm DOES NOT read
the simulator's ideal data syndrome. It uses ONLY measured check bits,
measured flags, and a lookup derived from the previously verified
single-gadget fault model.

Protocol (no preexisting data error; at most one listed fault):
  - Measure all four checks in a flagged round, twice.
  - If ANY flag or any difference in the four reported syndrome bits:
      measure a THIRD complete flagged round. With just one overall
      fault, this third round is fault-free. Decode with its measured
      bits and the first two rounds' flag information.
  - Otherwise decode using round 2's measured syndrome.

On accepting without a third round, a fault in the last round may leave
a one-module correctable residual. This is a 1-fault containment
property, NOT proof of unconditional zero-error output. Exact full
correction is obtained in the triggered third-round cases.

Faults include each nonidentity two-body Pauli after each CNOT/CZ,
Pauli preparation faults on the ancillas, and classical readout flips.
Exactly ONE fault anywhere in the first 2 x 4 gadgets (8 x 98 = 784),
plus the no-fault case and 15 preexisting one-module errors with zero
circuit faults, are verified.

No faults in the conditional third round, no preexisting+new faults,
no recovery-gate faults, no idle errors, and no physically compiled
SU(2) control gates are asserted. These are critical limitations.

Run: python measured_multiround_su2.py
"""
from collections import Counter
from itertools import product

from flagged_su2_measurement import (
    CHECKS, DATA_COUNT, FLAG, MEASURE, LOGICAL_X, LOGICAL_Z,
    class_after_simple_decoder, conjugate, commute_bit, fault_injections, from_word, label,
    multiply, recovery_classes, simulate, single, standard_decoder,
    syndrome_circuit, syndrome_data,
)


MASK = (1 << DATA_COUNT) - 1
NO_ERROR = (0, 0)


def gadget(data, check_index, injected=None):
    """Run exactly one flagged check with fresh ancillas.

    injected is one (description, Pauli error, first subsequent gate)
    from fault_injections. Output is obtained from the *ancilla* Pauli
    labels at the end of the circuit. The true data syndrome is not
    consulted here or by the recovery algorithm.
    """
    gates = syndrome_circuit(CHECKS[check_index], True, flag_after=(1, 3))
    state = data  # ancillas start in their ideal +X_m and +Z_f states.
    for step in range(len(gates) + 1):
        if injected is not None and injected[2] == step:
            state = multiply(state, injected[1])
        if step < len(gates):
            state = conjugate(state, gates[step])
    m_bit = (state[1] >> MEASURE) & 1  # X measurement flip
    f_bit = (state[0] >> FLAG) & 1     # Z measurement flip
    remaining_data = (state[0] & MASK, state[1] & MASK)
    return remaining_data, m_bit, f_bit


def full_round(data, injected=None):
    """Return updated data, four *recorded* syndrome bits, four flags."""
    results = []
    flags = []
    for j in range(4):
        fault = injected[1] if injected is not None and injected[0] == j else None
        data, m, f = gadget(data, j, fault)
        results.append(m)
        flags.append(f)
    return data, tuple(results), tuple(flags)


def flag_recovery_maps():
    """Precompute decisions from the original exact gadget enumeration.

    This is offline decoder construction, not a true-syndrome oracle
    at runtime. The lookup's KEY uses a third-round measured syndrome.
    """
    maps = {}
    for i, check in enumerate(CHECKS):
        corrections, ambiguities = recovery_classes(simulate(
            check, True, flag_after=(1, 3)))
        assert not ambiguities
        maps[i] = {
            key: from_word(pauli)
            for key, pauli in corrections.items()
            if key[0] == 1
        }
    return maps


def recover(data, round_one, round_two, third_round, maps):
    """Select correction from observed outcomes only.

    third_round is a callback which performs a *real*, noiseless
    measurement round on the current data when needed. This
    assumption is justified solely by the stated one-fault budget:
    reaching this callback implies that a fault already occurred
    in round one/two for a clean starting codeword.
    """
    data1, b1, f1 = round_one
    data2, b2, f2 = round_two
    assert data == data2
    decoder = standard_decoder()
    flagged_indices = [i for i in range(4) if f1[i] or f2[i]]
    triggered = bool(flagged_indices or b1 != b2)
    if triggered:
        data, observed_syndrome, observed_flags = third_round(data)
        # For the stated fault model a third-round fault is impossible
        # because the first fault already caused a trigger.
        assert observed_flags == (0, 0, 0, 0)
        if flagged_indices:
            assert len(flagged_indices) == 1
            which = flagged_indices[0]
            assert (1, observed_syndrome) in maps[which], (
                which, observed_syndrome)
            correction = maps[which][(1, observed_syndrome)]
        else:
            correction = decoder[observed_syndrome]
        source = "third_measured"
    else:
        correction = decoder[b2]
        source = "second_measured"
    return multiply(data, correction), triggered, source


def run_protocol(injection=None, preexisting=NO_ERROR, maps=None):
    """Injection = (round_index, gadget_index, fault_tuple), or None."""
    if maps is None:
        maps = flag_recovery_maps()
    fault1 = (injection[1], injection[2]) if injection is not None and injection[0] == 0 else None
    fault2 = (injection[1], injection[2]) if injection is not None and injection[0] == 1 else None
    r1 = full_round(preexisting, fault1)
    r2 = full_round(r1[0], fault2)
    third_calls = []
    def take_third(data):
        result = full_round(data)
        third_calls.append(result)
        return result
    after, triggered, source = recover(r2[0], r1, r2, take_third, maps)
    assert len(third_calls) == int(triggered)
    return after, {
        "measured_1": r1[1],
        "measured_2": r2[1],
        "flags_1": r1[2],
        "flags_2": r2[2],
        "triggered": triggered,
        "decision_source": source,
        "measured_3": third_calls[0][1] if triggered else None,
    }


def classify_residual(data):
    """Verifier-only use of TRUE syndrome. Never call in recover()."""
    decoder = standard_decoder()
    s = syndrome_data(data)
    remaining = multiply(data, decoder[s])
    logical = (
        commute_bit(remaining, from_word(LOGICAL_Z)),
        commute_bit(remaining, from_word(LOGICAL_X)),
    )
    return ("fully_corrected" if s == (0, 0, 0, 0) and logical == (0, 0)
            else "one_error_residual" if logical == (0, 0)
            else "logical_failure")


def enumerate_faults():
    for r, j in product(range(2), range(4)):
        gates = syndrome_circuit(CHECKS[j], True, flag_after=(1, 3))
        faults = list(fault_injections(gates, True))
        assert len(faults) == 98
        for fault in faults:
            yield (r, j, fault)


def test():
    decoder_maps = flag_recovery_maps()
    assert all(m for m in decoder_maps.values())

    # No fault anywhere: the complete two-round protocol terminates
    # with a no-error correction and no third round.
    after, trace = run_protocol(maps=decoder_maps)
    assert after == NO_ERROR and not trace["triggered"]

    # The two-round protocol without any gate/ancilla faults must also
    # correct an arbitrary *preexisting* single-module Pauli error.
    for i in range(DATA_COUNT):
        for op in "XYZ":
            after, trace = run_protocol(
                preexisting=single(i, op), maps=decoder_maps)
            assert classify_residual(after) == "fully_corrected"
            assert not trace["triggered"]
            assert trace["measured_1"] == trace["measured_2"]

    counts = Counter()
    failure_examples = []
    for fault in enumerate_faults():
        after, trace = run_protocol(injection=fault, maps=decoder_maps)
        status = classify_residual(after)
        counts[status] += 1
        counts["triggered" if trace["triggered"] else "accepted_two_rounds"] += 1
        if trace["triggered"]:
            # Here the third round was perfect, so the flagged decoder
            # should correct the data fully, not merely leave a single.
            if status != "fully_corrected":
                failure_examples.append((fault[0], fault[1], fault[2][0], status, trace))
        elif status not in ("fully_corrected", "one_error_residual"):
            failure_examples.append((fault[0], fault[1], fault[2][0], status, trace))

    assert sum(counts[s] for s in
               ("fully_corrected", "one_error_residual", "logical_failure")) == 784
    assert not failure_examples, failure_examples[:5]
    assert counts["logical_failure"] == 0
    assert counts["one_error_residual"] > 0
    assert counts["triggered"] > 0
    assert counts["accepted_two_rounds"] > 0
    print("Two flagged complete four-check rounds + adaptive third")
    print("0 circuit faults, preexisting one-module errors: 15/15 fully corrected")
    print("784 exact single-fault cases across 8 gadgets:", dict(counts))
    print("All triggered fault cases: fully corrected by actual third-round measurements")
    print("Untriggered fault cases: fully corrected OR one correctable residual")
    print("No perfect data syndrome is used by the runtime decoder")
    print("ALL MEASURED MULTIROUND SU2 CHECKS PASSED")



def stress_preexisting_plus_fault():
    """Diagnostic only: one preexisting module error PLUS one new fault.

    This is more demanding than the certified one-fault-from-clean
    experiment. Some flags and syndrome pairs may lack valid lookup
    entries, or the resulting residual may be a logical fault.
    Neither case is silently called success.
    """
    decoder_maps = flag_recovery_maps()
    counts = Counter()
    examples = []
    cases = tuple(enumerate_faults())
    seen_observations = {}
    intrinsically_ambiguous = set()
    witness = []
    ordinary_decoder = standard_decoder()
    for i in range(DATA_COUNT):
        for op in "XYZ":
            prior = single(i, op)
            for injection in cases:
                # Diagnostic-only reconstruction of the data Pauli and
                # measured transcript, independent of the runtime decoder.
                injected_one = ((injection[1], injection[2])
                                if injection[0] == 0 else None)
                injected_two = ((injection[1], injection[2])
                                if injection[0] == 1 else None)
                first = full_round(prior, injected_one)
                second = full_round(first[0], injected_two)
                trigger = (any(first[2]) or any(second[2]) or
                           first[1] != second[1])
                third_measured = (full_round(second[0])[1]
                                  if trigger else None)
                # Include even a TRUE final syndrome in the lookup key.
                # A collision with different logical classes here proves
                # insufficiency of the whole measured transcript EVEN
                # IF an ideal true syndrome were additionally supplied.
                key = (first[1], first[2], second[1], second[2],
                       third_measured, syndrome_data(second[0]))
                logical = class_after_simple_decoder(
                    second[0], ordinary_decoder)
                description = (i, op, injection[0], injection[1],
                               injection[2][0], label(second[0]))
                previous = seen_observations.get(key)
                if previous is None:
                    seen_observations[key] = (logical, description)
                elif previous[0] != logical:
                    intrinsically_ambiguous.add(key)
                    if len(witness) < 4:
                        witness.append((key, previous, (logical, description),
                                        "same full observed transcript and true syndrome"))
                try:
                    after, trace = run_protocol(
                        injection=injection,
                        preexisting=prior,
                        maps=decoder_maps,
                    )
                    outcome = classify_residual(after)
                except (AssertionError, KeyError, ValueError) as exc:
                    outcome = "decoder_undefined"
                    trace = {"reason": str(exc)[:180]}
                counts[outcome] += 1
                if outcome in ("logical_failure", "decoder_undefined") and len(examples) < 5:
                    examples.append((i, op, injection[0], injection[1],
                                     injection[2][0], outcome, trace))
    assert sum(counts.values()) == DATA_COUNT * 3 * len(cases) == 11760
    print("STRESS: one preexisting module error PLUS one circuit fault:")
    print(dict(counts))
    print("First uncorrectable/undefined examples:", examples)
    assert intrinsically_ambiguous and witness, (
        "The current stress fixture no longer establishes an "
        "intrinsic classical-transcript ambiguity; re-analyze it.")
    print("Distinct transcript+true-syndrome groups with conflicting "
          "logical error classes:", len(intrinsically_ambiguous))
    print("First ambiguity witnesses:", witness[:2])
    print("This is a DIAGNOSTIC, not a claim of fault tolerance.")
    return counts


if __name__ == "__main__":
    test()
    if "--stress" in __import__("sys").argv[1:]:
        stress_preexisting_plus_fault()
