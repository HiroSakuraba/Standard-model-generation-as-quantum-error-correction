"""Primary-source Figure 2(d) seven flagged errors, a reference regression.

Chao & Reichardt, Quantum Error Correction with Only Two Extra Qubits,
Phys Rev Lett 121, 050502 (2018), arXiv:1705.02329v1 p. 2 Figure 2(d).

Figure 2(c) is the published circuit; our earlier ideal-encoded
controlled-P flag gadgets use a DIFFERENT ancilla basis and do not
claim to be a literal transcription of that figure.

The paper lists four faults after gate b and four after gate c, with
IIZXI appearing twice; thus SEVEN distinct data errors. We check all
have distinct nonzero syndromes, explaining the paper's correction
lookup. We do NOT claim to have simulated the paper's physical gates,
or to establish full device-level fault tolerance from this catalog.
Run: python published_flag_figure2_catalog.py
"""
from flagged_su2_measurement import from_word, syndrome_data

# Manually transcribed after visually inspecting arXiv 1705.02329v1,
# p.2, Fig.2(d), which uses XZZXI for the first stabilizer.
PAPER_FAULT_TO_DATA = {
    ("b", "IZ"): "IIZXI",
    ("b", "XZ"): "IXZXI",
    ("b", "YZ"): "IYZXI",
    ("b", "ZZ"): "IZZXI",
    ("c", "IZ"): "IIIXI",
    ("c", "XZ"): "IIXXI",
    ("c", "YZ"): "IIYXI",
    ("c", "ZZ"): "IIZXI",
}


def test():
    assert len(PAPER_FAULT_TO_DATA) == 8
    distinct = set(PAPER_FAULT_TO_DATA.values())
    assert distinct == {
        "IIZXI", "IXZXI", "IYZXI", "IZZXI",
        "IIIXI", "IIXXI", "IIYXI",
    }
    assert len(distinct) == 7
    keys = {}
    for word in distinct:
        syndrome = syndrome_data(from_word(word))
        assert syndrome != (0, 0, 0, 0)
        assert syndrome not in keys, (
            "Published flagged errors are not distinguishable", word,
            keys.get(syndrome))
        keys[syndrome] = word
    assert len(keys) == 7
    print("Chao-Reichardt 2018 Fig.2(d): 8 fault rows, 7 distinct "
          "nontrivial data errors, 7 unique nonzero 4-bit syndromes.")
    print("Verified source catalog:", sorted(keys.items()))
    print("Not a gate-level reproduction of Fig.2(c) or hardware SU2 gadgets.")
    print("ALL PUBLISHED FLAG FIGURE2 CATALOG CHECKS PASSED")


if __name__ == "__main__":
    test()
