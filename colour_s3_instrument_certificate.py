"""Certify the FULL pointer-ready finite-S3 recovery instrument, not just F_e.

The code and electric syndrome each occupy rank four inside the 14D
physical Gauss sector; the rank-six orthogonal complement MUST be flagged.
All 14 inputs are tested without postselection. Ordered gR/Rg diagnostics
are deliberately separate from the original seven single-fault channel.

This is a *verification* tool, NOT a native decoder synthesis or a QCD
simulation. Projectors are used to compute diagnostic overlaps, never
silently credited as a physical correction step.
"""
from __future__ import annotations

import numpy as np

from colour_s3_recovery_pulse_compiler import (
    DIM, POINTER, TOTAL, target_three_outcome,
    single_fault_logical_recovery_fidelity,
)
from colour_s3_exact_electric_commutant import exact_generators
from colour_s3_order_six_active_electric import (
    closure_matrices, gauge_fixed_normalized_matrix,
)


def _logical_outcomes(outputs, code):
    """Unconditional channel ledger for four orthonormal logical inputs.

    Correct pointer outcomes 0/1 retaining the right code are credited.
    Flag 2 and off-code states are NOT renormalized into an apparent win.
    """
    split = outputs.reshape(DIM, POINTER, 4)
    kraus = [code.conj().T @ split[:, b, :] for b in (0, 1)]
    fe = float(sum(abs(np.trace(k) / 4) ** 2 for k in kraus))
    captured = float(sum(np.linalg.norm(k, "fro") ** 2 / 4 for k in kraus))
    flagged = float(np.linalg.norm(split[:, 2, :], "fro") ** 2 / 4)
    missed = float(1 - captured - flagged)
    if min(fe, captured, flagged, missed) < -2e-9:
        raise AssertionError("nonphysical probability ledger")
    return {
        "unconditional_logical_entanglement_F_e": fe,
        "recognized_pointer_and_code_probability": captured,
        "pointer2_flag_probability": flagged,
        "unflagged_outside_code_probability": max(0., missed),
        "probabilities_sum_to_one": float(captured + flagged + missed),
    }


def certify_three_outcome(U, label, p_electric=.37,
                          pointer_threshold=.99, transfer_threshold=.99):
    """Assess code(4), electric(4), and flag(6) on all 14 input directions.

    A passing *mathematical* instrument certificate requires each
    branch's intended pointer probability AND its branch transfer
    fidelity to clear thresholds. It does not certify locally native
    charged-carrier transport, preparation, measurement, or reset.

    The 42D U must be supplied explicitly (e.g. a product of admissible
    pulse Hamiltonians). Never substitute the ideal target silently.
    """
    U = np.asarray(U, dtype=complex)
    if U.shape != (TOTAL, TOTAL):
        raise ValueError("provide an explicit 42x42 system-pointer operator")
    if not np.all(np.isfinite(U)):
        raise ValueError("nonfinite system-pointer operator")
    if np.linalg.norm(U.conj().T @ U - np.eye(TOTAL), "fro") > 2e-8:
        raise ValueError("system-pointer operator is not unitary")
    for t in (pointer_threshold, transfer_threshold):
        if not np.isfinite(t) or not 0 <= t <= 1:
            raise ValueError("threshold must lie in [0,1]")
    if not np.isfinite(p_electric) or not 0 <= p_electric <= 1:
        raise ValueError("electric probability must lie in [0,1]")

    _, ideal, _, _, branches = target_three_outcome(label)
    ptr0 = np.eye(POINTER)[:, [0]]
    expected_labels = {"code": 0, "electric": 1, "flag": 2}
    code = branches["code"]
    sector = {"code": code, "electric": code, "flag": branches["flag"]}
    results = {}
    for name in ("code", "electric", "flag"):
        basis = branches[name]
        n = basis.shape[1]
        initial = np.kron(basis, ptr0)
        actual = U @ initial
        target = ideal @ initial
        split = actual.reshape(DIM, POINTER, n)
        ptr = np.sum(np.abs(split) ** 2, axis=(0, 2)) / n
        correct = expected_labels[name]
        in_expected_sector = float(
            np.linalg.norm(sector[name].conj().T @ split[:, correct, :],
                           "fro") ** 2 / n
        )
        branch_transfer = float(abs(np.vdot(target, actual) / n) ** 2)
        results[name] = {
            "input_rank": n,
            "intended_pointer_label": correct,
            "pointer_outcome_probabilities": [float(x) for x in ptr],
            "intended_pointer_probability": float(ptr[correct]),
            "intended_pointer_and_system_sector_probability": in_expected_sector,
            "branch_coherent_transfer_fidelity": branch_transfer,
            "branch_target_residual_fro_per_sqrt_rank": float(
                np.linalg.norm(target - actual, "fro") / np.sqrt(n)),
        }
        if abs(sum(ptr) - 1) > 2e-9:
            raise AssertionError("missing pointer probability")

    classes, ops, sizes, _ = __import__(
        "colour_s3_recovery_pulse_compiler",
        fromlist=["controls_for_group"]).controls_for_group(label)
    A, B, S, E = ops
    group = closure_matrices((A, S) if label == "AS" else (B, S))
    R = np.eye(DIM) - 2 * gauge_fixed_normalized_matrix(E, sizes)
    ordered = []
    for index, g in enumerate(group):
        gphys = gauge_fixed_normalized_matrix(g, sizes)
        for name, fault in (("gR", gphys @ R), ("Rg", R @ gphys)):
            final = U @ np.kron(fault @ code, ptr0)
            ordered.append({"exchange_index": index, "fault_order": name,
                            **_logical_outcomes(final, code)})

    seven = single_fault_logical_recovery_fidelity(U, label, p_electric)
    passed = all(
        rec["intended_pointer_probability"] >= pointer_threshold - 1e-12
        and rec["branch_coherent_transfer_fidelity"] >= transfer_threshold - 1e-12
        for rec in results.values()
    )
    return {
        "group": label,
        "physical_Gauss_sector_dimension": DIM,
        "pointer_dimension": POINTER,
        "tested_pointer_ready_input_directions": 14,
        "branch_partition": {"code": 4, "electric": 4, "flag": 6},
        "branches": results,
        "mathematical_three_outcome_instrument_pass": bool(passed),
        "required_pointer_probability_threshold": pointer_threshold,
        "required_coherent_transfer_threshold": transfer_threshold,
        "independent_seven_Kraus_single_fault_channel": seven,
        "ordered_two_fault_diagnostics": ordered,
        "all_syndrome_outcomes_unconditional": True,
        "projectors_used_only_as_measurement_diagnostics": True,
        "local_hardware_decoder_compiled": False,
        "native_device_calibration_and_noise_included": False,
        "scope": "finite S3 Gauss-14D toy system, NOT continuum SU3 or QCD",
    }


def main():
    """Fast deterministic identity/ideal controls; no optimizer hidden inside."""
    import json
    _, ideal, _, _, _ = target_three_outcome("AS")
    for name, candidate in (
        ("do_nothing", np.eye(TOTAL, dtype=complex)),
        ("uncompiled_global_ideal_upper_bound", ideal),
    ):
        outcome = certify_three_outcome(candidate, "AS")
        print(json.dumps({
            "policy": name,
            "passes_complete_three_outcome_transfer":
                outcome["mathematical_three_outcome_instrument_pass"],
            "logical_Fe": outcome[
                "independent_seven_Kraus_single_fault_channel"][
                    "unconditional_seven_fault_mixture_logical_F_e"],
            "branch_pointer_correct": {
                k: v["intended_pointer_probability"]
                for k, v in outcome["branches"].items()},
        }, sort_keys=True))


if __name__ == "__main__":
    main()
