"""Branch-aware recovery optimizer tests, without fictitious projections."""
import unittest

import numpy as np

from colour_s3_recovery_pulse_compiler import TOTAL
from colour_s3_branch_aware_optimizer import (
    BranchAwareCompiler, compile_branches, finite_gradient_audit,
)


class BranchAwareSearchTests(unittest.TestCase):
    def test_rank_partition_and_exact_ideal_upper_bound(self):
        for label in ("AS", "BS"):
            c = BranchAwareCompiler(label, 1)
            self.assertEqual(
                [(n, rank) for n, _, rank in c.spans],
                [("code", 4), ("electric", 4), ("flag", 6)])
            self.assertEqual(c.initial.shape, (TOTAL, 14))
            loss, adjoint, data = c._terminal(c.target)
            self.assertAlmostEqual(loss, 0., places=10)
            self.assertEqual(adjoint.shape, (TOTAL, 14))
            for score in data["branch_coherent_transfer"].values():
                self.assertAlmostEqual(score, 1., places=10)
            self.assertAlmostEqual(
                data["dimension_weighted_branch_transfer"], 1.,
                places=10)

    def test_identity_cannot_pass_full_instrument_certificate(self):
        for label in ("AS", "BS"):
            c = BranchAwareCompiler(label, 1)
            angles = np.zeros(c.n_angles)
            direct, independent = c.independently_audit(angles)
            self.assertAlmostEqual(
                direct["dimension_weighted_branch_transfer"],
                4 / 14, places=9)
            self.assertAlmostEqual(
                direct["branch_coherent_transfer"]["code"], 1., places=9)
            for name in ("electric", "flag"):
                self.assertAlmostEqual(
                    direct["branch_coherent_transfer"][name], 0., places=9)
            self.assertFalse(
                independent["mathematical_three_outcome_instrument_pass"])

    def test_independent_diagnostics_match_for_real_pulse_product(self):
        for label in ("AS", "BS"):
            c = BranchAwareCompiler(label, 1)
            a = np.random.default_rng(31).normal(
                0, .2, c.n_angles)
            direct, cert = c.independently_audit(a)
            for name in ("code", "electric", "flag"):
                self.assertAlmostEqual(
                    direct["branch_coherent_transfer"][name],
                    cert["branches"][name][
                        "branch_coherent_transfer_fidelity"], places=9)
                self.assertLessEqual(
                    direct["branch_coherent_transfer"][name],
                    cert["branches"][name][
                        "intended_pointer_probability"] + 1e-9)
            self.assertFalse(cert["local_hardware_decoder_compiled"])

    def test_finite_difference_gradients_with_weakest_branch_penalty(self):
        for label in ("AS", "BS"):
            for cycles in (1, 2):
                r = finite_gradient_audit(label, cycles)
                self.assertLess(r["max_absolute_error"], 4e-7)
                self.assertEqual(len(r["checks"]), 4)

    def test_invalid_configuration_rejected(self):
        for params in (
            {"min_transfer":1.01}, {"min_transfer":float("nan")},
            {"branch_penalty":-1.}, {"branch_penalty":float("inf")},
            {"effort_lambda":-1.}, {"effort_epsilon":0.},
        ):
            with self.assertRaises(ValueError):
                BranchAwareCompiler("AS", 1, **params)
        with self.assertRaises(ValueError):
            compile_branches(maxiter=0)
        with self.assertRaises(ValueError):
            compile_branches(maxiter=2, warmstart=np.zeros(13))

    def test_short_bounded_search_has_true_ungated_baseline(self):
        r = compile_branches("AS", cycles=1, maxiter=4, seed=67)
        self.assertLessEqual(
            r["selected_objective"],
            r["initial_no_gate_objective"]+1e-10)
        self.assertEqual(len(r["angles_rad"]), 13)
        self.assertEqual(len(r["attempts"]), 3)
        self.assertTrue(r["no_ideal_projector_applied_as_control_operation"])
        self.assertTrue(
            r["no_physical_device_model_or_hardware_success_claim"])
        self.assertFalse(r["no_uncompiled_feedback_gate_counted"] is False)
        self.assertEqual(
            r["all_three_branch_acceptance_passed"],
            r["complete_instrument_certificate"][
                "mathematical_three_outcome_instrument_pass"])
        self.assertGreaterEqual(
            r["optimistic_seven_fault_logical_Fe_with_ideal_code_projection"],
            0.)
        self.assertEqual(len(
            r["complete_instrument_certificate"][
                "ordered_two_fault_diagnostics"]), 12)


if __name__ == "__main__":
    unittest.main()
