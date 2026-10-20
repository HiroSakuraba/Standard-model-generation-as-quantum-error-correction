"""Deterministic checks for sector-capture and 13-to-104 continuation."""
import unittest
import numpy as np

from colour_s3_continuation_search import (
    SectorCaptureCompiler, extend_angles, verify_sector_gradient,
    search_depth, run_continuation,
)
from colour_s3_branch_aware_optimizer import BranchAwareCompiler


class SectorCaptureContinuationTests(unittest.TestCase):
    def test_ideal_target_captures_all_14_directions(self):
        for label in ("AS", "BS"):
            c = SectorCaptureCompiler(label)
            loss, back, diag = c._capture_and_adjoint(c.target)
            self.assertAlmostEqual(loss, 0., places=10)
            self.assertEqual(back.shape, (42, 14))
            self.assertAlmostEqual(diag[
                "dimension_weighted_capture_probability"], 1., places=10)
            self.assertEqual(
                set(diag["per_branch_correct_sector_and_pointer_probability"]),
                {"code", "electric", "flag"})

    def test_zero_pulse_identity_exposes_electric_and_flag_failure(self):
        for label in ("AS", "BS"):
            c = SectorCaptureCompiler(label)
            report = c.diagnostics(np.zeros(c.n_angles))
            b = report["per_branch_correct_sector_and_pointer_probability"]
            self.assertAlmostEqual(b["code"], 1., places=10)
            self.assertAlmostEqual(b["electric"], 0., places=10)
            self.assertAlmostEqual(b["flag"], 0., places=10)
            self.assertAlmostEqual(
                report["dimension_weighted_capture_probability"],
                4/14, places=10)

    def test_analytic_sector_capture_gradient(self):
        for label in ("AS", "BS"):
            for cycles in (1, 2):
                r = verify_sector_gradient(label, cycles)
                self.assertLess(r["max_absolute_error"], 4e-7)

    def test_sector_capture_is_upper_bound_on_coherent_transfer(self):
        c = SectorCaptureCompiler("AS", 2)
        a = np.random.default_rng(45).normal(0, .25, c.n_angles)
        capture = c.diagnostics(a)[
            "per_branch_correct_sector_and_pointer_probability"]
        coherent = c.coherent.transfer_diagnostics(a)[
            "branch_coherent_transfer"]
        for branch in ("code", "electric", "flag"):
            self.assertGreaterEqual(capture[branch] + 1e-10,
                                    coherent[branch])

    def test_identity_padding_preserves_entire_pulse_product(self):
        for label in ("AS", "BS"):
            small = BranchAwareCompiler(label, 1)
            large = BranchAwareCompiler(label, 4)
            a = np.random.default_rng(7).normal(0, .25, small.n_angles)
            padded = extend_angles(a, large.n_angles)
            small_unitary = small.base.propagate(
                a, input_states=np.eye(42))
            large_unitary = large.base.propagate(
                padded, input_states=np.eye(42))
            self.assertLess(np.linalg.norm(
                small_unitary-large_unitary, "fro"), 2e-10)
            for name in ("code", "electric", "flag"):
                self.assertAlmostEqual(
                    small.transfer_diagnostics(a)[
                        "branch_coherent_transfer"][name],
                    large.transfer_diagnostics(padded)[
                        "branch_coherent_transfer"][name], places=9)

    def test_small_continuation_never_loses_prior_worst_branch_score(self):
        result = run_continuation(
            depths=(1, 2), labels=("AS",),
            capture_iter=3, coherent_iter=4, seed=77)
        first, second = result["results"]
        self.assertEqual([first["pulse_count"],second["pulse_count"]],
                         [13, 26])
        self.assertGreaterEqual(
            second["worst_coherent_branch_transfer"] + 2e-9,
            first["worst_coherent_branch_transfer"])
        self.assertTrue(
            all(x["all_three_branches_audited_without_postselection"]
                for x in (first, second)))
        self.assertTrue(
            all(x["no_ideal_code_projection_inserted_as_pulse"]
                for x in (first, second)))
        self.assertFalse(
            first["full_physical_feedback_and_reset_implemented"])
        self.assertEqual(len(first["independent_instrument_certificate"][
            "ordered_two_fault_diagnostics"]), 12)
        self.assertEqual(len(second["optimization_trajectory"]), 4)
        self.assertEqual(len(second["candidate_metrics"]), 12)
        self.assertEqual(result["depth_cycles"], [1, 2])

    def test_parameter_validation(self):
        with self.assertRaises(ValueError):
            SectorCaptureCompiler(min_capture=1.5)
        with self.assertRaises(ValueError):
            SectorCaptureCompiler(balance_weight=-1.)
        with self.assertRaises(ValueError):
            SectorCaptureCompiler(effort_lambda=float("inf"))
        with self.assertRaises(ValueError):
            verify_sector_gradient(step=0.)
        with self.assertRaises(ValueError):
            extend_angles(np.zeros(13), 13)
        with self.assertRaises(ValueError):
            extend_angles(np.full(2, float("nan")), 13)
        with self.assertRaises(ValueError):
            search_depth("AS", 1, capture_iter=0)
        with self.assertRaises(ValueError):
            run_continuation(depths=(1, 1), labels=("AS",),
                             capture_iter=1, coherent_iter=1)


if __name__ == "__main__":
    unittest.main()
