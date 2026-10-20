"""Regressions for full three-outcome instrument certification."""
import unittest
import numpy as np

from colour_s3_recovery_pulse_compiler import (
    TOTAL, SpectralCompiler, target_three_outcome,
)
from colour_s3_instrument_certificate import certify_three_outcome


class ThreeOutcomeCertificateTests(unittest.TestCase):
    def test_global_target_passes_all_14_directions_both_subgroups(self):
        for label in ("AS", "BS"):
            _, ideal, _, _, _ = target_three_outcome(label)
            r = certify_three_outcome(ideal, label)
            self.assertTrue(r["mathematical_three_outcome_instrument_pass"])
            self.assertEqual(r["branch_partition"], {
                "code": 4, "electric": 4, "flag": 6})
            self.assertEqual(r["tested_pointer_ready_input_directions"], 14)
            self.assertTrue(r["projectors_used_only_as_measurement_diagnostics"])
            self.assertFalse(r["local_hardware_decoder_compiled"])
            for name, row in r["branches"].items():
                self.assertAlmostEqual(
                    row["intended_pointer_probability"], 1., places=9)
                self.assertAlmostEqual(
                    row["intended_pointer_and_system_sector_probability"],
                    1., places=9)
                self.assertAlmostEqual(
                    row["branch_coherent_transfer_fidelity"], 1., places=9)
                self.assertAlmostEqual(
                    sum(row["pointer_outcome_probabilities"]), 1., places=9)
            self.assertAlmostEqual(
                r["independent_seven_Kraus_single_fault_channel"][
                    "unconditional_seven_fault_mixture_logical_F_e"],
                1., places=9)

    def test_identity_baseline_fails_electric_and_flag_requirements(self):
        for label in ("AS", "BS"):
            r = certify_three_outcome(np.eye(TOTAL), label)
            self.assertFalse(r["mathematical_three_outcome_instrument_pass"])
            self.assertAlmostEqual(
                r["branches"]["electric"]["intended_pointer_probability"],
                0., places=10)
            self.assertAlmostEqual(
                r["branches"]["flag"]["intended_pointer_probability"],
                0., places=10)
            self.assertAlmostEqual(
                r["independent_seven_Kraus_single_fault_channel"][
                    "unconditional_seven_fault_mixture_logical_F_e"],
                .63 + .37 / 9, places=9)

    def test_exact_target_recovers_Rg_but_not_every_gR(self):
        for label in ("AS", "BS"):
            _, ideal, _, _, _ = target_three_outcome(label)
            data = certify_three_outcome(ideal, label)
            ordered = data["ordered_two_fault_diagnostics"]
            self.assertEqual(len(ordered), 12)
            for row in ordered:
                self.assertAlmostEqual(
                    row["probabilities_sum_to_one"], 1., places=9)
            for row in ordered:
                if row["fault_order"] == "Rg":
                    self.assertAlmostEqual(
                        row["unconditional_logical_entanglement_F_e"],
                        1., places=9)
            gR = [x for x in ordered if x["fault_order"] == "gR"]
            self.assertEqual(sum(
                abs(x["unconditional_logical_entanglement_F_e"] -
                    11 / 96) < 1e-8 for x in gR), 4)
            self.assertEqual(sum(
                abs(x["unconditional_logical_entanglement_F_e"] - 1.) < 1e-8
                for x in gR), 2)
            for x in gR:
                if x["unconditional_logical_entanglement_F_e"] < .5:
                    self.assertAlmostEqual(
                        x["pointer2_flag_probability"], 13 / 24, places=8)

    def test_admissible_pulse_product_is_audited_without_projection_credit(self):
        comp = SpectralCompiler("AS", 1)
        angles = np.random.default_rng(6719).normal(
            0, .15, len(comp.order))
        U = comp.propagate(angles, input_states=np.eye(TOTAL))
        record = certify_three_outcome(U, "AS")
        self.assertFalse(record["mathematical_three_outcome_instrument_pass"])
        self.assertEqual(len(record["ordered_two_fault_diagnostics"]), 12)
        self.assertTrue(record["all_syndrome_outcomes_unconditional"])
        self.assertFalse(record["native_device_calibration_and_noise_included"])
        for rec in record["ordered_two_fault_diagnostics"]:
            self.assertGreaterEqual(
                rec["unconditional_logical_entanglement_F_e"], -1e-10)
            self.assertLessEqual(
                rec["unconditional_logical_entanglement_F_e"], 1+1e-10)

    def test_rejects_invalid_operator_probabilities_and_thresholds(self):
        with self.assertRaises(ValueError):
            certify_three_outcome(np.eye(14), "AS")
        bad = np.eye(TOTAL); bad[0, 0] = 2
        with self.assertRaises(ValueError):
            certify_three_outcome(bad, "AS")
        with self.assertRaises(ValueError):
            certify_three_outcome(np.eye(TOTAL), "AS",
                                  p_electric=float("inf"))
        with self.assertRaises(ValueError):
            certify_three_outcome(np.eye(TOTAL), "AS",
                                  transfer_threshold=1.01)
        with self.assertRaises(ValueError):
            certify_three_outcome(np.eye(TOTAL), "AS",
                                  pointer_threshold=float("nan"))


if __name__ == "__main__":
    unittest.main()
