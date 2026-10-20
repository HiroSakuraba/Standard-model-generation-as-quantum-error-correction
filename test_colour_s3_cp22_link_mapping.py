"""Exact 6=2x3 link encoding, native CP22 cost and Gauss-path regression."""
import unittest

import numpy as np

from colour_gauge_link_s3 import REP
from colour_s3_cp22_link_mapping import (
    _parity, parity_cycle_order, cycle_fourier,
    standard_projector_diagonal_audit, pointer_electric_cp22_mapping,
    electric_link_local_phases, same_colour_cp22_mapping,
    onsite_swap_cp22_mapping, controlled_sum_via_cp22,
    link_selector_gauss_audit, report, targeted_cp22,
)


class CP22FiniteS3MappingTests(unittest.TestCase):
    def test_parity_cycle_link_factorization_and_central_projector(self):
        order=parity_cycle_order()
        self.assertEqual(len(set(order)),6)
        self.assertEqual([_parity(x) for x in order],[1,1,1,-1,-1,-1])
        f=cycle_fourier()
        self.assertLess(np.linalg.norm(f.conj().T@f-np.eye(3)),2e-12)
        E,D,U,error=standard_projector_diagonal_audit()
        self.assertLess(error,2e-12)
        self.assertEqual(int(round(np.trace(E).real)),4)
        self.assertLess(np.linalg.norm(U@E@U.conj().T-D),2e-12)

    def test_two_native_selective_cp22_exact_for_arbitrary_angles(self):
        for angle in (0.,.17,np.pi/7,np.pi/2,np.pi,2*np.pi):
            row=pointer_electric_cp22_mapping(angle)
            self.assertEqual(row["selective_CP22_count"],2)
            self.assertLess(row["unitary_residual"],2e-11)
            self.assertLess(np.linalg.norm(
                row["original_group_basis_unitary"].conj().T@
                row["original_group_basis_unitary"]-np.eye(18)),2e-11)
            self.assertTrue(row["link_parity_register_idle_during_entangler"])
            self.assertTrue(row["link_fourier_is_encoding_not_free_hardware_pulse"])

    def test_link_electric_no_conditional_cp22_but_requires_local_phases(self):
        for angle in (0.,.25,np.pi):
            row=electric_link_local_phases(angle)
            self.assertLess(row["residual"],1e-11)
            self.assertEqual(row["conditional_two_qutrit_CP22_count"],0)
            self.assertEqual(row["single_cycle_qutrit_level_phases"],2)

    def test_onsite_equal_colour_three_cp22_strict_gauge_final(self):
        for angle in (.01,.3,np.pi):
            row=same_colour_cp22_mapping(angle)
            self.assertEqual(row["selective_CP22_count"],3)
            self.assertLess(row["residual"],1e-11)
            for h in REP:
                action=np.kron(h,h)
                self.assertLess(np.linalg.norm(
                    row["unitary"]@action-action@row["unitary"]),2e-11)
            self.assertFalse(row[
                "individual_selective_CP22_factors_gauge_invariant"])
            # One colour-specific phase does NOT commute with all permutations.
            partial=targeted_cp22(0,0,angle)
            self.assertGreater(max(np.linalg.norm(
                partial@np.kron(h,h)-np.kron(h,h)@partial)
                for h in REP),.001)

    def test_onsite_swap_twelve_cp22_upper_bound(self):
        row=onsite_swap_cp22_mapping()
        self.assertLess(row["residual"],2e-10)
        self.assertEqual(row["selective_CP22_count"],12)
        self.assertEqual(row["single_qutrit_fourier_count"],6)
        for h in REP:
            U=np.kron(h,h)
            self.assertLess(np.linalg.norm(
                row["unitary"]@U-U@row["unitary"]),2e-10)
        # Every modular SUM gate is native-factor-decomposed in the
        # optimizer but need not individually commute with S3.
        plus=controlled_sum_via_cp22(0,1,+1)
        for a in range(3):
            for b in range(3):
                src=np.eye(9)[:,3*a+b]
                expected=np.eye(9)[:,3*a+(a+b)%3]
                self.assertLess(np.linalg.norm(plus@src-expected),2e-11)

    def test_one_selective_cp22_breaks_one_vertex_gauss_while_product_restores(self):
        row=link_selector_gauss_audit(.7)
        self.assertTrue(row["single_CP22_intermediate_symmetry_breaking"])
        self.assertTrue(row["complete_product_Gauss_invariant"])
        self.assertGreater(row["per_vertex"]["B"][
            "maximum_single_mode_projector_commutator"],.5)
        for v in ("A","B"):
            self.assertLess(row["per_vertex"][v][
                "maximum_complete_central_projector_commutator"],2e-11)
            self.assertAlmostEqual(row["physical_Gauss_acceptance"][v][
                "after_complete_electric_rotation"],1.,places=9)
        self.assertLess(row["physical_Gauss_acceptance"]["B"][
            "after_one_selective_phase"],.99999)

    def test_reject_bad_angles_and_validate_reproducible_summary(self):
        with self.assertRaises(ValueError):
            pointer_electric_cp22_mapping(float("nan"))
        with self.assertRaises(ValueError):
            same_colour_cp22_mapping(float("inf"))
        with self.assertRaises(ValueError):
            targeted_cp22(3,2,.5)
        with self.assertRaises(ValueError):
            controlled_sum_via_cp22(0,0)
        doc=report()
        self.assertEqual(len(doc["exact_gate_mappings"]),3)
        self.assertEqual(doc["onsite_qutrit_SWAP_ideal_CP22_upper_bound"],12)
        self.assertGreater(len(doc["explicitly_NOT_compiled"]),2)
        self.assertIn("not measured hardware gate counts",
                      " ".join(doc["limitations"]).lower())


if __name__=="__main__":
    unittest.main()
