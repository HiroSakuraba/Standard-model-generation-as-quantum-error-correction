"""Regress exact 4D active error-correction advantage for S3 exchange subgroups."""
import unittest
import numpy as np
import sympy as sp

from colour_s3_exact_electric_commutant import (
    exact_generators, physical_orbit_isometry,
)
from colour_s3_order_six_active_electric import (
    LOGICAL, closure_matrices, ordered_pair_centralizer_basis,
    exact_restricted_algebra_certificate, exact_sector_spectrum,
    gauge_fixed_normalized_matrix, construct_active_code,
    verify_five_plus_one_kraus_recovery,
    compound_exchange_after_electric_diagnostic,
    common_gauss_physical_crosscheck,
    study,
)


class OrderSixActiveCodeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.classes,cls.lookup,cls.ops=exact_generators()
        cls.sizes=[len(cl) for cl in cls.classes]

    def groups(self):
        A,B,S,E=self.ops
        for label,gens in (("AS",(A,S)),("BS",(B,S))):
            yield label,closure_matrices(gens),E

    def test_both_exchange_groups_have_exactly_six_elements(self):
        for label,group,E in self.groups():
            self.assertEqual(len(group),6)
            for g in group:
                self.assertEqual(g.shape,(14,14))
                self.assertEqual(g.T*g,sp.eye(14))

    def test_exact_group_centralizer_has_52_orbit_indicators(self):
        for label,group,E in self.groups():
            basis=ordered_pair_centralizer_basis(group)
            self.assertEqual(len(basis),52)
            self.assertEqual(sum(basis,sp.zeros(14)),sp.ones(14,14))
            for X in basis:
                for g in group:
                    self.assertEqual(X*g,g*X)

    def test_exact_electric_constraints_rank_41_joint_commutant_11(self):
        for label,group,E in self.groups():
            cert=exact_restricted_algebra_certificate(label,group,E)
            self.assertEqual(cert["electric_constraint_exact_rational_rank"],41)
            self.assertEqual(cert["joint_error_commutant_exact_dimension"],11)
            self.assertTrue(cert["no_four_logical_noiseless_subsystem"])

    def test_exact_fixed_sector_spectrum_supports_rank4_compression(self):
        for label,group,E in self.groups():
            P,r=exact_sector_spectrum(group,E)
            self.assertEqual(sp.trace(P),6)
            self.assertEqual(r["exact_compressed_electric_eigenvalues_and_multiplicities"],
                             {"0":1,"1/3":1,"2/3":2,"1":2})
            self.assertEqual(P*P,P)

    def test_four_logical_isometry_electric_scalar_compression(self):
        for label,group,E in self.groups():
            C,V,record=construct_active_code(label,group,E,self.sizes)
            self.assertEqual(C.shape,(14,4))
            self.assertEqual(V.shape,(14,4))
            self.assertLess(np.linalg.norm(C.T@C-np.eye(4)),2e-11)
            self.assertLess(np.linalg.norm(C.T@V),2e-11)
            self.assertLess(np.linalg.norm(V.T@V-np.eye(4)),2e-11)
            self.assertLess(max(record["numerical_Knill_Laflamme_residuals"].values()),
                            2e-11)

    def test_all_seven_individual_error_Kraus_branches_perfectly_decode(self):
        for label,group,E in self.groups():
            C,V,_=construct_active_code(label,group,E,self.sizes)
            rows=verify_five_plus_one_kraus_recovery(C,V,group,E,self.sizes)
            self.assertEqual(len(rows),7)
            self.assertEqual({x["error"] for x in rows},
                             {*(f"exchange_{i}" for i in range(6)),
                              "electric_phase"})
            for r in rows:
                self.assertAlmostEqual(
                    r["after_ideal_CPTP_recovery_entanglement_fidelity"],
                    1.,places=10)
                self.assertAlmostEqual(r["decoded_probability"],1.,places=10)
                self.assertAlmostEqual(r["unrecoverable_flagged_probability"],
                                       0.,places=10)
                self.assertLess(r["max_logical_Kraus_nonscalar_residual"],2e-11)

    def test_electric_error_has_nonorthogonal_pre_measurement_branches(self):
        for label,group,E in self.groups():
            C,V,_=construct_active_code(label,group,E,self.sizes)
            r=verify_five_plus_one_kraus_recovery(C,V,group,E,self.sizes)[-1]
            self.assertAlmostEqual(r["code_syndrome_probability"],1/9,places=11)
            self.assertAlmostEqual(r["electric_syndrome_probability"],8/9,places=11)

    def test_sequential_two_fault_error_not_silently_claimed_corrected(self):
        for label,group,E in self.groups():
            C,V,_=construct_active_code(label,group,E,self.sizes)
            rows=compound_exchange_after_electric_diagnostic(
                C,V,group,E,self.sizes)
            self.assertEqual(len(rows),6)
            self.assertAlmostEqual(
                rows[0]["recovered_F_e_without_postselection"],1.,places=10)
            self.assertLess(min(x[
                "recovered_F_e_without_postselection"] for x in rows),.2)
            self.assertGreater(max(x[
                "flagged_outside_syndrome_probability"] for x in rows),.5)

    def test_extended_gauge_sector_physical_embedding_486D(self):
        for label,group,E in self.groups():
            C,V,_=construct_active_code(label,group,E,self.sizes)
            test=common_gauss_physical_crosscheck(
                self.classes,self.sizes,C,V,label)
            self.assertTrue(test[
                "logical_and_syndrome_both_independent_vertex_Gauss_invariant"])
            self.assertLess(test[
                "physical_logical_to_syndrome_overlap_norm"],2e-11)

    def test_active_code_has_no_passive_four_state_or_subsystem_equivalent(self):
        r=study()
        for case in r["results"]:
            self.assertEqual(case["passive_common_scalar_sector_max_dimension"],2)
            self.assertEqual(case["exact_no_4D_passive_noiseless_subsystem"][
                "joint_error_commutant_exact_dimension"],11)
            self.assertEqual(case["active_four_logical_code"]["logical_dimension"],4)

    def test_theorem_applies_to_individual_faults_not_compound_words(self):
        r=study()
        self.assertTrue(r["active_advantage_against_passive_at_abstract_KL_level"])
        self.assertFalse(r["full_instrument_native_hardware_compiled"])
        self.assertTrue(any("sequential gR compound faults" in x for x in r["limitations"]))
        self.assertEqual(len(r["results"]),2)
        self.assertTrue(r["virtual_intermediate_physics"][
            "Hubbard_previous_PR55_to_PR58"].startswith("explicit high-energy"))

if __name__=="__main__":
    unittest.main()
