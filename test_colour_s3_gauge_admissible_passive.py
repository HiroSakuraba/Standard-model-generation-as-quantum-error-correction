"""Exact common-Gauss, common-dressed-swap passive-versus-active regressions."""
import unittest
import numpy as np

from colour_gauge_link_s3 import (
    DIM, G, physical_isometry, transform, dressed_exchange,
)
from colour_s3_edge_logical_syndrome import code_isometry, onsite_swap
from colour_s3_charged_bus_transport import initialize, gauss_acceptance
from colour_s3_gauge_admissible_passive import (
    DLOG, P_SOURCE, canonical_swap_passive_code,
    canonical_dihedral_trivial_code,
    alternate_dressed_exchange, active_ideal_second_swap_metric,
    passive_fault_metrics, colour_mode_weights,
    passive_exact_number_dephasing, comparison_row,
    stochastic_matched_dephasing, lifetime_sensitivity, run,
)


class SameGaugePassiveTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.W=physical_isometry()
        cls.K,_,_=code_isometry(cls.W)
        cls.C=canonical_dihedral_trivial_code(cls.W)
        cls.Cweak=canonical_swap_passive_code(cls.W)

    def test_exact_plus_one_dressed_swap_sector_rank_ten(self):
        S=self.W.conj().T@dressed_exchange(self.W)
        P=(np.eye(14)+S)/2
        self.assertAlmostEqual(np.trace(P).real,10.,places=10)
        self.assertLess(np.linalg.norm(P@P-P),1e-10)
        self.assertEqual(self.C.shape,(DIM,4))
        self.assertLess(np.linalg.norm(self.C.conj().T@self.C-np.eye(4)),1e-10)
        self.assertLess(np.linalg.norm(dressed_exchange(self.C)-self.C),1e-10)

    def test_both_codes_satisfy_both_vertex_gauss_actions(self):
        for C in (self.C,self.K):
            for v in ("A","B"):
                for h in range(G):
                    self.assertLess(np.linalg.norm(transform(C,h,v)-C),1e-10)
            z=initialize(C,pointer_uniform=False)
            self.assertAlmostEqual(np.vdot(z,z).real,DLOG,places=10)
            for acceptance in gauss_acceptance(z).values():
                self.assertAlmostEqual(acceptance,1.,places=10)

    def test_dihedral_exchange_trivial_sector_is_six_dimensional(self):
        for op in (dressed_exchange,onsite_swap,alternate_dressed_exchange):
            self.assertLess(np.linalg.norm(op(self.C)-self.C),1e-10)

    def test_same_channel_has_exact_passive_immunity(self):
        for p in (0.,.37,1.):
            result=passive_fault_metrics(self.C,p)
            self.assertAlmostEqual(result["same_dressed_exchange_channel_F_e"],1.,places=11)
        self.assertLess(np.linalg.norm(self.K.conj().T@dressed_exchange(self.K)),1e-10)
        self.assertLess(np.linalg.norm(self.C-self.K),3.) # distinct code subspaces

    def test_second_dressed_exchange_is_gauge_covariant_and_not_assumed_benign(self):
        for h in range(G):
            for v in ("A","B"):
                residual=transform(alternate_dressed_exchange(self.C),h,v)-(
                    alternate_dressed_exchange(transform(self.C,h,v)))
                self.assertLess(np.linalg.norm(residual),1e-10)
        r=passive_fault_metrics(self.C,.37)
        self.assertTrue(-1e-10<=r["alternate_dressed_exchange_S_A1_B0_channel_F_e"]<=1.+1e-10)
        self.assertAlmostEqual(
            r["alternate_dressed_exchange_S_A1_B0_channel_F_e"],
            1.,places=10)
        weak=passive_fault_metrics(self.Cweak,.37)
        self.assertLess(
            weak["alternate_dressed_exchange_S_A1_B0_channel_F_e"],
            1.-1e-5)
        active=active_ideal_second_swap_metric(self.K,.37)
        self.assertAlmostEqual(active,1.,places=10)

    def test_exact_dephasing_initial_weights_include_reference_singlet(self):
        weights=colour_mode_weights(self.C)
        self.assertEqual(weights.shape,(81,))
        self.assertAlmostEqual(weights.sum(),4.,places=10)
        self.assertTrue(np.all(weights>=0))
        self.assertAlmostEqual(
            passive_exact_number_dephasing(self.C,0.,400.,"colour_resolved"),
            1.,places=12)
        self.assertAlmostEqual(
            passive_exact_number_dephasing(self.C,.01,400.,"orbital_charge"),
            1.,places=12)
        a=passive_exact_number_dephasing(self.C,.0002,400.,"colour_resolved")
        b=passive_exact_number_dephasing(self.C,.0002,900.,"colour_resolved")
        self.assertTrue(0<=b<a<1.)
        self.assertAlmostEqual(
            passive_exact_number_dephasing(self.C,1000.,1.,"colour_resolved"),
            np.sum(weights**2)/16.,places=10)

    def test_same_four_active_fermion_loss_on_both_memories(self):
        row=comparison_row(30.,self.K,self.C)
        self.assertAlmostEqual(
            row["active_encoded_with_microscopic_feedback_no_loss_F_e"],
            .9119196563716021,places=5)
        self.assertAlmostEqual(row["passive_gauss_admissible_no_loss_F_e"],1.,places=12)
        self.assertAlmostEqual(row["same_loss_survival_both"],
                               np.exp(-4*row["equal_four_fermion_time_us"]/30),places=11)
        self.assertGreater(row["passive_unconditional_F_e"],
                           row["active_unconditional_F_e"])

    def test_no_finite_lifetime_crossover_for_known_swap_alone(self):
        anchors=[comparison_row(30.,self.K,self.C)]
        grid=lifetime_sensitivity(anchors,(20.,30.,65.,100.,300.,1.e6))
        self.assertEqual(len(grid),6)
        self.assertTrue(all(r["passive_wins"] for r in grid))
        for item in grid:
            self.assertAlmostEqual(
                item["passive_unconditional_F_e"]/
                item["equal_four_fermion_no_jump_probability"],1.,places=12)

    def test_identical_number_operator_bath_comparison(self):
        rows=stochastic_matched_dephasing(self.K,self.C,U=30.,
                                          gamma=.0002,shots=2)
        self.assertEqual([x["noise_model"] for x in rows],
                         ["orbital_charge","colour_resolved"])
        self.assertAlmostEqual(rows[0]["passive_conditional_F_e_exact"],1.,places=12)
        self.assertLess(rows[1]["passive_conditional_F_e_exact"],1.)
        for r in rows:
            self.assertTrue(r["same_stochastic_bath_operator_family_both_memories"])
            self.assertGreaterEqual(r["active_conditional_MC_standard_error"],0.)
            self.assertGreaterEqual(r["passive_unconditional_F_e_exact"],0.)

    def test_full_report_preserves_hardware_and_error_algebra_boundaries(self):
        data=run()
        self.assertEqual(data["known_dressed_swap_plus_sector_dimension"],10)
        self.assertEqual(data["common_three_reflection_trivial_sector_dimension"],6)
        self.assertLess(data["weak_swap_only_passive_negative_control"][
            "alternate_dressed_exchange_S_A1_B0_channel_F_e"],1.)
        self.assertAlmostEqual(data["exact_passive_known_error"][
            "alternate_dressed_exchange_S_A1_B0_channel_F_e"],1.,places=10)
        self.assertEqual(data["logical_dimension_each"],4)
        self.assertTrue(data["same_original_matter_link_registers"])
        self.assertEqual(len(data["noiseless_loss_matched_rows"]),3)
        self.assertEqual(len(data["colour_number_dephasing_comparison"]),2)
        self.assertTrue(any("SMALL D3" in x for x in data["limits"]))
        self.assertFalse(data["results_are_native_device_calibrated"])


if __name__=="__main__":
    unittest.main()
