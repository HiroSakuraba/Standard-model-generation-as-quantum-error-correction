"""Regression checks for temporal fault ordering, global syndrome and T1 cost."""
import unittest
import numpy as np

from colour_s3_exact_electric_commutant import exact_generators
from colour_s3_order_six_active_electric import (
    closure_matrices,construct_active_code,
)
from colour_s3_correlated_recovery_cadence import (
    COMPOUND_LOGICAL_LOSS, BAD_COMPOUND_F_E,
    BAD_COMPOUND_CAPTURE,branch_kraus,corrected_two_round,
    global_gauge_preserving_instrument,exact_order_diagnostic,
    matched_channel,passive_group_trivial_best,loss_survival,
    lifetime_cadence_comparison,study,
)


class OrderedCorrelatedFaultTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        classes,_,(A,B,S,E)=exact_generators()
        cls.cases=[]
        sizes=[len(x) for x in classes]
        for label,gens in (("AS",(A,S)),("BS",(B,S))):
            G=closure_matrices(gens)
            C,V,_=construct_active_code(label,G,E,sizes)
            cls.cases.append((label,C,V,G,E,sizes))

    def test_one_round_exact_good_and_bad_compound_fractions(self):
        self.assertAlmostEqual(COMPOUND_LOGICAL_LOSS,85/144,places=14)
        self.assertAlmostEqual(BAD_COMPOUND_F_E,11/96,places=14)
        self.assertAlmostEqual(BAD_COMPOUND_CAPTURE,11/24,places=14)
        for label,C,V,G,E,sizes in self.cases:
            rows=exact_order_diagnostic(C,V,G,E,sizes)
            self.assertEqual(len(rows),12)
            forward=[r for r in rows if "gR" in r["physical_fault_order"]]
            reverse=[r for r in rows if "Rg" in r["physical_fault_order"]]
            self.assertEqual(len(forward),6)
            self.assertEqual(len(reverse),6)
            self.assertEqual(sum(r["one_round"]["unconditional_F_e"]<.2
                                 for r in forward),4)
            self.assertTrue(all(abs(r["one_round"]["unconditional_F_e"]-1)<1e-10
                                for r in reverse))

    def test_no_postselection_for_bad_ordering(self):
        for label,C,V,G,E,sizes in self.cases:
            rows=exact_order_diagnostic(C,V,G,E,sizes)
            bad=[r for r in rows if r["one_round"]["unconditional_F_e"]<.2]
            self.assertEqual(len(bad),4)
            for r in bad:
                self.assertAlmostEqual(r["one_round"]["unconditional_F_e"],
                                       11/96,places=10)
                self.assertAlmostEqual(r["one_round"]["captured_probability"],
                                       11/24,places=10)
                self.assertAlmostEqual(r["one_round"]["flagged_probability"],
                                       13/24,places=10)

    def test_interleaving_two_recoveries_corrects_each_order(self):
        for label,C,V,G,E,sizes in self.cases:
            rows=exact_order_diagnostic(C,V,G,E,sizes)
            for r in rows:
                self.assertAlmostEqual(r["two_round"]["unconditional_F_e"],1.,
                                       places=10)
                self.assertAlmostEqual(r["two_round"]["captured_probability"],1.,
                                       places=10)

    def test_coherent_pointer_recovery_unitary_is_explicit_not_local(self):
        for label,C,V,G,E,sizes in self.cases:
            result=global_gauge_preserving_instrument(C,V)
            self.assertEqual(result["joint_unitary_dimension"],28)
            self.assertEqual(result["pointer_levels"],2)
            self.assertLess(result["full_logical_code_path_residual"],1e-10)
            self.assertLess(result["full_error_syndrome_path_residual"],1e-10)
            self.assertTrue(result["only_global_projector_implementation_shown"])

    def test_bernoulli_fault_model_is_cptp_and_exact_polynomial(self):
        for label,C,V,G,E,sizes in self.cases:
            for e,g,f in ((0.,0.,0.),(.37,.4,0.),
                          (.37,.4,.5),(.37,.4,1.),(1.,1.,1.)):
                result=matched_channel(C,V,G,E,sizes,e,g,f)
                self.assertAlmostEqual(result["all_20_events_probability_sum"],1.,
                                       places=12)
                self.assertAlmostEqual(
                    result["single_round_conditional_F_e"],
                    1-(85/144)*e*g*f,places=10)
                self.assertTrue(0<=result["single_round_flagged_probability"]<=1+1e-10)
                self.assertEqual(len(result["enumerated_event_rows"]),20)

    def test_order_asymmetry_is_not_just_error_count(self):
        label,C,V,G,E,sizes=self.cases[0]
        reverse=matched_channel(C,V,G,E,sizes,.37,.4,0.)
        forward=matched_channel(C,V,G,E,sizes,.37,.4,1.)
        self.assertAlmostEqual(reverse["single_round_conditional_F_e"],1.,
                               places=10)
        self.assertLess(forward["single_round_conditional_F_e"],.92)

    def test_passive_four_state_fixed_sector_optimal_electric_trace(self):
        for label,C,V,G,E,sizes in self.cases:
            result=passive_group_trivial_best(C,G,E,sizes)
            self.assertAlmostEqual(
                result["best_conditional_F_e_for_electric_event"],
                4/9,places=10)
            self.assertTrue(result["no_claim_globally_optimal_passive_gauge_code"])

    def test_two_round_lifetime_threshold_and_timing_fairness(self):
        label,C,V,G,E,sizes=self.cases[0]
        channel=matched_channel(C,V,G,E,sizes,.37,.4,1.)
        passive=passive_group_trivial_best(C,G,E,sizes)
        low=lifetime_cadence_comparison(channel,passive,.37,1.,10.)
        high=lifetime_cadence_comparison(channel,passive,.37,1.,100.)
        self.assertFalse(low["two_active_rounds_win_over_one_under_this_loss_only_model"])
        self.assertTrue(high["two_active_rounds_win_over_one_under_this_loss_only_model"])
        cross=low["break_even_T1_us_two_rounds_vs_one"]
        self.assertGreater(cross,10.)
        self.assertLess(cross,100.)
        self.assertAlmostEqual(loss_survival(1.,30.),np.exp(-4/30),places=12)

    def test_perfect_error_order_does_not_justify_extra_correction(self):
        label,C,V,G,E,sizes=self.cases[0]
        channel=matched_channel(C,V,G,E,sizes,.37,.4,0.)
        passive=passive_group_trivial_best(C,G,E,sizes)
        report=lifetime_cadence_comparison(channel,passive,.37,1.,30.)
        self.assertFalse(report["two_active_rounds_win_over_one_under_this_loss_only_model"])
        self.assertIsNone(report["break_even_T1_us_two_rounds_vs_one"])

    def test_probability_and_lifetime_inputs_reject_nonphysical_values(self):
        label,C,V,G,E,sizes=self.cases[0]
        for p in (-.2,1.01,float("nan")):
            with self.assertRaises(ValueError):
                matched_channel(C,V,G,E,sizes,p,.4,.5)
            with self.assertRaises(ValueError):
                matched_channel(C,V,G,E,sizes,.37,.4,p)
        with self.assertRaises(ValueError):
            loss_survival(1.,0.)
        with self.assertRaises(ValueError):
            loss_survival(-1.,30.)

    def test_both_codes_remain_finite_gauge_model_and_no_native_claim(self):
        result=study()
        self.assertEqual(len(result["comparisons"]),2)
        self.assertTrue(result["ideal_two_round_recovery_compound_faults"])
        self.assertFalse(result["hardware_native_recovery_compiled"])
        for case in result["comparisons"]:
            self.assertEqual(len(case["ordered_fault_unitary_audit"]),12)
            self.assertEqual(len(case["two_time_window_stochastic_scenarios"]),5)
            self.assertAlmostEqual(
                case["passive_G_fixed_four_dim_optimum"][
                    "best_conditional_F_e_for_electric_event"],4/9,places=10)

if __name__=="__main__":
    unittest.main()
