"""Direct channel-F_e optimization regressions, including no-worse-than-identity."""
import unittest
import numpy as np

from colour_s3_recovery_pulse_compiler import (
    target_three_outcome,SpectralCompiler,
    single_fault_logical_recovery_fidelity,
)
from colour_s3_direct_logical_compiler import (
    logical_fault_states,effort_weight_for_name,
    DirectLogicalCompiler,direct_gradient_check,
    compile_direct,benchmark,
)


class DirectLogicalFidelityTests(unittest.TestCase):
    def test_two_input_blocks_exactly_represent_seven_error_Kraus(self):
        for label in ("AS","BS"):
            d=logical_fault_states(label,.37)
            self.assertEqual(d["inputs"].shape,(42,8))
            self.assertEqual(d["readout_code_isometries"].shape,(2,42,4))
            self.assertEqual(d["individual_Kraus_errors_represented"],7)
            self.assertAlmostEqual(sum(d["weights"]),1.,places=12)
            self.assertAlmostEqual(d["identity_reference_logical_F_e"],
                                   .63+.37/9,places=12)

    def test_ideal_global_instrument_has_unconditional_logical_F_e_one(self):
        for label in ("AS","BS"):
            _,ideal,_,_,_=target_three_outcome(label)
            compiler=DirectLogicalCompiler(label,1,effort_lambda=0.)
            final=ideal@compiler.initial
            value,_,_=compiler._terminal_scores_and_adjoint(final)
            independent=single_fault_logical_recovery_fidelity(
                ideal,label,.37)
            self.assertAlmostEqual(value,1.,places=10)
            self.assertAlmostEqual(
                independent["unconditional_seven_fault_mixture_logical_F_e"],
                1.,places=10)

    def test_zero_angles_identity_is_exact_fair_no_gate_baseline(self):
        for label in ("AS","BS"):
            comp=DirectLogicalCompiler(label,1)
            z=np.zeros(comp.n_angles)
            val=comp.raw_fidelity(z)
            self.assertAlmostEqual(val,.63+.37/9,places=10)
            effort,deriv=comp.smooth_effort_and_gradient(z)
            self.assertAlmostEqual(effort,0.,places=12)
            self.assertTrue(np.allclose(deriv,0.,atol=1e-12))
            self.assertAlmostEqual(comp.objective_and_gradient(z)[0],
                                   1-val,places=10)

    def test_random_angles_match_independent_all_seven_Kraus_channel(self):
        for label in ("AS","BS"):
            comp=DirectLogicalCompiler(label,1)
            angles=np.random.default_rng(2719).normal(
                0,.2,comp.n_angles)
            report=comp.full_report(angles)
            self.assertAlmostEqual(
                report["actual_unconditional_logical_entanglement_F_e"],
                report["independent_all_seven_Kraus_recovery"][
                    "unconditional_seven_fault_mixture_logical_F_e"],
                places=10)
            self.assertTrue(report["no_postselection"])
            self.assertEqual(len(report[
                "independent_all_seven_Kraus_recovery"]["individual_faults"]),7)

    def test_analytic_logical_fidelity_gradient_is_verified(self):
        for label in ("AS","BS"):
            for lam in (0.,.0005):
                report=direct_gradient_check(
                    label,cycles=1,effort_lambda=lam)
                self.assertLess(report["max_error"],3e-7)
                self.assertEqual(len(report[
                    "sampled_coordinate_derivatives"]),5)

    def test_gradient_correct_at_longer_depth(self):
        comp=DirectLogicalCompiler("AS",2,effort_lambda=.0005)
        angles=np.random.default_rng(777).normal(
            0,.15,comp.n_angles)
        val,deriv=comp.objective_and_gradient(angles)
        for k in (0,7,12,13,25):
            v=np.zeros(comp.n_angles)
            v[k]=1e-6
            fd=(comp.objective_and_gradient(angles+v)[0]-
                comp.objective_and_gradient(angles-v)[0])/(2e-6)
            self.assertLess(abs(fd-deriv[k]),3e-7)

    def test_effort_penalty_is_nonnegative_and_weights_bus(self):
        comp=DirectLogicalCompiler("AS",1,effort_lambda=.001)
        a=np.zeros(comp.n_angles)
        dressed=comp.names.index("dressed_swap_link_A0_B0")
        onsite=next(i for i,n in enumerate(comp.names)
                    if n.startswith("onsite_swap_"))
        self.assertEqual(effort_weight_for_name(
            comp.names[dressed]),5.)
        a[dressed]=.3
        high,_=comp.smooth_effort_and_gradient(a)
        a[dressed]=0.;a[onsite]=.3
        low,_=comp.smooth_effort_and_gradient(a)
        self.assertAlmostEqual(high/low,5.,places=10)

    def test_compilation_guarantees_no_regression_below_identity(self):
        for label in ("AS","BS"):
            result=compile_direct(label,1,maxiter=10,seed=617)
            fid=result["logical_fidelity"]
            self.assertGreaterEqual(
                fid["actual_unconditional_logical_entanglement_F_e"]+1e-10,
                fid["do_nothing_baseline_logical_F_e"])
            self.assertEqual(len(result["angles_rad"]),13)
            self.assertTrue(np.all(np.isfinite(result["angles_rad"])))
            self.assertEqual(len(result["attempts"]),3)
            self.assertTrue(
                result["identity_policy_retained_even_if_optimizer_fails"])

    def test_zeroed_policy_has_no_fake_readout_clock_cost(self):
        result=compile_direct("AS",1,maxiter=1,effort_lambda=1000.)
        self.assertTrue(result["readout_skipped_for_literal_identity_policy"])
        self.assertEqual(result["native_primitive_cost"][
            "two_register_native_or_abstract_primitive_count"],0)
        self.assertAlmostEqual(result["native_primitive_cost"][
            "hypothetical_serial_duration_ns_including_readout"],0.,places=12)

    def test_invalid_probabilities_and_penalty_are_rejected(self):
        for v in (-.01,1.01,float("nan")):
            with self.assertRaises(ValueError):
                logical_fault_states("AS",v)
        for v in (-1.,float("nan")):
            with self.assertRaises(ValueError):
                DirectLogicalCompiler("AS",1,effort_lambda=v)
        with self.assertRaises(ValueError):
            DirectLogicalCompiler("AS",1,effort_epsilon=0.)
        with self.assertRaises(ValueError):
            compile_direct("AS",1,maxiter=0)

    def test_seed_is_repeatable(self):
        a=compile_direct("AS",1,maxiter=7,seed=42)
        b=compile_direct("AS",1,maxiter=7,seed=42)
        self.assertEqual(a["selected_initialization"],
                         b["selected_initialization"])
        self.assertLess(np.linalg.norm(
            np.array(a["angles_rad"])-np.array(b["angles_rad"])),1e-10)

    def test_warmstart_depth_transition_and_small_benchmark(self):
        b=benchmark(cycles=(1,2),maxiter=6,labels=("AS",))
        self.assertEqual(len(b["benchmark"]),2)
        self.assertEqual([r["pulse_count"] for r in b["benchmark"]],[13,26])
        self.assertEqual(len(b["benchmark"][1]["attempts"]),4)
        self.assertTrue(b["equality_of_two_block_and_full_seven_error_channel"])
        self.assertTrue(any("NOT the same" in t for t in b["claim_limits"]))
        self.assertTrue(all(x["logical_fidelity"][
            "actual_unconditional_logical_entanglement_F_e"]>=
            x["logical_fidelity"]["do_nothing_baseline_logical_F_e"]-1e-10
            for x in b["benchmark"]))


if __name__=="__main__":
    unittest.main()
