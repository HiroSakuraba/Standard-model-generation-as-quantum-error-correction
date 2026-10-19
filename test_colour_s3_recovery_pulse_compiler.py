"""Numerical 42D compiler and all-three-outcome recovery regression tests."""
import unittest
import numpy as np

from colour_s3_recovery_pulse_compiler import (
    DIM,POINTER,TOTAL,NINPUT,
    pointer_generators,controls_for_group,target_three_outcome,
    pulse_pattern,SpectralCompiler,compile_depth,
    native_operation_expansion,finite_gradient_error,
    single_fault_logical_recovery_fidelity,study,
)


class LocalFortyTwoDimensionalCompilerTests(unittest.TestCase):
    def test_all_thirteen_controls_are_exactly_hermitian(self):
        self.assertEqual(len(pointer_generators()),8)
        for label in ("AS","BS"):
            _,_,_,named=controls_for_group(label)
            self.assertEqual(len(named),13)
            for name,matrix in named:
                self.assertEqual(matrix.shape,(42,42))
                self.assertLess(np.linalg.norm(matrix-matrix.conj().T),1e-10)

    def test_target_has_complete_code_electric_and_flagged_paths(self):
        for label in ("AS","BS"):
            named,U,source,target,branches=target_three_outcome(label)
            self.assertEqual(U.shape,(42,42))
            self.assertEqual(source.shape,(42,14))
            self.assertEqual(target.shape,(42,14))
            self.assertEqual({k:v.shape[1] for k,v in branches.items()},
                             {"code":4,"electric":4,"flag":6})
            self.assertLess(np.linalg.norm(U.conj().T@U-np.eye(42)),1e-10)
            self.assertLess(np.linalg.norm(target.conj().T@target-
                                           np.eye(14)),1e-10)

    def test_variable_angle_gate_exponentials_are_unitary(self):
        for label in ("AS","BS"):
            compiler=SpectralCompiler(label,1)
            angles=np.random.default_rng(11).normal(
                0.,.2,len(compiler.order))
            for U in compiler.matrices(angles):
                self.assertLess(np.linalg.norm(
                    U.conj().T@U-np.eye(42)),1e-10)
            final=compiler.propagate(angles)
            self.assertLess(np.linalg.norm(
                final.conj().T@final-np.eye(14)),1e-10)

    def test_analytic_gradient_agrees_with_finite_difference(self):
        for label in ("AS","BS"):
            d=finite_gradient_error(label=label,cycles=1)
            self.assertLess(d["max_gradient_error"],1e-7)
            self.assertEqual(len(d["gradient_checks"]),5)

    def test_gradient_of_full_unitary_mode_is_also_correct(self):
        compiler=SpectralCompiler("AS",1,full_unitary=True)
        angles=np.random.default_rng(24).normal(
            0.,.08,len(compiler.order))
        loss,deriv=compiler.objective_and_gradient(angles)
        self.assertTrue(0<=loss<=1+1e-10)
        for k in (0,4,12):
            delta=np.zeros(len(angles));delta[k]=1e-6
            fd=(compiler.objective(angles+delta)-
                compiler.objective(angles-delta))/(2e-6)
            self.assertLess(abs(fd-deriv[k]),1e-7)

    def test_zero_angles_do_not_claim_to_correct_or_flag(self):
        for label in ("AS","BS"):
            compiler=SpectralCompiler(label,1)
            diags=compiler.diagnostics(np.zeros(13))
            self.assertLess(diags[
                "three_outcome_input_isometry_coherent_fidelity"],.5)
            # The clean input is already in the code, so even the
            # identity schedule correctly preserves its state/pointer.
            self.assertAlmostEqual(diags["per_branch"]["code"][
                "intended_pointer_probability"],1.,places=10)
            self.assertLess(diags["per_branch"]["electric"][
                "intended_pointer_probability"],1e-10)
            self.assertLess(diags["per_branch"]["flag"][
                "intended_pointer_probability"],1e-10)
            self.assertLess(diags[
                "entire_42D_unitary_fidelity_to_selected_extension"],1.)

    def test_bounded_optimizer_produces_actual_angles_and_no_regression(self):
        for label in ("AS","BS"):
            report=compile_depth(label,cycles=1,maxiter=12,seed=17)
            self.assertEqual(len(report["angles_rad"]),13)
            self.assertEqual(report["pulse_count"],13)
            self.assertTrue(np.all(np.abs(report["angles_rad"])<=np.pi+1e-10))
            self.assertGreaterEqual(report["fidelity_diagnostics"][
                "three_outcome_input_isometry_coherent_fidelity"]+1e-9,
                report["zero_angle_identity_baseline_fidelity"])
            self.assertEqual(len(report["attempts"]),2)
            self.assertFalse(report["fidelity_diagnostics"][
                "no_postselection_used"] is False)

    def test_repeatable_seed_and_cycle_order(self):
        a=compile_depth("AS",1,maxiter=5,seed=614)
        b=compile_depth("AS",1,maxiter=5,seed=614)
        self.assertEqual(a["pulse_names"],b["pulse_names"])
        self.assertLess(np.linalg.norm(np.array(a["angles_rad"])-
                                       np.array(b["angles_rad"])),1e-10)
        self.assertEqual(pulse_pattern(13,2),tuple(range(13))*2)
        with self.assertRaises(ValueError):
            pulse_pattern(13,0)

    def test_five_charged_bus_operations_per_active_dressed_pulse(self):
        compiler=SpectralCompiler("AS",1)
        angles=np.zeros(13)
        dressed=next(i for i,name in enumerate(
            compiler.named) if name[0].startswith("dressed_swap_"))
        angles[dressed]=.6
        info={
            "pulse_count":13,
            "pulse_names":[x[0] for x in compiler.named],
            "angles_rad":angles.tolist(),
        }
        steps=native_operation_expansion(info,transport_ns=100.)
        self.assertEqual(steps["two_register_native_or_abstract_primitive_count"],5)
        self.assertEqual(steps["charged_bus_swap_and_fractional_steps"],3)
        self.assertEqual(steps["charged_link_transport_steps"],2)
        self.assertAlmostEqual(
            steps["expanded_primitive_schedule"][2]["rotation_rad"],
            -1.2,places=12)
        self.assertGreater(steps[
            "hypothetical_serial_duration_ns_including_readout"],600)

    def test_native_estimates_reject_nonphysical_speed_or_negative_time(self):
        info={"pulse_count":13,"pulse_names":["electric_link_standard_projector"],
              "angles_rad":[.5]}
        with self.assertRaises(ValueError):
            native_operation_expansion(info,electric_mhz=0)
        with self.assertRaises(ValueError):
            native_operation_expansion(info,transport_ns=-1)
        with self.assertRaises(ValueError):
            native_operation_expansion(info,readout_ns=float("nan"))

    def test_ideal_42D_target_perfectly_recovers_all_seven_faults(self):
        for label in ("AS","BS"):
            _,target,_,_,_=target_three_outcome(label)
            r=single_fault_logical_recovery_fidelity(target,label,.37)
            self.assertAlmostEqual(r[
                "unconditional_seven_fault_mixture_logical_F_e"],1.,places=10)
            self.assertEqual(len(r["individual_faults"]),7)
            for fault in r["individual_faults"]:
                self.assertAlmostEqual(fault[
                    "unconditional_recovered_logical_entanglement_fidelity"],
                    1.,places=10)
                self.assertAlmostEqual(fault[
                    "logical_code_capture_probability"],1.,places=10)
                self.assertLess(fault[
                    "pointer_2_flagged_probability"],1e-10)
                self.assertLess(fault[
                    "unflagged_outside_logical_code_probability"],1e-10)

    def test_identity_pulse_is_not_falsely_treated_as_recovered_electric_fault(self):
        for label in ("AS","BS"):
            r=single_fault_logical_recovery_fidelity(
                np.eye(42),label,.37)
            self.assertAlmostEqual(r[
                "unconditional_seven_fault_mixture_logical_F_e"],
                .63+.37/9,places=10)
            bad=next(x for x in r["individual_faults"] if
                     x["error"]=="electric_phase")
            self.assertAlmostEqual(bad[
                "unconditional_recovered_logical_entanglement_fidelity"],
                1/9,places=10)
            with self.assertRaises(ValueError):
                single_fault_logical_recovery_fidelity(
                    np.eye(42),label,1.1)

    def test_warmstart_depth_never_silent_reset(self):
        first=compile_depth("AS",1,maxiter=5)
        second=compile_depth("AS",2,maxiter=5,
                             warmstart=first["angles_rad"])
        self.assertEqual(second["pulse_count"],26)
        self.assertEqual(len(second["attempts"]),3)
        self.assertGreaterEqual(
            second["fidelity_diagnostics"][
                "three_outcome_input_isometry_coherent_fidelity"]+1e-9,
            second["zero_angle_identity_baseline_fidelity"])

    def test_small_benchmark_does_not_conflate_reachability_with_success(self):
        report=study(cycles=(1,),maxiter=5,labels=("AS",))
        self.assertEqual(len(report["benchmark"]),1)
        sample=report["benchmark"][0]
        self.assertEqual(sample["objective_mode"],
            "coherent_14D_input_three_outcome_isometry")
        self.assertTrue(sample["explicit_compiled_schedule_exists"])
        self.assertTrue(sample["no_device_performance_claim"])
        self.assertFalse(report["benchmark"][0][
            "fidelity_diagnostics"]["no_postselection_used"] is False)
        self.assertTrue(any("NOT the same" in item
                            for item in report["limitations"]))

if __name__=="__main__":
    unittest.main()
