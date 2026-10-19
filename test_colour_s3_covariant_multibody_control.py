"""Tests of continuous strictly Gauss-invariant multibody Hamiltonian control."""
import unittest

import numpy as np

from colour_gauge_link_s3 import (
    DIM, G, physical_isometry, dressed_exchange, transform, projector
)
from colour_s3_edge_logical_syndrome import (
    onsite_swap, pointer_circuit, code_isometry
)
from colour_s3_two_register_syndrome import (
    prepared_pointer, circuit_gauge_diagnostics
)
from colour_s3_covariant_multibody_control import (
    PULSE_TIME, FEEDBACK_TIME, LOGICAL,
    dressed_hamiltonian, onsite_hamiltonian, physical_pulse,
    controlled_pulse, control_schedule, phase_readout,
    gauge_acceptance, continuous_gauge_diagnostics,
    gaussian_phase_noise, analog_gadget, pair_gadget,
    analog_feedback, pair_feedback, circuit_fidelity, study
)


class CovariantMultibodyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.W=physical_isometry()
        cls.K,cls.F,_=code_isometry(cls.W)

    def test_three_register_generator_and_exact_pi_gate(self):
        rng=np.random.default_rng(7)
        x=rng.normal(size=(DIM,3))+1j*rng.normal(size=(DIM,3))
        for name,op,H in (
            ("edge",dressed_exchange,dressed_hamiltonian),
            ("onsite",onsite_swap,onsite_hamiltonian)):
            self.assertLess(np.linalg.norm(H(H(x))),0.1+np.linalg.norm(x))
            self.assertLess(np.linalg.norm(physical_pulse(x,name,np.pi)-op(x)),2e-12)
            self.assertLess(np.linalg.norm(physical_pulse(x,name,.7).conj().T
                  @physical_pulse(x,name,.7)-x.conj().T@x),1e-10)
            self.assertLess(np.linalg.norm(physical_pulse(
                physical_pulse(x,name,.19),name,.27)
                -physical_pulse(x,name,.46)),2e-12)

    def test_hamiltonian_commutes_with_gauss_for_arbitrary_full_states(self):
        rng=np.random.default_rng(12)
        x=rng.normal(size=(DIM,2))+1j*rng.normal(size=(DIM,2))
        for name in ("edge","onsite"):
            for angle in (.0,.18,1.1,np.pi):
                for v in ("A","B"):
                    for h in range(G):
                        lhs=transform(physical_pulse(x,name,angle),h,v)
                        rhs=physical_pulse(transform(x,h,v),name,angle)
                        self.assertLess(np.linalg.norm(lhs-rhs),2e-11)

    def test_controlled_gauge_path_is_physical_at_quarter_times(self):
        result=continuous_gauge_diagnostics(self.K)
        self.assertLess(result["maximum_Gauss_covariance_residual"],1e-10)
        self.assertEqual(len(result["trace"]),21)
        for v in ("A","B"):
            self.assertGreater(result["minimum_Gauss_acceptance"][v],1-1e-10)
            prior=circuit_gauge_diagnostics(self.K)["minimum_intermediate_Gauss_accept"][v]
            self.assertLess(prior,.99)

    def test_four_pulse_pointer_matches_full_486dim_reference(self):
        rng=np.random.default_rng(19)
        x=rng.normal(size=(DIM,2))+1j*rng.normal(size=(DIM,2))
        z=prepared_pointer(x)
        for name,active in control_schedule():
            z=controlled_pulse(z,name,active)
        self.assertLess(np.linalg.norm(phase_readout(z)-pointer_circuit(x)),2e-11)
        self.assertLess(abs(np.vdot(z,z)-np.vdot(x,x)),2e-10)

    def test_pointer_control_is_four_register_not_three_register(self):
        # Unconditioned edge pulse spans A0-link-B0 (three).
        # Neutral-pointer conditioning also touches pointer (four).
        z=prepared_pointer(self.K)
        for kind,active in control_schedule():
            partial=controlled_pulse(z,kind,active,np.pi*.37)
            for v in ("A","B"):
                self.assertGreater(gauge_acceptance(partial)[v],1-1e-10)
            z=controlled_pulse(z,kind,active)
        # Error state remains gauge invariant, correct outcome is label 2.
        e=prepared_pointer(self.F)
        for kind,active in control_schedule():
            e=controlled_pulse(e,kind,active)
        self.assertLess(np.linalg.norm(phase_readout(e)[2]-self.F),1e-10)

    def test_phase_diffusion_random_unitary_preserves_norm(self):
        rng=np.random.default_rng(82)
        z=prepared_pointer(self.K)
        for dt in (0,1,16):
            y=gaussian_phase_noise(z,dt,.01,rng)
            self.assertLess(abs(np.vdot(y,y)-np.vdot(z,z)),1e-11)
        self.assertLess(np.linalg.norm(gaussian_phase_noise(z,4.,0.,rng)-z),1e-12)

    def test_noiseless_known_error_recovery_both_architectures(self):
        for architecture in ("covariant_multibody","pair_decomposed"):
            row=circuit_fidelity(self.K,architecture,0.,shots=2,seed=82)
            self.assertAlmostEqual(row["mean_entanglement_fidelity"],1.,places=10)
            self.assertEqual(row["notional_syndrome_time"],PULSE_TIME)
            self.assertEqual(row["notional_feedback_or_idle_time"],FEEDBACK_TIME)

    def test_nonzero_matched_phase_noise_is_bounded_and_reproducible(self):
        a=circuit_fidelity(self.K,"covariant_multibody",.02,shots=4,seed=90)
        b=circuit_fidelity(self.K,"pair_decomposed",.02,shots=4,seed=90)
        c=circuit_fidelity(self.K,"covariant_multibody",.02,shots=4,seed=90)
        self.assertTrue(0<=a["mean_entanglement_fidelity"]<=1)
        self.assertTrue(0<=b["mean_entanglement_fidelity"]<=1)
        self.assertEqual(a["individual_trajectory_fidelities"],
                         c["individual_trajectory_fidelities"])
        self.assertLess(a["mean_entanglement_fidelity"],1)
        self.assertLess(b["mean_entanglement_fidelity"],1)

    def test_study_reports_validated_resources_and_scopes(self):
        result=study(shots=2)
        self.assertEqual(result["logical_dimension"],4)
        self.assertEqual(result["physical_dimension"],486)
        self.assertEqual(result["resource_comparison"]["covariant_multibody"]
                         ["conditional_edge_support_registers"],4)
        self.assertEqual(result["resource_comparison"]["pair_decomposed"]
                         ["general_U9_or_U18_two_register_factors"],16)
        self.assertEqual(len(result["matched_noise_comparison"]),6)


if __name__=="__main__":
    unittest.main()
