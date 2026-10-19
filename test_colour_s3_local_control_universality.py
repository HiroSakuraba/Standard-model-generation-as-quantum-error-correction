"""Exact modular local-gauge Hamiltonian controllability and charged-bus tests."""
import unittest
import numpy as np
import sympy as sp

from colour_gauge_link_s3 import DIM,G,transform,physical_isometry
from colour_s3_exact_electric_commutant import (
    exact_generators,physical_orbit_isometry,
)
from colour_s3_order_six_active_electric import closure_matrices
from colour_s3_local_control_universality import (
    FIELD_PRIMES,same_colour_local_projector,same_colour_full_486,
    admissible_controls,rational_mod_p,finite_field_rank,
    exact_modular_lie_closure,algebra_commutant_prime,
    gauge_covariance_local_projector,controlled_link_electric_full_486,
    charged_mediator_fractional_dressed_swap,
    pointer_universality_argument,study,
)


class GaugeLocalControllabilityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.classes,cls.lookup,cls.ops=exact_generators()
        cls.W,cls.sizes=physical_orbit_isometry(cls.classes)

    def options(self):
        A,B,S,E=self.ops
        for name,g in (("AS",(A,S)),("BS",(B,S))):
            H,inventory=admissible_controls(self.classes,self.ops,name)
            yield name,H,inventory,closure_matrices(g)

    def test_opposite_vertex_equality_projectors_are_exact_idempotents(self):
        for pair in ((0,1),(2,3)):
            Q=same_colour_local_projector(self.classes,pair)
            self.assertEqual(Q*Q,Q)
            self.assertEqual(Q,Q.T)
            self.assertTrue(0<sp.trace(Q)<14)
        with self.assertRaises(ValueError):
            same_colour_local_projector(self.classes,(0,2))

    def test_equality_projectors_act_only_on_true_two_register_vertex(self):
        rng=np.random.default_rng(643)
        x=rng.normal(size=(DIM,3))+1j*rng.normal(size=(DIM,3))
        for pair in ((0,1),(2,3)):
            projected=same_colour_full_486(x,pair)
            self.assertLess(np.linalg.norm(
                same_colour_full_486(projected,pair)-projected),2e-11)
            self.assertLess(gauge_covariance_local_projector(pair),2e-11)
            compressed=self.W.conj().T@same_colour_full_486(self.W,pair)
            Q=same_colour_local_projector(self.classes,pair)
            normalized=np.array(Q.tolist(),float)
            self.assertLess(np.linalg.norm(compressed-normalized),2e-11)

    def test_each_generator_is_physically_hermitian(self):
        from colour_s3_order_six_active_electric import gauge_fixed_normalized_matrix
        for label,H,inventory,group in self.options():
            self.assertEqual([x[1] for x in inventory],[2,3,1,2])
            for h in H:
                m=gauge_fixed_normalized_matrix(h,self.sizes)
                self.assertLess(np.linalg.norm(m-m.conj().T),2e-11)

    def test_full_u14_control_lie_rank_exact_over_two_fields(self):
        for label,H,inventory,group in self.options():
            for p in FIELD_PRIMES:
                r=exact_modular_lie_closure(H,p)
                self.assertEqual(r["Lie_rank_mod_prime_lower_bound"],196)
                self.assertEqual(r["real_skew_rank"],91)
                self.assertEqual(r["imaginary_symmetric_rank"],105)
                self.assertTrue(r["full_u14_rational_rank_certified"])

    def test_omit_opposite_vertex_projector_is_not_full_control(self):
        for label,H,inventory,group in self.options():
            for p in FIELD_PRIMES:
                result=exact_modular_lie_closure(H[:3],p)
                self.assertEqual(result["Lie_rank_mod_prime_lower_bound"],21)
                self.assertFalse(result["full_u14_rational_rank_certified"])

    def test_same_vertex_control_is_inadequate_negative_control(self):
        classes=self.classes
        for label,H,inventory,group in self.options():
            wrong_pair=(0,1) if label=="AS" else (2,3)
            wrong=same_colour_local_projector(classes,wrong_pair)
            r=exact_modular_lie_closure(H[:3]+(wrong,),101)
            self.assertEqual(r["Lie_rank_mod_prime_lower_bound"],106)

    def test_scalar_commutant_independent_modular_certificate(self):
        E=self.ops[-1]
        for label,H,inventory,group in self.options():
            for p in FIELD_PRIMES:
                r=algebra_commutant_prime(group,E,H[-1],p)
                self.assertEqual(r["electric_plus_pair_equality_constraint_rank"],51)
                self.assertTrue(r["exact_scalar_only_commutant"])

    def test_charged_mediator_fractional_dressed_swap_continuum(self):
        for angle in (0.,np.pi/5,np.pi/2,np.pi):
            r=charged_mediator_fractional_dressed_swap(
                self.W[:,:3],angle)
            self.assertEqual(r["two_register_primitive_count"],5)
            self.assertLess(r[
                "original_three_register_S_fractional_pulse_reconstruction_Fro"],2e-10)

    def test_charged_mediator_stays_gauss_invariant_at_every_step(self):
        result=charged_mediator_fractional_dressed_swap(
            self.W[:,:2],np.pi/2,history=True)
        self.assertEqual(len(result["intermediate_Gauss_acceptance"]),5)
        self.assertTrue(all(min(a.values())>1-1e-10
                            for a in result["intermediate_Gauss_acceptance"]))

    def test_pointer_link_electric_gate_is_exactly_two_register_gauss_covar(self):
        rng=np.random.default_rng(4)
        z=rng.normal(size=(3,DIM,2))+1j*rng.normal(size=(3,DIM,2))
        for theta in (0.,.4,np.pi):
            out=controlled_link_electric_full_486(z,theta)
            self.assertAlmostEqual(
                np.linalg.norm(out),np.linalg.norm(z),places=10)
            for vertex in ("A","B"):
                for h in range(G):
                    lhs=controlled_link_electric_full_486(
                        np.stack([transform(x,h,vertex) for x in z]),theta)
                    rhs=np.stack([transform(x,h,vertex) for x in out])
                    self.assertLess(np.linalg.norm(lhs-rhs),2e-10)

    def test_pointer_universality_is_conditional_not_pulse_list(self):
        E=self.ops[-1]
        report=pointer_universality_argument(E)
        self.assertEqual(report["coupled_link_projector_rank_on_physical_system"],9)
        self.assertEqual(report["physical_plus_pointer_dimension"],42)
        self.assertEqual(report["traceless_dimension"],1763)
        self.assertTrue(report["full_unitary_reachable_in_ideal_continuous_control_model"])
        self.assertFalse(report["explicit_solved_pulse_sequence_for_PR64_global_recovery"])

    def test_end_to_end_evidence_and_honest_claim_boundary(self):
        result=study()
        self.assertEqual(len(result["subgroup_local_control_results"]),2)
        self.assertEqual(len(result["continuous_dressed_swap_five_pair_gate_exact_identity"]),4)
        self.assertFalse(result["claimed_finite_native_compiled_42D_PR64_recovery"])
        self.assertLess(result["neutral_pointer_link_electric_pair_covariance_Fro"],2e-10)
        self.assertTrue(any("NOT arbitrary continuous SU3" in x for x in
                            result["limitations"]))

if __name__=="__main__":
    unittest.main()
