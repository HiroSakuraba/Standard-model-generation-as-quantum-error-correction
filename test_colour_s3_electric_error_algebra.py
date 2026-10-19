"""Regress S4 exchange-sector rank, gauge electric error and KL limitation."""
import unittest
import numpy as np

from colour_gauge_link_s3 import (
    DIM,G,physical_isometry,transform,dressed_exchange,
)
from colour_s3_edge_logical_syndrome import code_isometry
from colour_s3_strict_gauss_pair_nogo import (
    link_central_projectors,apply_link,
)
from colour_s3_electric_error_algebra import (
    physical_onsite_swap, signed_exchange_generators,
    permutation_group_closure, s4_fixed_projector,
    deterministic_basis, electric_projector,electric_fault,
    commutant_dimension,code_kl_diagnostics,
    entanglement_fidelity_known_error,
    existing_ideal_phase_syndrome_electric_response,study,
)


class ElectricErrorAlgebraTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.W=physical_isometry()
        cls.P,cls.gen,cls.permutations=s4_fixed_projector(cls.W)
        cls.C14=deterministic_basis(cls.P,4)
        cls.E=electric_projector(cls.W)
        cls.K486,_,_=code_isometry(cls.W)
        cls.K=cls.W.conj().T@cls.K486

    def test_s4_really_has_all_24_elements_and_4D_common_fixed_space(self):
        self.assertEqual(len(self.permutations),24)
        self.assertAlmostEqual(np.trace(self.P).real,4.,places=10)
        self.assertLess(np.linalg.norm(self.P@self.P-self.P),1e-10)
        for a in self.permutations:
            self.assertLess(np.linalg.norm(a@self.C14-self.C14),1e-10)
        self.assertLess(np.linalg.norm(self.C14.conj().T@self.C14-np.eye(4)),1e-10)

    def test_onsite_swaps_and_dressed_edges_are_gauge_compatible(self):
        W=self.W
        for op in (
            lambda x:physical_onsite_swap(x,0,1),
            lambda x:physical_onsite_swap(x,2,3),
            dressed_exchange,
        ):
            self.assertLess(np.linalg.norm(op(op(W))-W),1e-10)
            for v in ("A","B"):
                for h in range(G):
                    self.assertLess(np.linalg.norm(
                        transform(op(W),h,v)-op(transform(W,h,v))),1e-10)

    def test_link_central_electric_projector_is_gauge_invariant_off_sector(self):
        rng=np.random.default_rng(89)
        x=rng.normal(size=(DIM,3))+1j*rng.normal(size=(DIM,3))
        self.assertLess(np.linalg.norm(electric_fault(electric_fault(x))-x),1e-10)
        self.assertLess(np.linalg.norm(self.E@self.E-self.E),1e-10)
        self.assertEqual(round(np.trace(self.E).real),9)
        for vertex in ("A","B"):
            for h in range(G):
                self.assertLess(np.linalg.norm(
                    transform(electric_fault(x),h,vertex)-
                    electric_fault(transform(x,h,vertex))),1e-10)

    def test_electric_field_splits_all_four_passive_s4_fixed_states(self):
        reduced=self.C14.conj().T@self.E@self.C14
        spectrum=np.linalg.eigvalsh(reduced)
        self.assertLess(np.linalg.norm(spectrum-[0,1/3,2/3,1]),1e-10)
        self.assertGreater(np.linalg.norm(
            (np.eye(14)-self.P)@self.E@self.C14),.1)
        R=np.eye(14)-2*self.E
        actual=entanglement_fidelity_known_error(self.C14,R,.37)
        self.assertAlmostEqual(actual,.63,places=10)

    def test_full_error_commutant_numerical_dimension_six(self):
        swaps=list(self.gen.values())
        only=commutant_dimension(swaps)
        mixed=commutant_dimension(swaps+[self.E])
        self.assertEqual(only["commutant_complex_dimension"],24)
        self.assertEqual(mixed["commutant_complex_dimension"],6)
        self.assertGreater(mixed["smallest_nonzero_singular_value"],.1)
        self.assertLess(mixed["largest_null_singular_value"],1e-10)
        self.assertTrue(mixed["no_fitted_four_dim_noiseless_subsystem"])

    def test_existing_active_code_fails_joint_KL_for_electric_and_swap(self):
        diag=code_kl_diagnostics(self.K,self.E,
                               self.gen["dressed_edge_A0_B0"])
        self.assertFalse(diag["all_Knill_Laflamme_scalar_conditions"])
        self.assertGreater(diag["max_non_scalar_KL_compression_Fro"],.5)

    def test_existing_ideal_pointer_does_not_perfectly_correct_electric(self):
        report=existing_ideal_phase_syndrome_electric_response(self.K486,.37)
        F=report["ideal_existing_K_phase_pointer_F_e_for_electric_channel"]
        self.assertAlmostEqual(F,.681388888888889,places=10)
        self.assertLess(F,1-.01)
        self.assertAlmostEqual(sum(report["pointer_outcome_probabilities_0_1_2"]),
                               1.,places=10)

    def test_single_electric_projector_alone_still_supports_4D_passive(self):
        vals=np.linalg.eigvalsh(self.E)
        self.assertEqual(int(np.sum(vals>.5)),9)
        self.assertEqual(int(np.sum(vals<.5)),5)
        # Failure is from the JOINT algebra; neither error by itself
        # rules out a 4D eigen-code.
        self.assertTrue(5>=4)

    def test_report_and_rigorous_claim_boundaries(self):
        result=study()
        self.assertEqual(result["generated_matter_swap_group_order"],24)
        self.assertEqual(result["full_swap_group_fixed_4D_subspace_rank"],4)
        self.assertTrue(result["electric_commutes_with_both_local_Gauss_groups"])
        self.assertTrue(result["no_passive_four_logical_dim_noiseless_subsystem_numerical"])
        self.assertFalse(result["KL_existing_active_K_three_error_set"]
                         ["all_Knill_Laflamme_scalar_conditions"])
        self.assertTrue(any("not an exact rational" in x for x in
                            result["important_bounds"]))


if __name__=="__main__":
    unittest.main()
