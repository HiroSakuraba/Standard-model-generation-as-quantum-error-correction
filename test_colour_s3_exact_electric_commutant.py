"""Regression checks for rational centralizer certificate and restricted KL code."""
import unittest
import numpy as np
import sympy as sp

from colour_s3_exact_electric_commutant import (
    four_colour_orbits, exact_generators, positional_s4,
    invariant_orbital_matrices, electric_centralizer_constraint,
    modular_rank, full_commutator_ranks, physical_orbit_isometry,
    crosscheck_existing_physical_model, restricted_swap_electric_code,
    exact_study,
)


class ExactElectricCommutantTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.classes,cls.lookup,cls.ops=exact_generators()
        cls.group=positional_s4(cls.classes,cls.lookup)
        cls.basis=invariant_orbital_matrices(cls.group)
        cls.D=electric_centralizer_constraint(cls.basis,cls.ops[-1])

    def test_colour_orbit_dimension_exactly_fourteen(self):
        self.assertEqual(len(self.classes),14)
        self.assertEqual(len(self.lookup),81)
        self.assertEqual([len(c) for c in self.classes].count(3),1)
        self.assertEqual([len(c) for c in self.classes].count(6),13)

    def test_three_matter_swaps_are_exact_integer_involutions(self):
        for x in self.ops[:3]:
            self.assertEqual(x*x,sp.eye(14))
            self.assertTrue(all(a.q==1 for a in x))
        self.assertNotEqual(self.ops[0]*self.ops[2],
                            self.ops[2]*self.ops[0])

    def test_electric_link_rational_projector_rank_nine(self):
        E=self.ops[-1]
        self.assertEqual(E*E,E)
        self.assertEqual(E.rank(),9)
        self.assertEqual(E.trace(),9)
        self.assertTrue(all((3*a).q==1 for a in E))

    def test_S4_orbit_centralizer_has_24_disjoint_indicators(self):
        self.assertEqual(len(self.group),24)
        self.assertEqual(len(self.basis),24)
        cover=sum(self.basis,sp.zeros(14))
        self.assertEqual(cover,sp.ones(14,14))
        for X in self.basis:
            for g in self.ops[:3]:
                self.assertEqual(X*g,g*X)

    def test_exact_rational_commutant_reduced_rank_eighteen(self):
        self.assertEqual(self.D.shape,(196,24))
        self.assertEqual(self.D.rank(),18)
        null=self.D.nullspace()
        self.assertEqual(len(null),6)
        for v in null:self.assertEqual(self.D*v,sp.zeros(196,1))

    def test_modular_full_784x196_rank_independently_190(self):
        ranks=full_commutator_ranks(self.ops)
        self.assertEqual(ranks,{"101":190,"65521":190})
        self.assertEqual(14*14-190,6)

    def test_exact_rational_representation_matches_original_486D_link(self):
        witness=crosscheck_existing_physical_model(self.classes,self.ops)
        self.assertLess(witness["max_residual"],2e-11)
        W,sizes=physical_orbit_isometry(self.classes)
        self.assertEqual(W.shape,(486,14))
        self.assertLess(np.linalg.norm(W.conj().T@W-np.eye(14)),1e-11)
        self.assertEqual(sorted(set(sizes)),[3.,6.])

    def test_restricted_active_four_state_code_satisfies_KL(self):
        r=restricted_swap_electric_code(self.ops,[len(x) for x in self.classes])
        self.assertEqual(r["S_plus_space_dimension"],10)
        self.assertEqual(r["active_syndrome_space_rank"],8)
        self.assertEqual(r["KL_correctable_active_code_dimension"],4)
        self.assertLess(max(r["exact_KL_residuals_numeric_construction"].values()),2e-11)

    def test_restricted_four_state_passive_code_also_exists(self):
        r=restricted_swap_electric_code(self.ops,[len(x) for x in self.classes])
        self.assertEqual(r["passive_joint_S_plus_E_one_space_dimension"],5)
        self.assertTrue(r["passive_4D_code_also_exists"])
        self.assertEqual(r["electric_compressed_on_S_plus_eigenspace"],
                         {"zero":2,"one_half":3,"one":5})

    def test_complete_exact_proof_and_scope(self):
        r=exact_study()
        self.assertEqual(r["joint_commutant_exact_complex_dimension"],6)
        self.assertEqual(r["reduced_commutator_exact_rational_rank"],18)
        self.assertEqual(r["electric_projector_exact_rank"],9)
        self.assertTrue(r["exact_no_four_dim_noiseless_subsystem_under_S4_plus_E"])
        self.assertTrue(any("does not rule out" in b for b in r["claim_limits"]))
        self.assertFalse(any("full SU3" in str(r.get("reason")) for _ in [0]))


if __name__=="__main__":
    unittest.main()
