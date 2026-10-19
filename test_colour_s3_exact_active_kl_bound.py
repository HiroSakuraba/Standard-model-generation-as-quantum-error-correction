"""Exact error-set classification and perfect 4D active KL obstruction tests."""
import unittest

import sympy as sp

from colour_s3_exact_electric_commutant import (
    exact_generators, crosscheck_existing_physical_model,
)
from colour_s3_exact_active_kl_bound import (
    K, EXPECTED_SUBGROUPS, CHARACTERS,
    all_S4_character_projectors, exact_electric_split_on_unique_swap_code,
    full_fault_operator_span_ranks,
    generated_group_order, generator_actions, common_eigenspace_dimension,
    classify_electric_subgroups, full_group_active_KL_no_go, exact_study,
    cycle_type,
)


class ExactActiveKLNoGoTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.classes,cls.lookup,cls.ops=exact_generators()
        cls.P,cls.multiplicities=all_S4_character_projectors()
        cls.electric=exact_electric_split_on_unique_swap_code(
            cls.P["trivial"],cls.ops[-1])
        cls.census=classify_electric_subgroups(cls.ops)

    def test_S4_character_table_has_all_five_irreps(self):
        self.assertEqual(len(CHARACTERS),5)
        self.assertEqual(cycle_type((0,1,2,3)),(1,1,1,1))
        self.assertEqual(cycle_type((1,0,2,3)),(2,1,1))
        self.assertEqual(cycle_type((1,2,0,3)),(3,1))
        self.assertEqual(cycle_type((1,2,3,0)),(4,))
        self.assertEqual(cycle_type((1,0,3,2)),(2,2))

    def test_exact_irrep_multiplicities_are_four_two_two(self):
        m=self.multiplicities
        self.assertEqual(m["trivial"]["multiplicity"],4)
        self.assertEqual(m["partition_22_2"]["multiplicity"],2)
        self.assertEqual(m["standard_3"]["multiplicity"],2)
        self.assertEqual(m["sign"]["multiplicity"],0)
        self.assertEqual(m["sign_twist_3"]["multiplicity"],0)
        self.assertEqual(sum(v["isotypic_physical_dimension"]
                             for v in m.values()),14)

    def test_all_character_projectors_are_mutually_orthogonal_exact(self):
        for label,P in self.P.items():
            self.assertEqual(P*P,P)
            for other,Q in self.P.items():
                if label!=other:self.assertEqual(P*Q,sp.zeros(14))
        self.assertEqual(sum(self.P.values(),sp.zeros(14)),sp.eye(14))

    def test_link_electric_is_non_scalar_on_unique_all_swap_code(self):
        e=self.electric
        self.assertEqual(e["compressed_electric_eigenvalues"],
                         ["0","1/3","2/3","1"])
        self.assertFalse(e["electric_is_scalar_on_unique_four_state_swap_sector"])
        self.assertFalse(e["complete_24_swaps_plus_electric_perfect_KL_four_state_code"])

    def test_number_of_errors_is_not_linear_span_dimension(self):
        exact=full_fault_operator_span_ranks(self.ops[-1])
        self.assertEqual(exact["listed_swap_Kraus_operators"],24)
        self.assertEqual(exact["swap_operator_linear_span_dimension"],14)
        self.assertEqual(exact["listed_operators_including_electric"],25)
        self.assertEqual(exact["swap_plus_electric_operator_linear_span_dimension"],15)

    def test_full_24_swaps_plus_R_exact_active_KL_no_go(self):
        result=full_group_active_KL_no_go(self.multiplicities,self.electric)
        self.assertTrue(result["exact_proof"])
        self.assertEqual(result["number_distinct_allowed_swap_kraus_operators"],24)
        self.assertEqual(result["irrep_blocks_permitting_rank4_map_from_logical_space"],
                         ["trivial"])
        self.assertFalse(result["perfect_four_state_KL_for_full_25_operator_set_exists"])

    def test_eight_exchange_subgroups_have_exact_orders(self):
        rows={r["generator_subset"]:r for r in self.census}
        self.assertEqual(len(rows),8)
        for label,(order,dimension) in EXPECTED_SUBGROUPS.items():
            self.assertEqual(rows[label]["whole_generated_exchange_group_order"],order)
            self.assertEqual(rows[label]["max_passive_common_eigenspace_dimension"],dimension)

    def test_proper_subgroups_with_no_four_dimensional_passive_DFS(self):
        rows={r["generator_subset"]:r for r in self.census}
        self.assertEqual(sorted(k for k,v in rows.items()
                         if not v["passive_four_state_DFS_exists"]),["ABS","AS","BS"])
        self.assertEqual(rows["AS"]["max_passive_common_eigenspace_dimension"],2)
        self.assertEqual(rows["BS"]["max_passive_common_eigenspace_dimension"],2)
        self.assertEqual(rows["ABS"]["max_passive_common_eigenspace_dimension"],1)

    def test_commuting_two_swap_subgroup_still_has_passive_ququart(self):
        rows={r["generator_subset"]:r for r in self.census}
        self.assertEqual(rows["AB"]["whole_generated_exchange_group_order"],4)
        self.assertEqual(rows["AB"]["max_passive_common_eigenspace_dimension"],4)
        self.assertTrue(rows["AB"]["passive_four_state_DFS_exists"])

    def test_exact_rational_model_matches_original_link_Gauss_embedding(self):
        residual=crosscheck_existing_physical_model(self.classes,self.ops)
        self.assertLess(residual["max_residual"],2e-11)

    def test_full_report_distinguishes_full_Kraus_list_and_generators(self):
        report=exact_study()
        self.assertEqual(report["candidate_subgroups_without_four_dimensional_common_passive_DFS"],
                         ["AS","BS","ABS"])
        self.assertEqual(len(report["all_eight_generator_subgroups_with_link_electric_phase"]),8)
        self.assertTrue(report["full_exchange_plus_electric_active_KL_no_go"]["exact_proof"])
        self.assertEqual(report["restricted_active_KL_existence_for_AS_BS_ABS"].split()[0],
                         "open")
        self.assertTrue(any("not merely three generator" in x
                            for x in report["scientific_claim_boundaries"]))


if __name__=="__main__":
    unittest.main()
