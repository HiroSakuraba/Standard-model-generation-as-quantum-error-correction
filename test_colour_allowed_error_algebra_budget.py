"""Regression: exact 65-syndrome colour Gram and allowed exchange algebra."""
import unittest
import numpy as np
from colour_allowed_error_algebra_budget import (
    N,COLORS,DIM,LOGICAL,allowed_graph_components,
    exact_single_colour_gram,orthogonal_modes,
    physical_gauge_commutator,exchange_obstructions,study,
)
class AllowedColourErrorTests(unittest.TestCase):
    def test_full_KL_Gram_rank_and_spectrum(self):
        C=exact_single_colour_gram()
        self.assertEqual(C.shape,(73,73))
        eigen=np.linalg.eigvalsh(C)
        self.assertEqual(int(np.sum(eigen>1e-10)),65)
        self.assertEqual(int(np.sum(abs(eigen)<1e-10)),8)
        self.assertLess(np.min(eigen),1e-10)
        self.assertEqual(int(np.sum(abs(eigen-3/16)<1e-10)),64)

    def test_eight_relative_site_modes_per_adj_component(self):
        C=exact_single_colour_gram()
        h,B=orthogonal_modes(C)
        self.assertEqual(h.shape,(8,9))
        self.assertEqual(B.shape,(73,64))
        self.assertLess(np.linalg.norm(B.T@C@B-np.eye(64)),1e-11)

    def test_same_vertex_global_gauge_vs_independent_gauges(self):
        p=physical_gauge_commutator()
        self.assertLess(p['max_common_vertex_colour_commutator_fro'],1e-12)
        self.assertGreater(p['min_independent_vertex_gauss_commutator_fro'],.01)

    def test_exchange_connected_graph_is_full_symmetric_group(self):
        self.assertEqual(allowed_graph_components([(i,i+1) for i in range(1,N)]),
                         [tuple(range(1,N+1))])
        self.assertEqual(allowed_graph_components([(1,2),(2,3),(7,8)]),
                         [(1,2,3),(4,),(5,),(6,),(7,8),(9,)])

    def test_single_exchange_detectable_but_two_locations_not_correctable(self):
        x=exchange_obstructions()
        self.assertLess(x['maximum_single_exchange_KLK_fro'],1e-10)
        self.assertGreater(x['two_unknown_swap_error_KL_non_scalar_fro'],1.)
        self.assertGreater(x['second_order_native_exchange_return_non_scalar_fro'],.1)
        self.assertFalse(x['passively_protected_from_all_allowed_exchanges'])
        self.assertEqual(x['complex_associative_exchange_algebra_dimension'],42**2)

    def test_ancilla_information_count_and_nonclaim(self):
        report=study()
        g=report['single_site_error_gram']
        self.assertEqual(g['rank_of_KL_gram'],65)
        self.assertEqual(g['orthogonal_error_images_dimension'],1365)
        self.assertEqual(g['minimal_orthogonal_syndrome_pointer_qutrits'],4)
        self.assertEqual(g['minimal_orthogonal_syndrome_pointer_qubits'],7)
        self.assertLess(g['orthonormalized_error_basis_residual_fro'],1e-10)

if __name__=='__main__':
    unittest.main()
