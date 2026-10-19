"""Regressions for strict Gauss-at-every-step two-register no-go."""
import unittest

import numpy as np

from colour_gauge_link_s3 import (
    DIM, G, transform, physical_isometry, dressed_exchange,
)
from colour_s3_edge_logical_syndrome import code_isometry
from colour_s3_strict_gauss_pair_nogo import (
    IRREPS, link_central_projectors, apply_link, pointer_link_irrep_gate,
    irrep_pointer_kraus, compressed_pair_matrix_units, generated_pair_algebra,
    main_result,
)


class StrictGaussPairTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.W=physical_isometry()
        cls.K,cls.F,_=code_isometry(cls.W)
        cls.P=link_central_projectors()

    def test_central_link_projectors_are_exact_and_gauge_covariant(self):
        self.assertEqual(list(self.P),list(IRREPS))
        self.assertEqual([round(np.trace(x).real) for x in self.P.values()],
                         [1,1,4])
        for name,P in self.P.items():
            self.assertLess(np.linalg.norm(P@P-P),1e-11)
            for other,Q in self.P.items():
                if name!=other:
                    self.assertLess(np.linalg.norm(P@Q),1e-11)
            z=apply_link(self.W,P)
            for v in ('A','B'):
                for h in range(G):
                    self.assertLess(np.linalg.norm(transform(z,h,v)-z),1e-10)

    def test_peter_weyl_physical_sector_ranks(self):
        ranks=[np.trace(self.W.conj().T@apply_link(self.W,P)).real
               for P in self.P.values()]
        self.assertLess(np.linalg.norm(np.array(ranks)-[4,1,9]),1e-10)
        self.assertAlmostEqual(sum(ranks),14,places=10)

    def test_dressed_swap_breaks_sector_conservation(self):
        for P in self.P.values():
            comm=apply_link(self.F,P)
            blocked=self.K.conj().T@comm
            self.assertGreater(np.linalg.norm(blocked),.4)
            # Full code K and F are orthogonal, but NOT per sector:
            self.assertLess(np.linalg.norm(self.K.conj().T@self.F),1e-10)

    def test_safe_two_register_pointer_is_complete_and_covariant(self):
        U=pointer_link_irrep_gate()
        self.assertEqual(U.shape,(18,18))
        self.assertLess(np.linalg.norm(U.conj().T@U-np.eye(18)),1e-11)
        rng=np.random.default_rng(11)
        x=rng.normal(size=(DIM,2))+1j*rng.normal(size=(DIM,2))
        got=irrep_pointer_kraus(x)
        target=np.stack([apply_link(x,self.P[k]) for k in IRREPS])
        self.assertLess(np.linalg.norm(got-target),1e-10)
        self.assertLess(abs(np.vdot(got,got)-np.vdot(x,x)),1e-10)
        for v in ('A','B'):
            for h in range(G):
                transformed=irrep_pointer_kraus(transform(x,h,v))
                reference=np.stack([transform(z,h,v) for z in got])
                self.assertLess(np.linalg.norm(transformed-reference),1e-10)

    def test_one_link_irrep_pointer_does_not_distinguish_error(self):
        def probabilities(x):
            return np.array([np.linalg.norm(z)**2/4
                             for z in irrep_pointer_kraus(x)])
        good,bad=probabilities(self.K),probabilities(self.F)
        self.assertLess(np.linalg.norm(good-bad),1e-10)
        self.assertLess(np.linalg.norm(good-np.array([5/24,1/8,2/3])),1e-10)

    def test_pair_matrix_algebra_closure_is_block_diagonal(self):
        pair=compressed_pair_matrix_units(self.W,0,1)
        self.assertEqual(pair.shape,(18**2,14**2))
        result=generated_pair_algebra(self.W)
        self.assertEqual(result['pair_linear_span_dimension'],29)
        self.assertEqual(result['associative_closure_ranks'],[29,98,98])
        self.assertEqual(result['associative_algebra_dimension'],4**2+1+9**2)

    def test_complete_report_keeps_scope(self):
        result=main_result(include_algebra=False)
        self.assertTrue(result['no_go']['conditions'])
        for value in result['healthy_vs_dressed_syndrome_block_overlap_fro'].values():
            self.assertGreater(value,.4)
        self.assertFalse(result['one_pointer_qutrit_irrep_readout']
                         ['resolves_known_dressed_swap'])


if __name__=="__main__":
    unittest.main()
