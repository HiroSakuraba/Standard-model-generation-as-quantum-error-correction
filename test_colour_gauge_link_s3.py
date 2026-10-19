"""Regression: finite S3-subgroup SU3 gauge link and local syndrome pointer."""
import unittest
import numpy as np
from colour_gauge_link_s3 import (
    PERMS,G,DIM,REP,compose,inverse,ID,physical_isometry,
    dressed_exchange,transform,projector,dft_ancilla_kraus,
    logical_swap_code,study,local_gauss_probabilities
)

class S3GaugeLinkTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.W=physical_isometry()

    def test_s3_faithful_SU3_representation(self):
        self.assertEqual(len(PERMS),6)
        for a in range(G):
            self.assertAlmostEqual(np.linalg.det(REP[a]),1.)
            self.assertEqual(compose(a,inverse(a)),ID)
            for b in range(G):
                self.assertLess(np.linalg.norm(REP[compose(a,b)]-REP[a]@REP[b]),1e-12)

    def test_vertex_transformations_commute_and_invariant_basis(self):
        W=self.W
        self.assertEqual(W.shape,(486,14))
        self.assertLess(np.linalg.norm(W.conj().T@W-np.eye(14)),1e-11)
        rng=np.random.default_rng(5)
        v=rng.standard_normal(DIM)+1j*rng.standard_normal(DIM)
        for a in range(G):
            for b in range(G):
                self.assertLess(np.linalg.norm(transform(transform(v,a,'A'),b,'B')
                    -transform(transform(v,b,'B'),a,'A')),1e-11)
        for vertex in ('A','B'):
            for h in range(G):
                self.assertLess(np.linalg.norm(transform(W,h,vertex)-W),1e-10)

    def test_dressed_swap_gauge_covariant_bare_fails(self):
        W=self.W
        dressed=dressed_exchange(W)
        bare=dressed_exchange(W,False)
        self.assertLess(np.linalg.norm(dressed_exchange(dressed)-W),1e-11)
        self.assertLess(np.linalg.norm(projector(dressed,'A')-dressed),1e-10)
        self.assertLess(np.linalg.norm(projector(dressed,'B')-dressed),1e-10)
        self.assertGreater(np.linalg.norm(projector(bare,'A')-bare),1)
        self.assertGreater(np.linalg.norm(projector(bare,'B')-bare),1)

    def test_vertex_six_outcome_pointer_CPTP(self):
        rng=np.random.default_rng(44)
        y=rng.standard_normal((DIM,2))+1j*rng.standard_normal((DIM,2))
        for vertex in ('A','B'):
            branch=dft_ancilla_kraus(y,vertex)
            self.assertEqual(branch.shape,(6,DIM,2))
            self.assertLess(abs(np.vdot(branch,branch)-np.vdot(y,y)),1e-10)
            self.assertLess(np.linalg.norm(branch[0]-projector(y,vertex)),1e-10)
            pure=dft_ancilla_kraus(self.W,vertex)
            self.assertLess(np.linalg.norm(pure[0]-self.W),1e-10)
            self.assertLess(np.linalg.norm(pure[1:]),1e-10)

    def test_gauss_pointer_cannot_detect_gauge_invariant_swap(self):
        W=self.W
        swapped=dressed_exchange(W)
        p=local_gauss_probabilities(swapped)
        self.assertAlmostEqual(p['A'],1.,places=10)
        self.assertAlmostEqual(p['B'],1.,places=10)
        bare=local_gauss_probabilities(dressed_exchange(W,False))
        self.assertAlmostEqual(bare['A'],.25,places=10)
        self.assertAlmostEqual(bare['B'],.25,places=10)

    def test_logical_syndrome_known_dressed_swap_recovery(self):
        K,sk,stats=logical_swap_code(self.W)
        self.assertEqual(stats['logical_code_dimension'],4)
        self.assertEqual(stats['swap_positive_multiplicity'],10)
        self.assertEqual(stats['swap_negative_multiplicity'],4)
        self.assertLess(np.linalg.norm(K.conj().T@sk),1e-10)
        result=study()
        self.assertAlmostEqual(result['known_swap_channel_ideal_nonlocal_decoder_F_e'],1.)
        self.assertEqual(result['controlled_group_branches_per_vertex_in_this_circuit'],6)
        self.assertAlmostEqual(result['illustrative_reported_gauss_accept_probabilities']['bare_swap_A'],.249)

if __name__=='__main__':
    unittest.main()
