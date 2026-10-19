"""Fast independent regression for the 42D nine-quark sign-twist code."""
import unittest
import numpy as np
from colour_sign_twist_theorem import (
    tableaux, adjacent_generators, signed_tableau_transpose, code_isometry,
    swaps, exchange_survival, exchange_channel_recovery, run, DIM
)


class SignTwistTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.T=tableaux()
        cls.S=adjacent_generators(cls.T)
        cls.J=signed_tableau_transpose(cls.T)
        cls.K=code_isometry(cls.J)
        cls.swap=swaps(cls.S)

    def test_42_dim_independent_specht(self):
        self.assertEqual(len(self.T),42)
        self.assertEqual(self.K.shape,(42,21))
        self.assertEqual(len(self.swap),36)

    def test_braid_orthogonal_and_sign_twist(self):
        self.assertLess(np.linalg.norm(self.J.T@self.J-np.eye(DIM)),1e-12)
        self.assertLess(np.linalg.norm(self.J+self.J.T),1e-12)
        self.assertLess(np.linalg.norm(self.J@self.J+np.eye(DIM)),1e-12)
        for i,s in enumerate(self.S):
            self.assertLess(np.linalg.norm(s@s-np.eye(DIM)),1e-11)
            self.assertLess(np.linalg.norm(self.J@s+s@self.J),1e-11)
            if i<7:
                b=self.S[i+1]
                self.assertLess(np.linalg.norm(s@b@s-b@s@b),1e-11)

    def test_every_exchange_maps_code_out(self):
        for M in self.swap.values():
            self.assertLess(np.linalg.norm(self.K.conj().T@M@self.K),1e-10)

    def test_three_cycle_acts_nontrivially(self):
        V=self.swap[(1,2)]@self.swap[(2,3)]
        sub=self.K.conj().T@V@self.K
        self.assertGreater(np.linalg.norm(sub-np.trace(sub)/21*np.eye(21)),1)

    def test_exchange_dynamics_not_harmless(self):
        S=self.swap[(1,2)]
        for t in (0.,np.pi/3,np.pi/2,np.pi):
            U=np.cos(t/2)*np.eye(DIM)-1j*np.sin(t/2)*S
            val=np.linalg.norm(self.K.conj().T@U@self.K)**2/21
            self.assertAlmostEqual(val,exchange_survival(t),places=11)

    def test_explicit_channel_recovery(self):
        for p in (0.,.03,.37,1.):
            result=exchange_channel_recovery(self.K,self.J,self.swap[(1,2)],p)
            self.assertAlmostEqual(result['entanglement_fidelity_after_explicit_recovery'],1.,places=11)

    def test_complete_report(self):
        report=run()
        self.assertLess(report['all_exchange_KLK_max_fro_residual'],1e-10)
        self.assertEqual(report['claims']['maximality_21'],
                         'only among zero-exchange-compression subspaces, NOT globally proved')


if __name__=='__main__':
    unittest.main()
