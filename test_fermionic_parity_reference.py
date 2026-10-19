"""Fast regression for fermionic parity-frame code versus unique RM(1,4) state."""
import unittest
import numpy as np
from fermionic_parity_reference import (
    construct_even_code,rank_certificate,audit,logical_dressed_X,
    total_parity,reference_Z
)


class ParityReferenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.checks,cls.z0,cls.z1,cls.K=construct_even_code()

    def test_two_dim_even_parity_code(self):
        self.assertEqual(self.K.shape,(65536,2))
        self.assertLess(np.linalg.norm(self.K.conj().T@self.K-np.eye(2)),1e-12)
        self.assertLess(np.linalg.norm(total_parity(self.K[:,0])-self.K[:,0]),1e-12)
        self.assertLess(np.linalg.norm(total_parity(self.K[:,1])-self.K[:,1]),1e-12)

    def test_dressed_even_logical_paulis(self):
        for j in (0,1):
            self.assertLess(np.linalg.norm(logical_dressed_X(self.K[:,j])-self.K[:,1-j]),1e-12)
            self.assertLess(np.linalg.norm(reference_Z(self.K[:,j])-(-1)**j*self.K[:,j]),1e-12)

    def test_extra_stabilizer_eliminates_encoded_qubit(self):
        plus=(self.K[:,0]+self.K[:,1])/np.sqrt(2)
        minus=(self.K[:,0]-self.K[:,1])/np.sqrt(2)
        self.assertLess(np.linalg.norm(logical_dressed_X(plus)-plus),1e-12)
        self.assertLess(np.linalg.norm(logical_dressed_X(minus)+minus),1e-12)
        ranks=rank_certificate(self.checks)
        self.assertEqual(ranks['relational_logical_qubits'],1)
        self.assertEqual(ranks['full_unpunctured_logical_qubits'],0)

    def test_parity_reference_retains_relational_coherence(self):
        r=audit()
        self.assertGreater(r['relational_logical_X_expectation'],.9)
        self.assertAlmostEqual(r['isolated_15_mode_parity_twirled_X_expectation'],0.,places=11)
        self.assertEqual(r['full_unpunctured_16_mode_rank'],16)


if __name__=='__main__':unittest.main()
