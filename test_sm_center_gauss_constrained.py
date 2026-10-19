"""Tests for physical Gauss-space finite quotient centre-echo negative controls."""
import unittest
import numpy as np
from sm_center_gauss_constrained import ORDERS,basis,operators,experiment,valid

class GaussCentreEchoTests(unittest.TestCase):
    def test_dimensions_and_vertex_constraints(self):
        for q in ORDERS:
            states=basis(q)
            self.assertEqual(len(states),16*(6//q))
            self.assertTrue(all(valid(*s,6//q) for s in states))
            self.assertEqual(len(set(states)),len(states))

    def test_matter_transport_preserved_and_hermiticity(self):
        for q in ORDERS:
            E,M,G,P,meta=operators(q)
            self.assertGreater(np.linalg.norm(M),1)
            self.assertLess(np.linalg.norm(M-M.conj().T),1e-11)
            self.assertLess(np.linalg.norm(M@P-P@M),1e-11)

    def test_global_form_quotient_no_nontrivial_centre_pulse(self):
        E,M,G,P,_=operators(6)
        self.assertLess(np.linalg.norm(P-np.eye(len(P))),1e-12)
        self.assertLess(np.linalg.norm(G),1e-12)

    def test_smaller_quotients_have_charge_sensitive_pulses(self):
        for q in (1,2,3):
            E,M,G,P,_=operators(q)
            self.assertGreater(np.linalg.norm(P-np.eye(len(P))),1)
            self.assertGreater(np.linalg.norm(G),1)

    def test_full_six_frame_twirling_and_finite_convergence(self):
        for q in ORDERS:
            result=experiment(q,cycles=(4,16))
            self.assertLess(result['averaged_magnetic_cancellation_residual_fro'],1e-10)
            a,b=result['finite_echo_errors_per_sqrt_dimension']
            if q==6:
                self.assertLess(max(a,b),1e-9)
            else:
                self.assertLess(b,a)


if __name__=='__main__':unittest.main()
