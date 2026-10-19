"""Exact qutrit phase-estimation measurement of edge-local gauge-invariant swap."""
import unittest
import numpy as np
from colour_gauge_link_s3 import (
    DIM,G,physical_isometry,dressed_exchange,transform,projector,
)
from colour_s3_edge_logical_syndrome import (
    W3,DLOG,onsite_swap,edge_cycle,projectors,pointer_circuit,
    code_isometry,recover_channel,run,
)

class EdgeLogicalSyndromeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.W=physical_isometry()
        cls.K,cls.F,cls.U=code_isometry(cls.W)

    def test_order_three_and_spectrum(self):
        self.assertLess(np.linalg.norm(edge_cycle(edge_cycle(edge_cycle(self.W)))-self.W),1e-10)
        eig=np.linalg.eigvals(self.U)
        self.assertEqual([int(sum(abs(eig-W3**m)<1e-8)) for m in range(3)],[6,4,4])
        self.assertLess(np.linalg.norm(self.K.conj().T@self.F),1e-10)

    def test_separate_edge_and_vertex_operations_gauge_covariant(self):
        for op in (onsite_swap,dressed_exchange,edge_cycle):
            for v in ('A','B'):
                for h in range(G):
                    self.assertLess(np.linalg.norm(
                        transform(op(self.W),h,v)-op(transform(self.W,h,v))),1e-10)

    def test_phase_pointer_is_explicit_complete_cptp_instrument(self):
        rng=np.random.default_rng(989)
        y=rng.normal(size=(DIM,4))+1j*rng.normal(size=(DIM,4))
        branches=pointer_circuit(y)
        self.assertEqual(branches.shape,(3,DIM,4))
        self.assertLess(np.linalg.norm(branches-projectors(y)),1e-10)
        self.assertLess(abs(np.vdot(branches,branches)-np.vdot(y,y)),1e-10)
        # Coherent superpositions are tagged linearly; ancilla readout is
        # not measuring the logical basis inside a given eigenphase.
        a=.6; b=.8j
        z=a*self.K+b*self.F
        out=pointer_circuit(z)
        self.assertLess(np.linalg.norm(out[1]-a*self.K),1e-10)
        self.assertLess(np.linalg.norm(out[2]-b*self.F),1e-10)
        self.assertLess(np.linalg.norm(out[0]),1e-10)

    def test_syndrome_labels_preserve_arbitrary_logical_superposition(self):
        rng=np.random.default_rng(201)
        c=rng.normal(size=4)+1j*rng.normal(size=4)
        c=c/np.linalg.norm(c)
        for code,expected in ((self.K,1),(self.F,2)):
            out=pointer_circuit(code@c)
            self.assertLess(np.linalg.norm(out[expected]-code@c),1e-10)
            for m in range(3):
                if m!=expected:self.assertLess(np.linalg.norm(out[m]),1e-10)
        self.assertTrue(0<abs(c[0])<1)

    def test_logical_phase_pointer_is_not_just_gauss_measurement(self):
        for y in (self.K,self.F):
            self.assertEqual(y.shape,(DIM,4))
            for v in ('A','B'):
                self.assertLess(np.linalg.norm(projector(y,v)-y),1e-10)
        self.assertLess(np.linalg.norm(pointer_circuit(self.K)[1]-self.K),1e-10)
        self.assertLess(np.linalg.norm(pointer_circuit(self.F)[2]-self.F),1e-10)

    def test_nonpostselected_known_swap_channel_recovered(self):
        for p in (0.,.17,.37,1.):
            data=recover_channel(self.K,p,assignment=.01)
            self.assertAlmostEqual(data['correctly_reported_F_e'],1.,places=10)
            self.assertAlmostEqual(data['flag_and_replace_F_e'],.99+.01/16,places=10)
            self.assertLess(np.linalg.norm(data['true_pointer_probabilities_0_1_2']-
                         np.array([0,1-p,p])),1e-10)

    def test_cost_and_scope_are_explicit(self):
        data=run()
        self.assertEqual(data['high_level_cost']['pointer_qutrits'],1)
        self.assertEqual(data['high_level_cost']['pointer_controlled_dressed_edge_swaps'],2)
        self.assertEqual(data['high_level_cost']['pointer_controlled_onsite_swaps'],2)
        self.assertEqual(data['logical_code_dimension'],4)
        self.assertEqual(data['eigenphases_1_w_w2_multiplicities'],[6,4,4])
        self.assertIn('missing',data['high_level_cost']['four_register_pointer_link_edge_primitive_native_compilation'])

if __name__=='__main__':unittest.main()
