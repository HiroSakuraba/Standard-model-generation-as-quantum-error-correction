"""Tests of *exact* pair-gate factorization and native-cost claim boundaries."""
import unittest
import numpy as np

from colour_gauge_link_s3 import (
    DIM, G, physical_isometry, dressed_exchange, transform, projector
)
from colour_s3_edge_logical_syndrome import (
    pointer_circuit, onsite_swap, code_isometry
)
from colour_s3_two_register_syndrome import (
    exchange_eigenbasis_q, link_b_matrix, pointer_phase_matrix,
    pair_gate, controlled_dressed_factors, controlled_onsite_factors,
    measured_syndrome_factors, correction_factors, prepared_pointer,
    apply_syndrome, apply_correction_branch, gauss_probability,
    circuit_gauge_diagnostics, local_weyl, sample_two_register_weyl,
    channel_fidelity, run
)

class TwoRegisterSyndromeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.W=physical_isometry()
        cls.K,cls.F,_=code_isometry(cls.W)

    def test_pair_gates_are_exact_unitaries(self):
        self.assertLess(np.linalg.norm(exchange_eigenbasis_q().conj().T
                     @exchange_eigenbasis_q()-np.eye(9)),1e-12)
        self.assertLess(np.linalg.norm(link_b_matrix().conj().T
                     @link_b_matrix()-np.eye(18)),1e-12)
        for gate in measured_syndrome_factors()+correction_factors():
            self.assertEqual(len(gate),3)
            self.assertEqual(gate[2].shape[0],gate[2].shape[1])
            self.assertLess(np.linalg.norm(gate[2].conj().T@gate[2]
                      -np.eye(gate[2].shape[0])),2e-12)
        self.assertEqual(len(measured_syndrome_factors()),16)
        self.assertEqual(len(correction_factors()),3)

    def test_controlled_dressed_swap_for_every_pointer_branch(self):
        rng=np.random.default_rng(3)
        z=rng.normal(size=(3,DIM,2))+1j*rng.normal(size=(3,DIM,2))
        for active in ((2,),(1,2),(0,)):
            y=z.copy()
            for _,axes,U in controlled_dressed_factors(active):
                y=pair_gate(y,axes,U)
            expected=np.stack([dressed_exchange(z[p]) if p in active
                               else z[p] for p in range(3)])
            self.assertLess(np.linalg.norm(y-expected),1e-10)

    def test_controlled_onsite_swap_for_every_pointer_branch(self):
        rng=np.random.default_rng(4)
        z=rng.normal(size=(3,DIM,2))+1j*rng.normal(size=(3,DIM,2))
        for active in ((2,),(1,2)):
            y=z.copy()
            for _,axes,U in controlled_onsite_factors(active):
                y=pair_gate(y,axes,U)
            expected=np.stack([onsite_swap(z[p]) if p in active
                               else z[p] for p in range(3)])
            self.assertLess(np.linalg.norm(y-expected),1e-10)

    def test_unconditioned_three_pair_feedback_matches_source(self):
        rng=np.random.default_rng(81)
        y=rng.normal(size=(DIM,2))+1j*rng.normal(size=(DIM,2))
        self.assertLess(np.linalg.norm(apply_correction_branch(y)
                            -dressed_exchange(y)),1e-10)

    def test_entire_four_controlled_gate_syndrome_circuit(self):
        rng=np.random.default_rng(41)
        y=rng.normal(size=(DIM,4))+1j*rng.normal(size=(DIM,4))
        actual=apply_syndrome(y)
        expected=pointer_circuit(y)
        self.assertLess(np.linalg.norm(actual-expected),1e-10)
        self.assertLess(abs(np.vdot(actual,actual)-np.vdot(y,y)),1e-10)

    def test_intermediate_gauss_violations_are_measured_and_resolved(self):
        result=circuit_gauge_diagnostics(self.K)
        self.assertEqual(len(result['trace']),18)
        for v in ('A','B'):
            self.assertAlmostEqual(result['final_Gauss_accept'][v],1.,places=10)
            self.assertGreaterEqual(result['minimum_intermediate_Gauss_accept'][v],0.)
            self.assertLessEqual(result['minimum_intermediate_Gauss_accept'][v],1.)
        # After each full high-level gauge-invariant controlled gate,
        # the physical input must be back in Gauss sector.
        for i in (5,8,13,16):
            row=result['trace'][i]
            for v in ('A','B'):
                self.assertAlmostEqual(row['Gauss_accept'][v],1.,places=10)
        # At least one non-gauge-covariant constituent breaks Gauss
        self.assertTrue(any(result['minimum_intermediate_Gauss_accept'][v]<.99
                            for v in ('A','B')))

    def test_exact_cptp_two_register_weyl_unitaries(self):
        rng=np.random.default_rng(82)
        for pair in ((0,4),(1,4),(2,4),(2,3)):
            d1=(3,6,3,3,3,3)[pair[0]]
            d2=(3,6,3,3,3,3)[pair[1]]
            gate=sample_two_register_weyl(pair,rng)
            self.assertEqual(gate.shape,(d1*d2,d1*d2))
            self.assertLess(np.linalg.norm(gate.conj().T@gate
                             -np.eye(d1*d2)),2e-12)

    def test_noise_zero_exact_recovery_and_nontrivial_noisy_trajectories(self):
        ideal=channel_fidelity(self.K,fault_probability=0.,shots=3)
        self.assertAlmostEqual(ideal['mean_recovered_entanglement_fidelity'],1.,places=10)
        noisy=channel_fidelity(self.K,fault_probability=.2,shots=10,seed=82)
        self.assertTrue(0<=noisy['mean_recovered_entanglement_fidelity']<=1)
        self.assertLess(noisy['mean_recovered_entanglement_fidelity'],1.)

    def test_full_report_respects_native_gate_caveats(self):
        data=run()
        self.assertEqual(data['two_register_pair_ops_syndrome'],16)
        self.assertEqual(data['two_register_pair_ops_feedback_if_error'],3)
        self.assertLess(data['full_pointer_unitary_ref_error_fro'],1e-10)
        self.assertEqual(data['native_two_qubit_or_two_qutrit_gate_count'],'NOT derived')

if __name__=="__main__":
    unittest.main()
