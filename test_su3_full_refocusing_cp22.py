"""Lightweight regression tests for the five-plaquette CP22 refocusing compiler.

Run from repository root:
  OPENBLAS_NUM_THREADS=1 python -m unittest -v test_su3_full_refocusing_cp22
"""
import unittest
import numpy as np

from su3_cubic_casimir_echo import construct
from su3_full_refocusing_cp22 import (
    DIM, MASK, N, PATH, THETA, apply_gate_schedule, compiled_site,
    compile_boundary, local_gate_regression, propagate, toy_risk, threshold_error,
)
from su3_qutrit_gate_compiler import index_groups, strang_segments


class FullFivePlaquetteTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.states, cls.K, cls.d, _, _ = construct()

    def test_local_gate_embedding_and_nine_block_source(self):
        rows = local_gate_regression(self.states)
        self.assertEqual(len(rows), N)
        self.assertEqual([r['cp22'] for r in rows], [4,40,40,40,4])

    def test_cp22_phase_is_correct_on_full_register(self):
        rng=np.random.default_rng(11)
        y=rng.normal(size=(DIM,2))+1j*rng.normal(size=(DIM,2))
        for c,controller_axis in (('L',1),('R',3)):
            actual=apply_gate_schedule(y,2,[('CP22',(c,.32))])
            expected=y.copy().reshape((3,)*N+(2,))
            select=[slice(None)]*(N+1)
            select[2]=select[controller_axis]=2
            expected[tuple(select)]*=np.exp(.32j)
            self.assertLess(np.max(abs(actual-expected.reshape(DIM,2))),1e-13)

    def test_boundary_requires_only_four_entanglers(self):
        for site in (0,4):
            for phase in (0,THETA,2*THETA):
                sequence=compile_boundary(site,.012,phase)
                self.assertEqual(sum(t=='CP22' for t,_ in sequence),4)

    def test_full_exact_matches_existing_direct_strang(self):
        y,counts=propagate(self.K,self.d,1,approximate=False)
        direct,_=strang_segments(self.K,self.states,self.d,index_groups(self.states),rounds=1)
        self.assertLess(np.linalg.norm(y-direct,'fro'),3e-8)
        self.assertEqual(counts['cp22'],10080)
        self.assertEqual(counts['magnetic_site_blocks'],360)
        self.assertEqual(counts['cp22_by_plaquette'],
                         {'0':320,'1':3200,'2':3200,'3':3200,'4':160})

    def test_approximate_budget_and_different_channel(self):
        exact,ce=propagate(self.K,self.d,1,approximate=False)
        approximate,ca=propagate(self.K,self.d,1,approximate=True)
        self.assertEqual(ca['cp22'],8160)
        self.assertLess(np.linalg.norm(approximate.conj().T@approximate-np.eye(2)),2e-9)
        self.assertGreater(np.linalg.norm(exact-approximate),1e-9)

    def test_replacement_model_is_cptp_and_monotone(self):
        counts={'cp22':10080,'local_qutrit':60000}
        f=.999
        base=.632749
        p=threshold_error(f,base,counts['cp22'])
        self.assertGreater(p,0)
        self.assertLess(p,1e-4)
        perfect=toy_risk(f,base,counts,0,0)
        noisy=toy_risk(f,base,counts,1e-4,1e-5)
        self.assertAlmostEqual(perfect['logical_replacement_F_e'],f)
        self.assertLess(noisy['logical_replacement_F_e'],perfect['logical_replacement_F_e'])
        self.assertGreaterEqual(noisy['logical_replacement_F_e'],.25)


if __name__=='__main__': unittest.main()
