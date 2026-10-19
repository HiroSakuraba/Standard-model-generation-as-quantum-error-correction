"""Regression tests for full five-ququart SU3 leakage/noise trajectories and pulse optimizer."""
import unittest
import numpy as np

from su3_cubic_casimir_echo import construct
from su3_channel_recovery import syndrome_decoder, source_data
from su3_full_refocusing_cp22 import apply_gate_schedule, compiled_site
from su3_five_ququart_trajectories import (
    SIZE,COMPUTATIONAL,IN_CODE,initial_entangled_map,
    apply_one,apply_cp22_mixing,cp_mixing,loss_kraus,
    sample_noise,readout_fidelity,apply_compiled_site,simulate
)
from su3_cp22_pulse_optimization import optimize_angle


class FullFiveQuquartTrajectoryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.states,cls.K,cls.d,_,_=construct()
        cls.V,_,cls.G=syndrome_decoder(*source_data())

    def test_computational_embedding(self):
        m=initial_entangled_map(self.K)
        self.assertEqual(m.shape,(1024,2))
        self.assertTrue(np.array_equal(m[COMPUTATIONAL],self.K))
        self.assertEqual(len(COMPUTATIONAL),243)
        self.assertLess(abs(np.vdot(m,m).real-2),1e-12)

    def test_cp22_coherent_leakage_is_unitary(self):
        for theta,eps in ((.1,0.),(.4,.01),(-.1,.07)):
            u=cp_mixing(theta,eps)
            self.assertLess(np.linalg.norm(u.conj().T@u-np.eye(2)),1e-12)
        m=np.zeros((SIZE,2),complex)
        m[2*4**4+2*4**3,0]=1
        m[3*4**4+1*4**3,1]=1
        result=apply_cp22_mixing(m,(0,1),.25,.1)
        self.assertLess(abs(np.vdot(m,m)-np.vdot(result,result)),1e-12)
        self.assertGreater(abs(result[3*4**4+1*4**3,0]),0)

    def test_four_level_kraus_completeness(self):
        for dt in (0.,25.,100000.):
            ops=loss_kraus(dt,30000.)
            total=sum((k.conj().T@k for k in ops),np.zeros((4,4),complex))
            self.assertLess(np.linalg.norm(total-np.eye(4)),5e-11)

    def test_ideal_recovery_and_noisy_six_sector_mislabel(self):
        m=initial_entangled_map(self.K)
        f,leak=readout_fidelity(m,self.V,0.)
        self.assertAlmostEqual(f,1.,places=11)
        self.assertAlmostEqual(leak,0.,places=11)
        f2,_=readout_fidelity(m,self.V,.5)
        self.assertAlmostEqual(f2,.4375,places=10)

    def test_loss_phase_unraveling_norm_preserved(self):
        m=initial_entangled_map(self.K)
        y,jumps=sample_noise(m,1000.,30000.,45000.,np.random.default_rng(44))
        self.assertTrue(np.isfinite(y).all())
        self.assertAlmostEqual(np.vdot(y,y).real,2.,places=7)
        self.assertGreaterEqual(jumps,0)

    def test_exact_embedded_site_matches_prior_three_level_gate(self):
        m=initial_entangled_map(self.K)
        for site in (0,2,4):
            seq=compiled_site(site,.013,0.,False)
            old=apply_gate_schedule(self.K,site,seq)
            got,counts=apply_compiled_site(m,site,.013,0.,False,0.)
            self.assertLess(np.linalg.norm(got[COMPUTATIONAL]-old),5e-10)
            self.assertLess(np.linalg.norm(got[~IN_CODE]),2e-10)
            self.assertEqual(counts['cp22'],4 if site in (0,4) else 40)

    def test_one_cycle_full_five_ququart_trajectory_executes(self):
        outcome=simulate(True,shots=1,seed=17,cycles=1,
                        leakage=.003,t1_ns=30000.,tphi_ns=45000.)
        self.assertEqual(outcome['serial_resources']['cp22'],1020)
        self.assertEqual(outcome['hilbert_dimension'],1024)
        self.assertTrue(0<=outcome['sample_mean_recovered_entanglement_fidelity']<=1)

    def test_low_budget_optimizer_reports_bounded_point(self):
        row=optimize_angle(np.pi,iterations=1,steps=8)
        self.assertEqual(row['target_angle_radians'],float(np.pi))
        self.assertTrue(32.<=row['finer_step_evaluation']['duration_ns']<=200.)
        self.assertTrue(210.<=row['finer_step_evaluation']['delta_pulse_mhz']<=420.)
        self.assertTrue(0<=row['finer_step_evaluation']['coherent_target_F_e']<=1)


if __name__=='__main__':unittest.main()
