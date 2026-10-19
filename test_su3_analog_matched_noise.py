"""Source-grounded tests: native SU3 analog, equal-time noise, two qutrit syndrome pointers."""
import unittest
from dataclasses import replace
import numpy as np
from su3_cubic_casimir_echo import construct, entanglement_fidelity
from su3_channel_recovery import source_data, syndrome_decoder
from su3_analog_matched_noise import (
    Scenario,DIM,POINTERS,POINTER_INDEX,syndrome_projectors,
    tag_syndrome_pointers,pointer_single_transition,pointer_confusion,
    syndrome_recovered_fidelity,loss_kraus,noise_segment,
    propagation_schedule,time_ns,run_one,study
)


class AnalogMatchedNoiseTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.states,cls.K,cls.d,cls.W,cls.H=construct()
        cls.V,_,gram=syndrome_decoder(*source_data())
        assert gram.rank()==5

    def test_six_projectors_resolve_identity_and_ranks(self):
        ps=syndrome_projectors(self.V)
        self.assertEqual(len(ps),6)
        self.assertEqual([round(np.trace(p).real) for p in ps],
                         [2,2,2,2,2,DIM-10])
        for i in range(6):
            for j in range(i):
                self.assertLess(np.linalg.norm(ps[i]@ps[j]),2e-10)
        self.assertLess(np.linalg.norm(sum(ps)-np.eye(DIM)),2e-10)

    def test_two_pointer_isometry_and_independent_tagging(self):
        ps=syndrome_projectors(self.V)
        rng=np.random.default_rng(12)
        m=rng.normal(size=(DIM,2))+1j*rng.normal(size=(DIM,2))
        tagged=tag_syndrome_pointers(m,ps)
        self.assertEqual(tagged.shape,(9,DIM,2))
        self.assertLess(abs(np.vdot(tagged,tagged)-np.vdot(m,m)),2e-9)
        for r,p in enumerate(ps):
            self.assertLess(np.linalg.norm(tagged[POINTER_INDEX[r]]-p@m),2e-10)
        for unused in set(range(9))-set(POINTER_INDEX):
            self.assertLess(np.linalg.norm(tagged[unused]),1e-12)
        self.assertEqual(len(set(POINTERS)),6)

    def test_pointer_noise_confusion_trace_preserving(self):
        perfect=pointer_single_transition(0.,30000.,0.)
        self.assertLess(np.linalg.norm(perfect-np.eye(3)),1e-12)
        imperfect=pointer_confusion(Scenario())
        self.assertEqual(imperfect.shape,(9,6))
        self.assertLess(np.max(abs(imperfect.sum(axis=0)-1)),1e-12)
        self.assertGreater(imperfect[POINTER_INDEX[0],0],
                           imperfect[POINTER_INDEX[5],5])

    def test_exact_syndrome_recovery_matches_parent_decoder(self):
        from su3_channel_recovery import recover_fidelity
        C=np.zeros((9,6))
        for i in range(6): C[POINTER_INDEX[i],i]=1
        expected=recover_fidelity([self.K],self.V,self.K)[0]
        actual,weight=syndrome_recovered_fidelity(self.K,self.V,C)
        self.assertAlmostEqual(actual,expected,places=11)
        self.assertAlmostEqual(actual,1.,places=11)
        self.assertAlmostEqual(weight,1.,places=11)
        noisy,_=syndrome_recovered_fidelity(self.K,self.V,pointer_confusion(Scenario()))
        self.assertLess(noisy,actual)

    def test_qutrit_pure_loss_completeness(self):
        for dt in (0,1000.,50000.):
            ks=loss_kraus(dt,30000.)
            self.assertEqual(len(ks),3)
            self.assertLess(np.linalg.norm(sum((k.conj().T@k for k in ks),
                    np.zeros((3,3),complex))-np.eye(3)),2e-11)

    def test_physical_noise_preserves_reference_norm(self):
        y=self.K.copy()
        rng=np.random.default_rng(27)
        unchanged,jumps=noise_segment(y,0.,30000.,45000.,rng)
        self.assertLess(np.linalg.norm(unchanged-y),1e-12)
        self.assertEqual(jumps,0)
        noisy,jumps=noise_segment(y,3000.,30000.,45000.,rng)
        self.assertLess(abs(np.vdot(noisy,noisy).real-2),1e-7)
        self.assertTrue(np.isfinite(noisy).all())
        self.assertGreaterEqual(jumps,0)

    def test_equal_time_schedule_and_nonlinear_magnetic_control(self):
        scen=Scenario()
        schedule=propagation_schedule(self.states,self.H,scen)
        self.assertEqual(len(schedule),72)
        self.assertAlmostEqual(sum(seg[2] for seg in schedule),
             time_ns(.6+32*.002,scen.scale_mhz),places=7)
        self.assertGreater(np.linalg.norm(schedule[1][0]-schedule[1][1]),1e-4)
        self.assertLess(np.linalg.norm(schedule[0][0]-schedule[0][1]),1e-12)

    def test_clean_channel_reproduces_two_matched_references(self):
        scen=Scenario(t1_ns=float('inf'),tphi_ns=float('inf'),
                pointer_t1_ns=float('inf'),pointer_read_ns=0.,
                pointer_assignment_error=0.,pointer_reset_ns=0.,shots=1)
        schedule=propagation_schedule(self.states,self.H,scen)
        conf=pointer_confusion(scen)
        for arm in ('analog','uncontrolled'):
            y=self.K.copy()
            for active,idle,_ in schedule:y=(active if arm=='analog' else idle)@y
            target=entanglement_fidelity([y],self.K,self.V)
            computed,weight,jumps=run_one(self.K.copy(),schedule,arm,
                    scen,self.V,conf,np.random.default_rng(81))
            self.assertAlmostEqual(computed,target,places=7)
            self.assertGreaterEqual(weight,0.)
            self.assertGreaterEqual(jumps,0)

    def test_study_reports_distinct_arms_and_limits(self):
        res=study(Scenario(shots=2,seed=123))
        self.assertEqual([r['arm'] for r in res['results']],
                         ['analog','uncontrolled'])
        self.assertEqual(res['meta']['control_pulse_windows'],32)
        self.assertLess(res['meta']['reference_difference_fro_max'],2e-9)
        self.assertEqual(res['meta']['syndrome_sector_ranks'],[2]*5+[233])
        for row in res['results']:
            self.assertEqual(len(row['noisy_pointer_channel_samples']),2)
            self.assertTrue(0<=row['noisy_F_e_mean']<=1)


if __name__=='__main__':
    unittest.main()
