"""Regression suite: microscopic Hubbard feedback, exact loss, full-Fock dephasing."""
import unittest
import numpy as np

from colour_gauge_link_s3 import physical_isometry
from colour_s3_edge_logical_syndrome import code_isometry
from colour_su3_hubbard_virtual_exchange import optimize_pulse
from colour_su3_hubbard_persistent_syndrome import (
    N_FOCK,initialize,gauss_acceptance,exact_full_channel
)
from colour_su3_hubbard_loss_dephase_feedback import (
    N_MODES,OCC,ORBITAL_OCC,number_phase_kick,time_unit_ns,
    flagged_uniform_loss,full_no_jump_branch,process_sample,run
)


class MicroscopicFeedbackNoiseTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        W=physical_isometry()
        cls.K,_,_=code_isometry(W)
        opt=optimize_pulse(30.)
        cls.U=30.
        cls.delta=opt["active_detuning"]
        cls.T=opt["time"]

    def test_495D_occupation_basis_stays_four_fermions(self):
        self.assertEqual(OCC.shape,(495,12))
        self.assertEqual(ORBITAL_OCC.shape,(495,4))
        self.assertEqual(N_MODES,12)
        self.assertTrue(np.all(OCC.sum(axis=1)==4))
        self.assertTrue(np.all(ORBITAL_OCC.sum(axis=1)==4))

    def test_stochastic_number_dephasing_is_unitary_and_seeded(self):
        state=initialize(self.K)
        for model in ("orbital_charge","colour_resolved"):
            noise_a=number_phase_kick(state,2.5,.01,
                         np.random.default_rng(12),model)
            noise_b=number_phase_kick(state,2.5,.01,
                         np.random.default_rng(12),model)
            self.assertLess(np.linalg.norm(noise_a-noise_b),1e-12)
            self.assertLess(abs(np.vdot(noise_a,noise_a)
                      -np.vdot(state,state)),1e-10)
            self.assertLess(np.linalg.norm(number_phase_kick(
                state,1.,0.,np.random.default_rng(13),model)-state),1e-10)

    def test_orbital_charge_dephasing_preserves_extended_gauss(self):
        state=initialize(self.K)
        noise=number_phase_kick(state,5.,.1,np.random.default_rng(22),
                                "orbital_charge")
        rates=gauss_acceptance(noise)
        self.assertGreater(min(rates.values()),1-1e-10)

    def test_colour_resolved_dephasing_breaks_gauss_in_general(self):
        state=initialize(self.K)
        noise=number_phase_kick(state,5.,.1,np.random.default_rng(22),
                                "colour_resolved")
        rates=gauss_acceptance(noise)
        self.assertTrue(any(v<.99 for v in rates.values()))

    def test_uniform_fermion_loss_flagged_erasure_is_cptp(self):
        for dt in (0.,1000.,11000.):
            data=flagged_uniform_loss(dt,30.)
            self.assertAlmostEqual(data["no_jump_probability"]+
                                   data["flagged_particle_loss_probability"],1.,places=12)
            self.assertAlmostEqual(data["no_jump_probability"],
                                   np.exp(-4*dt/30000.),places=12)
        self.assertEqual(flagged_uniform_loss(0.)["no_jump_probability"],1.)
        with self.assertRaises(ValueError):flagged_uniform_loss(-1.)
        with self.assertRaises(ValueError):flagged_uniform_loss(1.,-2.)

    def test_microscopic_feedback_no_loss_full_channel(self):
        item=full_no_jump_branch(self.K,self.U,self.delta,self.T,
                                 gamma=0.,substeps=1)
        original=exact_full_channel(self.K,self.U,self.delta,self.T)
        self.assertAlmostEqual(item["full_no_jump_trace"],1.,places=10)
        self.assertTrue(0<item["F_e_conditional_on_no_particle_loss"]<1.)
        # Real Hubbard feedback is no longer the ideal computational
        # subspace reference; it must not be silently substituted.
        self.assertGreater(abs(item["F_e_conditional_on_no_particle_loss"]-
                               original["unconditional_F_e_ideal_feedback"]),1e-5)

    def test_divided_pulse_matches_unsplit_when_dephasing_off(self):
        one=full_no_jump_branch(self.K,self.U,self.delta,self.T,
                                gamma=0.,substeps=1)
        two=full_no_jump_branch(self.K,self.U,self.delta,self.T,
                                gamma=0.,substeps=2)
        self.assertLess(abs(one["F_e_conditional_on_no_particle_loss"]-
                            two["F_e_conditional_on_no_particle_loss"]),2e-10)

    def test_loss_not_postselected_and_seeded_phase_trajectories(self):
        dry=process_sample(self.K,self.U,self.delta,self.T,
                           gamma=0.,shots=1,seed=14)
        loss=dry["loss_channel"]
        self.assertAlmostEqual(
            dry["mean_F_e_unconditional_including_flagged_loss"],
            dry["mean_F_e_conditional_on_no_loss"]*loss["no_jump_probability"],
            places=12)
        self.assertLess(dry["mean_F_e_unconditional_including_flagged_loss"],
                        dry["mean_F_e_conditional_on_no_loss"])
        noisy=process_sample(self.K,self.U,self.delta,self.T,
                            gamma=.0001,model="orbital_charge",
                            shots=2,seed=45,substeps=2,loss_enabled=False)
        again=process_sample(self.K,self.U,self.delta,self.T,
                            gamma=.0001,model="orbital_charge",
                            shots=2,seed=45,substeps=2,loss_enabled=False)
        self.assertEqual(noisy["individual_F_e_conditional"],
                         again["individual_F_e_conditional"])
        self.assertTrue(0<=noisy["mean_F_e_conditional_on_no_loss"]<=1.)
        self.assertGreater(noisy["mean_final_Gauss_acceptance_conditional"]["A"],1-1e-10)


if __name__=="__main__":
    unittest.main()
