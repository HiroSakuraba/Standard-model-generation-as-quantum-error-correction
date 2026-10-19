"""Regression tests for four-level transmon CP22 noise and SU3 analog/digital audit."""
import unittest
from dataclasses import replace
import numpy as np

from su3_transmon_noise_probe import (
    PairParameters, D, EYE, EYE_PAIR, TOTAL_EXCITATION, PROJECT,
    duffing_hamiltonian, smooth_pulse, local_phase_frame, ideal_cp22,
    loss_kraus, dephasing_kraus, single_site_noise, two_site_noise,
    channel_metrics, pulse_candidate
)
from su3_analog_digital_noise_comparison import (
    dimensional_mapping, syndrome_scenarios, analog_flagged_replacement
)


class DeviceNoiseTests(unittest.TestCase):
    def test_excitation_conserving_duffing_pulse_is_unitary(self):
        p=PairParameters(duration_ns=32,steps=12)
        h=duffing_hamiltonian(p,290.)
        self.assertLess(np.linalg.norm(h@TOTAL_EXCITATION-TOTAL_EXCITATION@h),1e-10)
        u=smooth_pulse(p)
        self.assertLess(np.linalg.norm(u.conj().T@u-EYE_PAIR),1e-10)

    def test_zero_exchange_only_local_phase(self):
        p=replace(PairParameters(),exchange_mhz=0,duration_ns=24,steps=12)
        u,phi=local_phase_frame(smooth_pulse(p))
        self.assertLess(np.linalg.norm(u-EYE_PAIR),1e-9)
        self.assertLess(abs(phi),1e-8)

    def test_loss_and_dephasing_kraus_complete(self):
        for time in (0,50.,1000.):
            for fn,args in ((loss_kraus,(time,30000.)),(dephasing_kraus,(time,45000.))):
                kraus=fn(*args)
                total=sum((a.conj().T@a for a in kraus),np.zeros((D,D),complex))
                self.assertLess(np.linalg.norm(total-EYE),2e-11)

    def test_physical_channel_identity_at_zero_duration(self):
        ks=two_site_noise(0.,30000.,45000.)
        m=channel_metrics(ideal_cp22(np.pi),np.pi,ks)
        self.assertAlmostEqual(m['unconditional_entanglement_fidelity'],1.,places=10)
        self.assertLess(m['final_fourth_level_leakage_probability'],1e-10)

    def test_finite_channels_reduce_ideal_process_fidelity(self):
        ks=two_site_noise(1000.,30000.,45000.)
        m=channel_metrics(ideal_cp22(np.pi),np.pi,ks)
        self.assertLess(m['unconditional_entanglement_fidelity'],1.)
        self.assertGreater(m['unconditional_entanglement_fidelity'],0.)
        self.assertLess(m['final_fourth_level_leakage_probability'],1e-9)

    def test_pulse_metrics_are_physical_unconditional(self):
        p=PairParameters(duration_ns=40.,steps=16)
        m=pulse_candidate(p,include_noise=True)
        self.assertGreaterEqual(m['coherent_leakage_probability'],0)
        self.assertLessEqual(m['coherent_target_F_e'],1.)
        self.assertLessEqual(m['post_pulse_T1_Tphi_target_F_e'],1.)

    def test_syndrome_lower_bound_and_control_scaling(self):
        m=dimensional_mapping(scale_mhz=.01,tau=.002)
        self.assertEqual(m['control_pulses_global'],32)
        self.assertEqual(m['addressed_face_clock_operations'],128)
        self.assertAlmostEqual(m['pulse_width_ns'],.002/(2*np.pi*.001*.01))
        self.assertEqual(syndrome_scenarios()[0]['minimum_distinct_syndrome_outcomes'],6)
        self.assertEqual(syndrome_scenarios()[0]['minimum_classical_bits'],3)

    def test_analog_toy_replacement_monotone(self):
        a=analog_flagged_replacement(.999,0)
        b=analog_flagged_replacement(.999,.01)
        self.assertGreater(a['logical_replacement_F_e'],b['logical_replacement_F_e'])


if __name__=='__main__':
    unittest.main()
