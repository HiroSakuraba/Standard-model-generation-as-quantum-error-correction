"""CPTP four-level pointer and matched-time finite-S3 benchmark regressions."""
import unittest
from dataclasses import replace

import numpy as np

from colour_s3_recovery_pulse_compiler import TOTAL, target_three_outcome
from colour_s3_open_system_instrument import (
    FULL, NoiseScenario, embed_pointer_ququart, pointer_leakage_unitary,
    channel_branch, apply_pointer, electric_phase_kraus,
    run_instrument, compare_policies, _readout_acceptance,
)
from colour_s3_branch_aware_optimizer import BranchAwareCompiler


class OpenSystemFiniteS3Tests(unittest.TestCase):
    def test_four_level_embedding_is_exact_and_unitary(self):
        for label in ("AS","BS"):
            _,ideal,_,_,_=target_three_outcome(label)
            u=embed_pointer_ququart(ideal)
            self.assertEqual(u.shape,(56,56))
            self.assertLess(np.linalg.norm(
                u.conj().T@u-np.eye(FULL),"fro"),2e-10)
            for i in range(14):
                index=4*i+3
                self.assertAlmostEqual(u[index,index].real,1.,places=12)
            with self.assertRaises(ValueError):
                embed_pointer_ququart(np.eye(14))

    def test_leakage_rotation_is_unitary_and_only_mixes_2_3(self):
        u=pointer_leakage_unitary(.15)
        self.assertLess(np.linalg.norm(
            u.conj().T@u-np.eye(4)),1e-12)
        self.assertGreater(abs(u[3,2]),0.)
        self.assertAlmostEqual(u[0,0],1.,places=12)
        self.assertAlmostEqual(u[1,1],1.,places=12)
        with self.assertRaises(ValueError):
            pointer_leakage_unitary(-.1)

    def test_complete_channel_kraus_trace_preservation(self):
        sc=NoiseScenario()
        E=np.diag([1.]*7+[0.]*7)
        R=np.eye(14)-2*E
        initial=np.zeros((56,4),complex)
        initial[0,0]=2.  # trace convention squared norm / 4 = 1
        x=channel_branch([initial],1000.,sc,R)
        self.assertAlmostEqual(sum(
            np.linalg.norm(k,"fro")**2/4 for k in x),1.,places=9)
        self.assertTrue(len(x)>1)
        self.assertLess(np.linalg.norm(
            sum((k.conj().T@k for k in electric_phase_kraus(
                1000.,sc.electric_Tphi_ns,R)),np.zeros((14,14),complex))
            -np.eye(14)),1e-10)
        with self.assertRaises(ValueError):
            channel_branch([initial],-1,sc,R)

    def test_noiseless_ideal_exactly_corrects_single_fault_mixture(self):
        clean=NoiseScenario(
            pointer_T1_ns=1e15, pointer_Tphi_ns=1e15,
            electric_Tphi_ns=1e15, active_particle_T1_ns=1e15,
            pointer_leak_angle_per_window=0.,
            readout_wrong_label_probability=0.,
            readout_ns=0.,reset_ns=0.)
        for label in ("AS","BS"):
            _,ideal,_,_,_=target_three_outcome(label)
            q=run_instrument(label,[ideal],0.,clean)
            self.assertAlmostEqual(
                q["unconditional_logical_F_e_including_flagged_erasure"],
                1.,places=9)
            self.assertAlmostEqual(
                q["mean_pointer_fourth_level_population"],0.,places=9)
            for row in q["fault_details"]:
                self.assertAlmostEqual(row["total_Kraus_trace"],1.,places=9)

    def test_do_nothing_reference_and_no_free_device_claim(self):
        clean=NoiseScenario(
            pointer_T1_ns=1e15,pointer_Tphi_ns=1e15,
            electric_Tphi_ns=1e15, active_particle_T1_ns=1e15,
            readout_wrong_label_probability=0.,
            pointer_leak_angle_per_window=0.,
            readout_ns=0.,reset_ns=0.)
        for label in ("AS","BS"):
            q=run_instrument(label,[],0.,clean)
            self.assertAlmostEqual(
                q["unconditional_logical_F_e_including_flagged_erasure"],
                .63+.37/9,places=9)
            self.assertTrue(q["no_device_specific_CP22_realization"])
            self.assertTrue(q["no_ideal_code_projection_applied_as_a_gate"])

    def test_measurement_label_confusion_is_normalized(self):
        self.assertEqual(_readout_acceptance(3,.05),0.)
        self.assertAlmostEqual(_readout_acceptance(0,.2),.9)
        self.assertAlmostEqual(_readout_acceptance(1,.2),.9)
        self.assertAlmostEqual(_readout_acceptance(2,.2),.2)
        for b in (0,1,2):
            self.assertTrue(0<=_readout_acceptance(b,.001)<=1)

    def test_matched_passive_pays_the_same_elapsed_particle_loss(self):
        comp=BranchAwareCompiler("AS",1)
        zeros=np.zeros(comp.n_angles)
        sc=NoiseScenario(readout_ns=300.,reset_ns=250.,
                         pointer_leak_angle_per_window=0.)
        r=compare_policies("AS",zeros,sc)
        self.assertAlmostEqual(
            r["candidate"]["total_elapsed_ns"],
            r["same_elapsed_time_passive"]["total_elapsed_ns"])
        self.assertAlmostEqual(
            r["candidate"]["four_active_particle_no_erasure_probability"],
            r["same_elapsed_time_passive"][
                "four_active_particle_no_erasure_probability"])
        self.assertTrue(r["candidate"]["no_postselection"])
        self.assertFalse(r["candidate_beats_matched_passive"])

    def test_positive_noise_and_readout_cost_cannot_improve_ideal_decode(self):
        _,ideal,_,_,_=target_three_outcome("AS")
        sc=NoiseScenario()
        a=run_instrument("AS",[ideal],0.,sc)
        b=run_instrument("AS",[ideal],0.,replace(sc,
            readout_ns=0.,reset_ns=0.,
            pointer_leak_angle_per_window=0.,
            readout_wrong_label_probability=0.))
        self.assertLessEqual(
            a["unconditional_logical_F_e_including_flagged_erasure"],
            b["unconditional_logical_F_e_including_flagged_erasure"]+1e-9)
        self.assertLess(a["four_active_particle_no_erasure_probability"],1)

    def test_nontrivial_pulse_product_is_audited_without_postselection(self):
        comp=BranchAwareCompiler("AS",1)
        a=np.zeros(comp.n_angles)
        a[4]=.5
        a[6]=.4
        r=compare_policies("AS",a,NoiseScenario())
        self.assertEqual(r["control_angles"],13)
        self.assertGreater(r["candidate"]["pulse_time_ns"],0.)
        self.assertAlmostEqual(
            r["candidate"]["total_elapsed_ns"],
            r["same_elapsed_time_passive"]["total_elapsed_ns"])
        for row in r["candidate"]["fault_details"]:
            self.assertAlmostEqual(row["total_Kraus_trace"],1.,places=8)
            self.assertLessEqual(row[
                "F_e_before_active_particle_erasure"],
                row["accepted_code_probability"]+1e-10)

    def test_invalid_noise_rejected(self):
        for changed in (
            {"electric_probability":1.3},
            {"pointer_T1_ns":0.},
            {"pointer_Tphi_ns":float("nan")},
            {"active_particle_T1_ns":-1.},
            {"readout_wrong_label_probability":1.2},
            {"pointer_leak_angle_per_window":-1.},
            {"reset_ns":float("inf")},
        ):
            with self.assertRaises(ValueError):
                NoiseScenario(**changed).validate()
        with self.assertRaises(ValueError):
            run_instrument("AS",[np.eye(TOTAL)],-1.)
        with self.assertRaises(ValueError):
            run_instrument("AS",[np.eye(TOTAL)],0.,nominal_pulse_count=2)


if __name__=="__main__":
    unittest.main()
