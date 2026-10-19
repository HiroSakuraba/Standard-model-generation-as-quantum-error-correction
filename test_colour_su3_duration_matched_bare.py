"""Reproducible direct four-state bare reference and lifetime break-even tests."""
import unittest

import numpy as np

from colour_su3_duration_matched_bare import (
    SOURCE_ERROR_PROBABILITY, DLOG, bare_colour_dephasing_fidelity,
    bare_reference, exactly_flagged_fermion_survival,
    required_T1_crossover_us, physical_duration_ns, analyze_anchor,
    duration_matched_grid, noisy_colour_bare_comparison, run,
)
from colour_gauge_link_s3 import physical_isometry, REP
from colour_s3_edge_logical_syndrome import code_isometry


class DurationMatchedBareTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        W=physical_isometry()
        cls.K,_,_=code_isometry(W)

    def test_duration_includes_microscopic_feedback_and_readout(self):
        t=physical_duration_ns(47.,hopping_mhz=10.,transport_ns=100.,
                               readout_ns=300.)
        self.assertAlmostEqual(t,15*47*1000/(2*np.pi*10)+900.,places=9)
        for args in [(-1.,10.,100.,300.),(47.,-1.,100.,300.),
                     (47.,10.,-10.,300.),(47.,10.,100.,-1.)]:
            with self.assertRaises(ValueError):
                physical_duration_ns(args[0],args[1],args[2],args[3])

    def test_analytic_colour_number_lindblad_against_exact_kraus(self):
        # One qubit in two occupied colour modes, Hermitian occupation
        # dephasing gives Z dephasing strength c=e^-2 gamma time.
        # Two independent qubits encode d=4; verify full 4x4 Kraus
        # entanglement fidelity and traceless logical-X channel explicitly.
        gamma=.0002
        elapsed=700.
        c=np.exp(-2*gamma*elapsed)
        I=np.eye(2,dtype=complex)
        Z=np.diag([1.,-1.])
        X=np.array([[0.,1.],[1.,0.]],complex)
        q=[np.sqrt((1+c)/2)*I,np.sqrt((1-c)/2)*Z]
        two=[np.kron(a,b) for a in q for b in q]
        Fe=sum(abs(np.trace(u)/4)**2 for u in two)
        FeX=sum(abs(np.trace(u@np.kron(X,I))/4)**2 for u in two)
        self.assertAlmostEqual(Fe,bare_colour_dephasing_fidelity(
            elapsed,gamma,"colour_resolved"),places=12)
        self.assertAlmostEqual(FeX,0.,places=12)
        self.assertAlmostEqual(bare_colour_dephasing_fidelity(
            elapsed,gamma,"orbital_charge"),1.,places=12)

    def test_bare_logical_fault_is_trace_zero_channel(self):
        x=bare_reference(12000.,T1_us=30.,gamma=.0002,
                         model="colour_resolved",source_p=.37)
        self.assertEqual(x["equal_logical_dimension"],4)
        self.assertEqual(x["active_fermions"],2)
        self.assertAlmostEqual(x["matched_probability_X_fault_unconditional_F_e"],
            .63*x["no_fault_unconditional_F_e"],places=12)
        self.assertLess(x["matched_probability_X_fault_unconditional_F_e"],
                        x["no_fault_unconditional_F_e"])

    def test_loss_reference_is_exact_two_versus_four_fermion_erasure(self):
        time=12000.
        bare=exactly_flagged_fermion_survival(2,time,30.)
        coded=exactly_flagged_fermion_survival(4,time,30.)
        self.assertAlmostEqual(bare**2,coded,places=12)
        self.assertAlmostEqual(bare,np.exp(-2*.012/.03),places=12)
        with self.assertRaises(ValueError):
            exactly_flagged_fermion_survival(2,time,0.)

    def test_no_crossover_against_ideal_unfaulted_bare_memory(self):
        for time in (100.,10000.,20000.):
            for F in (.5,.8,1.):
                self.assertIsNone(required_T1_crossover_us(F,1.,time,0.))

    def test_crossover_matches_analytic_channel_equality(self):
        F=.91; T=12500.
        threshold=required_T1_crossover_us(F,1.,T,.37)
        self.assertGreater(threshold,0.)
        ratio=np.exp(-2*T/(1000*threshold))
        self.assertAlmostEqual(F*ratio**2,.63*ratio,places=12)
        for lifetime in (.99*threshold,1.01*threshold):
            encoded=F*exactly_flagged_fermion_survival(4,T,lifetime)
            bare=.63*exactly_flagged_fermion_survival(2,T,lifetime)
            self.assertEqual(encoded>bare,lifetime>threshold)
        self.assertIsNone(required_T1_crossover_us(.6,1.,T,.37))
        self.assertEqual(required_T1_crossover_us(.8,1.,T,1.),0.)

    def test_bare_two_colour_subspace_is_not_an_S3_gauss_code(self):
        # Two-dimensional colour-basis storage is an external direct
        # hardware reference, not a gauge-invariant finite-S3 code.
        P=np.diag([1.,1.,0.])
        residual=max(np.linalg.norm(D@P-P@D) for D in REP)
        self.assertGreater(residual,0.5)
        sample=bare_reference(1000.)
        self.assertFalse(sample["bare_subspace_in_local_S3_Gauss_sector"])
        self.assertTrue(sample["bare_external_hardware_reference_not_gauge_admissible"])

    def test_merged_PR57_anchor_U30_is_recomputed(self):
        row=analyze_anchor(30.,self.K)
        self.assertAlmostEqual(row["encoded_microscopic_feedback_no_loss_F_e"],
                               .9119196563716021,places=5)
        self.assertAlmostEqual(row["encoded_unconditional_F_e_at_T1_30us"],
                               .1802217616023003,places=5)
        self.assertAlmostEqual(row["bare_no_fault_unconditional_F_e_at_T1_30us"],
                               .4445546,places=5)
        self.assertAlmostEqual(row["break_even_T1_us_vs_37pct_bare_X_fault_colourblind"],
                               65.761,places=2)
        self.assertLess(row["encoded_vs_bare_no_fault_difference_at_T1_30us"],0.)
        self.assertLess(row["encoded_vs_challenged_bare_difference_at_T1_30us"],0.)

    def test_sensitivity_grid_exposes_all_comparisons(self):
        anchors=[{"U_over_t":20.,"encoded_microscopic_feedback_no_loss_F_e":.77,
                  "optimized_15state_pulse_inverse_t":31.41}]
        grid=duration_matched_grid(anchors)
        self.assertEqual(len(grid),4*3*5)
        self.assertFalse(any(x["encoded_beats_clean_bare"] for x in grid))
        # A sufficiently faster device and longer T1 can make
        # nontrivial challenging channel correctable, without claiming
        # the same microscopic physical process for bare/code.
        self.assertTrue(any(x["encoded_beats_challenged_bare"] for x in grid))
        self.assertTrue(any(not x["encoded_beats_challenged_bare"] for x in grid))

    def test_colour_dephasing_baseline_stays_analytic(self):
        reports=noisy_colour_bare_comparison(self.K,U=30.,gamma=.0002,
                                              shots=2)
        self.assertEqual(len(reports),2)
        self.assertEqual(reports[0]["model"],"orbital_charge")
        self.assertEqual(reports[1]["model"],"colour_resolved")
        self.assertGreater(reports[0]["bare_no_fault_unconditional_F_e"],
                           reports[1]["bare_no_fault_unconditional_F_e"])
        for r in reports:
            self.assertTrue(0<=r["encoded_unconditional_loss_dephasing_F_e"]<=1)
            self.assertTrue(0<=r["bare_no_fault_unconditional_F_e"]<=1)
            self.assertTrue(r["same_gamma_is_not_equal_total_noise_power"])

if __name__=="__main__":
    unittest.main()
