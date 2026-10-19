"""Regression tests for full 15-state, two-site SU3 Hubbard exchange probe."""
import unittest
import numpy as np

from colour_su3_hubbard_virtual_exchange import (
    FOCK,W,SWAP,Q_ANTISYM,HOP,DOUBLON,IMBALANCE,
    one_body,hamiltonian,su3_total_generator,invariant_residual,
    effective_exchange,projected_evolution,gate_metrics,
    perturbative_schedule,optimize_pulse,max_transient_leakage,
    illustrative_loss_budget,run
)
from colour_singlet_code import gellmann


class MicroscopicHubbardTests(unittest.TestCase):
    def test_correct_fermionic_anticommutation_and_fock_dimensions(self):
        self.assertEqual(len(FOCK),15)
        self.assertEqual(W.shape,(15,9))
        self.assertEqual(round(np.trace(DOUBLON).real),6)
        self.assertLess(np.linalg.norm(W.conj().T@W-np.eye(9)),1e-12)
        self.assertLess(np.linalg.norm(HOP-HOP.conj().T),1e-12)
        self.assertLess(np.linalg.norm(SWAP@SWAP-np.eye(9)),1e-12)
        # both hopping and doublons included; one-body diagonal is n_colour
        for i in range(6):
            self.assertEqual(round(np.trace(one_body(i,i)).real),5)

    def test_full_same_vertex_SU3_colour_invariance(self):
        for delta in (0.,5.,10.):
            H=hamiltonian(20.,delta)
            self.assertLess(np.linalg.norm(H-H.conj().T),1e-12)
            for T in gellmann():
                G=su3_total_generator(T)
                self.assertLess(np.linalg.norm(G-G.conj().T),1e-12)
                self.assertLess(np.linalg.norm(H@G-G@H),2e-12)
            self.assertLess(invariant_residual(20.,delta),1e-11)

    def test_no_native_pointer_controlled_swap_is_assumed(self):
        # Within one particle per orbital, the site detuning is zero.
        self.assertLess(np.linalg.norm(W.conj().T@IMBALANCE@W),1e-12)
        # All actual exchange is through off-block virtual hopping.
        self.assertLess(np.linalg.norm(W.conj().T@HOP@W),1e-12)
        self.assertGreater(np.linalg.norm(DOUBLON@HOP@W),1.)

    def test_exact_schrieffer_wolff_antisymmetric_exchange(self):
        for U in (20.,30.,40.):
            for delta in (0.,U/np.sqrt(3)):
                J,H,err=effective_exchange(U,delta)
                self.assertLess(err,1e-12)
                self.assertLess(np.linalg.norm(H+J*Q_ANTISYM),1e-12)
                self.assertEqual(round(np.linalg.matrix_rank(H)),3)

    def test_pointer_ratio_two_pi_off_three_pi_on(self):
        for U in (20.,30.,40.):
            p=perturbative_schedule(U)
            self.assertAlmostEqual(p["on_over_off_exchange_ratio"],1.5,places=12)
            self.assertAlmostEqual(p["off_exchange_rate"]*p["time_in_inverse_t_angular_units"],
                                   2*np.pi,places=11)
            self.assertAlmostEqual(p["on_exchange_rate"]*p["time_in_inverse_t_angular_units"],
                                   3*np.pi,places=11)

    def test_exact_15_state_gate_includes_leakage_not_postselected(self):
        proposal=perturbative_schedule(20.)
        met=gate_metrics(20.,proposal["time_in_inverse_t_angular_units"],
                         proposal["on_pointer_detuning"])
        self.assertGreater(met["entanglement_fidelity"],.97)
        self.assertLess(met["entanglement_fidelity"],1.)
        self.assertGreater(met["on_pointer_leakage"],0.)
        self.assertLess(met["average_final_doublon_leakage"],.02)

    def test_bounded_optimizer_does_not_worsen_exact_gate(self):
        for U in (20.,30.):
            r=optimize_pulse(U)
            self.assertGreaterEqual(r["optimized_metrics"]["entanglement_fidelity"]+1e-8,
                                    r["baseline_metrics"]["entanglement_fidelity"])
            self.assertLess(r["optimized_metrics"]["entanglement_fidelity"],1.)
            self.assertGreater(r["time"],0.)
            self.assertLess(abs(r["active_detuning"]),U)

    def test_transient_doublons_and_illustrative_loss_budget(self):
        r=optimize_pulse(30.)
        T,d=r["time"],r["active_detuning"]
        trace=max_transient_leakage(30.,d,T,steps=31)
        self.assertGreater(trace["max_transient_doublon_occupation"],
                           trace["final_doublon_occupation"])
        self.assertLess(trace["max_transient_doublon_occupation"],.1)
        loss=illustrative_loss_budget(r["optimized_metrics"]["entanglement_fidelity"],T)
        self.assertGreater(loss["pulse_duration_ns"],700.)
        self.assertLess(loss["pulse_duration_ns"],800.)
        self.assertLess(loss["leakage_and_loss_survival_weighted_F_e"],
                        r["optimized_metrics"]["entanglement_fidelity"])

    def test_invalid_unphysical_hubbard_parameters(self):
        for U,delta,t in ((0.,0.,1.),(20.,20.,1.),(20.,0.,-1.)):
            with self.assertRaises(ValueError):
                hamiltonian(U,delta,t)
        with self.assertRaises(ValueError):
            projected_evolution(20.,0.,-1.)

    def test_full_reproducible_report_and_scientific_scope(self):
        data=run()
        self.assertEqual(len(data["benchmark_arms"]),3)
        self.assertEqual(data["virtual_doublon_hole_dimension"],6)
        self.assertEqual(data["one_per_orbital_colour_subspace_dimension"],9)
        self.assertTrue(any("ONE vertex" in s for s in data["limitations"]))
        for arm in data["benchmark_arms"]:
            self.assertGreater(arm["optimization"]["optimized_metrics"]["entanglement_fidelity"],.98)


if __name__=="__main__":
    unittest.main()
