"""Regression suite: isotropic SU3 quantum-conditional local charged exchange."""
import unittest
import numpy as np
from colour_gauge_link_s3 import G,DIM,REP,physical_isometry
from colour_s3_edge_logical_syndrome import code_isometry,pointer_circuit
from colour_s3_charged_bus_transport import (
    initialize,gauge_action,local_swap,gauss_acceptance,remove_ancilla_singlet,
)
from colour_s3_native_exchange_hamiltonian import (
    GENS,I9,SWAP9,H_EX,H_ANISO,verify_exchange_identity,
    su3_commutant_neutral_pointer,controlled_hamiltonian,
    local_two_qutrit_unitary,apply_local_hamiltonian,
    syndrome_hamiltonian_circuit,feedback_hamiltonian_circuit,
    full_SU3_pulse_covariance_witness,timing_scenario,fidelity,study,
)

class NativeExchangeHamiltonianTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.W=physical_isometry()
        cls.K,cls.F,_=code_isometry(cls.W)

    def test_eight_su3_generator_exchange_is_swap_projector(self):
        self.assertEqual(len(GENS),8)
        self.assertLess(verify_exchange_identity(),2e-14)
        self.assertLess(np.linalg.norm(H_EX@H_EX-H_EX),1e-12)
        for T in GENS:
            tot=np.kron(T,np.eye(3))+np.kron(np.eye(3),T)
            self.assertLess(np.linalg.norm(H_EX@tot-tot@H_EX),1e-12)

    def test_quantum_controlled_H_exact_fredkin(self):
        H=controlled_hamiltonian((1,2))
        self.assertEqual(H.shape,(27,27))
        self.assertLess(np.linalg.norm(H-H.conj().T),1e-12)
        U=np.kron(np.eye(3),I9)
        U-=np.kron(np.diag([0,1,1]),I9-SWAP9)
        # exp(-i pi P (I-S)/2) equals pointer-controlled SWAP
        from scipy.linalg import expm
        self.assertLess(np.linalg.norm(expm(-1j*np.pi*H)-U),2e-12)

    def test_full_su3_neutral_pointer_two_body_schur_obstruction(self):
        data=su3_commutant_neutral_pointer()
        self.assertEqual(data["neutral_pointer_fundamental_invariant_complex_dimension"],9)
        self.assertEqual(data["fundamental_quark_pair_invariant_complex_dimension"],2)
        self.assertFalse(data["generated_neutral_pointer_to_quark_entanglement"])

    def test_embedded_hamiltonian_matches_legacy_charged_bus_swap(self):
        rng=np.random.default_rng(8)
        x=rng.normal(size=(3,DIM,2,3,3,2))
        x=x+1j*rng.normal(size=x.shape)
        for site in ("A0","A1","B0"):
            for active in ((1,2),(2,),()):
                actual=apply_local_hamiltonian(x,site,active)
                expected=local_swap(x,site,active)
                self.assertLess(np.linalg.norm(actual-expected),3e-11)
                self.assertLess(abs(np.vdot(actual,actual)-np.vdot(x,x)),2e-10)

    def test_continuous_pulse_gauge_covariant_on_arbitrary_offsector_inputs(self):
        rng=np.random.default_rng(97)
        z=rng.normal(size=(3,DIM,2,3,3,1))+1j*rng.normal(
            size=(3,DIM,2,3,3,1))
        for site in ("A0","A1","B0"):
            for theta in (.11,1.23,np.pi):
                v=apply_local_hamiltonian(z,site,(1,2),theta)
                for vertex in ("A","B"):
                    for h in range(G):
                        left=gauge_action(v,h,vertex)
                        right=apply_local_hamiltonian(
                            gauge_action(z,h,vertex),site,(1,2),theta)
                        self.assertLess(np.linalg.norm(left-right),1e-10)

    def test_full_gauss_path_and_logical_pointer_exact(self):
        out,steps=syndrome_hamiltonian_circuit(self.K,return_snapshots=True)
        self.assertEqual(len(steps),16)
        obtained=remove_ancilla_singlet(out)
        self.assertLess(np.linalg.norm(obtained-pointer_circuit(self.K)),2e-10)
        for _,state in steps:
            self.assertGreater(min(gauss_acceptance(state).values()),1-1e-10)

    def test_anisotropic_control_breaks_full_su3_and_finite_gauge(self):
        clean=full_SU3_pulse_covariance_witness(0.)
        drift=full_SU3_pulse_covariance_witness(.01)
        self.assertLess(clean["SU3_generators_max_commutator_fro"],1e-11)
        self.assertLess(clean["S3_group_max_commutator_fro"],1e-11)
        self.assertGreater(drift["SU3_generators_max_commutator_fro"],.001)
        self.assertGreater(drift["S3_group_max_commutator_fro"],.001)

    def test_unconditional_recovered_channel_and_calibration_sensitivity(self):
        ideal=fidelity(self.K,0.)
        drift=fidelity(self.K,.05)
        self.assertAlmostEqual(ideal["unconditional_entanglement_fidelity"],1.,places=10)
        self.assertAlmostEqual(drift["total_output_probability"],1.,places=10)
        self.assertTrue(0<=drift["unconditional_entanglement_fidelity"]<1)

    def test_timing_is_hypothetical_and_complete_audit(self):
        r=timing_scenario(J_mhz=5.,transport_ns=100.)
        self.assertEqual(r["assumed_effective_control_pi_pulse_ns"],100.)
        self.assertEqual(r["syndrome_excluding_pointer_and_reset_ns"],1600.)
        data=study()
        self.assertTrue(data["Hamiltonian_resources"]
                        ["coherent_native_three_body_coupling_required"])
        self.assertEqual(data["model"].split(",")[0],"finite S3 charged-bus simulator")
        self.assertEqual(len(data["full_logical_channel_anisotropy_sweep"]),3)

if __name__=="__main__":
    unittest.main()
