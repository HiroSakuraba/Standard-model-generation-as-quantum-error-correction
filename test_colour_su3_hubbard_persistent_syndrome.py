"""Persistent doublon leakage and gauge constraints through full charged-bus syndrome."""
import unittest
import numpy as np
from colour_gauge_link_s3 import G,physical_isometry,dressed_exchange
from colour_s3_edge_logical_syndrome import code_isometry
from colour_su3_hubbard_virtual_exchange import projected_evolution,optimize_pulse
from colour_su3_hubbard_persistent_syndrome import (
    F_STATES,F_INDEX,N_FOCK,fock_index,subgroup_wedge_map,
    local_hamiltonian_blocks,local_fock_pulse,initialize,
    link_transport,apply_microscopic_gate,microscopic_schedule,
    gauge_transform,gauss_acceptance,occupation_leakage,
    exact_full_channel,run,
)

class PersistentHubbardTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        W=physical_isometry()
        cls.K,cls.F,_=code_isometry(W)

    def test_complete_four_orbital_fermion_Fock_space(self):
        self.assertEqual(N_FOCK,495)
        self.assertEqual(len(F_STATES),495)
        self.assertEqual(len(F_INDEX),495)
        self.assertEqual(len({fock_index(*t) for t in __import__('itertools').product(range(3),repeat=4)}),81)
        self.assertEqual(len(microscopic_schedule()),16)

    def test_two_orbital_block_propagation_is_exact_unitary(self):
        for site in ('A0','A1','B0'):
            U=local_fock_pulse(site,20.,20/np.sqrt(3),np.pi*20/2)
            self.assertLess(np.linalg.norm((U.getH()@U).toarray()-np.eye(495)),1e-10)
            self.assertGreater(U.nnz,495)
            self.assertEqual(sum(len(inds) for inds,_,_ in local_hamiltonian_blocks(site,20,0)),495)

    def test_shared_Fock_Hubbard_reduces_to_exact_two_orbital_gate(self):
        # For fixed two occupied spectator orbitals, compare same 9x9
        # computational projection against the exact 15-state PR55 model.
        opt=optimize_pulse(20.)
        T=opt['time'];delta=opt['active_detuning']
        indices=[fock_index(a,0,0,m) for a in range(3) for m in range(3)]
        block=local_fock_pulse('A0',20.,delta,T)[np.ix_(indices,indices)].toarray()
        reference,leak=projected_evolution(20.,delta,T)
        self.assertLess(np.linalg.norm(block-reference),1e-10)
        self.assertGreater(leak,0)

    def test_charge_transport_second_quantization_preserves_all_sectors(self):
        rng=np.random.default_rng(18)
        z=rng.normal(size=(3,6,3,2,495,3,1))+1j*rng.normal(size=(3,6,3,2,495,3,1))
        self.assertLess(np.linalg.norm(link_transport(link_transport(z))-z),1e-10)
        for g in range(G):
            idx,sign=subgroup_wedge_map(g,(3,))
            self.assertEqual(len(set(idx.tolist())),495)
            self.assertLess(max(abs(sign))-1,1e-12)
        for vertex in ('A','B'):
            for h in range(G):
                self.assertLess(np.linalg.norm(
                    gauge_transform(link_transport(z),h,vertex)
                    -link_transport(gauge_transform(z,h,vertex))),1e-10)

    def test_four_fermion_initial_state_fulfills_extended_Gauss(self):
        x=initialize(self.K)
        self.assertEqual(x.shape,(3,6,3,2,495,3,4))
        self.assertLess(abs(np.vdot(x,x)-4),1e-10)
        self.assertLess(abs(occupation_leakage(x)),1e-12)
        for v in gauss_acceptance(x).values():self.assertAlmostEqual(v,1.,places=10)

    def test_microscopic_Hubbard_pulse_is_Gauss_covariant_offsector(self):
        rng=np.random.default_rng(7)
        x=rng.normal(size=(3,6,3,2,495,3,1))+1j*rng.normal(size=(3,6,3,2,495,3,1))
        z=apply_microscopic_gate(x,'A0',(1,2),20.,10.,25.)
        for vertex in ('A','B'):
            for h in range(G):
                self.assertLess(np.linalg.norm(gauge_transform(z,h,vertex)-
                  apply_microscopic_gate(gauge_transform(x,h,vertex),'A0',(1,2),20.,10.,25.)),1e-10)

    def test_leakage_survives_transport_and_next_exchange(self):
        opt=optimize_pulse(20.)
        z=initialize(self.K)
        z1=apply_microscopic_gate(z,'A0',(1,2),20.,opt['active_detuning'],opt['time'])
        first=occupation_leakage(z1)
        self.assertGreater(first,.001)
        z2=link_transport(z1)
        self.assertAlmostEqual(occupation_leakage(z2),first,places=11)
        z3=apply_microscopic_gate(z2,'B0',(1,2),20.,opt['active_detuning'],opt['time'])
        self.assertGreater(occupation_leakage(z3),first)
        for ztest in (z1,z2,z3):
            self.assertGreater(min(gauss_acceptance(ztest).values()),1-1e-10)

    def test_complete_persistent_fock_logical_syndrome_is_not_single_gate(self):
        opt=optimize_pulse(20.)
        data=exact_full_channel(self.K,20.,opt['active_detuning'],opt['time'])
        self.assertAlmostEqual(data['total_probability'],1.,places=10)
        self.assertGreater(data['mean_final_occupation_leakage'],.001)
        self.assertLess(data['unconditional_F_e_ideal_feedback'],
                        opt['optimized_metrics']['entanglement_fidelity'])
        self.assertTrue(0<=data['unconditional_F_e_ideal_feedback']<1)
        self.assertAlmostEqual(sum(data['pointer_probabilities']),1.,places=10)

    def test_full_report_has_three_interaction_regimes_and_scope(self):
        data=run()
        self.assertEqual(len(data['arms']),3)
        self.assertEqual(data['active_fermionic_Fock_dimension'],495)
        self.assertTrue(data['all_four_active_orbital_occupations_persist_across_16_operations'])
        for arm in data['arms']:
            self.assertLess(arm['illustrative_no_loss_probability'],1.)
            self.assertLess(arm['illustrative_loss_weighted_syndrome_F_e'],
                            arm['complete_syndrome_ideal_feedback_reference']['unconditional_F_e_ideal_feedback'])

if __name__=='__main__':unittest.main()