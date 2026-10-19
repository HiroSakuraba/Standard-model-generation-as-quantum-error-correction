"""Regression suite for charged mobile mediator and gauge-covariant syndrome.

The mediator is paired with a stationary antifundamental reference in
a physical singlet: testing the charged bare mediator alone would violate
local Gauss-law preparation assumptions.
"""
import unittest
import numpy as np

from colour_gauge_link_s3 import (
    DIM,G,physical_isometry,transform,dressed_exchange
)
from colour_s3_edge_logical_syndrome import code_isometry,pointer_circuit
from colour_s3_charged_bus_transport import (
    initialize,gauge_action,transport,local_swap,involution_pulse,
    edge_bus,onsite_bus,run_syndrome,correction,
    remove_ancilla_singlet,exact_unconditional_fidelity,
    gauss_acceptance,schedule_ops,run,
)


class ChargedBusTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.W=physical_isometry()
        cls.K,cls.F,_=code_isometry(cls.W)

    def test_singlet_preparation_is_simultaneously_gauge_invariant(self):
        x=initialize(self.K)
        self.assertEqual(x.shape,(3,DIM,2,3,3,4))
        self.assertAlmostEqual(np.vdot(x,x).real,4.,places=10)
        for vertex in ("A","B"):
            for h in range(G):
                self.assertLess(np.linalg.norm(gauge_action(x,h,vertex)-x),1e-10)

    def test_transport_is_involutory_unitary_and_gauge_covariant(self):
        rng=np.random.default_rng(117)
        x=rng.normal(size=(3,DIM,2,3,3,1))+1j*rng.normal(
            size=(3,DIM,2,3,3,1))
        self.assertLess(np.linalg.norm(transport(transport(x))-x),1e-10)
        self.assertLess(abs(np.vdot(transport(x),transport(x))-np.vdot(x,x)),1e-10)
        for vertex in ("A","B"):
            for h in range(G):
                lhs=gauge_action(transport(x),h,vertex)
                rhs=transport(gauge_action(x,h,vertex))
                self.assertLess(np.linalg.norm(lhs-rhs),1e-10)

    def test_onsite_mobile_swap_respects_both_vertex_gauge_groups(self):
        rng=np.random.default_rng(19)
        x=rng.normal(size=(3,DIM,2,3,3,1))+1j*rng.normal(
            size=(3,DIM,2,3,3,1))
        for site in ("A0","A1","B0"):
            op=lambda z:local_swap(z,site,(1,2))
            self.assertLess(np.linalg.norm(op(op(x))-x),1e-10)
            for vertex in ("A","B"):
                for h in range(G):
                    self.assertLess(np.linalg.norm(gauge_action(op(x),h,vertex)
                        -op(gauge_action(x,h,vertex))),1e-10)

    def test_exact_link_dressed_swap_with_mediator_return(self):
        x=initialize(self.W[:, :2],pointer_uniform=False)
        y=edge_bus(x)
        self.assertLess(np.linalg.norm(remove_ancilla_singlet(y)[0]
                     -dressed_exchange(self.W[:,:2])),1e-10)
        self.assertLess(abs(np.vdot(y,y)-
                            np.vdot(remove_ancilla_singlet(y),
                                    remove_ancilla_singlet(y))),1e-10)
        self.assertLess(np.linalg.norm(edge_bus(y)-x),1e-10)

    def test_exact_pointer_syndrome_matches_original_code_in_full_space(self):
        output=run_syndrome(initialize(self.K))
        recovered=remove_ancilla_singlet(output)
        self.assertLess(np.linalg.norm(recovered-pointer_circuit(self.K)),1e-10)
        self.assertLess(abs(np.vdot(output,output)-np.vdot(recovered,recovered)),1e-10)
        self.assertLess(np.linalg.norm(recovered[1]-self.K),1e-10)
        self.assertLess(np.linalg.norm(recovered[0]),1e-10)
        errored=remove_ancilla_singlet(run_syndrome(initialize(self.F)))
        self.assertLess(np.linalg.norm(errored[2]-self.F),1e-10)

    def test_gauss_invariance_after_every_syndrome_operation(self):
        x=initialize(self.K)
        z,stages=run_syndrome(x,history=True)
        self.assertEqual(len(stages),16)
        self.assertEqual(sum(support==2 for _,support,_ in stages),4)
        for _,_,state in stages:
            for p in gauss_acceptance(state).values():
                self.assertAlmostEqual(p,1.,places=10)
        self.assertEqual(len(schedule_ops()),16)

    def test_full_channel_fidelity_no_ancilla_postselection(self):
        for p in (0.,.37,1.):
            result=exact_unconditional_fidelity(self.K,p,overrotation=0.)
            self.assertAlmostEqual(result['unconditional_entanglement_fidelity'],1.,places=10)
            self.assertAlmostEqual(result['total_output_probability'],1.,places=10)
        disturbed=exact_unconditional_fidelity(self.K,.37,overrotation=.04)
        self.assertGreaterEqual(disturbed['unconditional_entanglement_fidelity'],0.)
        self.assertLess(disturbed['unconditional_entanglement_fidelity'],1.)
        self.assertAlmostEqual(disturbed['total_output_probability'],1.,places=10)

    def test_full_report_explicit_scope_and_gauge_path(self):
        result=run()
        self.assertEqual(result['total_syndrome_elementary_count'],16)
        self.assertEqual(result['max_elementary_register_support'],3)
        self.assertEqual(result['two_register_link_transport_count'],4)
        self.assertEqual(result['three_register_pointer_conditioned_onsite_swaps'],12)
        self.assertLess(result['exact_high_level_PR48_syndrome_residual_fro'],1e-10)
        self.assertLess(result['maximum_elementary_gauge_covariance_residual_fro'],1e-10)


if __name__=="__main__":
    unittest.main()
