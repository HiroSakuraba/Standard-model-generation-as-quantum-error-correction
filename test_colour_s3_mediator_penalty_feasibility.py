"""Regression: mediator charge, exact finite-Rabi gadget and local gauge penalty."""
import unittest
import numpy as np

from colour_gauge_link_s3 import (
    DIM,G,physical_isometry,transform,projector,dressed_exchange
)
from colour_s3_edge_logical_syndrome import code_isometry
from colour_s3_two_register_syndrome import prepared_pointer
from colour_s3_mediator_penalty_feasibility import (
    charged_vertex_transform,charged_mediator_onsite_swap,
    charged_mediator_demo,gadget_amplitudes,encoded_two_level_gadget,
    gadget_recovery,penalty_single_pair_probe
)


class MediatorPenaltyFeasibilityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.W=physical_isometry()
        cls.K,cls.F,_=code_isometry(cls.W)

    def test_charged_mediator_local_swap_commutes_both_vertex_gauges(self):
        rng=np.random.default_rng(13)
        y=rng.normal(size=(DIM,3,2))+1j*rng.normal(size=(DIM,3,2))
        self.assertLess(np.linalg.norm(
            charged_mediator_onsite_swap(charged_mediator_onsite_swap(y))-y),1e-11)
        for vertex in ("A","B"):
            for h in range(G):
                self.assertLess(np.linalg.norm(
                    charged_vertex_transform(charged_mediator_onsite_swap(y),h,vertex)
                    -charged_mediator_onsite_swap(
                        charged_vertex_transform(y,h,vertex))),1e-10)

    def test_charged_mediator_only_solves_onsite_kinematics(self):
        r=charged_mediator_demo()
        self.assertEqual(r["mediator_register_dimension"],3)
        self.assertFalse(r["mediator_can_directly_read_logical_dressed_swap"])

    def test_rabi_two_level_channel_complete_and_longer_at_high_detuning(self):
        times=[]
        for D in (4.,8.,16.,32.):
            for l in (0,1,2):
                a,b,t=gadget_amplitudes(l,D)
                self.assertAlmostEqual(abs(a)**2+abs(b)**2,1.,places=11)
            times.append(t)
        self.assertTrue(all(times[j+1]>times[j] for j in range(3)))
        a,b,_=gadget_amplitudes(0,4.)
        self.assertLess(abs(a-1),1e-11)
        self.assertLess(abs(b),1e-11)

    def test_gadget_two_kraus_branches_gauge_covariant(self):
        rng=np.random.default_rng(8)
        x=rng.normal(size=(3,DIM,2))+1j*rng.normal(size=(3,DIM,2))
        for kind in ("edge","onsite"):
            a,b=encoded_two_level_gadget(x,kind,(1,2),8.)
            self.assertLess(abs(np.vdot(a,a)+np.vdot(b,b)-np.vdot(x,x)),2e-10)
            for v in ("A","B"):
                for h in range(G):
                    rot=np.stack([transform(y,h,v) for y in x])
                    aa,bb=encoded_two_level_gadget(rot,kind,(1,2),8.)
                    self.assertLess(np.linalg.norm(aa-np.stack([
                        transform(y,h,v) for y in a])),2e-10)
                    self.assertLess(np.linalg.norm(bb-np.stack([
                        transform(y,h,v) for y in b])),2e-10)

    def test_full_gadget_never_postselects_mediator_excitation(self):
        a=gadget_recovery(self.K,4.)
        b=gadget_recovery(self.K,32.)
        self.assertEqual(a["reset_environment_kraus_histories"],16)
        self.assertEqual(b["reset_environment_kraus_histories"],16)
        self.assertTrue(0<=a["unconditional_ideal_feedback_F_e"]<=1)
        self.assertTrue(0<=b["unconditional_ideal_feedback_F_e"]<=1)
        self.assertGreater(b["total_syndrome_time_four_windows"],
                           a["total_syndrome_time_four_windows"])
        self.assertGreater(b["unconditional_ideal_feedback_F_e"],
                           a["unconditional_ideal_feedback_F_e"])
        self.assertFalse(b["native_2body_gate_synthesis"])

    def test_one_pair_penalty_limit_has_nonzero_cost(self):
        rows=penalty_single_pair_probe(self.W,strengths=(0.,4.,16.))["results"]
        self.assertLess(rows[0]["single_pair_target_Fro_error_per_sqrt_14"],2e-9)
        for r in rows:
            for p in r["vertex_Gauss_accept"].values():
                self.assertGreaterEqual(p,-1e-10)
                self.assertLessEqual(p,1+1e-10)
            self.assertTrue(0<=r["joint_Gauss_return_probability"]<=1+1e-10)

    def test_physical_projector_is_preserved_by_vertex_penalty(self):
        x=np.random.default_rng(90).normal(size=(DIM,2))
        y=projector(projector(x,"A"),"B")
        self.assertLess(np.linalg.norm(self.W@self.W.conj().T@x-y),2e-10)


if __name__=="__main__":
    unittest.main()
