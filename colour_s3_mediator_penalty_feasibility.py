"""Three controlled-S3-gauge interaction routes: explicit feasibility bounds.

Route 1: an SU3-subgroup *charged mediator* qutrit at vertex A.
Its exchange with co-located matter is gauge covariant, but does NOT
implement a controlled dressed exchange or logical syndrome. It shows
what gauge charge the mediator must carry and which local operation is
actually allowed, without assuming a costless neutral pointer.

Route 2: perturbative two-level mediator gadget.
H=Delta |1><1|_m + g X_m (P_pointer + Q_dressed)
  + (g*g/Delta)(P_pointer + Q_dressed).
The exact two-level Rabi solution gives a controlled sign on the
pointer-and-dressed-antisymmetric sector in the Delta/g->infty limit.
Q_dressed=(I-S_dressed)/2 is already THREE-body. The mediator coupling
X_m Q_dressed is FOUR-register (mediator + edge), so the gadget does
NOT reduce the native many-body support. Finite detuning creates
ancilla excitation, imperfect phase, and runtime O(Delta/g^2).
All outcomes are traced without postselection.

Route 3: a gauge-penalty-assisted *single* pair factor from PR49.
Exactly propagate exp[-it(H_pair + lambda[(I-P_A)+(I-P_B)])].
Large penalty suppresses unphysical gauge transitions, but may also
distort the intended pair gate. This is ONE constituent gate, not a
complete hardware-synthesized 16-factor recovery circuit.

All three use the finite S3 two-vertex / 4D phase code. No native
transmon, continuous SU3, noisy physical syndrome, or nine-quark
21D/65-syndrome claim.
"""
from __future__ import annotations

from pathlib import Path
import json
import numpy as np
from scipy.linalg import logm, expm
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import expm_multiply

from colour_gauge_link_s3 import (
    DIM,G,MATTER,REP,physical_isometry,transform,projector,
    dressed_exchange,
)
from colour_s3_edge_logical_syndrome import code_isometry, onsite_swap
from colour_s3_two_register_syndrome import (
    exchange_eigenbasis_q, pair_gate, prepared_pointer, phase_readout,
)

LOGICAL=4


def charged_vertex_transform(x,h,vertex):
    """S3 gauge action on (link,matter4,charged-mediator-A,reference).

    Mediator is the fundamental signed-permutation D(h) at A; B acts
    only on original physical registers. Charge is transported only
    in models that explicitly supply the necessary gauge link.
    """
    y=np.asarray(x,complex)
    assert y.ndim==3 and y.shape[:2]==(DIM,3)
    local=np.stack([transform(y[:,m,:],h,vertex) for m in range(3)],axis=1)
    if vertex=="A":
        local=np.einsum("ab,ibk->iak",REP[h],local)
    return local


def charged_mediator_onsite_swap(x):
    """Two qutrits: exchange charged mediator at A with colour qutrit A0."""
    y=np.asarray(x,complex)
    assert y.ndim==3 and y.shape[:2]==(DIM,3)
    # 6 link, A0,A1,B0,B1, mediator, reference.
    arr=y.reshape((G,3,3,3,3,3,y.shape[-1]))
    return np.swapaxes(arr,1,5).reshape(y.shape)


def charged_mediator_demo():
    rng=np.random.default_rng(8)
    x=rng.normal(size=(DIM,3,2))+1j*rng.normal(size=(DIM,3,2))
    worst=0.
    for v in ("A","B"):
        for h in range(G):
            lhs=charged_vertex_transform(charged_mediator_onsite_swap(x),h,v)
            rhs=charged_mediator_onsite_swap(charged_vertex_transform(x,h,v))
            worst=max(worst,float(np.linalg.norm(lhs-rhs)))
    assert worst<1e-10
    assert np.linalg.norm(charged_mediator_onsite_swap(
        charged_mediator_onsite_swap(x))-x)<1e-10
    return {"mediator_register_dimension":3,
            "mediator_representation":"S3 signed-permutation fundamental at vertex A",
            "allowed_onsite_two_qutrit_exchange_covariance_residual":worst,
            "onsite_exchange_is_unitary_involution":True,
            "mediator_can_directly_read_logical_dressed_swap":False,
            "missing":"Charged mediator transfer across gauge link and gauge-invariant controlled logical readout are NOT synthesized"}


def gadget_amplitudes(level,detuning,coupling=1.):
    """Exact 2x2 Rabi amplitudes after duration t=pi*Delta/(2g²).

    Start ancilla ground. For eigenvalue l of L=P_pointer+Q_edge,
    H_l=[[g² l/Delta, g l],[g l,Delta+g² l/Delta]].
    Counterterm removes the undesired one-projector Stark shifts to
    second order, leaving -2g²/Delta P_pointer Q_edge.
    """
    if detuning<=0 or coupling<=0 or level not in (0,1,2):
        raise ValueError("parameters")
    t=np.pi*detuning/(2*coupling**2)
    w=np.hypot(detuning/2,coupling*level)
    scalar=coupling**2*level/detuning+detuning/2
    pref=np.exp(-1j*t*scalar)
    a=pref*(np.cos(t*w)+1j*(detuning/2)/w*np.sin(t*w))
    b=pref*(-1j*coupling*level/w*np.sin(t*w))
    assert abs(abs(a)**2+abs(b)**2-1)<2e-12
    return a,b,t


def encoded_two_level_gadget(z,kind,active,detuning):
    """Ancilla reset and trace: two CPTP Kraus branches, no postselection.

    All branches act through the exact dressed H_edge or onsite H_A
    and a pointer projector. Each branch commutes with vertex Gauss
    actions because its polynomial basis does.
    """
    y=np.asarray(z,complex)
    assert y.ndim==3 and y.shape[0:2]==(3,DIM)
    op=dressed_exchange if kind=="edge" else onsite_swap if kind=="onsite" else None
    if op is None:raise ValueError(kind)
    active=set(active)
    if not active.issubset({0,1,2}):raise ValueError(active)
    out0=np.empty_like(y)
    out1=np.empty_like(y)
    for p in range(3):
        q=(y[p]-op(y[p]))/2
        rest=y[p]-q
        bias=1 if p in active else 0
        a0,b0,_=gadget_amplitudes(bias,detuning)
        a1,b1,_=gadget_amplitudes(bias+1,detuning)
        out0[p]=a0*rest+a1*q
        out1[p]=b0*rest+b1*q
    assert abs(np.vdot(out0,out0).real+np.vdot(out1,out1).real
               -np.vdot(y,y).real)<5e-10
    return out0,out1


SCHEDULE=(("edge",(1,2)),("onsite",(1,2)),
          ("edge",(2,)),("onsite",(2,)))


def gadget_recovery(K,detuning,p_error=.37):
    """Full 2^4 mediator environment, 3 readout labels, two source Kraus.

    Exactly 16 reset-medicator Kraus histories, no trajectory sampling.
    Compute unconditional entanglement fidelity after *ideal* known
    dressed swap conditional correction. Reports true mediator reset
    failure probability and source-noise-independent runtime.
    """
    if not 0<=p_error<=1:raise ValueError(p_error)
    time=np.pi*detuning/2
    samples=[]
    max_excitations=[]
    fidelity=0.
    envtotal=0.
    for src in (np.sqrt(1-p_error)*K,np.sqrt(p_error)*dressed_exchange(K)):
        histories=[prepared_pointer(src)]
        for kind,active in SCHEDULE:
            updated=[]
            excited=0.
            for z in histories:
                a,b=encoded_two_level_gadget(z,kind,active,detuning)
                excited+=np.vdot(b,b).real/LOGICAL
                updated.extend((a,b))
            histories=updated
            max_excitations.append(float(excited))
        for z in histories:
            outcomes=phase_readout(z)
            for label,branch in enumerate(outcomes):
                corrected=dressed_exchange(branch) if label==2 else branch
                fidelity+=abs(np.trace(K.conj().T@corrected))**2/LOGICAL**2
                envtotal+=np.vdot(branch,branch).real/LOGICAL
    assert abs(envtotal-1)<5e-10
    assert 0<=fidelity<=1+1e-10
    return {"detuning_over_g":float(detuning),
            "coupling_g_arbitrary_units":1.,
            "time_per_syndrome_window":float(time),
            "total_syndrome_time_four_windows":float(4*time),
            "reset_environment_kraus_histories":16,
            "unconditional_ideal_feedback_F_e":float(fidelity),
            "maximum_unconditional_mediator_excitation_per_window":
                float(max(max_excitations)),
            "unconditional_mediator_excitation_per_window":max_excitations,
            "three_body_edge_Q_preexisting":True,
            "mediator_edge_interaction_register_support":4,
            "native_2body_gate_synthesis":False}


def penalty_single_pair_probe(W,strengths=(0.,2.,8.,32.),time=1.):
    """One PR49 Q_A0,B0 factor driven with local vertex Gauss penalty.

    H_pair=i log(Q_A0,B0); P_A, P_B are group-twirled local projectors.
    A lambda*(2I-P_A-P_B) penalty is (link + A matter2) plus
    (link + B matter2), each a THREE-register Hamiltonian.
    Single gate, not a complete syndrome sequence or native gate.
    """
    Q=exchange_eigenbasis_q()
    h=1j*logm(Q)
    h=(h+h.conj().T)/2
    assert np.linalg.norm(expm(-1j*h)-Q)<1e-9
    H=pair_gate(np.eye(DIM,dtype=complex),(2,4),h)
    H=(H+H.conj().T)/2
    pa=projector(np.eye(DIM,dtype=complex),"A")
    pb=projector(np.eye(DIM,dtype=complex),"B")
    assert np.linalg.norm(pa@pa-pa)<2e-10
    assert np.linalg.norm(pb@pb-pb)<2e-10
    penalty=2*np.eye(DIM)-pa-pb
    penalty=(penalty+penalty.conj().T)/2
    assert np.linalg.norm(penalty@W)<2e-10
    target=pair_gate(W,(2,4),Q)
    out=[]
    for lam in strengths:
        pulse=expm_multiply((-1j*time)*csr_matrix(H+lam*penalty),W)
        acceptance={v:float(np.linalg.norm(projector(pulse,v))**2/14)
                    for v in ("A","B")}
        target_error=float(np.linalg.norm(pulse-target)/np.sqrt(14))
        physical_return=float(np.linalg.norm(W.conj().T@pulse)**2/14)
        out.append({"penalty_lambda_relative_to_pair_drive":float(lam),
                    "vertex_Gauss_accept":acceptance,
                    "joint_Gauss_return_probability":physical_return,
                    "single_pair_target_Fro_error_per_sqrt_14":target_error})
    assert max(abs(x["single_pair_target_Fro_error_per_sqrt_14"])
               for x in out[:1])<2e-9
    return {"probe":"one general U9 pair factor Q_A0,B0 only, not completed PR49 syndrome",
            "local_gauss_penalty_support":"two sums of vertex-local THREE-register projectors",
            "penalty_natively_synthesized":False,
            "pulse_time":float(time),"results":out}


def run():
    W=physical_isometry()
    K,F,_=code_isometry(W)
    charge=charged_mediator_demo()
    gadget=[gadget_recovery(K,x) for x in (4.,8.,16.,32.)]
    penalty=penalty_single_pair_probe(W)
    assert gadget[-1]["unconditional_ideal_feedback_F_e"] > gadget[0]["unconditional_ideal_feedback_F_e"]
    assert all(0<=row["unconditional_ideal_feedback_F_e"]<=1+1e-10 for row in gadget)
    return {
        "code":"finite S3, 486-dimensional two-vertex physical register, 4D phase code",
        "charged_mediator":charge,
        "perturbative_mediator_gadget":gadget,
        "gauge_penalty_single_pair":penalty,
        "caveats":[
            "Charged mediator provides only a local covariant swap, not a complete logical readout",
            "Gadget coupling X_m*Q_edge itself requires FOUR-register mediator+edge interaction; no native 2-body reduction",
            "Gadget times are arbitrary units; increasing detuning decreases excitation but lengthens gates",
            "One-factor penalty probe cannot establish a complete gauge-penalized syndrome recovery",
            "Penalty projectors are formal 3-register vertex operators, not a calibrated physical Gauss penalty",
            "No environmental noise model, native pulse, T1/T2, link leakage, routing or device comparison",
            "Not continuous SU3 or original 21D nine-quark 65-syndrome problem"]}


def main():
    out=run()
    Path("colour_s3_mediator_penalty_results.json").write_text(
        json.dumps(out,indent=2,sort_keys=True)+"\n")
    print("S3_MEDIATOR_FEASIBILITY",json.dumps({
        "gadget": [{k:v for k,v in row.items() if
                    k!="unconditional_mediator_excitation_per_window"}
                   for row in out["perturbative_mediator_gadget"]],
        "charged_mediator":out["charged_mediator"],
        "penalty":out["gauge_penalty_single_pair"]
    },sort_keys=True))


if __name__=="__main__":main()
