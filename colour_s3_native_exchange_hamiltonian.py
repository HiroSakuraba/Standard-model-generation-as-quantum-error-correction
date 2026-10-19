"""SU(3)-isotropic quantum-controlled exchange Hamiltonian on charged bus.

We specify, rather than assume a Fredkin gate:
 T^a=lambda^a/2 (eight SU3 generators),
 C=sum_a T^a_M tensor T^a_q=(SWAP-I/3)/2,
 H_exchange=I/3-C=(I-SWAP)/2,
 H_control = (sum_{p in active}|p><p|) tensor H_exchange.

For an onsite charged mediator and matter qutrit at one gauge vertex,
H_control commutes with the FULL diagonal SU(3) action (hence S3)
at every pulse time, and e^{-i pi H_control} is an exact controlled
qutrit SWAP. The mobile mediator transport remains the separately
specified link+mediator involution from PR53.

SU3 obstruction: because a single fundamental is irreducible, a
gauge-invariant pointer-fundamental TWO-register coupling is
A_pointer tensor I_colour (Schur). A gauge-invariant colour-colour
2-register exchange cannot depend on the neutral pointer. Thus
their Lie closure cannot generate H_control. In this sense the
three-register, quantum-conditional interaction is irreducible
under the STRICT full-SU3, neutral-pointer, 2-register alphabet.
This is NOT a no-go for the REDUCIBLE finite S3 fundamental, charged
pointer/extra ancillas, effective many-body gadgets, or non-Gauss
intermediate control.

Implementation: exact 9x9 matrix exponential of H_exchange and of
a specified colour-anisotropic perturbation epsilon*T_M^3*T_q^3,
embedded into the complete 3x486x2x3x3xR bus tensor. The isotropic
full logical syndrome must exactly reproduce PR53; anisotropy can
break local Gauss covariance and degrade recovery. Counts and times
are scenario assumptions, not hardware measurements.
"""
from __future__ import annotations

from functools import lru_cache
from pathlib import Path
import json

import numpy as np
from scipy.linalg import expm
from colour_singlet_code import gellmann
from colour_gauge_link_s3 import (
    G, DIM, REP, physical_isometry, dressed_exchange,
)
from colour_s3_edge_logical_syndrome import code_isometry
from colour_s3_charged_bus_transport import (
    initialize, gauge_action, local_swap, involution_pulse,
    transport, gauss_acceptance, remove_ancilla_singlet,
    exact_unconditional_fidelity, SCHEDULE,
)

SWAP9 = np.zeros((9,9),complex)
for a in range(3):
    for b in range(3):
        SWAP9[3*b+a,3*a+b]=1
I9=np.eye(9,dtype=complex)
GENS=gellmann()
C_SU3=sum((np.kron(T,T) for T in GENS), np.zeros((9,9),complex))
H_EX=I9/3-C_SU3
# Off-axis SU3-exchange anisotropy, selected physically unsupported term.
H_ANISO=np.kron(GENS[6],GENS[6])


def verify_exchange_identity():
    residual=float(np.linalg.norm(H_EX-(I9-SWAP9)/2))
    assert residual<2e-14
    for T in GENS:
        total=np.kron(T,np.eye(3))+np.kron(np.eye(3),T)
        assert np.linalg.norm(H_EX@total-total@H_EX)<1e-12
    assert np.linalg.norm(H_EX@H_EX-H_EX)<2e-14
    return residual


def numerical_commutant_dim(representations):
    """Independent 9x9 full-matrix commutator nullity crosscheck.

    This floating-point rank check is NOT the proof; Schur–Weyl
    irreducibility provides the exact representation-theoretic result.
    """
    basis=np.eye(81,dtype=complex).reshape((81,9,9))
    rows=[]
    for generator in representations:
        diffs=np.einsum("kij,jl->kil",basis,generator)-np.einsum(
            "ij,kjl->kil",generator,basis)
        rows.append(diffs.reshape((81,81)).T)
    singular=np.linalg.svd(np.concatenate(rows),compute_uv=False)
    return int(np.sum(singular<1e-9))


def su3_commutant_neutral_pointer():
    """Construct a basis of SU3-invariant pointer⊗fundamental operators.

    The pointer is gauge-neutral dimension 3. The only *single-quark*
    SU3 invariants are scalars, so all two-register pointer-quark
    invariants have form A_3⊗I_3 (dimension 9 as a complex vector space).
    This identity is independent of a 9x9 numerical null-space search.
    """
    ops=[]
    for i in range(3):
        for j in range(3):
            E=np.zeros((3,3),complex)
            E[i,j]=1
            M=np.kron(E,np.eye(3))
            ops.append(M)
            for T in GENS:
                GT=np.kron(np.eye(3),T)
                assert np.linalg.norm(M@GT-GT@M)<1e-12
    # Independently measure both complex commutants in dimension nine.
    pointer_dim=numerical_commutant_dim([np.kron(np.eye(3),T) for T in GENS])
    matter_pair_dim=numerical_commutant_dim([np.kron(T,np.eye(3))+np.kron(np.eye(3),T) for T in GENS])
    assert pointer_dim==9 and matter_pair_dim==2
    # A gauge-invariant two-quark coupling is span I,SWAP.
    test=np.kron(GENS[6],np.eye(3))
    assert np.linalg.norm(SWAP9@test-test@SWAP9)>0.1
    return {"neutral_pointer_fundamental_invariant_complex_dimension":pointer_dim,
            "basis_form":"A_pointer tensor I_3",
            "fundamental_quark_pair_invariant_complex_dimension":matter_pair_dim,
            "pair_form":"a I_9 + b SWAP",
            "generated_neutral_pointer_to_quark_entanglement":False,
            "strict_scope":"irreducible SU3 fundamental, neutral pointer, SU3-invariant 1- and 2-register Hamiltonians only",
            "does_not_apply":"S3 restricted reducible representation, charged ancillas, 3-body terms, gauge-noninvariant intermediate operations"}


def controlled_hamiltonian(active=(1,2),anisotropy=0.):
    """Exact 27x27 pointer + two fundamentals H, one vertex."""
    if set(active)-{0,1,2} or not np.isfinite(anisotropy):
        raise ValueError("active/pulse anisotropy")
    P=np.diag([float(p in active) for p in range(3)])
    return np.kron(P,H_EX+anisotropy*H_ANISO)


@lru_cache(maxsize=64)
def local_two_qutrit_unitary(theta,anisotropy=0.):
    """Controlled-active block gate e^-i(theta*(H_EX+eps H_ANISO))."""
    if not np.isfinite(theta) or not np.isfinite(anisotropy):
        raise ValueError("finite pulse parameters required")
    return expm(-1j*theta*(H_EX+anisotropy*H_ANISO))


def apply_local_hamiltonian(z,site,active=(1,2),theta=np.pi,anisotropy=0.):
    """Embed specified 9x9 exchange H into the full charged-bus register.

    Local pair: the mobile mediator's colour and ONE matter colour at
    the mediator's present vertex. Pointer controls the pulse. If the
    mediator is at the wrong vertex or pointer inactive, apply identity.
    """
    x=np.asarray(z,complex)
    if x.ndim!=6 or x.shape[:5]!=(3,DIM,2,3,3):
        raise ValueError("shape of full charged register")
    if site not in ("A0","A1","B0","B1"):
        raise ValueError(site)
    if set(active)-{0,1,2}:
        raise ValueError(active)
    place=0 if site[0]=="A" else 1
    k={"A0":1,"A1":2,"B0":3,"B1":4}[site]
    U=local_two_qutrit_unitary(float(theta),float(anisotropy))
    dest=x.copy()
    arr=dest.reshape((3,G,3,3,3,3,2,3,3,x.shape[-1]))
    for pointer in set(active):
        block=arr[pointer,:,:,:,:,:,place,:,:,:]
        # axes now: link 0, A0 1, A1 2, B0 3, B1 4,
        # mediator colour 5, stationary reference 6, logical R 7.
        front=np.moveaxis(block,(k,5),(0,1))
        shape=front.shape
        evolved=(U@front.reshape(9,-1)).reshape(shape)
        arr[pointer,:,:,:,:,:,place,:,:,:]=np.moveaxis(evolved,(0,1),(k,5))
    return dest


def primitive_schedule():
    operations=[]
    for kind,active in SCHEDULE:
        if kind=="edge":
            operations.extend((
                ("A0",active),
                ("T",None),
                ("B0",active),
                ("T",None),
                ("A0",active)))
        else:
            operations.extend((("A0",active),("A1",active),("A0",active)))
    assert len(operations)==16
    return operations


def syndrome_hamiltonian_circuit(K,anisotropy=0.,transport_overrotation=0.,
                                  return_snapshots=False):
    """Full charged-bus qutrit Fourier syndrome from specified local H.

    The link T is the exact gauge-covariant two-register transport
    assumed in PR53. Every onsite controlled-SWAP is now realized
    by e^{-i pi P_pointer H_exchange} on the embedded 27D subsystem.
    """
    x=initialize(K)
    snapshots=[]
    for name,active in primitive_schedule():
        if name=="T":
            x=involution_pulse(x,transport,np.pi*(1+transport_overrotation))
        else:
            x=apply_local_hamiltonian(x,name,active,np.pi,anisotropy)
        if return_snapshots:snapshots.append((name,x.copy()))
    x=np.fft.fft(x,axis=0)/np.sqrt(3)
    return (x,snapshots) if return_snapshots else x


def feedback_hamiltonian_circuit(z,anisotropy=0.):
    """Label 2 selected after readout: 3 matter exchange + 2 transport."""
    x=z
    for name in ("A0","T","B0","T","A0"):
        if name=="T":
            x=transport(x)
        else:
            x=apply_local_hamiltonian(x,name,(2,),np.pi,anisotropy)
    return x


def fidelity(K,anisotropy=0.,p=.37):
    """Exact CPTP logical channel; trace every pointer/charged environment."""
    if not 0<=p<=1:raise ValueError(p)
    total=0.; fe=0.
    pointer=np.zeros(3)
    for src in (np.sqrt(1-p)*K,np.sqrt(p)*dressed_exchange(K)):
        z=feedback_hamiltonian_circuit(
            syndrome_hamiltonian_circuit(src,anisotropy),anisotropy)
        pointer+=np.array([np.vdot(b,b).real/4 for b in z])
        for label in range(3):
            for loc in range(2):
                for m in range(3):
                    for ref in range(3):
                        b=z[label,:,loc,m,ref,:]
                        total+=np.vdot(b,b).real/4
                        fe+=abs(np.trace(K.conj().T@b))**2/16
    assert abs(total-1)<3e-10
    assert -1e-8<=fe<=1+1e-8
    return {"anisotropic_colour_exchange_coefficient_relative_to_isotropic":float(anisotropy),
            "source_known_dressed_swap_probability":p,
            "unconditional_entanglement_fidelity":float(fe),
            "total_output_probability":float(total),
            "pointer_probabilities":pointer.tolist()}


def full_SU3_pulse_covariance_witness(epsilon):
    """Onsite generator SU3 and finite subgroup commutator witnesses."""
    H=H_EX+epsilon*H_ANISO
    su3=max(float(np.linalg.norm(H@(
        np.kron(T,np.eye(3))+np.kron(np.eye(3),T))-
        (np.kron(T,np.eye(3))+np.kron(np.eye(3),T))@H))
        for T in GENS)
    s3=max(float(np.linalg.norm(H@np.kron(D,D)-np.kron(D,D)@H))
           for D in REP)
    if abs(epsilon)<1e-14:
        assert su3<1e-12 and s3<1e-12
    else:
        assert su3>0 and s3>0
    return {"SU3_generators_max_commutator_fro":su3,
            "S3_group_max_commutator_fro":s3}


def timing_scenario(J_mhz=5.,transport_ns=100.):
    """Illustrative possible *effective* three-body coupling strength.

    H/h = J_MHz * P_pointer (I-SWAP)/2:
    theta=2pi J t, so pi-pulse time is 1/(2 J) us.
    No transmon or laboratory J for this three-body Hamiltonian exists.
    """
    if J_mhz<=0 or transport_ns<0:raise ValueError("positive coupling/time")
    control_ns=500./J_mhz
    return {"hypothetical_control_H_over_h_MHz":J_mhz,
            "assumed_effective_control_pi_pulse_ns":control_ns,
            "assumed_charged_transport_per_operation_ns":transport_ns,
            "twelve_controlled_local_swaps":12,
            "four_link_transports":4,
            "syndrome_excluding_pointer_and_reset_ns":
                 12*control_ns+4*transport_ns,
            "feedback_excluding_measurement_ns":
                 3*control_ns+2*transport_ns,
            "hardware_calibration":False}


def study():
    residual=verify_exchange_identity()
    no_go=su3_commutant_neutral_pointer()
    W=physical_isometry()
    K,F,_=code_isometry(W)
    raw=np.random.default_rng(810).normal(size=(3,DIM,2,3,3,1))
    raw=raw+1j*np.random.default_rng(811).normal(size=raw.shape)
    err=float(np.linalg.norm(apply_local_hamiltonian(
        raw,"A0",(1,2),np.pi)-local_swap(raw,"A0",(1,2))))
    assert err<2e-11
    identity=float(np.linalg.norm(apply_local_hamiltonian(
        raw,"B0",(),np.pi)-raw))
    assert identity<1e-12
    output,steps=syndrome_hamiltonian_circuit(K,return_snapshots=True)
    old_reference=remove_ancilla_singlet(output)
    from colour_s3_edge_logical_syndrome import pointer_circuit
    match=float(np.linalg.norm(old_reference-pointer_circuit(K)))
    assert match<1e-10
    stage_cov={v:1. for v in ("A","B")}
    for _,z in steps:
        gp=gauss_acceptance(z)
        for v in stage_cov:stage_cov[v]=min(stage_cov[v],gp[v])
    assert min(stage_cov.values())>1-1e-10
    witness=[full_SU3_pulse_covariance_witness(eps)
             for eps in (0.,.01,.05)]
    results=[fidelity(K,eps) for eps in (0.,.01,.05)]
    assert abs(results[0]["unconditional_entanglement_fidelity"]-1)<1e-10
    assert results[-1]["unconditional_entanglement_fidelity"]<1
    return {"model":"finite S3 charged-bus simulator, exact SU3 isotropic onsite conditional exchange Hamiltonian",
        "exchange_identity_residual":residual,
        "SU3_two_register_neutral_pointer_no_go":no_go,
        "local_pi_pulse_vs_Fredkin_residual":err,
        "inactive_pointer_branch_identity_residual":identity,
        "full_ideal_syndrome_reference_residual":match,
        "minimum_vertex_Gauss_accept_after_each_stage":stage_cov,
        "anisotropy_covariance_witnesses":witness,
        "full_logical_channel_anisotropy_sweep":results,
        "hypothetical_timing":timing_scenario(),
        "Hamiltonian_resources":{
            "coloured_exchange_generators":8,
            "control_generator":"P_pointer tensor [I/3 - sum_a T_m^a tensor T_site^a]",
            "max_support":3,
            "coherent_native_three_body_coupling_required":True,
            "Gauss_neutral_pointer_two_body_synthesis_without_new_resources":False,
            "requires_unsynthesized_two_body_link_transport":True
        },
        "limitations":[
            "This is a specified three-register Hamiltonian, not a realization of the three-body coupling from native two-body gates",
            "Pointer-dependent SU3 Heisenberg exchange amplitude is an engineering assumption, not measured hardware control",
            "Two-body obstruction is strict SU3 irreducible fundamental; S3 signed-permutation fundamental is reducible and requires separate reachability analysis",
            "Full system's link transport remains abstract two-register S3-controlled D(g) operation, not native decomposed",
            "Time uses hypothetical J_MHz and hypothetical transport duration, not a device calibration",
            "Anisotropy is synthetic coherent SU3-breaking colour interaction, not T1/T2 or four-level transmon leakage",
            "Finite S3 four-dimensional code, NOT continuous SU3 21D nine-quark or 65-syndrome recovery"
        ]}


def main():
    result=study()
    Path("colour_s3_native_exchange_hamiltonian_results.json").write_text(
        json.dumps(result,indent=2,sort_keys=True)+"\n")
    print("LOCAL_EXCHANGE_HAMILTONIAN",json.dumps({
        k:v for k,v in result.items()
        if k not in ("full_logical_channel_anisotropy_sweep","limitations")
    },sort_keys=True),flush=True)
    for r in result["full_logical_channel_anisotropy_sweep"]:
        print("EXCHANGE_ANISOTROPY_FIDELITY",json.dumps(r,sort_keys=True),flush=True)

if __name__=="__main__":
    main()
