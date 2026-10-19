"""Gauge-admissible passive 4D memory vs active finite-S3 logical syndrome.

Both start in the SAME two-vertex 486D gauge Hilbert space, have
four logical states, and experience the SAME *specified gauge-dressed*
matter A0-B0 SWAP channel {sqrt(1-p) I,sqrt(p) S}. The passive code is
a canonical orthonormal four-dimensional subspace of the 10D +1
eigenspace of S restricted to the 14D Gauss invariant sector.
Therefore S C = C exactly and a passive decoder does not need a
syndrome for this ONE particular physically allowed error.

For loss/dephasing and time comparisons, both use the SAME four active
fermionic orbitals A0,A1,B0 and charged mediator M plus the same
spectator colour B1, gauge link and mobile/reference singlet at rest,
as PRs #53-57. Passive storage sits at one particle per active orbital
and idles for the ACTUAL PR57 syndrome+microscopic-feedback duration.
Both suffer loss of the same 4 active fermions with equal T1. We do
NOT silently compare a two-fermion non-Gauss bare colour-qubit memory.

During passive idling the exact colour-mode number Lindblad channel
has commuting diagonal operators and can be evaluated analytically,
including mediator colour correlations with the charged reference.
The code's dephasing comparator uses existing stochastic full 495D
PR57 dynamics, so estimates have Monte Carlo error.

A SECOND physically allowed dressed exchange (A1,B0) is tested to
show how error-tailored immunity can fail for unanticipated faults.
The active 4D phase code is NOT presumed to correct that second swap.

Not device calibrated: gauge link transport and charged reference
preparation remain ideal, both 4D encodings' preparation is assumed
available, and full SU3/nine-quark 65-syndrome recovery is open.
"""
from __future__ import annotations

from pathlib import Path
import json
import numpy as np

from colour_gauge_link_s3 import (
    DIM, G, MATTER, physical_isometry, dressed_exchange, transform,
)
from colour_s3_edge_logical_syndrome import (
    code_isometry, onsite_swap, pointer_circuit, edge_cycle,
)
from colour_s3_charged_bus_transport import initialize, gauss_acceptance
from colour_su3_hubbard_virtual_exchange import optimize_pulse
from colour_su3_hubbard_loss_dephase_feedback import (
    process_sample, time_unit_ns,
)
from colour_su3_duration_matched_bare import (
    physical_duration_ns, exactly_flagged_fermion_survival,
)

DLOG=4
P_SOURCE=.37


def canonical_swap_passive_code(W):
    """Deterministically pick four gauge-invariant S=+1 eigenstates.

    The full + eigenspace has dimension 10, in a 14D Gauss sector.
    Use projected coordinate vectors and deterministic Gram-Schmidt
    with fixed traversal (no eigenvalue-degeneracy dependence).
    """
    if W.shape!=(DIM,14):
        raise ValueError("Expected the 14D physical Gauss isometry")
    S=W.conj().T@dressed_exchange(W)
    assert np.linalg.norm(S-S.conj().T)<1e-10
    assert np.linalg.norm(S@S-np.eye(14))<1e-10
    P=(np.eye(14)+S)/2
    rank=int(round(np.trace(P).real))
    assert rank==10
    basis=[]
    for i in range(14):
        v=P[:,i].copy()
        for u in basis:
            v-=u*np.vdot(u,v)
        # Second pass strengthens numerical orthogonality.
        for u in basis:
            v-=u*np.vdot(u,v)
        norm=np.linalg.norm(v)
        if norm>1e-9:
            basis.append(v/norm)
        if len(basis)==DLOG:
            break
    if len(basis)!=DLOG:
        raise AssertionError("Physical + sector cannot fit ququart")
    C=W@np.stack(basis,axis=1)
    assert np.linalg.norm(C.conj().T@C-np.eye(DLOG))<1e-10
    assert np.linalg.norm(dressed_exchange(C)-C)<1e-10
    for vertex in ("A","B"):
        for h in range(G):
            assert np.linalg.norm(transform(C,h,vertex)-C)<1e-10
    return C


def canonical_dihedral_trivial_code(W):
    """Best four-dimensional passive subspace for the WHOLE <A,S> group.

    The two gauge-compatible involutions A=SWAP(A0,A1) and
    S=dressed SWAP(A0,B0) obey A^2=S^2=(AS)^3=I: a D3=S3
    exchange group. On the physical 14D Gauss sector its fixed
    (trivial representation) subspace has rank SIX. Four dimensions
    can thus be stored without *any* active correction for every
    unitary in this small exchange group.

    The choice is error-algebra tailored and not a general DFS for
    all physically permitted gauge-compatible operations.
    """
    if W.shape!=(DIM,14):
        raise ValueError("Expected 14D Gauss isometry")
    U=W.conj().T@edge_cycle(W)
    I=np.eye(14)
    assert np.linalg.norm(U@U@U-I)<1e-10
    P=(I+U+U@U)/3
    assert np.linalg.norm(P-P.conj().T)<1e-10
    assert np.linalg.norm(P@P-P)<1e-10
    rank=int(round(np.trace(P).real))
    assert rank==6,rank
    vectors=[]
    for i in range(14):
        v=P[:,i].copy()
        for _ in range(2):
            for u in vectors:
                v-=u*np.vdot(u,v)
        norm=np.linalg.norm(v)
        if norm>1e-9:
            vectors.append(v/norm)
        if len(vectors)==DLOG:break
    assert len(vectors)==DLOG
    C=W@np.stack(vectors,axis=1)
    assert np.linalg.norm(C.conj().T@C-np.eye(DLOG))<1e-10
    assert np.linalg.norm(dressed_exchange(C)-C)<1e-10
    assert np.linalg.norm(onsite_swap(C)-C)<1e-10
    assert np.linalg.norm(alternate_dressed_exchange(C)-C)<1e-10
    for v in ("A","B"):
        for h in range(G):
            assert np.linalg.norm(transform(C,h,v)-C)<1e-10
    return C


def alternate_dressed_exchange(y):
    """Gauge-preserving across-link dressed SWAP(A1,B0)=A S(A0,B0) A."""
    return onsite_swap(dressed_exchange(onsite_swap(y)))


def passive_fault_metrics(C,p=P_SOURCE):
    """Exact logical process fidelity for common specified dressed errors."""
    if not np.isfinite(p) or not 0<=p<=1:
        raise ValueError("Source-error probability must be [0,1]")
    projected=C.conj().T@dressed_exchange(C)
    alternate=C.conj().T@alternate_dressed_exchange(C)
    identity_fidelity=1.
    same=float((1-p)+p*abs(np.trace(projected)/DLOG)**2)
    alt=float((1-p)+p*abs(np.trace(alternate)/DLOG)**2)
    assert abs(same-1)<1e-10
    return {
        "known_S_A0_B0_one_event_probability":p,
        "same_dressed_exchange_channel_F_e":same,
        "alternate_dressed_exchange_S_A1_B0_channel_F_e":alt,
        "alternate_error_code_projection_trace":
            [float(np.trace(alternate).real),float(np.trace(alternate).imag)],
        "alternate_error_is_same_as_target":False,
        "single_known_pair_exact_decoherence_free":True,
    }


def active_ideal_second_swap_metric(K,p=P_SOURCE):
    """What the *ideal* existing logical pointer + known-S correction
    does if the actual error was a different, gauge-allowed dressed SWAP.

    This is a diagnostic separate from the noisy microscopic PR57 path.
    """
    fidelity=(1-p)
    source=alternate_dressed_exchange(K)*np.sqrt(p)
    out=pointer_circuit(source)
    for m,branch in enumerate(out):
        corrected=dressed_exchange(branch) if m==2 else branch
        fidelity+=abs(np.trace(K.conj().T@corrected))**2/(DLOG**2)
    assert 0<=fidelity<=1+1e-10
    return float(fidelity)


def colour_mode_weights(C):
    """Exact logical trace weight over 81 four-active-colour patterns.

    Physical registers: link(g), matter(A0,A1,B0,B1); spectator B1
    and gauge link are traced. Mediator M is entangled with stationary
    reference colour in (|00>+|11>+|22>)/sqrt(3), so each mediator
    colour appears with probability 1/3. Initial particles occupy
    exactly one colour mode per active orbital A0,A1,B0,M.
    """
    if C.shape!=(DIM,DLOG):
        raise ValueError("Expected 486 x 4 code")
    tensor=C.reshape(G,3,3,3,3,DLOG)
    q=(abs(tensor)**2).sum(axis=(0,4,5))
    assert q.shape==(3,3,3)
    assert abs(float(q.sum())-DLOG)<1e-10
    weights=np.broadcast_to(q[:,:,:,None]/3,(3,3,3,3)).reshape(-1).copy()
    assert weights.shape==(81,)
    assert abs(weights.sum()-DLOG)<1e-10
    return weights


def passive_exact_number_dephasing(C,gamma,elapsed_inverse_t,
                                   model="orbital_charge"):
    """Analytic CPTP idling dephasing on all four active colour orbitals.

    Each Lindblad operator L=sqrt(2gamma)n_{site,colour} is diagonal,
    and for any two 4-orbital basis strings a,b,
    coherence multiplier = exp[-2gamma*T * Hamming(a,b)].
    Entanglement fidelity is a quadratic form of diagonal
    projector weights t_x=Tr(P_code|x><x|):
      F_e=(sum_xy t_x t_y damping_xy)/d^2.
    For colour-blind orbital total-number the jumps equal I on the
    fully singly occupied passive code, so F_e=1.
    """
    if model not in ("orbital_charge","colour_resolved"):
        raise ValueError(model)
    if not np.isfinite(gamma) or gamma<0 or not np.isfinite(elapsed_inverse_t) or elapsed_inverse_t<0:
        raise ValueError("nonnegative finite gamma/time")
    if model=="orbital_charge" or gamma==0 or elapsed_inverse_t==0:
        return 1.
    weights=colour_mode_weights(C)
    patterns=np.indices((3,3,3,3)).reshape(4,81).T
    mismatch=np.count_nonzero(patterns[:,None,:]!=patterns[None,:,:],axis=2)
    kernel=np.exp(-2*gamma*elapsed_inverse_t*mismatch)
    fidelity=float(weights@kernel@weights/DLOG**2)
    assert 0<=fidelity<=1+1e-10
    return fidelity


def comparison_row(U,K,C,source_p=P_SOURCE,T1_us=30.,
                   hopping_mhz=10.,transport_ns=100.,readout_ns=300.):
    """Recompute encoded 495D dynamics at established nominal controls.

    A change of hopping MHz here converts fixed optimized dimensionless
    pulses to durations without reoptimizing device-native bandwidth.
    """
    opt=optimize_pulse(U)
    duration=physical_duration_ns(opt["time"],hopping_mhz,
                                  transport_ns,readout_ns)
    sample=process_sample(K,U,opt["active_detuning"],opt["time"],
                          shots=1,gamma=0.,substeps=1,loss_enabled=False)
    fe=sample["mean_F_e_conditional_on_no_loss"]
    survival=exactly_flagged_fermion_survival(4,duration,T1_us)
    passive=passive_fault_metrics(C,source_p)
    return {
        "U_over_t":float(U),
        "equal_four_fermion_time_us":duration/1000,
        "same_active_fermion_T1_us":T1_us,
        "source_identical_dressed_exchange_probability":source_p,
        "active_encoded_with_microscopic_feedback_no_loss_F_e":fe,
        "passive_gauss_admissible_no_loss_F_e":passive[
            "same_dressed_exchange_channel_F_e"],
        "same_loss_survival_both":survival,
        "active_unconditional_F_e":float(survival*fe),
        "passive_unconditional_F_e":float(survival*passive[
            "same_dressed_exchange_channel_F_e"]),
        "passive_minus_active_unconditional_F_e":float(survival*(1-fe)),
        "active_source_error_physical_mechanism":
            "specified gauge-invariant S_dressed on A0/B0",
        "passive_source_error_physical_mechanism":
            "same specified gauge-invariant S_dressed on A0/B0",
    }


def stochastic_matched_dephasing(K,C,U=30.,gamma=.0002,shots=6):
    """Encoded PR57 trajectories vs exact passive idle Lindblad channel.

    Same 4 active fermions, duration, T1, source swap, dephasing
    operator family, rate per mode. Passive idles; active runs
    15 microscopic local Hubbard pulses and 6 high-level link transfers.
    """
    opt=optimize_pulse(U)
    duration=physical_duration_ns(opt["time"])
    tdim=duration/time_unit_ns()
    survival=exactly_flagged_fermion_survival(4,duration,30.)
    output=[]
    for i,model in enumerate(("orbital_charge","colour_resolved")):
        active=process_sample(K,U,opt["active_detuning"],opt["time"],
            model=model,gamma=gamma,shots=shots,seed=126+i,
            T1_us=30.,substeps=2)
        passive=passive_exact_number_dephasing(C,gamma,tdim,model)
        output.append({
            "noise_model":model,
            "gamma_per_occupation_Lindblad":gamma,
            "full_time_us":duration/1000,
            "identical_active_fermion_loss_survival":survival,
            "same_stochastic_bath_operator_family_both_memories":True,
            "passive_conditional_F_e_exact":passive,
            "passive_unconditional_F_e_exact":survival*passive,
            "active_conditional_F_e_estimate":
                active["mean_F_e_conditional_on_no_loss"],
            "active_unconditional_F_e_estimate":
                active["mean_F_e_unconditional_including_flagged_loss"],
            "active_conditional_MC_standard_error":
                active["Monte_Carlo_standard_error_conditional"],
            "active_MC_shots":shots,
            "passive_exact_minus_active_estimated_unconditional_F_e":
                survival*passive-active[
                    "mean_F_e_unconditional_including_flagged_loss"],
            "gamma_is_equal_per_mode_not_necessarily_equal_total_accumulated_noise":True,
        })
    return output


def lifetime_sensitivity(rows,T1_values=(30.,60.,100.,300.)):
    """Same 4-fermion erasure factor for BOTH architectures.
    Therefore no finite lifetime can invert a no-dephasing
    passive/active ranking when code is exactly S immune.
    """
    results=[]
    for r in rows:
        for T1 in T1_values:
            eta=exactly_flagged_fermion_survival(
                4,r["equal_four_fermion_time_us"]*1000,T1)
            code=r["active_encoded_with_microscopic_feedback_no_loss_F_e"]*eta
            bare=eta
            results.append({
                "U_over_t":r["U_over_t"],
                "T1_us":float(T1),
                "equal_four_fermion_no_jump_probability":eta,
                "passive_unconditional_F_e":bare,
                "active_unconditional_F_e":code,
                "passive_wins":bool(bare>code),
            })
    return results


def run():
    W=physical_isometry()
    K,_,_=code_isometry(W)
    C=canonical_dihedral_trivial_code(W)
    swap_only_C=canonical_swap_passive_code(W)
    p=passive_fault_metrics(C)
    swap_only_p=passive_fault_metrics(swap_only_C)
    alt_active=active_ideal_second_swap_metric(K)
    charged=initialize(C,pointer_uniform=False)
    extended=gauss_acceptance(charged)
    assert min(extended.values())>1-1e-10
    assert np.linalg.norm(dressed_exchange(C)-C)<1e-10
    rows=[comparison_row(U,K,C) for U in (20.,30.,40.)]
    assert all(r["passive_unconditional_F_e"]>
               r["active_unconditional_F_e"] for r in rows)
    sensitivity=lifetime_sensitivity(rows)
    assert all(r["passive_wins"] for r in sensitivity)
    noisy=stochastic_matched_dephasing(K,C)
    return {
        "gauge_group":"S3 subgroup of SU3, two independent vertices",
        "Gauss_physical_sector_dimension":14,
        "known_dressed_swap_plus_sector_dimension":10,
        "common_three_reflection_trivial_sector_dimension":6,
        "weak_swap_only_passive_negative_control":swap_only_p,
        "logical_dimension_each":4,
        "passive_subspace_definition":"canonical 4D subspace of the 6D trivial representation under the whole generated D3 exchange group, local-Gauss invariant",
        "active_subspace_definition":"4D omega eigenphase code of A_onsite S_dressed",
        "same_original_matter_link_registers":True,
        "same_4_active_fermion_orbitals":"A0,A1,B0,charged mediator M",
        "same_spectator_reference_link_initialization":True,
        "passive_reference_code_gauss_acceptance":extended,
        "exact_passive_known_error":p,
        "active_ideal_syndrome_alternate_error_F_e":alt_active,
        "noiseless_loss_matched_rows":rows,
        "T1_sensitivity":sensitivity,
        "colour_number_dephasing_comparison":noisy,
        "results_are_native_device_calibrated":False,
        "limits":[
            "The passive memory is a four-dimensional gauge-invariant code, NOT an unencoded fundamental colour qubit",
            "The main passive code is tailored to the SMALL D3 exchange group generated by onsite A and dressed S; its six-dimensional trivial representation supports four perfectly passive logical states, not all physically allowed exchanges",
            "A deliberately weaker S=+1-only four-state passive code fails against a distinct D3 reflection, but the stronger trivial-sector passive code succeeds; do not mistake the weaker counterexample for an optimal passive reference",
            "The shared specified fault is the same finite-S3 dressed exchange OPERATOR, not a hardware-derived microscopic stochastic fault or continuous SU3 lattice process",
            "Same four ACTIVE fermions suffer uniform flagged loss; spectator B1, charged stationary reference, finite group link and pointer are ideal/no loss",
            "Both state preparations and eventual decoding are assumed; passive eigenstate preparation costs are unmeasured",
            "Encoded microscopic Hamiltonian operation/error and link transports from merged PR57; idle passive store has no operations but equal elapsed time",
            "Colour-resolved bath can break Gauss even though passive memory starts physical; this model does not add a Gauss-enforcing penalty",
            "Active noise estimates have only six Monte Carlo trajectories per bath; passive diagonal idle dynamics are analytically exact",
            "No device routing, native gauge-link pulses, post-jump propagation, Tphi calibration, full SU3 or nine-quark 21D/65-syndrome construction",
            "A single-eigenvalue error-invariant passive subspace may cease to exist for a wider error algebra; no general passive protection claimed"
        ]}


def main():
    data=run()
    Path("colour_s3_gauge_admissible_passive_results.json").write_text(
        json.dumps(data,indent=2,sort_keys=True)+"\n")
    print("PASSIVE_GAUSS_SWAP",json.dumps({
        "passive_known_error":data["exact_passive_known_error"],
        "weaker_passive_negative_control":data["weak_swap_only_passive_negative_control"],
        "active_second_swap_F_e":data["active_ideal_syndrome_alternate_error_F_e"],
        "extended_Gauss":data["passive_reference_code_gauss_acceptance"],
    },sort_keys=True),flush=True)
    for row in data["noiseless_loss_matched_rows"]:
        print("SAME_FAULT_SAME_GAUGE_BARE_V_ACTIVE",json.dumps(row,sort_keys=True),
              flush=True)
    for row in data["colour_number_dephasing_comparison"]:
        print("SAME_NOISE_GAUSS_BARE_V_ACTIVE",json.dumps(row,sort_keys=True),
              flush=True)


if __name__=="__main__":
    main()
