"""Exact finite-S3 link encoding and selective-CP22 gate-cost witnesses.

THIS IS NOT A PHYSICAL DEVICE MODEL. The 6-level oriented S3 link can
be encoded as parity qubit (2) x cyclic qutrit (3). In its fixed local
Fourier-cycle encoding, the central standard-electric projector is
I_parity x diag(0,1,1). Thus a neutral-pointer-controlled electric
rotation uses exactly TWO selective two-qutrit |22> phases, modulo
local qutrit permutations and the link-basis encoding choice.

Also synthesize finite-S3 local same-colour projectors with THREE CP22
and an onsite two-qutrit SWAP with TWELVE selective CP22 (nonoptimal
upper bound via three modular SUMs).

Crucial: the individual CP22 factors do not generally commute with BOTH
local Gauss constraints even when their complete product does. No
claim of gauge-preserving intermediate pulses, native charged transport,
frequency-selective waveforms, device-calibrated arbitrary-angle CP22,
or completion of PR64's full syndrome instrument follows.
"""
from __future__ import annotations

import json
from pathlib import Path
import numpy as np

from colour_gauge_link_s3 import (
    PERMS, G, DIM, MATTER, compose, inverse, ID,
    REP, transform, physical_isometry, projector,
)
from colour_s3_strict_gauss_pair_nogo import (
    link_central_projectors, apply_link,
)

OMEGA = np.exp(2j*np.pi/3)


def _parity(index):
    perm = PERMS[index]
    return (-1)**sum(perm[a]>perm[b] for a in range(3)
                     for b in range(a+1,3))


def parity_cycle_order():
    """Group elements [e,r,r²,s,sr,sr²]; r a fixed 3-cycle, s odd."""
    r=PERMS.index((1,2,0))
    s=next(i for i in range(G) if _parity(i)==-1)
    cycle=(ID,r,compose(r,r))
    order=tuple(compose(p,g) for p in (ID,s) for g in cycle)
    assert len(set(order))==G
    return order


def cycle_fourier():
    """Rows t, columns k of the unitary 3-point character transform."""
    return np.array([[OMEGA**(-t*k)/np.sqrt(3)
                      for k in range(3)] for t in range(3)],complex)


def link_encoding_unitary():
    """Original |g> link basis -> |parity> tensor |Fourier cycle>.

    This is a CHOICE OF ENCODING. If hardware stores the group basis,
    physically executing the Fourier basis change is an ADDITIONAL cost,
    not a free conversion.
    """
    F=cycle_fourier()
    order=parity_cycle_order()
    U=np.zeros((G,G),complex)
    for parity in range(2):
        for mode in range(3):
            for cyclic in range(3):
                U[3*parity+mode,order[3*parity+cyclic]]=F[mode,cyclic]
    assert np.linalg.norm(U.conj().T@U-np.eye(G))<1e-12
    return U


def standard_projector_diagonal_audit():
    U=link_encoding_unitary()
    E=link_central_projectors()["standard"]
    diagonal=np.diag([0,1,1,0,1,1]).astype(complex)
    residual=float(np.linalg.norm(U@E@U.conj().T-diagonal))
    if residual>2e-12:
        raise AssertionError("S3 central electric projector does not diagonalize")
    return E,diagonal,U,residual


def _level_to_two(level):
    if level not in (0,1,2):
        raise ValueError("invalid qutrit basis value")
    P=np.eye(3, dtype=complex)
    if level!=2:
        P[[2,level]]=P[[level,2]]
    return P


def cp22_pair(angle):
    """One native-selective-phase IDEALIZATION, diag_{|22>}(e^-i angle)."""
    if not np.isfinite(angle):
        raise ValueError("finite conditional phase required")
    u=np.eye(9,dtype=complex)
    u[8,8]=np.exp(-1j*angle)
    return u


def targeted_cp22(control, target, angle):
    """Exactly one CP22 conjugated by local basis permutations."""
    local=np.kron(_level_to_two(control),_level_to_two(target))
    return local.conj().T@cp22_pair(angle)@local


def lift_pointer_cycle(gate9):
    """Two-qutrit gate on neutral pointer and cycle; parity qubit idle."""
    gate9=np.asarray(gate9,complex)
    if gate9.shape!=(9,9):
        raise ValueError("expected 9x9 pointer-cycle gate")
    out=np.zeros((18,18),complex)
    for p_out in range(3):
        for k_out in range(3):
            for p_in in range(3):
                for k_in in range(3):
                    element=gate9[p_out*3+k_out,p_in*3+k_in]
                    for parity in range(2):
                        out[p_out*6+parity*3+k_out,
                            p_in*6+parity*3+k_in]=element
    return out


def pointer_electric_cp22_mapping(angle):
    """Two ideal CP22 primitives reproduce exp(-i angle |1><1| Pstd)."""
    if not np.isfinite(angle):
        raise ValueError("finite phase required")
    E,_,U,_=standard_projector_diagonal_audit()
    native=lift_pointer_cycle(targeted_cp22(1,2,angle))@(
        lift_pointer_cycle(targeted_cp22(1,1,angle)))
    basis=np.kron(np.eye(3),U)
    original=basis.conj().T@native@basis
    pointer=np.diag([0.,1.,0.])
    projector_term=np.kron(pointer,E)
    ideal=np.eye(18,dtype=complex)+(np.exp(-1j*angle)-1)*projector_term
    residual=float(np.linalg.norm(original-ideal))
    if residual>2e-11:
        raise AssertionError("two CP22 synthesis differs from central gate")
    return {
        "angle_radians":float(angle),
        "original_group_basis_unitary":original,
        "encoded_parity_cycle_unitary":native,
        "ideal_central_unitary":ideal,
        "unitary_residual":residual,
        "selective_CP22_count":2,
        "registers":("pointer_qutrit","link_parity_qubit","link_cycle_qutrit"),
        "CP22_support":("pointer_qutrit","link_cycle_qutrit"),
        "link_parity_register_idle_during_entangler":True,
        "link_fourier_is_encoding_not_free_hardware_pulse":True,
    }


def electric_link_local_phases(angle):
    """Single-link E rotation: two local qutrit phase shifts after encoding."""
    if not np.isfinite(angle):
        raise ValueError("finite angle required")
    E,D,U,_=standard_projector_diagonal_audit()
    encoded=np.eye(6)+(np.exp(-1j*angle)-1)*D
    canonical=U.conj().T@encoded@U
    expected=np.eye(6)+(np.exp(-1j*angle)-1)*E
    residual=float(np.linalg.norm(canonical-expected))
    assert residual<1e-11
    return {"original_group_basis_unitary":canonical,
            "encoded_link_unitary":encoded,
            "conditional_two_qutrit_CP22_count":0,
            "single_cycle_qutrit_level_phases":2,
            "residual":residual}


def same_colour_cp22_mapping(angle):
    """At ONE vertex, δ(a=b) projector from three conditional CP22 gates."""
    if not np.isfinite(angle):
        raise ValueError("finite angle required")
    out=np.eye(9,dtype=complex)
    for i in range(3):
        out=targeted_cp22(i,i,angle)@out
    expected=np.diag([
        np.exp(-1j*angle) if a==b else 1.
        for a in range(3) for b in range(3)]).astype(complex)
    error=float(np.linalg.norm(out-expected))
    assert error<1e-11
    return {"unitary":out,"residual":error,
            "selective_CP22_count":3,
            "gauge_invariant_entire_operator_under_S3":True,
            "individual_selective_CP22_factors_gauge_invariant":False}


def _qutrit_fourier():
    return np.conj(cycle_fourier())


def controlled_sum_via_cp22(control,target,sign=1):
    """Generalized SUM, |a,b> -> |a,b+a> or reverse target registers.

    One modular controlled SUM is a local Fourier sandwich of four
    selective CP22 phases; sign=±1 selects add/subtract.
    """
    if (control,target) not in ((0,1),(1,0)) or sign not in (-1,1):
        raise ValueError("invalid directed qutrit SUM")
    F=_qutrit_fourier()
    rotation=np.kron(F,np.eye(3)) if target==0 else np.kron(np.eye(3),F)
    phase=np.eye(9,dtype=complex)
    for a in (1,2):
        for b in (1,2):
            # exp(+2πi sign ab/3) = exp(-i angle) in CP22 convention
            angle=-2*np.pi*sign*a*b/3
            phase=targeted_cp22(a,b,angle)@phase
    return rotation.conj().T@phase@rotation


def onsite_swap_cp22_mapping():
    """Explicit (nonoptimal) 12-CP22 upper bound for onsite SWAP."""
    a=controlled_sum_via_cp22(0,1,1)
    b=controlled_sum_via_cp22(1,0,-1)
    neg=np.zeros((3,3),complex)
    for i in range(3):
        neg[(-i)%3,i]=1
    compiled=np.kron(neg,np.eye(3))@a@b@a
    ideal=np.zeros((9,9),complex)
    for i in range(3):
        for j in range(3):
            ideal[3*j+i,3*i+j]=1
    residual=float(np.linalg.norm(compiled-ideal))
    if residual>2e-10:
        raise AssertionError("CP22-to-SUM-to-SWAP identity failed")
    return {"unitary":compiled,"ideal":ideal,
            "residual":residual,
            "selective_CP22_count":12,
            "single_qutrit_fourier_count":6,
            "single_qutrit_negation_count":1,
            "individual_CP22_factors_strictly_Gauss_preserving":False}


def link_selector_gauss_audit(angle=.7):
    """Commutators and full-register intermediate Gauss violation.

    The first selective cycle-mode phase, as a physical link operation,
    is not invariant under the independent RIGHT S3 vertex action.
    The complete central electric product is invariant at both vertices.
    """
    if not np.isfinite(angle):
        raise ValueError("finite angle")
    E,_,U,_=standard_projector_diagonal_audit()
    native_modes=[]
    for mode in (1,2):
        diag=np.diag([float(j==mode) for p in range(2)
                      for j in range(3)])
        native_modes.append(U.conj().T@diag@U)
    errors={}
    for vertex in ("A","B"):
        worst_partial=0.
        worst_complete=0.
        for h in range(G):
            transform_link=np.zeros((6,6),complex)
            for g in range(G):
                dest=compose(h,g) if vertex=="A" else compose(g,inverse(h))
                transform_link[dest,g]=1
            partial=native_modes[0]
            worst_partial=max(worst_partial,float(np.linalg.norm(
                partial@transform_link-transform_link@partial)))
            worst_complete=max(worst_complete,float(np.linalg.norm(
                E@transform_link-transform_link@E)))
        errors[vertex]={"maximum_single_mode_projector_commutator":
                            worst_partial,
                        "maximum_complete_central_projector_commutator":
                            worst_complete}
    W=physical_isometry()
    first=np.eye(6)+(np.exp(-1j*angle)-1)*native_modes[0]
    combined=np.eye(6)+(np.exp(-1j*angle)-1)*E
    partial_states=apply_link(W,first)
    final_states=apply_link(W,combined)
    gauss={}
    for vertex in ("A","B"):
        intermediate=projector(partial_states,vertex)
        complete=projector(final_states,vertex)
        gauss[vertex]={
            "after_one_selective_phase":float(
                np.linalg.norm(intermediate,"fro")**2/14),
            "after_complete_electric_rotation":float(
                np.linalg.norm(complete,"fro")**2/14),
        }
    assert min(v["after_complete_electric_rotation"]
               for v in gauss.values())>1-1e-10
    return {"per_vertex":errors,"physical_Gauss_acceptance":gauss,
            "single_CP22_intermediate_symmetry_breaking":True,
            "complete_product_Gauss_invariant":True,
            "Gauss_penalty_or_native_simultaneous_interaction_needed":True}


def report():
    angles=(np.pi/7,np.pi/2,np.pi)
    rows=[]
    for angle in angles:
        ptr=pointer_electric_cp22_mapping(angle)
        eq=same_colour_cp22_mapping(angle)
        electric=electric_link_local_phases(angle)
        rows.append({
            "angle":float(angle),
            "pointer_link_electric_CP22":ptr["selective_CP22_count"],
            "pointer_link_electric_exact_residual":ptr["unitary_residual"],
            "onsite_equality_CP22":eq["selective_CP22_count"],
            "onsite_equality_exact_residual":eq["residual"],
            "link_electric_local_qutrit_phases":
                electric["single_cycle_qutrit_level_phases"],
            "link_electric_residual":electric["residual"],
        })
    swap=onsite_swap_cp22_mapping()
    return {
        "study":"exact CP22-based gauge-register hardware mapping witness",
        "group":"finite signed-permutation S3; NOT continuum SU3",
        "link_encoding":"parity qubit x 3-cycle Fourier qutrit",
        "link_original_group_element_order":list(parity_cycle_order()),
        "link_standard_electric_rank":4,
        "exact_gate_mappings":rows,
        "onsite_qutrit_SWAP_ideal_CP22_upper_bound":
            swap["selective_CP22_count"],
        "onsite_SWAP_residual":swap["residual"],
        "local_Gauss_diagnostics":link_selector_gauss_audit(),
        "explicitly_NOT_compiled":[
            "Charged-mediator link transport and local charge reference preparation",
            "Device-specific CP22 pulse waveforms for all angles and spectator states",
            "Simultaneously gauge-preserving pulse paths for selective CP22 factors",
            "Physical pointer readout, feedback, reset and leaked-sector dynamics",
            "Complete 42D three-outcome syndrome instrument",
        ],
        "limitations":[
            "Counts are algebraic upper bounds of ideal selective CP22 primitives, NOT measured hardware gate counts.",
            "Link parity x Fourier-cycle basis is an encoding choice; preparing/using it incurs unspecified mapping/routing cost.",
            "Intermediate individually selective CP22 rotations need not preserve local S3 Gauss constraints.",
            "The full finite-S3 recovery alphabet includes charged transport and pointer drives that remain uncompiled.",
            "This is NOT the separate five-plaquette continuous-SU3 CP22 implementation.",
        ],
    }


def main():
    result=report()
    Path("colour_s3_cp22_link_mapping_results.json").write_text(
        json.dumps(result,indent=2,sort_keys=True)+"\n")
    print("CP22_LINK_MAPPING",json.dumps({
        "pointer_link_CP22":2,
        "same_colour_CP22":3,
        "onsite_SWAP_CP22":12,
        "gauss":result["local_Gauss_diagnostics"],
        "mapping_done":False,
    },sort_keys=True))
    return result


if __name__=="__main__":
    main()
