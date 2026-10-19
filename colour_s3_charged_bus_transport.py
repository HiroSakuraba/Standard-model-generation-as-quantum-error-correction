"""Charge-carrying mediator transport and gauge-covariant logical syndrome.

Finite S3 subgroup SU(3) gauge-link test only, *not* 21D nine-quark code.

System: four matter qutrits, one six-level oriented gauge link (486D),
neutral syndrome pointer qutrit; mobile fundamental-charge mediator with
location A/B (6D), plus a stationary anti-fundamental reference qutrit
at A. Because this S3 representation is real orthogonal, the local
mediator/reference singlet is (|00>+|11>+|22>)/sqrt(3).

Gauge transformation on mobile mediator acts as D(h_A) when at A or
D(h_B) when at B. Charge transport A->B uses D(g^-1), B->A D(g):
each T is a gauge-invariant TWO-register link+mobile-mediator gate.
Onsite swap of the mobile charged mediator and one matter qudit is a
gauge-invariant TWO-register gate; conditioning it on a neutral qutrit
pointer promotes it to a THREE-register controlled gate.

Identity: SWAP_m,A T SWAP_m,B T SWAP_m,A exactly implements the
gauge-link-dressed S_A,B while returning the mediator to A unchanged.
Three local swaps similarly implement an onsite A0<->A1 swap. Replace
both controlled S and controlled A in the earlier three-cycle phase
syndrome with these constructions: **16 gauge-covariant operations with
at most THREE-register support** (4 two-register transports and 12
three-register controlled local swaps), plus pointer preparation/DFT.
Measured feedback uses 5 gauge-covariant two-register operations.

This *escapes the neutral-ancilla pair-gate restriction* of PR #50
because the mediator is GAUGE CHARGED. It does NOT imply a native
two-body realization of controlled onsite Fredkin, local charge
preparation, remote transfer, nor hardware time/noise advantage.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from colour_gauge_link_s3 import (
    DIM, G, MATTER, REP, inverse, transform, physical_isometry,
    dressed_exchange, projector,
)
from colour_s3_edge_logical_syndrome import (
    code_isometry, pointer_circuit, W3,
)

# Pointer (3), gauge system (6 link * 3^4 matter = 486),
# mediator location (A=0,B=1), mediator colour (3), reference colour (3),
# arbitrary encoded logical-reference columns.
P = 3
MEDLOC = 2
COLOUR = 3
SOURCE_DIM = 4
GAUGE_DIMS = (G,3,3,3,3)
SCHEDULE=(("edge",(1,2)),("onsite",(1,2)),
          ("edge",(2,)),("onsite",(2,)))


def initialize(K, pointer_uniform=True):
    """Prepare a GAUSS-INVARIANT mediator/reference charged singlet at A."""
    x=np.asarray(K,complex)
    assert x.ndim==2 and x.shape[0]==DIM
    z=np.zeros((P,DIM,MEDLOC,COLOUR,COLOUR,x.shape[1]),complex)
    for c in range(COLOUR):
        if pointer_uniform:
            z[:,:,0,c,c,:]=x[None,:,:]/3
        else:
            z[0,:,0,c,c,:]=x/np.sqrt(3)
    assert abs(np.vdot(z,z)-np.vdot(x,x))<5e-10
    return z


def gauge_action(z,h,vertex):
    """Both local Gauss transformations including mobile and ref charges."""
    x=np.asarray(z,complex)
    assert x.ndim==6 and x.shape[:5]==(P,DIM,2,3,3)
    if vertex not in ("A","B"):raise ValueError(vertex)
    out=np.empty_like(x)
    D=REP[h]
    for p in range(P):
        for loc in range(2):
            data=transform(x[p,:,loc].reshape(DIM,-1),h,vertex)
            data=data.reshape(DIM,3,3,x.shape[-1])
            if (vertex=="A" and loc==0) or (vertex=="B" and loc==1):
                data=np.einsum("ab,ibcr->iacr",D,data)
            if vertex=="A":
                # Reference carries conjugate fundamental at vertex A.
                data=np.einsum("ab,imbr->imar",D.conj(),data)
            out[p,:,loc]=data
    return out


def transport(z):
    """T: A->B with D(g^-1); B->A with D(g); exact T^2=I."""
    x=np.asarray(z,complex)
    assert x.ndim==6 and x.shape[:5]==(P,DIM,2,3,3)
    r=x.shape[-1]
    src=x.reshape(P,G,MATTER,2,3,3,r)
    dst=np.empty_like(src)
    for g in range(G):
        dst[:,g,:,1,:,:,:]=np.einsum(
            "ab,pmbrk->pmark", REP[inverse(g)],src[:,g,:,0,:,:,:])
        dst[:,g,:,0,:,:,:]=np.einsum(
            "ab,pmbrk->pmark", REP[g],src[:,g,:,1,:,:,:])
    return dst.reshape(x.shape)


def local_swap(z,site,active=(0,1,2)):
    """Conditional mediator-matter swap at occupied vertex.

    High-level register support: syndrome pointer + mobile (loc/charge)
    + local matter qutrit. An inactive pointer branch is untouched.
    """
    if site not in ("A0","A1","B0","B1"):raise ValueError(site)
    x=np.asarray(z,complex)
    assert x.ndim==6 and x.shape[:5]==(P,DIM,2,3,3)
    out=x.copy()
    loc=0 if site[0]=="A" else 1
    matter_axis={"A0":1,"A1":2,"B0":3,"B1":4}[site]
    r=x.shape[-1]
    u=out.reshape(P,G,3,3,3,3,2,3,3,r)
    for p in active:
        if p not in (0,1,2):raise ValueError(active)
        # after indexing p and loc, axes are (link,A0,A1,B0,B1,med,ref,R)
        slice_of_charge=u[p,:,:,:,:,:,loc,:,:,:]
        u[p,:,:,:,:,:,loc,:,:,:]=np.swapaxes(
            slice_of_charge,matter_axis,5)
    return out


def involution_pulse(z,op,theta=np.pi):
    """Continuous extension, for coherent overrotations of charge gates."""
    f=np.exp(-1j*theta)
    return (1+f)*z/2 + (1-f)*op(z)/2


def edge_bus(z,active=(0,1,2),angle=np.pi):
    """Five exact GAUSS-INVARIANT operations; acts on active pointer arms.

    Unconditional T moves inactive pointer arms across the link AND
    back; those arms return identically. This avoids pointer-controlled
    transport and its larger register support.
    """
    y=involution_pulse(z,lambda a:local_swap(a,"A0",active),angle)
    y=involution_pulse(y,transport,angle)
    y=involution_pulse(y,lambda a:local_swap(a,"B0",active),angle)
    y=involution_pulse(y,transport,angle)
    y=involution_pulse(y,lambda a:local_swap(a,"A0",active),angle)
    return y


def onsite_bus(z,active=(0,1,2),angle=np.pi):
    """Three pointer-controlled mediator swaps implement SWAP(A0,A1)."""
    y=involution_pulse(z,lambda a:local_swap(a,"A0",active),angle)
    y=involution_pulse(y,lambda a:local_swap(a,"A1",active),angle)
    y=involution_pulse(y,lambda a:local_swap(a,"A0",active),angle)
    return y


def schedule_ops():
    """16 elementary charge-covariant operations, 12 conditional + 4 T."""
    operations=[]
    for kind,active in SCHEDULE:
        if kind=="edge":
            operations.extend([
                ("local_A0",lambda x,a=active:local_swap(x,"A0",a),3),
                ("transport",transport,2),
                ("local_B0",lambda x,a=active:local_swap(x,"B0",a),3),
                ("transport",transport,2),
                ("local_A0",lambda x,a=active:local_swap(x,"A0",a),3)])
        else:
            operations.extend([
                ("local_A0",lambda x,a=active:local_swap(x,"A0",a),3),
                ("local_A1",lambda x,a=active:local_swap(x,"A1",a),3),
                ("local_A0",lambda x,a=active:local_swap(x,"A0",a),3)])
    assert len(operations)==16
    assert sum(n==2 for _,_,n in operations)==4
    assert sum(n==3 for _,_,n in operations)==12
    return operations


def run_syndrome(z,overrotation=0.,history=False):
    """Pointer operations followed by exact three-point phase readout."""
    if not np.isfinite(overrotation):raise ValueError("finite error")
    snapshots=[]
    x=z
    theta=np.pi*(1+overrotation)
    for name,op,support in schedule_ops():
        x=involution_pulse(x,op,theta)
        if history:snapshots.append((name,support,x.copy()))
    result=np.fft.fft(x,axis=0)/np.sqrt(3)
    return (result,snapshots) if history else result


def correction(z,overrotation=0.):
    """After pointer readout, classically correct ONLY outcome label 2.

    Five two-register gauge-covariant operations (pointer is now a
    classical selector), though array layout keeps all labels.
    """
    return edge_bus(z,active=(2,),angle=np.pi*(1+overrotation))


def remove_ancilla_singlet(z):
    """Diagnostic only! Coherent contraction, NEVER used to compute noisy F_e."""
    x=np.asarray(z)
    return sum(x[:,:,0,c,c,:] for c in range(3))/np.sqrt(3)


def exact_unconditional_fidelity(K,p=.37,overrotation=0.):
    """All pointer/mediator/reference outcomes are traced, no postselection.

    Each branch carries a 486x4 operator mapping the original logical
    ququart. Sum |Tr(K^dag E_branch)|² / 4² over all environmental labels
    and both source noise Kraus operators.
    """
    if not (0<=p<=1):raise ValueError(p)
    fe=0.
    total=0.
    pointer_prob=np.zeros(3)
    for E in (np.sqrt(1-p)*K,np.sqrt(p)*dressed_exchange(K)):
        z=run_syndrome(initialize(E),overrotation)
        z=correction(z,overrotation)
        pointer_prob+=np.array([np.vdot(b,b).real/4 for b in z])
        for label in range(3):
            for loc in range(2):
                for c in range(3):
                    for ref in range(3):
                        b=z[label,:,loc,c,ref,:]
                        fe+=abs(np.trace(K.conj().T@b))**2/16
                        total+=np.vdot(b,b).real/4
    assert abs(total-1)<2e-10
    assert -2e-10<=fe<=1+2e-10
    return {"source_known_dressed_swap_probability":p,
            "coherent_per_primitive_fractional_overrotation":overrotation,
            "unconditional_entanglement_fidelity":float(fe),
            "total_output_probability":float(total),
            "pointer_probabilities_0_1_2":pointer_prob.tolist()}


def gauss_acceptance(z):
    """Project onto Gauss sectors of FULL mediator+reference system."""
    norm=np.vdot(z,z).real
    assert norm>0
    output={}
    for vertex in ("A","B"):
        proj=sum((gauge_action(z,h,vertex) for h in range(G)))/G
        output[vertex]=float(np.vdot(proj,proj).real/norm)
    return output


def run():
    W=physical_isometry()
    K,F,_=code_isometry(W)
    x=initialize(K)
    init=gauss_acceptance(x)
    assert min(init.values())>1-1e-10
    output,history=run_syndrome(x,history=True)
    # Source PR48 logical-syndrome acting on original physical system,
    # with original charged Bell-reference returned and disentangled.
    extracted=remove_ancilla_singlet(output)
    reference=pointer_circuit(K)
    difference=float(np.linalg.norm(extracted-reference))
    assert difference<1e-10
    assert abs(np.vdot(extracted,extracted)-np.vdot(output,output))<1e-10
    gauss_trace=[{"stage":"initial","accept":init}]
    for i,(name,support,state) in enumerate(history):
        p=gauss_acceptance(state)
        gauss_trace.append({"stage":i+1,"gate":name,
                            "register_support":support,"accept":p})
    final=gauss_acceptance(output)
    gauss_trace.append({"stage":"pointer_fourier","accept":final})
    minima={v:min(row["accept"][v] for row in gauss_trace)
            for v in ("A","B")}
    assert min(minima.values())>1-1e-10
    # Exact representation covariance on arbitrary non-Gauss inputs,
    # checking each distinct elementary operation at both vertices.
    rng=np.random.default_rng(19)
    arbitrary=rng.normal(size=x.shape)+1j*rng.normal(size=x.shape)
    unique_ops=(("transport",transport),
                ("A0",lambda z:local_swap(z,"A0",(1,2))),
                ("A1",lambda z:local_swap(z,"A1",(2,))),
                ("B0",lambda z:local_swap(z,"B0",(1,2))))
    worst=0.
    for name,op in unique_ops:
        for vertex in ("A","B"):
            for h in range(G):
                lhs=gauge_action(op(arbitrary),h,vertex)
                rhs=op(gauge_action(arbitrary,h,vertex))
                worst=max(worst,float(np.linalg.norm(lhs-rhs)))
    assert worst<1e-10
    results=[exact_unconditional_fidelity(K,overrotation=e)
             for e in (0.,.005,.02,.05)]
    assert abs(results[0]["unconditional_entanglement_fidelity"]-1)<1e-10
    return {
        "model":"S3 two-vertex four-matter-qutrit gauge link, mobile charged mediator and reference singlet",
        "original_physical_dimension":DIM,
        "mobile_mediator_dimension":6,
        "stationary_reference_dimension":3,
        "neutral_pointer_dimension":3,
        "code_dimension":4,
        "reference_gauge_invariant_preparation":"mobile fundamental + stationary conjugate fundamental in local maximally entangled singlet at vertex A",
        "two_register_link_transport_count":4,
        "three_register_pointer_conditioned_onsite_swaps":12,
        "total_syndrome_elementary_count":16,
        "max_elementary_register_support":3,
        "unconditional_feedback_two_register_operations":5,
        "exact_high_level_PR48_syndrome_residual_fro":difference,
        "maximum_elementary_gauge_covariance_residual_fro":worst,
        "intermediate_Gauss_acceptance":gauss_trace,
        "minimum_Gauss_acceptance":minima,
        "coherent_calibration_sensitivity":results,
        "science_boundary":[
            "Exact gauge-covariant 2-register charged transfer and 3-register controlled LOCAL exchanges, not native 2-body gate compilation",
            "Requires preparing a fundamental/conjugate singlet and a mobile 6-dimensional location+colour register, with realistic movement/routing unspecified",
            "No pulse strength, dimensional gate time, T1/T2, mediator/link leakage, reset or link loss",
            "Overrotation scan is coherent deterministic control error, not calibrated noise",
            "Only one known dressed-swap error on the finite S3 4D phase code, NOT SU3 nine-quark 21D 65-syndrome code",
            "Neutral-pointer strict two-register obstruction of PR50 remains valid; charged mediator changes the assumptions"
        ]}


def main():
    result=run()
    Path("colour_s3_charged_bus_transport_results.json").write_text(
        json.dumps(result,indent=2,sort_keys=True)+"\n")
    print("S3_CHARGED_BUS_RESULTS",json.dumps({
        k:v for k,v in result.items() if k not in
        ("intermediate_Gauss_acceptance","science_boundary")
    },sort_keys=True),flush=True)


if __name__=="__main__":
    main()
