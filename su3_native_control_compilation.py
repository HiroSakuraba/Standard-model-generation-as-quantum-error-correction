"""SU(3) 5-face Z3 refocusing hardware-control audit.

A) Local cubic-Casimir clocks, subset of addressed links, finite-width pulses.
B) Direct phase-programmed magnetic plaquette quadratures (conditional technology).

No hardware demonstrated, all times and strengths dimensionless.
"""
import json
import numpy as np
import sympy as sp
from scipy.sparse import csr_matrix, diags
from scipy.sparse.linalg import expm_multiply
from su3_cubic_casimir_echo import (construct, exact_clock_test,
    entanglement_fidelity, evolution, jump_operator)
from su3_channel_recovery import syndrome_decoder, source_data
from su3_four_five_plaquette import electric_links, one_site_actions

N=5;DIM=243;T=0.6;THETA=2*np.pi/3
SUBSETS=((4,),(3,4),(2,3,4),(1,2,3,4),(0,1,2,3,4))
ORDER=(0,1,2,1,0);WEIGHTS=(1/6,1/6,1/3,1/6,1/6)

def sympy_checks(states):
    omega=(-1+sp.I*sp.sqrt(3))/2
    assert sp.simplify(sum(omega**k for k in (0,1,2)))==0
    val={0:sp.Integer(0),1:sp.Rational(10,9),2:-sp.Rational(10,9)}
    f=lambda n:sp.Rational(9,10)*val[n]
    assert [f(k) for k in range(3)]==[0,1,-1]
    # In a gauge-fixed face-qutrit register the clock is a *single* Z(12) rotation:
    # diag(0,1,-1) on the qutrit, NOT a two-qudit gate.
    z12=sp.diag(0,1,-1)
    assert sp.diag(*[f(k) for k in range(3)])==z12
    clock_gate=sp.diag(1,omega,omega**2)
    assert all(sp.simplify(clock_gate[k,k]-omega**z12[k,k])==0 for k in range(3))
    for site in range(N):
        ntrans=0
        for state in states:
            for to,coeff in one_site_actions(state,site,'W').items():
                assert (f(to[site])-f(state[site])-1) % 3 == 0
                ntrans+=1
        assert ntrans>200
    a=(0,0,0,0,1);b=(0,0,0,0,2)
    assert electric_links(a)==electric_links(b)
    # The same electric eigenvalues imply identical phases for ANY f(E_0^2,...E_15^2).
    # Thus no quadratic-Casimir-only diagonal control can change <a|B4|b>.
    assert one_site_actions(a,4,'B')[b]==1
    return 'Exact modulo-3 Wilson charges and quadratic-only diagonal pulse obstruction PASS'

def phase_mask(states,subset):
    return np.array([sum(0 if q[i]==0 else (1 if q[i]==1 else -1) for i in subset) for q in states],dtype=float)

def toggled_matrices(d,W,subset,phase_error=None,crosstalk=0):
    """Direct complex Wilson-phase programming and equivalent Q H Q† states.
    Nonselected sites can have optional phase crosstalk coefficient.
    """
    phase_error=phase_error or {}
    factors=[(1+phase_error.get(i,0)) if i in subset else crosstalk for i in range(N)]
    mats=[]
    for k in (0,1,2):
        h=diags(d,format='csr',dtype=complex)
        for i,w in enumerate(W):
            phi=k*THETA*factors[i]
            h=h-(np.exp(1j*phi)*w+np.exp(-1j*phi)*w.getH())
        mats.append(h.tocsr())
    return mats

def conjugation_crosscheck(states,d,W,H,subset,matrices):
    s=phase_mask(states,subset);hi=H.tocoo()
    for k,hk in enumerate(matrices):
        z=csr_matrix((hi.data*np.exp(1j*k*THETA*(s[hi.row]-s[hi.col])),(hi.row,hi.col)),shape=H.shape)
        assert np.max(np.abs((hk-z).data),initial=0)<5e-13
        assert np.max(np.abs((hk-hk.getH()).data),initial=0)<5e-13
    expected=diags(d,format='csr')-sum((w+w.getH() for i,w in enumerate(W) if i not in subset),csr_matrix((DIM,DIM),dtype=complex))
    avg=sum(matrices,csr_matrix((DIM,DIM),dtype=complex))/3
    assert np.max(np.abs((avg-expected).data),initial=0)<2e-13
    # Linear combination of Hermitian quadratures, both individually gauge invariant
    return True

def miscalibrated_electric_clock(states,H,d,subset,fraction):
    # Electric cubic-Casimir pulse miscalibration is representation-dependent.
    # In particular 2 -> 0 differs by -2 levels rather than +1. This is
    # NOT equivalent to a constant phase error on every Wilson transition.
    s=phase_mask(states,subset)
    mat=H.tocoo(); h=[]
    for k in (0,1,2):
        p=np.exp(1j*k*THETA*(1+fraction)*(s[mat.row]-s[mat.col]))
        h.append(csr_matrix((mat.data*p,(mat.row,mat.col)),shape=H.shape))
    return h

def lab_cubic_pulses(states,H,K,subset,cycles,tau,phase_error=None,crosstalk=0):
    """4 instantaneous-clock-equivalent physical pulses per symmetric cycle.
    Full native H remains on during all pulse windows of width tau.
    Per-addressed-site control coefficient = theta/tau (in Hamiltonian units).
    """
    errs=phase_error or {}
    addr={i:((1+errs.get(i,0)) if i in subset else crosstalk) for i in range(N)}
    clock=np.array([sum(v*(0 if q[i]==0 else (1 if q[i]==1 else -1)) for i,v in addr.items()) for q in states],float)
    Dclock=diags(clock,format='csr')
    y=K.copy()
    for _ in range(cycles):
        for k,w in enumerate(WEIGHTS):
            y=expm_multiply(-1j*(w*T/cycles)*H,y,traceA=0)
            if k<4:
                sign=1 if k<2 else -1
                y=expm_multiply((-1j*tau)*(H+sign*THETA/tau*Dclock),y,traceA=0)
    return y

def ramped_drive(states,d,W,K,subset,cycles,tau,substeps=4):
    """Magnetic drive phase ramps between 0,theta,2theta,theta,0,0.
    During each ramp the complex W coupling stays on; time is added.
    Requires physical programmable Wilson quadratures (not just virtual frame bookkeeping).
    """
    assert substeps>=1
    y=K.copy();levels=np.array((0,1,2,1,0),float)
    for _ in range(cycles):
        for j,w in enumerate(WEIGHTS):
            h=toggled_matrices(d,W,subset)[int(levels[j])]
            y=expm_multiply((-1j*w*T/cycles)*h,y,traceA=0)
            if j<4:
                first,second=levels[j:j+2]
                for t in range(substeps):
                    phase=THETA*(first+(second-first)*(t+.5)/substeps)
                    h=diags(d,format='csr',dtype=complex)
                    for i,Wi in enumerate(W):
                        phi=phase if i in subset else 0.
                        h=h-(np.exp(1j*phi)*Wi+np.exp(-1j*phi)*Wi.getH())
                    y=expm_multiply((-1j*tau/substeps)*h.tocsr(),y,traceA=0)
    return y

def main():
    states,K,d,W,H=construct()
    V,_,G=syndrome_decoder(*source_data())
    assert G.rank()==5
    assert sympy_checks(states)
    rows=[]
    native=expm_multiply(-1j*T*H,K,traceA=0)
    rows.append({'type':'native_no_control','raw':entanglement_fidelity([native],K),'recovered':entanglement_fidelity([native],K,V)})
    ideal={}
    for subset in SUBSETS:
        hs=toggled_matrices(d,W,subset)
        assert conjugation_crosscheck(states,d,W,H,subset,hs)
        for cycles in (2,4,8):
            y=evolution(K,hs,T,cycles)
            row={'type':'subset_refocus','subset':list(subset),'cycles':cycles,'pulse_count':4*cycles,'link_phase_actions':4*cycles*len(subset), 'raw':entanglement_fidelity([y],K),'recovered':entanglement_fidelity([y],K,V)}
            rows.append(row)
            if cycles==8:ideal[subset]=row
            print('IDEAL',json.dumps(row,sort_keys=True),flush=True)
    assert ideal[(1,2,3,4)]['recovered']>0.99999
    assert ideal[(4,)]['recovered']<0.9
    for subset in ((1,2,3,4),(0,1,2,3,4)):
        for fraction in (.01,.05):
            hs=miscalibrated_electric_clock(states,H,d,subset,fraction)
            y=evolution(K,hs,T,8)
            row={'type':'electric_clock_angle_error','subset':list(subset),'fractional_error':fraction,'raw':entanglement_fidelity([y],K),'recovered':entanglement_fidelity([y],K,V)}
            rows.append(row)
            print('ELECTRIC_ANGLE_CAL',json.dumps(row,sort_keys=True),flush=True)
        for tau in (.0005,.002,.01):
            y=lab_cubic_pulses(states,H,K,subset,8,tau)
            row={'type':'finite_electric_pulse','subset':list(subset),'pulse_tau':tau,'elapsed_time':T+32*tau,'per_addressed_link_control_amplitude':THETA/tau, 'integrated_abs_control_link_radians':32*THETA*len(subset),'raw':entanglement_fidelity([y],K),'recovered':entanglement_fidelity([y],K,V)}
            rows.append(row)
            print('FINITE',json.dumps(row,sort_keys=True),flush=True)
        for error in (.01,.05):
            hs=toggled_matrices(d,W,subset,{i:error for i in subset})
            y=evolution(K,hs,T,8)
            row={'type':'magnetic_drive_phase_error','subset':list(subset),'fractional_error':error,'raw':entanglement_fidelity([y],K),'recovered':entanglement_fidelity([y],K,V)}
            rows.append(row)
            print('CAL',json.dumps(row,sort_keys=True),flush=True)
        for ct in (.02,.05):
            hs=toggled_matrices(d,W,subset,crosstalk=ct)
            y=evolution(K,hs,T,8)
            row={'type':'drive_phase_crosstalk','subset':list(subset),'fraction_of_phase_on_unaddressed_links':ct,'raw':entanglement_fidelity([y],K),'recovered':entanglement_fidelity([y],K,V)}
            rows.append(row)
            print('XTALK',json.dumps(row,sort_keys=True),flush=True)
    for tau in (.002,.01):
        subset=(1,2,3,4)
        y=ramped_drive(states,d,W,K,subset,8,tau,substeps=16)
        row={'type':'finite_drive_phase_ramp','subset':list(subset),'phase_ramp_tau':tau,'elapsed_time':T+32*tau,'steps_per_ramp':16,'ramp_slew_max_rad_per_dimensionless_time':THETA/tau,'raw':entanglement_fidelity([y],K),'recovered':entanglement_fidelity([y],K,V)}
        rows.append(row)
        print('DRIVE_RAMP',json.dumps(row,sort_keys=True),flush=True)
    # One unconditioned local directed jump after 4 ideal cycles and before 4.
    for subset in ((1,2,3,4),(0,1,2,3,4)):
        hs=toggled_matrices(d,W,subset)
        first=evolution(K,hs,T/2,4)
        for site in (3,4):
            maps=[evolution(k@first,hs,T/2,4) for k in jump_operator(states,.1,site)]
            row={'type':'single_directed_wilson_jump','subset':list(subset),'site':site,'jump_strength_p':.1,'raw':entanglement_fidelity(maps,K),'recovered':entanglement_fidelity(maps,K,V)}
            rows.append(row)
            print('JUMP',json.dumps(row,sort_keys=True),flush=True)
    data={'dimension':DIM,'time':T,'g2':1.,'model':'minimal-link SU3 square chain Chen-Muller-Yao with constructed controls','rows':rows}
    with open('su3_native_control_compilation_results.json','w') as f:json.dump(data,f,sort_keys=True,indent=2)
    print('ALL EXACT CLOCK/SELECTIVE CONTROL/SOURCE-BASED CHANNEL CHECKS PASSED',flush=True)

if __name__=='__main__':main()
