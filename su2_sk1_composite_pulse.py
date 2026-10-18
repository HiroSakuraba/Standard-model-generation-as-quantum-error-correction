"""SK1-inspired SU(2)-invariant composite exchange gate: full-spin benchmarks.

All constituent microscopic spin operators are co-located SU(2)-scalar
Heisenberg exchanges. The outer reference code is the product of two
four-spin singlet qubits in the 8-spin mediator from PR #15. Physical
cross-block sin^2 pulses use the full 8-spin Hamiltonian, not the
second-order approximation. Intra-block X/Z exchange rotations are
ideal, instantaneous *SU(2)-invariant* local control operations. Their
finite durations and elementary pulse faults remain unmodeled.

For each symmetric echo E(t)=L U(t) L U(t), with L=R_XA(pi)R_XD(pi)
and U(t) a simultaneous sin^2 cross-exchange pulse, calibrate times:
 E(t_target): nonlocal phase = pi   (CZ class).
 E(t_corr):   nonlocal phase = 4*pi (nominal 2*pi ZZ rotation).
These phases are the UNWRAPPED diagonal phase invariant, not modulo 2pi.

SK1-like correction rotates ZZ about the DATA logical Y axis:
  phi=acos(-1/8)
  U_SK1 = RY(-phi) E_corr RY(phi) *
          RY(phi) E_corr RY(-phi) * E_target .
The tilted noncommuting exchange axes cancel the linear amplitude error
of the logical entangling phase; the second-order correction is tested
with both physical cross couplers faulted separately and together.

We compare faulty unnormalized projected 4x4 matrices to the *fixed*
nominal output, allowing a GLOBAL phase alignment only. No per-fault
single-logical-Z recalibration is allowed. This does NOT imply the
finite-time gate is exactly Clifford, and is NOT a hardware threshold.

Critical cost: six smooth cross-exchange envelopes, total active
duration ~344/Delta (~9x the PR #16 38.23/Delta gate) plus local
exchange rotations whose finite duration is not included.

Run: OPENBLAS_NUM_THREADS=1 python su2_sk1_composite_pulse.py
"""
from functools import lru_cache
import math

import numpy as np
from scipy.linalg import eigh
from scipy.optimize import brentq
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import expm_multiply

from su2_shaped_pulse import setup, metric
from two_body_su2_mediator import swap_matrix

I4 = np.eye(4)
SIGN_ZZ = np.array([1, -1, -1, 1])
CAL_BINS = 512
AUDIT_BINS = 1024


def physical_local_generators():
    """Explicit, SU(2)-invariant, same-vertex TWO-SPIN exchanges."""
    eye = np.eye(256)
    def heis(i, j):
        return (swap_matrix(8, i, j)-0.5*eye)/2
    xa = 2*(heis(0, 3)-heis(0, 2))/np.sqrt(3)
    xd = 2*(heis(4, 7)-heis(4, 6))/np.sqrt(3)
    za = -2*heis(0, 1)-0.5*eye
    zd = -2*heis(4, 5)-0.5*eye
    return {"Xa": xa, "Xd": xd, "Za": za, "Zd": zd}


class PhysicalComposite:
    def __init__(self):
        self.h0f, self.h1f, self.h2f, self.embed, self.code, self.code_full = setup()
        self.h0, self.h1, self.h2 = [
            self.embed.T @ h @ self.embed
            for h in (self.h0f, self.h1f, self.h2f)
        ]
        self.locals_full = physical_local_generators()
        self.locals = {
            name:self.embed.T @ matrix @ self.embed
            for name,matrix in self.locals_full.items()
        }
        X = np.array([[0, 1], [1, 0]])
        Z = np.diag([1, -1])
        for name,expected in {
            "Xa":np.kron(X,np.eye(2)), "Xd":np.kron(np.eye(2),X),
            "Za":np.kron(Z,np.eye(2)), "Zd":np.kron(np.eye(2),Z)
        }.items():
            A = self.locals[name]
            assert np.linalg.norm(self.code.T @ A @ self.code - expected) < 1e-10
            assert np.linalg.norm(
                self.locals_full[name] @ self.embed -
                self.embed @ A) < 1e-10
            # No local pulse can leave either product-of-singlets code.
            assert np.linalg.norm(
                A@self.code-self.code@expected) < 1e-10
        assert np.linalg.norm(self.h0@self.code) < 1e-10
        self._diagonal = {}
        self._rotation = {}

    def rotation(self, name, angle, local_error=0.0):
        """An intra-block exchange rotation with coherent amplitude error."""
        key = (name, float(angle), float(local_error))
        if key not in self._rotation:
            values, vectors = eigh(self.locals[name])
            self._rotation[key] = (
                vectors * np.exp(-0.5j*angle*(1+local_error)*values)) @ vectors.T
        return self._rotation[key]

    def ry(self, angle, local_error=0.0):
        """R_Y on data, synthesized from separately faulty X/Z exchanges."""
        return (self.rotation("Zd", np.pi/2,local_error) @
                self.rotation("Xd", angle,local_error) @
                self.rotation("Zd", -np.pi/2,local_error))

    def echo_flip(self, local_error=0.0):
        return (self.rotation("Xa",np.pi,local_error) @
                self.rotation("Xd",np.pi,local_error))

    def eigensteps(self, fault1, fault2, bins):
        key = (float(fault1), float(fault2), int(bins))
        if key not in self._diagonal:
            steps = []
            for i in range(bins):
                env = np.sin(np.pi*(i+0.5)/bins)**2
                energies, vectors = eigh(
                    self.h0 + env*((1+fault1)*self.h1 +
                                   (1+fault2)*self.h2))
                steps.append((energies, vectors))
            self._diagonal[key] = steps
        return self._diagonal[key]

    def smooth(self, duration, fault1=0.0, fault2=0.0, bins=CAL_BINS):
        """Actual time-ordered 14D singlet-sector Hamiltonian evolution."""
        U = np.eye(14, dtype=complex)
        for energies, vectors in self.eigensteps(fault1, fault2, bins):
            U = (vectors * np.exp(-1j*duration*energies/bins)) @ vectors.T @ U
        return U

    def echo(self, duration, fault1=0.0, fault2=0.0, bins=CAL_BINS,
             local_error=0.0):
        pulse = self.smooth(duration, fault1, fault2, bins)
        flip = self.echo_flip(local_error)
        return flip @ pulse @ flip @ pulse

    def nonlocal_phase(self, U):
        W = self.code.conj().T @ U @ self.code
        phi = np.angle(np.diag(W))
        return float(np.angle(np.exp(1j*np.dot(phi, SIGN_ZZ))))

    def unwrapped_phase(self, duration, bins=CAL_BINS):
        phase = self.nonlocal_phase(self.echo(duration,bins=bins))
        # On the specified 18..79/Delta calibration range the observed
        # nonlocal phase winds approximately at rate 0.164/Delta.
        # Choose that winding branch explicitly; do not mistake
        # a 4*pi correction for a trivial zero-angle pulse.
        turns = round((0.164*duration - phase)/(2*np.pi))
        return phase+2*np.pi*turns

    def calibrate(self):
        target = brentq(
            lambda t: self.unwrapped_phase(t)-np.pi, 18., 21., xtol=1e-10)
        correct = brentq(
            lambda t: self.unwrapped_phase(t)-4*np.pi, 74., 79., xtol=1e-10)
        assert abs(self.unwrapped_phase(target)-np.pi) < 1e-9
        assert abs(self.unwrapped_phase(correct)-4*np.pi) < 1e-9
        return float(target), float(correct)

    def composite(self, t_target, t_correct,
                  fault1=0.0, fault2=0.0, bins=CAL_BINS,
                  local_error=0.0):
        phi = math.acos(-1/8)
        plus = self.ry(phi,local_error)
        minus = self.ry(-phi,local_error)
        target = self.echo(t_target,fault1,fault2,bins,local_error)
        correct = self.echo(t_correct,fault1,fault2,bins,local_error)
        # Two noncommuting tilted correction axes generate SK1-style
        # first-order amplitude cancellation, not a commuting echo.
        return (minus @ correct @ minus.conj().T @
                plus @ correct @ plus.conj().T @ target)

    def code_metrics(self, U):
        W = self.code.conj().T @ U @ self.code
        return metric(U,self.code), W


def aligned_error(code,U_fault,U_nom):
    """Nominal fixed calibration; NO fault-dependent Z-phase correction."""
    ref = code.conj().T@U_nom@code
    out = code.conj().T@U_fault@code
    phase = np.exp(-1j*np.angle(np.trace(ref.conj().T@out)))
    deviation = float(np.linalg.norm(phase*out-ref,ord=2))
    leakage = float(max(
        0., 1.-np.linalg.eigvalsh(out.conj().T@out).min()))
    return deviation, leakage


def independent_full_spin_check(model,t_target,t_correct):
    """Actual 256-state pulse propagation vs symmetry-reduced 14D.

    Uses separate sparse expm_multiply on all four physical codewords,
    including intra-register exchange rotations. Same 12-bin midpoint
    waveform in both representations; numerical agreement is not a
    continuum-time gate accuracy certificate.
    """
    bins=12
    psi=model.code_full.astype(complex).copy()

    def full_smooth(psi,duration):
        for j in range(bins):
            env=np.sin(np.pi*(j+.5)/bins)**2
            H=model.h0f+env*(model.h1f+model.h2f)
            X=csr_matrix(-1j*(duration/bins)*H)
            psi=expm_multiply(X,psi,traceA=-1j*(duration/bins)*np.trace(H))
        return psi

    full_cache={}
    def full_rotate(psi,name,angle):
        key=(name,float(angle))
        if key not in full_cache:
            vals,vecs=eigh(model.locals_full[name])
            full_cache[key]=(vecs*np.exp(-.5j*angle*vals))@vecs.T
        return full_cache[key]@psi

    def full_flip(psi):
        return full_rotate(full_rotate(psi,"Xd",np.pi),"Xa",np.pi)

    def full_echo(psi,duration):
        return full_flip(full_smooth(full_flip(full_smooth(psi,duration)),duration))

    def full_ry(psi,angle):
        psi=full_rotate(psi,"Zd",-np.pi/2)
        psi=full_rotate(psi,"Xd",angle)
        return full_rotate(psi,"Zd",np.pi/2)

    angle=math.acos(-1/8)
    psi=full_echo(psi,t_target)
    psi=full_ry(psi,-angle)
    psi=full_echo(psi,t_correct)
    psi=full_ry(psi,angle)
    psi=full_ry(psi,angle)
    psi=full_echo(psi,t_correct)
    psi=full_ry(psi,-angle)

    reduced=model.composite(t_target,t_correct,bins=bins)
    from_full=model.code_full.conj().T@psi
    from_reduced=model.code.conj().T@reduced@model.code
    difference=float(np.linalg.norm(from_full-from_reduced))
    assert difference < 1e-8, difference
    print("Independent FULL 256D vs 14D composite pulse check:",difference)
    return difference


def test():
    model=PhysicalComposite()
    t,c=model.calibrate()
    assert 19<t<20 and 76<c<77
    gate_time=2*(t+2*c)
    print(f"Pulse calibration: target t={t:.9f}, "
          f"2pi correction t={c:.9f}, total cross-exchange duration "
          f"{gate_time:.6f}/Delta plus uncounted local rotations")
    assert 340<gate_time<350

    nominal_single=model.echo(t)
    nominal_composite=model.composite(t,c)
    nominal_metric, W=model.code_metrics(nominal_composite)
    print("Nominal SK1-class gate:",nominal_metric)
    assert nominal_metric["deviation"]<1e-4
    assert nominal_metric["leakage"]<5e-6

    comparisons={}
    for name,indices in (
        ("first",(1,0)),("second",(0,1)),("both",(1,1))
    ):
        for error in (-.10,-.05,-.01,.01,.05,.10):
            e1,e2=error*indices[0],error*indices[1]
            bare=aligned_error(model.code,
                               model.echo(t,e1,e2),nominal_single)
            composite=aligned_error(model.code,
                model.composite(t,c,e1,e2),nominal_composite)
            comparisons[(name,error)]=(bare,composite)
            print(f"Fault {name} {error:+.0%}:"
                  f" single_dev={bare[0]:.8f}, "
                  f"SK1_dev={composite[0]:.8f}, "
                  f"SK1_leak={composite[1]:.9g}")
            assert composite[0]<bare[0]
            assert composite[1]<5e-6
            if name!="both":
                if abs(error)<=.01:
                    assert composite[0]<.001
                    assert composite[0]<bare[0]/10
                elif abs(error)<=.05:
                    assert composite[0]<.012
                else:
                    assert composite[0]<.04

    # Crucial negative control: the composite is robust to the
    # CROSS couplers above, not necessarily to calibration of its
    # additional intra-block exchange rotations.
    for e in (-.05,-.01,.01,.05):
        faulty = model.composite(t,c,local_error=e)
        err,leak = aligned_error(model.code,faulty,nominal_composite)
        print(f"All intra-block exchange rotations miscalibrated {e:+.0%}:",
              {"deviation":err,"leakage":leak})
        assert np.isfinite(err) and np.isfinite(leak)
        assert 0<=leak<=1

    # Check temporal-discretization stability INDEPENDENT of the grid
    # used for duration calibration (512). Low-bin interference can
    # create spurious extraordinarily small errors, so do not trust it.
    for bins in (256,512,AUDIT_BINS):
        U=model.composite(t,c,.1,0,bins)
        nominal=model.composite(t,c,0,0,bins)
        dev,leak=aligned_error(model.code,U,nominal)
        print(f"Composite +10% first-coupler, {bins} bins:",
              {"deviation":dev,"leakage":leak})
        assert .029<dev<.032
        assert leak<5e-6
    independently=independent_full_spin_check(model,t,c)
    assert independently<1e-8
    print("ALL SK1-COMPOSITE SU2 EXCHANGE TESTS PASSED")


if __name__=="__main__":
    test()
