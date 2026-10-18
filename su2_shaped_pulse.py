"""Numerically optimized smooth two-body SU(2) exchange gate with leakage audit.

Previously, two four-spin singlet qubits, using a block Casimir penalty,
mediated a perturbative ZZ gate from pairwise SU(2)-scalar exchanges.
This script improves its finite pulse while accounting for the complete
8-spin Hamiltonian and coherent one-coupler calibration errors.

Hamiltonian, in units Delta=1:
  H(t) = S_A^2 + S_D^2 + J1(t) S_A0.S_D0 + J2(t) S_A1.S_D1
  J1(t)=J2(t)=Jpeak*sin(pi*t/T)^2, 0 <= t <= T.

Optimize T alone for Jpeak=1 by bounded scalar optimization. Use an exact
numerical SU(2) total-spin-zero invariant-sector reduction from 256 to 14
states; check the full 256x256 spin matrices preserve that sector, and
independently compare a direct 256-state pulse propagation to 14-state
propagation. The four-dimensional code consists of a product of two
four-spin singlets.

Gate error: unnormalized 4x4 projected unitary minus a specified CZ
representative after global/individual logical Z-phase compensation.
This is an OPERATOR-NORM deviation, NOT average infidelity. Leakage is
the maximal probability for any code input to be outside the product
of two 4-spin singlet logical blocks at pulse end. Gauge-invariant
leakage out of the chosen code is NOT a Gauss law violation.

The optimized pulse is not immune to calibration faults, and no physical
native device, broadband controls, stochastic fault rate, fault-tolerant
multi-module syndrome extraction, or Standard Model mechanism is proven.

Requires numpy, scipy and sympy, plus two_body_su2_mediator.py.
Run: OPENBLAS_NUM_THREADS=1 python su2_shaped_pulse.py
"""
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import minimize_scalar
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import expm_multiply

from two_body_su2_mediator import basis4, swap_matrix

DELTA = 1.0
PEAK = 1.0
SEARCH_BINS = 512
ID4 = np.eye(4)
SIGN_ZZ = np.array([1, -1, -1, 1])


def setup():
    """Construct the full eight-spin Hamiltonian and its full S=0 block.

    The global SU(2) spin-zero multiplicity of 8 spin halves is 14.
    Build the 70-dimensional total-Sz=0 sector, diagonalize S^2, and
    embed its exact (to eigensolver precision) 14D kernel into 256D.
    """
    dim = 256
    eye = np.eye(dim)
    exchanges = {
        (i, j): (swap_matrix(8, i, j) - 0.5 * eye) / 2
        for i in range(8) for j in range(i+1, 8)
    }
    h0 = 6 * eye.copy()
    for block in (range(4), range(4, 8)):
        for i in block:
            for j in block:
                if i < j:
                    h0 += 2 * exchanges[(i, j)]
    h1, h2 = exchanges[(0, 4)], exchanges[(1, 5)]

    inner = basis4()
    code = np.zeros((dim, 4))
    for a in range(2):
        for b in range(2):
            for spin_a, c_a in inner[a].items():
                for spin_b, c_b in inner[b].items():
                    bits = spin_a + spin_b
                    index = sum(bit << (7-k) for k, bit in enumerate(bits))
                    code[index, 2*a+b] = float(c_a*c_b)
    assert np.linalg.norm(code.T @ code - ID4) < 1e-11

    sz0 = np.array([i for i in range(256) if i.bit_count() == 4])
    assert len(sz0) == 70
    total_s2 = 6*eye.copy()
    for exchange in exchanges.values():
        total_s2 += 2 * exchange
    eig, vec = eigh(total_s2[np.ix_(sz0, sz0)])
    assert np.max(np.abs(eig[:14])) < 1e-11
    assert eig[14] > 1.9
    embed = np.zeros((dim, 14))
    embed[sz0, :] = vec[:, :14]
    assert np.linalg.norm(embed.T @ embed - np.eye(14)) < 1e-11

    # These are the actual 256-state SU(2)-scalar Hamiltonian terms.
    # They never take a global color singlet into non-singlet states.
    for full_h in (h0, h1, h2):
        assert np.linalg.norm(full_h @ embed -
                              embed @ (embed.T @ full_h @ embed)) < 1e-10
    reduced_code = embed.T @ code
    assert np.linalg.norm(reduced_code.T @ reduced_code-ID4) < 1e-11
    assert np.linalg.norm(code - embed @ reduced_code) < 1e-10
    return h0, h1, h2, embed, reduced_code, code


def pulse(h0, h1, h2, peak, duration, bins,
          first_coupler_fault=0.0, profile="sin2"):
    """Time-ordered midpoint piecewise-constant evolution in S=0.

    The faults change one physical exchange strength. They are not
    assumed equivalent to a final Pauli error and preserve SU(2).
    """
    U = np.eye(14, dtype=complex)
    dt = duration / bins
    for i in range(bins):
        u = (i+0.5)/bins
        envelope = np.sin(np.pi*u)**2 if profile == "sin2" else 1.0
        coupling = peak*envelope
        h = h0+coupling*((1+first_coupler_fault)*h1+h2)
        energies, vectors = eigh(h)
        U = (vectors*np.exp(-1j*energies*dt)) @ vectors.T @ U
    return U


def metric(U, code):
    """Candidate phase-corrected CZ representative and exact leakage.

    Removing global/local logical Z phases means only the invariant
    nonlocal phase phi00-phi01-phi10+phi11 matters. We distribute its
    mismatch evenly over four diagonal phases to specify a candidate
    calibrated target; global optimality is NOT claimed.
    """
    W = code.conj().T @ U @ code
    phi = np.angle(np.diag(W))
    theta = np.angle(np.exp(1j*np.dot(SIGN_ZZ, phi)))
    mismatch = np.angle(np.exp(1j*(theta-np.pi)))
    target = np.diag(np.exp(1j*(phi-SIGN_ZZ*mismatch/4)))
    deviation = float(np.linalg.norm(W-target, ord=2))
    leakage = float(max(0, 1-np.linalg.eigvalsh(W.conj().T @ W).min()))
    offdiag = float(np.linalg.norm(W-np.diag(np.diag(W))))
    return {
        "deviation": deviation, "leakage": leakage,
        "entangling_phase_error": abs(float(mismatch)),
        "offdiag": offdiag,
    }


def test():
    h0_full, h1_full, h2_full, embed, code, physical_code = setup()
    h0, h1, h2 = [
        embed.T @ h @ embed for h in (h0_full, h1_full, h2_full)
    ]

    # Original PR #15 square-pulse baseline, now allowing calibrated
    # global/local Z phases (unlike its raw operator norm difference).
    nominal_time = 6*np.pi/.4**2
    baseline = metric(pulse(h0, h1, h2, .4, nominal_time, 1,
                            profile="square"), code)
    # Same PEAK exchange amplitude comparison, without retiming.
    peak_square = metric(pulse(h0, h1, h2, PEAK, 6*np.pi, 1,
                               profile="square"), code)

    objective = lambda t: metric(
        pulse(h0, h1, h2, PEAK, t, SEARCH_BINS), code)["deviation"]
    opt = minimize_scalar(
        objective, bounds=(30.0, 46.0), method="bounded",
        options={"xatol": 1e-7, "maxiter": 120},
    )
    assert opt.success
    duration = float(opt.x)
    results = {
        bins: metric(pulse(h0, h1, h2, PEAK, duration, bins), code)
        for bins in (256, 512, 1024)
    }
    print("PR15 square baseline J/Delta=0.4:", nominal_time, baseline)
    print("Square at same Jpeak/Delta=1:", 6*np.pi, peak_square)
    print("Local 1D optimized sin^2 pulse Jpeak/Delta=1, "
          f"T Delta={duration:.9f}")
    for bins, row in results.items():
        print(f"Midpoint bins={bins}:", row)

    accurate = results[1024]
    assert 30 < duration < 46
    assert duration < nominal_time / 2
    assert accurate["deviation"] < 1e-4
    assert accurate["leakage"] < 1e-6
    assert accurate["offdiag"] < 1e-9
    assert abs(results[1024]["leakage"] -
               results[512]["leakage"]) < 1e-7
    assert abs(results[1024]["deviation"] -
               results[512]["deviation"]) < 1e-4
    assert accurate["deviation"] < baseline["deviation"]/100

    # A pulse can have extremely low nominal leakage while retaining
    # first-order sensitivity to a GAUGE-INVARIANT miscalibrated exchange.
    for fault in (-.10, -.05, .05, .10):
        row = metric(pulse(h0, h1, h2, PEAK, duration, 512,
                           first_coupler_fault=fault), code)
        print(f"One elementary coupler miscalibrated {fault:+.0%}:", row)
        if abs(fault) > .09:
            assert row["deviation"] > .05

    # Independent, non-reduced 256D sparse exponential evolution
    # versus the projected total-singlet propagation. These must agree
    # for any gauge-invariant time-dependent Hamiltonian.
    bins = 24
    dt = duration / bins
    psi = physical_code.astype(complex).copy()
    for k in range(bins):
        env = np.sin(np.pi*(k+.5)/bins)**2
        H = h0_full+PEAK*env*(h1_full+h2_full)
        psi = expm_multiply(csr_matrix(-1j*dt*H), psi,
                            traceA=-1j*dt*np.trace(H))
    full_projection = physical_code.T @ psi
    reduced_projection = code.T @ pulse(
        h0, h1, h2, PEAK, duration, bins) @ code
    error = float(np.linalg.norm(full_projection-reduced_projection))
    print("Independent 256D vs 14D full-pulse check:", error)
    assert error < 1e-9
    print("ALL SMOOTH SU2 PULSE OPTIMIZATION CHECKS PASSED")


if __name__ == "__main__":
    test()
