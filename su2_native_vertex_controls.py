"""Source-motivated SU(2) four-mode matter/rishon vertex audit.

Four co-located color spin halves (one matter, three occupied rishon modes).
Tests native all-pair AFM exchange Gauss penalty, its symmetry-preserving
coupling defects, and an *assumed-tunable* positive-coupling exchange echo.
Exact 16-dimensional symbolic calculations plus 2D numerical channels.
No physical device, readout, fault-tolerant recovery, or Standard Model code.

Run: python su2_native_vertex_controls.py
Requires SymPy and NumPy (no SciPy).
"""
from itertools import combinations, product
import json
import numpy as np
import sympy as sp

PAIRS = tuple(combinations(range(4), 2))
I16, I2 = sp.eye(16), sp.eye(2)
X = sp.Matrix([[0, 1], [1, 0]])
Z = sp.diag(1, -1)


def swap(i, j):
    m = sp.zeros(16)
    for idx in range(16):
        bits = [(idx >> (3-k)) & 1 for k in range(4)]
        bits[i], bits[j] = bits[j], bits[i]
        dst = sum(b << (3-k) for k, b in enumerate(bits))
        m[dst, idx] = 1
    return m


SWAPS = {pair: swap(*pair) for pair in PAIRS}
H = {pair: (SWAPS[pair] - I16 / 2) / 2 for pair in PAIRS}
SUM_H = sum(H.values(), sp.zeros(16))
G2 = 3*I16 + 2*SUM_H


def valence_pairing(pairing):
    ket = sp.zeros(16, 1)
    for choices in product((0, 1), repeat=2):
        bits = [0]*4
        for (a, b), c in zip(pairing, choices):
            bits[a], bits[b] = c, 1-c
        index = sum(b << (3-i) for i, b in enumerate(bits))
        ket[index] = (-1)**sum(choices) / sp.Integer(2)
    return ket


A = valence_pairing(((0,1),(2,3)))
B = valence_pairing(((0,2),(1,3)))
Q = sp.Matrix.hstack(A, (2*B-A)/sp.sqrt(3))


def compress(op):
    return sp.simplify(Q.T.conjugate()*op*Q)


def arrays2(mat):
    return np.array(mat, dtype=complex)


def unitary(hermitian, time):
    vals, vecs = np.linalg.eigh(hermitian)
    return (vecs * np.exp(-1j*time*vals)) @ vecs.conj().T


def fidelity(branches):
    """Entanglement fidelity relative to identity for an equiprobable random unitary channel."""
    return float(np.mean([abs(np.trace(u)/2)**2 for u in branches]))


def exact_audit():
    assert Q.T.conjugate()*Q == I2
    assert G2*Q == sp.zeros(16,2)
    assert G2.eigenvals() == {sp.Integer(0):2,sp.Integer(2):9,sp.Integer(6):5}
    assert compress(SUM_H) == -sp.Rational(3,2)*I2

    # Independent pair interactions span six dimensions modulo identity.
    cols = [sp.Matrix(256,1,list(I16))] + [sp.Matrix(256,1,list(H[p])) for p in PAIRS]
    assert sp.Matrix.hstack(*cols).rank() == 7

    sigma_z = sp.diag(1,-1)
    sigma_plus = sp.Matrix([[0,1],[0,0]])
    eye = sp.eye(2)
    def local_sum(op):
        return sum((sp.kronecker_product(*[op if j==i else eye for j in range(4)]) for i in range(4)),sp.zeros(16))
    gz, gp = local_sum(sigma_z)/2, local_sum(sigma_plus)
    assert gz*Q == sp.zeros(16,2) and gp*Q == sp.zeros(16,2)
    for p in PAIRS:
        assert H[p]*gz==gz*H[p] and H[p]*gp==gp*H[p]

    # Native gauge protection needs a K4 complete all-pair exchange graph.
    # Uniform AFM J>0 yields exact singlet gap J, versus triplets.
    e0 = compress(SUM_H)
    assert e0 == -sp.Rational(3,2)*I2
    # The equal *star* AFM graph has a LOWER spin-one energy than singlets.
    star = H[(0,1)] + H[(0,2)] + H[(0,3)]
    assert compress(star) == -sp.Rational(3,4)*I2
    assert star.eigenvals()=={
        -sp.Rational(5,4):3, -sp.Rational(3,4):2,
        sp.Rational(1,4):6, sp.Rational(3,4):5}

    J = sp.symbols('J01 J02 J03 J12 J13 J23', real=True)
    full = sum((J[k]*H[PAIRS[k]] for k in range(6)),sp.zeros(16))
    sA = J[0] + J[5]  # 01,23
    sB = J[1] + J[4]  # 02,13
    sC = J[2] + J[3]  # 03,12
    predicted = -(sA+sB+sC)/4 * I2 + (-2*sA+sB+sC)/4 * Z + sp.sqrt(3)*(sC-sB)/4*X
    assert compress(full)==predicted
    assert compress(H[(0,1)])==sp.diag(-sp.Rational(3,4),sp.Rational(1,4))
    assert compress(H[(0,2)])-compress(H[(0,3)]) == -sp.sqrt(3)*X/2

    # Gauge-compatible logical X from the difference of two co-located exchanges.
    HX = 2*(H[(0,3)] - H[(0,2)])/sp.sqrt(3)
    assert compress(HX)==X
    assert HX*Q==Q*X  # no physical code leakage for ideal pulsing
    # Positive-coupling implementation: shift every coupling by common Jbase.
    # J02=Jbase-2Omega/sqrt3, J03=Jbase+2Omega/sqrt3.
    Omega, Jbase=sp.symbols('Omega Jbase',positive=True,real=True)
    Hpulse = Jbase*SUM_H + Omega*HX
    assert compress(Hpulse) == -sp.Rational(3,2)*Jbase*I2 + Omega*X
    assert sp.simplify(Hpulse*Q-Q*compress(Hpulse)) == sp.zeros(16,2)

    # Gauge charge monitoring sees no error even with a single-bond mismatch.
    assert compress(G2*H[(0,1)])==sp.zeros(2)
    print('Exact: K4 uniform pair exchange = (G^2-3I)/2, gap J; star has spin-one below singlet')
    print('Exact: 6 independent pair couplings; singlet degeneracy iff three complementary pair sums equal')
    print('Exact: single-bond detuning acts as logical Z and stays gauge singlet; symmetric X exchange control exists')
    return {'native_spin0':2,'native_spin1':9,'native_spin2':5,'star_lowest_spin':1,
            'Jgap_coeff':1,'delta_splitting_coeff':1}


def numerical_audit():
    # Scalar gauge-protection Hamiltonian on spin-zero, dropped as global phase.
    # An unknown +/- delta shift to h01 induces a logical Z detuning delta/2.
    hnoise = arrays2(compress(H[(0,1)]))
    h_other = arrays2(compress(H[(0,2)]))
    hx = arrays2(X)
    hold_t=3.0
    Omega=5.0
    tau=np.pi/(2*Omega)
    rows=[]
    for delta in [0.0,0.02,0.05,0.1,0.2]:
        bare=[unitary(sign*delta*hnoise,hold_t) for sign in (-1,1)]
        fe_bare=fidelity(bare)
        predicted=(1+np.cos(delta*hold_t))/2
        assert abs(fe_bare-predicted)<1e-13
        # Ideal instantaneous X echo acts as logical identity up to phase.
        perfect=[]
        for sign in (-1,1):
            w=unitary(sign*delta*hnoise,hold_t/2)
            perfect.append(hx@w@hx@w)
        assert abs(fidelity(perfect)-1)<1e-12
        # Two FINITE exchange pulses, both with the unknown drift present.
        static=[]
        for sign in (-1,1):
            w=unitary(sign*delta*hnoise,hold_t/2)
            pulse=unitary(Omega*hx+sign*delta*hnoise,tau)
            static.append(pulse@w@pulse@w)
        fe_echo=fidelity(static)
        # 4 equally likely noise scenarios: detuning sign may flip between halves.
        telegraph=[]
        for s1,s2 in product((-1,1),repeat=2):
            w1=unitary(s1*delta*hnoise,hold_t/2)
            w2=unitary(s2*delta*hnoise,hold_t/2)
            p1=unitary(Omega*hx+s1*delta*hnoise,tau)
            p2=unitary(Omega*hx+s2*delta*hnoise,tau)
            telegraph.append(p2@w2@p1@w1)
        fe_telegraph=fidelity(telegraph)
        # Independent h02 coupling drift contains an X component, which an
        # X echo cannot reverse; this shows the operation is not QEC.
        other_branches=[]
        for sign in (-1,1):
            w=unitary(sign*delta*h_other,hold_t/2)
            pulse=unitary(Omega*hx+sign*delta*h_other,tau)
            other_branches.append(pulse@w@pulse@w)
        fe_other=fidelity(other_branches)
        assert fe_echo<=1+1e-12 and fe_telegraph<=1+1e-12
        # Physical (nonnegative) simultaneous control strengths for Omega=5,
        # baseline Jbase=12, max 17.7735; test control assumptions explicitly.
        baseline=12.0
        coefficient_shift=2*Omega/np.sqrt(3)
        assert baseline-coefficient_shift>0
        row={'delta':delta,'time':hold_t,'F_bare':fe_bare,
             'F_finite_echo_static':fe_echo,'F_finite_echo_sign_flip':fe_telegraph,
             'F_finite_echo_other_bond':fe_other,
             'pulse_duration':tau,'total_with_echo':hold_t+2*tau}
        rows.append(row)
        print('delta=%.3f bare=%.9f static_echo=%.9f flip_echo=%.9f other_bond=%.9f' % (
            delta,fe_bare,fe_echo,fe_telegraph,fe_other))
    # A FULL 16-component test of the 2D pulse evolution vs its block form.
    HX = 2*(H[(0,3)] - H[(0,2)])/sp.sqrt(3)
    for sign in (-1,1):
        delta=0.1
        full_h=np.array(12*SUM_H + Omega*HX +sign*delta*H[(0,1)],dtype=complex)
        full_v=unitary(full_h,tau) @ np.array(Q,dtype=complex)
        logical_v=np.array(Q,dtype=complex) @ unitary(np.array(compress(12*SUM_H+Omega*HX+sign*delta*H[(0,1)]),dtype=complex),tau)
        assert np.linalg.norm(full_v-logical_v)<3e-12
    print('Full 16D vs compressed 2D finite-pulse evolution PASS')
    return rows


if __name__=='__main__':
    summary=exact_audit()
    rows=numerical_audit()
    print('ALL EXACT/NATIVE CONTROL AUDITS PASS')
    print('JSON_RESULTS='+json.dumps({'summary':summary,'rows':rows},sort_keys=True))
