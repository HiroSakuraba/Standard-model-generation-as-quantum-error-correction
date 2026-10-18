"""Exchange-only SU(2) two-local mediator: exact theory, full spin numerics.

A color-neutral two-level ancilla coupled to ONE fundamental SU(2) spin
has no gauge-invariant entangling two-body Hamiltonian (Schur lemma).
We instead encode the ancilla in the TWO spin-zero states of FOUR SU(2)
fundamentals at the same gauge vertex as the rich four-spin matter/rishon
singlet block. Both blocks separately carry total spin zero.

Allowed elementary interactions: tunable Heisenberg pair exchanges
S_i.S_j (exact SU(2) scalars) and a block-spin Casimir penalty.
H=Delta*(S_A^2+S_D^2)+J1*S_A0.S_D0+J2*S_A1.S_D1.

The cross coupling maps singlet(A)*singlet(D) exclusively to spin-1(A)
* spin-1(D), energy 4Delta. PVP=0, and exact second-order reduction:
H_eff= -P V^2 P / (4Delta);
coefficient of Z_A Z_D = -J1*J2 / (24 Delta).
This is an engineered perturbative entangler, not an EXACT controlled
gate. Compare against the FULL 256-state (eight-spin) Hamiltonian to
measure leakage and coherent deviation after nominal pi/4 ZZ phase.

Native engineering, finite-pulse errors, gauge-sector preparation,
fault-tolerant ancilla readout and Standard Model realization remain open.
"""
from itertools import product
import numpy as np
import sympy as sp

I2 = sp.eye(2)
Z2 = sp.diag(1, -1)
I4 = sp.eye(4)


def singlet_pair_state(pairing):
    state = {}
    for bits in product((0, 1), repeat=2):
        spins = [0]*4
        for (i, j), bit in zip(pairing, bits):
            spins[i], spins[j] = bit, 1-bit
        state[tuple(spins)] = sp.Rational((-1)**sum(bits), 2)
    return state


def overlap(left, right):
    return sp.simplify(sum(
        sp.conjugate(coeff)*right.get(bits, 0)
        for bits, coeff in left.items()))


def combine(left, right, a, b):
    return {k: sp.simplify(a*left.get(k, 0)+b*right.get(k, 0))
            for k in left.keys() | right.keys()
            if sp.simplify(a*left.get(k, 0)+b*right.get(k, 0)) != 0}


def basis4():
    a = singlet_pair_state(((0, 1), (2, 3)))
    b = singlet_pair_state(((0, 2), (1, 3)))
    assert overlap(a, b) == sp.Rational(1, 2)
    q1 = combine(a, b, -1/sp.sqrt(3), 2/sp.sqrt(3))
    assert overlap(a, q1) == 0 and overlap(q1, q1) == 1
    return a, q1


def tensor4(a, b):
    return {x+y: sp.simplify(c*d)
            for x, c in a.items() for y, d in b.items()}


def heisenberg(state, i, j):
    """Physical pair interaction S_i.S_j=(SWAP_ij - I/2)/2."""
    out = {}
    for bits, coeff in state.items():
        moved = list(bits)
        moved[i], moved[j] = moved[j], moved[i]
        moved = tuple(moved)
        out[moved] = out.get(moved, 0) + coeff/2
        out[bits] = out.get(bits, 0) - coeff/4
    return {k: sp.simplify(v) for k, v in out.items()
            if sp.simplify(v) != 0}


def add(a, b):
    return {k: sp.simplify(a.get(k, 0)+b.get(k, 0))
            for k in a.keys() | b.keys()
            if sp.simplify(a.get(k, 0)+b.get(k, 0)) != 0}


def neutral_ancilla_no_go():
    """All SU(2) invariants on (neutral qubit) x (one color doublet).

    Commutator constraints with 3 color Pauli generators on the complete
    two-body operator basis have rank 12 / 16: invariant commutant dim 4.
    The four surviving operators are ancilla-only, tensored with I_color.
    """
    basis = (sp.eye(2), sp.Matrix([[0,1],[1,0]]),
             sp.Matrix([[0,-sp.I],[sp.I,0]]), sp.diag(1,-1))
    two_body = [sp.kronecker_product(a, b) for a in basis for b in basis]
    constraints = []
    for generator in basis[1:]:
        G = sp.kronecker_product(sp.eye(2), generator)
        constraints.append(sp.Matrix.hstack(*[
            sp.Matrix(16,1,list(m*G-G*m)) for m in two_body]))
    commutators = sp.Matrix.vstack(*constraints)
    assert commutators.rank() == 12
    assert len(commutators.nullspace()) == 4
    print("Neutral ancilla + one color fundamental: 2-local invariant "
          "commutant dimension 4; no entangling interaction PASS")


def gauss_action(state, component):
    result = {}
    for bits, coeff in state.items():
        if component == "2z":
            eigen = sum(1-2*v for v in bits)
            if eigen:
                result[bits] = result.get(bits,0)+eigen*coeff
        else:
            for k in range(8):
                if bits[k] == 1:
                    moved = list(bits)
                    moved[k] = 0
                    moved = tuple(moved)
                    result[moved] = result.get(moved,0)+coeff
    return {k: sp.simplify(v) for k,v in result.items()
            if sp.simplify(v) != 0}


def code_and_second_order():
    code = tuple(tensor4(a,b) for a in basis4() for b in basis4())
    assert all(overlap(code[i],code[j]) == int(i==j)
               for i in range(4) for j in range(4))
    for bits in product((0,1),repeat=8):
        state = {bits: sp.Integer(1)}
        for i,j in ((0,4),(1,5)):
            for component in ("2z","plus"):
                lhs = gauss_action(heisenberg(state,i,j),component)
                rhs = heisenberg(gauss_action(state,component),i,j)
                assert not add(lhs,{k:-v for k,v in rhs.items()})
    print("All 2048 SU(2) Gauss commutators on eight-spin basis PASS")

    e1 = [heisenberg(psi,0,4) for psi in code]
    e2 = [heisenberg(psi,1,5) for psi in code]
    assert all(overlap(code[i],e) == 0 for i in range(4)
               for e in e1+e2)
    w1 = sp.Matrix(4,4,lambda i,j:overlap(e1[i],e1[j]))
    w2 = sp.Matrix(4,4,lambda i,j:overlap(e2[i],e2[j]))
    cross = sp.Matrix(4,4,lambda i,j:
                      overlap(e1[i],e2[j])+overlap(e2[i],e1[j]))
    assert w1 == w2 == sp.Rational(3,16)*I4
    pair = -sp.Rational(1,4)*I2-sp.Rational(1,2)*Z2
    assert cross == sp.Rational(2,3)*sp.kronecker_product(pair,pair)
    ZZ = sp.kronecker_product(Z2,Z2)
    J1,J2,Delta = sp.symbols("J1 J2 Delta",positive=True,real=True)
    heff = -(J1*J1*w1+J2*J2*w2+J1*J2*cross)/(4*Delta)
    assert sp.simplify(sp.trace(heff*ZZ)/4+J1*J2/(24*Delta)) == 0
    equal = sp.simplify(heff.subs({J1:1,J2:1,Delta:1}))
    Za = sp.kronecker_product(Z2,I2)
    Zd = sp.kronecker_product(I2,Z2)
    assert equal == -(sp.Rational(5,12)*I4+
                      sp.Rational(1,12)*(Za+Zd)+
                      sp.Rational(1,6)*ZZ)/4
    print("Exact Schrieffer-Wolff order-two: ZZ coefficient "
          "-J1*J2/(24 Delta), nonzero entangler PASS")
    return code,equal


def swap_matrix(n,i,j):
    dim=1<<n
    ids = np.arange(dim)
    bi=(ids>>(n-1-i))&1
    bj=(ids>>(n-1-j))&1
    moved=ids^((bi^bj)<<(n-1-i))^((bi^bj)<<(n-1-j))
    return np.eye(dim,dtype=float)[moved]


def full_eight_spin_dynamics(code,heff):
    """Exact diagonalization of H=Delta*Casimir + J*two cross Heisenbergs.

    Numerically compare full 256d unitary to the predicted order-two
    4d unitary for nominal ZZ angle pi/4 at J/Delta=0.4, 0.2, 0.1.
    Report operator norm error without correcting global/local phases
    and maximal code-subspace leakage probability (not infidelity).
    """
    dim=256
    identity = np.eye(dim)
    swaps={(i,j):(swap_matrix(8,i,j)-.5*identity)/2
           for i in range(8) for j in range(i+1,8)}
    penalty = 6*identity.copy()
    for block in (range(4),range(4,8)):
        for i in block:
            for j in block:
                if j>i:
                    penalty += 2*swaps[(i,j)]
    cross = swaps[(0,4)]+swaps[(1,5)]
    embedding = np.zeros((dim,4),dtype=complex)
    for c,psi in enumerate(code):
        for bits,amp in psi.items():
            idx=sum(bit<<(7-k) for k,bit in enumerate(bits))
            embedding[idx,c]=complex(sp.N(amp))
    assert np.allclose(embedding.conj().T@embedding,np.eye(4),atol=1e-12)
    assert np.linalg.norm(penalty@embedding)<1e-11
    assert np.linalg.norm(embedding.conj().T@cross@embedding)<1e-11
    h_eff=np.array(heff,float)
    results=[]
    for coupling in (0.4,0.2,0.1):
        H=penalty+coupling*cross
        energy,u=np.linalg.eigh(H)
        t=6*np.pi/coupling**2
        projected=(embedding.conj().T@u)@np.diag(np.exp(-1j*energy*t))@(u.conj().T@embedding)
        exact_energy,exact_u=np.linalg.eigh(coupling**2*h_eff)
        effective=exact_u@np.diag(np.exp(-1j*exact_energy*t))@exact_u.conj().T
        error=float(np.linalg.norm(projected-effective,ord=2))
        leakage=float(max(0,1-min(np.linalg.eigvalsh(projected.conj().T@projected))))
        results.append((coupling,t,error,leakage))
        assert error<1.0 and leakage<0.1
        print(f"J/Delta={coupling:.1f}, t*Delta={t:.3f}, "
              f"operator-deviation={error:.6f}, max-leakage={leakage:.6f}")
    assert results[0][2]>results[1][2]>results[2][2]
    assert results[0][3]>results[1][3]>results[2][3]
    print("Finite-pulse effect and leakage decrease with J/Delta, "
          "while required gate time increases: verified")


def test():
    neutral_ancilla_no_go()
    code,heff=code_and_second_order()
    full_eight_spin_dynamics(code,heff)
    print("ALL TWO-BODY SU2 MEDIATOR CHECKS PASSED")


if __name__=="__main__":
    test()
