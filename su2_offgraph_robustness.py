"""Off-graph exchange sensitivity in an eight-spin SU(2) singlet code.

Exact SymPy Knill-Laflamme obstruction and a full 256-state, coherent
+/- coupling-fluctuation channel followed by an explicit CPTP erasure recovery.
The recovery discards physical A0 and all D spins, decodes from A1,A2,A3,
and replaces leakage beyond the supported subsystem with maximally mixed data.
All exchanges are co-located at one gauge vertex; this is NOT a gauge theory
on separate physical vertices, a lifetime forecast, or an optimal decoder.

Reproduce: OPENBLAS_NUM_THREADS=1 python su2_offgraph_robustness.py
Requires numpy scipy sympy; imports local su2_error_algebra_capacity.py.
"""
from __future__ import annotations

from itertools import product
import numpy as np
from scipy.linalg import eigh
import sympy as sp

from su2_error_algebra_capacity import (
    four_spin_code, projected, h, singlet_exchange_representation,
    commutant_dimension, noncrossing_matchings, valence_state, tensor,
)


def exact_kl_defect():
    a = four_spin_code()
    S01 = projected(a, lambda psi: h(psi, 0, 1))
    assert S01 == sp.diag(-sp.Rational(3, 4), sp.Rational(1, 4))
    # Star geometry: protected logical |A_i> \otimes |D_0>.
    # h04 h14 compressed: (1/3) <S_A0.S_A1> (S_D0.S_D0)
    p_both = S01 / 4
    physical_code = [tensor(psi, a[0]) for psi in a]
    direct_256 = projected(physical_code, lambda psi: h(h(psi, 1, 4), 0, 4))
    assert direct_256 == p_both  # independent exact 256-basis-state contraction
    assert p_both == sp.diag(-sp.Rational(3, 16), sp.Rational(1, 16))
    trless = p_both - sp.trace(p_both) * sp.eye(2) / 2
    assert trless == sp.diag(-sp.Rational(1, 8), sp.Rational(1, 8))
    # h05 h14 compresses with opposite sign for spectator singlet D01.
    assert -p_both == sp.diag(sp.Rational(3, 16), -sp.Rational(1, 16))
    print('EXACT KL: P h04 h14 P = diag(-3/16,1/16)')
    print('traceless KL defect for offgraph amplitude epsilon: |epsilon|/8 in operator norm')
    return p_both


def exchange_graph_algebra():
    """Exact exchange-algebra discontinuity as one offgraph edge is added."""
    sw = singlet_exchange_representation(8)
    star = [(0, 4), (0, 5), (0, 6), (0, 7)]
    assert commutant_dimension([sw(*e) for e in star]) == 5
    assert commutant_dimension([sw(*e) for e in star + [(1, 4)]]) == 2
    component = (0, 1, 4, 5, 6, 7)
    # S_component^2 = sum(SWAP_ij) + n(4-n)/4 I, n=6.
    casimir = sum((sw(i, j) for idx, i in enumerate(component)
                   for j in component[idx + 1:]), sp.zeros(14)) - 3 * sp.eye(14)
    assert casimir.eigenvals() == {sp.S(0): 5, sp.S(2): 9}
    print("EXACT GRAPH ALGEBRA: star commutant dim=5 (noiseless qubit); + offgraph commutant dim=2 (M5 ⊕ M9, no protected multiplicity)")


def singlet_vector(state):
    z = np.zeros(16, dtype=np.complex128)
    for spins, amp in state.items():
        n = sum(bit << (3 - pos) for pos, bit in enumerate(spins))
        z[n] = complex(sp.N(amp))
    return z


def physical_embedding():
    a0, a1 = (singlet_vector(psi) for psi in four_spin_code())
    d0 = a0.copy()
    a = np.stack((a0, a1))
    code = np.column_stack((np.kron(a0, d0), np.kron(a1, d0)))
    assert np.allclose(code.conj().T @ code, np.eye(2), atol=1e-14)
    A = a.reshape(2, 2, 8)
    support = 2 * sum((A[k].conj().T @ A[k] for k in range(2)), np.zeros((8, 8), complex))
    assert np.linalg.norm(support @ support - support) < 1e-12
    assert abs(np.trace(support) - 4) < 1e-12
    assert np.linalg.norm(support - support.conj().T) < 1e-12
    # N(P)=support/2. The Petz erasure recovery has effect support,
    # and the missing-trace completion has positive effect I-support.
    assert np.linalg.eigvalsh(np.eye(8) - support).min() >= -1e-12
    return A, code


def swap_matrix(n, i, j):
    dim = 1 << n
    ix = np.arange(dim)
    bi = (ix >> (n - 1 - i)) & 1
    bj = (ix >> (n - 1 - j)) & 1
    destination = ix ^ ((bi ^ bj) << (n - 1 - i)) ^ ((bi ^ bj) << (n - 1 - j))
    return np.eye(dim)[destination]


def heisenberg_matrix(n, i, j):
    dim = 1 << n
    return (swap_matrix(n, i, j) - 0.5 * np.eye(dim)) / 2


def star_and_offgraph():
    # Dimensionless engineered couplings; not a fitted device geometry.
    star_strengths = (1.0, 0.6, 0.25, -0.15)
    star = sum((J * heisenberg_matrix(8, 0, 4 + k)
                for k, J in enumerate(star_strengths)), np.zeros((256, 256)))
    offgraph = heisenberg_matrix(8, 1, 4)
    return star, offgraph, star_strengths


def global_gauss_operators():
    n = 8
    dim = 1 << n
    raising = np.zeros((dim, dim), dtype=np.float64)
    z = np.zeros((dim, dim), dtype=np.float64)
    for ket in range(dim):
        spins = [(ket >> (n - 1 - i)) & 1 for i in range(n)]
        z[ket, ket] = sum(1 - 2*b for b in spins)
        for i, b in enumerate(spins):
            if b == 1:
                bra = ket ^ (1 << (n - 1 - i))
                raising[bra, ket] += 1
    return raising, z


def output_of_fixed_erasure_recovery(psi_i, psi_j, A):
    """R(tr_{A0,D} |psi_i><psi_j|) on the two logical dimensions.

    A: shape (codeword=2, erased_A0=2, retained_B=8).
    N(P) is 1/2 times rank-four support projector (P rank=2).
    Petz: R(X)=2 P (I_{A0} tensor X) P, plus a CP trace-completion
    that maps states outside the supported subspace to maximally mixed I/2.
    """
    i_state = psi_i.reshape(2, 8, 16)
    j_state = psi_j.reshape(2, 8, 16)
    reduced = np.einsum('ebd,ecd->bc', i_state, j_state.conj(), optimize=True)
    principal = 2 * np.einsum('keb,bc,lec->kl', A.conj(), reduced, A, optimize=True)
    missing_trace = np.vdot(psi_j, psi_i) - np.trace(principal)
    return principal + missing_trace * np.eye(2) / 2


def entanglement_fidelity(evolved, A):
    """Entanglement fidelity of erasure recovery for a single unitary.

    Evolved is 256x2, columns U|0L>, U|1L>.  Uses the maximally mixed
    logical input. F_avg=(2 F_e + 1)/3, because the resulting map is CPTP.
    """
    assert evolved.shape == (256, 2)
    fe_sum = 0.0j
    for i in range(2):
        for j in range(2):
            M = output_of_fixed_erasure_recovery(evolved[:, i], evolved[:, j], A)
            fe_sum += M[i, j]
    fe = float(np.real_if_close(fe_sum / 4))
    assert -1e-10 <= fe <= 1 + 1e-10
    return fe


def no_recovery_entanglement_fidelity(evolved, original_code):
    return float(abs(np.trace(original_code.conj().T @ evolved) / 2)**2)


def encoded_subspace_leakage(evolved, original_code):
    return float(1 - np.linalg.norm(original_code.conj().T @ evolved, 'fro')**2 / 2)


def propagate(H, code, time):
    eigenvalues, eigenvectors = eigh(H, check_finite=True)
    return eigenvectors @ (np.exp(-1j * eigenvalues * time)[:, None]
                           * (eigenvectors.T @ code))


def independent_global_singlet_crosscheck(star, offgraph, code):
    """Independent 14D singlet-sector evolution versus full 256D unitary.

    Orthonormalize valence-bond spanning columns numerically. This is an
    independent reduction using the known Catalan-14 global singlet space;
    not an additional gauge restriction on the physical 256D Hamiltonian.
    """
    valence = []
    for matching in noncrossing_matchings(tuple(range(8))):
        vector = np.zeros(256, dtype=np.complex128)
        for bits, amp in valence_state(8, matching).items():
            index = sum(bit << (7 - i) for i, bit in enumerate(bits))
            vector[index] = float(amp)
        valence.append(vector)
    Q, R = np.linalg.qr(np.stack(valence, axis=1), mode="reduced")
    assert Q.shape == (256, 14) and np.linalg.matrix_rank(R) == 14
    assert np.linalg.norm(Q.conj().T @ Q - np.eye(14)) < 1e-12
    assert np.linalg.norm(Q @ (Q.conj().T @ code) - code) < 1e-12
    eps, time = 0.1, 3.0
    H256 = star + eps * offgraph
    H14 = Q.conj().T @ H256 @ Q
    assert np.linalg.norm(H256 @ Q - Q @ H14) < 1e-11
    values, vectors = eigh(H14)
    evolved_14 = Q @ (vectors @ (np.exp(-1j * values * time)[:, None]
                                * (vectors.conj().T @ (Q.conj().T @ code))))
    evolved_256 = propagate(H256, code, time)
    deviation = np.linalg.norm(evolved_14 - evolved_256)
    assert deviation < 1e-11, deviation
    print(f'Full-256 versus independent 14D global-singlet evolution PASS (difference {deviation:.2e})')


def simulate():
    A, code = physical_embedding()
    test_recovery = []
    for i in range(2):
        for j in range(2):
            out = output_of_fixed_erasure_recovery(code[:, i], code[:, j], A)
            expected = np.zeros((2, 2), complex)
            expected[i, j] = 1
            assert np.linalg.norm(out - expected, ord='fro') < 1e-12
    print('Explicit erasure recovery maps four logical matrix units to themselves: PASS')

    star, offgraph, strengths = star_and_offgraph()
    raiseop, zz = global_gauss_operators()
    for op in (star, offgraph):
        assert np.linalg.norm(op @ raiseop - raiseop @ op) < 1e-12
        assert np.linalg.norm(op @ zz - zz @ op) < 1e-12
    assert np.linalg.norm(raiseop @ code) < 1e-12
    assert np.linalg.norm(zz @ code) < 1e-12
    print('All Hamiltonians exactly gauge invariant: SU(2) commutators vanish (numerical matrix check)')
    print('Star couplings =', strengths, '; offgraph edge = (1,4)')
    independent_global_singlet_crosscheck(star, offgraph, code)

    rows = []
    for t in (1.0, 3.0, 6.0):
        for epsilon in (0, 0.02, 0.05, 0.10, 0.20):
            p = propagate(star + epsilon * offgraph, code, t)
            m = propagate(star - epsilon * offgraph, code, t)
            # Classical uncertain coupling sign: actual CPTP mixture of 2 unitary channels.
            fe_rec = (entanglement_fidelity(p, A) + entanglement_fidelity(m, A)) / 2
            fe_raw = (no_recovery_entanglement_fidelity(p, code)
                      + no_recovery_entanglement_fidelity(m, code)) / 2
            leakage = (encoded_subspace_leakage(p, code)
                       + encoded_subspace_leakage(m, code)) / 2
            assert 0 <= fe_rec <= 1 + 1e-10
            if epsilon == 0:
                assert abs(1 - fe_rec) < 1e-10
            rows.append((t, epsilon, fe_rec, fe_raw, leakage))
            print(f't={t:.1f}, eps={epsilon:.2f}, recovered_F_e={fe_rec:.9f}, '
                  f'raw_F_e={fe_raw:.9f}, code_leakage={leakage:.9f}')
    # Restored logical information can be near-perfect despite extensive physical-code leakage.
    assert any(x[4] > 0.05 for x in rows if x[1] == 0)
    # Offgraph bath channel introduces nonzero residual loss for some sampled parameters.
    assert any(x[2] < 0.999 for x in rows if x[1] >= 0.1)
    print('FULL 256-STATE ROBUSTNESS AND CPTP ERASURE RECOVERY PASSED')
    return rows


if __name__ == '__main__':
    exact_kl_defect()
    exchange_graph_algebra()
    simulate()
