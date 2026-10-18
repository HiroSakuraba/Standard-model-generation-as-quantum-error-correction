"""Independent exact nine-state benchmark against Chen, Mueller & Yao.

Primary source: Chen, Mueller, Yao, Phys. Rev. D 114, 074502 (2026),
DOI 10.1103/1qt7-3ywl, arXiv:2601.10065v1.
Equation (3): Kogut-Susskind normalization; equations (10)-(11): square
chain with zero exterior plaquette flux; Appendix B Table VI (published; Table IV in arXiv v1): the magnetic
transition blocks. The three Table-IV boundary-zero blocks are explicitly
transcribed in this script, independent of Haar/loop expressions.

Run (requires SymPy, NumPy and both preceding repo scripts):
  OPENBLAS_NUM_THREADS=1 python su3_published_square_chain_benchmark.py
This is a comparison of two finite-link truncations, NOT glueball spectroscopy.
"""
from itertools import product
import sympy as sp
import numpy as np

from su3_complete_minimal_links import (
    NEW_LABELS, enumerate_basis, loop_pattern, inner,
)
from su3_two_plaquette_spin_network import EDGES

# Z3 plaquette electric-flux numbers in the fixed order of the previously
# established nine-state spin-network basis. 2 == -1 (antifundamental).
PLAQUETTE = {
    'vac': (0, 0),
    'L+': (1, 0), 'L-': (2, 0),
    'R+': (0, 1), 'R-': (0, 2),
    'O+': (1, 1), 'O-': (2, 2),
    'B+': (1, 2), 'B-': (2, 1),
}


def link_flux(left, right):
    """Link reps 0, +/-1 for the explicit graph in prior Gauss audit.

    Oriented edge order: 0(top L),1(top R),2(bottom L),3(bottom R),
    4(left side),5(shared side),6(right side). This is the paper's
    'link = active-face minus adjacent-face' rule, equations (5),(10),(11).
    """
    mods = (left, right, -left, -right, -left, left - right, right)
    return tuple(0 if n % 3 == 0 else (1 if n % 3 == 1 else -1) for n in mods)


def published_boundary_blocks():
    """Independent transcription: Appendix B Table VI (published; Table IV in arXiv v1), (PL,PR)=(0,n).

    Each 3x3 real symmetric block is ordered active plaquette qutrit 0,1,2.
    Exterior top and bottom plaquettes are fixed at 0 (paper p. 4).
    No Haar integrations or previous-code operator values are used.
    """
    t = sp.sqrt(3) / 3
    a = sp.Rational(1, 3)
    return {
        0: sp.Matrix([[0,1,1],[1,0,1],[1,1,0]]),
        1: sp.Matrix([[0,a,t],[a,0,t],[t,t,0]]),
        2: sp.Matrix([[0,t,a],[t,0,t],[a,t,0]]),
    }


def published_matrix(face, oriented=False):
    """Construct full 9x9 two-face Table-IV magnetic operator.

    For each plaquette, the other dynamical face and the frozen exterior
    face are its two neighbors. Only (0,n) or (n,0) blocks are needed.
    The directed operator implements q -> q+1 mod 3, equation (4).
    """
    n = len(NEW_LABELS)
    B = sp.zeros(n)
    for i, dst in enumerate(NEW_LABELS):
        dleft, dright = PLAQUETTE[dst]
        for j, src in enumerate(NEW_LABELS):
            sleft, sright = PLAQUETTE[src]
            if face == 'L':
                if dright != sright: continue
                active_dst, active_src, neighbor = dleft, sleft, sright
            elif face == 'R':
                if dleft != sleft: continue
                active_dst, active_src, neighbor = dright, sright, sleft
            else: raise ValueError(face)
            if oriented and active_dst != (active_src + 1) % 3:
                continue
            B[i,j] = published_boundary_blocks()[neighbor][active_dst,active_src]
    if oriented:
        assert B + B.T == published_matrix(face)
    else:
        assert B == B.T
    return B


def published_electric():
    """Paper equation (3), nl=2 and Sigma=1 for square plaquette chain."""
    electric = []
    for lab in NEW_LABELS:
        ql, qr = PLAQUETTE[lab]
        n_fundamental_links = 3*int(ql != 0) + 3*int(qr != 0) + int(ql != qr)
        electric.append(sp.Rational(4,3)*n_fundamental_links)
    return sp.diag(*electric)


def independently_normalized_operators():
    """Build prior Haar operators using physical normalized qutrit states."""
    G = inner()
    assert G == sp.diag(1,1,1,1,1,1,1,sp.Rational(4,3),sp.Rational(4,3))
    # Raw basis |B> has norm squared 4/3, so normalized states are
    # |B_norm> = sqrt(3)/2 |B_raw>.
    S = sp.diag(*([sp.Integer(1)]*7), sp.sqrt(3)/2, sp.sqrt(3)/2)
    assert S*G*S == sp.eye(9)
    wl = sp.simplify(S.inv()*G.inv()*inner('L+')*S)
    wr = sp.simplify(S.inv()*G.inv()*inner('R+')*S)
    return wl, wr


def compare():
    assert len(NEW_LABELS) == len(PLAQUETTE) == 9
    assert set(NEW_LABELS) == set(PLAQUETTE)
    assert len(set(PLAQUETTE.values())) == 9 == len(list(product(range(3), repeat=2)))
    valid, counts = enumerate_basis()
    assert len(valid) == 9 and sum(counts.values()) == 9
    # A direct bijection between the paper's two plaquette qutrit numbers
    # and all nine microscopic link-representation assignments.
    mapped = [link_flux(*PLAQUETTE[label]) for label in NEW_LABELS]
    assert len(set(mapped)) == 9 and set(mapped) == set(valid)
    assert all(mapped[i] == loop_pattern(NEW_LABELS[i]) for i in range(7))
    assert mapped[7] == (1,-1,-1,1,-1,-1,-1)
    assert mapped[8] == tuple(-v for v in mapped[7])
    print('EXACT BASIS BIJECTION: all 3^2=9 face qutrits match the complete 3^7 link cutoff Gauss sector')

    physical_wl, physical_wr = independently_normalized_operators()
    paper_wl = published_matrix('L', oriented=True)
    paper_wr = published_matrix('R', oriented=True)
    assert sp.simplify(physical_wl - paper_wl) == sp.zeros(9)
    assert sp.simplify(physical_wr - paper_wr) == sp.zeros(9)
    for face, physical, paper in [('L',physical_wl,paper_wl),('R',physical_wr,paper_wr)]:
        assert sp.simplify(physical+physical.T-published_matrix(face))==sp.zeros(9)
        # All 81 exact matrix entries, not just a common eigenvalue.
        print(f'EXACT PUBLISHED W_{face}: 81 of 81 directed matrix entries agree; Hermitian sum agrees')
        examples=sorted({str(x) for x in paper if x!=0})
        print(f'  directed nonzero values: {examples}; number of nonzeros: {sum(x!=0 for x in paper)}')

    D=published_electric()
    independent_configs = [link_flux(*PLAQUETTE[label]) for label in NEW_LABELS]
    from_links=sp.diag(*[sum(sp.Rational(4,3) for rep in cfg if rep) for cfg in independent_configs])
    assert D == from_links
    assert D == sp.diag(0,*[sp.Rational(16,3)]*4,*[sp.Integer(8)]*2,*[sp.Rational(28,3)]*2)
    print('EXACT ELECTRIC: all 9 Casimir diagonal energies agree; all 7 link assignments checked')

    paper_b=published_matrix('L')+published_matrix('R')
    prior_b=physical_wl+physical_wl.T+physical_wr+physical_wr.T
    assert sp.simplify(paper_b-prior_b)==sp.zeros(9)
    # Equation (3) for square plaquettes nl=2 and Sigma=1:
    # H = g^2 D + (12 I - B_L - B_R)/g^2.
    # 12/g^2 is the trivial plaquette magnetic diagonal shift, explicitly
    # dropped after equation (3) in the paper. The preexisting dimensionless
    # H(c)=D-c/2*(B_L+B_R) corresponds to c=2/(g^2)^2.
    I=sp.eye(9)
    for g2 in (sp.Rational(1,2), sp.Integer(1), sp.Rational(3,2), sp.Integer(2)):
        hpaper = g2*D + (12*I-paper_b)/g2
        c = 2 / (g2*g2)
        hold = D - c*prior_b/2
        assert sp.simplify(hpaper-(g2*hold+12*I/g2))==sp.zeros(9)
    print('EXACT HAMILTONIAN: full 9x9 electric + magnetic matches PRD Eq. (3) after overall scale and scalar shift')

    # Independent spectral comparison at the paper's g²=1 and 2.
    for g2 in (1.,2.):
        physical = np.array(g2*D + (12*I-paper_b)/g2,float)
        recomposed = np.array(g2*D + (12*I-prior_b)/g2,float)
        assert np.array_equal(physical,recomposed)
        eig = np.linalg.eigvalsh(physical)
        print(f'g²={g2:g}: E0={eig[0]:.10f}, gap={eig[1]-eig[0]:.10f} (dimensionless; numerical)')
    print('PUBLISHED SQUARE-CHAIN MATRIX-ENTRY BENCHMARK ALL PASS')

if __name__ == '__main__':
    compare()
