#!/usr/bin/env python3
"""Figures for paper/generation_codes.md. Run from the repository root:
    python3 paper/make_figures.py
Writes PDF (for the LaTeX build) and SVG (for the HTML build) into paper/figures/.
Every plotted number is computed here with the exact routines in charge_lattice_tools.py."""
import os, sys
from fractions import Fraction
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle
sys.path.insert(0, os.path.abspath('.'))
from charge_lattice_tools import (W16, B5, coupled, single, class_weight_distributions,
                                  direct_sum_distributions, exact_value, dot)

OUT = os.path.join('paper', 'figures')
plt.rcParams.update({'font.family': 'serif', 'font.size': 10, 'mathtext.fontset': 'cm',
                     'axes.spines.top': False, 'axes.spines.right': False})
INK, MUTED, POS, NEG, ACCENT, WARM = '#1b1f24', '#6b7280', '#2f7d6d', '#8a4f7d', '#c47f17', '#b3423a'
LABELS = [r'$\nu^c$', r'$u^c_3$', r'$u^c_2$', r'$d_1$', r'$u_1$', r'$u^c_1$', r'$d_2$', r'$u_2$',
          r'$d_3$', r'$u_3$', r'$e^c$', r'$d^c_1$', r'$d^c_2$', r'$d^c_3$', r'$\nu$', r'$e$']


def save(fig, name):
    for ext in ('pdf', 'svg'):
        fig.savefig(os.path.join(OUT, f'{name}.{ext}'), bbox_inches='tight')
    plt.close(fig)


def fig1():
    """(a) one-generation charge matrix; (b) two-family coupled B5 lattice."""
    fig, axes = plt.subplots(2, 1, figsize=(7.2, 4.6), gridspec_kw={'height_ratios': [5, 7.6]})
    ax = axes[0]
    M = np.array(W16).T
    for r in range(5):
        for c in range(16):
            v = M[r, c]
            ax.add_patch(Rectangle((c, 4 - r), 0.92, 0.92, color=POS if v > 0 else NEG))
            ax.text(c + 0.46, 4 - r + 0.46, '+' if v > 0 else '−', ha='center', va='center', color='white', fontsize=9)
    ax.add_patch(Rectangle((-0.06, -0.06), 1.04, 5.04, fill=False, ec=ACCENT, lw=2))
    ax.set_xlim(-0.2, 16); ax.set_ylim(-0.2, 5.1); ax.set_aspect('equal'); ax.axis('off')
    for c, lab in enumerate(LABELS):
        ax.text(c + 0.46, 5.25, lab, ha='center', va='bottom', fontsize=8)
    for r in range(5):
        ax.text(-0.35, 4 - r + 0.46, f'$Q_{r+1}$', ha='right', va='center', fontsize=8)
    ax.set_title('(a) One generation: five charges × sixteen fields', loc='left', fontsize=10, pad=24)
    ax = axes[1]
    for g, s in enumerate((1, -1)):
        for c, q in enumerate(W16):
            col = list(q) + [s * dot(b, q) for b in B5]
            x = c + g * 16.6
            for r, v in enumerate(col):
                y = 6 - r
                if r < 5:
                    color = POS if v > 0 else NEG
                else:
                    color = POS if v > 0 else (NEG if v < 0 else '#d9dde3')
                ax.add_patch(Rectangle((x * 0.47, y * 0.47), 0.43, 0.43, color=color))
                if r >= 5:
                    ax.text(x * 0.47 + 0.215, y * 0.47 + 0.215, f'{v:+d}' if v else '0', ha='center', va='center',
                            color='white' if v else INK, fontsize=5.5)
        ax.add_patch(Rectangle((g * 16.6 * 0.47 - 0.03, -0.03 + 0.0), 0.49, 7 * 0.47, fill=False, ec=ACCENT, lw=1.6))
        ax.text((g * 16.6 + 8) * 0.47, 7 * 0.47 + 0.12, 'electron family ($s=+1$)' if g == 0 else 'muon family ($s=-1$)',
                ha='center', va='bottom', fontsize=8.5)
    ax.plot([-0.1, 32.6 * 0.47], [2 * 0.47 - 0.02] * 2, color=INK, lw=0.6, ls=':')
    ax.text(-0.15, 4.5 * 0.47, 'shared\nSO(10)', ha='right', va='center', fontsize=7.5, color=MUTED)
    ax.text(-0.15, 1.0 * 0.47, r'$s\,Bq$', ha='right', va='center', fontsize=8, color=MUTED)
    ax.set_xlim(-0.9, 32.6 * 0.47); ax.set_ylim(-0.1, 7 * 0.47 + 0.45); ax.set_aspect('equal'); ax.axis('off')
    ax.set_title(r'(b) Two coupled generations, $B=B_5$: five shared rows and two differential rows', loc='left', fontsize=10)
    fig.tight_layout()
    save(fig, 'fig1_charge_matrix')


def fig2():
    """Decision path: the twins problem and its resolutions."""
    fig, ax = plt.subplots(figsize=(7.4, 2.6))
    ax.axis('off'); ax.set_xlim(0, 120); ax.set_ylim(0, 40)
    boxes = [
        (1, 'shared charges', r'$[[30,1,3]]_p$', 'one logical qudit', WARM),
        (25, 'family number', r'$[[30,2,2]]_p$', 'twins are logical', WARM),
        (49, 'smart badges', r'$[[30,2,3]]_p$', 'breaks SU(3)', ACCENT),
        (73, r'family $B{-}L$, $T_{3R}$', r'$[[30,2,3]]_{p\geq5}$', 'gauge compatible', POS),
        (97, 'real Yukawas', r'$k=1$', 'labels erased', WARM),
    ]
    for x, title, code, note, color in boxes:
        ax.add_patch(FancyBboxPatch((x, 11), 21, 20, boxstyle='round,pad=0.4,rounding_size=1.5',
                                    fc='white', ec=color, lw=1.6))
        ax.text(x + 10.5, 26.5, title, ha='center', va='center', fontsize=8, color=INK)
        ax.text(x + 10.5, 21, code, ha='center', va='center', fontsize=9, color=INK)
        ax.text(x + 10.5, 15, note, ha='center', va='center', fontsize=7.4, color=color)
    for x in (22.6, 46.6, 70.6, 94.6):
        ax.annotate('', xy=(x + 2.0, 21), xytext=(x, 21), arrowprops=dict(arrowstyle='-|>', color=INK, lw=0.9))
    ax.text(60, 4, 'Each step fixes the failure of the one before; the last step is imposed by nature.',
            ha='center', fontsize=8.2, color=MUTED)
    save(fig, 'fig2_path')


def fig3():
    """Logical error rates after post-selection, coupled B5 versus direct sum, p = 5."""
    p = 5
    c = coupled(B5, 2, p)
    Dc = class_weight_distributions(c)
    D1 = class_weight_distributions(single(p))
    Dd = direct_sum_distributions(D1, 15, 2)
    Dd = {k: (v + [0] * 31)[:31] for k, v in Dd.items()}
    nz = lambda cc: sum(1 for x in cc if x)

    def curves(D):
        acc = [sum(D[k][w] for k in D) for w in range(31)]
        one = [sum(D[k][w] for k in D if nz(k) == 1) for w in range(31)]
        mul = [sum(D[k][w] for k in D if nz(k) >= 2) for w in range(31)]
        eps = np.logspace(-4, -1, 61)
        s, m = [], []
        for e in eps:
            E = Fraction(float(e)).limit_denominator(10 ** 12)
            pa = exact_value(acc, p, 30, E)
            s.append(float(exact_value(one, p, 30, E) / pa))
            m.append(float(exact_value(mul, p, 30, E) / pa))
        return eps, np.array(s), np.array(m)
    eps, sc, mc = curves(Dc)
    _, sd, md = curves(Dd)
    fig, ax = plt.subplots(figsize=(5.6, 3.9))
    ax.loglog(eps, sc, color=POS, lw=2, label='coupled, one output wrong')
    ax.loglog(eps, sd, color=POS, lw=1.2, ls='--', label='independent, one output wrong')
    ax.loglog(eps, mc, color=WARM, lw=2, label='coupled, both outputs wrong')
    ax.loglog(eps, md, color=WARM, lw=1.2, ls='--', label='independent, both outputs wrong')
    ax.set_xlabel(r'physical error rate $\varepsilon$')
    ax.set_ylabel('probability given acceptance')
    ax.set_ylim(1e-26, 1)
    ax.text(1.15e-4, 3e-14, r'$\frac{9}{2}\,\varepsilon^4$', color=WARM, fontsize=10)
    ax.text(1.15e-4, 2e-22, r'$\frac{625}{256}\,\varepsilon^6$', color=WARM, fontsize=10)
    ax.text(1.15e-4, 2e-10, r'$\frac{25}{8}\,\varepsilon^3$  (both curves)', color=POS, fontsize=10)
    ax.legend(frameon=False, fontsize=8, loc='lower right')
    ax.set_title(r'Two generations at $p=5$: coupled $[[30,2,3]]_5$ vs. two $[[15,1,3]]_5$', fontsize=9.5, loc='left')
    save(fig, 'fig3_channels')


if __name__ == '__main__':
    os.makedirs(OUT, exist_ok=True)
    fig1(); fig2(); fig3()
    print('figures written to', OUT)
