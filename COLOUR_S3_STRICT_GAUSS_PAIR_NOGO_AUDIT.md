# Strict Gauss-preserving two-register operations: a finite-S3 syndrome no-go

**19 October 2026.** Continuation after PR #49. Applies only to the finite six-element S3 subgroup simulator and its **four-dimensional** phase code, not to the original nine-quark 21D SU(3) code.

## Problem and result

PR #49 found an exact decomposition of the logical phase-syndrome readout into 16 *general* two-register gates, but some intermediate gates nearly completely violated Gauss law. Could **different** two-register gates remain gauge invariant at *every* step and still measure the same known dressed-swap error?

**Not with the current physical registers and gauge-neutral pointers.** The obstruction is independent of circuit depth:

- A six-state S3 link transforms by **left regular action at A** and **right regular action at B**, independently.
- By the finite-group Peter–Weyl decomposition, the link has three irreducible representation labels: trivial (\(d=1\)), sign (\(d=1\)), and standard (\(d=2\)). Their link-space projector ranks are \(1,1,4\).
- On the **14D simultaneous vertex-Gauss space**, these link sectors have ranks **4,1,9**. The physical subspace may be represented as \(\bigoplus_r {\rm End}(\mathbb C^{m_r})\), with multiplicities \(m_r=(2,1,3)\) for the two-qutrit matter pair at each vertex.
- A one- or two-register gate that commutes with **both** local Gauss actions and operates on only one endpoint's matter plus the link must commute with the **other endpoint's** regular gauge action on the link. Its link irrep label is therefore conserved. Operators without the link conserve the label trivially. A neutral pointer does not change this. Thus **every** sequence of such gates preserves the three sectors.

The source dressed swap \(S\) mixes those sectors. It cannot be implemented exactly with this restricted invariant two-register alphabet, no matter how many steps are used.

### Stronger obstruction: exact logical-syndrome distinction

Let \(K\) be the four-dimensional \(\omega\)-eigenspace code of \(U=A_{\rm onsite}S_{\rm dressed}\) from PR #48, and \(F=SK\) its dressed-swap error space. Globally \(K^\dagger F=0\), so an unrestricted high-level measurement can distinguish them.

Nevertheless the individual sector overlaps are nonzero:

| Link representation r | Gauss-sector rank | \(\|K^\dagger P_r F\|_F\) | \(K\) and \(F\) sector probability |
|---|---:|---:|---:|
| trivial | 4 | \(1/2\) | \(5/24\) |
| sign | 1 | \(1/2\) | \(1/8\) |
| standard | 9 | \(1/\sqrt3\) | \(2/3\) |

These simple numerical values arise from the exact finite-matrix construction; tolerance-based regression checks compare them. If a **deterministic nondemolition** measurement had a "no error" effect \(E\) which accepted every K state and rejected every F state, while commuting with every \(P_r\), the two projected spaces \(P_rK\) and \(P_rF\) would have to be orthogonal separately. They are not. Thus **no such exact measurement is possible** using the stated restricted gate set and neutral pointer ancillas.

This is a model-specific **selection-rule** result, not a general quantum error-correction no-go. It does not rule out a noisy or incomplete syndrome, different codes, additional gauge-field links, gauge-charged ancillas, three-register operations, or intermediate gauge violations.

## Pair-supported gauge-invariant algebra: independent numerical check

On the physical 14D sector, let \(W\) be the exact isometry. For any operator \(O\) supported on at most two registers, compress to \(W^\dagger O W\). Twirl O by the two vertex gauge groups. The twirled operator **still has pair support**, is gauge invariant, and has the **same compression** because W is invariant.

The matrix-unit basis for all ten pairs spans 29 linearly independent 14×14 compressed matrices. Closing this span under multiplication yields ranks

\[
29 \longrightarrow 98 \longrightarrow 98,
\]

matching the upper bound

\[
4^2+1^2+9^2 = 98.
\]

This is consistent with the full **block-diagonal algebra \(M_4\oplus M_1\oplus M_9\)** available under the symmetry restriction and confirms that neither the dressed swap nor the PR #48 phase-projector readout belongs to that algebra. The rank computation is double-precision numerical evidence; the obstruction itself follows from exact irrep conservation.

## Positive control: what can two-register invariant hardware measure?

One qutrit pointer can measure the **link irrep label** with a single *general* two-register 18×18 unitary

\[
V_{pL}=\sum_{r=0}^{2} X_p^{\,r}\otimes P_r^{\rm link},
\]

where \(X_p\) cyclically shifts the pointer. Prepare the pointer in \(|0\rangle\). Its three readout Kraus maps are precisely \(P_r\). The script verifies complete trace preservation, commutation with both vertex Gauss transformations, and exactly zero gauge-sector departure on physical inputs.

**However**, both K and F have the same link-irrep probability distribution \((5/24,\,1/8,\,2/3)\). The pointer is a genuine gauge-invariant local irrep measurement, but **not** a syndrome that detects the known dressed swap. The two sector-conditioned spaces also overlap, so repeated sector readouts alone cannot solve this.

The pointer–link controlled irrep gate is still an arbitrary U(18) **two-register** primitive; no pulse sequence, native elementary entangler count, or implementation time is provided. More fundamentally, the measured link irrep is an actual superselection label under the strict primitive restriction; it is not a replacement for the desired logical eigenphase.

## Decision

For **this exact model and register inventory**, one cannot demand *simultaneously*:

1. Every elementary gate is two-register or one-register only;
2. Every intermediate gate commutes with both independent vertex Gauss symmetries;
3. All pointer registers are gauge neutral;
4. The error syndrome perfectly distinguishes K from its dressed-swap image without revealing logical information.

At least one condition must change. The most useful next constructions to test are a **gauge-invariant three-register** matter–link–matter interaction as an elementary operation, a **gauge-charged mediator** whose propagation respects Gauss law, or a **gauge-penalty control** that permits and quantifies temporary violations. Any resulting scheme still needs calibrated gate synthesis and a matched noise/storage baseline.

This no-go is strictly for the finite S3 subgroup, the four-dimensional phase code, and the specified registers; it **does not establish** a no-go for continuous SU(3), physical confinement, or the original nine-quark 21D/65-syndrome code.

## Reproduce

\`\`\`bash
python -m pip install numpy scipy
OPENBLAS_NUM_THREADS=1 python -m unittest -v test_colour_s3_strict_gauss_pair_nogo
OPENBLAS_NUM_THREADS=1 python colour_s3_strict_gauss_pair_nogo.py
\`\`\`

The dedicated CI workflow archives a JSON record of link sectors, overlap witnesses, covariance, pointer outcomes, numerical algebra closure and all claim limitations.
