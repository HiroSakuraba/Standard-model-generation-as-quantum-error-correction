# Walsh characters on the 16 SO(10) weights: exact rigidity and a disguised direct sum

**Research note — 17 October 2026.** This is a mathematical result about finite-field CSS codes. It does **not** establish new Standard Model symmetries, physical quantum error correction, or new inequivalent qudit codes. Run \`python3 walsh_character_lattice.py\` (requires NumPy and the existing \`charge_lattice_tools.py\`) to reproduce the finite enumerations and the exact code isomorphisms.

## 1. The source of the sixteen columns

The 16 chiral spinor sign vectors are

\[
S=\{q=(q_1,\ldots,q_5)\in\{\pm1\}^5:\prod_iq_i=1\}.
\]

Choosing the first four signs determines the last: \(q_5=q_1q_2q_3q_4\). Thus \(S\) is naturally the Boolean group \(G=(\mathbb Z_2)^4\), and its **16 Walsh characters** are the functions
\[
\chi_a(q)=\prod_{i=1}^4q_i^{a_i},\quad a\in\mathbb F_2^4.
\]

Indexing characters by four-bit masks, the five Cartan charge rows are
\[
A=\{1,2,4,8,15\},\qquad q_i=\chi_{a_i}.
\]

Over every field of odd characteristic, the character table is invertible:
\[
\sum_{q\in S}\chi_a(q)\chi_b(q)=16\,\delta_{a,b}.
\]

This is a restricted evaluation/character-code viewpoint, **not** the length-\(p^m\) full-evaluation Reed–Muller construction excluded in \`FULL_REED_MULLER_COMPARISON.md\`. Evaluation codes on chosen point sets are a standard coding-theory construction; see [Macaulay2 CodingTheory evaluation-code documentation](https://www.macaulay2.com/doc/Macaulay2/share/doc/Macaulay2/CodingTheory/html/___Evaluation__Code.html) and [Galindo, Hernando, and Ruano, *New quantum codes from evaluation and matrix-product codes* (2015)](https://doi.org/10.1016/j.ffa.2015.07.003).

## 2. Theorem — complete one-generation rigidity for **every odd prime**

Consider **any** function \(y:S\rightarrow\mathbb F_p\) that might be added as a sixth row to the original five Cartan charge rows while retaining cubic triorthogonality. The necessary mixed cubic identities are
\[
\sum_{q\in S} y(q)\,q_iq_j=0,\qquad 1\leq i\leq j\leq5.
\]

There are 11 *distinct* pair-product characters:
\[
\{0\}\cup\{a_i\oplus a_j:i<j\}.
\]
The ten distinct pairwise XOR masks are exactly all characters outside \(\{0\}\cup A\). Since the Walsh basis has 16 characters, orthogonality to those eleven forces
\[
\boxed{y\in\mathrm{span}_{\mathbb F_p}\{q_1,q_2,q_3,q_4,q_5\}.}
\]

Thus **no new independent row exists**, linear or nonlinear, while retaining the original five Cartan rows and their cubic conditions within a single generation. This is a proof for **all odd primes**, strengthening the earlier computer checks restricted to \(p=3,5,7\).

This theorem requires the full 16-column, single-family sign set and all five original rows. It says nothing about dropping some original rows, altering the point set, or using additional families to cancel cubic terms.

## 3. Coupled character rows: complete five-row classification

For two families, retain the shared five rows \((q_i,q_i)\) and propose nonlinear family-differential rows
\[
(\chi_{s_j},-\chi_{s_j}),\qquad s_j\in\mathbb F_2^4.
\]

The cubic identities with two shared rows cancel between families. The identities with three differential rows also cancel. The remaining condition with one shared and two differential rows gives
\[
2\sum_{q\in S}q_i\chi_{s_a}(q)\chi_{s_b}(q)=0
\quad\Longleftrightarrow\quad s_a\oplus s_b\notin A\quad(a\ne b).
\]
The condition is exact over **the integers**, not merely modulo \(p\); identical masks also pass because \(\sum q_i=0\).

For \(m\) distinct chosen masks, an undetectable same-species twin pair occurs when all \(m\) character values agree. Let
\[
r=\operatorname{rank}_{\mathbb F_2}
\{s_j\oplus s_1:j>1\}.
\]
There are exactly \(2^{4-r}-1\) unpunctured weights with this property. Hence eliminating all matched-twin weight-two logical errors requires \(r=4\) and therefore **at least five differential character rows**.

The exhaustive finite combinatorial enumeration over the **16-character set** gives:

| Number of rows \(m\) | Compatible unordered mask sets |
|---:|---:|
| 1 | 16 |
| 2 | 80 |
| 3 | 160 |
| 4 | 120 |
| 5 | **16** |
| 6 or more | **0** |

All 16 compatible five-mask sets are precisely
\[
\boxed{\,S_t=\{a\oplus t:a\in A\},\quad t\in\mathbb F_2^4.\,}
\]
Each one removes the matched-twin distance-two obstruction. The script verifies this complete finite classification, the rank formula, and the absence of further rows.

## 4. The surprise: every resulting distance-three code is a **direct sum in disguise**

Take a valid five-row set \(S_t\). Since
\[
\chi_{a_i\oplus t}(q)=\chi_t(q)q_i,
\]
the full charge vector for the two *original* families is
\[
h(q,\pm)=\begin{pmatrix}q\\\pm\chi_t(q)q\end{pmatrix}.
\]

Over any odd prime field, apply the invertible row transformation
\[
T=\frac12\begin{pmatrix}I_5&I_5\\I_5&-I_5\end{pmatrix}.
\]

For \(\chi_t(q)=+1\), this sends \(h(q,+)\) to \((q,0)\) and \(h(q,-)\) to \((0,q)\). For \(\chi_t(q)=-1\), **the assignments swap**. Thus a coordinate permutation taking
\[
(g,q)\longmapsto \big(g\oplus {\bf1}_{\chi_t(q)=-1},\,q\big)
\]
turns the complete two-family matrix into **two independent copies of the original five-row, 16-column matrix**. Crucially, \(\chi_t(q_0)=+1\) for the punctured all-plus singlet \(q_0\), so the punctures preserve the two logical slots.

**Therefore every compatible five-Walsh-character code is permutation-equivalent to two independent \([[15,1,3]]_p\) codes, for all odd primes.**

The implemented verification constructs all 16 cases for \(p=3,5\) and checks **every transformed matrix column exactly** against the direct-sum reference; the mathematical row-transform argument proves this for arbitrary odd primes. The stabilizer-group component sizes are \([15,15]\) in every case.

### A misleading cross-generation signature

For \(t=1\), the five differential characters have masks \((0,3,5,9,14)\). Although this is **exactly a direct sum after reshuffling the particle labels**, its 50 minimum-weight logical Z supports in the *original* family partition are distributed as

| Support uses first original family | Number of three-site logical supports |
|---:|---:|
| 3 coordinates | 3 |
| 2 coordinates | **22** |
| 1 coordinate | **22** |
| 0 coordinates | 3 |

For \(t=3\), the corresponding counts are \(6,19,19,6\). For \(t=0\), all 50 stay inside one original family (25+25). **These are all presentations of the same tensor-product code.**

This exposes a useful interpretation trap: *cross-family logical supports, by themselves, do not establish genuine code coupling*. Family labels can be mixed by a simple coordinate permutation. The original B3/B5 codes still have a separate, stronger **indecomposability** certificate from their stabilizer-support connected components; the present negative example does **not** invalidate that previous proof.

Likewise, for \(t=1\), 132 of the 150 four-site Z-stabilizer-supporting subsets happen to have a 2+2 split across the original generations (9 reside entirely in each original family). Their total reduced-projector purity multiset nevertheless matches the ordinary direct sum. Thus the four-site purity multiset alone is **not a complete equivalence invariant**.

## 5. Implications

1. **Rigidity is now proved**, rather than checked at a few primes: the one-generation five-Cartan-row lattice cannot be extended by any new cubic-orthogonal function on those same 16 points.
2. **Naive nonlinear character extensions are exhausted.** With two families, five pure differential characters are necessary to remove twin-pair errors, and all admissible five-character constructions are trivial direct sums after reassigning species coordinates.
3. **The successful B3/B5 construction is meaningfully different.** It uses *linear combinations* of the five gauge characters, not distinct Walsh characters, and its stabilizer indecomposability cannot be removed by coordinate relabelling.
4. For more interesting **nonlinear** coupled lattices, allow functions with mixed Walsh spectra or mixed-family multipliers beyond a single antisymmetric character; require **integer** cubic cancellation if interpreting rows as continuous-anomaly charges. Exclude trivial code factorizations **before** interpreting cross-family support counts as novel.
5. These character extensions are not exact Standard Model charges: they are generically nonconstant within non-Abelian gauge multiplets, and the earlier Yukawa obstruction still applies.

### Reproduce

\`\`\`bash
python -m pip install numpy
python walsh_character_lattice.py
\`\`\`

The program independently verifies the 16×16 Walsh orthogonality table, the eleven-dimensional pair-product span, the complete 16-mask admissibility classification, the weight-two twin-pair formula, integral cubic identities, 50 logical support patterns, code parameters, and explicit direct-sum equivalence for all 16 maximal choices at \(p=3,5\).

**Scope:** This is an exhaustive theorem for pure *character rows* of the form \((\chi_s,-\chi_s)\) and a single-generation theorem for functions added to the full Cartan lattice. It does not classify arbitrary nonlinear rows or all possible odd-prime triorthogonal codes.
