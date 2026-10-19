---
title: "Beyond One Generation"
subtitle: "Family-Coupled Quantum Codes from Standard Model Charges"
author: "Benjamin John Schulz"
date: "October 2026"
abstract: |
  Prakash and Ramya observed that cubic anomaly cancellation for one Standard Model generation, completed by a right-handed neutrino, is the triorthogonality condition that defines a [[15,1,3]]~p~ quantum code with a transversal non-Clifford gate. We ask what happens when the three generations are combined. Because the generations carry identical gauge charges, stacking them yields one logical qudit, and a family-number charge yields distance two. Family-coupled charge lattices give [[30,2,3]]~p~ and [[45,3,3]]~p~ triorthogonal codes; one such lattice uses only charges that commute with SU(3)×SU(2), namely family-dependent $B-L$ and right-handed isospin, and works for every prime $p\geq 5$. We prove that no new charge is compatible with a single generation over any odd prime, and that every lattice built from family-wise Cartan charges has distance at most three. Exact logical-channel enumerators show that coupling leaves the leading distillation error unchanged but produces correlated two-output errors at order $\varepsilon^4$, where independent codes first produce them at $\varepsilon^6$; we derive the routing rule that keeps these correlations out of a second distillation round. The coupled codes are not equivalent to any tensor product of smaller codes under local unitaries and qudit permutations. Nature does not supply the required charges: with one Higgs doublet, full-rank quark masses and non-block-diagonal mixing force every such family charge to be family-universal, and the charges needed for distance three act on lepton doublets, which neutrino mixing breaks at order one. The same structure explains, in coding language, why $E_8$ embeddings fail: anomaly cancellation by mirror fermions reduces the code distance to one.
---

# 1. Introduction

Quantum error correction has become a useful language for fundamental physics. Holographic codes made precise the sense in which a bulk region of anti-de Sitter space is encoded redundantly in its boundary [1, 2]. The Hayden–Preskill analysis of black-hole information is, at heart, a statement about how fast a scrambled code can be decoded [3]. In the other direction, Kitaev's toric code is a $\mathbb Z_2$ lattice gauge theory, and its protection comes from the structure of gauge constraints [4, 5].

Prakash and Ramya have added a new and unexpected entry to this dictionary [6]. Take one generation of Standard Model fermions together with a right-handed neutrino, sixteen Weyl fields in all, and arrange their charges under the Cartan generators of SO(10) as a matrix. Cubic anomaly cancellation says that the entrywise product of any three rows sums to zero. In coding theory this is triorthogonality, the property that lets a CSS code implement a non-Clifford gate transversally [7]. Puncturing the right-handed neutrino produces a [[15,1,3]]~p~ code for every prime $p$, the qudit analogue of the 15-qubit Reed–Muller code used in magic-state distillation [8, 9].

A single generation is only part of the story. Nature repeats the generation three times, and the repetition is one of the oldest unexplained features of the Standard Model. This paper asks whether the three generations, taken together, define a larger code, and what that code can or cannot do. The answer has three parts.

1. **Constructions.** With suitable family-dependent charges the generations couple into [[30,2,3]]~p~ and [[45,3,3]]~p~ triorthogonal codes. One version uses only charges that commute with the Standard Model gauge group.
2. **Limits.** A single generation admits no further charges (Theorem 1), and no lattice of this type exceeds distance three (Theorem 2). The coupling produces correlated logical errors but no new transversal gates.
3. **Physics.** The family charges that make the construction work are precisely those that the observed masses and mixings destroy (Proposition 3). The obstruction has the same structure as the one that defeats $E_8$ unification.

All results are exact. Codes are computed in finite-field arithmetic, anomaly sums are checked over the integers, and every number quoted here is asserted by a script in the accompanying repository (Appendix A).

**Outline.** Section 2 places the work in context. Section 3 sets up the correspondence between charge lattices and codes. Section 4 describes the obstruction to combining generations, and Section 5 resolves it. Section 6 restricts the construction to gauge-compatible charges. Section 7 proves the limits. Section 8 computes the logical channels and shows the coupled codes are genuinely new. Section 9 asks whether nature supplies the required charges. Section 10 summarizes. Figure 1 gives the route in one picture.

![**Figure 1.** The sequence of constructions. Each fixes the failure of the one before it; the final step is imposed by the observed flavour structure.](figures/fig2_path)

# 2. Context

**Codes and fundamental physics.** Error-correcting structure has appeared in quantum gravity [1, 2, 3] and in lattice gauge theory, where Gauss-law constraints have been used directly for error detection [10]. The present work differs in its object. The code here is not built from a spacetime geometry or a Hilbert space of gauge fields; it is read off from the charge table of the matter content.

**Anomalies and transversal gates.** A growing literature connects transversal non-Clifford gates to symmetry-protected topological phases and to higher-group symmetries and their anomalies in gauge theories [11, 12, 13]. The Prakash–Ramya correspondence is a direct, elementary instance of the same theme: the cubic anomaly is exactly what makes the transversal phase gate act logically (Section 3.3). The codes considered here are not geometrically local, so the Bravyi–König restrictions on transversal gates in local codes [14] do not apply.

**Triorthogonal and qudit distillation codes.** Triorthogonal codes were introduced by Bravyi and Haah [7]; small qubit examples are classified in [15], with an extended classification through length 54 in [16]. Prime-dimension distillation was developed in [9] and refined in [17, 18, 19]; cyclic triorthogonal codes in prime dimension are treated in [20]. We compare our codes with the published qudit families in Section 8.3.

**Flavour.** Hierarchical Yukawa couplings from broken family symmetries go back to Froggatt and Nielsen [21]; family-dependent gauged hypercharge has been proposed as an origin of flavour [22]. Section 9 uses these ideas to ask how badly the charges our codes need are broken.

**$E_8$.** Lisi's proposal places one generation inside $E_8$ [23]; Distler and Garibaldi proved that no such embedding of gravity and the Standard Model is chiral [24]. Section 9.4 translates their obstruction into a statement about code distance.

# 3. Charge lattices as codes

## 3.1 Notation

| Symbol | Meaning |
|:--|:--|
| $\mathcal Q$ | the 16 SO(10) spinor weights, Eq. (1) |
| $q_0=(1,1,1,1,1)$ | the weight of $\nu^c$, the conjugate right-handed neutrino |
| $G$ | integer charge matrix; rows are charges, column $h_j$ belongs to field $j$ |
| $R$, $P$ | punctured (logical) and physical columns; $G_R$, $G_P$ the corresponding submatrices |
| $p$ | an odd prime; all code arithmetic is in $\mathbb F_p$, $\omega=e^{2\pi i/p}$ |
| $[[n,k,d]]_p$ | $n$ physical qudits, $k$ logical qudits, distance $d$ |
| $s$, $f_g$ | family sign (two families) or family label (three families) |
| $B$ | integer matrix of differential charge rows |

: **Table 1.** Notation.

## 3.2 The code

The sixteen fields of one generation are labelled by sign vectors
$$
\mathcal Q=\big\{\,q\in\{+1,-1\}^5:\ \text{the number of }-1\text{ entries is even}\,\big\}. \tag{1}
$$
The first three coordinates carry colour and the last two weak isospin. The column of field $q$ in the charge matrix is $q$ itself (Figure 2a).

**Definition 1.** Let $G$ be an integer matrix with rows of length $N$ and let $R$ be a set of columns. The *charge-lattice code* $\mathcal C(G,R;p)$ is the CSS code on the qudits $P=\{1,\dots,N\}\setminus R$ with
$$
\text{X checks } X(aG_P)\ \text{for}\ aG_R=0,\qquad \text{Z checks } Z(z)\ \text{for}\ G_Pz=0, \tag{2}
$$
and logical basis states
$$
|c\rangle\ \propto \sum_{a\,:\,aG_R=c}|\,aG_P\,\rangle,\qquad c\in\mathbb F_p^{k},\quad k=\operatorname{rank}G_R. \tag{3}
$$

A $Z$-type error $e$ passes all X checks if and only if $G_Pe\in\operatorname{span}G_R$; it acts as the logical operator $Z^c$ when $G_Pe=G_Rc$, and as a stabilizer when $c=0$. The code distance is the smaller of the minimum weights of nontrivial $Z$- and $X$-type logical operators.

Equation (3) has a direct interpretation. Each basis state records, for every physical field, its charge under one combination $a$ of the generators, and the logical label $c$ is the charge carried by the punctured neutrinos. After a Fourier transform on every qudit, the X checks become operators $\exp(2\pi i\,a\cdot Q/p)$ that measure total charge modulo $p$. The code space is therefore a direct sum of charge sectors, and the logical label identifies the sector. We return to this in Section 9.3.

![**Figure 2.** (a) The charge matrix of one generation. The highlighted column, $\nu^c$, is punctured to form the logical qudit. (b) The two-family lattice with $B=B_5$: five shared rows and two differential rows whose sign depends on the family.](figures/fig1_charge_matrix){width=92%}

## 3.3 Where the anomaly enters

Apply $D=\omega^{x^3}$ to every physical qudit, $p\geq5$. On $|aG_P\rangle$ it produces the phase $\omega^{\Phi(a)}$ with $\Phi(a)=\sum_{j\in P}(a\cdot h_j)^3$. Expanding the cube gives the triple sums $\sum_j h_{ij}h_{kj}h_{lj}$. Cubic anomaly cancellation states that these vanish when summed over all sixteen fields, so
$$
\sum_{j\in P}(a\cdot h_j)^3 \;=\; -\sum_{j\in R}(a\cdot h_j)^3 \;=\; -\sum_i c_i^3 . \tag{4}
$$
The phase depends only on the logical label, and the transversal gate acts as $D^{-1}$ on each logical qudit. The anomaly is precisely the statement that the physical phases assemble into a logical one. For $p=3$ the same argument holds with ninth roots of unity because the cancellation is exact over the integers [6].

For one generation, $\mathcal C$ is a [[15,1,3]]~p~ code for every odd prime, with 4 independent X checks, 10 Z checks, X distance 7, and exactly 25 three-site supports of nontrivial $Z$ logicals.

# 4. The twins problem

The second generation alone gives the same code. Mass does not enter the charge table, and the muon family's columns coincide with the electron family's. To obtain two logical qudits we place two generations side by side: 32 columns, two punctured neutrinos, 30 physical qudits. Each field then has a twin in the other family with an identical column, and so do the two logical slots.

Two simple choices fail. With shared charges only, the punctured columns coincide and $k=1$: the code is [[30,1,3]]~p~. With an added family-number charge, $+1$ on one family and $-1$ on the other, the punctured columns separate and $k=2$. But a twin pair carrying $Z$ and $Z^{-1}$ has charge $(0,\dots,0,2)$, which is exactly the difference of the logical columns, so every twin pair is a weight-two logical operator and the code is [[30,2,2]]~p~.

The general statement rests on a property of the weights.

**Lemma 1.** For $p$ odd, $q_0$ is neither a multiple of one non-singlet weight nor a linear combination of two.

*Proof.* If $q=\lambda q_0$, every entry of $q$ has the same sign; $q=-q_0$ has odd parity and is excluded. Suppose $\alpha q_a+\beta q_b=q_0$ with $q_a\neq q_b$. On coordinates where $q_a,q_b$ agree with common sign $\sigma$, $\sigma(\alpha+\beta)=1$; where they disagree, $\tau(\alpha-\beta)=1$ with $\tau$ the sign of $q_a$. Hence $\sigma$ is constant on the agreeing set $A$ and $\tau$ on the disagreeing set $D$. Since both vectors have even parity, $|D|$ is even and $|A|$ is odd, so parity of $q_a$ forces $\sigma=1$. Then $\alpha+\beta=1$ and $\alpha-\beta=\pm1$, giving $(\alpha,\beta)=(1,0)$ or $(0,1)$. $\square$

Consider a two-family lattice in which field $q$ of family $s=\pm1$ has column $(q,\,sBq)$, with $B$ an integer $m\times 5$ matrix. By Lemma 1, a nontrivial $Z$ logical of weight at most two must cancel in the shared rows, which forces a twin pair with opposite coefficients. Its charge is $(0,2Bq)$, while the logical columns differ by $(0,2Bq_0)$. This proves:

**Lemma 2 (twin criterion).** Let $Bq_0\neq0$. The lattice $(q,sBq)$ defines a code of distance three if and only if, for every non-singlet weight $q$, the vector $Bq$ is either zero or not a scalar multiple of $Bq_0$.

(Distance at least three follows from the lemma; distance exactly three follows from Theorem 2 below.) A family-number charge fails because it is the same for every field. The extra charges must depend on the field, not only on the family.

# 5. Family-coupled lattices

## 5.1 Two generations

Lemma 2 reduces the construction to a small combinatorial problem. Two solutions are
$$
B_3=\begin{pmatrix}1&0&0&0&0\\0&1&0&0&1\\0&0&1&1&0\end{pmatrix},
\qquad
B_5=\begin{pmatrix}1&0&0&0&0\\0&1&1&1&1\end{pmatrix}. \tag{5}
$$
For $B_5$, $B_5q_0=(1,4)$. If $B_5q=\lambda(1,4)$ then $\lambda=q_1=\pm1$. With $k$ minus signs among $q_2,\dots,q_5$ the second entry is $4-2k$, so $q_1=1$ requires $2k\equiv0$ and $q_1=-1$ requires $2(4-k)\equiv0$. For $p\geq5$ these give $k=0$ or $k=4$, that is $q_0$ or the excluded all-minus vector. Hence $B_5$ works for every $p\geq5$. At $p=3$ the argument fails, and $B_3$ satisfies the criterion by a similar check.

| Lattice | $p$ | $n$ | $k$ | $d$ | X checks | Z checks | X dist. |
|:-------------------------------|:-----:|:-------:|:-----:|:--:|:------:|:-------:|:------:|
| One generation | odd | 15 | 1 | 3 | 4 | 10 | 7 |
| Shared charges | 5 | 30 | 1 | 3 | 4 | 25 | 14 |
| Family number | 5 | 30 | 2 | 2 | 4 | 24 | 14 |
| Coupled, $B_3$ | 3 | 30 | 2 | 3 | 6 | 22 | 7 |
| Coupled, $B_5$ | 5, 7 | 30 | 2 | 3 | 5 | 23 | 9 |
| Two independent copies | odd | 30 | 2 | 3 | 8 | 20 | 7 |
| Coupled, three families, $B_3$ | 3 | 45 | 3 | 3 | 8 | 34 | 7 |
| Coupled, three families, $B_5$ | 5 | 45 | 3 | 3 | 6 | 36 | 9 |
| Gauge-compatible, Eq. (7) | $\geq5$ | 30 / 45 | 2 / 3 | 3 | 5 / 6 | 23 / 36 | 7 |

: **Table 2.** The codes. "X dist." is the minimum weight of an X-type logical; the code distance $d$ is the smaller of the two distances. The gauge-compatible row is verified at $p=5,7$ (and $p=11$ for $d$); Proposition 1 proves $d=3$ for all $p\geq5$.

Both choices are minimal within this architecture. A single differential row never suffices: every nonzero $b\cdot q$ is then a multiple of $b\cdot q_0$, and $b\cdot q$ cannot vanish on all fifteen non-singlet weights, which span $\mathbb F_p^5$. For $p=3$, an exhaustive search of all 1,080 two-dimensional row spaces over $\mathbb F_3$ that separate the logical slots finds a bad twin pair in each, so three rows are necessary.

## 5.2 Three generations

For three families we use the labels $f_1=(1,0)$, $f_2=(-1,0)$, $f_3=(0,1)$ and give family $g$ the column $(q,\,f_{g,1}Bq,\,f_{g,2}Bq)$. The labels are affinely independent modulo every odd prime, since the determinant of their differences is $\pm2$, and the twin criterion becomes the requirement that $(f_g-f_h)\otimes Bq$ avoid the span of the logical columns. The same matrices satisfy it, giving [[45,3,3]]~p~. All 75 three-site logical supports lie within single families.

# 6. Gauge-compatible charges

The differential rows of Section 5 were chosen freely. A family-dependent charge is a symmetry of the Standard Model only if it commutes with the gauge interactions, which act identically on all families. Within each family the charge must therefore be constant on every colour triplet and every weak doublet. Solving these linear conditions over $\mathbb F_p$ (checked for $p=3,5,7,101$; the solution space is the same as over $\mathbb Q$) gives a two-dimensional space per family, spanned by
$$
S_3=(1,1,1,0,0)\cdot q=3(B-L),\qquad S_2=(0,0,0,1,1)\cdot q=-4\,T_{3R}. \tag{6}
$$
Hypercharge is $Y=S_3/6-S_2/4$. Every row of $B_3$ and $B_5$ assigns different values to members of the quark doublet, so those lattices are not invariant under colour SU(3).

Take instead
$$
B_{\mathrm{SM}}=\begin{pmatrix}1&1&1&0&0\\0&0&0&1&1\end{pmatrix},\qquad B_{\mathrm{SM}}\,q=(S_3,S_2). \tag{7}
$$

| Multiplet | $Q$ | $u^c$ | $d^c$ | $L$ | $e^c$ | $\nu^c$ |
|:--|:-:|:-:|:-:|:-:|:-:|:-:|
| $(S_3,S_2)$ | $(1,0)$ | $(-1,2)$ | $(-1,-2)$ | $(-3,0)$ | $(3,-2)$ | $(3,2)$ |

: **Table 3.** The gauge-compatible family charges on each multiplet.

**Proposition 1.** For every prime $p\geq5$, the lattices built from $B_{\mathrm{SM}}$ define [[30,2,3]]~p~ and [[45,3,3]]~p~ triorthogonal codes. For $p=3$ the distance is two.

*Proof.* By Lemma 2 we need $(S_3,S_2)\neq\lambda(3,2)$ for $\lambda\neq0$ on every non-singlet multiplet. From Table 3: $Q$ requires $1\equiv0$; $u^c$ requires $\lambda=1$ and $4\equiv0$; $d^c$ requires $\lambda=-1$ and $2\equiv0$; $L$ requires $\lambda=0$; $e^c$ requires $\lambda=-1$ and $6\equiv0$. None holds for $p\geq5$. At $p=3$ the $e^c$ condition holds. Triorthogonality holds because each family's rows are linear in $q$ and each family's sixteen weights cancel separately. $\square$

The X distance is 7 at $p=5$ and $p=7$. Both rows are necessary: with $S_3$ alone the lepton-doublet twins are logical, and with $S_2$ alone the right-handed twins are. A further consequence is that gauging a separate hypercharge for each family [22] does not by itself separate the logical qudits, because $\nu^c$ has zero hypercharge; those lattices give $k=1$.

# 7. Limits

## 7.1 Rigidity of a single generation

Because $q_5=q_1q_2q_3q_4$, the set $\mathcal Q$ is the group $(\mathbb Z_2)^4$. Its characters $\chi_s(q)=\prod_{i\le4}q_i^{s_i}$, $s\in\mathbb F_2^4$, satisfy
$$
\sum_{q\in\mathcal Q}\chi_s(q)\chi_t(q)=16\,\delta_{st}, \tag{8}
$$
and 16 is invertible modulo every odd prime. The five SO(10) charges are the characters with masks $A=\{1,2,4,8,15\}$.

**Theorem 1 (rigidity).** Let $p$ be an odd prime. Any function $y:\mathcal Q\to\mathbb F_p$ that can be adjoined as a row to the five SO(10) charges of one generation while preserving triorthogonality lies in their span.

*Proof.* Triorthogonality with two existing rows requires $\sum_q y(q)\,q_iq_j=0$ for all $i\le j$. For $i=j$ the product is the trivial character; for $i<j$ it is $\chi_{a_i\oplus a_j}$. The ten masks $a_i\oplus a_j$ are precisely the masks outside $\{0\}\cup A$. Thus $y$ is orthogonal to eleven of the sixteen characters and lies in the span of the remaining five. $\square$

## 7.2 A distance ceiling

The weights satisfy an integer identity,
$$
q(u^c_3)+q(d_3)-q(e)=q_0 . \tag{9}
$$

**Theorem 2 (distance ceiling).** Suppose each family's columns are $M_gq$ for some integer matrix $M_g$, and each punctured column $M_gq_0$ is nonzero, with the punctured columns linearly independent. Then the code distance is at most three.

*Proof.* Applying $M_g$ to Eq. (9) shows that the error $Z$, $Z$, $Z^{-1}$ on $u^c_3$, $d_3$, $e$ of family $g$ has charge $M_gq_0$, which is a nonzero element of $\operatorname{span}G_R$. $\square$

The hypothesis covers every lattice in Table 2 and any lattice of family-dependent Cartan charges.

## 7.3 Two further observations

**Apparent coupling can be illusory.** Consider nonlinear differential rows $(\chi_s,-\chi_s)$. The triorthogonality conditions reduce to $s_a\oplus s_b\notin A$, and protecting every twin pair requires five rows. Exactly sixteen sets qualify, $\{a\oplus t:a\in A\}$ for $t\in\mathbb F_2^4$. Each is two independent copies of the one-generation code after a change of variables: replacing shared and differential rows by their half-sum and half-difference and exchanging each field with its twin whenever $\chi_t(q)=-1$ yields a direct sum. Nevertheless, 44 of the 50 minimum-weight logical supports cross the original family boundary. Cross-family support counts therefore do not establish coupling; Section 8.2 gives an invariant that does.

**Rows that cancel only modulo $p$.** If one allows a row that is triorthogonal modulo $p$ but not over the integers, as for a discrete $\mathbb Z_p$ charge whose anomaly conditions hold only modulo $p$ [25], new options appear. Extending $B_5$, there are none at $p=5$ or $p=11$; there are 588 rows in 17 orbits at $p=7$ and 4,056 rows in 57 orbits at $p=13$. All give distance three, but with 8 rather than 50 three-site logical supports. At $p=7$ the leading distillation error falls from $\tfrac{25}{18}\varepsilon^3$ to $\tfrac{2}{9}\varepsilon^3$, a factor of 6.25, and the X distance is 9. These rows are not invariant under the gauge group and do not correspond to consistent $U(1)$ charges.

# 8. What the coupling does

## 8.1 Exact logical channels

We use the noise model of [6]: each input qudit independently suffers $Z^j$, $j=1,\dots,p-1$, each with probability $\varepsilon/(p-1)$. The X checks are measured and the output is kept only if all pass. For any target charge $s$,
$$
\sum_{e\,:\,G_Pe=s}x^{\mathrm{wt}(e)}=p^{-r}\sum_{a\in\mathbb F_p^r}\omega^{-a\cdot s}\big(1+(p-1)x\big)^{n-w(a)}(1-x)^{w(a)},\qquad w(a)=\mathrm{wt}(aG_P). \tag{10}
$$
With $s=G_Rc$ one has $a\cdot s=b\cdot c$ for $b=aG_R$. Grouping the $p^r$ row combinations by $(w(a),b)$ gives integer counts, and the character sum collapses because $a\mapsto ua$ permutes the nonzero classes. The result is the exact weight distribution of every logical class. It reproduces $35\varepsilon^3$ for the 15-qubit code and $6.25\varepsilon^3$ per output for the qutrit code of [6], and it agrees with brute-force enumeration of all errors of weight at most four.

| Code | one output in error | both outputs in error | per-output infidelity at $\varepsilon=10^{-3}$ |
|:--|:--|:--|:--|
| Coupled $B_3$, $p=3$ | $\tfrac{25}{2}\varepsilon^3+O(\varepsilon^4)$ | $5\,\varepsilon^4+O(\varepsilon^5)$ | $6.284\times10^{-9}$ |
| Two copies, $p=3$ | $\tfrac{25}{2}\varepsilon^3+O(\varepsilon^4)$ | $\tfrac{625}{16}\varepsilon^6+O(\varepsilon^7)$ | $6.276\times10^{-9}$ |
| Coupled $B_5$, $p=5$ | $\tfrac{25}{8}\varepsilon^3+O(\varepsilon^4)$ | $\tfrac92\,\varepsilon^4+O(\varepsilon^5)$ | $1.573\times10^{-9}$ |
| Two copies, $p=5$ | $\tfrac{25}{8}\varepsilon^3+O(\varepsilon^4)$ | $\tfrac{625}{256}\varepsilon^6+O(\varepsilon^7)$ | $1.568\times10^{-9}$ |
| Coupled, three families, $p=3$ | $\tfrac{75}{4}\varepsilon^3+O(\varepsilon^4)$ | $15\,\varepsilon^4+O(\varepsilon^5)$ | $6.292\times10^{-9}$ |
| Coupled, three families, $p=5$ | $\tfrac{75}{16}\varepsilon^3+O(\varepsilon^4)$ | $\tfrac{27}{2}\varepsilon^4+O(\varepsilon^5)$ | $1.579\times10^{-9}$ |

: **Table 4.** Leading terms of the exact logical channels, conditioned on acceptance and summed over outputs.

![**Figure 3.** Logical error after post-selection for two generations at $p=5$. The single-output curves of the coupled and independent codes coincide to the resolution of the plot; the correlated curves differ by two powers of $\varepsilon$.](figures/fig3_channels){width=78%}

Through order $\varepsilon^3$ the coupled codes and independent copies have the same acceptance probability and the same single-output error. They differ in correlation. Errors that corrupt both outputs appear at order $\varepsilon^4$ in the coupled code and at order $\varepsilon^6$ for independent copies. They arise entirely from weight-four patterns with two sites in each of two families: at $p=5$ there are 1,584 such logical patterns, of which 1,152 corrupt both outputs, and at $p=3$ there are 176, of which 80 do. Independent copies have none. At $\varepsilon=10^{-3}$ the per-output infidelity is 0.1% to 0.7% higher than for independent copies, and correlated failures are $10^5$ to $10^6$ times more likely, though still below $10^{-11}$.

## 8.2 Correlated outputs in a distillation factory

Fault-tolerance thresholds are usually proven under assumptions that limit correlations among faults [26, 27]. The correlations found here are of a different kind. They are not correlations in the physical noise; the code itself produces them in its outputs. Whether they matter depends on how outputs are consumed.

Consider two rounds of distillation in which each first-round output has error probability $\delta\approx c\,\varepsilon^3$ and the second round uses a code with distance three. If the inputs to a second-round block are independent, its leading failure probability is $O(\delta^3)=O(\varepsilon^9)$. Suppose instead that two of its inputs come from the same coupled first-round block and occupy two sites of a three-site logical support of the second-round code. Both are wrong with probability $O(\varepsilon^4)$, and one further error on the third site, with probability $O(\varepsilon^3)$, completes the logical. The failure probability becomes $O(\varepsilon^7)$, and two orders of suppression are lost.

The remedy is a routing rule: **outputs of the same coupled block must feed different second-round blocks.** With that rule the correlations do not enter at leading order, and the coupled codes behave like independent copies in a two-round protocol. A circuit-level analysis with noisy checks of weight 15 and 16 would be needed to settle the overhead, and we have not done it.

The coupled codes are also more expensive to measure. At $p=5$ a minimum-weight check basis has X checks of weight 15 and 16 against 8 for independent copies, for a total check weight of 170 against 144.

## 8.3 The coupled codes are new

Section 7.3 showed that a code can appear coupled while being a direct sum. The following criterion is invariant under any local change of basis.

**Proposition 2.** A stabilizer code with stabilizer group $S$ factorizes as $\mathcal C_A\otimes\mathcal C_B$ across a partition $A|B$ of its qudits if and only if $S=S_A\times S_B$, where $S_A$ denotes the subgroup supported on $A$.

*Proof.* The code projector is $\Pi=|S|^{-1}\sum_{s\in S}s$. If $\Pi=\Pi_A\otimes\Pi_B$, its Pauli expansion is supported on a product set $T_A\times T_B$. The identity lies in $T_B$, so $(t_A,1)\in S$ for all $t_A\in T_A$, and $S=S_A\times S_B$. The converse is immediate. $\square$

Factorizability across some partition is preserved by local unitaries and qudit permutations. For a CSS code, the partitions across which $S$ splits are unions of connected components of two linear matroids, the X-check space and the row space of $G_P$, and these components can be computed exactly. Every coupled code in Table 2, for two and three families at every prime tested, has a single component; a direct sum has one component per block.

**Corollary.** No coupled code in Table 2 is equivalent, under local unitaries and qudit permutations, to independent one-generation codes or to any tensor product of smaller codes.

An independent invariant gives the same conclusion. For four qudits $A$ the reduced code state has purity $\operatorname{Tr}\rho_A^2=p^{\dim S_A-4}$. The number of four-qudit sets with enhanced purity is 150 for two independent copies, against 184 ($p=3$) and 249 ($p=5$) for the coupled codes; for three families it is 225 against 327 and 522.

We have also excluded two published qudit families. The Prakash–Saha qutrit codes [[9m−k,k,2]]~3~ [18] have distance two and cannot match $(n,k)=(30,2)$ or $(45,3)$. Full-evaluation punctured Reed–Muller codes [19] require $n+k=p^m$, whereas here $n+k$ is 32 or 48. A systematic comparison with all known qudit triorthogonal codes remains to be done.

## 8.4 No new transversal gates

Consider all transversal diagonal gates $\bigotimes_jD^{w_j}$ with $D=\omega^{x^3}$ or $\omega^{x^2}$ and arbitrary powers $w_j\in\mathbb F_p$, $p\geq5$. Such a gate preserves the code if and only if
$$
\sum_jw_j\,(x\cdot h_j)(y\cdot h_j)(\sigma\cdot h_j)=0\qquad\text{for all }x,y\in\mathbb F_p^r\text{ and all }\sigma\text{ with }\sigma G_R=0, \tag{11}
$$
or the corresponding quadratic condition. The constraint is linear in $w$, so the full space of such gates is computable. The induced logical gate is the restriction of a cubic or quadratic form to logical labels, and it is entangling exactly when that form has mixed terms. For every code considered here (coupled, gauge-compatible, mod-7, direct sum, and family number) the mixed part vanishes. The coupling is visible in the errors and invisible in the gates. Non-diagonal and permutation-assisted transversal gates are not covered by this analysis.

# 9. Does nature supply the charges?

## 9.1 Exact family charges and the quark masses

**Proposition 3.** Consider $G\geq2$ quark families, one Higgs doublet, and an unbroken family charge in the gauge-compatible span, $a_iS_3+b_iS_2$ on family $i$. If the up and down mass matrices have full rank and the quark mixing matrix admits no nontrivial common block decomposition, then $a_i$ and $b_i$ are independent of $i$.

*Proof.* The doublet $Q_i$ carries charge $a_i$, and $u^c_j$, $d^c_j$ carry $-a_j+2b_j$ and $-a_j-2b_j$. An up-type Yukawa entry $(i,j)$ is allowed only if $a_i-a_j+2b_j+\eta_u=0$, where $\eta_u$ is the Higgs charge. Invariance implies that $Y_uY_u^\dagger$ and $Y_dY_d^\dagger$ commute with the charge on the doublets. Distinct $a_i$ would then make both block diagonal in a common basis, and the mixing matrix with them; hence all $a_i$ are equal. Each column $j$ of $Y_u$ is then entirely allowed or entirely forbidden according to $2b_j+\eta_u$. Full rank forbids a zero column, and 2 is invertible, so all $b_j$ are equal. $\square$

A deterministic enumeration of every charge assignment modulo $p$, for $G=2$ at $p=3,5,7$ and $G=3$ at $p=3,5$, confirms that exactly $p^2$ assignments survive, all family-universal. Family-universal charges make the punctured columns coincide, so $k=1$. With exact family charges, three generations encode one logical qudit.

The proposition is the familiar statement that generic Yukawa couplings leave no family-dependent abelian symmetry. Its content here is the consequence for the code. It concerns unbroken symmetries in the one-Higgs quark sector; models with several Higgs doublets, spurions, or non-abelian flavour groups require separate treatment.

## 9.2 How badly are the charges broken?

Exact invariance is excluded, so the relevant question is the size of the breaking, in the Froggatt–Nielsen sense [21]. A family charge acting on the doublets forces the mixing matrix to be block diagonal; its breaking is bounded below by the largest mixing element that no choice of basis removes.

| Doublet charge separates | Quark doublets (CKM) | Lepton doublets (PMNS) |
|:--|:-:|:-:|
| family 3 from families 1 and 2 | 0.041 | 0.66 |
| family 1 from families 2 and 3 | 0.23 | 0.55 |

: **Table 5.** Smallest symmetry breaking compatible with measured mixing, from the CKM global fit [28] and NuFIT 6.0 [29].

A charge acting only on right-handed fields is unaffected by mixing. With one Higgs doublet it allows one Yukawa column per sector, and its breaking is set by mass ratios: $m_c/m_t=0.004$, $m_s/m_b=0.019$ and $m_\mu/m_\tau=0.059$ at $M_Z$ [30].

These numbers decide the question. The right-handed charge $S_2$ is broken only at the percent level, but by Section 6 it gives distance two. Distance three requires $S_3$ as well, and $S_3$ acts on the lepton doublets with value $-3$, where neutrino mixing breaks it at order one. **No distance-three code in this family is protected by a small symmetry-breaking parameter in the observed flavour sector.**

## 9.3 Superselection

A second obstruction is independent of flavour. In the conjugate basis of Section 3.2, logical states lie in different sectors of total charge. For an exact gauge charge, total charge is fixed by the flux at infinity, and coherent superpositions of different total charge cannot be prepared or detected with local operations [31]. Such superpositions acquire operational meaning only relative to a charged reference system [32, 33]. A logical qudit of this code would then be stored in the relation between the matter and the reference, not in the matter alone. The codes describe the charge table; they are implemented on engineered hardware, not by elementary particles.

## 9.4 $E_8$

Under $E_8\supset SO(10)\times SO(6)$, the 128 spinorial roots of $E_8$ decompose as
$$
\mathbf{128}=(\mathbf{16},\mathbf 4)\oplus(\overline{\mathbf{16}},\bar{\mathbf 4}), \tag{12}
$$
four copies of a generation labelled by SU(4) weights, together with their mirrors. Distler and Garibaldi showed that embeddings of this kind cannot be chiral [24]. In the present language:

* The chiral block $(\mathbf{16},\mathbf 4)$ is not triorthogonal. The triple sum over the three SU(4) Cartan rows equals 64, which fails over the integers and modulo 3, 5 and 7.
* The full $\mathbf{128}$ is triorthogonal because every column appears with its negative. Each punctured column's mirror is then a physical column and a weight-one logical operator, and the code is [[124,4,1]]~p~.
* Using the SU(4) weights as family labels on the chiral block gives the family-number construction, [[60,4,2]]~p~. Cartan charges of a product group act additively, $(x,y)\cdot(q,w)=x\cdot q+y\cdot w$, and cannot produce the products $f_g\,(b\cdot q)$ that distance three requires.

$E_8$ cancels anomalies in the only way available to it, by pairing each field with its mirror, and in coding terms that pairing removes the protection.

# 10. Discussion

The argument can be summarized in one chain. Anomaly cancellation is triorthogonality, so a single generation is a magic-state code. The other generations give the same code because the charge table does not see mass. Combining them requires charges that distinguish the families field by field; the gauge group permits exactly two such charges per family, and with both the families couple into [[30,2,3]]~p~ and [[45,3,3]]~p~ codes for $p\geq5$. A single generation admits nothing further, distance three is a ceiling, and the coupling appears in correlated errors rather than in gates. The observed masses and mixings then remove every exact family charge, and the charge that distance three requires is broken at order one by neutrino mixing.

The result says something about the generations themselves. At the level of gauge charges the three generations are indistinguishable. At the level of masses, where particles propagate, they are mixed. A multi-family code requires a sharp and conserved answer to the question "which family is this?", and neither level provides one. The $E_8$ case reaches the same conclusion from the opposite direction: a large symmetry group cannot supply families that are both chiral and sharply labelled.

We list the open problems.

1. Compare the coupled codes systematically with known qudit triorthogonal codes, using canonical forms rather than invariants.
2. Extend the gate analysis to non-diagonal and permutation-assisted transversal gates, and to the ninth-root gates at $p=3$.
3. Explain why the nonlinear rows of Section 7.3 appear only for $p\equiv1\pmod 3$, and whether the floor of eight weight-three supports can be lowered.
4. Carry out a circuit-level analysis of a two-round factory, with noisy checks and the routing rule of Section 8.2.
5. Determine whether multi-Higgs or non-abelian flavour models can protect a family charge that acts on lepton doublets at a useful level.

# Acknowledgments {.unnumbered}

The exact computations, verification scripts and drafting were carried out with the assistance of Claude (Anthropic). Every claim is reproduced by the scripts listed in Appendix A.

# Appendix A. Reproducibility {.unnumbered}

All scripts are in the repository [github.com/HiroSakuraba/Standard-model-generation-as-quantum-error-correction](https://github.com/HiroSakuraba/Standard-model-generation-as-quantum-error-correction) and assert the numbers they print. They require Python 3.9 or later; all but the first three require NumPy.

| Result | Script | Runtime |
|:------------------------------|:--------------------------------------|:-------|
| One- and two-generation codes, minimality (§4, §5.1) | `matter_qec_study.py` | seconds |
| Three-generation codes and supports (§5.2) | `three_generation_flavor.py` | ~20 s |
| Yukawa enumeration (§9.1) | `exact_yukawa.py` | seconds |
| Gauge-compatible lattice, channels, equivalence, gates, mod-$p$ rows, breaking (§6–§9) | `open_steps.py` | ~4 min |
| Rigidity and disguised direct sums (§7.1, §7.3) | `walsh_character_lattice.py` | seconds |
| Four-site purity (§8.3) | `local_purity.py` | seconds |
| Prakash–Saha family (§8.3) | `published_qutrit_family.py` | seconds |
| Reed–Muller family (§8.3) | `reed_muller_family_exclusion.py` | seconds |
| $E_8$ (§9.4) | `e8_comparison.py` | ~1 s |
| Figures | `paper/make_figures.py` | seconds |

# References {.unnumbered}

1. A. Almheiri, X. Dong and D. Harlow, "Bulk locality and quantum error correction in AdS/CFT," JHEP 04 (2015) 163, [arXiv:1411.7041](https://arxiv.org/abs/1411.7041).
2. F. Pastawski, B. Yoshida, D. Harlow and J. Preskill, "Holographic quantum error-correcting codes: toy models for the bulk/boundary correspondence," JHEP 06 (2015) 149, [arXiv:1503.06237](https://arxiv.org/abs/1503.06237).
3. P. Hayden and J. Preskill, "Black holes as mirrors: quantum information in random subsystems," JHEP 09 (2007) 120, [arXiv:0708.4025](https://arxiv.org/abs/0708.4025).
4. A. Yu. Kitaev, "Fault-tolerant quantum computation by anyons," Ann. Phys. 303, 2 (2003), [arXiv:quant-ph/9707021](https://arxiv.org/abs/quant-ph/9707021).
5. E. Dennis, A. Kitaev, A. Landahl and J. Preskill, "Topological quantum memory," J. Math. Phys. 43, 4452 (2002), [arXiv:quant-ph/0110143](https://arxiv.org/abs/quant-ph/0110143).
6. S. Prakash and K. L. Ramya, "One generation of the Standard Model as a quantum error-correcting code," [arXiv:2610.21879](https://arxiv.org/abs/2610.21879) (2026).
7. S. Bravyi and J. Haah, "Magic-state distillation with low overhead," Phys. Rev. A 86, 052329 (2012), [arXiv:1209.2426](https://arxiv.org/abs/1209.2426).
8. S. Bravyi and A. Kitaev, "Universal quantum computation with ideal Clifford gates and noisy ancillas," Phys. Rev. A 71, 022316 (2005), [arXiv:quant-ph/0403025](https://arxiv.org/abs/quant-ph/0403025).
9. E. T. Campbell, H. Anwar and D. E. Browne, "Magic-state distillation in all prime dimensions using quantum Reed–Muller codes," Phys. Rev. X 2, 041021 (2012), [arXiv:1205.3104](https://arxiv.org/abs/1205.3104).
10. A. Rajput, A. Roggero and N. Wiebe, "Quantum error correction with gauge symmetries," npj Quantum Information (2023), [arXiv:2112.05186](https://arxiv.org/abs/2112.05186).
11. B. Yoshida, "Topological color code and symmetry-protected topological phases," [arXiv:1503.07208](https://arxiv.org/abs/1503.07208) (2015).
12. M. Barkeshli, Y.-A. Chen, P.-S. Hsin and N. Manjunath, "Higher-group symmetry in finite gauge theory and stabilizer codes," SciPost Phys. 16, 089 (2024), [arXiv:2211.11764](https://arxiv.org/abs/2211.11764).
13. P.-S. Hsin and R. Kobayashi, "Automorphism in gauge theories: higher symmetries and transversal non-Clifford logical gates," [arXiv:2511.15783](https://arxiv.org/abs/2511.15783) (2025).
14. S. Bravyi and R. König, "Classification of topologically protected gates for local stabilizer codes," Phys. Rev. Lett. 110, 170503 (2013), [arXiv:1206.1609](https://arxiv.org/abs/1206.1609).
15. S. Nezami and J. Haah, "Classification of small triorthogonal codes," Phys. Rev. A 106, 012437 (2022), [arXiv:2107.09684](https://arxiv.org/abs/2107.09684).
16. A. Wills, S. P. Jain and S. Singh, "Classification of generalised triorthogonal codes through length 54," [arXiv:2609.30860](https://arxiv.org/abs/2609.30860) (2026).
17. A. Krishna and J.-P. Tillich, "Towards low overhead magic state distillation," [arXiv:1811.08461](https://arxiv.org/abs/1811.08461) (2019).
18. S. Prakash and T. Saha, "Low overhead qutrit magic state distillation," Quantum 9, 1768 (2025), [doi:10.22331/q-2025-06-12-1768](https://doi.org/10.22331/q-2025-06-12-1768).
19. T. Saha and S. Prakash, "Sublogarithmic distillation in all prime dimensions using punctured Reed–Muller codes," Phys. Rev. A 113, 042405 (2026), [doi:10.1103/tzkz-9p2f](https://doi.org/10.1103/tzkz-9p2f).
20. S. Prakash, "Cyclic triorthogonal codes in prime dimension," [arXiv:2610.08012](https://arxiv.org/abs/2610.08012) (2026).
21. C. D. Froggatt and H. B. Nielsen, "Hierarchy of quark masses, Cabibbo angles and CP violation," Nucl. Phys. B 147, 277 (1979).
22. M. Fernández Navarro and S. F. King, "Tri-hypercharge: a separate gauged weak hypercharge for each fermion family as the origin of flavour," [arXiv:2305.07690](https://arxiv.org/abs/2305.07690) (2023).
23. A. G. Lisi, "An exceptionally simple theory of everything," [arXiv:0711.0770](https://arxiv.org/abs/0711.0770) (2007).
24. J. Distler and S. Garibaldi, "There is no 'Theory of Everything' inside E8," Commun. Math. Phys. (2010), [arXiv:0905.2658](https://arxiv.org/abs/0905.2658), [doi:10.1007/s00220-010-1006-y](https://doi.org/10.1007/s00220-010-1006-y).
25. L. E. Ibáñez and G. G. Ross, "Discrete gauge symmetry anomalies," Phys. Lett. B 260, 291 (1991).
26. P. Aliferis, D. Gottesman and J. Preskill, "Quantum accuracy threshold for concatenated distance-3 codes," Quantum Inf. Comput. 6, 97 (2006), [arXiv:quant-ph/0504218](https://arxiv.org/abs/quant-ph/0504218).
27. J. Preskill, "Sufficient condition on noise correlations for scalable quantum computing," Quantum Inf. Comput. 13, 181 (2013), [arXiv:1207.6131](https://arxiv.org/abs/1207.6131).
28. Particle Data Group, "CKM quark-mixing matrix," Review of Particle Physics (2020), [review PDF](https://pdg.web.cern.ch/2020/reviews/rpp2020-rev-ckm-matrix.pdf).
29. I. Esteban et al., "NuFit-6.0: updated global analysis of three-flavor neutrino oscillations," [arXiv:2410.05380](https://arxiv.org/abs/2410.05380) (2024).
30. G.-y. Huang and S. Zhou, "Precise values of running quark and lepton masses in the Standard Model," [arXiv:2009.04851](https://arxiv.org/abs/2009.04851) (2020).
31. G. C. Wick, A. S. Wightman and E. P. Wigner, "The intrinsic parity of elementary particles," Phys. Rev. 88, 101 (1952).
32. Y. Aharonov and L. Susskind, "Charge superselection rule," Phys. Rev. 155, 1428 (1967).
33. A. Kitaev, D. Mayers and J. Preskill, "Superselection rules and quantum protocols," Phys. Rev. A 69, 052326 (2004), [arXiv:quant-ph/0310088](https://arxiv.org/abs/quant-ph/0310088).
