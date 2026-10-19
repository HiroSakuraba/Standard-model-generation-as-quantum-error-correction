---
title: "Beyond One Generation"
subtitle: "Family-Coupled Quantum Codes from Standard Model Charges"
author: "Benjamin John Schulz"
date: "October 2026"
abstract: |
  Prakash and Ramya observed that the cubic anomaly-cancellation condition for one Standard Model generation, completed by a right-handed neutrino, is the triorthogonality condition that defines a [[15,1,3]]~p~ quantum code with a transversal non-Clifford gate. This paper asks what happens when the three generations are put together. Because the generations carry identical gauge charges, stacking them gives one logical qudit, and a family-number charge gives distance two. We construct family-coupled charge lattices that give [[30,2,3]]~p~ and [[45,3,3]]~p~ triorthogonal codes. One version uses only charges that commute with SU(3)×SU(2), namely family-dependent $B-L$ and right-handed isospin, and it works for every prime $p\geq 5$. We prove two limits. No new charge row is compatible with a single generation over any odd prime, and every lattice built from family-wise Cartan charges has distance at most three. We compute the logical error channels exactly. Coupling the families leaves the leading distillation error unchanged but produces correlated errors on two outputs at order $\varepsilon^4$, where independent codes first produce them at $\varepsilon^6$. The coupled codes are not equivalent to independent copies under any local unitary and qudit relabelling. Nature does not supply the needed charges. In a one-Higgs quark sector, full-rank masses and non-block-diagonal mixing force every family charge to be family-universal, and the charges that distance three requires sit on lepton doublets, which neutrino mixing breaks at order one. The same reasoning shows why $E_8$ fails: restoring anomaly cancellation with mirror fermions drops the code distance to one.
---

# 1. The question

There is a curious fact about the Standard Model that Shiroman Prakash and K. L. Ramya pointed out this month [1]. Take one generation of matter, the up and down quarks in three colours, the electron, the neutrino, and add a right-handed neutrino. That makes sixteen fields. Write down how each field responds to the forces, which is a table of numbers called charges. For the theory to make sense quantum mechanically, the charges must cancel a certain cubic sum. This is anomaly cancellation, and it is the reason a generation has exactly the fields it has.

Now forget the physics and look at the table as a quantum engineer would. The anomaly condition says that if you take any three rows of the table, multiply them entry by entry, and add up the result, you get zero. Quantum coding theorists have a name for a matrix with that property. They call it triorthogonal, and a triorthogonal matrix defines a quantum error-correcting code with a special property: a certain non-Clifford gate can be applied to the encoded information just by applying the same simple gate to every physical qudit. Codes of this kind are the engine of magic-state distillation [2, 3, 4], which is how a fault-tolerant quantum computer gets the one ingredient it cannot make cheaply.

So one generation of matter, with the right-handed neutrino set aside as the logical slot, is a 15-qudit code that protects one quantum digit. For qubits it is the familiar 15-qubit Reed–Muller code. For a prime $p$ it is a [[15,1,3]]~p~ code.

The obvious question is what happens with the other generations. Nature repeats the generation twice: charm, strange, muon and its neutrino; then top, bottom, tau and its neutrino. Do three generations make a bigger code? This paper answers that question, and the answer has three parts.

1. **Mathematically, yes.** With the right extra charges the families couple into [[30,2,3]]~p~ and [[45,3,3]]~p~ codes. These are not disguised copies of the one-generation code.
2. **There are hard limits.** A single generation admits no extra charges at all, and no lattice of this kind gets past distance three.
3. **Physically, no.** The extra charges that make the coupling work are exactly the ones the observed masses and mixings destroy.

Along the way we will see that the coupling has one visible signature, correlated errors, and that the reason the construction fails in nature is the same reason $E_8$ grand unification fails.

Figure 1 shows the route. Everything here is exact. The codes are computed in finite-field arithmetic, the anomaly sums are checked over the integers, and every number in this paper is asserted by a script in the repository that accompanies it (Appendix A).

![**Figure 1.** The path of this paper. Each construction fixes the failure of the one before it. The last arrow is drawn by nature.](figures/fig2_path)

# 2. One generation as a code

## 2.1 Sixteen fields as sixteen sign patterns

The cleanest way to write the charges of one generation uses SO(10), the group that holds a full generation in a single sixteen-dimensional spinor. Each of the sixteen fields is labelled by a vector of five signs,
$$
q = (q_1,q_2,q_3,q_4,q_5),\qquad q_i=\pm 1,
$$
with an even number of minus signs. There are exactly sixteen such vectors. The first three coordinates carry colour, the last two carry weak isospin. The all-plus vector $q_0=(1,1,1,1,1)$ is the right-handed neutrino, more precisely its conjugate $\nu^c$, which carries no charge under any Standard Model force.

The charge matrix $G$ has five rows, one for each coordinate, and sixteen columns, one for each field. Its column for field $q$ is just $q$ itself. Figure 2(a) shows it.

![**Figure 2.** (a) The charge matrix of one generation. The highlighted column is the right-handed neutrino, which becomes the logical slot. (b) The coupled lattice for two generations: five shared rows, and two differential rows that flip sign between families.](figures/fig1_charge_matrix)

## 2.2 From a charge matrix to a code

Here is the dictionary. Work modulo a prime $p$. Pick some columns to be logical slots; for one generation that is the $\nu^c$ column. Call the matrix restricted to those columns $G_R$ and the matrix on the remaining "physical" columns $G_P$. A row combination $a\in\mathbb F_p^r$ gives a vector $aG_P$ on the physical qudits and a value $aG_R$ on the logical slots.

* **X checks** are the row combinations that vanish on the logical slots: $X(aG_P)$ with $aG_R=0$.
* **Z checks** are the patterns $z$ of physical fields with zero total charge: $G_P z = 0$.
* **Z errors that do damage** are patterns $e$ whose total charge equals the charge of some logical slots: $G_P e = G_R c$ with $c\neq 0$. Such an error applies the logical operator $Z^c$ and passes every X check.
* The **distance** is the size of the smallest damaging error.

The codewords have a simple form,
$$
|c\rangle \;\propto \sum_{a\,:\,aG_R=c} |\,aG_P\,\rangle .
$$
Each basis state lists, for every physical field, its charge under one particular combination $a$ of the forces. The logical label $c$ is the charge that the missing neutrino would have carried. In a Fourier-rotated basis the X checks become measurements of the total charge modulo $p$, so the code space is a set of charge sectors and the logical digit says which sector you are in. Keep this picture in mind. It will matter in Section 8.

## 2.3 Where anomaly cancellation enters

Apply the gate $D=\omega^{x^3}$, with $\omega=e^{2\pi i/p}$, to every physical qudit, for $p\geq 5$. On the basis state $|aG_P\rangle$ it produces the phase
$$
\omega^{\sum_{j\in P}(a\cdot h_j)^3},
$$
where $h_j$ is the column of field $j$. Expand the cube and you get the triple sums $\sum_j h_{ij}h_{kj}h_{lj}$ over all fields. Anomaly cancellation says those sums vanish when taken over *all* sixteen fields. So the sum over the physical fields equals minus the contribution of the logical slots:
$$
\sum_{j\in P}(a\cdot h_j)^3 = -\sum_{j\in R}(a\cdot h_j)^3 = -\sum_i c_i^3 .
$$
The phase depends only on the logical label. The transversal gate acts on the encoded digit as $D^{-1}$, a non-Clifford gate. That is the whole trick: the anomaly is what makes the phase of the physical qudits add up to a phase of the logical one. For $p=3$ the same argument works with ninth roots of unity because the cancellation holds over the integers, and for $p=2$ one uses the bit form of the charges [1].

For one generation the result is [[15,1,3]]~p~ for every odd prime, with 4 independent X checks, 10 Z checks, and exactly 25 three-site supports for damaging errors.

# 3. Two generations and the twins problem

The second generation alone gives exactly the same code. The muon is two hundred times heavier than the electron, but mass does not appear in the charge table. Column for column, the muon family's table is the electron family's table. Only the labels change.

So put two generations side by side: 32 columns, two right-handed neutrinos as logical slots, 30 physical qudits. We would like two logical digits. Here is the problem. Every field has a twin in the other family with an identical column. The two logical slots are themselves twins.

**Shared charges only.** If both families use the same five rows, the two logical columns are identical, so the two slots hold one digit. The code is [[30,1,3]]~p~. Two identical hiding places hold one secret.

**Family number.** Add a sixth row: $+1$ on every field of the first family, $-1$ on the second. It is anomaly free, and now the logical columns differ, so there are two digits. But consider an error that applies $Z$ to a field in the first family and $Z^{-1}$ to its twin. The shared rows cancel. The family row contributes $1-(-1)=2$, which is exactly the difference of the two logical columns. The error is invisible and damaging. Every one of the fifteen twin pairs does this, so the code is [[30,2,2]]~p~.

The general lesson is in a lemma. We use one fact from [1], and we give it a short proof since everything rests on it.

**Lemma 1.** Over $\mathbb F_p$ with $p$ odd, $q_0$ is not a multiple of a single non-singlet weight, and not a combination of two of them.

*Proof.* If $q=\lambda q_0$ then all signs of $q$ agree, so $q=q_0$ or $q=-q_0$, and $-q_0$ has five minus signs, which is odd. For two weights $q_a\neq q_b$, write $\alpha q_a+\beta q_b=q_0$. On coordinates where $q_a$ and $q_b$ agree with sign $\sigma$ we need $\sigma(\alpha+\beta)=1$; where they disagree we need $\tau(\alpha-\beta)=1$, with $\tau$ the sign of $q_a$. So $\sigma$ is constant on the agreeing set $A$ and $\tau$ is constant on the disagreeing set $D$. Two even-parity vectors differ in an even number of places, so $|D|$ is even and $|A|$ is odd. Parity of $q_a$ then forces $\sigma=+1$, so $\alpha+\beta=1$ and $\alpha-\beta=\pm1$, which gives $(\alpha,\beta)=(1,0)$ or $(0,1)$. That is the one-weight case again. $\square$

Now take any two-family lattice with columns $(q,\;s\,Bq)$, where $s=\pm1$ labels the family and $B$ is an $m\times 5$ integer matrix of extra "differential" rows. Lemma 1 says that the shared part of a damaging error with at most two sites must cancel, which only happens for a matched twin pair with opposite coefficients. Such a pair has charge $(0,\,2Bq)$. The difference of the logical columns is $(0,\,2Bq_0)$. So:

**Lemma 2 (twin criterion).** A two-family lattice $(q,sBq)$ with $Bq_0\neq 0$ has distance three if and only if, for every non-singlet weight $q$, $Bq$ is either zero or not a multiple of $Bq_0$.

For family number, $B$ is effectively a constant, so every $Bq$ is a multiple of $Bq_0$, and every twin pair fails. To pass, the extra charges have to depend on *which field* it is, not just *which family*.

# 4. Coupled lattices

## 4.1 Two generations

Lemma 2 turns the problem into a small puzzle: find $B$ so that no non-singlet weight lines up with $q_0$. Here are two solutions.
$$
B_3=\begin{pmatrix}1&0&0&0&0\\0&1&0&0&1\\0&0&1&1&0\end{pmatrix},
\qquad
B_5=\begin{pmatrix}1&0&0&0&0\\0&1&1&1&1\end{pmatrix}.
$$
Take $B_5$. We have $B_5q_0=(1,4)$. If $B_5q=\lambda(1,4)$, the first entry forces $\lambda=q_1=\pm1$. Let $k$ be the number of minus signs among $q_2,\dots,q_5$, so the second entry is $4-2k$. If $q_1=1$ we need $2k\equiv 0$, so $k=0$; if $q_1=-1$ we need $2(4-k)\equiv0$, so $k=4$. For $p\geq 5$ these give $q_0$ and the forbidden all-minus vector. So $B_5$ works for every prime $p\geq5$. At $p=3$, $4\equiv1$ and the argument fails; $B_3$ handles $p=3$ by a similar check.

| Lattice | $p$ | $n$ | $k$ | $d$ | X checks | Z checks | X distance |
|:--|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| One generation | any odd | 15 | 1 | 3 | 4 | 10 | 7 |
| Shared charges | 5 | 30 | 1 | 3 | 4 | 25 | 14 |
| Family number | 5 | 30 | 2 | 2 | 4 | 24 | 14 |
| Coupled, $B_3$ | 3 | 30 | 2 | 3 | 6 | 22 | 7 |
| Coupled, $B_5$ | 5, 7 | 30 | 2 | 3 | 5 | 23 | 9 |
| Two independent copies | any odd | 30 | 2 | 3 | 8 | 20 | 7 |
| Coupled, three families, $B_3$ | 3 | 45 | 3 | 3 | 8 | 34 | 7 |
| Coupled, three families, $B_5$ | 5 | 45 | 3 | 3 | 6 | 36 | 9 |

: **Table 1.** The constructions. "X distance" is the smallest X-type logical; the code distance is the smaller of the two.

These are as small as this architecture allows. One differential row never works: with a single row every nonzero $b\cdot q$ is a multiple of $b\cdot q_0$, and $b\cdot q$ cannot vanish on all fifteen non-singlet weights because they span $\mathbb F_p^5$. So $B_5$ is minimal for $p\geq5$. For $p=3$ we enumerated every two-dimensional row space over $\mathbb F_3$ that separates the logical slots, all 1,080 of them, and each one leaves a bad twin pair. So $B_3$, with three rows, is minimal at $p=3$.

## 4.2 Three generations

For three families, replace the sign $s$ by a two-component label,
$$
f_1=(1,0),\qquad f_2=(-1,0),\qquad f_3=(0,1),
$$
and give family $g$ the column $\big(q,\;f_{g,1}Bq,\;f_{g,2}Bq\big)$. The three labels are affinely independent modulo every odd prime, since the determinant of their differences is $\pm2$. The twin criterion becomes: $(f_g-f_h)\otimes Bq$ must not lie in the span of the logical columns, which holds for the same $B$ matrices. The result is [[45,3,3]]~p~ (Table 1). Every three-site damaging error lies inside a single family, 25 per family, 75 in all.

# 5. Charges the forces allow

We have been treating the extra rows as free choices. They are not. A charge that depends on the family is a symmetry of the Standard Model only if it commutes with the forces, and the gluons and $W$ bosons are shared by all three families. A gluon can turn a red up quark into a green one in the first family. If the first family's red and green quarks carried different extra charge, that process would violate it.

So within each family, an allowed charge must take the same value on every member of a colour triplet and every member of a weak doublet. Solving those linear conditions over $\mathbb F_p$ (we did it for $p=3,5,7,101$; the answer is the same as over $\mathbb Q$) gives exactly a two-dimensional space per family:
$$
S_3 = (1,1,1,0,0)\cdot q = 3(B-L),\qquad S_2=(0,0,0,1,1)\cdot q = -4\,T_{3R}.
$$
These are baryon minus lepton number and right-handed isospin; ordinary hypercharge is $Y=S_3/6-S_2/4$. Our $B_3$ and $B_5$ rows fail the test: every one of them gives different values to the members of the quark doublet. So the smart badges of Section 4 break colour SU(3) itself.

Can the allowed charges do the job? Use $B=\begin{pmatrix}1&1&1&0&0\\0&0&0&1&1\end{pmatrix}$, so $Bq=(S_3,S_2)$. The values on each multiplet are:

| Multiplet | $Q$ | $u^c$ | $d^c$ | $L$ | $e^c$ | $\nu^c$ |
|:--|:-:|:-:|:-:|:-:|:-:|:-:|
| $(S_3,S_2)$ | $(1,0)$ | $(-1,2)$ | $(-1,-2)$ | $(-3,0)$ | $(3,-2)$ | $(3,2)$ |

: **Table 2.** Family $B-L$ and right-handed isospin on the multiplets.

Apply Lemma 2. A twin pair is damaging if $(S_3,S_2)=\lambda(3,2)$ for some $\lambda\ne0$. For $Q$ that needs $1\equiv0$. For $u^c$, $\lambda=1$ and $4\equiv0$. For $d^c$, $\lambda=-1$ and $2\equiv0$. For $L$, $\lambda=0$, excluded. For $e^c$, $\lambda=-1$ and $6\equiv0$. None of these can happen for $p\geq5$. So:

**Proposition 1.** Family-dependent $B-L$ and right-handed isospin, both of which commute with the Standard Model gauge group, give [[30,2,3]]~p~ and [[45,3,3]]~p~ triorthogonal codes for every prime $p\geq5$. At $p=3$ the $e^c$ twins are damaging and the distance drops to two.

The X distance is 7 at $p=5$ and $p=7$. Both rows are needed. $B-L$ alone leaves the lepton-doublet twins damaging, and right-handed isospin alone leaves the right-handed twins damaging; either way the distance is two.

One consequence follows at once. Models that gauge a separate hypercharge for each family [5] cannot separate the logical slots by themselves, because the right-handed neutrino has zero hypercharge. Those lattices give $k=1$. Family-dependent $B-L$ is required.

# 6. How far this can go

## 6.1 A single generation is rigid

Could we add a new charge to a single generation, a row that is not a combination of the five SO(10) charges, and still keep triorthogonality? The answer is no, and the proof is short once you see the right basis.

Since $q_5=q_1q_2q_3q_4$, the sixteen weights are the group $(\mathbb Z_2)^4$. Its sixteen characters are the products $\chi_s(q)=\prod_i q_i^{s_i}$ for $s\in\mathbb F_2^4$, and they are orthogonal: $\sum_q\chi_s(q)\chi_t(q)=16\,\delta_{st}$, which is invertible modulo every odd prime. The five charges are five of the characters, with masks $A=\{1,2,4,8,15\}$.

**Theorem 1 (rigidity).** For every odd prime $p$, any row $y:\{\text{16 weights}\}\to\mathbb F_p$ that can be added to the five SO(10) rows of one generation while keeping triorthogonality is a linear combination of those rows.

*Proof.* Triorthogonality with two old rows requires $\sum_q y(q)\,q_iq_j=0$ for all $i\le j$. The products $q_iq_j$ are characters: for $i=j$ the trivial one, for $i<j$ the character with mask $a_i\oplus a_j$. The ten masks $a_i\oplus a_j$ are exactly the ten masks outside $\{0\}\cup A$. So $y$ is orthogonal to eleven of the sixteen characters, and lies in the span of the remaining five, which are the original rows. $\square$

So the only freedom is in how the families are coupled.

## 6.2 Distance three is a ceiling

Every lattice so far assigns family $g$ a column of the form $M_gq$, a linear function of the SO(10) weight. That covers $B_3$, $B_5$, the gauge-compatible lattice, direct sums, and anything built from family-dependent Cartan charges. There is a three-term identity among the weights that holds over the integers:
$$
q(u^c_3)+q(d_3)-q(e) = q_0 .
$$
Apply $M_g$ to both sides. The error that applies $Z$, $Z$ and $Z^{-1}$ to $u^c_3$, $d_3$ and $e$ in family $g$ has charge $M_gq_0$, the charge of family $g$'s logical slot.

**Theorem 2 (ceiling).** If each family's columns are a linear function of the weight and each logical slot carries nonzero charge, the code distance is at most three.

That closes the obvious route to distance four.

## 6.3 What can and cannot be coupled

Two more results sharpen the picture.

**Coupling can be an illusion.** Try nonlinear differential rows built from characters, $(\chi_s,-\chi_s)$ on the two families. The conditions reduce to $s_a\oplus s_b\notin A$, and protecting all twins needs five such rows. There are exactly sixteen valid sets, $\{a\oplus t:a\in A\}$ for $t\in\mathbb F_2^4$. Every one of them is two independent copies of the one-generation code in disguise: replacing the shared and differential rows by their half-sum and half-difference, then swapping each field with its twin whenever $\chi_t(q)=-1$ turns the matrix into a direct sum. Yet in the original labels, 44 of the 50 smallest damaging errors cross between families. So counting cross-family errors does not prove that the families are coupled. Section 7.2 gives a test that does.

**Nonlinear rows mod 7.** Allow a row that is triorthogonal modulo $p$ but not over the integers. Extending $B_5$ this way, there are no solutions at $p=5$ or $p=11$, and solutions exist at $p=7$ (588 rows in 17 orbits) and $p=13$ (4,056 rows in 57 orbits). Every one of them still has distance three, but the number of three-site damaging errors drops from 50 to 8. At $p=7$ this lowers the leading distillation error from $\tfrac{25}{18}\varepsilon^3$ to $\tfrac{2}{9}\varepsilon^3$, a factor of 6.25, with X distance 9. These rows are not consistent U(1) charges and do not respect the gauge group. They are coding objects, not physics.

# 7. What the coupling does

## 7.1 Error channels, exactly

Feed noisy magic states into the code. In the symmetric model of [1], each physical qudit suffers $Z^j$, $j=1,\dots,p-1$, each with probability $\varepsilon/(p-1)$. Measure the X checks and keep the output only if they all pass. We want the exact probability of every logical outcome $Z^c$.

There is a clean formula. For any target charge $s$,
$$
\sum_{e\,:\,G_Pe=s} x^{\mathrm{wt}(e)} \;=\; p^{-r}\sum_{a\in\mathbb F_p^r}\omega^{-a\cdot s}\,\big(1+(p-1)x\big)^{\,n-w(a)}\,(1-x)^{\,w(a)},\qquad w(a)=\mathrm{wt}(aG_P).
$$
Set $s=G_Rc$, so that $a\cdot s=b\cdot c$ with $b=aG_R$. Sorting all $p^r$ row combinations by the pair $(w(a),b)$ gives integer counts, and the character sum collapses because $a\mapsto ua$ permutes the nonzero classes. The output is the exact number of errors of each weight in each logical class, for every code in Table 1. It reproduces the known $35\varepsilon^3$ for the 15-qubit code and the $6.25\varepsilon^3$ per output of [1] for qutrits, and it agrees with a brute-force enumeration of every error of weight at most four.

| Code | one output wrong | both outputs wrong | per-output infidelity, $\varepsilon=10^{-3}$ |
|:--|:--|:--|:--|
| Coupled $B_3$, $p=3$ | $\tfrac{25}{2}\varepsilon^3+\dots$ | $5\,\varepsilon^4+\dots$ | $6.284\times10^{-9}$ |
| Two copies, $p=3$ | $\tfrac{25}{2}\varepsilon^3+\dots$ | $\tfrac{625}{16}\varepsilon^6+\dots$ | $6.276\times10^{-9}$ |
| Coupled $B_5$, $p=5$ | $\tfrac{25}{8}\varepsilon^3+\dots$ | $\tfrac92\,\varepsilon^4+\dots$ | $1.573\times10^{-9}$ |
| Two copies, $p=5$ | $\tfrac{25}{8}\varepsilon^3+\dots$ | $\tfrac{625}{256}\varepsilon^6+\dots$ | $1.568\times10^{-9}$ |
| Coupled, three families, $p=3$ | $\tfrac{75}{4}\varepsilon^3+\dots$ | $15\,\varepsilon^4+\dots$ | $6.292\times10^{-9}$ |
| Coupled, three families, $p=5$ | $\tfrac{75}{16}\varepsilon^3+\dots$ | $\tfrac{27}{2}\varepsilon^4+\dots$ | $1.579\times10^{-9}$ |

: **Table 3.** Leading terms of the exact logical channels. Rates for "one output wrong" and "both outputs wrong" are summed over outputs.

![**Figure 3.** Logical error after post-selection for two generations at $p=5$. The single-output curves of the coupled and independent codes lie on top of each other. The correlated curves are separated by two powers of $\varepsilon$.](figures/fig3_channels)

Read Figure 3 this way. To leading order the coupled code is no better and no worse than two separate codes: the same acceptance through $\varepsilon^3$, the same single-output error. What is new is correlation. The coupled code produces errors that hit both outputs at order $\varepsilon^4$; separate codes need two independent failures and first do it at $\varepsilon^6$. All of this comes from errors with two sites in each of two families. At $p=5$ there are 1,584 such weight-four damaging patterns, and 1,152 of them hit both outputs. At $p=3$ the counts are 176 and 80. Separate codes have none.

The practical cost is small. At $\varepsilon=10^{-3}$ the per-output infidelity is 0.1% to 0.7% worse than for separate codes, and correlated failures are $10^5$ to $10^6$ times more likely but still below $10^{-11}$. Measuring the checks also costs more: the coupled code has fewer X checks but heavier ones, weights 15 and 16 against 8, for a total check weight of 170 against 144 at $p=5$.

## 7.2 The coupled codes are genuinely coupled

Section 6.3 showed that a code can look coupled and be a direct sum in disguise. Here is a test that cannot be fooled.

**Proposition 2.** A stabilizer code with stabilizer group $S$ factorizes as $\mathcal C_A\otimes\mathcal C_B$ across a division $A|B$ of the qudits if and only if $S=S_A\times S_B$, where $S_A$ is the subgroup supported on $A$.

*Proof.* The code projector is $\Pi=|S|^{-1}\sum_{s\in S}s$. If $\Pi=\Pi_A\otimes\Pi_B$, its Pauli expansion is supported on a product set $T_A\times T_B$. The identity lies in $T_B$, so every $(t_A,1)$ lies in $S$, which gives $S=S_A\times S_B$. The converse is immediate. $\square$

Whether a code factorizes across some division does not change under local unitaries or qudit relabellings. For a CSS code the splittings of $S$ are unions of connected components of two linear matroids, the X-check space and the row space of $G_P$, which we compute exactly. Every coupled code, for two and three families, at every prime tested, and the gauge-compatible lattice too, has a single component. Every direct sum has two.

**Corollary.** No coupled code in Table 1 is equivalent to two or three independent one-generation codes, or to any tensor product of smaller codes, under any local unitary and qudit relabelling.

A second, independent invariant agrees. For any four qudits $A$, the reduced code state has purity $\mathrm{Tr}\,\rho_A^2=p^{\dim S_A-4}$. Count the four-qudit sets where the purity is raised. Two independent copies have 150. The coupled codes have 184 at $p=3$ and 249 at $p=5$; with three families the numbers are 327 and 522 against 225.

We also compared the codes with published qudit families. They are not members of the Prakash–Saha qutrit family [[9m-k,k,2]]~3~ [6], which has distance two and cannot match $(n,k)=(30,2)$ or $(45,3)$, nor of the full-evaluation punctured Reed–Muller codes of [7], which need $n+k=p^m$; here $n+k$ is 32 or 48.

## 7.3 No new gates

Does the coupling buy an entangling logical gate? Consider every transversal diagonal gate $\bigotimes_j D^{w_j}$, with $D=\omega^{x^3}$ or $\omega^{x^2}$ and arbitrary powers $w_j$. The gate preserves the code if and only if
$$
\sum_j w_j\,(x\cdot h_j)(y\cdot h_j)(s\cdot h_j)=0\quad\text{for all } x,y \text{ and all } s \text{ with } sG_R=0,
$$
or the analogous quadratic condition. This is linear in $w$, so the full space of such gates can be computed. The logical gate is a cubic or quadratic form in $c$, and it entangles exactly when the form has mixed terms. For every code here, coupled, gauge-compatible, mod-7 nonlinear, direct sum and family number, the mixed part is zero. The coupling shows up in errors and never in gates.

# 8. Does nature cooperate?

## 8.1 Exact family charges do not survive the masses

The codes need family-dependent charges that are exactly conserved. The masses and mixings of the quarks decide whether that is possible.

**Proposition 3.** Consider $G\geq2$ quark families with one Higgs doublet and an unbroken family charge in the gauge-compatible span, $a_iS_3+b_iS_2$ on family $i$. If the up and down mass matrices have full rank and the quark mixing matrix has no nontrivial block-diagonal form, then $a_i$ and $b_i$ are the same for every family.

*Proof.* The doublet $Q_i$ carries $a_i$, and $u^c_j$ and $d^c_j$ carry $-a_j+2b_j$ and $-a_j-2b_j$. An up-type Yukawa entry $(i,j)$ is allowed only if $a_i-a_j+2b_j+\eta_u=0$, with $\eta_u$ the Higgs charge. Invariance makes $Y_uY_u^\dagger$ and $Y_dY_d^\dagger$ commute with the charge on the doublets, so if the $a_i$ differed both matrices would be block diagonal in the same basis, and the mixing matrix would be too. So all $a_i$ are equal. Then each column $j$ of $Y_u$ is either entirely allowed or entirely forbidden according to $2b_j+\eta_u$. Full rank forbids a zero column, and 2 is invertible, so all $b_j$ are equal. $\square$

We checked the same conclusion by deterministic enumeration of every charge assignment mod $p$ for $G=2$ at $p=3,5,7$ and $G=3$ at $p=3,5$: exactly $p^2$ assignments survive, all family-universal. With family-universal charges the logical columns coincide and $k=1$.

So the Standard Model with exact family charges gives back one logical qudit, no matter how many generations there are.

## 8.2 How badly broken?

If exact symmetry is impossible, the next question is how small the breaking can be. The answer depends on which fields carry the charge.

A family charge on the **doublets** forces the mixing matrix to be block diagonal. The breaking is at least the largest mixing element that no choice of basis can remove. A charge on the **right-handed fields** alone is not touched by mixing; with one Higgs it allows one Yukawa column per sector, and the breaking is set by mass ratios.

| Charge on the doublets separates | Quark doublets (CKM) | Lepton doublets (PMNS) |
|:--|:-:|:-:|
| family 3 from families 1 and 2 | 0.041 | 0.66 |
| family 1 from families 2 and 3 | 0.23 | 0.55 |

: **Table 4.** The largest mixing element that no basis choice removes, for a family charge on the doublets. Inputs: the CKM global fit [8] and NuFIT 6.0 [9].

For charges on the right-handed fields only, the breaking is set by the second-heaviest mass in each sector relative to the heaviest. With running masses at $M_Z$ [10] these are $m_c/m_t=0.004$, $m_s/m_b=0.019$ and $m_\mu/m_\tau=0.059$.

Now connect this to the codes. Right-handed isospin alone is broken only at the percent level, but by Section 5 it gives distance two. Distance three needs family $B-L$ as well, and $B-L$ sits on the lepton doublets with $S_3=-3$. Neutrino mixing breaks that at order one. **No distance-three code is protected by a small symmetry-breaking parameter in the observed world.**

## 8.3 Superselection

There is a deeper obstruction that does not care about flavour. Recall from Section 2.2 that in the conjugate basis the logical digit labels a charge sector. Physical gauge charges obey superselection: nature does not let you prepare coherent superpositions of states with different total gauge charge. A logical qudit that needs such superpositions cannot be stored in real particle charges. The codes in this paper are codes about the charge table, to be run on engineered hardware, not codes that elementary particles run.

## 8.4 A remark on $E_8$

Lisi's proposal [11] places one generation of matter inside the exceptional group $E_8$, and $E_8$ looks tailor-made for our question. Under $E_8\supset SO(10)\times SO(6)$, the 128 spinorial roots split as
$$
128 = (16,4)\oplus(\overline{16},\bar 4),
$$
four copies of a generation labelled by an SU(4) weight, plus four mirror copies. Distler and Garibaldi proved that no such embedding of gravity and the Standard Model is chiral [12]. In our language that theorem becomes a statement about distance.

* The chiral block $(16,4)$ alone is not triorthogonal. The triple sum over the three SU(4) rows is 64, so it fails over the integers and modulo 3, 5 and 7.
* The full 128 is triorthogonal, because every column now has its negative. Then each logical slot's mirror sits among the physical qudits as a single-site damaging error. The code is [[124,4,1]]~p~.
* Using the SU(4) weights as family labels on the chiral block gives the family-number code [[60,4,2]]~p~. The labels of a Lie algebra add, $(x,y)\cdot(q,w)=x\cdot q+y\cdot w$, so they can never produce the products $f_g\,(b\cdot q)$ that distance three needs.

Chirality and anomaly cancellation are the same thing here, and $E_8$ can only cancel anomalies the trivial way, by pairing every field with its mirror. That kills the code.

# 9. What we learned

Here is the logic in one chain. Anomaly cancellation is triorthogonality, so one generation is a 15-qudit magic code. The other generations give the same code, because the charge table does not see mass. Joining them needs charges that tell the families apart, and those charges must depend on the field, not just the family. The forces allow exactly two such charges per family, $B-L$ and right-handed isospin, and with both of them the families couple into new [[30,2,3]]~p~ and [[45,3,3]]~p~ codes for $p\geq5$. A single generation admits no further charges, distance three is a ceiling, and the coupling shows up only as correlated errors. Then nature closes the door: masses and mixings erase every exact family charge, and the part that distance three needs is broken at order one by neutrinos.

The deeper point is about levels of description. At the level of symmetry, where only charges exist, the three generations are perfect twins. At the level of mass, where particles actually travel, they are blended together. A multi-family code needs a sharp and permanent answer to the question "which family is this?" Neither level provides one. The $E_8$ story says the same thing from the other side: a single large symmetry cannot supply chiral, sharply labelled families either.

Several questions remain open. We have compared the coupled codes with direct sums and with two published families, not with every known qudit code. The gate search covers diagonal transversal gates; non-diagonal and permutation-assisted gates are untested. Nobody has explained why the nonlinear rows of Section 6.3 appear only for primes $p\equiv1 \pmod 3$. And a circuit-level simulation would show whether the correlated $\varepsilon^4$ errors matter inside a real distillation factory.

# Acknowledgments {.unnumbered}

The exact computations, verification scripts and drafting were carried out with the assistance of Claude (Anthropic). Every claim is reproduced by the scripts listed in Appendix A.

# Appendix A. Reproducibility {.unnumbered}

All scripts live in the repository [github.com/HiroSakuraba/Standard-model-generation-as-quantum-error-correction](https://github.com/HiroSakuraba/Standard-model-generation-as-quantum-error-correction) and assert the numbers they print. Python 3.9 or later; the last five need NumPy.

| Result | Script | Runtime |
|:--|:--|:--|
| One- and two-generation codes, minimality (§3, §4.1) | `matter_qec_study.py` | seconds |
| Three-generation codes and supports (§4.2) | `three_generation_flavor.py` | ~20 s |
| Exact Yukawa enumeration (§8.1) | `exact_yukawa.py` | seconds |
| Gauge-compatible lattice, channels, equivalence, gates, mod-$p$ rows, breaking table (§5–§8) | `open_steps.py` | ~4 min |
| Rigidity and disguised direct sums (§6.1, §6.3) | `walsh_character_lattice.py` | seconds |
| Four-site purity invariant (§7.2) | `local_purity.py` | seconds |
| Published qutrit family (§7.2) | `published_qutrit_family.py` | seconds |
| Reed–Muller family (§7.2) | `reed_muller_family_exclusion.py` | seconds |
| $E_8$ comparison (§8.4) | `e8_comparison.py` | ~1 s |
| Figures | `paper/make_figures.py` | seconds |

# References {.unnumbered}

1. S. Prakash and K. L. Ramya, "One generation of the Standard Model as a quantum error-correcting code," [arXiv:2610.21879](https://arxiv.org/abs/2610.21879) (2026).
2. S. Bravyi and A. Kitaev, "Universal quantum computation with ideal Clifford gates and noisy ancillas," Phys. Rev. A 71, 022316 (2005), [arXiv:quant-ph/0403025](https://arxiv.org/abs/quant-ph/0403025).
3. S. Bravyi and J. Haah, "Magic-state distillation with low overhead," Phys. Rev. A 86, 052329 (2012), [arXiv:1209.2426](https://arxiv.org/abs/1209.2426).
4. E. T. Campbell, H. Anwar and D. E. Browne, "Magic-state distillation in all prime dimensions using quantum Reed–Muller codes," Phys. Rev. X 2, 041021 (2012), [arXiv:1205.3104](https://arxiv.org/abs/1205.3104).
5. M. Fernández Navarro and S. F. King, "Tri-hypercharge: a separate gauged weak hypercharge for each fermion family as the origin of flavour," [arXiv:2305.07690](https://arxiv.org/abs/2305.07690) (2023).
6. S. Prakash and T. Saha, "Low overhead qutrit magic state distillation," Quantum 9, 1768 (2025), [doi:10.22331/q-2025-06-12-1768](https://doi.org/10.22331/q-2025-06-12-1768).
7. T. Saha and S. Prakash, "Sublogarithmic distillation in all prime dimensions using punctured Reed–Muller codes," Phys. Rev. A 113, 042405 (2026), [doi:10.1103/tzkz-9p2f](https://doi.org/10.1103/tzkz-9p2f).
8. Particle Data Group, "CKM quark-mixing matrix," Review of Particle Physics (2020), [review PDF](https://pdg.web.cern.ch/2020/reviews/rpp2020-rev-ckm-matrix.pdf).
9. I. Esteban et al., "NuFit-6.0: updated global analysis of three-flavor neutrino oscillations," [arXiv:2410.05380](https://arxiv.org/abs/2410.05380) (2024).
10. G.-y. Huang and S. Zhou, "Precise values of running quark and lepton masses in the Standard Model," [arXiv:2009.04851](https://arxiv.org/abs/2009.04851) (2020).
11. A. G. Lisi, "An exceptionally simple theory of everything," [arXiv:0711.0770](https://arxiv.org/abs/0711.0770) (2007).
12. J. Distler and S. Garibaldi, "There is no 'Theory of Everything' inside E8," Commun. Math. Phys. (2010), [arXiv:0905.2658](https://arxiv.org/abs/0905.2658), [doi:10.1007/s00220-010-1006-y](https://doi.org/10.1007/s00220-010-1006-y).
