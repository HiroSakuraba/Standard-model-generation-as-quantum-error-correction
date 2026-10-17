# A local-unitary invariant of the family-coupled codes

**Status:** reproducible calculation and proof for the specific charge-lattice CSS constructions at p=3 and p=5. Not a complete classification against known qudit codes. This note continues `EXACT_SYMMETRY_AND_CLASSIFICATION.md` and is verified by `local_purity.py`.

## Result

A stronger, operationally interpretable way to distinguish the coupled codes from copies of the 15-qudit code is to calculate the purity of every four-qudit reduction of the **normalized code-space projector**.

For the p-level code with stabilizer group S and code projector P, put `rho = P / p^k`. For any set A of four physical qudits let `rho_A = Tr_(complement A) rho`. Then

```
Tr(rho_A^2) = p^(dim S_A - 4),
S_A = {s in S : support(s) subset of A}.
```

This is a **local-unitary (LU) invariant**: arbitrary local changes of basis, and a permutation of the qudits, do not change the multiset of four-site purities.

For the explicit 2- and 3-family examples, every four-site subset has `dim S_A` either 0 or 1. Call a subset a *purity hotspot* when `dim S_A=1`; its purity is p^-3 rather than the baseline p^-4.

| Prime | Families | Independent 15-qudit copies | Coupled construction | Extra 2+2 hotspots per family pair |
|---:|---:|---:|---:|---:|
| 3 | 2 | 150 | **184** | 34 |
| 3 | 3 | 225 | **327** | 34 |
| 5 | 2 | 150 | **249** | 99 |
| 5 | 3 | 225 | **522** | 99 |

The count of 2+2 hotspots is a **Z-stabilizer count**, NOT a Z-logical-error count. It must not be confused with the 88 (p=3) and 396 (p=5) mixed 2+2 logical supports found in the earlier study. These are distinct types of four-site relations: null charge vs nonzero puncture charge.

### Proof of the purity formula

The normalized CSS stabilizer-code projector is

```
rho = p^(-n) sum_(s in S) s.
```

Taking the partial trace outside A kills every Pauli term that acts nontrivially there and produces

```
rho_A = p^(-|A|) sum_(s in S_A) s|_A .
```

Orthogonality of distinct qudit Pauli operators under the Hilbert–Schmidt inner product gives

```
Tr(rho_A^2) = |S_A| / p^|A| = p^(dim S_A - |A|).
```

The purity is unchanged by conjugation with `U_A = tensor_(i in A) U_i`. Any permutation merely relabels A.

### Efficient exact computation

For the family's Cartan-linear constructions, every X-stabilizer restricted to one family is the evaluation of some linear functional `a·q` on the 15 non-singlet SO(10) spinor weights `q`. Exhaustive evaluation of all nonzero `a in F_p^5` gives minimum nonzero family support **7**, for both p=3 and p=5. Therefore a nonidentity X stabilizer cannot be supported on at most four sites.

For any four-site subset A the **entire** locally supported stabilizer subgroup is Z-type. The Z stabilizers are the kernel of the punctured charge matrix G_P, so

```
dim S_A = dim ker (G_P restricted to columns A)
        = 4 - rank_Fp(G_P|_A).
```

The script enumerates all 4-subsets, computes these tiny ranks exactly modulo p, checks `dim S_A <=1`, and groups the hotspots by distribution across families. The total four-subset counts are `binomial(30,4) = 27,405` and `binomial(45,4) = 148,995` per model.

```bash
python -m pip install numpy
python local_purity.py
# ALL LOCAL-PURITY CHECKS PASSED
```

### Consequence and non-consequence

Different four-site purity multisets give a direct proof that these particular coupled codes are inequivalent to their independent-copy baselines under **arbitrary single-qudit local unitaries plus permutations**, independent of the previously used stabilizer-factorization proof.

They do **not** prove novelty relative to all previously constructed indecomposable p-ary triorthogonal codes. Four-site purity is only one invariant; its equality would not establish equivalence.

No claim is made about gauge superselection sectors, physical recovery operations, actual quark/lepton noise, or an experimentally realized encoding.

## External literature comparison: scope audit

1. [Shiroman Prakash, *Cyclic triorthogonal codes in prime dimension*, arXiv:2610.08012](https://arxiv.org/abs/2610.08012), formulates its spectral-support classification for cyclic lengths **coprime** to p. Here n=30 or 45 and the checked p=3 or 5 always divide n. Thus that particular length-coprime spectral classification does **not** directly cover the present parameters; this says nothing about equivalence to some other cyclic construction.
2. [Tanay Saha and Shiroman Prakash, *Sublogarithmic distillation in all prime dimensions using punctured Reed-Muller codes*, Physical Review A 113, 042405 (2026)](https://doi.org/10.1103/tzkz-9p2f) is genuinely relevant to odd-prime qudit comparisons. No explicit candidate matrix from that work was compared here; no equivalence or inequivalence claim is made.
3. [Adam Wills, Shubham P. Jain, and Shraddha Singh, *Classification of Generalised Triorthogonal Codes through Length 54*, arXiv:2609.30860](https://arxiv.org/abs/2609.30860) covers the qubit T-state/generalized-binary setting. Despite including lengths 30 and 45, it is **not** a catalogue of p=3 or p=5 codes.

**Next comparison milestone:** obtain a specifically documented, explicit odd-prime competitor generator matrix with matching (n,k,p), normalize its X and Z stabilizer conventions, and run these LU-invariant diagnostics plus an equivalence algorithm. Until then, publication language should say *coupled, indecomposable constructions; equivalence to the wider odd-prime code literature untested*.
