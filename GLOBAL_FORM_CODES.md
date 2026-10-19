# The global form of the gauge group is a choice of stabilizer code

**Status, 19 October 2026:** exploratory note, not peer reviewed. Verified by `python3 global_form_codes.py` (numpy; under a second; ends with `ALL GLOBAL-FORM CODE CHECKS PASSED`).

**Question.** The Standard Model gauge group is (SU(3) × SU(2) × U(1)_Y)/Γ with Γ = 1, ℤ2, ℤ3 or ℤ6. No experiment so far distinguishes them. The difference is which line operators exist (Tong, [arXiv:1705.01853](https://arxiv.org/abs/1705.01853)). Can the choice be written as a code?

**Answer.** Yes. The possible line spectra are the twelve stabilizer bases of one six-level qudit. The four values of Γ are the orbits of these bases under θ-angle shifts, and those shifts act as Clifford phase gates.

## From lines to classes

Tong labels a line by electric charges (z_e2, z_e3, q) and magnetic charges (z_m2, z_m3, g), with the Dirac condition (his eq. 3.1)

  3 z_e2 z_m2 + 2 z_e3 z_m3 − 6qg ∈ 6ℤ.

Standard Model matter screens Wilson lines, so only classes matter. The script checks that:

- every SM field, and the Higgs, has electric class a = (z_e2 − q mod 2, z_e3 − q mod 3) = (0, 0);
- the magnetic labels local with all matter are exactly those with 6g ≡ 3 z_m2 + 2 z_m3 (mod 6); write b = (z_m2, z_m3);
- on classes, mutual locality of two dyonic lines x = (a, b) and x′ = (a′, b′) becomes

  ⟨x, x′⟩ = 3(a₂b₂′ − b₂a₂′) + 2(a₃b₃′ − b₃a₃′) ≡ 0 (mod 6).

So lines live in ℤ6 × ℤ6 with a symplectic form.

## Twelve line spectra

A consistent theory has a maximal mutually local set of lines: a Lagrangian subgroup L, isotropic and of order 6. There are exactly 12. Sort them by the image of L in the magnetic labels, which is the quotient Γ:

| Γ | Line spectra | Allowed Wilson-line classes |
|---|---|---|
| 1 | 1 | 6 |
| ℤ2 | 2 | 3 |
| ℤ3 | 3 | 2 |
| ℤ6 | 6 | 1 |

The θ-angle shifts of SU(2) and SU(3) act through the Witten effect a → a + b. They permute the 12 spectra, and the orbits have sizes 1, 2, 3, 6, one for each Γ. So the choice of Γ is a discrete choice, while the spectrum within a fixed Γ moves with θ, as Tong describes. This is the same counting as Aharony, Seiberg and Tachikawa's SU(6) case ([arXiv:1305.0318](https://arxiv.org/abs/1305.0318)), because the ℤ6 data are identical.

## The qudit dictionary

Use ℤ6 = ℤ2 × ℤ3 and send a line to a Weyl operator on a qubit ⊗ qutrit:

  (a₂, b₂, a₃, b₃) ↦ Z^{a₂} X^{b₂} ⊗ Z₃^{a₃} X₃^{b₃}.

| Gauge theory | Qudit |
|---|---|
| line operator | Weyl operator |
| mutual locality | commutation (all 36 × 36 pairs checked) |
| global form, i.e. line spectrum L | stabilizer group; its six joint eigenstates form a basis |
| θ₂ → θ₂ + 2π | qubit phase gate S |
| θ₃ → θ₃ + 2π | qutrit phase gate diag(1, 1, ω) |
| electric ↔ magnetic exchange | Fourier gate H ⊗ F₃; maps Γ = 1 to Γ = ℤ6 at θ = 0 |

The script checks every row: each of the 12 spectra gives six one-dimensional joint eigenspaces, and each gate maps the stabilizer basis of one spectrum onto the stabilizer basis of its image. Electric–magnetic exchange is not a symmetry of the Standard Model. It appears here only as a map between codes.

## On a spatial 3-torus

On T³ each global form has flux sectors labelled by L³: 6³ = 216 sectors in every case. Different forms pick different sectors. The number shared with SU(3) × SU(2) × U(1) is |L ∩ L₁|³:

| Γ | Sectors | Shared with Γ = 1 |
|---|---|---|
| 1 | 216 | 216 |
| ℤ2 | 216 | 27 |
| ℤ3 | 216 | 8 |
| ℤ6 | 216 | 1 |

In code language, there are three six-level qudits, one per pair of cycles. The global form is the basis they are read out in.

## Link to the centre echo

In the theory with Γ = ℤ_q, the allowed link classes form a group of order 6/q. The ℤ6 echo of [`SM_CENTER_ECHO.md`](SM_CENTER_ECHO.md) takes 6/q distinct values on them, which is 6, 3, 2 or 1. For Γ = ℤ6 it is the identity. Whether the echo is a working control is therefore a test of the global form, the same question as whether a fractionally charged colour singlet exists.

## Scope

- **Abstraction.** The finite ℤ6 × ℤ6 classes are exact consequences of Tong's eq. 3.1 and the SM charges. The continuous U(1) charges and θ̃ are not modelled. Tong shows that some shifts of θ₂ and θ₃ are equivalent to shifts of θ̃, which does not change the orbit count.
- **Physics.** The qudit is a bookkeeping device, the finite part of the symmetry-TFT picture. No new physical prediction is claimed.
