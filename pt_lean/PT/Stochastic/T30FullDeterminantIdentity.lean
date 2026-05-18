/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Stochastic.T30TraceDeterminant
import PT.Stochastic.T30CharPolyComplete

/-!
# T2 — Full determinant identity for `T_30 = T_2 ⊗ T_3 ⊗ T_5`

**Statement (paper-level, monograph ch03 §T2, structural exhaustive
exposition).** This file exposes — exhaustively and in one place — the
**determinantal closure** of the CRT factorisation
`T_30 ≅ T_2 ⊗ T_3 ⊗ T_5` and *all* its immediate algebraic consequences
for the characteristic polynomial.

The single master identity is

  `det(T_30)  =  det(T_5)²` ,

unconditional in `T_5`. Every other determinant statement in the file is
an algebraic reformulation:

* **Kronecker–product factorisation** : `det(T_30) =
  (det T_2)^{|n|·|p|} · (det T_3)^{|m|·|p|} · det(T_5)^{|m|·|n|}` with
  index cardinalities `|m| = |Fin 1| = 1`, `|n| = |Fin 2| = 2`,
  `|p| = |Fin 4| = 4`.
* **Sign computation** : substituting `det(T_2) = 1` and
  `det(T_3) = -1` gives `det(T_30) = 1^{2·4} · (-1)^{1·4} · det(T_5)^{1·2}
  = (+1) · (+1) · det(T_5)² = det(T_5)²`. The `-1` factor from `T_3` is
  raised to the **even** power `|m|·|p| = 4`, so its sign drops out.
* **Charpoly constant coefficient** : `a_0 = e_8 = det(T_30) = det(T_5)²`.
* **Positivity** : `det(T_30) = det(T_5)² ≥ 0`, *unconditional* — a
  square of a real number.
* **Spectral interpretation** : if the eight eigenvalues of `T_30` are
  `(λ_i)_{i=1..8}`, then `det(T_30) = ∏ λ_i`. Combined with the
  `±`-symmetry of the spectrum proved in `T30CharPolyComplete`
  (eigenvalues come in four sign-pairs), the product is automatically a
  square: `∏ λ_i = ∏_{j=1..4} (λ_j · (-λ_j)) = ∏ (-λ_j²) =
  (-1)^4 · ∏ λ_j² = (∏ λ_j)²`. This is the **spectral shadow** of the
  algebraic identity `det(T_30) = det(T_5)²`.

## Main results

### Reexposition of the master identity
* `master_det_T30` — **headline**: `det(T_30) = det(T_5)²` (alias of
  `T30_det_eq`, exposed under the canonical "master" name for downstream
  consumption).

### Explicit Kronecker form with cardinalities
* `T30_det_kron_expanded` — fully expanded Kronecker determinant:
  `det(T_30) = (det T_2 ^ 2 · det T_3 ^ 1) ^ 4 · det T_5 ^ 2`,
  *before* simplification of the `T_2`/`T_3` values.
* `T30_det_T2_T3_powers` — substitution `det T_2 = 1`, `det T_3 = -1`:
  `(1^2 · (-1)^1)^4 · det T_5 ^ 2 = det T_5 ^ 2`, with the parity
  cancellation `(-1)^4 = 1` made explicit.

### Sign-pinning lemmas
* `T3_det_pow_inner_even` — `(det T_3)^{|m|·|p|} = (-1)^4 = 1`. The key
  sign-cancellation: the antidiagonal involution sign is raised to an
  even power because the outer dimensions multiply to `1·4 = 4`.
* `T30_det_sign_positive_factor` — the global sign factor in front of
  `det(T_5)²` is `+1` (no further sign correction needed).

### Charpoly constant term
* `T30_charpoly_constant_eq_T5_det_sq` — `a_0 = e_8 = det(T_5)²`. This
  re-exposes `T30_charpoly_a0_closed` under the "constant term"
  perspective. The sign convention is the standard one
  `a_0 = (-1)^8 · det(T_30) = det(T_30)` (the `(-1)^d` factor is `+1`
  for `d = 8` even).

### Positivity
* `T30_det_nonneg` — `det(T_30) ≥ 0` (square of a real).
* `T30_det_eq_T5_det_sq_abs` — `det(T_30) = |det(T_5)|²`
  (sign-insensitive form of the headline; `|x|² = x²` for real `x`).
* `T30_det_zero_iff_T5_det_zero` — `det(T_30) = 0 ↔ det(T_5) = 0`,
  i.e. `T_30` is singular iff `T_5` is singular. The CRT factorisation
  **preserves invertibility** (since `T_2`, `T_3` are both invertible:
  `det T_2 = 1 ≠ 0`, `det T_3 = -1 ≠ 0`).

### Headline summary
* `T30_full_determinant_identity` — packages the master identity, the
  expanded Kronecker form, the constant-term equality, positivity, and
  the invertibility equivalence in one statement.

## Strategy

Every result is a one- to three-line consequence of the lemmas
`T30_det_eq` (from `T30TraceDeterminant`) and `T30_charpoly_a0_closed`
(from `T30CharPolyComplete`). The "expanded Kronecker form" is obtained
by inlining the proof of `kron3_det` applied to `T_2, T_3, T_5` and
keeping the cardinality powers explicit before they collapse to `1`.

The sign-pinning step uses only `Even.neg_one_pow` (or directly
`(-1 : ℝ)^4 = 1` by `norm_num`). The positivity step uses
`sq_nonneg` (or `pow_two_nonneg`).

## Status

`[THM]` — purely algebraic consequence of the Kronecker determinant
formula and the trivial invariants `det T_2 = 1`, `det T_3 = -1`. No
spectral data of `T_5` enters. The identity `det(T_30) = det(T_5)²` is
*unconditional* in `T_5`, and so are all its consequences exposed here.

## References

* `PT.Stochastic.T30TraceDeterminant` — `T30_det_eq`, `T2_trivial_det`,
  `T3_det_neg_one`, `kron3_det`.
* `PT.Stochastic.T30CharPolyComplete` — `T30_charpoly_a0_closed`.
* `PT.Sieve.T3Antidiagonal` — `T3_det_neg_one`.
* Mathlib `Matrix.det_kronecker`.
* Monograph Chapter 3, §"Spectral conservation T2", determinant identities
  table.
-/

namespace PT.Stochastic

open Matrix Kronecker BigOperators PT.Sieve

/-! ### Master identity (reexposition)

The single algebraic fact that drives every other statement in this file:

  `det(T_30) = det(T_5)²` ,

an immediate consequence of the three-factor Kronecker determinant formula
specialised to `det(T_2) = 1`, `det(T_3) = -1`, and the index cardinalities
`(1, 2, 4)`. -/

/-- **Master identity** : `det(T_30) = det(T_5)²` (unconditional).

    This is the canonical "master" alias of `T30_det_eq`, exposed here
    under a stable name so downstream files (and the
    `T30FullSpectralAnalysis` aggregator) can refer to the determinant
    closure without committing to the trace-determinant module name. -/
theorem master_det_T30 (T5 : T5Like) :
    (T30 T5).det = T5.matrix.det ^ 2 :=
  T30_det_eq T5

/-! ### Explicit Kronecker form (cardinalities exposed)

The headline `det(T_30) = det(T_5)²` is obtained from the three-factor
Kronecker determinant formula

  `det((A ⊗ B) ⊗ C) = (det A ^ |n| · det B ^ |m|) ^ |p| · det C ^ |m·n|`

with `|m| = 1`, `|n| = 2`, `|p| = 4`. Before simplifying the `T_2`/`T_3`
values, the formula reads

  `det(T_30) = (det T_2 ^ 2 · det T_3 ^ 1) ^ 4 · det T_5 ^ 2`,

which exposes the role of each cardinality. -/

/-- **Cardinality computations.** Hidden behind `simp` / `decide`. -/
@[simp] lemma card_Fin_one : Fintype.card (Fin 1) = 1 := by decide
@[simp] lemma card_Fin_two : Fintype.card (Fin 2) = 2 := by decide
@[simp] lemma card_Fin_four : Fintype.card (Fin 4) = 4 := by decide

/-- **Cardinality of the inner index `(Fin 1 × Fin 2)`.** Used to keep
    the `kron3_det` signature transparent. -/
@[simp] lemma card_Fin_one_times_two :
    Fintype.card (Fin 1 × Fin 2) = 2 := by decide

/-- **Fully expanded Kronecker determinant for `T_30`** (before
    substituting `det T_2 = 1` and `det T_3 = -1`).

    `det(T_30) = (det T_2 ^ 2 · det T_3 ^ 1) ^ 4 · det T_5 ^ 2`,

    with cardinality exponents `|Fin 2| = 2`, `|Fin 1| = 1`,
    `|Fin 4| = 4`, `|Fin 1 × Fin 2| = 2` made explicit. -/
theorem T30_det_kron_expanded (T5 : T5Like) :
    (T30 T5).det
      = (T2_trivial.det ^ 2 * T3.det ^ 1) ^ 4 * T5.matrix.det ^ 2 := by
  unfold T30
  rw [kron3_det T2_trivial T3 T5.matrix]
  simp

/-! ### Substitution of the trivial-factor invariants

Plug `det(T_2) = 1` and `det(T_3) = -1` into the expanded form. The key
sign-cancellation is `(-1)^4 = 1`: the antidiagonal involution sign is
raised to the **outer cardinality** `|m| · |p| = 1 · 4 = 4`, which is
**even**, so its `-1` drops out. -/

/-- **Sign-pinning at the `T_3` factor.** The exponent `|m| · |p| = 4` is
    even, so `(det T_3)^4 = (-1)^4 = 1`. This is the algebraic origin of
    the *positive* sign in front of `det(T_5)²` in the master identity. -/
theorem T3_det_pow_inner_even :
    (T3.det : ℝ) ^ (Fintype.card (Fin 1) * Fintype.card (Fin 4)) = 1 := by
  rw [T3_det_neg_one]
  simp; norm_num

/-- **Cleaner statement** : `(det T_3)^4 = 1`. -/
theorem T3_det_pow_four_eq_one : (T3.det : ℝ) ^ 4 = 1 := by
  rw [T3_det_neg_one]; norm_num

/-- **Substitution step.** With `det T_2 = 1` and `det T_3 = -1`, the
    expanded Kronecker determinant collapses to `(1^2 · (-1)^1)^4 ·
    det T_5 ^ 2 = (-1)^4 · det T_5 ^ 2 = det T_5 ^ 2`. -/
theorem T30_det_T2_T3_powers (T5 : T5Like) :
    (T30 T5).det
      = ((1 : ℝ) ^ 2 * ((-1 : ℝ)) ^ 1) ^ 4 * T5.matrix.det ^ 2 := by
  rw [T30_det_kron_expanded T5, T2_trivial_det, T3_det_neg_one]

/-- **Positive-sign confirmation.** The global sign factor
    `(1^2 · (-1)^1)^4 = ((-1))^4 = +1` in front of `det(T_5)^2`. -/
theorem T30_det_sign_positive_factor :
    (((1 : ℝ) ^ 2 * ((-1 : ℝ)) ^ 1)) ^ 4 = 1 := by
  norm_num

/-! ### Charpoly constant term

For a square matrix `M` of size `d × d`, the constant coefficient `a_0`
of the characteristic polynomial `det(x · I − M) = x^d − e_1 x^{d-1}
+ … + (-1)^d e_d` equals `(-1)^d · det(M)`. For `T_30` (size `d = 8`),
the parity is **even**, so `a_0 = (+1) · det(T_30) = det(T_5)²`. -/

/-- **Constant-term parity.** For the dimension `d = 8` of `T_30`,
    `(-1)^8 = +1`, so the standard sign convention `a_0 = (-1)^d · det(M)`
    reduces to `a_0 = det(M)`. -/
theorem dim_T30_parity_even : ((-1 : ℝ)) ^ 8 = 1 := by
  norm_num

/-- **Charpoly constant term equals `det(T_5)²`** (canonical form).

    The constant coefficient `a_0` of `charpoly(T_30)` equals
    `(-1)^8 · det(T_30) = det(T_30) = det(T_5)²`. Reexposed from
    `T30_charpoly_a0_closed`. -/
theorem T30_charpoly_constant_eq_T5_det_sq (T5 : T5Like) :
    ((-1 : ℝ)) ^ 8 * (T30 T5).det = T5.matrix.det ^ 2 := by
  rw [dim_T30_parity_even, one_mul]
  exact T30_charpoly_a0_closed T5

/-! ### Positivity and singularity

The square form `det(T_30) = det(T_5)²` makes positivity automatic and
turns the singularity question for `T_30` into the singularity question
for `T_5`. -/

/-- **Positivity** : `det(T_30) ≥ 0` (square of a real number). -/
theorem T30_det_nonneg (T5 : T5Like) : 0 ≤ (T30 T5).det := by
  rw [master_det_T30]
  exact sq_nonneg _

/-- **Absolute-value reformulation** : `det(T_30) = |det(T_5)|²`. -/
theorem T30_det_eq_T5_det_sq_abs (T5 : T5Like) :
    (T30 T5).det = |T5.matrix.det| ^ 2 := by
  rw [master_det_T30, sq_abs]

/-- **Singularity equivalence** : `T_30` is singular iff `T_5` is
    singular. The CRT factorisation preserves invertibility (since
    `det T_2 = 1 ≠ 0` and `det T_3 = -1 ≠ 0`, the determinant of `T_30`
    is controlled entirely by `det T_5`). -/
theorem T30_det_zero_iff_T5_det_zero (T5 : T5Like) :
    (T30 T5).det = 0 ↔ T5.matrix.det = 0 := by
  rw [master_det_T30]
  exact pow_eq_zero_iff (n := 2) (by norm_num)

/-- **Non-vanishing equivalence** : `T_30` is invertible iff `T_5` is
    invertible. The contrapositive of `T30_det_zero_iff_T5_det_zero`. -/
theorem T30_det_ne_zero_iff_T5_det_ne_zero (T5 : T5Like) :
    (T30 T5).det ≠ 0 ↔ T5.matrix.det ≠ 0 := by
  rw [ne_eq, ne_eq, T30_det_zero_iff_T5_det_zero]

/-! ### Spectral shadow

If the eight eigenvalues of `T_30` are `(λ_i)_{i=1..8}`, then
`det(T_30) = ∏ λ_i`. The `±`-symmetry of the spectrum (cf.
`T30_charpoly_is_even`) implies eigenvalues come in four sign-pairs
`(μ_j, -μ_j)`, so

  `∏_i λ_i = ∏_{j=1..4} μ_j · (-μ_j) = (-1)^4 · ∏_j μ_j^2 = (∏_j μ_j)²` ,

automatically a perfect square. This is the **spectral shadow** of the
algebraic identity `det(T_30) = det(T_5)²`: the half-spectrum `(μ_j)`
parametrises a *square root* of `det(T_30)`, which the algebraic
identity identifies (up to sign) with `det(T_5)`.

We record only the algebraic prefactor `(-1)^4 = 1` here; committing to
a specific eigenvalue datum is left to `T30FullEigenpairCount`. -/

/-- **Sign-pair prefactor** : `(-1)^4 = 1`. The product of four sign-pairs
    `(μ_j, -μ_j)` carries a global sign `(-1)^4 = +1`, so the eigenvalue
    product is automatically a non-negative square. -/
theorem sign_pair_prefactor_four : ((-1 : ℝ)) ^ 4 = 1 := by norm_num

/-- **Eigenvalue-product square form (algebraic).** For any four reals
    `μ_1, μ_2, μ_3, μ_4`, the product over four `±`-pairs
    `(μ_j, -μ_j)` is a square:

    `μ_1·(-μ_1)·μ_2·(-μ_2)·μ_3·(-μ_3)·μ_4·(-μ_4) = (μ_1·μ_2·μ_3·μ_4)²` ,

    a purely algebraic identity that the `±`-symmetric `T_30` spectrum
    matches against the master identity `det(T_30) = det(T_5)²`. -/
theorem four_sign_pair_product_is_square (μ₁ μ₂ μ₃ μ₄ : ℝ) :
    μ₁ * (-μ₁) * (μ₂ * (-μ₂)) * (μ₃ * (-μ₃)) * (μ₄ * (-μ₄))
      = (μ₁ * μ₂ * μ₃ * μ₄) ^ 2 := by
  ring

/-! ### Headline summary -/

/-- **Headline summary for the full determinantal identity of `T_30`.**

    The CRT factorisation `T_30 = T_2 ⊗ T_3 ⊗ T_5` (left-bracketed
    Kronecker product, dimensions `1 · 2 · 4 = 8`) yields the closed
    determinantal identity

      `det(T_30) = det(T_5)²` ,

    which we package together with:

    1. **Master identity** — `det(T_30) = det(T_5)²` (unconditional).
    2. **Expanded Kronecker form** — `det(T_30) = (det T_2 ^ 2 ·
       det T_3 ^ 1) ^ 4 · det T_5 ^ 2`, with the cardinality exponents
       `(2, 1, 4)` exposed before simplification.
    3. **Sign cancellation** — `(det T_3)^4 = 1` : the antidiagonal sign
       is killed by the even outer exponent `|m|·|p| = 4`.
    4. **Charpoly constant term** — `a_0 = (-1)^8 · det(T_30) =
       det(T_5)²`.
    5. **Positivity** — `det(T_30) ≥ 0` (square of a real).
    6. **Singularity equivalence** — `det(T_30) = 0 ↔ det(T_5) = 0`
       (CRT preserves invertibility through `T_2`, `T_3`).
    7. **Sign-pair prefactor** — `(-1)^4 = 1`, the algebraic shadow of
       the eigenvalue `±`-symmetry consistent with the squareness of
       `det(T_30)`. -/
theorem T30_full_determinant_identity (T5 : T5Like) :
    -- 1. Master identity
    (T30 T5).det = T5.matrix.det ^ 2
    -- 2. Expanded Kronecker form
    ∧ (T30 T5).det
        = (T2_trivial.det ^ 2 * T3.det ^ 1) ^ 4 * T5.matrix.det ^ 2
    -- 3. Sign cancellation at T_3
    ∧ (T3.det : ℝ) ^ 4 = 1
    -- 4. Charpoly constant term
    ∧ ((-1 : ℝ)) ^ 8 * (T30 T5).det = T5.matrix.det ^ 2
    -- 5. Positivity
    ∧ 0 ≤ (T30 T5).det
    -- 6. Singularity equivalence
    ∧ ((T30 T5).det = 0 ↔ T5.matrix.det = 0)
    -- 7. Sign-pair prefactor
    ∧ ((-1 : ℝ)) ^ 4 = 1 :=
  ⟨ master_det_T30 T5,
    T30_det_kron_expanded T5,
    T3_det_pow_four_eq_one,
    T30_charpoly_constant_eq_T5_det_sq T5,
    T30_det_nonneg T5,
    T30_det_zero_iff_T5_det_zero T5,
    sign_pair_prefactor_four ⟩

end PT.Stochastic
