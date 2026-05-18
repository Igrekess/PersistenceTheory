/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Holonomy.CouplingReconstruction
import PT.Holonomy.CouplingReconstructionBounds
import PT.Holonomy.SinSqProductBounds
import Mathlib.Tactic

/-!
# `∏ sin²θ_p` — Extended chain over `{3, 5, 7, 11, 13}` (Ch09 echo-prime extension)

This file extends `SinSqProductBounds` to the first **two echo primes**
`p = 11` and `p = 13` (the "echo block" `{p ≥ 11}` per the convention
recorded in `feedback_pt_casimir_nomenclature`). It exhibits the strictly
decreasing chain of partial products

  `P_5 < P_4 < P_3 < P_2 < P_1 < 1`

and gives rational brackets for the two new factors `sin²θ_11`, `sin²θ_13`
together with the two new partial products `P_4`, `P_5`.

* `P_4 := P_3 · sin²θ_11` ∈ `(1018/10⁶, 1021/10⁶) ≈ 1.02 × 10⁻³`.
* `P_5 := P_4 · sin²θ_13` ∈ `(127/10⁶, 129/10⁶) ≈ 1.28 × 10⁻⁴`.

Each new sin² factor is verified to lie strictly between `0` and `1`, hence
multiplying the previous partial product by such a factor strictly shrinks it.

## Reference

Monograph Chapter 9 §"Produit partiel des phases cycliques" — echo-prime
extension. `feedback_pt_casimir_nomenclature` for the term "echo primes".
-/

namespace PT.Holonomy

/-! ### Numerical bounds for `sin²θ_11` and `sin²θ_13`

The new sin²-factors are pure rationals (Pythagorean identity applied to
`δ_p = (1 - q_+^p)/p` at `q_+ = 13/15`). We bracket each to three decimals,
matching the float values `0.13895…` and `0.12569…`. -/

/-- `sin²θ_11` lies in `(0.138, 0.140)` (decimal value `≈ 0.13895`). -/
theorem sinSq_11_bracket : 138 / 1000 < sinSqQ 11 ∧ sinSqQ 11 < 140 / 1000 := by
  unfold sinSqQ deltaQ qPT
  refine ⟨?_, ?_⟩ <;> norm_num

/-- `sin²θ_13` lies in `(0.125, 0.127)` (decimal value `≈ 0.12569`). -/
theorem sinSq_13_bracket : 125 / 1000 < sinSqQ 13 ∧ sinSqQ 13 < 127 / 1000 := by
  unfold sinSqQ deltaQ qPT
  refine ⟨?_, ?_⟩ <;> norm_num

/-! ### Positivity and sub-unitarity for the new factors -/

/-- `sin²θ_11 > 0`. -/
theorem sinSq_11_pos : 0 < sinSqQ 11 :=
  lt_trans (by norm_num : (0 : ℚ) < 138 / 1000) sinSq_11_bracket.1

/-- `sin²θ_13 > 0`. -/
theorem sinSq_13_pos : 0 < sinSqQ 13 :=
  lt_trans (by norm_num : (0 : ℚ) < 125 / 1000) sinSq_13_bracket.1

/-- `sin²θ_11 < 1`. -/
theorem sinSq_11_lt_one : sinSqQ 11 < 1 :=
  lt_trans sinSq_11_bracket.2 (by norm_num : (140 / 1000 : ℚ) < 1)

/-- `sin²θ_13 < 1`. -/
theorem sinSq_13_lt_one : sinSqQ 13 < 1 :=
  lt_trans sinSq_13_bracket.2 (by norm_num : (127 / 1000 : ℚ) < 1)

/-! ### Partial-product definitions -/

/-- The four-factor partial product `P_4 = sin²θ_3 · sin²θ_5 · sin²θ_7 · sin²θ_11`. -/
def P4 : ℚ := P3 * sinSqQ 11

/-- The five-factor partial product
    `P_5 = sin²θ_3 · sin²θ_5 · sin²θ_7 · sin²θ_11 · sin²θ_13`. -/
def P5 : ℚ := P4 * sinSqQ 13

/-- `P_4` expanded over the four factors. -/
theorem P4_eq_prod :
    P4 = sinSqQ 3 * sinSqQ 5 * sinSqQ 7 * sinSqQ 11 := by
  unfold P4 P3; rfl

/-- `P_5` expanded over the five factors. -/
theorem P5_eq_prod :
    P5 = sinSqQ 3 * sinSqQ 5 * sinSqQ 7 * sinSqQ 11 * sinSqQ 13 := by
  unfold P5 P4 P3; rfl

/-! ### Positivity -/

theorem P4_pos : 0 < P4 := mul_pos P3_pos sinSq_11_pos
theorem P5_pos : 0 < P5 := mul_pos P4_pos sinSq_13_pos

/-! ### Brackets for `P_4` and `P_5` -/

/-- `P_4 ∈ (1018/10⁶, 1021/10⁶) ≈ 1.02 × 10⁻³`. -/
theorem P4_bracket : 1018 / 1000000 < P4 ∧ P4 < 1021 / 1000000 := by
  unfold P4 P3 sinSqQ deltaQ qPT
  refine ⟨?_, ?_⟩ <;> norm_num

/-- `P_5 ∈ (127/10⁶, 129/10⁶) ≈ 1.28 × 10⁻⁴`. -/
theorem P5_bracket : 127 / 1000000 < P5 ∧ P5 < 129 / 1000000 := by
  unfold P5 P4 P3 sinSqQ deltaQ qPT
  refine ⟨?_, ?_⟩ <;> norm_num

/-! ### Strict decreasing chain `P_5 < P_4 < P_3 < P_2 < P_1 < 1` -/

/-- `P_4 < P_3` (multiplying by `sin²θ_11 < 1` strictly shrinks). -/
theorem P4_lt_P3 : P4 < P3 := by
  unfold P4
  have h11_lt : sinSqQ 11 < 1 := sinSq_11_lt_one
  have h3pos : 0 < P3 := P3_pos
  nlinarith

/-- `P_5 < P_4` (multiplying by `sin²θ_13 < 1` strictly shrinks). -/
theorem P5_lt_P4 : P5 < P4 := by
  unfold P5
  have h13_lt : sinSqQ 13 < 1 := sinSq_13_lt_one
  have h4pos : 0 < P4 := P4_pos
  nlinarith

/-- **Extended strict chain.** `0 < P_5 < P_4 < P_3 < P_2 < P_1 < 1`. -/
theorem partial_product_chain_extended :
    0 < P5 ∧ P5 < P4 ∧ P4 < P3 ∧ P3 < P2 ∧ P2 < P1 ∧ P1 < 1 :=
  ⟨P5_pos, P5_lt_P4, P4_lt_P3, P3_lt_P2, P2_lt_P1, P1_lt_one⟩

/-! ### Ratio observation: the echo block strongly suppresses the product

The two echo factors `sin²θ_11 · sin²θ_13` introduce a multiplicative
suppression on the order of `1.7 × 10⁻²` relative to `P_3`. Concretely
`P_5 / P_3 < 1/50`, i.e. adding `{11, 13}` shrinks the product by more
than a factor 50. -/

/-- The echo factors shrink the product significantly:
    `P_5 < P_3 / 50`. -/
theorem P5_lt_P3_div_fifty : P5 < P3 / 50 := by
  have h5 : P5 < 129 / 1000000 := P5_bracket.2
  have h3 : 7335 / 1000000 < P3 := P3_bracket.1
  have h3div : 7335 / 1000000 / 50 < P3 / 50 := by linarith
  have : (129 : ℚ) / 1000000 < 7335 / 1000000 / 50 := by norm_num
  linarith

/-! ### Numerical headline -/

/-- **Headline (extended partial products, active cascade + echo primes 11, 13).**
    All quantities at the PT fixed point `q = 13/15`:

    * `P_1 = sin²θ_3 ∈ (0.219, 0.220)`,
    * `P_2 = P_1 · sin²θ_5 ∈ (0.042, 0.043)`,
    * `P_3 = P_2 · sin²θ_7 = α_bare ∈ (7335/10⁶, 7340/10⁶)`,
    * `P_4 = P_3 · sin²θ_11 ∈ (1018/10⁶, 1021/10⁶) ≈ 1.02 × 10⁻³`,
    * `P_5 = P_4 · sin²θ_13 ∈ (127/10⁶, 129/10⁶) ≈ 1.28 × 10⁻⁴`,
    * strictly decreasing: `P_5 < P_4 < P_3 < P_2 < P_1 < 1`,
    * echo-block suppression: `P_5 < P_3 / 50`. -/
theorem partial_product_summary_extended :
    -- brackets for the two new factors
    138 / 1000 < sinSqQ 11 ∧ sinSqQ 11 < 140 / 1000
    ∧ 125 / 1000 < sinSqQ 13 ∧ sinSqQ 13 < 127 / 1000
    -- brackets for the two new partial products
    ∧ 1018 / 1000000 < P4 ∧ P4 < 1021 / 1000000
    ∧ 127 / 1000000 < P5 ∧ P5 < 129 / 1000000
    -- full extended chain
    ∧ 0 < P5 ∧ P5 < P4 ∧ P4 < P3 ∧ P3 < P2 ∧ P2 < P1 ∧ P1 < 1
    -- echo suppression
    ∧ P5 < P3 / 50 :=
  ⟨sinSq_11_bracket.1, sinSq_11_bracket.2,
   sinSq_13_bracket.1, sinSq_13_bracket.2,
   P4_bracket.1, P4_bracket.2,
   P5_bracket.1, P5_bracket.2,
   P5_pos, P5_lt_P4, P4_lt_P3, P3_lt_P2, P2_lt_P1, P1_lt_one,
   P5_lt_P3_div_fifty⟩

end PT.Holonomy
