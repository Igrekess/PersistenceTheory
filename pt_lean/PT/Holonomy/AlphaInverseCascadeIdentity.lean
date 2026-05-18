/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Holonomy.CouplingReconstruction
import PT.Holonomy.CouplingReconstructionBounds
import PT.Holonomy.InverseSinSq
import PT.Holonomy.InverseSinSqProduct
import PT.Holonomy.SinSqProductChain
import Mathlib.Tactic

/-!
# `α_bare⁻¹` Cascade Identity (structural exhaustion)

This file consolidates the **product form** of the inverse fine-structure
coupling at the PT fixed point `q_+ = 13/15` and exposes its structural
variants exhaustively. The starting identity is

  `α_bare⁻¹ = (1/sin²θ_3) · (1/sin²θ_5) · (1/sin²θ_7)`

already established in `InverseSinSqProduct` as
`invProductActive_eq_alphaBareInv`. We add:

1. **Primary identity (named alias).**
   `alphaBareInv_eq_invSinSqProd : α_bare⁻¹ = ∏_{p∈{3,5,7}} 1/sin²θ_p`.

2. **5-prime extension.**
   `oneOverP5_eq_invSinSqChain : 1/P_5 = ∏_{p∈{3,5,7,11,13}} 1/sin²θ_p`.

3. **Recombination variants.** The 5-factor inverse product decomposes
   as `α_bare⁻¹ · (1/sin²θ_11) · (1/sin²θ_13)`, hence:
   * `oneOverP5_eq_alphaBareInv_times_echo`
   * `alphaBareInv_eq_invSinSq_5_times_7_times_sinSq_3` (3-factor → 2-factor
     by division by `(1/sin²θ_3)`, i.e. multiplication by `sin²θ_3`).

4. **Brackets** for the new derived quantities:
   * ratio `(1/P_5) / α_bare⁻¹ = (1/sin²θ_11)(1/sin²θ_13) ≈ 57.26`,
   * mean per active prime `α_bare⁻¹ / 3 ≈ 45.43`,
   * 2-factor head `(1/sin²θ_5)(1/sin²θ_7) = α_bare⁻¹ · sin²θ_3 ≈ 29.87`.

5. **Headline summary** collecting every invariant of the cascade.

All proofs are exact rational arithmetic (`field_simp` for identities,
`norm_num` for brackets). No analytic ingredient.

## Reference

Monograph Chapter 9, §"Reconstruction des couplages" (BA5 Step 4 — inverse
form). Builds on `CouplingReconstruction`, `CouplingReconstructionBounds`,
`InverseSinSq`, `InverseSinSqProduct`, `SinSqProductChain`.
-/

namespace PT.Holonomy

/-! ### Primary identity (active 3-prime cascade) -/

/-- **Primary identity.**
    `α_bare⁻¹ = (1/sin²θ_3) · (1/sin²θ_5) · (1/sin²θ_7)`. -/
theorem alphaBareInv_eq_invSinSqProd :
    alphaBareInvQ = invSinSqQ 3 * invSinSqQ 5 * invSinSqQ 7 := by
  -- `invProductActive` was defined exactly this way; reverse the existing identity.
  have h := invProductActive_eq_alphaBareInv
  -- `h : invProductActive = alphaBareInvQ`
  -- `invProductActive = invSinSqQ 3 * invSinSqQ 5 * invSinSqQ 7` by `rfl`.
  unfold invProductActive at h
  exact h.symm

/-! ### 5-prime extension over the echo block -/

/-- **5-prime inverse identity.**
    `1 / P_5 = (1/sin²θ_3) · (1/sin²θ_5) · (1/sin²θ_7) · (1/sin²θ_11) · (1/sin²θ_13)`. -/
theorem oneOverP5_eq_invSinSqChain :
    1 / P5 = invSinSqQ 3 * invSinSqQ 5 * invSinSqQ 7
              * invSinSqQ 11 * invSinSqQ 13 := by
  have h := invProductChain_eq_inv_P5
  unfold invProductChain at h
  exact h.symm

/-! ### Recombination variants -/

/-- **Echo factorisation.**
    The 5-prime inverse product equals `α_bare⁻¹` multiplied by the two
    echo inverse factors:
    `1 / P_5 = α_bare⁻¹ · (1/sin²θ_11) · (1/sin²θ_13)`. -/
theorem oneOverP5_eq_alphaBareInv_times_echo :
    1 / P5 = alphaBareInvQ * invSinSqQ 11 * invSinSqQ 13 := by
  rw [oneOverP5_eq_invSinSqChain, alphaBareInv_eq_invSinSqProd]

/-- **Ratio identity.** Dividing the 5-prime product by the 3-prime product
    leaves exactly the two echo factors:
    `(1 / P_5) / α_bare⁻¹ = (1/sin²θ_11) · (1/sin²θ_13)`. -/
theorem oneOverP5_div_alphaBareInv :
    (1 / P5) / alphaBareInvQ = invSinSqQ 11 * invSinSqQ 13 := by
  rw [oneOverP5_eq_alphaBareInv_times_echo]
  have h : alphaBareInvQ ≠ 0 := ne_of_gt alphaBareInvQ_pos
  field_simp

/-- **3-factor → 2-factor reduction.**
    Multiplying `α_bare⁻¹` by `sin²θ_3` cancels the leading inverse factor:
    `α_bare⁻¹ · sin²θ_3 = (1/sin²θ_5) · (1/sin²θ_7)`. -/
theorem alphaBareInv_mul_sinSq3_eq_invSinSq_5_7 :
    alphaBareInvQ * sinSqQ 3 = invSinSqQ 5 * invSinSqQ 7 := by
  rw [alphaBareInv_eq_invSinSqProd]
  unfold invSinSqQ
  have h3 : sinSqQ 3 ≠ 0 := ne_of_gt sinSq_3_pos
  field_simp

/-- **3-factor → 1-factor reduction.**
    Multiplying `α_bare⁻¹` by `sin²θ_3 · sin²θ_5` cancels two inverse factors:
    `α_bare⁻¹ · sin²θ_3 · sin²θ_5 = 1/sin²θ_7`. -/
theorem alphaBareInv_mul_sinSq3_sinSq5_eq_invSinSq_7 :
    alphaBareInvQ * sinSqQ 3 * sinSqQ 5 = invSinSqQ 7 := by
  rw [alphaBareInv_eq_invSinSqProd]
  unfold invSinSqQ
  have h3 : sinSqQ 3 ≠ 0 := ne_of_gt sinSq_3_pos
  have h5 : sinSqQ 5 ≠ 0 := ne_of_gt sinSq_5_pos
  field_simp

/-- **Single-prime extraction.** Dividing `α_bare⁻¹` by any of its three
    factors leaves the product of the other two:
    `α_bare⁻¹ / (1/sin²θ_7) = (1/sin²θ_3) · (1/sin²θ_5)`. -/
theorem alphaBareInv_div_invSinSq_7 :
    alphaBareInvQ / invSinSqQ 7 = invSinSqQ 3 * invSinSqQ 5 := by
  rw [alphaBareInv_eq_invSinSqProd]
  have h7 : invSinSqQ 7 ≠ 0 := ne_of_gt invSinSq_7_pos
  field_simp

/-! ### Positivity of derived quantities -/

/-- The ratio `(1/P_5) / α_bare⁻¹` is strictly positive. -/
theorem oneOverP5_div_alphaBareInv_pos :
    0 < (1 / P5) / alphaBareInvQ := by
  rw [oneOverP5_div_alphaBareInv]
  exact mul_pos invSinSq_11_pos invSinSq_13_pos

/-- The mean of `α_bare⁻¹` over the three active primes is strictly positive. -/
theorem alphaBareInv_div_three_pos :
    0 < alphaBareInvQ / 3 := by
  have hpos : 0 < alphaBareInvQ := alphaBareInvQ_pos
  positivity

/-- The 2-factor head `(1/sin²θ_5)(1/sin²θ_7)` is strictly positive. -/
theorem invSinSq_5_mul_7_pos : 0 < invSinSqQ 5 * invSinSqQ 7 :=
  mul_pos invSinSq_5_pos invSinSq_7_pos

/-! ### Decimal brackets for derived ratios

Calibrated by `#eval` (rational → float) on the closed-form expressions:

* `(1/sin²θ_11) · (1/sin²θ_13) ≈ 57.2598`,
* `α_bare⁻¹ / 3 ≈ 45.4261`,
* `(1/sin²θ_5) · (1/sin²θ_7) ≈ 29.8661`. -/

/-- **Ratio bracket.** `(1/sin²θ_11)·(1/sin²θ_13) ∈ (57.25, 57.27)`. -/
theorem invSinSq_11_mul_13_bracket :
    5725 / 100 < invSinSqQ 11 * invSinSqQ 13
    ∧ invSinSqQ 11 * invSinSqQ 13 < 5727 / 100 := by
  unfold invSinSqQ sinSqQ deltaQ qPT
  refine ⟨?_, ?_⟩ <;> norm_num

/-- **Ratio bracket (named via the recombination identity).**
    `(1/P_5)/α_bare⁻¹ ∈ (57.25, 57.27)`. -/
theorem oneOverP5_div_alphaBareInv_bracket :
    5725 / 100 < (1 / P5) / alphaBareInvQ
    ∧ (1 / P5) / alphaBareInvQ < 5727 / 100 := by
  rw [oneOverP5_div_alphaBareInv]
  exact invSinSq_11_mul_13_bracket

/-- **Mean bracket.** `α_bare⁻¹ / 3 ∈ (45.42, 45.43)`
    (exact value `≈ 45.4261`). -/
theorem alphaBareInv_div_three_bracket :
    4542 / 100 < alphaBareInvQ / 3
    ∧ alphaBareInvQ / 3 < 4543 / 100 := by
  unfold alphaBareInvQ alphaBareQ sinSqQ deltaQ qPT
  refine ⟨?_, ?_⟩ <;> norm_num

/-- **2-factor head bracket.** `(1/sin²θ_5)·(1/sin²θ_7) ∈ (29.86, 29.87)`. -/
theorem invSinSq_5_mul_7_bracket :
    2986 / 100 < invSinSqQ 5 * invSinSqQ 7
    ∧ invSinSqQ 5 * invSinSqQ 7 < 2987 / 100 := by
  unfold invSinSqQ sinSqQ deltaQ qPT
  refine ⟨?_, ?_⟩ <;> norm_num

/-- **2-factor head bracket (named via the reduction identity).**
    `α_bare⁻¹ · sin²θ_3 ∈ (29.86, 29.87)`. -/
theorem alphaBareInv_mul_sinSq3_bracket :
    2986 / 100 < alphaBareInvQ * sinSqQ 3
    ∧ alphaBareInvQ * sinSqQ 3 < 2987 / 100 := by
  rw [alphaBareInv_mul_sinSq3_eq_invSinSq_5_7]
  exact invSinSq_5_mul_7_bracket

/-! ### Headline -/

/-- **Headline (α_bare⁻¹ cascade identity, full structural exhaustion).**

    At the PT fixed point `q_+ = 13/15`:

    * **Primary** `α_bare⁻¹ = (1/sin²θ_3)·(1/sin²θ_5)·(1/sin²θ_7)`,
      bracketed `136 < α_bare⁻¹ < 137` (exact `≈ 136.278`).
    * **5-prime extension** `1/P_5 = ∏_{p∈{3,5,7,11,13}} 1/sin²θ_p`,
      bracketed `7800 < 1/P_5 < 7810` (exact `≈ 7803.27`).
    * **Echo factorisation** `1/P_5 = α_bare⁻¹ · (1/sin²θ_11)·(1/sin²θ_13)`.
    * **Echo ratio** `(1/P_5)/α_bare⁻¹ = (1/sin²θ_11)·(1/sin²θ_13)
      ∈ (57.25, 57.27)` (exact `≈ 57.260`).
    * **Mean per active prime** `α_bare⁻¹ / 3 ∈ (45.42, 45.43)`
      (exact `≈ 45.426`).
    * **3→2 reduction** `α_bare⁻¹ · sin²θ_3 = (1/sin²θ_5)·(1/sin²θ_7)
      ∈ (29.86, 29.87)` (exact `≈ 29.866`).
    * **3→1 reduction** `α_bare⁻¹ · sin²θ_3 · sin²θ_5 = 1/sin²θ_7`.
    * Strict monotone enlargement `α_bare⁻¹ < 1/P_5` (echo factors > 1). -/
theorem alphaBareInv_cascade_summary :
    -- primary identity (3-prime)
    alphaBareInvQ = invSinSqQ 3 * invSinSqQ 5 * invSinSqQ 7
    -- 5-prime extension
    ∧ 1 / P5 = invSinSqQ 3 * invSinSqQ 5 * invSinSqQ 7
                * invSinSqQ 11 * invSinSqQ 13
    -- echo factorisation
    ∧ 1 / P5 = alphaBareInvQ * invSinSqQ 11 * invSinSqQ 13
    -- echo ratio identity
    ∧ (1 / P5) / alphaBareInvQ = invSinSqQ 11 * invSinSqQ 13
    -- 3→2 reduction
    ∧ alphaBareInvQ * sinSqQ 3 = invSinSqQ 5 * invSinSqQ 7
    -- 3→1 reduction
    ∧ alphaBareInvQ * sinSqQ 3 * sinSqQ 5 = invSinSqQ 7
    -- single-prime extraction
    ∧ alphaBareInvQ / invSinSqQ 7 = invSinSqQ 3 * invSinSqQ 5
    -- positivity
    ∧ 0 < alphaBareInvQ ∧ 0 < (1 : ℚ) / P5
    ∧ 0 < (1 / P5) / alphaBareInvQ
    ∧ 0 < alphaBareInvQ / 3
    -- brackets for the primary quantities (inherited)
    ∧ 136 < alphaBareInvQ ∧ alphaBareInvQ < 137
    ∧ 7800 < (1 : ℚ) / P5 ∧ (1 : ℚ) / P5 < 7810
    -- brackets for the new derived quantities
    ∧ 5725 / 100 < (1 / P5) / alphaBareInvQ
    ∧ (1 / P5) / alphaBareInvQ < 5727 / 100
    ∧ 4542 / 100 < alphaBareInvQ / 3
    ∧ alphaBareInvQ / 3 < 4543 / 100
    ∧ 2986 / 100 < alphaBareInvQ * sinSqQ 3
    ∧ alphaBareInvQ * sinSqQ 3 < 2987 / 100
    -- strict monotone enlargement
    ∧ alphaBareInvQ < 1 / P5 := by
  refine ⟨alphaBareInv_eq_invSinSqProd,
          oneOverP5_eq_invSinSqChain,
          oneOverP5_eq_alphaBareInv_times_echo,
          oneOverP5_div_alphaBareInv,
          alphaBareInv_mul_sinSq3_eq_invSinSq_5_7,
          alphaBareInv_mul_sinSq3_sinSq5_eq_invSinSq_7,
          alphaBareInv_div_invSinSq_7,
          alphaBareInvQ_pos,
          one_div_pos.mpr P5_pos,
          oneOverP5_div_alphaBareInv_pos,
          alphaBareInv_div_three_pos,
          invProductActive_bracket_tight.1.trans_eq
            invProductActive_eq_alphaBareInv,
          ?_,
          ?_,
          ?_,
          oneOverP5_div_alphaBareInv_bracket.1,
          oneOverP5_div_alphaBareInv_bracket.2,
          alphaBareInv_div_three_bracket.1,
          alphaBareInv_div_three_bracket.2,
          alphaBareInv_mul_sinSq3_bracket.1,
          alphaBareInv_mul_sinSq3_bracket.2,
          ?_⟩
  · -- alphaBareInvQ < 137
    exact alphaBareInvQ_lt_137
  · -- 7800 < 1/P5
    rw [← invProductChain_eq_inv_P5]
    exact invProductChain_bracket_tight.1
  · -- 1/P5 < 7810
    rw [← invProductChain_eq_inv_P5]
    exact invProductChain_bracket_tight.2
  · -- alphaBareInvQ < 1/P5
    rw [← invProductActive_eq_alphaBareInv, ← invProductChain_eq_inv_P5]
    exact invProductActive_lt_invProductChain

end PT.Holonomy
