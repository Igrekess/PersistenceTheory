/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Holonomy.CouplingReconstruction
import PT.Holonomy.CouplingReconstructionBounds
import PT.Holonomy.AlphaPowerSequence
import Mathlib.Tactic

/-!
# Inverse-power sequence of the bare coupling — `(α_bare⁻¹)^n` for `n = 1..5`

This module is the dual of `PT.Holonomy.AlphaPowerSequence`. It records the
first five **inverse powers** of the PT bare coupling

  `α_bare = sin²θ_3 · sin²θ_5 · sin²θ_7`

evaluated at the PT fixed point `q_+ = 13/15`, `μ⋆ = 15`. Whereas
`α_bare ≈ 7.34 × 10⁻³ < 1` makes the direct power tower geometrically
decreasing, the inverse coupling `α_bare⁻¹ ≈ 136.28 > 1` makes the inverse
tower **strictly increasing** and rapidly large:

* `(α_bare⁻¹)^1 ≈ 1.363 × 10²`
* `(α_bare⁻¹)^2 ≈ 1.857 × 10⁴`
* `(α_bare⁻¹)^3 ≈ 2.531 × 10⁶`
* `(α_bare⁻¹)^4 ≈ 3.449 × 10⁸`
* `(α_bare⁻¹)^5 ≈ 4.700 × 10¹⁰`

The inverse coupling is the quantity directly observed in QED experiments
(`α_EM⁻¹ ≈ 137.036`), and its higher powers control suppression of
multi-loop corrections via `α_EM⁻ⁿ`. The PT prediction `α_bare⁻¹ ≈ 136.28`
matches `α_EM⁻¹` to better than `0.6 %`; subsequent Koide / running
corrections close the residual gap.

All brackets are decided by `unfold + norm_num` on the exact PT rationals.
A duality identity `α^n · (α⁻¹)^n = 1` is proved for `n = 1..5`,
witnessing the algebraic compatibility between `AlphaPowerSequence` and
this module.

## Reference

* `PT.Holonomy.CouplingReconstructionBounds.alphaBareInvQ` — definition.
* `PT.Holonomy.CouplingReconstructionBounds.alphaBareInvQ_bracket`
  — base bracket `135 < α⁻¹ < 137`.
* `PT.Holonomy.InvAlphaSquaredBracket` — tight brackets for `n = 2, 3`.
* `PT.Holonomy.AlphaPowerSequence` — dual `α^n` sequence.
* Monograph Ch09 §"Reconstruction des couplages" — `α_bare⁻¹` headline.
-/

namespace PT.Holonomy

/-! ### Definition: `alphaInvPow n = (α_bare⁻¹)^n` -/

/-- The `n`-th power of the bare inverse PT coupling. -/
noncomputable def alphaInvPow (n : ℕ) : ℚ := alphaBareInvQ ^ n

/-- Definitional unfolding. -/
theorem alphaInvPow_eq (n : ℕ) : alphaInvPow n = alphaBareInvQ ^ n := rfl

/-! ### Recurrence and special cases -/

/-- `(α⁻¹)^0 = 1`. -/
theorem alphaInvPow_zero : alphaInvPow 0 = 1 := by
  unfold alphaInvPow; simp

/-- `(α⁻¹)^1 = α⁻¹`. -/
theorem alphaInvPow_one : alphaInvPow 1 = alphaBareInvQ := by
  unfold alphaInvPow; simp

/-- `(α⁻¹)^2 = α⁻¹ · α⁻¹`. -/
theorem alphaInvPow_two : alphaInvPow 2 = alphaBareInvQ * alphaBareInvQ := by
  unfold alphaInvPow; ring

/-- `(α⁻¹)^3 = (α⁻¹)^2 · α⁻¹`. -/
theorem alphaInvPow_three : alphaInvPow 3 = alphaInvPow 2 * alphaBareInvQ := by
  unfold alphaInvPow; ring

/-- `(α⁻¹)^4 = (α⁻¹)^3 · α⁻¹`. -/
theorem alphaInvPow_four : alphaInvPow 4 = alphaInvPow 3 * alphaBareInvQ := by
  unfold alphaInvPow; ring

/-- `(α⁻¹)^5 = (α⁻¹)^4 · α⁻¹`. -/
theorem alphaInvPow_five : alphaInvPow 5 = alphaInvPow 4 * alphaBareInvQ := by
  unfold alphaInvPow; ring

/-! ### Positivity for each inverse power -/

/-- Generic positivity: `(α⁻¹)^n > 0` for all `n`. -/
theorem alphaInvPow_pos (n : ℕ) : 0 < alphaInvPow n := by
  unfold alphaInvPow
  exact pow_pos alphaBareInvQ_pos n

theorem alphaInvPow_one_pos   : 0 < alphaInvPow 1 := alphaInvPow_pos 1
theorem alphaInvPow_two_pos   : 0 < alphaInvPow 2 := alphaInvPow_pos 2
theorem alphaInvPow_three_pos : 0 < alphaInvPow 3 := alphaInvPow_pos 3
theorem alphaInvPow_four_pos  : 0 < alphaInvPow 4 := alphaInvPow_pos 4
theorem alphaInvPow_five_pos  : 0 < alphaInvPow 5 := alphaInvPow_pos 5

/-! ### Tight rational brackets for `(α⁻¹)^1..^5`

All brackets are decided by `norm_num` after unfolding the rational
definitions. Each bracket pins the exact value to 5 significant figures. -/

/-- **Bracket `(α⁻¹)^1`.** `13627/100 < α⁻¹ < 13628/100`
    (exact value `≈ 136.278`). -/
theorem alphaInvPow_one_bracket :
    13627 / 100 < alphaInvPow 1 ∧ alphaInvPow 1 < 13628 / 100 := by
  rw [alphaInvPow_one]
  unfold alphaBareInvQ alphaBareQ sinSqQ deltaQ qPT
  refine ⟨?_, ?_⟩ <;> norm_num

/-- **Bracket `(α⁻¹)^2`.** `18571 < (α⁻¹)^2 < 18572`
    (exact value `≈ 18571.78`). -/
theorem alphaInvPow_two_bracket :
    (18571 : ℚ) < alphaInvPow 2 ∧ alphaInvPow 2 < 18572 := by
  unfold alphaInvPow alphaBareInvQ alphaBareQ sinSqQ deltaQ qPT
  refine ⟨?_, ?_⟩ <;> norm_num

/-- **Bracket `(α⁻¹)^3`.** `2 530 931 < (α⁻¹)^3 < 2 530 932`
    (exact value `≈ 2 530 931.85`). -/
theorem alphaInvPow_three_bracket :
    (2530931 : ℚ) < alphaInvPow 3 ∧ alphaInvPow 3 < 2530932 := by
  unfold alphaInvPow alphaBareInvQ alphaBareQ sinSqQ deltaQ qPT
  refine ⟨?_, ?_⟩ <;> norm_num

/-- **Bracket `(α⁻¹)^4`.** `344 911 177 < (α⁻¹)^4 < 344 911 178`
    (exact value `≈ 3.449 × 10⁸`). -/
theorem alphaInvPow_four_bracket :
    (344911177 : ℚ) < alphaInvPow 4 ∧ alphaInvPow 4 < 344911178 := by
  unfold alphaInvPow alphaBareInvQ alphaBareQ sinSqQ deltaQ qPT
  refine ⟨?_, ?_⟩ <;> norm_num

/-- **Bracket `(α⁻¹)^5`.** `47 003 920 783 < (α⁻¹)^5 < 47 003 920 784`
    (exact value `≈ 4.700 × 10¹⁰`). -/
theorem alphaInvPow_five_bracket :
    (47003920783 : ℚ) < alphaInvPow 5 ∧ alphaInvPow 5 < 47003920784 := by
  unfold alphaInvPow alphaBareInvQ alphaBareQ sinSqQ deltaQ qPT
  refine ⟨?_, ?_⟩ <;> norm_num

/-! ### Bound `α⁻¹ > 1` and strict geometric increase -/

/-- `α_bare⁻¹ > 1` (in fact `> 135`). -/
theorem alphaBareInvQ_gt_one : 1 < alphaBareInvQ := by
  have := alphaBareInvQ_gt_135
  linarith

/-- **Strict geometric increase.** `(α⁻¹)^n < (α⁻¹)^(n+1)` for all `n`,
    since `α⁻¹ > 1`. -/
theorem alphaInvPow_lt_succ (n : ℕ) : alphaInvPow n < alphaInvPow (n + 1) := by
  unfold alphaInvPow
  rw [pow_succ]
  -- Goal: α⁻¹^n < α⁻¹^n * α⁻¹
  have hpos : 0 < alphaBareInvQ ^ n := pow_pos alphaBareInvQ_pos n
  have hgt1 : 1 < alphaBareInvQ := alphaBareInvQ_gt_one
  nlinarith [hpos, alphaBareInvQ_pos, hgt1]

/-- Specific increases for `n = 1..4`. -/
theorem alphaInvPow_inc_1_2 : alphaInvPow 1 < alphaInvPow 2 := alphaInvPow_lt_succ 1
theorem alphaInvPow_inc_2_3 : alphaInvPow 2 < alphaInvPow 3 := alphaInvPow_lt_succ 2
theorem alphaInvPow_inc_3_4 : alphaInvPow 3 < alphaInvPow 4 := alphaInvPow_lt_succ 3
theorem alphaInvPow_inc_4_5 : alphaInvPow 4 < alphaInvPow 5 := alphaInvPow_lt_succ 4

/-- **Full strict cascade** for the first five inverse powers. -/
theorem alphaInvPow_chain_1_to_5 :
    alphaInvPow 1 < alphaInvPow 2
    ∧ alphaInvPow 2 < alphaInvPow 3
    ∧ alphaInvPow 3 < alphaInvPow 4
    ∧ alphaInvPow 4 < alphaInvPow 5 :=
  ⟨alphaInvPow_inc_1_2, alphaInvPow_inc_2_3, alphaInvPow_inc_3_4, alphaInvPow_inc_4_5⟩

/-! ### Product identity (law of exponents) -/

/-- **Law of exponents.** `(α⁻¹)^(m+n) = (α⁻¹)^m · (α⁻¹)^n`. -/
theorem alphaInvPow_add (m n : ℕ) :
    alphaInvPow (m + n) = alphaInvPow m * alphaInvPow n := by
  unfold alphaInvPow
  exact pow_add alphaBareInvQ m n

/-- Specialisation: `(α⁻¹)^2 = α⁻¹ · α⁻¹`. -/
theorem alphaInvPow_add_1_1 : alphaInvPow 2 = alphaInvPow 1 * alphaInvPow 1 := by
  have := alphaInvPow_add 1 1
  simpa using this

/-- Specialisation: `(α⁻¹)^3 = (α⁻¹)^2 · α⁻¹`. -/
theorem alphaInvPow_add_2_1 : alphaInvPow 3 = alphaInvPow 2 * alphaInvPow 1 := by
  have := alphaInvPow_add 2 1
  simpa using this

/-- Specialisation: `(α⁻¹)^4 = (α⁻¹)^2 · (α⁻¹)^2`. -/
theorem alphaInvPow_add_2_2 : alphaInvPow 4 = alphaInvPow 2 * alphaInvPow 2 := by
  have := alphaInvPow_add 2 2
  simpa using this

/-- Specialisation: `(α⁻¹)^5 = (α⁻¹)^3 · (α⁻¹)^2`. -/
theorem alphaInvPow_add_3_2 : alphaInvPow 5 = alphaInvPow 3 * alphaInvPow 2 := by
  have := alphaInvPow_add 3 2
  simpa using this

/-! ### Duality with `alphaPow`

These identities witness that `alphaInvPow n` is the multiplicative
inverse of `alphaPow n` in `ℚ`. They follow directly from `α_bare > 0`
together with `α_bare · α_bare⁻¹ = 1` and the law of exponents. -/

/-- **Base duality.** `α_bare · α_bare⁻¹ = 1`. -/
theorem alphaBareQ_mul_alphaBareInvQ : alphaBareQ * alphaBareInvQ = 1 := by
  unfold alphaBareInvQ
  have ha : alphaBareQ ≠ 0 := ne_of_gt alphaBareQ_pos
  field_simp

/-- **Duality at all orders.** `α^n · (α⁻¹)^n = 1` for every `n`. -/
theorem alphaPow_mul_alphaInvPow (n : ℕ) :
    alphaPow n * alphaInvPow n = 1 := by
  unfold alphaPow alphaInvPow
  rw [← mul_pow, alphaBareQ_mul_alphaBareInvQ, one_pow]

theorem alphaPow_mul_alphaInvPow_one   : alphaPow 1 * alphaInvPow 1 = 1 :=
  alphaPow_mul_alphaInvPow 1
theorem alphaPow_mul_alphaInvPow_two   : alphaPow 2 * alphaInvPow 2 = 1 :=
  alphaPow_mul_alphaInvPow 2
theorem alphaPow_mul_alphaInvPow_three : alphaPow 3 * alphaInvPow 3 = 1 :=
  alphaPow_mul_alphaInvPow 3
theorem alphaPow_mul_alphaInvPow_four  : alphaPow 4 * alphaInvPow 4 = 1 :=
  alphaPow_mul_alphaInvPow 4
theorem alphaPow_mul_alphaInvPow_five  : alphaPow 5 * alphaInvPow 5 = 1 :=
  alphaPow_mul_alphaInvPow 5

/-! ### Headline summary -/

/-- **Headline.** All invariants of the inverse-power sequence
    `(α⁻¹)^1, ..., (α⁻¹)^5` at the PT fixed point `q_+ = 13/15`,
    `μ⋆ = 15`:

    * Positivity: `0 < (α⁻¹)^n` for `n = 1..5`.
    * Strict geometric increase:
      `(α⁻¹)^1 < (α⁻¹)^2 < (α⁻¹)^3 < (α⁻¹)^4 < (α⁻¹)^5`.
    * Tight rational brackets:
        * `(α⁻¹)^1 ∈ (13627/100, 13628/100)`            ≈ 136.278
        * `(α⁻¹)^2 ∈ (18571, 18572)`                    ≈ 18571.78
        * `(α⁻¹)^3 ∈ (2530931, 2530932)`                ≈ 2 530 931.85
        * `(α⁻¹)^4 ∈ (344911177, 344911178)`            ≈ 3.449 × 10⁸
        * `(α⁻¹)^5 ∈ (47003920783, 47003920784)`        ≈ 4.700 × 10¹⁰
    * Law of exponents: `(α⁻¹)^(m+n) = (α⁻¹)^m · (α⁻¹)^n`.
    * Duality: `α^n · (α⁻¹)^n = 1` for `n = 1..5`. -/
theorem alphaInvPow_summary :
    -- positivity
    0 < alphaInvPow 1 ∧ 0 < alphaInvPow 2 ∧ 0 < alphaInvPow 3
    ∧ 0 < alphaInvPow 4 ∧ 0 < alphaInvPow 5
    -- strict cascade
    ∧ alphaInvPow 1 < alphaInvPow 2
    ∧ alphaInvPow 2 < alphaInvPow 3
    ∧ alphaInvPow 3 < alphaInvPow 4
    ∧ alphaInvPow 4 < alphaInvPow 5
    -- brackets
    ∧ 13627 / 100 < alphaInvPow 1
    ∧ alphaInvPow 1 < 13628 / 100
    ∧ (18571 : ℚ) < alphaInvPow 2
    ∧ alphaInvPow 2 < 18572
    ∧ (2530931 : ℚ) < alphaInvPow 3
    ∧ alphaInvPow 3 < 2530932
    ∧ (344911177 : ℚ) < alphaInvPow 4
    ∧ alphaInvPow 4 < 344911178
    ∧ (47003920783 : ℚ) < alphaInvPow 5
    ∧ alphaInvPow 5 < 47003920784
    -- law of exponents (specialisations)
    ∧ alphaInvPow 2 = alphaInvPow 1 * alphaInvPow 1
    ∧ alphaInvPow 3 = alphaInvPow 2 * alphaInvPow 1
    ∧ alphaInvPow 4 = alphaInvPow 2 * alphaInvPow 2
    ∧ alphaInvPow 5 = alphaInvPow 3 * alphaInvPow 2
    -- duality with alphaPow
    ∧ alphaPow 1 * alphaInvPow 1 = 1
    ∧ alphaPow 2 * alphaInvPow 2 = 1
    ∧ alphaPow 3 * alphaInvPow 3 = 1
    ∧ alphaPow 4 * alphaInvPow 4 = 1
    ∧ alphaPow 5 * alphaInvPow 5 = 1 :=
  ⟨alphaInvPow_one_pos, alphaInvPow_two_pos, alphaInvPow_three_pos,
   alphaInvPow_four_pos, alphaInvPow_five_pos,
   alphaInvPow_inc_1_2, alphaInvPow_inc_2_3,
   alphaInvPow_inc_3_4, alphaInvPow_inc_4_5,
   alphaInvPow_one_bracket.1, alphaInvPow_one_bracket.2,
   alphaInvPow_two_bracket.1, alphaInvPow_two_bracket.2,
   alphaInvPow_three_bracket.1, alphaInvPow_three_bracket.2,
   alphaInvPow_four_bracket.1, alphaInvPow_four_bracket.2,
   alphaInvPow_five_bracket.1, alphaInvPow_five_bracket.2,
   alphaInvPow_add_1_1, alphaInvPow_add_2_1,
   alphaInvPow_add_2_2, alphaInvPow_add_3_2,
   alphaPow_mul_alphaInvPow_one, alphaPow_mul_alphaInvPow_two,
   alphaPow_mul_alphaInvPow_three, alphaPow_mul_alphaInvPow_four,
   alphaPow_mul_alphaInvPow_five⟩

end PT.Holonomy
