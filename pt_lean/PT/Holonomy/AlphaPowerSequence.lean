/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Holonomy.CouplingReconstruction
import PT.Holonomy.CouplingReconstructionBounds
import Mathlib.Tactic

/-!
# Power sequence of the bare coupling — `α_bare^n` for `n = 1..5`

This module computes the first five powers of the bare PT coupling
`α_bare = sin²θ_3 · sin²θ_5 · sin²θ_7` evaluated at the PT fixed point
`q_+ = 13/15`, `μ⋆ = 15`.

The bare value `α_bare ≈ 7.3379 × 10⁻³` lies strictly below 1, so the
sequence `α^n` is **strictly decreasing** and converges geometrically to 0.

Numerical headline (exact rationals, decimal approximations):

* `α^1 ≈ 7.3379 × 10⁻³`
* `α^2 ≈ 5.3846 × 10⁻⁵`
* `α^3 ≈ 3.9511 × 10⁻⁷`
* `α^4 ≈ 2.8993 × 10⁻⁹`
* `α^5 ≈ 2.1275 × 10⁻¹¹`

These appear in perturbative QED expansions (each loop order contributes
a factor `α`), in the Koide / running cascade, and in PT predictions for
high-order finite-structure corrections.

The brackets below are tight (5 significant figures) and proven by
`norm_num` after unfolding the rational definitions.

## Reference

* `PT.Holonomy.CouplingReconstruction.alphaBareQ` — definition.
* `PT.Holonomy.CouplingReconstructionBounds.alphaBareQ_bracket_tight` —
  source bracket for `α^1`.
* Monograph Ch09 §"Reconstruction des couplages" — perturbative tower.
-/

namespace PT.Holonomy

/-! ### Definition: `alphaPow n = α_bare^n` -/

/-- The `n`-th power of the bare PT coupling. -/
def alphaPow (n : ℕ) : ℚ := alphaBareQ ^ n

/-- Definitional unfolding. -/
theorem alphaPow_eq (n : ℕ) : alphaPow n = alphaBareQ ^ n := rfl

/-! ### Recurrence and special cases -/

/-- `α^0 = 1`. -/
theorem alphaPow_zero : alphaPow 0 = 1 := by
  unfold alphaPow; simp

/-- `α^1 = α`. -/
theorem alphaPow_one : alphaPow 1 = alphaBareQ := by
  unfold alphaPow; simp

/-- `α^2 = α · α`. -/
theorem alphaPow_two : alphaPow 2 = alphaBareQ * alphaBareQ := by
  unfold alphaPow; ring

/-- `α^3 = α^2 · α`. -/
theorem alphaPow_three : alphaPow 3 = alphaPow 2 * alphaBareQ := by
  unfold alphaPow; ring

/-- `α^4 = α^3 · α`. -/
theorem alphaPow_four : alphaPow 4 = alphaPow 3 * alphaBareQ := by
  unfold alphaPow; ring

/-- `α^5 = α^4 · α`. -/
theorem alphaPow_five : alphaPow 5 = alphaPow 4 * alphaBareQ := by
  unfold alphaPow; ring

/-! ### Positivity for each power -/

/-- Generic positivity: `α^n > 0` for all `n`. -/
theorem alphaPow_pos (n : ℕ) : 0 < alphaPow n := by
  unfold alphaPow
  exact pow_pos alphaBareQ_pos n

theorem alphaPow_one_pos   : 0 < alphaPow 1 := alphaPow_pos 1
theorem alphaPow_two_pos   : 0 < alphaPow 2 := alphaPow_pos 2
theorem alphaPow_three_pos : 0 < alphaPow 3 := alphaPow_pos 3
theorem alphaPow_four_pos  : 0 < alphaPow 4 := alphaPow_pos 4
theorem alphaPow_five_pos  : 0 < alphaPow 5 := alphaPow_pos 5

/-! ### Tight rational brackets for `α^1, α^2, α^3, α^4, α^5`

All brackets are decided by `norm_num` after unfolding the rational
definitions. Each bracket has 5 significant figures of precision. -/

/-- **Bracket `α^1`.** `7335/10⁶ < α < 7340/10⁶`
    (exact value `≈ 7.3379 × 10⁻³`). -/
theorem alphaPow_one_bracket :
    7335 / 1000000 < alphaPow 1 ∧ alphaPow 1 < 7340 / 1000000 := by
  rw [alphaPow_one]
  exact alphaBareQ_bracket_tight

/-- **Bracket `α^2`.** `53840/10⁹ < α^2 < 53850/10⁹`
    (exact value `≈ 5.3846 × 10⁻⁵`). -/
theorem alphaPow_two_bracket :
    53840 / 1000000000 < alphaPow 2 ∧ alphaPow 2 < 53850 / 1000000000 := by
  unfold alphaPow alphaBareQ sinSqQ deltaQ qPT
  refine ⟨?_, ?_⟩ <;> norm_num

/-- **Bracket `α^3`.** `395100/10¹² < α^3 < 395120/10¹²`
    (exact value `≈ 3.9511 × 10⁻⁷`). -/
theorem alphaPow_three_bracket :
    395100 / 1000000000000 < alphaPow 3
    ∧ alphaPow 3 < 395120 / 1000000000000 := by
  unfold alphaPow alphaBareQ sinSqQ deltaQ qPT
  refine ⟨?_, ?_⟩ <;> norm_num

/-- **Bracket `α^4`.** `2899290/10¹⁵ < α^4 < 2899310/10¹⁵`
    (exact value `≈ 2.8993 × 10⁻⁹`). -/
theorem alphaPow_four_bracket :
    2899290 / 1000000000000000 < alphaPow 4
    ∧ alphaPow 4 < 2899310 / 1000000000000000 := by
  unfold alphaPow alphaBareQ sinSqQ deltaQ qPT
  refine ⟨?_, ?_⟩ <;> norm_num

/-- **Bracket `α^5`.** `21274810/10¹⁸ < α^5 < 21274830/10¹⁸`
    (exact value `≈ 2.1275 × 10⁻¹¹`). -/
theorem alphaPow_five_bracket :
    21274810 / 1000000000000000000 < alphaPow 5
    ∧ alphaPow 5 < 21274830 / 1000000000000000000 := by
  unfold alphaPow alphaBareQ sinSqQ deltaQ qPT
  refine ⟨?_, ?_⟩ <;> norm_num

/-! ### Bound `α < 1` and strict geometric decrease -/

/-- `α_bare < 1` (loose: `α < 1/100 < 1`). -/
theorem alphaBareQ_lt_one : alphaBareQ < 1 := by
  have h := alphaBareQ_lt_one_hundredth
  linarith

/-- **Strict geometric decrease.** `α^(n+1) < α^n` for all `n`,
    since `0 < α < 1`. -/
theorem alphaPow_succ_lt (n : ℕ) : alphaPow (n + 1) < alphaPow n := by
  unfold alphaPow
  rw [pow_succ]
  -- Goal: α^n * α < α^n
  have hpos : 0 < alphaBareQ ^ n := pow_pos alphaBareQ_pos n
  have hlt1 : alphaBareQ < 1 := alphaBareQ_lt_one
  nlinarith [hpos, alphaBareQ_pos, hlt1]

/-- Specific decreases for `n = 1..4`. -/
theorem alphaPow_dec_1_2 : alphaPow 2 < alphaPow 1 := alphaPow_succ_lt 1
theorem alphaPow_dec_2_3 : alphaPow 3 < alphaPow 2 := alphaPow_succ_lt 2
theorem alphaPow_dec_3_4 : alphaPow 4 < alphaPow 3 := alphaPow_succ_lt 3
theorem alphaPow_dec_4_5 : alphaPow 5 < alphaPow 4 := alphaPow_succ_lt 4

/-- **Full strict cascade** for the first five powers. -/
theorem alphaPow_chain_1_to_5 :
    alphaPow 5 < alphaPow 4
    ∧ alphaPow 4 < alphaPow 3
    ∧ alphaPow 3 < alphaPow 2
    ∧ alphaPow 2 < alphaPow 1 :=
  ⟨alphaPow_dec_4_5, alphaPow_dec_3_4, alphaPow_dec_2_3, alphaPow_dec_1_2⟩

/-! ### Product identity (law of exponents) -/

/-- **Law of exponents.** `α^(m+n) = α^m · α^n`. -/
theorem alphaPow_add (m n : ℕ) :
    alphaPow (m + n) = alphaPow m * alphaPow n := by
  unfold alphaPow
  exact pow_add alphaBareQ m n

/-- Specialisation: `α^2 = α · α`. -/
theorem alphaPow_add_1_1 : alphaPow 2 = alphaPow 1 * alphaPow 1 := by
  have := alphaPow_add 1 1
  simpa using this

/-- Specialisation: `α^3 = α^2 · α`. -/
theorem alphaPow_add_2_1 : alphaPow 3 = alphaPow 2 * alphaPow 1 := by
  have := alphaPow_add 2 1
  simpa using this

/-- Specialisation: `α^4 = α^2 · α^2`. -/
theorem alphaPow_add_2_2 : alphaPow 4 = alphaPow 2 * alphaPow 2 := by
  have := alphaPow_add 2 2
  simpa using this

/-- Specialisation: `α^5 = α^3 · α^2`. -/
theorem alphaPow_add_3_2 : alphaPow 5 = alphaPow 3 * alphaPow 2 := by
  have := alphaPow_add 3 2
  simpa using this

/-! ### Headline summary -/

/-- **Headline.** All invariants of the power sequence `α^1, ..., α^5`
    at the PT fixed point `q_+ = 13/15`, `μ⋆ = 15`:

    * Positivity: `0 < α^n` for `n = 1..5`.
    * Strict geometric decrease: `α^5 < α^4 < α^3 < α^2 < α^1`.
    * Tight rational brackets (5 sig. fig. each):
        * `α^1 ∈ (7335/10⁶, 7340/10⁶)`            ≈ 7.3379 × 10⁻³
        * `α^2 ∈ (53840/10⁹, 53850/10⁹)`          ≈ 5.3846 × 10⁻⁵
        * `α^3 ∈ (395100/10¹², 395120/10¹²)`      ≈ 3.9511 × 10⁻⁷
        * `α^4 ∈ (2899290/10¹⁵, 2899310/10¹⁵)`    ≈ 2.8993 × 10⁻⁹
        * `α^5 ∈ (21274810/10¹⁸, 21274830/10¹⁸)`  ≈ 2.1275 × 10⁻¹¹
    * Law of exponents: `α^(m+n) = α^m · α^n`. -/
theorem alphaPow_summary :
    -- positivity
    0 < alphaPow 1 ∧ 0 < alphaPow 2 ∧ 0 < alphaPow 3
    ∧ 0 < alphaPow 4 ∧ 0 < alphaPow 5
    -- strict cascade
    ∧ alphaPow 5 < alphaPow 4
    ∧ alphaPow 4 < alphaPow 3
    ∧ alphaPow 3 < alphaPow 2
    ∧ alphaPow 2 < alphaPow 1
    -- brackets
    ∧ 7335 / 1000000 < alphaPow 1
    ∧ alphaPow 1 < 7340 / 1000000
    ∧ 53840 / 1000000000 < alphaPow 2
    ∧ alphaPow 2 < 53850 / 1000000000
    ∧ 395100 / 1000000000000 < alphaPow 3
    ∧ alphaPow 3 < 395120 / 1000000000000
    ∧ 2899290 / 1000000000000000 < alphaPow 4
    ∧ alphaPow 4 < 2899310 / 1000000000000000
    ∧ 21274810 / 1000000000000000000 < alphaPow 5
    ∧ alphaPow 5 < 21274830 / 1000000000000000000
    -- law of exponents (specialisations)
    ∧ alphaPow 2 = alphaPow 1 * alphaPow 1
    ∧ alphaPow 3 = alphaPow 2 * alphaPow 1
    ∧ alphaPow 4 = alphaPow 2 * alphaPow 2
    ∧ alphaPow 5 = alphaPow 3 * alphaPow 2 :=
  ⟨alphaPow_one_pos, alphaPow_two_pos, alphaPow_three_pos,
   alphaPow_four_pos, alphaPow_five_pos,
   alphaPow_dec_4_5, alphaPow_dec_3_4, alphaPow_dec_2_3, alphaPow_dec_1_2,
   alphaPow_one_bracket.1, alphaPow_one_bracket.2,
   alphaPow_two_bracket.1, alphaPow_two_bracket.2,
   alphaPow_three_bracket.1, alphaPow_three_bracket.2,
   alphaPow_four_bracket.1, alphaPow_four_bracket.2,
   alphaPow_five_bracket.1, alphaPow_five_bracket.2,
   alphaPow_add_1_1, alphaPow_add_2_1, alphaPow_add_2_2, alphaPow_add_3_2⟩

end PT.Holonomy
