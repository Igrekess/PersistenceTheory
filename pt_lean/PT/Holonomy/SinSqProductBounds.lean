/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Holonomy.CouplingReconstruction
import PT.Holonomy.CouplingReconstructionBounds
import Mathlib.Tactic

/-!
# `∏ sin²θ_p` — Partial product bounds on the active cascade (Ch09 extension)

This file gives **partial-product brackets** for the active cascade
`{p_1, p_2, p_3} = {3, 5, 7}`. Each partial product is bounded between two
explicit rationals.

* `P_1 := sin²θ_3` ∈ `(0.219, 0.220)` (already in `CouplingReconstruction`).
* `P_2 := sin²θ_3 · sin²θ_5` ∈ `(0.042, 0.043)`.
* `P_3 := sin²θ_3 · sin²θ_5 · sin²θ_7 = α_bare ≈ 7.34 × 10⁻³`.

Plus comparisons: `P_3 < P_2 < P_1 < 1`.

## Reference

Monograph Chapter 9 §"Produit partiel des phases cycliques".
-/

namespace PT.Holonomy

/-! ### Partial-product definitions -/

/-- The single-factor partial product `P_1 = sin²θ_3`. -/
def P1 : ℚ := sinSqQ 3

/-- The two-factor partial product `P_2 = sin²θ_3 · sin²θ_5`. -/
def P2 : ℚ := sinSqQ 3 * sinSqQ 5

/-- The three-factor partial product `P_3 = sin²θ_3 · sin²θ_5 · sin²θ_7
    = α_bare`. -/
def P3 : ℚ := sinSqQ 3 * sinSqQ 5 * sinSqQ 7

@[simp] theorem P3_eq_alphaBareQ : P3 = alphaBareQ := rfl

/-! ### Positivity -/

theorem P1_pos : 0 < P1 := sinSq_3_pos
theorem P2_pos : 0 < P2 := mul_pos sinSq_3_pos sinSq_5_pos
theorem P3_pos : 0 < P3 := alphaBareQ_pos

/-! ### Brackets -/

/-- `0.219 < P_1 < 0.220`. -/
theorem P1_bracket : 219 / 1000 < P1 ∧ P1 < 220 / 1000 := sinSq_3_bracket

/-- `0.042 < P_2 < 0.043`. -/
theorem P2_bracket : 42 / 1000 < P2 ∧ P2 < 43 / 1000 := by
  unfold P2 sinSqQ deltaQ qPT
  refine ⟨?_, ?_⟩ <;> norm_num

/-- `7335/10⁶ < P_3 < 7340/10⁶` (already in `alphaBareQ_bracket_tight`). -/
theorem P3_bracket : 7335 / 1000000 < P3 ∧ P3 < 7340 / 1000000 := by
  rw [P3_eq_alphaBareQ]; exact alphaBareQ_bracket_tight

/-! ### Strict decreasing chain `P_3 < P_2 < P_1 < 1` -/

/-- `P_3 < P_2` (multiplying by `sin²θ_7 < 1` strictly shrinks). -/
theorem P3_lt_P2 : P3 < P2 := by
  unfold P3 P2
  have h7_lt : sinSqQ 7 < 1 := sinSq_7_lt_one
  have h35_pos : 0 < sinSqQ 3 * sinSqQ 5 := mul_pos sinSq_3_pos sinSq_5_pos
  nlinarith

/-- `P_2 < P_1` (multiplying by `sin²θ_5 < 1` strictly shrinks). -/
theorem P2_lt_P1 : P2 < P1 := by
  unfold P2 P1
  have h5_lt : sinSqQ 5 < 1 := sinSq_5_lt_one
  nlinarith [sinSq_3_pos]

/-- `P_1 < 1`. -/
theorem P1_lt_one : P1 < 1 := sinSq_3_lt_one

/-- **Strict chain.** `0 < P_3 < P_2 < P_1 < 1`. -/
theorem partial_product_chain :
    0 < P3 ∧ P3 < P2 ∧ P2 < P1 ∧ P1 < 1 :=
  ⟨P3_pos, P3_lt_P2, P2_lt_P1, P1_lt_one⟩

/-! ### Numerical headline -/

/-- **Headline (partial products of the active cascade).** All four
    quantities at the PT fixed point `q = 13/15`:

    * `P_1 = sin²θ_3 ∈ (0.219, 0.220)`,
    * `P_2 = P_1 · sin²θ_5 ∈ (0.042, 0.043)`,
    * `P_3 = P_2 · sin²θ_7 = α_bare ∈ (7335/10⁶, 7340/10⁶) ≈ 7.34 · 10⁻³`,
    * strictly decreasing: `P_3 < P_2 < P_1 < 1`. -/
theorem partial_product_summary :
    -- brackets
    219 / 1000 < P1 ∧ P1 < 220 / 1000
    ∧ 42 / 1000 < P2 ∧ P2 < 43 / 1000
    ∧ 7335 / 1000000 < P3 ∧ P3 < 7340 / 1000000
    -- chain
    ∧ 0 < P3 ∧ P3 < P2 ∧ P2 < P1 ∧ P1 < 1
    -- P_3 is α_bare
    ∧ P3 = alphaBareQ :=
  ⟨P1_bracket.1, P1_bracket.2,
   P2_bracket.1, P2_bracket.2,
   P3_bracket.1, P3_bracket.2,
   P3_pos, P3_lt_P2, P2_lt_P1, P1_lt_one,
   P3_eq_alphaBareQ⟩

end PT.Holonomy
