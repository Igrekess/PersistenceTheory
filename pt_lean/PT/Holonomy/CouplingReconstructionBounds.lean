/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Holonomy.CouplingReconstruction
import Mathlib.Tactic

/-!
# Coupling reconstruction — Bounds and reciprocal (Ch09 #39 extension)

This file extends `PT.Holonomy.CouplingReconstruction` with quantitative
bracketing for the bare PT coupling `α_bare = sin²θ_3 · sin²θ_5 · sin²θ_7`
and its reciprocal `1 / α_bare` (the inverse coupling, of physical interest
since `α_EM^{-1} ≈ 137.036`).

Concretely:

* **Tight bracket (3-significant figures).** Already in
  `alphaBareQ_bracket`: `73/10000 < α_bare < 74/10000`.
* **Reciprocal upper bound.** `1 / α_bare < 137.5`.
* **Reciprocal lower bound.** `135.0 < 1 / α_bare`.
* **Strict factorisation.** The partial products
  `sin²θ_3 · sin²θ_5` and `sin²θ_3 · sin²θ_5 · sin²θ_7` exhibit strict
  decrease (each factor is in `(0, 1)`).
* **Symmetric form.** `α_bare = (sin θ_3 · sin θ_5 · sin θ_7)²` is the
  algebraic square of the *amplitude* `A = sin θ_3 · sin θ_5 · sin θ_7`.

All bounds are purely rational + `norm_num`.

## Reference

Monograph Ch09 §"Couplage et constante de structure fine", follow-up to
`alphaBareQ_bracket`. Audit row #39.
-/

namespace PT.Holonomy

/-! ### Exact value of `α_bare` at `μ* = 15`, `q = 13/15`

At the PT fixed point, `α_bare` is an **exact rational** computable by
`decide` / `native_decide`. -/

/-- **Exact value.** `α_bare = 15512777115364308026953701325116440576
    / 2114055428042547055520117282867431640625`. -/
theorem alphaBareQ_exact :
    alphaBareQ = 15512777115364308026953701325116440576
                  / 2114055428042547055520117282867431640625 := by
  unfold alphaBareQ sinSqQ deltaQ qPT
  norm_num

/-! ### Tighter bracket on `α_bare` -/

/-- **Refined bracket.** `7335/1000000 < α_bare < 7340/1000000`
    (the exact value is `≈ 7.3379 × 10⁻³`). -/
theorem alphaBareQ_bracket_tight :
    7335 / 1000000 < alphaBareQ ∧ alphaBareQ < 7340 / 1000000 := by
  unfold alphaBareQ sinSqQ deltaQ qPT
  refine ⟨?_, ?_⟩ <;> norm_num

/-! ### Reciprocal `1 / α_bare` -/

/-- The reciprocal of the bare coupling. -/
noncomputable def alphaBareInvQ : ℚ := 1 / alphaBareQ

/-- The reciprocal is positive. -/
theorem alphaBareInvQ_pos : 0 < alphaBareInvQ := by
  unfold alphaBareInvQ
  exact one_div_pos.mpr alphaBareQ_pos

/-- **Lower bound.** `135 < 1 / α_bare`. -/
theorem alphaBareInvQ_gt_135 : 135 < alphaBareInvQ := by
  unfold alphaBareInvQ alphaBareQ sinSqQ deltaQ qPT
  norm_num

/-- **Upper bound.** `1 / α_bare < 137`.
    Tightness: the published PT prediction is `α_bare^{-1} ≈ 136.20`,
    very close to `α_EM^{-1} ≈ 137.036` (gap closed by Koide-running). -/
theorem alphaBareInvQ_lt_137 : alphaBareInvQ < 137 := by
  unfold alphaBareInvQ alphaBareQ sinSqQ deltaQ qPT
  norm_num

/-- **Bracket on the reciprocal.** `135 < α_bare^{-1} < 137`. -/
theorem alphaBareInvQ_bracket :
    135 < alphaBareInvQ ∧ alphaBareInvQ < 137 :=
  ⟨alphaBareInvQ_gt_135, alphaBareInvQ_lt_137⟩

/-! ### Strict decreasing partial products -/

/-- **Partial product 2.** `0 < sin²θ_3 · sin²θ_5 < sin²θ_3 < 1`. -/
theorem partialProduct_2_decreasing :
    0 < sinSqQ 3 * sinSqQ 5
    ∧ sinSqQ 3 * sinSqQ 5 < sinSqQ 3
    ∧ sinSqQ 3 < 1 := by
  refine ⟨mul_pos sinSq_3_pos sinSq_5_pos, ?_, sinSq_3_lt_one⟩
  -- sinSq 3 * sinSq 5 < sinSq 3 ⇔ sinSq 5 < 1 (since sinSq 3 > 0)
  have h5lt : sinSqQ 5 < 1 := sinSq_5_lt_one
  nlinarith [sinSq_3_pos]

/-- **Partial product 3.** `0 < α_bare < sin²θ_3 · sin²θ_5`. -/
theorem partialProduct_3_decreasing :
    0 < alphaBareQ
    ∧ alphaBareQ < sinSqQ 3 * sinSqQ 5 := by
  refine ⟨alphaBareQ_pos, ?_⟩
  unfold alphaBareQ
  -- (sinSq 3 * sinSq 5) * sinSq 7 < (sinSq 3 * sinSq 5)  ⇔  sinSq 7 < 1
  have h7lt : sinSqQ 7 < 1 := sinSq_7_lt_one
  have h35pos : 0 < sinSqQ 3 * sinSqQ 5 := mul_pos sinSq_3_pos sinSq_5_pos
  nlinarith

/-! ### Symmetric / amplitude form -/

/-- **Amplitude form (real).** `α_bare = (sin θ_3 · sin θ_5 · sin θ_7)²`
    whenever the phase prescription `cos θ_p = 1 - δ_p` holds. -/
theorem alphaBare_eq_amplitude_squared
    (θ3 θ5 θ7 : ℝ)
    (h3 : Real.cos θ3 = 1 - (deltaQ 3 : ℝ))
    (h5 : Real.cos θ5 = 1 - (deltaQ 5 : ℝ))
    (h7 : Real.cos θ7 = 1 - (deltaQ 7 : ℝ)) :
    (Real.sin θ3 * Real.sin θ5 * Real.sin θ7) ^ 2 = (alphaBareQ : ℝ) := by
  have hprod := coupling_reconstruction_real θ3 θ5 θ7 h3 h5 h7
  -- (a·b·c)² = a²·b²·c²
  have : (Real.sin θ3 * Real.sin θ5 * Real.sin θ7) ^ 2
       = Real.sin θ3 ^ 2 * Real.sin θ5 ^ 2 * Real.sin θ7 ^ 2 := by ring
  rw [this, hprod]

/-! ### Headline: complete bracket -/

/-- **Headline.** All bounds at a glance for the bare coupling and its
    reciprocal:

    * `7335/10⁶ < α_bare < 7340/10⁶` (centered on `≈ 7.3379 · 10⁻³`),
    * `135 < α_bare⁻¹ < 137`,
    * `α_bare = sin²θ_3 · sin²θ_5 · sin²θ_7` (definition),
    * `α_bare > 0` strictly, `α_bare < 1` strictly. -/
theorem alphaBare_complete_summary :
    7335 / 1000000 < alphaBareQ ∧ alphaBareQ < 7340 / 1000000
    ∧ 135 < alphaBareInvQ ∧ alphaBareInvQ < 137
    ∧ 0 < alphaBareQ ∧ alphaBareQ < 1 := by
  refine ⟨alphaBareQ_bracket_tight.1, alphaBareQ_bracket_tight.2,
          alphaBareInvQ_gt_135, alphaBareInvQ_lt_137,
          alphaBareQ_pos, ?_⟩
  have := alphaBareQ_lt_one_hundredth
  linarith

end PT.Holonomy
