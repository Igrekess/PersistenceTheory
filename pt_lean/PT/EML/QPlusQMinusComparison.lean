/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.EML.QSheffer
import PT.EML.QParameterMonotonicity
import Mathlib.Analysis.SpecialFunctions.Exp
import Mathlib.Tactic

/-!
# Comparison of `q_+(μ)` and `q_-(μ)` on `(0, ∞)` (App O addendum)

This file isolates the **global comparison** of the two PT branch parameters
on the positive reals:

  `q_+(μ) := 1 - 2/μ`     (statistical / geometric)
  `q_-(μ) := exp(-1/μ)`   (thermal / Gibbs)

The pivotal observation is the convex inequality `exp(x) ≥ 1 + x` (Mathlib's
`Real.add_one_le_exp`). Specialised at `x = -1/μ` with `μ > 0`, it yields

  `q_-(μ) = exp(-1/μ) ≥ 1 - 1/μ > 1 - 2/μ = q_+(μ)`,

so **`q_+(μ) < q_-(μ)` everywhere on `(0, ∞)`** — the two branches do *not*
cross. Furthermore the gap is quantitatively bounded below:

  `q_-(μ) - q_+(μ) ≥ 1/μ`  for every `μ > 0`.

Both parameters share the same limit `1` as `μ → ∞`, but `q_-` approaches `1`
strictly faster (the gap closes like `1/μ`).

## Reference

Appendix O of the monograph, follow-up to `QSheffer.lean` and
`QParameterMonotonicity.lean`. Uses the elementary convex bound
`exp(x) ≥ 1 + x` (Mathlib `Real.add_one_le_exp`).
-/

namespace PT.EML

open Real

/-! ### Quantitative gap `q_- - q_+ ≥ 1/μ` -/

/-- **Quantitative gap.** For every `μ > 0`,
    `q_-(μ) - q_+(μ) ≥ 1/μ`.

    Proof: `q_-(μ) = exp(-1/μ) ≥ 1 + (-1/μ) = 1 - 1/μ` by the convex bound
    `exp(x) ≥ 1 + x`. Subtracting `q_+(μ) = 1 - 2/μ` gives
    `q_-(μ) - q_+(μ) ≥ (1 - 1/μ) - (1 - 2/μ) = 1/μ`. -/
theorem qMinus_sub_qPlus_ge_inv (μ : ℝ) (_hμ : 0 < μ) :
    1 / μ ≤ qMinus μ - qPlus μ := by
  unfold qMinus qPlus
  -- Step 1: exp(-1/μ) ≥ 1 - 1/μ (Mathlib gives `x + 1 ≤ exp x`)
  have h_exp : -(1 / μ) + 1 ≤ Real.exp (-(1 / μ)) :=
    Real.add_one_le_exp (-(1 / μ))
  -- Tie `2/μ` and `1/μ` together so `linarith` sees one atom.
  have h2μ : (2 : ℝ) / μ = 2 * (1 / μ) := by ring
  rw [h2μ]
  linarith

/-! ### Global strict ordering `q_+ < q_-` on `(0, ∞)` -/

/-- **Main theorem — no crossing.** For every `μ > 0`,
    `q_+(μ) < q_-(μ)`.

    Proof: the quantitative gap `q_-(μ) - q_+(μ) ≥ 1/μ > 0`. -/
theorem qPlus_lt_qMinus (μ : ℝ) (hμ : 0 < μ) : qPlus μ < qMinus μ := by
  have hgap : 1 / μ ≤ qMinus μ - qPlus μ := qMinus_sub_qPlus_ge_inv μ hμ
  have hpos : 0 < 1 / μ := by positivity
  linarith

/-- **No crossing (negation form).** The equation `q_+(μ) = q_-(μ)` has no
    solution in `(0, ∞)`. -/
theorem qPlus_ne_qMinus (μ : ℝ) (hμ : 0 < μ) : qPlus μ ≠ qMinus μ := by
  intro h
  have hlt : qPlus μ < qMinus μ := qPlus_lt_qMinus μ hμ
  exact (lt_irrefl _ (h ▸ hlt))

/-- **No crossing — universal form.** There is no `μ > 0` such that
    `q_+(μ) = q_-(μ)`. -/
theorem no_crossing_on_pos :
    ¬ ∃ μ : ℝ, 0 < μ ∧ qPlus μ = qMinus μ := by
  rintro ⟨μ, hμ, heq⟩
  exact qPlus_ne_qMinus μ hμ heq

/-! ### Specialisation at the PT fixed point `μ = 15` -/

/-- **Gap at `μ = 15`.** `q_-(15) - q_+(15) ≥ 1/15`. -/
theorem qMinus_sub_qPlus_at_15_ge : 1 / 15 ≤ qMinus 15 - qPlus 15 :=
  qMinus_sub_qPlus_ge_inv 15 (by norm_num)

/-- **Explicit bracket at `μ = 15`.**
    `q_+(15) = 13/15 < 14/15 ≤ q_-(15)`. -/
theorem qPlus_lt_14_15_le_qMinus_at_15 :
    qPlus 15 < 14 / 15 ∧ 14 / 15 ≤ qMinus 15 := by
  refine ⟨?_, ?_⟩
  · -- q_+(15) = 13/15 < 14/15
    unfold qPlus; norm_num
  · -- q_-(15) = exp(-1/15) ≥ 1 + (-1/15) = 14/15
    unfold qMinus
    have h : -((1 : ℝ) / 15) + 1 ≤ Real.exp (-((1 : ℝ) / 15)) :=
      Real.add_one_le_exp (-((1 : ℝ) / 15))
    have habs : -((1 : ℝ) / 15) + 1 = 14 / 15 := by norm_num
    linarith

/-! ### Limit behaviour as `μ → ∞`

    Both `q_+(μ)` and `q_-(μ)` approach `1` from below. We record this in
    elementary algebraic form (no `Filter.Tendsto`): for every `ε > 0`
    there is `N > 0` such that for all `μ > N`, `|q_±(μ) - 1| < ε`. -/

/-- **`q_+ → 1` (elementary form).** For every `ε > 0` there is `N > 0`
    such that for all `μ > N`, `1 - ε < qPlus μ < 1`. -/
theorem qPlus_tendsto_one_elem :
    ∀ ε : ℝ, 0 < ε → ∃ N : ℝ, 0 < N ∧
      ∀ μ : ℝ, N < μ → 1 - ε < qPlus μ ∧ qPlus μ < 1 := by
  intro ε hε
  refine ⟨max 1 (2 / ε), ?_, ?_⟩
  · exact lt_max_of_lt_left (by norm_num : (0 : ℝ) < 1)
  intro μ hμN
  have hμ_pos : 0 < μ :=
    lt_of_lt_of_le (by norm_num : (0 : ℝ) < 1)
      (le_trans (le_max_left (1 : ℝ) (2 / ε)) (le_of_lt hμN))
  have hμ_2 : 2 / ε < μ :=
    lt_of_le_of_lt (le_max_right (1 : ℝ) (2 / ε)) hμN
  refine ⟨?_, qPlus_lt_one μ hμ_pos⟩
  -- Need 1 - ε < 1 - 2/μ  ⇔  2/μ < ε.
  unfold qPlus
  -- From 2/ε < μ and ε > 0, μ > 0: 2/μ < ε.
  have h2μ : 2 / μ < ε := by
    rw [div_lt_iff₀ hμ_pos]
    have : 2 / ε * ε < μ * ε := by
      exact (mul_lt_mul_iff_of_pos_right hε).mpr hμ_2
    have hε_ne : ε ≠ 0 := ne_of_gt hε
    rw [div_mul_cancel₀ 2 hε_ne] at this
    linarith
  linarith

/-- **`q_- → 1` (elementary form).** For every `ε > 0` there is `N > 0`
    such that for all `μ > N`, `1 - ε < qMinus μ < 1`. -/
theorem qMinus_tendsto_one_elem :
    ∀ ε : ℝ, 0 < ε → ∃ N : ℝ, 0 < N ∧
      ∀ μ : ℝ, N < μ → 1 - ε < qMinus μ ∧ qMinus μ < 1 := by
  intro ε hε
  -- Use the same `N` as for `q_+`: by `q_+ < q_-`, the lower bound
  -- `1 - ε < q_+` already implies `1 - ε < q_-`.
  obtain ⟨N, hN_pos, hN_qp⟩ := qPlus_tendsto_one_elem ε hε
  refine ⟨N, hN_pos, ?_⟩
  intro μ hμN
  have hμ_pos : 0 < μ := lt_trans hN_pos hμN
  have hqp := hN_qp μ hμN
  have hpqm : qPlus μ < qMinus μ := qPlus_lt_qMinus μ hμ_pos
  refine ⟨lt_trans hqp.1 hpqm, qMinus_lt_one μ hμ_pos⟩

/-! ### Headline summary -/

/-- **Headline (`q_+` vs `q_-` global comparison).** For every `μ > 0`:

    * **No crossing**: `q_+(μ) < q_-(μ)`.
    * **Quantitative gap**: `q_-(μ) - q_+(μ) ≥ 1/μ`.
    * **Specialisation at `μ = 15`**: `q_+(15) < 14/15 ≤ q_-(15)`.
    * **Common limit**: both `q_+(μ)` and `q_-(μ)` tend to `1` as `μ → ∞`
      (elementary `ε`–`N` form). -/
theorem qPlus_qMinus_comparison_summary :
    (∀ μ : ℝ, 0 < μ → qPlus μ < qMinus μ)
    ∧ (∀ μ : ℝ, 0 < μ → 1 / μ ≤ qMinus μ - qPlus μ)
    ∧ (qPlus 15 < 14 / 15 ∧ 14 / 15 ≤ qMinus 15)
    ∧ (¬ ∃ μ : ℝ, 0 < μ ∧ qPlus μ = qMinus μ) :=
  ⟨qPlus_lt_qMinus, qMinus_sub_qPlus_ge_inv,
   qPlus_lt_14_15_le_qMinus_at_15, no_crossing_on_pos⟩

end PT.EML
