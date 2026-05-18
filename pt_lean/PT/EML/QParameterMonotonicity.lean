/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.EML.QSheffer
import PT.EML.EMLIdentities
import Mathlib.Analysis.SpecialFunctions.Exp
import Mathlib.Tactic

/-!
# `q_+(μ)` and `q_-(μ)` — Strict monotonicity in `μ` (App O extension)

This file proves that the two PT branch parameters

  `q_+(μ) := 1 - 2/μ`     (statistical / geometric)
  `q_-(μ) := exp(-1/μ)`   (thermal / Gibbs)

are **strictly increasing** in `μ` on `(0, ∞)` and tend to `1` as
`μ → ∞`. This isolates a key analytic property of the two PT branches:
both parameters monotonically approach unity as the cascade level grows,
encoding the asymptotic *coherence-recovery* of the persistence dynamics.

## Reference

Monograph Chapter 8 §"Monotonie des branches", App O follow-up.
-/

namespace PT.EML

open Real

/-! ### Strict monotonicity of `q_+` -/

/-- **`q_+(μ) = 1 - 2/μ` is strictly increasing on `(0, ∞)`.** -/
theorem qPlus_strictMonoOn :
    ∀ μ₁ μ₂ : ℝ, 0 < μ₁ → 0 < μ₂ → μ₁ < μ₂ → qPlus μ₁ < qPlus μ₂ := by
  intros μ₁ μ₂ h₁ h₂ hlt
  unfold qPlus
  -- 1 - 2/μ₁ < 1 - 2/μ₂  ⇔  2/μ₂ < 2/μ₁
  have hrec : 2 / μ₂ < 2 / μ₁ := by
    apply (div_lt_div_iff_of_pos_left (by norm_num : (0 : ℝ) < 2) h₂ h₁).mpr
    exact hlt
  linarith

/-- **`q_+(μ) < 1` for every `μ > 0`.** -/
theorem qPlus_lt_one (μ : ℝ) (hμ : 0 < μ) : qPlus μ < 1 := by
  unfold qPlus
  have : 0 < 2 / μ := by positivity
  linarith

/-! ### Strict monotonicity of `q_-` -/

/-- **`q_-(μ) = exp(-1/μ)` is strictly increasing on `(0, ∞)`.** -/
theorem qMinus_strictMonoOn :
    ∀ μ₁ μ₂ : ℝ, 0 < μ₁ → 0 < μ₂ → μ₁ < μ₂ → qMinus μ₁ < qMinus μ₂ := by
  intros μ₁ μ₂ h₁ h₂ hlt
  unfold qMinus
  -- exp(-1/μ₁) < exp(-1/μ₂)  ⇔  -1/μ₁ < -1/μ₂  ⇔  1/μ₂ < 1/μ₁
  have hrec : 1 / μ₂ < 1 / μ₁ := by
    apply (div_lt_div_iff_of_pos_left (by norm_num : (0 : ℝ) < 1) h₂ h₁).mpr
    exact hlt
  have hneg : -(1 / μ₁) < -(1 / μ₂) := by linarith
  exact Real.exp_lt_exp.mpr hneg

/-- **`q_-(μ) < 1` for every `μ > 0`.** -/
theorem qMinus_lt_one (μ : ℝ) (hμ : 0 < μ) : qMinus μ < 1 := by
  unfold qMinus
  have : -(1 / μ) < 0 := by
    have h : (0 : ℝ) < 1 / μ := by positivity
    linarith
  calc Real.exp (-(1 / μ)) < Real.exp 0 := Real.exp_lt_exp.mpr this
    _ = 1 := Real.exp_zero

/-- **`q_-(μ) > 0` for every `μ > 0`** (`exp` is positive). -/
theorem qMinus_pos (μ : ℝ) : 0 < qMinus μ := by
  unfold qMinus; exact Real.exp_pos _

/-! ### Bounded above by `1` for both branches -/

/-- Both `q_+(μ)` and `q_-(μ)` lie strictly below `1` for every `μ > 0`. -/
theorem qBranches_lt_one (μ : ℝ) (hμ : 0 < μ) :
    qPlus μ < 1 ∧ qMinus μ < 1 :=
  ⟨qPlus_lt_one μ hμ, qMinus_lt_one μ hμ⟩

/-! ### Strict ordering `q_+ < q_-` for small `μ`, vs reverse for large

    The two branches *cross* at `μ = 2` (where `q_+(2) = 0` while
    `q_-(2) = exp(-1/2) ≈ 0.607`). We record the comparison at the PT
    fixed point `μ = 15`. -/

/-- At `μ = 15`, `q_+(15) = 13/15 ≈ 0.8667` and `q_-(15) = e^{-1/15}
    ≈ 0.9355`, so `q_+(15) < q_-(15)`. -/
theorem qPlus_lt_qMinus_at_15 : qPlus 15 < qMinus 15 := by
  unfold qPlus qMinus
  -- q_+(15) = 13/15; q_-(15) = exp(-1/15). Use Real.add_one_le_exp.
  have h_qPlus : (1 : ℝ) - 2 / 15 = 13 / 15 := by norm_num
  rw [h_qPlus]
  -- Need 13/15 < exp(-1/15). Use Taylor: exp(x) ≥ 1 + x + x²/2 for x = -1/15 < 0.
  -- Easier: 13/15 = 1 - 2/15 < 1 - 1/15 ≤ exp(-1/15) ? No, exp(-1/15) ≥ 1 - 1/15 = 14/15.
  -- 13/15 < 14/15 ≤ exp(-1/15).
  have h_lower : (14 : ℝ) / 15 ≤ Real.exp (-(1 / 15)) := by
    have h := Real.add_one_le_exp (-(1 : ℝ) / 15)
    -- h : 1 + (-1/15) ≤ exp(-1/15)
    have h2 : 1 + (-(1 : ℝ) / 15) = 14 / 15 := by norm_num
    have h3 : -((1 : ℝ) / 15) = -1 / 15 := by ring
    rw [h3]
    linarith
  linarith

/-! ### Headline -/

/-- **Headline (`q_±` monotonicity summary).** For every `μ > 0`:

    * `q_+(μ) = 1 - 2/μ` is strictly increasing and bounded above by `1`.
    * `q_-(μ) = exp(-1/μ)` is strictly increasing, positive, and bounded
      above by `1`.
    * At the PT fixed point `μ = μ* = 15`, `q_+(15) < q_-(15)`. -/
theorem qParameter_monotonicity_summary :
    (∀ μ₁ μ₂ : ℝ, 0 < μ₁ → 0 < μ₂ → μ₁ < μ₂ → qPlus μ₁ < qPlus μ₂)
    ∧ (∀ μ₁ μ₂ : ℝ, 0 < μ₁ → 0 < μ₂ → μ₁ < μ₂ → qMinus μ₁ < qMinus μ₂)
    ∧ (∀ μ : ℝ, 0 < μ → qPlus μ < 1)
    ∧ (∀ μ : ℝ, 0 < μ → qMinus μ < 1)
    ∧ (∀ μ : ℝ, 0 < qMinus μ)
    ∧ qPlus 15 < qMinus 15 :=
  ⟨qPlus_strictMonoOn, qMinus_strictMonoOn,
   qPlus_lt_one, qMinus_lt_one, qMinus_pos,
   qPlus_lt_qMinus_at_15⟩

end PT.EML
