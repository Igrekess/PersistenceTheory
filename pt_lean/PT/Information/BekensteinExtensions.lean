/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Information.BekensteinBound
import PT.Information.GFTSpecialisations
import Mathlib.Tactic

/-!
# Bekenstein Bound — Saturation and corner equalities (Ch04 #18 extension)

This file extends `PT.Information.BekensteinBound` with the **saturation
analysis**:

* **Delta saturates the bound.** When `P = δ_{r₀}`, the Bekenstein bound
  is achieved with equality: `D_KL(δ_{r₀} ‖ U_m) = log m`.
* **Uniform attains the lower bound 0.** When `P = U_m` (uniform on `s`
  of cardinality `m`), `D_KL = 0` (trivially) and `H = log m` (entropy
  saturated).
* **Strict inequality away from delta.** If `P` is strictly between
  delta and uniform (every coordinate strictly less than 1), then
  `D_KL < log m` strictly (since `H > 0`).

These are direct algebraic corollaries of `BekensteinBound.bekenstein_bound`
and `GFTSpecialisations.GFT_at_delta`.

## Reference

Monograph Ch04 §"Saturation Bekenstein", `\label{cor:bekenstein-saturation}`.
M4 article.
-/

namespace PT.Information

open Real Finset

/-! ### Delta saturates Bekenstein -/

/-- **Delta saturates the Bekenstein bound.** `D_KL(δ_{r₀} ‖ U_m) = log m`. -/
theorem bekenstein_saturated_at_delta {α : Type*} [DecidableEq α]
    (s : Finset α) (m : ℝ) (hm : 0 < m) (r₀ : α) (hr : r₀ ∈ s) :
    klToUniform s m (deltaDist r₀) = Real.log m :=
  klToUniform_delta s m hm r₀ hr

/-! ### Uniform corner: D_KL = 0, H = log m -/

/-- The KL divergence of the uniform distribution `1/m` to itself is `0`.
    Pointwise: `(1/m) · log (m · (1/m)) = (1/m) · log 1 = 0`. -/
theorem klToUniform_uniform {α : Type*} (s : Finset α) (m : ℝ) (hm : 0 < m) :
    klToUniform s m (fun _ => (1 : ℝ) / m) = 0 := by
  unfold klToUniform
  apply Finset.sum_eq_zero
  intro r _
  have hmne : (m : ℝ) ≠ 0 := ne_of_gt hm
  -- Goal: (1/m) * log(m * (1/m)) = 0
  have hone : m * ((1 : ℝ) / m) = 1 := by field_simp
  show (1 / m) * Real.log (m * (1 / m)) = 0
  rw [hone, Real.log_one, mul_zero]

/-! ### Strict inequality away from saturation -/

/-- **Strict Bekenstein.** If the entropy is strictly positive
    (`H > 0`), then the Bekenstein inequality is strict:
    `D_KL < log m`. -/
theorem bekenstein_strict_of_entropy_pos
    {α : Type*} (s : Finset α) (m : ℝ) (hm : 0 < m)
    (p : α → ℝ)
    (hp_nonneg : ∀ r ∈ s, 0 ≤ p r)
    (hp_sum : ∑ r ∈ s, p r = 1)
    (hH_pos : 0 < shannonH s p) :
    klToUniform s m p < Real.log m := by
  have hid := GFT_identity s m hm p hp_nonneg hp_sum
  -- klToUniform + shannonH = log m, and shannonH > 0, so klToUniform < log m
  linarith

/-! ### Headline (saturation summary) -/

/-- **Headline (saturation summary, `m ≥ 1` form).** For a finite set of
    cardinality `m ≥ 1`, the Bekenstein bound `D_KL(P ‖ U_m) ≤ log m` is

    * **saturated** at the delta distributions `P = δ_{r₀}`
      (`D_KL = log m`),
    * **attained at 0** at the uniform distribution `P = U_m`
      (`D_KL = 0`),

    and is **strict** in between whenever the entropy is strictly positive.
    The `m ≥ 1` hypothesis is what guarantees `log m ≥ 0` so the comparison
    `0 ≤ log m` is non-vacuous. -/
theorem bekenstein_saturation_summary {α : Type*} [DecidableEq α]
    (s : Finset α) (m : ℝ) (hm : 1 ≤ m) (r₀ : α) (hr : r₀ ∈ s) :
    -- Saturation at delta:
    klToUniform s m (deltaDist r₀) = Real.log m
    -- Vanishing at uniform:
    ∧ klToUniform s m (fun _ => (1 : ℝ) / m) = 0
    -- Both ≤ log m (Bekenstein):
    ∧ klToUniform s m (deltaDist r₀) ≤ Real.log m
    ∧ klToUniform s m (fun _ => (1 : ℝ) / m) ≤ Real.log m := by
  have hm_pos : 0 < m := lt_of_lt_of_le zero_lt_one hm
  have hlog_nn : 0 ≤ Real.log m := Real.log_nonneg hm
  refine ⟨bekenstein_saturated_at_delta s m hm_pos r₀ hr,
          klToUniform_uniform s m hm_pos,
          ?_, ?_⟩
  · rw [bekenstein_saturated_at_delta s m hm_pos r₀ hr]
  · rw [klToUniform_uniform s m hm_pos]
    exact hlog_nn

end PT.Information
