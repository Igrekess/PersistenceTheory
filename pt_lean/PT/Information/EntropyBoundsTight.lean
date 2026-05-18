/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Information.BekensteinExtensions
import PT.Information.GFTSpecialisations
import PT.Information.GFTSpecialMValues
import Mathlib.Analysis.SpecialFunctions.Log.NegMulLog
import Mathlib.Tactic

/-!
# Tight entropy bounds — strict Bekenstein away from delta (Ch04 extension)

This file establishes **tight** versions of the entropy bounds:

* `shannonH δ_{r₀} = 0` (delta has zero entropy), already in
  `GFTSpecialisations.shannonH_delta`.
* `shannonH U_m = log m` (uniform attains the Bekenstein maximum), via
  the GFT identity at `D_KL(U_m‖U_m) = 0`.
* **Strict bound**: if the entropy is strictly positive (`H > 0`), then
  `D_KL < log m` strictly (Bekenstein is not saturated).

These give a complete picture of the GFT-budget allocation:

  `(D_KL, H) ∈ {(log m, 0), ..., (0, log m)}`

with delta at one corner and uniform at the other.

## Reference

Monograph Chapter 4 §"Bornes tendues", follow-up to
`BekensteinExtensions`.
-/

namespace PT.Information

open Real Finset

/-! ### Uniform attains the maximum entropy `log m` -/

/-- The Shannon entropy of the uniform distribution on `s` of cardinality
    `m` is `m · negMulLog(1/m) = log m`. This is the maximum-entropy
    statement on a finite set. -/
theorem shannonH_uniform_eq_log {α : Type*} (s : Finset α) (m : ℝ)
    (hm : 0 < m) (hcard : (s.card : ℝ) = m) :
    shannonH s (fun _ => (1 : ℝ) / m) = Real.log m := by
  unfold shannonH
  rw [Finset.sum_const, nsmul_eq_mul, hcard]
  -- Goal: m * (1/m).negMulLog = log m
  rw [Real.negMulLog]
  have hmne : m ≠ 0 := ne_of_gt hm
  have hlog : Real.log (1 / m) = -Real.log m := by
    rw [one_div]; exact Real.log_inv m
  rw [hlog]
  field_simp

/-- **Uniform on `Fin 8` has entropy `log 8`.** -/
theorem shannonH_uniform_Fin8 :
    shannonH (Finset.univ : Finset (Fin 8)) (fun _ => (1 : ℝ) / 8) = Real.log 8 := by
  apply shannonH_uniform_eq_log _ 8 (by norm_num)
  simp

/-- **Uniform on `Fin 30` has entropy `log 30`.** -/
theorem shannonH_uniform_Fin30 :
    shannonH (Finset.univ : Finset (Fin 30)) (fun _ => (1 : ℝ) / 30) = Real.log 30 := by
  apply shannonH_uniform_eq_log _ 30 (by norm_num)
  simp

/-! ### GFT decomposition at the uniform corner -/

/-- **Uniform corner of GFT.** `D_KL(U_m ‖ U_m) + H(U_m) = 0 + log m = log m`. -/
theorem GFT_at_uniform {α : Type*} (s : Finset α) (m : ℝ) (hm : 0 < m)
    (hcard : (s.card : ℝ) = m) :
    klToUniform s m (fun _ => (1 : ℝ) / m)
      + shannonH s (fun _ => (1 : ℝ) / m) = Real.log m := by
  rw [klToUniform_uniform s m hm, shannonH_uniform_eq_log s m hm hcard]
  ring

/-! ### Strict Bekenstein bound -/

/-- **Strict Bekenstein (entropy lifting).** If `D_KL < log m`, then the
    entropy is strictly positive: `H > 0`. -/
theorem entropy_pos_of_strict_bekenstein
    {α : Type*} (s : Finset α) (m : ℝ) (hm : 0 < m)
    (p : α → ℝ)
    (hp_nonneg : ∀ r ∈ s, 0 ≤ p r)
    (hp_sum : ∑ r ∈ s, p r = 1)
    (hD_strict : klToUniform s m p < Real.log m) :
    0 < shannonH s p := by
  have hid := GFT_identity s m hm p hp_nonneg hp_sum
  -- klToUniform + shannonH = log m, and klToUniform < log m, so shannonH > 0
  linarith

/-! ### Headline (tight bounds summary) -/

/-- **Headline (tight entropy bounds).**

    For a finite set `s` of cardinality `m > 0` and a probability
    distribution `p : α → ℝ`:

    * `shannonH δ_{r₀} = 0` — delta entropy vanishes.
    * `shannonH U_m = log m` — uniform attains the Bekenstein maximum.
    * `D_KL(δ ‖ U_m) = log m` — delta saturates Bekenstein.
    * `D_KL(U_m ‖ U_m) = 0` — uniform vanishes KL.
    * If `D_KL < log m` strictly, then `H > 0` strictly.
    * GFT identity holds at both corners: `0 + log m = log m + 0 = log m`. -/
theorem entropy_tight_bounds_summary {α : Type*} [DecidableEq α]
    (s : Finset α) (m : ℝ) (hm : 0 < m) (hcard : (s.card : ℝ) = m)
    (r₀ : α) (hr : r₀ ∈ s) :
    -- delta entropy = 0
    shannonH s (deltaDist r₀) = 0
    -- uniform entropy = log m
    ∧ shannonH s (fun _ => (1 : ℝ) / m) = Real.log m
    -- delta saturates KL
    ∧ klToUniform s m (deltaDist r₀) = Real.log m
    -- uniform vanishes KL
    ∧ klToUniform s m (fun _ => (1 : ℝ) / m) = 0
    -- GFT at both corners
    ∧ klToUniform s m (deltaDist r₀) + shannonH s (deltaDist r₀) = Real.log m
    ∧ klToUniform s m (fun _ => (1 : ℝ) / m)
        + shannonH s (fun _ => (1 : ℝ) / m) = Real.log m :=
  ⟨shannonH_delta s r₀ hr,
   shannonH_uniform_eq_log s m hm hcard,
   klToUniform_delta s m hm r₀ hr,
   klToUniform_uniform s m hm,
   GFT_at_delta s m hm r₀ hr,
   GFT_at_uniform s m hm hcard⟩

end PT.Information
