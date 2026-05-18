/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.FixedPoint.T7MuStar
import Mathlib.Tactic

/-!
# T7 — Global uniqueness of `μ* = 15` (Appendix, item #34)

`PT.FixedPoint.T7MuStar` proves combinatorially that `μ* = 15` is the unique
fixed point of the persistence-active map `F_pers` on the finite range
`μ ∈ [8, 20]`. This file lifts that result to **global uniqueness on
`μ ≥ 8`**: for every integer `μ ≥ 8`, `F_pers μ = μ` implies `μ = 15`.

The strategy is elementary:

* For `μ ≥ 7`, every conditional inside `F_pers` triggers, so
  `F_pers μ = 3 + 5 + 7 = 15` is **constant**.
* Hence for `μ ≥ 16` (i.e. strictly above `μ*`), `F_pers μ = 15 < μ`, ruling
  out any fixed point.
* Combined with the `[8, 20]` enumeration from `T7MuStar`, this closes the
  uniqueness statement over the full half-line `μ ≥ 8`.

This is the form expected by item #34 of `AUDIT_FORMALISABLE.md` ("Unique
fixed point — μ* = 15 globally unique; no alternatives"). The complementary
restriction `μ ≥ 8` (excluding the degenerate one-prime fixed point `μ = 3`)
is exactly the PT side condition `|P_active| ≥ 2`.

## Reference

Item #34 of `AUDIT_FORMALISABLE.md`; monograph Chapter 7, headline
`unique_fixedpoint`.
-/

namespace PT.FixedPoint

/-! ### `F_pers` is constant on `μ ≥ 7` -/

/-- For every `μ ≥ 7`, `F_pers μ = 15`. All three thresholds (`3, 5, 7`)
    are active. -/
theorem Fpers_eq_fifteen_of_ge_seven (μ : ℕ) (h : 7 ≤ μ) : Fpers μ = 15 := by
  unfold Fpers
  have h3 : 3 ≤ μ := by omega
  have h5 : 5 ≤ μ := by omega
  rw [if_pos h3, if_pos h5, if_pos h]

/-- For every `μ ≥ 16`, `F_pers μ = 15 ≠ μ`. -/
theorem Fpers_ne_self_of_ge_sixteen (μ : ℕ) (h : 16 ≤ μ) : Fpers μ ≠ μ := by
  rw [Fpers_eq_fifteen_of_ge_seven μ (by omega)]
  omega

/-! ### Global uniqueness on `μ ≥ 8` -/

/-- **T7 global uniqueness (item #34).** For every integer `μ ≥ 8`,
    `F_pers μ = μ` implies `μ = 15`. Equivalently: `μ* = 15` is the unique
    persistence-active fixed point on the non-degenerate range. -/
theorem T7_muStar_unique_global (μ : ℕ) (h : 8 ≤ μ) (hfix : Fpers μ = μ) :
    μ = 15 := by
  by_cases h20 : μ ≤ 20
  · exact (T7_unique_fixedpoint_small μ h h20).mp hfix
  · have h16 : 16 ≤ μ := by omega
    exact absurd hfix (Fpers_ne_self_of_ge_sixteen μ h16)

/-- **T7 global uniqueness — `μ*` form.** Direct statement in terms of
    `muStar`. -/
theorem T7_muStar_unique_global_muStar (μ : ℕ) (h : 8 ≤ μ) (hfix : Fpers μ = μ) :
    μ = muStar := by
  rw [muStar_eq]
  exact T7_muStar_unique_global μ h hfix

/-- **T7 global existence + uniqueness, iff form.**
    For `μ ≥ 8`, `μ` is a fixed point of `F_pers` iff `μ = 15`. -/
theorem T7_fixedpoint_iff (μ : ℕ) (h : 8 ≤ μ) : Fpers μ = μ ↔ μ = 15 := by
  refine ⟨T7_muStar_unique_global μ h, ?_⟩
  rintro rfl
  exact Fpers_at_15

/-! ### Bounds on `F_pers` outside the fixed-point regime -/

/-- For `μ ∈ [8, 14]` (between the two-prime activation threshold and `μ*`),
    `F_pers μ` overshoots `μ`. Specifically `F_pers μ = 15 > μ`. -/
theorem Fpers_gt_self_on_8_to_14 (μ : ℕ) (h₁ : 8 ≤ μ) (h₂ : μ ≤ 14) :
    μ < Fpers μ := by
  rw [Fpers_eq_fifteen_of_ge_seven μ (by omega)]
  omega

/-- For `μ ≥ 16`, `F_pers μ` undershoots `μ`. -/
theorem Fpers_lt_self_of_ge_sixteen (μ : ℕ) (h : 16 ≤ μ) : Fpers μ < μ := by
  rw [Fpers_eq_fifteen_of_ge_seven μ (by omega)]
  omega

/-- The fixed-point equation `F_pers μ = μ` has **exactly two solutions** on
    `μ ≥ 2`: the trivial one-prime cascade `μ = 3` (where only `p = 3` is
    active and `3 + 0 + 0 = 3`) and the non-trivial PT fixed point `μ = 15`. -/
theorem T7_all_fixed_points (μ : ℕ) (h : 2 ≤ μ) (hfix : Fpers μ = μ) :
    μ = 3 ∨ μ = 15 := by
  by_cases h8 : 8 ≤ μ
  · exact Or.inr (T7_muStar_unique_global μ h8 hfix)
  · -- μ ∈ [2, 7]; enumerate
    interval_cases μ
    all_goals first | (left; rfl) | (exfalso; revert hfix; unfold Fpers; decide)

end PT.FixedPoint
