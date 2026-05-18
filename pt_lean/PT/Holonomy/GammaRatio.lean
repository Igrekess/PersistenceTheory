/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Holonomy.ActivePrimeCriterion
import PT.Holonomy.GammaTablesExtended
import PT.Holonomy.GammaMonotonicity
import PT.Holonomy.GammaProduct
import Mathlib.Tactic

/-!
# γ_p / γ_{p'} — Ratios of anomalous dimensions across the PT cascade

**Statement (paper-level, Ch06 §"Cascade arithmétique de persistance").**
The anomalous dimensions `γ_p` of Persistence Theory at the fixed point
`μ* = 15`, `q = 13/15`, are organised in a strictly decreasing cascade
indexed by the odd primes. The *ratios* `γ_p / γ_{p'}` (with `p < p'`)
quantify the *spacing* of the cascade and provide a natural arithmetic
signature of activity:

* **Active/active ratios** (`p, p' ∈ {3, 5, 7}`) are very close to `1`
  (the active primes form a *tight cluster* above the threshold `s = 1/2`).
  Concretely `γ_3/γ_5 ≈ 1.16`, `γ_3/γ_7 ≈ 1.36`, `γ_5/γ_7 ≈ 1.17`.

* **Active/inactive ratios** (e.g. `γ_7/γ_{11} ≈ 1.40`, `γ_7/γ_{13} ≈ 1.67`)
  open up significantly: the gap across the activity threshold leaves a
  measurable algebraic trace.

* **Global decay**: `γ_p / γ_{p'} > 1` whenever `p < p'` (immediate from
  strict monotonicity, `GammaMonotonicity`).

This file formalises both the exact rational ratios and the decimal
brackets (proven by `unfold ... ; norm_num`), and records the simple
global "ratio > 1" theorem as a closure result.

## Reference

Monograph Chapter 6, §"Cascade arithmétique de persistance",
extended table `tab:gamma_p_extended`.
-/

namespace PT.Holonomy

/-! ### Definition -/

/-- The PT cascade ratio `γ_p / γ_{p'}` between two anomalous dimensions
    at the fixed point `μ* = 15`, `q = qPT = 13/15`. -/
def gammaRatio (p p' : ℕ) : ℚ := gammaQ p / gammaQ p'

/-! ### Active/active ratios (tight cluster above the threshold) -/

/-- `γ_3 / γ_5 ≈ 1.16` (active/active, close cluster). -/
theorem gammaRatio_3_5_bracket :
    (1.15 : ℚ) < gammaRatio 3 5 ∧ gammaRatio 3 5 < (1.16 : ℚ) := by
  unfold gammaRatio gammaQ deltaQ qPT muStar
  refine ⟨?_, ?_⟩ <;> norm_num

/-- `γ_3 / γ_7 ≈ 1.36` (active/active across two steps). -/
theorem gammaRatio_3_7_bracket :
    (1.35 : ℚ) < gammaRatio 3 7 ∧ gammaRatio 3 7 < (1.36 : ℚ) := by
  unfold gammaRatio gammaQ deltaQ qPT muStar
  refine ⟨?_, ?_⟩ <;> norm_num

/-- `γ_5 / γ_7 ≈ 1.17` (active/active, second step). -/
theorem gammaRatio_5_7_bracket :
    (1.16 : ℚ) < gammaRatio 5 7 ∧ gammaRatio 5 7 < (1.17 : ℚ) := by
  unfold gammaRatio gammaQ deltaQ qPT muStar
  refine ⟨?_, ?_⟩ <;> norm_num

/-! ### Active/inactive ratios (gap across the threshold) -/

/-- `γ_7 / γ_{11} ≈ 1.40` (active/inactive: signature of the threshold
    crossing between `p = 7` and `p = 11`). -/
theorem gammaRatio_7_11_bracket :
    (1.39 : ℚ) < gammaRatio 7 11 ∧ gammaRatio 7 11 < (1.40 : ℚ) := by
  unfold gammaRatio gammaQ deltaQ qPT muStar
  refine ⟨?_, ?_⟩ <;> norm_num

/-- `γ_7 / γ_{13} ≈ 1.67` (active/inactive, two steps past the threshold). -/
theorem gammaRatio_7_13_bracket :
    (1.67 : ℚ) < gammaRatio 7 13 ∧ gammaRatio 7 13 < (1.68 : ℚ) := by
  unfold gammaRatio gammaQ deltaQ qPT muStar
  refine ⟨?_, ?_⟩ <;> norm_num

/-! ### Inactive/inactive ratios (deeper cascade) -/

/-- `γ_{11} / γ_{13} ≈ 1.20` (inactive/inactive, both below the threshold). -/
theorem gammaRatio_11_13_bracket :
    (1.19 : ℚ) < gammaRatio 11 13 ∧ gammaRatio 11 13 < (1.20 : ℚ) := by
  unfold gammaRatio gammaQ deltaQ qPT muStar
  refine ⟨?_, ?_⟩ <;> norm_num

/-- `γ_{13} / γ_{17} ≈ 1.46`. -/
theorem gammaRatio_13_17_bracket :
    (1.45 : ℚ) < gammaRatio 13 17 ∧ gammaRatio 13 17 < (1.46 : ℚ) := by
  unfold gammaRatio gammaQ deltaQ qPT muStar
  refine ⟨?_, ?_⟩ <;> norm_num

/-! ### Global decay: every cascade ratio exceeds 1 -/

/-- `gammaRatio p p' > 1` follows from `γ_p > γ_{p'}` and `γ_{p'} > 0`.
    Helper used to derive the active/active and active/inactive
    "ratio > 1" closure facts uniformly. -/
private theorem gammaRatio_gt_one_of_gt
    {p p' : ℕ} (h : gammaQ p > gammaQ p') (hpos : gammaQ p' > 0) :
    gammaRatio p p' > 1 := by
  unfold gammaRatio
  rw [gt_iff_lt, lt_div_iff₀ hpos, one_mul]
  exact h

-- `γ_5 > 0`, `γ_7 > 0`, `γ_{11} > 0`, `γ_{13} > 0` viennent de `GammaProduct`.
-- Seul `γ_{17} > 0` reste local ici.

/-- `γ_{17} > 0`. -/
theorem gammaQ_17_pos : gammaQ 17 > 0 := by
  unfold gammaQ deltaQ qPT muStar
  norm_num

/-- Closure: every consecutive cascade ratio exceeds `1`. -/
theorem gammaRatio_consecutive_gt_one :
    gammaRatio 3 5 > 1 ∧ gammaRatio 5 7 > 1
    ∧ gammaRatio 7 11 > 1 ∧ gammaRatio 11 13 > 1
    ∧ gammaRatio 13 17 > 1 :=
  ⟨gammaRatio_gt_one_of_gt gamma_3_gt_gamma_5  gammaQ_5_pos,
   gammaRatio_gt_one_of_gt gamma_5_gt_gamma_7  gammaQ_7_pos,
   gammaRatio_gt_one_of_gt gamma_7_gt_gamma_11 gammaQ_11_pos,
   gammaRatio_gt_one_of_gt gamma_11_gt_gamma_13 gammaQ_13_pos,
   gammaRatio_gt_one_of_gt gammaQ_13_gt_17     gammaQ_17_pos⟩

/-! ### Threshold-crossing signature

The *largest* active/active ratio (`γ_3/γ_7 ≈ 1.36`) is strictly smaller
than the *smallest* active/inactive ratio (`γ_7/γ_{11} ≈ 1.40`). This
"ratio gap" is a clean arithmetic signature of the threshold crossing:
ratios within the active cluster `{3,5,7}` are uniformly below the ratios
that bridge into the inactive sector `{11, 13, …}`. -/

/-- **Active/inactive ratio gap.** The arithmetic signature of the
    threshold crossing: `γ_3/γ_7 < γ_7/γ_{11}`. -/
theorem gammaRatio_threshold_gap :
    gammaRatio 3 7 < gammaRatio 7 11 := by
  unfold gammaRatio gammaQ deltaQ qPT muStar
  norm_num

/-! ### Headline -/

/-- **Headline (cascade ratios).** At `q = qPT = 13/15`, `μ* = 15`:

    * Active/active ratios live in a tight cluster:
      `γ_3/γ_5 ∈ (1.15, 1.16)`, `γ_5/γ_7 ∈ (1.16, 1.17)`,
      `γ_3/γ_7 ∈ (1.35, 1.36)`.
    * Active/inactive ratios open up at the threshold crossing:
      `γ_7/γ_{11} ∈ (1.39, 1.40)`, `γ_7/γ_{13} ∈ (1.67, 1.68)`.
    * Every consecutive cascade ratio exceeds `1` (the cascade is strictly
      decreasing).
    * The "ratio gap" `γ_3/γ_7 < γ_7/γ_{11}` is the clean arithmetic
      signature of the activity threshold at `s = 1/2`. -/
theorem gammaRatio_headline :
    ((1.15 : ℚ) < gammaRatio 3 5 ∧ gammaRatio 3 5 < (1.16 : ℚ))
    ∧ ((1.16 : ℚ) < gammaRatio 5 7 ∧ gammaRatio 5 7 < (1.17 : ℚ))
    ∧ ((1.35 : ℚ) < gammaRatio 3 7 ∧ gammaRatio 3 7 < (1.36 : ℚ))
    ∧ ((1.39 : ℚ) < gammaRatio 7 11 ∧ gammaRatio 7 11 < (1.40 : ℚ))
    ∧ ((1.67 : ℚ) < gammaRatio 7 13 ∧ gammaRatio 7 13 < (1.68 : ℚ))
    ∧ gammaRatio 3 7 < gammaRatio 7 11 :=
  ⟨gammaRatio_3_5_bracket, gammaRatio_5_7_bracket, gammaRatio_3_7_bracket,
   gammaRatio_7_11_bracket, gammaRatio_7_13_bracket,
   gammaRatio_threshold_gap⟩

end PT.Holonomy
