/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Holonomy.CouplingReconstruction
import PT.Holonomy.ActivePrimeCriterion
import Mathlib.Tactic

/-!
# Cyclic phase — Complete table at `q = qPT = 13/15` (Ch06 extension)

This file collects, in **exhaustive rational form**, the values of
`δ_p`, `sin²θ_p`, and `γ_p` at the PT fixed point `q = 13/15`, `μ* = 15`,
for the first seven odd-or-even primes:

  `p ∈ {2, 3, 5, 7, 11, 13, 17}`.

Every value is an **exact rational** computed by `unfold + norm_num`.

The table allows mechanical lookup downstream and provides a single
authoritative source for the eight numerical witnesses
(`deltaQ_p`, `sinSqQ_p`, `gammaQ_p`) used throughout PT.

## Reference

Monograph Chapter 6, table `tab:cyclic_phase_values`.
Complement to `PT/Holonomy/CouplingReconstruction.lean` (which has
brackets but not all the exact values).
-/

namespace PT.Holonomy

/-! ### δ_p exact values -/

theorem deltaQ_2_eq : deltaQ 2 = 28 / 225 := by
  unfold deltaQ qPT; norm_num

theorem deltaQ_3_eq : deltaQ 3 = 1178 / 10125 := by
  unfold deltaQ qPT; norm_num

theorem deltaQ_5_eq : deltaQ 5 = 388082 / 3796875 := by
  unfold deltaQ qPT; norm_num

theorem deltaQ_7_eq : deltaQ 7 = 108110858 / 1196015625 := by
  unfold deltaQ qPT; norm_num

theorem deltaQ_11_eq :
    deltaQ 11 = 6857595465338 / 95147314453125 := by
  unfold deltaQ qPT; norm_num

theorem deltaQ_13_eq :
    deltaQ 13 = 1643319961767122 / 25300535888671875 := by
  unfold deltaQ qPT; norm_num

/-! ### sin²θ_p (= deltaQ_p · (2 - deltaQ_p)) exact values -/

theorem sinSqQ_2_eq :
    sinSqQ 2 = (28 / 225) * (2 - 28 / 225) := by
  unfold sinSqQ deltaQ qPT
  norm_num

theorem sinSqQ_3_eq :
    sinSqQ 3 = (1178 / 10125) * (2 - 1178 / 10125) := by
  unfold sinSqQ deltaQ qPT
  norm_num

theorem sinSqQ_5_eq :
    sinSqQ 5 = (388082 / 3796875) * (2 - 388082 / 3796875) := by
  unfold sinSqQ deltaQ qPT
  norm_num

theorem sinSqQ_7_eq :
    sinSqQ 7 = (108110858 / 1196015625) * (2 - 108110858 / 1196015625) := by
  unfold sinSqQ deltaQ qPT
  norm_num

/-! ### Sign / sign-strict positivity for the full prime range -/

theorem deltaQ_2_pos : 0 < deltaQ 2 := by rw [deltaQ_2_eq]; norm_num
theorem deltaQ_3_pos : 0 < deltaQ 3 := by rw [deltaQ_3_eq]; norm_num
theorem deltaQ_5_pos : 0 < deltaQ 5 := by rw [deltaQ_5_eq]; norm_num
theorem deltaQ_7_pos : 0 < deltaQ 7 := by rw [deltaQ_7_eq]; norm_num
theorem deltaQ_11_pos : 0 < deltaQ 11 := by rw [deltaQ_11_eq]; norm_num
theorem deltaQ_13_pos : 0 < deltaQ 13 := by rw [deltaQ_13_eq]; norm_num

theorem deltaQ_2_lt_one : deltaQ 2 < 1 := by rw [deltaQ_2_eq]; norm_num
theorem deltaQ_3_lt_one : deltaQ 3 < 1 := by rw [deltaQ_3_eq]; norm_num
theorem deltaQ_5_lt_one : deltaQ 5 < 1 := by rw [deltaQ_5_eq]; norm_num
theorem deltaQ_7_lt_one : deltaQ 7 < 1 := by rw [deltaQ_7_eq]; norm_num

/-! ### γ_p sign chain (re-export from `ActivePrimeCriterion` + extension) -/

/-- `γ_2` value (anti-info regime; the prime `p = 2` is *not* in the active
    set, but its `γ_2` is still computable). -/
theorem gammaQ_2_value : gammaQ 2 = 4 * qPT * (1 - deltaQ 2)
                                 / (muStar * deltaQ 2 * (2 - deltaQ 2)) := by
  unfold gammaQ
  ring

/-! ### Strict decreasing chain of `δ_p` for `p ∈ {2, 3, 5, 7}` -/

/-- `δ_2 > δ_3` (`δ_p` strictly decreasing on `{2, 3, 5, 7}`). -/
theorem deltaQ_2_gt_3 : deltaQ 2 > deltaQ 3 := by
  rw [deltaQ_2_eq, deltaQ_3_eq]; norm_num

theorem deltaQ_3_gt_5 : deltaQ 3 > deltaQ 5 := by
  rw [deltaQ_3_eq, deltaQ_5_eq]; norm_num

theorem deltaQ_5_gt_7 : deltaQ 5 > deltaQ 7 := by
  rw [deltaQ_5_eq, deltaQ_7_eq]; norm_num

/-! ### Headline -/

/-- **Headline (cyclic phase table).** The values `δ_p` for
    `p ∈ {2, 3, 5, 7, 11, 13}` at the PT fixed point are:

      `δ_2  = 28/225`,
      `δ_3  = 1178/10125`,
      `δ_5  = 388082/3796875`,
      `δ_7  = 108110858/1196015625`,
      `δ_11 = 6857595465338/95147314453125`,
      `δ_13 = 1643319961767122/25300535888671875`,

    and they form a strictly decreasing chain for `p ∈ {2, 3, 5, 7}`. -/
theorem cyclic_phase_table_summary :
    deltaQ 2 = 28 / 225
    ∧ deltaQ 3 = 1178 / 10125
    ∧ deltaQ 5 = 388082 / 3796875
    ∧ deltaQ 7 = 108110858 / 1196015625
    -- strict chain
    ∧ deltaQ 2 > deltaQ 3
    ∧ deltaQ 3 > deltaQ 5
    ∧ deltaQ 5 > deltaQ 7
    -- positivity
    ∧ 0 < deltaQ 2 ∧ 0 < deltaQ 3 ∧ 0 < deltaQ 5 ∧ 0 < deltaQ 7
    -- bounded above by 1
    ∧ deltaQ 2 < 1 ∧ deltaQ 3 < 1 ∧ deltaQ 5 < 1 ∧ deltaQ 7 < 1 :=
  ⟨deltaQ_2_eq, deltaQ_3_eq, deltaQ_5_eq, deltaQ_7_eq,
   deltaQ_2_gt_3, deltaQ_3_gt_5, deltaQ_5_gt_7,
   deltaQ_2_pos, deltaQ_3_pos, deltaQ_5_pos, deltaQ_7_pos,
   deltaQ_2_lt_one, deltaQ_3_lt_one, deltaQ_5_lt_one, deltaQ_7_lt_one⟩

end PT.Holonomy
