/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Holonomy.ActivePrimeCriterion
import PT.Holonomy.ActivePrimeMargins
import Mathlib.Tactic

/-!
# γ_p — Extended exact table for `p ∈ {3, 5, 7, 11, 13, 17, 19, 23}`

This file extends `PT.Holonomy.ActivePrimeMargins` with the **exact rational
values** of `γ_p` for the first eight odd primes at the PT fixed point
`μ* = 15`, `q = 13/15`:

  `γ_3  = 4536129 / 5616704`
  `γ_5  = 486792684365 / 699097512194`
  `γ_7  = 2827519972576117 / 4748396022746468`
  `γ_11 = 133886362559711681370529793 / 314484242174863298491777064`
  `γ_13 = 7165181795537222421135938678509
          / 20113312710444408001435070506154`
  `γ_17 = 17930424521713546612756886650523475297101
          / 73248985204208930372556937795746032669534`
  `γ_19 = 855713881248444794557952426258479962643584577
          / 4253582256587939600561313383855525888781591824`
  `γ_23 = 1826993706763809744070478198227115452180065072741739889
          / 13653478210144391635849057056556448578741295121252037604`

The table provides a single authoritative source for downstream comparisons,
extending `ActivePrimeMargins.gamma_*_margin_exact` (which gives only the
`γ_p - s` differences).

## Reference

Monograph Chapter 6, extended table `tab:gamma_p_extended`.
-/

namespace PT.Holonomy

/-! ### γ_p exact values (active primes {3, 5, 7}) -/

theorem gammaQ_3_exact : gammaQ 3 = 4536129 / 5616704 := by
  unfold gammaQ deltaQ qPT muStar
  norm_num

theorem gammaQ_5_exact : gammaQ 5 = 486792684365 / 699097512194 := by
  unfold gammaQ deltaQ qPT muStar
  norm_num

theorem gammaQ_7_exact :
    gammaQ 7 = 2827519972576117 / 4748396022746468 := by
  unfold gammaQ deltaQ qPT muStar
  norm_num

/-! ### γ_p exact values (inactive primes {11, 13, 17}) -/

theorem gammaQ_11_exact :
    gammaQ 11
      = 133886362559711681370529793
        / 314484242174863298491777064 := by
  unfold gammaQ deltaQ qPT muStar
  norm_num

theorem gammaQ_13_exact :
    gammaQ 13
      = 7165181795537222421135938678509
        / 20113312710444408001435070506154 := by
  unfold gammaQ deltaQ qPT muStar
  norm_num

theorem gammaQ_17_exact :
    gammaQ 17
      = 17930424521713546612756886650523475297101
        / 73248985204208930372556937795746032669534 := by
  unfold gammaQ deltaQ qPT muStar
  norm_num

/-! ### γ_p exact values (deeper inactive primes {19, 23}) -/

theorem gammaQ_19_exact :
    gammaQ 19
      = 855713881248444794557952426258479962643584577
        / 4253582256587939600561313383855525888781591824 := by
  unfold gammaQ deltaQ qPT muStar
  norm_num

theorem gammaQ_23_exact :
    gammaQ 23
      = 1826993706763809744070478198227115452180065072741739889
        / 13653478210144391635849057056556448578741295121252037604 := by
  unfold gammaQ deltaQ qPT muStar
  norm_num

/-! ### Sign chain `> 1/2` (active) vs `< 1/2` (inactive) -/

theorem gammaQ_19_inactive : gammaQ 19 < sPT := by
  unfold gammaQ deltaQ qPT muStar sPT
  norm_num

theorem gammaQ_23_inactive : gammaQ 23 < sPT := by
  unfold gammaQ deltaQ qPT muStar sPT
  norm_num

/-! ### Strict decreasing chain across all 8 -/

theorem gammaQ_13_gt_17 : gammaQ 13 > gammaQ 17 := by
  unfold gammaQ deltaQ qPT muStar
  norm_num

theorem gammaQ_17_gt_19 : gammaQ 17 > gammaQ 19 := by
  unfold gammaQ deltaQ qPT muStar
  norm_num

theorem gammaQ_19_gt_23 : gammaQ 19 > gammaQ 23 := by
  unfold gammaQ deltaQ qPT muStar
  norm_num

/-! ### Headline -/

/-- **Headline (extended γ table).** At `q = qPT = 13/15`, `μ* = 15`:

    * Active primes `{3, 5, 7}` satisfy `γ_p > s = 1/2`.
    * Inactive primes `{11, 13, 17, 19, 23}` satisfy `γ_p < s = 1/2`.
    * The full strict ordering `γ_3 > γ_5 > γ_7 > 1/2 > γ_11 > γ_13 > γ_17
      > γ_19 > γ_23` holds.

    Combined with the smaller `gamma_*_margin_exact` rational equalities, the
    full numerical content of the PT γ-cascade is laid out exhaustively up
    to `p = 23`. -/
theorem gamma_extended_table_summary :
    gammaQ 3 > sPT ∧ gammaQ 5 > sPT ∧ gammaQ 7 > sPT
    ∧ sPT > gammaQ 11 ∧ sPT > gammaQ 13 ∧ sPT > gammaQ 17
    ∧ sPT > gammaQ 19 ∧ sPT > gammaQ 23
    ∧ gammaQ 13 > gammaQ 17
    ∧ gammaQ 17 > gammaQ 19
    ∧ gammaQ 19 > gammaQ 23 :=
  ⟨gamma_3_active, gamma_5_active, gamma_7_active,
   gamma_11_inactive, gamma_13_inactive, gamma_17_inactive,
   gammaQ_19_inactive, gammaQ_23_inactive,
   gammaQ_13_gt_17, gammaQ_17_gt_19, gammaQ_19_gt_23⟩

end PT.Holonomy
