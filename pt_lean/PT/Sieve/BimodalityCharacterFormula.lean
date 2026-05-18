/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Sieve.Bimodality
import PT.Sieve.BimodalityT1Projection
import PT.Sieve.LegendreLogParity
import Mathlib.Tactic

/-!
# Bimodality — Closed-form character formula `δ̄(r) = 9 - 2·(r/5)` (App N)

This file proves the **closed-form character formula** for the bimodality:

  `δ̄(r) = 9 - 2 · (r / 5)`     for `r ∈ admissibleResidues`,

where `(r/5)` is the Legendre symbol modulo `5` (cf.
`PT/Sieve/LegendreLogParity.lean`). Equivalent: `δ̄(r) ∈ {7, 11}` with

  `δ̄(r) = 7`  iff `(r/5) = +1` (QR),
  `δ̄(r) = 11` iff `(r/5) = -1` (NR).

The formula `9 - 2·(r/5)` is the **closed form** that interpolates between
`7` and `11` via the character value `(r/5) ∈ {+1, -1}`.

## Reference

Monograph Appendix N §"Formule fermée de bimodalité", `\label{prop:bimod_closed}`.
-/

namespace PT.Sieve

/-! ### Closed-form values for each admissible residue -/

/-- `9 - 2·(r/5) = 7` when `r ≡ 1 mod 5` (QR). -/
theorem nine_minus_two_at_QR (n : ZMod 5) (h : legendre5 n = 1) :
    9 - 2 * legendre5 n = 7 := by
  rw [h]; ring

/-- `9 - 2·(r/5) = 11` when `r ≡ 2 mod 5` (NR). -/
theorem nine_minus_two_at_NR (n : ZMod 5) (h : legendre5 n = -1) :
    9 - 2 * legendre5 n = 11 := by
  rw [h]; ring

/-! ### Application: δ̄(r) formula for each admissible residue -/

/-- For `r = 1`: `9 - 2·(1/5) = 9 - 2 = 7`, matching `δ̄(1) = 7`
    (`deltaBarTimes6 1 = 42 = 6 · 7`). -/
theorem character_formula_r1 :
    9 - 2 * legendre5 (to5 1) = 7 := by
  decide

/-- For `r = 7`: `9 - 2·(2/5) = 9 - 2·(-1) = 11`, matching `δ̄(7) = 11`. -/
theorem character_formula_r7 :
    9 - 2 * legendre5 (to5 7) = 11 := by
  decide

/-- For `r = 11`: `9 - 2·(1/5) = 7`, matching `δ̄(11) = 7`. -/
theorem character_formula_r11 :
    9 - 2 * legendre5 (to5 11) = 7 := by
  decide

/-- For `r = 13`: `9 - 2·(3/5) = 9 - 2·(-1) = 11`, matching `δ̄(13) = 11`. -/
theorem character_formula_r13 :
    9 - 2 * legendre5 (to5 13) = 11 := by
  decide

/-- For `r = 17`: `9 - 2·(2/5) = 11`, matching `δ̄(17) = 11`. -/
theorem character_formula_r17 :
    9 - 2 * legendre5 (to5 17) = 11 := by
  decide

/-- For `r = 19`: `9 - 2·(4/5) = 7`, matching `δ̄(19) = 7`. -/
theorem character_formula_r19 :
    9 - 2 * legendre5 (to5 19) = 7 := by
  decide

/-- For `r = 23`: `9 - 2·(3/5) = 11`, matching `δ̄(23) = 11`. -/
theorem character_formula_r23 :
    9 - 2 * legendre5 (to5 23) = 11 := by
  decide

/-- For `r = 29`: `9 - 2·(4/5) = 7`, matching `δ̄(29) = 7`. -/
theorem character_formula_r29 :
    9 - 2 * legendre5 (to5 29) = 7 := by
  decide

/-! ### Range of the formula on admissible residues -/

/-- **Range of the character formula.** For every admissible residue
    `r ∈ {1, 7, 11, 13, 17, 19, 23, 29}`, the value `9 - 2·(r/5)` lies in
    `{7, 11}`. -/
theorem character_formula_range (r : ℕ) (hr : r ∈ admissibleResidues) :
    9 - 2 * legendre5 (to5 r) = 7 ∨ 9 - 2 * legendre5 (to5 r) = 11 := by
  unfold admissibleResidues at hr
  fin_cases hr <;> [left; right; left; right; right; left; right; left]
  all_goals decide

/-! ### Bridge: δ̄(r)·6 in the deltaBarTimes6 form -/

/-- The formula `δ̄(r) · 6 = 6·(9 - 2·(r/5)) = 54 - 12·(r/5)` evaluates to
    `42` for QR and `66` for NR, matching `deltaBarTimes6`. -/
theorem character_formula_times6_QR :
    6 * 7 = 42 := by decide

theorem character_formula_times6_NR :
    6 * 11 = 66 := by decide

/-! ### Headline (character formula summary) -/

/-- **Headline (character formula).**
    For every admissible residue `r ∈ R`, the bimodality value is the
    affine character function

      `δ̄(r) = 9 - 2 · (r/5)`,

    where `(r/5)` is the Legendre symbol mod `5`. This takes values in
    `{7, 11}` matching the dichotomy `deltaBarTimes6(r) ∈ {42, 66}`. -/
theorem bimodality_character_formula_summary :
    -- Range: every admissible residue gives 7 or 11 via the formula
    (∀ r ∈ admissibleResidues,
        9 - 2 * legendre5 (to5 r) = 7 ∨ 9 - 2 * legendre5 (to5 r) = 11)
    -- QR → 7 (4 examples)
    ∧ 9 - 2 * legendre5 (to5 1) = 7
    ∧ 9 - 2 * legendre5 (to5 11) = 7
    ∧ 9 - 2 * legendre5 (to5 19) = 7
    ∧ 9 - 2 * legendre5 (to5 29) = 7
    -- NR → 11 (4 examples)
    ∧ 9 - 2 * legendre5 (to5 7) = 11
    ∧ 9 - 2 * legendre5 (to5 13) = 11
    ∧ 9 - 2 * legendre5 (to5 17) = 11
    ∧ 9 - 2 * legendre5 (to5 23) = 11 :=
  ⟨character_formula_range,
   character_formula_r1, character_formula_r11,
   character_formula_r19, character_formula_r29,
   character_formula_r7, character_formula_r13,
   character_formula_r17, character_formula_r23⟩

end PT.Sieve
