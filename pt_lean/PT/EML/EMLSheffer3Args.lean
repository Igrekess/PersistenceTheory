/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.EML.QSheffer
import PT.EML.EMLIdentities
import PT.EML.EMLAlgebra
import Mathlib.Tactic

/-!
# EML — 3-argument extension and Sheffer composition (App O extension)

This file defines and studies a **ternary** extension of the EML operator,
along with composition identities for 3-argument forms:

* `eml3(x, y, z) := eml(eml(x, y), z) = exp(exp x - log y) - log z`
* `eml3'(x, y, z) := eml(x, eml(y, z)) = exp x - log(exp y - log z)`

These two ternary forms are **different** (EML is non-associative, see
`PT.EML.EMLAlgebra`), so `eml3(x, y, z) ≠ eml3'(x, y, z)` in general.

We compute specific values and exhibit the non-equality witness.

## Reference

Monograph Appendix O §"Compositions ternaires EML",
follow-up to `PT.EML.EMLAlgebra`.
-/

namespace PT.EML

open Real

/-! ### Two non-equivalent 3-arg compositions -/

/-- Left-associative ternary EML: `eml3(x, y, z) = eml(eml(x, y), z)`. -/
noncomputable def eml3 (x y z : ℝ) : ℝ := eml (eml x y) z

/-- Right-associative ternary EML: `eml3'(x, y, z) = eml(x, eml(y, z))`. -/
noncomputable def eml3' (x y z : ℝ) : ℝ := eml x (eml y z)

/-! ### Unit collapses -/

/-- `eml3(x, 1, 1) = exp(exp x)`. -/
theorem eml3_arg23_one (x : ℝ) :
    eml3 x 1 1 = Real.exp (Real.exp x) := by
  unfold eml3
  rw [eml_arg2_one, eml_arg2_one]

/-- `eml3'(0, 0, 1) = 1 - log(1 - 0) = 1`.
    (`eml 0 1 = exp 0 - log 1 = 1`, so `eml 0 (eml 0 1) = eml 0 1 = 1`.) -/
theorem eml3'_zero_zero_one :
    eml3' 0 0 1 = 1 := by
  unfold eml3' eml
  rw [Real.log_one, sub_zero, Real.exp_zero, Real.log_one, sub_zero]

/-! ### Non-associativity witness -/

/-- **Non-associativity witness.** `eml3(1, 1, 1) = exp(exp 1)` whereas
    `eml3'(1, 1, 1) = exp 1 - log(exp 1) = exp 1 - 1`. These differ
    since `exp(exp 1) > exp 1 > exp 1 - 1`. -/
theorem eml3_ne_eml3'_at_111 : eml3 1 1 1 ≠ eml3' 1 1 1 := by
  unfold eml3 eml3'
  -- eml3 1 1 1 = eml (eml 1 1) 1 = eml (exp 1) 1 = exp(exp 1)
  -- eml3' 1 1 1 = eml 1 (eml 1 1) = eml 1 (exp 1) = exp 1 - log(exp 1) = exp 1 - 1
  intro habs
  rw [eml_arg2_one, eml_arg2_one] at habs
  unfold eml at habs
  rw [Real.log_exp] at habs
  -- habs : exp(exp 1) = exp 1 - 1
  -- But exp 1 > 2 (Real.add_one_lt_exp), so exp(exp 1) > exp 2 > 7,
  -- while exp 1 - 1 < exp 1 - 0 = exp 1 ≈ 2.718.
  have hexp1_gt2 : Real.exp 1 > 2 := by
    have h := Real.add_one_lt_exp (x := 1) (by norm_num)
    linarith
  have hexpexp_gt2 : Real.exp (Real.exp 1) > Real.exp 2 := by
    exact Real.exp_lt_exp.mpr hexp1_gt2
  have hexp2_gt_e1 : Real.exp 2 > Real.exp 1 := by
    exact Real.exp_lt_exp.mpr (by norm_num)
  -- exp(exp 1) > exp 2 > exp 1, so exp(exp 1) > exp 1, so exp(exp 1) > exp 1 - 1
  linarith

/-! ### Symmetry-breaking comparison `eml3 ≠ eml3'` -/

/-- The two ternary compositions are **not** equal as functions of `(x, y, z)`. -/
theorem eml3_ne_eml3' : eml3 ≠ eml3' := by
  intro habs
  have := congrFun (congrFun (congrFun habs 1) 1) 1
  exact eml3_ne_eml3'_at_111 this

/-! ### Headline -/

/-- **Headline (3-arg EML composition).** Two distinct ternary EML
    compositions:

    * `eml3(x, y, z) := eml(eml(x, y), z)` (left-associative)
    * `eml3'(x, y, z) := eml(x, eml(y, z))` (right-associative)

    are explicitly non-equal: `eml3(1, 1, 1) = exp(exp 1) ≈ 15.15` whereas
    `eml3'(1, 1, 1) = exp 1 - 1 ≈ 1.72`. This concretises the
    non-associativity of `eml` from `EMLAlgebra` to the ternary level. -/
theorem eml_ternary_summary :
    eml3 1 1 1 = Real.exp (Real.exp 1)
    ∧ eml3' 0 0 1 = 1
    ∧ eml3 1 1 1 ≠ eml3' 1 1 1
    ∧ eml3 ≠ eml3' := by
  refine ⟨?_, eml3'_zero_zero_one, eml3_ne_eml3'_at_111, eml3_ne_eml3'⟩
  unfold eml3
  rw [eml_arg2_one, eml_arg2_one]

end PT.EML
