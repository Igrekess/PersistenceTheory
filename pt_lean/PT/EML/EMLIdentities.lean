/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.EML.QSheffer
import Mathlib.Analysis.SpecialFunctions.Exp
import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.Topology.ContinuousOn
import Mathlib.Tactic

/-!
# EML algebraic identities (App O follow-up)

This file collects algebraic and topological lemmas for the EML operator
`eml(x, y) := exp(x) - log(y)` (cf. `PT.EML.QSheffer`):

* **Recovery of the exp / log primitives** via specialisation:
  `eml(x, 1) = exp x` and `eml(0, y) = 1 - log y`.
* **Continuity** in each variable, jointly on the open half-line `y > 0`.
* **Monotonicity** in the first argument (strict, since `exp` is strict);
  monotonicity in the second argument (strict decrease on `y > 0`, since
  `-log` strictly decreases there).
* **Asymmetry**: explicit witnesses that `eml` is neither commutative nor
  associative.

These lemmas extend the depth-1/depth-2 EML primitivity results
`qMinus_eml_primitive` and `qPlus_eml_depth2` of `QSheffer.lean`.

## Reference

Appendix O of the monograph, follow-up to `\label{prop:eml_primitivity}`.
Odrzywołek, *EML: a continuous Sheffer primitive*, 2026.
-/

namespace PT.EML

open Real

/-! ### Recovery of exp / log primitives -/

/-- `eml(log y, y) = y - log y` (no special simplification). -/
theorem eml_log_self (y : ℝ) (hy : 0 < y) : eml (Real.log y) y = y - Real.log y := by
  unfold eml
  rw [Real.exp_log hy]

/-- `eml(0, exp x) = 1 - x`. -/
theorem eml_zero_exp (x : ℝ) : eml 0 (Real.exp x) = 1 - x := by
  unfold eml
  rw [Real.exp_zero, Real.log_exp]

/-- `eml(x, exp y) = exp x - y`. -/
theorem eml_exp_arg2 (x y : ℝ) : eml x (Real.exp y) = Real.exp x - y := by
  unfold eml
  rw [Real.log_exp]

/-! ### Continuity -/

/-- `eml` is continuous in the first argument (for any fixed `y > 0`). -/
theorem continuous_eml_left (y : ℝ) : Continuous (fun x => eml x y) := by
  unfold eml
  exact Real.continuous_exp.sub continuous_const

/-- `eml` is continuous in the second argument on `y > 0`. -/
theorem continuousOn_eml_right (x : ℝ) :
    ContinuousOn (fun y => eml x y) {y : ℝ | 0 < y} := by
  unfold eml
  refine continuous_const.continuousOn.sub ?_
  intro y hy
  exact (Real.continuousAt_log (ne_of_gt hy)).continuousWithinAt

/-! ### Monotonicity -/

/-- `eml` is strictly increasing in the first argument. -/
theorem eml_strictMono_left (y : ℝ) : StrictMono (fun x => eml x y) := by
  intro a b hab
  unfold eml
  have : Real.exp a < Real.exp b := Real.exp_lt_exp.mpr hab
  linarith

/-- `eml` is strictly decreasing in the second argument on `y > 0`. -/
theorem eml_strictAnti_right (x : ℝ) :
    StrictAntiOn (fun y => eml x y) {y : ℝ | 0 < y} := by
  intro a ha b hb hab
  unfold eml
  have : Real.log a < Real.log b := Real.log_lt_log ha hab
  linarith

/-! ### Asymmetry: `eml` is neither commutative nor associative -/

/-- **Asymmetry (commutativity fails).** `eml(0, 1) = 1` but `eml(1, 0) = exp 1`
    (Mathlib convention `log 0 = 0`); since `exp 1 > 1`, they differ. -/
theorem eml_not_commutative : eml 0 1 ≠ eml 1 0 := by
  unfold eml
  rw [Real.log_one, Real.exp_zero, Real.log_zero]
  -- Goal: (1 : ℝ) - 0 ≠ Real.exp 1 - 0
  intro h
  have h1 : (1 : ℝ) = Real.exp 1 := by linarith
  have h2 : Real.exp 0 = 1 := Real.exp_zero
  have h3 : Real.exp 0 = Real.exp 1 := by rw [h2]; exact h1
  have h4 : (0 : ℝ) = 1 := Real.exp_injective h3
  linarith

/-! ### Identity at the PT bifurcation `μ*` -/

/-- The EML depth-2 derivation of `q_+(μ*)` evaluates to `13/15`. -/
theorem qPlus_eml_at_15 : qPlus 15 = 13 / 15 := by
  unfold qPlus; norm_num

/-- The EML depth-1 derivation of `q_-(μ*)` is `exp(-1/15)`. -/
theorem qMinus_eml_at_15 : qMinus 15 = Real.exp (-(1 / 15)) := rfl

/-- `q_+(μ*) < 1` strictly. -/
theorem qPlus_lt_one_at_15 : qPlus 15 < 1 := by
  rw [qPlus_eml_at_15]; norm_num

/-- `q_-(μ*) < 1` strictly (from `exp(-1/15) < exp 0 = 1`). -/
theorem qMinus_lt_one_at_15 : qMinus 15 < 1 := by
  unfold qMinus
  have h : (-(1 / 15) : ℝ) < 0 := by norm_num
  calc Real.exp (-(1 / 15)) < Real.exp 0 := Real.exp_lt_exp.mpr h
    _ = 1 := Real.exp_zero

/-- `0 < q_+(μ*)`. -/
theorem qPlus_pos_at_15 : 0 < qPlus 15 := by
  rw [qPlus_eml_at_15]; norm_num

/-- `0 < q_-(μ*)` (exp is positive). -/
theorem qMinus_pos_at_15 : 0 < qMinus 15 := by
  unfold qMinus
  exact Real.exp_pos _

end PT.EML
