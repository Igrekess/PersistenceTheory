/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import Mathlib.Analysis.SpecialFunctions.Exp
import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.Tactic

/-!
# EML primitivity of `q_+` and `q_-` (App O)

**Statement (paper-level, App O §§"Primitivity of $q_-$", "Derivation of $q_+$").**

The EML (Exp–Minus–Log) operator of Odrzywołek (2026) is the binary
operation `eml(x, y) := exp(x) - log(y)` on `ℝ`. It is a continuous Sheffer
primitive: together with the constant `1`, it generates every elementary
function.

Within Persistence Theory, the two canonical PT branch parameters at level
`μ` are:

* `q_-(μ) := exp(-1/μ)` — the *thermal* parameter (Gibbs / soliton);
* `q_+(μ) := 1 - 2/μ` — the *statistical* parameter (geometric / max-entropy).

App O records:

1. **`q_- = eml(-1/μ, 1)` is a depth-1 EML primitive** (Prop O.1).
2. **`q_+ = eml(0, eml(2/μ, 1))` is a depth-2 EML construction** (Prop O.3).
3. **Bifurcation asymmetry:** `q_-` is "native" (depth 1), `q_+` is "derived"
   (depth ≥ 2). The numerical witness at the PT fixed point `μ = 15`
   confirms `q_-(15) ≠ q_+(15)`.

The depth lower bound for `q_+` (Prop O.2: any depth-1 EML tree fails) is a
structural claim about EML syntactic trees and is **left as a comment** —
it would require formalising the syntactic algebra of EML expressions, out
of scope for Vague 4.

## Reference

Appendix O of the monograph (`app:eml`), Propositions O.1 and O.3.
Odrzywołek, *EML: a continuous Sheffer primitive*, 2026.
-/

namespace PT.EML

open Real

/-- The EML operator `eml(x, y) := exp(x) - log(y)`. -/
noncomputable def eml (x y : ℝ) : ℝ := Real.exp x - Real.log y

/-- `eml(x, 1) = exp(x)`: the "log half" collapses on the unit. -/
@[simp] theorem eml_arg2_one (x : ℝ) : eml x 1 = Real.exp x := by
  unfold eml
  rw [Real.log_one, sub_zero]

/-- `eml(0, y) = 1 - log y`: the "exp half" collapses to `1`. -/
@[simp] theorem eml_arg1_zero (y : ℝ) : eml 0 y = 1 - Real.log y := by
  unfold eml
  rw [Real.exp_zero]

/-! ### The thermal parameter `q_-(μ) = e^{-1/μ}` -/

/-- The PT thermal branch parameter `q_-(μ) := exp(-1/μ)`. -/
noncomputable def qMinus (μ : ℝ) : ℝ := Real.exp (-(1 / μ))

/-- **Prop O.1 — `q_- = eml(-1/μ, 1)`: depth-1 EML primitive.**
    Total nodes = 3 (one EML, leaves `-1/μ` and `1`); depth = 1. -/
theorem qMinus_eml_primitive (μ : ℝ) :
    qMinus μ = eml (-(1 / μ)) 1 := by
  unfold qMinus
  rw [eml_arg2_one]

/-! ### The statistical parameter `q_+(μ) = 1 - 2/μ` -/

/-- The PT statistical branch parameter `q_+(μ) := 1 - 2/μ`. -/
noncomputable def qPlus (μ : ℝ) : ℝ := 1 - 2 / μ

/-- **Prop O.3 — `q_+(μ) = eml(0, eml(2/μ, 1))`: depth-2 EML construction.**

    Two-line verification:
    `eml(2/μ, 1) = exp(2/μ) - log 1 = exp(2/μ)`, and
    `eml(0, exp(2/μ)) = 1 - log(exp(2/μ)) = 1 - 2/μ = q_+(μ)`. -/
theorem qPlus_eml_depth2 (μ : ℝ) :
    qPlus μ = eml 0 (eml (2 / μ) 1) := by
  unfold qPlus
  rw [eml_arg2_one, eml_arg1_zero, Real.log_exp]

/-! ### Bifurcation asymmetry (numerical witness at the PT fixed point) -/

/-- **Bifurcation asymmetry at `μ = 15`.** `q_-(15) ≠ q_+(15)`.
    The proof uses the strict-convexity bound `exp(x) ≥ 1 + x`, applied
    at `x = -1/15`: `q_-(15) ≥ 14/15 > 13/15 = q_+(15)`. -/
theorem qMinus_ne_qPlus_at_15 : qMinus 15 ≠ qPlus 15 := by
  intro h
  -- q_+(15) = 1 - 2/15 = 13/15
  have hp : qPlus 15 = 13 / 15 := by unfold qPlus; norm_num
  -- q_-(15) = exp(-1/15) ≥ 1 + (-1/15) = 14/15 (Real.add_one_le_exp)
  have hm_lower : (14 : ℝ) / 15 ≤ qMinus 15 := by
    unfold qMinus
    have := Real.add_one_le_exp (-((1 : ℝ) / 15))
    have habs : 1 + -((1 : ℝ) / 15) = 14 / 15 := by norm_num
    linarith [this]
  rw [h, hp] at hm_lower
  norm_num at hm_lower

end PT.EML
