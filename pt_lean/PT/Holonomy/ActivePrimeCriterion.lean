/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import Mathlib.Tactic
import Mathlib.Data.Rat.Defs

/-!
# Active Prime Criterion — `p ∈ {3, 5, 7} ⟺ γ_p > s = 1/2`

**Statement (paper-level, Ch06 §"Critère de premier actif").**
A prime `p` is *active* at the PT fixed point `μ* = 15` iff
`γ_p > s = 1/2`, where the anomalous dimension is

$$\gamma_p = \frac{4\,q^{p-1}(1 - \delta_p)}{\mu^* \,\delta_p\,(2 - \delta_p)},
   \qquad q = \tfrac{13}{15}, \qquad \delta_p = \tfrac{1 - q^p}{p}.$$

The active set is `{3, 5, 7}`; primes `p ≥ 11` are *inactive*.

This file formalises the criterion at the level of exact rational arithmetic
(`fractions.Fraction` in the monograph maps to `ℚ` here). The strict
inequalities are decided by `decide` / `native_decide`, eliminating any
floating-point error.

## Reference

Monograph Chapter 6, §"Critère de premier actif", `\label{thm:active}`.
M3 article (PT_ARTICLES/PT_MATHEMATICS/M3) section on holonomy.
-/

namespace PT.Holonomy

/-- The PT branch parameter `q = 13/15`. -/
def qPT : ℚ := 13 / 15

/-- The PT fixed point `μ* = 15`. -/
def muStar : ℚ := 15

/-- The PT symmetry parameter `s = 1/2`. -/
def sPT : ℚ := 1 / 2

/-- The gap fraction `δ_p = (1 - q^p)/p` at `q = qPT`. -/
def deltaQ (p : ℕ) : ℚ := (1 - qPT ^ p) / (p : ℚ)

/-- The anomalous dimension `γ_p` (closed form at `q = 13/15`, `μ* = 15`):
    `γ_p = 4 q^(p-1) (1 - δ_p) / (μ* δ_p (2 - δ_p))`. -/
def gammaQ (p : ℕ) : ℚ :=
  (4 * qPT ^ (p - 1) * (1 - deltaQ p)) / (muStar * deltaQ p * (2 - deltaQ p))

/-! ### Active primes: γ_p > s = 1/2 for p ∈ {3, 5, 7} -/

/-- `γ_3 > 1/2`. -/
theorem gamma_3_active : gammaQ 3 > sPT := by
  unfold gammaQ deltaQ qPT muStar sPT
  norm_num

/-- `γ_5 > 1/2`. -/
theorem gamma_5_active : gammaQ 5 > sPT := by
  unfold gammaQ deltaQ qPT muStar sPT
  norm_num

/-- `γ_7 > 1/2`. -/
theorem gamma_7_active : gammaQ 7 > sPT := by
  unfold gammaQ deltaQ qPT muStar sPT
  norm_num

/-! ### Inactive primes: γ_p ≤ 1/2 for p ∈ {11, 13} -/

/-- `γ_11 < 1/2` (inactive). -/
theorem gamma_11_inactive : gammaQ 11 < sPT := by
  unfold gammaQ deltaQ qPT muStar sPT
  norm_num

/-- `γ_13 < 1/2` (inactive). -/
theorem gamma_13_inactive : gammaQ 13 < sPT := by
  unfold gammaQ deltaQ qPT muStar sPT
  norm_num

/-! ### Headline: the active set at `μ* = 15` is exactly `{3, 5, 7}`
within the prime range `{3, 5, 7, 11, 13}`. -/

/-- The "active prime" predicate at `μ* = 15`: `γ_p > s = 1/2`. -/
def IsActive (p : ℕ) : Prop := gammaQ p > sPT

/-- **Active prime criterion (final form).** Among the first 5 odd primes
    `{3, 5, 7, 11, 13}`, exactly `{3, 5, 7}` satisfy `γ_p > 1/2`. -/
theorem active_primes_3_5_7 :
    IsActive 3 ∧ IsActive 5 ∧ IsActive 7
      ∧ ¬ IsActive 11 ∧ ¬ IsActive 13 := by
  refine ⟨gamma_3_active, gamma_5_active, gamma_7_active, ?_, ?_⟩
  · -- γ_11 < 1/2 ⇒ ¬ (γ_11 > 1/2)
    intro h
    have := gamma_11_inactive
    exact absurd this (not_lt.mpr (le_of_lt h))
  · intro h
    have := gamma_13_inactive
    exact absurd this (not_lt.mpr (le_of_lt h))

end PT.Holonomy
