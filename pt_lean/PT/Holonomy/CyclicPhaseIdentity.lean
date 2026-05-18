/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Basic
import Mathlib.Tactic

/-!
# Cyclic Phase Identity — `sin²θ_p = δ_p (2 - δ_p)`

**Statement (paper-level, Ch06 §6.2).** The cyclic-phase angle `θ_p` of
Persistence Theory is *defined* by `cos θ_p = 1 - δ_p`, where `δ_p ∈ (0,1)`
is the gap fraction. Consequently, by the Pythagorean identity,

$$\sin^2 \theta_p \;=\; 1 - \cos^2 \theta_p
   \;=\; 1 - (1 - \delta_p)^2 \;=\; \delta_p (2 - \delta_p).$$

This file records the algebraic identity in two flavours:

* `sin_sq_of_cos_eq_one_sub` — purely scalar: for any real `δ` and `θ` with
  `cos θ = 1 - δ`, we have `sin² θ = δ (2 - δ)`.
* `delta_two_minus_delta_eq` — pure-algebra version: `1 - (1-δ)² = δ(2-δ)`.

## Reference

Chapter 6, §6.2 of the monograph (`thm:cyclic-phase`). Identity is exact and
unconditional; it does not depend on the value of `δ_p`.
-/

namespace PT.Holonomy

open Real

/-- Pure-algebra identity: `1 - (1 - δ)² = δ (2 - δ)`. -/
theorem one_sub_one_sub_sq (δ : ℝ) : 1 - (1 - δ)^2 = δ * (2 - δ) := by ring

/-- **Cyclic phase identity (Pythagorean form).**
    For any real number `θ` and any `δ` with `cos θ = 1 - δ`,
    `sin² θ = δ (2 - δ)`. -/
theorem sin_sq_of_cos_eq_one_sub (θ δ : ℝ) (h : Real.cos θ = 1 - δ) :
    Real.sin θ ^ 2 = δ * (2 - δ) := by
  have hpyth : Real.sin θ ^ 2 + Real.cos θ ^ 2 = 1 := Real.sin_sq_add_cos_sq θ
  have hcos2 : Real.cos θ ^ 2 = (1 - δ)^2 := by rw [h]
  have : Real.sin θ ^ 2 = 1 - (1 - δ)^2 := by linarith
  rw [this, one_sub_one_sub_sq]

/-- The same identity rephrased via the gap-fraction parameterisation
    `δ_p(q) = (1 - q^p)/p` for a generic real `q` and natural prime exponent
    `p`. The formula is `sin² θ_p = δ_p (2 - δ_p)` whenever the angle is
    set by `cos θ_p = 1 - δ_p`. -/
theorem cyclic_phase_via_delta (θ : ℝ) (q : ℝ) (p : ℕ)
    (h : Real.cos θ = 1 - (1 - q^p) / p) :
    Real.sin θ ^ 2 = ((1 - q^p) / p) * (2 - (1 - q^p) / p) :=
  sin_sq_of_cos_eq_one_sub θ ((1 - q^p) / p) h

end PT.Holonomy
