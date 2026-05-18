/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Holonomy.ActivePrimeCriterion
import PT.Holonomy.CyclicPhaseIdentity
import PT.Holonomy.GammaMonotonicity
import Mathlib.Tactic
import Mathlib.Data.Rat.Defs

/-!
# Coupling Reconstruction — `α_bare = ∏_{p ∈ {3,5,7}} sin²θ_p`

**Statement (paper-level, Ch09 §"Reconstruction des couplages"
and BA5 Step 4).** At the PT fixed point `μ* = 15` with `q_+ = 13/15`,
the *bare* sieve coupling is the product of the cyclic-phase squared sines
over the active primes:

$$\alpha_{\rm bare} \;=\; \prod_{p \in \{3, 5, 7\}} \sin^2 \theta_p(q_+).$$

Where the cyclic phase squared sine is given (Pythagorean identity,
`PT.Holonomy.CyclicPhaseIdentity.sin_sq_of_cos_eq_one_sub`) by

$$\sin^2 \theta_p \;=\; \delta_p (2 - \delta_p), \qquad
   \delta_p \;=\; \frac{1 - q_+^p}{p}.$$

The published reference values
(monograph `app_e_constants.tex`, PT_Geophysics `PROMPT_AUDIT_CRITIQUE.md`):

* `sin²θ_3 ≈ 0.2192`
* `sin²θ_5 ≈ 0.1940`
* `sin²θ_7 ≈ 0.1726`
* `α_bare = ∏ ≈ 7.342 × 10⁻³` (subsequent Koide–running corrections give
  the physical `α_EM ≈ 1/137.036`).

This module records:

1. The **algebraic definitions** of `sinSqQ p` (cyclic-phase squared sine
   at `q_+ = 13/15`) and `alphaBareQ` (product over active primes).
2. The **link to `CyclicPhaseIdentity`**: `sinSqQ p = δ_p(2 - δ_p)` exact.
3. **Positivity and ordering**: `0 < sinSqQ p < 1` for each active prime,
   strict decrease `sinSqQ 3 > sinSqQ 5 > sinSqQ 7` (parallel to
   `gamma_monotonicity` but at the level of `sin²`).
4. **Numerical bounds** matching the published values (rational interval
   bracketing each to 10⁻³).
5. **Bare coupling positivity** and **alpha bound**: `α_bare > 0`,
   `α_bare < 1/100` (loose), and the tight bracket
   `7.3 × 10⁻³ < α_bare < 7.4 × 10⁻³` decided by `norm_num`.

We work in `ℚ` for exact computation. The product formula is therefore
a purely arithmetic identity — no measure-theoretic or analytic ingredient.

## Reference

* Monograph Chapter 9, §"Reconstruction des couplages".
* `PT-WEBSITE/.../theorems/en/BA5.mdx`, Step 4 — Product over active primes.
* `app_e_constants.tex`, `\alpha_{\rm EM}` computation.
-/

namespace PT.Holonomy

/-! ### Cyclic-phase squared sine `sin²θ_p` at `q_+ = 13/15` -/

/-- The cyclic-phase squared sine evaluated at the PT branch parameter
    `q_+ = 13/15`:

    `sinSqQ p := δ_p (2 - δ_p)`, where `δ_p = (1 - q_+^p)/p` is
    `deltaQ p` defined in `ActivePrimeCriterion`.

    By the Pythagorean identity (cf. `sin_sq_of_cos_eq_one_sub`), this
    coincides with the real-valued `sin² θ_p` whenever
    `cos θ_p = 1 - δ_p`. -/
def sinSqQ (p : ℕ) : ℚ := deltaQ p * (2 - deltaQ p)

/-- Reformulation: `sinSqQ p = 1 - (1 - δ_p)²`. -/
theorem sinSqQ_eq_one_sub_sq (p : ℕ) :
    sinSqQ p = 1 - (1 - deltaQ p) ^ 2 := by
  unfold sinSqQ
  ring

/-! ### Numerical values of `sin²θ_p` for active primes

We verify exact rational brackets matching the published values:
`sin²θ_3 ≈ 0.2192`, `sin²θ_5 ≈ 0.1940`, `sin²θ_7 ≈ 0.1726`. -/

/-- `sin²θ_3` lies in `(0.219, 0.220)` (matches the published value `0.2192`). -/
theorem sinSq_3_bracket : 219 / 1000 < sinSqQ 3 ∧ sinSqQ 3 < 220 / 1000 := by
  unfold sinSqQ deltaQ qPT
  refine ⟨?_, ?_⟩ <;> norm_num

/-- `sin²θ_5` lies in `(0.193, 0.195)` (matches the published value `0.1940`). -/
theorem sinSq_5_bracket : 193 / 1000 < sinSqQ 5 ∧ sinSqQ 5 < 195 / 1000 := by
  unfold sinSqQ deltaQ qPT
  refine ⟨?_, ?_⟩ <;> norm_num

/-- `sin²θ_7` lies in `(0.172, 0.174)` (matches the published value `0.1726`). -/
theorem sinSq_7_bracket : 172 / 1000 < sinSqQ 7 ∧ sinSqQ 7 < 174 / 1000 := by
  unfold sinSqQ deltaQ qPT
  refine ⟨?_, ?_⟩ <;> norm_num

/-! ### Positivity and ordering -/

/-- `sin²θ_3 > 0`. -/
theorem sinSq_3_pos : 0 < sinSqQ 3 :=
  lt_trans (by norm_num : (0 : ℚ) < 219 / 1000) sinSq_3_bracket.1

/-- `sin²θ_5 > 0`. -/
theorem sinSq_5_pos : 0 < sinSqQ 5 :=
  lt_trans (by norm_num : (0 : ℚ) < 193 / 1000) sinSq_5_bracket.1

/-- `sin²θ_7 > 0`. -/
theorem sinSq_7_pos : 0 < sinSqQ 7 :=
  lt_trans (by norm_num : (0 : ℚ) < 172 / 1000) sinSq_7_bracket.1

/-- `sin²θ_3 < 1`. -/
theorem sinSq_3_lt_one : sinSqQ 3 < 1 :=
  lt_trans sinSq_3_bracket.2 (by norm_num : (220 / 1000 : ℚ) < 1)

/-- `sin²θ_5 < 1`. -/
theorem sinSq_5_lt_one : sinSqQ 5 < 1 :=
  lt_trans sinSq_5_bracket.2 (by norm_num : (195 / 1000 : ℚ) < 1)

/-- `sin²θ_7 < 1`. -/
theorem sinSq_7_lt_one : sinSqQ 7 < 1 :=
  lt_trans sinSq_7_bracket.2 (by norm_num : (174 / 1000 : ℚ) < 1)

/-- Strict decrease `sin²θ_3 > sin²θ_5`. -/
theorem sinSq_3_gt_sinSq_5 : sinSqQ 3 > sinSqQ 5 := by
  unfold sinSqQ deltaQ qPT
  norm_num

/-- Strict decrease `sin²θ_5 > sin²θ_7`. -/
theorem sinSq_5_gt_sinSq_7 : sinSqQ 5 > sinSqQ 7 := by
  unfold sinSqQ deltaQ qPT
  norm_num

/-- Full strict cascade `sin²θ_3 > sin²θ_5 > sin²θ_7`. -/
theorem sinSq_chain_active :
    sinSqQ 3 > sinSqQ 5 ∧ sinSqQ 5 > sinSqQ 7 :=
  ⟨sinSq_3_gt_sinSq_5, sinSq_5_gt_sinSq_7⟩

/-! ### Bare coupling `α_bare = ∏_{p ∈ {3,5,7}} sin²θ_p` -/

/-- The bare coupling at the PT fixed point: product of cyclic-phase squared
    sines over the active primes. -/
def alphaBareQ : ℚ := sinSqQ 3 * sinSqQ 5 * sinSqQ 7

/-- **Coupling reconstruction (definitional unfolding).**
    The bare coupling factorises into the product over active primes. -/
theorem alphaBareQ_eq_prod :
    alphaBareQ = sinSqQ 3 * sinSqQ 5 * sinSqQ 7 := rfl

/-- The bare coupling is strictly positive. -/
theorem alphaBareQ_pos : 0 < alphaBareQ := by
  unfold alphaBareQ
  exact mul_pos (mul_pos sinSq_3_pos sinSq_5_pos) sinSq_7_pos

/-- The bare coupling is strictly below `1/100` (loose bound from
    `sin²θ_3 < 1/4`, `sin²θ_5 < 1/4`, `sin²θ_7 < 1/4`). -/
theorem alphaBareQ_lt_one_hundredth : alphaBareQ < 1 / 100 := by
  unfold alphaBareQ sinSqQ deltaQ qPT
  norm_num

/-- **Tight bracket for the bare coupling.**
    `7.3 × 10⁻³ < α_bare < 7.4 × 10⁻³`, matching the published
    `α_bare ≈ 7.342 × 10⁻³`. The physical fine-structure constant
    `α_EM ≈ 1/137.036 ≈ 7.297 × 10⁻³` is recovered after Koide-running
    corrections (outside this kernel). -/
theorem alphaBareQ_bracket :
    73 / 10000 < alphaBareQ ∧ alphaBareQ < 74 / 10000 := by
  unfold alphaBareQ sinSqQ deltaQ qPT
  refine ⟨?_, ?_⟩ <;> norm_num

/-! ### Link to the cyclic phase identity (`CyclicPhaseIdentity`)

For each active prime `p ∈ {3, 5, 7}`, `sinSqQ p` is *by construction*
the rational form of `sin² θ_p` whenever `cos θ_p = 1 - δ_p`.
We record the explicit numerical equality `sinSqQ p = (Q-version of) sin²θ_p`. -/

/-- **Bridge to `CyclicPhaseIdentity` (real-valued).**
    For any real `θ_3` satisfying `cos θ_3 = 1 - (δ_3 : ℝ)`, we have
    `sin² θ_3 = (sinSqQ 3 : ℝ)`. -/
theorem sin_sq_cyclic_3 (θ : ℝ) (h : Real.cos θ = 1 - (deltaQ 3 : ℝ)) :
    Real.sin θ ^ 2 = (sinSqQ 3 : ℝ) := by
  have hpy := sin_sq_of_cos_eq_one_sub θ (deltaQ 3 : ℝ) h
  rw [hpy]
  unfold sinSqQ
  push_cast
  ring

/-- **Bridge to `CyclicPhaseIdentity` (real-valued, prime 5).** -/
theorem sin_sq_cyclic_5 (θ : ℝ) (h : Real.cos θ = 1 - (deltaQ 5 : ℝ)) :
    Real.sin θ ^ 2 = (sinSqQ 5 : ℝ) := by
  have hpy := sin_sq_of_cos_eq_one_sub θ (deltaQ 5 : ℝ) h
  rw [hpy]
  unfold sinSqQ
  push_cast
  ring

/-- **Bridge to `CyclicPhaseIdentity` (real-valued, prime 7).** -/
theorem sin_sq_cyclic_7 (θ : ℝ) (h : Real.cos θ = 1 - (deltaQ 7 : ℝ)) :
    Real.sin θ ^ 2 = (sinSqQ 7 : ℝ) := by
  have hpy := sin_sq_of_cos_eq_one_sub θ (deltaQ 7 : ℝ) h
  rw [hpy]
  unfold sinSqQ
  push_cast
  ring

/-! ### Headline: coupling reconstruction (real-valued product form) -/

/-- **Coupling reconstruction (real-valued).**
    Given three angles `θ_3, θ_5, θ_7` with `cos θ_p = 1 - δ_p`
    (the PT phase prescription), the product of squared sines equals
    the bare coupling:

    `(sin θ_3 · sin θ_5 · sin θ_7)² = α_bare`.

    Stated as the product of three `sin² θ_p`. -/
theorem coupling_reconstruction_real
    (θ3 θ5 θ7 : ℝ)
    (h3 : Real.cos θ3 = 1 - (deltaQ 3 : ℝ))
    (h5 : Real.cos θ5 = 1 - (deltaQ 5 : ℝ))
    (h7 : Real.cos θ7 = 1 - (deltaQ 7 : ℝ)) :
    Real.sin θ3 ^ 2 * Real.sin θ5 ^ 2 * Real.sin θ7 ^ 2 = (alphaBareQ : ℝ) := by
  rw [sin_sq_cyclic_3 θ3 h3, sin_sq_cyclic_5 θ5 h5, sin_sq_cyclic_7 θ7 h7]
  unfold alphaBareQ
  push_cast
  ring

end PT.Holonomy
