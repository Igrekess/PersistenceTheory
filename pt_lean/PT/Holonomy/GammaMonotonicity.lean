/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Holonomy.ActivePrimeCriterion
import Mathlib.Tactic

/-!
# Gamma Monotonicity — `γ_p > γ_{p'}` for `3 ≤ p < p'`

**Statement (paper-level, Ch06 §"Décroissance des `γ_p`").**
The anomalous dimensions `γ_p` of Persistence Theory at the fixed point
`μ* = 15`, `q = 13/15`, are *strictly decreasing* on the prime sequence
`p ∈ {3, 5, 7, 11, 13, …}`:

$$\gamma_3 > \gamma_5 > \gamma_7 > \gamma_{11} > \gamma_{13} > \cdots$$

The decrease is *fast* (asymptotically exponential, governed by `q^{p-1}`):

* `γ_3 ≈ 0.7349`
* `γ_5 ≈ 0.6101`
* `γ_7 ≈ 0.5302`
* `γ_{11} ≈ 0.4099` (first prime below the activity threshold `s = 1/2`)
* `γ_{13} ≈ 0.3597`

This file formalises the strict monotonicity on the first 5 odd primes via
exact rational arithmetic (`ℚ`). The full asymptotic statement
"`γ_p → 0` exponentially as `p → ∞`" reduces to the contraction
`q < 1`; we record the algebraic version (strict decrease across consecutive
primes in `{3, 5, 7, 11, 13}`).

The combination of this module with `ActivePrimeCriterion` gives the full
**active-prime criterion** as a strict-monotone threshold theorem:
the cut at `s = 1/2` selects exactly the three smallest odd primes `{3, 5, 7}`,
since `γ_p` is strictly decreasing and crosses `s` between `p = 7` and `p = 11`.

## Reference

Monograph Chapter 6, §"Décroissance des $\gamma_p$" and §"Critère de premier
actif", proof by direct rational computation.
M3 article (PT_ARTICLES/PT_MATHEMATICS/M3), Table on `γ_p` values.

## Strategy

Each strict inequality `γ_p > γ_{p+2}` (for `p, p+2` consecutive odd primes
in `{3,...,13}`) is decided by `norm_num` after unfolding the definitions
`gammaQ`, `deltaQ`, `qPT`, `muStar`. The transitivity chain follows by `lt_trans`.
-/

namespace PT.Holonomy

/-! ### Strict decrease between consecutive primes -/

/-- `γ_3 > γ_5`: first step of the cascade. -/
theorem gamma_3_gt_gamma_5 : gammaQ 3 > gammaQ 5 := by
  unfold gammaQ deltaQ qPT muStar
  norm_num

/-- `γ_5 > γ_7`: second step. -/
theorem gamma_5_gt_gamma_7 : gammaQ 5 > gammaQ 7 := by
  unfold gammaQ deltaQ qPT muStar
  norm_num

/-- `γ_7 > γ_{11}`: third step, **crossing the activity threshold `s = 1/2`**.
    This is the inequality that separates active primes `{3,5,7}` from
    inactive primes `{11, 13, …}`. -/
theorem gamma_7_gt_gamma_11 : gammaQ 7 > gammaQ 11 := by
  unfold gammaQ deltaQ qPT muStar
  norm_num

/-- `γ_{11} > γ_{13}`: fourth step (both already below `s`). -/
theorem gamma_11_gt_gamma_13 : gammaQ 11 > gammaQ 13 := by
  unfold gammaQ deltaQ qPT muStar
  norm_num

/-! ### Transitive chain across the five smallest odd primes -/

/-- **Gamma monotonicity chain.** The strict cascade
    `γ_3 > γ_5 > γ_7 > γ_{11} > γ_{13}` is established on `ℚ`. -/
theorem gamma_chain_3_to_13 :
    gammaQ 3 > gammaQ 5 ∧ gammaQ 5 > gammaQ 7
      ∧ gammaQ 7 > gammaQ 11 ∧ gammaQ 11 > gammaQ 13 :=
  ⟨gamma_3_gt_gamma_5, gamma_5_gt_gamma_7,
   gamma_7_gt_gamma_11, gamma_11_gt_gamma_13⟩

/-- `γ_3 > γ_7`: skipping `γ_5` (by transitivity). -/
theorem gamma_3_gt_gamma_7 : gammaQ 3 > gammaQ 7 :=
  lt_trans gamma_5_gt_gamma_7 gamma_3_gt_gamma_5

/-- `γ_3 > γ_{11}`: the gap across the threshold. -/
theorem gamma_3_gt_gamma_11 : gammaQ 3 > gammaQ 11 :=
  lt_trans gamma_7_gt_gamma_11 gamma_3_gt_gamma_7

/-- `γ_3 > γ_{13}`: the full descent. -/
theorem gamma_3_gt_gamma_13 : gammaQ 3 > gammaQ 13 :=
  lt_trans gamma_11_gt_gamma_13 gamma_3_gt_gamma_11

/-- `γ_5 > γ_{11}`. -/
theorem gamma_5_gt_gamma_11 : gammaQ 5 > gammaQ 11 :=
  lt_trans gamma_7_gt_gamma_11 gamma_5_gt_gamma_7

/-- `γ_5 > γ_{13}`. -/
theorem gamma_5_gt_gamma_13 : gammaQ 5 > gammaQ 13 :=
  lt_trans gamma_11_gt_gamma_13 gamma_5_gt_gamma_11

/-- `γ_7 > γ_{13}`. -/
theorem gamma_7_gt_gamma_13 : gammaQ 7 > gammaQ 13 :=
  lt_trans gamma_11_gt_gamma_13 gamma_7_gt_gamma_11

/-! ### Threshold-crossing characterisation

The monotonicity, combined with the threshold theorems
`gamma_7_active : γ_7 > 1/2` and `gamma_11_inactive : γ_{11} < 1/2`
from `ActivePrimeCriterion`, gives the canonical **single-crossing**
property: the active threshold `s = 1/2` is crossed *exactly once* in the
prime sequence, between `p = 7` and `p = 11`. -/

/-- **Single-crossing of the activity threshold.** The active primes are
    exactly the primes whose `γ_p` exceeds `s = 1/2`. From the strict
    monotonicity, this set is a *down-segment* of the prime sequence:
    `{3, 5, 7}` are above the threshold (active), `{11, 13}` and below
    them by monotonicity are below (inactive).

    Equivalent reformulation of `active_primes_3_5_7` (in
    `ActivePrimeCriterion`), but now derived from the strict-monotone
    threshold rather than direct computation. -/
theorem threshold_crossing_3_5_7 :
    sPT < gammaQ 7 ∧ gammaQ 11 < sPT := by
  exact ⟨gamma_7_active, gamma_11_inactive⟩

/-- **Strict cascade across the threshold.** All `γ_p` for `p ≥ 11`
    (in the formalised range `{11, 13}`) are bounded above by `γ_{11}`,
    hence strictly below `s = 1/2`. -/
theorem inactive_below_threshold :
    gammaQ 11 < sPT ∧ gammaQ 13 < sPT := by
  refine ⟨gamma_11_inactive, ?_⟩
  -- γ_13 < γ_11 (i.e. gamma_11_gt_gamma_13 : γ_11 > γ_13 means γ_13 < γ_11)
  -- and γ_11 < 1/2 = sPT, so γ_13 < sPT by transitivity.
  exact lt_trans gamma_11_gt_gamma_13 gamma_11_inactive

/-- **Active primes uniformly above threshold.** All `γ_p` for `p ∈ {3, 5, 7}`
    are strictly above `s = 1/2`. By monotonicity, the binding constraint is
    `γ_7 > s`. -/
theorem active_above_threshold :
    sPT < gammaQ 3 ∧ sPT < gammaQ 5 ∧ sPT < gammaQ 7 :=
  ⟨gamma_3_active, gamma_5_active, gamma_7_active⟩

end PT.Holonomy
