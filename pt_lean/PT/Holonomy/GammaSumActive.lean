/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Holonomy.ActivePrimeCriterion
import PT.Holonomy.ActivePrimeMargins
import PT.Holonomy.GammaTablesExtended
import PT.Holonomy.GammaProduct
import Mathlib.Tactic

/-!
# Sum of anomalous dimensions `Σ γ_p` over the active cascade `{3, 5, 7}`

**Statement (paper-level, Ch06 §"Cascade arithmétique" / Ch09 supplement).**
Where `PT.Holonomy.GammaProduct` records the *product* `Γ_active = γ_3 γ_5 γ_7`,
this module records the **additive companion**

$$\Sigma_{\rm active} \;:=\; \sum_{p \in \{3, 5, 7\}} \gamma_p
   \;=\; \gamma_3 + \gamma_5 + \gamma_7
   \;\approx\; 0.8076 + 0.6963 + 0.5954 \;\approx\; 2.0994.$$

The sum cleanly separates from the trivial bounds:

* `Σ_active > 3·s = 3/2` — each `γ_p > 1/2 = s`, so the sum exceeds `3s`.
* `Σ_active < 3`         — each `γ_p < 1`, so the sum stays below `3`.
* `Σ_active > 2 = 1/s`   — strictly above the inverse-symmetry value `1/s`.

The product `Σ_active · μ*` with `μ* = 15` yields a sharp rational
bracket `(31.49, 31.50)`, recording the "γ-sum × fixed-point" invariant
discussed in the Ch09 cascade supplement.

Finally, the 5-prime extended sum `Σ_ext := γ_3+γ_5+γ_7+γ_11+γ_13` is
bracketed at `(2.881, 2.882)` and strictly exceeds `Σ_active` because
`γ_11, γ_13 > 0`.

All proofs are pure rational arithmetic — `unfold` + `norm_num`.

## Reference

* Monograph Chapter 6, §"Cascade arithmétique".
* Monograph Chapter 9, §"Reconstruction des couplages" supplement on
  additive holonomy invariants.
* Companion product file `PT.Holonomy.GammaProduct`.
-/

namespace PT.Holonomy

/-! ### Definition: sum over the active cascade `{3, 5, 7}` -/

/-- The PT anomalous-dimension *sum* over the active primes `{3, 5, 7}`
    at `q_+ = 13/15`, `μ* = 15`:

    `gammaSumActive := γ_3 + γ_5 + γ_7`. -/
def gammaSumActive : ℚ := gammaQ 3 + gammaQ 5 + gammaQ 7

/-- Definitional unfolding: the active sum is literally the triple sum. -/
theorem gammaSumActive_eq_sum :
    gammaSumActive = gammaQ 3 + gammaQ 5 + gammaQ 7 := rfl

/-! ### Tight rational bracket on the active sum -/

/-- **Tight bracket for `Σ_active`.**
    `2.099 < γ_3 + γ_5 + γ_7 < 2.100`. Numerical value `≈ 2.0994`. -/
theorem gammaSumActive_bracket :
    2099 / 1000 < gammaSumActive ∧ gammaSumActive < 2100 / 1000 := by
  unfold gammaSumActive gammaQ deltaQ qPT muStar
  refine ⟨?_, ?_⟩ <;> norm_num

/-! ### Positivity and basic comparisons -/

/-- **Positivity** `0 < Σ_active`. -/
theorem gammaSumActive_pos : 0 < gammaSumActive := by
  unfold gammaSumActive
  exact add_pos (add_pos gammaQ_3_pos gammaQ_5_pos) gammaQ_7_pos

/-- **Lower bound from active criterion.** `Σ_active > 3 · s = 3/2`.

    Direct sum of the three strict inequalities `γ_p > s = 1/2`
    (`gamma_3_active`, `gamma_5_active`, `gamma_7_active`). -/
theorem gammaSumActive_gt_three_s : gammaSumActive > 3 * sPT := by
  unfold gammaSumActive sPT
  -- 3 * (1/2) = 3/2; show γ_3 + γ_5 + γ_7 > 3/2
  have h1 : gammaQ 3 > (1 : ℚ) / 2 := gamma_3_active
  have h2 : gammaQ 5 > (1 : ℚ) / 2 := gamma_5_active
  have h3 : gammaQ 7 > (1 : ℚ) / 2 := gamma_7_active
  have hsum :
      gammaQ 3 + gammaQ 5 + gammaQ 7 > 1 / 2 + 1 / 2 + 1 / 2 := by
    have := add_lt_add (add_lt_add h1 h2) h3
    linarith
  linarith

/-- **Upper bound from `γ_p < 1` cap.** `Σ_active < 3`.

    Direct sum of the three strict inequalities `γ_p < 1`. -/
theorem gammaSumActive_lt_three : gammaSumActive < 3 := by
  unfold gammaSumActive
  have h1 : gammaQ 3 < 1 := gammaQ_3_lt_one
  have h2 : gammaQ 5 < 1 := gammaQ_5_lt_one
  have h3 : gammaQ 7 < 1 := gammaQ_7_lt_one
  linarith

/-- **Strict comparison with `1/s = 2`.** `Σ_active > 2`.

    Since `s = 1/2`, the value `1/s = 2` is a natural threshold; the
    cascade sits strictly above it. Follows from the bracket. -/
theorem gammaSumActive_gt_two : gammaSumActive > 2 :=
  lt_trans (by norm_num : (2 : ℚ) < 2099 / 1000) gammaSumActive_bracket.1

/-- **Strict comparison with `1/s = 2`** (named form, `1/sPT` on the RHS). -/
theorem gammaSumActive_gt_inv_s : gammaSumActive > 1 / sPT := by
  have h : (1 : ℚ) / sPT = 2 := by
    unfold sPT; norm_num
  rw [h]
  exact gammaSumActive_gt_two

/-! ### Link with the fixed point `μ* = 15` -/

/-- **Fixed-point invariant `Σ_active · μ*`.**

    Tight bracket `31.49 < (γ_3 + γ_5 + γ_7) · 15 < 31.50`.
    Numerical value `≈ 31.4910`. -/
theorem gammaSumActive_mul_muStar_bracket :
    3149 / 100 < gammaSumActive * muStar
    ∧ gammaSumActive * muStar < 3150 / 100 := by
  unfold gammaSumActive gammaQ deltaQ qPT muStar
  refine ⟨?_, ?_⟩ <;> norm_num

/-- **Positivity of `Σ_active · μ*`.** -/
theorem gammaSumActive_mul_muStar_pos : 0 < gammaSumActive * muStar := by
  have h1 : 0 < gammaSumActive := gammaSumActive_pos
  have h2 : (0 : ℚ) < muStar := by unfold muStar; norm_num
  exact mul_pos h1 h2

/-! ### Extension to the 5-prime tail `{3, 5, 7, 11, 13}` -/

/-- The 5-prime extended sum `γ_3 + γ_5 + γ_7 + γ_11 + γ_13`. -/
def gammaSumExtended : ℚ :=
  gammaQ 3 + gammaQ 5 + gammaQ 7 + gammaQ 11 + gammaQ 13

/-- The extended sum factors through the active sum. -/
theorem gammaSumExtended_eq :
    gammaSumExtended = gammaSumActive + gammaQ 11 + gammaQ 13 := by
  unfold gammaSumExtended gammaSumActive
  ring

/-- **Positivity of the extended sum.** -/
theorem gammaSumExtended_pos : 0 < gammaSumExtended := by
  unfold gammaSumExtended
  exact add_pos (add_pos (add_pos (add_pos
    gammaQ_3_pos gammaQ_5_pos) gammaQ_7_pos) gammaQ_11_pos) gammaQ_13_pos

/-- **Tight bracket for the extended sum.**
    `2.881 < γ_3 + γ_5 + γ_7 + γ_11 + γ_13 < 2.882`.
    Numerical value `≈ 2.8814`. -/
theorem gammaSumExtended_bracket :
    2881 / 1000 < gammaSumExtended ∧ gammaSumExtended < 2882 / 1000 := by
  unfold gammaSumExtended gammaQ deltaQ qPT muStar
  refine ⟨?_, ?_⟩ <;> norm_num

/-- **Extended > active.** Adding two positive terms `γ_11, γ_13 > 0`
    strictly increases the sum. -/
theorem gammaSumExtended_gt_active :
    gammaSumExtended > gammaSumActive := by
  rw [gammaSumExtended_eq]
  -- gammaSumActive + γ_11 + γ_13 > gammaSumActive
  have h11 : 0 < gammaQ 11 := gammaQ_11_pos
  have h13 : 0 < gammaQ 13 := gammaQ_13_pos
  linarith

/-- **Extended still bounded above by `5`.** Each of the five `γ_p < 1`. -/
theorem gammaSumExtended_lt_five : gammaSumExtended < 5 := by
  unfold gammaSumExtended
  have h1 : gammaQ 3 < 1 := gammaQ_3_lt_one
  have h2 : gammaQ 5 < 1 := gammaQ_5_lt_one
  have h3 : gammaQ 7 < 1 := gammaQ_7_lt_one
  have h4 : gammaQ 11 < 1 := gammaQ_11_lt_one
  have h5 : gammaQ 13 < 1 := gammaQ_13_lt_one
  linarith

/-! ### Headline -/

/-- **Headline (γ-sum cascade).** At `q_+ = 13/15`, `μ* = 15`, `s = 1/2`:

    1. `Σ_active = γ_3 + γ_5 + γ_7` is positive and lies in the tight
       bracket `(2.099, 2.100)`.
    2. `Σ_active > 3·s = 3/2` (active criterion lifts the sum above `3/2`).
    3. `Σ_active < 3` (the `γ_p < 1` cap keeps the sum below `3`).
    4. `Σ_active > 1/s = 2` (strict separation from the inverse-symmetry
       threshold).
    5. `Σ_active · μ*` lies in `(31.49, 31.50)` — the "γ-sum × fixed-point"
       invariant of the active cascade.
    6. The 5-prime extended sum `Σ_ext = Σ_active + γ_11 + γ_13` lies in
       `(2.881, 2.882)`, strictly exceeds `Σ_active`, and still stays
       below `5`. -/
theorem gamma_sum_summary :
    0 < gammaSumActive
    ∧ 2099 / 1000 < gammaSumActive ∧ gammaSumActive < 2100 / 1000
    ∧ gammaSumActive > 3 * sPT
    ∧ gammaSumActive < 3
    ∧ gammaSumActive > 1 / sPT
    ∧ 3149 / 100 < gammaSumActive * muStar
    ∧ gammaSumActive * muStar < 3150 / 100
    ∧ 2881 / 1000 < gammaSumExtended
    ∧ gammaSumExtended < 2882 / 1000
    ∧ gammaSumExtended > gammaSumActive
    ∧ gammaSumExtended < 5 :=
  ⟨gammaSumActive_pos,
   gammaSumActive_bracket.1,
   gammaSumActive_bracket.2,
   gammaSumActive_gt_three_s,
   gammaSumActive_lt_three,
   gammaSumActive_gt_inv_s,
   gammaSumActive_mul_muStar_bracket.1,
   gammaSumActive_mul_muStar_bracket.2,
   gammaSumExtended_bracket.1,
   gammaSumExtended_bracket.2,
   gammaSumExtended_gt_active,
   gammaSumExtended_lt_five⟩

end PT.Holonomy
