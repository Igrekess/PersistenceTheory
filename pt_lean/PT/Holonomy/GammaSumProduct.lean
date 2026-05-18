/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Holonomy.ActivePrimeCriterion
import PT.Holonomy.ActivePrimeMargins
import PT.Holonomy.GammaTablesExtended
import PT.Holonomy.GammaProduct
import PT.Holonomy.GammaSumActive
import Mathlib.Tactic

/-!
# Combined sum–product invariant `Σ_active · Π_active` over `{3, 5, 7}`

**Statement (paper-level, Ch06 §"Cascade arithmétique" / Ch09 supplement).**
Where `PT.Holonomy.GammaProduct` records the *product*
`Π_active = γ_3 γ_5 γ_7` and `PT.Holonomy.GammaSumActive` records the
*sum* `Σ_active = γ_3 + γ_5 + γ_7`, this module records the
**combined invariant**

$$\Theta_{\rm active} \;:=\; \Sigma_{\rm active} \cdot \Pi_{\rm active}
   \;=\; (\gamma_3 + \gamma_5 + \gamma_7) \cdot (\gamma_3 \gamma_5 \gamma_7)
   \;\approx\; 2.0994 \times 0.3348 \;\approx\; 0.7030.$$

The dual additive companion `Σ / Π ≈ 6.27` measures the spread of the
cascade in units of its product, and the AM-GM-type inequality

$$(\Sigma_{\rm active} / 3)^{3} \;>\; \Pi_{\rm active}$$

records the structural strict inequality between arithmetic and geometric
mean — strict because the three `γ_p` are unequal.

This module records:

1. The **definition** `gammaSumProductActive = Σ_active · Π_active`.
2. **Positivity** `0 < Θ_active`.
3. **Tight bracket** `0.7030 < Θ_active < 0.7031`.
4. **Structural lower bound** `Θ_active > (3/2) · Π_active`, which is
   `Σ_active > 3s = 3/2` multiplied by the positive `Π_active`.
5. **Ratio** `Σ / Π` with tight bracket `(6.269, 6.270)`.
6. **AM-GM strict inequality** `(Σ/3)^3 > Π` in rational form, verified
   numerically via the bracket `0.007 < (Σ/3)^3 - Π < 0.008`.
7. **Headline** combining all invariants.

All proofs are pure rational arithmetic — `unfold` + `norm_num`.

## Reference

* Monograph Chapter 6, §"Cascade arithmétique", table `tab:gamma_invariants`.
* Companion files `PT.Holonomy.GammaProduct`, `PT.Holonomy.GammaSumActive`.
-/

namespace PT.Holonomy

/-! ### Definition: combined sum × product over the active cascade -/

/-- The PT combined γ-invariant over the active primes `{3, 5, 7}` at
    `q_+ = 13/15`, `μ* = 15`:

    `gammaSumProductActive := Σ_active · Π_active`. -/
def gammaSumProductActive : ℚ := gammaSumActive * gammaProductActive

/-- Definitional unfolding: the combined invariant is literally the sum
    times the product. -/
theorem gammaSumProductActive_eq :
    gammaSumProductActive = gammaSumActive * gammaProductActive := rfl

/-! ### Positivity -/

/-- **Positivity** `0 < Θ_active = Σ · Π`.

    Both factors are positive (active criterion). -/
theorem gammaSumProductActive_pos : 0 < gammaSumProductActive := by
  unfold gammaSumProductActive
  exact mul_pos gammaSumActive_pos gammaProductActive_pos

/-! ### Tight rational bracket -/

/-- **Tight bracket for `Θ_active = Σ · Π`.**
    `0.7030 < Σ_active · Π_active < 0.7031`. Numerical value `≈ 0.70301`. -/
theorem gammaSumProductActive_bracket :
    7030 / 10000 < gammaSumProductActive
    ∧ gammaSumProductActive < 7031 / 10000 := by
  unfold gammaSumProductActive gammaSumActive gammaProductActive
        gammaQ deltaQ qPT muStar
  refine ⟨?_, ?_⟩ <;> norm_num

/-! ### Structural lower bound from `Σ > 3·s = 3/2` -/

/-- **Structural lower bound.** `Θ_active > (3/2) · Π_active`.

    Direct consequence of `Σ_active > 3·s = 3/2`
    (`gammaSumActive_gt_three_s`) multiplied by the positive
    `Π_active`. Numerically `≈ 0.7030 > 0.5023`. -/
theorem gammaSumProductActive_gt_threeHalves_prod :
    gammaSumProductActive > (3 / 2) * gammaProductActive := by
  unfold gammaSumProductActive
  have hsum : gammaSumActive > 3 * sPT := gammaSumActive_gt_three_s
  have hsum' : gammaSumActive > 3 / 2 := by
    have h : (3 : ℚ) * sPT = 3 / 2 := by unfold sPT; norm_num
    linarith
  have hpos : 0 < gammaProductActive := gammaProductActive_pos
  -- Σ * Π > (3/2) * Π
  have := mul_lt_mul_of_pos_right hsum' hpos
  linarith

/-- **Structural lower bound** in compact `3 · s · Π` form.
    `Θ_active > 3 · s · Π`. -/
theorem gammaSumProductActive_gt_three_s_prod :
    gammaSumProductActive > 3 * sPT * gammaProductActive := by
  have h : (3 : ℚ) * sPT * gammaProductActive
            = (3 / 2) * gammaProductActive := by
    unfold sPT; ring
  rw [h]
  exact gammaSumProductActive_gt_threeHalves_prod

/-! ### Upper bound from `Σ < 3` -/

/-- **Upper bound from the cap `Σ_active < 3`.**
    `Θ_active < 3 · Π_active`.

    Direct consequence of `γ_p < 1` (so `Σ < 3`) multiplied by the
    positive `Π_active`. -/
theorem gammaSumProductActive_lt_three_prod :
    gammaSumProductActive < 3 * gammaProductActive := by
  unfold gammaSumProductActive
  have hsum : gammaSumActive < 3 := gammaSumActive_lt_three
  have hpos : 0 < gammaProductActive := gammaProductActive_pos
  -- Σ * Π < 3 * Π
  have := mul_lt_mul_of_pos_right hsum hpos
  linarith

/-! ### The ratio `Σ / Π` -/

/-- The PT γ-ratio over the active primes — spread of the cascade in
    units of its product:

    `gammaSumOverProductActive := Σ_active / Π_active`. -/
def gammaSumOverProductActive : ℚ := gammaSumActive / gammaProductActive

/-- **Positivity of the ratio.** `0 < Σ_active / Π_active`. -/
theorem gammaSumOverProductActive_pos : 0 < gammaSumOverProductActive := by
  unfold gammaSumOverProductActive
  exact div_pos gammaSumActive_pos gammaProductActive_pos

/-- **Tight bracket for the ratio.**
    `6.269 < Σ_active / Π_active < 6.270`. Numerical value `≈ 6.2698`. -/
theorem gammaSumOverProductActive_bracket :
    6269 / 1000 < gammaSumOverProductActive
    ∧ gammaSumOverProductActive < 6270 / 1000 := by
  unfold gammaSumOverProductActive gammaSumActive gammaProductActive
        gammaQ deltaQ qPT muStar
  refine ⟨?_, ?_⟩ <;> norm_num

/-- **Ratio exceeds `6`.** `Σ_active / Π_active > 6`.

    Follows from the tight bracket; sharp consequence of the cascade's
    arithmetic spread (each γ_p stays close to `s = 1/2` while their
    product collapses). -/
theorem gammaSumOverProductActive_gt_six :
    gammaSumOverProductActive > 6 :=
  lt_trans (by norm_num : (6 : ℚ) < 6269 / 1000)
           gammaSumOverProductActive_bracket.1

/-! ### AM-GM strict inequality `(Σ/3)^3 > Π` -/

/-- **AM-GM strict inequality, rational form.**
    `((Σ_active / 3))^3 > Π_active`.

    Direct rational reformulation of the AM-GM inequality
    `(γ_3 + γ_5 + γ_7) / 3 ≥ (γ_3 γ_5 γ_7)^{1/3}` without invoking
    real cube roots: cubing both sides preserves the inequality since
    both sides are positive, and the inequality is **strict** because
    the three `γ_p` are pairwise distinct. -/
theorem gammaSumActive_cubed_gt_27_prod :
    (gammaSumActive / 3)^3 > gammaProductActive := by
  unfold gammaSumActive gammaProductActive gammaQ deltaQ qPT muStar
  norm_num

/-- **AM-GM gap bracket.**
    `0.007 < (Σ_active / 3)^3 - Π_active < 0.008`.

    Records the magnitude of the AM-GM gap; numerical value `≈ 0.00787`. -/
theorem gammaSumActive_amgm_gap_bracket :
    7 / 1000 < (gammaSumActive / 3)^3 - gammaProductActive
    ∧ (gammaSumActive / 3)^3 - gammaProductActive < 8 / 1000 := by
  unfold gammaSumActive gammaProductActive gammaQ deltaQ qPT muStar
  refine ⟨?_, ?_⟩ <;> norm_num

/-- **`(Σ/3)^3` tight bracket.**
    `0.3427 < (Σ_active / 3)^3 < 0.3428`. -/
theorem gammaSumActive_cubed_bracket :
    3427 / 10000 < (gammaSumActive / 3)^3
    ∧ (gammaSumActive / 3)^3 < 3428 / 10000 := by
  unfold gammaSumActive gammaQ deltaQ qPT muStar
  refine ⟨?_, ?_⟩ <;> norm_num

/-! ### Headline -/

/-- **Headline (combined sum × product γ-invariant).**
    At `q_+ = 13/15`, `μ* = 15`, `s = 1/2`:

    1. `Θ_active = Σ_active · Π_active` is positive and lies in the tight
       bracket `(0.7030, 0.7031)`.
    2. Structural lower bound: `Θ_active > (3/2) · Π_active = 3·s·Π`
       (from `Σ > 3s`).
    3. Structural upper bound: `Θ_active < 3 · Π_active` (from `Σ < 3`).
    4. The ratio `Σ_active / Π_active` lies in the tight bracket
       `(6.269, 6.270)` and strictly exceeds `6`.
    5. AM-GM strict inequality: `((Σ/3))^3 > Π`, with explicit gap
       bracket `0.007 < (Σ/3)^3 - Π < 0.008`. -/
theorem gamma_sum_product_summary :
    0 < gammaSumProductActive
    ∧ 7030 / 10000 < gammaSumProductActive
    ∧ gammaSumProductActive < 7031 / 10000
    ∧ gammaSumProductActive > (3 / 2) * gammaProductActive
    ∧ gammaSumProductActive > 3 * sPT * gammaProductActive
    ∧ gammaSumProductActive < 3 * gammaProductActive
    ∧ 0 < gammaSumOverProductActive
    ∧ 6269 / 1000 < gammaSumOverProductActive
    ∧ gammaSumOverProductActive < 6270 / 1000
    ∧ gammaSumOverProductActive > 6
    ∧ (gammaSumActive / 3)^3 > gammaProductActive
    ∧ 7 / 1000 < (gammaSumActive / 3)^3 - gammaProductActive
    ∧ (gammaSumActive / 3)^3 - gammaProductActive < 8 / 1000 :=
  ⟨gammaSumProductActive_pos,
   gammaSumProductActive_bracket.1,
   gammaSumProductActive_bracket.2,
   gammaSumProductActive_gt_threeHalves_prod,
   gammaSumProductActive_gt_three_s_prod,
   gammaSumProductActive_lt_three_prod,
   gammaSumOverProductActive_pos,
   gammaSumOverProductActive_bracket.1,
   gammaSumOverProductActive_bracket.2,
   gammaSumOverProductActive_gt_six,
   gammaSumActive_cubed_gt_27_prod,
   gammaSumActive_amgm_gap_bracket.1,
   gammaSumActive_amgm_gap_bracket.2⟩

end PT.Holonomy
