/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Holonomy.AlphaPowerSequence
import PT.Holonomy.GammaSumActive
import PT.Holonomy.InvAlphaSquaredBracket
import Mathlib.Tactic

/-!
# Mixed invariant `α_bare · Σ γ_p` and its power-tower `α_bare^n · Σ γ_p`

This module records the **mixed holonomy invariant** combining the
*bare* coupling `α_bare = ∏_{p ∈ {3,5,7}} sin²θ_p` (cf.
`PT.Holonomy.CouplingReconstruction`) with the *additive* γ-sum
`Σ_active = γ_3 + γ_5 + γ_7` (cf. `PT.Holonomy.GammaSumActive`).

This product appears in PT cascade formulae whenever a multiplicative
coupling factor multiplies an additive holonomy contribution — typically
in single-loop perturbative corrections where each loop contributes a
factor `α_bare` and the holonomy weight is the γ-sum.

Numerical headline (exact rationals at `q_+ = 13/15`, `μ⋆ = 15`):

* `α_bare · Σγ                 ≈ 1.5405 × 10⁻²`     — *headline*.
* `α_bare^2 · Σγ               ≈ 1.1304 × 10⁻⁴`.
* `α_bare^3 · Σγ               ≈ 8.2950 × 10⁻⁷`.
* `(α_bare · Σγ) / (α_bare⁻¹)² = α_bare^3 · Σγ`     — algebraic identity.

The tower decreases strictly geometrically (ratio `α_bare ≈ 7.34 × 10⁻³`).

All brackets are decided by `norm_num` after unfolding the rational
definitions. Decay is purely arithmetic.

## Reference

* `PT.Holonomy.CouplingReconstruction.alphaBareQ` — multiplicative factor.
* `PT.Holonomy.GammaSumActive.gammaSumActive` — additive factor.
* `PT.Holonomy.AlphaPowerSequence.alphaPow` — power tower.
* `PT.Holonomy.InvAlphaSquaredBracket.invAlphaSquared` — square inverse.
-/

namespace PT.Holonomy

/-! ### Definition: `α_bare · Σ γ_p` -/

/-- The **mixed holonomy invariant** `α_bare · (γ_3 + γ_5 + γ_7)`. -/
def alphaTimesGammaSum : ℚ := alphaBareQ * gammaSumActive

/-- Definitional unfolding. -/
theorem alphaTimesGammaSum_eq :
    alphaTimesGammaSum = alphaBareQ * gammaSumActive := rfl

/-! ### Positivity and bracket -/

/-- **Positivity.** `0 < α_bare · Σγ`. -/
theorem alphaTimesGammaSum_pos : 0 < alphaTimesGammaSum := by
  unfold alphaTimesGammaSum
  exact mul_pos alphaBareQ_pos gammaSumActive_pos

/-- **Decimal bracket.** `1.540 × 10⁻² < α_bare · Σγ < 1.541 × 10⁻²`
    (exact value `≈ 1.5405 × 10⁻²`). -/
theorem alphaTimesGammaSum_bracket :
    1540 / 100000 < alphaTimesGammaSum
    ∧ alphaTimesGammaSum < 1541 / 100000 := by
  simp only [alphaTimesGammaSum, alphaBareQ, sinSqQ,
    gammaSumActive, gammaQ, deltaQ, qPT, muStar]
  refine ⟨?_, ?_⟩ <;> norm_num

/-- **Tight bracket.** `15405/10⁶ < α_bare · Σγ < 15406/10⁶`. -/
theorem alphaTimesGammaSum_bracket_tight :
    15405 / 1000000 < alphaTimesGammaSum
    ∧ alphaTimesGammaSum < 15406 / 1000000 := by
  simp only [alphaTimesGammaSum, alphaBareQ, sinSqQ,
    gammaSumActive, gammaQ, deltaQ, qPT, muStar]
  refine ⟨?_, ?_⟩ <;> norm_num

/-- The mixed invariant is bounded above by `1/50` (a clean loose cap). -/
theorem alphaTimesGammaSum_lt_one_fiftieth :
    alphaTimesGammaSum < 1 / 50 := by
  have h := alphaTimesGammaSum_bracket.2
  linarith

/-! ### Power tower: `α_bare^n · Σ γ_p` -/

/-- The **`n`-th power mixed invariant** `α_bare^n · Σγ`. -/
def alphaPowTimesGammaSum (n : ℕ) : ℚ := alphaPow n * gammaSumActive

/-- Definitional unfolding. -/
theorem alphaPowTimesGammaSum_eq (n : ℕ) :
    alphaPowTimesGammaSum n = alphaPow n * gammaSumActive := rfl

/-- At `n = 1`, the tower reduces to the headline mixed invariant. -/
theorem alphaPowTimesGammaSum_one :
    alphaPowTimesGammaSum 1 = alphaTimesGammaSum := by
  unfold alphaPowTimesGammaSum alphaTimesGammaSum
  rw [alphaPow_one]

/-- At `n = 0`, the tower is just the γ-sum. -/
theorem alphaPowTimesGammaSum_zero :
    alphaPowTimesGammaSum 0 = gammaSumActive := by
  unfold alphaPowTimesGammaSum
  rw [alphaPow_zero]; ring

/-! ### Generic positivity -/

/-- **Generic positivity.** `0 < α_bare^n · Σγ` for all `n`. -/
theorem alphaPowTimesGammaSum_pos (n : ℕ) : 0 < alphaPowTimesGammaSum n := by
  unfold alphaPowTimesGammaSum
  exact mul_pos (alphaPow_pos n) gammaSumActive_pos

theorem alphaPowTimesGammaSum_one_pos   : 0 < alphaPowTimesGammaSum 1 :=
  alphaPowTimesGammaSum_pos 1
theorem alphaPowTimesGammaSum_two_pos   : 0 < alphaPowTimesGammaSum 2 :=
  alphaPowTimesGammaSum_pos 2
theorem alphaPowTimesGammaSum_three_pos : 0 < alphaPowTimesGammaSum 3 :=
  alphaPowTimesGammaSum_pos 3

/-! ### Tight brackets for `n = 1, 2, 3` -/

/-- **Bracket `n = 1`.** `15405/10⁶ < α_bare · Σγ < 15406/10⁶`
    (exact `≈ 1.5405 × 10⁻²`). -/
theorem alphaPowTimesGammaSum_one_bracket :
    15405 / 1000000 < alphaPowTimesGammaSum 1
    ∧ alphaPowTimesGammaSum 1 < 15406 / 1000000 := by
  rw [alphaPowTimesGammaSum_one]
  exact alphaTimesGammaSum_bracket_tight

/-- **Bracket `n = 2`.** `11304/10⁸ < α_bare² · Σγ < 11305/10⁸`
    (exact `≈ 1.1304 × 10⁻⁴`). -/
theorem alphaPowTimesGammaSum_two_bracket :
    11304 / 100000000 < alphaPowTimesGammaSum 2
    ∧ alphaPowTimesGammaSum 2 < 11305 / 100000000 := by
  simp only [alphaPowTimesGammaSum, alphaPow, alphaBareQ, sinSqQ,
    gammaSumActive, gammaQ, deltaQ, qPT, muStar]
  refine ⟨?_, ?_⟩ <;> norm_num

/-- **Bracket `n = 3`.** `82949/10¹¹ < α_bare³ · Σγ < 82950/10¹¹`
    (exact `≈ 8.2950 × 10⁻⁷`). -/
theorem alphaPowTimesGammaSum_three_bracket :
    82949 / 100000000000 < alphaPowTimesGammaSum 3
    ∧ alphaPowTimesGammaSum 3 < 82950 / 100000000000 := by
  simp only [alphaPowTimesGammaSum, alphaPow, alphaBareQ, sinSqQ,
    gammaSumActive, gammaQ, deltaQ, qPT, muStar]
  refine ⟨?_, ?_⟩ <;> norm_num

/-! ### Strict geometric decrease -/

/-- **Strict geometric decrease.** `α^(n+1) · Σγ < α^n · Σγ` for all `n`,
    since `0 < α < 1` and `Σγ > 0`. -/
theorem alphaPowTimesGammaSum_succ_lt (n : ℕ) :
    alphaPowTimesGammaSum (n + 1) < alphaPowTimesGammaSum n := by
  unfold alphaPowTimesGammaSum
  have hdec : alphaPow (n + 1) < alphaPow n := alphaPow_succ_lt n
  have hpos : 0 < gammaSumActive := gammaSumActive_pos
  exact mul_lt_mul_of_pos_right hdec hpos

theorem alphaPowTimesGammaSum_dec_1_2 :
    alphaPowTimesGammaSum 2 < alphaPowTimesGammaSum 1 :=
  alphaPowTimesGammaSum_succ_lt 1

theorem alphaPowTimesGammaSum_dec_2_3 :
    alphaPowTimesGammaSum 3 < alphaPowTimesGammaSum 2 :=
  alphaPowTimesGammaSum_succ_lt 2

/-- **Strict cascade** for `n = 1, 2, 3`. -/
theorem alphaPowTimesGammaSum_chain_1_to_3 :
    alphaPowTimesGammaSum 3 < alphaPowTimesGammaSum 2
    ∧ alphaPowTimesGammaSum 2 < alphaPowTimesGammaSum 1 :=
  ⟨alphaPowTimesGammaSum_dec_2_3, alphaPowTimesGammaSum_dec_1_2⟩

/-! ### Link with `invAlphaSquared`

The square inverse `(α_bare⁻¹)² = 1/α_bare²` interacts cleanly with the
mixed invariant: dividing `α_bare · Σγ` by `(α_bare⁻¹)²` collapses to
`α_bare^3 · Σγ`. -/

/-- **Algebraic identity.**
    `(α_bare · Σγ) / (α_bare⁻¹)² = α_bare³ · Σγ`. -/
theorem alphaTimesGammaSum_div_invAlphaSq :
    alphaTimesGammaSum / invAlphaSquared = alphaPowTimesGammaSum 3 := by
  unfold alphaTimesGammaSum alphaPowTimesGammaSum alphaPow
  rw [invAlphaSquared_eq_inv_alpha_sq]
  have ha : alphaBareQ ≠ 0 := ne_of_gt alphaBareQ_pos
  field_simp

/-- **Positivity of the ratio.** -/
theorem alphaTimesGammaSum_div_invAlphaSq_pos :
    0 < alphaTimesGammaSum / invAlphaSquared := by
  rw [alphaTimesGammaSum_div_invAlphaSq]
  exact alphaPowTimesGammaSum_three_pos

/-- **Bracket for the ratio.** Same as the `n = 3` bracket. -/
theorem alphaTimesGammaSum_div_invAlphaSq_bracket :
    82949 / 100000000000 < alphaTimesGammaSum / invAlphaSquared
    ∧ alphaTimesGammaSum / invAlphaSquared < 82950 / 100000000000 := by
  rw [alphaTimesGammaSum_div_invAlphaSq]
  exact alphaPowTimesGammaSum_three_bracket

/-! ### Headline summary -/

/-- **Headline.** All invariants at a glance for `α_bare · Σγ` and its
    power tower at the PT fixed point `q_+ = 13/15`, `μ⋆ = 15`:

    * Positivity: `0 < α_bare^n · Σγ` for `n = 1, 2, 3`.
    * Tight rational brackets:
        * `α_bare · Σγ            ∈ (15405/10⁶, 15406/10⁶)`  ≈ 1.5405 × 10⁻²
        * `α_bare² · Σγ           ∈ (11304/10⁸, 11305/10⁸)`  ≈ 1.1304 × 10⁻⁴
        * `α_bare³ · Σγ           ∈ (82949/10¹¹, 82950/10¹¹)` ≈ 8.2950 × 10⁻⁷
    * Strict geometric decrease: `α_bare³ · Σγ < α_bare² · Σγ < α_bare · Σγ`.
    * Algebraic identity: `(α_bare · Σγ) / (α_bare⁻¹)² = α_bare³ · Σγ`.
    * Loose cap: `α_bare · Σγ < 1/50`. -/
theorem alphaTimesGammaSum_summary :
    -- positivity
    0 < alphaPowTimesGammaSum 1
    ∧ 0 < alphaPowTimesGammaSum 2
    ∧ 0 < alphaPowTimesGammaSum 3
    -- n = 1 bracket
    ∧ 15405 / 1000000 < alphaPowTimesGammaSum 1
    ∧ alphaPowTimesGammaSum 1 < 15406 / 1000000
    -- n = 2 bracket
    ∧ 11304 / 100000000 < alphaPowTimesGammaSum 2
    ∧ alphaPowTimesGammaSum 2 < 11305 / 100000000
    -- n = 3 bracket
    ∧ 82949 / 100000000000 < alphaPowTimesGammaSum 3
    ∧ alphaPowTimesGammaSum 3 < 82950 / 100000000000
    -- strict cascade
    ∧ alphaPowTimesGammaSum 3 < alphaPowTimesGammaSum 2
    ∧ alphaPowTimesGammaSum 2 < alphaPowTimesGammaSum 1
    -- algebraic link with invAlphaSquared
    ∧ alphaTimesGammaSum / invAlphaSquared = alphaPowTimesGammaSum 3
    -- loose cap
    ∧ alphaTimesGammaSum < 1 / 50 :=
  ⟨alphaPowTimesGammaSum_one_pos,
   alphaPowTimesGammaSum_two_pos,
   alphaPowTimesGammaSum_three_pos,
   alphaPowTimesGammaSum_one_bracket.1, alphaPowTimesGammaSum_one_bracket.2,
   alphaPowTimesGammaSum_two_bracket.1, alphaPowTimesGammaSum_two_bracket.2,
   alphaPowTimesGammaSum_three_bracket.1, alphaPowTimesGammaSum_three_bracket.2,
   alphaPowTimesGammaSum_dec_2_3, alphaPowTimesGammaSum_dec_1_2,
   alphaTimesGammaSum_div_invAlphaSq,
   alphaTimesGammaSum_lt_one_fiftieth⟩

end PT.Holonomy
