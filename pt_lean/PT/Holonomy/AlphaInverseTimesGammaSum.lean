/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Holonomy.AlphaInversePowerSequence
import PT.Holonomy.AlphaTimesGammaSum
import PT.Holonomy.CouplingReconstructionBounds
import PT.Holonomy.GammaSumActive
import Mathlib.Tactic

/-!
# Mixed invariant `α_bare⁻¹ · Σ γ_p` and its inverse power-tower

This module is the **inverse-coupling dual** of
`PT.Holonomy.AlphaTimesGammaSum`. It records the mixed invariant
combining the *inverse* bare coupling
`α_bare⁻¹ = 1 / (sin²θ_3 · sin²θ_5 · sin²θ_7)` (cf.
`PT.Holonomy.CouplingReconstructionBounds`) with the *additive* γ-sum
`Σ_active = γ_3 + γ_5 + γ_7` (cf. `PT.Holonomy.GammaSumActive`).

This product appears in PT cascade formulae whenever an *inverse*
coupling factor multiplies an additive holonomy contribution — typically
in the experimentally observable form `α_EM⁻¹ ≈ 137`, where each loop
order amplifies the γ-sum by a factor `α_bare⁻¹` rather than suppressing
it.

Numerical headline (exact rationals at `q_+ = 13/15`, `μ⋆ = 15`):

* `α_bare⁻¹ · Σγ                ≈ 2.8610 × 10²`     — *headline*.
* `(α_bare⁻¹)² · Σγ             ≈ 3.8990 × 10⁴`.
* `(α_bare⁻¹)³ · Σγ             ≈ 5.3134 × 10⁶`.
* `(α · Σγ) · (α⁻¹ · Σγ) = (Σγ)²`               — duality identity.
* `(α · Σγ) · (α⁻¹ · Σγ) / (Σγ)² = 1`           — normalised duality.

The inverse tower **increases** strictly (ratio `α_bare⁻¹ ≈ 136.28 > 1`)
in contrast with the direct `α^n · Σγ` tower which decreases.

All brackets are decided by `unfold + norm_num` after unfolding the
rational definitions.

## Reference

* `PT.Holonomy.CouplingReconstructionBounds.alphaBareInvQ` — inverse factor.
* `PT.Holonomy.GammaSumActive.gammaSumActive` — additive factor.
* `PT.Holonomy.AlphaInversePowerSequence.alphaInvPow` — inverse power tower.
* `PT.Holonomy.AlphaTimesGammaSum` — direct dual.
-/

namespace PT.Holonomy

/-! ### Definition: `α_bare⁻¹ · Σ γ_p` -/

/-- The **mixed inverse-coupling holonomy invariant**
    `α_bare⁻¹ · (γ_3 + γ_5 + γ_7)`. -/
noncomputable def alphaInvTimesGammaSum : ℚ := alphaBareInvQ * gammaSumActive

/-- Definitional unfolding. -/
theorem alphaInvTimesGammaSum_eq :
    alphaInvTimesGammaSum = alphaBareInvQ * gammaSumActive := rfl

/-! ### Positivity and bracket -/

/-- **Positivity.** `0 < α_bare⁻¹ · Σγ`. -/
theorem alphaInvTimesGammaSum_pos : 0 < alphaInvTimesGammaSum := by
  unfold alphaInvTimesGammaSum
  exact mul_pos alphaBareInvQ_pos gammaSumActive_pos

/-- **Decimal bracket.** `286 < α_bare⁻¹ · Σγ < 287`
    (exact value `≈ 286.10`). -/
theorem alphaInvTimesGammaSum_bracket :
    (286 : ℚ) < alphaInvTimesGammaSum
    ∧ alphaInvTimesGammaSum < 287 := by
  unfold alphaInvTimesGammaSum alphaBareInvQ alphaBareQ sinSqQ
    gammaSumActive gammaQ deltaQ qPT muStar
  refine ⟨?_, ?_⟩ <;> norm_num

/-- **Tight bracket.** `28610/100 < α_bare⁻¹ · Σγ < 28611/100`
    (exact `≈ 286.1025`). -/
theorem alphaInvTimesGammaSum_bracket_tight :
    28610 / 100 < alphaInvTimesGammaSum
    ∧ alphaInvTimesGammaSum < 28611 / 100 := by
  unfold alphaInvTimesGammaSum alphaBareInvQ alphaBareQ sinSqQ
    gammaSumActive gammaQ deltaQ qPT muStar
  refine ⟨?_, ?_⟩ <;> norm_num

/-- The mixed inverse invariant is bounded below by `285` (a clean cap). -/
theorem alphaInvTimesGammaSum_gt_285 :
    (285 : ℚ) < alphaInvTimesGammaSum := by
  have h := alphaInvTimesGammaSum_bracket.1
  linarith

/-! ### Power tower: `(α_bare⁻¹)^n · Σ γ_p` -/

/-- The **`n`-th inverse-power mixed invariant** `(α_bare⁻¹)^n · Σγ`. -/
noncomputable def alphaInvPowTimesGammaSum (n : ℕ) : ℚ :=
  alphaInvPow n * gammaSumActive

/-- Definitional unfolding. -/
theorem alphaInvPowTimesGammaSum_eq (n : ℕ) :
    alphaInvPowTimesGammaSum n = alphaInvPow n * gammaSumActive := rfl

/-- At `n = 1`, the tower reduces to the headline mixed invariant. -/
theorem alphaInvPowTimesGammaSum_one :
    alphaInvPowTimesGammaSum 1 = alphaInvTimesGammaSum := by
  unfold alphaInvPowTimesGammaSum alphaInvTimesGammaSum
  rw [alphaInvPow_one]

/-- At `n = 0`, the tower is just the γ-sum. -/
theorem alphaInvPowTimesGammaSum_zero :
    alphaInvPowTimesGammaSum 0 = gammaSumActive := by
  unfold alphaInvPowTimesGammaSum
  rw [alphaInvPow_zero]; ring

/-! ### Generic positivity -/

/-- **Generic positivity.** `0 < (α⁻¹)^n · Σγ` for all `n`. -/
theorem alphaInvPowTimesGammaSum_pos (n : ℕ) :
    0 < alphaInvPowTimesGammaSum n := by
  unfold alphaInvPowTimesGammaSum
  exact mul_pos (alphaInvPow_pos n) gammaSumActive_pos

theorem alphaInvPowTimesGammaSum_one_pos   : 0 < alphaInvPowTimesGammaSum 1 :=
  alphaInvPowTimesGammaSum_pos 1
theorem alphaInvPowTimesGammaSum_two_pos   : 0 < alphaInvPowTimesGammaSum 2 :=
  alphaInvPowTimesGammaSum_pos 2
theorem alphaInvPowTimesGammaSum_three_pos : 0 < alphaInvPowTimesGammaSum 3 :=
  alphaInvPowTimesGammaSum_pos 3

/-! ### Tight brackets for `n = 1, 2, 3` -/

/-- **Bracket `n = 1`.** `28610/100 < (α⁻¹) · Σγ < 28611/100`
    (exact `≈ 286.1025`). -/
theorem alphaInvPowTimesGammaSum_one_bracket :
    28610 / 100 < alphaInvPowTimesGammaSum 1
    ∧ alphaInvPowTimesGammaSum 1 < 28611 / 100 := by
  rw [alphaInvPowTimesGammaSum_one]
  exact alphaInvTimesGammaSum_bracket_tight

/-- **Bracket `n = 2`.** `38989 < (α⁻¹)² · Σγ < 38990`
    (exact `≈ 38989.57`). -/
theorem alphaInvPowTimesGammaSum_two_bracket :
    (38989 : ℚ) < alphaInvPowTimesGammaSum 2
    ∧ alphaInvPowTimesGammaSum 2 < 38990 := by
  unfold alphaInvPowTimesGammaSum alphaInvPow alphaBareInvQ alphaBareQ
    sinSqQ gammaSumActive gammaQ deltaQ qPT muStar
  refine ⟨?_, ?_⟩ <;> norm_num

/-- **Bracket `n = 3`.** `5313434 < (α⁻¹)³ · Σγ < 5313435`
    (exact `≈ 5 313 434.12`). -/
theorem alphaInvPowTimesGammaSum_three_bracket :
    (5313434 : ℚ) < alphaInvPowTimesGammaSum 3
    ∧ alphaInvPowTimesGammaSum 3 < 5313435 := by
  unfold alphaInvPowTimesGammaSum alphaInvPow alphaBareInvQ alphaBareQ
    sinSqQ gammaSumActive gammaQ deltaQ qPT muStar
  refine ⟨?_, ?_⟩ <;> norm_num

/-! ### Strict geometric increase -/

/-- **Strict geometric increase.** `(α⁻¹)^n · Σγ < (α⁻¹)^(n+1) · Σγ`
    for all `n`, since `α⁻¹ > 1` and `Σγ > 0`. -/
theorem alphaInvPowTimesGammaSum_lt_succ (n : ℕ) :
    alphaInvPowTimesGammaSum n < alphaInvPowTimesGammaSum (n + 1) := by
  unfold alphaInvPowTimesGammaSum
  have hinc : alphaInvPow n < alphaInvPow (n + 1) := alphaInvPow_lt_succ n
  have hpos : 0 < gammaSumActive := gammaSumActive_pos
  exact mul_lt_mul_of_pos_right hinc hpos

theorem alphaInvPowTimesGammaSum_inc_1_2 :
    alphaInvPowTimesGammaSum 1 < alphaInvPowTimesGammaSum 2 :=
  alphaInvPowTimesGammaSum_lt_succ 1

theorem alphaInvPowTimesGammaSum_inc_2_3 :
    alphaInvPowTimesGammaSum 2 < alphaInvPowTimesGammaSum 3 :=
  alphaInvPowTimesGammaSum_lt_succ 2

/-- **Strict cascade** for `n = 1, 2, 3`. -/
theorem alphaInvPowTimesGammaSum_chain_1_to_3 :
    alphaInvPowTimesGammaSum 1 < alphaInvPowTimesGammaSum 2
    ∧ alphaInvPowTimesGammaSum 2 < alphaInvPowTimesGammaSum 3 :=
  ⟨alphaInvPowTimesGammaSum_inc_1_2, alphaInvPowTimesGammaSum_inc_2_3⟩

/-! ### Duality with `alphaTimesGammaSum`

The product of the direct invariant `α · Σγ` and the inverse invariant
`α⁻¹ · Σγ` collapses to `(Σγ)²`, since `α · α⁻¹ = 1`. This is the
arithmetic dual of `alphaPow_mul_alphaInvPow` evaluated at `n = 1`. -/

/-- **Duality identity.**
    `(α · Σγ) · (α⁻¹ · Σγ) = (Σγ)²`. -/
theorem alphaTimesGammaSum_mul_alphaInvTimesGammaSum :
    alphaTimesGammaSum * alphaInvTimesGammaSum = gammaSumActive ^ 2 := by
  unfold alphaTimesGammaSum alphaInvTimesGammaSum
  have h : alphaBareQ * alphaBareInvQ = 1 := alphaBareQ_mul_alphaBareInvQ
  -- Rearrange: (α · Σγ) · (α⁻¹ · Σγ) = (α · α⁻¹) · Σγ · Σγ = Σγ^2
  have hsq : gammaSumActive ^ 2 = gammaSumActive * gammaSumActive := sq gammaSumActive
  rw [hsq]
  calc
    alphaBareQ * gammaSumActive * (alphaBareInvQ * gammaSumActive)
        = (alphaBareQ * alphaBareInvQ) * (gammaSumActive * gammaSumActive) := by ring
    _   = 1 * (gammaSumActive * gammaSumActive) := by rw [h]
    _   = gammaSumActive * gammaSumActive := one_mul _

/-- **Normalised duality identity.**
    `(α · Σγ) · (α⁻¹ · Σγ) / (Σγ)² = 1`. -/
theorem alphaTimesGammaSum_mul_alphaInvTimesGammaSum_div_gammaSq :
    alphaTimesGammaSum * alphaInvTimesGammaSum / gammaSumActive ^ 2 = 1 := by
  rw [alphaTimesGammaSum_mul_alphaInvTimesGammaSum]
  have hne : gammaSumActive ^ 2 ≠ 0 := pow_ne_zero _ (ne_of_gt gammaSumActive_pos)
  exact div_self hne

/-! ### Headline summary -/

/-- **Headline.** All invariants at a glance for `α_bare⁻¹ · Σγ` and its
    inverse power tower at the PT fixed point `q_+ = 13/15`, `μ⋆ = 15`:

    * Positivity: `0 < (α⁻¹)^n · Σγ` for `n = 1, 2, 3`.
    * Tight rational brackets:
        * `α⁻¹ · Σγ        ∈ (28610/100, 28611/100)`     ≈ 286.10
        * `(α⁻¹)² · Σγ     ∈ (38989, 38990)`             ≈ 3.899 × 10⁴
        * `(α⁻¹)³ · Σγ     ∈ (5313434, 5313435)`         ≈ 5.313 × 10⁶
    * Strict geometric increase:
        `α⁻¹ · Σγ < (α⁻¹)² · Σγ < (α⁻¹)³ · Σγ`.
    * Duality identity: `(α · Σγ) · (α⁻¹ · Σγ) = (Σγ)²`.
    * Normalised duality: `(α · Σγ) · (α⁻¹ · Σγ) / (Σγ)² = 1`.
    * Loose cap: `α⁻¹ · Σγ > 285`. -/
theorem alphaInvTimesGammaSum_summary :
    -- positivity
    0 < alphaInvPowTimesGammaSum 1
    ∧ 0 < alphaInvPowTimesGammaSum 2
    ∧ 0 < alphaInvPowTimesGammaSum 3
    -- n = 1 bracket
    ∧ 28610 / 100 < alphaInvPowTimesGammaSum 1
    ∧ alphaInvPowTimesGammaSum 1 < 28611 / 100
    -- n = 2 bracket
    ∧ (38989 : ℚ) < alphaInvPowTimesGammaSum 2
    ∧ alphaInvPowTimesGammaSum 2 < 38990
    -- n = 3 bracket
    ∧ (5313434 : ℚ) < alphaInvPowTimesGammaSum 3
    ∧ alphaInvPowTimesGammaSum 3 < 5313435
    -- strict cascade
    ∧ alphaInvPowTimesGammaSum 1 < alphaInvPowTimesGammaSum 2
    ∧ alphaInvPowTimesGammaSum 2 < alphaInvPowTimesGammaSum 3
    -- duality
    ∧ alphaTimesGammaSum * alphaInvTimesGammaSum = gammaSumActive ^ 2
    ∧ alphaTimesGammaSum * alphaInvTimesGammaSum / gammaSumActive ^ 2 = 1
    -- loose cap
    ∧ (285 : ℚ) < alphaInvTimesGammaSum :=
  ⟨alphaInvPowTimesGammaSum_one_pos,
   alphaInvPowTimesGammaSum_two_pos,
   alphaInvPowTimesGammaSum_three_pos,
   alphaInvPowTimesGammaSum_one_bracket.1, alphaInvPowTimesGammaSum_one_bracket.2,
   alphaInvPowTimesGammaSum_two_bracket.1, alphaInvPowTimesGammaSum_two_bracket.2,
   alphaInvPowTimesGammaSum_three_bracket.1, alphaInvPowTimesGammaSum_three_bracket.2,
   alphaInvPowTimesGammaSum_inc_1_2, alphaInvPowTimesGammaSum_inc_2_3,
   alphaTimesGammaSum_mul_alphaInvTimesGammaSum,
   alphaTimesGammaSum_mul_alphaInvTimesGammaSum_div_gammaSq,
   alphaInvTimesGammaSum_gt_285⟩

end PT.Holonomy
