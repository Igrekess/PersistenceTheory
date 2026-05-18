/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Holonomy.CouplingReconstruction
import PT.Holonomy.CouplingReconstructionBounds
import PT.Holonomy.InverseSinSqProduct
import PT.Holonomy.ActivePrimeCriterion
import Mathlib.Tactic

/-!
# `(1/α_bare)²` and higher inverse powers (Ch09 extension)

This module records the **square** and **cube** of the bare inverse
fine-structure coupling `α_bare⁻¹` at the PT fixed point
`q_+ = 13/15`, `μ⋆ = 15`, plus a comparison against the canonical scaling
`μ⋆ · α_bare⁻¹` that appears in PT mass / coupling running formulae.

Numerical values (exact rationals at the fixed point):

* `α_bare⁻¹                ≈ 136.278`     (recalled from
  `alphaBareInvQ_bracket`).
* `(α_bare⁻¹)²            ≈ 18571.78`     — *headline of this module*.
* `(α_bare⁻¹)³            ≈ 2 530 931.85`.
* `μ⋆ · α_bare⁻¹          ≈ 2044.18`.
* `(α_bare⁻¹)² / (μ⋆ · α_bare⁻¹) = α_bare⁻¹ / μ⋆ ≈ 9.0852`.

The algebraic content is straightforward:

  `(α_bare⁻¹)² = 1 / α_bare²`     and     `(α_bare⁻¹)³ = 1 / α_bare³`,

so the value question reduces to bracketing `α_bare⁻¹` itself (already
done in `CouplingReconstructionBounds.alphaBareInvQ_bracket`) and to a
single `norm_num` pass on the rational squaring/cubing.

All brackets are exact rational arithmetic; positivity is inherited from
`alphaBareInvQ_pos`.

## Reference

* Monograph Ch09 §"Reconstruction des couplages" — `α_bare⁻¹` headline.
* `PT.Holonomy.CouplingReconstructionBounds.alphaBareInvQ_bracket`.
* `PT.Holonomy.InverseSinSqProduct.invProductActive_eq_alphaBareInv`.
-/

namespace PT.Holonomy

/-! ### Definition: `invAlphaSquared = (1 / α_bare)²` -/

/-- The **squared inverse bare coupling** `(α_bare⁻¹)²`. Pure rational at
    the PT fixed point. -/
noncomputable def invAlphaSquared : ℚ := alphaBareInvQ * alphaBareInvQ

/-- Definitional unfolding. -/
theorem invAlphaSquared_eq :
    invAlphaSquared = alphaBareInvQ * alphaBareInvQ := rfl

/-- **Algebraic identity (square).**
    `(α_bare⁻¹)² = 1 / α_bare²`. -/
theorem invAlphaSquared_eq_inv_alpha_sq :
    invAlphaSquared = 1 / (alphaBareQ * alphaBareQ) := by
  unfold invAlphaSquared alphaBareInvQ
  have ha : alphaBareQ ≠ 0 := ne_of_gt alphaBareQ_pos
  field_simp

/-- **Positivity.** `0 < (α_bare⁻¹)²`. -/
theorem invAlphaSquared_pos : 0 < invAlphaSquared := by
  unfold invAlphaSquared
  exact mul_pos alphaBareInvQ_pos alphaBareInvQ_pos

/-! ### Brackets for `invAlphaSquared` -/

/-- **Loose decimal bracket.** `18500 < (α_bare⁻¹)² < 18600`
    (exact value `≈ 18571.78`). -/
theorem invAlphaSquared_bracket :
    18500 < invAlphaSquared ∧ invAlphaSquared < 18600 := by
  unfold invAlphaSquared alphaBareInvQ alphaBareQ sinSqQ deltaQ qPT
  refine ⟨?_, ?_⟩ <;> norm_num

/-- **Tight decimal bracket.** `18570 < (α_bare⁻¹)² < 18575`
    (exact value `≈ 18571.78`). -/
theorem invAlphaSquared_bracket_tight :
    18570 < invAlphaSquared ∧ invAlphaSquared < 18575 := by
  unfold invAlphaSquared alphaBareInvQ alphaBareQ sinSqQ deltaQ qPT
  refine ⟨?_, ?_⟩ <;> norm_num

/-! ### Comparison with `μ⋆ · α_bare⁻¹` -/

/-- The **PT-running scale** `μ⋆ · α_bare⁻¹ = 15 · α_bare⁻¹`.
    Identically equal to `μ⋆ · invProductActive`. -/
noncomputable def muStarTimesAlphaInv : ℚ := muStar * alphaBareInvQ

/-- Rewriting through `invProductActive`. -/
theorem muStarTimesAlphaInv_eq_muStar_invProductActive :
    muStarTimesAlphaInv = muStar * invProductActive := by
  unfold muStarTimesAlphaInv
  rw [invProductActive_eq_alphaBareInv]

/-- **Positivity.** -/
theorem muStarTimesAlphaInv_pos : 0 < muStarTimesAlphaInv := by
  unfold muStarTimesAlphaInv muStar
  exact mul_pos (by norm_num : (0 : ℚ) < 15) alphaBareInvQ_pos

/-- **Bracket.** `2044 < μ⋆ · α_bare⁻¹ < 2045` (exact `≈ 2044.18`). -/
theorem muStarTimesAlphaInv_bracket :
    2044 < muStarTimesAlphaInv ∧ muStarTimesAlphaInv < 2045 := by
  unfold muStarTimesAlphaInv alphaBareInvQ alphaBareQ sinSqQ deltaQ qPT muStar
  refine ⟨?_, ?_⟩ <;> norm_num

/-- The **scale ratio** `(α_bare⁻¹)² / (μ⋆ · α_bare⁻¹) = α_bare⁻¹ / μ⋆`. -/
noncomputable def invAlphaSquaredOverScale : ℚ :=
  invAlphaSquared / muStarTimesAlphaInv

/-- **Identity.** `(α_bare⁻¹)² / (μ⋆ · α_bare⁻¹) = α_bare⁻¹ / μ⋆`. -/
theorem invAlphaSquaredOverScale_eq_alphaInv_div_muStar :
    invAlphaSquaredOverScale = alphaBareInvQ / muStar := by
  unfold invAlphaSquaredOverScale invAlphaSquared muStarTimesAlphaInv
  have hα : alphaBareInvQ ≠ 0 := ne_of_gt alphaBareInvQ_pos
  have hμ : (muStar : ℚ) ≠ 0 := by unfold muStar; norm_num
  field_simp

/-- **Positivity.** -/
theorem invAlphaSquaredOverScale_pos : 0 < invAlphaSquaredOverScale := by
  rw [invAlphaSquaredOverScale_eq_alphaInv_div_muStar]
  unfold muStar
  exact div_pos alphaBareInvQ_pos (by norm_num : (0 : ℚ) < 15)

/-- **Bracket.** `9.08 < (α_bare⁻¹)² / (μ⋆ · α_bare⁻¹) < 9.09`
    (exact value `≈ 9.0852`). -/
theorem invAlphaSquaredOverScale_bracket :
    908 / 100 < invAlphaSquaredOverScale
    ∧ invAlphaSquaredOverScale < 909 / 100 := by
  rw [invAlphaSquaredOverScale_eq_alphaInv_div_muStar]
  unfold alphaBareInvQ alphaBareQ sinSqQ deltaQ qPT muStar
  refine ⟨?_, ?_⟩ <;> norm_num

/-! ### Cube `(α_bare⁻¹)³` -/

/-- The **cubed inverse bare coupling** `(α_bare⁻¹)³`. -/
noncomputable def invAlphaCubed : ℚ := alphaBareInvQ ^ 3

/-- **Algebraic identity (cube).** `(α_bare⁻¹)³ = 1 / α_bare³`. -/
theorem invAlphaCubed_eq_inv_alpha_cubed :
    invAlphaCubed = 1 / alphaBareQ ^ 3 := by
  unfold invAlphaCubed alphaBareInvQ
  have ha : alphaBareQ ≠ 0 := ne_of_gt alphaBareQ_pos
  field_simp

/-- **Positivity.** -/
theorem invAlphaCubed_pos : 0 < invAlphaCubed := by
  unfold invAlphaCubed
  exact pow_pos alphaBareInvQ_pos 3

/-- **Identity linking square and cube.**
    `(α_bare⁻¹)³ = (α_bare⁻¹)² · α_bare⁻¹`. -/
theorem invAlphaCubed_eq_sq_mul_inv :
    invAlphaCubed = invAlphaSquared * alphaBareInvQ := by
  unfold invAlphaCubed invAlphaSquared
  ring

/-- **Loose bracket.** `2 500 000 < (α_bare⁻¹)³ < 2 600 000`
    (exact `≈ 2 530 931.85`). -/
theorem invAlphaCubed_bracket :
    2500000 < invAlphaCubed ∧ invAlphaCubed < 2600000 := by
  unfold invAlphaCubed alphaBareInvQ alphaBareQ sinSqQ deltaQ qPT
  refine ⟨?_, ?_⟩ <;> norm_num

/-- **Tight bracket.** `2 530 000 < (α_bare⁻¹)³ < 2 531 000`
    (exact `≈ 2 530 931.85`). -/
theorem invAlphaCubed_bracket_tight :
    2530000 < invAlphaCubed ∧ invAlphaCubed < 2531000 := by
  unfold invAlphaCubed alphaBareInvQ alphaBareQ sinSqQ deltaQ qPT
  refine ⟨?_, ?_⟩ <;> norm_num

/-! ### Strict monotone chain `α_bare⁻¹ < (α_bare⁻¹)² < (α_bare⁻¹)³` -/

/-- The inverse-power tower is strictly increasing (since `α_bare⁻¹ > 1`).
    `α_bare⁻¹ < (α_bare⁻¹)² < (α_bare⁻¹)³`. -/
theorem invAlphaPower_chain :
    alphaBareInvQ < invAlphaSquared ∧ invAlphaSquared < invAlphaCubed := by
  have hpos : 0 < alphaBareInvQ := alphaBareInvQ_pos
  have hgt1 : 1 < alphaBareInvQ := by
    have := alphaBareInvQ_gt_135
    linarith
  refine ⟨?_, ?_⟩
  · -- α⁻¹ < (α⁻¹)² ↔ 1 < α⁻¹ (multiply by α⁻¹ > 0).
    unfold invAlphaSquared
    have := (lt_mul_iff_one_lt_right hpos).mpr hgt1
    -- α⁻¹ < α⁻¹ * α⁻¹
    linarith
  · -- (α⁻¹)² < (α⁻¹)³  ↔  1 < α⁻¹ (multiply by (α⁻¹)² > 0).
    rw [invAlphaCubed_eq_sq_mul_inv]
    have hsq : 0 < invAlphaSquared := invAlphaSquared_pos
    have := (lt_mul_iff_one_lt_right hsq).mpr hgt1
    linarith

/-! ### Headline -/

/-- **Headline.** All invariants at a glance for the squared / cubed
    inverse bare coupling, at the PT fixed point `q_+ = 13/15`,
    `μ⋆ = 15`:

    * `(α_bare⁻¹)² ∈ (18500, 18600)`, tight `(18570, 18575)`, exact ≈ 18571.78.
    * `(α_bare⁻¹)³ ∈ (2 500 000, 2 600 000)`, tight `(2 530 000, 2 531 000)`,
      exact ≈ 2 530 931.85.
    * `μ⋆ · α_bare⁻¹ ∈ (2044, 2045)`, exact ≈ 2044.18.
    * `(α_bare⁻¹)² / (μ⋆ · α_bare⁻¹) = α_bare⁻¹ / μ⋆ ∈ (9.08, 9.09)`,
      exact ≈ 9.0852.
    * Strict monotonicity: `α_bare⁻¹ < (α_bare⁻¹)² < (α_bare⁻¹)³`.
    * Algebraic identities: `(α_bare⁻¹)² = 1/α_bare²`,
      `(α_bare⁻¹)³ = 1/α_bare³`. -/
theorem invAlphaSquared_summary :
    -- positivity
    0 < invAlphaSquared
    ∧ 0 < invAlphaCubed
    ∧ 0 < muStarTimesAlphaInv
    ∧ 0 < invAlphaSquaredOverScale
    -- square: loose + tight
    ∧ 18500 < invAlphaSquared ∧ invAlphaSquared < 18600
    ∧ 18570 < invAlphaSquared ∧ invAlphaSquared < 18575
    -- cube: loose + tight
    ∧ 2500000 < invAlphaCubed ∧ invAlphaCubed < 2600000
    ∧ 2530000 < invAlphaCubed ∧ invAlphaCubed < 2531000
    -- scale ratio
    ∧ 2044 < muStarTimesAlphaInv ∧ muStarTimesAlphaInv < 2045
    ∧ 908 / 100 < invAlphaSquaredOverScale
    ∧ invAlphaSquaredOverScale < 909 / 100
    -- algebraic identities
    ∧ invAlphaSquared = 1 / (alphaBareQ * alphaBareQ)
    ∧ invAlphaCubed = 1 / alphaBareQ ^ 3
    ∧ invAlphaSquaredOverScale = alphaBareInvQ / muStar
    -- monotone chain
    ∧ alphaBareInvQ < invAlphaSquared
    ∧ invAlphaSquared < invAlphaCubed :=
  ⟨invAlphaSquared_pos,
   invAlphaCubed_pos,
   muStarTimesAlphaInv_pos,
   invAlphaSquaredOverScale_pos,
   invAlphaSquared_bracket.1, invAlphaSquared_bracket.2,
   invAlphaSquared_bracket_tight.1, invAlphaSquared_bracket_tight.2,
   invAlphaCubed_bracket.1, invAlphaCubed_bracket.2,
   invAlphaCubed_bracket_tight.1, invAlphaCubed_bracket_tight.2,
   muStarTimesAlphaInv_bracket.1, muStarTimesAlphaInv_bracket.2,
   invAlphaSquaredOverScale_bracket.1, invAlphaSquaredOverScale_bracket.2,
   invAlphaSquared_eq_inv_alpha_sq,
   invAlphaCubed_eq_inv_alpha_cubed,
   invAlphaSquaredOverScale_eq_alphaInv_div_muStar,
   invAlphaPower_chain.1,
   invAlphaPower_chain.2⟩

end PT.Holonomy
