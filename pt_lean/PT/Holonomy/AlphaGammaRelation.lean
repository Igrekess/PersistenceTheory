/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Holonomy.CouplingReconstruction
import PT.Holonomy.CouplingReconstructionBounds
import PT.Holonomy.GammaProduct
import PT.Holonomy.GammaSumActive
import PT.Holonomy.InverseSinSqProduct
import Mathlib.Tactic

/-!
# Algebraic relations between `α_bare`, `Π γ_p` and `Σ γ_p` on the active cascade

**Statement (paper-level, Ch09 §"Reconstruction des couplages",
Ch06 §"Cascade arithmétique" — joint identities).** The PT fixed point
`μ⋆ = 15`, `q_+ = 13/15`, `s = 1/2` produces three canonical holonomy
invariants on the active set `{3, 5, 7}`:

* the **bare coupling** `α_bare = ∏_p sin²θ_p`
  (`PT.Holonomy.CouplingReconstruction.alphaBareQ`),
* the **γ-product** `Γ_active = ∏_p γ_p`
  (`PT.Holonomy.GammaProduct.gammaProductActive`),
* the **γ-sum**     `Σ_active = ∑_p γ_p`
  (`PT.Holonomy.GammaSumActive.gammaSumActive`).

The closed forms

  `γ_p = 4 q^{p-1}(1 - δ_p) / (μ⋆ · sin²θ_p)`,
  `sin²θ_p = δ_p (2 - δ_p)`,

immediately yield the **factor-by-factor identity**

  `sin²θ_p · γ_p = 4 q^{p-1}(1 - δ_p) / μ⋆`,

so the **algebraic product** `α_bare · Γ_active` collapses to the rational

  `α_bare · Γ_active = ∏_{p ∈ {3,5,7}} 4 q^{p-1}(1 - δ_p) / μ⋆`.

This module records the four standard mixed invariants:

1. `alphaGammaProductActive := α_bare · Γ_active` — bracket
   `0.002456 < · < 0.002458`, exact value `≈ 0.002457`.
2. `alphaInvOverGammaProductActive := (1 / α_bare) / Γ_active` —
   bracket `406.9 < · < 407.0`, exact value `≈ 406.966`.
3. `gammaSumOverAlphaActive := Σ_active / α_bare` — bracket
   `286.0 < · < 286.2`, exact value `≈ 286.103`.
4. `alphaGammaCubedActive := α_bare · Γ_active · s^3 = (α_bare · Γ_active) / 8`
   — bracket `306/10^6 < · < 308/10^6`, exact value `≈ 307·10⁻⁶`.

Plus two cross-invariants of common interest:

5. `gammaSumTimesProductActive := Σ_active · Γ_active` — bracket
   `0.702 < · < 0.704` (exact `≈ 0.7030`).
6. `alphaPlusGammaProductActive := α_bare + Γ_active` — bracket
   `0.342 < · < 0.343` (exact `≈ 0.3422`).

All identities are pure rational arithmetic; brackets decided by `norm_num`.

## Reference

* Monograph Ch06 §"Cascade arithmétique" and Ch09 §"Reconstruction des
  couplages" (joint α–γ identities).
* `PT.Holonomy.GammaProduct.gammaProductActive_over_alphaBare_bracket`
  records the dual ratio `Γ_active / α_bare ∈ (45.6, 45.7)`.
* `PT.Holonomy.InverseSinSqProduct.invProductActive_eq_alphaBareInv` records
  the inverse-coupling identity used in (2).
-/

namespace PT.Holonomy

/-! ### (1) Multiplicative pairing `α_bare · Γ_active` -/

/-- The **mixed product** `α_bare · Γ_active`: bare coupling times the
    γ-product over the active cascade. -/
def alphaGammaProductActive : ℚ := alphaBareQ * gammaProductActive

/-- Definitional unfolding. -/
theorem alphaGammaProductActive_eq :
    alphaGammaProductActive = alphaBareQ * gammaProductActive := rfl

/-- **Positivity.** `0 < α_bare · Γ_active`. -/
theorem alphaGammaProductActive_pos : 0 < alphaGammaProductActive := by
  unfold alphaGammaProductActive
  exact mul_pos alphaBareQ_pos gammaProductActive_pos

/-- **Tight bracket.** `0.002456 < α_bare · Γ_active < 0.002458`
    (exact `≈ 0.002457`). -/
theorem alphaGammaProductActive_bracket :
    2456 / 1000000 < alphaGammaProductActive
    ∧ alphaGammaProductActive < 2458 / 1000000 := by
  unfold alphaGammaProductActive alphaBareQ gammaProductActive
        sinSqQ gammaQ deltaQ qPT muStar
  refine ⟨?_, ?_⟩ <;> norm_num

/-! ### (2) Inverse pairing `α_bare⁻¹ / Γ_active` -/

/-- The **inverse pairing** `(1 / α_bare) / Γ_active`: scale at which
    the γ-product is comparable to the inverse coupling. -/
noncomputable def alphaInvOverGammaProductActive : ℚ :=
  (1 / alphaBareQ) / gammaProductActive

/-- **Identity in terms of the mixed product.**
    `(1 / α_bare) / Γ_active = 1 / (α_bare · Γ_active)
                              = 1 / alphaGammaProductActive`. -/
theorem alphaInvOverGammaProductActive_eq_inv :
    alphaInvOverGammaProductActive = 1 / alphaGammaProductActive := by
  unfold alphaInvOverGammaProductActive alphaGammaProductActive
  have ha : alphaBareQ ≠ 0 := ne_of_gt alphaBareQ_pos
  have hg : gammaProductActive ≠ 0 := ne_of_gt gammaProductActive_pos
  field_simp

/-- **Positivity.** -/
theorem alphaInvOverGammaProductActive_pos :
    0 < alphaInvOverGammaProductActive := by
  rw [alphaInvOverGammaProductActive_eq_inv]
  exact one_div_pos.mpr alphaGammaProductActive_pos

/-- **Tight bracket.** `406.9 < (1/α_bare)/Γ_active < 407.0`
    (exact `≈ 406.966`). -/
theorem alphaInvOverGammaProductActive_bracket :
    4069 / 10 < alphaInvOverGammaProductActive
    ∧ alphaInvOverGammaProductActive < 4070 / 10 := by
  unfold alphaInvOverGammaProductActive alphaBareQ gammaProductActive
        sinSqQ gammaQ deltaQ qPT muStar
  refine ⟨?_, ?_⟩ <;> norm_num

/-! ### (3) Additive-over-multiplicative `Σ_active / α_bare` -/

/-- The **sum-over-coupling ratio** `Σ_active / α_bare`. Pure rational
    number once `α_bare ≠ 0`. -/
noncomputable def gammaSumOverAlphaActive : ℚ := gammaSumActive / alphaBareQ

/-- **Positivity.** -/
theorem gammaSumOverAlphaActive_pos : 0 < gammaSumOverAlphaActive := by
  unfold gammaSumOverAlphaActive
  exact div_pos gammaSumActive_pos alphaBareQ_pos

/-- **Tight bracket.** `286.0 < Σ_active / α_bare < 286.2`
    (exact `≈ 286.103`). -/
theorem gammaSumOverAlphaActive_bracket :
    2860 / 10 < gammaSumOverAlphaActive
    ∧ gammaSumOverAlphaActive < 2862 / 10 := by
  unfold gammaSumOverAlphaActive gammaSumActive alphaBareQ
        sinSqQ gammaQ deltaQ qPT muStar
  refine ⟨?_, ?_⟩ <;> norm_num

/-! ### (4) Cubed-symmetry pairing `α_bare · Γ_active · s³` -/

/-- The **`s³`-scaled mixed product**
    `α_bare · Γ_active · s^3 = (α_bare · Γ_active) / 8`. Records the
    *cubed-symmetry* invariant, i.e. the response of the mixed product
    under three insertions of the PT symmetry parameter `s = 1/2`. -/
def alphaGammaCubedActive : ℚ := alphaBareQ * gammaProductActive * sPT ^ 3

/-- **Identity.** `α_bare · Γ_active · s^3 = alphaGammaProductActive / 8`. -/
theorem alphaGammaCubedActive_eq_div_eight :
    alphaGammaCubedActive = alphaGammaProductActive / 8 := by
  unfold alphaGammaCubedActive alphaGammaProductActive sPT
  ring

/-- **Positivity.** -/
theorem alphaGammaCubedActive_pos : 0 < alphaGammaCubedActive := by
  rw [alphaGammaCubedActive_eq_div_eight]
  exact div_pos alphaGammaProductActive_pos (by norm_num : (0 : ℚ) < 8)

/-- **Tight bracket.** `306 / 10⁶ < α_bare · Γ_active · s^3 < 308 / 10⁶`
    (exact `≈ 307·10⁻⁶`). -/
theorem alphaGammaCubedActive_bracket :
    306 / 1000000 < alphaGammaCubedActive
    ∧ alphaGammaCubedActive < 308 / 1000000 := by
  unfold alphaGammaCubedActive alphaBareQ gammaProductActive sPT
        sinSqQ gammaQ deltaQ qPT muStar
  refine ⟨?_, ?_⟩ <;> norm_num

/-! ### (5) Cross invariant `Σ_active · Γ_active` -/

/-- The **sum-times-product** invariant `Σ_active · Γ_active`. -/
def gammaSumTimesProductActive : ℚ := gammaSumActive * gammaProductActive

/-- **Positivity.** -/
theorem gammaSumTimesProductActive_pos : 0 < gammaSumTimesProductActive := by
  unfold gammaSumTimesProductActive
  exact mul_pos gammaSumActive_pos gammaProductActive_pos

/-- **Tight bracket.** `0.702 < Σ_active · Γ_active < 0.704`
    (exact `≈ 0.7030`). -/
theorem gammaSumTimesProductActive_bracket :
    702 / 1000 < gammaSumTimesProductActive
    ∧ gammaSumTimesProductActive < 704 / 1000 := by
  unfold gammaSumTimesProductActive gammaSumActive gammaProductActive
        gammaQ deltaQ qPT muStar
  refine ⟨?_, ?_⟩ <;> norm_num

/-! ### (6) Additive combination `α_bare + Γ_active` -/

/-- The **additive combination** `α_bare + Γ_active`. -/
def alphaPlusGammaProductActive : ℚ := alphaBareQ + gammaProductActive

/-- **Positivity.** -/
theorem alphaPlusGammaProductActive_pos : 0 < alphaPlusGammaProductActive := by
  unfold alphaPlusGammaProductActive
  exact add_pos alphaBareQ_pos gammaProductActive_pos

/-- **Tight bracket.** `0.342 < α_bare + Γ_active < 0.343`
    (exact `≈ 0.3422`). The sum is dominated by `Γ_active` (since
    `α_bare ≈ 7.34·10⁻³` is two orders of magnitude smaller). -/
theorem alphaPlusGammaProductActive_bracket :
    342 / 1000 < alphaPlusGammaProductActive
    ∧ alphaPlusGammaProductActive < 343 / 1000 := by
  unfold alphaPlusGammaProductActive alphaBareQ gammaProductActive
        sinSqQ gammaQ deltaQ qPT muStar
  refine ⟨?_, ?_⟩ <;> norm_num

/-- The additive combination is strictly closer to `Γ_active` than to
    `α_bare`: `Γ_active < α_bare + Γ_active`. -/
theorem alphaPlusGammaProductActive_gt_gammaProductActive :
    gammaProductActive < alphaPlusGammaProductActive := by
  unfold alphaPlusGammaProductActive
  have := alphaBareQ_pos
  linarith

/-! ### Algebraic factor-by-factor identity

For each active prime, the product `sin²θ_p · γ_p` reduces to
`4 q^{p-1}(1 - δ_p)/μ⋆` (closed form of `γ_p`). This is the local
ingredient behind `alphaGammaProductActive`. -/

-- We give per-prime closed-form factor identities. A "generic" statement
-- over arbitrary `p` would require non-vanishing of `δ_p` and `2 - δ_p`,
-- which is parameter-dependent; the three concrete cases used downstream
-- discharge by `norm_num` directly.

/-- **Factor identity at `p = 3`.** -/
theorem sinSq_mul_gamma_3 :
    sinSqQ 3 * gammaQ 3 = 4 * qPT ^ 2 * (1 - deltaQ 3) / muStar := by
  unfold sinSqQ gammaQ deltaQ qPT muStar
  norm_num

/-- **Factor identity at `p = 5`.** -/
theorem sinSq_mul_gamma_5 :
    sinSqQ 5 * gammaQ 5 = 4 * qPT ^ 4 * (1 - deltaQ 5) / muStar := by
  unfold sinSqQ gammaQ deltaQ qPT muStar
  norm_num

/-- **Factor identity at `p = 7`.** -/
theorem sinSq_mul_gamma_7 :
    sinSqQ 7 * gammaQ 7 = 4 * qPT ^ 6 * (1 - deltaQ 7) / muStar := by
  unfold sinSqQ gammaQ deltaQ qPT muStar
  norm_num

/-- **Aggregate identity.** The active mixed product expands as a triple
    product of `4 q^{p-1}(1 - δ_p)/μ⋆` over `p ∈ {3, 5, 7}`. -/
theorem alphaGammaProductActive_eq_factored :
    alphaGammaProductActive
      = (4 * qPT ^ 2 * (1 - deltaQ 3) / muStar)
      * (4 * qPT ^ 4 * (1 - deltaQ 5) / muStar)
      * (4 * qPT ^ 6 * (1 - deltaQ 7) / muStar) := by
  unfold alphaGammaProductActive alphaBareQ gammaProductActive
        sinSqQ gammaQ deltaQ qPT muStar
  norm_num

/-! ### Link with `invProductActive = α_bare⁻¹` -/

/-- **Reciprocal interpretation.**
    `α_bare⁻¹ / Γ_active = invProductActive / Γ_active`. -/
theorem alphaInvOverGammaProductActive_eq_invProductActive_div :
    alphaInvOverGammaProductActive = invProductActive / gammaProductActive := by
  unfold alphaInvOverGammaProductActive
  rw [invProductActive_eq_inv_alphaBare]

/-! ### Headline -/

/-- **Headline (α–γ relation summary).** At the PT fixed point
    `q_+ = 13/15`, `μ⋆ = 15`, `s = 1/2`, the four canonical mixed
    invariants on the active cascade `{3, 5, 7}` satisfy:

    1. `α_bare · Γ_active ∈ (0.002456, 0.002458)`, strictly positive.
    2. `(1/α_bare) / Γ_active = 1 / (α_bare · Γ_active) ∈ (40.69·10, 40.70·10)`,
       i.e. `(406.9, 407.0)`.
    3. `Σ_active / α_bare ∈ (286.0, 286.2)`.
    4. `α_bare · Γ_active · s^3 = (α_bare · Γ_active)/8 ∈ (306·10⁻⁶, 308·10⁻⁶)`.
    5. `Σ_active · Γ_active ∈ (0.702, 0.704)`.
    6. `α_bare + Γ_active ∈ (0.342, 0.343)`, strictly above `Γ_active`. -/
theorem alpha_gamma_relation_summary :
    0 < alphaGammaProductActive
    ∧ 2456 / 1000000 < alphaGammaProductActive
    ∧ alphaGammaProductActive < 2458 / 1000000
    ∧ alphaInvOverGammaProductActive = 1 / alphaGammaProductActive
    ∧ 4069 / 10 < alphaInvOverGammaProductActive
    ∧ alphaInvOverGammaProductActive < 4070 / 10
    ∧ 2860 / 10 < gammaSumOverAlphaActive
    ∧ gammaSumOverAlphaActive < 2862 / 10
    ∧ alphaGammaCubedActive = alphaGammaProductActive / 8
    ∧ 306 / 1000000 < alphaGammaCubedActive
    ∧ alphaGammaCubedActive < 308 / 1000000
    ∧ 702 / 1000 < gammaSumTimesProductActive
    ∧ gammaSumTimesProductActive < 704 / 1000
    ∧ 342 / 1000 < alphaPlusGammaProductActive
    ∧ alphaPlusGammaProductActive < 343 / 1000
    ∧ gammaProductActive < alphaPlusGammaProductActive :=
  ⟨alphaGammaProductActive_pos,
   alphaGammaProductActive_bracket.1, alphaGammaProductActive_bracket.2,
   alphaInvOverGammaProductActive_eq_inv,
   alphaInvOverGammaProductActive_bracket.1,
   alphaInvOverGammaProductActive_bracket.2,
   gammaSumOverAlphaActive_bracket.1, gammaSumOverAlphaActive_bracket.2,
   alphaGammaCubedActive_eq_div_eight,
   alphaGammaCubedActive_bracket.1, alphaGammaCubedActive_bracket.2,
   gammaSumTimesProductActive_bracket.1,
   gammaSumTimesProductActive_bracket.2,
   alphaPlusGammaProductActive_bracket.1,
   alphaPlusGammaProductActive_bracket.2,
   alphaPlusGammaProductActive_gt_gammaProductActive⟩

end PT.Holonomy
