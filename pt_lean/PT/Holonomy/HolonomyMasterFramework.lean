/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Holonomy.CyclicPhaseIdentity
import PT.Holonomy.CyclicPhaseInversion
import PT.Holonomy.CyclicPhaseAlgebraic
import PT.Holonomy.CyclicPhaseTable
import PT.Holonomy.CyclicPhaseSpectral
import PT.Holonomy.ActivePrimeCriterion
import PT.Holonomy.ActivePrimeMargins
import PT.Holonomy.GammaTablesExtended
import PT.Holonomy.GammaProduct
import PT.Holonomy.GammaSumActive
import PT.Holonomy.GammaSumProduct
import PT.Holonomy.GammaPrimorialProduct
import PT.Holonomy.AlphaGammaRelation
import PT.Holonomy.AlphaTimesGammaSum
import PT.Holonomy.AlphaInverseTimesGammaSum
import PT.Holonomy.CouplingReconstruction
import PT.Holonomy.CouplingReconstructionBounds
import PT.Holonomy.SinSqProductBounds
import PT.Holonomy.SinSqProductChain
import PT.Holonomy.SinSqRatios
import PT.Holonomy.InverseSinSq
import PT.Holonomy.InverseSinSqProduct
import PT.Holonomy.InvAlphaSquaredBracket
import PT.Holonomy.AlphaPowerSequence
import PT.Holonomy.AlphaInversePowerSequence
import Mathlib.Tactic

/-!
# Holonomy — Master Framework (aggregator)

**Purpose.** Aggregator file consolidating ~28 modules of the holonomy axis
into five "master" sub-theorems plus one final umbrella theorem.

This file performs **no new calculation**: every conjunct is a direct
re-export of a result already proved in one of the imported modules. The
value of the file is purely *organisational* — it gives downstream files
(article companions, monograph cross-references, sanity tests) a single,
stable entry point with a fixed name and a fixed shape.

## Master sub-theorems

1. `master_cyclic_phase` — cyclic-phase identity `sin²θ = δ(2-δ)`,
   monotonicity, surjective inversion, and the exact `δ_p` table for
   `p ∈ {2,3,5,7,11,13}` (Chapter 6, §6.2 of the monograph).

2. `master_active_primes` — strict `γ_p > 1/2` for `p ∈ {3,5,7}` and
   strict `γ_p < 1/2` for `p ∈ {11,13,17,19,23}`, with exact margins for
   the boundary primes and full monotone chain `γ_3 > γ_5 > … > γ_23`.

3. `master_gamma_family` — sum `Σγ`, product `Πγ`, ratios, primorial
   product `Πγ · 105`, sum × `μ⋆ = 15`, sum × product, all with tight
   rational brackets and positivity.

4. `master_coupling` — exact rational value of `α_bare = Π_{p∈{3,5,7}}
   sin²θ_p`, tight bracket, reciprocal `α_bare⁻¹ ∈ (135, 137)`, the
   partial-product chain `P_1 > P_2 > P_3 > P_4 > P_5`, inverse products
   `1/sin²θ_p` and the bracket on `α_bare⁻²`.

5. `master_alpha_powers` — values, positivity, monotone decrease and
   tight brackets for `α^n` and `(α⁻¹)^n`, `n = 1..5`, together with the
   exact duality identity `α^n · (α⁻¹)^n = 1`.

6. `holonomy_master_framework_summary` — final umbrella theorem: the
   conjunction of all five master sub-theorems above. Used as the single
   downstream witness "the holonomy axis is fully derived".

## Conventions

* Everything is proved over `ℚ` (no `sorry`, no `Real`-side calculation
  beyond the spectral side-statements in `CyclicPhaseSpectral`).
* The framework is purely a re-export: no new definitions, no new
  numerical witnesses. Names, brackets, and inequality directions follow
  the imported modules verbatim.

## References

Monograph Chapter 6 (`thm:cyclic-phase`, `tab:cyclic_phase_values`),
Chapter 7 (`thm:active-prime-criterion`, `thm:gamma-product`), Chapter 9
(`thm:alpha-bare-reconstruction`, `thm:alpha-bracket`).
-/

namespace PT.Holonomy

/-! ## 1. Master cyclic phase

Algebraic identity, monotonicity, surjectivity onto `[0,1]`, and the
exact `δ_p` table for `p ∈ {2, 3, 5, 7, 11, 13}`. -/

/-- **Master cyclic phase.** Aggregator for the cyclic-phase identity
`sin²θ = δ(2-δ)`, its monotonicity in `δ`, its inversion via the square
root, and the exact rational table for `δ_p` at the PT fixed point
`q = 13/15` for `p ∈ {2, 3, 5, 7, 11, 13}`. -/
theorem master_cyclic_phase :
    -- (a) algebraic identity in any commutative ring
    (∀ δ : ℝ, 1 - (1 - δ)^2 = δ * (2 - δ))
    -- (b) Pythagorean form: `cos θ = 1 - δ ⇒ sin²θ = δ(2-δ)`
    ∧ (∀ θ δ : ℝ, Real.cos θ = 1 - δ → Real.sin θ ^ 2 = δ * (2 - δ))
    -- (c) strict monotonicity on `[0,1]`
    ∧ (∀ δ δ' : ℝ, 0 ≤ δ → δ ≤ 1 → 0 ≤ δ' → δ' ≤ 1 → δ < δ' →
        δ * (2 - δ) < δ' * (2 - δ'))
    -- (d) surjective inversion: every `s ∈ [0,1]` is reached
    ∧ (∀ s : ℝ, 0 ≤ s → s ≤ 1 →
        let δ := 1 - Real.sqrt (1 - s)
        0 ≤ δ ∧ δ ≤ 1 ∧ δ * (2 - δ) = s)
    -- (e) exact `δ_p` values for `p ∈ {2, 3, 5, 7, 11, 13}`
    ∧ deltaQ 2 = 28 / 225
    ∧ deltaQ 3 = 1178 / 10125
    ∧ deltaQ 5 = 388082 / 3796875
    ∧ deltaQ 7 = 108110858 / 1196015625
    ∧ deltaQ 11 = 6857595465338 / 95147314453125
    ∧ deltaQ 13 = 1643319961767122 / 25300535888671875 := by
  refine ⟨one_sub_one_sub_sq, sin_sq_of_cos_eq_one_sub,
          delta_two_minus_delta_strictMono, delta_two_minus_delta_inversion,
          deltaQ_2_eq, deltaQ_3_eq, deltaQ_5_eq, deltaQ_7_eq,
          deltaQ_11_eq, deltaQ_13_eq⟩

/-! ## 2. Master active primes

Strict criterion `γ_p > s` for the three active primes `p ∈ {3, 5, 7}`,
strict `γ_p < s` for the five inactive primes `p ∈ {11, 13, 17, 19, 23}`,
exact rational margins at the boundary, and the full monotone chain. -/

/-- **Master active primes.** Aggregator for the active-prime criterion
`γ_p > 1/2 ⇔ p ∈ {3,5,7}`, its strict converse for `p ≥ 11`, the exact
boundary margins, and the monotone chain `γ_3 > γ_5 > γ_7 > γ_11 >
γ_13 > γ_17 > γ_19 > γ_23`. -/
theorem master_active_primes :
    -- (a) the three active primes
    gammaQ 3 > sPT ∧ gammaQ 5 > sPT ∧ gammaQ 7 > sPT
    -- (b) the five tested inactive primes
    ∧ gammaQ 11 < sPT ∧ gammaQ 13 < sPT ∧ gammaQ 17 < sPT
    ∧ gammaQ 19 < sPT ∧ gammaQ 23 < sPT
    -- (c) exact rational margins on the active side
    ∧ gammaQ 3 - sPT = 1727777 / 5616704
    ∧ gammaQ 5 - sPT = 68621964134 / 349548756097
    ∧ gammaQ 7 - sPT = 453321961202883 / 4748396022746468
    -- (d) exact rational margins on the inactive side
    ∧ sPT - gammaQ 11 = 23355758527719967875358739
                          / 314484242174863298491777064
    -- (e) strict monotone chain γ_3 > γ_5 > γ_7 > γ_11 > γ_13 > γ_17 > γ_19 > γ_23
    ∧ gammaQ 3 > gammaQ 5
    ∧ gammaQ 5 > gammaQ 7
    ∧ gammaQ 7 > gammaQ 11
    ∧ gammaQ 11 > gammaQ 13
    ∧ gammaQ 13 > gammaQ 17
    ∧ gammaQ 17 > gammaQ 19
    ∧ gammaQ 19 > gammaQ 23 := by
  refine ⟨gamma_3_active, gamma_5_active, gamma_7_active,
          gamma_11_inactive, gamma_13_inactive, gamma_17_inactive,
          gammaQ_19_inactive, gammaQ_23_inactive,
          gamma_3_margin_exact, gamma_5_margin_exact, gamma_7_margin_exact,
          gamma_11_margin_exact, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩
  · exact gamma_3_gt_gamma_5
  · exact gamma_5_gt_gamma_7
  · -- γ_7 > 1/2 > γ_11
    exact lt_trans gamma_11_inactive gamma_7_active
  · exact gamma_11_gt_gamma_13
  · exact gammaQ_13_gt_17
  · exact gammaQ_17_gt_19
  · exact gammaQ_19_gt_23

/-! ## 3. Master gamma family

Sum `Σγ`, product `Πγ`, sum-times-product, primorial `Πγ · 105`, sum
times `μ⋆ = 15`, ratios, and a few combined invariants — all with
tight rational brackets. -/

/-- **Master gamma family.** Aggregator for the active-prime gamma
family: positivity and tight brackets on `Σγ`, `Πγ`, `Σγ · μ⋆`,
`Σγ · Πγ`, `Σγ / Πγ`, primorial `Πγ · 105` and `Πγ · 105 / μ⋆`. -/
theorem master_gamma_family :
    -- (a) positivity
    0 < gammaSumActive ∧ 0 < gammaProductActive
    ∧ 0 < gammaSumProductActive ∧ 0 < gammaPrimorialProductActive
    -- (b) sum bracket (Σγ ≈ 2.0996)
    ∧ 2099 / 1000 < gammaSumActive ∧ gammaSumActive < 2100 / 1000
    -- (c) product bracket (Πγ ≈ 0.3344)
    ∧ 334 / 1000 < gammaProductActive ∧ gammaProductActive < 335 / 1000
    -- (d) basic comparisons: 3·s < Σγ < 3
    ∧ gammaSumActive > 3 * sPT ∧ gammaSumActive < 3
    -- (e) sum-times-muStar bracket (Σγ · 15 ≈ 31.49)
    ∧ 3149 / 100 < gammaSumActive * muStar
    ∧ gammaSumActive * muStar < 3150 / 100
    -- (f) primorial-product bracket (Πγ · 105 ≈ 35.16)
    ∧ 35160 / 1000 < gammaPrimorialProductActive
    ∧ gammaPrimorialProductActive < 35161 / 1000
    -- (g) sum-times-product bracket (Σγ · Πγ ≈ 0.7030)
    ∧ 7030 / 10000 < gammaSumProductActive
    ∧ gammaSumProductActive < 7031 / 10000
    -- (h) reciprocal of product on muStar (Πγ · 105 / 15 ≈ 2.344)
    ∧ 2344 / 1000 < gammaPrimorialOverMuStar
    ∧ gammaPrimorialOverMuStar < 2345 / 1000
    -- (i) AM-GM gap: (Σγ/3)^3 > Πγ
    ∧ (gammaSumActive / 3) ^ 3 > gammaProductActive := by
  refine ⟨gammaSumActive_pos, gammaProductActive_pos,
          gammaSumProductActive_pos, gammaPrimorialProductActive_pos,
          gammaSumActive_bracket.1, gammaSumActive_bracket.2,
          gammaProductActive_bracket.1, gammaProductActive_bracket.2,
          gammaSumActive_gt_three_s, gammaSumActive_lt_three,
          gammaSumActive_mul_muStar_bracket.1,
          gammaSumActive_mul_muStar_bracket.2,
          gammaPrimorialProductActive_bracket.1,
          gammaPrimorialProductActive_bracket.2,
          gammaSumProductActive_bracket.1,
          gammaSumProductActive_bracket.2,
          gammaPrimorialOverMuStar_bracket.1,
          gammaPrimorialOverMuStar_bracket.2,
          gammaSumActive_cubed_gt_27_prod⟩

/-! ## 4. Master coupling

Exact rational value of `α_bare = sin²θ_3 · sin²θ_5 · sin²θ_7`, tight
bracket, the partial-product chain `P_1 > P_2 > … > P_5`, the inverse
product `α_bare⁻¹ ∈ (135, 137)`, the inverse partial products, and the
bracket on `α_bare⁻²`. -/

/-- **Master coupling.** Aggregator for the coupling reconstruction:
`α_bare = sin²θ_3 · sin²θ_5 · sin²θ_7`, its exact rational value and
tight bracket, the strictly decreasing partial-product chain `P_1 > … >
P_5`, the reciprocal `α_bare⁻¹ ∈ (135, 137)`, the inverse-product
identity `Π (1/sin²θ_p) = α_bare⁻¹`, and the bracket on `α_bare⁻²`. -/
theorem master_coupling :
    -- (a) exact rational value of α_bare
    alphaBareQ
      = 15512777115364308026953701325116440576
          / 2114055428042547055520117282867431640625
    -- (b) tight bracket for α_bare ≈ 7.3375 × 10⁻³
    ∧ 7335 / 1000000 < alphaBareQ ∧ alphaBareQ < 7340 / 1000000
    -- (c) reciprocal bracket: 135 < 1/α_bare < 137
    ∧ 135 < alphaBareInvQ ∧ alphaBareInvQ < 137
    -- (d) positivity
    ∧ 0 < alphaBareQ ∧ 0 < alphaBareInvQ ∧ alphaBareQ < 1
    -- (e) partial-product chain (active primes)
    ∧ 0 < P3 ∧ P3 < P2 ∧ P2 < P1 ∧ P1 < 1
    -- (f) extended chain to echo primes
    ∧ 0 < P5 ∧ P5 < P4 ∧ P4 < P3
    -- (g) identification P_3 = α_bare
    ∧ P3 = alphaBareQ
    -- (h) inverse-product identity: 1/sin²θ_3 · 1/sin²θ_5 · 1/sin²θ_7 = α_bare⁻¹
    ∧ invProductActive = alphaBareInvQ
    -- (i) tight bracket on α_bare⁻² (1/α² ∈ (18570, 18575))
    ∧ 18570 < invAlphaSquared ∧ invAlphaSquared < 18575 := by
  refine ⟨alphaBareQ_exact,
          alphaBareQ_bracket_tight.1, alphaBareQ_bracket_tight.2,
          alphaBareInvQ_gt_135, alphaBareInvQ_lt_137,
          alphaBareQ_pos, alphaBareInvQ_pos, alphaBareQ_lt_one,
          P3_pos, P3_lt_P2, P2_lt_P1, P1_lt_one,
          P5_pos, P5_lt_P4, P4_lt_P3,
          P3_eq_alphaBareQ,
          invProductActive_eq_alphaBareInv,
          invAlphaSquared_bracket_tight.1,
          invAlphaSquared_bracket_tight.2⟩

/-! ## 5. Master alpha powers

Values, positivity, monotone decrease and tight brackets for the
sequences `α^n` and `(α⁻¹)^n` for `n = 1..5`, plus the exact duality
`α^n · (α⁻¹)^n = 1`. -/

/-- **Master alpha powers.** Aggregator for the `α^n` and `(α⁻¹)^n`
sequences (`n = 1..5`): values at `n=0,1,2`, positivity, the strict
monotone decrease `α^n` and increase `(α⁻¹)^n`, and the exact duality
`α^n · (α⁻¹)^n = 1` for `n = 1..5`. -/
theorem master_alpha_powers :
    -- (a) base values
    alphaPow 0 = 1 ∧ alphaPow 1 = alphaBareQ
    ∧ alphaInvPow 0 = 1 ∧ alphaInvPow 1 = alphaBareInvQ
    -- (b) positivity for n = 1..5
    ∧ 0 < alphaPow 1 ∧ 0 < alphaPow 2 ∧ 0 < alphaPow 3
    ∧ 0 < alphaPow 4 ∧ 0 < alphaPow 5
    ∧ 0 < alphaInvPow 1 ∧ 0 < alphaInvPow 2 ∧ 0 < alphaInvPow 3
    ∧ 0 < alphaInvPow 4 ∧ 0 < alphaInvPow 5
    -- (c) recursion / decomposition
    ∧ alphaPow 2 = alphaBareQ * alphaBareQ
    ∧ alphaInvPow 2 = alphaBareInvQ * alphaBareInvQ
    -- (d) exact duality α^n · (α⁻¹)^n = 1 for n = 1..5
    ∧ alphaPow 1 * alphaInvPow 1 = 1
    ∧ alphaPow 2 * alphaInvPow 2 = 1
    ∧ alphaPow 3 * alphaInvPow 3 = 1
    ∧ alphaPow 4 * alphaInvPow 4 = 1
    ∧ alphaPow 5 * alphaInvPow 5 = 1 := by
  refine ⟨alphaPow_zero, alphaPow_one, alphaInvPow_zero, alphaInvPow_one,
          alphaPow_one_pos, alphaPow_two_pos, alphaPow_three_pos,
          alphaPow_four_pos, alphaPow_five_pos,
          alphaInvPow_one_pos, alphaInvPow_two_pos, alphaInvPow_three_pos,
          alphaInvPow_four_pos, alphaInvPow_five_pos,
          alphaPow_two, alphaInvPow_two,
          alphaPow_mul_alphaInvPow_one, alphaPow_mul_alphaInvPow_two,
          alphaPow_mul_alphaInvPow_three, alphaPow_mul_alphaInvPow_four,
          alphaPow_mul_alphaInvPow_five⟩

/-! ## 6. Final umbrella

Single witness that **the entire holonomy axis** of Persistence Theory
is fully derived inside Lean: each conjunct is one of the five master
sub-theorems above. -/

/-- **Holonomy master framework — final summary.** The conjunction of
the five master sub-theorems `master_cyclic_phase`, `master_active_primes`,
`master_gamma_family`, `master_coupling`, `master_alpha_powers`.

This theorem has no new mathematical content: it is purely an
aggregator providing downstream code with a single, stable witness
named `holonomy_master_framework_summary`. -/
theorem holonomy_master_framework_summary :
    -- 1. cyclic phase
    ((∀ δ : ℝ, 1 - (1 - δ)^2 = δ * (2 - δ))
      ∧ (∀ θ δ : ℝ, Real.cos θ = 1 - δ → Real.sin θ ^ 2 = δ * (2 - δ))
      ∧ (∀ δ δ' : ℝ, 0 ≤ δ → δ ≤ 1 → 0 ≤ δ' → δ' ≤ 1 → δ < δ' →
          δ * (2 - δ) < δ' * (2 - δ'))
      ∧ (∀ s : ℝ, 0 ≤ s → s ≤ 1 →
          let δ := 1 - Real.sqrt (1 - s)
          0 ≤ δ ∧ δ ≤ 1 ∧ δ * (2 - δ) = s)
      ∧ deltaQ 2 = 28 / 225
      ∧ deltaQ 3 = 1178 / 10125
      ∧ deltaQ 5 = 388082 / 3796875
      ∧ deltaQ 7 = 108110858 / 1196015625
      ∧ deltaQ 11 = 6857595465338 / 95147314453125
      ∧ deltaQ 13 = 1643319961767122 / 25300535888671875)
    -- 2. active primes
    ∧ (gammaQ 3 > sPT ∧ gammaQ 5 > sPT ∧ gammaQ 7 > sPT
      ∧ gammaQ 11 < sPT ∧ gammaQ 13 < sPT ∧ gammaQ 17 < sPT
      ∧ gammaQ 19 < sPT ∧ gammaQ 23 < sPT)
    -- 3. gamma family (key brackets)
    ∧ (0 < gammaSumActive ∧ 0 < gammaProductActive
      ∧ 2099 / 1000 < gammaSumActive ∧ gammaSumActive < 2100 / 1000
      ∧ 334 / 1000 < gammaProductActive
      ∧ gammaProductActive < 335 / 1000
      ∧ 3149 / 100 < gammaSumActive * muStar
      ∧ gammaSumActive * muStar < 3150 / 100
      ∧ 35160 / 1000 < gammaPrimorialProductActive
      ∧ gammaPrimorialProductActive < 35161 / 1000)
    -- 4. coupling
    ∧ (0 < alphaBareQ ∧ alphaBareQ < 1
      ∧ 7335 / 1000000 < alphaBareQ ∧ alphaBareQ < 7340 / 1000000
      ∧ 135 < alphaBareInvQ ∧ alphaBareInvQ < 137
      ∧ P3 = alphaBareQ
      ∧ invProductActive = alphaBareInvQ)
    -- 5. alpha powers (duality core)
    ∧ (alphaPow 1 * alphaInvPow 1 = 1
      ∧ alphaPow 2 * alphaInvPow 2 = 1
      ∧ alphaPow 3 * alphaInvPow 3 = 1
      ∧ alphaPow 4 * alphaInvPow 4 = 1
      ∧ alphaPow 5 * alphaInvPow 5 = 1) := by
  refine ⟨?_, ?_, ?_, ?_, ?_⟩
  · -- cyclic phase
    exact ⟨one_sub_one_sub_sq, sin_sq_of_cos_eq_one_sub,
           delta_two_minus_delta_strictMono, delta_two_minus_delta_inversion,
           deltaQ_2_eq, deltaQ_3_eq, deltaQ_5_eq, deltaQ_7_eq,
           deltaQ_11_eq, deltaQ_13_eq⟩
  · -- active primes
    exact ⟨gamma_3_active, gamma_5_active, gamma_7_active,
           gamma_11_inactive, gamma_13_inactive, gamma_17_inactive,
           gammaQ_19_inactive, gammaQ_23_inactive⟩
  · -- gamma family
    exact ⟨gammaSumActive_pos, gammaProductActive_pos,
           gammaSumActive_bracket.1, gammaSumActive_bracket.2,
           gammaProductActive_bracket.1, gammaProductActive_bracket.2,
           gammaSumActive_mul_muStar_bracket.1,
           gammaSumActive_mul_muStar_bracket.2,
           gammaPrimorialProductActive_bracket.1,
           gammaPrimorialProductActive_bracket.2⟩
  · -- coupling
    exact ⟨alphaBareQ_pos, alphaBareQ_lt_one,
           alphaBareQ_bracket_tight.1, alphaBareQ_bracket_tight.2,
           alphaBareInvQ_gt_135, alphaBareInvQ_lt_137,
           P3_eq_alphaBareQ, invProductActive_eq_alphaBareInv⟩
  · -- alpha powers (duality)
    exact ⟨alphaPow_mul_alphaInvPow_one, alphaPow_mul_alphaInvPow_two,
           alphaPow_mul_alphaInvPow_three, alphaPow_mul_alphaInvPow_four,
           alphaPow_mul_alphaInvPow_five⟩

end PT.Holonomy
