/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Holonomy.GammaProduct
import PT.Holonomy.GammaSumActive
import Mathlib.Tactic

/-!
# Joint invariant `Γ_active · ∏ p` over the active cascade `{3, 5, 7}`

**Statement (paper-level, Ch06 §"Cascade arithmétique" / Ch09).**
Where `PT.Holonomy.GammaProduct` records the **anomalous-dimension product**
`Γ_active = γ_3 γ_5 γ_7 ∈ (0.334, 0.335)` and
`PT.Holonomy.GammaSumActive` records the **anomalous-dimension sum**
`Σ_active = γ_3 + γ_5 + γ_7 ∈ (2.099, 2.100)`, this module records the
**joint product** with the *active prime product*

$$P_{\rm active} \;:=\; \prod_{p \in \{3, 5, 7\}} p \;=\; 3 \cdot 5 \cdot 7
   \;=\; 105.$$

The two arithmetic sides — one *spectral* (`Γ`), one *number-theoretic* (`P`)
— combine into

$$\Gamma_{\rm active} \cdot P_{\rm active}
   \;=\; (\gamma_3 \gamma_5 \gamma_7) \cdot (3 \cdot 5 \cdot 7)
   \;\approx\; 0.3348 \cdot 105 \;\approx\; 35.16.$$

This file records four observables:

1. **Definition** `activePrimeProduct = 3·5·7 = 105` and
   `gammaPrimorialProductActive = Γ_active · 105`.
2. **Tight bracket** `(35.160, 35.161)`.
3. **μ⋆-quotient** `(Γ_active · 105) / μ⋆ ∈ (2.344, 2.345)` where `μ⋆ = 15`.
4. **Σ–Π ratio** `(Γ · 105) / (Σ · 15) ∈ (1.116, 1.117)`, the
   dimensionless "geometric vs arithmetic" balance.

Two combinatorial connectors are also recorded:

* **Primorial bridge** `2 · activePrimeProduct = primorial 3 = 210`
  (`PT.Sieve.N3aFactorisation.nat_factorise_primorial_4` reads
  `210 = 2·3·5·7`; we use only the active prime triple here).
* **N3a factorisation bridge** `activePrimeProduct = 105 = 3 · 5 · 7` —
  the explicit decomposition, matching
  `PT.Sieve.N3aFactorisation.nat_factorise_primorial_3` style.

All proofs are pure rational arithmetic — `unfold` + `norm_num`.

## Reference

* Monograph Chapter 6, §"Cascade arithmétique", cross-reference to
  `PT.Holonomy.GammaProduct.gamma_product_summary` and
  `PT.Holonomy.GammaSumActive.gamma_sum_summary`.
* Monograph Chapter 9, §"Reconstruction des couplages" (joint
  spectral × arithmetic invariants).
-/

namespace PT.Holonomy

/-! ### Definition: active prime product `3 · 5 · 7 = 105` -/

/-- The product of the **active primes** `{3, 5, 7}` of the PT cascade,
    promoted to `ℚ`:

    `activePrimeProduct := 3 · 5 · 7 = 105`.

    Companion of the number-theoretic identity `3 · 5 · 7 = 105`
    (`PT.Sieve.N3aFactorisation` style, here on `ℚ`). -/
def activePrimeProduct : ℚ := 3 * 5 * 7

/-- The active prime product is exactly `105`. -/
theorem activePrimeProduct_eq_105 : activePrimeProduct = 105 := by
  unfold activePrimeProduct; norm_num

/-- The active prime product is **positive**. -/
theorem activePrimeProduct_pos : 0 < activePrimeProduct := by
  unfold activePrimeProduct; norm_num

/-- Twice the active prime product equals the **fourth primorial** `210`.
    `2 · (3·5·7) = 210 = primorial 4` (cross-reference to
    `PT.Sieve.N3aFactorisation.nat_factorise_primorial_4 : 210 = 2·3·5·7`). -/
theorem two_mul_activePrimeProduct_eq_primorial4 :
    2 * activePrimeProduct = 210 := by
  unfold activePrimeProduct; norm_num

/-! ### Definition: joint spectral × arithmetic invariant -/

/-- The **joint invariant** `Γ_active · P_active`:

    `gammaPrimorialProductActive := gammaProductActive · activePrimeProduct
                                 = (γ_3 γ_5 γ_7) · (3 · 5 · 7)`.

    Combines the spectral product (anomalous dimensions) with the
    arithmetic product (active primes themselves). -/
def gammaPrimorialProductActive : ℚ :=
  gammaProductActive * activePrimeProduct

/-- **Definitional identity.** The joint invariant unfolds to
    `Γ_active · 105`. -/
theorem gammaPrimorialProductActive_eq_def :
    gammaPrimorialProductActive = gammaProductActive * activePrimeProduct :=
  rfl

/-- **Explicit-105 form.** The joint invariant equals `Γ_active · 105`. -/
theorem gammaPrimorialProductActive_eq_mul_105 :
    gammaPrimorialProductActive = gammaProductActive * 105 := by
  unfold gammaPrimorialProductActive
  rw [activePrimeProduct_eq_105]

/-! ### Positivity -/

/-- **Positivity** `0 < Γ_active · 105`. Follows from
    `0 < Γ_active` and `0 < 105`. -/
theorem gammaPrimorialProductActive_pos :
    0 < gammaPrimorialProductActive := by
  unfold gammaPrimorialProductActive
  exact mul_pos gammaProductActive_pos activePrimeProduct_pos

/-! ### Tight rational bracket `(35.160, 35.161)` -/

/-- **Tight bracket for the joint invariant.**
    `35.160 < Γ_active · 105 < 35.161`. Numerical value `≈ 35.1608`. -/
theorem gammaPrimorialProductActive_bracket :
    35160 / 1000 < gammaPrimorialProductActive
    ∧ gammaPrimorialProductActive < 35161 / 1000 := by
  unfold gammaPrimorialProductActive activePrimeProduct
        gammaProductActive gammaQ deltaQ qPT muStar
  refine ⟨?_, ?_⟩ <;> norm_num

/-- **Looser headline bracket.** `35 < Γ_active · 105 < 36`. -/
theorem gammaPrimorialProductActive_loose_bracket :
    35 < gammaPrimorialProductActive
    ∧ gammaPrimorialProductActive < 36 := by
  refine ⟨?_, ?_⟩
  · exact lt_trans (by norm_num : (35 : ℚ) < 35160 / 1000)
      gammaPrimorialProductActive_bracket.1
  · exact lt_trans gammaPrimorialProductActive_bracket.2
      (by norm_num : (35161 / 1000 : ℚ) < 36)

/-! ### μ⋆-quotient: `(Γ · 105) / 15 ∈ (2.344, 2.345)` -/

/-- The **μ⋆-quotient** of the joint invariant. -/
def gammaPrimorialOverMuStar : ℚ := gammaPrimorialProductActive / muStar

/-- The μ⋆-quotient unfolds to `Γ_active · 105 / 15`. -/
theorem gammaPrimorialOverMuStar_eq :
    gammaPrimorialOverMuStar = gammaProductActive * 105 / 15 := by
  unfold gammaPrimorialOverMuStar gammaPrimorialProductActive
        activePrimeProduct muStar
  ring

/-- **Tight bracket for the μ⋆-quotient.**
    `2.344 < (Γ_active · 105) / μ⋆ < 2.345`. Numerical value `≈ 2.3441`. -/
theorem gammaPrimorialOverMuStar_bracket :
    2344 / 1000 < gammaPrimorialOverMuStar
    ∧ gammaPrimorialOverMuStar < 2345 / 1000 := by
  unfold gammaPrimorialOverMuStar gammaPrimorialProductActive
        activePrimeProduct gammaProductActive gammaQ deltaQ qPT muStar
  refine ⟨?_, ?_⟩ <;> norm_num

/-- **Positivity of the μ⋆-quotient.** -/
theorem gammaPrimorialOverMuStar_pos : 0 < gammaPrimorialOverMuStar :=
  lt_trans (by norm_num : (0 : ℚ) < 2344 / 1000)
    gammaPrimorialOverMuStar_bracket.1

/-! ### Σ–Π ratio: `(Γ · 105) / (Σ · 15) ∈ (1.116, 1.117)` -/

/-- The **Σ–Π ratio** comparing the spectral product (`Γ · P`) to the
    spectral sum (`Σ · μ⋆`):

    `gammaPrimorialSumRatio := (Γ_active · 105) / (Σ_active · 15)`. -/
def gammaPrimorialSumRatio : ℚ :=
  gammaPrimorialProductActive / (gammaSumActive * muStar)

/-- **Tight bracket for the Σ–Π ratio.**
    `1.116 < (Γ · 105) / (Σ · 15) < 1.117`. Numerical value `≈ 1.1165`. -/
theorem gammaPrimorialSumRatio_bracket :
    1116 / 1000 < gammaPrimorialSumRatio
    ∧ gammaPrimorialSumRatio < 1117 / 1000 := by
  unfold gammaPrimorialSumRatio gammaPrimorialProductActive
        activePrimeProduct gammaProductActive gammaSumActive
        gammaQ deltaQ qPT muStar
  refine ⟨?_, ?_⟩ <;> norm_num

/-- **Σ–Π ratio is positive.** -/
theorem gammaPrimorialSumRatio_pos : 0 < gammaPrimorialSumRatio :=
  lt_trans (by norm_num : (0 : ℚ) < 1116 / 1000)
    gammaPrimorialSumRatio_bracket.1

/-- **Σ–Π ratio is greater than `1`.** The geometric side
    `Γ · P_active` strictly exceeds the spectral side `Σ · μ⋆`. -/
theorem gammaPrimorialSumRatio_gt_one : 1 < gammaPrimorialSumRatio :=
  lt_trans (by norm_num : (1 : ℚ) < 1116 / 1000)
    gammaPrimorialSumRatio_bracket.1

/-- **Σ–Π ratio is less than `2`.** -/
theorem gammaPrimorialSumRatio_lt_two : gammaPrimorialSumRatio < 2 :=
  lt_trans gammaPrimorialSumRatio_bracket.2
    (by norm_num : (1117 / 1000 : ℚ) < 2)

/-! ### Algebraic identity Γ·P = sum of triple products

The joint invariant expands as a fully symmetric algebraic combination of
`γ_p` with the active prime weights. -/

/-- **Algebraic identity (expanded form).** The joint invariant equals
    the sum over the diagonal `γ_p γ_q γ_r · p q r`, here with the single
    triple `(3, 5, 7)`. -/
theorem gammaPrimorialProductActive_expanded :
    gammaPrimorialProductActive
      = gammaQ 3 * gammaQ 5 * gammaQ 7 * (3 * 5 * 7) := by
  unfold gammaPrimorialProductActive activePrimeProduct gammaProductActive
  ring

/-- **Algebraic identity (factored form by prime).** Each spectral
    `γ_p` is paired with its arithmetic prime `p`:

    `(γ_3 · 3) · (γ_5 · 5) · (γ_7 · 7) = Γ · 105`. -/
theorem gammaPrimorialProductActive_factored :
    gammaPrimorialProductActive
      = (gammaQ 3 * 3) * (gammaQ 5 * 5) * (gammaQ 7 * 7) := by
  unfold gammaPrimorialProductActive activePrimeProduct gammaProductActive
  ring

/-! ### Headline -/

/-- **Headline (γ × active-prime joint invariant).** At `q_+ = 13/15`,
    `μ⋆ = 15`, on the active cascade `{3, 5, 7}`:

    1. `activePrimeProduct = 3 · 5 · 7 = 105` is positive, and
       `2 · activePrimeProduct = 210` matches the fourth primorial.
    2. The joint invariant
       `Γ_active · P_active = gammaPrimorialProductActive`
       is positive and lies in the tight bracket `(35.160, 35.161)`.
    3. The μ⋆-quotient `(Γ_active · 105) / μ⋆` lies in the tight bracket
       `(2.344, 2.345)`.
    4. The Σ–Π ratio `(Γ · 105) / (Σ · 15)` lies in `(1.116, 1.117)`
       and is strictly between `1` and `2`. -/
theorem gamma_primorial_product_summary :
    activePrimeProduct = 105
    ∧ 0 < activePrimeProduct
    ∧ 2 * activePrimeProduct = 210
    ∧ 0 < gammaPrimorialProductActive
    ∧ 35160 / 1000 < gammaPrimorialProductActive
    ∧ gammaPrimorialProductActive < 35161 / 1000
    ∧ 35 < gammaPrimorialProductActive
    ∧ gammaPrimorialProductActive < 36
    ∧ 0 < gammaPrimorialOverMuStar
    ∧ 2344 / 1000 < gammaPrimorialOverMuStar
    ∧ gammaPrimorialOverMuStar < 2345 / 1000
    ∧ 0 < gammaPrimorialSumRatio
    ∧ 1116 / 1000 < gammaPrimorialSumRatio
    ∧ gammaPrimorialSumRatio < 1117 / 1000
    ∧ 1 < gammaPrimorialSumRatio
    ∧ gammaPrimorialSumRatio < 2 :=
  ⟨activePrimeProduct_eq_105,
   activePrimeProduct_pos,
   two_mul_activePrimeProduct_eq_primorial4,
   gammaPrimorialProductActive_pos,
   gammaPrimorialProductActive_bracket.1,
   gammaPrimorialProductActive_bracket.2,
   gammaPrimorialProductActive_loose_bracket.1,
   gammaPrimorialProductActive_loose_bracket.2,
   gammaPrimorialOverMuStar_pos,
   gammaPrimorialOverMuStar_bracket.1,
   gammaPrimorialOverMuStar_bracket.2,
   gammaPrimorialSumRatio_pos,
   gammaPrimorialSumRatio_bracket.1,
   gammaPrimorialSumRatio_bracket.2,
   gammaPrimorialSumRatio_gt_one,
   gammaPrimorialSumRatio_lt_two⟩

end PT.Holonomy
