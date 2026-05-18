/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import Mathlib.Analysis.SpecialFunctions.Log.NegMulLog
import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.Tactic

/-!
# GFT Fundamental Identity — `log m = D_KL(P ‖ U_m) + H(P)`

**Statement (paper-level, Ch04 §4.2).** Let `P = (p_r)_{r ∈ R}` be a discrete
probability distribution on a finite set `R` of cardinality `m`, and let
`U_m` denote the uniform distribution on `R`. Then

$$\log m \;=\; D_{KL}(P \,\Vert\, U_m) \;+\; H(P),$$

an *exact* algebraic identity (residue identically zero). Here
`D_KL(P ‖ U_m) = ∑_r p_r log(m · p_r)`, `H(P) = -∑_r p_r log p_r`, and the
identity is base-independent — we state it in the natural log.

Proof: termwise expansion
`p_r log(m p_r) = p_r log m + p_r log p_r`, summed and using `∑ p_r = 1`.

We work in `Finset`s for compactness; the choice of base (`log₂` vs. `log`)
only multiplies both sides by a positive constant.

## Reference

Monograph Chapter 4, §4.2, `\label{thm:GFT-ID}`. M4 article.
-/

namespace PT.Information

open Real Finset

/-- The Shannon entropy on a finite distribution `p : α → ℝ`
    (in natural log). -/
noncomputable def shannonH {α : Type*} (s : Finset α) (p : α → ℝ) : ℝ :=
  ∑ r ∈ s, Real.negMulLog (p r)

/-- The (unnormalised) KL divergence to the uniform distribution on a
    finite set of cardinality `m`: `∑_r p_r * log (m * p_r)`. -/
noncomputable def klToUniform {α : Type*} (s : Finset α) (m : ℝ) (p : α → ℝ) : ℝ :=
  ∑ r ∈ s, p r * Real.log (m * p r)

/-- **GFT identity (algebraic kernel, termwise form).**
    For `m > 0`, `0 < p_r`, the equality
    `p_r · log(m · p_r) = p_r · log m + p_r · log p_r` holds. -/
lemma kl_pointwise {m p : ℝ} (hm : 0 < m) (hp : 0 < p) :
    p * Real.log (m * p) = p * Real.log m + p * Real.log p := by
  rw [Real.log_mul (ne_of_gt hm) (ne_of_gt hp)]
  ring

/-- **GFT identity (algebraic kernel, termwise form, allows `p = 0`).**
    For `m > 0`, `0 ≤ p`, the equality
    `p · log(m · p) = p · log m + p · log p` holds (with the convention
    `0 * log 0 = 0`). -/
lemma kl_pointwise_nonneg {m p : ℝ} (hm : 0 < m) (hp : 0 ≤ p) :
    p * Real.log (m * p) = p * Real.log m + p * Real.log p := by
  rcases eq_or_lt_of_le hp with hp0 | hp_pos
  · -- p = 0: both sides are 0 (via 0 * _ = 0).
    have : p = 0 := hp0.symm
    rw [this]; simp
  · exact kl_pointwise hm hp_pos

/-- **GFT Fundamental Identity (exact algebraic partition).**
    Let `s : Finset α` enumerate the residue classes of cardinality `m`,
    let `p : α → ℝ` be a probability vector (`∑ p_r = 1`, `p_r ≥ 0`).
    Then for any `m > 0`,

    `klToUniform s m p + shannonH s p = log m`. -/
theorem GFT_identity {α : Type*} (s : Finset α) (m : ℝ) (hm : 0 < m)
    (p : α → ℝ) (hp_nonneg : ∀ r ∈ s, 0 ≤ p r)
    (hp_sum : ∑ r ∈ s, p r = 1) :
    klToUniform s m p + shannonH s p = Real.log m := by
  unfold klToUniform shannonH
  -- termwise: p_r log(m p_r) + (- p_r log p_r) = p_r log m
  have hterm : ∀ r ∈ s,
      p r * Real.log (m * p r) + Real.negMulLog (p r) = p r * Real.log m := by
    intro r hr
    have hp_r : 0 ≤ p r := hp_nonneg r hr
    have := kl_pointwise_nonneg hm hp_r
    unfold Real.negMulLog
    linarith
  rw [← Finset.sum_add_distrib]
  rw [Finset.sum_congr rfl hterm]
  -- ∑ p_r * log m = log m · ∑ p_r = log m
  rw [← Finset.sum_mul, hp_sum, one_mul]

/-- **GFT Identity, rewritten as a partition `log m = D_KL + H`.** -/
theorem GFT_partition {α : Type*} (s : Finset α) (m : ℝ) (hm : 0 < m)
    (p : α → ℝ) (hp_nonneg : ∀ r ∈ s, 0 ≤ p r)
    (hp_sum : ∑ r ∈ s, p r = 1) :
    Real.log m = klToUniform s m p + shannonH s p :=
  (GFT_identity s m hm p hp_nonneg hp_sum).symm

/-- **GFT Identity, specialised to a natural-number-valued modulus `m ≥ 1`.** -/
theorem GFT_identity_nat {α : Type*} (s : Finset α) (m : ℕ) (hm : 1 ≤ m)
    (p : α → ℝ) (hp_nonneg : ∀ r ∈ s, 0 ≤ p r)
    (hp_sum : ∑ r ∈ s, p r = 1) :
    klToUniform s (m : ℝ) p + shannonH s p = Real.log (m : ℝ) :=
  GFT_identity s (m : ℝ) (by exact_mod_cast Nat.lt_of_lt_of_le Nat.zero_lt_one hm) p hp_nonneg hp_sum

end PT.Information
