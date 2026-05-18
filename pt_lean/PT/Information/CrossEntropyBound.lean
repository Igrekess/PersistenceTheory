/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Information.GFTIdentity
import PT.Information.GFTSpecialisations
import PT.Information.BekensteinExtensions
import PT.Information.GFTSpecialMValues
import Mathlib.Tactic

/-!
# Cross-entropy `H(p, q) := -∑ p log q` — algebraic identities and bounds

This file introduces the **cross-entropy** of a (finite-support) probability
distribution `p` relative to a strictly positive reference distribution `q`:

$$H(p, q) \;:=\; -\sum_{r \in s} p_r \, \log q_r.$$

We collect the unconditional algebraic facts that follow from this definition:

* **Self-identity.** `H(p, p) = H(p)` (the Shannon entropy of `p`).
* **Uniform reference.** When `q ≡ 1/m`, `H(p, U_m) = log m` for any
  probability vector `p` of total mass `1`. In particular for `p = δ_a`
  one again gets `log m`.
* **Bridge to KL.** `D_{KL}(p \,\Vert\, U_m) = H(p, U_m) - H(p) = \log m - H(p)`.
  This is just a re-packaging of `GFT_identity`.
* **Parametric Gibbs.** Under the *hypothesis* that the KL divergence is
  non-negative (which the literature establishes via convexity of `log`),
  we deduce `H(p, q) ≥ H(p)`. We do **not** prove Gibbs internally — it
  is supplied as a parameter, keeping the file purely algebraic.
* **PT specialisations.** Instances at `m ∈ {2, 8, 30}` (the trivial
  binary partition, the admissible residue cardinality, and the third
  primorial).

## References

Monograph Ch04 §4.3 (cross-entropy & uniform reference); follow-up to
`PT.Information.GFTIdentity` and `PT.Information.BekensteinExtensions`.
-/

namespace PT.Information

open Real Finset

/-! ## Definition -/

/-- **Cross-entropy** of a distribution `p` relative to a reference `q`
    on a finite support `s`:
    `H(p, q) := -∑_{r ∈ s} p_r * log q_r`. -/
noncomputable def crossEntropy {α : Type*} (s : Finset α) (p q : α → ℝ) : ℝ :=
  -(∑ r ∈ s, p r * Real.log (q r))

/-! ## Self-cross = Shannon entropy -/

/-- **Self-identity.** `H(p, p) = H(p)`. The proof is termwise:
    `-(p_r * log p_r) = negMulLog p_r`. -/
theorem crossEntropy_self {α : Type*} (s : Finset α) (p : α → ℝ) :
    crossEntropy s p p = shannonH s p := by
  unfold crossEntropy shannonH
  -- Rewrite negMulLog p_r = -(p_r * log p_r) and pull the negation across the sum.
  rw [← Finset.sum_neg_distrib]
  refine Finset.sum_congr rfl ?_
  intro r _
  rw [Real.negMulLog]
  ring

/-! ## Uniform reference: `crossEntropy s p (1/m) = log m` -/

/-- **Cross-entropy against the uniform reference** is `log m`, for any
    probability vector `p` of total mass `1`. -/
theorem crossEntropy_uniform {α : Type*} (s : Finset α) (m : ℝ) (_hm : 0 < m)
    (p : α → ℝ) (hp_sum : ∑ r ∈ s, p r = 1) :
    crossEntropy s p (fun _ => (1 : ℝ) / m) = Real.log m := by
  unfold crossEntropy
  -- log (1/m) = -log m
  have hlog : Real.log ((1 : ℝ) / m) = -Real.log m := by
    rw [one_div, Real.log_inv]
  -- ∑ p_r * log(1/m) = (∑ p_r) * log(1/m) = log(1/m)
  rw [← Finset.sum_mul, hp_sum, one_mul, hlog]
  ring

/-- **Cross-entropy of a delta against the uniform** is `log m`. -/
theorem crossEntropy_delta_uniform {α : Type*} [DecidableEq α]
    (s : Finset α) (m : ℝ) (hm : 0 < m) (r₀ : α) (hr : r₀ ∈ s) :
    crossEntropy s (deltaDist r₀) (fun _ => (1 : ℝ) / m) = Real.log m :=
  crossEntropy_uniform s m hm (deltaDist r₀) (deltaDist_sum s r₀ hr)

/-! ## GFT identity, repackaged through cross-entropy -/

/-- **Bridge to KL via cross-entropy.** For a probability vector `p`,
    `D_KL(p ‖ U_m) = H(p, U_m) - H(p) = log m - H(p)`.

    This is just a re-arrangement of `GFT_identity`. -/
theorem klToUniform_eq_crossEntropy_sub_shannonH
    {α : Type*} (s : Finset α) (m : ℝ) (hm : 0 < m)
    (p : α → ℝ) (hp_nonneg : ∀ r ∈ s, 0 ≤ p r)
    (hp_sum : ∑ r ∈ s, p r = 1) :
    klToUniform s m p
      = crossEntropy s p (fun _ => (1 : ℝ) / m) - shannonH s p := by
  have hcross := crossEntropy_uniform s m hm p hp_sum
  have hid := GFT_identity s m hm p hp_nonneg hp_sum
  -- klToUniform + shannonH = log m  and  crossEntropy = log m
  linarith

/-- **Bridge to KL, explicit `log m - H(p)` form.** -/
theorem klToUniform_eq_logm_sub_shannonH
    {α : Type*} (s : Finset α) (m : ℝ) (hm : 0 < m)
    (p : α → ℝ) (hp_nonneg : ∀ r ∈ s, 0 ≤ p r)
    (hp_sum : ∑ r ∈ s, p r = 1) :
    klToUniform s m p = Real.log m - shannonH s p := by
  have hid := GFT_identity s m hm p hp_nonneg hp_sum
  linarith

/-! ## Parametric Gibbs inequality

We do **not** prove `D_KL ≥ 0` inside this file (that requires convexity
of `log` and is supplied elsewhere). The statement below shows that, as
soon as one *hypothesises* `D_KL(p ‖ q) ≥ 0`, the cross-entropy bound
`H(p, q) ≥ H(p)` is immediate for the uniform reference. -/

/-- **Parametric Gibbs (uniform reference).** Assuming the KL divergence
    of `p` to the uniform `U_m` is non-negative, the cross-entropy
    `H(p, U_m)` is at least the Shannon entropy `H(p)`. -/
theorem crossEntropy_uniform_ge_shannonH_of_kl_nonneg
    {α : Type*} (s : Finset α) (m : ℝ) (hm : 0 < m)
    (p : α → ℝ) (hp_nonneg : ∀ r ∈ s, 0 ≤ p r)
    (hp_sum : ∑ r ∈ s, p r = 1)
    (h_kl_nn : 0 ≤ klToUniform s m p) :
    shannonH s p ≤ crossEntropy s p (fun _ => (1 : ℝ) / m) := by
  have hcross := crossEntropy_uniform s m hm p hp_sum
  have hid := GFT_identity s m hm p hp_nonneg hp_sum
  -- crossEntropy = log m, and klToUniform + shannonH = log m with klToUniform ≥ 0
  linarith

/-! ## PT specialisations: `m ∈ {2, 8, 30}` -/

/-- **Cross-entropy at `m = 2`, uniform reference, delta input.** -/
theorem crossEntropy_m2_delta_uniform :
    crossEntropy (Finset.univ : Finset (Fin 2)) (deltaDist 0)
      (fun _ => (1 : ℝ) / 2) = Real.log 2 :=
  crossEntropy_delta_uniform _ 2 (by norm_num) 0 (by decide)

/-- **Cross-entropy at `m = 8`, uniform reference, delta input.** -/
theorem crossEntropy_m8_delta_uniform :
    crossEntropy (Finset.univ : Finset (Fin 8)) (deltaDist 0)
      (fun _ => (1 : ℝ) / 8) = Real.log 8 :=
  crossEntropy_delta_uniform _ 8 (by norm_num) 0 (by decide)

/-- **Cross-entropy at `m = 30`, uniform reference, delta input.** -/
theorem crossEntropy_m30_delta_uniform :
    crossEntropy (Finset.univ : Finset (Fin 30)) (deltaDist 0)
      (fun _ => (1 : ℝ) / 30) = Real.log 30 :=
  crossEntropy_delta_uniform _ 30 (by norm_num) 0 (by decide)

/-- **Cross-entropy at `m = 2`, uniform-vs-uniform** equals `log 2`. -/
theorem crossEntropy_m2_uniform_uniform :
    crossEntropy (Finset.univ : Finset (Fin 2))
      (fun _ => (1 : ℝ) / 2) (fun _ => (1 : ℝ) / 2) = Real.log 2 := by
  apply crossEntropy_uniform _ 2 (by norm_num)
  simp [Finset.sum_const, Finset.card_univ, Fintype.card_fin]

/-- **Cross-entropy at `m = 8`, uniform-vs-uniform** equals `log 8`. -/
theorem crossEntropy_m8_uniform_uniform :
    crossEntropy (Finset.univ : Finset (Fin 8))
      (fun _ => (1 : ℝ) / 8) (fun _ => (1 : ℝ) / 8) = Real.log 8 := by
  apply crossEntropy_uniform _ 8 (by norm_num)
  simp [Finset.sum_const, Finset.card_univ, Fintype.card_fin]

/-- **Cross-entropy at `m = 30`, uniform-vs-uniform** equals `log 30`. -/
theorem crossEntropy_m30_uniform_uniform :
    crossEntropy (Finset.univ : Finset (Fin 30))
      (fun _ => (1 : ℝ) / 30) (fun _ => (1 : ℝ) / 30) = Real.log 30 := by
  apply crossEntropy_uniform _ 30 (by norm_num)
  simp [Finset.sum_const, Finset.card_univ, Fintype.card_fin]

/-! ## Headline summary -/

/-- **Headline.** Algebraic identities of cross-entropy w.r.t. the
    uniform reference:

    * `H(p, p) = H(p)` (self-cross is the Shannon entropy).
    * `H(p, U_m) = log m` for any probability vector `p` (in particular
      both for `p = δ_a` and for `p = U_m`).
    * `D_KL(p ‖ U_m) = H(p, U_m) - H(p) = log m - H(p)`. -/
theorem crossEntropy_summary {α : Type*} [DecidableEq α]
    (s : Finset α) (m : ℝ) (hm : 0 < m) (r₀ : α) (hr : r₀ ∈ s)
    (p : α → ℝ) (hp_nonneg : ∀ r ∈ s, 0 ≤ p r)
    (hp_sum : ∑ r ∈ s, p r = 1) :
    -- self-cross
    crossEntropy s p p = shannonH s p
    -- uniform reference, any p
    ∧ crossEntropy s p (fun _ => (1 : ℝ) / m) = Real.log m
    -- delta corner
    ∧ crossEntropy s (deltaDist r₀) (fun _ => (1 : ℝ) / m) = Real.log m
    -- bridge to KL
    ∧ klToUniform s m p = Real.log m - shannonH s p :=
  ⟨crossEntropy_self s p,
   crossEntropy_uniform s m hm p hp_sum,
   crossEntropy_delta_uniform s m hm r₀ hr,
   klToUniform_eq_logm_sub_shannonH s m hm p hp_nonneg hp_sum⟩

end PT.Information
