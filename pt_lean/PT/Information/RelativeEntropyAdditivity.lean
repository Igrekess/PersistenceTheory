/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Information.GFTIdentity
import PT.Information.GFTSpecialisations
import PT.Information.BekensteinExtensions
import PT.Information.MutualInformationBasic
import PT.Information.KLAdditivityProduct
import Mathlib.Algebra.BigOperators.Group.Finset.Sigma
import Mathlib.Tactic

/-!
# Additivity of relative entropy on tensor products (Ch04 §"Additivité KL")

The standard KL divergence to the uniform is **additive** on independent
product distributions. Concretely, if `P = (p_a)_{a ∈ s}` and
`Q = (q_b)_{b ∈ t}` are two probability vectors and `U_m, U_n` are the
uniform distributions on `s, t` of cardinalities `m, n`, then

  `D_KL(P ⊗ Q ‖ U_m ⊗ U_n) = D_KL(P ‖ U_m) + D_KL(Q ‖ U_n).`

Equivalently in our `klToUniform` API (where `U_m ⊗ U_n` is the uniform
on `s ×ˢ t` of cardinality `m · n`, see
`MutualInformationBasic.joint_uniform_eq_uniform_prod`):

  `klToUniform (s ×ˢ t) (m · n) (joint p q)
      = klToUniform s m p + klToUniform t n q.`

The pointwise algebraic kernel is

  `(p·q) · log((m·n)·(p·q)) = (p·q) · (log(m·p) + log(n·q))
                            = q · (p · log(m·p)) + p · (q · log(n·q)),`

and summing over `s ×ˢ t` collapses one of the two indices via
`∑ p = ∑ q = 1`.

We then specialise:

* **Delta corner.** `D_KL(δ ⊗ δ ‖ U_m ⊗ U_n) = log m + log n = log(m·n)`
  (matches `KLAdditivityProduct.klDelta_add_klDelta_eq_log_mul`).
* **Uniform corner.** `D_KL(U_m ⊗ U_n ‖ U_m ⊗ U_n) = 0`.
* **PT-canonical numerics.** `(m, n) = (2, 8)` gives `log 16`,
  `(m, n) = (8, 30)` gives `log 240`.

## Reference

Monograph Chapter 4, `\label{prop:KL_additivity}`. M4 article.
-/

namespace PT.Information

open Real Finset

/-! ### Pointwise algebraic kernel on the joint -/

/-- **Joint pointwise kernel.** For `m, n > 0` and `0 ≤ p, q`,
    `(p · q) · log((m·n) · (p · q))
        = q · (p · log(m · p)) + p · (q · log(n · q))`. -/
lemma kl_joint_pointwise_nonneg
    {m n p q : ℝ} (hm : 0 < m) (hn : 0 < n)
    (hp : 0 ≤ p) (hq : 0 ≤ q) :
    (p * q) * Real.log ((m * n) * (p * q))
      = q * (p * Real.log (m * p)) + p * (q * Real.log (n * q)) := by
  -- Use the existing `kl_pointwise_nonneg` lemmas: for each side,
  -- `p · log(m·p) = p·log m + p·log p`, etc.
  have hmn_pos : 0 < m * n := mul_pos hm hn
  have hpq_nn : 0 ≤ p * q := mul_nonneg hp hq
  have h_lhs := kl_pointwise_nonneg hmn_pos hpq_nn
  -- h_lhs : (p*q) * log((m*n) * (p*q)) = (p*q) * log(m*n) + (p*q) * log(p*q)
  have h_p := kl_pointwise_nonneg hm hp
  -- h_p : p * log(m * p) = p * log m + p * log p
  have h_q := kl_pointwise_nonneg hn hq
  -- h_q : q * log(n * q) = q * log n + q * log q
  -- Expand log(m*n) and log(p*q) when possible.
  rw [h_lhs]
  rw [Real.log_mul (ne_of_gt hm) (ne_of_gt hn)]
  -- Now we need (p*q) * log(p*q) = q * (p * log p) + p * (q * log q),
  -- and we handle zero cases by splitting on whether p = 0 or q = 0.
  -- We pull RHS apart using h_p and h_q.
  rw [h_p, h_q]
  -- LHS now: (p*q) * (log m + log n) + (p*q) * log(p*q)
  -- RHS now: q * (p * log m + p * log p) + p * (q * log n + q * log q)
  rcases eq_or_lt_of_le hp with hp0 | hp_pos
  · -- p = 0
    have : p = 0 := hp0.symm
    subst this
    simp
  · rcases eq_or_lt_of_le hq with hq0 | hq_pos
    · -- q = 0
      have : q = 0 := hq0.symm
      subst this
      simp
    · -- p > 0 and q > 0: log(p*q) = log p + log q
      rw [Real.log_mul (ne_of_gt hp_pos) (ne_of_gt hq_pos)]
      ring

/-! ### Additivity of `klToUniform` on the joint -/

/-- **KL additivity on tensor products.**
    For probability vectors `p, q` (non-negative, summing to `1`) and
    moduli `m, n > 0`:

      `klToUniform (s ×ˢ t) (m · n) (joint p q)
          = klToUniform s m p + klToUniform t n q.`

    This is the relative-entropy analogue of `shannonH_joint`. -/
theorem klToUniform_joint
    {α β : Type*} (s : Finset α) (t : Finset β)
    (m n : ℝ) (hm : 0 < m) (hn : 0 < n)
    (p : α → ℝ) (q : β → ℝ)
    (hp_nonneg : ∀ a ∈ s, 0 ≤ p a) (hq_nonneg : ∀ b ∈ t, 0 ≤ q b)
    (hp_sum : ∑ a ∈ s, p a = 1) (hq_sum : ∑ b ∈ t, q b = 1) :
    klToUniform (s ×ˢ t) (m * n) (joint p q)
      = klToUniform s m p + klToUniform t n q := by
  unfold klToUniform joint
  -- ∑_{ab ∈ s ×ˢ t} (p ab.1 * q ab.2) * log((m·n) * (p ab.1 * q ab.2))
  --   = ∑_a ∑_b (...)
  rw [Finset.sum_product]
  -- Rewrite each inner term using the joint pointwise kernel
  have h_inner : ∀ a ∈ s,
      ∑ b ∈ t, (p a * q b) * Real.log ((m * n) * (p a * q b))
        = (∑ b ∈ t, q b) * (p a * Real.log (m * p a))
          + p a * (∑ b ∈ t, q b * Real.log (n * q b)) := by
    intro a ha
    have hp_a : 0 ≤ p a := hp_nonneg a ha
    -- termwise expansion
    have h_term : ∀ b ∈ t,
        (p a * q b) * Real.log ((m * n) * (p a * q b))
          = q b * (p a * Real.log (m * p a))
            + p a * (q b * Real.log (n * q b)) := by
      intro b hb
      exact kl_joint_pointwise_nonneg hm hn hp_a (hq_nonneg b hb)
    rw [Finset.sum_congr rfl h_term]
    rw [Finset.sum_add_distrib]
    congr 1
    · -- ∑_b q b * (p a * log (m * p a)) = (∑ q) * (p a * log (m * p a))
      rw [← Finset.sum_mul]
    · -- ∑_b p a * (q b * log (n * q b)) = p a * ∑_b (q b * log (n * q b))
      rw [← Finset.mul_sum]
  rw [Finset.sum_congr rfl h_inner]
  -- ∑_a [ (∑ q) * (p a * log (m * p a)) + p a * (∑_b q b * log (n * q b)) ]
  rw [Finset.sum_add_distrib]
  -- First piece: (∑ q) * ∑_a p a * log (m * p a)
  -- Second piece: (∑_a p a) * ∑_b q b * log (n * q b)
  rw [hq_sum]
  simp_rw [one_mul]
  rw [← Finset.sum_mul, hp_sum, one_mul]

/-! ### Delta corner -/

/-- **Delta corner additivity** (re-derived through `klToUniform_joint`).
    `D_KL(δ_(a₀, b₀) ‖ U_m ⊗ U_n) = log m + log n = log(m · n)`. -/
theorem klToUniform_joint_delta
    {α β : Type*} [DecidableEq α] [DecidableEq β]
    (s : Finset α) (t : Finset β)
    (m n : ℝ) (hm : 0 < m) (hn : 0 < n)
    (a₀ : α) (ha : a₀ ∈ s) (b₀ : β) (hb : b₀ ∈ t) :
    klToUniform (s ×ˢ t) (m * n) (joint (deltaDist a₀) (deltaDist b₀))
      = Real.log m + Real.log n := by
  rw [klToUniform_joint s t m n hm hn (deltaDist a₀) (deltaDist b₀)
        (fun r _ => deltaDist_nonneg a₀ r)
        (fun r _ => deltaDist_nonneg b₀ r)
        (deltaDist_sum s a₀ ha) (deltaDist_sum t b₀ hb)]
  rw [klToUniform_delta s m hm a₀ ha, klToUniform_delta t n hn b₀ hb]

/-- **Delta corner additivity, log-product form.** Matches
    `klDelta_add_klDelta_eq_log_mul` and is equal to `log(m·n)`. -/
theorem klToUniform_joint_delta_eq_log_mul
    {α β : Type*} [DecidableEq α] [DecidableEq β]
    (s : Finset α) (t : Finset β)
    (m n : ℝ) (hm : 0 < m) (hn : 0 < n)
    (a₀ : α) (ha : a₀ ∈ s) (b₀ : β) (hb : b₀ ∈ t) :
    klToUniform (s ×ˢ t) (m * n) (joint (deltaDist a₀) (deltaDist b₀))
      = Real.log (m * n) := by
  rw [klToUniform_joint_delta s t m n hm hn a₀ ha b₀ hb]
  exact (Real.log_mul (ne_of_gt hm) (ne_of_gt hn)).symm

/-! ### Uniform corner -/

/-- **Uniform corner additivity.** Both sides are zero:
    `D_KL(U_m ⊗ U_n ‖ U_m ⊗ U_n) = D_KL(U_m ‖ U_m) + D_KL(U_n ‖ U_n) = 0`. -/
theorem klToUniform_joint_uniform
    {α β : Type*} (s : Finset α) (t : Finset β)
    (m n : ℝ) (hm : 0 < m) (hn : 0 < n)
    (hcard_s : (s.card : ℝ) = m) (hcard_t : (t.card : ℝ) = n) :
    klToUniform (s ×ˢ t) (m * n)
        (joint (fun _ : α => (1 : ℝ) / m) (fun _ : β => (1 : ℝ) / n)) = 0 := by
  -- Use additivity, then the fact that each uniform-to-uniform KL is 0.
  have hp_sum : ∑ a ∈ s, ((1 : ℝ) / m) = 1 := by
    rw [Finset.sum_const, nsmul_eq_mul, hcard_s, mul_one_div]
    exact div_self (ne_of_gt hm)
  have hq_sum : ∑ b ∈ t, ((1 : ℝ) / n) = 1 := by
    rw [Finset.sum_const, nsmul_eq_mul, hcard_t, mul_one_div]
    exact div_self (ne_of_gt hn)
  rw [klToUniform_joint s t m n hm hn (fun _ => 1 / m) (fun _ => 1 / n)
        (fun _ _ => by positivity) (fun _ _ => by positivity)
        hp_sum hq_sum]
  rw [klToUniform_uniform s m hm, klToUniform_uniform t n hn]
  ring

/-! ### Cross-check with `MutualInformationBasic.joint_uniform_eq_uniform_prod` -/

/-- **Uniform corner, "true" uniform form.** Plugging the actual uniform
    `1/(m·n)` on `s ×ˢ t` into `klToUniform` gives `0`, consistently with
    the additive form above. -/
theorem klToUniform_uniform_prod
    {α β : Type*} (s : Finset α) (t : Finset β)
    (m n : ℝ) (hm : 0 < m) (hn : 0 < n) :
    klToUniform (s ×ˢ t) (m * n) (fun _ : α × β => (1 : ℝ) / (m * n)) = 0 :=
  klToUniform_uniform (s ×ˢ t) (m * n) (mul_pos hm hn)

/-! ### PT-canonical numerical instances -/

/-- **(m, n) = (2, 8).** Delta-corner additivity yields
    `log 2 + log 8 = log 16`. -/
theorem klToUniform_joint_delta_2_8 :
    klToUniform ((Finset.univ : Finset (Fin 2)) ×ˢ (Finset.univ : Finset (Fin 8)))
        (2 * 8) (joint (deltaDist 0) (deltaDist 0))
      = Real.log 16 := by
  have h := klToUniform_joint_delta_eq_log_mul
    (Finset.univ : Finset (Fin 2)) (Finset.univ : Finset (Fin 8))
    2 8 (by norm_num) (by norm_num) 0 (by decide) 0 (by decide)
  convert h using 2
  norm_num

/-- **(m, n) = (8, 30).** Delta-corner additivity yields
    `log 8 + log 30 = log 240`. -/
theorem klToUniform_joint_delta_8_30 :
    klToUniform ((Finset.univ : Finset (Fin 8)) ×ˢ (Finset.univ : Finset (Fin 30)))
        (8 * 30) (joint (deltaDist 0) (deltaDist 0))
      = Real.log 240 := by
  have h := klToUniform_joint_delta_eq_log_mul
    (Finset.univ : Finset (Fin 8)) (Finset.univ : Finset (Fin 30))
    8 30 (by norm_num) (by norm_num) 0 (by decide) 0 (by decide)
  convert h using 2
  norm_num

/-! ### Headline -/

/-- **Headline (KL additivity on tensor products).** For two finite
    probability distributions `p, q` on `s, t` and moduli `m, n > 0`:

    * **General additivity.**
      `D_KL(P ⊗ Q ‖ U_m ⊗ U_n) = D_KL(P ‖ U_m) + D_KL(Q ‖ U_n)`.
    * **Delta corner.** `D_KL(δ ⊗ δ ‖ U_m ⊗ U_n) = log(m · n)`.
    * **Uniform corner.** `D_KL(U_m ⊗ U_n ‖ U_m ⊗ U_n) = 0`.
    * **PT numerics.** `(2, 8) ⇒ log 16`, `(8, 30) ⇒ log 240`. -/
theorem KL_tensor_additivity_summary :
    -- General KL additivity
    (∀ {α β : Type*} (s : Finset α) (t : Finset β)
        (m n : ℝ), 0 < m → 0 < n →
      ∀ (p : α → ℝ) (q : β → ℝ),
        (∀ a ∈ s, 0 ≤ p a) → (∀ b ∈ t, 0 ≤ q b) →
        (∑ a ∈ s, p a = 1) → (∑ b ∈ t, q b = 1) →
          klToUniform (s ×ˢ t) (m * n) (joint p q)
            = klToUniform s m p + klToUniform t n q)
    -- (2, 8) instance
    ∧ klToUniform ((Finset.univ : Finset (Fin 2)) ×ˢ (Finset.univ : Finset (Fin 8)))
        (2 * 8) (joint (deltaDist 0) (deltaDist 0))
        = Real.log 16
    -- (8, 30) instance
    ∧ klToUniform ((Finset.univ : Finset (Fin 8)) ×ˢ (Finset.univ : Finset (Fin 30)))
        (8 * 30) (joint (deltaDist 0) (deltaDist 0))
        = Real.log 240 :=
  ⟨@klToUniform_joint,
   klToUniform_joint_delta_2_8,
   klToUniform_joint_delta_8_30⟩

end PT.Information
