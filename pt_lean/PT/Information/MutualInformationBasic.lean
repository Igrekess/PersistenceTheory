/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Information.GFTIdentity
import PT.Information.GFTSpecialisations
import PT.Information.BekensteinExtensions
import PT.Information.EntropyBoundsTight
import PT.Information.EntropyAdditivityCorners
import Mathlib.Algebra.BigOperators.Group.Finset.Sigma
import Mathlib.Tactic

/-!
# Basic Mutual Information `I(X; Y) := H(X) + H(Y) - H(X, Y)`

For two independent finite probability distributions
`p : α → ℝ` on `s : Finset α` and `q : β → ℝ` on `t : Finset β`, the
product distribution

  `joint p q (a, b) := p a · q b`

normalises on `s ×ˢ t` whenever `∑ p = ∑ q = 1`. The Shannon entropy of
the product equals the sum of the marginal entropies:

  `H(p ⊗ q) = H(p) + H(q)`

(algebraic kernel : `negMulLog (x · y) = y · negMulLog x + x · negMulLog y`,
already in Mathlib). The mutual information is then

  `I(p, q) := H(p) + H(q) − H(joint p q)`

and vanishes on independent product distributions:

  `I(p, q) = 0`  whenever the joint factorises as `p a · q b`.

In particular:

* **Delta corner.** `I(δ_a, δ_b) = 0` (all three entropies vanish).
* **Uniform corner.** `I(U_m, U_n) = 0` since `U_m ⊗ U_n = U_{m n}` and
  `H(U_{m n}) = log(m n) = log m + log n`.

## Reference

Monograph Chapter 4 §"Information mutuelle", elementary algebraic
corollary of the GFT identity at independence. M4 article.
-/

namespace PT.Information

open Real Finset

/-! ### Joint product distribution -/

/-- The **product / joint** distribution `p ⊗ q` on `α × β`:
    `(p ⊗ q)(a, b) = p a · q b`. -/
noncomputable def joint {α β : Type*} (p : α → ℝ) (q : β → ℝ) :
    α × β → ℝ := fun ab => p ab.1 * q ab.2

@[simp] lemma joint_apply {α β : Type*} (p : α → ℝ) (q : β → ℝ) (a : α) (b : β) :
    joint p q (a, b) = p a * q b := rfl

/-- Non-negativity is preserved by the product distribution. -/
theorem joint_nonneg {α β : Type*}
    (s : Finset α) (t : Finset β) (p : α → ℝ) (q : β → ℝ)
    (hp : ∀ a ∈ s, 0 ≤ p a) (hq : ∀ b ∈ t, 0 ≤ q b) :
    ∀ ab ∈ s ×ˢ t, 0 ≤ joint p q ab := by
  intro ab hab
  rcases ab with ⟨a, b⟩
  simp only [Finset.mem_product] at hab
  exact mul_nonneg (hp a hab.1) (hq b hab.2)

/-- The joint distribution normalises:
    `∑_{(a, b)} p a · q b = (∑ p) · (∑ q) = 1`. -/
theorem joint_sum {α β : Type*}
    (s : Finset α) (t : Finset β) (p : α → ℝ) (q : β → ℝ)
    (hp_sum : ∑ a ∈ s, p a = 1) (hq_sum : ∑ b ∈ t, q b = 1) :
    ∑ ab ∈ s ×ˢ t, joint p q ab = 1 := by
  unfold joint
  -- ∑_{(a,b) ∈ s ×ˢ t} p a · q b = ∑_a ∑_b p a · q b
  rw [Finset.sum_product]
  -- = ∑_a p a · ∑_b q b = (∑ p)·(∑ q)
  simp_rw [← Finset.mul_sum]
  rw [← Finset.sum_mul, hp_sum, hq_sum, one_mul]

/-! ### Shannon entropy of a product distribution -/

/-- **Entropy of a product factorises additively.**
    `H(p ⊗ q) = H(p) + H(q)`, assuming `∑ p = ∑ q = 1`.

    Algebraic kernel: `negMulLog (x · y) = y · negMulLog x + x · negMulLog y`
    (Mathlib). Summing over `(a, b) ∈ s ×ˢ t`:

    `∑ negMulLog(p a · q b)
      = ∑_{a,b} q b · negMulLog(p a) + ∑_{a,b} p a · negMulLog(q b)
      = (∑ q) · ∑_a negMulLog(p a) + (∑ p) · ∑_b negMulLog(q b)
      = H(p) + H(q)`. -/
theorem shannonH_joint
    {α β : Type*} (s : Finset α) (t : Finset β)
    (p : α → ℝ) (q : β → ℝ)
    (hp_sum : ∑ a ∈ s, p a = 1) (hq_sum : ∑ b ∈ t, q b = 1) :
    shannonH (s ×ˢ t) (joint p q) = shannonH s p + shannonH t q := by
  unfold shannonH joint
  -- ∑_{ab ∈ s ×ˢ t} negMulLog(p ab.1 · q ab.2)
  --   = ∑_a ∑_b negMulLog(p a · q b)
  rw [Finset.sum_product]
  -- expand the inner via negMulLog_mul
  have h_inner : ∀ a ∈ s,
      ∑ b ∈ t, Real.negMulLog (p a * q b)
        = (∑ b ∈ t, q b) * Real.negMulLog (p a)
          + p a * (∑ b ∈ t, Real.negMulLog (q b)) := by
    intro a _
    simp_rw [Real.negMulLog_mul]
    rw [Finset.sum_add_distrib]
    congr 1
    · rw [← Finset.sum_mul]
    · rw [← Finset.mul_sum]
  rw [Finset.sum_congr rfl h_inner]
  -- ∑_a [ (∑ q) · negMulLog (p a) + p a · ∑_b negMulLog (q b) ]
  --   = (∑ q) · ∑_a negMulLog(p a) + (∑_a p a) · ∑_b negMulLog(q b)
  rw [Finset.sum_add_distrib]
  rw [hq_sum]
  simp_rw [one_mul]
  rw [← Finset.sum_mul, hp_sum, one_mul]

/-! ### Mutual information -/

/-- **Mutual information** of two finite distributions (basic / set-theoretic
    form): `I(X; Y) := H(X) + H(Y) − H(X, Y)`, where `H(X, Y)` is the joint
    entropy on `α × β` and `X, Y` are taken independent (joint = product). -/
noncomputable def mutualInfo
    {α β : Type*} (s : Finset α) (t : Finset β)
    (p : α → ℝ) (q : β → ℝ) : ℝ :=
  shannonH s p + shannonH t q - shannonH (s ×ˢ t) (joint p q)

/-- **Independence ⇒ zero mutual information.**
    If `∑ p = ∑ q = 1`, then `I(p, q) = 0` for the product joint. -/
theorem mutualInfo_indep_eq_zero
    {α β : Type*} (s : Finset α) (t : Finset β)
    (p : α → ℝ) (q : β → ℝ)
    (hp_sum : ∑ a ∈ s, p a = 1) (hq_sum : ∑ b ∈ t, q b = 1) :
    mutualInfo s t p q = 0 := by
  unfold mutualInfo
  rw [shannonH_joint s t p q hp_sum hq_sum]
  ring

/-! ### Delta corner: `I(δ_a, δ_b) = 0` -/

/-- The product of two delta distributions is itself a delta on the pair:
    `δ_a ⊗ δ_b = δ_{(a, b)}` (pointwise on `s ×ˢ t`). -/
theorem joint_deltaDist_eq_deltaDist
    {α β : Type*} [DecidableEq α] [DecidableEq β]
    (a₀ : α) (b₀ : β) :
    joint (deltaDist a₀) (deltaDist b₀) = deltaDist (a₀, b₀) := by
  funext ab
  rcases ab with ⟨a, b⟩
  unfold joint deltaDist
  by_cases ha : a = a₀
  · by_cases hb : b = b₀
    · simp [ha, hb]
    · simp [ha, hb]
  · by_cases hb : b = b₀
    · simp [ha, hb]
    · simp [ha, hb]

/-- **Joint entropy of two deltas vanishes.** Direct corollary of
    `shannonH_delta` applied to the pair `δ_{(a₀, b₀)}`. -/
theorem shannonH_joint_delta
    {α β : Type*} [DecidableEq α] [DecidableEq β]
    (s : Finset α) (t : Finset β)
    (a₀ : α) (ha : a₀ ∈ s) (b₀ : β) (hb : b₀ ∈ t) :
    shannonH (s ×ˢ t) (joint (deltaDist a₀) (deltaDist b₀)) = 0 := by
  rw [joint_deltaDist_eq_deltaDist]
  apply shannonH_delta
  rw [Finset.mem_product]
  exact ⟨ha, hb⟩

/-- **Mutual information at the delta corner vanishes.** `I(δ_a, δ_b) = 0`
    since the three constituent entropies are all `0`. -/
theorem mutualInfo_delta
    {α β : Type*} [DecidableEq α] [DecidableEq β]
    (s : Finset α) (t : Finset β)
    (a₀ : α) (ha : a₀ ∈ s) (b₀ : β) (hb : b₀ ∈ t) :
    mutualInfo s t (deltaDist a₀) (deltaDist b₀) = 0 := by
  unfold mutualInfo
  rw [shannonH_delta s a₀ ha, shannonH_delta t b₀ hb,
      shannonH_joint_delta s t a₀ ha b₀ hb]
  ring

/-! ### Uniform corner: `I(U_m, U_n) = 0` -/

/-- The product of two uniform distributions: pointwise equality
    `U_m(a) · U_n(b) = 1/(m · n) = U_{m n}(a, b)`. -/
theorem joint_uniform_eq_uniform_prod
    {α β : Type*} (m n : ℝ) (hm : m ≠ 0) (hn : n ≠ 0) :
    joint (α := α) (β := β) (fun _ => (1 : ℝ) / m) (fun _ => (1 : ℝ) / n)
      = fun _ => (1 : ℝ) / (m * n) := by
  funext ab
  unfold joint
  field_simp

/-- **Mutual information at the uniform corner vanishes.**
    `I(U_m, U_n) = H(U_m) + H(U_n) − H(U_{m n})
                 = log m + log n − log(m n) = 0`. -/
theorem mutualInfo_uniform
    {α β : Type*} (s : Finset α) (t : Finset β)
    (m n : ℝ) (hm : 0 < m) (hn : 0 < n)
    (hcard_s : (s.card : ℝ) = m) (hcard_t : (t.card : ℝ) = n) :
    mutualInfo s t (fun _ : α => (1 : ℝ) / m) (fun _ : β => (1 : ℝ) / n) = 0 := by
  -- Apply the independence formula with `∑ U_m = 1`, `∑ U_n = 1`.
  apply mutualInfo_indep_eq_zero
  · -- ∑_{a ∈ s} 1/m = s.card / m = m / m = 1
    rw [Finset.sum_const, nsmul_eq_mul, hcard_s, mul_one_div]
    exact div_self (ne_of_gt hm)
  · rw [Finset.sum_const, nsmul_eq_mul, hcard_t, mul_one_div]
    exact div_self (ne_of_gt hn)

/-! ### PT-canonical numerical instances -/

/-- **`I(U_2, U_8) = 0`.** Product of uniforms on `Fin 2 × Fin 8`. -/
theorem mutualInfo_uniform_2_8 :
    mutualInfo (Finset.univ : Finset (Fin 2)) (Finset.univ : Finset (Fin 8))
      (fun _ => (1 : ℝ) / 2) (fun _ => (1 : ℝ) / 8) = 0 := by
  apply mutualInfo_uniform _ _ 2 8 (by norm_num) (by norm_num) (by simp) (by simp)

/-- **`I(δ_0, δ_0) = 0` on `Fin 8 × Fin 30`.** -/
theorem mutualInfo_delta_8_30 :
    mutualInfo (Finset.univ : Finset (Fin 8)) (Finset.univ : Finset (Fin 30))
      (deltaDist 0) (deltaDist 0) = 0 := by
  apply mutualInfo_delta _ _ (0 : Fin 8) (by decide) (0 : Fin 30) (by decide)

/-! ### Headline -/

/-- **Headline (basic mutual information).** For any two finite
    probability distributions `(p, q)` on `(s, t)` with `∑ p = ∑ q = 1`:

    * **Normalisation.** `∑ joint p q = 1` on `s ×ˢ t`.
    * **Additivity of entropy under independence.**
      `H(joint p q) = H(p) + H(q)`.
    * **Mutual information vanishes at independence.**
      `I(p, q) := H(p) + H(q) − H(joint p q) = 0`.
    * **Delta corner.** `I(δ_a, δ_b) = 0`.
    * **Uniform corner.** `I(U_m, U_n) = 0` (and `U_m ⊗ U_n = U_{m n}`). -/
theorem mutualInfo_basic_summary
    {α β : Type*} [DecidableEq α] [DecidableEq β]
    (s : Finset α) (t : Finset β)
    (p : α → ℝ) (q : β → ℝ)
    (hp_sum : ∑ a ∈ s, p a = 1) (hq_sum : ∑ b ∈ t, q b = 1)
    (m n : ℝ) (hm : 0 < m) (hn : 0 < n)
    (hcard_s : (s.card : ℝ) = m) (hcard_t : (t.card : ℝ) = n)
    (a₀ : α) (ha : a₀ ∈ s) (b₀ : β) (hb : b₀ ∈ t) :
    -- Normalisation
    (∑ ab ∈ s ×ˢ t, joint p q ab) = 1
    -- Additivity of entropy under independence
    ∧ shannonH (s ×ˢ t) (joint p q) = shannonH s p + shannonH t q
    -- Mutual information vanishes on the product joint
    ∧ mutualInfo s t p q = 0
    -- Delta corner: I = 0
    ∧ mutualInfo s t (deltaDist a₀) (deltaDist b₀) = 0
    -- Uniform corner: I = 0
    ∧ mutualInfo s t (fun _ : α => (1 : ℝ) / m) (fun _ : β => (1 : ℝ) / n) = 0 :=
  ⟨joint_sum s t p q hp_sum hq_sum,
   shannonH_joint s t p q hp_sum hq_sum,
   mutualInfo_indep_eq_zero s t p q hp_sum hq_sum,
   mutualInfo_delta s t a₀ ha b₀ hb,
   mutualInfo_uniform s t m n hm hn hcard_s hcard_t⟩

end PT.Information
