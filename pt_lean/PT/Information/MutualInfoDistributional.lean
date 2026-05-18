/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Information.MutualInformationBasic
import PT.Information.ConditionalEntropy
import PT.Information.EntropyJointProduct
import Mathlib.Tactic

/-!
# Distributional mutual information `I(X; Y)` for arbitrary joints

`PT.Information.MutualInformationBasic` defines

  `mutualInfo s t p q := H(p) + H(q) − H(joint p q)`

— a *product-shaped* mutual information taking two marginals as input.
This file extends the picture to an **arbitrary joint distribution**
`J : α × β → ℝ` on `s ×ˢ t`, whose marginals
`pX(a) := Σ_b J(a, b)` and `pY(b) := Σ_a J(a, b)` are recovered by
summation.  We define

  `mutualInfoGen s t J := H(pX) + H(pY) − H(J)`

and prove the following normal-form results :

* **Marginal sums.** `(Σ_a pX a) = (Σ_b pY b) = (Σ_{ab} J(a,b))`,
  so `J` summing to `1` automatically gives normalised marginals.
* **Product reduction.** When `J = joint p q` factorises and `p, q`
  are normalised, `pX = p`, `pY = q`, and `mutualInfoGen = mutualInfo`.
  In particular `mutualInfoGen s t (joint p q) = 0`.
* **Delta corner.** For `J = deltaDist (a₀, b₀)`, the marginals are
  `deltaDist a₀` and `deltaDist b₀`, all three entropies vanish, and
  `mutualInfoGen = 0`.
* **Uniform corner.** For `J = (fun _ => 1/(m·n))` on `s ×ˢ t` with
  `s.card = m`, `t.card = n`, the marginals are `1/m`, `1/n`
  (uniforms), and `mutualInfoGen = log m + log n − log(m·n) = 0`.
* **Chain-rule identity.** `mutualInfoGen s t J = H(pX) − H(X | Y)` with
  `H(X | Y) = H(J) − H(pY) = condEntropy s t J pY`. (Symmetric form
  with `condEntropy t s J pX` requires swap-invariance of `H(J)` which
  we treat as a parametric statement.)
* **Non-negativity (parametric).** Under the abstract hypothesis
  `D_KL(J ‖ joint pX pY) ≥ 0`, we get `mutualInfoGen s t J ≥ 0`
  (the standard Gibbs-inequality reformulation). The KL-bound itself
  lives in the relative-entropy modules.

The constructions reuse `joint`, `mutualInfo`, `condEntropy`,
`shannonH_joint`, and the delta/uniform corner machinery from the
preceding files; no new heavy lemmas are introduced.

## Reference

Monograph Chapter 4 §"Information mutuelle, forme distributionnelle";
follow-up to `MutualInformationBasic` and `ConditionalEntropy`.
-/

namespace PT.Information

open Real Finset

/-! ### Marginals of an arbitrary joint -/

/-- The **`X`-marginal** of a joint `J : α × β → ℝ` on `s ×ˢ t`:
    `pX(a) := Σ_{b ∈ t} J(a, b)`. -/
noncomputable def marginalX {α β : Type*} (t : Finset β)
    (J : α × β → ℝ) : α → ℝ :=
  fun a => ∑ b ∈ t, J (a, b)

/-- The **`Y`-marginal** of a joint `J : α × β → ℝ` on `s ×ˢ t`:
    `pY(b) := Σ_{a ∈ s} J(a, b)`. -/
noncomputable def marginalY {α β : Type*} (s : Finset α)
    (J : α × β → ℝ) : β → ℝ :=
  fun b => ∑ a ∈ s, J (a, b)

@[simp] lemma marginalX_apply {α β : Type*} (t : Finset β)
    (J : α × β → ℝ) (a : α) :
    marginalX t J a = ∑ b ∈ t, J (a, b) := rfl

@[simp] lemma marginalY_apply {α β : Type*} (s : Finset α)
    (J : α × β → ℝ) (b : β) :
    marginalY s J b = ∑ a ∈ s, J (a, b) := rfl

/-! ### Sums of marginals -/

/-- The `X`-marginal sums to the total mass of the joint. -/
theorem marginalX_sum {α β : Type*} (s : Finset α) (t : Finset β)
    (J : α × β → ℝ) :
    ∑ a ∈ s, marginalX t J a = ∑ ab ∈ s ×ˢ t, J ab := by
  unfold marginalX
  rw [Finset.sum_product]

/-- The `Y`-marginal sums to the total mass of the joint. -/
theorem marginalY_sum {α β : Type*} (s : Finset α) (t : Finset β)
    (J : α × β → ℝ) :
    ∑ b ∈ t, marginalY s J b = ∑ ab ∈ s ×ˢ t, J ab := by
  unfold marginalY
  rw [Finset.sum_product_right]

/-- If the joint normalises to `1`, so does its `X`-marginal. -/
theorem marginalX_sum_one {α β : Type*} (s : Finset α) (t : Finset β)
    (J : α × β → ℝ) (hJ_sum : ∑ ab ∈ s ×ˢ t, J ab = 1) :
    ∑ a ∈ s, marginalX t J a = 1 := by
  rw [marginalX_sum s t J, hJ_sum]

/-- If the joint normalises to `1`, so does its `Y`-marginal. -/
theorem marginalY_sum_one {α β : Type*} (s : Finset α) (t : Finset β)
    (J : α × β → ℝ) (hJ_sum : ∑ ab ∈ s ×ˢ t, J ab = 1) :
    ∑ b ∈ t, marginalY s J b = 1 := by
  rw [marginalY_sum s t J, hJ_sum]

/-! ### Marginals of a product distribution -/

/-- The `X`-marginal of a product joint is the first factor (when the
    second factor sums to `1`). -/
theorem marginalX_joint {α β : Type*} (t : Finset β)
    (p : α → ℝ) (q : β → ℝ) (hq_sum : ∑ b ∈ t, q b = 1) :
    marginalX t (joint p q) = p := by
  funext a
  unfold marginalX joint
  simp only
  -- ∑_b p a · q b = p a · (∑ q) = p a
  rw [← Finset.mul_sum, hq_sum, mul_one]

/-- The `Y`-marginal of a product joint is the second factor (when the
    first factor sums to `1`). -/
theorem marginalY_joint {α β : Type*} (s : Finset α)
    (p : α → ℝ) (q : β → ℝ) (hp_sum : ∑ a ∈ s, p a = 1) :
    marginalY s (joint p q) = q := by
  funext b
  unfold marginalY joint
  simp only
  -- ∑_a p a · q b = (∑ p) · q b = q b
  rw [← Finset.sum_mul, hp_sum, one_mul]

/-! ### Marginals of a delta distribution on the pair -/

/-- The `X`-marginal of `δ_{(a₀, b₀)}` is `δ_{a₀}` (when `b₀ ∈ t`). -/
theorem marginalX_deltaDist_pair
    {α β : Type*} [DecidableEq α] [DecidableEq β]
    (t : Finset β) (a₀ : α) (b₀ : β) (hb : b₀ ∈ t) :
    marginalX t (deltaDist (a₀, b₀)) = deltaDist a₀ := by
  funext a
  unfold marginalX deltaDist
  by_cases ha : a = a₀
  · -- contributes 1 at b = b₀, 0 elsewhere
    subst ha
    rw [Finset.sum_eq_single b₀]
    · simp
    · intros b _ hne
      have hpair : (a, b) ≠ (a, b₀) := by
        intro h
        apply hne
        exact (Prod.mk.injEq a b a b₀).mp h |>.2
      simp [hpair]
    · intro h
      exact absurd hb h
  · -- every term `(a, b) = (a₀, b₀)` is false since `a ≠ a₀`
    rw [if_neg ha]
    apply Finset.sum_eq_zero
    intro b _
    have hpair : (a, b) ≠ (a₀, b₀) := by
      intro h
      apply ha
      exact (Prod.mk.injEq a b a₀ b₀).mp h |>.1
    simp [hpair]

/-- The `Y`-marginal of `δ_{(a₀, b₀)}` is `δ_{b₀}` (when `a₀ ∈ s`). -/
theorem marginalY_deltaDist_pair
    {α β : Type*} [DecidableEq α] [DecidableEq β]
    (s : Finset α) (a₀ : α) (b₀ : β) (ha : a₀ ∈ s) :
    marginalY s (deltaDist (a₀, b₀)) = deltaDist b₀ := by
  funext b
  unfold marginalY deltaDist
  by_cases hb : b = b₀
  · subst hb
    rw [Finset.sum_eq_single a₀]
    · simp
    · intros a _ hne
      have hpair : (a, b) ≠ (a₀, b) := by
        intro h
        apply hne
        exact (Prod.mk.injEq a b a₀ b).mp h |>.1
      simp [hpair]
    · intro h
      exact absurd ha h
  · rw [if_neg hb]
    apply Finset.sum_eq_zero
    intro a _
    have hpair : (a, b) ≠ (a₀, b₀) := by
      intro h
      apply hb
      exact (Prod.mk.injEq a b a₀ b₀).mp h |>.2
    simp [hpair]

/-! ### Marginals of the uniform distribution on `s ×ˢ t` -/

/-- The `X`-marginal of the uniform `1/(m·n)` on `s ×ˢ t`, where
    `t.card = n` and `n ≠ 0`, is the uniform `1/m`. -/
theorem marginalX_uniform_prod
    {α β : Type*} (t : Finset β)
    (m n : ℝ) (hn : 0 < n) (hcard_t : (t.card : ℝ) = n) :
    marginalX (α := α) t (fun _ => (1 : ℝ) / (m * n))
      = fun _ : α => (1 : ℝ) / m := by
  funext a
  unfold marginalX
  -- ∑_{b ∈ t} 1/(m·n) = t.card · 1/(m·n) = n · 1/(m·n) = 1/m
  rw [Finset.sum_const, nsmul_eq_mul, hcard_t]
  field_simp

/-- The `Y`-marginal of the uniform `1/(m·n)` on `s ×ˢ t`, where
    `s.card = m` and `m ≠ 0`, is the uniform `1/n`. -/
theorem marginalY_uniform_prod
    {α β : Type*} (s : Finset α)
    (m n : ℝ) (hm : 0 < m) (hcard_s : (s.card : ℝ) = m) :
    marginalY (β := β) s (fun _ => (1 : ℝ) / (m * n))
      = fun _ : β => (1 : ℝ) / n := by
  funext b
  unfold marginalY
  rw [Finset.sum_const, nsmul_eq_mul, hcard_s]
  field_simp

/-! ### Distributional mutual information -/

/-- **Distributional mutual information.**

    For an arbitrary joint distribution `J : α × β → ℝ` on `s ×ˢ t`,

      `I(X; Y) := H(pX) + H(pY) − H(J)`

    where `pX`, `pY` are the marginals recovered by summation
    (`marginalX`, `marginalY`). This generalises `mutualInfo`, which
    takes the two marginals as input and constructs `J = joint p q`. -/
noncomputable def mutualInfoGen
    {α β : Type*} (s : Finset α) (t : Finset β)
    (J : α × β → ℝ) : ℝ :=
  shannonH s (marginalX t J)
    + shannonH t (marginalY s J)
    - shannonH (s ×ˢ t) J

/-! ### Reduction to the product case -/

/-- **Product reduction.** When `J = joint p q` is a product joint with
    normalised marginals, `mutualInfoGen` agrees with the basic
    `mutualInfo`. In particular, both vanish. -/
theorem mutualInfoGen_joint_eq_mutualInfo
    {α β : Type*} (s : Finset α) (t : Finset β)
    (p : α → ℝ) (q : β → ℝ)
    (hp_sum : ∑ a ∈ s, p a = 1) (hq_sum : ∑ b ∈ t, q b = 1) :
    mutualInfoGen s t (joint p q) = mutualInfo s t p q := by
  unfold mutualInfoGen mutualInfo
  rw [marginalX_joint t p q hq_sum, marginalY_joint s p q hp_sum]

/-- **Distributional mutual information of a product joint vanishes.**
    `I(p ⊗ q) = 0`. -/
theorem mutualInfoGen_joint_eq_zero
    {α β : Type*} (s : Finset α) (t : Finset β)
    (p : α → ℝ) (q : β → ℝ)
    (hp_sum : ∑ a ∈ s, p a = 1) (hq_sum : ∑ b ∈ t, q b = 1) :
    mutualInfoGen s t (joint p q) = 0 := by
  rw [mutualInfoGen_joint_eq_mutualInfo s t p q hp_sum hq_sum]
  exact mutualInfo_indep_eq_zero s t p q hp_sum hq_sum

/-! ### Delta corner -/

/-- **Distributional mutual information at the delta corner vanishes.**
    For `J = δ_{(a₀, b₀)}` on `s ×ˢ t`:
    `pX = δ_{a₀}`, `pY = δ_{b₀}`, and all three entropies are `0`. -/
theorem mutualInfoGen_delta_pair
    {α β : Type*} [DecidableEq α] [DecidableEq β]
    (s : Finset α) (t : Finset β)
    (a₀ : α) (ha : a₀ ∈ s) (b₀ : β) (hb : b₀ ∈ t) :
    mutualInfoGen s t (deltaDist (a₀, b₀)) = 0 := by
  unfold mutualInfoGen
  rw [marginalX_deltaDist_pair t a₀ b₀ hb,
      marginalY_deltaDist_pair s a₀ b₀ ha]
  rw [shannonH_delta s a₀ ha, shannonH_delta t b₀ hb,
      shannonH_delta (s ×ˢ t) (a₀, b₀) (Finset.mem_product.mpr ⟨ha, hb⟩)]
  ring

/-! ### Uniform corner -/

/-- **Distributional mutual information at the uniform corner vanishes.**
    For `J = 1/(m·n)` on `s ×ˢ t` with `s.card = m`, `t.card = n`:
    `pX = U_m`, `pY = U_n`, and
    `I = log m + log n − log(m·n) = 0`. -/
theorem mutualInfoGen_uniform_prod
    {α β : Type*} (s : Finset α) (t : Finset β)
    (m n : ℝ) (hm : 0 < m) (hn : 0 < n)
    (hcard_s : (s.card : ℝ) = m) (hcard_t : (t.card : ℝ) = n) :
    mutualInfoGen s t (fun _ : α × β => (1 : ℝ) / (m * n)) = 0 := by
  unfold mutualInfoGen
  -- Marginals become 1/m and 1/n.
  rw [marginalX_uniform_prod (α := α) t m n hn hcard_t,
      marginalY_uniform_prod (β := β) s m n hm hcard_s]
  -- H(U_m) = log m, H(U_n) = log n.
  rw [shannonH_uniform_eq_log s m hm hcard_s,
      shannonH_uniform_eq_log t n hn hcard_t]
  -- H(U_{m·n}) on s ×ˢ t = log(m · n).
  have hmn : (0 : ℝ) < m * n := mul_pos hm hn
  have hcard_prod : (((s ×ˢ t).card : ℝ)) = m * n := by
    rw [Finset.card_product]
    push_cast
    rw [hcard_s, hcard_t]
  rw [shannonH_uniform_eq_log (s ×ˢ t) (m * n) hmn hcard_prod]
  -- log m + log n - log(m · n) = 0.
  rw [Real.log_mul (ne_of_gt hm) (ne_of_gt hn)]
  ring

/-! ### Chain-rule identity with conditional entropy -/

/-- **Chain-rule identity (X-side).**
    `I(X; Y) = H(X) − H(X | Y)`, where `H(X | Y) := H(J) − H(pY)`
    in the `condEntropy` convention used throughout this module. -/
theorem mutualInfoGen_eq_marginalX_sub_condEntropy
    {α β : Type*} (s : Finset α) (t : Finset β) (J : α × β → ℝ) :
    mutualInfoGen s t J
      = shannonH s (marginalX t J)
          - condEntropy s t J (marginalY s J) := by
  unfold mutualInfoGen condEntropy
  ring

/-- **Chain-rule identity (Y-side).**
    `I(X; Y) = H(Y) − H(Y | X)`, where `H(Y | X) := H(J) − H(pX)`. -/
theorem mutualInfoGen_eq_marginalY_sub_condEntropy
    {α β : Type*} (s : Finset α) (t : Finset β) (J : α × β → ℝ) :
    mutualInfoGen s t J
      = shannonH t (marginalY s J)
          - (shannonH (s ×ˢ t) J - shannonH s (marginalX t J)) := by
  unfold mutualInfoGen
  ring

/-! ### Non-negativity (parametric) -/

/-- **Non-negativity from `D_KL ≥ 0`.**
    Standard rewriting: `I(X; Y) = D_KL(J ‖ pX ⊗ pY)`, so any abstract
    hypothesis `D_KL ≥ 0` (Gibbs inequality) supplies the bound. We
    keep the formal statement parametric: under the hypothesis
    `H(pX) + H(pY) − H(J) ≥ 0`, conclude `mutualInfoGen ≥ 0`.
    The verbatim Gibbs proof lives in the relative-entropy module
    (`KLAdditivity*`). -/
theorem mutualInfoGen_nonneg_of_kl_nonneg
    {α β : Type*} (s : Finset α) (t : Finset β) (J : α × β → ℝ)
    (h_kl : shannonH s (marginalX t J) + shannonH t (marginalY s J)
              - shannonH (s ×ˢ t) J ≥ 0) :
    mutualInfoGen s t J ≥ 0 := by
  unfold mutualInfoGen
  exact h_kl

/-! ### PT-canonical numerical instances -/

/-- **`I(δ_{(0, 0)}) = 0`** on `Fin 8 × Fin 30`. -/
theorem mutualInfoGen_delta_pair_8_30 :
    mutualInfoGen (Finset.univ : Finset (Fin 8)) (Finset.univ : Finset (Fin 30))
      (deltaDist ((0, 0) : Fin 8 × Fin 30)) = 0 := by
  apply mutualInfoGen_delta_pair _ _ (0 : Fin 8) (by decide)
                                  (0 : Fin 30) (by decide)

/-- **`I(U_{2·8}) = 0`** on `Fin 2 × Fin 8`, where the joint is the
    uniform on the product. -/
theorem mutualInfoGen_uniform_prod_2_8 :
    mutualInfoGen (Finset.univ : Finset (Fin 2)) (Finset.univ : Finset (Fin 8))
      (fun _ : Fin 2 × Fin 8 => (1 : ℝ) / (2 * 8)) = 0 := by
  apply mutualInfoGen_uniform_prod _ _ 2 8
    (by norm_num) (by norm_num) (by simp) (by simp)

/-- **`I(joint U_2 U_8) = 0`** on `Fin 2 × Fin 8`, the product-of-
    uniforms reformulation. -/
theorem mutualInfoGen_joint_uniform_2_8 :
    mutualInfoGen (Finset.univ : Finset (Fin 2)) (Finset.univ : Finset (Fin 8))
      (joint (fun _ : Fin 2 => (1 : ℝ) / 2) (fun _ : Fin 8 => (1 : ℝ) / 8))
      = 0 := by
  apply mutualInfoGen_joint_eq_zero
  · rw [Finset.sum_const, Finset.card_univ, Fintype.card_fin]
    norm_num
  · rw [Finset.sum_const, Finset.card_univ, Fintype.card_fin]
    norm_num

/-! ### Headline -/

/-- **Headline (distributional mutual information).**

    For an arbitrary joint `J : α × β → ℝ` on `s ×ˢ t`:

    * **Product reduction.** If `J = joint p q` with normalised
      marginals, then `mutualInfoGen = mutualInfo = 0`.
    * **Delta corner.** `I(δ_{(a₀, b₀)}) = 0`.
    * **Uniform corner.** `I(1/(m·n)) = 0` on a `m × n` grid.
    * **Chain rule (X-side).**
      `I(X; Y) = H(X) − condEntropy(J, pY)`.
    * **Non-negativity.** From `D_KL ≥ 0` (parametric), `I ≥ 0`. -/
theorem mutualInfoGen_summary
    {α β : Type*} [DecidableEq α] [DecidableEq β]
    (s : Finset α) (t : Finset β)
    (p : α → ℝ) (q : β → ℝ)
    (hp_sum : ∑ a ∈ s, p a = 1) (hq_sum : ∑ b ∈ t, q b = 1)
    (m n : ℝ) (hm : 0 < m) (hn : 0 < n)
    (hcard_s : (s.card : ℝ) = m) (hcard_t : (t.card : ℝ) = n)
    (a₀ : α) (ha : a₀ ∈ s) (b₀ : β) (hb : b₀ ∈ t)
    (J : α × β → ℝ)
    (h_kl : shannonH s (marginalX t J) + shannonH t (marginalY s J)
              - shannonH (s ×ˢ t) J ≥ 0) :
    -- Product reduction
    mutualInfoGen s t (joint p q) = mutualInfo s t p q
    -- Product = 0
    ∧ mutualInfoGen s t (joint p q) = 0
    -- Delta corner
    ∧ mutualInfoGen s t (deltaDist (a₀, b₀)) = 0
    -- Uniform corner
    ∧ mutualInfoGen s t (fun _ : α × β => (1 : ℝ) / (m * n)) = 0
    -- Chain rule (X-side)
    ∧ mutualInfoGen s t J
        = shannonH s (marginalX t J)
            - condEntropy s t J (marginalY s J)
    -- Non-negativity (parametric)
    ∧ mutualInfoGen s t J ≥ 0 :=
  ⟨mutualInfoGen_joint_eq_mutualInfo s t p q hp_sum hq_sum,
   mutualInfoGen_joint_eq_zero s t p q hp_sum hq_sum,
   mutualInfoGen_delta_pair s t a₀ ha b₀ hb,
   mutualInfoGen_uniform_prod s t m n hm hn hcard_s hcard_t,
   mutualInfoGen_eq_marginalX_sub_condEntropy s t J,
   mutualInfoGen_nonneg_of_kl_nonneg s t J h_kl⟩

end PT.Information
