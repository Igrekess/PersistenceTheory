/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Information.GFTIdentity
import PT.Information.GFTSpecialisations
import PT.Information.EntropyBoundsTight
import PT.Information.MutualInformationBasic
import Mathlib.Tactic

/-!
# Conditional entropy `H(X | Y) := H(X, Y) − H(Y)`

For two finite probability distributions `p : α → ℝ` on `s : Finset α`
and `q : β → ℝ` on `t : Finset β`, the **conditional entropy** of `X`
given `Y` (relative to a joint distribution `jointPQ : α × β → ℝ`) is
defined by

  `condEntropy s t jointPQ q := shannonH (s ×ˢ t) jointPQ − shannonH t q`.

This is the standard definition `H(X | Y) = H(X, Y) − H(Y)`. The chain
rule `H(X, Y) = H(Y) + H(X | Y)` is then a trivial rearrangement.

When the joint distribution factorises as a product `jointPQ = joint p q`
with `∑ p = ∑ q = 1`, then `H(X, Y) = H(X) + H(Y)` (see
`shannonH_joint`), so

  `condEntropy s t (joint p q) q = shannonH s p`

i.e. *under independence, `H(X | Y) = H(X)`* — the marginal entropy of
`X` is unaffected by the (independent) `Y`.

The link with mutual information at independence is
`mutualInfo p q = shannonH p − condEntropy = 0` (both sides vanish).

## Corner cases

* **Delta corner.** `H(δ_a | δ_b) = 0`.
* **Uniform corner.** `H(U_m | U_n) = log m` (the conditional entropy of
  a uniform on `m` symbols, given any independent variable, is `log m`).

## Reference

Monograph Chapter 4 §"Entropie conditionnelle", elementary algebraic
corollary of `shannonH_joint`. Follow-up to `MutualInformationBasic`.
-/

namespace PT.Information

open Real Finset

/-! ### Definition -/

/-- **Conditional entropy** `H(X | Y) := H(X, Y) − H(Y)`, where `H(X, Y)`
    is the joint entropy on `s ×ˢ t` (with joint distribution `jointPQ`)
    and `H(Y)` is the marginal entropy of `Y` on `t`. -/
noncomputable def condEntropy
    {α β : Type*} (s : Finset α) (t : Finset β)
    (jointPQ : α × β → ℝ) (q : β → ℝ) : ℝ :=
  shannonH (s ×ˢ t) jointPQ - shannonH t q

/-! ### Chain rule -/

/-- **Chain rule (rearranged definition).**
    `H(X, Y) = H(Y) + H(X | Y)`. This is immediate from the definition
    of `condEntropy`. -/
theorem shannonH_joint_eq_add_condEntropy
    {α β : Type*} (s : Finset α) (t : Finset β)
    (jointPQ : α × β → ℝ) (q : β → ℝ) :
    shannonH (s ×ˢ t) jointPQ
      = shannonH t q + condEntropy s t jointPQ q := by
  unfold condEntropy
  ring

/-! ### Independence: `H(X | Y) = H(X)` for product joints -/

/-- **Independence ⇒ `H(X | Y) = H(X)`.**
    If `jointPQ = joint p q` is the product distribution and
    `∑ p = ∑ q = 1`, then the conditional entropy reduces to the
    marginal entropy of `X`. Direct corollary of `shannonH_joint`. -/
theorem condEntropy_indep
    {α β : Type*} (s : Finset α) (t : Finset β)
    (p : α → ℝ) (q : β → ℝ)
    (hp_sum : ∑ a ∈ s, p a = 1) (hq_sum : ∑ b ∈ t, q b = 1) :
    condEntropy s t (joint p q) q = shannonH s p := by
  unfold condEntropy
  rw [shannonH_joint s t p q hp_sum hq_sum]
  ring

/-! ### Link with mutual information -/

/-- **Mutual information ↔ conditional entropy (independence).**
    At independence (`jointPQ = joint p q`, `∑ p = ∑ q = 1`),
    `I(p ; q) = H(p) − H(X | Y)`. Both sides are zero. -/
theorem mutualInfo_eq_shannonH_sub_condEntropy_indep
    {α β : Type*} (s : Finset α) (t : Finset β)
    (p : α → ℝ) (q : β → ℝ)
    (hp_sum : ∑ a ∈ s, p a = 1) (hq_sum : ∑ b ∈ t, q b = 1) :
    mutualInfo s t p q
      = shannonH s p - condEntropy s t (joint p q) q := by
  rw [condEntropy_indep s t p q hp_sum hq_sum,
      mutualInfo_indep_eq_zero s t p q hp_sum hq_sum]
  ring

/-! ### Delta corner: `H(δ_a | δ_b) = 0` -/

/-- **Conditional entropy at the delta corner vanishes.**
    `H(δ_a | δ_b) = H(δ_{(a,b)}) − H(δ_b) = 0 − 0 = 0`. -/
theorem condEntropy_delta
    {α β : Type*} [DecidableEq α] [DecidableEq β]
    (s : Finset α) (t : Finset β)
    (a₀ : α) (ha : a₀ ∈ s) (b₀ : β) (hb : b₀ ∈ t) :
    condEntropy s t (joint (deltaDist a₀) (deltaDist b₀)) (deltaDist b₀) = 0 := by
  unfold condEntropy
  rw [shannonH_joint_delta s t a₀ ha b₀ hb, shannonH_delta t b₀ hb]
  ring

/-! ### Uniform corner: `H(U_m | U_n) = log m` -/

/-- **Conditional entropy at the uniform corner.**
    `H(U_m | U_n) = H(U_m ⊗ U_n) − H(U_n) = log(m·n) − log n = log m`,
    obtained algebraically via `shannonH_joint` + `shannonH_uniform_eq_log`. -/
theorem condEntropy_uniform
    {α β : Type*} (s : Finset α) (t : Finset β)
    (m n : ℝ) (hm : 0 < m) (hn : 0 < n)
    (hcard_s : (s.card : ℝ) = m) (hcard_t : (t.card : ℝ) = n) :
    condEntropy s t
        (joint (fun _ : α => (1 : ℝ) / m) (fun _ : β => (1 : ℝ) / n))
        (fun _ : β => (1 : ℝ) / n) = Real.log m := by
  -- Use the independence reduction: condEntropy = H(U_m) = log m.
  have hp_sum : ∑ _a ∈ s, (1 : ℝ) / m = 1 := by
    rw [Finset.sum_const, nsmul_eq_mul, hcard_s, mul_one_div]
    exact div_self (ne_of_gt hm)
  have hq_sum : ∑ _b ∈ t, (1 : ℝ) / n = 1 := by
    rw [Finset.sum_const, nsmul_eq_mul, hcard_t, mul_one_div]
    exact div_self (ne_of_gt hn)
  rw [condEntropy_indep s t _ _ hp_sum hq_sum]
  exact shannonH_uniform_eq_log s m hm hcard_s

/-! ### PT-canonical numerical instances -/

/-- **`H(U_2 | U_8) = log 2`** on `Fin 2 × Fin 8`. -/
theorem condEntropy_uniform_2_8 :
    condEntropy (Finset.univ : Finset (Fin 2)) (Finset.univ : Finset (Fin 8))
      (joint (fun _ => (1 : ℝ) / 2) (fun _ => (1 : ℝ) / 8))
      (fun _ => (1 : ℝ) / 8) = Real.log 2 := by
  apply condEntropy_uniform _ _ 2 8 (by norm_num) (by norm_num) (by simp) (by simp)

/-- **`H(δ_0 | δ_0) = 0`** on `Fin 8 × Fin 30`. -/
theorem condEntropy_delta_8_30 :
    condEntropy (Finset.univ : Finset (Fin 8)) (Finset.univ : Finset (Fin 30))
      (joint (deltaDist (0 : Fin 8)) (deltaDist (0 : Fin 30)))
      (deltaDist (0 : Fin 30)) = 0 := by
  apply condEntropy_delta _ _ (0 : Fin 8) (by decide) (0 : Fin 30) (by decide)

/-! ### Headline -/

/-- **Headline (conditional entropy).** For any two finite probability
    distributions `(p, q)` on `(s, t)` with `∑ p = ∑ q = 1`:

    * **Chain rule.** `H(X, Y) = H(Y) + H(X | Y)`.
    * **Independence.** `H(X | Y) = H(X)` for the product joint.
    * **Link with mutual information at independence.**
      `I(p ; q) = H(p) − H(X | Y)` (both sides are zero).
    * **Delta corner.** `H(δ_a | δ_b) = 0`.
    * **Uniform corner.** `H(U_m | U_n) = log m`. -/
theorem condEntropy_summary
    {α β : Type*} [DecidableEq α] [DecidableEq β]
    (s : Finset α) (t : Finset β)
    (p : α → ℝ) (q : β → ℝ)
    (hp_sum : ∑ a ∈ s, p a = 1) (hq_sum : ∑ b ∈ t, q b = 1)
    (m n : ℝ) (hm : 0 < m) (hn : 0 < n)
    (hcard_s : (s.card : ℝ) = m) (hcard_t : (t.card : ℝ) = n)
    (a₀ : α) (ha : a₀ ∈ s) (b₀ : β) (hb : b₀ ∈ t) :
    -- Chain rule
    shannonH (s ×ˢ t) (joint p q)
        = shannonH t q + condEntropy s t (joint p q) q
    -- Independence: H(X | Y) = H(X)
    ∧ condEntropy s t (joint p q) q = shannonH s p
    -- Mutual information ↔ conditional entropy at independence
    ∧ mutualInfo s t p q = shannonH s p - condEntropy s t (joint p q) q
    -- Delta corner: H(X | Y) = 0
    ∧ condEntropy s t (joint (deltaDist a₀) (deltaDist b₀)) (deltaDist b₀) = 0
    -- Uniform corner: H(X | Y) = log m
    ∧ condEntropy s t
        (joint (fun _ : α => (1 : ℝ) / m) (fun _ : β => (1 : ℝ) / n))
        (fun _ : β => (1 : ℝ) / n) = Real.log m :=
  ⟨shannonH_joint_eq_add_condEntropy s t (joint p q) q,
   condEntropy_indep s t p q hp_sum hq_sum,
   mutualInfo_eq_shannonH_sub_condEntropy_indep s t p q hp_sum hq_sum,
   condEntropy_delta s t a₀ ha b₀ hb,
   condEntropy_uniform s t m n hm hn hcard_s hcard_t⟩

end PT.Information
