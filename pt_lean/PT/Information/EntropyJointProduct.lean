/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Information.GFTIdentity
import PT.Information.GFTSpecialisations
import PT.Information.GFTSpecialMValues
import PT.Information.BekensteinExtensions
import PT.Information.EntropyBoundsTight
import PT.Information.EntropyMonotonicity
import PT.Information.MutualInformationBasic
import Mathlib.Analysis.SpecialFunctions.Log.NegMulLog
import Mathlib.Tactic

/-!
# Joint entropy and the `log m` increment from independent uniforms

This file consolidates the **algebraic identity**

  `H(p ⊗ U_m) − H(p) = log m`

— adjoining an independent uniform on `m` outcomes shifts the Shannon
entropy of any distribution by exactly `log m` — and derives a small
chain of corollaries: double / triple uniform cascades, delta-uniform
products, the three-factor chain `H(p ⊗ U_m ⊗ U_n)`, PT-canonical
numerical instances (`m, n ∈ {2, 8, 30}`), and a parametric strict
Bekenstein bound for the product distribution.

The core identity `H(p ⊗ U_m) = H(p) + log m` is already in
`EntropyMonotonicity.shannonH_joint_uniform_add_log`; we re-expose
it under a name that makes the increment manifest, and extend it.

## Headline results

* `shannonH_joint_uniform_increment` — the headline identity
  `H(p ⊗ U_m) − H(p) = log m` (subtractive form).
* `shannonH_joint_uniform_uniform` — `H(U_m ⊗ U_n) = log m + log n`.
* `shannonH_joint_uniform_uniform_log_mul` — `H(U_m ⊗ U_n) = log (m·n)`.
* `shannonH_joint_delta_uniform` — `H(δ ⊗ U_m) = log m`.
* `shannonH_joint_delta_delta_uniform` — `H(δ ⊗ δ ⊗ U_m) = log m`.
* `shannonH_joint_chain_three` — chain identity
  `H((p ⊗ U_m) ⊗ U_n) = H(p) + log m + log n`.
* `shannonH_joint_chain_three_right` — right-associated chain
  `H(p ⊗ (U_m ⊗ U_n)) = H(p) + log m + log n`.
* `shannonH_joint_uniform_strict_lt` — parametric strict bound
  `H(p ⊗ U_m) < log(s.card · m)`, assuming `H(p) < log s.card`.
* `shannonH_joint_increment_summary` — single headline bundle.

## Reference

Monograph Chapter 4 §"Additivité de l'entropie sous indépendance" and
follow-up to `EntropyMonotonicity.shannonH_joint_uniform_add_log`.
-/

namespace PT.Information

open Real Finset

/-! ### 1. Headline identity: `H(p ⊗ U_m) − H(p) = log m` -/

/-- **Headline (uniform increment identity).**
    Adjoining an independent uniform on `m` outcomes increases the
    Shannon entropy by exactly `log m`:

    `H(p ⊗ U_m) − H(p) = log m`.

    This is the subtractive form of
    `EntropyMonotonicity.shannonH_joint_uniform_add_log`. -/
theorem shannonH_joint_uniform_increment
    {α β : Type*} (s : Finset α) (t : Finset β)
    (p : α → ℝ) (hp_sum : ∑ a ∈ s, p a = 1)
    (m : ℝ) (hm : 0 < m) (hcard_t : (t.card : ℝ) = m) :
    shannonH (s ×ˢ t) (joint p (fun _ : β => (1 : ℝ) / m)) - shannonH s p
      = Real.log m := by
  have h := shannonH_joint_uniform_add_log s t p hp_sum m hm hcard_t
  linarith

/-- **Re-exposition (additive form).** `H(p ⊗ U_m) = H(p) + log m`. -/
theorem shannonH_joint_uniform_add
    {α β : Type*} (s : Finset α) (t : Finset β)
    (p : α → ℝ) (hp_sum : ∑ a ∈ s, p a = 1)
    (m : ℝ) (hm : 0 < m) (hcard_t : (t.card : ℝ) = m) :
    shannonH (s ×ˢ t) (joint p (fun _ : β => (1 : ℝ) / m))
      = shannonH s p + Real.log m :=
  shannonH_joint_uniform_add_log s t p hp_sum m hm hcard_t

/-! ### 2. Double uniform cascade: `H(U_m ⊗ U_n) = log m + log n` -/

/-- **Double uniform product.** `H(U_m ⊗ U_n) = log m + log n`.
    Direct corollary of `shannonH_joint_uniform_add_log` applied to
    the uniform on `s`. -/
theorem shannonH_joint_uniform_uniform
    {α β : Type*} (s : Finset α) (t : Finset β)
    (m n : ℝ) (hm : 0 < m) (hn : 0 < n)
    (hcard_s : (s.card : ℝ) = m) (hcard_t : (t.card : ℝ) = n) :
    shannonH (s ×ˢ t)
        (joint (fun _ : α => (1 : ℝ) / m) (fun _ : β => (1 : ℝ) / n))
      = Real.log m + Real.log n := by
  -- The uniform on `s` sums to 1.
  have hp_sum : ∑ _a ∈ s, (1 : ℝ) / m = 1 := by
    rw [Finset.sum_const, nsmul_eq_mul, hcard_s, mul_one_div]
    exact div_self (ne_of_gt hm)
  rw [shannonH_joint_uniform_add_log s t (fun _ => (1 : ℝ) / m) hp_sum n hn hcard_t]
  rw [shannonH_uniform_eq_log s m hm hcard_s]

/-- **Double uniform product (log-of-product form).**
    `H(U_m ⊗ U_n) = log (m · n)`. -/
theorem shannonH_joint_uniform_uniform_log_mul
    {α β : Type*} (s : Finset α) (t : Finset β)
    (m n : ℝ) (hm : 0 < m) (hn : 0 < n)
    (hcard_s : (s.card : ℝ) = m) (hcard_t : (t.card : ℝ) = n) :
    shannonH (s ×ˢ t)
        (joint (fun _ : α => (1 : ℝ) / m) (fun _ : β => (1 : ℝ) / n))
      = Real.log (m * n) := by
  rw [shannonH_joint_uniform_uniform s t m n hm hn hcard_s hcard_t]
  exact (Real.log_mul (ne_of_gt hm) (ne_of_gt hn)).symm

/-! ### 3. Delta + uniform: `H(δ ⊗ U_m) = log m` -/

/-- **Delta times uniform.** `H(δ_{a₀} ⊗ U_m) = 0 + log m = log m`. -/
theorem shannonH_joint_delta_uniform
    {α β : Type*} [DecidableEq α]
    (s : Finset α) (t : Finset β)
    (a₀ : α) (ha : a₀ ∈ s)
    (m : ℝ) (hm : 0 < m) (hcard_t : (t.card : ℝ) = m) :
    shannonH (s ×ˢ t) (joint (deltaDist a₀) (fun _ : β => (1 : ℝ) / m))
      = Real.log m := by
  have h := shannonH_joint_uniform_add_log s t (deltaDist a₀)
              (deltaDist_sum s a₀ ha) m hm hcard_t
  rw [shannonH_delta s a₀ ha] at h
  linarith

/-! ### 4. Three-factor chains -/

/-- **Three-factor chain (left-associated).**
    `H((p ⊗ U_m) ⊗ U_n) = H(p) + log m + log n`. -/
theorem shannonH_joint_chain_three
    {α β γ : Type*}
    (s : Finset α) (t : Finset β) (u : Finset γ)
    (p : α → ℝ) (hp_sum : ∑ a ∈ s, p a = 1)
    (m n : ℝ) (hm : 0 < m) (hn : 0 < n)
    (hcard_t : (t.card : ℝ) = m) (hcard_u : (u.card : ℝ) = n) :
    shannonH ((s ×ˢ t) ×ˢ u)
        (joint (joint p (fun _ : β => (1 : ℝ) / m)) (fun _ : γ => (1 : ℝ) / n))
      = shannonH s p + Real.log m + Real.log n := by
  -- First: the joint `p ⊗ U_m` sums to 1.
  have hUm_sum : ∑ _b ∈ t, (1 : ℝ) / m = 1 := by
    rw [Finset.sum_const, nsmul_eq_mul, hcard_t, mul_one_div]
    exact div_self (ne_of_gt hm)
  have hjoint_sum : ∑ ab ∈ s ×ˢ t,
      joint p (fun _ : β => (1 : ℝ) / m) ab = 1 :=
    joint_sum s t p (fun _ : β => (1 : ℝ) / m) hp_sum hUm_sum
  -- Adjoin U_n: H((p ⊗ U_m) ⊗ U_n) = H(p ⊗ U_m) + log n.
  rw [shannonH_joint_uniform_add_log (s ×ˢ t) u
        (joint p (fun _ : β => (1 : ℝ) / m)) hjoint_sum n hn hcard_u]
  -- Then expand H(p ⊗ U_m) = H(p) + log m.
  rw [shannonH_joint_uniform_add_log s t p hp_sum m hm hcard_t]

/-- **Three-factor chain (right-associated).**
    `H(p ⊗ (U_m ⊗ U_n)) = H(p) + log m + log n`. -/
theorem shannonH_joint_chain_three_right
    {α β γ : Type*}
    (s : Finset α) (t : Finset β) (u : Finset γ)
    (p : α → ℝ) (hp_sum : ∑ a ∈ s, p a = 1)
    (m n : ℝ) (hm : 0 < m) (hn : 0 < n)
    (hcard_t : (t.card : ℝ) = m) (hcard_u : (u.card : ℝ) = n) :
    shannonH (s ×ˢ (t ×ˢ u))
        (joint p (joint (fun _ : β => (1 : ℝ) / m) (fun _ : γ => (1 : ℝ) / n)))
      = shannonH s p + Real.log m + Real.log n := by
  -- The uniforms on `t` and `u` sum to 1.
  have hUm_sum : ∑ _b ∈ t, (1 : ℝ) / m = 1 := by
    rw [Finset.sum_const, nsmul_eq_mul, hcard_t, mul_one_div]
    exact div_self (ne_of_gt hm)
  have hUn_sum : ∑ _c ∈ u, (1 : ℝ) / n = 1 := by
    rw [Finset.sum_const, nsmul_eq_mul, hcard_u, mul_one_div]
    exact div_self (ne_of_gt hn)
  -- Step 1: additivity p vs (U_m ⊗ U_n).
  rw [shannonH_joint s (t ×ˢ u) p
        (joint (fun _ : β => (1 : ℝ) / m) (fun _ : γ => (1 : ℝ) / n)) hp_sum
        (joint_sum t u (fun _ : β => (1 : ℝ) / m) (fun _ : γ => (1 : ℝ) / n)
          hUm_sum hUn_sum)]
  -- Step 2: H(U_m ⊗ U_n) = log m + log n.
  rw [shannonH_joint_uniform_uniform t u m n hm hn hcard_t hcard_u]
  ring

/-! ### 5. Triple delta-cascade: `H(δ ⊗ δ ⊗ U_m) = log m` -/

/-- **Two deltas + uniform.** `H((δ_{a₀} ⊗ δ_{b₀}) ⊗ U_m) = log m`.
    Both deltas contribute zero entropy; only the uniform survives. -/
theorem shannonH_joint_delta_delta_uniform
    {α β γ : Type*} [DecidableEq α] [DecidableEq β]
    (s : Finset α) (t : Finset β) (u : Finset γ)
    (a₀ : α) (ha : a₀ ∈ s) (b₀ : β) (hb : b₀ ∈ t)
    (m : ℝ) (hm : 0 < m) (hcard_u : (u.card : ℝ) = m) :
    shannonH ((s ×ˢ t) ×ˢ u)
        (joint (joint (deltaDist a₀) (deltaDist b₀)) (fun _ : γ => (1 : ℝ) / m))
      = Real.log m := by
  -- Use the 3-factor chain with p = δ_{a₀} and m's role taken by 1 (we
  -- inline). Simpler: directly invoke chain_three with p = δ_{a₀}, then
  -- both H(δ_a) and H(δ_b) collapse to 0.
  -- Adjoin Um to the joint of two deltas:
  have h1 :
      shannonH ((s ×ˢ t) ×ˢ u)
        (joint (joint (deltaDist a₀) (deltaDist b₀)) (fun _ : γ => (1 : ℝ) / m))
        = shannonH (s ×ˢ t) (joint (deltaDist a₀) (deltaDist b₀))
            + Real.log m := by
    apply shannonH_joint_uniform_add_log
    · exact joint_sum s t (deltaDist a₀) (deltaDist b₀)
        (deltaDist_sum s a₀ ha) (deltaDist_sum t b₀ hb)
    · exact hm
    · exact hcard_u
  rw [h1, shannonH_joint_delta s t a₀ ha b₀ hb, zero_add]

/-! ### 6. PT-canonical numerical instances `m, n ∈ {2, 8, 30}` -/

/-- **`H(U_2 ⊗ U_8) = log 2 + log 8 = log 16`** on `Fin 2 × Fin 8`. -/
theorem shannonH_joint_uniform_2_8 :
    shannonH ((Finset.univ : Finset (Fin 2)) ×ˢ (Finset.univ : Finset (Fin 8)))
        (joint (fun _ => (1 : ℝ) / 2) (fun _ => (1 : ℝ) / 8))
      = Real.log 2 + Real.log 8 := by
  apply shannonH_joint_uniform_uniform _ _ 2 8
    (by norm_num) (by norm_num) (by simp) (by simp)

/-- **`H(U_2 ⊗ U_8) = log 16`** (compact form). -/
theorem shannonH_joint_uniform_2_8_log_mul :
    shannonH ((Finset.univ : Finset (Fin 2)) ×ˢ (Finset.univ : Finset (Fin 8)))
        (joint (fun _ => (1 : ℝ) / 2) (fun _ => (1 : ℝ) / 8))
      = Real.log (2 * 8) := by
  apply shannonH_joint_uniform_uniform_log_mul _ _ 2 8
    (by norm_num) (by norm_num) (by simp) (by simp)

/-- **`H(U_8 ⊗ U_30) = log 8 + log 30 = log 240`** on `Fin 8 × Fin 30`. -/
theorem shannonH_joint_uniform_8_30 :
    shannonH ((Finset.univ : Finset (Fin 8)) ×ˢ (Finset.univ : Finset (Fin 30)))
        (joint (fun _ => (1 : ℝ) / 8) (fun _ => (1 : ℝ) / 30))
      = Real.log 8 + Real.log 30 := by
  apply shannonH_joint_uniform_uniform _ _ 8 30
    (by norm_num) (by norm_num) (by simp) (by simp)

/-- **`H(δ_0 ⊗ U_30) = log 30`** on `Fin 8 × Fin 30`. -/
theorem shannonH_joint_delta_uniform_8_30 :
    shannonH ((Finset.univ : Finset (Fin 8)) ×ˢ (Finset.univ : Finset (Fin 30)))
        (joint (deltaDist (0 : Fin 8)) (fun _ : Fin 30 => (1 : ℝ) / 30))
      = Real.log 30 := by
  apply shannonH_joint_delta_uniform _ _ (0 : Fin 8) (by decide) 30
    (by norm_num) (by simp)

/-- **`H(δ_0 ⊗ δ_0 ⊗ U_30) = log 30`** on `Fin 2 × Fin 8 × Fin 30`. -/
theorem shannonH_joint_delta_delta_uniform_2_8_30 :
    shannonH (((Finset.univ : Finset (Fin 2)) ×ˢ (Finset.univ : Finset (Fin 8)))
              ×ˢ (Finset.univ : Finset (Fin 30)))
        (joint (joint (deltaDist (0 : Fin 2)) (deltaDist (0 : Fin 8)))
               (fun _ : Fin 30 => (1 : ℝ) / 30))
      = Real.log 30 := by
  apply shannonH_joint_delta_delta_uniform
    _ _ _ (0 : Fin 2) (by decide) (0 : Fin 8) (by decide) 30
    (by norm_num) (by simp)

/-! ### 7. Strict Bekenstein for the product distribution -/

/-- **Strict upper bound on `H(p ⊗ U_m)`** (parametric form).
    If the base entropy is strictly below its uniform maximum,
    `H(p) < log s.card`, then the product entropy is strictly below its
    uniform maximum `log(s.card · m)`:

    `H(p ⊗ U_m) < log(s.card · m)`.

    The hypothesis `H(p) < log s.card` is the parametric strict form
    of Bekenstein on `s`, consistent with the surrounding files. -/
theorem shannonH_joint_uniform_strict_lt
    {α β : Type*} (s : Finset α) (t : Finset β)
    (p : α → ℝ) (hp_sum : ∑ a ∈ s, p a = 1)
    (m : ℝ) (hm : 0 < m) (hcard_t : (t.card : ℝ) = m)
    (hs_card_pos : 0 < (s.card : ℝ))
    (hH_lt : shannonH s p < Real.log (s.card : ℝ)) :
    shannonH (s ×ˢ t) (joint p (fun _ : β => (1 : ℝ) / m))
      < Real.log ((s.card : ℝ) * m) := by
  -- H(p ⊗ U_m) = H(p) + log m
  rw [shannonH_joint_uniform_add_log s t p hp_sum m hm hcard_t]
  -- log(s.card · m) = log s.card + log m
  rw [Real.log_mul (ne_of_gt hs_card_pos) (ne_of_gt hm)]
  -- H(p) + log m < log s.card + log m  ↔  H(p) < log s.card
  linarith

/-! ### 8. Headline bundle -/

/-- **Headline (independent-uniform increment, full bundle).**

    For any probability distribution `p : α → ℝ` on `s : Finset α`
    summing to `1`, and any pair `(m, n)` of cardinalities `> 0` with
    `t.card = m`, `u.card = n`:

    * **Increment identity.** `H(p ⊗ U_m) − H(p) = log m`.
    * **Additive form.** `H(p ⊗ U_m) = H(p) + log m`.
    * **Double uniform.** `H(U_m ⊗ U_n) = log m + log n = log(m · n)`.
    * **Delta + uniform.** `H(δ_{a₀} ⊗ U_m) = log m`.
    * **Three-factor chain (left).**
      `H((p ⊗ U_m) ⊗ U_n) = H(p) + log m + log n`.
    * **Three-factor chain (right).**
      `H(p ⊗ (U_m ⊗ U_n)) = H(p) + log m + log n`.
    * **Triple delta-cascade.**
      `H((δ_{a₀} ⊗ δ_{b₀}) ⊗ U_n) = log n`.
-/
theorem shannonH_joint_increment_summary
    {α β γ : Type*} [DecidableEq α] [DecidableEq β]
    (s : Finset α) (t : Finset β) (u : Finset γ)
    (p : α → ℝ) (hp_sum : ∑ a ∈ s, p a = 1)
    (m n : ℝ) (hm : 0 < m) (hn : 0 < n)
    (hcard_s : (s.card : ℝ) = m)
    (hcard_t : (t.card : ℝ) = m)
    (hcard_u : (u.card : ℝ) = n)
    (a₀ : α) (ha : a₀ ∈ s) (b₀ : β) (hb : b₀ ∈ t) :
    -- 1. Increment identity (subtractive)
    shannonH (s ×ˢ t) (joint p (fun _ : β => (1 : ℝ) / m)) - shannonH s p
        = Real.log m
    -- 2. Additive form
    ∧ shannonH (s ×ˢ t) (joint p (fun _ : β => (1 : ℝ) / m))
        = shannonH s p + Real.log m
    -- 3a. Double uniform: log m + log n
    ∧ shannonH (s ×ˢ u)
        (joint (fun _ : α => (1 : ℝ) / m) (fun _ : γ => (1 : ℝ) / n))
        = Real.log m + Real.log n
    -- 3b. Double uniform: log (m · n)
    ∧ shannonH (s ×ˢ u)
        (joint (fun _ : α => (1 : ℝ) / m) (fun _ : γ => (1 : ℝ) / n))
        = Real.log (m * n)
    -- 4. Delta + uniform
    ∧ shannonH (s ×ˢ t)
        (joint (deltaDist a₀) (fun _ : β => (1 : ℝ) / m)) = Real.log m
    -- 5. Three-factor chain (left)
    ∧ shannonH ((s ×ˢ t) ×ˢ u)
        (joint (joint p (fun _ : β => (1 : ℝ) / m)) (fun _ : γ => (1 : ℝ) / n))
        = shannonH s p + Real.log m + Real.log n
    -- 6. Three-factor chain (right)
    ∧ shannonH (s ×ˢ (t ×ˢ u))
        (joint p (joint (fun _ : β => (1 : ℝ) / m) (fun _ : γ => (1 : ℝ) / n)))
        = shannonH s p + Real.log m + Real.log n
    -- 7. Triple delta-cascade
    ∧ shannonH ((s ×ˢ t) ×ˢ u)
        (joint (joint (deltaDist a₀) (deltaDist b₀)) (fun _ : γ => (1 : ℝ) / n))
        = Real.log n := by
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩
  · exact shannonH_joint_uniform_increment s t p hp_sum m hm hcard_t
  · exact shannonH_joint_uniform_add s t p hp_sum m hm hcard_t
  · exact shannonH_joint_uniform_uniform s u m n hm hn hcard_s hcard_u
  · exact shannonH_joint_uniform_uniform_log_mul s u m n hm hn hcard_s hcard_u
  · exact shannonH_joint_delta_uniform s t a₀ ha m hm hcard_t
  · exact shannonH_joint_chain_three s t u p hp_sum m n hm hn hcard_t hcard_u
  · exact shannonH_joint_chain_three_right s t u p hp_sum m n hm hn hcard_t hcard_u
  · exact shannonH_joint_delta_delta_uniform s t u a₀ ha b₀ hb n hn hcard_u

end PT.Information
