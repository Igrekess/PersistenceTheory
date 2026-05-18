/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Information.GFTIdentity
import PT.Information.GFTSpecialisations
import PT.Information.BekensteinExtensions
import PT.Information.MutualInformationBasic
import PT.Information.RelativeEntropyAdditivity
import Mathlib.Algebra.BigOperators.Group.Finset.Sigma
import Mathlib.Tactic

/-!
# KL Additivity via the GFT identity (alternative route through `H(joint)`)

This file gives a **second, independent derivation** of the KL additivity
theorem `klToUniform_joint` already proved in
`PT.Information.RelativeEntropyAdditivity`.

The original proof in `RelativeEntropyAdditivity.lean` proceeds by an
explicit pointwise expansion of `(p·q) · log((m·n)·(p·q))` (kernel
`kl_joint_pointwise_nonneg`), summed over `s ×ˢ t`. Here we obtain the
same identity as a *direct algebraic consequence* of:

* the **GFT identity** `klToUniform s m p + shannonH s p = log m`
  (`GFTIdentity.GFT_identity`),
* **additivity of Shannon entropy on a product**
  `shannonH (s ×ˢ t) (joint p q) = shannonH s p + shannonH t q`
  (`MutualInformationBasic.shannonH_joint`),
* **normalisation of the joint** `∑ joint p q = 1`
  (`MutualInformationBasic.joint_sum`),
* and the elementary `log_mul` identity.

The argument is "morally": at independence the mutual information vanishes,
the entropy `H` is additive, and the GFT identity then forces `D_KL` to be
additive too (with the convention `U_{mn} = U_m ⊗ U_n`).

Concretely, applying GFT three times:

```
log(m·n) = D_KL(joint p q ‖ U_{mn}) + H(joint p q)
         = D_KL(joint p q ‖ U_{mn}) + H(p) + H(q)             -- shannonH_joint
log m    = D_KL(p ‖ U_m) + H(p)
log n    = D_KL(q ‖ U_n) + H(q)
```

Adding the last two and substituting `log(m·n) = log m + log n`:

```
D_KL(joint p q ‖ U_{mn}) = log m + log n − H(p) − H(q)
                         = (log m − H(p)) + (log n − H(q))
                         = D_KL(p ‖ U_m) + D_KL(q ‖ U_n).
```

The two proofs (this one and `klToUniform_joint`) give the same equality;
the comparison theorem `klToUniform_joint_via_GFT_eq_klToUniform_joint`
records this redundancy as a trivial equality between the two routes.

## Reference

Monograph Chapter 4, `\label{prop:KL_additivity}`, alternative derivation
via the mutual-information vanishing. M4 article.
-/

namespace PT.Information

open Real Finset

/-! ### Intermediate: GFT in rearranged "D_KL = log m − H" form -/

/-- **GFT rearranged.** From `klToUniform s m p + shannonH s p = log m`
    we read off the closed form `D_KL(P ‖ U_m) = log m − H(P)`. -/
theorem klToUniform_eq_log_sub_H
    {α : Type*} (s : Finset α) (m : ℝ) (hm : 0 < m)
    (p : α → ℝ) (hp_nonneg : ∀ r ∈ s, 0 ≤ p r)
    (hp_sum : ∑ r ∈ s, p r = 1) :
    klToUniform s m p = Real.log m - shannonH s p := by
  have h := GFT_identity s m hm p hp_nonneg hp_sum
  linarith

/-- **GFT rearranged, applied on the product set `s ×ˢ t` at modulus `m·n`.**
    `klToUniform (s ×ˢ t) (m·n) (joint p q) = log(m·n) − H(joint p q)`. -/
theorem klToUniform_joint_eq_log_mul_sub_H
    {α β : Type*} (s : Finset α) (t : Finset β)
    (m n : ℝ) (hm : 0 < m) (hn : 0 < n)
    (p : α → ℝ) (q : β → ℝ)
    (hp_nonneg : ∀ a ∈ s, 0 ≤ p a) (hq_nonneg : ∀ b ∈ t, 0 ≤ q b)
    (hp_sum : ∑ a ∈ s, p a = 1) (hq_sum : ∑ b ∈ t, q b = 1) :
    klToUniform (s ×ˢ t) (m * n) (joint p q)
      = Real.log (m * n) - shannonH (s ×ˢ t) (joint p q) := by
  apply klToUniform_eq_log_sub_H (s ×ˢ t) (m * n) (mul_pos hm hn) (joint p q)
  · exact joint_nonneg s t p q hp_nonneg hq_nonneg
  · exact joint_sum s t p q hp_sum hq_sum

/-! ### Main theorem: KL additivity via the GFT identity -/

/-- **KL additivity on tensor products — derived via the GFT identity.**

    Same statement as `klToUniform_joint` (proved in
    `RelativeEntropyAdditivity.lean` by a direct pointwise expansion); here
    the proof goes through:

    1. `GFT` on `(s ×ˢ t, m · n)`: `D_KL = log(m·n) − H(joint)`,
    2. `shannonH_joint`: `H(joint p q) = H(p) + H(q)`,
    3. `log_mul`: `log(m · n) = log m + log n`,
    4. `GFT` on `(s, m)` and `(t, n)`: `D_KL = log m − H(p)` and `log n − H(q)`. -/
theorem klToUniform_joint_via_GFT
    {α β : Type*} (s : Finset α) (t : Finset β)
    (m n : ℝ) (hm : 0 < m) (hn : 0 < n)
    (p : α → ℝ) (q : β → ℝ)
    (hp_nonneg : ∀ a ∈ s, 0 ≤ p a) (hq_nonneg : ∀ b ∈ t, 0 ≤ q b)
    (hp_sum : ∑ a ∈ s, p a = 1) (hq_sum : ∑ b ∈ t, q b = 1) :
    klToUniform (s ×ˢ t) (m * n) (joint p q)
      = klToUniform s m p + klToUniform t n q := by
  -- Step 1: rewrite LHS via GFT on the product.
  rw [klToUniform_joint_eq_log_mul_sub_H s t m n hm hn p q
        hp_nonneg hq_nonneg hp_sum hq_sum]
  -- Step 2: `H(joint p q) = H(p) + H(q)`.
  rw [shannonH_joint s t p q hp_sum hq_sum]
  -- Step 3: `log(m * n) = log m + log n`.
  rw [Real.log_mul (ne_of_gt hm) (ne_of_gt hn)]
  -- Step 4: rewrite each RHS factor via GFT.
  rw [klToUniform_eq_log_sub_H s m hm p hp_nonneg hp_sum,
      klToUniform_eq_log_sub_H t n hn q hq_nonneg hq_sum]
  ring

/-! ### Cross-check: the two proofs coincide -/

/-- **The two derivations of KL additivity coincide.** Equality between the
    pointwise-kernel proof (`klToUniform_joint`) and the GFT-route proof
    (`klToUniform_joint_via_GFT`). Both are theorems about the same equality
    of real numbers, hence trivially equal — but this lemma makes the
    redundancy explicit and serves as a sanity check that the two routes
    yield genuinely the same statement. -/
theorem klToUniform_joint_via_GFT_eq_klToUniform_joint
    {α β : Type*} (s : Finset α) (t : Finset β)
    (m n : ℝ) (hm : 0 < m) (hn : 0 < n)
    (p : α → ℝ) (q : β → ℝ)
    (hp_nonneg : ∀ a ∈ s, 0 ≤ p a) (hq_nonneg : ∀ b ∈ t, 0 ≤ q b)
    (hp_sum : ∑ a ∈ s, p a = 1) (hq_sum : ∑ b ∈ t, q b = 1) :
    klToUniform_joint_via_GFT s t m n hm hn p q
        hp_nonneg hq_nonneg hp_sum hq_sum
      = klToUniform_joint s t m n hm hn p q
          hp_nonneg hq_nonneg hp_sum hq_sum := by
  -- Both sides are equality proofs (the same real-number equality),
  -- so by proof irrelevance they coincide.
  rfl

/-! ### Corner instances re-derived through the GFT route -/

/-- **Delta corner via GFT route.** `D_KL(δ_(a₀, b₀) ‖ U_{m·n}) = log(m·n)`. -/
theorem klToUniform_joint_delta_via_GFT
    {α β : Type*} [DecidableEq α] [DecidableEq β]
    (s : Finset α) (t : Finset β)
    (m n : ℝ) (hm : 0 < m) (hn : 0 < n)
    (a₀ : α) (ha : a₀ ∈ s) (b₀ : β) (hb : b₀ ∈ t) :
    klToUniform (s ×ˢ t) (m * n) (joint (deltaDist a₀) (deltaDist b₀))
      = Real.log m + Real.log n := by
  rw [klToUniform_joint_via_GFT s t m n hm hn
        (deltaDist a₀) (deltaDist b₀)
        (fun r _ => deltaDist_nonneg a₀ r)
        (fun r _ => deltaDist_nonneg b₀ r)
        (deltaDist_sum s a₀ ha) (deltaDist_sum t b₀ hb)]
  rw [klToUniform_delta s m hm a₀ ha, klToUniform_delta t n hn b₀ hb]

/-- **Uniform corner via GFT route.** Both sides are zero. -/
theorem klToUniform_joint_uniform_via_GFT
    {α β : Type*} (s : Finset α) (t : Finset β)
    (m n : ℝ) (hm : 0 < m) (hn : 0 < n)
    (hcard_s : (s.card : ℝ) = m) (hcard_t : (t.card : ℝ) = n) :
    klToUniform (s ×ˢ t) (m * n)
        (joint (fun _ : α => (1 : ℝ) / m) (fun _ : β => (1 : ℝ) / n)) = 0 := by
  have hp_sum : ∑ _a ∈ s, ((1 : ℝ) / m) = 1 := by
    rw [Finset.sum_const, nsmul_eq_mul, hcard_s, mul_one_div]
    exact div_self (ne_of_gt hm)
  have hq_sum : ∑ _b ∈ t, ((1 : ℝ) / n) = 1 := by
    rw [Finset.sum_const, nsmul_eq_mul, hcard_t, mul_one_div]
    exact div_self (ne_of_gt hn)
  rw [klToUniform_joint_via_GFT s t m n hm hn
        (fun _ => 1 / m) (fun _ => 1 / n)
        (fun _ _ => by positivity) (fun _ _ => by positivity)
        hp_sum hq_sum]
  rw [klToUniform_uniform s m hm, klToUniform_uniform t n hn]
  ring

/-! ### Headline (alternative-route summary) -/

/-- **Headline (KL additivity via the GFT identity).**

    The same `D_KL(P ⊗ Q ‖ U_m ⊗ U_n) = D_KL(P ‖ U_m) + D_KL(Q ‖ U_n)`
    is obtained by a different route — through the GFT identity and the
    additivity of Shannon entropy on a product — and we record the trivial
    equality with the pointwise-kernel route from
    `RelativeEntropyAdditivity.lean`. -/
theorem KL_additivity_via_GFT_summary
    {α β : Type*} (s : Finset α) (t : Finset β)
    (m n : ℝ) (hm : 0 < m) (hn : 0 < n)
    (p : α → ℝ) (q : β → ℝ)
    (hp_nonneg : ∀ a ∈ s, 0 ≤ p a) (hq_nonneg : ∀ b ∈ t, 0 ≤ q b)
    (hp_sum : ∑ a ∈ s, p a = 1) (hq_sum : ∑ b ∈ t, q b = 1) :
    -- GFT rearrangement on the product set:
    klToUniform (s ×ˢ t) (m * n) (joint p q)
        = Real.log (m * n) - shannonH (s ×ˢ t) (joint p q)
    -- Main additivity (GFT route):
    ∧ klToUniform (s ×ˢ t) (m * n) (joint p q)
        = klToUniform s m p + klToUniform t n q
    -- Same identity as the pointwise-kernel route:
    ∧ klToUniform_joint_via_GFT s t m n hm hn p q
          hp_nonneg hq_nonneg hp_sum hq_sum
        = klToUniform_joint s t m n hm hn p q
            hp_nonneg hq_nonneg hp_sum hq_sum :=
  ⟨klToUniform_joint_eq_log_mul_sub_H s t m n hm hn p q
      hp_nonneg hq_nonneg hp_sum hq_sum,
   klToUniform_joint_via_GFT s t m n hm hn p q
      hp_nonneg hq_nonneg hp_sum hq_sum,
   klToUniform_joint_via_GFT_eq_klToUniform_joint s t m n hm hn p q
      hp_nonneg hq_nonneg hp_sum hq_sum⟩

end PT.Information
