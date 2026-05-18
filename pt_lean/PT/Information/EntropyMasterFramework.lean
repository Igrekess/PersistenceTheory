/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Information.GFTIdentity
import PT.Information.GFTSpecialisations
import PT.Information.BekensteinBound
import PT.Information.BekensteinExtensions
import PT.Information.EntropyBoundsTight
import PT.Information.EntropyMonotonicity
import PT.Information.ShannonEntropyConcavity
import PT.Information.MutualInformationBasic
import PT.Information.RelativeEntropyAdditivity
import PT.Information.CrossEntropyBound
import PT.Information.ConditionalEntropy
import Mathlib.Tactic

/-!
# Entropy Master Framework — Unified PT view of GFT, Bekenstein, additivity,
  concavity, mutual information, conditional entropy, and the chain rule.

This file is a **pure synthesis aggregator**: every theorem stated here is
*re-exposed* from an existing module of `PT.Information`. No new computation
is performed. The purpose is to bundle the principal information-theoretic
results of the PT corpus into a single, high-level master view that can be
cited as one block in downstream PT theorems (Bekenstein-saturation arguments,
Kähler-Fisher entropy potentials, RH bridges, etc.).

## Structure

The master framework is organised into **six pillars**, each captured by a
self-contained sub-bundle, and then assembled into the top-level
`entropy_master_framework_summary`:

* **Pillar 1 — Master GFT identity.** `master_GFT_identity`
  re-exposes `GFT_identity` (the algebraic partition
  `log m = D_KL(P ‖ U_m) + H(P)` of the GFT budget) as the central PT
  theorem of the framework.
* **Pillar 2 — Master Bekenstein.** `master_bekenstein_full` bundles
  the Bekenstein upper bound `D_KL ≤ log m`, its saturation at deltas
  `D_KL(δ) = log m`, its vanishing at uniform `D_KL(U_m) = 0`, the
  entropy bound `H ≤ log m`, and saturation at uniform `H(U_m) = log m`.
* **Pillar 3 — Master tensor additivity.** `master_tensor_additivity`
  bundles Shannon entropy additivity `H(p ⊗ q) = H(p) + H(q)`, KL
  additivity `D_KL(p ⊗ q ‖ U_m ⊗ U_n) = D_KL(p ‖ U_m) + D_KL(q ‖ U_n)`,
  and `I(p, q) = 0` at independence.
* **Pillar 4 — Master concavity.** `master_concavity` bundles concavity
  of `negMulLog`, concavity of `shannonH` (two-distribution Jensen),
  and the Jensen application to the Bekenstein bound.
* **Pillar 5 — Master chain rule.** `master_chain_rule` re-exposes
  `H(X, Y) = H(Y) + H(X | Y)`, plus `H(X | Y) = H(X)` at independence.
* **Pillar 6 — Master cross-entropy / Gibbs bridge.**
  `master_cross_entropy_bridge` re-exposes
  `H(p, U_m) = log m` and `D_KL(p ‖ U_m) = log m - H(p)`.

The final theorem `entropy_master_framework_summary` aggregates the six
pillars into a single conjunction. Because every component is already
established in its source module with `0 sorry`, this file inherits the
same guarantee.

## Reference

Monograph Ch04. M4 article. Companion modules `PT.Information.*`.
-/

namespace PT.Information

open Real Finset

/-! ## Pillar 1 — Master GFT identity -/

/-- **Master GFT identity (re-exposure of `GFT_identity`).**
    The fundamental exact algebraic partition of the GFT budget:

    $$\log m \;=\; D_{KL}(P \,\Vert\, U_m) \;+\; H(P).$$

    For any probability distribution `p` on a finite set `s` of cardinality
    `m`, the sum of the KL divergence to the uniform and the Shannon
    entropy equals `log m` exactly. -/
theorem master_GFT_identity {α : Type*} (s : Finset α) (m : ℝ) (hm : 0 < m)
    (p : α → ℝ) (hp_nonneg : ∀ r ∈ s, 0 ≤ p r)
    (hp_sum : ∑ r ∈ s, p r = 1) :
    klToUniform s m p + shannonH s p = Real.log m :=
  GFT_identity s m hm p hp_nonneg hp_sum

/-! ## Pillar 2 — Master Bekenstein (bound, saturation, vanishing) -/

/-- **Master Bekenstein (full bundle).**

    For any probability distribution `p : α → ℝ` on a finite set `s` of
    cardinality `m > 0` with `s.card = m`:

    1. **Bekenstein bound.** `D_KL(P ‖ U_m) ≤ log m`.
    2. **Saturation at delta.** `D_KL(δ_{r₀} ‖ U_m) = log m`.
    3. **Vanishing at uniform.** `D_KL(U_m ‖ U_m) = 0`.
    4. **Entropy bound (Bekenstein dual).** `H(U_m) = log m`, so the
       Shannon entropy reaches its maximum `log m` on the uniform.
    5. **Entropy vanishing at delta.** `H(δ_{r₀}) = 0`.
    6. **Both corners satisfy the GFT identity.** -/
theorem master_bekenstein_full {α : Type*} [DecidableEq α]
    (s : Finset α) (m : ℝ) (hm : 0 < m) (hcard : (s.card : ℝ) = m)
    (r₀ : α) (hr : r₀ ∈ s)
    (p : α → ℝ)
    (hp_nonneg : ∀ r ∈ s, 0 ≤ p r) (hp_le_one : ∀ r ∈ s, p r ≤ 1)
    (hp_sum : ∑ r ∈ s, p r = 1) :
    -- 1. Bekenstein bound
    klToUniform s m p ≤ Real.log m
    -- 2. Saturation at delta
    ∧ klToUniform s m (deltaDist r₀) = Real.log m
    -- 3. Vanishing at uniform
    ∧ klToUniform s m (fun _ => (1 : ℝ) / m) = 0
    -- 4. Entropy saturation at uniform
    ∧ shannonH s (fun _ => (1 : ℝ) / m) = Real.log m
    -- 5. Entropy vanishing at delta
    ∧ shannonH s (deltaDist r₀) = 0
    -- 6. GFT partition at both corners
    ∧ klToUniform s m (deltaDist r₀) + shannonH s (deltaDist r₀) = Real.log m
    ∧ klToUniform s m (fun _ => (1 : ℝ) / m)
        + shannonH s (fun _ => (1 : ℝ) / m) = Real.log m :=
  ⟨bekenstein_bound s m hm p hp_nonneg hp_le_one hp_sum,
   bekenstein_saturated_at_delta s m hm r₀ hr,
   klToUniform_uniform s m hm,
   shannonH_uniform_eq_log s m hm hcard,
   shannonH_delta s r₀ hr,
   GFT_at_delta s m hm r₀ hr,
   GFT_at_uniform s m hm hcard⟩

/-! ## Pillar 3 — Master tensor additivity (Shannon, KL, mutual info) -/

/-- **Master tensor additivity (full bundle).**

    For two probability distributions `p` on `s` (cardinality `m`) and
    `q` on `t` (cardinality `n`) with `m, n > 0`:

    1. **Shannon additivity.** `H(p ⊗ q) = H(p) + H(q)`.
    2. **KL additivity.** `D_KL(p ⊗ q ‖ U_m ⊗ U_n)
                          = D_KL(p ‖ U_m) + D_KL(q ‖ U_n)`.
    3. **Mutual information vanishes at independence.** `I(p, q) = 0`.
    4. **Joint normalisation.** `∑ joint p q = 1`. -/
theorem master_tensor_additivity
    {α β : Type*} (s : Finset α) (t : Finset β)
    (m n : ℝ) (hm : 0 < m) (hn : 0 < n)
    (p : α → ℝ) (q : β → ℝ)
    (hp_nonneg : ∀ a ∈ s, 0 ≤ p a) (hq_nonneg : ∀ b ∈ t, 0 ≤ q b)
    (hp_sum : ∑ a ∈ s, p a = 1) (hq_sum : ∑ b ∈ t, q b = 1) :
    -- 1. Shannon additivity
    shannonH (s ×ˢ t) (joint p q) = shannonH s p + shannonH t q
    -- 2. KL additivity
    ∧ klToUniform (s ×ˢ t) (m * n) (joint p q)
        = klToUniform s m p + klToUniform t n q
    -- 3. Mutual information vanishes at independence
    ∧ mutualInfo s t p q = 0
    -- 4. Joint normalisation
    ∧ (∑ ab ∈ s ×ˢ t, joint p q ab) = 1 :=
  ⟨shannonH_joint s t p q hp_sum hq_sum,
   klToUniform_joint s t m n hm hn p q hp_nonneg hq_nonneg hp_sum hq_sum,
   mutualInfo_indep_eq_zero s t p q hp_sum hq_sum,
   joint_sum s t p q hp_sum hq_sum⟩

/-! ## Pillar 4 — Master concavity (negMulLog, shannonH, Jensen Bekenstein) -/

/-- **Master concavity (full bundle).**

    1. **Pointwise concavity of `negMulLog`.** The map `x ↦ -x log x`
       is concave on `[0, ∞)`.
    2. **Two-distribution concavity of `shannonH`.** For non-negative
       distributions `p, q : α → ℝ` and weights `a, b ≥ 0` with
       `a + b = 1`,
       `a · H(p) + b · H(q) ≤ H(a · p + b · q)`.
    3. **Bekenstein via Jensen.** `H(p) ≤ log m` follows from Jensen
       applied to `negMulLog` with uniform weights. -/
theorem master_concavity
    {α : Type*} (s : Finset α) (m : ℝ) (hm : 0 < m)
    (hcard : (s.card : ℝ) = m)
    (p q : α → ℝ)
    (hp_nonneg : ∀ r ∈ s, 0 ≤ p r) (hq_nonneg : ∀ r ∈ s, 0 ≤ q r)
    (hp_sum : ∑ r ∈ s, p r = 1)
    {a b : ℝ} (ha : 0 ≤ a) (hb : 0 ≤ b) (hab : a + b = 1) :
    -- 1. Pointwise concavity of negMulLog
    ConcaveOn ℝ (Set.Ici (0 : ℝ)) Real.negMulLog
    -- 2. Two-distribution concavity of shannonH
    ∧ a * shannonH s p + b * shannonH s q
        ≤ shannonH s (fun r => a * p r + b * q r)
    -- 3. Bekenstein via Jensen
    ∧ shannonH s p ≤ Real.log m :=
  ⟨concaveOn_negMulLog_PT,
   shannonH_concave_two s p q hp_nonneg hq_nonneg ha hb hab,
   shannonH_le_log_card_of_jensen s m hm hcard p hp_nonneg hp_sum⟩

/-! ## Pillar 5 — Master chain rule (joint = marginal + conditional) -/

/-- **Master chain rule (full bundle).**

    For two probability distributions `p` on `s` and `q` on `t` with
    `∑ p = ∑ q = 1`:

    1. **Chain rule.** `H(X, Y) = H(Y) + H(X | Y)`.
    2. **Independence ⇒ `H(X | Y) = H(X)`.** -/
theorem master_chain_rule
    {α β : Type*} (s : Finset α) (t : Finset β)
    (p : α → ℝ) (q : β → ℝ)
    (hp_sum : ∑ a ∈ s, p a = 1) (hq_sum : ∑ b ∈ t, q b = 1) :
    -- 1. Chain rule
    shannonH (s ×ˢ t) (joint p q)
        = shannonH t q + condEntropy s t (joint p q) q
    -- 2. Independence reduction
    ∧ condEntropy s t (joint p q) q = shannonH s p :=
  ⟨shannonH_joint_eq_add_condEntropy s t (joint p q) q,
   condEntropy_indep s t p q hp_sum hq_sum⟩

/-! ## Pillar 6 — Master cross-entropy / Gibbs bridge -/

/-- **Master cross-entropy bridge (full bundle).**

    1. **Self-cross is Shannon.** `H(p, p) = H(p)`.
    2. **Uniform reference.** `H(p, U_m) = log m`.
    3. **GFT bridge.** `D_KL(p ‖ U_m) = log m - H(p)`. -/
theorem master_cross_entropy_bridge
    {α : Type*} (s : Finset α) (m : ℝ) (hm : 0 < m)
    (p : α → ℝ) (hp_nonneg : ∀ r ∈ s, 0 ≤ p r)
    (hp_sum : ∑ r ∈ s, p r = 1) :
    -- 1. Self-cross
    crossEntropy s p p = shannonH s p
    -- 2. Uniform reference
    ∧ crossEntropy s p (fun _ => (1 : ℝ) / m) = Real.log m
    -- 3. GFT bridge
    ∧ klToUniform s m p = Real.log m - shannonH s p :=
  ⟨crossEntropy_self s p,
   crossEntropy_uniform s m hm p hp_sum,
   klToUniform_eq_logm_sub_shannonH s m hm p hp_nonneg hp_sum⟩

/-! ## Top-level master framework -/

/-- **Entropy Master Framework — top-level synthesis.**

    For any two probability distributions `p` on `s` (cardinality `m`) and
    `q` on `t` (cardinality `n`), with `m, n > 0`, `0 ≤ p ≤ 1`,
    `∑ p = ∑ q = 1`, and a chosen residue class `r₀ ∈ s`, the six pillars
    of the PT information framework hold jointly:

    * **(GFT)** The fundamental partition `log m = D_KL + H` is exact.
    * **(Bekenstein)** `D_KL ≤ log m` with saturation at deltas, vanishing
      at uniform, and the dual `H ≤ log m`.
    * **(Tensor additivity)** `H` and `D_KL` are additive under independent
      products; `I = 0` at independence.
    * **(Concavity)** `negMulLog` is concave; `H` inherits two-distribution
      concavity; Jensen recovers the Bekenstein bound.
    * **(Chain rule)** `H(X, Y) = H(Y) + H(X|Y)`; under independence
      `H(X|Y) = H(X)`.
    * **(Cross-entropy bridge)** `H(p, p) = H(p)` and
      `D_KL(p ‖ U_m) = log m - H(p)`.

    This is a **pure aggregator** of theorems already established in
    `PT.Information.*`. -/
theorem entropy_master_framework_summary
    {α β : Type*} [DecidableEq α]
    (s : Finset α) (t : Finset β)
    (m n : ℝ) (hm : 0 < m) (hn : 0 < n)
    (hcard_s : (s.card : ℝ) = m) (_hcard_t : (t.card : ℝ) = n)
    (r₀ : α) (hr : r₀ ∈ s)
    (p : α → ℝ) (q : β → ℝ)
    (hp_nonneg : ∀ r ∈ s, 0 ≤ p r) (hp_le_one : ∀ r ∈ s, p r ≤ 1)
    (hp_sum : ∑ r ∈ s, p r = 1)
    (hq_nonneg : ∀ b ∈ t, 0 ≤ q b) (hq_sum : ∑ b ∈ t, q b = 1) :
    -- Pillar 1: GFT identity
    klToUniform s m p + shannonH s p = Real.log m
    -- Pillar 2a: Bekenstein bound
    ∧ klToUniform s m p ≤ Real.log m
    -- Pillar 2b: delta saturates
    ∧ klToUniform s m (deltaDist r₀) = Real.log m
    -- Pillar 2c: uniform vanishes KL
    ∧ klToUniform s m (fun _ => (1 : ℝ) / m) = 0
    -- Pillar 2d: uniform attains entropy max
    ∧ shannonH s (fun _ => (1 : ℝ) / m) = Real.log m
    -- Pillar 2e: delta has zero entropy
    ∧ shannonH s (deltaDist r₀) = 0
    -- Pillar 3a: Shannon additivity
    ∧ shannonH (s ×ˢ t) (joint p q) = shannonH s p + shannonH t q
    -- Pillar 3b: KL additivity
    ∧ klToUniform (s ×ˢ t) (m * n) (joint p q)
        = klToUniform s m p + klToUniform t n q
    -- Pillar 3c: I = 0 at independence
    ∧ mutualInfo s t p q = 0
    -- Pillar 4a: negMulLog concave
    ∧ ConcaveOn ℝ (Set.Ici (0 : ℝ)) Real.negMulLog
    -- Pillar 4b: Bekenstein via Jensen
    ∧ shannonH s p ≤ Real.log m
    -- Pillar 5a: chain rule
    ∧ shannonH (s ×ˢ t) (joint p q)
        = shannonH t q + condEntropy s t (joint p q) q
    -- Pillar 5b: independence reduces conditional to marginal
    ∧ condEntropy s t (joint p q) q = shannonH s p
    -- Pillar 6a: self-cross = Shannon
    ∧ crossEntropy s p p = shannonH s p
    -- Pillar 6b: uniform cross-entropy
    ∧ crossEntropy s p (fun _ => (1 : ℝ) / m) = Real.log m
    -- Pillar 6c: GFT bridge via cross-entropy
    ∧ klToUniform s m p = Real.log m - shannonH s p :=
  ⟨GFT_identity s m hm p hp_nonneg hp_sum,
   bekenstein_bound s m hm p hp_nonneg hp_le_one hp_sum,
   bekenstein_saturated_at_delta s m hm r₀ hr,
   klToUniform_uniform s m hm,
   shannonH_uniform_eq_log s m hm hcard_s,
   shannonH_delta s r₀ hr,
   shannonH_joint s t p q hp_sum hq_sum,
   klToUniform_joint s t m n hm hn p q hp_nonneg hq_nonneg hp_sum hq_sum,
   mutualInfo_indep_eq_zero s t p q hp_sum hq_sum,
   concaveOn_negMulLog_PT,
   shannonH_le_log_card_of_jensen s m hm hcard_s p hp_nonneg hp_sum,
   shannonH_joint_eq_add_condEntropy s t (joint p q) q,
   condEntropy_indep s t p q hp_sum hq_sum,
   crossEntropy_self s p,
   crossEntropy_uniform s m hm p hp_sum,
   klToUniform_eq_logm_sub_shannonH s m hm p hp_nonneg hp_sum⟩

end PT.Information
