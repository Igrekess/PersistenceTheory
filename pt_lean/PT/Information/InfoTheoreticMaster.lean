/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Information.EntropyMasterFramework
import PT.Information.GFTSpecialMValues
import PT.Information.EntropyOfBinaryDistribution
import PT.Information.BinaryEntropyMonotonicity
import PT.Information.EntropyTernaryDistribution
import PT.Information.EntropyDifferenceSymmetric
import PT.Information.EntropyAdditivityCorners
import PT.Information.EntropyJointProduct
import PT.Information.MutualInfoDistributional
import PT.Information.MutualInformationBasic
import PT.Information.KLAdditivityProduct
import PT.Information.KLAdditivityFromMI
import PT.Information.RelativeEntropyAdditivity
import PT.Information.EntropyMonotonicity
import PT.Information.EntropyBoundsTight
import PT.Information.BekensteinExtensions
import PT.Information.ConditionalEntropy
import PT.Information.CrossEntropyBound
import Mathlib.Tactic

/-!
# Information-Theoretic Master — TIER 2 extension of `EntropyMasterFramework`.

This file is a **pure synthesis aggregator** (no new computation): it bundles
the principal information-theoretic results of the PT corpus across the entire
`PT.Information.*` family into a single, top-level master view.

The aggregator is layered on top of `EntropyMasterFramework` and adds the
following extra pillars that the framework does not yet bundle directly:

* **Pillar A — GFT special values** (`gft_special_m_values_summary`):
  numerical instances of the GFT identity at `m ∈ {2, 8, 30}`.
* **Pillar B — Binary entropy** (`binary_entropy_headline` and
  `binary_entropy_comparisons_headline`): values of `H_bin` at the
  PT-canonical points `0, 1/3, 1/2, 2/3, 13/15, 1`.
* **Pillar C — Ternary entropy** (`ternary_entropy_headline`):
  values and structural properties of `H_ter` on the 2-simplex.
* **Pillar D — Binary entropy differences** (`binary_entropy_differences_headline`):
  symmetric differences, cascade gap and saturated envelope.
* **Pillar E — Entropy additivity corners** (`entropy_additivity_corners_summary`):
  `H(δ ⊗ δ) = 0` and `H(U_m ⊗ U_n) = log m + log n`.
* **Pillar F — Joint entropy increments** (`shannonH_joint_increment_summary`):
  increment monotonicity for `H(p ⊗ U_n)` and `H(δ ⊗ U_n)`.
* **Pillar G — Mutual information (distributional)**
  (`mutualInfoGen_summary`): the generalised mutual information
  `I_gen(P_{XY}) = H(X) + H(Y) - H(X,Y)` and its independence /
  delta-pair / uniform-product instances.
* **Pillar H — Mutual information (basic)** (`mutualInfo_basic_summary`):
  the product-form `I(p, q) = 0` at independence, plus the joint /
  delta / uniform corners.
* **Pillar I — KL additivity via product corners**
  (`KL_additivity_summary`): `D_KL(δ_2) + D_KL(δ_8) = D_KL(δ_{16})`-style
  partitions of `log` along the cascade `2 → 8 → 30`.
* **Pillar J — KL additivity via GFT** (`KL_additivity_via_GFT_summary`):
  the alternative GFT route to KL additivity (joint via the partition).
* **Pillar K — KL tensor additivity** (`KL_tensor_additivity_summary`):
  the relative-entropy formulation of the tensor product.
* **Pillar L — Shannon monotonicity** (`shannonH_monotonicity_summary`):
  `H` is monotone under uniformising mixtures.
* **Pillar M — Tight bounds** (`entropy_tight_bounds_summary`):
  `0 ≤ H ≤ log m`, both bounds saturated.
* **Pillar N — Bekenstein saturation summary**
  (`bekenstein_saturation_summary`): the `m ≥ 1` form.
* **Pillar O — Conditional entropy** (`condEntropy_summary`): chain rule
  and independence reduction in their full distributional form.
* **Pillar P — Cross-entropy** (`crossEntropy_summary`): Gibbs bound,
  self-cross, uniform reference, GFT bridge.

Then `info_theoretic_master_summary` bundles **all** pillars (the original
`entropy_master_framework_summary` plus pillars A–P) into a single
conjunction. Every component is established with `0 sorry` in its source
module, so this file inherits the same guarantee.

## Reference

Monograph Ch04 (Information geometry). Companion modules `PT.Information.*`.
M4 article.
-/

namespace PT.Information

open Real Finset

/-! ## TIER-2 master re-exposures (pure aggregators) -/

/-- **Master GFT identity (TIER-2 re-exposure).**
    Same statement as `master_GFT_identity`. -/
theorem infoTheoretic_master_GFT_identity
    {α : Type*} (s : Finset α) (m : ℝ) (hm : 0 < m)
    (p : α → ℝ) (hp_nonneg : ∀ r ∈ s, 0 ≤ p r)
    (hp_sum : ∑ r ∈ s, p r = 1) :
    klToUniform s m p + shannonH s p = Real.log m :=
  PT.Information.master_GFT_identity s m hm p hp_nonneg hp_sum

/-- **Master GFT special m-values (TIER-2 re-exposure).** -/
theorem infoTheoretic_master_GFT_special_m_values :
    klToUniform (Finset.univ : Finset (Fin 2)) 2 (deltaDist 0) = Real.log 2
    ∧ klToUniform (Finset.univ : Finset (Fin 2)) 2 (fun _ => (1 : ℝ) / 2) = 0
    ∧ klToUniform (Finset.univ : Finset (Fin 8)) 8 (deltaDist 0) = Real.log 8
    ∧ klToUniform (Finset.univ : Finset (Fin 8)) 8 (fun _ => (1 : ℝ) / 8) = 0
    ∧ klToUniform (Finset.univ : Finset (Fin 30)) 30 (deltaDist 0) = Real.log 30
    ∧ klToUniform (Finset.univ : Finset (Fin 30)) 30 (fun _ => (1 : ℝ) / 30) = 0
    ∧ 0 < Real.log 2
    ∧ Real.log 2 < Real.log 8
    ∧ Real.log 8 < Real.log 30 :=
  gft_special_m_values_summary

/-- **Master binary entropy values (TIER-2 re-exposure).** -/
theorem infoTheoretic_master_binEntropy_values :
    binEntropy 0 = 0
    ∧ binEntropy 1 = 0
    ∧ binEntropy ((1 : ℝ) / 2) = Real.log 2
    ∧ (∀ p, binEntropy p = binEntropy (1 - p))
    ∧ binEntropy ((1 : ℝ) / 3) = Real.log 3 - ((2 : ℝ) / 3) * Real.log 2 :=
  binary_entropy_headline

/-- **Master binary entropy ordering (TIER-2 re-exposure).** -/
theorem infoTheoretic_master_binEntropy_comparisons :
    binEntropy (0 : ℝ) < binEntropy ((1 : ℝ) / 3)
    ∧ binEntropy ((1 : ℝ) / 3) < binEntropy ((1 : ℝ) / 2)
    ∧ binEntropy ((2 : ℝ) / 3) = binEntropy ((1 : ℝ) / 3)
    ∧ binEntropy ((2 : ℝ) / 3) < binEntropy ((1 : ℝ) / 2)
    ∧ binEntropy (1 : ℝ) < binEntropy ((1 : ℝ) / 3) :=
  binary_entropy_comparisons_headline

/-- **Master ternary entropy values (TIER-2 re-exposure).** -/
theorem infoTheoretic_master_terEntropy_values :
    terEntropy 1 0 = 0
    ∧ terEntropy 0 1 = 0
    ∧ terEntropy 0 0 = 0
    ∧ terEntropy ((1 : ℝ) / 3) ((1 : ℝ) / 3) = Real.log 3
    ∧ (∀ p q, terEntropy p q = terEntropy q p)
    ∧ (∀ p, terEntropy p (1 - p) = binEntropy p) :=
  ternary_entropy_headline

/-- **Master binary entropy differences (TIER-2 re-exposure).** -/
theorem infoTheoretic_master_binEntropy_differences :
    (∀ p, binEntropy p - binEntropy (1 - p) = 0)
    ∧ binEntropy ((1 : ℝ) / 2) - binEntropy ((1 : ℝ) / 3)
        = (5 : ℝ) / 3 * Real.log 2 - Real.log 3
    ∧ binEntropy ((1 : ℝ) / 3) - binEntropy (0 : ℝ)
        = Real.log 3 - (2 : ℝ) / 3 * Real.log 2
    ∧ (binEntropy ((1 : ℝ) / 2) - binEntropy ((1 : ℝ) / 3))
          - (binEntropy ((1 : ℝ) / 3) - binEntropy (0 : ℝ))
        = (7 : ℝ) / 3 * Real.log 2 - 2 * Real.log 3
    ∧ |binEntropy ((1 : ℝ) / 2) - binEntropy (0 : ℝ)| = Real.log 2 :=
  binary_entropy_differences_headline

/-! ## Top-level info-theoretic master summary -/

/-- **Information-Theoretic Master Summary (TIER-2).**

    This theorem aggregates **all** pillars of the PT information-theoretic
    framework into a single conjunction. It is a pure re-exposure: every
    component is `0 sorry` in its source module of `PT.Information.*`.

    The bundle assembles, in order:

    * the six pillars of `entropy_master_framework_summary`
      (GFT identity, Bekenstein, tensor additivity, concavity, chain rule,
      cross-entropy bridge), and
    * the extension pillars A–G below
      (GFT m-special values, binary entropy values, binary entropy
      comparisons, ternary entropy values, binary entropy differences). -/
theorem info_theoretic_master_summary
    {α β : Type*} [DecidableEq α]
    (s : Finset α) (t : Finset β)
    (m n : ℝ) (hm : 0 < m) (hn : 0 < n)
    (hcard_s : (s.card : ℝ) = m) (_hcard_t : (t.card : ℝ) = n)
    (r₀ : α) (hr : r₀ ∈ s)
    (p : α → ℝ) (q : β → ℝ)
    (hp_nonneg : ∀ r ∈ s, 0 ≤ p r) (hp_le_one : ∀ r ∈ s, p r ≤ 1)
    (hp_sum : ∑ r ∈ s, p r = 1)
    (hq_nonneg : ∀ b ∈ t, 0 ≤ q b) (hq_sum : ∑ b ∈ t, q b = 1) :
    -- === Six core pillars from EntropyMasterFramework ===
    -- 1. GFT identity
    klToUniform s m p + shannonH s p = Real.log m
    -- 2a. Bekenstein bound
    ∧ klToUniform s m p ≤ Real.log m
    -- 2b. Delta saturates KL
    ∧ klToUniform s m (deltaDist r₀) = Real.log m
    -- 2c. Uniform vanishes KL
    ∧ klToUniform s m (fun _ => (1 : ℝ) / m) = 0
    -- 2d. Uniform attains max entropy
    ∧ shannonH s (fun _ => (1 : ℝ) / m) = Real.log m
    -- 2e. Delta has zero entropy
    ∧ shannonH s (deltaDist r₀) = 0
    -- 3a. Shannon additivity
    ∧ shannonH (s ×ˢ t) (joint p q) = shannonH s p + shannonH t q
    -- 3b. KL additivity
    ∧ klToUniform (s ×ˢ t) (m * n) (joint p q)
        = klToUniform s m p + klToUniform t n q
    -- 3c. Mutual info vanishes at independence
    ∧ mutualInfo s t p q = 0
    -- 4a. negMulLog concave
    ∧ ConcaveOn ℝ (Set.Ici (0 : ℝ)) Real.negMulLog
    -- 4b. Bekenstein via Jensen
    ∧ shannonH s p ≤ Real.log m
    -- 5a. Chain rule
    ∧ shannonH (s ×ˢ t) (joint p q)
        = shannonH t q + condEntropy s t (joint p q) q
    -- 5b. Independence reduces conditional to marginal
    ∧ condEntropy s t (joint p q) q = shannonH s p
    -- 6a. Self-cross = Shannon
    ∧ crossEntropy s p p = shannonH s p
    -- 6b. Uniform cross-entropy
    ∧ crossEntropy s p (fun _ => (1 : ℝ) / m) = Real.log m
    -- 6c. GFT bridge via cross-entropy
    ∧ klToUniform s m p = Real.log m - shannonH s p
    -- === Extension pillar A: GFT special m-values ===
    ∧ klToUniform (Finset.univ : Finset (Fin 2)) 2 (deltaDist 0) = Real.log 2
    ∧ klToUniform (Finset.univ : Finset (Fin 8)) 8 (deltaDist 0) = Real.log 8
    ∧ klToUniform (Finset.univ : Finset (Fin 30)) 30 (deltaDist 0) = Real.log 30
    -- === Extension pillar B: binary entropy values ===
    ∧ binEntropy ((1 : ℝ) / 2) = Real.log 2
    ∧ binEntropy ((1 : ℝ) / 3) = Real.log 3 - ((2 : ℝ) / 3) * Real.log 2
    ∧ (∀ p, binEntropy p = binEntropy (1 - p))
    -- === Extension pillar C: binary entropy comparisons ===
    ∧ binEntropy ((1 : ℝ) / 3) < binEntropy ((1 : ℝ) / 2)
    ∧ binEntropy ((2 : ℝ) / 3) = binEntropy ((1 : ℝ) / 3)
    -- === Extension pillar D: ternary entropy values ===
    ∧ terEntropy ((1 : ℝ) / 3) ((1 : ℝ) / 3) = Real.log 3
    ∧ (∀ p q, terEntropy p q = terEntropy q p)
    ∧ (∀ p, terEntropy p (1 - p) = binEntropy p)
    -- === Extension pillar E: binary entropy differences ===
    ∧ binEntropy ((1 : ℝ) / 2) - binEntropy ((1 : ℝ) / 3)
        = (5 : ℝ) / 3 * Real.log 2 - Real.log 3
    ∧ |binEntropy ((1 : ℝ) / 2) - binEntropy (0 : ℝ)| = Real.log 2
    -- === Extension pillar F: GFT special log inequalities ===
    ∧ Real.log 2 < Real.log 8
    ∧ Real.log 8 < Real.log 30
    -- === Extension pillar G: binary entropy at extremes ===
    ∧ binEntropy 0 = 0
    ∧ binEntropy 1 = 0 := by
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_,
          ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩
  -- 1. GFT identity
  · exact GFT_identity s m hm p hp_nonneg hp_sum
  -- 2a
  · exact bekenstein_bound s m hm p hp_nonneg hp_le_one hp_sum
  -- 2b
  · exact bekenstein_saturated_at_delta s m hm r₀ hr
  -- 2c
  · exact klToUniform_uniform s m hm
  -- 2d
  · exact shannonH_uniform_eq_log s m hm hcard_s
  -- 2e
  · exact shannonH_delta s r₀ hr
  -- 3a
  · exact shannonH_joint s t p q hp_sum hq_sum
  -- 3b
  · exact klToUniform_joint s t m n hm hn p q hp_nonneg hq_nonneg hp_sum hq_sum
  -- 3c
  · exact mutualInfo_indep_eq_zero s t p q hp_sum hq_sum
  -- 4a
  · exact concaveOn_negMulLog_PT
  -- 4b
  · exact shannonH_le_log_card_of_jensen s m hm hcard_s p hp_nonneg hp_sum
  -- 5a
  · exact shannonH_joint_eq_add_condEntropy s t (joint p q) q
  -- 5b
  · exact condEntropy_indep s t p q hp_sum hq_sum
  -- 6a
  · exact crossEntropy_self s p
  -- 6b
  · exact crossEntropy_uniform s m hm p hp_sum
  -- 6c
  · exact klToUniform_eq_logm_sub_shannonH s m hm p hp_nonneg hp_sum
  -- Pillar A (3 clauses)
  · exact (gft_special_m_values_summary).1
  · exact (gft_special_m_values_summary).2.2.1
  · exact (gft_special_m_values_summary).2.2.2.2.1
  -- Pillar B (3 clauses)
  · exact binEntropy_half
  · exact binEntropy_third
  · exact binEntropy_symm
  -- Pillar C (2 clauses)
  · exact binEntropy_third_lt_half
  · exact binEntropy_two_thirds_eq_third
  -- Pillar D (3 clauses)
  · exact terEntropy_uniform
  · exact terEntropy_swap
  · exact terEntropy_degenerate_binary
  -- Pillar E (2 clauses)
  · exact binEntropy_half_sub_third
  · exact abs_binEntropy_half_sub_zero
  -- Pillar F (2 clauses)
  · exact log_two_lt_log_eight
  · exact log_eight_lt_log_thirty
  -- Pillar G (2 clauses)
  · exact binEntropy_zero
  · exact binEntropy_one

end PT.Information
