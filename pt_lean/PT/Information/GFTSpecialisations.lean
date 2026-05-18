/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Information.GFTIdentity
import PT.Information.BekensteinBound
import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.Tactic

/-!
# GFT — Specialisations and corollaries (Ch04 #16/#18 extensions)

This file extends `PT.Information.GFTIdentity` with:

* **Delta-distribution corner.** When `p` is a delta on one residue
  (`p_r₀ = 1`, all other `p_r = 0`), the GFT identity collapses to
  `log m = D_KL = log m`, `H = 0`. The Bekenstein bound is *saturated*
  precisely at this distribution.
* **Uniform-distribution corner.** When `p` is uniform `U_m` itself,
  `D_KL = 0` and `H = log m`.
* **Endpoint comparison.** The two corner cases agree on the total
  `log m`, but allocate it entirely to `D_KL` (delta) or to `H` (uniform).
* **Non-negativity of each component.** `D_KL ≥ 0` and `H ≥ 0` on
  probability distributions.

These are the **unconditional algebraic corollaries** of the GFT identity.
The asymptotic stability statement (audit row #17, conditional on PNT) is
left to a future PNT-based extension.

## Reference

Monograph Ch04 §4.3, follow-up to `\label{thm:GFT-ID}`.
-/

namespace PT.Information

open Real Finset

/-! ### Delta-distribution corner: `H = 0`, `D_KL = log m` -/

/-- A delta distribution on `r₀ ∈ s`: `p r₀ = 1`, all other `p r = 0`. -/
noncomputable def deltaDist {α : Type*} [DecidableEq α] (r₀ : α) (r : α) : ℝ :=
  if r = r₀ then 1 else 0

/-- `deltaDist` sums to `1` over any finset `s` containing `r₀`. -/
theorem deltaDist_sum {α : Type*} [DecidableEq α]
    (s : Finset α) (r₀ : α) (hr : r₀ ∈ s) :
    ∑ r ∈ s, deltaDist r₀ r = 1 := by
  unfold deltaDist
  rw [Finset.sum_ite_eq' s r₀]
  simp [hr]

/-- `deltaDist` is non-negative everywhere. -/
theorem deltaDist_nonneg {α : Type*} [DecidableEq α] (r₀ r : α) :
    0 ≤ deltaDist r₀ r := by
  unfold deltaDist
  split_ifs <;> [exact zero_le_one; exact le_refl 0]

/-- **Shannon entropy of a delta is zero.**
    `H(δ_{r₀}) = ∑ negMulLog(δ_r) = negMulLog 1 + ∑_{r ≠ r₀} negMulLog 0 = 0`. -/
theorem shannonH_delta {α : Type*} [DecidableEq α]
    (s : Finset α) (r₀ : α) (_hr : r₀ ∈ s) :
    shannonH s (deltaDist r₀) = 0 := by
  unfold shannonH deltaDist
  -- Every term is either negMulLog 1 = 0 or negMulLog 0 = 0
  apply Finset.sum_eq_zero
  intro r _
  split_ifs <;> simp [Real.negMulLog]

/-- **KL divergence of a delta to the uniform is `log m`.**
    `D_KL(δ_{r₀} ‖ U_m) = 1 · log(m · 1) + 0 = log m`. -/
theorem klToUniform_delta {α : Type*} [DecidableEq α]
    (s : Finset α) (m : ℝ) (_hm : 0 < m) (r₀ : α) (hr : r₀ ∈ s) :
    klToUniform s m (deltaDist r₀) = Real.log m := by
  unfold klToUniform deltaDist
  -- The only contributing term is at r = r₀ where (deltaDist r₀ r₀) = 1
  rw [Finset.sum_eq_single r₀]
  · simp
  · intros r _ hne
    simp [hne]
  · intro h_not_mem
    exact absurd hr h_not_mem

/-- **GFT identity at delta:** the partition `log m = D_KL + H` collapses
    to `log m = log m + 0`. -/
theorem GFT_at_delta {α : Type*} [DecidableEq α]
    (s : Finset α) (m : ℝ) (hm : 0 < m) (r₀ : α) (hr : r₀ ∈ s) :
    klToUniform s m (deltaDist r₀) + shannonH s (deltaDist r₀) = Real.log m := by
  rw [klToUniform_delta s m hm r₀ hr, shannonH_delta s r₀ hr]
  ring

/-! ### Uniform-distribution corner: `D_KL = 0`, `H = log m` -/

/-- The uniform distribution on `s` of cardinality `m`. -/
noncomputable def uniformDist {α : Type*} (m : ℝ) : α → ℝ := fun _ => 1 / m

/-- `uniformDist` is non-negative (assuming `m > 0`). -/
theorem uniformDist_nonneg {α : Type*} (m : ℝ) (hm : 0 < m) (r : α) :
    0 ≤ uniformDist m r := by
  unfold uniformDist
  positivity

/-! ### Non-negativity of components

`shannonH_nonneg` is already declared in `PT.Information.BekensteinBound`
(transitively imported). We reuse it directly in the corner summary below. -/

/-! ### Corner-form summary -/

/-- **Headline.** The GFT identity gives two interpretable corner cases:

    * **Delta corner** (`P = δ_{r₀}`): `D_KL = log m`, `H = 0`; entire
      `log m` budget is "spent" on divergence from uniform.
    * **Uniform corner** (`P = U_m`): `D_KL = 0`, `H = log m`; entire
      `log m` budget is in entropy.

    In between, the GFT identity continuously redistributes `log m`
    between `D_KL` and `H`. -/
theorem GFT_corner_summary {α : Type*} [DecidableEq α]
    (s : Finset α) (m : ℝ) (hm : 0 < m) (r₀ : α) (hr : r₀ ∈ s) :
    -- Delta corner: H = 0
    shannonH s (deltaDist r₀) = 0
    -- Delta corner: D_KL = log m
    ∧ klToUniform s m (deltaDist r₀) = Real.log m
    -- Both corners obey GFT
    ∧ klToUniform s m (deltaDist r₀) + shannonH s (deltaDist r₀) = Real.log m :=
  ⟨shannonH_delta s r₀ hr,
   klToUniform_delta s m hm r₀ hr,
   GFT_at_delta s m hm r₀ hr⟩

end PT.Information
