/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Holonomy.InverseSinSq
import PT.Holonomy.CouplingReconstructionBounds
import PT.Holonomy.SinSqProductChain
import Mathlib.Tactic

/-!
# Product of inverses `∏ 1 / sin²θ_p` (active cascade + echo extension)

This module computes the **multiplicative** counterpart of the inverse
cyclic-phase squared sines defined in `InverseSinSq`. Whereas the published
PT inverse-coupling formula `α_EM^{-1} ≈ 137` is itself a **product** (not a
sum) over the active primes — namely the reciprocal of the bare coupling
`α_bare = ∏_{p ∈ {3,5,7}} sin²θ_p` — this file makes that identity explicit
and records:

1. The **active inverse product**
   `invProductActive = (1/sin²θ_3) · (1/sin²θ_5) · (1/sin²θ_7)`,
   which is **algebraically equal** to `1 / α_bare = α_bare⁻¹`.
2. The corresponding **bracket** `135 < invProductActive < 137`
   (inherited from `alphaBareInvQ_bracket`) plus a tighter bracket
   `136 < invProductActive < 137` decided by `norm_num`, matching the
   published value `α_bare⁻¹ ≈ 136.278`.
3. The **5-factor extension** `invProductChain` over the active *and* echo
   primes `{3, 5, 7, 11, 13}`, equal to `1 / P_5` (with `P_5` from
   `SinSqProductChain`), bracketed at `7800 < invProductChain < 7810`
   (published value `≈ 7803.27`).

All identities are **exact rational arithmetic** at the PT fixed point
`q_+ = 13/15`. No analytic or measure-theoretic ingredient.

## Reference

* Monograph Chapter 9 §"Reconstruction des couplages", BA5 Step 4
  (inverse-coupling form).
* `app_e_constants.tex`, value `α_bare^{-1} ≈ 136.278`.
-/

namespace PT.Holonomy

/-! ### Active 3-prime inverse product -/

/-- The product over the active primes of the inverse cyclic-phase
    squared sines at `q_+ = 13/15`:

    `invProductActive := (1/sin²θ_3) · (1/sin²θ_5) · (1/sin²θ_7)`.

    By construction this equals `1 / α_bare = α_bare⁻¹`. -/
def invProductActive : ℚ := invSinSqQ 3 * invSinSqQ 5 * invSinSqQ 7

/-- **Algebraic identity (active cascade).**
    `invProductActive = 1 / α_bare`. -/
theorem invProductActive_eq_inv_alphaBare :
    invProductActive = 1 / alphaBareQ := by
  unfold invProductActive invSinSqQ alphaBareQ
  have h3 : sinSqQ 3 ≠ 0 := ne_of_gt sinSq_3_pos
  have h5 : sinSqQ 5 ≠ 0 := ne_of_gt sinSq_5_pos
  have h7 : sinSqQ 7 ≠ 0 := ne_of_gt sinSq_7_pos
  field_simp

/-- **Algebraic identity (active cascade, reciprocal form).**
    `invProductActive = α_bare⁻¹`. -/
theorem invProductActive_eq_alphaBareInv :
    invProductActive = alphaBareInvQ := by
  unfold alphaBareInvQ
  exact invProductActive_eq_inv_alphaBare

/-- The active inverse product is strictly positive. -/
theorem invProductActive_pos : 0 < invProductActive := by
  unfold invProductActive
  exact mul_pos (mul_pos invSinSq_3_pos invSinSq_5_pos) invSinSq_7_pos

/-! ### Brackets (active cascade) -/

/-- **Coarse bracket (inherited from `alphaBareInvQ_bracket`).**
    `135 < invProductActive < 137`. -/
theorem invProductActive_bracket :
    135 < invProductActive ∧ invProductActive < 137 := by
  rw [invProductActive_eq_alphaBareInv]
  exact alphaBareInvQ_bracket

/-- **Tight bracket.** `136 < invProductActive < 137`
    (the exact value is `≈ 136.278`). -/
theorem invProductActive_bracket_tight :
    136 < invProductActive ∧ invProductActive < 137 := by
  unfold invProductActive invSinSqQ sinSqQ deltaQ qPT
  refine ⟨?_, ?_⟩ <;> norm_num

/-- **Three-decimal bracket.** `13627/100 < invProductActive < 13629/100`
    centered on the exact value `≈ 136.278`. -/
theorem invProductActive_bracket_three_decimals :
    13627 / 100 < invProductActive ∧ invProductActive < 13629 / 100 := by
  unfold invProductActive invSinSqQ sinSqQ deltaQ qPT
  refine ⟨?_, ?_⟩ <;> norm_num

/-! ### 5-factor extension over the echo primes -/

/-- The product over the active primes **and** the first two echo primes
    `{3, 5, 7, 11, 13}` of the inverse cyclic-phase squared sines at
    `q_+ = 13/15`:

    `invProductChain :=
       (1/sin²θ_3) · (1/sin²θ_5) · (1/sin²θ_7) · (1/sin²θ_11) · (1/sin²θ_13)`.

    Equals `1 / P_5` where `P_5 = sin²θ_3 · … · sin²θ_13` is from
    `SinSqProductChain`. -/
def invProductChain : ℚ :=
  invSinSqQ 3 * invSinSqQ 5 * invSinSqQ 7 * invSinSqQ 11 * invSinSqQ 13

/-- **Algebraic identity (5-factor chain).**
    `invProductChain = 1 / P_5`. -/
theorem invProductChain_eq_inv_P5 :
    invProductChain = 1 / P5 := by
  unfold invProductChain invSinSqQ
  rw [P5_eq_prod]
  have h3 : sinSqQ 3 ≠ 0 := ne_of_gt sinSq_3_pos
  have h5 : sinSqQ 5 ≠ 0 := ne_of_gt sinSq_5_pos
  have h7 : sinSqQ 7 ≠ 0 := ne_of_gt sinSq_7_pos
  have h11 : sinSqQ 11 ≠ 0 := ne_of_gt sinSq_11_pos
  have h13 : sinSqQ 13 ≠ 0 := ne_of_gt sinSq_13_pos
  field_simp

/-- The 5-factor inverse product is strictly positive. -/
theorem invProductChain_pos : 0 < invProductChain := by
  unfold invProductChain
  exact mul_pos (mul_pos (mul_pos (mul_pos
    invSinSq_3_pos invSinSq_5_pos) invSinSq_7_pos)
    invSinSq_11_pos) invSinSq_13_pos

/-- **Echo extension is strictly larger than the active inverse product**
    (each echo factor `1/sin²θ_p > 1`):
    `invProductActive < invProductChain`. -/
theorem invProductActive_lt_invProductChain :
    invProductActive < invProductChain := by
  -- Strategy: write invProductChain = invProductActive * invSinSqQ 11 * invSinSqQ 13.
  -- Both echo factors are > 1 and invProductActive > 0, so the product strictly grows.
  have hchain :
      invProductChain = invProductActive * invSinSqQ 11 * invSinSqQ 13 := by
    unfold invProductActive invProductChain
    ring
  rw [hchain]
  have hpos : 0 < invProductActive := invProductActive_pos
  have h11 : 1 < invSinSqQ 11 := invSinSq_11_gt_one
  have h13 : 1 < invSinSqQ 13 := invSinSq_13_gt_one
  -- Step 1: invProductActive < invProductActive * invSinSqQ 11
  have step1 : invProductActive < invProductActive * invSinSqQ 11 := by
    have := (lt_mul_iff_one_lt_right hpos).mpr h11
    linarith
  -- Step 2: invProductActive * invSinSqQ 11 < invProductActive * invSinSqQ 11 * invSinSqQ 13
  have h11pos : 0 < invSinSqQ 11 := invSinSq_11_pos
  have hmid_pos : 0 < invProductActive * invSinSqQ 11 := mul_pos hpos h11pos
  have step2 : invProductActive * invSinSqQ 11 <
      invProductActive * invSinSqQ 11 * invSinSqQ 13 := by
    have := (lt_mul_iff_one_lt_right hmid_pos).mpr h13
    linarith
  linarith

/-! ### Brackets (5-factor chain) -/

/-- **Coarse bracket.** `7000 < invProductChain < 8000`
    (the exact value is `≈ 7803.27`). -/
theorem invProductChain_bracket :
    7000 < invProductChain ∧ invProductChain < 8000 := by
  unfold invProductChain invSinSqQ sinSqQ deltaQ qPT
  refine ⟨?_, ?_⟩ <;> norm_num

/-- **Tight bracket.** `7800 < invProductChain < 7810`
    (the exact value is `≈ 7803.27`). -/
theorem invProductChain_bracket_tight :
    7800 < invProductChain ∧ invProductChain < 7810 := by
  unfold invProductChain invSinSqQ sinSqQ deltaQ qPT
  refine ⟨?_, ?_⟩ <;> norm_num

/-! ### Headline: complete inverse-product summary -/

/-- **Headline (inverse products, active cascade + echo extension).**
    All quantities at the PT fixed point `q_+ = 13/15`:

    * `invProductActive := ∏_{p ∈ {3,5,7}} 1/sin²θ_p`,
      bracketed `136 < invProductActive < 137` (exact `≈ 136.278`),
      and equal to `α_bare⁻¹`.
    * `invProductChain := ∏_{p ∈ {3,5,7,11,13}} 1/sin²θ_p`,
      bracketed `7800 < invProductChain < 7810` (exact `≈ 7803.27`),
      and equal to `1 / P_5`.
    * Strict echo enlargement:
      `invProductActive < invProductChain` (the echo factors
      `1/sin²θ_11`, `1/sin²θ_13` are both `> 1`). -/
theorem invProduct_summary :
    -- algebraic identities
    invProductActive = alphaBareInvQ
    ∧ invProductChain = 1 / P5
    -- positivity
    ∧ 0 < invProductActive ∧ 0 < invProductChain
    -- coarse brackets
    ∧ 135 < invProductActive ∧ invProductActive < 137
    ∧ 7000 < invProductChain ∧ invProductChain < 8000
    -- tight brackets
    ∧ 136 < invProductActive ∧ invProductActive < 137
    ∧ 7800 < invProductChain ∧ invProductChain < 7810
    -- echo enlargement
    ∧ invProductActive < invProductChain :=
  ⟨invProductActive_eq_alphaBareInv,
   invProductChain_eq_inv_P5,
   invProductActive_pos, invProductChain_pos,
   invProductActive_bracket.1, invProductActive_bracket.2,
   invProductChain_bracket.1, invProductChain_bracket.2,
   invProductActive_bracket_tight.1, invProductActive_bracket_tight.2,
   invProductChain_bracket_tight.1, invProductChain_bracket_tight.2,
   invProductActive_lt_invProductChain⟩

end PT.Holonomy
