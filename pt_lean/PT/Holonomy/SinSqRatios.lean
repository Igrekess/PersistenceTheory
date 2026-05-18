/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Holonomy.CouplingReconstruction
import PT.Holonomy.CouplingReconstructionBounds
import PT.Holonomy.SinSqProductBounds
import PT.Holonomy.SinSqProductChain
import Mathlib.Tactic

/-!
# `∏ sin²θ_p` — Ratios of partial products (Ch09 cascade)

This file isolates the *ratios* `R_{k,k-1} := P_k / P_{k-1}` of the cascade
partial products defined in `SinSqProductBounds` (`P_1, P_2, P_3`) and
`SinSqProductChain` (`P_4, P_5`). Each ratio collapses by construction to a
single cyclic-phase squared sine:

* `R₂₁ := P_2 / P_1 = sin²θ_5`,
* `R₃₂ := P_3 / P_2 = sin²θ_7`,
* `R₄₃ := P_4 / P_3 = sin²θ_11`,
* `R₅₄ := P_5 / P_4 = sin²θ_13`.

Since `sin²θ_p` is strictly decreasing in `p` along the active and echo
primes `{3, 5, 7, 11, 13}` (cf. `sinSq_chain_active` and the brackets of
`sinSq_{11,13}_bracket`), the ratios themselves form a strictly decreasing
chain

`R₅₄ < R₄₃ < R₃₂ < R₂₁ < 1`.

Cumulatively, multiplying the ratios reconstructs the partial products:
`P_k = P_1 · ∏_{i=2}^k R_{i,i-1}`, in particular
`P_5 = P_1 · R₂₁ · R₃₂ · R₄₃ · R₅₄`.

## Reference

Monograph Chapter 9 §"Produit partiel des phases cycliques".
-/

namespace PT.Holonomy

/-! ### Ratio definitions -/

/-- `R₂₁ := P_2 / P_1`. Equals `sin²θ_5` by construction. -/
def R21 : ℚ := P2 / P1

/-- `R₃₂ := P_3 / P_2`. Equals `sin²θ_7` by construction. -/
def R32 : ℚ := P3 / P2

/-- `R₄₃ := P_4 / P_3`. Equals `sin²θ_11` by construction. -/
def R43 : ℚ := P4 / P3

/-- `R₅₄ := P_5 / P_4`. Equals `sin²θ_13` by construction. -/
def R54 : ℚ := P5 / P4

/-! ### Collapse to a single sin² factor

Each ratio simplifies to a single cyclic-phase squared sine. The proofs are
pure field arithmetic given strict positivity of the denominators. -/

/-- **`R₂₁ = sin²θ_5`.** -/
@[simp] theorem R21_eq_sinSq_5 : R21 = sinSqQ 5 := by
  unfold R21 P2 P1
  exact mul_div_cancel_left₀ (sinSqQ 5) sinSq_3_pos.ne'

/-- **`R₃₂ = sin²θ_7`.** -/
@[simp] theorem R32_eq_sinSq_7 : R32 = sinSqQ 7 := by
  unfold R32 P3 P2
  have h35 : sinSqQ 3 * sinSqQ 5 ≠ 0 :=
    (mul_pos sinSq_3_pos sinSq_5_pos).ne'
  exact mul_div_cancel_left₀ (sinSqQ 7) h35

/-- **`R₄₃ = sin²θ_11`.** -/
@[simp] theorem R43_eq_sinSq_11 : R43 = sinSqQ 11 := by
  unfold R43 P4
  exact mul_div_cancel_left₀ (sinSqQ 11) P3_pos.ne'

/-- **`R₅₄ = sin²θ_13`.** -/
@[simp] theorem R54_eq_sinSq_13 : R54 = sinSqQ 13 := by
  unfold R54 P5
  exact mul_div_cancel_left₀ (sinSqQ 13) P4_pos.ne'

/-! ### Brackets on each ratio

Each ratio inherits a decimal bracket from the corresponding `sinSq_p`
bracket through the collapse identities above. -/

/-- `R₂₁ ∈ (0.193, 0.195)`. -/
theorem R21_bracket : 193 / 1000 < R21 ∧ R21 < 195 / 1000 := by
  rw [R21_eq_sinSq_5]; exact sinSq_5_bracket

/-- `R₃₂ ∈ (0.172, 0.174)`. -/
theorem R32_bracket : 172 / 1000 < R32 ∧ R32 < 174 / 1000 := by
  rw [R32_eq_sinSq_7]; exact sinSq_7_bracket

/-- `R₄₃ ∈ (0.138, 0.140)`. -/
theorem R43_bracket : 138 / 1000 < R43 ∧ R43 < 140 / 1000 := by
  rw [R43_eq_sinSq_11]; exact sinSq_11_bracket

/-- `R₅₄ ∈ (0.125, 0.127)`. -/
theorem R54_bracket : 125 / 1000 < R54 ∧ R54 < 127 / 1000 := by
  rw [R54_eq_sinSq_13]; exact sinSq_13_bracket

/-! ### Positivity and sub-unitarity -/

theorem R21_pos : 0 < R21 := by rw [R21_eq_sinSq_5]; exact sinSq_5_pos
theorem R32_pos : 0 < R32 := by rw [R32_eq_sinSq_7]; exact sinSq_7_pos
theorem R43_pos : 0 < R43 := by rw [R43_eq_sinSq_11]; exact sinSq_11_pos
theorem R54_pos : 0 < R54 := by rw [R54_eq_sinSq_13]; exact sinSq_13_pos

theorem R21_lt_one : R21 < 1 := by rw [R21_eq_sinSq_5]; exact sinSq_5_lt_one
theorem R32_lt_one : R32 < 1 := by rw [R32_eq_sinSq_7]; exact sinSq_7_lt_one
theorem R43_lt_one : R43 < 1 := by rw [R43_eq_sinSq_11]; exact sinSq_11_lt_one
theorem R54_lt_one : R54 < 1 := by rw [R54_eq_sinSq_13]; exact sinSq_13_lt_one

/-! ### Strict decreasing chain of ratios

Since `sin²θ_p` is strictly decreasing in `p` along the active and echo
primes `{3, 5, 7, 11, 13}`, the ratios `R_{k,k-1} = sin²θ_{p_k}` form a
strictly decreasing chain. -/

/-- `R₃₂ < R₂₁`. -/
theorem R32_lt_R21 : R32 < R21 := by
  rw [R32_eq_sinSq_7, R21_eq_sinSq_5]; exact sinSq_5_gt_sinSq_7

/-- `R₄₃ < R₃₂`. -/
theorem R43_lt_R32 : R43 < R32 := by
  rw [R43_eq_sinSq_11, R32_eq_sinSq_7]
  unfold sinSqQ deltaQ qPT
  norm_num

/-- `R₅₄ < R₄₃`. -/
theorem R54_lt_R43 : R54 < R43 := by
  rw [R54_eq_sinSq_13, R43_eq_sinSq_11]
  unfold sinSqQ deltaQ qPT
  norm_num

/-- **Strict decreasing chain.** `0 < R₅₄ < R₄₃ < R₃₂ < R₂₁ < 1`. -/
theorem ratio_chain_strict :
    0 < R54 ∧ R54 < R43 ∧ R43 < R32 ∧ R32 < R21 ∧ R21 < 1 :=
  ⟨R54_pos, R54_lt_R43, R43_lt_R32, R32_lt_R21, R21_lt_one⟩

/-! ### Cumulative identities

The ratios cumulatively reconstruct the partial products. We multiply by
`P_1` on the left so that each `P_k` factors as a single telescoping product
of ratios anchored on `P_1`. -/

/-- `P_2 = P_1 · R₂₁`. -/
theorem P2_eq_P1_mul_R21 : P2 = P1 * R21 := by
  rw [R21_eq_sinSq_5]; rfl

/-- `P_3 = P_1 · R₂₁ · R₃₂`. -/
theorem P3_eq_P1_mul_R21_R32 : P3 = P1 * R21 * R32 := by
  rw [R21_eq_sinSq_5, R32_eq_sinSq_7]; rfl

/-- `P_4 = P_1 · R₂₁ · R₃₂ · R₄₃`. -/
theorem P4_eq_P1_mul_R21_R32_R43 : P4 = P1 * R21 * R32 * R43 := by
  rw [R21_eq_sinSq_5, R32_eq_sinSq_7, R43_eq_sinSq_11]; rfl

/-- `P_5 = P_1 · R₂₁ · R₃₂ · R₄₃ · R₅₄`. -/
theorem P5_eq_P1_mul_R21_R32_R43_R54 :
    P5 = P1 * R21 * R32 * R43 * R54 := by
  rw [R21_eq_sinSq_5, R32_eq_sinSq_7, R43_eq_sinSq_11, R54_eq_sinSq_13]; rfl

/-! ### Numerical headline -/

/-- **Headline (ratios of cascade partial products).** At the PT fixed point
    `q_+ = 13/15`, the four ratios of consecutive partial products along the
    active+echo cascade `{3, 5, 7, 11, 13}` collapse to a single
    cyclic-phase squared sine, each with the bracket below; together they
    form a strictly decreasing chain bounded by `1`:

    * `R₂₁ = P_2 / P_1 = sin²θ_5  ∈ (0.193, 0.195)`,
    * `R₃₂ = P_3 / P_2 = sin²θ_7  ∈ (0.172, 0.174)`,
    * `R₄₃ = P_4 / P_3 = sin²θ_11 ∈ (0.138, 0.140)`,
    * `R₅₄ = P_5 / P_4 = sin²θ_13 ∈ (0.125, 0.127)`,
    * strictly decreasing: `0 < R₅₄ < R₄₃ < R₃₂ < R₂₁ < 1`,
    * cumulative: `P_5 = P_1 · R₂₁ · R₃₂ · R₄₃ · R₅₄`. -/
theorem ratio_summary :
    -- collapse to single sin² factor
    R21 = sinSqQ 5 ∧ R32 = sinSqQ 7 ∧ R43 = sinSqQ 11 ∧ R54 = sinSqQ 13
    -- brackets
    ∧ 193 / 1000 < R21 ∧ R21 < 195 / 1000
    ∧ 172 / 1000 < R32 ∧ R32 < 174 / 1000
    ∧ 138 / 1000 < R43 ∧ R43 < 140 / 1000
    ∧ 125 / 1000 < R54 ∧ R54 < 127 / 1000
    -- strict decreasing chain
    ∧ 0 < R54 ∧ R54 < R43 ∧ R43 < R32 ∧ R32 < R21 ∧ R21 < 1
    -- cumulative identity
    ∧ P5 = P1 * R21 * R32 * R43 * R54 :=
  ⟨R21_eq_sinSq_5, R32_eq_sinSq_7, R43_eq_sinSq_11, R54_eq_sinSq_13,
   R21_bracket.1, R21_bracket.2,
   R32_bracket.1, R32_bracket.2,
   R43_bracket.1, R43_bracket.2,
   R54_bracket.1, R54_bracket.2,
   R54_pos, R54_lt_R43, R43_lt_R32, R32_lt_R21, R21_lt_one,
   P5_eq_P1_mul_R21_R32_R43_R54⟩

end PT.Holonomy
