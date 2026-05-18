/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Holonomy.CouplingReconstruction
import PT.Holonomy.SinSqProductChain
import Mathlib.Tactic

/-!
# Inverse cyclic-phase squared sines `1 / sin²θ_p`

The fine-structure inverse coupling `α_EM^{-1} ≈ 137` is reconstructed from
**sums of inverse squared sines** over the active and echo primes
(monograph Chapter 9 §"Reconstruction des couplages", BA5 Step 4).
This module records the elementary properties of the rational quantities
`1 / sin²θ_p` at the PT fixed point `q_+ = 13/15`:

* positivity (every factor `> 0`),
* sub-unitarity inversion (`sin²θ_p < 1` implies `1 / sin²θ_p > 1`),
* exact rational brackets to three decimals,
* the strict increasing cascade
  `1/sin²θ_3 < 1/sin²θ_5 < 1/sin²θ_7 < 1/sin²θ_11 < 1/sin²θ_13`,
* numerical headline matching the published values

  | `p` | `1/sin²θ_p` | bracket           |
  |-----|-------------|-------------------|
  | 3   | `≈ 4.5630`  | `(4.55, 4.57)`    |
  | 5   | `≈ 5.1553`  | `(5.15, 5.16)`    |
  | 7   | `≈ 5.7933`  | `(5.78, 5.80)`    |
  | 11  | `≈ 7.1967`  | `(7.19, 7.21)`    |
  | 13  | `≈ 7.9564`  | `(7.95, 7.97)`    |

These inverses are the rational summands used downstream by the bracket
`α_EM^{-1} ∈ (135, 137)` of `CouplingReconstructionBounds`.

## Reference

* Monograph Chapter 9, §"Reconstruction des couplages".
* `PT-WEBSITE/.../theorems/en/BA5.mdx`, Step 4 — inverse-coupling form.
-/

namespace PT.Holonomy

/-! ### Definition -/

/-- The inverse cyclic-phase squared sine at the PT branch parameter
    `q_+ = 13/15`:

    `invSinSqQ p := 1 / sin²θ_p`, where `sin²θ_p` is `sinSqQ p` defined in
    `CouplingReconstruction`. -/
def invSinSqQ (p : ℕ) : ℚ := 1 / sinSqQ p

/-! ### Positivity (`> 0`) for active and echo primes -/

/-- `1 / sin²θ_3 > 0`. -/
theorem invSinSq_3_pos : 0 < invSinSqQ 3 := by
  unfold invSinSqQ
  exact one_div_pos.mpr sinSq_3_pos

/-- `1 / sin²θ_5 > 0`. -/
theorem invSinSq_5_pos : 0 < invSinSqQ 5 := by
  unfold invSinSqQ
  exact one_div_pos.mpr sinSq_5_pos

/-- `1 / sin²θ_7 > 0`. -/
theorem invSinSq_7_pos : 0 < invSinSqQ 7 := by
  unfold invSinSqQ
  exact one_div_pos.mpr sinSq_7_pos

/-- `1 / sin²θ_11 > 0`. -/
theorem invSinSq_11_pos : 0 < invSinSqQ 11 := by
  unfold invSinSqQ
  exact one_div_pos.mpr sinSq_11_pos

/-- `1 / sin²θ_13 > 0`. -/
theorem invSinSq_13_pos : 0 < invSinSqQ 13 := by
  unfold invSinSqQ
  exact one_div_pos.mpr sinSq_13_pos

/-! ### Greater than one (`> 1`), since `sin²θ_p < 1` -/

/-- `1 / sin²θ_3 > 1`. -/
theorem invSinSq_3_gt_one : 1 < invSinSqQ 3 := by
  unfold invSinSqQ
  exact one_lt_one_div sinSq_3_pos sinSq_3_lt_one

/-- `1 / sin²θ_5 > 1`. -/
theorem invSinSq_5_gt_one : 1 < invSinSqQ 5 := by
  unfold invSinSqQ
  exact one_lt_one_div sinSq_5_pos sinSq_5_lt_one

/-- `1 / sin²θ_7 > 1`. -/
theorem invSinSq_7_gt_one : 1 < invSinSqQ 7 := by
  unfold invSinSqQ
  exact one_lt_one_div sinSq_7_pos sinSq_7_lt_one

/-- `1 / sin²θ_11 > 1`. -/
theorem invSinSq_11_gt_one : 1 < invSinSqQ 11 := by
  unfold invSinSqQ
  exact one_lt_one_div sinSq_11_pos sinSq_11_lt_one

/-- `1 / sin²θ_13 > 1`. -/
theorem invSinSq_13_gt_one : 1 < invSinSqQ 13 := by
  unfold invSinSqQ
  exact one_lt_one_div sinSq_13_pos sinSq_13_lt_one

/-! ### Decimal brackets (three-decimal precision) -/

/-- `1/sin²θ_3 ∈ (4.55, 4.57)` (decimal value `≈ 4.5630`). -/
theorem invSinSq_3_bracket :
    455 / 100 < invSinSqQ 3 ∧ invSinSqQ 3 < 457 / 100 := by
  unfold invSinSqQ sinSqQ deltaQ qPT
  refine ⟨?_, ?_⟩ <;> norm_num

/-- `1/sin²θ_5 ∈ (5.15, 5.16)` (decimal value `≈ 5.1553`). -/
theorem invSinSq_5_bracket :
    515 / 100 < invSinSqQ 5 ∧ invSinSqQ 5 < 516 / 100 := by
  unfold invSinSqQ sinSqQ deltaQ qPT
  refine ⟨?_, ?_⟩ <;> norm_num

/-- `1/sin²θ_7 ∈ (5.78, 5.80)` (decimal value `≈ 5.7933`). -/
theorem invSinSq_7_bracket :
    578 / 100 < invSinSqQ 7 ∧ invSinSqQ 7 < 580 / 100 := by
  unfold invSinSqQ sinSqQ deltaQ qPT
  refine ⟨?_, ?_⟩ <;> norm_num

/-- `1/sin²θ_11 ∈ (7.19, 7.21)` (decimal value `≈ 7.1967`). -/
theorem invSinSq_11_bracket :
    719 / 100 < invSinSqQ 11 ∧ invSinSqQ 11 < 721 / 100 := by
  unfold invSinSqQ sinSqQ deltaQ qPT
  refine ⟨?_, ?_⟩ <;> norm_num

/-- `1/sin²θ_13 ∈ (7.95, 7.97)` (decimal value `≈ 7.9564`). -/
theorem invSinSq_13_bracket :
    795 / 100 < invSinSqQ 13 ∧ invSinSqQ 13 < 797 / 100 := by
  unfold invSinSqQ sinSqQ deltaQ qPT
  refine ⟨?_, ?_⟩ <;> norm_num

/-! ### Strict increasing cascade

The inverse cascade reverses the monotonicity of `sinSqQ` (cf.
`sinSq_chain_active` and `SinSqProductChain.partial_product_chain_extended`),
since `x ↦ 1/x` is strictly decreasing on `(0, ∞)`. -/

/-- `1/sin²θ_3 < 1/sin²θ_5`. -/
theorem invSinSq_3_lt_invSinSq_5 : invSinSqQ 3 < invSinSqQ 5 := by
  have h1 : invSinSqQ 3 < 457 / 100 := invSinSq_3_bracket.2
  have h2 : (457 / 100 : ℚ) < 515 / 100 := by norm_num
  have h3 : (515 / 100 : ℚ) < invSinSqQ 5 := invSinSq_5_bracket.1
  linarith

/-- `1/sin²θ_5 < 1/sin²θ_7`. -/
theorem invSinSq_5_lt_invSinSq_7 : invSinSqQ 5 < invSinSqQ 7 := by
  have h1 : invSinSqQ 5 < 516 / 100 := invSinSq_5_bracket.2
  have h2 : (516 / 100 : ℚ) < 578 / 100 := by norm_num
  have h3 : (578 / 100 : ℚ) < invSinSqQ 7 := invSinSq_7_bracket.1
  linarith

/-- `1/sin²θ_7 < 1/sin²θ_11`. -/
theorem invSinSq_7_lt_invSinSq_11 : invSinSqQ 7 < invSinSqQ 11 := by
  have h1 : invSinSqQ 7 < 580 / 100 := invSinSq_7_bracket.2
  have h2 : (580 / 100 : ℚ) < 719 / 100 := by norm_num
  have h3 : (719 / 100 : ℚ) < invSinSqQ 11 := invSinSq_11_bracket.1
  linarith

/-- `1/sin²θ_11 < 1/sin²θ_13`. -/
theorem invSinSq_11_lt_invSinSq_13 : invSinSqQ 11 < invSinSqQ 13 := by
  have h1 : invSinSqQ 11 < 721 / 100 := invSinSq_11_bracket.2
  have h2 : (721 / 100 : ℚ) < 795 / 100 := by norm_num
  have h3 : (795 / 100 : ℚ) < invSinSqQ 13 := invSinSq_13_bracket.1
  linarith

/-- **Full strict increasing cascade** for the inverse squared sines over
    the active and echo primes:
    `1/sin²θ_3 < 1/sin²θ_5 < 1/sin²θ_7 < 1/sin²θ_11 < 1/sin²θ_13`. -/
theorem invSinSq_chain_strict :
    invSinSqQ 3 < invSinSqQ 5 ∧ invSinSqQ 5 < invSinSqQ 7
    ∧ invSinSqQ 7 < invSinSqQ 11 ∧ invSinSqQ 11 < invSinSqQ 13 :=
  ⟨invSinSq_3_lt_invSinSq_5, invSinSq_5_lt_invSinSq_7,
   invSinSq_7_lt_invSinSq_11, invSinSq_11_lt_invSinSq_13⟩

/-! ### Numerical headline -/

/-- **Headline (inverse cyclic-phase squared sines, active + echo cascade).**
    All quantities at the PT fixed point `q_+ = 13/15`:

    * `1/sin²θ_3 ∈ (4.55, 4.57)`,
    * `1/sin²θ_5 ∈ (5.15, 5.16)`,
    * `1/sin²θ_7 ∈ (5.78, 5.80)`,
    * `1/sin²θ_11 ∈ (7.19, 7.21)`,
    * `1/sin²θ_13 ∈ (7.95, 7.97)`,
    * strictly increasing:
      `1/sin²θ_3 < 1/sin²θ_5 < 1/sin²θ_7 < 1/sin²θ_11 < 1/sin²θ_13`,
    * each inverse exceeds `1` (since `sin²θ_p < 1`). -/
theorem invSinSq_summary :
    -- decimal brackets
    455 / 100 < invSinSqQ 3 ∧ invSinSqQ 3 < 457 / 100
    ∧ 515 / 100 < invSinSqQ 5 ∧ invSinSqQ 5 < 516 / 100
    ∧ 578 / 100 < invSinSqQ 7 ∧ invSinSqQ 7 < 580 / 100
    ∧ 719 / 100 < invSinSqQ 11 ∧ invSinSqQ 11 < 721 / 100
    ∧ 795 / 100 < invSinSqQ 13 ∧ invSinSqQ 13 < 797 / 100
    -- strict increasing cascade
    ∧ invSinSqQ 3 < invSinSqQ 5 ∧ invSinSqQ 5 < invSinSqQ 7
    ∧ invSinSqQ 7 < invSinSqQ 11 ∧ invSinSqQ 11 < invSinSqQ 13
    -- each inverse exceeds 1
    ∧ 1 < invSinSqQ 3 ∧ 1 < invSinSqQ 5 ∧ 1 < invSinSqQ 7
    ∧ 1 < invSinSqQ 11 ∧ 1 < invSinSqQ 13 :=
  ⟨invSinSq_3_bracket.1, invSinSq_3_bracket.2,
   invSinSq_5_bracket.1, invSinSq_5_bracket.2,
   invSinSq_7_bracket.1, invSinSq_7_bracket.2,
   invSinSq_11_bracket.1, invSinSq_11_bracket.2,
   invSinSq_13_bracket.1, invSinSq_13_bracket.2,
   invSinSq_3_lt_invSinSq_5, invSinSq_5_lt_invSinSq_7,
   invSinSq_7_lt_invSinSq_11, invSinSq_11_lt_invSinSq_13,
   invSinSq_3_gt_one, invSinSq_5_gt_one, invSinSq_7_gt_one,
   invSinSq_11_gt_one, invSinSq_13_gt_one⟩

end PT.Holonomy
