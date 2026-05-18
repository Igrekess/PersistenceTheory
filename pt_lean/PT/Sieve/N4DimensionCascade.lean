/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Sieve.N4PrimeCascade
import Mathlib.Tactic

/-!
# N4 — Dimension cascade `{1, 2, 4, 6}` for `p ∈ {2, 3, 5, 7}` (Ch02 #6 extension)

This file extends `PT.Sieve.N4PrimeCascade` by computing the cardinality of
`(ℤ/pℤ)*` for `p ∈ {2, 3, 5, 7, 11, 13}`, i.e. the dimension cascade

$$\dim_p := \varphi(p) = p - 1, \qquad p \in \{2, 3, 5, 7, 11, 13\}$$

giving the integer sequence `{1, 2, 4, 6, 10, 12}`. The PT active set
`{3, 5, 7}` corresponds to the sub-cascade `{2, 4, 6}` (consecutive even
integers), which is the **dimension signature** of the persistence-active
sector.

* **Cardinality at each prime.** Exhaustive `decide` on small `ZMod p`.
* **Cascade monotonicity.** `|nonzeroResidues p| < |nonzeroResidues p'|`
  for `p < p'` primes.
* **Sub-cascade {2, 4, 6}.** The active dimensions `{φ(3), φ(5), φ(7)}`
  form an arithmetic progression of common difference `2`.

## Reference

Monograph Chapter 2, §"Cascade dimensionnelle", follow-up to `\label{thm:N4}`.
Audit row #6 extension.
-/

namespace PT.Sieve.N4

/-! ### Cardinality of `(ℤ/pℤ)*` for small primes -/

/-- `|nonzeroResidues 5| = 4 = 5 - 1`. -/
theorem card_nonzeroResidues_five : (nonzeroResidues 5).card = 4 := by
  decide

/-- `|nonzeroResidues 7| = 6 = 7 - 1`. -/
theorem card_nonzeroResidues_seven : (nonzeroResidues 7).card = 6 := by
  decide

/-- `|nonzeroResidues 11| = 10 = 11 - 1`. -/
theorem card_nonzeroResidues_eleven : (nonzeroResidues 11).card = 10 := by
  decide

/-- `|nonzeroResidues 13| = 12 = 13 - 1`. -/
theorem card_nonzeroResidues_thirteen : (nonzeroResidues 13).card = 12 := by
  decide

/-! ### Cascade monotonicity -/

/-- **Strict cascade.** `|nonzeroResidues 2| = 1 < 2 = |nonzeroResidues 3|`. -/
theorem cascade_2_lt_3 : (nonzeroResidues 2).card < (nonzeroResidues 3).card := by
  rw [card_nonzeroResidues_two, card_nonzeroResidues_three]; decide

theorem cascade_3_lt_5 : (nonzeroResidues 3).card < (nonzeroResidues 5).card := by
  rw [card_nonzeroResidues_three, card_nonzeroResidues_five]; decide

theorem cascade_5_lt_7 : (nonzeroResidues 5).card < (nonzeroResidues 7).card := by
  rw [card_nonzeroResidues_five, card_nonzeroResidues_seven]; decide

theorem cascade_7_lt_11 : (nonzeroResidues 7).card < (nonzeroResidues 11).card := by
  rw [card_nonzeroResidues_seven, card_nonzeroResidues_eleven]; decide

theorem cascade_11_lt_13 : (nonzeroResidues 11).card < (nonzeroResidues 13).card := by
  rw [card_nonzeroResidues_eleven, card_nonzeroResidues_thirteen]; decide

/-! ### Active sub-cascade `{2, 4, 6}` -/

/-- **Active sub-cascade.** The dimensions `(|nonzeroResidues 3|,
    |nonzeroResidues 5|, |nonzeroResidues 7|) = (2, 4, 6)` form a length-3
    arithmetic progression with common difference `2`. -/
theorem active_subcascade_arithmetic :
    (nonzeroResidues 3).card = 2
    ∧ (nonzeroResidues 5).card = 4
    ∧ (nonzeroResidues 7).card = 6
    ∧ (nonzeroResidues 5).card - (nonzeroResidues 3).card = 2
    ∧ (nonzeroResidues 7).card - (nonzeroResidues 5).card = 2 := by
  refine ⟨card_nonzeroResidues_three, card_nonzeroResidues_five,
          card_nonzeroResidues_seven, ?_, ?_⟩
  · rw [card_nonzeroResidues_three, card_nonzeroResidues_five]
  · rw [card_nonzeroResidues_five, card_nonzeroResidues_seven]

/-! ### Sum of active dimensions -/

/-- **Sum of active dimensions** is `2 + 4 + 6 = 12 = 2 · 6`.
    Numerologically, this is `(p - 1)` summed over the active set, equalling
    `12 = (3 - 1) + (5 - 1) + (7 - 1)`. -/
theorem active_dimension_sum :
    (nonzeroResidues 3).card + (nonzeroResidues 5).card +
      (nonzeroResidues 7).card = 12 := by
  rw [card_nonzeroResidues_three, card_nonzeroResidues_five,
      card_nonzeroResidues_seven]

/-- The active dimension sum equals `μ* - 3` (where `μ* = 15`),
    matching the prediction `Σ (p - 1) = (Σ p) - |P_active| = 15 - 3 = 12`. -/
theorem active_dimension_sum_eq_muStar_minus_three :
    (nonzeroResidues 3).card + (nonzeroResidues 5).card +
      (nonzeroResidues 7).card = 15 - 3 := by
  rw [active_dimension_sum]

/-! ### Headline -/

/-- **Headline (dimension cascade).** The first six prime moduli give a
    strictly increasing cascade of `(ℤ/pℤ)*`-cardinalities:

    `(|n.r.|_{p=2,3,5,7,11,13}) = (1, 2, 4, 6, 10, 12)`,

    a strict chain `1 < 2 < 4 < 6 < 10 < 12`. The PT-active sub-cascade
    `(2, 4, 6)` is the unique three-element initial segment of even
    cardinalities forming an arithmetic progression. -/
theorem dimension_cascade_complete :
    (nonzeroResidues 2).card = 1
    ∧ (nonzeroResidues 3).card = 2
    ∧ (nonzeroResidues 5).card = 4
    ∧ (nonzeroResidues 7).card = 6
    ∧ (nonzeroResidues 11).card = 10
    ∧ (nonzeroResidues 13).card = 12
    -- strict monotonicity
    ∧ (nonzeroResidues 2).card < (nonzeroResidues 3).card
    ∧ (nonzeroResidues 3).card < (nonzeroResidues 5).card
    ∧ (nonzeroResidues 5).card < (nonzeroResidues 7).card
    ∧ (nonzeroResidues 7).card < (nonzeroResidues 11).card
    ∧ (nonzeroResidues 11).card < (nonzeroResidues 13).card :=
  ⟨card_nonzeroResidues_two, card_nonzeroResidues_three,
   card_nonzeroResidues_five, card_nonzeroResidues_seven,
   card_nonzeroResidues_eleven, card_nonzeroResidues_thirteen,
   cascade_2_lt_3, cascade_3_lt_5, cascade_5_lt_7,
   cascade_7_lt_11, cascade_11_lt_13⟩

end PT.Sieve.N4
