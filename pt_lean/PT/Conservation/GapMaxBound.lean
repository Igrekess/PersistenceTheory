/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Conservation.ConservationID
import PT.Conservation.ConservationIDExtensions
import PT.Conservation.CumulativeBoundsExtended
import Mathlib.Tactic

/-!
# Gap sequence — Upper bounds on `g_n` (Ch03 extension)

This file proves **upper bounds** on the PT prime-gap sequence
`g_n := ptPrimeExt(n+1) - ptPrimeExt(n)`, complementing the lower
bounds of `GapBoundedBelow.lean`.

The bounds proven here are *small-`N` finite witnesses* in the spirit
of Bertrand's postulate (`g_n < p_n`), restricted to `n ∈ {1, …, 10}`
where `ptPrimeExt` is concretely defined.

## Content

* **Individual gap values** for `n ∈ {1, …, 10}`:
  `(g_1, …, g_10) = (1, 2, 2, 4, 2, 4, 2, 4, 6, 2)`.
* **Bertrand-like upper bound** `g_n ≤ ptPrimeExt n` for `n ∈ {1, …, 10}`
  (weak form of `g_n < p_n`).
* **Global maximum** `g_n ≤ 6` over `n ∈ {1, …, 10}`.
* **Maximum attained at `n = 9`**: `g_9 = 6` and `g_9 > g_k` for all
  `k ∈ {1, …, 8} ∪ {10}`.
* **Local sum bounds**: `g_n + g_{n+1} ≤ 10` for `n ∈ {1, …, 9}`.
* **Ratio bound**: `g_n ≤ ∑_{k=1}^n g_k` (trivial, every gap is bounded
  by the cumulative sum since gaps are nonnegative).

The asymptotic Bertrand statement `g_n < p_n` is **not** in this file
— only the small-`N` finite witnesses.

## Reference

Monograph Chapter 3 §"Bornes supérieures sur les gaps", follow-up to
`CumulativeBoundsExtended.lean`.
-/

namespace PT.Conservation

open Finset

/-! ### Individual gap values on `n ∈ {1, …, 10}`

These are restatements of the `gapExt_*` lemmas from
`CumulativeBoundsExtended.lean`, gathered here for convenience and
self-containment. The full tuple
`(g_1, …, g_10) = (1, 2, 2, 4, 2, 4, 2, 4, 6, 2)` is bundled in
`gap_values_tuple` below. -/

/-- The first ten gaps of `ptPrimeExt`, packaged as a single conjunction. -/
theorem gap_values_tuple :
    gap ptPrimeExt 1  = 1
    ∧ gap ptPrimeExt 2  = 2
    ∧ gap ptPrimeExt 3  = 2
    ∧ gap ptPrimeExt 4  = 4
    ∧ gap ptPrimeExt 5  = 2
    ∧ gap ptPrimeExt 6  = 4
    ∧ gap ptPrimeExt 7  = 2
    ∧ gap ptPrimeExt 8  = 4
    ∧ gap ptPrimeExt 9  = 6
    ∧ gap ptPrimeExt 10 = 2 := by
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩ <;> (unfold gap; decide)

/-! ### Bertrand-like upper bound `g_n ≤ ptPrimeExt n` -/

/-- `g_1 = 1 ≤ p_1 = 2`. -/
theorem gapExt_one_le_ptPrime :
    gap ptPrimeExt 1 ≤ ptPrimeExt 1 := by unfold gap; decide

/-- `g_2 = 2 ≤ p_2 = 3`. -/
theorem gapExt_two_le_ptPrime :
    gap ptPrimeExt 2 ≤ ptPrimeExt 2 := by unfold gap; decide

/-- `g_3 = 2 ≤ p_3 = 5`. -/
theorem gapExt_three_le_ptPrime :
    gap ptPrimeExt 3 ≤ ptPrimeExt 3 := by unfold gap; decide

/-- `g_4 = 4 ≤ p_4 = 7`. -/
theorem gapExt_four_le_ptPrime :
    gap ptPrimeExt 4 ≤ ptPrimeExt 4 := by unfold gap; decide

/-- `g_5 = 2 ≤ p_5 = 11`. -/
theorem gapExt_five_le_ptPrime :
    gap ptPrimeExt 5 ≤ ptPrimeExt 5 := by unfold gap; decide

/-- `g_6 = 4 ≤ p_6 = 13`. -/
theorem gapExt_six_le_ptPrime :
    gap ptPrimeExt 6 ≤ ptPrimeExt 6 := by unfold gap; decide

/-- `g_7 = 2 ≤ p_7 = 17`. -/
theorem gapExt_seven_le_ptPrime :
    gap ptPrimeExt 7 ≤ ptPrimeExt 7 := by unfold gap; decide

/-- `g_8 = 4 ≤ p_8 = 19`. -/
theorem gapExt_eight_le_ptPrime :
    gap ptPrimeExt 8 ≤ ptPrimeExt 8 := by unfold gap; decide

/-- `g_9 = 6 ≤ p_9 = 23`. -/
theorem gapExt_nine_le_ptPrime :
    gap ptPrimeExt 9 ≤ ptPrimeExt 9 := by unfold gap; decide

/-- `g_10 = 2 ≤ p_10 = 29`. -/
theorem gapExt_ten_le_ptPrime :
    gap ptPrimeExt 10 ≤ ptPrimeExt 10 := by unfold gap; decide

/-- **Bertrand-like upper bound (weak form), small-`N` witnesses.**
    For `n ∈ {1, …, 10}`, `g_n ≤ p_n`. -/
theorem gap_le_ptPrime_small_n :
    gap ptPrimeExt 1  ≤ ptPrimeExt 1
    ∧ gap ptPrimeExt 2  ≤ ptPrimeExt 2
    ∧ gap ptPrimeExt 3  ≤ ptPrimeExt 3
    ∧ gap ptPrimeExt 4  ≤ ptPrimeExt 4
    ∧ gap ptPrimeExt 5  ≤ ptPrimeExt 5
    ∧ gap ptPrimeExt 6  ≤ ptPrimeExt 6
    ∧ gap ptPrimeExt 7  ≤ ptPrimeExt 7
    ∧ gap ptPrimeExt 8  ≤ ptPrimeExt 8
    ∧ gap ptPrimeExt 9  ≤ ptPrimeExt 9
    ∧ gap ptPrimeExt 10 ≤ ptPrimeExt 10 :=
  ⟨gapExt_one_le_ptPrime,   gapExt_two_le_ptPrime,
   gapExt_three_le_ptPrime, gapExt_four_le_ptPrime,
   gapExt_five_le_ptPrime,  gapExt_six_le_ptPrime,
   gapExt_seven_le_ptPrime, gapExt_eight_le_ptPrime,
   gapExt_nine_le_ptPrime,  gapExt_ten_le_ptPrime⟩

/-! ### Global maximum `g_n ≤ 6` on `n ∈ {1, …, 10}` -/

/-- `g_1 ≤ 6`. -/
theorem gapExt_one_le_six   : gap ptPrimeExt 1  ≤ 6 := by unfold gap; decide
/-- `g_2 ≤ 6`. -/
theorem gapExt_two_le_six   : gap ptPrimeExt 2  ≤ 6 := by unfold gap; decide
/-- `g_3 ≤ 6`. -/
theorem gapExt_three_le_six : gap ptPrimeExt 3  ≤ 6 := by unfold gap; decide
/-- `g_4 ≤ 6`. -/
theorem gapExt_four_le_six  : gap ptPrimeExt 4  ≤ 6 := by unfold gap; decide
/-- `g_5 ≤ 6`. -/
theorem gapExt_five_le_six  : gap ptPrimeExt 5  ≤ 6 := by unfold gap; decide
/-- `g_6 ≤ 6`. -/
theorem gapExt_six_le_six   : gap ptPrimeExt 6  ≤ 6 := by unfold gap; decide
/-- `g_7 ≤ 6`. -/
theorem gapExt_seven_le_six : gap ptPrimeExt 7  ≤ 6 := by unfold gap; decide
/-- `g_8 ≤ 6`. -/
theorem gapExt_eight_le_six : gap ptPrimeExt 8  ≤ 6 := by unfold gap; decide
/-- `g_9 ≤ 6`. -/
theorem gapExt_nine_le_six  : gap ptPrimeExt 9  ≤ 6 := by unfold gap; decide
/-- `g_10 ≤ 6`. -/
theorem gapExt_ten_le_six   : gap ptPrimeExt 10 ≤ 6 := by unfold gap; decide

/-- **Global maximum bound (`g_n ≤ 6` for `n ∈ {1, …, 10}`).** -/
theorem gap_max_le_six :
    gap ptPrimeExt 1  ≤ 6
    ∧ gap ptPrimeExt 2  ≤ 6
    ∧ gap ptPrimeExt 3  ≤ 6
    ∧ gap ptPrimeExt 4  ≤ 6
    ∧ gap ptPrimeExt 5  ≤ 6
    ∧ gap ptPrimeExt 6  ≤ 6
    ∧ gap ptPrimeExt 7  ≤ 6
    ∧ gap ptPrimeExt 8  ≤ 6
    ∧ gap ptPrimeExt 9  ≤ 6
    ∧ gap ptPrimeExt 10 ≤ 6 :=
  ⟨gapExt_one_le_six,   gapExt_two_le_six,
   gapExt_three_le_six, gapExt_four_le_six,
   gapExt_five_le_six,  gapExt_six_le_six,
   gapExt_seven_le_six, gapExt_eight_le_six,
   gapExt_nine_le_six,  gapExt_ten_le_six⟩

/-! ### Maximum attained strictly at `n = 9` -/

/-- `g_9 > g_1`. -/
theorem gapExt_nine_gt_one   : gap ptPrimeExt 9 > gap ptPrimeExt 1  := by
  unfold gap; decide
/-- `g_9 > g_2`. -/
theorem gapExt_nine_gt_two   : gap ptPrimeExt 9 > gap ptPrimeExt 2  := by
  unfold gap; decide
/-- `g_9 > g_3`. -/
theorem gapExt_nine_gt_three : gap ptPrimeExt 9 > gap ptPrimeExt 3  := by
  unfold gap; decide
/-- `g_9 > g_4`. -/
theorem gapExt_nine_gt_four  : gap ptPrimeExt 9 > gap ptPrimeExt 4  := by
  unfold gap; decide
/-- `g_9 > g_5`. -/
theorem gapExt_nine_gt_five  : gap ptPrimeExt 9 > gap ptPrimeExt 5  := by
  unfold gap; decide
/-- `g_9 > g_6`. -/
theorem gapExt_nine_gt_six   : gap ptPrimeExt 9 > gap ptPrimeExt 6  := by
  unfold gap; decide
/-- `g_9 > g_7`. -/
theorem gapExt_nine_gt_seven : gap ptPrimeExt 9 > gap ptPrimeExt 7  := by
  unfold gap; decide
/-- `g_9 > g_8`. -/
theorem gapExt_nine_gt_eight : gap ptPrimeExt 9 > gap ptPrimeExt 8  := by
  unfold gap; decide
/-- `g_9 > g_10`. -/
theorem gapExt_nine_gt_ten   : gap ptPrimeExt 9 > gap ptPrimeExt 10 := by
  unfold gap; decide

/-- **Strict maximum at `n = 9`.** The gap `g_9 = 6` is strictly larger
    than every other gap `g_k` for `k ∈ {1, …, 8} ∪ {10}`. -/
theorem gap_strict_max_at_nine :
    gap ptPrimeExt 9 = 6
    ∧ gap ptPrimeExt 9 > gap ptPrimeExt 1
    ∧ gap ptPrimeExt 9 > gap ptPrimeExt 2
    ∧ gap ptPrimeExt 9 > gap ptPrimeExt 3
    ∧ gap ptPrimeExt 9 > gap ptPrimeExt 4
    ∧ gap ptPrimeExt 9 > gap ptPrimeExt 5
    ∧ gap ptPrimeExt 9 > gap ptPrimeExt 6
    ∧ gap ptPrimeExt 9 > gap ptPrimeExt 7
    ∧ gap ptPrimeExt 9 > gap ptPrimeExt 8
    ∧ gap ptPrimeExt 9 > gap ptPrimeExt 10 :=
  ⟨gapExt_nine,
   gapExt_nine_gt_one,   gapExt_nine_gt_two,
   gapExt_nine_gt_three, gapExt_nine_gt_four,
   gapExt_nine_gt_five,  gapExt_nine_gt_six,
   gapExt_nine_gt_seven, gapExt_nine_gt_eight,
   gapExt_nine_gt_ten⟩

/-! ### Local sum bounds `g_n + g_{n+1} ≤ 10` -/

/-- `g_1 + g_2 = 1 + 2 = 3 ≤ 10`. -/
theorem gapExt_sum_1_2 : gap ptPrimeExt 1 + gap ptPrimeExt 2 ≤ 10 := by
  unfold gap; decide

/-- `g_2 + g_3 = 2 + 2 = 4 ≤ 10`. -/
theorem gapExt_sum_2_3 : gap ptPrimeExt 2 + gap ptPrimeExt 3 ≤ 10 := by
  unfold gap; decide

/-- `g_3 + g_4 = 2 + 4 = 6 ≤ 10`. -/
theorem gapExt_sum_3_4 : gap ptPrimeExt 3 + gap ptPrimeExt 4 ≤ 10 := by
  unfold gap; decide

/-- `g_4 + g_5 = 4 + 2 = 6 ≤ 10`. -/
theorem gapExt_sum_4_5 : gap ptPrimeExt 4 + gap ptPrimeExt 5 ≤ 10 := by
  unfold gap; decide

/-- `g_5 + g_6 = 2 + 4 = 6 ≤ 10`. -/
theorem gapExt_sum_5_6 : gap ptPrimeExt 5 + gap ptPrimeExt 6 ≤ 10 := by
  unfold gap; decide

/-- `g_6 + g_7 = 4 + 2 = 6 ≤ 10`. -/
theorem gapExt_sum_6_7 : gap ptPrimeExt 6 + gap ptPrimeExt 7 ≤ 10 := by
  unfold gap; decide

/-- `g_7 + g_8 = 2 + 4 = 6 ≤ 10`. -/
theorem gapExt_sum_7_8 : gap ptPrimeExt 7 + gap ptPrimeExt 8 ≤ 10 := by
  unfold gap; decide

/-- `g_8 + g_9 = 4 + 6 = 10 ≤ 10`. -/
theorem gapExt_sum_8_9 : gap ptPrimeExt 8 + gap ptPrimeExt 9 ≤ 10 := by
  unfold gap; decide

/-- `g_9 + g_10 = 6 + 2 = 8 ≤ 10`. -/
theorem gapExt_sum_9_10 : gap ptPrimeExt 9 + gap ptPrimeExt 10 ≤ 10 := by
  unfold gap; decide

/-- **Local sum bound.** For every consecutive pair `(g_n, g_{n+1})` with
    `n ∈ {1, …, 9}`, the sum is `≤ 10`. The bound is *attained* at
    `n = 8` since `g_8 + g_9 = 4 + 6 = 10`. -/
theorem gap_local_sum_le_ten :
    gap ptPrimeExt 1 + gap ptPrimeExt 2  ≤ 10
    ∧ gap ptPrimeExt 2 + gap ptPrimeExt 3  ≤ 10
    ∧ gap ptPrimeExt 3 + gap ptPrimeExt 4  ≤ 10
    ∧ gap ptPrimeExt 4 + gap ptPrimeExt 5  ≤ 10
    ∧ gap ptPrimeExt 5 + gap ptPrimeExt 6  ≤ 10
    ∧ gap ptPrimeExt 6 + gap ptPrimeExt 7  ≤ 10
    ∧ gap ptPrimeExt 7 + gap ptPrimeExt 8  ≤ 10
    ∧ gap ptPrimeExt 8 + gap ptPrimeExt 9  ≤ 10
    ∧ gap ptPrimeExt 9 + gap ptPrimeExt 10 ≤ 10 :=
  ⟨gapExt_sum_1_2, gapExt_sum_2_3, gapExt_sum_3_4,
   gapExt_sum_4_5, gapExt_sum_5_6, gapExt_sum_6_7,
   gapExt_sum_7_8, gapExt_sum_8_9, gapExt_sum_9_10⟩

/-! ### Ratio bound `g_n ≤ ∑_{k=1}^n g_k` -/

/-- `g_1 ≤ ∑_{k=1}^1 g_k = 1`. -/
theorem gapExt_one_le_cumul :
    gap ptPrimeExt 1 ≤ ∑ n ∈ Ico 1 2, gap ptPrimeExt n := by decide

/-- `g_2 ≤ ∑_{k=1}^2 g_k = 3`. -/
theorem gapExt_two_le_cumul :
    gap ptPrimeExt 2 ≤ ∑ n ∈ Ico 1 3, gap ptPrimeExt n := by decide

/-- `g_3 ≤ ∑_{k=1}^3 g_k = 5`. -/
theorem gapExt_three_le_cumul :
    gap ptPrimeExt 3 ≤ ∑ n ∈ Ico 1 4, gap ptPrimeExt n := by decide

/-- `g_4 ≤ ∑_{k=1}^4 g_k = 9`. -/
theorem gapExt_four_le_cumul :
    gap ptPrimeExt 4 ≤ ∑ n ∈ Ico 1 5, gap ptPrimeExt n := by decide

/-- `g_5 ≤ ∑_{k=1}^5 g_k = 11`. -/
theorem gapExt_five_le_cumul :
    gap ptPrimeExt 5 ≤ ∑ n ∈ Ico 1 6, gap ptPrimeExt n := by decide

/-- `g_6 ≤ ∑_{k=1}^6 g_k = 15`. -/
theorem gapExt_six_le_cumul :
    gap ptPrimeExt 6 ≤ ∑ n ∈ Ico 1 7, gap ptPrimeExt n := by decide

/-- `g_7 ≤ ∑_{k=1}^7 g_k = 17`. -/
theorem gapExt_seven_le_cumul :
    gap ptPrimeExt 7 ≤ ∑ n ∈ Ico 1 8, gap ptPrimeExt n := by decide

/-- `g_8 ≤ ∑_{k=1}^8 g_k = 21`. -/
theorem gapExt_eight_le_cumul :
    gap ptPrimeExt 8 ≤ ∑ n ∈ Ico 1 9, gap ptPrimeExt n := by decide

/-- `g_9 ≤ ∑_{k=1}^9 g_k = 27`. -/
theorem gapExt_nine_le_cumul :
    gap ptPrimeExt 9 ≤ ∑ n ∈ Ico 1 10, gap ptPrimeExt n := by decide

/-- `g_10 ≤ ∑_{k=1}^{10} g_k = 29`. -/
theorem gapExt_ten_le_cumul :
    gap ptPrimeExt 10 ≤ ∑ n ∈ Ico 1 11, gap ptPrimeExt n := by decide

/-- **Ratio bound.** Each individual gap `g_n` is bounded by the cumulative
    sum `∑_{k=1}^n g_k`, for `n ∈ {1, …, 10}`. This is a trivial but useful
    consequence of `g_k ≥ 0`. -/
theorem gap_le_cumul_small_n :
    gap ptPrimeExt 1  ≤ ∑ n ∈ Ico 1 2,  gap ptPrimeExt n
    ∧ gap ptPrimeExt 2  ≤ ∑ n ∈ Ico 1 3,  gap ptPrimeExt n
    ∧ gap ptPrimeExt 3  ≤ ∑ n ∈ Ico 1 4,  gap ptPrimeExt n
    ∧ gap ptPrimeExt 4  ≤ ∑ n ∈ Ico 1 5,  gap ptPrimeExt n
    ∧ gap ptPrimeExt 5  ≤ ∑ n ∈ Ico 1 6,  gap ptPrimeExt n
    ∧ gap ptPrimeExt 6  ≤ ∑ n ∈ Ico 1 7,  gap ptPrimeExt n
    ∧ gap ptPrimeExt 7  ≤ ∑ n ∈ Ico 1 8,  gap ptPrimeExt n
    ∧ gap ptPrimeExt 8  ≤ ∑ n ∈ Ico 1 9,  gap ptPrimeExt n
    ∧ gap ptPrimeExt 9  ≤ ∑ n ∈ Ico 1 10, gap ptPrimeExt n
    ∧ gap ptPrimeExt 10 ≤ ∑ n ∈ Ico 1 11, gap ptPrimeExt n :=
  ⟨gapExt_one_le_cumul,   gapExt_two_le_cumul,
   gapExt_three_le_cumul, gapExt_four_le_cumul,
   gapExt_five_le_cumul,  gapExt_six_le_cumul,
   gapExt_seven_le_cumul, gapExt_eight_le_cumul,
   gapExt_nine_le_cumul,  gapExt_ten_le_cumul⟩

/-! ### Headline summary -/

/-- **Headline (gap upper bounds, `n ∈ {1, …, 10}`).** For the eleven-prime
    PT sequence `ptPrimeExt`:

    * Values: `(g_1, …, g_10) = (1, 2, 2, 4, 2, 4, 2, 4, 6, 2)`.
    * Bertrand-like: `g_n ≤ p_n` for every `n ∈ {1, …, 10}`.
    * Global max: `g_n ≤ 6` for every `n ∈ {1, …, 10}`.
    * Strict max at `n = 9`: `g_9 = 6` and `g_9 > g_k` for all
      `k ∈ {1, …, 8} ∪ {10}`.
    * Local sum bound: `g_n + g_{n+1} ≤ 10` for `n ∈ {1, …, 9}`. -/
theorem gap_upper_bound_headline :
    -- (i) exact values
    (gap ptPrimeExt 1 = 1 ∧ gap ptPrimeExt 2 = 2
      ∧ gap ptPrimeExt 3 = 2 ∧ gap ptPrimeExt 4 = 4
      ∧ gap ptPrimeExt 5 = 2 ∧ gap ptPrimeExt 6 = 4
      ∧ gap ptPrimeExt 7 = 2 ∧ gap ptPrimeExt 8 = 4
      ∧ gap ptPrimeExt 9 = 6 ∧ gap ptPrimeExt 10 = 2)
    -- (ii) Bertrand-like: g_n ≤ p_n
    ∧ (gap ptPrimeExt 1 ≤ ptPrimeExt 1
        ∧ gap ptPrimeExt 2 ≤ ptPrimeExt 2
        ∧ gap ptPrimeExt 3 ≤ ptPrimeExt 3
        ∧ gap ptPrimeExt 4 ≤ ptPrimeExt 4
        ∧ gap ptPrimeExt 5 ≤ ptPrimeExt 5
        ∧ gap ptPrimeExt 6 ≤ ptPrimeExt 6
        ∧ gap ptPrimeExt 7 ≤ ptPrimeExt 7
        ∧ gap ptPrimeExt 8 ≤ ptPrimeExt 8
        ∧ gap ptPrimeExt 9 ≤ ptPrimeExt 9
        ∧ gap ptPrimeExt 10 ≤ ptPrimeExt 10)
    -- (iii) global max ≤ 6
    ∧ (gap ptPrimeExt 1 ≤ 6 ∧ gap ptPrimeExt 2 ≤ 6
        ∧ gap ptPrimeExt 3 ≤ 6 ∧ gap ptPrimeExt 4 ≤ 6
        ∧ gap ptPrimeExt 5 ≤ 6 ∧ gap ptPrimeExt 6 ≤ 6
        ∧ gap ptPrimeExt 7 ≤ 6 ∧ gap ptPrimeExt 8 ≤ 6
        ∧ gap ptPrimeExt 9 ≤ 6 ∧ gap ptPrimeExt 10 ≤ 6)
    -- (iv) strict max at n = 9
    ∧ (gap ptPrimeExt 9 = 6
        ∧ gap ptPrimeExt 9 > gap ptPrimeExt 1
        ∧ gap ptPrimeExt 9 > gap ptPrimeExt 2
        ∧ gap ptPrimeExt 9 > gap ptPrimeExt 3
        ∧ gap ptPrimeExt 9 > gap ptPrimeExt 4
        ∧ gap ptPrimeExt 9 > gap ptPrimeExt 5
        ∧ gap ptPrimeExt 9 > gap ptPrimeExt 6
        ∧ gap ptPrimeExt 9 > gap ptPrimeExt 7
        ∧ gap ptPrimeExt 9 > gap ptPrimeExt 8
        ∧ gap ptPrimeExt 9 > gap ptPrimeExt 10)
    -- (v) local sum bound ≤ 10
    ∧ (gap ptPrimeExt 1 + gap ptPrimeExt 2  ≤ 10
        ∧ gap ptPrimeExt 2 + gap ptPrimeExt 3  ≤ 10
        ∧ gap ptPrimeExt 3 + gap ptPrimeExt 4  ≤ 10
        ∧ gap ptPrimeExt 4 + gap ptPrimeExt 5  ≤ 10
        ∧ gap ptPrimeExt 5 + gap ptPrimeExt 6  ≤ 10
        ∧ gap ptPrimeExt 6 + gap ptPrimeExt 7  ≤ 10
        ∧ gap ptPrimeExt 7 + gap ptPrimeExt 8  ≤ 10
        ∧ gap ptPrimeExt 8 + gap ptPrimeExt 9  ≤ 10
        ∧ gap ptPrimeExt 9 + gap ptPrimeExt 10 ≤ 10) :=
  ⟨gap_values_tuple,
   gap_le_ptPrime_small_n,
   gap_max_le_six,
   gap_strict_max_at_nine,
   gap_local_sum_le_ten⟩

end PT.Conservation
