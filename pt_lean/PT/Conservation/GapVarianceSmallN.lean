/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Conservation.ConservationID
import PT.Conservation.ConservationIDExtensions
import Mathlib.Algebra.BigOperators.Intervals
import Mathlib.Data.Rat.Defs
import Mathlib.Tactic

/-!
# Gap mean and variance on small `N` — `N ∈ {2, 3, 4}`

This file computes, in `ℚ`, the **arithmetic mean** and **(population)
variance** of the first `N` PT prime gaps, for `N ∈ {2, 3, 4}`.

Recall (from `PT.Conservation.ConservationIDExtensions`):
`(p_1, …, p_5) = (2, 3, 5, 7, 11)`, hence
`g_1 = 1, g_2 = 2, g_3 = 2, g_4 = 4`.

We define

* `meanGapsTo N     := (∑_{n=1}^{N} g_n) / N`,
* `sqDevSumGapsTo N := ∑_{n=1}^{N} (g_n − meanGapsTo N)²`,
* `varianceGapsTo N := sqDevSumGapsTo N / N` (population variance),
* `varianceGapsToBessel N := sqDevSumGapsTo N / (N − 1)`
  (Bessel-corrected, only valid for `N ≥ 2`).

Concrete values proved here:

| N | mean | population variance | Bessel variance |
|---|------|---------------------|-----------------|
| 2 | 3/2  | 1/4                 | 1/2             |
| 3 | 5/3  | 2/9                 | 1/3             |
| 4 | 9/4  | 19/16               | 19/12           |

All values are pure finite arithmetic in `ℚ`; we discharge them by
unrolling the `Finset.Ico` sums with `Finset.sum_Ico_succ_top`, then
unfolding `gap`/`ptPrime` and closing with `norm_num` (or `ring` for
the deviation sums).

## Reference

* `PT.Conservation.ConservationID` — `gap`.
* `PT.Conservation.ConservationIDExtensions` — `ptPrime`,
  `gap_one, gap_two, gap_three, gap_four`.
* `PT.Conservation.ConservationActivePrimes` — `activeGapsSum_eq_five`.
* Monograph Ch03 §3.1.
-/

namespace PT.Conservation

open Finset

/-! ### Definitions -/

/-- Mean of the first `N` PT prime gaps, computed in `ℚ`. -/
def meanGapsTo (N : ℕ) : ℚ :=
  (∑ n ∈ Ico 1 (N + 1), (gap ptPrime n : ℚ)) / (N : ℚ)

/-- Sum of squared deviations of the first `N` PT prime gaps from
    their mean, computed in `ℚ`. -/
def sqDevSumGapsTo (N : ℕ) : ℚ :=
  ∑ n ∈ Ico 1 (N + 1), ((gap ptPrime n : ℚ) - meanGapsTo N) ^ 2

/-- Population variance of the first `N` PT prime gaps. -/
def varianceGapsTo (N : ℕ) : ℚ :=
  sqDevSumGapsTo N / (N : ℚ)

/-- Bessel-corrected (sample) variance of the first `N` PT prime gaps.
    Meaningful only for `N ≥ 2`. -/
def varianceGapsToBessel (N : ℕ) : ℚ :=
  sqDevSumGapsTo N / ((N : ℚ) - 1)

/-! ### Helper — cast lemmas for the first four gaps in `ℚ` -/

theorem gap_one_rat : ((gap ptPrime 1 : ℤ) : ℚ) = 1 := by
  rw [gap_one]; norm_num

theorem gap_two_rat : ((gap ptPrime 2 : ℤ) : ℚ) = 2 := by
  rw [gap_two]; norm_num

theorem gap_three_rat : ((gap ptPrime 3 : ℤ) : ℚ) = 2 := by
  rw [gap_three]; norm_num

theorem gap_four_rat : ((gap ptPrime 4 : ℤ) : ℚ) = 4 := by
  rw [gap_four]; norm_num

/-! ### Helper — closed forms for the gap sums in `ℚ` -/

/-- Closed form for the gap sum on `Ico 1 3` (the `N = 2` numerator): `1 + 2 = 3`. -/
theorem sum_gap_ico_1_3_rat :
    (∑ n ∈ Ico 1 3, (gap ptPrime n : ℚ)) = 3 := by
  rw [show (3 : ℕ) = 2 + 1 from rfl,
      Finset.sum_Ico_succ_top (by omega : 1 ≤ 2),
      show (2 : ℕ) = 1 + 1 from rfl,
      Finset.sum_Ico_succ_top (by omega : 1 ≤ 1)]
  simp [gap_one_rat, gap_two_rat]
  norm_num

/-- Closed form for the gap sum on `Ico 1 4` (the `N = 3` numerator): `1 + 2 + 2 = 5`. -/
theorem sum_gap_ico_1_4_rat :
    (∑ n ∈ Ico 1 4, (gap ptPrime n : ℚ)) = 5 := by
  rw [show (4 : ℕ) = 3 + 1 from rfl,
      Finset.sum_Ico_succ_top (by omega : 1 ≤ 3),
      sum_gap_ico_1_3_rat, gap_three_rat]
  norm_num

/-- Closed form for the gap sum on `Ico 1 5` (the `N = 4` numerator): `1 + 2 + 2 + 4 = 9`. -/
theorem sum_gap_ico_1_5_rat :
    (∑ n ∈ Ico 1 5, (gap ptPrime n : ℚ)) = 9 := by
  rw [show (5 : ℕ) = 4 + 1 from rfl,
      Finset.sum_Ico_succ_top (by omega : 1 ≤ 4),
      sum_gap_ico_1_4_rat, gap_four_rat]
  norm_num

/-! ### `N = 2` — mean `3/2`, variance `1/4` -/

/-- Mean at `N = 2`: `(1 + 2) / 2 = 3/2`. -/
theorem meanGapsTo_two : meanGapsTo 2 = 3 / 2 := by
  unfold meanGapsTo
  rw [show (2 : ℕ) + 1 = 3 from rfl, sum_gap_ico_1_3_rat]
  norm_num

/-- Squared-deviation sum at `N = 2`:
    `(1 − 3/2)² + (2 − 3/2)² = 1/2`. -/
theorem sqDevSumGapsTo_two : sqDevSumGapsTo 2 = 1 / 2 := by
  unfold sqDevSumGapsTo
  rw [meanGapsTo_two,
      show (2 : ℕ) + 1 = 3 from rfl,
      show (3 : ℕ) = 2 + 1 from rfl,
      Finset.sum_Ico_succ_top (by omega : 1 ≤ 2),
      show (2 : ℕ) = 1 + 1 from rfl,
      Finset.sum_Ico_succ_top (by omega : 1 ≤ 1)]
  simp [gap_one_rat, gap_two_rat]
  norm_num

/-- Population variance at `N = 2`: `(1/2) / 2 = 1/4`. -/
theorem varianceGapsTo_two : varianceGapsTo 2 = 1 / 4 := by
  unfold varianceGapsTo
  rw [sqDevSumGapsTo_two]
  norm_num

/-- Bessel-corrected variance at `N = 2`: `(1/2) / 1 = 1/2`. -/
theorem varianceGapsToBessel_two : varianceGapsToBessel 2 = 1 / 2 := by
  unfold varianceGapsToBessel
  rw [sqDevSumGapsTo_two]
  norm_num

/-! ### `N = 3` — mean `5/3`, variance `2/9` -/

/-- Mean at `N = 3`: `(1 + 2 + 2) / 3 = 5/3`. -/
theorem meanGapsTo_three : meanGapsTo 3 = 5 / 3 := by
  unfold meanGapsTo
  rw [show (3 : ℕ) + 1 = 4 from rfl, sum_gap_ico_1_4_rat]
  norm_num

/-- Squared-deviation sum at `N = 3`:
    `(1 − 5/3)² + (2 − 5/3)² + (2 − 5/3)² = 4/9 + 1/9 + 1/9 = 2/3`. -/
theorem sqDevSumGapsTo_three : sqDevSumGapsTo 3 = 2 / 3 := by
  unfold sqDevSumGapsTo
  rw [meanGapsTo_three,
      show (3 : ℕ) + 1 = 4 from rfl,
      show (4 : ℕ) = 3 + 1 from rfl,
      Finset.sum_Ico_succ_top (by omega : 1 ≤ 3),
      show (3 : ℕ) = 2 + 1 from rfl,
      Finset.sum_Ico_succ_top (by omega : 1 ≤ 2),
      show (2 : ℕ) = 1 + 1 from rfl,
      Finset.sum_Ico_succ_top (by omega : 1 ≤ 1)]
  simp [gap_one_rat, gap_two_rat, gap_three_rat]
  norm_num

/-- Population variance at `N = 3`: `(2/3) / 3 = 2/9`. -/
theorem varianceGapsTo_three : varianceGapsTo 3 = 2 / 9 := by
  unfold varianceGapsTo
  rw [sqDevSumGapsTo_three]
  norm_num

/-- Bessel-corrected variance at `N = 3`: `(2/3) / 2 = 1/3`. -/
theorem varianceGapsToBessel_three : varianceGapsToBessel 3 = 1 / 3 := by
  unfold varianceGapsToBessel
  rw [sqDevSumGapsTo_three]
  norm_num

/-! ### `N = 4` — mean `9/4`, variance `19/16` -/

/-- Mean at `N = 4`: `(1 + 2 + 2 + 4) / 4 = 9/4`. -/
theorem meanGapsTo_four : meanGapsTo 4 = 9 / 4 := by
  unfold meanGapsTo
  rw [show (4 : ℕ) + 1 = 5 from rfl, sum_gap_ico_1_5_rat]
  norm_num

/-- Squared-deviation sum at `N = 4`:
    `(1 − 9/4)² + (2 − 9/4)² + (2 − 9/4)² + (4 − 9/4)²
     = 25/16 + 1/16 + 1/16 + 49/16 = 19/4`. -/
theorem sqDevSumGapsTo_four : sqDevSumGapsTo 4 = 19 / 4 := by
  unfold sqDevSumGapsTo
  rw [meanGapsTo_four,
      show (4 : ℕ) + 1 = 5 from rfl,
      show (5 : ℕ) = 4 + 1 from rfl,
      Finset.sum_Ico_succ_top (by omega : 1 ≤ 4),
      show (4 : ℕ) = 3 + 1 from rfl,
      Finset.sum_Ico_succ_top (by omega : 1 ≤ 3),
      show (3 : ℕ) = 2 + 1 from rfl,
      Finset.sum_Ico_succ_top (by omega : 1 ≤ 2),
      show (2 : ℕ) = 1 + 1 from rfl,
      Finset.sum_Ico_succ_top (by omega : 1 ≤ 1)]
  simp [gap_one_rat, gap_two_rat, gap_three_rat, gap_four_rat]
  norm_num

/-- Population variance at `N = 4`: `(19/4) / 4 = 19/16`. -/
theorem varianceGapsTo_four : varianceGapsTo 4 = 19 / 16 := by
  unfold varianceGapsTo
  rw [sqDevSumGapsTo_four]
  norm_num

/-- Bessel-corrected variance at `N = 4`: `(19/4) / 3 = 19/12`. -/
theorem varianceGapsToBessel_four : varianceGapsToBessel 4 = 19 / 12 := by
  unfold varianceGapsToBessel
  rw [sqDevSumGapsTo_four]
  norm_num

/-! ### Cross-checks: consistency with the conservation identity

The numerator of `meanGapsTo N` is exactly the conservation cumulative
sum from `ConservationIDExtensions`, namely `p_{N+1} − 2`. The closed
forms `sum_gap_ico_1_{3,4,5}_rat` above act as the `ℚ`-cast version of
`conservation_N{2,3,4}`. -/

/-- Numerator of the mean at `N = 2` equals `p_3 − 2 = 3`. -/
theorem meanGapsTo_two_numerator :
    (∑ n ∈ Ico 1 3, (gap ptPrime n : ℚ)) = 3 :=
  sum_gap_ico_1_3_rat

/-- Numerator of the mean at `N = 3` equals `p_4 − 2 = 5`. -/
theorem meanGapsTo_three_numerator :
    (∑ n ∈ Ico 1 4, (gap ptPrime n : ℚ)) = 5 :=
  sum_gap_ico_1_4_rat

/-- Numerator of the mean at `N = 4` equals `p_5 − 2 = 9`. -/
theorem meanGapsTo_four_numerator :
    (∑ n ∈ Ico 1 5, (gap ptPrime n : ℚ)) = 9 :=
  sum_gap_ico_1_5_rat

/-! ### Headline summary -/

/-- **Headline (small-`N` gap mean/variance summary).** Putting it all
    together:

    * `N = 2`: `mean = 3/2`, `var = 1/4`
    * `N = 3`: `mean = 5/3`, `var = 2/9`
    * `N = 4`: `mean = 9/4`, `var = 19/16` -/
theorem gap_mean_variance_small_N_summary :
    meanGapsTo 2 = 3 / 2 ∧ varianceGapsTo 2 = 1 / 4
    ∧ meanGapsTo 3 = 5 / 3 ∧ varianceGapsTo 3 = 2 / 9
    ∧ meanGapsTo 4 = 9 / 4 ∧ varianceGapsTo 4 = 19 / 16 :=
  ⟨meanGapsTo_two, varianceGapsTo_two,
   meanGapsTo_three, varianceGapsTo_three,
   meanGapsTo_four, varianceGapsTo_four⟩

end PT.Conservation
