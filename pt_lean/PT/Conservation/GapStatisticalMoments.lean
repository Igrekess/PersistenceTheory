/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Conservation.ConservationID
import PT.Conservation.ConservationIDExtensions
import PT.Conservation.GapVarianceSmallN
import Mathlib.Algebra.BigOperators.Intervals
import Mathlib.Data.Rat.Defs
import Mathlib.Tactic

/-!
# Higher central moments of PT prime gaps — `N ∈ {3, 4}`

This file extends `PT.Conservation.GapVarianceSmallN` (mean and variance)
to the **third and fourth central moments** of the first `N` PT prime
gaps, together with the dimensionless **raw kurtosis ratio**.

For each `N`, we define

* `centralMoment3 N := (∑_{n=1}^{N} (g_n − meanGapsTo N)³) / N`,
* `centralMoment4 N := (∑_{n=1}^{N} (g_n − meanGapsTo N)⁴) / N`,
* `kurtosisRaw   N := centralMoment4 N / (varianceGapsTo N)²`.

The classical **skewness** ratio
`m_3 / σ³ = centralMoment3 / (varianceGapsTo)^(3/2)`
is *not* rational, so we expose only the un-normalised numerator
`centralMoment3 N` (the *scaled* skewness).

Concrete values proved here:

| N | mean | var   | m_3    | m_4     | kurtosisRaw |
|---|------|-------|--------|---------|-------------|
| 3 | 5/3  | 2/9   | -2/27  | 2/27    | 3/2         |
| 4 | 9/4  | 19/16 | 27/32  | 757/256 | 757/361     |

For `N = 3`, gaps are `(1, 2, 2)`, mean `5/3`, deviations
`(-2/3, 1/3, 1/3)`. For `N = 4`, gaps are `(1, 2, 2, 4)`,
mean `9/4`, deviations `(-5/4, -1/4, -1/4, 7/4)`.

## Reference

* `PT.Conservation.ConservationID` — `gap`.
* `PT.Conservation.ConservationIDExtensions` — `ptPrime`,
  `gap_one, gap_two, gap_three, gap_four`.
* `PT.Conservation.GapVarianceSmallN` — `meanGapsTo`, `varianceGapsTo`,
  `gap_*_rat`.
* Monograph Ch03 §3.1.
-/

namespace PT.Conservation

open Finset

/-! ### Definitions -/

/-- Third central moment of the first `N` PT prime gaps, in `ℚ`:
    `m_3(N) := (∑_{n=1}^{N} (g_n − meanGapsTo N)³) / N`. -/
def centralMoment3 (N : ℕ) : ℚ :=
  (∑ n ∈ Ico 1 (N + 1), ((gap ptPrime n : ℚ) - meanGapsTo N) ^ 3) / (N : ℚ)

/-- Fourth central moment of the first `N` PT prime gaps, in `ℚ`:
    `m_4(N) := (∑_{n=1}^{N} (g_n − meanGapsTo N)⁴) / N`. -/
def centralMoment4 (N : ℕ) : ℚ :=
  (∑ n ∈ Ico 1 (N + 1), ((gap ptPrime n : ℚ) - meanGapsTo N) ^ 4) / (N : ℚ)

/-- Raw (un-shifted) kurtosis ratio of the first `N` PT prime gaps:
    `m_4(N) / (varianceGapsTo N)²`. This is rational. -/
def kurtosisRaw (N : ℕ) : ℚ :=
  centralMoment4 N / (varianceGapsTo N) ^ 2

/-! ### `N = 3` — third and fourth central moments

Gaps `(1, 2, 2)`, mean `5/3`, deviations `(-2/3, 1/3, 1/3)`.

* `m_3(3) = ((-2/3)³ + (1/3)³ + (1/3)³) / 3 = (-6/27) / 3 = -2/27`.
* `m_4(3) = ((-2/3)⁴ + (1/3)⁴ + (1/3)⁴) / 3 = (18/81) / 3 = 2/27`.
-/

/-- Third central moment at `N = 3`: `-2/27`. -/
theorem centralMoment3_three : centralMoment3 3 = -2 / 27 := by
  unfold centralMoment3
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

/-- Fourth central moment at `N = 3`: `2/27`. -/
theorem centralMoment4_three : centralMoment4 3 = 2 / 27 := by
  unfold centralMoment4
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

/-- Raw kurtosis ratio at `N = 3`: `(2/27) / (2/9)² = (2/27) / (4/81) = 3/2`. -/
theorem kurtosisRaw_three : kurtosisRaw 3 = 3 / 2 := by
  unfold kurtosisRaw
  rw [centralMoment4_three, varianceGapsTo_three]
  norm_num

/-! ### `N = 4` — third and fourth central moments

Gaps `(1, 2, 2, 4)`, mean `9/4`, deviations `(-5/4, -1/4, -1/4, 7/4)`.

* `m_3(4) = ((-5/4)³ + (-1/4)³ + (-1/4)³ + (7/4)³) / 4`
          `= (-125/64 − 1/64 − 1/64 + 343/64) / 4`
          `= (216/64) / 4 = 27/32`.
* `m_4(4) = ((-5/4)⁴ + (-1/4)⁴ + (-1/4)⁴ + (7/4)⁴) / 4`
          `= (625/256 + 1/256 + 1/256 + 2401/256) / 4`
          `= (3028/256) / 4 = 757/256`.
-/

/-- Third central moment at `N = 4`: `27/32`. -/
theorem centralMoment3_four : centralMoment3 4 = 27 / 32 := by
  unfold centralMoment3
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

/-- Fourth central moment at `N = 4`: `757/256`. -/
theorem centralMoment4_four : centralMoment4 4 = 757 / 256 := by
  unfold centralMoment4
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

/-- Raw kurtosis ratio at `N = 4`:
    `(757/256) / (19/16)² = (757/256) / (361/256) = 757/361`. -/
theorem kurtosisRaw_four : kurtosisRaw 4 = 757 / 361 := by
  unfold kurtosisRaw
  rw [centralMoment4_four, varianceGapsTo_four]
  norm_num

/-! ### Sign of the third central moment

The third central moment is **negative at `N = 3`** (left-skew: the
unique larger gap `g_1 = 1` does not yet appear large enough to flip
the sign) and becomes **positive at `N = 4`** once `g_4 = 4` enters
and dominates the cubed deviation. -/

/-- `m_3(3) < 0` — left-skew at `N = 3`. -/
theorem centralMoment3_three_neg : centralMoment3 3 < 0 := by
  rw [centralMoment3_three]; norm_num

/-- `m_3(4) > 0` — right-skew at `N = 4`. -/
theorem centralMoment3_four_pos : 0 < centralMoment3 4 := by
  rw [centralMoment3_four]; norm_num

/-! ### Headline summary -/

/-- **Headline (small-`N` higher-moment summary).**

    * `N = 3`: `m_3 = -2/27`, `m_4 = 2/27`, `kurtosisRaw = 3/2`.
    * `N = 4`: `m_3 = 27/32`, `m_4 = 757/256`, `kurtosisRaw = 757/361`. -/
theorem gap_higher_moments_small_N_summary :
    centralMoment3 3 = -2 / 27 ∧ centralMoment4 3 = 2 / 27
    ∧ kurtosisRaw 3 = 3 / 2
    ∧ centralMoment3 4 = 27 / 32 ∧ centralMoment4 4 = 757 / 256
    ∧ kurtosisRaw 4 = 757 / 361 :=
  ⟨centralMoment3_three, centralMoment4_three, kurtosisRaw_three,
   centralMoment3_four, centralMoment4_four, kurtosisRaw_four⟩

end PT.Conservation
