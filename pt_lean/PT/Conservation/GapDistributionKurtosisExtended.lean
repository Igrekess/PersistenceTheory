/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Conservation.ConservationID
import PT.Conservation.ConservationIDExtensions
import PT.Conservation.CumulativeBoundsExtended
import PT.Conservation.GapMomentsNExtended
import PT.Conservation.GapDistributionVariance
import PT.Conservation.GapDistributionMeanExtended
import PT.Conservation.GapDistributionVarianceExtended
import PT.Conservation.GapDistributionSkewExtended
import PT.Conservation.GapStatisticalMoments
import Mathlib.Algebra.BigOperators.Intervals
import Mathlib.Data.Rat.Defs
import Mathlib.Tactic

/-!
# Distributional kurtosis (fourth central moment) of the PT gap
# distribution — extension `N ∈ {5, …, 10}`

This file is the **fourth-moment** companion to
`PT.Conservation.GapDistributionSkewExtended`. On the eleven-prime
sequence `ptPrimeExt`, it computes the **distributional fourth central
moment** of the index `n` under the probability distribution
`n ↦ g_n / S_N`:

  `gapsDistCentralMoment4ExtQ N
      := (Σ_{n=1}^{N} (n − μ_N)⁴ · g_n) / S_N`,

where `μ_N := gapsDistMeanExtQ N` and `S_N := Σ_{n=1}^{N} g_n`.

Unlike the third central moment, this quantity is **always strictly
positive** (it is a sum of fourth powers with positive weights `g_n > 0`
and at least one non-zero deviation, divided by the positive cumulative
sum `S_N`).

The **dimensionless kurtosis ratio**

  `kurtosisExtQ N := gapsDistCentralMoment4ExtQ N / (gapsDistVarianceExtQ N)²`

measures the relative tail weight of the index distribution. For a
Gaussian distribution this ratio equals `3`; values below `3` are
*platykurtic* (light tails, more uniform than a Gaussian), and values
above `3` are *leptokurtic* (heavy tails).

With gaps `(g_1, …, g_{10}) = (1, 2, 2, 4, 2, 4, 2, 4, 6, 2)`, the six
exact rationals are

| `N`  | `S_N` | `μ_N`    | `gapsDistCentralMoment4ExtQ N` | `kurtosisExtQ N`  |
|------|-------|----------|--------------------------------|-------------------|
| `5`  | `11`  | `37/11`  | `70754 / 14641`                | `35377/16562`     |
| `6`  | `15`  | `61/15`  | `208034 / 16875`               | `312051/153458`   |
| `7`  | `17`  | `75/17`  | `1604890 / 83521`              | `802445/392498`   |
| `8`  | `21`  | `107/21` | `2540354 / 64827`              | `3810531/1940450` |
| `9`  | `27`  | `161/27` | `12445406 / 177147`            | `18668109/9919058`|
| `10` | `29`  | `181/29` | `61313616 / 707281`            | `425789/222784`   |

Numerically the six kurtosis ratios are

  `K_5 ≈ 2.136, K_6 ≈ 2.033, K_7 ≈ 2.044,`
  `K_8 ≈ 1.964, K_9 ≈ 1.882, K_{10} ≈ 1.911`.

All six are **strictly less than `3`**: the gap-weighted index
distribution on `{1, …, N}` is **platykurtic** across the extended
range — the tails are systematically lighter than a Gaussian of equal
variance, reflecting the bounded support of the index distribution.

All values are exact rationals, discharged by unfolding the `Finset.Ico`
sums step by step and closing with `norm_num`.

## Reference

* `PT.Conservation.CumulativeBoundsExtended` — `ptPrimeExt`, gaps,
  cumulative sums `S_N`.
* `PT.Conservation.GapDistributionMeanExtended` — `gapsDistMeanExtQ N`
  and the six exact mean values used here.
* `PT.Conservation.GapDistributionVarianceExtended` — `gapsDistVarianceExtQ N`
  and the rational gap lemmas `gapExt_one_rat, …, gapExt_ten_rat` reused
  below.
* `PT.Conservation.GapDistributionSkewExtended` — the third central
  moment companion file.
* `PT.Conservation.GapStatisticalMoments` — `centralMoment4` on
  `ptPrime` (small-`N`, counting measure).
-/

namespace PT.Conservation

open Finset

/-! ### Distributional fourth central moment on the extended range -/

/-- **Extended distributional fourth central moment.** Fourth central
    moment of the index `n` under the probability distribution
    `n ↦ g_n / S_N` on the eleven-prime sequence `ptPrimeExt`:

      `gapsDistCentralMoment4ExtQ N
          = (Σ_{n=1}^{N} (n − μ_N)⁴ · g_n) / S_N`,

    where `μ_N = gapsDistMeanExtQ N`. Always strictly positive on the
    extended range since the deviations are not all zero. -/
def gapsDistCentralMoment4ExtQ (N : ℕ) : ℚ :=
  (∑ n ∈ Ico (1 : ℕ) (N + 1),
      ((n : ℚ) - gapsDistMeanExtQ N) ^ 4 * (gap ptPrimeExt n : ℚ))
    / cumGapExtQ N

/-! ### Six exact values `N ∈ {5, …, 10}`

Each fourth central moment is discharged by unfolding the `Ico` sum
step by step (reducing it to a sum of `N` explicit terms) and closing
with `norm_num`, after rewriting `μ_N = gapsDistMeanExtQ N` and
`S_N = cumGapExtQ N` to their exact rational values.
-/

/-- **Distributional fourth central moment at `N = 5`.**
    `70754 / 14641` with `14641 = 11⁴`. -/
theorem gapsDistCentralMoment4ExtQ_five :
    gapsDistCentralMoment4ExtQ 5 = 70754 / 14641 := by
  unfold gapsDistCentralMoment4ExtQ
  rw [gapsDistMeanExtQ_five, cumGapExtQ_five,
      show (5 : ℕ) + 1 = 6 from rfl,
      show (6 : ℕ) = 5 + 1 from rfl,
      Finset.sum_Ico_succ_top (by omega : 1 ≤ 5),
      show (5 : ℕ) = 4 + 1 from rfl,
      Finset.sum_Ico_succ_top (by omega : 1 ≤ 4),
      show (4 : ℕ) = 3 + 1 from rfl,
      Finset.sum_Ico_succ_top (by omega : 1 ≤ 3),
      show (3 : ℕ) = 2 + 1 from rfl,
      Finset.sum_Ico_succ_top (by omega : 1 ≤ 2),
      show (2 : ℕ) = 1 + 1 from rfl,
      Finset.sum_Ico_succ_top (by omega : 1 ≤ 1)]
  simp [gapExt_one_rat, gapExt_two_rat, gapExt_three_rat,
        gapExt_four_rat, gapExt_five_rat]
  norm_num

/-- **Distributional fourth central moment at `N = 6`.** `208034 / 16875`. -/
theorem gapsDistCentralMoment4ExtQ_six :
    gapsDistCentralMoment4ExtQ 6 = 208034 / 16875 := by
  unfold gapsDistCentralMoment4ExtQ
  rw [gapsDistMeanExtQ_six, cumGapExtQ_six,
      show (6 : ℕ) + 1 = 7 from rfl,
      show (7 : ℕ) = 6 + 1 from rfl,
      Finset.sum_Ico_succ_top (by omega : 1 ≤ 6),
      show (6 : ℕ) = 5 + 1 from rfl,
      Finset.sum_Ico_succ_top (by omega : 1 ≤ 5),
      show (5 : ℕ) = 4 + 1 from rfl,
      Finset.sum_Ico_succ_top (by omega : 1 ≤ 4),
      show (4 : ℕ) = 3 + 1 from rfl,
      Finset.sum_Ico_succ_top (by omega : 1 ≤ 3),
      show (3 : ℕ) = 2 + 1 from rfl,
      Finset.sum_Ico_succ_top (by omega : 1 ≤ 2),
      show (2 : ℕ) = 1 + 1 from rfl,
      Finset.sum_Ico_succ_top (by omega : 1 ≤ 1)]
  simp [gapExt_one_rat, gapExt_two_rat, gapExt_three_rat,
        gapExt_four_rat, gapExt_five_rat, gapExt_six_rat]
  norm_num

/-- **Distributional fourth central moment at `N = 7`.** `1604890 / 83521`
    with `83521 = 17⁴`. -/
theorem gapsDistCentralMoment4ExtQ_seven :
    gapsDistCentralMoment4ExtQ 7 = 1604890 / 83521 := by
  unfold gapsDistCentralMoment4ExtQ
  rw [gapsDistMeanExtQ_seven, cumGapExtQ_seven,
      show (7 : ℕ) + 1 = 8 from rfl,
      show (8 : ℕ) = 7 + 1 from rfl,
      Finset.sum_Ico_succ_top (by omega : 1 ≤ 7),
      show (7 : ℕ) = 6 + 1 from rfl,
      Finset.sum_Ico_succ_top (by omega : 1 ≤ 6),
      show (6 : ℕ) = 5 + 1 from rfl,
      Finset.sum_Ico_succ_top (by omega : 1 ≤ 5),
      show (5 : ℕ) = 4 + 1 from rfl,
      Finset.sum_Ico_succ_top (by omega : 1 ≤ 4),
      show (4 : ℕ) = 3 + 1 from rfl,
      Finset.sum_Ico_succ_top (by omega : 1 ≤ 3),
      show (3 : ℕ) = 2 + 1 from rfl,
      Finset.sum_Ico_succ_top (by omega : 1 ≤ 2),
      show (2 : ℕ) = 1 + 1 from rfl,
      Finset.sum_Ico_succ_top (by omega : 1 ≤ 1)]
  simp [gapExt_one_rat, gapExt_two_rat, gapExt_three_rat,
        gapExt_four_rat, gapExt_five_rat, gapExt_six_rat,
        gapExt_seven_rat]
  norm_num

/-- **Distributional fourth central moment at `N = 8`.** `2540354 / 64827`. -/
theorem gapsDistCentralMoment4ExtQ_eight :
    gapsDistCentralMoment4ExtQ 8 = 2540354 / 64827 := by
  unfold gapsDistCentralMoment4ExtQ
  rw [gapsDistMeanExtQ_eight, cumGapExtQ_eight,
      show (8 : ℕ) + 1 = 9 from rfl,
      show (9 : ℕ) = 8 + 1 from rfl,
      Finset.sum_Ico_succ_top (by omega : 1 ≤ 8),
      show (8 : ℕ) = 7 + 1 from rfl,
      Finset.sum_Ico_succ_top (by omega : 1 ≤ 7),
      show (7 : ℕ) = 6 + 1 from rfl,
      Finset.sum_Ico_succ_top (by omega : 1 ≤ 6),
      show (6 : ℕ) = 5 + 1 from rfl,
      Finset.sum_Ico_succ_top (by omega : 1 ≤ 5),
      show (5 : ℕ) = 4 + 1 from rfl,
      Finset.sum_Ico_succ_top (by omega : 1 ≤ 4),
      show (4 : ℕ) = 3 + 1 from rfl,
      Finset.sum_Ico_succ_top (by omega : 1 ≤ 3),
      show (3 : ℕ) = 2 + 1 from rfl,
      Finset.sum_Ico_succ_top (by omega : 1 ≤ 2),
      show (2 : ℕ) = 1 + 1 from rfl,
      Finset.sum_Ico_succ_top (by omega : 1 ≤ 1)]
  simp [gapExt_one_rat, gapExt_two_rat, gapExt_three_rat,
        gapExt_four_rat, gapExt_five_rat, gapExt_six_rat,
        gapExt_seven_rat, gapExt_eight_rat]
  norm_num

/-- **Distributional fourth central moment at `N = 9`.** `12445406 / 177147`. -/
theorem gapsDistCentralMoment4ExtQ_nine :
    gapsDistCentralMoment4ExtQ 9 = 12445406 / 177147 := by
  unfold gapsDistCentralMoment4ExtQ
  rw [gapsDistMeanExtQ_nine, cumGapExtQ_nine,
      show (9 : ℕ) + 1 = 10 from rfl,
      show (10 : ℕ) = 9 + 1 from rfl,
      Finset.sum_Ico_succ_top (by omega : 1 ≤ 9),
      show (9 : ℕ) = 8 + 1 from rfl,
      Finset.sum_Ico_succ_top (by omega : 1 ≤ 8),
      show (8 : ℕ) = 7 + 1 from rfl,
      Finset.sum_Ico_succ_top (by omega : 1 ≤ 7),
      show (7 : ℕ) = 6 + 1 from rfl,
      Finset.sum_Ico_succ_top (by omega : 1 ≤ 6),
      show (6 : ℕ) = 5 + 1 from rfl,
      Finset.sum_Ico_succ_top (by omega : 1 ≤ 5),
      show (5 : ℕ) = 4 + 1 from rfl,
      Finset.sum_Ico_succ_top (by omega : 1 ≤ 4),
      show (4 : ℕ) = 3 + 1 from rfl,
      Finset.sum_Ico_succ_top (by omega : 1 ≤ 3),
      show (3 : ℕ) = 2 + 1 from rfl,
      Finset.sum_Ico_succ_top (by omega : 1 ≤ 2),
      show (2 : ℕ) = 1 + 1 from rfl,
      Finset.sum_Ico_succ_top (by omega : 1 ≤ 1)]
  simp [gapExt_one_rat, gapExt_two_rat, gapExt_three_rat,
        gapExt_four_rat, gapExt_five_rat, gapExt_six_rat,
        gapExt_seven_rat, gapExt_eight_rat, gapExt_nine_rat]
  norm_num

/-- **Distributional fourth central moment at `N = 10`.** `61313616 / 707281`
    with `707281 = 29⁴`. -/
theorem gapsDistCentralMoment4ExtQ_ten :
    gapsDistCentralMoment4ExtQ 10 = 61313616 / 707281 := by
  unfold gapsDistCentralMoment4ExtQ
  rw [gapsDistMeanExtQ_ten, cumGapExtQ_ten,
      show (10 : ℕ) + 1 = 11 from rfl,
      show (11 : ℕ) = 10 + 1 from rfl,
      Finset.sum_Ico_succ_top (by omega : 1 ≤ 10),
      show (10 : ℕ) = 9 + 1 from rfl,
      Finset.sum_Ico_succ_top (by omega : 1 ≤ 9),
      show (9 : ℕ) = 8 + 1 from rfl,
      Finset.sum_Ico_succ_top (by omega : 1 ≤ 8),
      show (8 : ℕ) = 7 + 1 from rfl,
      Finset.sum_Ico_succ_top (by omega : 1 ≤ 7),
      show (7 : ℕ) = 6 + 1 from rfl,
      Finset.sum_Ico_succ_top (by omega : 1 ≤ 6),
      show (6 : ℕ) = 5 + 1 from rfl,
      Finset.sum_Ico_succ_top (by omega : 1 ≤ 5),
      show (5 : ℕ) = 4 + 1 from rfl,
      Finset.sum_Ico_succ_top (by omega : 1 ≤ 4),
      show (4 : ℕ) = 3 + 1 from rfl,
      Finset.sum_Ico_succ_top (by omega : 1 ≤ 3),
      show (3 : ℕ) = 2 + 1 from rfl,
      Finset.sum_Ico_succ_top (by omega : 1 ≤ 2),
      show (2 : ℕ) = 1 + 1 from rfl,
      Finset.sum_Ico_succ_top (by omega : 1 ≤ 1)]
  simp [gapExt_one_rat, gapExt_two_rat, gapExt_three_rat,
        gapExt_four_rat, gapExt_five_rat, gapExt_six_rat,
        gapExt_seven_rat, gapExt_eight_rat, gapExt_nine_rat,
        gapExt_ten_rat]
  norm_num

/-! ### Strict positivity of the fourth central moment

The fourth central moment is a sum of fourth powers with positive
weights `g_n > 0`. On the extended range `N ∈ {5, …, 10}` at least one
deviation `(n − μ_N)` is non-zero (the mean `μ_N` is irrational w.r.t.
the integer indices in the sense that it never coincides with all
indices simultaneously), so the sum is strictly positive. Dividing by
the positive `S_N > 0` keeps the sign. We discharge each case by
rewriting to the explicit rational value. -/

/-- `gapsDistCentralMoment4ExtQ 5 > 0`. -/
theorem gapsDistCentralMoment4ExtQ_five_pos :
    gapsDistCentralMoment4ExtQ 5 > 0 := by
  rw [gapsDistCentralMoment4ExtQ_five]; norm_num

/-- `gapsDistCentralMoment4ExtQ 6 > 0`. -/
theorem gapsDistCentralMoment4ExtQ_six_pos :
    gapsDistCentralMoment4ExtQ 6 > 0 := by
  rw [gapsDistCentralMoment4ExtQ_six]; norm_num

/-- `gapsDistCentralMoment4ExtQ 7 > 0`. -/
theorem gapsDistCentralMoment4ExtQ_seven_pos :
    gapsDistCentralMoment4ExtQ 7 > 0 := by
  rw [gapsDistCentralMoment4ExtQ_seven]; norm_num

/-- `gapsDistCentralMoment4ExtQ 8 > 0`. -/
theorem gapsDistCentralMoment4ExtQ_eight_pos :
    gapsDistCentralMoment4ExtQ 8 > 0 := by
  rw [gapsDistCentralMoment4ExtQ_eight]; norm_num

/-- `gapsDistCentralMoment4ExtQ 9 > 0`. -/
theorem gapsDistCentralMoment4ExtQ_nine_pos :
    gapsDistCentralMoment4ExtQ 9 > 0 := by
  rw [gapsDistCentralMoment4ExtQ_nine]; norm_num

/-- `gapsDistCentralMoment4ExtQ 10 > 0`. -/
theorem gapsDistCentralMoment4ExtQ_ten_pos :
    gapsDistCentralMoment4ExtQ 10 > 0 := by
  rw [gapsDistCentralMoment4ExtQ_ten]; norm_num

/-! ### Dimensionless kurtosis ratio

The kurtosis ratio `K_N := m_4(N) / σ²(N)²` is the standard
scale-invariant measure of the tail weight of the distribution.
For a Gaussian, `K = 3`. -/

/-- **Kurtosis ratio.** Dimensionless ratio
    `m_4(N) / (σ²(N))²` measuring tail weight of the gap-weighted
    index distribution. -/
def kurtosisExtQ (N : ℕ) : ℚ :=
  gapsDistCentralMoment4ExtQ N / (gapsDistVarianceExtQ N) ^ 2

/-- **Kurtosis at `N = 5`.** `35377 / 16562 ≈ 2.136`. -/
theorem kurtosisExtQ_five : kurtosisExtQ 5 = 35377 / 16562 := by
  unfold kurtosisExtQ
  rw [gapsDistCentralMoment4ExtQ_five, gapsDistVarianceExtQ_five]
  norm_num

/-- **Kurtosis at `N = 6`.** `312051 / 153458 ≈ 2.033`. -/
theorem kurtosisExtQ_six : kurtosisExtQ 6 = 312051 / 153458 := by
  unfold kurtosisExtQ
  rw [gapsDistCentralMoment4ExtQ_six, gapsDistVarianceExtQ_six]
  norm_num

/-- **Kurtosis at `N = 7`.** `802445 / 392498 ≈ 2.044`. -/
theorem kurtosisExtQ_seven : kurtosisExtQ 7 = 802445 / 392498 := by
  unfold kurtosisExtQ
  rw [gapsDistCentralMoment4ExtQ_seven, gapsDistVarianceExtQ_seven]
  norm_num

/-- **Kurtosis at `N = 8`.** `3810531 / 1940450 ≈ 1.964`. -/
theorem kurtosisExtQ_eight : kurtosisExtQ 8 = 3810531 / 1940450 := by
  unfold kurtosisExtQ
  rw [gapsDistCentralMoment4ExtQ_eight, gapsDistVarianceExtQ_eight]
  norm_num

/-- **Kurtosis at `N = 9`.** `18668109 / 9919058 ≈ 1.882`. -/
theorem kurtosisExtQ_nine : kurtosisExtQ 9 = 18668109 / 9919058 := by
  unfold kurtosisExtQ
  rw [gapsDistCentralMoment4ExtQ_nine, gapsDistVarianceExtQ_nine]
  norm_num

/-- **Kurtosis at `N = 10`.** `425789 / 222784 ≈ 1.911`. -/
theorem kurtosisExtQ_ten : kurtosisExtQ 10 = 425789 / 222784 := by
  unfold kurtosisExtQ
  rw [gapsDistCentralMoment4ExtQ_ten, gapsDistVarianceExtQ_ten]
  norm_num

/-! ### Comparison with the Gaussian reference value `K = 3`

All six kurtosis ratios are **strictly less than `3`**: the gap-weighted
index distribution on `{1, …, N}` is **platykurtic** across the extended
range — the tails are lighter than a Gaussian of equal variance,
reflecting the bounded support `{1, …, N}` of the index distribution. -/

/-- `kurtosisExtQ 5 < 3` — platykurtic at `N = 5`. -/
theorem kurtosisExtQ_five_lt_three : kurtosisExtQ 5 < 3 := by
  rw [kurtosisExtQ_five]; norm_num

/-- `kurtosisExtQ 6 < 3` — platykurtic at `N = 6`. -/
theorem kurtosisExtQ_six_lt_three : kurtosisExtQ 6 < 3 := by
  rw [kurtosisExtQ_six]; norm_num

/-- `kurtosisExtQ 7 < 3` — platykurtic at `N = 7`. -/
theorem kurtosisExtQ_seven_lt_three : kurtosisExtQ 7 < 3 := by
  rw [kurtosisExtQ_seven]; norm_num

/-- `kurtosisExtQ 8 < 3` — platykurtic at `N = 8`. -/
theorem kurtosisExtQ_eight_lt_three : kurtosisExtQ 8 < 3 := by
  rw [kurtosisExtQ_eight]; norm_num

/-- `kurtosisExtQ 9 < 3` — platykurtic at `N = 9`. -/
theorem kurtosisExtQ_nine_lt_three : kurtosisExtQ 9 < 3 := by
  rw [kurtosisExtQ_nine]; norm_num

/-- `kurtosisExtQ 10 < 3` — platykurtic at `N = 10`. -/
theorem kurtosisExtQ_ten_lt_three : kurtosisExtQ 10 < 3 := by
  rw [kurtosisExtQ_ten]; norm_num

/-! ### Headline summary -/

/-- **Headline (extended distributional kurtosis, `N ∈ {5, …, 10}`).**

    For the eleven-prime PT sequence `ptPrimeExt` with gaps
    `(1, 2, 2, 4, 2, 4, 2, 4, 6, 2)`, the six distributional fourth
    central moments

      `gapsDistCentralMoment4ExtQ N
          = (Σ (n − μ_N)⁴ · g_n) / (Σ g_n)`

    are exactly

      `70754/14641, 208034/16875, 1604890/83521,
       2540354/64827, 12445406/177147, 61313616/707281`.

    All six are **strictly positive**, and the corresponding
    dimensionless kurtosis ratios `K_N := m_4 / σ⁴` are all
    **strictly less than `3`** (the Gaussian reference): the gap-weighted
    index distribution is consistently **platykurtic** across the
    extended range. -/
theorem gapsDistCentralMoment4Ext_headline :
    -- six exact fourth-moment values
    (gapsDistCentralMoment4ExtQ 5  = 70754    / 14641)
    ∧ (gapsDistCentralMoment4ExtQ 6  = 208034   / 16875)
    ∧ (gapsDistCentralMoment4ExtQ 7  = 1604890  / 83521)
    ∧ (gapsDistCentralMoment4ExtQ 8  = 2540354  / 64827)
    ∧ (gapsDistCentralMoment4ExtQ 9  = 12445406 / 177147)
    ∧ (gapsDistCentralMoment4ExtQ 10 = 61313616 / 707281)
    -- strict positivity on the extended range
    ∧ (gapsDistCentralMoment4ExtQ 5  > 0)
    ∧ (gapsDistCentralMoment4ExtQ 6  > 0)
    ∧ (gapsDistCentralMoment4ExtQ 7  > 0)
    ∧ (gapsDistCentralMoment4ExtQ 8  > 0)
    ∧ (gapsDistCentralMoment4ExtQ 9  > 0)
    ∧ (gapsDistCentralMoment4ExtQ 10 > 0)
    -- six exact kurtosis ratios
    ∧ (kurtosisExtQ 5  = 35377    / 16562)
    ∧ (kurtosisExtQ 6  = 312051   / 153458)
    ∧ (kurtosisExtQ 7  = 802445   / 392498)
    ∧ (kurtosisExtQ 8  = 3810531  / 1940450)
    ∧ (kurtosisExtQ 9  = 18668109 / 9919058)
    ∧ (kurtosisExtQ 10 = 425789   / 222784)
    -- uniform platykurtosis across the extended range
    ∧ (kurtosisExtQ 5  < 3)
    ∧ (kurtosisExtQ 6  < 3)
    ∧ (kurtosisExtQ 7  < 3)
    ∧ (kurtosisExtQ 8  < 3)
    ∧ (kurtosisExtQ 9  < 3)
    ∧ (kurtosisExtQ 10 < 3) :=
  ⟨gapsDistCentralMoment4ExtQ_five,  gapsDistCentralMoment4ExtQ_six,
   gapsDistCentralMoment4ExtQ_seven, gapsDistCentralMoment4ExtQ_eight,
   gapsDistCentralMoment4ExtQ_nine,  gapsDistCentralMoment4ExtQ_ten,
   gapsDistCentralMoment4ExtQ_five_pos,   gapsDistCentralMoment4ExtQ_six_pos,
   gapsDistCentralMoment4ExtQ_seven_pos,  gapsDistCentralMoment4ExtQ_eight_pos,
   gapsDistCentralMoment4ExtQ_nine_pos,   gapsDistCentralMoment4ExtQ_ten_pos,
   kurtosisExtQ_five,  kurtosisExtQ_six,  kurtosisExtQ_seven,
   kurtosisExtQ_eight, kurtosisExtQ_nine, kurtosisExtQ_ten,
   kurtosisExtQ_five_lt_three,   kurtosisExtQ_six_lt_three,
   kurtosisExtQ_seven_lt_three,  kurtosisExtQ_eight_lt_three,
   kurtosisExtQ_nine_lt_three,   kurtosisExtQ_ten_lt_three⟩

end PT.Conservation
