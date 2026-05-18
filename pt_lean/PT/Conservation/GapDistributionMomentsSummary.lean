/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Conservation.GapDistributionVariance
import PT.Conservation.GapDistributionMeanExtended
import PT.Conservation.GapDistributionVarianceExtended
import PT.Conservation.GapDistributionSkewExtended
import PT.Conservation.GapDistributionKurtosisExtended
import Mathlib.Data.Rat.Defs
import Mathlib.Tactic

/-!
# Summary table — first four distributional moments of the PT gap distribution

This file is a **pure aggregation module**: it does not introduce a single
new computation. It reorganises the four families of moment values that
are scattered across

* `PT.Conservation.GapDistributionVariance`
  — `gapsDistMeanQ 4 = 3`, `gapsDistVarianceQ 4 = 10/9`
  (anchor case on the five-prime sequence `ptPrime`),
* `PT.Conservation.GapDistributionMeanExtended`
  — `gapsDistMeanExtQ N` for `N ∈ {5, …, 10}`,
* `PT.Conservation.GapDistributionVarianceExtended`
  — `gapsDistVarianceExtQ N` for `N ∈ {5, …, 10}`,
* `PT.Conservation.GapDistributionSkewExtended`
  — `gapsDistCentralMoment3ExtQ N` for `N ∈ {5, …, 10}`,
* `PT.Conservation.GapDistributionKurtosisExtended`
  — `gapsDistCentralMoment4ExtQ N` for `N ∈ {5, …, 10}`,

into a single conjunctive theorem that records all 24 extended values
plus the anchor pair at `N = 4`, packaged as four-tuples
`(mean, variance, m3, m4)` per `N`.

## Important — two different notions of "moment"

The values aggregated here are **distributional**: averages over the
**index set** `{1, …, N}` weighted by the probability law
`n ↦ g_n / S_N`. They are *not* the classical population moments of
the gap **values** `(g_1, …, g_N)`. The latter are exposed by
`PT.Conservation.GapStatisticalMoments` under the names
`centralMoment3 4 = 27/32`, `centralMoment4 4 = 757/256`, etc., and are
**genuinely different** (see e.g. `gapsDistVarianceQ_four = 10/9`
vs `varianceGapsTo 4 = 19/16`).

## Layout

* Table 1: the six four-tuples for `N ∈ {5, …, 10}` on `ptPrimeExt`.
* Table 2: the anchor pair `(mean, variance) = (3, 10/9)` at `N = 4`
  on `ptPrime` (the third and fourth distributional moments at `N = 4`
  are not separately defined in the corpus — only the **value-based**
  central moments are, in `GapStatisticalMoments`, and those are not
  the same object).
* `gapsDistMomentsExt_summary`: the 24-conjunct headline.
* `gapsDistMoments_full_summary`: the headline plus the `N = 4` anchor.
-/

namespace PT.Conservation

open Finset

/-! ### Per-`N` four-tuple of distributional moments -/

/-- **Four-tuple of distributional moments at `N`** on the extended
    sequence `ptPrimeExt`:

      `(mean, variance, third central moment, fourth central moment)`. -/
def gapsDistMomentsExt_tuple (N : ℕ) : ℚ × ℚ × ℚ × ℚ :=
  ( gapsDistMeanExtQ N,
    gapsDistVarianceExtQ N,
    gapsDistCentralMoment3ExtQ N,
    gapsDistCentralMoment4ExtQ N )

/-! ### Six exact four-tuples for `N ∈ {5, …, 10}` -/

/-- **Moments four-tuple at `N = 5`** — `(37/11, 182/121, -1038/1331, 70754/14641)`. -/
theorem gapsDistMomentsExt_tuple_five :
    gapsDistMomentsExt_tuple 5 =
      ( 37 / 11, 182 / 121, -1038 / 1331, 70754 / 14641 ) := by
  unfold gapsDistMomentsExt_tuple
  rw [gapsDistMeanExtQ_five, gapsDistVarianceExtQ_five,
      gapsDistCentralMoment3ExtQ_five, gapsDistCentralMoment4ExtQ_five]

/-- **Moments four-tuple at `N = 6`** — `(61/15, 554/225, -4138/3375, 208034/16875)`. -/
theorem gapsDistMomentsExt_tuple_six :
    gapsDistMomentsExt_tuple 6 =
      ( 61 / 15, 554 / 225, -4138 / 3375, 208034 / 16875 ) := by
  unfold gapsDistMomentsExt_tuple
  rw [gapsDistMeanExtQ_six, gapsDistVarianceExtQ_six,
      gapsDistCentralMoment3ExtQ_six, gapsDistCentralMoment4ExtQ_six]

/-- **Moments four-tuple at `N = 7`** — `(75/17, 886/289, -6522/4913, 1604890/83521)`. -/
theorem gapsDistMomentsExt_tuple_seven :
    gapsDistMomentsExt_tuple 7 =
      ( 75 / 17, 886 / 289, -6522 / 4913, 1604890 / 83521 ) := by
  unfold gapsDistMomentsExt_tuple
  rw [gapsDistMeanExtQ_seven, gapsDistVarianceExtQ_seven,
      gapsDistCentralMoment3ExtQ_seven, gapsDistCentralMoment4ExtQ_seven]

/-- **Moments four-tuple at `N = 8`** — `(107/21, 1970/441, -16238/9261, 2540354/64827)`. -/
theorem gapsDistMomentsExt_tuple_eight :
    gapsDistMomentsExt_tuple 8 =
      ( 107 / 21, 1970 / 441, -16238 / 9261, 2540354 / 64827 ) := by
  unfold gapsDistMomentsExt_tuple
  rw [gapsDistMeanExtQ_eight, gapsDistVarianceExtQ_eight,
      gapsDistCentralMoment3ExtQ_eight, gapsDistCentralMoment4ExtQ_eight]

/-- **Moments four-tuple at `N = 9`** — `(161/27, 4454/729, -92342/19683, 12445406/177147)`. -/
theorem gapsDistMomentsExt_tuple_nine :
    gapsDistMomentsExt_tuple 9 =
      ( 161 / 27, 4454 / 729, -92342 / 19683, 12445406 / 177147 ) := by
  unfold gapsDistMomentsExt_tuple
  rw [gapsDistMeanExtQ_nine, gapsDistVarianceExtQ_nine,
      gapsDistCentralMoment3ExtQ_nine, gapsDistCentralMoment4ExtQ_nine]

/-- **Moments four-tuple at `N = 10`** —
    `(181/29, 5664/841, -133584/24389, 61313616/707281)`. -/
theorem gapsDistMomentsExt_tuple_ten :
    gapsDistMomentsExt_tuple 10 =
      ( 181 / 29, 5664 / 841, -133584 / 24389, 61313616 / 707281 ) := by
  unfold gapsDistMomentsExt_tuple
  rw [gapsDistMeanExtQ_ten, gapsDistVarianceExtQ_ten,
      gapsDistCentralMoment3ExtQ_ten, gapsDistCentralMoment4ExtQ_ten]

/-! ### Aggregate headline — the 24 extended values, one conjunction -/

/-- **Aggregate headline (extended distributional moments).**

    All 24 distributional moment values for `N ∈ {5, …, 10}` on the
    eleven-prime sequence `ptPrimeExt`, organised by `N` as
    `(μ_N, σ²_N, m_{3,N}, m_{4,N})`:

    | `N`  | `μ`        | `σ²`       | `m_3`           | `m_4`              |
    |------|------------|------------|-----------------|--------------------|
    | `5`  | `37/11`    | `182/121`  | `-1038/1331`    | `70754/14641`      |
    | `6`  | `61/15`    | `554/225`  | `-4138/3375`    | `208034/16875`     |
    | `7`  | `75/17`    | `886/289`  | `-6522/4913`    | `1604890/83521`    |
    | `8`  | `107/21`   | `1970/441` | `-16238/9261`   | `2540354/64827`    |
    | `9`  | `161/27`   | `4454/729` | `-92342/19683`  | `12445406/177147`  |
    | `10` | `181/29`   | `5664/841` | `-133584/24389` | `61313616/707281`  |

    This is the single-statement view of the four moment tables exposed by
    `GapDistributionMeanExtended`, `GapDistributionVarianceExtended`,
    `GapDistributionSkewExtended`, and `GapDistributionKurtosisExtended`. -/
theorem gapsDistMomentsExt_summary :
    (gapsDistMomentsExt_tuple 5  =
      ( 37 / 11, 182 / 121, -1038 / 1331, 70754 / 14641 ))
    ∧ (gapsDistMomentsExt_tuple 6  =
      ( 61 / 15, 554 / 225, -4138 / 3375, 208034 / 16875 ))
    ∧ (gapsDistMomentsExt_tuple 7  =
      ( 75 / 17, 886 / 289, -6522 / 4913, 1604890 / 83521 ))
    ∧ (gapsDistMomentsExt_tuple 8  =
      ( 107 / 21, 1970 / 441, -16238 / 9261, 2540354 / 64827 ))
    ∧ (gapsDistMomentsExt_tuple 9  =
      ( 161 / 27, 4454 / 729, -92342 / 19683, 12445406 / 177147 ))
    ∧ (gapsDistMomentsExt_tuple 10 =
      ( 181 / 29, 5664 / 841, -133584 / 24389, 61313616 / 707281 )) :=
  ⟨gapsDistMomentsExt_tuple_five,  gapsDistMomentsExt_tuple_six,
   gapsDistMomentsExt_tuple_seven, gapsDistMomentsExt_tuple_eight,
   gapsDistMomentsExt_tuple_nine,  gapsDistMomentsExt_tuple_ten⟩

/-! ### Sign / shape qualitative summary

A compact restatement of what the 24 numbers above *say* about the
shape of the gap-weighted index distribution on the extended range. -/

/-- **Qualitative shape headline (extended range).**

    The gap-weighted index distribution `n ↦ g_n / S_N` on `ptPrimeExt`
    is, for every `N ∈ {5, …, 10}`:

    * **right-biased**: `μ_N > (N+1)/2`;
    * **dispersed**: `σ²_N > 0`;
    * **left-skewed in shape**: `m_{3,N} < 0`;
    * **non-degenerate in the fourth moment**: `m_{4,N} > 0`. -/
theorem gapsDistMomentsExt_shape_summary :
    -- right-bias of the mean
    (gapsDistMeanExtQ 5  > (5  + 1 : ℚ) / 2)
    ∧ (gapsDistMeanExtQ 6  > (6  + 1 : ℚ) / 2)
    ∧ (gapsDistMeanExtQ 7  > (7  + 1 : ℚ) / 2)
    ∧ (gapsDistMeanExtQ 8  > (8  + 1 : ℚ) / 2)
    ∧ (gapsDistMeanExtQ 9  > (9  + 1 : ℚ) / 2)
    ∧ (gapsDistMeanExtQ 10 > (10 + 1 : ℚ) / 2)
    -- left-skew (negative third central moment)
    ∧ (gapsDistCentralMoment3ExtQ 5  < 0)
    ∧ (gapsDistCentralMoment3ExtQ 6  < 0)
    ∧ (gapsDistCentralMoment3ExtQ 7  < 0)
    ∧ (gapsDistCentralMoment3ExtQ 8  < 0)
    ∧ (gapsDistCentralMoment3ExtQ 9  < 0)
    ∧ (gapsDistCentralMoment3ExtQ 10 < 0)
    -- strict positivity of the fourth central moment
    ∧ (gapsDistCentralMoment4ExtQ 5  > 0)
    ∧ (gapsDistCentralMoment4ExtQ 6  > 0)
    ∧ (gapsDistCentralMoment4ExtQ 7  > 0)
    ∧ (gapsDistCentralMoment4ExtQ 8  > 0)
    ∧ (gapsDistCentralMoment4ExtQ 9  > 0)
    ∧ (gapsDistCentralMoment4ExtQ 10 > 0) :=
  ⟨gapsDistMeanExtQ_five_gt_uniform,  gapsDistMeanExtQ_six_gt_uniform,
   gapsDistMeanExtQ_seven_gt_uniform, gapsDistMeanExtQ_eight_gt_uniform,
   gapsDistMeanExtQ_nine_gt_uniform,  gapsDistMeanExtQ_ten_gt_uniform,
   gapsDistCentralMoment3ExtQ_five_neg,   gapsDistCentralMoment3ExtQ_six_neg,
   gapsDistCentralMoment3ExtQ_seven_neg,  gapsDistCentralMoment3ExtQ_eight_neg,
   gapsDistCentralMoment3ExtQ_nine_neg,   gapsDistCentralMoment3ExtQ_ten_neg,
   gapsDistCentralMoment4ExtQ_five_pos,   gapsDistCentralMoment4ExtQ_six_pos,
   gapsDistCentralMoment4ExtQ_seven_pos,  gapsDistCentralMoment4ExtQ_eight_pos,
   gapsDistCentralMoment4ExtQ_nine_pos,   gapsDistCentralMoment4ExtQ_ten_pos⟩

/-! ### Anchor case `N = 4` (on `ptPrime`)

The third and fourth distributional central moments at `N = 4` are not
separately formalised in the corpus (only `gapsDistMeanQ 4` and
`gapsDistVarianceQ 4` are exposed in `GapDistributionVariance`). The
classical population central moments of the gap **values**
`centralMoment3 4 = 27/32` and `centralMoment4 4 = 757/256` from
`GapStatisticalMoments` are a **different** statistical object and
should not be confused with the distributional values here.

We therefore record at `N = 4` the two definite distributional values,
plus the contrast against the classical value-based variance. -/

/-- **Anchor pair at `N = 4`** — `(mean, variance) = (3, 10/9)`,
    together with the contrast against the classical value-based
    variance `varianceGapsTo 4 = 19/16`. -/
theorem gapsDistMoments_four_anchor :
    gapsDistMeanQ 4 = 3
    ∧ gapsDistVarianceQ 4 = 10 / 9
    ∧ varianceGapsTo 4 = 19 / 16
    ∧ gapsDistVarianceQ 4 < varianceGapsTo 4 :=
  gapsDist_mean_variance_four_summary

/-! ### Full headline — anchor `N = 4` + extended `N ∈ {5, …, 10}` -/

/-- **Full headline.** The complete distributional-moments summary:

    * **Anchor (`N = 4`, on `ptPrime`)**: `μ_4 = 3`, `σ²_4 = 10/9`,
      strictly below the classical population variance `19/16`.
    * **Extended table (`N ∈ {5, …, 10}`, on `ptPrimeExt`)**:
      the six four-tuples `(μ, σ², m_3, m_4)` listed in
      `gapsDistMomentsExt_summary`.

    This single statement aggregates the four moment tables that live in
    `GapDistributionVariance`, `GapDistribution{Mean,Variance,Skew,
    Kurtosis}Extended` into one citeable headline. -/
theorem gapsDistMoments_full_summary :
    -- anchor at N = 4 on ptPrime
    (gapsDistMeanQ 4 = 3)
    ∧ (gapsDistVarianceQ 4 = 10 / 9)
    ∧ (varianceGapsTo 4 = 19 / 16)
    ∧ (gapsDistVarianceQ 4 < varianceGapsTo 4)
    -- six extended four-tuples on ptPrimeExt
    ∧ (gapsDistMomentsExt_tuple 5  =
        ( 37 / 11, 182 / 121, -1038 / 1331, 70754 / 14641 ))
    ∧ (gapsDistMomentsExt_tuple 6  =
        ( 61 / 15, 554 / 225, -4138 / 3375, 208034 / 16875 ))
    ∧ (gapsDistMomentsExt_tuple 7  =
        ( 75 / 17, 886 / 289, -6522 / 4913, 1604890 / 83521 ))
    ∧ (gapsDistMomentsExt_tuple 8  =
        ( 107 / 21, 1970 / 441, -16238 / 9261, 2540354 / 64827 ))
    ∧ (gapsDistMomentsExt_tuple 9  =
        ( 161 / 27, 4454 / 729, -92342 / 19683, 12445406 / 177147 ))
    ∧ (gapsDistMomentsExt_tuple 10 =
        ( 181 / 29, 5664 / 841, -133584 / 24389, 61313616 / 707281 )) :=
  ⟨gapsDistMeanQ_four,
   gapsDistVarianceQ_four,
   varianceGapsTo_four,
   gapsDistVarianceQ_four_lt_varianceGapsTo_four,
   gapsDistMomentsExt_tuple_five,  gapsDistMomentsExt_tuple_six,
   gapsDistMomentsExt_tuple_seven, gapsDistMomentsExt_tuple_eight,
   gapsDistMomentsExt_tuple_nine,  gapsDistMomentsExt_tuple_ten⟩

end PT.Conservation
