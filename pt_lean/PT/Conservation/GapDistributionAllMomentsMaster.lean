/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Conservation.GapDistributionVariance
import PT.Conservation.GapDistributionMeanExtended
import PT.Conservation.GapDistributionVarianceExtended
import PT.Conservation.GapDistributionSkewExtended
import PT.Conservation.GapDistributionKurtosisExtended
import PT.Conservation.GapDistributionFifthMoment
import PT.Conservation.GapDistributionSixthMoment
import PT.Conservation.GapDistributionMomentsSummary
import PT.Conservation.GapStatisticalMoments
import Mathlib.Data.Rat.Defs
import Mathlib.Tactic

/-!
# Master aggregator — all six distributional central moments of the PT
# gap distribution on `N ∈ {4, 5, …, 10}`

This file is a **pure aggregator**. It introduces no new computation and
no new mathematical content: every conjunct below is discharged by an
already-proved lemma in one of the six per-moment files listed in the
imports. Its only purpose is to record, in a single citeable place, the
exhaustive distributional-moment table of the PT gap distribution
together with the **sign / parity** pattern observed on the extended
range.

## Six families, one table

For the PT prime sequence (small range `N = 4` on `ptPrime`, extended
range `N ∈ {5, …, 10}` on `ptPrimeExt`), the **distributional central
moments** of the index `n` under the probability law `n ↦ g_n / S_N`
are:

* **Moment 1** (mean, `gapsDistMeanQ` / `gapsDistMeanExtQ`) — strictly
  above the uniform-index mean `(N+1)/2` (right-bias).
* **Moment 2** (variance, `gapsDistVarianceQ` /
  `gapsDistVarianceExtQ`) — strictly positive and strictly increasing
  in `N` on the extended range.
* **Moment 3** (`gapsDistCentralMoment3ExtQ`) — uniformly **negative**
  (left skew).
* **Moment 4** (`gapsDistCentralMoment4ExtQ`) — uniformly **positive**
  (even power), platykurtic (`kurtosisExtQ N < 3`).
* **Moment 5** (`gapsDistCentralMoment5ExtQ`) — uniformly **negative**
  (left-skew amplification by an odd power).
* **Moment 6** (`gapsDistCentralMoment6ExtQ`) — uniformly **positive**
  (even power).

## The parity-of-moment-order pattern

On `N ∈ {5, …, 10}`:

* **Even-order central moments** (orders `2`, `4`, `6`) are **always
  strictly positive**.
* **Odd-order central moments** (orders `3`, `5`) are **always strictly
  negative**.

This is a direct **structural diagnosis** of the PT gap distribution —
the odd-order moments share a common sign (left skew) and the
even-order moments share the opposite sign (mass dispersion).

## Layout

* `master_moment_1` — moment 1 (means): right-bias on `{4, 5, …, 10}`.
* `master_moment_2` — moment 2 (variances): strict monotonicity on
  `{4, 5, …, 10}`.
* `master_moment_3` — moment 3 on `{5, …, 10}`: all negative.
* `master_moment_4` — moment 4 on `{5, …, 10}`: all positive,
  platykurtic.
* `master_moment_5` — moment 5 on `{5, …, 10}`: all negative.
* `master_moment_6` — moment 6 on `{5, …, 10}`: all positive.
* `master_moment_parity_pattern` — even / odd dichotomy on `{5, …, 10}`.
* `gap_distribution_all_moments_master` — final aggregator.
-/

namespace PT.Conservation

/-! ### Master moment 1 (means + right-bias) -/

/-- **Master moment 1** — distributional means.

    Records the seven exact mean values
    `gapsDistMeanQ 4 = 3` and `gapsDistMeanExtQ N` for `N ∈ {5, …, 10}`,
    together with the **right-bias witnesses** `μ_N > (N+1)/2` on the
    same range (the `N = 4` witness `3 > 5/2` is by `norm_num` on the
    exact value `gapsDistMeanQ 4 = 3`). -/
theorem master_moment_1 :
    -- seven exact means
    (gapsDistMeanQ    4  = 3)
    ∧ (gapsDistMeanExtQ 5  = 37  / 11)
    ∧ (gapsDistMeanExtQ 6  = 61  / 15)
    ∧ (gapsDistMeanExtQ 7  = 75  / 17)
    ∧ (gapsDistMeanExtQ 8  = 107 / 21)
    ∧ (gapsDistMeanExtQ 9  = 161 / 27)
    ∧ (gapsDistMeanExtQ 10 = 181 / 29)
    -- seven right-bias witnesses μ_N > (N+1)/2
    ∧ (gapsDistMeanQ    4  > (4  + 1 : ℚ) / 2)
    ∧ (gapsDistMeanExtQ 5  > (5  + 1 : ℚ) / 2)
    ∧ (gapsDistMeanExtQ 6  > (6  + 1 : ℚ) / 2)
    ∧ (gapsDistMeanExtQ 7  > (7  + 1 : ℚ) / 2)
    ∧ (gapsDistMeanExtQ 8  > (8  + 1 : ℚ) / 2)
    ∧ (gapsDistMeanExtQ 9  > (9  + 1 : ℚ) / 2)
    ∧ (gapsDistMeanExtQ 10 > (10 + 1 : ℚ) / 2) := by
  refine ⟨gapsDistMeanQ_four, gapsDistMeanExtQ_five, gapsDistMeanExtQ_six,
          gapsDistMeanExtQ_seven, gapsDistMeanExtQ_eight,
          gapsDistMeanExtQ_nine, gapsDistMeanExtQ_ten,
          ?_, gapsDistMeanExtQ_five_gt_uniform,
          gapsDistMeanExtQ_six_gt_uniform, gapsDistMeanExtQ_seven_gt_uniform,
          gapsDistMeanExtQ_eight_gt_uniform, gapsDistMeanExtQ_nine_gt_uniform,
          gapsDistMeanExtQ_ten_gt_uniform⟩
  rw [gapsDistMeanQ_four]; norm_num

/-! ### Master moment 2 (variances + strict monotonicity) -/

/-- **Master moment 2** — distributional variances.

    Records the seven exact variances `gapsDistVarianceQ 4 = 10/9` and
    `gapsDistVarianceExtQ N` for `N ∈ {5, …, 10}`, together with the
    **strict monotonicity** chain
    `σ²_4 < σ²_5 < σ²_6 < σ²_7 < σ²_8 < σ²_9 < σ²_{10}`. -/
theorem master_moment_2 :
    -- seven exact variances
    (gapsDistVarianceQ    4  = 10 / 9)
    ∧ (gapsDistVarianceExtQ 5  = 182  / 121)
    ∧ (gapsDistVarianceExtQ 6  = 554  / 225)
    ∧ (gapsDistVarianceExtQ 7  = 886  / 289)
    ∧ (gapsDistVarianceExtQ 8  = 1970 / 441)
    ∧ (gapsDistVarianceExtQ 9  = 4454 / 729)
    ∧ (gapsDistVarianceExtQ 10 = 5664 / 841)
    -- strict monotonicity from N = 4 through N = 10
    ∧ (gapsDistVarianceQ    4 < gapsDistVarianceExtQ 5)
    ∧ (gapsDistVarianceExtQ 5 < gapsDistVarianceExtQ 6)
    ∧ (gapsDistVarianceExtQ 6 < gapsDistVarianceExtQ 7)
    ∧ (gapsDistVarianceExtQ 7 < gapsDistVarianceExtQ 8)
    ∧ (gapsDistVarianceExtQ 8 < gapsDistVarianceExtQ 9)
    ∧ (gapsDistVarianceExtQ 9 < gapsDistVarianceExtQ 10) :=
  ⟨gapsDistVarianceQ_four, gapsDistVarianceExtQ_five,
   gapsDistVarianceExtQ_six, gapsDistVarianceExtQ_seven,
   gapsDistVarianceExtQ_eight, gapsDistVarianceExtQ_nine,
   gapsDistVarianceExtQ_ten,
   gapsDistVarianceQ_four_lt_gapsDistVarianceExtQ_five,
   gapsDistVarianceExtQ_five_lt_six, gapsDistVarianceExtQ_six_lt_seven,
   gapsDistVarianceExtQ_seven_lt_eight, gapsDistVarianceExtQ_eight_lt_nine,
   gapsDistVarianceExtQ_nine_lt_ten⟩

/-! ### Master moment 3 (third central moment, left skew) -/

/-- **Master moment 3** — distributional third central moments on
    `N ∈ {5, …, 10}`.

    Records the six exact values together with the uniform **left-skew**
    diagnosis `m_{3,N} < 0`. -/
theorem master_moment_3 :
    -- six exact third central moments
    (gapsDistCentralMoment3ExtQ 5  = -1038   / 1331)
    ∧ (gapsDistCentralMoment3ExtQ 6  = -4138   / 3375)
    ∧ (gapsDistCentralMoment3ExtQ 7  = -6522   / 4913)
    ∧ (gapsDistCentralMoment3ExtQ 8  = -16238  / 9261)
    ∧ (gapsDistCentralMoment3ExtQ 9  = -92342  / 19683)
    ∧ (gapsDistCentralMoment3ExtQ 10 = -133584 / 24389)
    -- uniform left skew on the extended range
    ∧ (gapsDistCentralMoment3ExtQ 5  < 0)
    ∧ (gapsDistCentralMoment3ExtQ 6  < 0)
    ∧ (gapsDistCentralMoment3ExtQ 7  < 0)
    ∧ (gapsDistCentralMoment3ExtQ 8  < 0)
    ∧ (gapsDistCentralMoment3ExtQ 9  < 0)
    ∧ (gapsDistCentralMoment3ExtQ 10 < 0) :=
  ⟨gapsDistCentralMoment3ExtQ_five,  gapsDistCentralMoment3ExtQ_six,
   gapsDistCentralMoment3ExtQ_seven, gapsDistCentralMoment3ExtQ_eight,
   gapsDistCentralMoment3ExtQ_nine,  gapsDistCentralMoment3ExtQ_ten,
   gapsDistCentralMoment3ExtQ_five_neg,  gapsDistCentralMoment3ExtQ_six_neg,
   gapsDistCentralMoment3ExtQ_seven_neg, gapsDistCentralMoment3ExtQ_eight_neg,
   gapsDistCentralMoment3ExtQ_nine_neg,  gapsDistCentralMoment3ExtQ_ten_neg⟩

/-! ### Master moment 4 (fourth central moment + platykurtic) -/

/-- **Master moment 4** — distributional fourth central moments and
    standardised kurtosis on `N ∈ {5, …, 10}`.

    Records the six exact fourth central moments, their strict
    positivity (even-power), and the **platykurtic** diagnosis
    `kurtosisExtQ N < 3` for every `N` in the range. -/
theorem master_moment_4 :
    -- six exact fourth central moments
    (gapsDistCentralMoment4ExtQ 5  = 70754    / 14641)
    ∧ (gapsDistCentralMoment4ExtQ 6  = 208034   / 16875)
    ∧ (gapsDistCentralMoment4ExtQ 7  = 1604890  / 83521)
    ∧ (gapsDistCentralMoment4ExtQ 8  = 2540354  / 64827)
    ∧ (gapsDistCentralMoment4ExtQ 9  = 12445406 / 177147)
    ∧ (gapsDistCentralMoment4ExtQ 10 = 61313616 / 707281)
    -- strict positivity (even power)
    ∧ (gapsDistCentralMoment4ExtQ 5  > 0)
    ∧ (gapsDistCentralMoment4ExtQ 6  > 0)
    ∧ (gapsDistCentralMoment4ExtQ 7  > 0)
    ∧ (gapsDistCentralMoment4ExtQ 8  > 0)
    ∧ (gapsDistCentralMoment4ExtQ 9  > 0)
    ∧ (gapsDistCentralMoment4ExtQ 10 > 0)
    -- platykurtic on the extended range
    ∧ (kurtosisExtQ 5  < 3)
    ∧ (kurtosisExtQ 6  < 3)
    ∧ (kurtosisExtQ 7  < 3)
    ∧ (kurtosisExtQ 8  < 3)
    ∧ (kurtosisExtQ 9  < 3)
    ∧ (kurtosisExtQ 10 < 3) :=
  ⟨gapsDistCentralMoment4ExtQ_five,  gapsDistCentralMoment4ExtQ_six,
   gapsDistCentralMoment4ExtQ_seven, gapsDistCentralMoment4ExtQ_eight,
   gapsDistCentralMoment4ExtQ_nine,  gapsDistCentralMoment4ExtQ_ten,
   gapsDistCentralMoment4ExtQ_five_pos,  gapsDistCentralMoment4ExtQ_six_pos,
   gapsDistCentralMoment4ExtQ_seven_pos, gapsDistCentralMoment4ExtQ_eight_pos,
   gapsDistCentralMoment4ExtQ_nine_pos,  gapsDistCentralMoment4ExtQ_ten_pos,
   kurtosisExtQ_five_lt_three,  kurtosisExtQ_six_lt_three,
   kurtosisExtQ_seven_lt_three, kurtosisExtQ_eight_lt_three,
   kurtosisExtQ_nine_lt_three,  kurtosisExtQ_ten_lt_three⟩

/-! ### Master moment 5 (fifth central moment, amplified left skew) -/

/-- **Master moment 5** — distributional fifth central moments on
    `N ∈ {5, …, 10}`.

    Records the six exact values together with their uniform strict
    negativity: the odd-power fifth moment preserves the left-skew sign
    seen at the third central moment. -/
theorem master_moment_5 :
    -- six exact fifth central moments
    (gapsDistCentralMoment5ExtQ 5  = -868710     / 161051)
    ∧ (gapsDistCentralMoment5ExtQ 6  = -12146546   / 759375)
    ∧ (gapsDistCentralMoment5ExtQ 7  = -30390210   / 1419857)
    ∧ (gapsDistCentralMoment5ExtQ 8  = -180349510  / 4084101)
    ∧ (gapsDistCentralMoment5ExtQ 9  = -2044934590 / 14348907)
    ∧ (gapsDistCentralMoment5ExtQ 10 = -3622740720 / 20511149)
    -- uniform strict negativity (odd-power amplified left skew)
    ∧ (gapsDistCentralMoment5ExtQ 5  < 0)
    ∧ (gapsDistCentralMoment5ExtQ 6  < 0)
    ∧ (gapsDistCentralMoment5ExtQ 7  < 0)
    ∧ (gapsDistCentralMoment5ExtQ 8  < 0)
    ∧ (gapsDistCentralMoment5ExtQ 9  < 0)
    ∧ (gapsDistCentralMoment5ExtQ 10 < 0) :=
  ⟨gapsDistCentralMoment5ExtQ_five,  gapsDistCentralMoment5ExtQ_six,
   gapsDistCentralMoment5ExtQ_seven, gapsDistCentralMoment5ExtQ_eight,
   gapsDistCentralMoment5ExtQ_nine,  gapsDistCentralMoment5ExtQ_ten,
   gapsDistCentralMoment5ExtQ_five_neg,  gapsDistCentralMoment5ExtQ_six_neg,
   gapsDistCentralMoment5ExtQ_seven_neg, gapsDistCentralMoment5ExtQ_eight_neg,
   gapsDistCentralMoment5ExtQ_nine_neg,  gapsDistCentralMoment5ExtQ_ten_neg⟩

/-! ### Master moment 6 (sixth central moment, positive tail mass) -/

/-- **Master moment 6** — distributional sixth central moments on
    `N ∈ {5, …, 10}`.

    Records the six exact values together with their uniform strict
    positivity (even sixth power). -/
theorem master_moment_6 :
    -- six exact sixth central moments
    (gapsDistCentralMoment6ExtQ 5  = 36381842     / 1771561)
    ∧ (gapsDistCentralMoment6ExtQ 6  = 182362814    / 2278125)
    ∧ (gapsDistCentralMoment6ExtQ 7  = 3765640426   / 24137569)
    ∧ (gapsDistCentralMoment6ExtQ 8  = 37379593750  / 85766121)
    ∧ (gapsDistCentralMoment6ExtQ 9  = 420003215194 / 387420489)
    ∧ (gapsDistCentralMoment6ExtQ 10 = 894548024304 / 594823321)
    -- uniform strict positivity (even-power)
    ∧ (gapsDistCentralMoment6ExtQ 5  > 0)
    ∧ (gapsDistCentralMoment6ExtQ 6  > 0)
    ∧ (gapsDistCentralMoment6ExtQ 7  > 0)
    ∧ (gapsDistCentralMoment6ExtQ 8  > 0)
    ∧ (gapsDistCentralMoment6ExtQ 9  > 0)
    ∧ (gapsDistCentralMoment6ExtQ 10 > 0) :=
  ⟨gapsDistCentralMoment6ExtQ_five,  gapsDistCentralMoment6ExtQ_six,
   gapsDistCentralMoment6ExtQ_seven, gapsDistCentralMoment6ExtQ_eight,
   gapsDistCentralMoment6ExtQ_nine,  gapsDistCentralMoment6ExtQ_ten,
   gapsDistCentralMoment6ExtQ_five_pos,  gapsDistCentralMoment6ExtQ_six_pos,
   gapsDistCentralMoment6ExtQ_seven_pos, gapsDistCentralMoment6ExtQ_eight_pos,
   gapsDistCentralMoment6ExtQ_nine_pos,  gapsDistCentralMoment6ExtQ_ten_pos⟩

/-! ### Master parity-of-order dichotomy -/

/-- **Master parity-of-moment-order pattern.**

    On the extended range `N ∈ {5, …, 10}`, the distributional central
    moments of the PT gap-weighted index distribution split by parity
    of the moment order:

    * **Even orders (`2`, `4`, `6`)** are uniformly **strictly positive**.
    * **Odd orders (`3`, `5`)** are uniformly **strictly negative**
      (uniform left skew across orders).

    This is a `5 × 6 = 30` conjunction of sign witnesses, all already
    proved in the per-moment files. -/
theorem master_moment_parity_pattern :
    -- even-order moments: order 2 (variances) > 0
    (gapsDistVarianceExtQ 5  > 0)
    ∧ (gapsDistVarianceExtQ 6  > 0)
    ∧ (gapsDistVarianceExtQ 7  > 0)
    ∧ (gapsDistVarianceExtQ 8  > 0)
    ∧ (gapsDistVarianceExtQ 9  > 0)
    ∧ (gapsDistVarianceExtQ 10 > 0)
    -- even-order moments: order 4 > 0
    ∧ (gapsDistCentralMoment4ExtQ 5  > 0)
    ∧ (gapsDistCentralMoment4ExtQ 6  > 0)
    ∧ (gapsDistCentralMoment4ExtQ 7  > 0)
    ∧ (gapsDistCentralMoment4ExtQ 8  > 0)
    ∧ (gapsDistCentralMoment4ExtQ 9  > 0)
    ∧ (gapsDistCentralMoment4ExtQ 10 > 0)
    -- even-order moments: order 6 > 0
    ∧ (gapsDistCentralMoment6ExtQ 5  > 0)
    ∧ (gapsDistCentralMoment6ExtQ 6  > 0)
    ∧ (gapsDistCentralMoment6ExtQ 7  > 0)
    ∧ (gapsDistCentralMoment6ExtQ 8  > 0)
    ∧ (gapsDistCentralMoment6ExtQ 9  > 0)
    ∧ (gapsDistCentralMoment6ExtQ 10 > 0)
    -- odd-order moments: order 3 < 0
    ∧ (gapsDistCentralMoment3ExtQ 5  < 0)
    ∧ (gapsDistCentralMoment3ExtQ 6  < 0)
    ∧ (gapsDistCentralMoment3ExtQ 7  < 0)
    ∧ (gapsDistCentralMoment3ExtQ 8  < 0)
    ∧ (gapsDistCentralMoment3ExtQ 9  < 0)
    ∧ (gapsDistCentralMoment3ExtQ 10 < 0)
    -- odd-order moments: order 5 < 0
    ∧ (gapsDistCentralMoment5ExtQ 5  < 0)
    ∧ (gapsDistCentralMoment5ExtQ 6  < 0)
    ∧ (gapsDistCentralMoment5ExtQ 7  < 0)
    ∧ (gapsDistCentralMoment5ExtQ 8  < 0)
    ∧ (gapsDistCentralMoment5ExtQ 9  < 0)
    ∧ (gapsDistCentralMoment5ExtQ 10 < 0) := by
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_,
          gapsDistCentralMoment4ExtQ_five_pos,
          gapsDistCentralMoment4ExtQ_six_pos,
          gapsDistCentralMoment4ExtQ_seven_pos,
          gapsDistCentralMoment4ExtQ_eight_pos,
          gapsDistCentralMoment4ExtQ_nine_pos,
          gapsDistCentralMoment4ExtQ_ten_pos,
          gapsDistCentralMoment6ExtQ_five_pos,
          gapsDistCentralMoment6ExtQ_six_pos,
          gapsDistCentralMoment6ExtQ_seven_pos,
          gapsDistCentralMoment6ExtQ_eight_pos,
          gapsDistCentralMoment6ExtQ_nine_pos,
          gapsDistCentralMoment6ExtQ_ten_pos,
          gapsDistCentralMoment3ExtQ_five_neg,
          gapsDistCentralMoment3ExtQ_six_neg,
          gapsDistCentralMoment3ExtQ_seven_neg,
          gapsDistCentralMoment3ExtQ_eight_neg,
          gapsDistCentralMoment3ExtQ_nine_neg,
          gapsDistCentralMoment3ExtQ_ten_neg,
          gapsDistCentralMoment5ExtQ_five_neg,
          gapsDistCentralMoment5ExtQ_six_neg,
          gapsDistCentralMoment5ExtQ_seven_neg,
          gapsDistCentralMoment5ExtQ_eight_neg,
          gapsDistCentralMoment5ExtQ_nine_neg,
          gapsDistCentralMoment5ExtQ_ten_neg⟩
  · rw [gapsDistVarianceExtQ_five];  norm_num
  · rw [gapsDistVarianceExtQ_six];   norm_num
  · rw [gapsDistVarianceExtQ_seven]; norm_num
  · rw [gapsDistVarianceExtQ_eight]; norm_num
  · rw [gapsDistVarianceExtQ_nine];  norm_num
  · rw [gapsDistVarianceExtQ_ten];   norm_num

/-! ### Final master aggregator -/

/-- **Master final theorem — all six distributional central moments of
    the PT gap distribution on `{4, 5, …, 10}`.**

    Single citeable headline assembling the six per-moment masters and
    the parity-of-order dichotomy:

    * `master_moment_1` — means on `{4, 5, …, 10}` with right-bias
      `μ_N > (N+1)/2`.
    * `master_moment_2` — variances on `{4, 5, …, 10}` with strict
      monotonicity.
    * `master_moment_3` — third central moments on `{5, …, 10}`,
      uniformly negative (left skew).
    * `master_moment_4` — fourth central moments on `{5, …, 10}`,
      uniformly positive, platykurtic.
    * `master_moment_5` — fifth central moments on `{5, …, 10}`,
      uniformly negative (amplified left skew).
    * `master_moment_6` — sixth central moments on `{5, …, 10}`,
      uniformly positive.
    * `master_moment_parity_pattern` — even orders positive, odd orders
      negative.

    This is an aggregator: every conjunct is already a theorem of one of
    the per-moment files. The point is to expose the complete diagnostic
    table of the PT gap distribution as one citeable object. -/
theorem gap_distribution_all_moments_master :
    -- moment 1: means + right-bias
    ((gapsDistMeanQ    4  = 3)
      ∧ (gapsDistMeanExtQ 5  = 37  / 11)
      ∧ (gapsDistMeanExtQ 6  = 61  / 15)
      ∧ (gapsDistMeanExtQ 7  = 75  / 17)
      ∧ (gapsDistMeanExtQ 8  = 107 / 21)
      ∧ (gapsDistMeanExtQ 9  = 161 / 27)
      ∧ (gapsDistMeanExtQ 10 = 181 / 29)
      ∧ (gapsDistMeanQ    4  > (4  + 1 : ℚ) / 2)
      ∧ (gapsDistMeanExtQ 5  > (5  + 1 : ℚ) / 2)
      ∧ (gapsDistMeanExtQ 6  > (6  + 1 : ℚ) / 2)
      ∧ (gapsDistMeanExtQ 7  > (7  + 1 : ℚ) / 2)
      ∧ (gapsDistMeanExtQ 8  > (8  + 1 : ℚ) / 2)
      ∧ (gapsDistMeanExtQ 9  > (9  + 1 : ℚ) / 2)
      ∧ (gapsDistMeanExtQ 10 > (10 + 1 : ℚ) / 2))
    -- moment 2: variances + monotonicity
    ∧ ((gapsDistVarianceQ    4  = 10 / 9)
      ∧ (gapsDistVarianceExtQ 5  = 182  / 121)
      ∧ (gapsDistVarianceExtQ 6  = 554  / 225)
      ∧ (gapsDistVarianceExtQ 7  = 886  / 289)
      ∧ (gapsDistVarianceExtQ 8  = 1970 / 441)
      ∧ (gapsDistVarianceExtQ 9  = 4454 / 729)
      ∧ (gapsDistVarianceExtQ 10 = 5664 / 841)
      ∧ (gapsDistVarianceQ    4 < gapsDistVarianceExtQ 5)
      ∧ (gapsDistVarianceExtQ 5 < gapsDistVarianceExtQ 6)
      ∧ (gapsDistVarianceExtQ 6 < gapsDistVarianceExtQ 7)
      ∧ (gapsDistVarianceExtQ 7 < gapsDistVarianceExtQ 8)
      ∧ (gapsDistVarianceExtQ 8 < gapsDistVarianceExtQ 9)
      ∧ (gapsDistVarianceExtQ 9 < gapsDistVarianceExtQ 10))
    -- moment 3
    ∧ ((gapsDistCentralMoment3ExtQ 5  = -1038   / 1331)
      ∧ (gapsDistCentralMoment3ExtQ 6  = -4138   / 3375)
      ∧ (gapsDistCentralMoment3ExtQ 7  = -6522   / 4913)
      ∧ (gapsDistCentralMoment3ExtQ 8  = -16238  / 9261)
      ∧ (gapsDistCentralMoment3ExtQ 9  = -92342  / 19683)
      ∧ (gapsDistCentralMoment3ExtQ 10 = -133584 / 24389)
      ∧ (gapsDistCentralMoment3ExtQ 5  < 0)
      ∧ (gapsDistCentralMoment3ExtQ 6  < 0)
      ∧ (gapsDistCentralMoment3ExtQ 7  < 0)
      ∧ (gapsDistCentralMoment3ExtQ 8  < 0)
      ∧ (gapsDistCentralMoment3ExtQ 9  < 0)
      ∧ (gapsDistCentralMoment3ExtQ 10 < 0))
    -- moment 4 + platykurtic
    ∧ ((gapsDistCentralMoment4ExtQ 5  = 70754    / 14641)
      ∧ (gapsDistCentralMoment4ExtQ 6  = 208034   / 16875)
      ∧ (gapsDistCentralMoment4ExtQ 7  = 1604890  / 83521)
      ∧ (gapsDistCentralMoment4ExtQ 8  = 2540354  / 64827)
      ∧ (gapsDistCentralMoment4ExtQ 9  = 12445406 / 177147)
      ∧ (gapsDistCentralMoment4ExtQ 10 = 61313616 / 707281)
      ∧ (gapsDistCentralMoment4ExtQ 5  > 0)
      ∧ (gapsDistCentralMoment4ExtQ 6  > 0)
      ∧ (gapsDistCentralMoment4ExtQ 7  > 0)
      ∧ (gapsDistCentralMoment4ExtQ 8  > 0)
      ∧ (gapsDistCentralMoment4ExtQ 9  > 0)
      ∧ (gapsDistCentralMoment4ExtQ 10 > 0)
      ∧ (kurtosisExtQ 5  < 3)
      ∧ (kurtosisExtQ 6  < 3)
      ∧ (kurtosisExtQ 7  < 3)
      ∧ (kurtosisExtQ 8  < 3)
      ∧ (kurtosisExtQ 9  < 3)
      ∧ (kurtosisExtQ 10 < 3))
    -- moment 5
    ∧ ((gapsDistCentralMoment5ExtQ 5  = -868710     / 161051)
      ∧ (gapsDistCentralMoment5ExtQ 6  = -12146546   / 759375)
      ∧ (gapsDistCentralMoment5ExtQ 7  = -30390210   / 1419857)
      ∧ (gapsDistCentralMoment5ExtQ 8  = -180349510  / 4084101)
      ∧ (gapsDistCentralMoment5ExtQ 9  = -2044934590 / 14348907)
      ∧ (gapsDistCentralMoment5ExtQ 10 = -3622740720 / 20511149)
      ∧ (gapsDistCentralMoment5ExtQ 5  < 0)
      ∧ (gapsDistCentralMoment5ExtQ 6  < 0)
      ∧ (gapsDistCentralMoment5ExtQ 7  < 0)
      ∧ (gapsDistCentralMoment5ExtQ 8  < 0)
      ∧ (gapsDistCentralMoment5ExtQ 9  < 0)
      ∧ (gapsDistCentralMoment5ExtQ 10 < 0))
    -- moment 6
    ∧ ((gapsDistCentralMoment6ExtQ 5  = 36381842     / 1771561)
      ∧ (gapsDistCentralMoment6ExtQ 6  = 182362814    / 2278125)
      ∧ (gapsDistCentralMoment6ExtQ 7  = 3765640426   / 24137569)
      ∧ (gapsDistCentralMoment6ExtQ 8  = 37379593750  / 85766121)
      ∧ (gapsDistCentralMoment6ExtQ 9  = 420003215194 / 387420489)
      ∧ (gapsDistCentralMoment6ExtQ 10 = 894548024304 / 594823321)
      ∧ (gapsDistCentralMoment6ExtQ 5  > 0)
      ∧ (gapsDistCentralMoment6ExtQ 6  > 0)
      ∧ (gapsDistCentralMoment6ExtQ 7  > 0)
      ∧ (gapsDistCentralMoment6ExtQ 8  > 0)
      ∧ (gapsDistCentralMoment6ExtQ 9  > 0)
      ∧ (gapsDistCentralMoment6ExtQ 10 > 0)) :=
  ⟨master_moment_1, master_moment_2, master_moment_3,
   master_moment_4, master_moment_5, master_moment_6⟩

end PT.Conservation
