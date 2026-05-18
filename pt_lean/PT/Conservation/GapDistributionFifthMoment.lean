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
import PT.Conservation.GapDistributionKurtosisExtended
import PT.Conservation.GapStatisticalMoments
import Mathlib.Algebra.BigOperators.Intervals
import Mathlib.Data.Rat.Defs
import Mathlib.Tactic

/-!
# Distributional fifth central moment of the PT gap distribution —
# extension `N ∈ {5, …, 10}`

This file is the **fifth-moment** companion to
`PT.Conservation.GapDistributionSkewExtended` (third central moment)
and `PT.Conservation.GapDistributionKurtosisExtended` (fourth central
moment). On the eleven-prime sequence `ptPrimeExt`, it computes the
**distributional fifth central moment** of the index `n` under the
probability distribution `n ↦ g_n / S_N`:

  `gapsDistCentralMoment5ExtQ N
      := (Σ_{n=1}^{N} (n − μ_N)⁵ · g_n) / S_N`,

where `μ_N := gapsDistMeanExtQ N` and `S_N := Σ_{n=1}^{N} g_n`.

Like the third central moment (an odd power), the **sign** of the fifth
central moment encodes the **asymmetry** of the index distribution:
negative values indicate a left tail (the centred index distribution
has more mass weighted to the left of the mean once amplified by the
fifth power), positive values indicate a right tail. Compared with the
third central moment, the fifth power **amplifies the tails** even
further, giving a more pronounced view of the same asymmetry.

With gaps `(g_1, …, g_{10}) = (1, 2, 2, 4, 2, 4, 2, 4, 6, 2)`, the six
exact rationals are

| `N`  | `S_N` | `μ_N`    | `gapsDistCentralMoment5ExtQ N`     |
|------|-------|----------|------------------------------------|
| `5`  | `11`  | `37/11`  | `-868710 / 161051`                 |
| `6`  | `15`  | `61/15`  | `-12146546 / 759375`               |
| `7`  | `17`  | `75/17`  | `-30390210 / 1419857`              |
| `8`  | `21`  | `107/21` | `-180349510 / 4084101`             |
| `9`  | `27`  | `161/27` | `-2044934590 / 14348907`           |
| `10` | `29`  | `181/29` | `-3622740720 / 20511149`           |

The denominators factor as `S_N⁵`: `161051 = 11⁵`, `759375 = 15⁵`,
`1419857 = 17⁵`, `4084101 = 21⁵`, `14348907 = 27⁵`, `20511149 = 29⁵`.

**All six values are strictly negative** — exactly like the third
central moment counterparts in `GapDistributionSkewExtended` — so the
**left skew** of the gap-weighted index distribution is uniformly
confirmed by the fifth-moment diagnostic on the extended range
`N ∈ {5, …, 10}`. The fifth-power amplification of tails does not
change the sign of the asymmetry.

All values are exact rationals, discharged by unfolding the `Finset.Ico`
sums step by step and closing with `norm_num`.

## Reference

* `PT.Conservation.CumulativeBoundsExtended` — `ptPrimeExt`, gaps,
  cumulative sums `S_N`.
* `PT.Conservation.GapDistributionMeanExtended` — `gapsDistMeanExtQ N`
  and the six exact mean values used here.
* `PT.Conservation.GapDistributionVarianceExtended` — the rational gap
  lemmas `gapExt_one_rat, …, gapExt_ten_rat` reused below.
* `PT.Conservation.GapDistributionSkewExtended` — the third central
  moment companion file (also strictly negative across `N ∈ {5,…,10}`).
* `PT.Conservation.GapDistributionKurtosisExtended` — the fourth
  central moment companion file (strictly positive, platykurtic).
-/

namespace PT.Conservation

open Finset

/-! ### Distributional fifth central moment on the extended range -/

/-- **Extended distributional fifth central moment.** Fifth central
    moment of the index `n` under the probability distribution
    `n ↦ g_n / S_N` on the eleven-prime sequence `ptPrimeExt`:

      `gapsDistCentralMoment5ExtQ N
          = (Σ_{n=1}^{N} (n − μ_N)⁵ · g_n) / S_N`,

    where `μ_N = gapsDistMeanExtQ N`. The sign of this quantity is the
    sign of the (fifth-power amplified) distributional skewness of the
    index. -/
def gapsDistCentralMoment5ExtQ (N : ℕ) : ℚ :=
  (∑ n ∈ Ico (1 : ℕ) (N + 1),
      ((n : ℚ) - gapsDistMeanExtQ N) ^ 5 * (gap ptPrimeExt n : ℚ))
    / cumGapExtQ N

/-! ### Six exact values `N ∈ {5, …, 10}`

Each fifth central moment is discharged by unfolding the `Ico` sum
step by step (reducing it to a sum of `N` explicit terms) and closing
with `norm_num`, after rewriting `μ_N = gapsDistMeanExtQ N` and
`S_N = cumGapExtQ N` to their exact rational values.
-/

/-- **Distributional fifth central moment at `N = 5`.**
    `-868710 / 161051` with `161051 = 11⁵`. -/
theorem gapsDistCentralMoment5ExtQ_five :
    gapsDistCentralMoment5ExtQ 5 = -868710 / 161051 := by
  unfold gapsDistCentralMoment5ExtQ
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

/-- **Distributional fifth central moment at `N = 6`.**
    `-12146546 / 759375` with `759375 = 15⁵`. -/
theorem gapsDistCentralMoment5ExtQ_six :
    gapsDistCentralMoment5ExtQ 6 = -12146546 / 759375 := by
  unfold gapsDistCentralMoment5ExtQ
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

/-- **Distributional fifth central moment at `N = 7`.**
    `-30390210 / 1419857` with `1419857 = 17⁵`. -/
theorem gapsDistCentralMoment5ExtQ_seven :
    gapsDistCentralMoment5ExtQ 7 = -30390210 / 1419857 := by
  unfold gapsDistCentralMoment5ExtQ
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

/-- **Distributional fifth central moment at `N = 8`.**
    `-180349510 / 4084101` with `4084101 = 21⁵`. -/
theorem gapsDistCentralMoment5ExtQ_eight :
    gapsDistCentralMoment5ExtQ 8 = -180349510 / 4084101 := by
  unfold gapsDistCentralMoment5ExtQ
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

/-- **Distributional fifth central moment at `N = 9`.**
    `-2044934590 / 14348907` with `14348907 = 27⁵`. -/
theorem gapsDistCentralMoment5ExtQ_nine :
    gapsDistCentralMoment5ExtQ 9 = -2044934590 / 14348907 := by
  unfold gapsDistCentralMoment5ExtQ
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

/-- **Distributional fifth central moment at `N = 10`.**
    `-3622740720 / 20511149` with `20511149 = 29⁵`. -/
theorem gapsDistCentralMoment5ExtQ_ten :
    gapsDistCentralMoment5ExtQ 10 = -3622740720 / 20511149 := by
  unfold gapsDistCentralMoment5ExtQ
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

/-! ### Strict negativity (uniform left skew on the extended range)

Each fifth central moment is strictly negative. Like the third central
moment counterpart (`GapDistributionSkewExtended`), this expresses a
**left skew** of the gap-weighted index distribution: the few large
gaps (in particular `g_9 = 6`) anchor the mean `μ_N` close to the right
boundary, so the residual asymmetry of the centred index distribution
points to the left tail — and the fifth-power amplification preserves
the sign. We discharge each case by rewriting to the explicit rational
value. -/

/-- `gapsDistCentralMoment5ExtQ 5 < 0` — left skew at `N = 5`. -/
theorem gapsDistCentralMoment5ExtQ_five_neg :
    gapsDistCentralMoment5ExtQ 5 < 0 := by
  rw [gapsDistCentralMoment5ExtQ_five]; norm_num

/-- `gapsDistCentralMoment5ExtQ 6 < 0` — left skew at `N = 6`. -/
theorem gapsDistCentralMoment5ExtQ_six_neg :
    gapsDistCentralMoment5ExtQ 6 < 0 := by
  rw [gapsDistCentralMoment5ExtQ_six]; norm_num

/-- `gapsDistCentralMoment5ExtQ 7 < 0` — left skew at `N = 7`. -/
theorem gapsDistCentralMoment5ExtQ_seven_neg :
    gapsDistCentralMoment5ExtQ 7 < 0 := by
  rw [gapsDistCentralMoment5ExtQ_seven]; norm_num

/-- `gapsDistCentralMoment5ExtQ 8 < 0` — left skew at `N = 8`. -/
theorem gapsDistCentralMoment5ExtQ_eight_neg :
    gapsDistCentralMoment5ExtQ 8 < 0 := by
  rw [gapsDistCentralMoment5ExtQ_eight]; norm_num

/-- `gapsDistCentralMoment5ExtQ 9 < 0` — left skew at `N = 9`. -/
theorem gapsDistCentralMoment5ExtQ_nine_neg :
    gapsDistCentralMoment5ExtQ 9 < 0 := by
  rw [gapsDistCentralMoment5ExtQ_nine]; norm_num

/-- `gapsDistCentralMoment5ExtQ 10 < 0` — left skew at `N = 10`. -/
theorem gapsDistCentralMoment5ExtQ_ten_neg :
    gapsDistCentralMoment5ExtQ 10 < 0 := by
  rw [gapsDistCentralMoment5ExtQ_ten]; norm_num

/-! ### Headline summary -/

/-- **Headline (extended distributional fifth central moment,
    `N ∈ {5, …, 10}`).**

    For the eleven-prime PT sequence `ptPrimeExt` with gaps
    `(1, 2, 2, 4, 2, 4, 2, 4, 6, 2)`, the six distributional fifth
    central moments

      `gapsDistCentralMoment5ExtQ N
          = (Σ (n − μ_N)⁵ · g_n) / (Σ g_n)`

    are exactly

      `-868710/161051,   -12146546/759375, -30390210/1419857,
       -180349510/4084101, -2044934590/14348907,
       -3622740720/20511149`,

    with denominators `S_N⁵` (`11⁵, 15⁵, 17⁵, 21⁵, 27⁵, 29⁵`).

    **All six are strictly negative**: the fifth-power amplification of
    the tails confirms the uniform **left skew** of the gap-weighted
    index distribution on the extended range — the same qualitative
    diagnosis as the third central moment in
    `GapDistributionSkewExtended`. -/
theorem gapsDistCentralMoment5Ext_headline :
    -- six exact fifth-moment values
    (gapsDistCentralMoment5ExtQ 5  = -868710     / 161051)
    ∧ (gapsDistCentralMoment5ExtQ 6  = -12146546   / 759375)
    ∧ (gapsDistCentralMoment5ExtQ 7  = -30390210   / 1419857)
    ∧ (gapsDistCentralMoment5ExtQ 8  = -180349510  / 4084101)
    ∧ (gapsDistCentralMoment5ExtQ 9  = -2044934590 / 14348907)
    ∧ (gapsDistCentralMoment5ExtQ 10 = -3622740720 / 20511149)
    -- strict negativity (uniform left skew) on the extended range
    ∧ (gapsDistCentralMoment5ExtQ 5  < 0)
    ∧ (gapsDistCentralMoment5ExtQ 6  < 0)
    ∧ (gapsDistCentralMoment5ExtQ 7  < 0)
    ∧ (gapsDistCentralMoment5ExtQ 8  < 0)
    ∧ (gapsDistCentralMoment5ExtQ 9  < 0)
    ∧ (gapsDistCentralMoment5ExtQ 10 < 0) :=
  ⟨gapsDistCentralMoment5ExtQ_five,  gapsDistCentralMoment5ExtQ_six,
   gapsDistCentralMoment5ExtQ_seven, gapsDistCentralMoment5ExtQ_eight,
   gapsDistCentralMoment5ExtQ_nine,  gapsDistCentralMoment5ExtQ_ten,
   gapsDistCentralMoment5ExtQ_five_neg,   gapsDistCentralMoment5ExtQ_six_neg,
   gapsDistCentralMoment5ExtQ_seven_neg,  gapsDistCentralMoment5ExtQ_eight_neg,
   gapsDistCentralMoment5ExtQ_nine_neg,   gapsDistCentralMoment5ExtQ_ten_neg⟩

end PT.Conservation
