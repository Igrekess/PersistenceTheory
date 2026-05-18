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
import PT.Conservation.GapStatisticalMoments
import Mathlib.Algebra.BigOperators.Intervals
import Mathlib.Data.Rat.Defs
import Mathlib.Tactic

/-!
# Distributional skewness (third central moment) of the PT gap
# distribution — extension `N ∈ {5, …, 10}`

`PT.Conservation.GapStatisticalMoments` introduces, on the small-`N`
sequence `ptPrime` (first five primes), the **unweighted third central
moment**

  `centralMoment3 N := (Σ_{n=1}^{N} (g_n − meanGapsTo N)³) / N`,

i.e. the scaled skewness of the *gap values* `g_n` viewed as a sample
of size `N` under the **counting** measure.

This file is the **distributional** analogue, on the extended sequence
`ptPrimeExt`: instead of the counting measure on `{1, …, N}`, we put on
each index `n` the probability `g_n / S_N`, and compute the third
central moment of the **index** under that distribution:

  `gapsDistCentralMoment3ExtQ N
      := (Σ_{n=1}^{N} (n − μ_N)³ · g_n) / S_N`,

where `μ_N := gapsDistMeanExtQ N` and `S_N := Σ_{n=1}^{N} g_n`. The
sign of this quantity classifies the asymmetry of the index
distribution: negative means a **left tail** (mass weighted towards
large `n` and skew towards small `n`), positive means a right tail.

With gaps `(g_1, …, g_{10}) = (1, 2, 2, 4, 2, 4, 2, 4, 6, 2)`, the
six exact rationals are

| `N`  | `S_N` | `μ_N`    | `gapsDistCentralMoment3ExtQ N` |
|------|-------|----------|--------------------------------|
| `5`  | `11`  | `37/11`  | `-1038 / 1331`                 |
| `6`  | `15`  | `61/15`  | `-4138 / 3375`                 |
| `7`  | `17`  | `75/17`  | `-6522 / 4913`                 |
| `8`  | `21`  | `107/21` | `-16238 / 9261`                |
| `9`  | `27`  | `161/27` | `-92342 / 19683`               |
| `10` | `29`  | `181/29` | `-133584 / 24389`              |

**All six values are strictly negative**: under the gap-weighted
distribution on the index range `{1, …, N}`, the centre of mass
`μ_N > (N+1)/2` (see `GapDistributionMeanExtended`), and the residual
asymmetry of the index distribution is a **left skew** — small indices
sit further from the mean than large ones because the few large gaps
(in particular `g_9 = 6`) anchor the mean very close to the right
boundary, so the residual asymmetry of the centred index distribution
points to the left tail.

All values are exact rationals, discharged by unfolding the `Finset.Ico`
sums step by step and closing with `norm_num`.

## Reference

* `PT.Conservation.CumulativeBoundsExtended` — `ptPrimeExt`, gaps,
  cumulative sums `S_N`.
* `PT.Conservation.GapDistributionMeanExtended` — `gapsDistMeanExtQ N`
  and the six exact mean values used here.
* `PT.Conservation.GapDistributionVarianceExtended` — the second central
  moment counterpart and the rational gap lemmas
  `gapExt_one_rat, …, gapExt_ten_rat` reused below.
* `PT.Conservation.GapStatisticalMoments` — `centralMoment3` on
  `ptPrime` (small-`N`, counting measure).
-/

namespace PT.Conservation

open Finset

/-! ### Distributional third central moment on the extended range -/

/-- **Extended distributional third central moment.** Third central
    moment of the index `n` under the probability distribution
    `n ↦ g_n / S_N` on the eleven-prime sequence `ptPrimeExt`:

      `gapsDistCentralMoment3ExtQ N
          = (Σ_{n=1}^{N} (n − μ_N)³ · g_n) / S_N`,

    where `μ_N = gapsDistMeanExtQ N`. The sign of this quantity is the
    sign of the distributional skewness of the index. -/
def gapsDistCentralMoment3ExtQ (N : ℕ) : ℚ :=
  (∑ n ∈ Ico (1 : ℕ) (N + 1),
      ((n : ℚ) - gapsDistMeanExtQ N) ^ 3 * (gap ptPrimeExt n : ℚ))
    / cumGapExtQ N

/-! ### Six exact values `N ∈ {5, …, 10}`

Each skewness is discharged by unfolding the `Ico` sum step by step
(reducing it to a sum of `N` explicit terms) and closing with `norm_num`,
after rewriting `μ_N = gapsDistMeanExtQ N` and `S_N = cumGapExtQ N` to
their exact rational values.
-/

/-- **Distributional skewness at `N = 5`.** `-1038 / 1331 = -1038 / 11³`. -/
theorem gapsDistCentralMoment3ExtQ_five :
    gapsDistCentralMoment3ExtQ 5 = -1038 / 1331 := by
  unfold gapsDistCentralMoment3ExtQ
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

/-- **Distributional skewness at `N = 6`.** `-4138 / 3375 = -4138 / 15³`. -/
theorem gapsDistCentralMoment3ExtQ_six :
    gapsDistCentralMoment3ExtQ 6 = -4138 / 3375 := by
  unfold gapsDistCentralMoment3ExtQ
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

/-- **Distributional skewness at `N = 7`.** `-6522 / 4913 = -6522 / 17³`. -/
theorem gapsDistCentralMoment3ExtQ_seven :
    gapsDistCentralMoment3ExtQ 7 = -6522 / 4913 := by
  unfold gapsDistCentralMoment3ExtQ
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

/-- **Distributional skewness at `N = 8`.** `-16238 / 9261 = -16238 / 21³`. -/
theorem gapsDistCentralMoment3ExtQ_eight :
    gapsDistCentralMoment3ExtQ 8 = -16238 / 9261 := by
  unfold gapsDistCentralMoment3ExtQ
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

/-- **Distributional skewness at `N = 9`.** `-92342 / 19683 = -92342 / 27³`. -/
theorem gapsDistCentralMoment3ExtQ_nine :
    gapsDistCentralMoment3ExtQ 9 = -92342 / 19683 := by
  unfold gapsDistCentralMoment3ExtQ
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

/-- **Distributional skewness at `N = 10`.** `-133584 / 24389 = -133584 / 29³`. -/
theorem gapsDistCentralMoment3ExtQ_ten :
    gapsDistCentralMoment3ExtQ 10 = -133584 / 24389 := by
  unfold gapsDistCentralMoment3ExtQ
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

/-! ### Sign of the distributional third central moment

All six values on the extended range are **strictly negative**: the
index distribution `n ↦ g_n / S_N` on `{1, …, N}` is **left-skewed**
once the heavy gap `g_9 = 6` (and its predecessors at index `4, 6, 8`)
pull the mean towards the right of the index range. The residual
asymmetry of the centred distribution then lives on the left tail.

This is the **opposite** sign from the small-`N` unweighted
`centralMoment3 4 = 27/32 > 0`: the gap-weighted distribution flips
the sign of the residual asymmetry compared to the counting measure on
gap *values*. -/

/-- `gapsDistCentralMoment3ExtQ 5 < 0` — left-skew at `N = 5`. -/
theorem gapsDistCentralMoment3ExtQ_five_neg :
    gapsDistCentralMoment3ExtQ 5 < 0 := by
  rw [gapsDistCentralMoment3ExtQ_five]; norm_num

/-- `gapsDistCentralMoment3ExtQ 6 < 0` — left-skew at `N = 6`. -/
theorem gapsDistCentralMoment3ExtQ_six_neg :
    gapsDistCentralMoment3ExtQ 6 < 0 := by
  rw [gapsDistCentralMoment3ExtQ_six]; norm_num

/-- `gapsDistCentralMoment3ExtQ 7 < 0` — left-skew at `N = 7`. -/
theorem gapsDistCentralMoment3ExtQ_seven_neg :
    gapsDistCentralMoment3ExtQ 7 < 0 := by
  rw [gapsDistCentralMoment3ExtQ_seven]; norm_num

/-- `gapsDistCentralMoment3ExtQ 8 < 0` — left-skew at `N = 8`. -/
theorem gapsDistCentralMoment3ExtQ_eight_neg :
    gapsDistCentralMoment3ExtQ 8 < 0 := by
  rw [gapsDistCentralMoment3ExtQ_eight]; norm_num

/-- `gapsDistCentralMoment3ExtQ 9 < 0` — left-skew at `N = 9`. -/
theorem gapsDistCentralMoment3ExtQ_nine_neg :
    gapsDistCentralMoment3ExtQ 9 < 0 := by
  rw [gapsDistCentralMoment3ExtQ_nine]; norm_num

/-- `gapsDistCentralMoment3ExtQ 10 < 0` — left-skew at `N = 10`. -/
theorem gapsDistCentralMoment3ExtQ_ten_neg :
    gapsDistCentralMoment3ExtQ 10 < 0 := by
  rw [gapsDistCentralMoment3ExtQ_ten]; norm_num

/-! ### Headline summary -/

/-- **Headline (extended distributional skewness, `N ∈ {5, …, 10}`).**

    For the eleven-prime PT sequence `ptPrimeExt` with gaps
    `(1, 2, 2, 4, 2, 4, 2, 4, 6, 2)`, the six distributional third
    central moments

      `gapsDistCentralMoment3ExtQ N
          = (Σ (n − μ_N)³ · g_n) / (Σ g_n)`

    are exactly

      `-1038/1331, -4138/3375, -6522/4913,
       -16238/9261, -92342/19683, -133584/24389`.

    All six are **strictly negative**: the gap-weighted index
    distribution on `{1, …, N}` is consistently left-skewed across the
    extended range. -/
theorem gapsDistCentralMoment3Ext_headline :
    -- six exact values
    (gapsDistCentralMoment3ExtQ 5  = -1038   / 1331)
    ∧ (gapsDistCentralMoment3ExtQ 6  = -4138   / 3375)
    ∧ (gapsDistCentralMoment3ExtQ 7  = -6522   / 4913)
    ∧ (gapsDistCentralMoment3ExtQ 8  = -16238  / 9261)
    ∧ (gapsDistCentralMoment3ExtQ 9  = -92342  / 19683)
    ∧ (gapsDistCentralMoment3ExtQ 10 = -133584 / 24389)
    -- uniform left-skew on the extended range
    ∧ (gapsDistCentralMoment3ExtQ 5  < 0)
    ∧ (gapsDistCentralMoment3ExtQ 6  < 0)
    ∧ (gapsDistCentralMoment3ExtQ 7  < 0)
    ∧ (gapsDistCentralMoment3ExtQ 8  < 0)
    ∧ (gapsDistCentralMoment3ExtQ 9  < 0)
    ∧ (gapsDistCentralMoment3ExtQ 10 < 0) :=
  ⟨gapsDistCentralMoment3ExtQ_five,  gapsDistCentralMoment3ExtQ_six,
   gapsDistCentralMoment3ExtQ_seven, gapsDistCentralMoment3ExtQ_eight,
   gapsDistCentralMoment3ExtQ_nine,  gapsDistCentralMoment3ExtQ_ten,
   gapsDistCentralMoment3ExtQ_five_neg,   gapsDistCentralMoment3ExtQ_six_neg,
   gapsDistCentralMoment3ExtQ_seven_neg,  gapsDistCentralMoment3ExtQ_eight_neg,
   gapsDistCentralMoment3ExtQ_nine_neg,   gapsDistCentralMoment3ExtQ_ten_neg⟩

end PT.Conservation
