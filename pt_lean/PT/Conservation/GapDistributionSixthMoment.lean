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
import PT.Conservation.GapDistributionFifthMoment
import PT.Conservation.GapStatisticalMoments
import Mathlib.Algebra.BigOperators.Intervals
import Mathlib.Data.Rat.Defs
import Mathlib.Tactic

/-!
# Distributional sixth central moment of the PT gap distribution —
# extension `N ∈ {5, …, 10}`

This file is the **sixth-moment** companion to
`PT.Conservation.GapDistributionFifthMoment` (fifth central moment),
`PT.Conservation.GapDistributionKurtosisExtended` (fourth central
moment), `PT.Conservation.GapDistributionSkewExtended` (third central
moment) and the variance/mean files. On the eleven-prime sequence
`ptPrimeExt`, it computes the **distributional sixth central moment**
of the index `n` under the probability distribution `n ↦ g_n / S_N`:

  `gapsDistCentralMoment6ExtQ N
      := (Σ_{n=1}^{N} (n − μ_N)⁶ · g_n) / S_N`,

where `μ_N := gapsDistMeanExtQ N` and `S_N := Σ_{n=1}^{N} g_n`.

Because the sixth power is **even**, the sixth central moment is a
**uniformly non-negative** functional of the centred index
distribution. Together with the (negative) third and fifth central
moments, it gives a finer view of the tail mass: the sixth power
amplifies the tails even more strongly than the fourth, so the
sixth-to-fourth ratio is a higher-resolution tail diagnostic for the
gap-weighted index distribution.

With gaps `(g_1, …, g_{10}) = (1, 2, 2, 4, 2, 4, 2, 4, 6, 2)`, the six
exact rationals are

| `N`  | `S_N` | `μ_N`    | `gapsDistCentralMoment6ExtQ N`         |
|------|-------|----------|----------------------------------------|
| `5`  | `11`  | `37/11`  | `36381842 / 1771561`                   |
| `6`  | `15`  | `61/15`  | `182362814 / 2278125`                  |
| `7`  | `17`  | `75/17`  | `3765640426 / 24137569`                |
| `8`  | `21`  | `107/21` | `37379593750 / 85766121`               |
| `9`  | `27`  | `161/27` | `420003215194 / 387420489`             |
| `10` | `29`  | `181/29` | `894548024304 / 594823321`             |

The denominators are reduced rationals derived from `S_N⁶`:
`1771561 = 11⁶`, `2278125 = 15⁶/5`, `24137569 = 17⁶`,
`85766121 = 21⁶`, `387420489 = 27⁶`, `594823321 = 29⁶`.

**All six values are strictly positive** — as forced by the even
sixth power applied to the centred index distribution — confirming
the non-degeneracy of the gap-weighted index distribution on the
extended range `N ∈ {5, …, 10}` at the sixth-moment level.

All values are exact rationals, discharged by unfolding the `Finset.Ico`
sums step by step and closing with `norm_num`.

## Reference

* `PT.Conservation.CumulativeBoundsExtended` — `ptPrimeExt`, gaps,
  cumulative sums `S_N`.
* `PT.Conservation.GapDistributionMeanExtended` — `gapsDistMeanExtQ N`
  and the six exact mean values used here.
* `PT.Conservation.GapDistributionVarianceExtended` — the rational gap
  lemmas `gapExt_one_rat, …, gapExt_ten_rat` reused below.
* `PT.Conservation.GapDistributionSkewExtended` — third central
  moment companion (strictly negative across `N ∈ {5,…,10}`).
* `PT.Conservation.GapDistributionKurtosisExtended` — fourth central
  moment companion (strictly positive, platykurtic).
* `PT.Conservation.GapDistributionFifthMoment` — fifth central
  moment companion (strictly negative; tail-amplified left skew).
-/

namespace PT.Conservation

open Finset

/-! ### Distributional sixth central moment on the extended range -/

/-- **Extended distributional sixth central moment.** Sixth central
    moment of the index `n` under the probability distribution
    `n ↦ g_n / S_N` on the eleven-prime sequence `ptPrimeExt`:

      `gapsDistCentralMoment6ExtQ N
          = (Σ_{n=1}^{N} (n − μ_N)⁶ · g_n) / S_N`,

    where `μ_N = gapsDistMeanExtQ N`. Being an even-power moment, this
    quantity is always non-negative; strict positivity is established
    individually for each `N ∈ {5, …, 10}` below. -/
def gapsDistCentralMoment6ExtQ (N : ℕ) : ℚ :=
  (∑ n ∈ Ico (1 : ℕ) (N + 1),
      ((n : ℚ) - gapsDistMeanExtQ N) ^ 6 * (gap ptPrimeExt n : ℚ))
    / cumGapExtQ N

/-! ### Six exact values `N ∈ {5, …, 10}`

Each sixth central moment is discharged by unfolding the `Ico` sum
step by step (reducing it to a sum of `N` explicit terms) and closing
with `norm_num`, after rewriting `μ_N = gapsDistMeanExtQ N` and
`S_N = cumGapExtQ N` to their exact rational values.
-/

/-- **Distributional sixth central moment at `N = 5`.**
    `36381842 / 1771561` with `1771561 = 11⁶`. -/
theorem gapsDistCentralMoment6ExtQ_five :
    gapsDistCentralMoment6ExtQ 5 = 36381842 / 1771561 := by
  unfold gapsDistCentralMoment6ExtQ
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

/-- **Distributional sixth central moment at `N = 6`.**
    `182362814 / 2278125`. -/
theorem gapsDistCentralMoment6ExtQ_six :
    gapsDistCentralMoment6ExtQ 6 = 182362814 / 2278125 := by
  unfold gapsDistCentralMoment6ExtQ
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

/-- **Distributional sixth central moment at `N = 7`.**
    `3765640426 / 24137569` with `24137569 = 17⁶`. -/
theorem gapsDistCentralMoment6ExtQ_seven :
    gapsDistCentralMoment6ExtQ 7 = 3765640426 / 24137569 := by
  unfold gapsDistCentralMoment6ExtQ
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

/-- **Distributional sixth central moment at `N = 8`.**
    `37379593750 / 85766121` with `85766121 = 21⁶`. -/
theorem gapsDistCentralMoment6ExtQ_eight :
    gapsDistCentralMoment6ExtQ 8 = 37379593750 / 85766121 := by
  unfold gapsDistCentralMoment6ExtQ
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

/-- **Distributional sixth central moment at `N = 9`.**
    `420003215194 / 387420489` with `387420489 = 27⁶`. -/
theorem gapsDistCentralMoment6ExtQ_nine :
    gapsDistCentralMoment6ExtQ 9 = 420003215194 / 387420489 := by
  unfold gapsDistCentralMoment6ExtQ
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

/-- **Distributional sixth central moment at `N = 10`.**
    `894548024304 / 594823321` with `594823321 = 29⁶`. -/
theorem gapsDistCentralMoment6ExtQ_ten :
    gapsDistCentralMoment6ExtQ 10 = 894548024304 / 594823321 := by
  unfold gapsDistCentralMoment6ExtQ
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

/-! ### Strict positivity (uniform on the extended range)

Each sixth central moment is strictly positive, as forced by the even
sixth power applied to the centred index distribution: the centred
indices `(n − μ_N)` are not all zero (the supports are non-trivial),
hence raising to an even power and weighting by the positive gap
counts gives a strictly positive sum. We discharge each case by
rewriting to the explicit rational value. -/

/-- `gapsDistCentralMoment6ExtQ 5 > 0`. -/
theorem gapsDistCentralMoment6ExtQ_five_pos :
    gapsDistCentralMoment6ExtQ 5 > 0 := by
  rw [gapsDistCentralMoment6ExtQ_five]; norm_num

/-- `gapsDistCentralMoment6ExtQ 6 > 0`. -/
theorem gapsDistCentralMoment6ExtQ_six_pos :
    gapsDistCentralMoment6ExtQ 6 > 0 := by
  rw [gapsDistCentralMoment6ExtQ_six]; norm_num

/-- `gapsDistCentralMoment6ExtQ 7 > 0`. -/
theorem gapsDistCentralMoment6ExtQ_seven_pos :
    gapsDistCentralMoment6ExtQ 7 > 0 := by
  rw [gapsDistCentralMoment6ExtQ_seven]; norm_num

/-- `gapsDistCentralMoment6ExtQ 8 > 0`. -/
theorem gapsDistCentralMoment6ExtQ_eight_pos :
    gapsDistCentralMoment6ExtQ 8 > 0 := by
  rw [gapsDistCentralMoment6ExtQ_eight]; norm_num

/-- `gapsDistCentralMoment6ExtQ 9 > 0`. -/
theorem gapsDistCentralMoment6ExtQ_nine_pos :
    gapsDistCentralMoment6ExtQ 9 > 0 := by
  rw [gapsDistCentralMoment6ExtQ_nine]; norm_num

/-- `gapsDistCentralMoment6ExtQ 10 > 0`. -/
theorem gapsDistCentralMoment6ExtQ_ten_pos :
    gapsDistCentralMoment6ExtQ 10 > 0 := by
  rw [gapsDistCentralMoment6ExtQ_ten]; norm_num

/-! ### Headline summary -/

/-- **Headline (extended distributional sixth central moment,
    `N ∈ {5, …, 10}`).**

    For the eleven-prime PT sequence `ptPrimeExt` with gaps
    `(1, 2, 2, 4, 2, 4, 2, 4, 6, 2)`, the six distributional sixth
    central moments

      `gapsDistCentralMoment6ExtQ N
          = (Σ (n − μ_N)⁶ · g_n) / (Σ g_n)`

    are exactly

      `36381842/1771561, 182362814/2278125, 3765640426/24137569,
       37379593750/85766121, 420003215194/387420489,
       894548024304/594823321`.

    **All six are strictly positive**: as forced by the even sixth
    power, the gap-weighted index distribution has non-degenerate
    tail mass at the sixth-moment level on the extended range. -/
theorem gapsDistCentralMoment6Ext_headline :
    -- six exact sixth-moment values
    (gapsDistCentralMoment6ExtQ 5  = 36381842     / 1771561)
    ∧ (gapsDistCentralMoment6ExtQ 6  = 182362814    / 2278125)
    ∧ (gapsDistCentralMoment6ExtQ 7  = 3765640426   / 24137569)
    ∧ (gapsDistCentralMoment6ExtQ 8  = 37379593750  / 85766121)
    ∧ (gapsDistCentralMoment6ExtQ 9  = 420003215194 / 387420489)
    ∧ (gapsDistCentralMoment6ExtQ 10 = 894548024304 / 594823321)
    -- strict positivity on the extended range
    ∧ (gapsDistCentralMoment6ExtQ 5  > 0)
    ∧ (gapsDistCentralMoment6ExtQ 6  > 0)
    ∧ (gapsDistCentralMoment6ExtQ 7  > 0)
    ∧ (gapsDistCentralMoment6ExtQ 8  > 0)
    ∧ (gapsDistCentralMoment6ExtQ 9  > 0)
    ∧ (gapsDistCentralMoment6ExtQ 10 > 0) :=
  ⟨gapsDistCentralMoment6ExtQ_five,  gapsDistCentralMoment6ExtQ_six,
   gapsDistCentralMoment6ExtQ_seven, gapsDistCentralMoment6ExtQ_eight,
   gapsDistCentralMoment6ExtQ_nine,  gapsDistCentralMoment6ExtQ_ten,
   gapsDistCentralMoment6ExtQ_five_pos,   gapsDistCentralMoment6ExtQ_six_pos,
   gapsDistCentralMoment6ExtQ_seven_pos,  gapsDistCentralMoment6ExtQ_eight_pos,
   gapsDistCentralMoment6ExtQ_nine_pos,   gapsDistCentralMoment6ExtQ_ten_pos⟩

end PT.Conservation
