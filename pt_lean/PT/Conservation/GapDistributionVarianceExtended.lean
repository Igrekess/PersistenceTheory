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
import Mathlib.Algebra.BigOperators.Intervals
import Mathlib.Data.Rat.Defs
import Mathlib.Tactic

/-!
# Distributional variance of the PT gap distribution — extension `N ∈ {5, …, 10}`

`PT.Conservation.GapDistributionVariance` introduces, on the small-`N`
sequence `ptPrime` (first five primes), the **distributional variance**

  `gapsDistVarianceQ N := (Σ_{n=1}^{N} (n − μ_N)² · g_n) / (Σ_{n=1}^{N} g_n)`,

where `μ_N := gapsDistMeanQ N`, and discharges the value
`gapsDistVarianceQ 4 = 10/9`.

This file extends that construction to the **eleven-prime sequence**
`ptPrimeExt` defined in `PT.Conservation.CumulativeBoundsExtended`,
covering the next six cases `N ∈ {5, 6, 7, 8, 9, 10}`.

With gaps `(g_1, …, g_{10}) = (1, 2, 2, 4, 2, 4, 2, 4, 6, 2)`, the
distributional means `μ_N = (Σ n · g_n)/S_N` and second weighted moments
`Σ n² · g_n` give the following exact rational variances
(from `Σ (n − μ_N)² g_n = (Σ n² g_n) · S_N − (Σ n g_n)²) / S_N²`):

| `N`  | `S_N` | `Σ n · g_n` | `Σ n² · g_n` | `gapsDistVarianceExtQ N` |
|------|-------|-------------|---------------|--------------------------|
| `5`  | `11`  | `37`        | `141`         | `182  / 121 ≈ 1.504`     |
| `6`  | `15`  | `61`        | `285`         | `554  / 225 ≈ 2.462`     |
| `7`  | `17`  | `75`        | `383`         | `886  / 289 ≈ 3.066`     |
| `8`  | `21`  | `107`       | `639`         | `1970 / 441 ≈ 4.467`     |
| `9`  | `27`  | `161`       | `1125`        | `4454 / 729 ≈ 6.110`     |
| `10` | `29`  | `181`       | `1325`        | `5664 / 841 ≈ 6.735`     |

The sequence `gapsDistVarianceExtQ N` is **strictly increasing in `N`**
on this range: as more high-`n` gaps enter, the index variance under
the gap weighting grows monotonically — a discrete witness of the
spread of the right-biased distribution.

All values are exact rationals, discharged by unfolding the `Finset.Ico`
sums and closing with `norm_num`.

## Reference

* `PT.Conservation.CumulativeBoundsExtended` — `ptPrimeExt`, gaps,
  cumulative sums `S_N`.
* `PT.Conservation.GapDistributionMeanExtended` — `gapsDistMeanExtQ N`
  and the six exact mean values.
* `PT.Conservation.GapMomentsNExtended` — integer linear-weighted sum
  `linWeightedSumGapsExt N`.
* `PT.Conservation.GapDistributionVariance` — `gapsDistVarianceQ N` on
  the small-`N` sequence (`N ≤ 4`).
-/

namespace PT.Conservation

open Finset

/-! ### Distributional variance on the extended range -/

/-- **Extended distributional variance.** Variance of the index `n`
    under the probability distribution `n ↦ g_n / S_N` on the
    eleven-prime sequence `ptPrimeExt`:

      `gapsDistVarianceExtQ N
        = (Σ_{n=1}^{N} (n − μ_N)² · g_n) / S_N`,

    where `μ_N = gapsDistMeanExtQ N`. -/
def gapsDistVarianceExtQ (N : ℕ) : ℚ :=
  (∑ n ∈ Ico (1 : ℕ) (N + 1),
      ((n : ℚ) - gapsDistMeanExtQ N) ^ 2 * (gap ptPrimeExt n : ℚ))
    / cumGapExtQ N

/-! ### Rational gap values on the extended range -/

/-- `g_1 = 1` as a rational. -/
theorem gapExt_one_rat : ((gap ptPrimeExt 1 : ℤ) : ℚ) = 1 := by
  rw [gapExt_one]; norm_num

/-- `g_2 = 2` as a rational. -/
theorem gapExt_two_rat : ((gap ptPrimeExt 2 : ℤ) : ℚ) = 2 := by
  rw [gapExt_two]; norm_num

/-- `g_3 = 2` as a rational. -/
theorem gapExt_three_rat : ((gap ptPrimeExt 3 : ℤ) : ℚ) = 2 := by
  rw [gapExt_three]; norm_num

/-- `g_4 = 4` as a rational. -/
theorem gapExt_four_rat : ((gap ptPrimeExt 4 : ℤ) : ℚ) = 4 := by
  rw [gapExt_four]; norm_num

/-- `g_5 = 2` as a rational. -/
theorem gapExt_five_rat : ((gap ptPrimeExt 5 : ℤ) : ℚ) = 2 := by
  rw [gapExt_five]; norm_num

/-- `g_6 = 4` as a rational. -/
theorem gapExt_six_rat : ((gap ptPrimeExt 6 : ℤ) : ℚ) = 4 := by
  rw [gapExt_six]; norm_num

/-- `g_7 = 2` as a rational. -/
theorem gapExt_seven_rat : ((gap ptPrimeExt 7 : ℤ) : ℚ) = 2 := by
  rw [gapExt_seven]; norm_num

/-- `g_8 = 4` as a rational. -/
theorem gapExt_eight_rat : ((gap ptPrimeExt 8 : ℤ) : ℚ) = 4 := by
  rw [gapExt_eight]; norm_num

/-- `g_9 = 6` as a rational. -/
theorem gapExt_nine_rat : ((gap ptPrimeExt 9 : ℤ) : ℚ) = 6 := by
  rw [gapExt_nine]; norm_num

/-- `g_{10} = 2` as a rational. -/
theorem gapExt_ten_rat : ((gap ptPrimeExt 10 : ℤ) : ℚ) = 2 := by
  rw [gapExt_ten]; norm_num

/-! ### Six exact values `N ∈ {5, …, 10}`

Each variance is discharged by unfolding the `Ico` sum step by step and
closing with `norm_num`. The closed-form check used to derive the
target values is

  `gapsDistVarianceExtQ N
    = ((Σ n² · g_n) · S_N − (Σ n · g_n)²) / S_N²`,

so e.g. for `N = 5`:
  `(141·11 − 37²) / 11² = (1551 − 1369) / 121 = 182 / 121`.
-/

/-- **Distributional variance at `N = 5`.** `182 / 121`. -/
theorem gapsDistVarianceExtQ_five : gapsDistVarianceExtQ 5 = 182 / 121 := by
  unfold gapsDistVarianceExtQ
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

/-- **Distributional variance at `N = 6`.** `554 / 225`. -/
theorem gapsDistVarianceExtQ_six : gapsDistVarianceExtQ 6 = 554 / 225 := by
  unfold gapsDistVarianceExtQ
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

/-- **Distributional variance at `N = 7`.** `886 / 289`. -/
theorem gapsDistVarianceExtQ_seven : gapsDistVarianceExtQ 7 = 886 / 289 := by
  unfold gapsDistVarianceExtQ
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

/-- **Distributional variance at `N = 8`.** `1970 / 441`. -/
theorem gapsDistVarianceExtQ_eight : gapsDistVarianceExtQ 8 = 1970 / 441 := by
  unfold gapsDistVarianceExtQ
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

/-- **Distributional variance at `N = 9`.** `4454 / 729`. -/
theorem gapsDistVarianceExtQ_nine : gapsDistVarianceExtQ 9 = 4454 / 729 := by
  unfold gapsDistVarianceExtQ
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

/-- **Distributional variance at `N = 10`.** `5664 / 841`. -/
theorem gapsDistVarianceExtQ_ten : gapsDistVarianceExtQ 10 = 5664 / 841 := by
  unfold gapsDistVarianceExtQ
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

/-! ### Strict monotonicity on `{5, …, 10}`

The distributional variance is strictly increasing in `N` on the extended
range: numerically `182/121 < 554/225 < 886/289 < 1970/441 < 4454/729 <
5664/841`. -/

/-- `gapsDistVarianceExtQ 5 < gapsDistVarianceExtQ 6`. -/
theorem gapsDistVarianceExtQ_five_lt_six :
    gapsDistVarianceExtQ 5 < gapsDistVarianceExtQ 6 := by
  rw [gapsDistVarianceExtQ_five, gapsDistVarianceExtQ_six]; norm_num

/-- `gapsDistVarianceExtQ 6 < gapsDistVarianceExtQ 7`. -/
theorem gapsDistVarianceExtQ_six_lt_seven :
    gapsDistVarianceExtQ 6 < gapsDistVarianceExtQ 7 := by
  rw [gapsDistVarianceExtQ_six, gapsDistVarianceExtQ_seven]; norm_num

/-- `gapsDistVarianceExtQ 7 < gapsDistVarianceExtQ 8`. -/
theorem gapsDistVarianceExtQ_seven_lt_eight :
    gapsDistVarianceExtQ 7 < gapsDistVarianceExtQ 8 := by
  rw [gapsDistVarianceExtQ_seven, gapsDistVarianceExtQ_eight]; norm_num

/-- `gapsDistVarianceExtQ 8 < gapsDistVarianceExtQ 9`. -/
theorem gapsDistVarianceExtQ_eight_lt_nine :
    gapsDistVarianceExtQ 8 < gapsDistVarianceExtQ 9 := by
  rw [gapsDistVarianceExtQ_eight, gapsDistVarianceExtQ_nine]; norm_num

/-- `gapsDistVarianceExtQ 9 < gapsDistVarianceExtQ 10`. -/
theorem gapsDistVarianceExtQ_nine_lt_ten :
    gapsDistVarianceExtQ 9 < gapsDistVarianceExtQ 10 := by
  rw [gapsDistVarianceExtQ_nine, gapsDistVarianceExtQ_ten]; norm_num

/-! ### Continuity with the small-`N` value

At `N = 4`, the eleven-prime sequence `ptPrimeExt` and the small-`N`
sequence `ptPrime` agree, so the distributional variance is unambiguous
there: `10 / 9`. The first extended value `gapsDistVarianceExtQ 5` is
strictly larger than this anchor. -/

/-- `gapsDistVarianceQ 4 = 10/9 < 182/121 = gapsDistVarianceExtQ 5`. -/
theorem gapsDistVarianceQ_four_lt_gapsDistVarianceExtQ_five :
    gapsDistVarianceQ 4 < gapsDistVarianceExtQ 5 := by
  rw [gapsDistVarianceQ_four, gapsDistVarianceExtQ_five]; norm_num

/-! ### Headline summary -/

/-- **Headline (extended distributional variances, `N ∈ {5, …, 10}`).**

    For the eleven-prime PT sequence `ptPrimeExt` with gaps
    `(1, 2, 2, 4, 2, 4, 2, 4, 6, 2)`, the six distributional variances

      `gapsDistVarianceExtQ N = (Σ (n − μ_N)² · g_n) / (Σ g_n)`

    are exactly

      `182/121, 554/225, 886/289, 1970/441, 4454/729, 5664/841`.

    The sequence is strictly increasing in `N` on this range, and the
    first extended value `182/121 ≈ 1.504` is itself strictly above the
    small-`N` anchor `gapsDistVarianceQ 4 = 10/9 ≈ 1.111`. -/
theorem gapsDistVarianceExt_headline :
    -- six exact values
    (gapsDistVarianceExtQ 5  = 182  / 121)
    ∧ (gapsDistVarianceExtQ 6  = 554  / 225)
    ∧ (gapsDistVarianceExtQ 7  = 886  / 289)
    ∧ (gapsDistVarianceExtQ 8  = 1970 / 441)
    ∧ (gapsDistVarianceExtQ 9  = 4454 / 729)
    ∧ (gapsDistVarianceExtQ 10 = 5664 / 841)
    -- strict monotonicity on the extended range
    ∧ (gapsDistVarianceExtQ 5 < gapsDistVarianceExtQ 6)
    ∧ (gapsDistVarianceExtQ 6 < gapsDistVarianceExtQ 7)
    ∧ (gapsDistVarianceExtQ 7 < gapsDistVarianceExtQ 8)
    ∧ (gapsDistVarianceExtQ 8 < gapsDistVarianceExtQ 9)
    ∧ (gapsDistVarianceExtQ 9 < gapsDistVarianceExtQ 10)
    -- continuity with the small-`N` anchor
    ∧ (gapsDistVarianceQ 4 < gapsDistVarianceExtQ 5) :=
  ⟨gapsDistVarianceExtQ_five,  gapsDistVarianceExtQ_six,
   gapsDistVarianceExtQ_seven, gapsDistVarianceExtQ_eight,
   gapsDistVarianceExtQ_nine,  gapsDistVarianceExtQ_ten,
   gapsDistVarianceExtQ_five_lt_six,   gapsDistVarianceExtQ_six_lt_seven,
   gapsDistVarianceExtQ_seven_lt_eight, gapsDistVarianceExtQ_eight_lt_nine,
   gapsDistVarianceExtQ_nine_lt_ten,
   gapsDistVarianceQ_four_lt_gapsDistVarianceExtQ_five⟩

end PT.Conservation
