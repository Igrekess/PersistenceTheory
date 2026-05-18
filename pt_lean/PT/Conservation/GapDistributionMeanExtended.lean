/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Conservation.ConservationID
import PT.Conservation.ConservationIDExtensions
import PT.Conservation.CumulativeBoundsExtended
import PT.Conservation.GapMomentsNExtended
import PT.Conservation.GapDistributionVariance
import Mathlib.Algebra.BigOperators.Intervals
import Mathlib.Data.Rat.Defs
import Mathlib.Tactic

/-!
# Distributional mean of the PT gap distribution — extension `N ∈ {5, …, 10}`

`PT.Conservation.GapDistributionVariance` introduces, on the small-`N`
sequence `ptPrime` (first five primes), the **distributional mean**

  `gapsDistMeanQ N := (Σ_{n=1}^{N} n · g_n) / (Σ_{n=1}^{N} g_n)`,

and discharges the worked-out value `gapsDistMeanQ 4 = 3`.

This file extends that construction to the **eleven-prime sequence**
`ptPrimeExt` defined in `PT.Conservation.CumulativeBoundsExtended`,
covering the next six cases `N ∈ {5, 6, 7, 8, 9, 10}`.

With gaps `(g_1, …, g_{10}) = (1, 2, 2, 4, 2, 4, 2, 4, 6, 2)`, the
cumulative sums and weighted index sums on the extended range are

| `N`  | `S_N = Σ g_n` | `Σ n · g_n` | `gapsDistMeanExtQ N` |
|------|---------------|-------------|----------------------|
| `5`  | `11`          | `37`        | `37/11 ≈ 3.364`      |
| `6`  | `15`          | `61`        | `61/15 ≈ 4.067`      |
| `7`  | `17`          | `75`        | `75/17 ≈ 4.412`      |
| `8`  | `21`          | `107`       | `107/21 ≈ 5.095`     |
| `9`  | `27`          | `161`       | `161/27 ≈ 5.963`     |
| `10` | `29`          | `181`       | `181/29 ≈ 6.241`     |

The sequence `gapsDistMeanExtQ N` is **strictly increasing in `N`** on
this range: the distributional mean tracks the centre of mass of the
heavy-tail-shifted gap distribution.

For comparison, the index-uniform mean on `{1, …, N}` is `(N+1)/2`:

| `N`  | `(N+1)/2`        | `gapsDistMeanExtQ N`  | direction |
|------|------------------|-----------------------|-----------|
| `5`  | `3`              | `37/11 ≈ 3.364`       | `> 3`     |
| `6`  | `7/2  = 3.5`     | `61/15 ≈ 4.067`       | `> 7/2`   |
| `7`  | `4`              | `75/17 ≈ 4.412`       | `> 4`     |
| `8`  | `9/2  = 4.5`     | `107/21 ≈ 5.095`      | `> 9/2`   |
| `9`  | `5`              | `161/27 ≈ 5.963`      | `> 5`     |
| `10` | `11/2 = 5.5`     | `181/29 ≈ 6.241`      | `> 11/2`  |

In every case, `gapsDistMeanExtQ N > (N+1)/2`: the gap distribution is
biased towards the **right** of the index range (large gaps cluster at
high `n`, in line with the slow growth of primes).

All values are exact rationals, discharged by unfolding the `Finset.Ico`
sums and closing with `decide` / `norm_num`.

## Reference

* `PT.Conservation.CumulativeBoundsExtended` — `ptPrimeExt`, gaps,
  cumulative sums `S_N`.
* `PT.Conservation.GapMomentsNExtended` — `linWeightedSumGapsExt N`
  (the integer numerator `Σ n · g_n`).
* `PT.Conservation.GapDistributionVariance` — `gapsDistMeanQ N`
  on the small-`N` sequence.
-/

namespace PT.Conservation

open Finset

/-! ### Rational cumulative gap sum on the extended range -/

/-- The cumulative gap sum on `Ico 1 (N + 1)` for the extended sequence
    `ptPrimeExt`, in `ℚ`. This is the rational counterpart of the
    integer-valued `∑ n ∈ Ico 1 (N+1), gap ptPrimeExt n` defined via
    `CumulativeBoundsExtended`. -/
def cumGapExtQ (N : ℕ) : ℚ :=
  ∑ n ∈ Ico 1 (N + 1), (gap ptPrimeExt n : ℚ)

/-- **Integer bridge.** The rational cumulative gap sum equals the rational
    cast of the integer-valued cumulative sum from
    `PT.Conservation.CumulativeBoundsExtended`. -/
theorem cumGapExtQ_eq_int_cast (N : ℕ) :
    cumGapExtQ N = ((∑ n ∈ Ico 1 (N + 1), gap ptPrimeExt n : ℤ) : ℚ) := by
  unfold cumGapExtQ
  push_cast
  rfl

/-- `cumGapExtQ 5 = 11`. -/
theorem cumGapExtQ_five : cumGapExtQ 5 = 11 := by
  rw [cumGapExtQ_eq_int_cast, cumulativeExt_N5]; rfl

/-- `cumGapExtQ 6 = 15`. -/
theorem cumGapExtQ_six : cumGapExtQ 6 = 15 := by
  rw [cumGapExtQ_eq_int_cast, cumulativeExt_N6]; rfl

/-- `cumGapExtQ 7 = 17`. -/
theorem cumGapExtQ_seven : cumGapExtQ 7 = 17 := by
  rw [cumGapExtQ_eq_int_cast, cumulativeExt_N7]; rfl

/-- `cumGapExtQ 8 = 21`. -/
theorem cumGapExtQ_eight : cumGapExtQ 8 = 21 := by
  rw [cumGapExtQ_eq_int_cast, cumulativeExt_N8]; rfl

/-- `cumGapExtQ 9 = 27`. -/
theorem cumGapExtQ_nine : cumGapExtQ 9 = 27 := by
  rw [cumGapExtQ_eq_int_cast, cumulativeExt_N9]; rfl

/-- `cumGapExtQ 10 = 29`. -/
theorem cumGapExtQ_ten : cumGapExtQ 10 = 29 := by
  rw [cumGapExtQ_eq_int_cast, cumulativeExt_N10]; rfl

/-! ### Distributional mean on the extended range -/

/-- **Extended distributional mean.** Mean of the index `n` under the
    probability distribution `n ↦ g_n / S_N` on the eleven-prime
    sequence `ptPrimeExt`:

      `gapsDistMeanExtQ N = (Σ_{n=1}^{N} n · g_n) / S_N`. -/
def gapsDistMeanExtQ (N : ℕ) : ℚ :=
  (∑ n ∈ Ico (1 : ℕ) (N + 1),
      (n : ℚ) * (gap ptPrimeExt n : ℚ)) / cumGapExtQ N

/-! ### Numerator bridge with `linWeightedSumGapsExt`

The numerator of `gapsDistMeanExtQ N` is exactly the rational image of the
integer first moment `linWeightedSumGapsExt N` from
`PT.Conservation.GapMomentsNExtended`. -/

/-- **Numerator bridge.** The integer first-moment `linWeightedSumGapsExt N`,
    cast in `ℚ`, equals the numerator of `gapsDistMeanExtQ N`. -/
theorem gapsDistMeanExtQ_numerator_eq_linWeighted (N : ℕ) :
    (∑ n ∈ Ico (1 : ℕ) (N + 1), (n : ℚ) * (gap ptPrimeExt n : ℚ))
      = (linWeightedSumGapsExt N : ℚ) := by
  unfold linWeightedSumGapsExt
  push_cast
  rfl

/-! ### Six exact values `N ∈ {5, …, 10}` -/

/-- **Distributional mean at `N = 5`.** `(1·1 + 2·2 + 3·2 + 4·4 + 5·2) / 11
    = 37 / 11`. -/
theorem gapsDistMeanExtQ_five : gapsDistMeanExtQ 5 = 37 / 11 := by
  unfold gapsDistMeanExtQ
  rw [gapsDistMeanExtQ_numerator_eq_linWeighted, linWeightedSumGapsExt_five,
      cumGapExtQ_five]
  norm_num

/-- **Distributional mean at `N = 6`.** `61 / 15`. -/
theorem gapsDistMeanExtQ_six : gapsDistMeanExtQ 6 = 61 / 15 := by
  unfold gapsDistMeanExtQ
  rw [gapsDistMeanExtQ_numerator_eq_linWeighted, linWeightedSumGapsExt_six,
      cumGapExtQ_six]
  norm_num

/-- **Distributional mean at `N = 7`.** `75 / 17`. -/
theorem gapsDistMeanExtQ_seven : gapsDistMeanExtQ 7 = 75 / 17 := by
  unfold gapsDistMeanExtQ
  rw [gapsDistMeanExtQ_numerator_eq_linWeighted, linWeightedSumGapsExt_seven,
      cumGapExtQ_seven]
  norm_num

/-- **Distributional mean at `N = 8`.** `107 / 21`. -/
theorem gapsDistMeanExtQ_eight : gapsDistMeanExtQ 8 = 107 / 21 := by
  unfold gapsDistMeanExtQ
  rw [gapsDistMeanExtQ_numerator_eq_linWeighted, linWeightedSumGapsExt_eight,
      cumGapExtQ_eight]
  norm_num

/-- **Distributional mean at `N = 9`.** `161 / 27`. -/
theorem gapsDistMeanExtQ_nine : gapsDistMeanExtQ 9 = 161 / 27 := by
  unfold gapsDistMeanExtQ
  rw [gapsDistMeanExtQ_numerator_eq_linWeighted, linWeightedSumGapsExt_nine,
      cumGapExtQ_nine]
  norm_num

/-- **Distributional mean at `N = 10`.** `181 / 29`. -/
theorem gapsDistMeanExtQ_ten : gapsDistMeanExtQ 10 = 181 / 29 := by
  unfold gapsDistMeanExtQ
  rw [gapsDistMeanExtQ_numerator_eq_linWeighted, linWeightedSumGapsExt_ten,
      cumGapExtQ_ten]
  norm_num

/-! ### Comparison with the uniform-index mean `(N+1)/2` -/

/-- **Right-bias at `N = 5`.** `37/11 > 3 = (5+1)/2`. -/
theorem gapsDistMeanExtQ_five_gt_uniform :
    gapsDistMeanExtQ 5 > (5 + 1 : ℚ) / 2 := by
  rw [gapsDistMeanExtQ_five]; norm_num

/-- **Right-bias at `N = 6`.** `61/15 > 7/2 = (6+1)/2`. -/
theorem gapsDistMeanExtQ_six_gt_uniform :
    gapsDistMeanExtQ 6 > (6 + 1 : ℚ) / 2 := by
  rw [gapsDistMeanExtQ_six]; norm_num

/-- **Right-bias at `N = 7`.** `75/17 > 4 = (7+1)/2`. -/
theorem gapsDistMeanExtQ_seven_gt_uniform :
    gapsDistMeanExtQ 7 > (7 + 1 : ℚ) / 2 := by
  rw [gapsDistMeanExtQ_seven]; norm_num

/-- **Right-bias at `N = 8`.** `107/21 > 9/2 = (8+1)/2`. -/
theorem gapsDistMeanExtQ_eight_gt_uniform :
    gapsDistMeanExtQ 8 > (8 + 1 : ℚ) / 2 := by
  rw [gapsDistMeanExtQ_eight]; norm_num

/-- **Right-bias at `N = 9`.** `161/27 > 5 = (9+1)/2`. -/
theorem gapsDistMeanExtQ_nine_gt_uniform :
    gapsDistMeanExtQ 9 > (9 + 1 : ℚ) / 2 := by
  rw [gapsDistMeanExtQ_nine]; norm_num

/-- **Right-bias at `N = 10`.** `181/29 > 11/2 = (10+1)/2`. -/
theorem gapsDistMeanExtQ_ten_gt_uniform :
    gapsDistMeanExtQ 10 > (10 + 1 : ℚ) / 2 := by
  rw [gapsDistMeanExtQ_ten]; norm_num

/-! ### Strict monotonicity on `{5, …, 10}`

The distributional mean is strictly increasing in `N` on the extended
range. This is a finite numerical witness, not a structural result.
-/

/-- `gapsDistMeanExtQ 5 < gapsDistMeanExtQ 6`. -/
theorem gapsDistMeanExtQ_five_lt_six :
    gapsDistMeanExtQ 5 < gapsDistMeanExtQ 6 := by
  rw [gapsDistMeanExtQ_five, gapsDistMeanExtQ_six]; norm_num

/-- `gapsDistMeanExtQ 6 < gapsDistMeanExtQ 7`. -/
theorem gapsDistMeanExtQ_six_lt_seven :
    gapsDistMeanExtQ 6 < gapsDistMeanExtQ 7 := by
  rw [gapsDistMeanExtQ_six, gapsDistMeanExtQ_seven]; norm_num

/-- `gapsDistMeanExtQ 7 < gapsDistMeanExtQ 8`. -/
theorem gapsDistMeanExtQ_seven_lt_eight :
    gapsDistMeanExtQ 7 < gapsDistMeanExtQ 8 := by
  rw [gapsDistMeanExtQ_seven, gapsDistMeanExtQ_eight]; norm_num

/-- `gapsDistMeanExtQ 8 < gapsDistMeanExtQ 9`. -/
theorem gapsDistMeanExtQ_eight_lt_nine :
    gapsDistMeanExtQ 8 < gapsDistMeanExtQ 9 := by
  rw [gapsDistMeanExtQ_eight, gapsDistMeanExtQ_nine]; norm_num

/-- `gapsDistMeanExtQ 9 < gapsDistMeanExtQ 10`. -/
theorem gapsDistMeanExtQ_nine_lt_ten :
    gapsDistMeanExtQ 9 < gapsDistMeanExtQ 10 := by
  rw [gapsDistMeanExtQ_nine, gapsDistMeanExtQ_ten]; norm_num

/-! ### Headline summary -/

/-- **Headline (extended distributional means, `N ∈ {5, …, 10}`).**

    For the eleven-prime PT sequence `ptPrimeExt` with gaps
    `(1, 2, 2, 4, 2, 4, 2, 4, 6, 2)`, the six distributional means

      `gapsDistMeanExtQ N = (Σ n · g_n) / (Σ g_n)`

    are exactly

      `37/11, 61/15, 75/17, 107/21, 161/27, 181/29`.

    All six lie strictly above the index-uniform mean `(N+1)/2`,
    confirming the right-bias of the PT gap distribution, and the
    sequence is strictly increasing in `N` on this range. -/
theorem gapsDistMeanExt_headline :
    -- six exact values
    (gapsDistMeanExtQ 5  = 37 / 11)
    ∧ (gapsDistMeanExtQ 6  = 61 / 15)
    ∧ (gapsDistMeanExtQ 7  = 75 / 17)
    ∧ (gapsDistMeanExtQ 8  = 107 / 21)
    ∧ (gapsDistMeanExtQ 9  = 161 / 27)
    ∧ (gapsDistMeanExtQ 10 = 181 / 29)
    -- six right-bias witnesses vs uniform mean
    ∧ (gapsDistMeanExtQ 5  > (5  + 1 : ℚ) / 2)
    ∧ (gapsDistMeanExtQ 6  > (6  + 1 : ℚ) / 2)
    ∧ (gapsDistMeanExtQ 7  > (7  + 1 : ℚ) / 2)
    ∧ (gapsDistMeanExtQ 8  > (8  + 1 : ℚ) / 2)
    ∧ (gapsDistMeanExtQ 9  > (9  + 1 : ℚ) / 2)
    ∧ (gapsDistMeanExtQ 10 > (10 + 1 : ℚ) / 2)
    -- strict monotonicity on the extended range
    ∧ (gapsDistMeanExtQ 5 < gapsDistMeanExtQ 6)
    ∧ (gapsDistMeanExtQ 6 < gapsDistMeanExtQ 7)
    ∧ (gapsDistMeanExtQ 7 < gapsDistMeanExtQ 8)
    ∧ (gapsDistMeanExtQ 8 < gapsDistMeanExtQ 9)
    ∧ (gapsDistMeanExtQ 9 < gapsDistMeanExtQ 10) :=
  ⟨gapsDistMeanExtQ_five,  gapsDistMeanExtQ_six,
   gapsDistMeanExtQ_seven, gapsDistMeanExtQ_eight,
   gapsDistMeanExtQ_nine,  gapsDistMeanExtQ_ten,
   gapsDistMeanExtQ_five_gt_uniform,  gapsDistMeanExtQ_six_gt_uniform,
   gapsDistMeanExtQ_seven_gt_uniform, gapsDistMeanExtQ_eight_gt_uniform,
   gapsDistMeanExtQ_nine_gt_uniform,  gapsDistMeanExtQ_ten_gt_uniform,
   gapsDistMeanExtQ_five_lt_six,   gapsDistMeanExtQ_six_lt_seven,
   gapsDistMeanExtQ_seven_lt_eight, gapsDistMeanExtQ_eight_lt_nine,
   gapsDistMeanExtQ_nine_lt_ten⟩

end PT.Conservation
