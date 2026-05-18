/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Conservation.ConservationID
import PT.Conservation.ConservationIDExtensions
import PT.Conservation.GapVarianceSmallN
import PT.Conservation.GapEntropyBound
import Mathlib.Algebra.BigOperators.Intervals
import Mathlib.Data.Rat.Defs
import Mathlib.Tactic

/-!
# Variance of the **normalised PT gap distribution** — small `N`

`PT.Conservation.GapVarianceSmallN` computes the **classical (population)
variance** of the gap *values* `(g_1, …, g_N)`:

  `meanGapsTo N     = (1/N) · Σ g_n`,
  `varianceGapsTo N = (1/N) · Σ (g_n − meanGapsTo N)²`.

In contrast, `PT.Conservation.GapEntropyBound` normalises the gaps into a
probability distribution on the index set `{1, …, N}`:

  `gapsDist N n := g_n / S_N`     with     `S_N := Σ g_n`.

The corresponding statistical functionals are taken **over the indices
`n`** weighted by the distribution `gapsDist N`:

  `gapsDistMeanQ N      := Σ n · (g_n / S_N)`,
  `gapsDistVarianceQ N  := Σ (n − gapsDistMeanQ N)² · (g_n / S_N)`.

This file gives definitions and discharges the small-`N` cases. The key
worked-out instance is `N = 4` (gaps `(1, 2, 2, 4)`, `S_4 = 9`):

  `gapsDistMeanQ 4 = (1·1 + 2·2 + 3·2 + 4·4) / 9 = 27/9 = 3`,
  `gapsDistVarianceQ 4 = ((1-3)² · 1 + (2-3)² · 2 + (3-3)² · 2 + (4-3)² · 4) / 9
                       = (4 + 2 + 0 + 4) / 9 = 10/9`.

In particular,
`gapsDistVarianceQ 4 = 10/9 ≠ 19/16 = varianceGapsTo 4`,
which shows that the two notions of "variance" (over *values*
vs. over *indices weighted by the value distribution*) genuinely
differ on the PT gap sequence.

All computations are exact rationals; we discharge them by unrolling
`Finset.Ico` sums and closing with `norm_num`.

## References

* `PT.Conservation.GapVarianceSmallN` — classical variance over values.
* `PT.Conservation.GapEntropyBound` — `gapsDist` (real-valued).
-/

namespace PT.Conservation

open Finset

/-! ### Rational cumulative gap sum -/

/-- The cumulative gap sum on `Ico 1 (N + 1)`, in `ℚ`.
    This is the rational counterpart of `cumGapReal`. -/
def cumGapQ (N : ℕ) : ℚ :=
  ∑ n ∈ Ico 1 (N + 1), (gap ptPrime n : ℚ)

theorem cumGapQ_two : cumGapQ 2 = 3 := by
  unfold cumGapQ
  rw [show (2 : ℕ) + 1 = 3 from rfl, sum_gap_ico_1_3_rat]

theorem cumGapQ_three : cumGapQ 3 = 5 := by
  unfold cumGapQ
  rw [show (3 : ℕ) + 1 = 4 from rfl, sum_gap_ico_1_4_rat]

theorem cumGapQ_four : cumGapQ 4 = 9 := by
  unfold cumGapQ
  rw [show (4 : ℕ) + 1 = 5 from rfl, sum_gap_ico_1_5_rat]

/-! ### Distributional mean of the PT gap distribution -/

/-- **Distributional mean.** The mean of the index `n` under the
    probability distribution `n ↦ g_n / S_N`, computed in `ℚ`:

      `gapsDistMeanQ N = (Σ n · g_n) / S_N`. -/
def gapsDistMeanQ (N : ℕ) : ℚ :=
  (∑ n ∈ Ico (1 : ℕ) (N + 1), (n : ℚ) * (gap ptPrime n : ℚ)) / cumGapQ N

/-- **Distributional variance.** Variance of the index `n` under the
    distribution `n ↦ g_n / S_N`:

      `gapsDistVarianceQ N = (Σ (n − μ)² · g_n) / S_N`,

    where `μ = gapsDistMeanQ N`. -/
def gapsDistVarianceQ (N : ℕ) : ℚ :=
  (∑ n ∈ Ico (1 : ℕ) (N + 1),
      ((n : ℚ) - gapsDistMeanQ N) ^ 2 * (gap ptPrime n : ℚ))
    / cumGapQ N

/-! ### Helper — closed form for `Σ n · g_n` on `Ico 1 5` -/

/-- Closed form `∑_{n=1}^{4} n · g_n = 1·1 + 2·2 + 3·2 + 4·4 = 27`. -/
theorem sum_n_gap_ico_1_5_rat :
    (∑ n ∈ Ico (1 : ℕ) 5, (n : ℚ) * (gap ptPrime n : ℚ)) = 27 := by
  rw [show (5 : ℕ) = 4 + 1 from rfl,
      Finset.sum_Ico_succ_top (by omega : 1 ≤ 4),
      show (4 : ℕ) = 3 + 1 from rfl,
      Finset.sum_Ico_succ_top (by omega : 1 ≤ 3),
      show (3 : ℕ) = 2 + 1 from rfl,
      Finset.sum_Ico_succ_top (by omega : 1 ≤ 2),
      show (2 : ℕ) = 1 + 1 from rfl,
      Finset.sum_Ico_succ_top (by omega : 1 ≤ 1)]
  simp [gap_one_rat, gap_two_rat, gap_three_rat, gap_four_rat]
  norm_num

/-! ### `N = 4` — distributional mean equals `3` -/

/-- **Distributional mean at `N = 4`.**

    `gapsDistMeanQ 4 = (1·1 + 2·2 + 3·2 + 4·4) / 9 = 27/9 = 3`. -/
theorem gapsDistMeanQ_four : gapsDistMeanQ 4 = 3 := by
  unfold gapsDistMeanQ
  rw [show (4 : ℕ) + 1 = 5 from rfl, sum_n_gap_ico_1_5_rat, cumGapQ_four]
  norm_num

/-! ### `N = 4` — distributional variance equals `10/9` -/

/-- **Distributional variance at `N = 4`.**

    With `μ = 3`, the weighted deviations are
    `(1-3)²·1 + (2-3)²·2 + (3-3)²·2 + (4-3)²·4 = 4 + 2 + 0 + 4 = 10`,
    hence `gapsDistVarianceQ 4 = 10 / 9`. -/
theorem gapsDistVarianceQ_four : gapsDistVarianceQ 4 = 10 / 9 := by
  unfold gapsDistVarianceQ
  rw [gapsDistMeanQ_four, cumGapQ_four,
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

/-! ### Comparison with the classical variance -/

/-- **Inequality.** The distributional variance and the classical
    (population) variance of the PT gaps **differ** at `N = 4`:

      `gapsDistVarianceQ 4 = 10/9 ≠ 19/16 = varianceGapsTo 4`. -/
theorem gapsDistVarianceQ_four_ne_varianceGapsTo_four :
    gapsDistVarianceQ 4 ≠ varianceGapsTo 4 := by
  rw [gapsDistVarianceQ_four, varianceGapsTo_four]
  norm_num

/-- **Strict ordering.** Numerically, `10/9 < 19/16`, i.e.
    the distributional variance of the index is *smaller* than the
    classical population variance of the gap values at `N = 4`. -/
theorem gapsDistVarianceQ_four_lt_varianceGapsTo_four :
    gapsDistVarianceQ 4 < varianceGapsTo 4 := by
  rw [gapsDistVarianceQ_four, varianceGapsTo_four]
  norm_num

/-! ### Headline summary -/

/-- **Headline (distributional mean & variance summary, `N = 4`).**

    For the PT prime gap sequence `(g_1, g_2, g_3, g_4) = (1, 2, 2, 4)`
    with cumulative sum `S_4 = 9`, the index distribution
    `n ↦ g_n / S_4` has

      mean `μ = 3`,  variance `σ² = 10/9 ≈ 1.111`,

    which is **strictly smaller** than the classical population variance
    of the gap *values* `varianceGapsTo 4 = 19/16 ≈ 1.188`. -/
theorem gapsDist_mean_variance_four_summary :
    gapsDistMeanQ 4 = 3
    ∧ gapsDistVarianceQ 4 = 10 / 9
    ∧ varianceGapsTo 4 = 19 / 16
    ∧ gapsDistVarianceQ 4 < varianceGapsTo 4 :=
  ⟨gapsDistMeanQ_four,
   gapsDistVarianceQ_four,
   varianceGapsTo_four,
   gapsDistVarianceQ_four_lt_varianceGapsTo_four⟩

end PT.Conservation
