/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Conservation.ConservationID
import PT.Conservation.CumulativeBoundsExtended
import PT.Conservation.PrimeGapMoments
import Mathlib.Algebra.BigOperators.Intervals
import Mathlib.Tactic

/-!
# Prime gap moments — extension to `N ∈ {5, …, 10}`

This file extends `PT.Conservation.PrimeGapMoments` (which computes weighted
moments of the prime gaps on `ptPrime`, the first five PT primes) to the
**eleven-prime sequence** `ptPrimeExt` (see `CumulativeBoundsExtended`).

The four moments studied here are, for `N ≥ 1`:

* **Alternating sum** `altSumGapsExt N := ∑_{n=1}^{N} (-1)^{n+1} g_n`,
* **Linear-weighted (first moment)** `linWeightedSumGapsExt N := ∑_{n=1}^{N} n · g_n`,
* **Sum of squares** `sqSumGapsExt N := ∑_{n=1}^{N} g_n²`,
* **Sum of cubes** `cubeSumGapsExt N := ∑_{n=1}^{N} g_n³`,

where `g_n := gap ptPrimeExt n`.

The first ten extended gaps are `(g_1, …, g_{10}) = (1, 2, 2, 4, 2, 4, 2, 4, 6, 2)`,
giving the headline table:

| N  | altSum | linWeighted | sqSum | cubeSum |
|----|--------|-------------|-------|---------|
|  5 |  -1    |    37       |  29   |   89    |
|  6 |  -5    |    61       |  45   |  153    |
|  7 |  -3    |    75       |  49   |  161    |
|  8 |  -7    |   107       |  65   |  225    |
|  9 |  -1    |   161       | 101   |  441    |
| 10 |  -3    |   181       | 105   |  449    |

All values are decidable finite arithmetic in `ℤ`, discharged by `decide`.

## Reference

* `PT.Conservation.ConservationID` — `gap`, `sum_gap_telescope`.
* `PT.Conservation.CumulativeBoundsExtended` — `ptPrimeExt`,
  `gapExt_one … gapExt_ten`, `cumulativeExt_N5 … cumulativeExt_N10`.
* `PT.Conservation.PrimeGapMoments` — the four moments on the base
  five-prime sequence `ptPrime` (this file is the eleven-prime sibling).
* Monograph Chapter 3, §3.1 (extended range).
-/

namespace PT.Conservation

open Finset

/-! ### Definitions on `ptPrimeExt` -/

/-- Alternating sum of the first `N` extended PT prime gaps:
    `∑_{n=1}^{N} (−1)^{n+1} g_n`. -/
def altSumGapsExt (N : ℕ) : ℤ :=
  ∑ n ∈ Ico 1 (N + 1), (-1 : ℤ) ^ (n + 1) * gap ptPrimeExt n

/-- Linear-weighted (first-moment) sum of the first `N` extended PT prime gaps:
    `∑_{n=1}^{N} n · g_n`. -/
def linWeightedSumGapsExt (N : ℕ) : ℤ :=
  ∑ n ∈ Ico 1 (N + 1), (n : ℤ) * gap ptPrimeExt n

/-- Sum of squares of the first `N` extended PT prime gaps:
    `∑_{n=1}^{N} g_n²`. -/
def sqSumGapsExt (N : ℕ) : ℤ :=
  ∑ n ∈ Ico 1 (N + 1), (gap ptPrimeExt n) ^ 2

/-- Sum of cubes of the first `N` extended PT prime gaps:
    `∑_{n=1}^{N} g_n³`. -/
def cubeSumGapsExt (N : ℕ) : ℤ :=
  ∑ n ∈ Ico 1 (N + 1), (gap ptPrimeExt n) ^ 3

/-! ### Alternating sum — values on `N ∈ {5, …, 10}` -/

/-- **Alternating sum at `N = 5`.**
    `1 − 2 + 2 − 4 + 2 = −1`. -/
theorem altSumGapsExt_five : altSumGapsExt 5 = -1 := by
  unfold altSumGapsExt; decide

/-- **Alternating sum at `N = 6`.**
    `1 − 2 + 2 − 4 + 2 − 4 = −5`. -/
theorem altSumGapsExt_six : altSumGapsExt 6 = -5 := by
  unfold altSumGapsExt; decide

/-- **Alternating sum at `N = 7`.**
    `1 − 2 + 2 − 4 + 2 − 4 + 2 = −3`. -/
theorem altSumGapsExt_seven : altSumGapsExt 7 = -3 := by
  unfold altSumGapsExt; decide

/-- **Alternating sum at `N = 8`.**
    `1 − 2 + 2 − 4 + 2 − 4 + 2 − 4 = −7`. -/
theorem altSumGapsExt_eight : altSumGapsExt 8 = -7 := by
  unfold altSumGapsExt; decide

/-- **Alternating sum at `N = 9`.**
    `1 − 2 + 2 − 4 + 2 − 4 + 2 − 4 + 6 = −1`. -/
theorem altSumGapsExt_nine : altSumGapsExt 9 = -1 := by
  unfold altSumGapsExt; decide

/-- **Alternating sum at `N = 10`.**
    `1 − 2 + 2 − 4 + 2 − 4 + 2 − 4 + 6 − 2 = −3`. -/
theorem altSumGapsExt_ten : altSumGapsExt 10 = -3 := by
  unfold altSumGapsExt; decide

/-! ### Linear-weighted (first-moment) sum — values on `N ∈ {5, …, 10}` -/

/-- **Linear-weighted sum at `N = 5`.**
    `1·1 + 2·2 + 3·2 + 4·4 + 5·2 = 1 + 4 + 6 + 16 + 10 = 37`. -/
theorem linWeightedSumGapsExt_five : linWeightedSumGapsExt 5 = 37 := by
  unfold linWeightedSumGapsExt; decide

/-- **Linear-weighted sum at `N = 6`.**
    `37 + 6·4 = 61`. -/
theorem linWeightedSumGapsExt_six : linWeightedSumGapsExt 6 = 61 := by
  unfold linWeightedSumGapsExt; decide

/-- **Linear-weighted sum at `N = 7`.**
    `61 + 7·2 = 75`. -/
theorem linWeightedSumGapsExt_seven : linWeightedSumGapsExt 7 = 75 := by
  unfold linWeightedSumGapsExt; decide

/-- **Linear-weighted sum at `N = 8`.**
    `75 + 8·4 = 107`. -/
theorem linWeightedSumGapsExt_eight : linWeightedSumGapsExt 8 = 107 := by
  unfold linWeightedSumGapsExt; decide

/-- **Linear-weighted sum at `N = 9`.**
    `107 + 9·6 = 161`. -/
theorem linWeightedSumGapsExt_nine : linWeightedSumGapsExt 9 = 161 := by
  unfold linWeightedSumGapsExt; decide

/-- **Linear-weighted sum at `N = 10`.**
    `161 + 10·2 = 181`. -/
theorem linWeightedSumGapsExt_ten : linWeightedSumGapsExt 10 = 181 := by
  unfold linWeightedSumGapsExt; decide

/-! ### Sum of squares — values on `N ∈ {5, …, 10}` -/

/-- **Sum of squares at `N = 5`.** `1 + 4 + 4 + 16 + 4 = 29`. -/
theorem sqSumGapsExt_five : sqSumGapsExt 5 = 29 := by
  unfold sqSumGapsExt; decide

/-- **Sum of squares at `N = 6`.** `29 + 16 = 45`. -/
theorem sqSumGapsExt_six : sqSumGapsExt 6 = 45 := by
  unfold sqSumGapsExt; decide

/-- **Sum of squares at `N = 7`.** `45 + 4 = 49`. -/
theorem sqSumGapsExt_seven : sqSumGapsExt 7 = 49 := by
  unfold sqSumGapsExt; decide

/-- **Sum of squares at `N = 8`.** `49 + 16 = 65`. -/
theorem sqSumGapsExt_eight : sqSumGapsExt 8 = 65 := by
  unfold sqSumGapsExt; decide

/-- **Sum of squares at `N = 9`.** `65 + 36 = 101`. -/
theorem sqSumGapsExt_nine : sqSumGapsExt 9 = 101 := by
  unfold sqSumGapsExt; decide

/-- **Sum of squares at `N = 10`.** `101 + 4 = 105`. -/
theorem sqSumGapsExt_ten : sqSumGapsExt 10 = 105 := by
  unfold sqSumGapsExt; decide

/-! ### Sum of cubes — values on `N ∈ {5, …, 10}` -/

/-- **Sum of cubes at `N = 5`.** `1 + 8 + 8 + 64 + 8 = 89`. -/
theorem cubeSumGapsExt_five : cubeSumGapsExt 5 = 89 := by
  unfold cubeSumGapsExt; decide

/-- **Sum of cubes at `N = 6`.** `89 + 64 = 153`. -/
theorem cubeSumGapsExt_six : cubeSumGapsExt 6 = 153 := by
  unfold cubeSumGapsExt; decide

/-- **Sum of cubes at `N = 7`.** `153 + 8 = 161`. -/
theorem cubeSumGapsExt_seven : cubeSumGapsExt 7 = 161 := by
  unfold cubeSumGapsExt; decide

/-- **Sum of cubes at `N = 8`.** `161 + 64 = 225`. -/
theorem cubeSumGapsExt_eight : cubeSumGapsExt 8 = 225 := by
  unfold cubeSumGapsExt; decide

/-- **Sum of cubes at `N = 9`.** `225 + 216 = 441`. -/
theorem cubeSumGapsExt_nine : cubeSumGapsExt 9 = 441 := by
  unfold cubeSumGapsExt; decide

/-- **Sum of cubes at `N = 10`.** `441 + 8 = 449`. -/
theorem cubeSumGapsExt_ten : cubeSumGapsExt 10 = 449 := by
  unfold cubeSumGapsExt; decide

/-! ### Cross-checks: Cauchy–Schwarz witnesses -/

/-- **Cauchy–Schwarz witness at `N = 5`.** `(∑ g_n)² = 11² = 121 ≤ 5 · 29 = 145`. -/
theorem cauchy_schwarz_witness_ext_five :
    (∑ n ∈ Ico 1 6, gap ptPrimeExt n) ^ 2 ≤ 5 * sqSumGapsExt 5 := by
  rw [cumulativeExt_N5, sqSumGapsExt_five]
  decide

/-- **Cauchy–Schwarz witness at `N = 10`.** `(∑ g_n)² = 29² = 841 ≤ 10 · 105 = 1050`. -/
theorem cauchy_schwarz_witness_ext_ten :
    (∑ n ∈ Ico 1 11, gap ptPrimeExt n) ^ 2 ≤ 10 * sqSumGapsExt 10 := by
  rw [cumulativeExt_N10, sqSumGapsExt_ten]
  decide

/-! ### First-moment vs telescoping sum at `N = 10` -/

/-- **First-moment vs telescoping sum at `N = 10`.**
    The difference between the linear-weighted sum and the plain cumulative
    gap sum is `181 − 29 = 152`. -/
theorem linWeightedSumGapsExt_minus_cumulative_ten :
    linWeightedSumGapsExt 10 - (∑ n ∈ Ico 1 11, gap ptPrimeExt n) = 152 := by
  rw [linWeightedSumGapsExt_ten, cumulativeExt_N10]; rfl

/-! ### Headline summary -/

/-- **Headline (extended-range moments table, `N ∈ {5, …, 10}`).**
    The four weighted moments of the first `N` PT prime gaps, for each
    `N ∈ {5, 6, 7, 8, 9, 10}`, packaged as a single conjunction.

    | N  | altSum | linWeighted | sqSum | cubeSum |
    |----|--------|-------------|-------|---------|
    |  5 |  -1    |    37       |  29   |   89    |
    |  6 |  -5    |    61       |  45   |  153    |
    |  7 |  -3    |    75       |  49   |  161    |
    |  8 |  -7    |   107       |  65   |  225    |
    |  9 |  -1    |   161       | 101   |  441    |
    | 10 |  -3    |   181       | 105   |  449    | -/
theorem prime_gap_moments_ext_headline :
    -- N = 5
    altSumGapsExt 5 = -1
    ∧ linWeightedSumGapsExt 5 = 37
    ∧ sqSumGapsExt 5 = 29
    ∧ cubeSumGapsExt 5 = 89
    -- N = 6
    ∧ altSumGapsExt 6 = -5
    ∧ linWeightedSumGapsExt 6 = 61
    ∧ sqSumGapsExt 6 = 45
    ∧ cubeSumGapsExt 6 = 153
    -- N = 7
    ∧ altSumGapsExt 7 = -3
    ∧ linWeightedSumGapsExt 7 = 75
    ∧ sqSumGapsExt 7 = 49
    ∧ cubeSumGapsExt 7 = 161
    -- N = 8
    ∧ altSumGapsExt 8 = -7
    ∧ linWeightedSumGapsExt 8 = 107
    ∧ sqSumGapsExt 8 = 65
    ∧ cubeSumGapsExt 8 = 225
    -- N = 9
    ∧ altSumGapsExt 9 = -1
    ∧ linWeightedSumGapsExt 9 = 161
    ∧ sqSumGapsExt 9 = 101
    ∧ cubeSumGapsExt 9 = 441
    -- N = 10
    ∧ altSumGapsExt 10 = -3
    ∧ linWeightedSumGapsExt 10 = 181
    ∧ sqSumGapsExt 10 = 105
    ∧ cubeSumGapsExt 10 = 449 :=
  ⟨altSumGapsExt_five,   linWeightedSumGapsExt_five,   sqSumGapsExt_five,   cubeSumGapsExt_five,
   altSumGapsExt_six,    linWeightedSumGapsExt_six,    sqSumGapsExt_six,    cubeSumGapsExt_six,
   altSumGapsExt_seven,  linWeightedSumGapsExt_seven,  sqSumGapsExt_seven,  cubeSumGapsExt_seven,
   altSumGapsExt_eight,  linWeightedSumGapsExt_eight,  sqSumGapsExt_eight,  cubeSumGapsExt_eight,
   altSumGapsExt_nine,   linWeightedSumGapsExt_nine,   sqSumGapsExt_nine,   cubeSumGapsExt_nine,
   altSumGapsExt_ten,    linWeightedSumGapsExt_ten,    sqSumGapsExt_ten,    cubeSumGapsExt_ten⟩

end PT.Conservation
