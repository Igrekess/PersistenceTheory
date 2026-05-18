/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Conservation.ConservationID
import PT.Conservation.ConservationIDExtensions
import Mathlib.Algebra.BigOperators.Intervals
import Mathlib.Tactic

/-!
# Prime gap moments — weighted sums on the PT prime sequence

This file computes **weighted moments** of the prime gaps
`g_n := p_{n+1} - p_n` on the PT prime sequence `ptPrime`
(see `PT.Conservation.ConservationIDExtensions`), for small `N`.

The four moments studied here are, for `N ≥ 1`:

* **Alternating sum** `altSumGaps N := ∑_{n=1}^{N} (-1)^{n+1} g_n`,
* **Linear-weighted (first moment)** `linWeightedSumGaps N := ∑_{n=1}^{N} n · g_n`,
* **Sum of squares** `sqSumGaps N := ∑_{n=1}^{N} g_n²`,
* **Sum of cubes** `cubeSumGaps N := ∑_{n=1}^{N} g_n³`.

With `(g_1, g_2, g_3, g_4) = (1, 2, 2, 4)` we obtain the concrete values

| moment           | `N = 4` value                         |
|------------------|---------------------------------------|
| `altSumGaps`     | `1 − 2 + 2 − 4 = −3`                  |
| `linWeightedSum` | `1·1 + 2·2 + 3·2 + 4·4 = 27`          |
| `sqSumGaps`      | `1 + 4 + 4 + 16 = 25`                 |
| `cubeSumGaps`    | `1 + 8 + 8 + 64 = 81`                 |

We also record a **first-moment identity** that compares the linear-weighted
sum to the plain telescoping sum (which, by the conservation identity, equals
`p_{N+1} − 2 = 9` at `N = 4`).

All values are decidable finite arithmetic in `ℤ`, discharged by `decide`.

## Reference

* `PT.Conservation.ConservationID` — `gap`, `sum_gap_telescope`.
* `PT.Conservation.ConservationIDExtensions` — `ptPrime`,
  `gap_one, gap_two, gap_three, gap_four`, `conservation_N4`.
* Monograph Chapter 3, §3.1.
-/

namespace PT.Conservation

open Finset

/-! ### Definitions -/

/-- Alternating sum of the first `N` PT prime gaps:
    `∑_{n=1}^{N} (−1)^{n+1} g_n`. -/
def altSumGaps (N : ℕ) : ℤ :=
  ∑ n ∈ Ico 1 (N + 1), (-1 : ℤ) ^ (n + 1) * gap ptPrime n

/-- Linear-weighted (first-moment) sum of the first `N` PT prime gaps:
    `∑_{n=1}^{N} n · g_n`. -/
def linWeightedSumGaps (N : ℕ) : ℤ :=
  ∑ n ∈ Ico 1 (N + 1), (n : ℤ) * gap ptPrime n

/-- Sum of squares of the first `N` PT prime gaps:
    `∑_{n=1}^{N} g_n²`. -/
def sqSumGaps (N : ℕ) : ℤ :=
  ∑ n ∈ Ico 1 (N + 1), (gap ptPrime n) ^ 2

/-- Sum of cubes of the first `N` PT prime gaps:
    `∑_{n=1}^{N} g_n³`. -/
def cubeSumGaps (N : ℕ) : ℤ :=
  ∑ n ∈ Ico 1 (N + 1), (gap ptPrime n) ^ 3

/-! ### Alternating sum — small-`N` values -/

/-- `altSumGaps 1 = g_1 = 1`. -/
theorem altSumGaps_one : altSumGaps 1 = 1 := by
  unfold altSumGaps; decide

/-- `altSumGaps 2 = g_1 − g_2 = 1 − 2 = −1`. -/
theorem altSumGaps_two : altSumGaps 2 = -1 := by
  unfold altSumGaps; decide

/-- `altSumGaps 3 = g_1 − g_2 + g_3 = 1 − 2 + 2 = 1`. -/
theorem altSumGaps_three : altSumGaps 3 = 1 := by
  unfold altSumGaps; decide

/-- **Alternating sum at `N = 4`.**
    `g_1 − g_2 + g_3 − g_4 = 1 − 2 + 2 − 4 = −3`. -/
theorem altSumGaps_four : altSumGaps 4 = -3 := by
  unfold altSumGaps; decide

/-! ### Linear-weighted (first-moment) sum — small-`N` values -/

/-- `linWeightedSumGaps 1 = 1 · g_1 = 1`. -/
theorem linWeightedSumGaps_one : linWeightedSumGaps 1 = 1 := by
  unfold linWeightedSumGaps; decide

/-- `linWeightedSumGaps 2 = 1·1 + 2·2 = 5`. -/
theorem linWeightedSumGaps_two : linWeightedSumGaps 2 = 5 := by
  unfold linWeightedSumGaps; decide

/-- `linWeightedSumGaps 3 = 1·1 + 2·2 + 3·2 = 11`. -/
theorem linWeightedSumGaps_three : linWeightedSumGaps 3 = 11 := by
  unfold linWeightedSumGaps; decide

/-- **Linear-weighted (first-moment) sum at `N = 4`.**
    `1·1 + 2·2 + 3·2 + 4·4 = 1 + 4 + 6 + 16 = 27`. -/
theorem linWeightedSumGaps_four : linWeightedSumGaps 4 = 27 := by
  unfold linWeightedSumGaps; decide

/-! ### Sum of squares — small-`N` values -/

/-- `sqSumGaps 1 = 1`. -/
theorem sqSumGaps_one : sqSumGaps 1 = 1 := by
  unfold sqSumGaps; decide

/-- `sqSumGaps 2 = 1 + 4 = 5`. -/
theorem sqSumGaps_two : sqSumGaps 2 = 5 := by
  unfold sqSumGaps; decide

/-- `sqSumGaps 3 = 1 + 4 + 4 = 9`. -/
theorem sqSumGaps_three : sqSumGaps 3 = 9 := by
  unfold sqSumGaps; decide

/-- **Sum of squares at `N = 4`.** `1 + 4 + 4 + 16 = 25`. -/
theorem sqSumGaps_four : sqSumGaps 4 = 25 := by
  unfold sqSumGaps; decide

/-! ### Sum of cubes — small-`N` values -/

/-- `cubeSumGaps 1 = 1`. -/
theorem cubeSumGaps_one : cubeSumGaps 1 = 1 := by
  unfold cubeSumGaps; decide

/-- `cubeSumGaps 2 = 1 + 8 = 9`. -/
theorem cubeSumGaps_two : cubeSumGaps 2 = 9 := by
  unfold cubeSumGaps; decide

/-- `cubeSumGaps 3 = 1 + 8 + 8 = 17`. -/
theorem cubeSumGaps_three : cubeSumGaps 3 = 17 := by
  unfold cubeSumGaps; decide

/-- **Sum of cubes at `N = 4`.** `1 + 8 + 8 + 64 = 81`. -/
theorem cubeSumGaps_four : cubeSumGaps 4 = 81 := by
  unfold cubeSumGaps; decide

/-! ### First-moment identity vs the conservation telescoping sum

The conservation identity `conservation_N4` gives
`∑_{n=1}^{4} g_n = p_5 − 2 = 9`. The linear-weighted sum
`linWeightedSumGaps 4 = 27 = 3 · 9` happens, at this small `N`, to equal
`3` times the telescoping sum. We record both the difference and the
ratio identity at `N = 4`.
-/

/-- **First-moment vs telescoping sum at `N = 4`.**
    The difference between the linear-weighted sum and the plain
    cumulative gap sum is `27 − 9 = 18`. -/
theorem linWeightedSumGaps_minus_cumulative_four :
    linWeightedSumGaps 4 - (∑ n ∈ Ico 1 5, gap ptPrime n) = 18 := by
  rw [linWeightedSumGaps_four, conservation_N4]; rfl

/-- **First-moment ratio identity at `N = 4`.**
    `∑_{n=1}^{4} n · g_n = 3 · ∑_{n=1}^{4} g_n`. Both sides equal `27`. -/
theorem linWeightedSumGaps_eq_three_cumulative_four :
    linWeightedSumGaps 4 = 3 * (∑ n ∈ Ico 1 5, gap ptPrime n) := by
  rw [linWeightedSumGaps_four, conservation_N4]; rfl

/-! ### Cross-checks: sum of squares vs square of sum

For positive integer values `g_1, …, g_N` we have, by Cauchy–Schwarz,
`(∑ g_n)² ≤ N · ∑ g_n²`. At `N = 4` this reads `9² = 81 ≤ 4 · 25 = 100`,
which we record as a finite numerical witness. -/

/-- **Cauchy–Schwarz witness at `N = 4`.** `9² ≤ 4 · 25`. -/
theorem cauchy_schwarz_witness_four :
    (∑ n ∈ Ico 1 5, gap ptPrime n) ^ 2 ≤ 4 * sqSumGaps 4 := by
  rw [conservation_N4, sqSumGaps_four]
  decide

/-! ### Headline summary -/

/-- **Headline (`N = 4` moments summary).**
    The four weighted moments of the first four PT prime gaps:

    * `altSumGaps 4       = −3`
    * `linWeightedSumGaps 4 = 27`
    * `sqSumGaps 4        = 25`
    * `cubeSumGaps 4      = 81` -/
theorem prime_gap_moments_N4_summary :
    altSumGaps 4 = -3
    ∧ linWeightedSumGaps 4 = 27
    ∧ sqSumGaps 4 = 25
    ∧ cubeSumGaps 4 = 81 :=
  ⟨altSumGaps_four, linWeightedSumGaps_four,
   sqSumGaps_four, cubeSumGaps_four⟩

end PT.Conservation
