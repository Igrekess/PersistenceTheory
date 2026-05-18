/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Conservation.ConservationIDExtensions
import Mathlib.Algebra.BigOperators.Intervals
import Mathlib.Tactic

/-!
# Gap parity decomposition for the PT prime sequence (Ch03 extension)

The gap sequence on the PT prime sequence
`ptPrime = (2, 3, 5, 7, 11, ...)` is

  `g_1 = 1, g_2 = 2, g_3 = 2, g_4 = 4, ...`

so the **first** gap `g_1 = 1` is the unique odd gap among
`{g_1, g_2, g_3, g_4}`; all other small gaps are even (a manifestation of
the fact that consecutive odd primes differ by an even amount).

This file decomposes the cumulative sum `∑ g_n` into its **odd-index** and
**even-index** parts and gives explicit values for small `N`.

## Reference

Monograph Chapter 3 §"Parité des gaps", follow-up to
`ConservationIDExtensions`.
-/

namespace PT.Conservation

open Finset

/-! ### Odd/even index decomposition for small `N` -/

/-- The sum of odd-indexed gaps (`g_1, g_3`) up to `N = 4`:
    `g_1 + g_3 = 1 + 2 = 3`. -/
theorem oddIndex_sum_N4 :
    gap ptPrime 1 + gap ptPrime 3 = 3 := by
  unfold gap; decide

/-- The sum of even-indexed gaps (`g_2, g_4`) up to `N = 4`:
    `g_2 + g_4 = 2 + 4 = 6`. -/
theorem evenIndex_sum_N4 :
    gap ptPrime 2 + gap ptPrime 4 = 6 := by
  unfold gap; decide

/-- The complete cumulative sum up to `N = 4` is `odd-sum + even-sum
    = 3 + 6 = 9 = p_5 - 2`. -/
theorem cumulative_N4_decomposed :
    gap ptPrime 1 + gap ptPrime 2 + gap ptPrime 3 + gap ptPrime 4
      = (gap ptPrime 1 + gap ptPrime 3) + (gap ptPrime 2 + gap ptPrime 4) := by
  ring

/-! ### Parity of individual gaps -/

/-- `g_1 = 1` is **odd**. -/
theorem gap_one_odd : gap ptPrime 1 % 2 = 1 := by
  rw [gap_one]
  decide

/-- `g_2 = 2` is **even**. -/
theorem gap_two_even : gap ptPrime 2 % 2 = 0 := by
  rw [gap_two]
  decide

/-- `g_3 = 2` is **even**. -/
theorem gap_three_even : gap ptPrime 3 % 2 = 0 := by
  rw [gap_three]
  decide

/-- `g_4 = 4` is **even**. -/
theorem gap_four_even : gap ptPrime 4 % 2 = 0 := by
  rw [gap_four]
  decide

/-! ### Sum-of-gaps parity -/

/-- The cumulative sum `g_1 + g_2 + g_3 + g_4 = 9` is **odd**. -/
theorem cumulative_N4_odd :
    (gap ptPrime 1 + gap ptPrime 2 + gap ptPrime 3 + gap ptPrime 4) % 2 = 1 := by
  rw [gap_one, gap_two, gap_three, gap_four]
  decide

/-- The cumulative sum `g_1 + g_2 + g_3 = 5` is **odd**. -/
theorem cumulative_N3_odd :
    (gap ptPrime 1 + gap ptPrime 2 + gap ptPrime 3) % 2 = 1 := by
  rw [gap_one, gap_two, gap_three]
  decide

/-! ### The unique-odd-gap structural fact -/

/-- **Structural fact.** Among the first four gaps `(g_1, g_2, g_3, g_4)`,
    only `g_1 = 1` is odd. -/
theorem only_g1_is_odd :
    gap ptPrime 1 % 2 = 1
    ∧ gap ptPrime 2 % 2 = 0
    ∧ gap ptPrime 3 % 2 = 0
    ∧ gap ptPrime 4 % 2 = 0 :=
  ⟨gap_one_odd, gap_two_even, gap_three_even, gap_four_even⟩

/-! ### Headline -/

/-- **Headline (gap-parity decomposition).** For the PT prime sequence,
    the first four gaps `(1, 2, 2, 4)` exhibit:

    * `g_1 = 1` is the unique odd gap (parity-anomaly at `(p_1, p_2) = (2, 3)`).
    * `g_2, g_3, g_4 = 2, 2, 4` are all even.
    * Cumulative `∑_{n ≤ 4} g_n = 9 = 3 + 6`, decomposed as odd-indexed
      (`g_1 + g_3 = 3`) + even-indexed (`g_2 + g_4 = 6`).
    * Both cumulative sums at `N = 3` and `N = 4` are odd
      (forced by the single odd `g_1`). -/
theorem gap_parity_decomposition_summary :
    -- parities
    gap ptPrime 1 % 2 = 1
    ∧ gap ptPrime 2 % 2 = 0
    ∧ gap ptPrime 3 % 2 = 0
    ∧ gap ptPrime 4 % 2 = 0
    -- index-decomposition
    ∧ gap ptPrime 1 + gap ptPrime 3 = 3
    ∧ gap ptPrime 2 + gap ptPrime 4 = 6
    -- cumulative-parity
    ∧ (gap ptPrime 1 + gap ptPrime 2 + gap ptPrime 3) % 2 = 1
    ∧ (gap ptPrime 1 + gap ptPrime 2 + gap ptPrime 3 + gap ptPrime 4) % 2 = 1 :=
  ⟨gap_one_odd, gap_two_even, gap_three_even, gap_four_even,
   oddIndex_sum_N4, evenIndex_sum_N4,
   cumulative_N3_odd, cumulative_N4_odd⟩

end PT.Conservation
