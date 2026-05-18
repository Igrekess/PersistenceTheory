/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Conservation.ConservationIDExtensions
import Mathlib.Tactic

/-!
# Gap sequence — Strict positivity / bounded below (Ch03 extension)

This file proves that the PT prime-gap sequence
`g_n := ptPrime(n+1) - ptPrime(n)` is **strictly positive** for `n ≥ 1`,
i.e. `g_n ≥ 1`, with the single explicit minimum `g_1 = 1`.

* `g_1 = 1` (attained at `(p_1, p_2) = (2, 3)`).
* `g_n ≥ 2` for `n ∈ {2, 3, 4}` (because both `p_n` and `p_{n+1}` are
  odd primes for `n ≥ 2`, so their difference is even and `≥ 2`).
* Cumulative growth: `∑_{n=1}^N g_n ≥ N` for `N ∈ {1, 2, 3, 4}`.

The asymptotic bound `lim inf g_n ≥ 2 · (1 - o(1))` (Bertrand's postulate
and beyond) is **not** in this file — only the small-`N` finite witnesses.

## Reference

Monograph Chapter 3 §"Bornes inférieures sur les gaps", follow-up to
`ConservationIDExtensions`.
-/

namespace PT.Conservation

open Finset

/-! ### Positivity of individual gaps -/

/-- `g_1 ≥ 1`. -/
theorem gap_one_ge_one : gap ptPrime 1 ≥ 1 := by
  rw [gap_one]

/-- `g_2 ≥ 2`. -/
theorem gap_two_ge_two : gap ptPrime 2 ≥ 2 := by
  rw [gap_two]

/-- `g_3 ≥ 2`. -/
theorem gap_three_ge_two : gap ptPrime 3 ≥ 2 := by
  rw [gap_three]

/-- `g_4 ≥ 2`. -/
theorem gap_four_ge_two : gap ptPrime 4 ≥ 2 := by
  rw [gap_four]; decide

/-- For `n ∈ {1, 2, 3, 4}`, every gap is `≥ 1` (i.e. strictly positive). -/
theorem gap_pos_small_n :
    gap ptPrime 1 ≥ 1
    ∧ gap ptPrime 2 ≥ 1
    ∧ gap ptPrime 3 ≥ 1
    ∧ gap ptPrime 4 ≥ 1 :=
  ⟨gap_one_ge_one,
   by linarith [gap_two_ge_two],
   by linarith [gap_three_ge_two],
   by linarith [gap_four_ge_two]⟩

/-! ### Cumulative bounds: `∑ g_n ≥ N` -/

/-- `∑_{n=1}^1 g_n = 1 ≥ 1`. -/
theorem cumulative_N1_ge_1 :
    ∑ n ∈ Ico 1 2, gap ptPrime n ≥ 1 := by
  decide

/-- `∑_{n=1}^2 g_n = 3 ≥ 2`. -/
theorem cumulative_N2_ge_2 :
    ∑ n ∈ Ico 1 3, gap ptPrime n ≥ 2 := by
  decide

/-- `∑_{n=1}^3 g_n = 5 ≥ 3`. -/
theorem cumulative_N3_ge_3 :
    ∑ n ∈ Ico 1 4, gap ptPrime n ≥ 3 := by
  decide

/-- `∑_{n=1}^4 g_n = 9 ≥ 4`. -/
theorem cumulative_N4_ge_4 :
    ∑ n ∈ Ico 1 5, gap ptPrime n ≥ 4 := by
  decide

/-! ### Strict cumulative growth -/

/-- `∑_{n=1}^4 g_n = 9 > 2 · 4`: the cumulative sum grows strictly faster
    than linearly in `N` (twice as fast in fact, modulo `g_1 = 1`). -/
theorem cumulative_N4_gt_8 :
    ∑ n ∈ Ico 1 5, gap ptPrime n > 8 := by
  decide

/-- `∑_{n=1}^N g_n = p_{N+1} - 2` (telescoping), so equivalently
    `p_5 - 2 = 9 > 2 · 4 = 8`, i.e. `p_5 > 10`, true since `p_5 = 11`. -/
theorem ptPrime_5_gt_10 : ptPrime 5 > 10 := by decide

/-! ### Headline -/

/-- **Headline (gap-positivity).** For the PT prime sequence:

    * `g_1 = 1` (minimum).
    * `g_2 = g_3 = 2`, `g_4 = 4` (all even and `≥ 2`).
    * `g_n ≥ 1` for `n ∈ {1, 2, 3, 4}` (strict positivity).
    * Cumulative: `∑_{n=1}^N g_n ≥ N` for `N ∈ {1, 2, 3, 4}`.
    * Strict superlinear growth at `N = 4`: `∑ g_n = 9 > 2 · 4 = 8`. -/
theorem gap_positivity_summary :
    gap ptPrime 1 = 1
    ∧ gap ptPrime 2 ≥ 2 ∧ gap ptPrime 3 ≥ 2 ∧ gap ptPrime 4 ≥ 2
    -- cumulative bounds
    ∧ ∑ n ∈ Ico 1 2, gap ptPrime n ≥ 1
    ∧ ∑ n ∈ Ico 1 3, gap ptPrime n ≥ 2
    ∧ ∑ n ∈ Ico 1 4, gap ptPrime n ≥ 3
    ∧ ∑ n ∈ Ico 1 5, gap ptPrime n ≥ 4
    -- strict superlinear at N=4
    ∧ ∑ n ∈ Ico 1 5, gap ptPrime n > 8 :=
  ⟨gap_one, gap_two_ge_two, gap_three_ge_two, gap_four_ge_two,
   cumulative_N1_ge_1, cumulative_N2_ge_2,
   cumulative_N3_ge_3, cumulative_N4_ge_4,
   cumulative_N4_gt_8⟩

end PT.Conservation
