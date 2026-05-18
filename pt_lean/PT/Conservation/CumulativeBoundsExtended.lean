/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Conservation.ConservationID
import PT.Conservation.ConservationIDExtensions
import PT.Conservation.ConservationIDPrimorial
import PT.Conservation.GapBoundedBelow
import Mathlib.Tactic

/-!
# Cumulative gap bounds — extension to `N ∈ {5, …, 10}`

This file extends `PT.Conservation.GapBoundedBelow` to the next six prime
indices. The base file only goes up to `N = 4` because the symbolic
`ptPrime` defined in `ConservationIDExtensions` only lists `(2, 3, 5, 7, 11)`.

Here we introduce **a new sequence** `ptPrimeExt : ℕ → ℤ` covering the
first eleven primes:

  `(2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31)`   for indices `1, …, 11`.

This is strictly an extension; we do **not** modify the existing
`ptPrime`. The two sequences agree on indices `{1, …, 5}`.

From this we derive:

* The individual gaps `g_n` for `n ∈ {1, …, 10}` (gaps `(1, 2, 2, 4, 2, 4,
  2, 4, 6, 2)`).
* The cumulative sums `∑_{n=1}^N g_n` for `N ∈ {5, …, 10}`
  (values `(11, 15, 17, 21, 27, 29) = p_{N+1} - 2`).
* The **superlinear bound** `∑ g_n ≥ N + 1` (witnessed in the extended
  range — already strict at `N = 4` since `9 ≥ 5`).
* The **primorial-relative bound** `∑ g_n ≤ primorial(N+1) / N`
  (very loose; trivially decidable from the values above).

These are all decidable finite instances. The deeper Bertrand /
Mertens-type asymptotics are out of scope.

## Reference

Monograph Chapter 3 §"Bornes inférieures sur les gaps" (extension),
follow-up to `GapBoundedBelow.lean`.
-/

namespace PT.Conservation

open Finset

/-! ### Extended prime sequence -/

/-- The PT prime sequence (integer-valued) for the **first eleven primes**.
    This is a strict extension of `ptPrime` (which only covers `1 … 5`):
    the two sequences agree on `{1, 2, 3, 4, 5}` and `ptPrimeExt` adds the
    six further values `13, 17, 19, 23, 29, 31`. -/
def ptPrimeExt : ℕ → ℤ
  | 1  => 2
  | 2  => 3
  | 3  => 5
  | 4  => 7
  | 5  => 11
  | 6  => 13
  | 7  => 17
  | 8  => 19
  | 9  => 23
  | 10 => 29
  | 11 => 31
  | _  => 0  -- unused beyond first 11 in this file

@[simp] theorem ptPrimeExt_one    : ptPrimeExt 1  = 2  := rfl
@[simp] theorem ptPrimeExt_two    : ptPrimeExt 2  = 3  := rfl
@[simp] theorem ptPrimeExt_three  : ptPrimeExt 3  = 5  := rfl
@[simp] theorem ptPrimeExt_four   : ptPrimeExt 4  = 7  := rfl
@[simp] theorem ptPrimeExt_five   : ptPrimeExt 5  = 11 := rfl
@[simp] theorem ptPrimeExt_six    : ptPrimeExt 6  = 13 := rfl
@[simp] theorem ptPrimeExt_seven  : ptPrimeExt 7  = 17 := rfl
@[simp] theorem ptPrimeExt_eight  : ptPrimeExt 8  = 19 := rfl
@[simp] theorem ptPrimeExt_nine   : ptPrimeExt 9  = 23 := rfl
@[simp] theorem ptPrimeExt_ten    : ptPrimeExt 10 = 29 := rfl
@[simp] theorem ptPrimeExt_eleven : ptPrimeExt 11 = 31 := rfl

/-- The extended sequence agrees with `ptPrime` on `{1, …, 5}`. -/
theorem ptPrimeExt_eq_ptPrime_small :
    ptPrimeExt 1 = ptPrime 1
    ∧ ptPrimeExt 2 = ptPrime 2
    ∧ ptPrimeExt 3 = ptPrime 3
    ∧ ptPrimeExt 4 = ptPrime 4
    ∧ ptPrimeExt 5 = ptPrime 5 := by
  refine ⟨rfl, rfl, rfl, rfl, rfl⟩

/-! ### Individual gaps on the extended range -/

/-- `g_1 = p_2 - p_1 = 1`. -/
theorem gapExt_one   : gap ptPrimeExt 1  = 1 := by unfold gap; decide
/-- `g_2 = p_3 - p_2 = 2`. -/
theorem gapExt_two   : gap ptPrimeExt 2  = 2 := by unfold gap; decide
/-- `g_3 = p_4 - p_3 = 2`. -/
theorem gapExt_three : gap ptPrimeExt 3  = 2 := by unfold gap; decide
/-- `g_4 = p_5 - p_4 = 4`. -/
theorem gapExt_four  : gap ptPrimeExt 4  = 4 := by unfold gap; decide
/-- `g_5 = p_6 - p_5 = 2`. -/
theorem gapExt_five  : gap ptPrimeExt 5  = 2 := by unfold gap; decide
/-- `g_6 = p_7 - p_6 = 4`. -/
theorem gapExt_six   : gap ptPrimeExt 6  = 4 := by unfold gap; decide
/-- `g_7 = p_8 - p_7 = 2`. -/
theorem gapExt_seven : gap ptPrimeExt 7  = 2 := by unfold gap; decide
/-- `g_8 = p_9 - p_8 = 4`. -/
theorem gapExt_eight : gap ptPrimeExt 8  = 4 := by unfold gap; decide
/-- `g_9 = p_10 - p_9 = 6`. -/
theorem gapExt_nine  : gap ptPrimeExt 9  = 6 := by unfold gap; decide
/-- `g_10 = p_11 - p_10 = 2`. -/
theorem gapExt_ten   : gap ptPrimeExt 10 = 2 := by unfold gap; decide

/-! ### Cumulative sums on the extended range -/

/-- `∑_{n=1}^5 g_n = p_6 - 2 = 11`. -/
theorem cumulativeExt_N5 :
    ∑ n ∈ Ico 1 6, gap ptPrimeExt n = 11 := by decide

/-- `∑_{n=1}^6 g_n = p_7 - 2 = 15`. -/
theorem cumulativeExt_N6 :
    ∑ n ∈ Ico 1 7, gap ptPrimeExt n = 15 := by decide

/-- `∑_{n=1}^7 g_n = p_8 - 2 = 17`. -/
theorem cumulativeExt_N7 :
    ∑ n ∈ Ico 1 8, gap ptPrimeExt n = 17 := by decide

/-- `∑_{n=1}^8 g_n = p_9 - 2 = 21`. -/
theorem cumulativeExt_N8 :
    ∑ n ∈ Ico 1 9, gap ptPrimeExt n = 21 := by decide

/-- `∑_{n=1}^9 g_n = p_10 - 2 = 27`. -/
theorem cumulativeExt_N9 :
    ∑ n ∈ Ico 1 10, gap ptPrimeExt n = 27 := by decide

/-- `∑_{n=1}^{10} g_n = p_11 - 2 = 29`. -/
theorem cumulativeExt_N10 :
    ∑ n ∈ Ico 1 11, gap ptPrimeExt n = 29 := by decide

/-! ### Cross-check via the generic telescoping identity -/

/-- Generic telescoping cross-check at `N = 5`: `∑ g_n = p_6 - p_1 = 11`. -/
theorem cumulativeExt_N5_via_generic :
    ∑ n ∈ Ico 1 6, gap ptPrimeExt n = ptPrimeExt 6 - ptPrimeExt 1 := by
  simpa using sum_gap_telescope ptPrimeExt 5

/-- Generic telescoping cross-check at `N = 10`: `∑ g_n = p_11 - p_1 = 29`. -/
theorem cumulativeExt_N10_via_generic :
    ∑ n ∈ Ico 1 11, gap ptPrimeExt n = ptPrimeExt 11 - ptPrimeExt 1 := by
  simpa using sum_gap_telescope ptPrimeExt 10

/-! ### Superlinear lower bound `∑ g_n ≥ N + 1` (extended range) -/

/-- `∑_{n=1}^5 g_n = 11 ≥ 6`. -/
theorem cumulativeExt_N5_ge :
    ∑ n ∈ Ico 1 6, gap ptPrimeExt n ≥ 6 := by decide

/-- `∑_{n=1}^6 g_n = 15 ≥ 7`. -/
theorem cumulativeExt_N6_ge :
    ∑ n ∈ Ico 1 7, gap ptPrimeExt n ≥ 7 := by decide

/-- `∑_{n=1}^7 g_n = 17 ≥ 8`. -/
theorem cumulativeExt_N7_ge :
    ∑ n ∈ Ico 1 8, gap ptPrimeExt n ≥ 8 := by decide

/-- `∑_{n=1}^8 g_n = 21 ≥ 9`. -/
theorem cumulativeExt_N8_ge :
    ∑ n ∈ Ico 1 9, gap ptPrimeExt n ≥ 9 := by decide

/-- `∑_{n=1}^9 g_n = 27 ≥ 10`. -/
theorem cumulativeExt_N9_ge :
    ∑ n ∈ Ico 1 10, gap ptPrimeExt n ≥ 10 := by decide

/-- `∑_{n=1}^{10} g_n = 29 ≥ 11`. -/
theorem cumulativeExt_N10_ge :
    ∑ n ∈ Ico 1 11, gap ptPrimeExt n ≥ 11 := by decide

/-! ### Extended primorials -/

/-- The first eleven primorials, as integers.
    `primorialExt k = ∏_{i = 1}^{k} p_i`. -/
def primorialExt : ℕ → ℤ
  | 0  => 1
  | 1  => 2
  | 2  => 6
  | 3  => 30
  | 4  => 210
  | 5  => 2310
  | 6  => 30030
  | 7  => 510510
  | 8  => 9699690
  | 9  => 223092870
  | 10 => 6469693230
  | 11 => 200560490130
  | _  => 0  -- unused

@[simp] theorem primorialExt_6  : primorialExt 6  = 30030        := rfl
@[simp] theorem primorialExt_7  : primorialExt 7  = 510510       := rfl
@[simp] theorem primorialExt_8  : primorialExt 8  = 9699690      := rfl
@[simp] theorem primorialExt_9  : primorialExt 9  = 223092870    := rfl
@[simp] theorem primorialExt_10 : primorialExt 10 = 6469693230   := rfl
@[simp] theorem primorialExt_11 : primorialExt 11 = 200560490130 := rfl

/-! ### Primorial-relative upper bound `∑ g_n ≤ primorial(N+1) / N` -/

/-- `∑_{n=1}^5 g_n = 11 ≤ primorialExt 6 / 5 = 6006`. -/
theorem cumulativeExt_N5_le_primorial :
    ∑ n ∈ Ico 1 6, gap ptPrimeExt n ≤ primorialExt 6 / 5 := by decide

/-- `∑_{n=1}^6 g_n = 15 ≤ primorialExt 7 / 6 = 85085`. -/
theorem cumulativeExt_N6_le_primorial :
    ∑ n ∈ Ico 1 7, gap ptPrimeExt n ≤ primorialExt 7 / 6 := by decide

/-- `∑_{n=1}^7 g_n = 17 ≤ primorialExt 8 / 7 = 1385670`. -/
theorem cumulativeExt_N7_le_primorial :
    ∑ n ∈ Ico 1 8, gap ptPrimeExt n ≤ primorialExt 8 / 7 := by decide

/-- `∑_{n=1}^8 g_n = 21 ≤ primorialExt 9 / 8`. -/
theorem cumulativeExt_N8_le_primorial :
    ∑ n ∈ Ico 1 9, gap ptPrimeExt n ≤ primorialExt 9 / 8 := by decide

/-- `∑_{n=1}^9 g_n = 27 ≤ primorialExt 10 / 9`. -/
theorem cumulativeExt_N9_le_primorial :
    ∑ n ∈ Ico 1 10, gap ptPrimeExt n ≤ primorialExt 10 / 9 := by decide

/-- `∑_{n=1}^{10} g_n = 29 ≤ primorialExt 11 / 10`. -/
theorem cumulativeExt_N10_le_primorial :
    ∑ n ∈ Ico 1 11, gap ptPrimeExt n ≤ primorialExt 11 / 10 := by decide

/-! ### Headline summaries -/

/-- **Headline (extended cumulative bounds).** For the eleven-prime PT
    sequence `ptPrimeExt`, the cumulative gap sums on `N ∈ {5, …, 10}`
    are exactly `(11, 15, 17, 21, 27, 29)`, and each is bounded below
    by `N + 1` and above by `primorialExt(N+1) / N`. -/
theorem cumulativeExt_headline :
    -- exact values
    (∑ n ∈ Ico 1 6,  gap ptPrimeExt n = 11)
    ∧ (∑ n ∈ Ico 1 7,  gap ptPrimeExt n = 15)
    ∧ (∑ n ∈ Ico 1 8,  gap ptPrimeExt n = 17)
    ∧ (∑ n ∈ Ico 1 9,  gap ptPrimeExt n = 21)
    ∧ (∑ n ∈ Ico 1 10, gap ptPrimeExt n = 27)
    ∧ (∑ n ∈ Ico 1 11, gap ptPrimeExt n = 29)
    -- superlinear lower bound `≥ N + 1`
    ∧ (∑ n ∈ Ico 1 6,  gap ptPrimeExt n ≥ 6)
    ∧ (∑ n ∈ Ico 1 7,  gap ptPrimeExt n ≥ 7)
    ∧ (∑ n ∈ Ico 1 8,  gap ptPrimeExt n ≥ 8)
    ∧ (∑ n ∈ Ico 1 9,  gap ptPrimeExt n ≥ 9)
    ∧ (∑ n ∈ Ico 1 10, gap ptPrimeExt n ≥ 10)
    ∧ (∑ n ∈ Ico 1 11, gap ptPrimeExt n ≥ 11)
    -- primorial-relative upper bound `≤ primorial(N+1) / N`
    ∧ (∑ n ∈ Ico 1 6,  gap ptPrimeExt n ≤ primorialExt 6  / 5)
    ∧ (∑ n ∈ Ico 1 7,  gap ptPrimeExt n ≤ primorialExt 7  / 6)
    ∧ (∑ n ∈ Ico 1 8,  gap ptPrimeExt n ≤ primorialExt 8  / 7)
    ∧ (∑ n ∈ Ico 1 9,  gap ptPrimeExt n ≤ primorialExt 9  / 8)
    ∧ (∑ n ∈ Ico 1 10, gap ptPrimeExt n ≤ primorialExt 10 / 9)
    ∧ (∑ n ∈ Ico 1 11, gap ptPrimeExt n ≤ primorialExt 11 / 10) :=
  ⟨cumulativeExt_N5,  cumulativeExt_N6,  cumulativeExt_N7,
   cumulativeExt_N8,  cumulativeExt_N9,  cumulativeExt_N10,
   cumulativeExt_N5_ge,  cumulativeExt_N6_ge,  cumulativeExt_N7_ge,
   cumulativeExt_N8_ge,  cumulativeExt_N9_ge,  cumulativeExt_N10_ge,
   cumulativeExt_N5_le_primorial,  cumulativeExt_N6_le_primorial,
   cumulativeExt_N7_le_primorial,  cumulativeExt_N8_le_primorial,
   cumulativeExt_N9_le_primorial,  cumulativeExt_N10_le_primorial⟩

end PT.Conservation
