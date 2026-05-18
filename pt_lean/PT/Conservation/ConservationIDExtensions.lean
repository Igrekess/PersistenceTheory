/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Conservation.ConservationID
import Mathlib.Algebra.BigOperators.Intervals
import Mathlib.Tactic

/-!
# Conservation Identity — algebraic extensions

This file collects extension lemmas for the conservation identity
`PT.Conservation.ConservationID`, namely:

* Telescoping with the *left-shifted* sum convention
  `∑_{n=0}^{N-1} g_n = p N - p 0`.
* Recovery of a single gap from successive cumulative sums.
* Monotonicity of the cumulative sum when `p` is monotone.
* Specialisation to PT prime indexing: small-`N` values when
  `p_1 = 2, p_2 = 3, p_3 = 5, p_4 = 7`.

All proofs are short telescoping / `induction` arguments resting on
`sum_gap_telescope` (already kernel-verified).

## Reference

Monograph Chapter 3, §3.1, follow-up to `\label{thm:conservation-id}`.
Audit follow-up to Vague 4 Phase 1 EASY.
-/

namespace PT.Conservation

open Finset

/-- **Left-shifted conservation identity.**
    For any sequence `p : ℕ → ℤ` and any `N`,
    `∑_{n=0}^{N-1} (p(n+1) - p n) = p N - p 0`. -/
theorem sum_gap_telescope_zero (p : ℕ → ℤ) (N : ℕ) :
    ∑ n ∈ range N, gap p n = p N - p 0 := by
  induction N with
  | zero => simp [gap]
  | succ k ih =>
      rw [Finset.sum_range_succ, ih]
      show p k - p 0 + gap p k = p (k + 1) - p 0
      unfold gap
      ring

/-- **Single-gap recovery.** The `(N)`-th gap is the difference of the two
    successive cumulative sums on `Ico 1 (N+2)` and `Ico 1 (N+1)`. -/
theorem gap_eq_sum_diff (p : ℕ → ℤ) (N : ℕ) :
    gap p N = (∑ n ∈ range (N + 1), gap p n)
            - (∑ n ∈ range N, gap p n) := by
  rw [Finset.sum_range_succ]
  ring

/-! ### Monotonicity -/

/-- If `p` is monotone (each step `p(n+1) ≥ p n`), then every gap is nonnegative. -/
theorem gap_nonneg_of_mono (p : ℕ → ℤ) (hmono : ∀ n, p n ≤ p (n + 1)) (n : ℕ) :
    0 ≤ gap p n := by
  unfold gap
  linarith [hmono n]

/-- If `p` is monotone, the cumulative sum is non-decreasing in `N`. -/
theorem cumulative_sum_mono (p : ℕ → ℤ) (hmono : ∀ n, p n ≤ p (n + 1)) :
    ∀ N, ∑ n ∈ range N, gap p n ≤ ∑ n ∈ range (N + 1), gap p n := by
  intro N
  rw [Finset.sum_range_succ]
  linarith [gap_nonneg_of_mono p hmono N]

/-- If `p` is monotone, the cumulative sum on `Ico 1 (N+1)` is non-negative. -/
theorem cumulative_sum_nonneg (p : ℕ → ℤ) (hmono : ∀ n, p n ≤ p (n + 1)) (N : ℕ) :
    0 ≤ ∑ n ∈ Ico 1 (N + 1), gap p n := by
  induction N with
  | zero => simp
  | succ k ih =>
      rw [show k + 1 + 1 = (k + 1) + 1 from rfl,
          Finset.sum_Ico_succ_top (by omega : 1 ≤ k + 1)]
      linarith [gap_nonneg_of_mono p hmono (k + 1)]

/-! ### PT specialisation: small `N` values -/

/-- The PT prime sequence (integer-valued) for the first four primes. -/
def ptPrime : ℕ → ℤ
  | 1 => 2
  | 2 => 3
  | 3 => 5
  | 4 => 7
  | 5 => 11
  | _ => 0  -- unused beyond first 5 in this file

@[simp] theorem ptPrime_one : ptPrime 1 = 2 := rfl
@[simp] theorem ptPrime_two : ptPrime 2 = 3 := rfl
@[simp] theorem ptPrime_three : ptPrime 3 = 5 := rfl
@[simp] theorem ptPrime_four : ptPrime 4 = 7 := rfl
@[simp] theorem ptPrime_five : ptPrime 5 = 11 := rfl

/-- First gap: `g_1 = p_2 - p_1 = 1`. -/
theorem gap_one : gap ptPrime 1 = 1 := by
  unfold gap; decide

/-- Second gap: `g_2 = p_3 - p_2 = 2`. -/
theorem gap_two : gap ptPrime 2 = 2 := by
  unfold gap; decide

/-- Third gap: `g_3 = p_4 - p_3 = 2`. -/
theorem gap_three : gap ptPrime 3 = 2 := by
  unfold gap; decide

/-- Fourth gap: `g_4 = p_5 - p_4 = 4`. -/
theorem gap_four : gap ptPrime 4 = 4 := by
  unfold gap; decide

/-- **Conservation identity, `N = 1` instance.** `g_1 = p_2 - 2 = 1`. -/
theorem conservation_N1 : ∑ n ∈ Ico 1 2, gap ptPrime n = 1 := by
  decide

/-- **Conservation identity, `N = 2` instance.** `g_1 + g_2 = p_3 - 2 = 3`. -/
theorem conservation_N2 : ∑ n ∈ Ico 1 3, gap ptPrime n = 3 := by
  decide

/-- **Conservation identity, `N = 3` instance.** `g_1 + g_2 + g_3 = p_4 - 2 = 5`. -/
theorem conservation_N3 : ∑ n ∈ Ico 1 4, gap ptPrime n = 5 := by
  decide

/-- **Conservation identity, `N = 4` instance.** `g_1 + … + g_4 = p_5 - 2 = 9`. -/
theorem conservation_N4 : ∑ n ∈ Ico 1 5, gap ptPrime n = 9 := by
  decide

/-- **Cross-check (`N = 4`).** The four-prime cumulative sum equals `p_5 - p_1 = 9`,
    as predicted by the generic telescoping identity. -/
theorem conservation_N4_via_generic :
    ∑ n ∈ Ico 1 5, gap ptPrime n = ptPrime 5 - ptPrime 1 := by
  have := sum_gap_telescope ptPrime 4
  -- The Ico-form telescope yields `p 5 - p 1`.
  simpa using this

end PT.Conservation
