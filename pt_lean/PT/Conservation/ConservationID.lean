/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import Mathlib.Algebra.BigOperators.Intervals
import Mathlib.Tactic

/-!
# Conservation Identity — `∑ g_n = p_{N+1} - 2`

**Statement (paper-level, Ch03 §3.1).** Let `p : ℕ → ℕ` enumerate the prime
sequence with `p 1 = 2`, and let `g_n := p (n+1) - p n` denote the `n`-th
prime gap. Then for every `N ≥ 1`,

$$\sum_{n=1}^{N} g_n = p_{N+1} - p_1 = p_{N+1} - 2.$$

This is a pure telescoping identity: it does not depend on the primality of
`p`, only on `p 1 = 2` and the monotonicity convention. We formalise it
generically for any sequence `p : ℕ → ℤ` and then specialise.

## Reference

Monograph Chapter 3, §3.1, `\label{thm:conservation-id}`. M1 article
(PT_ARTICLES/PT_MATHEMATICS/M1), Conservation identity.
-/

namespace PT.Conservation

open Finset

/-- The `n`-th gap of an integer-valued sequence `p`: `g_n := p(n+1) - p(n)`. -/
def gap (p : ℕ → ℤ) (n : ℕ) : ℤ := p (n + 1) - p n

/-- **Conservation identity (generic telescoping form).**
    For any sequence `p : ℕ → ℤ` and any `N`,
    `∑_{n=1}^{N} (p(n+1) - p n) = p(N+1) - p 1`. -/
theorem sum_gap_telescope (p : ℕ → ℤ) (N : ℕ) :
    ∑ n ∈ Ico 1 (N + 1), gap p n = p (N + 1) - p 1 := by
  induction N with
  | zero => simp [gap]
  | succ k ih =>
      rw [show k + 1 + 1 = (k + 1) + 1 from rfl,
          Finset.sum_Ico_succ_top (by omega : 1 ≤ k + 1), ih]
      show p (k + 1) - p 1 + gap p (k + 1) = p (k + 1 + 1) - p 1
      unfold gap
      ring

/-- **Conservation identity (PT form, integer prime convention).**
    If `p : ℕ → ℤ` satisfies `p 1 = 2`, then for any `N`,
    `∑_{n=1}^{N} g_n = p(N+1) - 2`. -/
theorem conservation_id (p : ℕ → ℤ) (hp1 : p 1 = 2) (N : ℕ) :
    ∑ n ∈ Ico 1 (N + 1), gap p n = p (N + 1) - 2 := by
  rw [sum_gap_telescope, hp1]

end PT.Conservation
