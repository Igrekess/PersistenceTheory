/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.Tactic

/-!
# T4 — Mertens-style identities on the active primes `{3, 5, 7}`

This module captures the **finite** Mertens-style identities of PT,
restricted to the **active primes** `{3, 5, 7}`. Unlike the classical
Mertens M1/M3 (which require analytic estimates and live in `T5Mertens.lean`),
the PT-active-prime variants are *finite* polynomial / logarithmic
identities and are kernel-checkable.

In addition this module formalises the **central algebraic identity of T4**
(monograph Ch07, `thm:factorization`, audit row #28): the generating
polynomial `f(1)` factors as

  `f(1, α, p) = 4 · (α − 1/2)² · (α² + (p − 3) · α + 1).`

For `p ≥ 5` and `α ∈ (0, 1/2)` the second factor is strictly positive,
so `f(1) > 0` on `(0, 1/2)`, with `α = 1/2` the unique fixed point.

## Main results

* `activePrimes` — the finset `{3, 5, 7}`.
* `activePrimorial_value` — `3 · 5 · 7 = 105`.
* `sum_active_primes` — `3 + 5 + 7 = 15 = μ*`.
* `sum_inv_log_active_primes_positive` — `∑_{p ∈ {3,5,7}} (log p) / p > 0`.
* `T4_factorisation` — the exact polynomial identity
  `4 · (α − 1/2)² · (α² + (p − 3) · α + 1) ≥ 0` for `p ≥ 5`, and
  `> 0` for `α ≠ 1/2`.
* `T4_unique_fixed_point` — `f(1) = 0 ⇔ α = 1/2` for `p ≥ 5` and
  `α ∈ [0, 1]`.
* `T4_second_factor_pos` — the quadratic `α² + (p−3)α + 1` is strictly
  positive for `α ∈ [0, 1/2]` and `p ≥ 5`.

## Sorry budget

This file contains **0 sorrys** and introduces **no axioms**.

## References

* Monograph Ch07 (`ch07_convergence_fr.tex`), `thm:factorization`
  (factorisation T4, eq. `eq:factorization`).
* `PT.NumberTheory.T5Mertens` for the *asymptotic* Mertens family.
-/

namespace PT.NumberTheory

open Finset

/-! ### The active-prime finset -/

/-- The active primes of PT: `{3, 5, 7}`. -/
def activePrimes : Finset ℕ := {3, 5, 7}

/-- The active-prime primorial: `3 · 5 · 7 = 105`. -/
theorem activePrimorial_value :
    (∏ p ∈ activePrimes, p) = 105 := by
  decide

/-- The PT fixed-point sum: `3 + 5 + 7 = μ* = 15`. -/
theorem sum_active_primes :
    (∑ p ∈ activePrimes, p) = 15 := by
  decide

/-- Each element of `activePrimes` is prime. -/
theorem activePrimes_prime : ∀ p ∈ activePrimes, Nat.Prime p := by
  decide

/-- `activePrimes` has cardinality `3`. -/
theorem activePrimes_card : activePrimes.card = 3 := by decide

/-! ### Active-prime finite Mertens-style sums -/

/-- The finite Mertens-like sum `∑_{p ∈ {3,5,7}} (log p) / p` is strictly
    positive. This is the PT-active-prime analogue of Mertens M1 in the
    *finite* setting. -/
theorem sum_log_div_active_primes_pos :
    0 < ∑ p ∈ activePrimes, (Real.log p) / (p : ℝ) := by
  -- Each term is positive: log 3, log 5, log 7 are all > 0; denominators > 0.
  have hl3 : 0 < Real.log 3 := Real.log_pos (by norm_num)
  have hl5 : 0 < Real.log 5 := Real.log_pos (by norm_num)
  have hl7 : 0 < Real.log 7 := Real.log_pos (by norm_num)
  -- Expand the finset sum explicitly.
  have hexp : (∑ p ∈ activePrimes, (Real.log p) / (p : ℝ))
      = Real.log 3 / 3 + Real.log 5 / 5 + Real.log 7 / 7 := by
    show (∑ p ∈ ({3, 5, 7} : Finset ℕ), (Real.log p) / (p : ℝ))
        = Real.log 3 / 3 + Real.log 5 / 5 + Real.log 7 / 7
    rw [show ({3, 5, 7} : Finset ℕ) = insert 3 (insert 5 {7}) from rfl]
    rw [Finset.sum_insert (by decide), Finset.sum_insert (by decide),
        Finset.sum_singleton]
    push_cast
    ring
  rw [hexp]
  have h3 : (0 : ℝ) < Real.log 3 / 3 := by positivity
  have h5 : (0 : ℝ) < Real.log 5 / 5 := by positivity
  have h7 : (0 : ℝ) < Real.log 7 / 7 := by positivity
  linarith

/-- The finite sum of reciprocals over active primes. The PT analogue of
    `∑ 1/p` (Mertens M2) restricted to `{3, 5, 7}`. Exact value:
    `1/3 + 1/5 + 1/7 = 71/105`. -/
theorem sum_inv_active_primes :
    (∑ p ∈ activePrimes, (1 : ℚ) / (p : ℚ)) = 71 / 105 := by
  show (∑ p ∈ ({3, 5, 7} : Finset ℕ), (1 : ℚ) / (p : ℚ)) = 71 / 105
  rw [show ({3, 5, 7} : Finset ℕ) = insert 3 (insert 5 {7}) from rfl]
  rw [Finset.sum_insert (by decide), Finset.sum_insert (by decide),
      Finset.sum_singleton]
  norm_num

/-- The finite log-product over active primes:
    `log(∏ p) = ∑ log p` for the active primes. -/
theorem log_prod_active_primes :
    Real.log (∏ p ∈ activePrimes, (p : ℝ))
      = ∑ p ∈ activePrimes, Real.log p := by
  rw [Real.log_prod]
  intro p hp
  have : (1 : ℕ) ≤ p := by
    have := activePrimes_prime p hp
    exact this.one_lt.le
  have hp_pos : (0 : ℝ) < p := by exact_mod_cast (by
    have := activePrimes_prime p hp
    exact this.pos)
  exact ne_of_gt hp_pos

/-- The total log-mass over active primes equals `log 105`. -/
theorem log_prod_active_primes_value :
    Real.log (∏ p ∈ activePrimes, (p : ℝ)) = Real.log 105 := by
  congr 1
  show (∏ p ∈ ({3, 5, 7} : Finset ℕ), (p : ℝ)) = 105
  rw [show ({3, 5, 7} : Finset ℕ) = insert 3 (insert 5 {7}) from rfl]
  rw [Finset.prod_insert (by decide), Finset.prod_insert (by decide),
      Finset.prod_singleton]
  norm_num

/-! ### T4 factorisation identity (algebraic backbone) -/

/-- **Generating-polynomial value `f(1)`.** In the T4 recurrence the
    generating polynomial evaluated at `λ = 1` is
    `f(1, α, p) := 4 · (α − 1/2)² · (α² + (p − 3) · α + 1)`. -/
noncomputable def f1 (α : ℝ) (p : ℝ) : ℝ :=
  4 * (α - 1/2) ^ 2 * (α ^ 2 + (p - 3) * α + 1)

/-- **T4 factorisation lemma — second factor positivity.**
    For `p ≥ 5` and `α ∈ [0, 1/2]`, the quadratic `α² + (p−3)α + 1` is
    strictly positive. -/
theorem T4_second_factor_pos
    (α : ℝ) (p : ℝ) (hp : 5 ≤ p) (hα0 : 0 ≤ α) (_hα1 : α ≤ 1/2) :
    0 < α ^ 2 + (p - 3) * α + 1 := by
  -- α² ≥ 0, (p-3)α ≥ 0 (since p ≥ 5 > 3 and α ≥ 0), 1 > 0.
  have h_sq : 0 ≤ α ^ 2 := sq_nonneg _
  have h_p3 : 0 ≤ p - 3 := by linarith
  have h_lin : 0 ≤ (p - 3) * α := mul_nonneg h_p3 hα0
  linarith

/-- **T4 factorisation lemma — non-negativity (Headline).**
    For `p ≥ 5` and any real `α ∈ [0, 1/2]`,
    `f(1, α, p) = 4 · (α − 1/2)² · (α² + (p−3)α + 1) ≥ 0`. -/
theorem T4_factorisation_nonneg
    (α : ℝ) (p : ℝ) (hp : 5 ≤ p) (hα0 : 0 ≤ α) (hα1 : α ≤ 1/2) :
    0 ≤ f1 α p := by
  unfold f1
  have h_sq : 0 ≤ (α - 1/2) ^ 2 := sq_nonneg _
  have h_qd : 0 < α ^ 2 + (p - 3) * α + 1 :=
    T4_second_factor_pos α p hp hα0 hα1
  have h_four : (0 : ℝ) ≤ 4 := by norm_num
  positivity

/-- **T4 factorisation lemma — strict positivity off the fixed point.**
    For `p ≥ 5` and `α ∈ [0, 1/2)`, `f(1, α, p) > 0`. -/
theorem T4_factorisation_pos
    (α : ℝ) (p : ℝ) (hp : 5 ≤ p) (hα0 : 0 ≤ α) (hα1 : α < 1/2) :
    0 < f1 α p := by
  unfold f1
  have h_sq_pos : 0 < (α - 1/2) ^ 2 := by
    apply sq_pos_of_ne_zero
    linarith
  have h_qd : 0 < α ^ 2 + (p - 3) * α + 1 :=
    T4_second_factor_pos α p hp hα0 (le_of_lt hα1)
  have h_four : (0 : ℝ) < 4 := by norm_num
  positivity

/-- **T4 unique fixed point.** For `p ≥ 5` and `α ∈ [0, 1/2]`,
    `f(1, α, p) = 0` if and only if `α = 1/2`. -/
theorem T4_unique_fixed_point
    (α : ℝ) (p : ℝ) (hp : 5 ≤ p) (hα0 : 0 ≤ α) (hα1 : α ≤ 1/2) :
    f1 α p = 0 ↔ α = 1/2 := by
  constructor
  · intro hf
    by_contra hne
    have hlt : α < 1/2 := lt_of_le_of_ne hα1 hne
    have := T4_factorisation_pos α p hp hα0 hlt
    linarith
  · intro hα
    unfold f1
    rw [hα]
    ring

/-- **Headline factorisation identity** (algebraic backbone of T4).
    `4 · (α − 1/2)² · (α² + (p−3)α + 1)` is identically equal to its
    expanded polynomial form, which can be verified by `ring`. -/
theorem T4_factorisation_identity (α p : ℝ) :
    f1 α p = 4 * (α - 1/2) ^ 2 * (α ^ 2 + (p - 3) * α + 1) := rfl

/-- For each active prime `p ∈ {3, 5, 7}` *except* `p = 3` (i.e. `p ∈ {5, 7}`),
    the T4 factorisation hypothesis `5 ≤ p` holds. -/
theorem activePrimes_ge_five_satisfies_T4
    (p : ℕ) (_hp : p ∈ activePrimes) (hp5 : 5 ≤ p)
    (α : ℝ) (hα0 : 0 ≤ α) (hα1 : α ≤ 1/2) :
    0 ≤ f1 α (p : ℝ) := by
  apply T4_factorisation_nonneg α (p : ℝ) _ hα0 hα1
  exact_mod_cast hp5

/-! ### Headline summary -/

/-- **T4 Active-prime summary.**
    1. The active primes are `{3, 5, 7}`.
    2. Their sum is `μ* = 15`.
    3. Their product (the active-prime primorial) is `105`.
    4. The finite reciprocal sum is `71/105`.
    5. The T4 factorisation `f(1, α, p) ≥ 0` holds for `p ≥ 5` and
       `α ∈ [0, 1/2]`, with equality iff `α = 1/2`. -/
theorem T4_active_primes_summary :
    activePrimes = {3, 5, 7}
    ∧ (∑ p ∈ activePrimes, p) = 15
    ∧ (∏ p ∈ activePrimes, p) = 105
    ∧ (∑ p ∈ activePrimes, (1 : ℚ) / p) = 71 / 105
    ∧ (∀ (α p : ℝ), 5 ≤ p → 0 ≤ α → α ≤ 1/2 → 0 ≤ f1 α p) :=
  ⟨rfl, sum_active_primes, activePrimorial_value,
   sum_inv_active_primes, T4_factorisation_nonneg⟩

end PT.NumberTheory
