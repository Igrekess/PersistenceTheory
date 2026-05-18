/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Sieve.T3Antidiagonal

/-!
# s = 1/2 — Foundational Theorem of Persistence Theory

The unique stationary distribution of the sieve-level transfer matrix `T₃`
(see `PT.Sieve.T3Antidiagonal`) is `π = (1/2, 1/2)`. The corresponding
symmetry parameter

$$s := \pi(0) = \pi(1) = \tfrac{1}{2}$$

is the **fundamental axiom-as-theorem** of Persistence Theory: every other
PT result is downstream of `s = 1/2`.

This file proves that `s = 1/2` is *forced* by `T₃`, hence ultimately forced by
the elementary combinatorics of consecutive 6-rough residues (Theorem T1).

## Reference

Section 4 of: Y. Senez, M1 article (PT_ARTICLES/PT_MATHEMATICS/M1).
Monograph Chapter 1, §1.3 (the canonical statement of `s = 1/2`).
-/

namespace PT.Stochastic

open Matrix PT.Sieve

/-- A distribution `π : Fin 2 → ℝ` is **stationary for `T₃`** if it is a
    probability vector left-fixed by `T₃`:

    * `π 0 + π 1 = 1`
    * `π i ≥ 0`
    * `πᵀ · T₃ = πᵀ`, i.e. `∑ j, π j * T₃ j i = π i`. -/
structure IsStationary (π : Fin 2 → ℝ) : Prop where
  prob_sum : π 0 + π 1 = 1
  nonneg   : ∀ i, 0 ≤ π i
  fixed    : ∀ i, ∑ j, π j * T3 j i = π i

/-- The candidate stationary distribution `π = (1/2, 1/2)`. -/
noncomputable def piHalf : Fin 2 → ℝ := ![1/2, 1/2]

@[simp] lemma piHalf_zero : piHalf 0 = 1/2 := rfl
@[simp] lemma piHalf_one  : piHalf 1 = 1/2 := rfl

/-- `(1/2, 1/2)` is a stationary distribution for `T₃`. -/
theorem piHalf_isStationary : IsStationary piHalf where
  prob_sum := by norm_num [piHalf]
  nonneg := by
    intro i; fin_cases i <;> simp [piHalf]
  fixed := by
    intro i
    fin_cases i <;>
      simp [piHalf, T3, Fin.sum_univ_two]

/-- **Uniqueness of the stationary distribution.**
    Any stationary distribution of `T₃` equals `(1/2, 1/2)`. -/
theorem T3_unique_stationary (π : Fin 2 → ℝ) (h : IsStationary π) :
    π = piHalf := by
  -- The fixed-point equation πᵀ T₃ = πᵀ on index `0`:
  --   π 0 * T3 0 0 + π 1 * T3 1 0 = π 0
  -- i.e.                  π 1 = π 0
  have h0 := h.fixed 0
  have eq_at_0 : π 1 = π 0 := by
    have := h0
    simp [T3, Fin.sum_univ_two] at this
    linarith
  -- Combined with π 0 + π 1 = 1 we get π 0 = π 1 = 1/2
  have hsum := h.prob_sum
  have hp0 : π 0 = 1/2 := by linarith
  have hp1 : π 1 = 1/2 := by linarith
  ext i
  fin_cases i <;> simp [piHalf, hp0, hp1]

/-! ### The symmetry parameter `s = 1/2` -/

/-- The PT symmetry parameter `s`. -/
noncomputable def s : ℝ := 1 / 2

@[simp] lemma s_def : s = 1 / 2 := rfl

/-- `π 0 = s`. -/
theorem s_eq_piHalf_zero : piHalf 0 = s := by simp [s]

/-- `π 1 = s`. -/
theorem s_eq_piHalf_one : piHalf 1 = s := by simp [s]

/-- **Foundational theorem of PT.**
    The value `s = 1/2` is forced by `T₃`: for any stationary distribution `π`
    of `T₃`, `π i = s` for every `i`. -/
theorem s_is_T3_stationary_value (π : Fin 2 → ℝ) (h : IsStationary π) (i : Fin 2) :
    π i = s := by
  have := T3_unique_stationary π h
  rw [this]
  fin_cases i <;> simp [s]

/-- `s = 1/2` is the unique fixed value of any `T₃`-stationary distribution. -/
theorem s_eq_one_half : s = 1 / 2 := rfl

end PT.Stochastic
