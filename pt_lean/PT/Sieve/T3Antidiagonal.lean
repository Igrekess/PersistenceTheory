/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import Mathlib.Tactic
import Mathlib.LinearAlgebra.Matrix.Notation
import Mathlib.LinearAlgebra.Matrix.Trace
import Mathlib.LinearAlgebra.Matrix.Determinant.Basic
import Mathlib.Data.Matrix.Mul
import Mathlib.Data.Real.Basic

/-!
# T3 — The Sieve-Level Transfer Matrix

The transfer matrix of consecutive 6-rough residues modulo 3, with index
convention `0 ↦ 1 mod 3` and `1 ↦ 2 mod 3`, is the anti-diagonal permutation
matrix

$$T_3 = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}.$$

It is doubly stochastic, has trace `0`, determinant `−1`, eigenvalues `{+1, −1}`,
and is the unique trace-`0` doubly-stochastic 2×2 real matrix.

This file formalises these basic linear-algebraic facts. The stationary
distribution `π = (1/2, 1/2)` is treated in `PT.Stochastic.SHalf`.

## Reference

Theorem T1 → T3 link, Section 3 of:
Y. Senez, M1 article (PT_ARTICLES/PT_MATHEMATICS/M1).
-/

namespace PT.Sieve

open Matrix

/-- The sieve-level transfer matrix `T₃` on the two non-zero residues modulo 3,
    indexed by `0 ↦ 1 mod 3`, `1 ↦ 2 mod 3`. -/
def T3 : Matrix (Fin 2) (Fin 2) ℝ := !![0, 1; 1, 0]

/-! ### Pointwise entries -/

@[simp] lemma T3_zero_zero : T3 0 0 = 0 := rfl
@[simp] lemma T3_zero_one  : T3 0 1 = 1 := rfl
@[simp] lemma T3_one_zero  : T3 1 0 = 1 := rfl
@[simp] lemma T3_one_one   : T3 1 1 = 0 := rfl

/-! ### Antidiagonal structure -/

/-- `T₃` is the anti-diagonal permutation matrix: diagonal entries vanish. -/
theorem T3_is_antidiagonal :
    T3 0 0 = 0 ∧ T3 1 1 = 0 ∧ T3 0 1 = 1 ∧ T3 1 0 = 1 := by
  refine ⟨rfl, rfl, rfl, rfl⟩

/-! ### Doubly stochastic -/

/-- Each row of `T₃` sums to 1. -/
theorem T3_row_sum (i : Fin 2) : ∑ j, T3 i j = 1 := by
  fin_cases i <;> simp [T3, Fin.sum_univ_two]

/-- Each column of `T₃` sums to 1. -/
theorem T3_col_sum (j : Fin 2) : ∑ i, T3 i j = 1 := by
  fin_cases j <;> simp [T3, Fin.sum_univ_two]

/-- All entries of `T₃` are non-negative. -/
theorem T3_nonneg (i j : Fin 2) : 0 ≤ T3 i j := by
  fin_cases i <;> fin_cases j <;> simp [T3]

/-- `T₃` is doubly stochastic: non-negative entries, every row and every column
    summing to 1. -/
theorem T3_doubly_stochastic :
    (∀ i j, 0 ≤ T3 i j) ∧ (∀ i, ∑ j, T3 i j = 1) ∧ (∀ j, ∑ i, T3 i j = 1) :=
  ⟨T3_nonneg, T3_row_sum, T3_col_sum⟩

/-! ### Trace and determinant -/

/-- The trace of `T₃` is zero. -/
theorem T3_trace_zero : T3.trace = 0 := by
  simp [Matrix.trace, T3, Fin.sum_univ_two]

/-- The determinant of `T₃` is `−1`. -/
theorem T3_det_neg_one : T3.det = -1 := by
  simp [T3, Matrix.det_fin_two]

/-! ### Eigenvalues `{+1, -1}` -/

/-- `(1, 1)` is an eigenvector of `T₃` with eigenvalue `+1`. -/
theorem T3_eigen_plus_one :
    T3.mulVec ![1, 1] = (1 : ℝ) • ![1, 1] := by
  ext i
  fin_cases i <;>
    simp [T3, Matrix.mulVec, dotProduct]

/-- `(1, -1)` is an eigenvector of `T₃` with eigenvalue `−1`. -/
theorem T3_eigen_minus_one :
    T3.mulVec ![1, -1] = (-1 : ℝ) • ![1, -1] := by
  ext i
  fin_cases i <;>
    simp [T3, Matrix.mulVec, dotProduct]

/-- The eigenvalues of `T₃` are exactly `{+1, −1}`: the matrix satisfies
    `T₃ ^ 2 = I`, so its minimal polynomial divides `X^2 - 1`. Combined with
    trace `0` and determinant `−1`, the characteristic polynomial is `X^2 - 1`.
    Numerically: -/
theorem T3_sq_eq_one : T3 * T3 = (1 : Matrix (Fin 2) (Fin 2) ℝ) := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [T3, Matrix.mul_apply, Fin.sum_univ_two]

/-- The characteristic polynomial of `T₃` (encoded as the pair
    `(trace, det) = (0, -1)`) gives `t^2 - 1`, hence eigenvalues `±1`. -/
theorem T3_char_poly_data : T3.trace = 0 ∧ T3.det = -1 :=
  ⟨T3_trace_zero, T3_det_neg_one⟩

/-! ### Uniqueness among trace-0 doubly stochastic 2×2 matrices -/

/-- Among 2×2 doubly-stochastic real matrices with trace zero, `T₃` is the
    unique one.

    A 2×2 matrix `M` is doubly stochastic with trace `0` iff its entries satisfy
    `M 0 0 + M 0 1 = 1`, `M 1 0 + M 1 1 = 1`, `M 0 0 = M 1 1 = 0`, hence
    `M 0 1 = M 1 0 = 1`. -/
theorem T3_unique_trace_zero_doubly_stochastic
    (M : Matrix (Fin 2) (Fin 2) ℝ)
    (hrows : ∀ i, ∑ j, M i j = 1)
    (_hcols : ∀ j, ∑ i, M i j = 1)
    (htr   : M.trace = 0)
    (hnn   : ∀ i j, 0 ≤ M i j) :
    M = T3 := by
  -- Extract entry-level constraints (Fin 2 sums expand to two terms)
  have hr0 : M 0 0 + M 0 1 = 1 := by
    have h := hrows 0
    simp [Fin.sum_univ_two] at h
    linarith
  have hr1 : M 1 0 + M 1 1 = 1 := by
    have h := hrows 1
    simp [Fin.sum_univ_two] at h
    linarith
  have htrace : M 0 0 + M 1 1 = 0 := by
    have h := htr
    simp [Matrix.trace, Fin.sum_univ_two] at h
    linarith
  -- Non-negativity forces diagonal entries to zero
  have hM00 : M 0 0 = 0 := le_antisymm
    (by linarith [hnn 1 1]) (hnn 0 0)
  have hM11 : M 1 1 = 0 := le_antisymm
    (by linarith [hnn 0 0]) (hnn 1 1)
  -- Hence off-diagonals are 1
  have hM01 : M 0 1 = 1 := by linarith
  have hM10 : M 1 0 = 1 := by linarith
  -- Conclude entry-wise
  ext i j
  fin_cases i <;> fin_cases j <;> simp [T3, hM00, hM01, hM10, hM11]

end PT.Sieve
