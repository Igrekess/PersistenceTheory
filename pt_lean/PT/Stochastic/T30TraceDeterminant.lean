/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Stochastic.T2T30Normalisation

/-!
# T2 — Trace and determinant invariants of `T_30 = T_2 ⊗ T_3 ⊗ T_5`

**Statement (paper-level, monograph ch03 §T2).**
The CRT factorisation `T_30 ≅ T_2 ⊗ T_3 ⊗ T_5` (cf.
`PT.Stochastic.T2T30Normalisation`) makes the trace and determinant of `T_30`
*multiplicative* across the three prime factors. Combined with the elementary
invariants of the individual factors

* `trace(T_2) = 1`, `det(T_2) = 1` (the trivial `1 × 1` Perron block);
* `trace(T_3) = 0`, `det(T_3) = -1` (the antidiagonal involution
  `!![0,1;1,0]`, cf. `PT.Sieve.T3Antidiagonal`),

this yields immediately

  `trace(T_30) = trace(T_2) · trace(T_3) · trace(T_5) = 1 · 0 · trace(T_5) = 0`,

so **`trace(T_30) = 0` unconditionally**, regardless of the underlying `T_5`.
The determinant identity is the standard Kronecker formula
`det(A ⊗ B) = det(A)^m · det(B)^n` for `A` of size `n` and `B` of size `m`,
specialised to the three-factor case.

## Main results

### Invariants of the individual factors
* `T2_trivial_trace` — `trace(T_2) = 1`.
* `T2_trivial_det` — `det(T_2) = 1`.
* (`T3_trace_zero`, `T3_det_neg_one` are imported from `PT.Sieve.T3Antidiagonal`.)

### Multiplicativity of the trace under Kronecker product
* `kron3_trace` — `trace((A ⊗ B) ⊗ C) = trace(A) · trace(B) · trace(C)`.

### Multiplicativity of the determinant under Kronecker product
* `kron3_det` — closed form for `det((A ⊗ B) ⊗ C)` in terms of the three
  factor determinants and the cardinalities of the index sets.

### Headline specialisations to `T_30`
* `T30_trace_eq` — `trace(T_30) = trace(T_2) · trace(T_3) · trace(T_5)`.
* `T30_trace_zero` — **`trace(T_30) = 0`** (the headline, holds for *any*
  `T5Like` because the factor `trace(T_3) = 0` annihilates the product).
* `T30_det_eq` — closed form for `det(T_30)` in terms of `det(T_5)`,
  derived from the Kronecker determinant formula.

## Strategy

All proofs reduce to two applications of the two-factor Mathlib lemmas:

* `Matrix.trace_kronecker : trace(A ⊗ₖ B) = trace(A) · trace(B)`;
* `Matrix.det_kronecker   : det(A ⊗ₖ B) = det(A)^(card n) · det(B)^(card m)`.

Applied to the left-bracketed `kron3 A B C = (A ⊗ₖ B) ⊗ₖ C`, each lemma
delivers the three-factor identity in one rewrite step. The `T_30`
specialisations then plug in `trace(T_2) = 1`, `trace(T_3) = 0`,
`det(T_2) = 1`, `det(T_3) = -1` and let `simp`/`ring` close.

## Status

`[THM]` — these are elementary linear-algebraic invariants of the CRT
factorisation. The trace identity `trace(T_30) = 0` is *unconditional* in
`T_5` (no spectral data of `T_5` is needed), making it the cleanest sanity
check of the three-factor Kronecker structure.

## References

* Monograph Chapter 3, §"Spectral conservation T2", invariants table.
* `PT.Sieve.T3Antidiagonal` — `T3_trace_zero`, `T3_det_neg_one`.
* `PT.Stochastic.T2T30Normalisation` — `T2_trivial`, `T30`, `kron3`.
* Mathlib `Matrix.trace_kronecker`, `Matrix.det_kronecker`.
-/

namespace PT.Stochastic

open Matrix Kronecker BigOperators PT.Sieve

/-! ### Trace and determinant of the trivial `T_2` block -/

/-- **`trace(T_2) = 1`.** Since `T_2 = 1 : Matrix (Fin 1) (Fin 1) ℝ` is the
    `1 × 1` identity, its trace is the single diagonal entry `1`. -/
@[simp] theorem T2_trivial_trace : T2_trivial.trace = 1 := by
  simp [T2_trivial, Matrix.trace]

/-- **`det(T_2) = 1`.** The determinant of the `1 × 1` identity matrix is `1`. -/
@[simp] theorem T2_trivial_det : T2_trivial.det = 1 := by
  simp [T2_trivial]

/-! ### Trace multiplicativity under three-factor Kronecker product -/

/-- **Trace of a three-factor Kronecker product.**
    `trace((A ⊗ B) ⊗ C) = trace(A) · trace(B) · trace(C)`.

    Direct consequence of `Matrix.trace_kronecker` applied twice
    (once for the inner `A ⊗ B`, once for the outer `· ⊗ C`). -/
theorem kron3_trace
    {m n p : Type*} [Fintype m] [Fintype n] [Fintype p]
    (A : Matrix m m ℝ) (B : Matrix n n ℝ) (C : Matrix p p ℝ) :
    (kron3 A B C).trace = A.trace * B.trace * C.trace := by
  unfold kron3
  rw [Matrix.trace_kronecker, Matrix.trace_kronecker]

/-! ### Determinant multiplicativity under three-factor Kronecker product -/

/-- **Determinant of a three-factor Kronecker product (closed form).**
    For square matrices `A : Matrix m m`, `B : Matrix n n`, `C : Matrix p p`,
    the left-bracketed Kronecker product `(A ⊗ B) ⊗ C` has determinant

    `det((A ⊗ B) ⊗ C)
      = (det A ^ card n) ^ card p · (det B ^ card m) ^ card p · det C ^ (card m · card n)`.

    This follows by two applications of `Matrix.det_kronecker`: the inner
    product gives `det(A ⊗ B) = det A ^ card n · det B ^ card m`, and the
    outer product raises this to the `card p`-th power (and multiplies by
    `det C ^ (card (m × n))`). -/
theorem kron3_det
    {m n p : Type*} [Fintype m] [Fintype n] [Fintype p]
    [DecidableEq m] [DecidableEq n] [DecidableEq p]
    (A : Matrix m m ℝ) (B : Matrix n n ℝ) (C : Matrix p p ℝ) :
    (kron3 A B C).det
      = (A.det ^ Fintype.card n * B.det ^ Fintype.card m) ^ Fintype.card p
        * C.det ^ Fintype.card (m × n) := by
  unfold kron3
  rw [Matrix.det_kronecker, Matrix.det_kronecker, mul_pow]

/-! ### Headline specialisations to `T_30 = T_2 ⊗ T_3 ⊗ T_5` -/

/-- **Trace of `T_30` factorises across the three factors.** -/
theorem T30_trace_eq (T5 : T5Like) :
    (T30 T5).trace = T2_trivial.trace * T3.trace * T5.matrix.trace := by
  unfold T30
  exact kron3_trace T2_trivial T3 T5.matrix

/-- **Headline: `trace(T_30) = 0`.**

    Because the middle factor `T_3 = !![0,1;1,0]` has vanishing trace
    (`PT.Sieve.T3Antidiagonal.T3_trace_zero`), the multiplicative
    factorisation `trace(T_30) = trace(T_2) · trace(T_3) · trace(T_5)`
    forces `trace(T_30) = 0`, **unconditionally** in the choice of `T_5`.
    This is the cleanest sanity check of the CRT factorisation `T_30 ≅ T_2 ⊗
    T_3 ⊗ T_5`: a single zero factor annihilates the product. -/
@[simp] theorem T30_trace_zero (T5 : T5Like) : (T30 T5).trace = 0 := by
  rw [T30_trace_eq, T3_trace_zero]
  ring

/-- **Determinant of `T_30` (closed form in `det(T_5)`).**

    Plug `det(T_2) = 1` and `det(T_3) = -1` into the three-factor Kronecker
    determinant formula `kron3_det`. The index cardinalities are
    `card (Fin 1) = 1`, `card (Fin 2) = 2`, `card (Fin 4) = 4`, so

    `det(T_30) = (1^2 · (-1)^1)^4 · det(T_5)^(1·2) = 1 · det(T_5)^2`. -/
theorem T30_det_eq (T5 : T5Like) :
    (T30 T5).det = T5.matrix.det ^ 2 := by
  unfold T30
  rw [kron3_det T2_trivial T3 T5.matrix, T2_trivial_det, T3_det_neg_one]
  norm_num

/-! ### Two-factor specialisations (intermediate sanity checks) -/

/-- **Trace of `T_2 ⊗ T_3` is zero.** Immediate from `trace_kronecker` and
    `trace(T_3) = 0`. -/
@[simp] theorem T2_kron_T3_trace_zero :
    (T2_trivial ⊗ₖ T3).trace = 0 := by
  rw [Matrix.trace_kronecker, T3_trace_zero]
  ring

/-- **Determinant of `T_2 ⊗ T_3` is `-1`.** Computed via `det_kronecker`,
    with `det(T_2) = 1`, `det(T_3) = -1`, and the cardinalities
    `card (Fin 1) = 1`, `card (Fin 2) = 2`. -/
@[simp] theorem T2_kron_T3_det :
    (T2_trivial ⊗ₖ T3).det = -1 := by
  rw [Matrix.det_kronecker, T2_trivial_det, T3_det_neg_one]
  simp

/-! ### Summary -/

/-- **Headline summary for the invariants of `T_30`.**

    The CRT factorisation `T_30 = T_2 ⊗ T_3 ⊗ T_5` (left-bracketed Kronecker
    product) yields:

    1. `trace(T_30) = trace(T_2) · trace(T_3) · trace(T_5)` (multiplicativity
       of the trace under the three-factor Kronecker product).
    2. `trace(T_30) = 0` (the middle factor `trace(T_3) = 0` annihilates
       the product; **unconditional in `T_5`**).
    3. `det(T_30) = det(T_5)^2` (the determinant factorisation, plugging
       `det(T_2) = 1`, `det(T_3) = -1`, and the index cardinalities
       `1, 2, 4`).
-/
theorem T30_invariants_summary (T5 : T5Like) :
    (T30 T5).trace = T2_trivial.trace * T3.trace * T5.matrix.trace
    ∧ (T30 T5).trace = 0
    ∧ (T30 T5).det = T5.matrix.det ^ 2 :=
  ⟨T30_trace_eq T5, T30_trace_zero T5, T30_det_eq T5⟩

end PT.Stochastic
