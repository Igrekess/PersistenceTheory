/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Stochastic.T30TraceDeterminant
import PT.Stochastic.SpectralDominance

/-!
# T2 — Trace formula for `T_30^n` via Kronecker-power factorisation

**Statement (paper-level, monograph ch03 §T2, follow-up to
`T30TraceDeterminant`).** The headline `trace(T_30) = 0` of
`PT.Stochastic.T30TraceDeterminant` extends to every natural power `n`
of `T_30`. The key tool is the *Kronecker power* identity

  `(A ⊗ₖ B)^n  =  A^n ⊗ₖ B^n`     (for square matrices, `CommSemiring`)

obtained by induction from Mathlib's `Matrix.mul_kronecker_mul`. Combined
with the three-factor structure `T_30 = T_2 ⊗ T_3 ⊗ T_5`, this gives

  `(T_30)^n  =  T_2^n ⊗ T_3^n ⊗ T_5^n`

and, by `kron3_trace`,

  `trace(T_30^n)  =  trace(T_2^n) · trace(T_3^n) · trace(T_5^n)`.

Specialising:

* `trace(T_2^n) = 1` (since `T_2 = 1 : Matrix (Fin 1) (Fin 1) ℝ`).
* `trace(T_3^{2k}) = 2`, `trace(T_3^{2k+1}) = 0` (from `T3_pow_even`,
  `T3_pow_odd` and the explicit trace of the `2 × 2` identity / antidiagonal).
* Hence `trace(T_30^n) = 2 · trace(T_5^n)` for `n` even, and
  `trace(T_30^n) = 0` for `n` odd.

The odd case is the **unconditional** continuation of `T30_trace_zero`
to every odd power.

## Main results

### Kronecker power factorisation
* `kronecker_pow` — `(A ⊗ₖ B)^n = A^n ⊗ₖ B^n` (square matrices, induction).
* `kron3_pow` — `(kron3 A B C)^n = kron3 (A^n) (B^n) (C^n)`.

### Trace of `(T_30)^n` reduces to factor traces
* `T30_pow_eq_kron3_pow` — `(T_30)^n = kron3 (T_2^n) (T_3^n) (T_5^n)`.
* `T30_pow_trace_eq` — `trace((T_30)^n) =
  trace(T_2^n) · trace(T_3^n) · trace(T_5^n)`.

### Individual factor traces
* `T2_trivial_pow_trace` — `trace(T_2^n) = 1` for every `n`.
* `T3_pow_even_trace` — `trace(T_3^{2k}) = 2`.
* `T3_pow_odd_trace`  — `trace(T_3^{2k+1}) = 0`.

### Headline parity formula
* `T30_pow_trace_odd` — `trace((T_30)^{2k+1}) = 0` (every odd power
  has vanishing trace, generalising `T30_trace_zero`).
* `T30_pow_trace_even` — `trace((T_30)^{2k}) = 2 · trace(T_5^{2k})`
  (every even power factors `T_5`).
* `T30_pow_trace_dichotomy` — packaged even/odd dichotomy in one
  statement.

### Explicit small cases
* `T30_pow_one_trace`  — `trace(T_30^1) = 0` (recovers `T30_trace_zero`).
* `T30_pow_two_trace`  — `trace(T_30^2) = 2 · trace(T_5^2)`.
* `T30_pow_three_trace` — `trace(T_30^3) = 0`.

## Strategy

The Kronecker power identity `kronecker_pow` proceeds by induction on `n`,
using `Matrix.mul_kronecker_mul` at the inductive step (the case `n = 0`
is `kroneckerMap_one_one`-style, captured directly via `pow_zero`). The
three-factor version `kron3_pow` is then two applications of
`kronecker_pow`. The trace formula reduces to `kron3_trace` applied to
`(T_2^n, T_3^n, T_5^n)`.

The factor traces are elementary:

* `T_2 = (1 : Matrix (Fin 1) (Fin 1) ℝ)`, so `T_2^n = 1` and the trace
  of the `1 × 1` identity is `1`.
* `T_3^(2k) = I` has trace `2`; `T_3^(2k+1) = T_3` has trace `0`
  (re-using `T3_pow_even`, `T3_pow_odd`, `T3_trace_zero`).

The general parity formula follows from `Nat.even_or_odd'`.

## Status

`[THM]` — purely linear-algebraic; no PT numerical input beyond what is
already encoded in `T_3` (an explicit `2 × 2` matrix) and in `T_2` (the
trivial scalar). The structural three-factor reduction is kernel-verified
by Mathlib's Kronecker product machinery (`mul_kronecker_mul`,
`trace_kronecker`).

## References

* `PT.Stochastic.T30TraceDeterminant` — `T30_trace_zero` (the `n = 1`
  case), `kron3_trace`.
* `PT.Stochastic.SpectralDominance` — `T3_pow_even`, `T3_pow_odd`.
* `PT.Sieve.T3Antidiagonal` — `T3_trace_zero`.
* Mathlib `Matrix.mul_kronecker_mul`, `Matrix.trace_kronecker`.
-/

namespace PT.Stochastic

open Matrix Kronecker BigOperators PT.Sieve

/-! ### Kronecker power factorisation `(A ⊗ₖ B)^n = A^n ⊗ₖ B^n` -/

/-- **Kronecker product of identity matrices is the identity.**
    Auxiliary base case for `kronecker_pow`. -/
private lemma one_kronecker_one
    {m n : Type*} [Fintype m] [Fintype n] [DecidableEq m] [DecidableEq n] :
    (1 : Matrix m m ℝ) ⊗ₖ (1 : Matrix n n ℝ) = (1 : Matrix (m × n) (m × n) ℝ) := by
  ext ⟨i₁, i₂⟩ ⟨j₁, j₂⟩
  by_cases h₁ : i₁ = j₁
  · by_cases h₂ : i₂ = j₂
    · subst h₁; subst h₂
      simp [Matrix.kroneckerMap, Matrix.one_apply]
    · simp [Matrix.kroneckerMap, Matrix.one_apply, h₂]
  · simp [Matrix.kroneckerMap, Matrix.one_apply, h₁]

/-- **Kronecker power identity.** For square matrices `A : Matrix m m ℝ`,
    `B : Matrix n n ℝ` and every natural number `k`,

    `(A ⊗ₖ B)^k  =  A^k ⊗ₖ B^k`.

    Proof by induction on `k`, base case `kronecker(1, 1) = 1`,
    inductive step `Matrix.mul_kronecker_mul`. -/
theorem kronecker_pow
    {m n : Type*} [Fintype m] [Fintype n] [DecidableEq m] [DecidableEq n]
    (A : Matrix m m ℝ) (B : Matrix n n ℝ) (k : ℕ) :
    (A ⊗ₖ B) ^ k = (A ^ k) ⊗ₖ (B ^ k) := by
  induction k with
  | zero =>
    rw [pow_zero, pow_zero, pow_zero]
    exact (one_kronecker_one).symm
  | succ k ih =>
    rw [pow_succ, ih, pow_succ, pow_succ, ← Matrix.mul_kronecker_mul]

/-- **Three-factor Kronecker power identity.**
    `(kron3 A B C)^k = kron3 (A^k) (B^k) (C^k)`. Two applications of
    `kronecker_pow`. -/
theorem kron3_pow
    {m n p : Type*} [Fintype m] [Fintype n] [Fintype p]
    [DecidableEq m] [DecidableEq n] [DecidableEq p]
    (A : Matrix m m ℝ) (B : Matrix n n ℝ) (C : Matrix p p ℝ) (k : ℕ) :
    (kron3 A B C) ^ k = kron3 (A ^ k) (B ^ k) (C ^ k) := by
  unfold kron3
  rw [kronecker_pow, kronecker_pow]

/-! ### `(T_30)^n` factorises as `(T_2^n) ⊗ (T_3^n) ⊗ (T_5^n)` -/

/-- **Kronecker power of `T_30`.** Direct corollary of `kron3_pow`:
    `(T_30)^n = T_2^n ⊗ T_3^n ⊗ T_5^n` (left-bracketed). -/
theorem T30_pow_eq_kron3_pow (T5 : T5Like) (n : ℕ) :
    (T30 T5) ^ n = kron3 (T2_trivial ^ n) (T3 ^ n) (T5.matrix ^ n) := by
  unfold T30
  exact kron3_pow T2_trivial T3 T5.matrix n

/-- **Trace of `(T_30)^n` factors across the three factor traces.** -/
theorem T30_pow_trace_eq (T5 : T5Like) (n : ℕ) :
    ((T30 T5) ^ n).trace
      = (T2_trivial ^ n).trace * (T3 ^ n).trace * (T5.matrix ^ n).trace := by
  rw [T30_pow_eq_kron3_pow]
  exact kron3_trace _ _ _

/-! ### Trace of `T_2^n` -/

/-- **`T_2^n = T_2` for every `n ≥ 1`, and `T_2^0 = 1 = T_2`.**
    Since `T_2 = (1 : Matrix (Fin 1) (Fin 1) ℝ)`, every power equals `1`. -/
@[simp] theorem T2_trivial_pow (n : ℕ) : T2_trivial ^ n = T2_trivial := by
  unfold T2_trivial
  exact one_pow n

/-- **Trace of `T_2^n` is `1`.** -/
@[simp] theorem T2_trivial_pow_trace (n : ℕ) : (T2_trivial ^ n).trace = 1 := by
  rw [T2_trivial_pow]
  exact T2_trivial_trace

/-! ### Trace of `T_3^n` (even / odd parity)

The orbit `{T_3^n : n ∈ ℕ}` alternates `{I, T_3, I, T_3, …}`. The trace
of the `2 × 2` identity is `2`; the trace of `T_3` is `0`. Hence the
trace of `T_3^n` is `2` on even `n` and `0` on odd `n`. -/

/-- **Trace of `(I : Matrix (Fin 2) (Fin 2) ℝ)` is `2`.** -/
@[simp] theorem T3_one_trace : (1 : Matrix (Fin 2) (Fin 2) ℝ).trace = 2 := by
  simp [Matrix.trace]

/-- **Trace of `T_3^(2k)` is `2`.** Using `T3_pow_even`: `T_3^(2k) = I`,
    whose trace is `2`. -/
theorem T3_pow_even_trace (k : ℕ) : (T3 ^ (2 * k)).trace = 2 := by
  rw [T3_pow_even]
  exact T3_one_trace

/-- **Trace of `T_3^(2k+1)` is `0`.** Using `T3_pow_odd`:
    `T_3^(2k+1) = T_3`, whose trace is `0`. -/
theorem T3_pow_odd_trace (k : ℕ) : (T3 ^ (2 * k + 1)).trace = 0 := by
  rw [T3_pow_odd]
  exact T3_trace_zero

/-! ### Headline parity formula for `trace((T_30)^n)` -/

/-- **Headline (odd case).** For every odd power `n = 2k + 1`,

    `trace((T_30)^{2k+1})  =  0`.

    This is the unconditional extension of `T30_trace_zero` to every odd
    `n`: the middle factor `trace(T_3^{2k+1}) = 0` annihilates the
    three-factor product, regardless of `T_5`. -/
theorem T30_pow_trace_odd (T5 : T5Like) (k : ℕ) :
    ((T30 T5) ^ (2 * k + 1)).trace = 0 := by
  rw [T30_pow_trace_eq, T3_pow_odd_trace]
  ring

/-- **Headline (even case).** For every even power `n = 2k`,

    `trace((T_30)^{2k})  =  2 · trace(T_5^{2k})`.

    The factors contribute `trace(T_2^{2k}) = 1`, `trace(T_3^{2k}) = 2`,
    and `trace(T_5^{2k})`. -/
theorem T30_pow_trace_even (T5 : T5Like) (k : ℕ) :
    ((T30 T5) ^ (2 * k)).trace = 2 * (T5.matrix ^ (2 * k)).trace := by
  rw [T30_pow_trace_eq, T2_trivial_pow_trace, T3_pow_even_trace]
  ring

/-- **Dichotomy.** For every `n`, exactly one of the two parity formulas
    holds:

    * if `n = 2k`, then `trace((T_30)^n) = 2 · trace(T_5^n)`;
    * if `n = 2k+1`, then `trace((T_30)^n) = 0`.

    A packaged restatement of `T30_pow_trace_even` and
    `T30_pow_trace_odd`, via `Nat.even_or_odd'`. -/
theorem T30_pow_trace_dichotomy (T5 : T5Like) (n : ℕ) :
    (∃ k, n = 2 * k ∧ ((T30 T5) ^ n).trace = 2 * (T5.matrix ^ n).trace)
    ∨ (∃ k, n = 2 * k + 1 ∧ ((T30 T5) ^ n).trace = 0) := by
  obtain ⟨k, hk | hk⟩ := Nat.even_or_odd' n
  · left
    refine ⟨k, hk, ?_⟩
    rw [hk]
    exact T30_pow_trace_even T5 k
  · right
    refine ⟨k, hk, ?_⟩
    rw [hk]
    exact T30_pow_trace_odd T5 k

/-! ### Explicit small cases -/

/-- **`trace((T_30)^1) = 0`.** Recovers `T30_trace_zero` from the
    parity formula at `k = 0`. -/
theorem T30_pow_one_trace (T5 : T5Like) : ((T30 T5) ^ 1).trace = 0 := by
  have h := T30_pow_trace_odd T5 0
  have h1 : 2 * 0 + 1 = 1 := by norm_num
  rw [h1] at h
  exact h

/-- **`trace((T_30)^2) = 2 · trace(T_5^2)`.** Specialisation of
    `T30_pow_trace_even` at `k = 1`. -/
theorem T30_pow_two_trace (T5 : T5Like) :
    ((T30 T5) ^ 2).trace = 2 * (T5.matrix ^ 2).trace := by
  have h := T30_pow_trace_even T5 1
  simpa using h

/-- **`trace((T_30)^3) = 0`.** Specialisation of `T30_pow_trace_odd`
    at `k = 1`. -/
theorem T30_pow_three_trace (T5 : T5Like) : ((T30 T5) ^ 3).trace = 0 := by
  have h := T30_pow_trace_odd T5 1
  simpa using h

/-- **`trace((T_30)^4) = 2 · trace(T_5^4)`.** Specialisation at `k = 2`. -/
theorem T30_pow_four_trace (T5 : T5Like) :
    ((T30 T5) ^ 4).trace = 2 * (T5.matrix ^ 4).trace := by
  have h := T30_pow_trace_even T5 2
  simpa using h

/-! ### Summary -/

/-- **Headline summary for `trace((T_30)^n)`.**

    The CRT factorisation `T_30 = T_2 ⊗ T_3 ⊗ T_5` propagates through
    matrix powers via `kron3_pow`. The trace formula collapses two of
    the three factors:

    1. `T_2^n` always has trace `1` (trivial `1 × 1` block);
    2. `T_3^n` has trace `2` on even `n` and `0` on odd `n`
       (involution `T_3^2 = I`);
    3. `T_5^n` is left unevaluated (no spectral data of `T_5` is needed
       for the parity headline).

    Hence:

    * `n` odd  → `trace((T_30)^n) = 0` (the extension of `T30_trace_zero`
      to every odd power);
    * `n` even → `trace((T_30)^n) = 2 · trace(T_5^n)`. -/
theorem T30_pow_trace_summary (T5 : T5Like) (n : ℕ) :
    -- (1) Three-factor decomposition
    ((T30 T5) ^ n).trace
      = (T2_trivial ^ n).trace * (T3 ^ n).trace * (T5.matrix ^ n).trace
    -- (2) Parity dichotomy
    ∧ ((∃ k, n = 2 * k ∧ ((T30 T5) ^ n).trace = 2 * (T5.matrix ^ n).trace)
       ∨ (∃ k, n = 2 * k + 1 ∧ ((T30 T5) ^ n).trace = 0)) :=
  ⟨T30_pow_trace_eq T5 n, T30_pow_trace_dichotomy T5 n⟩

end PT.Stochastic
