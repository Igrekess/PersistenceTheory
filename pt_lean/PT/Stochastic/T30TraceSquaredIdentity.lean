/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Stochastic.T30TracePowerSequence

/-!
# T2 — Algebraic identities `tr(T_30^a) · tr(T_30^b) = 0` and consequences

**Statement (paper-level, monograph ch03 §T2, follow-up to
`T30TracePowerSequence`).** Since `tr(T_30) = 0` (cf.
`PT.Stochastic.T30TraceDeterminant.T30_trace_zero`) and, more generally,
`tr(T_30^(2k+1)) = 0` for every `k` (cf.
`PT.Stochastic.T30TraceFormulaExtended.T30_pow_trace_odd`), every product
of the form `tr(T_30^a) · tr(T_30^b)` in which at least one of the two
exponents is odd vanishes identically.

This file packages this *vanishing identity* and its immediate
consequences:

1. the squared / iterated identities `tr(T_30^2) · tr(T_30) = 0`,
   `tr(T_30^k) · tr(T_30) = 0`, and the symmetric variants;
2. the **bilinear odd-annihilation** `tr(T_30^(2k+1)) · tr(T_30^m) = 0`;
3. the **Newton-Girard residue** for the eigenvalues `(λ_i)` of `T_30`:
   `(Σ λ_i)^2 = (tr T_30)^2 = 0` while `Σ λ_i^2 = tr(T_30^2) =
   2 · tr(T_5^2)`, so the *power-sum residue*
   `Σ λ_i^2 − (Σ λ_i)^2 = 2 · tr(T_5^2)` is the entire variance of the
   spectrum;
4. the **scalar matrix identity** `T_30 · tr(T_30^2) = T_30 ·
   (2 · tr(T_5^2))` (a trivial consequence of `T30_pow_two_trace`,
   recorded for downstream matrix-level use).

All proofs are one- or two-line reductions to `T30_trace_zero` and to
the parity formulas of `T30TraceFormulaExtended`. **No new arithmetic
input is required** — the file is a structural reorganisation of the
existing trace identities into a *bilinear-vanishing* form that is
convenient to cite when manipulating products of traces (e.g. in
characteristic-polynomial calculations or moment-method estimates).

## Main results

### Direct squared identity
* `T30_trace_sq_mul_trace_zero` — `tr(T_30^2) · tr(T_30) = 0`.
* `T30_trace_mul_trace_sq_zero` — `tr(T_30) · tr(T_30^2) = 0`
  (symmetric form).

### Bilinear odd-annihilation
* `T30_pow_trace_odd_mul_any` — for every odd exponent `2k+1` and every
  `m`, `tr(T_30^(2k+1)) · tr(T_30^m) = 0`.
* `T30_any_mul_pow_trace_odd` — symmetric form.
* `T30_trace_mul_pow_zero` — for every `m`, `tr(T_30) · tr(T_30^m) = 0`.
* `T30_pow_mul_trace_zero` — symmetric form.

### Newton-Girard residue (power-sum interpretation)
* `T30_power_sum_residue` — `tr(T_30^2) − (tr T_30)^2 = 2 · tr(T_5^2)`.
  In eigenvalue language, `Σ λ_i^2 − (Σ λ_i)^2 = 2 · tr(T_5^2)`.

### Scalar matrix identity
* `T30_smul_trace_sq` — `tr(T_30^2) • T_30 = (2 · tr(T_5^2)) • T_30`
  (trivial consequence of `T30_pow_two_trace`).

### Headline summary
* `T30_trace_squared_invariants_summary` — packages the four facts
  above in one statement.

## Strategy

Every product identity reduces via `T30_trace_zero` (for `m = 1`) or
`T30_pow_trace_odd` (for general odd exponents) to `_ * 0 = 0` or
`0 * _ = 0`. The Newton-Girard residue is one rewrite using
`T30_pow_two_trace` and `T30_trace_zero`. The scalar matrix identity is
a `congrArg` on the scalar.

## Status

`[THM]` — purely algebraic restatement of the existing trace vanishing
results. No additional spectral hypothesis on `T_5`; the residue
identity carries the value `2 · tr(T_5^2)` symbolically.

## References

* `PT.Stochastic.T30TraceDeterminant` — `T30_trace_zero`.
* `PT.Stochastic.T30TraceFormulaExtended` — `T30_pow_trace_odd`,
  `T30_pow_trace_even`, `T30_pow_two_trace`.
* `PT.Stochastic.T30TracePowerSequence` — `T30_pow1_trace`,
  `T30_pow2_trace`.
-/

namespace PT.Stochastic

open Matrix BigOperators

/-! ### Direct squared identity `tr(T_30^2) · tr(T_30) = 0` -/

/-- **`tr(T_30^2) · tr(T_30) = 0`.**

    Immediate from `T30_trace_zero`: the right-hand factor is `0`, so
    the product is `0`. This is the cleanest *bilinear* form of the
    headline `tr(T_30) = 0`. -/
theorem T30_trace_sq_mul_trace_zero (T5 : T5Like) :
    ((T30 T5) ^ 2).trace * (T30 T5).trace = 0 := by
  rw [T30_trace_zero]
  ring

/-- **Symmetric form: `tr(T_30) · tr(T_30^2) = 0`.** -/
theorem T30_trace_mul_trace_sq_zero (T5 : T5Like) :
    (T30 T5).trace * ((T30 T5) ^ 2).trace = 0 := by
  rw [T30_trace_zero]
  ring

/-! ### Bilinear odd-annihilation -/

/-- **Bilinear odd-annihilation.** For every odd exponent `2k + 1` and
    every natural exponent `m`,

    `tr((T_30)^(2k+1)) · tr((T_30)^m) = 0`.

    Direct consequence of `T30_pow_trace_odd` (the left factor vanishes).
    Covers `T30_trace_sq_mul_trace_zero` at `k = 0`, `m = 2`. -/
theorem T30_pow_trace_odd_mul_any (T5 : T5Like) (k m : ℕ) :
    ((T30 T5) ^ (2 * k + 1)).trace * ((T30 T5) ^ m).trace = 0 := by
  rw [T30_pow_trace_odd]
  ring

/-- **Symmetric form: any · odd.** -/
theorem T30_any_mul_pow_trace_odd (T5 : T5Like) (m k : ℕ) :
    ((T30 T5) ^ m).trace * ((T30 T5) ^ (2 * k + 1)).trace = 0 := by
  rw [T30_pow_trace_odd]
  ring

/-- **`tr(T_30) · tr(T_30^m) = 0` for every `m`.**

    Special case of `T30_pow_trace_odd_mul_any` at `k = 0`, with the
    rewrite `2 * 0 + 1 = 1` made implicit. -/
theorem T30_trace_mul_pow_zero (T5 : T5Like) (m : ℕ) :
    (T30 T5).trace * ((T30 T5) ^ m).trace = 0 := by
  rw [T30_trace_zero]
  ring

/-- **Symmetric form: `tr(T_30^m) · tr(T_30) = 0` for every `m`.** -/
theorem T30_pow_mul_trace_zero (T5 : T5Like) (m : ℕ) :
    ((T30 T5) ^ m).trace * (T30 T5).trace = 0 := by
  rw [T30_trace_zero]
  ring

/-- **Triple-odd product vanishes.** Any product involving two odd
    powers also vanishes (only one zero factor is needed, but recording
    the iterated form is convenient). -/
theorem T30_triple_odd_mul (T5 : T5Like) (j k m : ℕ) :
    ((T30 T5) ^ (2 * j + 1)).trace
      * ((T30 T5) ^ (2 * k + 1)).trace
      * ((T30 T5) ^ m).trace = 0 := by
  rw [T30_pow_trace_odd, T30_pow_trace_odd]
  ring

/-! ### Newton-Girard residue (power-sum interpretation)

If `(λ_i)_{i = 1..30}` are the eigenvalues of `T_30` (counted with
algebraic multiplicity), Newton's identities relate the power sums

  `p_1  =  Σ λ_i        =  tr(T_30)`,
  `p_2  =  Σ λ_i^2      =  tr(T_30^2)`,

to the elementary symmetric polynomials. The **squared first power sum**
`(p_1)^2 = (tr T_30)^2 = 0` and the **second power sum** `p_2 =
tr(T_30^2) = 2 · tr(T_5^2)` (by `T30_pow_two_trace`). The *residue*

  `p_2 − (p_1)^2  =  Σ λ_i^2 − (Σ λ_i)^2  =  2 · tr(T_5^2)`

is therefore the entire **variance of the spectrum** of `T_30` (the
second elementary symmetric polynomial in disguise), reduced to a single
trace on the `T_5` block. -/

/-- **Newton-Girard residue for `T_30`.**

    `tr(T_30^2) − (tr T_30)^2 = 2 · tr(T_5^2)`.

    In eigenvalue terms (for the multiset `(λ_i)` of eigenvalues of
    `T_30`): `Σ λ_i^2 − (Σ λ_i)^2 = 2 · tr(T_5^2)`. This is the
    *variance-like* invariant of the `T_30` spectrum, expressed entirely
    in the `T_5` block. -/
theorem T30_power_sum_residue (T5 : T5Like) :
    ((T30 T5) ^ 2).trace - ((T30 T5).trace) ^ 2 = 2 * (T5.matrix ^ 2).trace := by
  rw [T30_pow_two_trace, T30_trace_zero]
  ring

/-- **Equivalent form: `(tr T_30)^2 = 0`.**

    The squared first power sum vanishes, so the residue identity
    `T30_power_sum_residue` is `tr(T_30^2) = 2 · tr(T_5^2)`. -/
theorem T30_trace_sq_eq_zero (T5 : T5Like) :
    ((T30 T5).trace) ^ 2 = 0 := by
  rw [T30_trace_zero]
  ring

/-- **Variance-like residue, multiplicative form.**

    `tr(T_30^2) = (tr T_30)^2 + 2 · tr(T_5^2)`. A trivial rearrangement
    of `T30_power_sum_residue`, useful when one wants `tr(T_30^2)` on
    the left-hand side. -/
theorem T30_pow_two_trace_decomp (T5 : T5Like) :
    ((T30 T5) ^ 2).trace = ((T30 T5).trace) ^ 2 + 2 * (T5.matrix ^ 2).trace := by
  rw [T30_trace_zero, T30_pow_two_trace]
  ring

/-! ### Scalar matrix identity

The scalar `tr(T_30^2) = 2 · tr(T_5^2)` (from `T30_pow_two_trace`) can
be applied as a scalar multiple of `T_30`. The two `smul` expressions
agree pointwise. -/

/-- **Scalar matrix identity.**

    `tr(T_30^2) • T_30 = (2 · tr(T_5^2)) • T_30`.

    Direct rewrite via `T30_pow_two_trace`. Useful when one wants to
    factor `T_30` against its second moment in a matrix expression. -/
theorem T30_smul_trace_sq (T5 : T5Like) :
    ((T30 T5) ^ 2).trace • (T30 T5) = (2 * (T5.matrix ^ 2).trace) • (T30 T5) := by
  rw [T30_pow_two_trace]

/-! ### Headline summary -/

/-- **Headline summary for the bilinear trace identities.**

    The vanishing `tr(T_30) = 0` propagates to a *bilinear-vanishing*
    family of identities, plus a residue interpretation:

    1. **Direct squared identity** — `tr(T_30^2) · tr(T_30) = 0`.
    2. **Bilinear odd-annihilation** — for every odd exponent `2k+1`
       and every `m`, `tr(T_30^(2k+1)) · tr(T_30^m) = 0`.
    3. **Newton-Girard residue** — `tr(T_30^2) − (tr T_30)^2 =
       2 · tr(T_5^2)`, i.e. the variance of the `T_30` spectrum reduces
       entirely to the `T_5` block.
    4. **Scalar matrix identity** — `tr(T_30^2) • T_30 =
       (2 · tr(T_5^2)) • T_30`. -/
theorem T30_trace_squared_invariants_summary (T5 : T5Like) :
    -- (1) Direct squared identity
    ((T30 T5) ^ 2).trace * (T30 T5).trace = 0
    -- (2) Bilinear odd-annihilation (k = 0 form, any m)
    ∧ (∀ m : ℕ, (T30 T5).trace * ((T30 T5) ^ m).trace = 0)
    -- (3) Newton-Girard residue
    ∧ ((T30 T5) ^ 2).trace - ((T30 T5).trace) ^ 2 = 2 * (T5.matrix ^ 2).trace
    -- (4) Scalar matrix identity
    ∧ ((T30 T5) ^ 2).trace • (T30 T5) = (2 * (T5.matrix ^ 2).trace) • (T30 T5) :=
  ⟨T30_trace_sq_mul_trace_zero T5,
   T30_trace_mul_pow_zero T5,
   T30_power_sum_residue T5,
   T30_smul_trace_sq T5⟩

end PT.Stochastic
