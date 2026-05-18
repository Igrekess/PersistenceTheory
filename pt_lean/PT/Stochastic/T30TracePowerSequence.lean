/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Stochastic.T30TraceFormulaExtended

/-!
# T2 — Explicit trace sequence `(tr T_30^n)_{n = 1..6}`

**Statement (paper-level, monograph ch03 §T2, follow-up to
`T30TraceFormulaExtended`).** The parity formulas

  `tr(T_30^{2k+1}) = 0`,    `tr(T_30^{2k})   = 2 · tr(T_5^{2k})`

(theorems `T30_pow_trace_odd` and `T30_pow_trace_even`) determine the
trace of every power of `T_30` from the corresponding power of `T_5`.
This file unrolls the first six values explicitly and packages them as
a headline sequence and a closed-form sum.

The sequence is parametric in the auxiliary structure `T5Like` (the
abstract `T_5` block) — no spectral or normalisation hypothesis on
`T_5` is used.

## Main results

### Individual values
* `T30_pow1_trace` — `tr(T_30^1) = 0`.
* `T30_pow2_trace` — `tr(T_30^2) = 2 · tr(T_5^2)`.
* `T30_pow3_trace` — `tr(T_30^3) = 0`.
* `T30_pow4_trace` — `tr(T_30^4) = 2 · tr(T_5^4)`.
* `T30_pow5_trace` — `tr(T_30^5) = 0`.
* `T30_pow6_trace` — `tr(T_30^6) = 2 · tr(T_5^6)`.

### Headline sequence
* `T30_trace_sequence_six` — packs the six values into one statement.

### Closed-form sum
* `T30_trace_sum_six` — telescoping of the three odd terms:
  `Σ_{n=1..6} tr(T_30^n) = 2 · (tr(T_5^2) + tr(T_5^4) + tr(T_5^6))`.

### Auxiliary inequality
* `abs_tr_le_self_of_nonneg` — trivial `|x| ≤ x` upper-bound when
  `x ≥ 0`, included for downstream convenience.

## Status

`[THM]` — every result reduces to one of the two parity formulas of
`T30TraceFormulaExtended` evaluated at `k ∈ {0, 1, 2, 3}`. No new
spectral input on `T_5`.

## References

* `PT.Stochastic.T30TraceFormulaExtended` — parity formulas
  `T30_pow_trace_odd` (k arbitrary) and `T30_pow_trace_even`.
-/

namespace PT.Stochastic

open Matrix BigOperators

/-! ### The six explicit values `(tr T_30^n)_{n = 1..6}` -/

/-- **`tr(T_30^1) = 0`.** Specialisation of `T30_pow_trace_odd` at
    `k = 0` (so `2k + 1 = 1`). -/
theorem T30_pow1_trace (T5 : T5Like) : ((T30 T5) ^ 1).trace = 0 := by
  have h := T30_pow_trace_odd T5 0
  simpa using h

/-- **`tr(T_30^2) = 2 · tr(T_5^2)`.** Specialisation of
    `T30_pow_trace_even` at `k = 1`. -/
theorem T30_pow2_trace (T5 : T5Like) :
    ((T30 T5) ^ 2).trace = 2 * (T5.matrix ^ 2).trace := by
  have h := T30_pow_trace_even T5 1
  simpa using h

/-- **`tr(T_30^3) = 0`.** Specialisation of `T30_pow_trace_odd` at
    `k = 1` (so `2k + 1 = 3`). -/
theorem T30_pow3_trace (T5 : T5Like) : ((T30 T5) ^ 3).trace = 0 := by
  have h := T30_pow_trace_odd T5 1
  simpa using h

/-- **`tr(T_30^4) = 2 · tr(T_5^4)`.** Specialisation of
    `T30_pow_trace_even` at `k = 2`. -/
theorem T30_pow4_trace (T5 : T5Like) :
    ((T30 T5) ^ 4).trace = 2 * (T5.matrix ^ 4).trace := by
  have h := T30_pow_trace_even T5 2
  simpa using h

/-- **`tr(T_30^5) = 0`.** Specialisation of `T30_pow_trace_odd` at
    `k = 2` (so `2k + 1 = 5`). -/
theorem T30_pow5_trace (T5 : T5Like) : ((T30 T5) ^ 5).trace = 0 := by
  have h := T30_pow_trace_odd T5 2
  simpa using h

/-- **`tr(T_30^6) = 2 · tr(T_5^6)`.** Specialisation of
    `T30_pow_trace_even` at `k = 3`. -/
theorem T30_pow6_trace (T5 : T5Like) :
    ((T30 T5) ^ 6).trace = 2 * (T5.matrix ^ 6).trace := by
  have h := T30_pow_trace_even T5 3
  simpa using h

/-! ### Headline sequence `(tr T_30^n)_{n = 1..6}` -/

/-- **Explicit trace sequence for `n = 1, …, 6`.** Headline conjunction
    bundling the six individual results. Reads off the alternation
    `0, 2·tr(T_5^2), 0, 2·tr(T_5^4), 0, 2·tr(T_5^6)`. -/
theorem T30_trace_sequence_six (T5 : T5Like) :
    ((T30 T5) ^ 1).trace = 0
    ∧ ((T30 T5) ^ 2).trace = 2 * (T5.matrix ^ 2).trace
    ∧ ((T30 T5) ^ 3).trace = 0
    ∧ ((T30 T5) ^ 4).trace = 2 * (T5.matrix ^ 4).trace
    ∧ ((T30 T5) ^ 5).trace = 0
    ∧ ((T30 T5) ^ 6).trace = 2 * (T5.matrix ^ 6).trace :=
  ⟨T30_pow1_trace T5, T30_pow2_trace T5, T30_pow3_trace T5,
   T30_pow4_trace T5, T30_pow5_trace T5, T30_pow6_trace T5⟩

/-! ### Closed-form sum of the six traces

The three odd terms vanish, so the sum collapses to `2 ·` the sum of
the three even-power `T_5` traces. -/

/-- **Sum of the six traces.** The three odd terms are zero, leaving
    `Σ_{n=1..6} tr(T_30^n) = 2 · (tr(T_5^2) + tr(T_5^4) + tr(T_5^6))`. -/
theorem T30_trace_sum_six (T5 : T5Like) :
    ((T30 T5) ^ 1).trace + ((T30 T5) ^ 2).trace + ((T30 T5) ^ 3).trace
      + ((T30 T5) ^ 4).trace + ((T30 T5) ^ 5).trace + ((T30 T5) ^ 6).trace
      = 2 * ((T5.matrix ^ 2).trace + (T5.matrix ^ 4).trace
              + (T5.matrix ^ 6).trace) := by
  rw [T30_pow1_trace, T30_pow2_trace, T30_pow3_trace,
      T30_pow4_trace, T30_pow5_trace, T30_pow6_trace]
  ring

/-! ### Auxiliary inequality

A trivial `|x| ≤ x` upper bound under nonnegativity, useful for
downstream estimates on `(T5.matrix^k).trace` when the trace is known
to be nonnegative (e.g. for symmetric positive-semidefinite blocks). -/

/-- **Trivial absolute-value bound under nonnegativity.** If `0 ≤ x`
    then `|x| ≤ x` (in fact `|x| = x`). -/
theorem abs_tr_le_self_of_nonneg (x : ℝ) (hx : 0 ≤ x) : |x| ≤ x := by
  rw [abs_of_nonneg hx]

end PT.Stochastic
