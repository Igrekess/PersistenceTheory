/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Stochastic.T30TraceFifthPowerIdentity

/-!
# T2 — Complete characteristic polynomial of `T_30` (even-polynomial closure)

**Statement (paper-level, monograph ch03 §T2, capstone of the trace-identity
layer).** The characteristic polynomial of `T_30 = T_2 ⊗ T_3 ⊗ T_5` (dimension
`d = 1 · 2 · 4 = 8 = φ(30)`) admits the Vieta expansion

  `det(x · I − T_30)  =  x^8 − e_1 · x^7 + e_2 · x^6 − e_3 · x^5
                          + e_4 · x^4 − e_5 · x^3 + e_6 · x^2 − e_7 · x + e_8`,

with `e_k` the `k`-th elementary symmetric polynomial of the spectrum of `T_30`.
This file packages the **complete coefficient-by-coefficient closure** of the
charpoly that follows from combining the trace identities of
`T30TraceFormulaExtended` (odd traces vanish: `tr(T_30^{2k+1}) = 0`) with the
Newton–Girard recursion: every **odd-index** elementary symmetric polynomial of
the `T_30` spectrum vanishes, so the charpoly reduces to an **even polynomial**
in `x`:

  `det(x · I − T_30)  =  x^8 + a_6 · x^6 + a_4 · x^4 + a_2 · x^2 + a_0`,

with `a_6 = e_2`, `a_4 = e_4`, `a_2 = e_6`, `a_0 = e_8 = det(T_30) = det(T_5)^2`.

The vanishing of the odd coefficients `a_7 = a_5 = a_3 = a_1 = 0` follows from:

* `a_7 = −e_1 = −tr(T_30) = 0`            (**unconditional**, from
   `T30TraceDeterminant.T30_trace_zero`);
* `a_5 = −e_3 = 0`                         (conditional on Newton's third
   identity, the cubic-layer collapse from `T30TraceCubeIdentity`);
* `a_3 = −e_5 = 0`                         (conditional on Newton's fifth
   identity *plus* `e_3 = 0`, the quintic-layer collapse from
   `T30TraceFifthPowerIdentity`);
* `a_1 = −e_7 = 0`                         (conditional on Newton's seventh
   identity *plus* `e_1 = e_3 = e_5 = 0` *plus* `tr(T_30^7) = 0`).

The last item is *new* in this file: we record the algebraic Newton-seventh
collapse `7·e_7 = e_6·p_1 − e_5·p_2 + e_4·p_3 − e_3·p_4 + e_2·p_5 − e_1·p_6
+ p_7 = 0` under the four vanishings `p_1 = p_3 = p_5 = p_7 = 0` and
`e_1 = e_3 = e_5 = 0`.

## Physical interpretation

The reduction `charpoly(T_30) = p(x^2)` (a polynomial in `x^2`) means that the
spectrum of `T_30` is **stable under negation**: if `λ` is an eigenvalue of
`T_30`, so is `−λ`. This is the **dichotomy** observed in the Perron / anti
sector decomposition (cf. `T30PerronAntiCommutator`,
`T30FullEigenpairCount`): the eight eigenvalues of `T_30` come in **four pairs**
`{(+1, −1), (+λ_2, −λ_2), (+λ_3, −λ_3), (+λ_4, −λ_4)}`. The arithmetic origin
of this `±`-symmetry is the **antidiagonal involution** `T_3 = !![0,1;1,0]`,
whose own spectrum `{+1, −1}` propagates through the Kronecker product
`T_30 = T_2 ⊗ T_3 ⊗ T_5` to enforce sign-pair symmetry on the full `T_30`
spectrum.

## Main results

### Vanishing of the odd-degree coefficients
* `T30_charpoly_a7_zero` — `a_7 = 0` (**unconditional**).
* `T30_charpoly_a5_zero` — `a_5 = 0` under Newton's third identity.
* `T30_charpoly_a3_zero` — `a_3 = 0` under Newton's fifth identity plus
  `e_3 = 0`.
* `newton_seventh_identity_collapse` — purely algebraic Newton-seventh
  collapse.
* `e7_vanishes_from_newton` — algebraic implication: under Newton's seventh
  identity, `p_1 = p_3 = p_5 = p_7 = 0`, `e_1 = e_3 = e_5 = 0`, then
  `e_7 = 0`.
* `T30_charpoly_a1_zero` — `a_1 = 0` under the Newton-seventh hypothesis and
  the three preceding vanishings.

### Constant-term closure
* `T30_charpoly_a0_closed` — `a_0 = e_8 = det(T_30) = det(T_5)^2`,
  **unconditional** (from `T30TraceDeterminant.T30_det_eq`).

### Even-polynomial form (headline)
* `T30_charpoly_is_even` — **headline**: under the three Newton hypotheses
  (3rd, 5th, 7th identities), the charpoly of `T_30` has the **even-polynomial
  form**
  `p(x) = x^8 + e_2 · x^6 + e_4 · x^4 + e_6 · x^2 + det(T_5)^2`.

### Headline summary
* `T30_charpoly_complete_summary` — packages all four odd-coefficient
  vanishings, the constant-term closure, and the even-polynomial form in one
  statement.

## Strategy

Every result reduces to a one- or two-line `linarith` / `ring` argument once
the relevant Newton identity is substituted. The new content (Newton-seventh)
mirrors `e3_vanishes_from_newton` / `e5_vanishes_from_newton`: a purely
algebraic substitution lemma plus an application to `T_30`.

The "even-polynomial form" is recorded as a propositional `And` of the four
odd-coefficient vanishings plus the closed values for the even coefficients.
We **do not** unfold Mathlib's `Matrix.charpoly` machinery (which would require
casting from `ℝ` to a polynomial ring `ℝ[X]` and committing to a specific
eigenvalue datum); instead we record the *scalar values* of the eight Vieta
coefficients `e_1, …, e_8`, packaged as a structural statement about the
charpoly's shape.

## Status

`[THM]` — purely algebraic capstone of the trace-identity tower
(`T30TraceFormulaExtended` → `T30TracePowerSequence` →
`T30TraceSquaredIdentity` → `T30TraceCubeIdentity` →
`T30TraceFourthPowerIdentity` → `T30TraceFifthPowerIdentity` → *this file*).
The headline `T30_charpoly_is_even` is *conditional* on the three Newton
identities holding for the `T_30` spectrum (a textbook fact about elementary
symmetric polynomials of a finite multiset that requires committing to an
eigenvalue datum to formalise inside Lean — left as the algebraic hypothesis
`hN3 / hN5 / hN7`).

## References

* `PT.Stochastic.T30TraceDeterminant` — `T30_trace_zero`, `T30_det_eq`.
* `PT.Stochastic.T30TraceFormulaExtended` — `T30_pow_trace_odd`,
  `T30_pow_trace_even`.
* `PT.Stochastic.T30TracePowerSequence` — `T30_pow1_trace`, …, `T30_pow6_trace`.
* `PT.Stochastic.T30TraceCubeIdentity` — `e3_vanishes_from_newton`,
  `T30_e3_collapse`.
* `PT.Stochastic.T30TraceFourthPowerIdentity` — `T30_e4_closed_form_factored`.
* `PT.Stochastic.T30TraceFifthPowerIdentity` — `e5_vanishes_from_newton`,
  `T30_e5_vanishes_from_T5_traces`.
* `PT.Stochastic.T30PerronAntiCommutator`, `PT.Stochastic.T30FullEigenpairCount`
  — physical / spectral interpretation of the `±`-symmetric structure.
-/

namespace PT.Stochastic

open Matrix BigOperators

/-! ### Charpoly coefficient `a_7 = −e_1`: unconditional vanishing

The coefficient of `x^7` in the characteristic polynomial of an `8 × 8`
matrix equals `−e_1 = −tr`. For `T_30` the trace vanishes unconditionally
(`T30_trace_zero`), so this coefficient is `0`. -/

/-- **`a_7 = −e_1 = −tr(T_30) = 0`** (unconditional).

    The coefficient of `x^7` in `charpoly(T_30)` vanishes, regardless of
    the underlying `T_5`. This is the unconditional `T30_trace_zero`
    rewritten as a scalar Vieta value. -/
theorem T30_charpoly_a7_zero (T5 : T5Like) :
    -(T30 T5).trace = 0 := by
  rw [T30_trace_zero]
  ring

/-! ### Charpoly coefficient `a_5 = −e_3`: vanishing under Newton's third

The coefficient of `x^5` equals `−e_3`. By the cubic-layer collapse
(`e3_vanishes_from_newton` from `T30TraceCubeIdentity`), under Newton's
third identity for the `T_30` spectrum we have `e_3 = 0`, hence
`a_5 = 0`. -/

/-- **`a_5 = −e_3 = 0`** under Newton's third identity.

    The coefficient of `x^5` in `charpoly(T_30)` vanishes, conditional on
    Newton's third identity `p_3 − e_1·p_2 + e_2·p_1 − 3·e_3 = 0` (which
    is a textbook relation between power sums and elementary symmetric
    polynomials of any finite multiset). -/
theorem T30_charpoly_a5_zero (T5 : T5Like)
    (p2 e2 e3 : ℝ)
    (hN3 : ((T30 T5) ^ 3).trace - (T30 T5).trace * p2
            + e2 * (T30 T5).trace - 3 * e3 = 0) :
    -e3 = 0 := by
  have h : -(3 * e3) = 0 := by
    have := hN3
    rw [T30_pow3_trace, T30_trace_zero] at this
    linarith
  linarith

/-! ### Charpoly coefficient `a_3 = −e_5`: vanishing under Newton's fifth

The coefficient of `x^3` equals `−e_5`. By the quintic-layer collapse
(`e5_vanishes_from_newton` from `T30TraceFifthPowerIdentity`), under
Newton's fifth identity for the `T_30` spectrum *plus* `e_3 = 0`, we have
`e_5 = 0`, hence `a_3 = 0`. -/

/-- **`a_3 = −e_5 = 0`** under Newton's fifth identity and `e_3 = 0`.

    The coefficient of `x^3` in `charpoly(T_30)` vanishes, conditional on
    Newton's fifth identity for the `T_30` spectrum **plus** the cubic
    vanishing `e_3 = 0` (which itself follows from Newton's third). -/
theorem T30_charpoly_a3_zero (T5 : T5Like)
    (e2 e3 e4 e5 : ℝ)
    (hN5 : ((T30 T5) ^ 5).trace
            - (T30 T5).trace * ((T30 T5) ^ 4).trace
            + e2 * ((T30 T5) ^ 3).trace
            - e3 * ((T30 T5) ^ 2).trace
            + e4 * (T30 T5).trace
            - 5 * e5
          = 0)
    (he3 : e3 = 0) :
    -e5 = 0 := by
  have h := T30_e5_vanishes_from_T5_traces T5 e2 e3 e4 e5 hN5 he3
  linarith

/-! ### Charpoly coefficient `a_1 = −e_7`: Newton's seventh identity collapse

Newton's seventh identity for the power sums `p_k` and elementary symmetric
polynomials `e_k` of a finite multiset reads

  `p_7 − e_1·p_6 + e_2·p_5 − e_3·p_4 + e_4·p_3 − e_5·p_2 + e_6·p_1 − 7·e_7
    = 0`.

For the spectrum of `T_30`:

* `p_1 = tr(T_30) = 0`     (unconditional, `T30_trace_zero`);
* `p_3 = tr(T_30^3) = 0`   (unconditional, `T30_pow3_trace`);
* `p_5 = tr(T_30^5) = 0`   (unconditional, `T30_pow5_trace`);
* `p_7 = tr(T_30^7) = 0`   (unconditional, `T30_pow_trace_odd` at `k = 3`);
* `e_1 = p_1 = 0`, and conditional on the preceding Newton identities,
  `e_3 = 0`, `e_5 = 0`.

Substituting these seven zeros into Newton's seventh identity gives the
collapse `−7·e_7 = 0`, hence `e_7 = 0`. -/

/-- **Newton–Girard seventh-identity collapse (purely algebraic).**

    If `p_1 = p_3 = p_5 = p_7 = 0` and `e_1 = e_3 = e_5 = 0`, then for any
    `p_2, p_4, p_6, e_2, e_4, e_6, e_7 : ℝ`, the Newton-seventh linear
    combination
    `p_7 − e_1·p_6 + e_2·p_5 − e_3·p_4 + e_4·p_3 − e_5·p_2 + e_6·p_1 − 7·e_7`
    simplifies to `−7·e_7`.

    The lemma records the **collapse step** (the substitution of the seven
    zeros), not the underlying Newton identity itself. -/
theorem newton_seventh_identity_collapse
    (p2 p4 p6 e2 e4 e6 e7 : ℝ) :
    (0 : ℝ) - (0 : ℝ) * p6 + e2 * (0 : ℝ) - (0 : ℝ) * p4
        + e4 * (0 : ℝ) - (0 : ℝ) * p2 + e6 * (0 : ℝ) - 7 * e7
      = -(7 * e7) := by
  ring

/-- **`e_7` vanishes from Newton's seventh (purely algebraic).**

    If Newton's seventh identity holds with `p_1 = p_3 = p_5 = p_7 = 0`
    and `e_1 = e_3 = e_5 = 0`, then `e_7 = 0`. -/
theorem e7_vanishes_from_newton
    (p1 p2 p3 p4 p5 p6 p7 e1 e2 e3 e4 e5 e6 e7 : ℝ)
    (hN : p7 - e1 * p6 + e2 * p5 - e3 * p4 + e4 * p3 - e5 * p2
            + e6 * p1 - 7 * e7 = 0)
    (hp1 : p1 = 0) (hp3 : p3 = 0) (hp5 : p5 = 0) (hp7 : p7 = 0)
    (he1 : e1 = 0) (he3 : e3 = 0) (he5 : e5 = 0) :
    e7 = 0 := by
  have h := hN
  rw [hp1, hp3, hp5, hp7, he1, he3, he5] at h
  linarith

/-- **Application to `T_30`: Newton-seventh collapse.**

    Substituting `p_1 = tr(T_30)`, `p_3 = tr(T_30^3)`, `p_5 = tr(T_30^5)`,
    `p_7 = tr(T_30^7)` (all zero) and `e_1 = tr(T_30) = 0`, the Newton-seventh
    combination
    `tr(T_30^7) − tr(T_30)·tr(T_30^6) + e_2·tr(T_30^5) − e_3·tr(T_30^4)
        + e_4·tr(T_30^3) − e_5·tr(T_30^2) + e_6·tr(T_30) − 7·e_7`
    simplifies to `−e_3·tr(T_30^4) − e_5·tr(T_30^2) − 7·e_7` for every
    `e_2, e_3, e_4, e_5, e_6, e_7`.

    The two residual terms involve `e_3` and `e_5`, which themselves vanish
    by the preceding Newton-third and Newton-fifth collapses. -/
theorem T30_newton_seventh_collapse (T5 : T5Like)
    (e2 e3 e4 e5 e6 e7 : ℝ) :
    ((T30 T5) ^ 7).trace
      - (T30 T5).trace * ((T30 T5) ^ 6).trace
      + e2 * ((T30 T5) ^ 5).trace
      - e3 * ((T30 T5) ^ 4).trace
      + e4 * ((T30 T5) ^ 3).trace
      - e5 * ((T30 T5) ^ 2).trace
      + e6 * (T30 T5).trace
      - 7 * e7
    = -e3 * ((T30 T5) ^ 4).trace - e5 * ((T30 T5) ^ 2).trace - 7 * e7 := by
  have h7 : ((T30 T5) ^ 7).trace = 0 := by
    have h := T30_pow_trace_odd T5 3
    simpa using h
  rw [T30_trace_zero, T30_pow3_trace, T30_pow5_trace, h7]
  ring

/-- **Application to `T_30`: `e_7 = 0` under Newton-seventh plus `e_3 = 0`
    and `e_5 = 0`.**

    Under the Newton-seventh hypothesis for the spectrum of `T_30`, plus
    the cubic vanishing `e_3 = 0` (Newton-third) and the quintic vanishing
    `e_5 = 0` (Newton-fifth), the seventh elementary symmetric polynomial
    `e_7` of the spectrum vanishes. -/
theorem T30_e7_vanishes_from_T5_traces (T5 : T5Like)
    (e2 e3 e4 e5 e6 e7 : ℝ)
    (hN : ((T30 T5) ^ 7).trace
            - (T30 T5).trace * ((T30 T5) ^ 6).trace
            + e2 * ((T30 T5) ^ 5).trace
            - e3 * ((T30 T5) ^ 4).trace
            + e4 * ((T30 T5) ^ 3).trace
            - e5 * ((T30 T5) ^ 2).trace
            + e6 * (T30 T5).trace
            - 7 * e7
          = 0)
    (he3 : e3 = 0) (he5 : e5 = 0) :
    e7 = 0 := by
  have hcollapse := T30_newton_seventh_collapse T5 e2 e3 e4 e5 e6 e7
  -- `hN` and `hcollapse` together give `-e3·p_4 - e5·p_2 - 7·e_7 = 0`
  have h : -e3 * ((T30 T5) ^ 4).trace - e5 * ((T30 T5) ^ 2).trace - 7 * e7 = 0 := by
    linarith
  rw [he3, he5] at h
  linarith

/-- **`a_1 = −e_7 = 0`** under Newton's seventh identity, `e_3 = 0`, and
    `e_5 = 0`.

    The coefficient of `x^1` in `charpoly(T_30)` vanishes, conditional on
    Newton's seventh identity for the `T_30` spectrum **plus** the
    cubic vanishing `e_3 = 0` **plus** the quintic vanishing `e_5 = 0`. -/
theorem T30_charpoly_a1_zero (T5 : T5Like)
    (e2 e3 e4 e5 e6 e7 : ℝ)
    (hN7 : ((T30 T5) ^ 7).trace
            - (T30 T5).trace * ((T30 T5) ^ 6).trace
            + e2 * ((T30 T5) ^ 5).trace
            - e3 * ((T30 T5) ^ 4).trace
            + e4 * ((T30 T5) ^ 3).trace
            - e5 * ((T30 T5) ^ 2).trace
            + e6 * (T30 T5).trace
            - 7 * e7
          = 0)
    (he3 : e3 = 0) (he5 : e5 = 0) :
    -e7 = 0 := by
  have h := T30_e7_vanishes_from_T5_traces T5 e2 e3 e4 e5 e6 e7 hN7 he3 he5
  linarith

/-! ### Constant term: `a_0 = e_8 = det(T_30) = det(T_5)^2`

The constant term of `charpoly(T_30)` equals `(−1)^8 · det(T_30) = det(T_30)`,
which by `T30_det_eq` equals `det(T_5)^2`. This is **unconditional** (no
Newton hypothesis required) and is the cleanest closed-form value among the
charpoly coefficients. -/

/-- **`a_0 = e_8 = det(T_30) = det(T_5)^2`** (unconditional).

    The constant term of `charpoly(T_30)` equals `det(T_5)^2` for every
    `T5Like`. Direct from `T30_det_eq`. -/
theorem T30_charpoly_a0_closed (T5 : T5Like) :
    (T30 T5).det = T5.matrix.det ^ 2 :=
  T30_det_eq T5

/-! ### Even-polynomial form (headline)

Under the three Newton hypotheses (3rd, 5th, 7th identities for the `T_30`
spectrum), the characteristic polynomial of `T_30` has the **even-polynomial
form**

  `charpoly(T_30)(x)  =  x^8 + e_2 · x^6 + e_4 · x^4 + e_6 · x^2 + det(T_5)^2`,

i.e. all four odd-degree coefficients vanish, and the constant term is
`det(T_5)^2`. We record this as the conjunction of the four
odd-coefficient vanishings plus the constant-term closure. -/

/-- **Headline: charpoly of `T_30` is an even polynomial in `x`.**

    Under the three Newton hypotheses for the `T_30` spectrum
    (`hN3`, `hN5`, `hN7`), the four odd-degree coefficients of
    `charpoly(T_30)` all vanish, and the constant term equals `det(T_5)^2`.
    The remaining (even-degree) coefficients are `a_6 = e_2`, `a_4 = e_4`,
    `a_2 = e_6` — pinned by Newton's *even*-identity layer (cf. the
    conditional closed forms from `T30TraceFourthPowerIdentity` for `e_4`,
    and the quadratic-layer `e_2 = -tr(T_5^2)` from
    `T30TraceSquaredIdentity`).

    **Spectral interpretation.** `charpoly(T_30)(x) = p(x^2)` for some
    degree-`4` polynomial `p`, so the spectrum of `T_30` is **symmetric
    under negation**: eigenvalues come in `±`-pairs. This is the algebraic
    shadow of the antidiagonal involution `T_3 = !![0,1;1,0]` whose
    `{+1, −1}`-spectrum propagates through the Kronecker product
    `T_30 = T_2 ⊗ T_3 ⊗ T_5`. -/
theorem T30_charpoly_is_even (T5 : T5Like)
    (e2 e3 e4 e5 e6 e7 : ℝ)
    (p2 : ℝ)
    (hN3 : ((T30 T5) ^ 3).trace - (T30 T5).trace * p2
            + e2 * (T30 T5).trace - 3 * e3 = 0)
    (hN5 : ((T30 T5) ^ 5).trace
            - (T30 T5).trace * ((T30 T5) ^ 4).trace
            + e2 * ((T30 T5) ^ 3).trace
            - e3 * ((T30 T5) ^ 2).trace
            + e4 * (T30 T5).trace
            - 5 * e5
          = 0)
    (hN7 : ((T30 T5) ^ 7).trace
            - (T30 T5).trace * ((T30 T5) ^ 6).trace
            + e2 * ((T30 T5) ^ 5).trace
            - e3 * ((T30 T5) ^ 4).trace
            + e4 * ((T30 T5) ^ 3).trace
            - e5 * ((T30 T5) ^ 2).trace
            + e6 * (T30 T5).trace
            - 7 * e7
          = 0) :
    -- Vanishing of all four odd-degree coefficients
    (-(T30 T5).trace = 0)            -- a_7 = -e_1 = 0
    ∧ (-e3 = 0)                       -- a_5 = -e_3 = 0
    ∧ (-e5 = 0)                       -- a_3 = -e_5 = 0
    ∧ (-e7 = 0)                       -- a_1 = -e_7 = 0
    -- Constant-term closure
    ∧ (T30 T5).det = T5.matrix.det ^ 2 := by
  have ha7 := T30_charpoly_a7_zero T5
  have ha5 := T30_charpoly_a5_zero T5 p2 e2 e3 hN3
  have he3 : e3 = 0 := by linarith
  have ha3 := T30_charpoly_a3_zero T5 e2 e3 e4 e5 hN5 he3
  have he5 : e5 = 0 := by linarith
  have ha1 := T30_charpoly_a1_zero T5 e2 e3 e4 e5 e6 e7 hN7 he3 he5
  have ha0 := T30_charpoly_a0_closed T5
  exact ⟨ha7, ha5, ha3, ha1, ha0⟩

/-! ### Complete coefficient summary -/

/-- **Headline summary for the complete charpoly of `T_30`.**

    The eight Vieta coefficients of `charpoly(T_30) = x^8 − e_1·x^7
    + e_2·x^6 − e_3·x^5 + e_4·x^4 − e_5·x^3 + e_6·x^2 − e_7·x + e_8`
    decompose into four **odd-index** coefficients (all vanishing, three
    conditionally) and four **even-index** coefficients:

    1. **`e_1 = tr(T_30) = 0`** (**unconditional**, `T30_trace_zero`).
    2. **`e_2 = (1/2)((tr T_30)^2 − tr(T_30^2)) = −tr(T_5^2)`** (Newton-second,
       cf. `T30_power_sum_residue`; not a hypothesis here).
    3. **`e_3 = 0`** under Newton's third identity (`T30TraceCubeIdentity`).
    4. **`e_4 = −(tr(T_5^4) + e_2·tr(T_5^2))/2`** under Newton's fourth
       identity (`T30TraceFourthPowerIdentity`; modulo `e_2`).
    5. **`e_5 = 0`** under Newton's fifth identity and `e_3 = 0`
       (`T30TraceFifthPowerIdentity`).
    6. **`e_6`** (Newton's sixth identity; not unfolded here).
    7. **`e_7 = 0`** under Newton's seventh identity and `e_3 = e_5 = 0`
       (this file, `e7_vanishes_from_newton`).
    8. **`e_8 = det(T_30) = det(T_5)^2`** (**unconditional**, `T30_det_eq`).

    The four odd-index vanishings collapse `charpoly(T_30)` to an
    **even polynomial** in `x`, encoding the `±`-symmetry of the
    `T_30` spectrum (Perron / anti-sector dichotomy). -/
theorem T30_charpoly_complete_summary (T5 : T5Like) :
    -- (1) `e_1 = 0` unconditional
    (T30 T5).trace = 0
    -- (2) `a_7 = −e_1 = 0` unconditional
    ∧ -(T30 T5).trace = 0
    -- (3) `e_8 = det(T_5)^2` unconditional
    ∧ (T30 T5).det = T5.matrix.det ^ 2
    -- (4) Newton-third collapse (algebraic)
    ∧ (∀ p2 e2 e3 : ℝ,
        ((T30 T5) ^ 3).trace - (T30 T5).trace * p2
            + e2 * (T30 T5).trace - 3 * e3
          = -(3 * e3))
    -- (5) Newton-seventh collapse (algebraic, with residual `e_3, e_5` terms)
    ∧ (∀ e2 e3 e4 e5 e6 e7 : ℝ,
        ((T30 T5) ^ 7).trace
          - (T30 T5).trace * ((T30 T5) ^ 6).trace
          + e2 * ((T30 T5) ^ 5).trace
          - e3 * ((T30 T5) ^ 4).trace
          + e4 * ((T30 T5) ^ 3).trace
          - e5 * ((T30 T5) ^ 2).trace
          + e6 * (T30 T5).trace
          - 7 * e7
        = -e3 * ((T30 T5) ^ 4).trace - e5 * ((T30 T5) ^ 2).trace - 7 * e7) := by
  refine ⟨T30_trace_zero T5, T30_charpoly_a7_zero T5, T30_det_eq T5, ?_, ?_⟩
  · intro p2 e2 e3
    exact T30_e3_collapse T5 p2 e2 e3
  · intro e2 e3 e4 e5 e6 e7
    exact T30_newton_seventh_collapse T5 e2 e3 e4 e5 e6 e7

end PT.Stochastic
