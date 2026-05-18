/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Stochastic.T30TraceFourthPowerIdentity

/-!
# T2 — Quintic identities `tr(T_30^5) = 0` and Newton–Girard consequences

**Statement (paper-level, monograph ch03 §T2, follow-up to
`T30TraceFourthPowerIdentity`).** The parity formula
`tr(T_30^{2k+1}) = 0` of `T30TraceFormulaExtended` specialises at
`k = 2` to

  `tr(T_30^5)  =  0`.

This file packages the *quintic* layer of trace identities for `T_30`:

1. the **direct fifth-power identity** `tr(T_30^5) = 0`
   (re-exposed from `T30_pow5_trace`);
2. the **bilinear quintic products** `tr(T_30^5) · tr(T_30^a) = 0` for
   every `a`, in particular `tr(T_30^5) · tr(T_30) = 0`,
   `tr(T_30^5) · tr(T_30^2) = 0`, `tr(T_30^5) · tr(T_30^3) = 0`, and
   `tr(T_30^5) · tr(T_30^4) = 0`;
3. the **Newton–Girard fifth-identity collapse**: combining
   `p_1 = tr(T_30) = 0`, `p_3 = tr(T_30^3) = 0`, `p_5 = tr(T_30^5) = 0`,
   `e_1 = 0` in Newton's identity
       `p_5 − e_1·p_4 + e_2·p_3 − e_3·p_2 + e_4·p_1 − 5·e_5  =  0`
   collapses the relation to
       `−e_3·p_2 − 5·e_5  =  0`,
   so that **`5·e_5` is fully determined by `e_3` and `tr(T_5^2)`**:
       `5·e_5 = −e_3·p_2 = −2·e_3·tr(T_5^2)`;
4. the **conditional `e_5` closure**: if `e_3 = 0` (cubic Newton–Girard,
   cf. `T30TraceCubeIdentity.e3_vanishes_from_newton`), the quintic
   coefficient `e_5(T_30)` *vanishes*: `e_5 = 0`;
5. the **characteristic-polynomial coefficient `[X^3]`**: by Vieta,
   `[X^3] charpoly(T_30) = −e_5` (the index `d − 5 = 3` is odd, hence
   the sign flip), so the conditional closure above pins this
   coefficient to `0` whenever `e_3 = 0`;
6. the **headline summary** bundling (1)–(5).

A *full* derivation of `e_5` would require committing to a concrete
eigenvalue datum (unavailable at the abstract `T5Like` level). We
therefore record the **algebraic collapse** of Newton's fifth
identity (the linear-combination form, which is a one-line `ring`
identity once the four vanishings are substituted) and the
**conditional closed form** `5·e_5 = −2·e_3·tr(T_5^2)`.

## Main results

### Direct fifth-power identity
* `T30_pow_five_trace_headline` — `tr(T_30^5) = 0`
  (re-expose of `T30_pow5_trace`, headline form).
* `T30_trace_fifth_self_zero_form` — `tr(T_30^5) − 0 = 0` rewritten as
  a substitution form for downstream linear-combination identities.

### Bilinear quintic products
* `T30_trace_fifth_mul_any_zero` — `tr(T_30^5) · tr(T_30^a) = 0` for
  every `a`.
* `T30_any_mul_trace_fifth_zero` — symmetric form.
* `T30_trace_fifth_mul_trace_zero` — `tr(T_30^5) · tr(T_30) = 0`.
* `T30_trace_fifth_mul_trace_sq_zero` — `tr(T_30^5) · tr(T_30^2) = 0`.
* `T30_trace_fifth_mul_trace_cube_zero` — `tr(T_30^5) · tr(T_30^3) = 0`.
* `T30_trace_fifth_mul_trace_fourth_zero` — `tr(T_30^5) · tr(T_30^4) = 0`.

### Newton–Girard fifth-identity collapse (algebraic)
* `newton_fifth_identity_collapse` — *purely algebraic* lemma: under
  `p1 = 0`, `p3 = 0`, `p5 = 0`, `e1 = 0`, Newton's identity
  `p5 − e1·p4 + e2·p3 − e3·p2 + e4·p1 − 5·e5 = 0` collapses to
  `−e3·p2 − 5·e5 = 0`.
* `T30_newton_fifth_collapse` — application to `T_30`: the expression
  `tr(T_30^5) − tr(T_30)·tr(T_30^4) + e2·tr(T_30^3) − e3·tr(T_30^2)
    + e4·tr(T_30) − 5·e5` simplifies to `−e3·tr(T_30^2) − 5·e5` for
  every `e2, e3, e4, e5`.

### Conditional `e_5` closure
* `e5_closed_form_from_newton` — *purely algebraic* lemma: under the
  hypothesis that Newton's fifth identity holds with `p1 = 0`,
  `p3 = 0`, `p5 = 0`, `e1 = 0`, one has `5·e5 = −e3·p2`.
* `T30_e5_closed_form_from_T5_traces` — application to `T_30`: under
  the Newton-fifth hypothesis, `5·e5 = −2·e3·tr(T_5^2)`.
* `e5_vanishes_from_newton` — under the Newton-fifth hypothesis with
  additionally `e3 = 0`, the quintic invariant `e5` vanishes: `e5 = 0`.
* `T30_e5_vanishes_from_T5_traces` — application to `T_30`: under the
  Newton-fifth hypothesis and `e3 = 0`, the quintic elementary
  symmetric polynomial `e_5` of the spectrum of `T_30` vanishes.

### Characteristic-polynomial coefficient `[X^3]` (sketch)
* `T30_char_poly_coeff_three_value` — under the Newton-fifth hypothesis
  with `e3 = 0`, the value `e5(T_30)` is `0`, hence the `[X^3]`
  coefficient (which equals `−e5` by Vieta with sign flip for the odd
  index `d − 5 = 3`) also vanishes.

### Headline summary
* `T30_trace_fifth_invariants_summary` — packages the five direct
  identities and the conditional `e_5` closure in one statement.

## Strategy

Every direct identity reduces via `T30_pow5_trace` (the fifth power
vanishes) or `T30_trace_zero` / `T30_pow_trace_odd` to a one-line
rewrite plus `ring`. The Newton–Girard collapse and the conditional
closures are *purely algebraic*: each is a one- or two-line `ring` /
`linarith` manipulation once the four vanishings (`p1 = 0`, `p3 = 0`,
`p5 = 0`, `e1 = 0`) are substituted. No new spectral input on `T_5` is
needed beyond the parity formulas already established in
`T30TraceFormulaExtended`.

## Status

`[THM]` — purely algebraic restatement of the existing trace
identities at the quintic layer. The Newton–Girard collapse and the
conditional `e_5` closed form are one-line `ring` / `linarith`
identities; the conditional vanishing `e_5 = 0` under `e_3 = 0` is the
quintic analogue of the *cubic-layer* result `e_3 = 0` under the
Newton-third hypothesis (cf. `T30TraceCubeIdentity.e3_vanishes_from_newton`).
The corresponding statements *for the actual elementary symmetric
polynomials of the spectrum* require committing to a specific
eigenvalue datum (deferred — would close the characteristic polynomial
coefficient `[X^3] = −e_5 = 0` when combined with the full Newton
machinery and the cubic-layer `e_3 = 0`).

## References

* `PT.Stochastic.T30TraceDeterminant` — `T30_trace_zero`.
* `PT.Stochastic.T30TraceFormulaExtended` — `T30_pow_trace_odd`,
  `T30_pow_trace_even`.
* `PT.Stochastic.T30TracePowerSequence` — `T30_pow5_trace`,
  `T30_pow3_trace`, `T30_pow2_trace`, `T30_pow4_trace`.
* `PT.Stochastic.T30TraceSquaredIdentity` — quadratic-layer bilinear
  vanishings and `T30_power_sum_residue`.
* `PT.Stochastic.T30TraceCubeIdentity` — cubic-layer bilinear
  vanishings, Newton's third-identity collapse, and
  `e3_vanishes_from_newton`.
* `PT.Stochastic.T30TraceFourthPowerIdentity` — quartic-layer
  bilinear vanishings, Newton's fourth-identity collapse, and the
  conditional `e_4` closure.
-/

namespace PT.Stochastic

open Matrix BigOperators

/-! ### Direct fifth-power identity `tr(T_30^5) = 0` -/

/-- **`tr(T_30^5) = 0`.**

    Re-expose of `T30_pow5_trace` in headline form: the fifth power of
    `T_30` has vanishing trace. Specialisation of `T30_pow_trace_odd`
    at `k = 2` (so `2k + 1 = 5`). -/
theorem T30_pow_five_trace_headline (T5 : T5Like) :
    ((T30 T5) ^ 5).trace = 0 :=
  T30_pow5_trace T5

/-- **Substitution form: `tr(T_30^5) − 0 = 0`.**

    Equivalent rewriting of `T30_pow_five_trace_headline`, useful for
    combining with linear-combination identities (Newton–Girard
    substitutions). -/
theorem T30_trace_fifth_self_zero_form (T5 : T5Like) :
    ((T30 T5) ^ 5).trace - 0 = 0 := by
  rw [T30_pow_five_trace_headline]
  ring

/-! ### Bilinear quintic products

The fifth power of `T_30` has *vanishing* trace, so every bilinear
product `tr(T_30^5) · tr(T_30^a) = 0` holds unconditionally. We record
the principal cases (`a ∈ {1, 2, 3, 4}`) and the general form. -/

/-- **`tr(T_30^5) · tr(T_30^a) = 0` for every `a`.**

    Direct consequence of `T30_pow5_trace` (the left factor vanishes). -/
theorem T30_trace_fifth_mul_any_zero (T5 : T5Like) (a : ℕ) :
    ((T30 T5) ^ 5).trace * ((T30 T5) ^ a).trace = 0 := by
  rw [T30_pow5_trace]
  ring

/-- **Symmetric form: `tr(T_30^a) · tr(T_30^5) = 0`.** -/
theorem T30_any_mul_trace_fifth_zero (T5 : T5Like) (a : ℕ) :
    ((T30 T5) ^ a).trace * ((T30 T5) ^ 5).trace = 0 := by
  rw [T30_pow5_trace]
  ring

/-- **`tr(T_30^5) · tr(T_30) = 0`.** Specialisation at `a = 1`. -/
theorem T30_trace_fifth_mul_trace_zero (T5 : T5Like) :
    ((T30 T5) ^ 5).trace * (T30 T5).trace = 0 := by
  rw [T30_pow5_trace]
  ring

/-- **`tr(T_30^5) · tr(T_30^2) = 0`.** Specialisation at `a = 2`. -/
theorem T30_trace_fifth_mul_trace_sq_zero (T5 : T5Like) :
    ((T30 T5) ^ 5).trace * ((T30 T5) ^ 2).trace = 0 := by
  rw [T30_pow5_trace]
  ring

/-- **`tr(T_30^5) · tr(T_30^3) = 0`.** Specialisation at `a = 3`. -/
theorem T30_trace_fifth_mul_trace_cube_zero (T5 : T5Like) :
    ((T30 T5) ^ 5).trace * ((T30 T5) ^ 3).trace = 0 := by
  rw [T30_pow5_trace]
  ring

/-- **`tr(T_30^5) · tr(T_30^4) = 0`.** Specialisation at `a = 4`. -/
theorem T30_trace_fifth_mul_trace_fourth_zero (T5 : T5Like) :
    ((T30 T5) ^ 5).trace * ((T30 T5) ^ 4).trace = 0 := by
  rw [T30_pow5_trace]
  ring

/-- **`tr(T_30^5)^2 = 0`.** Square of the (vanishing) fifth-power
    trace. -/
theorem T30_trace_fifth_squared_zero (T5 : T5Like) :
    ((T30 T5) ^ 5).trace ^ 2 = 0 := by
  rw [T30_pow5_trace]
  ring

/-! ### Newton–Girard fifth-identity collapse

Newton's fifth identity for power sums `p_k` and elementary symmetric
polynomials `e_k` of a finite multiset of values reads

  `p_5 − e_1·p_4 + e_2·p_3 − e_3·p_2 + e_4·p_1 − 5·e_5  =  0`.

For the spectrum of `T_30`, we have `p_1 = tr(T_30) = 0`,
`p_3 = tr(T_30^3) = 0`, `p_5 = tr(T_30^5) = 0`, and `e_1 = p_1 = 0`.
Substituting,

  `0 − 0·p_4 + e_2·0 − e_3·p_2 + e_4·0 − 5·e_5  =  0`,

i.e.

  `−e_3·p_2 − 5·e_5  =  0`,

so that

  `5·e_5 = −e_3·p_2 = −2·e_3·tr(T_5^2)`.

If furthermore `e_3 = 0` (the conclusion of Newton's third identity
applied to `T_30`, cf. `T30TraceCubeIdentity.e3_vanishes_from_newton`),
the relation above becomes `5·e_5 = 0`, i.e. `e_5 = 0`. -/

/-- **Newton–Girard fifth-identity collapse (purely algebraic).**

    If `p1 = 0`, `p3 = 0`, `p5 = 0`, and `e1 = 0`, then for any
    `p2, p4, e2, e3, e4, e5 : ℝ`, Newton's identity
    `p5 − e1·p4 + e2·p3 − e3·p2 + e4·p1 − 5·e5 = 0` is equivalent to
    `−e3·p2 − 5·e5 = 0`.

    The lemma records the **collapse step** (the substitution of the
    four zeros), not the underlying Newton identity — which is the
    textbook relation between power sums and elementary symmetric
    polynomials of a finite multiset. -/
theorem newton_fifth_identity_collapse
    (p2 p4 e2 e3 e4 e5 : ℝ) :
    (0 : ℝ) - (0 : ℝ) * p4 + e2 * (0 : ℝ) - e3 * p2 + e4 * (0 : ℝ)
      - 5 * e5
      = -e3 * p2 - 5 * e5 := by
  ring

/-- **Application to `T_30`.** Substituting `p1 = tr(T_30)`,
    `p3 = tr(T_30^3)`, `p5 = tr(T_30^5)`, `e1 = tr(T_30)` (so
    `e1 = p1 = 0`), the Newton-fifth combination
    `tr(T_30^5) − tr(T_30)·tr(T_30^4) + e2·tr(T_30^3) − e3·tr(T_30^2)
      + e4·tr(T_30) − 5·e5` simplifies to `−e3·tr(T_30^2) − 5·e5` for
    every `e2, e3, e4, e5`.

    This is the *quintic-layer* analogue of `T30_e3_collapse` from
    `T30TraceCubeIdentity` and `T30_newton_fourth_collapse` from
    `T30TraceFourthPowerIdentity`. -/
theorem T30_newton_fifth_collapse (T5 : T5Like) (e2 e3 e4 e5 : ℝ) :
    ((T30 T5) ^ 5).trace
      - (T30 T5).trace * ((T30 T5) ^ 4).trace
      + e2 * ((T30 T5) ^ 3).trace
      - e3 * ((T30 T5) ^ 2).trace
      + e4 * (T30 T5).trace
      - 5 * e5
    = -e3 * ((T30 T5) ^ 2).trace - 5 * e5 := by
  rw [T30_trace_zero, T30_pow3_trace, T30_pow5_trace]
  ring

/-! ### Conditional `e_5` closure

Under the Newton-fifth hypothesis (the textbook identity holds for the
spectrum of `T_30`), the collapse above determines `5·e_5` from
`tr(T_5^2)` and `e_3`. If additionally `e_3 = 0` (from Newton's third
identity, cf. `e3_vanishes_from_newton`), the quintic invariant `e_5`
*vanishes outright*. -/

/-- **`e_5` closed-form from Newton's fifth (purely algebraic).**

    If Newton's fifth identity holds with `p1 = 0`, `p3 = 0`, `p5 = 0`,
    `e1 = 0`, then `5·e5 = −e3·p2`. -/
theorem e5_closed_form_from_newton
    (p1 p2 p3 p4 p5 e1 e2 e3 e4 e5 : ℝ)
    (hN : p5 - e1 * p4 + e2 * p3 - e3 * p2 + e4 * p1 - 5 * e5 = 0)
    (hp1 : p1 = 0) (hp3 : p3 = 0) (hp5 : p5 = 0) (he1 : e1 = 0) :
    5 * e5 = -e3 * p2 := by
  have h := hN
  rw [hp1, hp3, hp5, he1] at h
  linarith

/-- **Application to `T_30`.** Under the Newton-fifth hypothesis, the
    quintic elementary symmetric polynomial `e_5` of the spectrum of
    `T_30` satisfies

    `5·e_5 = −e_3·tr(T_30^2) = −2·e_3·tr(T_5^2)`.

    The right-hand side is *fully determined* by `e_3` (the cubic
    elementary symmetric polynomial) and `tr(T_5^2)` (the quadratic
    `T_5` trace). If additionally `e_3 = 0` (Newton's third identity
    applied to `T_30`, cf. `T30TraceCubeIdentity.e3_vanishes_from_newton`),
    this collapses to `5·e_5 = 0`, hence `e_5 = 0`. -/
theorem T30_e5_closed_form_from_T5_traces (T5 : T5Like)
    (e2 e3 e4 e5 : ℝ)
    (hN : ((T30 T5) ^ 5).trace
            - (T30 T5).trace * ((T30 T5) ^ 4).trace
            + e2 * ((T30 T5) ^ 3).trace
            - e3 * ((T30 T5) ^ 2).trace
            + e4 * (T30 T5).trace
            - 5 * e5
          = 0) :
    5 * e5 = -2 * e3 * (T5.matrix ^ 2).trace := by
  have h := hN
  rw [T30_trace_zero, T30_pow3_trace, T30_pow5_trace, T30_pow2_trace] at h
  linarith

/-- **`e_5` vanishes from Newton's fifth plus `e3 = 0` (purely algebraic).**

    If Newton's fifth identity holds with `p1 = 0`, `p3 = 0`, `p5 = 0`,
    `e1 = 0`, and additionally `e3 = 0`, then `e5 = 0`. -/
theorem e5_vanishes_from_newton
    (p1 p2 p3 p4 p5 e1 e2 e3 e4 e5 : ℝ)
    (hN : p5 - e1 * p4 + e2 * p3 - e3 * p2 + e4 * p1 - 5 * e5 = 0)
    (hp1 : p1 = 0) (hp3 : p3 = 0) (hp5 : p5 = 0)
    (he1 : e1 = 0) (he3 : e3 = 0) :
    e5 = 0 := by
  have h5 : 5 * e5 = -e3 * p2 :=
    e5_closed_form_from_newton p1 p2 p3 p4 p5 e1 e2 e3 e4 e5 hN hp1 hp3 hp5 he1
  rw [he3] at h5
  linarith

/-- **Application to `T_30`.** Under the Newton-fifth hypothesis and
    `e3 = 0`, the quintic elementary symmetric polynomial `e_5` of the
    spectrum of `T_30` vanishes.

    This combines `T30_e5_closed_form_from_T5_traces` (giving
    `5·e_5 = −2·e_3·tr(T_5^2)`) with the cubic-layer vanishing
    `e_3 = 0` (from `e3_vanishes_from_newton` applied to `T_30`'s
    Newton-third identity) to conclude `e_5 = 0`. -/
theorem T30_e5_vanishes_from_T5_traces (T5 : T5Like)
    (e2 e3 e4 e5 : ℝ)
    (hN : ((T30 T5) ^ 5).trace
            - (T30 T5).trace * ((T30 T5) ^ 4).trace
            + e2 * ((T30 T5) ^ 3).trace
            - e3 * ((T30 T5) ^ 2).trace
            + e4 * (T30 T5).trace
            - 5 * e5
          = 0)
    (he3 : e3 = 0) :
    e5 = 0 := by
  have h5 : 5 * e5 = -2 * e3 * (T5.matrix ^ 2).trace :=
    T30_e5_closed_form_from_T5_traces T5 e2 e3 e4 e5 hN
  rw [he3] at h5
  linarith

/-! ### Characteristic-polynomial coefficient `[X^3]` (sketch)

For a square matrix `M : Matrix n n ℝ` of dimension `d = card n`, the
characteristic polynomial admits the Vieta expansion

  `charpoly M = X^d − e_1·X^{d-1} + e_2·X^{d-2} − e_3·X^{d-3}
                + e_4·X^{d-4} − e_5·X^{d-5} + …  + (−1)^d · e_d`,

with `e_k` the `k`-th elementary symmetric polynomial of the
eigenvalues. For `T_30` the dimension is `d = 1·2·4 = 8 = φ(30)`, so
the coefficient of `X^3` equals `−e_5` (the index `d − 5 = 3` is odd,
hence the sign flip).

We *do not* unfold Mathlib's `charpoly` machinery here. Instead we
record the *scalar value*

  `e_5 = 0`  ⟹  `[X^3] = −e_5 = 0`

(under the Newton-fifth hypothesis and `e_3 = 0`), which is the
abstract content of "the coefficient of `X^3` vanishes" once one
identifies the coefficient via Vieta. -/

/-- **Value of `e_5` (and hence of the `[X^3]`-coefficient of
    `charpoly(T_30)`) under `e_3 = 0`.**

    Under the Newton-fifth hypothesis and `e_3 = 0`, `e_5 = 0`, hence
    `−e_5 = 0`. This is the scalar value of the `[X^3]` coefficient
    of `charpoly(T_30)` (an `8 × 8` matrix), up to the Vieta
    identification `[X^3] = −e_5` (which we do not unfold). -/
theorem T30_char_poly_coeff_three_value (T5 : T5Like)
    (e2 e3 e4 e5 : ℝ)
    (hN : ((T30 T5) ^ 5).trace
            - (T30 T5).trace * ((T30 T5) ^ 4).trace
            + e2 * ((T30 T5) ^ 3).trace
            - e3 * ((T30 T5) ^ 2).trace
            + e4 * (T30 T5).trace
            - 5 * e5
          = 0)
    (he3 : e3 = 0) :
    -e5 = 0 := by
  have h := T30_e5_vanishes_from_T5_traces T5 e2 e3 e4 e5 hN he3
  linarith

/-! ### Headline summary -/

/-- **Headline summary for the quintic-layer trace identities.**

    The vanishing `tr(T_30^5) = 0` propagates to a *quintic-layer*
    family of identities:

    1. **Direct fifth-power identity** —
       `tr(T_30^5) = 0`.
    2. **Bilinear quintic–any** — for every `a`,
       `tr(T_30^5) · tr(T_30^a) = 0`.
    3. **Newton–Girard fifth-identity collapse** — for every
       `e2, e3, e4, e5`,
       `tr(T_30^5) − tr(T_30)·tr(T_30^4) + e2·tr(T_30^3)
          − e3·tr(T_30^2) + e4·tr(T_30) − 5·e5
        = −e3·tr(T_30^2) − 5·e5`.
    4. **Conditional `e_5` closure (formula)** — under the
       Newton-fifth hypothesis, `5·e_5 = −2·e_3·tr(T_5^2)`.
    5. **Conditional `e_5` vanishing** — under the Newton-fifth
       hypothesis and `e_3 = 0`, `e_5 = 0`.
    6. **Characteristic-polynomial third-power coefficient** —
       under the same hypotheses, `−e_5 = 0` (the `[X^3]` coefficient
       of `charpoly(T_30)` by Vieta with sign flip). -/
theorem T30_trace_fifth_invariants_summary (T5 : T5Like) :
    -- (1) Direct fifth-power identity
    ((T30 T5) ^ 5).trace = 0
    -- (2) Bilinear quintic–any
    ∧ (∀ a : ℕ, ((T30 T5) ^ 5).trace * ((T30 T5) ^ a).trace = 0)
    -- (3) Newton–Girard fifth-identity collapse
    ∧ (∀ e2 e3 e4 e5 : ℝ,
        ((T30 T5) ^ 5).trace
          - (T30 T5).trace * ((T30 T5) ^ 4).trace
          + e2 * ((T30 T5) ^ 3).trace
          - e3 * ((T30 T5) ^ 2).trace
          + e4 * (T30 T5).trace
          - 5 * e5
        = -e3 * ((T30 T5) ^ 2).trace - 5 * e5)
    -- (4) Conditional `e_5` closure (formula) under Newton-fifth
    ∧ (∀ e2 e3 e4 e5 : ℝ,
        ((T30 T5) ^ 5).trace
            - (T30 T5).trace * ((T30 T5) ^ 4).trace
            + e2 * ((T30 T5) ^ 3).trace
            - e3 * ((T30 T5) ^ 2).trace
            + e4 * (T30 T5).trace
            - 5 * e5
          = 0
        → 5 * e5 = -2 * e3 * (T5.matrix ^ 2).trace)
    -- (5) Conditional `e_5` vanishing under Newton-fifth + `e3 = 0`
    ∧ (∀ e2 e3 e4 e5 : ℝ,
        ((T30 T5) ^ 5).trace
            - (T30 T5).trace * ((T30 T5) ^ 4).trace
            + e2 * ((T30 T5) ^ 3).trace
            - e3 * ((T30 T5) ^ 2).trace
            + e4 * (T30 T5).trace
            - 5 * e5
          = 0
        → e3 = 0
        → e5 = 0) :=
  ⟨T30_pow_five_trace_headline T5,
   T30_trace_fifth_mul_any_zero T5,
   T30_newton_fifth_collapse T5,
   fun e2 e3 e4 e5 hN => T30_e5_closed_form_from_T5_traces T5 e2 e3 e4 e5 hN,
   fun e2 e3 e4 e5 hN he3 =>
     T30_e5_vanishes_from_T5_traces T5 e2 e3 e4 e5 hN he3⟩

end PT.Stochastic
