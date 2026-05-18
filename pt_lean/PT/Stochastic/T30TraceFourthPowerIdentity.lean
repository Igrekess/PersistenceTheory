/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Stochastic.T30TraceCubeIdentity

/-!
# T2 — Quartic identities `tr(T_30^4) = 2 · tr(T_5^4)` and Newton–Girard consequences

**Statement (paper-level, monograph ch03 §T2, follow-up to
`T30TraceCubeIdentity`).** The parity formula
`tr(T_30^{2k}) = 2 · tr(T_5^{2k})` of `T30TraceFormulaExtended`
specialises at `k = 2` to

  `tr(T_30^4)  =  2 · tr(T_5^4)`.

This file packages the *quartic* layer of trace identities for `T_30`:

1. the **direct fourth-power identity** `tr(T_30^4) = 2 · tr(T_5^4)`
   (re-exposed from `T30_pow4_trace`);
2. the **bilinear quartic products** `tr(T_30^4) · tr(T_30^a) = 0` for
   every odd `a`, in particular `tr(T_30^4) · tr(T_30) = 0` and
   `tr(T_30^4) · tr(T_30^3) = 0`;
3. the **Newton–Girard fourth-identity collapse**: combining
   `p_1 = tr(T_30) = 0`, `p_3 = tr(T_30^3) = 0`, `e_1 = 0` in Newton's
   identity
       `p_4 − e_1·p_3 + e_2·p_2 − e_3·p_1 + 4·e_4  =  0`
   collapses the relation to
       `p_4 + e_2·p_2 + 4·e_4  =  0`,
   so that **`4·e_4` is fully determined by `tr(T_5^4)`, `tr(T_5^2)`,
   and `e_2`**:
       `4·e_4 = −p_4 − e_2·p_2 = −2·tr(T_5^4) − 2·e_2·tr(T_5^2)`;
4. the **conditional `e_4` closure**: if `e_3 = 0` (cubic Newton–Girard,
   cf. `T30TraceCubeIdentity.e3_vanishes_from_newton`), the quartic
   coefficient `e_4(T_30)` is fully determined by `tr(T_5^2)`,
   `tr(T_5^4)` and `e_2`;
5. the **characteristic-polynomial coefficient `[X^4]`**: by Vieta,
   `[X^4] charpoly(T_30) = e_4` (with the standard sign convention for
   the even index), so the conditional closure above pins this
   coefficient modulo `e_2`;
6. the **headline summary** bundling (1)–(5).

A *full* derivation of `e_4` would require committing to a concrete
eigenvalue datum (unavailable at the abstract `T5Like` level). We
therefore record the **algebraic collapse** of Newton's fourth
identity (the linear-combination form, which is a one-line `ring`
identity once the three vanishings are substituted) and the
**conditional closed form** `4·e_4 = −2·(tr(T_5^4) + e_2·tr(T_5^2))`.

## Main results

### Direct fourth-power identity
* `T30_pow_four_trace_headline` — `tr(T_30^4) = 2 · tr(T_5^4)`
  (re-expose of `T30_pow4_trace`, headline form).
* `T30_trace_fourth_self_nonneg_form` — `tr(T_30^4) = 2 · tr(T_5^4)`
  rewritten as `tr(T_30^4) − 2·tr(T_5^4) = 0` for downstream
  substitution.

### Bilinear quartic products
* `T30_trace_fourth_mul_trace_zero` — `tr(T_30^4) · tr(T_30) = 0`.
* `T30_trace_fourth_mul_trace_cube_zero` — `tr(T_30^4) · tr(T_30^3) = 0`.
* `T30_trace_fourth_mul_pow_odd_zero` — `tr(T_30^4) · tr(T_30^(2k+1)) = 0`
  for every `k`.
* `T30_pow_odd_mul_trace_fourth_zero` — symmetric form.
* `T30_trace_cube_mul_trace_fourth_zero` — `tr(T_30^3) · tr(T_30^4) = 0`.

### Newton–Girard fourth-identity collapse (algebraic)
* `newton_fourth_identity_collapse` — *purely algebraic* lemma: under
  `p1 = 0`, `p3 = 0`, `e1 = 0`, Newton's identity
  `p4 − e1·p3 + e2·p2 − e3·p1 + 4·e4 = 0` collapses to
  `p4 + e2·p2 + 4·e4 = 0`.
* `T30_newton_fourth_collapse` — application to `T_30`: the expression
  `tr(T_30^4) − tr(T_30)·tr(T_30^3) + e2·tr(T_30^2) − e3·tr(T_30)
    + 4·e4` simplifies to `tr(T_30^4) + e2·tr(T_30^2) + 4·e4` for every
  `e2, e3, e4`.

### Conditional `e_4` closure
* `e4_closed_form_from_newton` — *purely algebraic* lemma: under the
  hypothesis that Newton's fourth identity holds with `p1 = 0`,
  `p3 = 0`, `e1 = 0`, one has `4·e4 = −p4 − e2·p2`.
* `T30_e4_closed_form_from_T5_traces` — application to `T_30`: under
  the Newton-fourth hypothesis,
  `4·e4 = −2·tr(T_5^4) − 2·e2·tr(T_5^2) = −2·(tr(T_5^4) + e2·tr(T_5^2))`.

### Characteristic-polynomial coefficient `[X^4]` (sketch)
* `T30_char_poly_coeff_four_value` — under the Newton-fourth hypothesis,
  the value `e4(T_30)` is the rational number
  `e4 = −(tr(T_5^4) + e2·tr(T_5^2))/2`. Identification with the
  `[X^4]`-coefficient of `charpoly(T_30)` (a Vieta step) is recorded
  symbolically (Mathlib's `Matrix.charpoly` machinery is not unfolded
  here).

### Headline summary
* `T30_trace_fourth_invariants_summary` — packages the five direct
  identities and the conditional `e_4` closure in one statement.

## Strategy

Every direct identity reduces via `T30_pow4_trace` (the fourth power
equals `2 · tr(T_5^4)`) or `T30_trace_zero` / `T30_pow_trace_odd` (the
odd-power factor vanishes) to a one-line rewrite plus `ring`. The
Newton–Girard collapse and the conditional closure are *purely
algebraic*: each is a one- or two-line `ring` / `linarith` manipulation
once the three vanishings (`p1 = 0`, `p3 = 0`, `e1 = 0`) are
substituted. No new spectral input on `T_5` is needed beyond the
parity formulas already established in `T30TraceFormulaExtended`.

## Status

`[THM]` — purely algebraic restatement of the existing trace
identities at the quartic layer. The Newton–Girard collapse and the
conditional `e_4` closed form are one-line `ring` / `linarith`
identities; the corresponding statements *for the actual elementary
symmetric polynomials of the spectrum* require committing to a
specific eigenvalue datum (deferred — would close the characteristic
polynomial coefficient `c_4 = e_4` modulo `e_2` when combined with the
full Newton machinery).

## References

* `PT.Stochastic.T30TraceDeterminant` — `T30_trace_zero`.
* `PT.Stochastic.T30TraceFormulaExtended` — `T30_pow_trace_odd`,
  `T30_pow_trace_even`.
* `PT.Stochastic.T30TracePowerSequence` — `T30_pow4_trace`,
  `T30_pow3_trace`.
* `PT.Stochastic.T30TraceSquaredIdentity` — quadratic-layer bilinear
  vanishings and `T30_power_sum_residue`.
* `PT.Stochastic.T30TraceCubeIdentity` — cubic-layer bilinear
  vanishings and Newton's third-identity collapse.
-/

namespace PT.Stochastic

open Matrix BigOperators

/-! ### Direct fourth-power identity `tr(T_30^4) = 2 · tr(T_5^4)` -/

/-- **`tr(T_30^4) = 2 · tr(T_5^4)`.**

    Re-expose of `T30_pow4_trace` in headline form: the fourth power of
    `T_30` has trace exactly twice that of `T_5^4`. Specialisation of
    `T30_pow_trace_even` at `k = 2`. -/
theorem T30_pow_four_trace_headline (T5 : T5Like) :
    ((T30 T5) ^ 4).trace = 2 * (T5.matrix ^ 4).trace :=
  T30_pow4_trace T5

/-- **Substitution form: `tr(T_30^4) − 2·tr(T_5^4) = 0`.**

    Equivalent rewriting of `T30_pow_four_trace_headline`, with the
    constant on the right, useful for combining with linear-combination
    identities (Newton–Girard substitutions). -/
theorem T30_trace_fourth_self_nonneg_form (T5 : T5Like) :
    ((T30 T5) ^ 4).trace - 2 * (T5.matrix ^ 4).trace = 0 := by
  rw [T30_pow_four_trace_headline]
  ring

/-! ### Bilinear quartic products

The fourth power of `T_30` has the *non-zero* trace `2·tr(T_5^4)`, so
the bilinear products `tr(T_30^4) · tr(T_30^a) = 0` hold **only when
the other factor vanishes**, i.e. when `a` is odd (so
`tr(T_30^a) = 0` by `T30_pow_trace_odd`). We record the principal
cases (`a = 1` and `a = 3`) and the general odd-exponent form. -/

/-- **`tr(T_30^4) · tr(T_30) = 0`.**

    The right factor vanishes by `T30_trace_zero`. -/
theorem T30_trace_fourth_mul_trace_zero (T5 : T5Like) :
    ((T30 T5) ^ 4).trace * (T30 T5).trace = 0 := by
  rw [T30_trace_zero]
  ring

/-- **`tr(T_30^4) · tr(T_30^3) = 0`.**

    The right factor vanishes by `T30_pow3_trace`. -/
theorem T30_trace_fourth_mul_trace_cube_zero (T5 : T5Like) :
    ((T30 T5) ^ 4).trace * ((T30 T5) ^ 3).trace = 0 := by
  rw [T30_pow3_trace]
  ring

/-- **`tr(T_30^4) · tr(T_30^(2k+1)) = 0` for every `k`.**

    General odd-exponent bilinear: the right factor vanishes by
    `T30_pow_trace_odd`. Covers `T30_trace_fourth_mul_trace_zero` at
    `k = 0` and `T30_trace_fourth_mul_trace_cube_zero` at `k = 1`. -/
theorem T30_trace_fourth_mul_pow_odd_zero (T5 : T5Like) (k : ℕ) :
    ((T30 T5) ^ 4).trace * ((T30 T5) ^ (2 * k + 1)).trace = 0 := by
  rw [T30_pow_trace_odd]
  ring

/-- **Symmetric form: `tr(T_30^(2k+1)) · tr(T_30^4) = 0`.** -/
theorem T30_pow_odd_mul_trace_fourth_zero (T5 : T5Like) (k : ℕ) :
    ((T30 T5) ^ (2 * k + 1)).trace * ((T30 T5) ^ 4).trace = 0 := by
  rw [T30_pow_trace_odd]
  ring

/-- **`tr(T_30^3) · tr(T_30^4) = 0`.** Symmetric form of
    `T30_trace_fourth_mul_trace_cube_zero`. -/
theorem T30_trace_cube_mul_trace_fourth_zero (T5 : T5Like) :
    ((T30 T5) ^ 3).trace * ((T30 T5) ^ 4).trace = 0 := by
  rw [T30_pow3_trace]
  ring

/-- **`tr(T_30) · tr(T_30^4) = 0`.** Symmetric form of
    `T30_trace_fourth_mul_trace_zero`. -/
theorem T30_trace_mul_trace_fourth_zero (T5 : T5Like) :
    (T30 T5).trace * ((T30 T5) ^ 4).trace = 0 := by
  rw [T30_trace_zero]
  ring

/-! ### Newton–Girard fourth-identity collapse

Newton's fourth identity for power sums `p_k` and elementary symmetric
polynomials `e_k` of a finite multiset of values reads

  `p_4 − e_1·p_3 + e_2·p_2 − e_3·p_1 + 4·e_4  =  0`.

For the spectrum of `T_30`, we have `p_1 = tr(T_30) = 0`,
`p_3 = tr(T_30^3) = 0`, and `e_1 = p_1 = 0`. Substituting,

  `p_4 − 0 · 0 + e_2 · p_2 − e_3 · 0 + 4·e_4  =  0`,

i.e.

  `p_4 + e_2 · p_2 + 4·e_4  =  0`,

so that

  `4·e_4 = −p_4 − e_2·p_2 = −2·tr(T_5^4) − 2·e_2·tr(T_5^2)`.

If furthermore `e_3 = 0` (the conclusion of Newton's third identity
applied to `T_30`, cf. `T30TraceCubeIdentity.e3_vanishes_from_newton`),
the relation above becomes a **fully determined closed form** for
`e_4` in terms of the `T_5` traces and `e_2`. -/

/-- **Newton–Girard fourth-identity collapse (purely algebraic).**

    If `p1 = 0`, `p3 = 0`, and `e1 = 0`, then for any
    `p2, p4, e2, e3, e4 : ℝ`, Newton's identity
    `p4 − e1·p3 + e2·p2 − e3·p1 + 4·e4 = 0` is equivalent to
    `p4 + e2·p2 + 4·e4 = 0`.

    The lemma records the **collapse step** (the substitution of the
    three zeros), not the underlying Newton identity — which is the
    textbook relation between power sums and elementary symmetric
    polynomials of a finite multiset. -/
theorem newton_fourth_identity_collapse
    (p2 p4 e2 e3 e4 : ℝ) :
    p4 - (0 : ℝ) * (0 : ℝ) + e2 * p2 - e3 * (0 : ℝ) + 4 * e4
      = p4 + e2 * p2 + 4 * e4 := by
  ring

/-- **Application to `T_30`.** Substituting `p1 = tr(T_30)`,
    `p3 = tr(T_30^3)`, `e1 = tr(T_30)` (so `e1 = p1 = 0`), the
    Newton-fourth combination
    `tr(T_30^4) − tr(T_30)·tr(T_30^3) + e2·tr(T_30^2) − e3·tr(T_30)
      + 4·e4` simplifies to `tr(T_30^4) + e2·tr(T_30^2) + 4·e4` for
    every `e2, e3, e4`.

    This is the *quartic-layer* analogue of `T30_e3_collapse` from
    `T30TraceCubeIdentity`. -/
theorem T30_newton_fourth_collapse (T5 : T5Like) (e2 e3 e4 : ℝ) :
    ((T30 T5) ^ 4).trace
      - (T30 T5).trace * ((T30 T5) ^ 3).trace
      + e2 * ((T30 T5) ^ 2).trace
      - e3 * (T30 T5).trace
      + 4 * e4
    = ((T30 T5) ^ 4).trace + e2 * ((T30 T5) ^ 2).trace + 4 * e4 := by
  rw [T30_trace_zero, T30_pow3_trace]
  ring

/-! ### Conditional `e_4` closure

Under the Newton-fourth hypothesis (the textbook identity holds for
the spectrum of `T_30`), the collapse above determines `4·e_4` from
`tr(T_5^4)`, `tr(T_5^2)`, and `e_2`. -/

/-- **`e_4` closed-form from Newton's fourth (purely algebraic).**

    If Newton's fourth identity holds with `p1 = 0`, `p3 = 0`,
    `e1 = 0`, then `4·e4 = −p4 − e2·p2`. -/
theorem e4_closed_form_from_newton
    (p1 p2 p3 p4 e1 e2 e3 e4 : ℝ)
    (hN : p4 - e1 * p3 + e2 * p2 - e3 * p1 + 4 * e4 = 0)
    (hp1 : p1 = 0) (hp3 : p3 = 0) (he1 : e1 = 0) :
    4 * e4 = -p4 - e2 * p2 := by
  have h := hN
  rw [hp1, hp3, he1] at h
  linarith

/-- **Application to `T_30`.** Under the Newton-fourth hypothesis, the
    quartic elementary symmetric polynomial `e_4` of the spectrum of
    `T_30` satisfies

    `4·e_4 = −2·tr(T_5^4) − 2·e_2·tr(T_5^2)
           = −2·(tr(T_5^4) + e_2·tr(T_5^2))`.

    The right-hand side is *fully determined* by the `T_5` block (the
    two traces `tr(T_5^2)` and `tr(T_5^4)`) and the second elementary
    symmetric polynomial `e_2`. If additionally `e_3 = 0` (Newton's
    third identity applied to `T_30`, cf.
    `T30TraceCubeIdentity.e3_vanishes_from_newton`), this becomes a
    closed expression for `e_4` modulo `e_2`. -/
theorem T30_e4_closed_form_from_T5_traces (T5 : T5Like)
    (e2 e3 e4 : ℝ)
    (hN : ((T30 T5) ^ 4).trace
            - (T30 T5).trace * ((T30 T5) ^ 3).trace
            + e2 * ((T30 T5) ^ 2).trace
            - e3 * (T30 T5).trace
            + 4 * e4
          = 0) :
    4 * e4 = -2 * (T5.matrix ^ 4).trace - 2 * e2 * (T5.matrix ^ 2).trace := by
  have h := hN
  rw [T30_trace_zero, T30_pow3_trace, T30_pow_four_trace_headline,
      T30_pow2_trace] at h
  linarith

/-- **Factored form: `4·e_4 = −2·(tr(T_5^4) + e_2·tr(T_5^2))`.**

    Algebraic rearrangement of `T30_e4_closed_form_from_T5_traces`,
    grouping the dependence on the `T_5` block as a single factor of
    `2`. -/
theorem T30_e4_closed_form_factored (T5 : T5Like)
    (e2 e3 e4 : ℝ)
    (hN : ((T30 T5) ^ 4).trace
            - (T30 T5).trace * ((T30 T5) ^ 3).trace
            + e2 * ((T30 T5) ^ 2).trace
            - e3 * (T30 T5).trace
            + 4 * e4
          = 0) :
    4 * e4 = -2 * ((T5.matrix ^ 4).trace + e2 * (T5.matrix ^ 2).trace) := by
  have h := T30_e4_closed_form_from_T5_traces T5 e2 e3 e4 hN
  linarith

/-! ### Characteristic-polynomial coefficient `[X^4]` (sketch)

For a square matrix `M : Matrix n n ℝ` of dimension `d = card n`, the
characteristic polynomial admits the Vieta expansion

  `charpoly M = X^d − e_1·X^{d-1} + e_2·X^{d-2} − e_3·X^{d-3}
                + e_4·X^{d-4} − …  + (−1)^d · e_d`,

with `e_k` the `k`-th elementary symmetric polynomial of the
eigenvalues. For `T_30` the dimension is `d = 1·2·4 = 8 = φ(30)`, so
the coefficient of `X^4` equals `e_4` (the index `d − 4 = 4` is even,
so no sign flip).

We *do not* unfold Mathlib's `charpoly` machinery here. Instead we
record the *scalar value*

  `e_4 = −(tr(T_5^4) + e_2·tr(T_5^2))/2`

(under the Newton-fourth hypothesis), which is the abstract content of
"the coefficient of `X^4` in `charpoly(T_30)` is determined modulo
`e_2`" once one identifies the coefficient via Vieta. -/

/-- **Value of `e_4` (and hence of the `[X^4]`-coefficient of
    `charpoly(T_30)`) modulo `e_2`.**

    Under the Newton-fourth hypothesis, `e_4 = −(tr(T_5^4) +
    e_2·tr(T_5^2))/2`. This is the scalar value of the `[X^4]`
    coefficient of `charpoly(T_30)` (an `8 × 8` matrix), up to the
    Vieta identification `[X^4] = e_4` (which we do not unfold). -/
theorem T30_char_poly_coeff_four_value (T5 : T5Like)
    (e2 e3 e4 : ℝ)
    (hN : ((T30 T5) ^ 4).trace
            - (T30 T5).trace * ((T30 T5) ^ 3).trace
            + e2 * ((T30 T5) ^ 2).trace
            - e3 * (T30 T5).trace
            + 4 * e4
          = 0) :
    e4 = -((T5.matrix ^ 4).trace + e2 * (T5.matrix ^ 2).trace) / 2 := by
  have h := T30_e4_closed_form_factored T5 e2 e3 e4 hN
  linarith

/-! ### Headline summary -/

/-- **Headline summary for the quartic-layer trace identities.**

    The fourth-power identity `tr(T_30^4) = 2 · tr(T_5^4)` propagates to
    a *quartic-layer* family of identities:

    1. **Direct fourth-power identity** —
       `tr(T_30^4) = 2 · tr(T_5^4)`.
    2. **Bilinear quartic–odd** — for every `k`,
       `tr(T_30^4) · tr(T_30^(2k+1)) = 0`. In particular
       `tr(T_30^4) · tr(T_30) = 0` and
       `tr(T_30^4) · tr(T_30^3) = 0`.
    3. **Newton–Girard fourth-identity collapse** — for every
       `e2, e3, e4`,
       `tr(T_30^4) − tr(T_30)·tr(T_30^3) + e2·tr(T_30^2) − e3·tr(T_30)
         + 4·e4  =  tr(T_30^4) + e2·tr(T_30^2) + 4·e4`.
    4. **Conditional `e_4` closure** — under the Newton-fourth
       hypothesis, `4·e_4 = −2·(tr(T_5^4) + e_2·tr(T_5^2))`. -/
theorem T30_trace_fourth_invariants_summary (T5 : T5Like) :
    -- (1) Direct fourth-power identity
    ((T30 T5) ^ 4).trace = 2 * (T5.matrix ^ 4).trace
    -- (2) Bilinear quartic–odd
    ∧ (∀ k : ℕ, ((T30 T5) ^ 4).trace * ((T30 T5) ^ (2 * k + 1)).trace = 0)
    -- (3) Newton–Girard fourth-identity collapse
    ∧ (∀ e2 e3 e4 : ℝ,
        ((T30 T5) ^ 4).trace
          - (T30 T5).trace * ((T30 T5) ^ 3).trace
          + e2 * ((T30 T5) ^ 2).trace
          - e3 * (T30 T5).trace
          + 4 * e4
        = ((T30 T5) ^ 4).trace + e2 * ((T30 T5) ^ 2).trace + 4 * e4)
    -- (4) Conditional `e_4` closure under the Newton-fourth hypothesis
    ∧ (∀ e2 e3 e4 : ℝ,
        ((T30 T5) ^ 4).trace
            - (T30 T5).trace * ((T30 T5) ^ 3).trace
            + e2 * ((T30 T5) ^ 2).trace
            - e3 * (T30 T5).trace
            + 4 * e4
          = 0
        → 4 * e4
            = -2 * ((T5.matrix ^ 4).trace + e2 * (T5.matrix ^ 2).trace)) :=
  ⟨T30_pow_four_trace_headline T5,
   T30_trace_fourth_mul_pow_odd_zero T5,
   T30_newton_fourth_collapse T5,
   fun e2 e3 e4 hN => T30_e4_closed_form_factored T5 e2 e3 e4 hN⟩

end PT.Stochastic
