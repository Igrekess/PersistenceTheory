/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Stochastic.T30TraceSquaredIdentity

/-!
# T2 — Cubic identities `tr(T_30^3) = 0` and Newton–Girard consequences

**Statement (paper-level, monograph ch03 §T2, follow-up to
`T30TraceSquaredIdentity`).** The parity formula
`tr(T_30^{2k+1}) = 0` of `T30TraceFormulaExtended` specialises at
`k = 1` to

  `tr(T_30^3)  =  0`.

This file packages the *cubic* layer of trace identities for `T_30`:

1. the **direct cube identity** `tr(T_30^3) = 0`;
2. the **bilinear / trilinear cube products** `tr(T_30^a) · tr(T_30^3) = 0`
   for arbitrary `a`, and the symmetric / iterated variants;
3. the **Newton–Girard `e_1 = 0` consequence**: the coefficient of
   `λ^{φ(30) - 1} = λ^7` in the characteristic polynomial of `T_30`
   equals `−tr(T_30) = 0`. We expose this as the scalar identity
   `−(tr T_30) = 0`, which is the *seventh-power coefficient* of the
   characteristic polynomial of any matrix (by Vieta / Newton);
4. the **Newton–Girard residue at the cubic level**: combining
   `p_1 = tr(T_30) = 0` and `p_3 = tr(T_30^3) = 0` in Newton's identity
   `p_3 − e_1·p_2 + e_2·p_1 − 3·e_3 = 0` collapses the relation to
   `3·e_3 = 0`, i.e. **the elementary symmetric `e_3` of the
   eigenvalues of `T_30` vanishes**. We expose the *algebraic form* of
   this identity (the linear-combination form `p_3 − e_1·p_2 + e_2·p_1
   − 3·e_3 = 0` with `e_1 = 0`, `p_1 = 0`, `p_3 = 0` ⟹ `3·e_3 = 0`)
   as a structural lemma: the closed-form expression `3·e_3 =
   −e_1·p_2 + e_2·p_1 + p_3 = 0` is recorded as a scalar equation
   parametric in `e_1, e_2, p_2` (those quantities being indeterminate
   here at the abstract `T5Like` level — the *implication* `e_1 = 0
   ∧ p_3 = 0 → 3·e_3 = 0` modulo Newton is what we record);
5. the **headline summary** bundling (1)–(4).

A full formalisation of Newton's identities for the spectrum of `T_30`
would require committing to a specific eigenvalue datum (which is not
available at the `T5Like` abstract level, where `T_5` is left
parametric). We therefore record the *algebraic identities* on the
traces and the **conditional implication** for `e_3` (under the
algebraic form of Newton's third identity).

## Main results

### Direct cube identity
* `T30_pow_three_trace_zero` — `tr(T_30^3) = 0`
  (re-expose of `T30_pow3_trace`, in headline form).
* `T30_trace_cube_self_zero` — `(tr T_30)^3 = 0` (cube of the trace,
  not trace of the cube — both vanish).

### Bilinear / trilinear cube products
* `T30_trace_cube_mul_any_zero` — `tr(T_30^3) · tr(T_30^a) = 0` for
  every `a`.
* `T30_any_mul_trace_cube_zero` — symmetric form.
* `T30_trace_cube_mul_trace_zero` — `tr(T_30^3) · tr(T_30) = 0`.
* `T30_trace_cube_mul_trace_sq_zero` — `tr(T_30^3) · tr(T_30^2) = 0`.
* `T30_triple_cube_product_zero` — `tr(T_30^3)^3 = 0`.

### Newton–Girard at the cubic level (algebraic implication)
* `newton_third_identity_collapse` — *purely algebraic* lemma:
  if `p1 = 0`, `p3 = 0`, `e1 = 0`, then Newton's identity
  `p3 − e1·p2 + e2·p1 − 3·e3 = 0` collapses to `3·e3 = 0`.
* `T30_e3_collapse` — application to `T_30`: with `p1 = tr(T_30)`,
  `p3 = tr(T_30^3)` and `e1 := tr(T_30)`, the linear combination
  `tr(T_30^3) − tr(T_30)·p2 + e2·tr(T_30) − 3·e3` simplifies to
  `−3·e3` for every choice of `p2, e2, e3`.

### Characteristic polynomial coefficient (sketch)
* `T30_char_poly_coeff_seven_zero` — the *seventh-power* coefficient
  of the characteristic polynomial of `T_30` equals `−(tr T_30)`, hence
  equals `0`. We record the value `−(tr T_30) = 0` directly; the
  identification with `c_7` is a textbook Vieta step (Mathlib provides
  `Matrix.charpoly_natDegree`, `Matrix.charpoly_coeff_neg_trace`).

### Headline summary
* `T30_trace_cube_invariants_summary` — packages the four direct
  identities of (1)–(2) plus the Newton–Girard collapse of (3) in one
  statement.

## Strategy

Every direct identity reduces via `T30_pow3_trace` (the cube vanishes)
or `T30_trace_zero` to `_ * 0 = 0` / `0 * _ = 0` / `0^k = 0`. The
Newton–Girard collapse is *purely algebraic*: a one-line `ring`
manipulation. The implication for `e_3` then follows by substitution.

## Status

`[THM]` — purely algebraic restatement of the existing trace vanishing
results at the cubic layer. The Newton–Girard collapse is a one-line
`ring` identity; the corresponding statement *for the actual
elementary symmetric polynomials of the spectrum* requires committing
to a specific eigenvalue datum (deferred — would close the
characteristic polynomial coefficient `c_5 = −e_3 = 0` when combined
with the full Newton machinery).

## References

* `PT.Stochastic.T30TraceDeterminant` — `T30_trace_zero`.
* `PT.Stochastic.T30TraceFormulaExtended` — `T30_pow_trace_odd`,
  `T30_pow_trace_even`.
* `PT.Stochastic.T30TracePowerSequence` — `T30_pow3_trace`.
* `PT.Stochastic.T30TraceSquaredIdentity` — squared-level bilinear
  vanishings.
-/

namespace PT.Stochastic

open Matrix BigOperators

/-! ### Direct cube identity `tr(T_30^3) = 0` -/

/-- **`tr(T_30^3) = 0`.**

    Re-expose of `T30_pow3_trace` in headline form: the cube of `T_30`
    has vanishing trace. Specialisation of `T30_pow_trace_odd` at
    `k = 1` (so `2k + 1 = 3`). -/
theorem T30_pow_three_trace_zero (T5 : T5Like) :
    ((T30 T5) ^ 3).trace = 0 :=
  T30_pow3_trace T5

/-- **`(tr T_30)^3 = 0`.**

    The *cube of the trace* (as opposed to the *trace of the cube*) is
    also zero, since `tr(T_30) = 0`. Both quantities vanish, but for
    distinct reasons: this one is `0^3 = 0`; the previous one is the
    parity formula at `k = 1`. -/
theorem T30_trace_cube_self_zero (T5 : T5Like) :
    ((T30 T5).trace) ^ 3 = 0 := by
  rw [T30_trace_zero]
  ring

/-! ### Bilinear / trilinear cube products -/

/-- **`tr(T_30^3) · tr(T_30^a) = 0` for every `a`.**

    Direct consequence of `T30_pow3_trace` (the left factor vanishes).
    Special case of `T30_pow_trace_odd_mul_any` at `k = 1`. -/
theorem T30_trace_cube_mul_any_zero (T5 : T5Like) (a : ℕ) :
    ((T30 T5) ^ 3).trace * ((T30 T5) ^ a).trace = 0 := by
  rw [T30_pow3_trace]
  ring

/-- **Symmetric form: `tr(T_30^a) · tr(T_30^3) = 0`.** -/
theorem T30_any_mul_trace_cube_zero (T5 : T5Like) (a : ℕ) :
    ((T30 T5) ^ a).trace * ((T30 T5) ^ 3).trace = 0 := by
  rw [T30_pow3_trace]
  ring

/-- **`tr(T_30^3) · tr(T_30) = 0`.** Specialisation at `a = 1`. -/
theorem T30_trace_cube_mul_trace_zero (T5 : T5Like) :
    ((T30 T5) ^ 3).trace * (T30 T5).trace = 0 := by
  rw [T30_pow3_trace]
  ring

/-- **`tr(T_30^3) · tr(T_30^2) = 0`.** Specialisation at `a = 2`. -/
theorem T30_trace_cube_mul_trace_sq_zero (T5 : T5Like) :
    ((T30 T5) ^ 3).trace * ((T30 T5) ^ 2).trace = 0 := by
  rw [T30_pow3_trace]
  ring

/-- **`tr(T_30^3)^3 = 0`.** The cube of `tr(T_30^3)` (a `0^3`). -/
theorem T30_triple_cube_product_zero (T5 : T5Like) :
    ((T30 T5) ^ 3).trace ^ 3 = 0 := by
  rw [T30_pow3_trace]
  ring

/-- **Triple cube-cube-any product vanishes.** Three-factor identity:
    any product involving two copies of `tr(T_30^3)` (or one copy of
    `tr(T_30^3)` and one copy of `tr(T_30)`) vanishes. -/
theorem T30_cube_cube_any_zero (T5 : T5Like) (a : ℕ) :
    ((T30 T5) ^ 3).trace * ((T30 T5) ^ 3).trace * ((T30 T5) ^ a).trace = 0 := by
  rw [T30_pow3_trace]
  ring

/-! ### Newton–Girard collapse at the cubic level

Newton's third identity for the power sums `p_k = Σ λ_i^k` and the
elementary symmetric polynomials `e_k = Σ_{i_1 < … < i_k} λ_{i_1}·…·λ_{i_k}`
of a multiset of values `(λ_i)` reads

  `p_3 − e_1·p_2 + e_2·p_1 − 3·e_3  =  0`.

For the spectrum of `T_30`, we have `p_1 = tr(T_30) = 0`,
`p_3 = tr(T_30^3) = 0`, and `e_1 = p_1 = 0`. Substituting,

  `0 − 0·p_2 + e_2·0 − 3·e_3  =  0`     ⟹     `3·e_3 = 0`     ⟹     `e_3 = 0`.

We expose this collapse *algebraically* (without committing to a
specific eigenvalue datum, which is unavailable at the abstract
`T5Like` level): the linear combination simplifies for *any* choice of
`p_2`, `e_2`, `e_3` once `p_1 = 0`, `p_3 = 0`, `e_1 = 0`. -/

/-- **Newton–Girard third-identity collapse (purely algebraic).**

    If `p1 = 0`, `p3 = 0`, and `e1 = 0`, then for any `p2, e2, e3 : ℝ`,
    Newton's identity `p3 − e1·p2 + e2·p1 − 3·e3 = 0` collapses to
    `−3·e3 = 0`, hence `e3 = 0`.

    The lemma records the **collapse step** (the substitution of the
    three zeros), not the underlying Newton identity itself — which is
    the textbook relation between power sums and elementary symmetric
    polynomials of a finite multiset. -/
theorem newton_third_identity_collapse
    (p2 e2 e3 : ℝ) :
    (0 : ℝ) - (0 : ℝ) * p2 + e2 * (0 : ℝ) - 3 * e3 = -(3 * e3) := by
  ring

/-- **Newton–Girard third-identity collapse, `e3 = 0` form.**

    If we additionally assume Newton's identity itself
    (`p3 − e1·p2 + e2·p1 − 3·e3 = 0`) and the three vanishings
    `p1 = 0`, `p3 = 0`, `e1 = 0`, then `e3 = 0`.

    This packages the *implication* that, in the spectrum of `T_30`,
    the elementary symmetric polynomial `e_3 = Σ_{i<j<k} λ_i·λ_j·λ_k`
    vanishes — modulo Newton's third identity (which is a textbook
    relation we leave as the algebraic hypothesis `hN`). -/
theorem e3_vanishes_from_newton
    (p1 p2 p3 e1 e2 e3 : ℝ)
    (hN : p3 - e1 * p2 + e2 * p1 - 3 * e3 = 0)
    (hp1 : p1 = 0) (hp3 : p3 = 0) (he1 : e1 = 0) :
    e3 = 0 := by
  have : -(3 * e3) = 0 := by
    have := hN
    rw [hp1, hp3, he1] at this
    linarith
  linarith

/-- **Application to `T_30`.** Substituting `p1 = tr(T_30)`,
    `p3 = tr(T_30^3)`, `e1 = tr(T_30)` (so `e1 = p1 = 0`), the
    expression `tr(T_30^3) − tr(T_30)·p2 + e2·tr(T_30) − 3·e3`
    simplifies to `−3·e3` for every `p2, e2, e3`.

    This is the *cubic-layer* analogue of `T30_power_sum_residue` from
    `T30TraceSquaredIdentity`: at the quadratic layer the residue
    `p2 − p1^2` was `2 · tr(T_5^2)`; at the cubic layer the residue
    `p3 − 3·e1·e2 + 3·e3` (Newton's third identity, rewritten using
    `p1 = 0`) collapses to `−3·e3`, isolating `e3` as the only
    non-trivial third-order spectral invariant of `T_30`. -/
theorem T30_e3_collapse (T5 : T5Like) (p2 e2 e3 : ℝ) :
    ((T30 T5) ^ 3).trace - (T30 T5).trace * p2 + e2 * (T30 T5).trace
      - 3 * e3 = -(3 * e3) := by
  rw [T30_pow3_trace, T30_trace_zero]
  ring

/-! ### Characteristic-polynomial coefficient interpretation

For a square matrix `M : Matrix n n ℝ` of dimension `d = card n`, the
characteristic polynomial `charpoly M = det (X • 1 − M)` has the form

  `X^d − (tr M) · X^{d-1} + … + (−1)^d · det M`.

In particular, the coefficient of `X^{d-1}` equals `−tr M`. For `T_30`
the dimension is `d = 1 · 2 · 4 = 8 = φ(30)`, so the coefficient of
`X^7` equals `−tr(T_30) = 0`.

We *do not* unfold Mathlib's `charpoly` machinery here (which would
require importing `Mathlib.LinearAlgebra.Matrix.Charpoly.Basic` and
casting from `ℝ` to a polynomial ring). Instead we record the *scalar
value* `−(tr T_30) = 0`, which is the abstract content of "the
coefficient of `X^7` vanishes" once one identifies the coefficient via
Vieta. -/

/-- **Coefficient-of-`X^7` value.** The scalar `−(tr T_30) = 0`. By
    the Vieta formula for the characteristic polynomial of an `8 × 8`
    matrix, this is the coefficient of `X^7` in `charpoly (T_30)`. -/
theorem T30_char_poly_coeff_seven_zero (T5 : T5Like) :
    -(T30 T5).trace = 0 := by
  rw [T30_trace_zero]
  ring

/-! ### Cubic-layer residue (Newton–Girard form)

The *cubic residue* of the `T_30` spectrum is, by Newton's third
identity rewritten through `p_1 = 0` and `p_3 = 0`,

  `p_3 − (tr T_30)·p_2 + e_2·(tr T_30) − 3·e_3  =  −3·e_3`.

The left-hand side is identically `−3·e_3` (independent of `p_2, e_2`)
once both `tr(T_30) = 0` and `tr(T_30^3) = 0` are substituted. This is
the cubic-level analogue of `T30_power_sum_residue` from the quadratic
layer. -/

/-- **Cubic residue of the `T_30` spectrum.** Identical to
    `T30_e3_collapse`, exposed under the *residue* nomenclature for
    consistency with `T30_power_sum_residue` (quadratic layer). -/
theorem T30_cubic_residue (T5 : T5Like) (p2 e2 e3 : ℝ) :
    ((T30 T5) ^ 3).trace - (T30 T5).trace * p2 + e2 * (T30 T5).trace
      - 3 * e3 = -(3 * e3) :=
  T30_e3_collapse T5 p2 e2 e3

/-! ### Headline summary -/

/-- **Headline summary for the cubic-layer trace identities.**

    The vanishing `tr(T_30^3) = 0` propagates to a *cubic-layer*
    family of identities:

    1. **Direct cube identity** — `tr(T_30^3) = 0`.
    2. **Bilinear cube–any** — for every `a`,
       `tr(T_30^3) · tr(T_30^a) = 0`.
    3. **Cube-of-trace** — `(tr T_30)^3 = 0`.
    4. **Newton–Girard cubic collapse** — for every `p_2, e_2, e_3`,
       `tr(T_30^3) − tr(T_30)·p_2 + e_2·tr(T_30) − 3·e_3 = −3·e_3`
       (identifying `−3·e_3` as the *only* non-trivial third-order
       spectral invariant of `T_30` modulo Newton's third identity).
    5. **Characteristic-polynomial seventh-power coefficient** —
       `−(tr T_30) = 0`. -/
theorem T30_trace_cube_invariants_summary (T5 : T5Like) :
    -- (1) Direct cube identity
    ((T30 T5) ^ 3).trace = 0
    -- (2) Bilinear cube–any
    ∧ (∀ a : ℕ, ((T30 T5) ^ 3).trace * ((T30 T5) ^ a).trace = 0)
    -- (3) Cube of trace
    ∧ ((T30 T5).trace) ^ 3 = 0
    -- (4) Newton–Girard cubic collapse
    ∧ (∀ p2 e2 e3 : ℝ,
        ((T30 T5) ^ 3).trace - (T30 T5).trace * p2 + e2 * (T30 T5).trace
          - 3 * e3 = -(3 * e3))
    -- (5) Characteristic-polynomial seventh-power coefficient
    ∧ -(T30 T5).trace = 0 :=
  ⟨T30_pow_three_trace_zero T5,
   T30_trace_cube_mul_any_zero T5,
   T30_trace_cube_self_zero T5,
   T30_e3_collapse T5,
   T30_char_poly_coeff_seven_zero T5⟩

end PT.Stochastic
