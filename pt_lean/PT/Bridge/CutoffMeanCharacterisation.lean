/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Bridge.ShoreJohnsonG1Spectral

/-!
# Characteristic scale of the PT spectral cutoff

This module is part of the axiomatic refactoring of
`HiggsCutoffUniqueness.lean` so that the value-at-one axiom

  `PT_cutoff_at_one_eq_exp_minus_one_over_Nb : f 1 = exp(-1/N_b)`

is replaced by a more **structural** scale axiom (`τ = N_b`) together
with a kernel-verified **theorem** deriving the original value-at-one
statement.

## Conceptual content

By the chain
* `SJG1_spectral_cutoff_factorises` (`ShoreJohnsonG1Spectral.lean`, axiom),
* `cauchy_bilinear_from_SJG1` (derived theorem),
* `cauchy_mult_exp` (`CauchyMultiplicativeExp.lean`, kernel-verified
  theorem, 0 PT axiom),

any PT spectral cutoff `f` has the exponential form `f(x) = exp(c · x)`
for some `c ∈ ℝ`. The remaining question is the value of `c`.

The PT axiom states that the characteristic scale of the cutoff is
`N_b = 2`, the cardinality of the `q_+/q_-` bifurcation, fixed by
`Tr_F(I_F) = 2` in the finite spectral triple `ST_F`. In terms of
`c`, this means `c = -1/N_b`.

## What this module accomplishes

* **Axiom `PT_cutoff_exponent_eq_neg_one_over_Nb`** — the structural
  axiom that the exponent `c` in `f(x) = exp(c · x)` equals `-1/N_b`.
  This is the **purest form** of the scale axiom: it acts on `c`
  directly (a real number) rather than on `f(1)` (a function value).

* **Theorem `cutoff_at_one_from_scale_axiom`** — derives the original
  statement `f(1) = exp(-1/N_b)`, previously posited as an axiom.

* **Theorem `cutoff_pointwise_from_scale_axiom`** — derives the
  stronger pointwise statement `f(x) = exp(-x/N_b)` for all `x`.

## Epistemic progress

Before this module, `HiggsCutoffUniqueness.lean` declared an
*ad-hoc value-at-one* axiom. After this module, the dependency chain
reads:

| Statement | Status | Module |
|---|---|---|
| `PT_cutoff_exponent_eq_neg_one_over_Nb` | **axiom** (purest form of scale axiom) | this file |
| `cutoff_at_one_from_scale_axiom` | **theorem** | this file |
| `cutoff_pointwise_from_scale_axiom` | **theorem** | this file |
| `cutoff_PT_unique_eq_cutoffPT` | theorem given multilinear + decay + scale axioms | `HiggsCutoffUniqueness.lean` |

The PT axiomatic basis is unchanged in number but improved in
conceptual character: the scale axiom now captures the structural
content "the scale is `N_b`" rather than the bilinear value
"`f(1) = exp(-1/N_b)`".

## What this module does NOT accomplish

A full closure would derive the scale axiom from `Tr_F(I_F) = 2 =
N_b` inside a formalised finite spectral triple `ST_F`. That
formalisation is provided in `FiniteSpectralTriple.lean` and
`ScaleFromFiniteSpectralTriple.lean`.

## References

* PT monograph: appendix Y §Y.13.6, `prop:E7_unicity_f` (P4 condition).
-/

namespace PT.Bridge

open Real

/-! ## The structural scale axiom -/

/-- **Axiom (Characteristic scale = `N_b`).**

    For any PT spectral cutoff `f`, the exponent `c` in the exponential
    form `f(x) = exp(c · x)` (provided by
    `cauchy_mult_exp` via `SJG1_spectral_cutoff_factorises`) equals
    `-1/N_b`, where `N_b = 2` is the cardinality of the `q_+/q_-`
    bifurcation.

    **Justification.** The finite spectral triple of the bifurcation,
    `ST_F = (ℂ², ℂ², m·σ_x, σ_x)` (cf.\ monograph ch37b §17), has
    `Tr_F(I_F) = 2 = N_b` levels. This fixes the unique natural scale
    `τ = N_b` for the cutoff function, equivalently `c = -1/τ =
    -1/N_b` in the exponential parametrisation.

    Equivalent formulations:
    * **Mean** `∫₀^∞ u·f(u) du / ∫₀^∞ f(u) du = N_b` (the canonical
      mean of an exponential `(1/τ) exp(-u/τ)` is `τ`).
    * **e-folding scale** `f(N_b) = e^{-1} · f(0)`.
    * **Value at one** `f(1) = exp(-1/N_b)` (derived as a theorem below).

    This axiom captures the *structural* content of `Tr_F(I_F) = N_b`
    without requiring a full Lean formalisation of `ST_F` (a
    formalisation that is itself work in progress; see
    `FiniteSpectralTriple.lean` for the current state).

    **Why "exponent" rather than "value at one"?** The exponent `c`
    is a coordinate-free characterisation of the scale, whereas
    `f(1)` is sensitive to the choice of unit on the spectral
    coordinate `u`. Axiomatising the exponent is the *epistemically
    pure* form of the scale fixing. -/
axiom PT_cutoff_exponent_eq_neg_one_over_Nb
    (f : ℝ → ℝ) (hf : IsPTSpectralCutoff f)
    (c : ℝ) (hc : ∀ x, f x = Real.exp (c * x)) :
    c = -1 / N_b

/-! ## Derived theorems -/

/-- **Value-at-one, derived.** Combining the scale axiom with the
    exponential form (from `cauchy_mult_exp` + multilinear SJ G1)
    yields the value `f(1) = exp(-1/N_b)`, previously posited as the
    axiom `PT_cutoff_at_one_eq_exp_minus_one_over_Nb`. -/
theorem cutoff_at_one_from_scale_axiom
    (f : ℝ → ℝ) (hf : IsPTSpectralCutoff f) :
    f 1 = Real.exp (-1 / N_b) := by
  -- Step 1: get exponential form `f(x) = exp(c · x)` from
  -- `cauchy_mult_exp` via the multilinear SJ G1.
  obtain ⟨c, hc⟩ := cutoff_to_cauchy_mult_exp f hf
  -- Step 2: the scale axiom says `c = -1/N_b`.
  have hc_eq : c = -1 / N_b := PT_cutoff_exponent_eq_neg_one_over_Nb f hf c hc
  -- Step 3: evaluate at `x = 1`.
  have h1 : f 1 = Real.exp (c * 1) := hc 1
  rw [hc_eq] at h1
  simpa using h1

/-- **Pointwise form, derived.** Same chain, but giving the full
    pointwise identification `f(x) = exp(-x/N_b)`. -/
theorem cutoff_pointwise_from_scale_axiom
    (f : ℝ → ℝ) (hf : IsPTSpectralCutoff f) :
    ∀ x, f x = Real.exp (-x / N_b) := by
  obtain ⟨c, hc⟩ := cutoff_to_cauchy_mult_exp f hf
  have hc_eq : c = -1 / N_b := PT_cutoff_exponent_eq_neg_one_over_Nb f hf c hc
  intro x
  rw [hc x, hc_eq]
  -- `exp((-1/N_b) * x) = exp(-x/N_b)`.
  congr 1
  ring

/-! ## Bridge to `cutoffPT` -/

/-- **Pointwise identification with `cutoffPT`.** The PT spectral
    cutoff `f` equals the canonical `cutoffPT(u) = exp(-u/N_b)`
    pointwise, derived from the multilinear SJ G1 axiom + the
    scale axiom (no use of the bilinear ad-hoc axiom or the
    value-at-one ad-hoc axiom). -/
theorem cutoff_eq_cutoffPT_from_structural_axioms
    (f : ℝ → ℝ) (hf : IsPTSpectralCutoff f) :
    ∀ x, f x = cutoffPT x := by
  intro x
  unfold cutoffPT
  exact cutoff_pointwise_from_scale_axiom f hf x

end PT.Bridge
