/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Bridge.CutoffMeanCharacterisation

/-!
# Uniqueness of the PT spectral cutoff `f_PT(u) = exp(-u/N_b)`

This module assembles the **PT structural axioms** (multilinear SJ G1
factorisation, scale `N_b`) with the kernel-verified Cauchy theorem
`cauchy_mult_exp` (`PT/Bridge/CauchyMultiplicativeExp.lean`) to
**conclude formally** that the PT spectral cutoff is uniquely
`f_PT(u) = exp(-u/N_b)`.

## Axiomatic history

* The initial formulation posited the cutoff uniqueness via **three**
  ad-hoc axioms: `CRT_ShoreJohnson_cutoff_factorises` (bilinear
  `f(x+y) = f(x)·f(y)`), `PT_cutoff_decay` (`f(1) < f(0)`),
  `PT_cutoff_at_one_eq_exp_minus_one_over_Nb` (`f(1) = exp(-1/N_b)`).
* The current refactor factors the bilinear axiom through the
  **multilinear** Shore--Johnson G1 spectral axiom
  (`PT/Bridge/ShoreJohnsonG1Spectral.lean`) and replaces the
  value-at-one axiom by the **structural scale axiom**
  (`PT/Bridge/CutoffMeanCharacterisation.lean`):

  | Before (ad-hoc) | After (structural) |
  |---|---|
  | `axiom CRT_ShoreJohnson_cutoff_factorises` (bilinear) | derived via multilinear `SJG1_spectral_cutoff_factorises` |
  | `axiom PT_cutoff_at_one_eq_exp_minus_one_over_Nb` | derived via `PT_cutoff_exponent_eq_neg_one_over_Nb` |
  | `axiom PT_cutoff_decay` | retained (definitional UV regulator) |

  **Net axiomatic reduction:** the proof of
  `cutoff_PT_unique_eq_cutoffPT` no longer depends on
  `PT_cutoff_decay` (it now flows entirely through the scale axiom).
  So the principal theorem rests on **2 PT axioms** instead of 3.

## Status of the conclusion

* Theorem `cutoff_PT_unique_eq_cutoffPT` is **kernel-verified**
  (0 sorry) given the two structural axioms
  `SJG1_spectral_cutoff_factorises`
  and `PT_cutoff_exponent_eq_neg_one_over_Nb`.
* The combination of (a) kernel-verified theorem + (b) PT structural
  axioms gives **K4 [DER strict modulo SJ-cutoff]**.
* Removing both remaining axioms requires a deeper formalisation of
  `ST_F` and the spectral action in Lean (work in progress).

## References

* PT monograph: appendix Y §Y.13.6.
* Cauchy core: `PT/Bridge/CauchyMultiplicativeExp.lean`.
* Multilinear SJ G1: `PT/Bridge/ShoreJohnsonG1Spectral.lean`.
* Scale axiom: `PT/Bridge/CutoffMeanCharacterisation.lean`.
-/

namespace PT.Bridge

open Real

/-! ## Definition of a PT spectral cutoff is in `CauchyMultiplicativeExp` -/

-- The structure `IsPTSpectralCutoff` is defined in `ShoreJohnsonG1Spectral.lean`'s
-- dependency chain (via `CauchyMultiplicativeExp` re-export). It is restated
-- below for documentation continuity but not re-declared.

/-! ## The retained decay axiom (P1)

After the refactor, the only ad-hoc axiom of this module that is
*not* factored through a more structural form is `PT_cutoff_decay`.
This is because decay is the *definitional* property of any UV
regulator function in a spectral action approach
(Chamseddine--Connes 1996, eq.\ 2.5); it is not a substantive
content of Shore--Johnson G1.

The principal theorem `cutoff_PT_unique_eq_cutoffPT` does **not**
depend on this axiom after the refactor. We retain it because it
is used by the corollary `PT_spectral_cutoff_decay`. -/

/-- **Axiom (P1) — UV regulator decay.** A PT spectral cutoff strictly
    decreases between `0` and `1` (in particular `f(1) < f(0)`).

    This is a definitional requirement of any UV regulator function in
    spectral action approaches (Chamseddine--Connes 1996, eq.~2.5). -/
axiom PT_cutoff_decay
    (f : ℝ → ℝ) (hf : IsPTSpectralCutoff f) :
    f 1 < f 0

/-! ## Compatibility re-exports (post-refactor)

The two previously ad-hoc axioms are now **theorems**, derived from
the multilinear SJ G1 axiom and the structural scale axiom. For
backwards compatibility with downstream modules that use the original
names, we re-export them as `theorem` (no longer `axiom`). -/

/-- **Theorem (formerly ad-hoc axiom).** Any PT spectral cutoff `f`
    satisfies the bilinear multiplicative Cauchy equation
    `f(x + y) = f(x) · f(y)`.

    **Post-refactor status.** This was the ad-hoc bilinear axiom
    `CRT_ShoreJohnson_cutoff_factorises` of the original formulation.
    It is now derived from the more fundamental multilinear axiom
    `SJG1_spectral_cutoff_factorises` (`ShoreJohnsonG1Spectral.lean`)
    by specialisation at `k = 2`. -/
theorem CRT_ShoreJohnson_cutoff_factorises
    (f : ℝ → ℝ) (hf : IsPTSpectralCutoff f) :
    ∀ x y, f (x + y) = f x * f y :=
  cauchy_bilinear_from_SJG1 f hf

/-- **Theorem (formerly ad-hoc axiom).** A PT spectral cutoff satisfies
    `f(1) = exp(-1/N_b)`.

    **Post-refactor status.** This was the ad-hoc value-at-one axiom
    `PT_cutoff_at_one_eq_exp_minus_one_over_Nb` of the original
    formulation. It is now derived from the structural scale axiom
    `PT_cutoff_exponent_eq_neg_one_over_Nb`
    (`CutoffMeanCharacterisation.lean`) combined with the multilinear
    SJ G1 axiom via `cauchy_mult_exp`. -/
theorem PT_cutoff_at_one_eq_exp_minus_one_over_Nb
    (f : ℝ → ℝ) (hf : IsPTSpectralCutoff f) :
    f 1 = Real.exp (-1 / N_b) :=
  cutoff_at_one_from_scale_axiom f hf

/-! ## Main uniqueness theorem

The principal theorem now rests on **2 PT axioms** instead of 3:
* `SJG1_spectral_cutoff_factorises` (multilinear, structural)
* `PT_cutoff_exponent_eq_neg_one_over_Nb` (scale, structural)

The decay axiom `PT_cutoff_decay` is no longer in the dependency chain
of `cutoff_PT_unique_eq_cutoffPT` (it is retained only for the corollary). -/

/-- **Main theorem.** Under the two structural PT axioms
    (multilinear SJ G1 + scale `= N_b`), every PT spectral cutoff `f`
    equals `cutoffPT(u) = exp(-u/N_b)` pointwise.

    **Proof.** Direct via
    `cutoff_eq_cutoffPT_from_structural_axioms` from
    `CutoffMeanCharacterisation.lean`. The chain is:
    1. Multilinear SJ G1 → bilinear Cauchy (k = 2 specialisation).
    2. Bilinear Cauchy + continuity + positivity → `f(x) = exp(c·x)`
       (kernel-verified `cauchy_mult_exp`).
    3. Scale axiom `c = -1/N_b` → pointwise identification with `cutoffPT`. -/
theorem cutoff_PT_unique_eq_cutoffPT
    (f : ℝ → ℝ) (hf : IsPTSpectralCutoff f) :
    ∀ x, f x = cutoffPT x :=
  cutoff_eq_cutoffPT_from_structural_axioms f hf

/-! ## Corollaries -/

/-- **Decay corollary.** The PT cutoff decays strictly when the argument
    increases: `x < y ⟹ cutoffPT y < cutoffPT x`. -/
theorem cutoffPT_decay (x y : ℝ) (h : x < y) : cutoffPT y < cutoffPT x := by
  unfold cutoffPT
  apply Real.exp_lt_exp.mpr
  -- Goal: `-y / N_b < -x / N_b`. We pass through `-(y/N_b) < -(x/N_b)`
  -- which follows from `x/N_b < y/N_b` (since `N_b > 0`) by negation.
  have hNb_pos : 0 < N_b := N_b_pos
  have h1 : x / N_b < y / N_b :=
    (div_lt_div_iff_of_pos_right hNb_pos).mpr h
  rw [neg_div, neg_div]
  linarith

/-- **Decay corollary applied to the original cutoff.** Any PT spectral
    cutoff inherits the strict decay of `cutoffPT`. -/
theorem PT_spectral_cutoff_decay
    (f : ℝ → ℝ) (hf : IsPTSpectralCutoff f) (x y : ℝ) (h : x < y) :
    f y < f x := by
  rw [cutoff_PT_unique_eq_cutoffPT f hf x,
      cutoff_PT_unique_eq_cutoffPT f hf y]
  exact cutoffPT_decay x y h

end PT.Bridge
