/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Bridge.FiniteSpectralTriple
import PT.Bridge.CutoffMeanCharacterisation

/-!
# Scale of the PT spectral cutoff from `ST_F`

This module factors the structural axiom

  `PT_cutoff_exponent_eq_neg_one_over_Nb : c = -1/N_b`
  (`PT/Bridge/CutoffMeanCharacterisation.lean`)

through the **dimension** of the finite spectral triple `ST_F`
formalised in `PT/Bridge/FiniteSpectralTriple.lean`. The new axiom
expresses the scale of the cutoff as `c = -1/Tr_F(I_F)`, where
`Tr_F(I_F)` is a **kernel-verified theorem** giving the dimension of
`H_F = ℂ²`.

## What this module accomplishes

* **Axiom `cutoff_scale_from_dim_FST`** — the structural axiom that
  the exponent `c` of any PT spectral cutoff equals
  `-1 / (Re Tr_F(I_F))`, where `Tr_F(I_F)` is the dimension of the
  bifurcation Hilbert space. The dimension itself is no longer an
  opaque parameter `N_b` but the **kernel-verified theorem**
  `Tr_F_identity_real` of the finite spectral triple.

* **Theorem `PT_cutoff_exponent_eq_neg_one_over_Nb_from_FST`** —
  derives the original axiom statement from the new one via the
  identity `(Re Tr_F(I_F)) = N_b` (`Tr_F_identity_real`).

## What this module does NOT accomplish

This refactor **does not eliminate** the structural axiom. It
**re-expresses** the axiom in terms of the dimension of `ST_F`,
which is itself a kernel-verified theorem. This is a *conceptual*
progress (the axiom now explicitly invokes the spectral-triple
dimension rather than an opaque `N_b`) but it does not derive the
identification "cutoff scale = dimension of H_F" from first
principles.

A full derivation of `cutoff_scale_from_dim_FST` from more
fundamental principles would require:

1. **Spectral action principle** (Chamseddine-Connes 1996): the
   variational origin of the cutoff is `S = Tr f(D²/Λ²)` extremised
   over `(f, Λ)`. This requires formalising `Tr f(D²/Λ²)` as a Lean
   object (work in progress).

2. **Shore-Johnson G1 extension to spectral inference**: the cutoff
   must factorise on CRT-independent eigenvalues, which forces the
   scale to coincide with the natural scale of the finite sector
   (i.e., the dimension). This is the conceptual content of a
   future extension.

Without (1) and (2), the new axiom remains structurally axiomatic.
What we gain here is **explicit linkage** between the cutoff scale
and the kernel-verified dimension `Tr_F(I_F) = 2`, replacing an
opaque parameter with a structurally identified quantity.

## Epistemic status

**Before this refactor**: K4 depends on
`PT_cutoff_exponent_eq_neg_one_over_Nb` (opaque scale `c = -1/N_b`).

**After this refactor**: K4 depends on `cutoff_scale_from_dim_FST`
(structural scale `c = -1/dim(H_F)`), where `dim(H_F) = 2 = N_b` is
a kernel-verified theorem of `ST_F`.

Net axiomatic basis: **1 axiom** in the cutoff-mean characterisation
chain, but now this axiom is expressed structurally (via the FST
dimension) rather than via an opaque parameter. Quantitative count
is unchanged; qualitative character is improved.

## References

* `PT/Bridge/FiniteSpectralTriple.lean` (`Tr_F_identity` =
  kernel-verified).
* Connes-Marcolli, *Noncommutative Geometry, Quantum Fields and
  Motives* (AMS Colloquium 55), §1.18 (Chamseddine-Connes spectral
  action).
-/

namespace PT.Bridge

open Real Matrix

/-! ## The structural scale axiom (FST version) -/

/-- **Axiom (scale of PT cutoff from `ST_F` dimension).**

    For any PT spectral cutoff `f` with exponential form
    `f(x) = exp(c · x)` (provided by `cauchy_mult_exp` via the
    multilinear Shore-Johnson G1 axiom), the exponent `c` equals
    minus the reciprocal of the **dimension** of the bifurcation
    Hilbert space `H_F = ℂ²`, encoded by
    `(Re Tr_F(I_F))` where `I_F` is the identity on `H_F`.

    **Justification.** The Chamseddine-Connes spectral action
    `Tr_F(f(D_F²/Λ²))` on `ST_F = (ℂ, ℂ², m·σ_x, σ_x)` reduces to
    `(dim H_F) · f(m²/Λ²) = (Re Tr_F(I_F)) · f(m²/Λ²)` since
    `D_F² = m²·I` (`D_F_PT_sq` in `FiniteSpectralTriple.lean`).
    The Shore-Johnson G1 axiom applied to this spectral inference
    then forces the cutoff scale to coincide with the **structural
    scale** of the finite sector, which is `dim H_F = (Re Tr_F(I_F))
    = N_b = 2` (`Tr_F_identity_real`).

    This axiom captures the *structural* content of the Chamseddine-
    Connes scale identification, expressed in terms of the
    kernel-verified dimension of `H_F` (no opaque parameter `N_b`).

    A full Lean derivation of this axiom from more fundamental
    spectral-action principles would require formalising
    `Tr f(D²/Λ²)` for arbitrary `f` and the Shore-Johnson G1 spectral
    extension. Until then, this is taken axiomatically. -/
axiom cutoff_scale_from_dim_FST
    (f : ℝ → ℝ) (hf : IsPTSpectralCutoff f)
    (c : ℝ) (hc : ∀ x, f x = Real.exp (c * x)) :
    c = -1 / (Matrix.trace (1 : Matrix (Fin N_b_PT_nat) (Fin N_b_PT_nat) ℂ)).re

/-! ## Derived theorems -/

/-- **Theorem (formerly an axiom in `CutoffMeanCharacterisation.lean`).**

    The original axiom `PT_cutoff_exponent_eq_neg_one_over_Nb`,
    asserting `c = -1/N_b`, is now derived from the FST scale axiom
    via `Tr_F_identity_real : (Re Tr_F(I_F)) = N_b`. -/
theorem PT_cutoff_exponent_eq_neg_one_over_Nb_from_FST
    (f : ℝ → ℝ) (hf : IsPTSpectralCutoff f)
    (c : ℝ) (hc : ∀ x, f x = Real.exp (c * x)) :
    c = -1 / N_b := by
  have h_fst := cutoff_scale_from_dim_FST f hf c hc
  -- `(Re Tr_F(I_F)) = N_b` (kernel-verified in FiniteSpectralTriple).
  rw [Tr_F_identity_real] at h_fst
  exact h_fst

/-! ## Bridge: equivalence of the old and new axioms

The two axioms `PT_cutoff_exponent_eq_neg_one_over_Nb` (opaque scale)
and `cutoff_scale_from_dim_FST` (structural scale via `Tr_F(I_F)`)
are formally equivalent given `Tr_F_identity_real`. The new axiom
is preferred because it explicitly links the scale to a
kernel-verified property of `ST_F` rather than to an opaque
parameter. -/

/-- **Equivalence of the two scale axioms.** Given the kernel-verified
    identity `Tr_F_identity_real`, the two axioms are equivalent.
    This documents that the refactor of `CutoffMeanCharacterisation`
    via `ScaleFromFiniteSpectralTriple` is sound: the same content is
    expressed in two notations. -/
theorem scale_axiom_equivalence
    (_f : ℝ → ℝ) (_hf : IsPTSpectralCutoff _f) (c : ℝ)
    (_hc : ∀ x, _f x = Real.exp (c * x)) :
    (c = -1 / N_b) ↔
    (c = -1 / (Matrix.trace (1 : Matrix (Fin N_b_PT_nat) (Fin N_b_PT_nat) ℂ)).re) := by
  rw [Tr_F_identity_real]

end PT.Bridge
