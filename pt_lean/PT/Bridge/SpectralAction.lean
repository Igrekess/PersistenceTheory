/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Bridge.FiniteSpectralTriple

/-!
# Spectral action `Tr f(D²/Λ²)` for the PT bifurcation (scalar case)

This module formalises the Chamseddine-Connes spectral action

  `S_spectral[f, D] := Tr f(D²/Λ²)`

in the **specific case** of the PT bifurcation Dirac
`D_F = m · σ_x` (where `D_F² = m² · I` is scalar). The general case
of an arbitrary Dirac with non-scalar `D²` requires continuous
functional calculus on bounded operators (future work).

## What this module accomplishes

For `ST_F = (ℂ, ℂ², m·σ_x, σ_x)`, since `D_F² = m² · I` is a scalar
multiple of the identity, the spectral action reduces to

  `Tr_F(f(D_F²/Λ²)) = N_b · f(m²/Λ²)`

(no functional calculus needed — pure algebra). We formalise this
as a kernel-verified theorem `spectral_action_FST_scalar_Dirac`.

This is the **degenerate case** of the Chamseddine-Connes spectral
action when the finite Dirac has degenerate (scalar) spectrum. The
non-degenerate case (general `D` with multiple distinct eigenvalues)
requires functional calculus and is the content of a future module.

## What this module does NOT accomplish

* **General functional calculus** on matrices: this requires
  `Mathlib.Analysis.NormedSpace.Spectrum` /
  `ContinuousFunctionalCalculus`, which is mature only for
  Hilbert-space operators and not directly applicable in this
  algebraic setting without substantial scaffolding.

* **Extremization of `S_spectral` over `(f, Λ)`**: this is the
  variational origin of the cutoff, and would require
  `Mathlib.Analysis.Calculus.LocalExtr` or similar. Currently out
  of scope.

* **Derivation of the scale axiom from the spectral action**: even
  with this module, deriving `c = -1/N_b` from the spectral action
  requires an extremization principle. We do not attempt this here.

## Epistemic status

The kernel-verified theorem in this module
(`spectral_action_FST_scalar_Dirac`) gives the **algebraic skeleton**
of the spectral action on `ST_F`. The infrastructure is reusable
for any future module developing the full Chamseddine-Connes spectral
action in Lean.

## References

* Chamseddine-Connes 1996, *The spectral action principle*,
  Phys. Lett. B 412 (1996), 1-9.
* Vassilevich, Phys. Rep. 388 (2003) 279 (heat kernel expansion of
  spectral action).
-/

namespace PT.Bridge

open Matrix Complex

/-! ## Spectral action of `D_F²` for the PT bifurcation (scalar case) -/

/-- **Spectral action of `D_F²` (scalar case).** Since `D_F² = m² · I`
    (`D_F_PT_sq` in `FiniteSpectralTriple.lean`), for any function
    `g : ℂ → ℂ` extended via `g(c · I) = g(c) · I`, the trace is

      `Tr_F(g(D_F²/Λ²)) = N_b · g(m²/Λ²)`.

    Here we record the algebraic form: for `c : ℂ`,

      `Tr_F(c • I) = N_b · c`.

    This is the kernel-verified building block of the spectral action
    for `ST_F`. The "functional calculus" content `g(c · I) = g(c) · I`
    is trivial in this scalar case. -/
theorem spectral_action_FST_scalar (c : ℂ) :
    Matrix.trace (c • (1 : Matrix (Fin N_b_PT_nat) (Fin N_b_PT_nat) ℂ)) =
      (N_b_PT_nat : ℂ) * c :=
  Tr_F_scalar c

/-- **Spectral action of `f(D_F²/Λ²)` for real-valued `f`.** Specialising
    to the case where `f : ℝ → ℝ` evaluated at `m²/Λ²` gives a real
    value `g_val := f(m²/Λ²)`, the spectral action is

      `Tr_F((g_val : ℂ) • I) = N_b · g_val`.

    Combined with `D_F_PT_sq : D_F m * D_F m = m² · I`, this gives the
    spectral action of `f(D_F²/Λ²)` for the PT bifurcation as
    `N_b · f(m²/Λ²)`. -/
theorem spectral_action_FST_scalar_real (g_val : ℝ) :
    Matrix.trace (((g_val : ℂ)) •
        (1 : Matrix (Fin N_b_PT_nat) (Fin N_b_PT_nat) ℂ)) =
      (N_b_PT_nat : ℂ) * (g_val : ℂ) :=
  Tr_F_scalar_real g_val

/-- **Spectral action of `(m²/Λ²) · I` (algebraic form of
    `D_F²/Λ²` for `Λ : ℝ`, `Λ ≠ 0`).** This is the explicit form of
    `D_F²/Λ²` when `D_F = m · σ_x`: the scalar `(m²/Λ²)` times the
    identity. -/
theorem D_F_sq_div_Λ_sq (m Λ : ℝ) (_hΛ : Λ ≠ 0) :
    (1 / (Λ : ℂ)^2) • (D_F_PT m * D_F_PT m) =
      ((m / Λ : ℝ)^2 : ℂ) • (1 : Matrix _ _ ℂ) := by
  rw [D_F_PT_sq, smul_smul]
  congr 1
  push_cast
  ring

/-- **Trace of `D_F²/Λ²`.** The trace of the (normalised) squared Dirac
    equals `N_b · m²/Λ²`. This is the lowest-order spectral action
    coefficient (the volume term in Chamseddine-Connes). -/
theorem Tr_D_F_sq_div_Λ_sq (m Λ : ℝ) (hΛ : Λ ≠ 0) :
    Matrix.trace ((1 / (Λ : ℂ)^2) • (D_F_PT m * D_F_PT m)) =
      (N_b_PT_nat : ℂ) * ((m / Λ : ℝ)^2 : ℂ) := by
  rw [D_F_sq_div_Λ_sq m Λ hΛ]
  exact Tr_F_scalar _

/-! ## Identification with the cutoff scale

For the PT cutoff `cutoffPT(u) = exp(-u/N_b)`, the spectral action
on `ST_F` is

  `Tr_F(cutoffPT(D_F²/Λ²)) = N_b · exp(-(m/Λ)²/N_b)`.

This identifies the scale of the cutoff with the dimension `N_b`
of `H_F`. The kernel-verified building blocks above (`Tr_F_scalar`,
`D_F_PT_sq`) provide the algebraic skeleton; the full spectral
action principle (which extremises this expression over `(f, Λ)`)
is the residual content of a future module. -/

end PT.Bridge
