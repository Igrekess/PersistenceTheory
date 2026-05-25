/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Bridge.CauchyMultiplicativeExp
import Mathlib.LinearAlgebra.Matrix.Trace
import Mathlib.LinearAlgebra.Matrix.Notation
import Mathlib.Data.Complex.Basic

/-!
# Finite spectral triple `ST_F` of the PT bifurcation

This module formalises in Lean the finite spectral triple

  `ST_F := (A_F = ℂ, H_F = ℂ², D_F = m·σ_x, γ_F = σ_x)`

that encodes the `q_+/q_-` bifurcation in PT, and proves the
structural identity

  `Tr_F(I_F) = 2 = N_b`

(`N_b` = cardinality of the bifurcation, defined in
`PT/Bridge/CauchyMultiplicativeExp.lean`).

## Conceptual content

In the Connes–Marcolli formulation of noncommutative geometry, a finite
spectral triple `(A, H, D, γ)` consists of an involutive algebra `A`
acting on a finite-dimensional Hilbert space `H`, an off-diagonal Dirac
operator `D : H → H` (self-adjoint, anticommuting with the `ℤ/2`-grading
`γ`), and the grading itself `γ : H → H` (`γ² = I`).

For PT, the bifurcation `q_+ / q_-` is the canonical finite-dimensional
sector that encodes the Higgs (cf. monograph ch37b §18). The
corresponding spectral triple is:

* `A_F = ℂ` — the trivial commutative algebra (the bifurcation acts
  by scalar multiplication on each branch).
* `H_F = ℂ²` — two-dimensional Hilbert space, one direction per
  branch (`q_+, q_-`).
* `D_F = m · σ_x` — off-diagonal Dirac, with `m = Δq = q_- - q_+`
  the bifurcation mass.
* `γ_F = σ_x` — chirality (off-diagonal Pauli).

## What this module accomplishes

* **Concrete definitions** `H_F_PT` (`Fin 2 → ℂ`), `D_F_PT (m)`,
  `γ_F_PT`, all built from Pauli `σ_x`. The dimension is
  `N_b_PT_nat = 2`.

* **Structural trace identity** `Tr_F(I_F) = N_b_PT_nat = 2` is
  proven as a kernel-verified theorem (`Tr_F_identity`).

* Basic kernel-verified algebraic facts about `σ_x`: it is
  self-adjoint over ℂ, squares to the identity, and has trace zero.
  These are the building blocks `D_F` and `γ_F` inherit.

## What this module does NOT accomplish

The derivation of the scale axiom
`PT_cutoff_exponent_eq_neg_one_over_Nb` from `Tr_F(I_F) = N_b`
requires connecting `Tr_F` to the integral `∫₀^∞ u · f(u) du` of the
cutoff via a spectral-action argument (`Tr_F f(D_F²/Λ²)`). A partial
derivation is provided in `ScaleFromFiniteSpectralTriple.lean` and
`SpectralAction.lean`; the full derivation requires developing the
spectral action `Tr f(D²/Λ²)` for arbitrary `f` (substantially
harder, left to future work).

The present module formalises `ST_F` and proves `Tr_F(I_F) = 2`
without claiming the full derivation of the scale axiom. The
**infrastructure** for `ST_F` in Lean is now available, reusable
for any later work in PT-NCG or in Connes-Marcolli SM.

## References

* Monograph ch37b §17-18 (bifurcation spectral structure).
* Connes-Marcolli, *Noncommutative Geometry, Quantum Fields
  and Motives* (AMS Colloquium Publications 55), §1.18.
* Chamseddine-Connes 1996, *The spectral action principle*,
  Phys. Lett. B 412.
-/

namespace PT.Bridge

open Matrix Complex
open scoped BigOperators

/-! ## Dimension of the finite sector -/

/-- **Bifurcation cardinality as a natural number.** The PT bifurcation
    `q_+ / q_-` has exactly two branches; this is the dimension of the
    finite Hilbert space `H_F`. We use `ℕ` here so that we can use it as
    the index type of `Fin N_b_PT_nat`, and provide a one-line bridge to
    the `ℝ`-valued `N_b` of `CauchyMultiplicativeExp` below. -/
abbrev N_b_PT_nat : ℕ := 2

@[simp]
theorem N_b_PT_nat_eq : N_b_PT_nat = 2 := rfl

/-- **Bridge to the real-valued `N_b`.** The `ℝ`-cast of `N_b_PT_nat`
    equals the real number `N_b = 2` used by `CauchyMultiplicativeExp`
    and `CutoffMeanCharacterisation`. This makes the structural
    identification "bifurcation cardinality = cutoff scale" explicit. -/
theorem N_b_PT_nat_cast_eq_N_b : ((N_b_PT_nat : ℕ) : ℝ) = N_b := by
  unfold N_b_PT_nat N_b
  norm_num

/-! ## Hilbert space `H_F = ℂ²` -/

/-- **PT finite Hilbert space.** `H_F_PT := Fin 2 → ℂ`, isomorphic to
    `ℂ²`. One direction for each branch of the `q_+/q_-` bifurcation.

    The standard `Fintype` and `AddCommGroup` / `Module ℂ` instances are
    automatic from `Fin 2 → ℂ`. -/
abbrev H_F_PT : Type := Fin 2 → ℂ

example : Fintype (Fin N_b_PT_nat) := inferInstance

/-! ## Pauli `σ_x` matrix

`σ_x = !![0, 1; 1, 0]` is the workhorse of the bifurcation: it is the
off-diagonal generator that mixes the two branches `q_+, q_-`. It serves
both as the chirality `γ_F` and (multiplied by the mass `m = Δq`) as the
Dirac `D_F`. -/

/-- **Pauli `σ_x` as a `2 × 2` complex matrix.** Off-diagonal involution
    that swaps the two branches. -/
def sigmaX : Matrix (Fin N_b_PT_nat) (Fin N_b_PT_nat) ℂ :=
  !![0, 1; 1, 0]

/-! ### Pointwise entries (definitional `rfl`). -/

@[simp] theorem sigmaX_zero_zero : sigmaX 0 0 = 0 := rfl
@[simp] theorem sigmaX_zero_one  : sigmaX 0 1 = 1 := rfl
@[simp] theorem sigmaX_one_zero  : sigmaX 1 0 = 1 := rfl
@[simp] theorem sigmaX_one_one   : sigmaX 1 1 = 0 := rfl

/-- **`σ_x` is self-adjoint** (conjugate-transpose equals itself).
    Algebraic identity: `(σ_x)†_{ij} = conj((σ_x)_{ji}) = (σ_x)_{ij}`
    since all entries are real. -/
theorem sigmaX_self_adjoint : sigmaX.conjTranspose = sigmaX := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [Matrix.conjTranspose_apply, sigmaX]

/-- **`σ_x² = I`** as a `2 × 2` matrix identity. The chirality squares
    to the identity, as required of a `ℤ/2` grading. -/
theorem sigmaX_sq : sigmaX * sigmaX = 1 := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [Matrix.mul_apply, Fin.sum_univ_two]

/-- **`σ_x` has trace zero.** Both diagonal entries are `0`. -/
theorem sigmaX_trace : Matrix.trace sigmaX = 0 := by
  trans (sigmaX 0 0 + sigmaX 1 1)
  · exact Matrix.trace_fin_two sigmaX
  · simp

/-! ## Chirality `γ_F = σ_x` -/

/-- **PT chirality** in the bifurcation sector: `γ_F := σ_x`. Concretely,
    the chirality interchanges the `q_+` and `q_-` branches.

    (The cleaner Z₂-grading `γ = diag(+1, -1) = σ_z` is also possible
    and leads to the same `Tr_F(I_F) = 2`. We follow the off-diagonal
    convention here for parity with the bifurcation commutator
    identity.) -/
def gammaF_PT : Matrix (Fin N_b_PT_nat) (Fin N_b_PT_nat) ℂ :=
  sigmaX

/-- `γ_F² = I`. Direct corollary of `sigmaX_sq`. -/
theorem gammaF_PT_sq : gammaF_PT * gammaF_PT = 1 :=
  sigmaX_sq

/-- `γ_F` is self-adjoint. -/
theorem gammaF_PT_self_adjoint : gammaF_PT.conjTranspose = gammaF_PT :=
  sigmaX_self_adjoint

/-! ## Dirac operator `D_F = m · σ_x` -/

/-- **PT finite Dirac operator** at mass parameter `m : ℝ`. The off-diagonal
    Dirac is the bifurcation mass `m = Δq = q_- - q_+` times Pauli `σ_x`.

    This is the canonical form: `D_F` is purely off-diagonal (mixes the
    two branches) and reduces to the Higgs mass term in the spectral
    action. -/
noncomputable def D_F_PT (m : ℝ) :
    Matrix (Fin N_b_PT_nat) (Fin N_b_PT_nat) ℂ :=
  (m : ℂ) • sigmaX

/-- `D_F` is self-adjoint. Combines real scalar (its conjugate equals
    itself) with `σ_x` self-adjointness. -/
theorem D_F_PT_self_adjoint (m : ℝ) : (D_F_PT m).conjTranspose = D_F_PT m := by
  unfold D_F_PT
  rw [Matrix.conjTranspose_smul]
  rw [sigmaX_self_adjoint]
  congr 1
  -- Need: star ((m : ℂ)) = (m : ℂ). True because m is real.
  exact Complex.conj_ofReal m

/-- `D_F² = m² · I`. The squared Dirac is a scalar (`m²`) times identity;
    its spectrum is `{m²}` with multiplicity `N_b = 2`.

    This is the key algebraic fact that makes the PT Higgs cutoff
    structural: every state has eigenvalue `m²` under `D_F²`, so the
    spectral action becomes `Tr f(D_F²/Λ²) = 2 · f(m²/Λ²)`. -/
theorem D_F_PT_sq (m : ℝ) :
    D_F_PT m * D_F_PT m = ((m : ℂ)^2) • (1 : Matrix _ _ ℂ) := by
  unfold D_F_PT
  rw [Matrix.smul_mul, Matrix.mul_smul, smul_smul, sigmaX_sq]
  congr 1
  ring

/-! ## Trace identity `Tr_F(I_F) = N_b = 2`

This is the central structural identity of `ST_F`: the trace of the
identity equals the dimension of `H_F`, which equals the cardinality
`N_b = 2` of the bifurcation. -/

/-- **`Tr_F(I_F) = N_b` (main trace identity).**

    The trace of the identity matrix on `H_F = ℂ²` equals the natural
    number `N_b_PT_nat = 2`, the cardinality of the bifurcation. This
    is the structural identity that links the dimension of the finite
    sector to the scale `N_b` of the PT spectral cutoff.

    Proof: standard `Matrix.trace_one` + `Fintype.card (Fin 2) = 2`. -/
theorem Tr_F_identity :
    Matrix.trace (1 : Matrix (Fin N_b_PT_nat) (Fin N_b_PT_nat) ℂ) =
      (N_b_PT_nat : ℂ) := by
  rw [Matrix.trace_one, Fintype.card_fin]

/-- **Real-valued version.** The trace of `I_F` cast to `ℝ` via `.re`
    (it is a real number since it equals `(N_b_PT_nat : ℕ)` cast to ℂ)
    equals the real-valued `N_b` of `CauchyMultiplicativeExp`. -/
theorem Tr_F_identity_real :
    (Matrix.trace (1 : Matrix (Fin N_b_PT_nat) (Fin N_b_PT_nat) ℂ)).re = N_b := by
  rw [Tr_F_identity]
  -- `(N_b_PT_nat : ℂ).re = (N_b_PT_nat : ℝ) = N_b`.
  rw [Complex.natCast_re]
  exact N_b_PT_nat_cast_eq_N_b

/-! ## Spectral action of `D_F²`: degenerate sum

Since `D_F² = m² · I`, the spectral action `Tr_F f(D_F²/Λ²)` is
simply `N_b · f(m²/Λ²)`, irrespective of the choice of `f`. This is
the "degenerate trace" formula that lets the PT bifurcation lift the
SJ-G1 spectral cutoff factorisation to the multiplicative Cauchy
equation. -/

/-- **Spectral action of `D_F²` for the PT bifurcation.** For any function
    `g : ℂ → ℂ` extended via the trace and the fact that
    `D_F² = m² · I`, we have

      `Tr_F(g(D_F²/Λ²)) = N_b · g(m²/Λ²)`.

    Here we record the algebraic specialisation `Tr_F((m²/Λ²) · I) =
    N_b · (m²/Λ²)`, which is the relevant case for the spectral action
    after diagonalisation.

    This wraps `Matrix.trace_smul` + `Tr_F_identity`. -/
theorem Tr_F_scalar (c : ℂ) :
    Matrix.trace (c • (1 : Matrix (Fin N_b_PT_nat) (Fin N_b_PT_nat) ℂ)) =
      (N_b_PT_nat : ℂ) * c := by
  rw [Matrix.trace_smul, Tr_F_identity]
  ring

/-- **Spectral action of `f(D_F²/Λ²)` for the PT bifurcation, real
    parameter version.** When `g : ℝ → ℝ` is applied to `D_F²/Λ² = (m²/Λ²) · I`
    via functional calculus (which on a scalar matrix reduces to
    `g(m²/Λ²) · I`), the trace equals `N_b · g(m²/Λ²)`.

    For now we state this only at the algebraic specialisation
    `Tr_F((g_val : ℂ) • I) = N_b · g_val`. Full functional calculus
    on matrices is future work (see `SpectralAction.lean` for the
    scalar case). -/
theorem Tr_F_scalar_real (g_val : ℝ) :
    Matrix.trace (((g_val : ℂ)) •
        (1 : Matrix (Fin N_b_PT_nat) (Fin N_b_PT_nat) ℂ)) =
      (N_b_PT_nat : ℂ) * (g_val : ℂ) :=
  Tr_F_scalar _

/-! ## Bridge identity between `N_b_PT_nat` and `N_b`

The final ingredient needed to use `ST_F` in the scale-from-dim
derivation (`ScaleFromFiniteSpectralTriple.lean`) is the
identification `N_b_PT_nat = N_b` as real numbers. -/

/-- **Bridge identity.** The cardinality of the bifurcation (as ℝ) equals
    the structural scale `N_b` of `CauchyMultiplicativeExp`. -/
theorem N_b_PT_eq_N_b : ((N_b_PT_nat : ℕ) : ℝ) = N_b :=
  N_b_PT_nat_cast_eq_N_b

end PT.Bridge
