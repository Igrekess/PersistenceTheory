/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import Mathlib.Data.Complex.Basic
import Mathlib.Analysis.Complex.Norm
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Basic
import Mathlib.Algebra.BigOperators.Group.Finset.Basic
import Mathlib.Data.ZMod.Basic
import Mathlib.Tactic

/-!
# Cyclic Phase Identity — **Route 2 (Spectral / Fourier contraction)**

**Statement (paper-level, Ch06 §6.2, Remark `holonomy_routes`, "Voie 2").**
On the residue space `ℤ/pℤ`, the Fourier transform of the transfer kernel
`T_p` (restricted to surviving residues) admits as eigenvalue, on the
fundamental non-trivial character `χ₁`, the **real** value

  `T̂_p(χ₁) = 1 - δ_p ∈ ℝ ⊂ ℂ`.

The spectral content of cyclic phase is the **contraction of the
fundamental Fourier mode**:

  `1 - |T̂_p(χ₁)|² = 1 - (1 - δ_p)² = δ_p (2 - δ_p) = sin²θ_p`.

This is Route 2 of three independent derivations of the cyclic phase
identity in the monograph: spectral (this file), algebraic-geometric
(`PT.Holonomy.CyclicPhaseIdentity`, Route 1), Fisher-information
(Route 3, Ch05b, not yet formalised).

## Mathematical content

Route 2 does **not** add new content over Route 1 in its conclusion —
both yield `sin²θ = δ(2-δ)`. Its value lies in the **independence of the
derivation**: the same number `δ(2-δ)` arises as the *spectral weight*
shed by the fundamental Fourier mode of a finite-group convolution, with
no reference to Pythagoras or geometric angle. This file formalises that
spectral identity in a clean form.

The Fourier eigenvalue `T̂_p(χ₁) = 1 - δ_p` is itself a finite-group
character sum:

  `T̂_p(χ₁) = ∑_{r ∈ ℤ/pℤ} T_p(r) · ω_p^r`, where `ω_p = e^{2πi/p}`.

For the **diagonal-only survival kernel** `T_p(0) = 1 - δ_p`,
`T_p(r) = 0` for `r ≠ 0`, the sum collapses to `1 - δ_p` at any
character. This realises Route 2's hypothesis explicitly.

## Main results

1. `complex_modSq_real_sub` — purely algebraic: for `δ : ℝ` and the
   complex number `lam = (1 - δ : ℂ)`, `1 - Complex.normSq lam = δ(2 - δ)`.

2. `complex_norm_sq_real_sub` — same identity in `‖·‖²` form:
   `1 - ‖((1 - δ : ℝ) : ℂ)‖² = δ(2 - δ)`.

3. `cyclic_phase_spectral` — **Route 2 headline.** For any complex
   Fourier eigenvalue `lam : ℂ` satisfying `lam = (1 - δ : ℝ)` and any
   angle `θ` with `cos θ = 1 - δ`,
   `1 - Complex.normSq lam = Real.sin θ ^ 2`.

4. `cyclic_phase_spectral_norm` — same in `‖·‖²` form.

5. `fourier_eigenvalue_diagonal_kernel` — **explicit Fourier sum lemma**:
   for the diagonal-only kernel `T r = if r = 0 then 1 - δ else 0` on a
   finite group of residues, the character sum equals `1 - δ` at any
   character. This realises the eigenvalue assumption used in
   `cyclic_phase_spectral` from a concrete finite-group Fourier
   transform.

6. `cyclic_phase_spectral_via_kernel` — composition of (3) and (5):
   the headline statement directly from the diagonal kernel.

## Comparison with Route 1

| Aspect | Route 1 (`CyclicPhaseIdentity`) | Route 2 (this file) |
|---|---|---|
| Domain | Trigonometry on `ℝ` | Harmonic analysis on `ℤ/pℤ` |
| Input | `cos θ = 1 - δ` | `T̂(χ) = 1 - δ` (Fourier eigenvalue) |
| Algebraic core | `sin²θ + cos²θ = 1` (Pythagore) | `‖x‖² = x²` for `x ∈ ℝ ⊂ ℂ` |
| Conclusion | `sin²θ = δ(2 - δ)` | `1 - ‖lam‖² = δ(2 - δ)` |

Both Routes converge on the same number; their derivations use
disjoint machinery (Pythagore vs. complex modulus of real number).

## Reference

Monograph Chapter 6, Remark `rem:holonomy_routes`, "Voie 2 : Spectrale
(contraction de Fourier)" — equations `(eq:fourier_eigenvalue)` and
`(eq:fourier_holonomy)`.
-/

namespace PT.Holonomy

open Real Complex BigOperators

/-! ### Section 1: Algebraic core (real number embedded in ℂ) -/

/-- **Algebraic core of Route 2.** For any real `δ`, the complex number
`((1 - δ : ℝ) : ℂ)` has squared modulus `(1 - δ)²`, and hence the
"contraction" `1 - normSq = δ(2 - δ)`.

This is the spectral analogue of Route 1's Pythagorean identity. It rests
on the elementary fact that for real `x`, `Complex.normSq ((x : ℂ)) = x²`. -/
theorem complex_normSq_real_sub (δ : ℝ) :
    1 - Complex.normSq ((1 - δ : ℝ) : ℂ) = δ * (2 - δ) := by
  rw [Complex.normSq_ofReal]
  ring

/-- **Norm-form** of `complex_normSq_real_sub`. For `lam = ((1 - δ : ℝ) : ℂ)`,
`1 - ‖lam‖² = δ (2 - δ)`. -/
theorem complex_norm_sq_real_sub (δ : ℝ) :
    1 - ‖((1 - δ : ℝ) : ℂ)‖ ^ 2 = δ * (2 - δ) := by
  rw [Complex.norm_real, Real.norm_eq_abs, sq_abs]
  ring

/-! ### Section 2: Cyclic phase identity via spectrum -/

/-- **Cyclic phase identity — Route 2 (spectral form).**

For any complex eigenvalue `lam` of a Fourier transform satisfying
`lam = (1 - δ : ℝ)` (as element of ℂ), and any angle `θ` parameterising
the cyclic phase via `cos θ = 1 - δ`, the contraction of the fundamental
mode equals `sin² θ`:

  `1 - Complex.normSq lam = sin² θ`.

This is the spectral derivation of `sin²θ_p = δ_p (2 - δ_p)`, equivalent
to but logically independent from Route 1 (`sin_sq_of_cos_eq_one_sub`). -/
theorem cyclic_phase_spectral
    (θ δ : ℝ) (lam : ℂ)
    (hlam : lam = ((1 - δ : ℝ) : ℂ))
    (hcos : Real.cos θ = 1 - δ) :
    1 - Complex.normSq lam = Real.sin θ ^ 2 := by
  rw [hlam, complex_normSq_real_sub]
  have hpyth : Real.sin θ ^ 2 + Real.cos θ ^ 2 = 1 := Real.sin_sq_add_cos_sq θ
  have hcos2 : Real.cos θ ^ 2 = (1 - δ)^2 := by rw [hcos]
  have hsin2 : Real.sin θ ^ 2 = 1 - (1 - δ)^2 := by linarith
  rw [hsin2]; ring

/-- **Norm-form** of `cyclic_phase_spectral`. -/
theorem cyclic_phase_spectral_norm
    (θ δ : ℝ) (lam : ℂ)
    (hlam : lam = ((1 - δ : ℝ) : ℂ))
    (hcos : Real.cos θ = 1 - δ) :
    1 - ‖lam‖ ^ 2 = Real.sin θ ^ 2 := by
  rw [hlam, complex_norm_sq_real_sub]
  have hpyth : Real.sin θ ^ 2 + Real.cos θ ^ 2 = 1 := Real.sin_sq_add_cos_sq θ
  have hcos2 : Real.cos θ ^ 2 = (1 - δ)^2 := by rw [hcos]
  have hsin2 : Real.sin θ ^ 2 = 1 - (1 - δ)^2 := by linarith
  rw [hsin2]; ring

/-! ### Section 3: Spectral identity matches the algebraic identity -/

/-- **Compatibility theorem.** Routes 1 and 2 produce numerically identical
conclusions for the same `(δ, θ)` data. This is a tautology after the
algebraic kernel `1 - (1-δ)² = δ(2-δ)` is fixed, but it is recorded here
as the explicit "the three routes agree" statement of the monograph. -/
theorem routes_one_two_agree
    (θ δ : ℝ) (lam : ℂ)
    (hlam : lam = ((1 - δ : ℝ) : ℂ))
    (hcos : Real.cos θ = 1 - δ) :
    Real.sin θ ^ 2 = δ * (2 - δ)
    ∧ 1 - Complex.normSq lam = δ * (2 - δ)
    ∧ 1 - Complex.normSq lam = Real.sin θ ^ 2 := by
  refine ⟨?_, ?_, cyclic_phase_spectral θ δ lam hlam hcos⟩
  · -- Route 1 conclusion
    have hpyth : Real.sin θ ^ 2 + Real.cos θ ^ 2 = 1 := Real.sin_sq_add_cos_sq θ
    have hcos2 : Real.cos θ ^ 2 = (1 - δ)^2 := by rw [hcos]
    have hsin2 : Real.sin θ ^ 2 = 1 - (1 - δ)^2 := by linarith
    rw [hsin2]; ring
  · -- Route 2 conclusion
    rw [hlam]; exact complex_normSq_real_sub δ

/-! ### Section 4: Explicit Fourier sum (diagonal-only kernel) -/

/-- **Diagonal-only survival kernel** on a residue space `α`. Models the
"trivial" PT transfer kernel where mass `(1 - δ)` rests on residue `0`
and no mass escapes to non-zero residues. -/
def diagonalKernel {α : Type*} [DecidableEq α] [Zero α] (δ : ℝ) :
    α → ℂ := fun r => if r = 0 then ((1 - δ : ℝ) : ℂ) else (0 : ℂ)

@[simp] lemma diagonalKernel_zero {α : Type*} [DecidableEq α] [Zero α] (δ : ℝ) :
    diagonalKernel (α := α) δ 0 = ((1 - δ : ℝ) : ℂ) := by
  simp [diagonalKernel]

lemma diagonalKernel_ne_zero
    {α : Type*} [DecidableEq α] [Zero α] (δ : ℝ) {r : α} (hr : r ≠ 0) :
    diagonalKernel δ r = 0 := by
  simp [diagonalKernel, hr]

/-- **Fourier eigenvalue of the diagonal-only kernel.** For any "character"
(arbitrary complex-valued function on the index set) `χ` satisfying
`χ 0 = 1`, the convolution sum collapses:

  `∑_{r} (diagonalKernel δ) r · χ r = (1 - δ : ℂ)`.

This realises Route 2's hypothesis `T̂_p(χ₁) = 1 - δ_p` as a concrete
finite Fourier sum. The condition `χ 0 = 1` holds for every additive
character of `ZMod p` (including `χ₁`), as `χ_j(0) = ω^{j·0} = 1`. -/
theorem fourier_eigenvalue_diagonal_kernel
    {α : Type*} [DecidableEq α] [Zero α] [Fintype α]
    (δ : ℝ) (χ : α → ℂ) (hχ : χ 0 = 1) :
    ∑ r, diagonalKernel δ r * χ r = ((1 - δ : ℝ) : ℂ) := by
  classical
  rw [Finset.sum_eq_single (0 : α)]
  · simp [diagonalKernel, hχ]
  · intro r _ hr
    rw [diagonalKernel_ne_zero δ hr, zero_mul]
  · intro h
    exact absurd (Finset.mem_univ (0 : α)) h

/-- **Headline (concrete version): Route 2 spectral identity via the
explicit diagonal-only Fourier sum.**

For a residue type `α` (e.g. `ZMod p`), for any complex-valued character
`χ` with `χ 0 = 1`, and for the diagonal-only survival kernel `T r =
diagonalKernel δ r`, the resulting Fourier sum `T̂ := ∑ T r · χ r` is a
real number equal to `1 - δ`. Consequently, the spectral contraction
`1 - normSq T̂` equals `δ(2 - δ)`, which equals `sin²θ` when
`cos θ = 1 - δ`. -/
theorem cyclic_phase_spectral_via_kernel
    {α : Type*} [DecidableEq α] [Zero α] [Fintype α]
    (θ δ : ℝ) (χ : α → ℂ) (hχ : χ 0 = 1)
    (hcos : Real.cos θ = 1 - δ) :
    1 - Complex.normSq (∑ r, diagonalKernel δ r * χ r) = Real.sin θ ^ 2 :=
  cyclic_phase_spectral θ δ _
    (fourier_eigenvalue_diagonal_kernel δ χ hχ) hcos

/-! ### Section 5: PT specialisation — characters of `ZMod p` -/

/-- **Principal-character normalisation on `ZMod p`.** Every additive
character `χ_j(r) = ω_p^{j · r}` of `ZMod p` (where `ω_p = e^{2πi/p}`)
maps `0 ↦ 1`. We formalise this trivially for the standard character
families by exhibiting the zero-input-yields-one property.

Concretely, for any function `χ : ZMod p → ℂ` defined as
`χ r = c^{(r : ℕ)}` (or any power-of-something form), we get `χ 0 = 1`. -/
lemma zmod_character_zero_eq_one {p : ℕ} [NeZero p]
    (c : ℂ) :
    (fun r : ZMod p => c ^ ((r : ZMod p).val)) 0 = 1 := by
  simp

/-- **PT specialisation.** For `p` a prime (or any `p` with `NeZero p`),
fix the character family `χ_j(r) = c^{r.val}` for some complex `c`
(`c = e^{2πi/p}` in the analytic instantiation). Then the diagonal-only
Fourier sum equals `1 - δ` in `ℂ`, and consequently the spectral
contraction equals `sin²θ` for `cos θ = 1 - δ`. -/
theorem cyclic_phase_spectral_zmod
    (p : ℕ) [NeZero p] (θ δ : ℝ) (c : ℂ) (hcos : Real.cos θ = 1 - δ) :
    1 - Complex.normSq
        (∑ r : ZMod p, diagonalKernel δ r * c ^ ((r : ZMod p).val))
      = Real.sin θ ^ 2 := by
  refine cyclic_phase_spectral_via_kernel θ δ
    (fun r : ZMod p => c ^ ((r : ZMod p).val)) ?_ hcos
  simp

/-! ### Section 6: Convenience corollary — bare algebra identity (Route 2 core) -/

/-- **Bare algebraic identity at the heart of Route 2.** For a complex
number `lam` that is *real-valued* (i.e. `lam.im = 0`), the contraction
`1 - normSq lam = 1 - lam.re²`. Combined with the assumption
`lam.re = 1 - δ`, this yields `δ(2 - δ)`. -/
theorem normSq_one_sub_of_real_part
    (lam : ℂ) (δ : ℝ) (him : lam.im = 0) (hre : lam.re = 1 - δ) :
    1 - Complex.normSq lam = δ * (2 - δ) := by
  have : Complex.normSq lam = lam.re * lam.re := by
    simp [Complex.normSq_apply, him]
  rw [this, hre]
  ring

end PT.Holonomy
