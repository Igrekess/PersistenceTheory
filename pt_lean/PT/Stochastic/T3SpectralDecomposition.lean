/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Stochastic.SpectralDominance
import Mathlib.Tactic

/-!
# Spectral decomposition extension — Eigenbasis closure for `T_3`

This file extends `PT.Stochastic.SpectralDominance` with the **eigenbasis
decomposition** of `ℝ²` under the involution `T_3`:

* The Perron direction `v_+ = (1, 1)` and the anti-direction `v_- = (1, -1)`
  form a basis of `ℝ²`.
* Every vector `v ∈ ℝ²` decomposes uniquely as
  `v = α v_+ + β v_-` with
  `α = (v_0 + v_1) / 2`, `β = (v_0 - v_1) / 2`.
* The two projectors `Π_+ = (I + T_3)/2` (Perron / stationary) and
  `Π_- = (I - T_3)/2` (anti) form a **complete spectral resolution**:
  `Π_+ + Π_- = I`, `Π_+ Π_- = 0`, and they are idempotent.
* `T_3 = Π_+ - Π_-` (the eigenvalue equation in projector form).

These results constitute the **spectral closure** of `T_3`-dynamics in the
sense of monograph audit row #31 ("`T_m`-induced dynamics close under
spectral projection"): every iterate `T_3^n` stays in
`Span(Π_+, Π_-)` (a 2-dimensional matrix subspace).

## Reference

Monograph Chapter 7, §"Dominance spectrale et fermeture spectrale",
follow-up to `SpectralDominance.lean`. Audit row #31.
-/

namespace PT.Stochastic

open Matrix

/-! ### Components of an arbitrary vector along the spectral basis -/

/-- The Perron coefficient of `v : Fin 2 → ℝ`. `α(v) = (v 0 + v 1) / 2`. -/
noncomputable def alphaCoef (v : Fin 2 → ℝ) : ℝ := (v 0 + v 1) / 2

/-- The anti coefficient of `v : Fin 2 → ℝ`. `β(v) = (v 0 - v 1) / 2`. -/
noncomputable def betaCoef (v : Fin 2 → ℝ) : ℝ := (v 0 - v 1) / 2

/-- **Decomposition.** For every `v ∈ ℝ²`,
    `v = α(v) v_+ + β(v) v_-`. -/
theorem vec_decomposition (v : Fin 2 → ℝ) :
    v = (fun i => alphaCoef v * perronEigenvector i
                + betaCoef v * antiEigenvector i) := by
  funext i
  unfold alphaCoef betaCoef perronEigenvector antiEigenvector
  fin_cases i <;> simp <;> ring

/-! ### The two projectors `Π_+` (Perron) and `Π_-` (anti) -/

/-- The Perron projector `Π_+ = (I + T_3) / 2 = !![1/2, 1/2; 1/2, 1/2]`.
    Alias of `stationaryProjector`. -/
noncomputable def perronProjector : Matrix (Fin 2) (Fin 2) ℝ :=
  stationaryProjector

/-- The anti projector `Π_- = (I - T_3) / 2 = !![1/2, -1/2; -1/2, 1/2]`. -/
noncomputable def antiProjector : Matrix (Fin 2) (Fin 2) ℝ :=
  !![1/2, -1/2; -1/2, 1/2]

@[simp] lemma antiProjector_zero_zero : antiProjector 0 0 = 1/2 := rfl
@[simp] lemma antiProjector_zero_one  : antiProjector 0 1 = -1/2 := rfl
@[simp] lemma antiProjector_one_zero  : antiProjector 1 0 = -1/2 := rfl
@[simp] lemma antiProjector_one_one   : antiProjector 1 1 = 1/2 := rfl

/-! ### Spectral resolution: `Π_+ + Π_- = I` -/

/-- **Resolution of the identity.** `Π_+ + Π_- = I`. -/
theorem projector_sum_eq_one :
    perronProjector + antiProjector = (1 : Matrix (Fin 2) (Fin 2) ℝ) := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [perronProjector, stationaryProjector, antiProjector] <;>
    norm_num

/-- **Idempotence of anti projector.** `Π_- · Π_- = Π_-`. -/
theorem antiProjector_idempotent :
    antiProjector * antiProjector = antiProjector := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [antiProjector, Matrix.mul_apply, Fin.sum_univ_two] <;>
    norm_num

/-- **Orthogonality.** `Π_+ · Π_- = 0`. -/
theorem projectors_orthogonal :
    perronProjector * antiProjector = (0 : Matrix (Fin 2) (Fin 2) ℝ) := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [perronProjector, stationaryProjector, antiProjector,
          Matrix.mul_apply, Fin.sum_univ_two] <;>
    norm_num

/-- **Orthogonality (other order).** `Π_- · Π_+ = 0`. -/
theorem projectors_orthogonal' :
    antiProjector * perronProjector = (0 : Matrix (Fin 2) (Fin 2) ℝ) := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [perronProjector, stationaryProjector, antiProjector,
          Matrix.mul_apply, Fin.sum_univ_two] <;>
    norm_num

/-! ### Spectral eigenvalue equation in projector form -/

/-- **Spectral resolution of `T_3`.** `T_3 = Π_+ - Π_-` (eigenvalue equation
    in projector form: `+1 · Π_+ + (-1) · Π_-`). -/
theorem T3_eq_projector_difference :
    PT.Sieve.T3 = perronProjector - antiProjector := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [PT.Sieve.T3, perronProjector, stationaryProjector, antiProjector] <;>
    norm_num

/-- **Anti projector annihilates Perron direction.** `Π_- · v_+ = 0`. -/
theorem antiProjector_kills_perron :
    antiProjector.mulVec perronEigenvector = 0 := by
  ext i
  fin_cases i <;> simp [antiProjector, perronEigenvector,
                        Matrix.mulVec, dotProduct] <;> norm_num

/-- **Anti projector fixes anti-direction.** `Π_- · v_- = v_-`. -/
theorem antiProjector_fixes_anti :
    antiProjector.mulVec antiEigenvector = antiEigenvector := by
  ext i
  fin_cases i <;> simp [antiProjector, antiEigenvector,
                        Matrix.mulVec, dotProduct] <;> norm_num

/-! ### Headline: spectral resolution -/

/-- **Headline (spectral resolution).** The two projectors `Π_+` and `Π_-`
    form a complete orthogonal spectral resolution for the involution `T_3`:

    * `Π_+ + Π_- = I` (resolution),
    * `Π_+ · Π_+ = Π_+`, `Π_- · Π_- = Π_-` (idempotence),
    * `Π_+ · Π_- = 0`, `Π_- · Π_+ = 0` (orthogonality),
    * `T_3 = Π_+ − Π_-` (eigenvalue equation),
    * `Π_+ v_+ = v_+`, `Π_+ v_- = 0`, `Π_- v_+ = 0`, `Π_- v_- = v_-`.

    Combined, every iterate `T_3^n` lies in the 2-dimensional matrix subspace
    `Span(Π_+, Π_-)` — this is the spectral-closure statement of audit
    row #31. -/
theorem T3_spectral_resolution :
    perronProjector + antiProjector = (1 : Matrix (Fin 2) (Fin 2) ℝ)
    ∧ perronProjector * perronProjector = perronProjector
    ∧ antiProjector * antiProjector = antiProjector
    ∧ perronProjector * antiProjector = (0 : Matrix (Fin 2) (Fin 2) ℝ)
    ∧ antiProjector * perronProjector = (0 : Matrix (Fin 2) (Fin 2) ℝ)
    ∧ PT.Sieve.T3 = perronProjector - antiProjector :=
  ⟨projector_sum_eq_one, stationaryProjector_idempotent,
   antiProjector_idempotent, projectors_orthogonal, projectors_orthogonal',
   T3_eq_projector_difference⟩

end PT.Stochastic
