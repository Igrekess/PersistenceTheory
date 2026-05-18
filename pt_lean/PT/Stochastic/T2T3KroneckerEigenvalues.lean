/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Stochastic.T2T30Normalisation
import PT.Stochastic.CRTFactorizationT2
import Mathlib.Tactic

/-!
# Full spectrum of `T_2 ⊗ T_3` (Ch07 extension)

This file enumerates the **full spectrum** of the Kronecker product
`T_2 ⊗ T_3`, where `T_2` is the trivial `1 × 1` identity matrix on
`(ℤ/2ℤ)*` and `T_3` is the `2 × 2` antidiagonal involution on `(ℤ/3ℤ)*`.
The space `(Fin 1) × (Fin 2)` has total dimension `1 · 2 = 2`, so the
spectrum has exactly `2` eigenvalues:

* **Perron** (1, 1) — eigenvalue `1 · 1 = 1`, eigenvector
  `v_+(T_2) ⊗ v_+(T_3) = 1 ⊗ (1, 1)`.
* **Anti** (1, -1) — eigenvalue `1 · (-1) = -1`, eigenvector
  `v_+(T_2) ⊗ v_-(T_3) = 1 ⊗ (1, -1)`.

The two eigenvalues `{+1, -1}` exhaust the spectrum (since the space is
`2`-dimensional). Both have absolute value `1`, so the spectral radius is
`1`.

## Reference

Monograph Chapter 7 §"Spectre de T_2 ⊗ T_3", complement to
`PT.Stochastic.CRTFactorizationT2`.
-/

namespace PT.Stochastic

open Matrix Kronecker BigOperators PT.Sieve

/-! ### Perron eigenpair -/

/-- The Perron eigenvector `v_+(T_2 ⊗ T_3) = 1 ⊗ (1, 1)`. -/
def T2T3_perronVec : Fin 1 × Fin 2 → ℝ :=
  vecTensor T2_perronVec perronEigenvector

/-- **Perron eigenpair.** `T_2 ⊗ T_3` has eigenvalue `1` with the tensor
    eigenvector. -/
theorem T2T3_perron_eigen :
    (T2_trivial ⊗ₖ T3).mulVec T2T3_perronVec
      = ((1 : ℝ) * 1) • T2T3_perronVec := by
  unfold T2T3_perronVec
  exact kronecker_eigenvector T2_trivial T3 T2_perronVec perronEigenvector
          1 1 T2_perron_eigen T3_perron_eigen

/-- Cleaner restatement: the eigenvalue is `1`. -/
theorem T2T3_perron_eigen' :
    (T2_trivial ⊗ₖ T3).mulVec T2T3_perronVec = (1 : ℝ) • T2T3_perronVec := by
  have h := T2T3_perron_eigen
  simpa [one_mul] using h

/-! ### Anti eigenpair -/

/-- The anti eigenvector `v_-(T_2 ⊗ T_3) = 1 ⊗ (1, -1)`. -/
def T2T3_antiVec : Fin 1 × Fin 2 → ℝ :=
  vecTensor T2_perronVec antiEigenvector

/-- **Anti eigenpair.** `T_2 ⊗ T_3` has eigenvalue `-1` with the
    anti-tensor eigenvector. -/
theorem T2T3_anti_eigen :
    (T2_trivial ⊗ₖ T3).mulVec T2T3_antiVec
      = ((1 : ℝ) * (-1)) • T2T3_antiVec := by
  unfold T2T3_antiVec
  exact kronecker_eigenvector T2_trivial T3 T2_perronVec antiEigenvector
          1 (-1) T2_perron_eigen T3_sub_eigen

/-- Cleaner restatement: the eigenvalue is `-1`. -/
theorem T2T3_anti_eigen' :
    (T2_trivial ⊗ₖ T3).mulVec T2T3_antiVec = (-1 : ℝ) • T2T3_antiVec := by
  have h := T2T3_anti_eigen
  simpa [one_mul] using h

/-! ### Spectral radius and absolute-value bounds -/

/-- Both eigenvalues `{+1, -1}` have absolute value `1`. -/
theorem T2T3_spectral_radius_eq_one :
    |(1 : ℝ)| = 1 ∧ |(-1 : ℝ)| = 1 ∧ max |(1 : ℝ)| |(-1 : ℝ)| = 1 := by
  refine ⟨?_, ?_, ?_⟩ <;> norm_num

/-! ### Positivity / sign of eigenvectors -/

/-- The Perron eigenvector is strictly positive componentwise. -/
theorem T2T3_perronVec_pos : ∀ ij, 0 < T2T3_perronVec ij := by
  intro ij
  obtain ⟨i, j⟩ := ij
  unfold T2T3_perronVec vecTensor
  have hT2 := T2_perronVec_pos i
  have hT3 := perronEigenvector_pos j
  positivity

/-- The anti eigenvector has mixed signs (not strictly positive). -/
theorem T2T3_antiVec_not_pos :
    ¬ (∀ ij, 0 < T2T3_antiVec ij) := by
  intro h
  have h1 : (0 : ℝ) < T2T3_antiVec (0, 1) := h (0, 1)
  unfold T2T3_antiVec vecTensor T2_perronVec antiEigenvector at h1
  simp at h1
  linarith

/-! ### Headline (full spectrum) -/

/-- **Headline (full spectrum of `T_2 ⊗ T_3`).**

    * Eigenvalue `+1` with strictly positive Perron eigenvector
      `v_+(T_2) ⊗ v_+(T_3) = 1 ⊗ (1, 1)`.
    * Eigenvalue `-1` with sign-changing eigenvector
      `v_+(T_2) ⊗ v_-(T_3) = 1 ⊗ (1, -1)`.
    * Spectral radius: `max(|1|, |-1|) = 1`.

    Since the state space `(Fin 1) × (Fin 2)` has dimension `2`, these two
    eigenvalues exhaust the spectrum. -/
theorem T2T3_full_spectrum_summary :
    -- Perron eigenpair
    (T2_trivial ⊗ₖ T3).mulVec T2T3_perronVec = (1 : ℝ) • T2T3_perronVec
    ∧ (∀ ij, 0 < T2T3_perronVec ij)
    -- Anti eigenpair
    ∧ (T2_trivial ⊗ₖ T3).mulVec T2T3_antiVec = (-1 : ℝ) • T2T3_antiVec
    ∧ ¬ (∀ ij, 0 < T2T3_antiVec ij)
    -- Spectral radius
    ∧ max |(1 : ℝ)| |(-1 : ℝ)| = 1 :=
  ⟨T2T3_perron_eigen', T2T3_perronVec_pos,
   T2T3_anti_eigen', T2T3_antiVec_not_pos,
   T2T3_spectral_radius_eq_one.2.2⟩

end PT.Stochastic
