/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Stochastic.T2T3KroneckerEigenvalues
import PT.Stochastic.T2T3StationaryUniqueness
import Mathlib.Tactic

/-!
# Uniqueness (up to scaling) of the Perron eigenvector of `T_2 ⊗ T_3`

This file establishes the **one-dimensionality of the Perron eigenspace**
of the Kronecker product `T_2 ⊗ T_3`, where `T_2 = (1)` is the trivial
`1 × 1` identity matrix on `(ℤ/2ℤ)*` and `T_3 = !![0,1;1,0]` is the
antidiagonal involution on `(ℤ/3ℤ)*`.

Since `T_2` acts trivially on the singleton factor `Fin 1` and `T_3`
acts as a swap on `Fin 2`, the eigenvalue equation
`(T_2 ⊗ T_3).mulVec v = (1 : ℝ) • v` reduces to two scalar equations

  `v (0, 1) = v (0, 0)`  and  `v (0, 0) = v (0, 1)`

(the second is redundant). Combined, they force `v (0, 0) = v (0, 1)`,
i.e. every Perron eigenvector is constant on the two-element state space
`{(0, 0), (0, 1)}`. Equivalently, every Perron eigenvector is a scalar
multiple of `T2T3_perronVec = (1, 1)`.

## Main results

* `IsPerronEigenvector v` — `v` satisfies `(T_2 ⊗ T_3).mulVec v = 1 • v`.
* `T2T3_perronVec_isPerronEigenvector` — `T2T3_perronVec` satisfies it.
* `T2T3_perron_symmetry` — for every `v` with `IsPerronEigenvector v`,
  `v (0, 0) = v (0, 1)`. **Algebraic characterisation**.
* `T2T3_perron_uniqueness` — for every `v` with `IsPerronEigenvector v`,
  `v = v (0, 0) • T2T3_perronVec`. **Uniqueness up to scaling**.
* `T2T3_perron_uniqueness_nonzero` — variant: every nonzero Perron
  eigenvector is a (nonzero) scalar multiple of `T2T3_perronVec`.
* `T2T3_perron_eigenspace_one_dim` — characterisation: the Perron
  eigenspace equals `{c • T2T3_perronVec : c ∈ ℝ}`.
* `T2T3_perron_positive_unique` — **Perron-Frobenius corollary**: every
  componentwise nonnegative Perron eigenvector is a nonnegative multiple
  of `T2T3_perronVec`; if moreover it is nonzero, the multiplier is
  strictly positive.

## Reference

Monograph Chapter 7 §"Unicité du vecteur de Perron de `T_2 ⊗ T_3`",
companion to `T2T3StationaryUniqueness.T2T3_unique_stationary`.
-/

namespace PT.Stochastic

open Matrix Kronecker BigOperators PT.Sieve

/-! ### Predicate: Perron eigenvector for `T_2 ⊗ T_3` -/

/-- **Perron eigenvector predicate.** `v : Fin 1 × Fin 2 → ℝ` is a Perron
    eigenvector of `T_2 ⊗ T_3` iff `(T_2 ⊗ T_3).mulVec v = (1 : ℝ) • v`.
    The zero vector trivially satisfies this; the non-trivial content is
    the existence/uniqueness of *non-zero* eigenvectors. -/
def IsPerronEigenvector (v : Fin 1 × Fin 2 → ℝ) : Prop :=
  (T2_trivial ⊗ₖ T3).mulVec v = (1 : ℝ) • v

/-- `T2T3_perronVec` is a Perron eigenvector. Direct restatement of
    `T2T3_perron_eigen'` in predicate form. -/
theorem T2T3_perronVec_isPerronEigenvector :
    IsPerronEigenvector T2T3_perronVec :=
  T2T3_perron_eigen'

/-! ### Explicit values of `T2T3_perronVec` -/

@[simp] lemma T2T3_perronVec_zero_zero : T2T3_perronVec (0, 0) = 1 := by
  unfold T2T3_perronVec vecTensor T2_perronVec
  simp [perronEigenvector]

@[simp] lemma T2T3_perronVec_zero_one : T2T3_perronVec (0, 1) = 1 := by
  unfold T2T3_perronVec vecTensor T2_perronVec
  simp [perronEigenvector]

/-! ### Algebraic characterisation: symmetry forced by the eigenvalue equation -/

/-- **Symmetry forced by `+1` eigenvalue.** If `v` is a Perron eigenvector
    of `T_2 ⊗ T_3`, then `v (0, 0) = v (0, 1)`.

    *Proof.* The mulVec equation evaluated at `(0, 0)` reads
    `v (0, 1) = v (0, 0)` (from `T2T3_mulVec_apply_zero` and the fact
    that `(1 : ℝ) • v` at `(0, 0)` equals `v (0, 0)`). -/
theorem T2T3_perron_symmetry (v : Fin 1 × Fin 2 → ℝ)
    (hv : IsPerronEigenvector v) : v (0, 0) = v (0, 1) := by
  -- Evaluate the eigenvalue equation at (0, 0).
  have h0 : ((T2_trivial ⊗ₖ T3).mulVec v) (0, 0) = ((1 : ℝ) • v) (0, 0) := by
    rw [hv]
  rw [T2T3_mulVec_apply_zero] at h0
  -- h0 : v (0, 1) = (1 • v) (0, 0)
  have : ((1 : ℝ) • v) (0, 0) = v (0, 0) := by
    simp
  rw [this] at h0
  -- h0 : v (0, 1) = v (0, 0)
  exact h0.symm

/-! ### Uniqueness up to scaling -/

/-- **Uniqueness up to scaling.** Every Perron eigenvector `v` of
    `T_2 ⊗ T_3` is `v (0, 0)` times the canonical Perron vector
    `T2T3_perronVec = (1, 1)`.

    Since `T2T3_perronVec (0, 0) = 1`, the scaling factor is exactly the
    `(0, 0)`-coordinate of `v`. -/
theorem T2T3_perron_uniqueness (v : Fin 1 × Fin 2 → ℝ)
    (hv : IsPerronEigenvector v) : v = (v (0, 0)) • T2T3_perronVec := by
  have hsym := T2T3_perron_symmetry v hv
  funext ij
  obtain ⟨i, j⟩ := ij
  fin_cases i
  fin_cases j
  · -- (0, 0)
    show v (0, 0) = ((v (0, 0)) • T2T3_perronVec) (0, 0)
    simp
  · -- (0, 1) — use the symmetry v(0,0) = v(0,1).
    show v (0, 1) = ((v (0, 0)) • T2T3_perronVec) (0, 1)
    rw [← hsym]
    simp

/-- **Uniqueness up to scaling (nonzero variant).** Every nonzero Perron
    eigenvector is a (necessarily nonzero) scalar multiple of
    `T2T3_perronVec`. -/
theorem T2T3_perron_uniqueness_nonzero (v : Fin 1 × Fin 2 → ℝ)
    (hv : IsPerronEigenvector v) (hne : v ≠ 0) :
    ∃ c : ℝ, c ≠ 0 ∧ v = c • T2T3_perronVec := by
  refine ⟨v (0, 0), ?_, T2T3_perron_uniqueness v hv⟩
  -- The scalar v(0,0) cannot be zero, else v = 0 (since v(0,1) = v(0,0)).
  intro hc
  apply hne
  have hrep := T2T3_perron_uniqueness v hv
  rw [hc, zero_smul] at hrep
  exact hrep

/-! ### Full characterisation of the Perron eigenspace -/

/-- **Perron eigenspace = line generated by `T2T3_perronVec`.**
    `v` is a Perron eigenvector iff it is some scalar multiple of
    `T2T3_perronVec`. -/
theorem T2T3_perron_eigenspace_one_dim (v : Fin 1 × Fin 2 → ℝ) :
    IsPerronEigenvector v ↔ ∃ c : ℝ, v = c • T2T3_perronVec := by
  constructor
  · intro hv
    exact ⟨v (0, 0), T2T3_perron_uniqueness v hv⟩
  · rintro ⟨c, rfl⟩
    -- Need (T_2 ⊗ T_3).mulVec (c • T2T3_perronVec) = 1 • (c • T2T3_perronVec).
    unfold IsPerronEigenvector
    rw [Matrix.mulVec_smul]
    rw [T2T3_perron_eigen']
    -- Goal: c • (1 • T2T3_perronVec) = 1 • (c • T2T3_perronVec)
    rw [smul_comm]

/-! ### Perron-Frobenius corollary: positivity -/

/-- **Componentwise-nonneg Perron eigenvectors are nonneg multiples
    of `T2T3_perronVec`.** If `v` is a Perron eigenvector with
    `v ij ≥ 0` for all `ij`, then `v (0, 0) ≥ 0` and `v = v (0, 0) • T2T3_perronVec`. -/
theorem T2T3_perron_nonneg_multiple (v : Fin 1 × Fin 2 → ℝ)
    (hv : IsPerronEigenvector v) (hnn : ∀ ij, 0 ≤ v ij) :
    0 ≤ v (0, 0) ∧ v = (v (0, 0)) • T2T3_perronVec :=
  ⟨hnn (0, 0), T2T3_perron_uniqueness v hv⟩

/-- **Positive Perron-Frobenius corollary.** Every componentwise-nonneg
    nonzero Perron eigenvector is a *strictly positive* multiple of
    `T2T3_perronVec`. -/
theorem T2T3_perron_positive_unique (v : Fin 1 × Fin 2 → ℝ)
    (hv : IsPerronEigenvector v) (hnn : ∀ ij, 0 ≤ v ij) (hne : v ≠ 0) :
    ∃ c : ℝ, 0 < c ∧ v = c • T2T3_perronVec := by
  refine ⟨v (0, 0), ?_, T2T3_perron_uniqueness v hv⟩
  -- v(0,0) ≥ 0 by hnn; we need strict positivity. Suppose v(0,0) = 0;
  -- then by symmetry v(0,1) = 0, so v ≡ 0, contradicting hne.
  have hnn0 : 0 ≤ v (0, 0) := hnn (0, 0)
  rcases lt_or_eq_of_le hnn0 with hpos | hzero
  · exact hpos
  · exfalso
    apply hne
    have hsym := T2T3_perron_symmetry v hv
    funext ij
    obtain ⟨i, j⟩ := ij
    fin_cases i
    fin_cases j
    · show v (0, 0) = (0 : Fin 1 × Fin 2 → ℝ) (0, 0)
      simp [← hzero]
    · show v (0, 1) = (0 : Fin 1 × Fin 2 → ℝ) (0, 1)
      simp [← hsym, ← hzero]

/-! ### Headline summary -/

/-- **Headline (uniqueness of the Perron eigenvector of `T_2 ⊗ T_3`).**

    1. `T2T3_perronVec` is a Perron eigenvector.
    2. Every Perron eigenvector `v` satisfies the symmetry `v (0, 0) = v (0, 1)`.
    3. Every Perron eigenvector is a scalar multiple of `T2T3_perronVec`,
       with explicit scalar `v (0, 0)`.
    4. Every componentwise-nonneg, nonzero Perron eigenvector is a
       strictly positive multiple of `T2T3_perronVec`.

    Together, these show that the Perron eigenspace of `T_2 ⊗ T_3` is
    `1`-dimensional, and that its positive cone reduces to a single ray
    `ℝ_{>0} · T2T3_perronVec`. This is the Perron-Frobenius statement
    in the trivial-Kronecker case, and it complements
    `T2T3_unique_stationary` (uniqueness at the level of probability
    distributions) by giving uniqueness at the level of eigenvectors. -/
theorem T2T3_perron_uniqueness_summary :
    IsPerronEigenvector T2T3_perronVec
    ∧ (∀ v : Fin 1 × Fin 2 → ℝ, IsPerronEigenvector v →
         v (0, 0) = v (0, 1))
    ∧ (∀ v : Fin 1 × Fin 2 → ℝ, IsPerronEigenvector v →
         v = (v (0, 0)) • T2T3_perronVec)
    ∧ (∀ v : Fin 1 × Fin 2 → ℝ, IsPerronEigenvector v →
         (∀ ij, 0 ≤ v ij) → v ≠ 0 →
         ∃ c : ℝ, 0 < c ∧ v = c • T2T3_perronVec) :=
  ⟨T2T3_perronVec_isPerronEigenvector,
   T2T3_perron_symmetry,
   T2T3_perron_uniqueness,
   T2T3_perron_positive_unique⟩

end PT.Stochastic
