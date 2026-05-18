/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Stochastic.T2T3KroneckerEigenvalues
import PT.Stochastic.SHalf
import Mathlib.Tactic

/-!
# Unicité de la distribution stationnaire pour `T_2 ⊗ T_3`

Ce fichier démontre que la distribution stationnaire de la matrice de
Kronecker `T_2 ⊗ T_3` (où `T_2 = (1)` est le bloc trivial `1 × 1` de
`T2T30Normalisation` et `T_3 = !![0,1;1,0]` est l'involution antidiagonale
de `T3Antidiagonal`) est **unique**, égale à

  `π_{T_2 ⊗ T_3} (0, 0) = 1/2,   π_{T_2 ⊗ T_3} (0, 1) = 1/2`,

et coïncide, modulo le facteur trivial `T_2`, avec la distribution
stationnaire `(1/2, 1/2)` de `T_3` (cf. `SHalf.T3_unique_stationary`).

La preuve est élémentaire car l'espace d'états `Fin 1 × Fin 2` est de
dimension `2` et `T_2 = 1` agit trivialement sur le facteur `Fin 1` :
le système se réduit donc à `T_3` sur la composante `Fin 2`.

## Résultats principaux

* `T2T3_stationary` — la distribution candidate `(1/2, 1/2)`.
* `IsT2T3Stationary` — la propriété de distribution stationnaire pour
  `T_2 ⊗ T_3` : `mulVec`-invariance + somme `= 1` + positivité.
* `T2T3_stationary_isStationary` — `(1/2, 1/2)` est bien une distribution
  stationnaire.
* `T2T3_stationary_sum_one` et `T2T3_stationary_nonneg` — propriétés
  élémentaires (somme = 1 et positivité).
* `T2T3_unique_stationary` — **unicité** : toute distribution stationnaire
  pour `T_2 ⊗ T_3` est égale à `T2T3_stationary`.

## Stratégie

La structure de `T_2 ⊗ T_3` comme `Matrix (Fin 1 × Fin 2) (Fin 1 × Fin 2) ℝ`
donne, après évaluation, deux équations sur les deux composantes :

* `π (0, 0) = π (0, 1)` (la mulVec-invariance en `(0, 0)`),
* `π (0, 1) = π (0, 0)` (la mulVec-invariance en `(0, 1)` — redondante).

Combinées avec `π (0, 0) + π (0, 1) = 1`, elles forcent `π = (1/2, 1/2)`.

## Référence

Monographe ch07 §"Distribution stationnaire de `T_2 ⊗ T_3`", suite à
`PT.Stochastic.SHalf.T3_unique_stationary`.
-/

namespace PT.Stochastic

open Matrix Kronecker BigOperators PT.Sieve

/-! ### Distribution stationnaire candidate -/

/-- **Distribution stationnaire candidate de `T_2 ⊗ T_3`.**
    `π (0, 0) = 1/2, π (0, 1) = 1/2`. -/
noncomputable def T2T3_stationary : Fin 1 × Fin 2 → ℝ :=
  fun _ => 1 / 2

@[simp] lemma T2T3_stationary_zero_zero : T2T3_stationary (0, 0) = 1 / 2 := rfl
@[simp] lemma T2T3_stationary_zero_one  : T2T3_stationary (0, 1) = 1 / 2 := rfl

/-! ### Propriété de distribution stationnaire -/

/-- **Distribution stationnaire pour `T_2 ⊗ T_3`.** Une distribution
    `π : Fin 1 × Fin 2 → ℝ` est *stationnaire* pour `T_2 ⊗ T_3` si :

    * elle est laissée fixe par l'action de mulVec :
      `(T_2 ⊗ T_3).mulVec π = π` ;
    * elle est une distribution de probabilité : `∑ π = 1` et `π ≥ 0`. -/
structure IsT2T3Stationary (π : Fin 1 × Fin 2 → ℝ) : Prop where
  fixed    : (T2_trivial ⊗ₖ T3).mulVec π = π
  prob_sum : ∑ ij, π ij = 1
  nonneg   : ∀ ij, 0 ≤ π ij

/-! ### Entrées explicites du produit de Kronecker -/

/-- Évaluation entrée-par-entrée de `T_2 ⊗ T_3`. -/
@[simp] lemma T2T3_apply (i₁ k₁ : Fin 1) (j l : Fin 2) :
    (T2_trivial ⊗ₖ T3) (i₁, j) (k₁, l) = T2_trivial i₁ k₁ * T3 j l := by
  rfl

/-- Calcul explicite de `(T_2 ⊗ T_3).mulVec π` en `(0, 0)`. -/
lemma T2T3_mulVec_apply_zero (π : Fin 1 × Fin 2 → ℝ) :
    ((T2_trivial ⊗ₖ T3).mulVec π) (0, 0) = π (0, 1) := by
  simp [Matrix.mulVec, dotProduct, Fintype.sum_prod_type,
        Fin.sum_univ_two, T2_trivial, Matrix.one_apply_eq, T3]

/-- Calcul explicite de `(T_2 ⊗ T_3).mulVec π` en `(0, 1)`. -/
lemma T2T3_mulVec_apply_one (π : Fin 1 × Fin 2 → ℝ) :
    ((T2_trivial ⊗ₖ T3).mulVec π) (0, 1) = π (0, 0) := by
  simp [Matrix.mulVec, dotProduct, Fintype.sum_prod_type,
        Fin.sum_univ_two, T2_trivial, Matrix.one_apply_eq, T3]

/-! ### `(1/2, 1/2)` est stationnaire -/

/-- **Invariance par mulVec.** `(T_2 ⊗ T_3).mulVec T2T3_stationary = T2T3_stationary`. -/
theorem T2T3_stationary_fixed :
    (T2_trivial ⊗ₖ T3).mulVec T2T3_stationary = T2T3_stationary := by
  funext ij
  obtain ⟨i, j⟩ := ij
  fin_cases i
  fin_cases j
  · -- (0, 0)
    show ((T2_trivial ⊗ₖ T3).mulVec T2T3_stationary) (0, 0)
          = T2T3_stationary (0, 0)
    rw [T2T3_mulVec_apply_zero]
    rfl
  · -- (0, 1)
    show ((T2_trivial ⊗ₖ T3).mulVec T2T3_stationary) (0, 1)
          = T2T3_stationary (0, 1)
    rw [T2T3_mulVec_apply_one]
    rfl

/-- **Somme = 1.** `∑ ij, T2T3_stationary ij = 1`. -/
theorem T2T3_stationary_sum_one :
    ∑ ij, T2T3_stationary ij = (1 : ℝ) := by
  rw [Fintype.sum_prod_type, Fin.sum_univ_one, Fin.sum_univ_two]
  norm_num [T2T3_stationary]

/-- **Positivité.** `T2T3_stationary ij ≥ 0` pour tout `ij`. -/
theorem T2T3_stationary_nonneg : ∀ ij, 0 ≤ T2T3_stationary ij := by
  intro ij
  unfold T2T3_stationary
  norm_num

/-- **`(1/2, 1/2)` est une distribution stationnaire pour `T_2 ⊗ T_3`.** -/
theorem T2T3_stationary_isStationary : IsT2T3Stationary T2T3_stationary where
  fixed    := T2T3_stationary_fixed
  prob_sum := T2T3_stationary_sum_one
  nonneg   := T2T3_stationary_nonneg

/-! ### Unicité -/

/-- **Unicité de la distribution stationnaire.** Toute distribution
    stationnaire `π` pour `T_2 ⊗ T_3` (au sens de `IsT2T3Stationary`)
    est égale à `T2T3_stationary = (1/2, 1/2)`.

    *Stratégie.* Les deux équations mulVec-fixe donnent
    `π (0, 1) = π (0, 0)` (et l'autre est redondante). Combinée avec la
    contrainte de probabilité `π (0, 0) + π (0, 1) = 1`, on obtient
    `π (0, 0) = π (0, 1) = 1/2`. -/
theorem T2T3_unique_stationary (π : Fin 1 × Fin 2 → ℝ)
    (h : IsT2T3Stationary π) : π = T2T3_stationary := by
  -- Mulvec-fixity at (0, 0): π(0, 1) = π(0, 0)
  have hfix := h.fixed
  have hfix0 : ((T2_trivial ⊗ₖ T3).mulVec π) (0, 0) = π (0, 0) := by
    rw [hfix]
  rw [T2T3_mulVec_apply_zero] at hfix0
  -- hfix0 : π (0, 1) = π (0, 0)
  -- Probability sum: π (0, 0) + π (0, 1) = 1
  have hsum := h.prob_sum
  rw [Fintype.sum_prod_type, Fin.sum_univ_one, Fin.sum_univ_two] at hsum
  -- hsum : π (0, 0) + π (0, 1) = 1
  have hp0 : π (0, 0) = 1 / 2 := by linarith
  have hp1 : π (0, 1) = 1 / 2 := by linarith
  -- Conclude pointwise.
  funext ij
  obtain ⟨i, j⟩ := ij
  fin_cases i
  fin_cases j
  · show π (0, 0) = T2T3_stationary (0, 0); simpa using hp0
  · show π (0, 1) = T2T3_stationary (0, 1); simpa using hp1

/-! ### Réduction au cas `T_3` -/

/-- **Réduction au cas `T_3`.** Si `π : Fin 1 × Fin 2 → ℝ` est
    `T_2 ⊗ T_3`-stationnaire, alors la « tranche » `i ↦ π (0, i)` est
    `T_3`-stationnaire au sens de `IsStationary` de `SHalf`. -/
theorem T2T3_slice_isStationary (π : Fin 1 × Fin 2 → ℝ)
    (h : IsT2T3Stationary π) :
    IsStationary (fun i => π (0, i)) := by
  refine ⟨?_, ?_, ?_⟩
  · -- prob_sum : π(0,0) + π(0,1) = 1
    have hsum := h.prob_sum
    rw [Fintype.sum_prod_type, Fin.sum_univ_one, Fin.sum_univ_two] at hsum
    exact hsum
  · -- nonneg
    intro i
    exact h.nonneg (0, i)
  · -- fixed : ∑ j, π(0, j) * T3 j i = π (0, i)
    have hfix := h.fixed
    intro i
    fin_cases i
    · -- i = 0 : ∑ j, π(0, j) * T3 j 0 = π(0, 0)
      have hfix0 : ((T2_trivial ⊗ₖ T3).mulVec π) (0, 0) = π (0, 0) := by
        rw [hfix]
      rw [T2T3_mulVec_apply_zero] at hfix0
      -- hfix0 : π (0, 1) = π (0, 0)
      show ∑ j, π (0, j) * T3 j 0 = π (0, 0)
      simp [T3, Fin.sum_univ_two]
      linarith
    · -- i = 1 : ∑ j, π(0, j) * T3 j 1 = π(0, 1)
      have hfix1 : ((T2_trivial ⊗ₖ T3).mulVec π) (0, 1) = π (0, 1) := by
        rw [hfix]
      rw [T2T3_mulVec_apply_one] at hfix1
      -- hfix1 : π (0, 0) = π (0, 1)
      show ∑ j, π (0, j) * T3 j 1 = π (0, 1)
      simp [T3, Fin.sum_univ_two]
      linarith

/-- **Unicité via réduction au cas `T_3`.** Variante de
    `T2T3_unique_stationary` passant par la stationnarité de `T_3`
    (qui factorise via `SHalf.T3_unique_stationary`). -/
theorem T2T3_unique_stationary_via_T3 (π : Fin 1 × Fin 2 → ℝ)
    (h : IsT2T3Stationary π) : π = T2T3_stationary := by
  have hslice := T2T3_slice_isStationary π h
  have hπ : (fun i => π (0, i)) = piHalf := T3_unique_stationary _ hslice
  funext ij
  obtain ⟨i, j⟩ := ij
  fin_cases i
  fin_cases j
  · have hj : π (0, 0) = piHalf 0 := congrFun hπ 0
    show π (0, 0) = T2T3_stationary (0, 0)
    simpa [piHalf] using hj
  · have hj : π (0, 1) = piHalf 1 := congrFun hπ 1
    show π (0, 1) = T2T3_stationary (0, 1)
    simpa [piHalf] using hj

/-! ### Résumé -/

/-- **Résumé (unicité de la stationnaire de `T_2 ⊗ T_3`).**

    * `T2T3_stationary` est une distribution stationnaire :
      `(T_2 ⊗ T_3).mulVec T2T3_stationary = T2T3_stationary`,
      `∑ T2T3_stationary = 1`, `T2T3_stationary ≥ 0`.
    * Toute distribution stationnaire est égale à `T2T3_stationary`.
    * La « tranche » `i ↦ π(0, i)` d'une distribution stationnaire
      `π` est `T_3`-stationnaire, ce qui ramène l'unicité à
      `SHalf.T3_unique_stationary`. -/
theorem T2T3_stationary_summary :
    -- Existence
    IsT2T3Stationary T2T3_stationary
    -- Unicité
    ∧ (∀ π, IsT2T3Stationary π → π = T2T3_stationary)
    -- Valeurs explicites
    ∧ T2T3_stationary (0, 0) = 1 / 2
    ∧ T2T3_stationary (0, 1) = 1 / 2 :=
  ⟨T2T3_stationary_isStationary, T2T3_unique_stationary, rfl, rfl⟩

end PT.Stochastic
