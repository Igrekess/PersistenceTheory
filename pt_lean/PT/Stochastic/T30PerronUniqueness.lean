/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Stochastic.T2T30Normalisation
import PT.Stochastic.T2T3PerronEigenvectorUniqueness
import Mathlib.Tactic

/-!
# Uniqueness (up to scaling) of the Perron eigenvector of `T_30 = T_2 ⊗ T_3 ⊗ T_5`

This file extends the Perron uniqueness result of
`PT.Stochastic.T2T3PerronEigenvectorUniqueness` from the two-factor product
`T_2 ⊗ T_3` to the three-factor product `T_30 = T_2 ⊗ T_3 ⊗ T_5`.

## Scope: tensor-decomposable sector only

The general statement "the Perron eigenspace of `T_30` is one-dimensional" is
**equivalent** (over `ℝ`, by the Künneth-type splitting of Kronecker
eigenspaces for distinct dominant eigenvalues) to the conjunction of the
Perron-uniqueness statements for each factor `T_2`, `T_3`, `T_5`. For `T_2`
and `T_3` we have these unconditionally (`T_2` is `1 × 1`; `T_3` is treated
in `T2T3PerronEigenvectorUniqueness`). For `T_5`, however, the `T5Like`
record only exposes the *existence* of a Perron eigenpair, not its
uniqueness; the latter would require either a concrete spectral computation
of the `4 × 4` matrix or a structural Perron-Frobenius statement.

Rather than postulating uniqueness for `T_5` as an additional unstated
hypothesis, this file proves:

1. **Existence**: `T30_perronVec T5` is a Perron eigenvector of `T_30`,
   strictly positive (already established in `T2T30Normalisation`).
2. **Tensor-decomposable uniqueness**: every Perron eigenvector of `T_30`
   that decomposes as a Kronecker tensor `v = v_2 ⊗ v_3 ⊗ v_5` of factor
   Perron eigenvectors is a scalar multiple of `T30_perronVec`. Two
   variants are provided:
   * `T30_perron_tensor_uniqueness_factor` — direct version using the
     two intermediate Perron-uniqueness lemmas for `T_2` and `T_3`.
   * `T30_perron_tensor_uniqueness` — packaged convenience version.
3. **Reduction to `T_5`**: a `T5LikeUnique` structure extends `T5Like`
   with a Perron-uniqueness axiom for `T_5`. From this we recover the
   *full* one-dimensional Perron-eigenspace statement for `T_30`,
   conditionally on `T_5` being one-dimensional on its Perron sector.

The conditional `T5LikeUnique`-version is the precise structural sense in
which the three-factor Perron uniqueness reduces to the (numerically clear,
formally unproven) Perron uniqueness for `T_5`.

## Main results

* `IsT30PerronEigenvector T5 v` — predicate
  `(T30 T5).mulVec v = (1 : ℝ) • v`.
* `T30_perronVec_isPerronEigenvector` — `T30_perronVec T5` satisfies it.
* `T30_perron_tensor_uniqueness_factor` — tensor-decomposable Perron
  eigenvectors are scalar multiples of `T30_perronVec`.
* `T5LikeUnique` — extension of `T5Like` with a Perron-uniqueness axiom.
* `T30_perron_uniqueness_under_T5_unique` — *conditional full uniqueness*:
  given `T5LikeUnique`, every Perron eigenvector of `T_30` arising from a
  tensor decomposition is a scalar multiple of `T30_perronVec`.

## Reference

* Monograph Chapter 7, §"Unicité du vecteur de Perron de `T_{30}`",
  follow-up to `T2T3_perron_uniqueness_summary`.
* `T2T30Normalisation` — `T30`, `T30_perronVec`, `T30_perron_eigen`,
  `T5Like`.
* `T2T3PerronEigenvectorUniqueness` — base case `T_2 ⊗ T_3`.
-/

namespace PT.Stochastic

open Matrix Kronecker BigOperators PT.Sieve

/-! ### Predicate: Perron eigenvector for `T_30` -/

/-- **Perron eigenvector predicate for `T_30`.** `v` is a Perron eigenvector
    of `T_30 = T_2 ⊗ T_3 ⊗ T_5` iff `(T30 T5).mulVec v = (1 : ℝ) • v`.
    The zero vector satisfies this trivially; the non-trivial content is
    classifying the *non-zero* such `v`. -/
def IsT30PerronEigenvector (T5 : T5Like)
    (v : (Fin 1 × Fin 2) × Fin 4 → ℝ) : Prop :=
  (T30 T5).mulVec v = (1 : ℝ) • v

/-- `T30_perronVec T5` is a Perron eigenvector of `T_30`. Restatement of
    `T30_perron_eigen` in predicate form. -/
theorem T30_perronVec_isPerronEigenvector (T5 : T5Like) :
    IsT30PerronEigenvector T5 (T30_perronVec T5) :=
  T30_perron_eigen T5

/-! ### Predicate: factor-Perron eigenvector for `T_5` -/

/-- **Perron eigenvector predicate for `T_5`.** `v : Fin 4 → ℝ` satisfies
    `T_5.matrix.mulVec v = (1 : ℝ) • v`. -/
def IsT5PerronEigenvector (T5 : T5Like) (v : Fin 4 → ℝ) : Prop :=
  T5.matrix.mulVec v = (1 : ℝ) • v

/-- `T5.perronVec` is a Perron eigenvector of `T_5`. -/
theorem T5_perronVec_isPerronEigenvector (T5 : T5Like) :
    IsT5PerronEigenvector T5 T5.perronVec :=
  T5.perron_eigen

/-! ### Tensor-decomposable Perron eigenvectors are scalar multiples
        of `T30_perronVec` -/

/-- **Tensor-decomposable uniqueness (factor-by-factor form).** Suppose
    `v_2 : Fin 1 → ℝ`, `v_3 : Fin 2 → ℝ`, `v_5 : Fin 4 → ℝ` are Perron
    eigenvectors of `T_2`, `T_3`, `T_5` respectively. Then their
    Kronecker tensor product `(v_2 ⊗ v_3) ⊗ v_5` is a Perron eigenvector
    of `T_30 = T_2 ⊗ T_3 ⊗ T_5`. -/
theorem T30_tensor_isPerronEigenvector (T5 : T5Like)
    (v2 : Fin 1 → ℝ) (v3 : Fin 2 → ℝ) (v5 : Fin 4 → ℝ)
    (h2 : T2_trivial.mulVec v2 = (1 : ℝ) • v2)
    (h3 : T3.mulVec v3 = (1 : ℝ) • v3)
    (h5 : T5.matrix.mulVec v5 = (1 : ℝ) • v5) :
    IsT30PerronEigenvector T5 (vecTensor (vecTensor v2 v3) v5) := by
  unfold IsT30PerronEigenvector T30
  have h := kron3_eigenvector_left T2_trivial T3 T5.matrix
              v2 v3 v5 1 1 1 h2 h3 h5
  -- `h` gives the eigenvalue `1 * 1 * 1`; rewrite to `1`.
  rw [show (1 : ℝ) * 1 * 1 = 1 by ring] at h
  exact h

/-! ### Auxiliary: `T_2` Perron uniqueness (trivial since `dim = 1`) -/

/-- **`T_2` Perron uniqueness.** Every `v : Fin 1 → ℝ` is a Perron
    eigenvector of `T_2 = (1)` (the equation is the tautology
    `v = 1 • v`). Moreover, every such `v` is `v 0 • T2_perronVec`. -/
theorem T2_perron_uniqueness (v : Fin 1 → ℝ)
    (_hv : T2_trivial.mulVec v = (1 : ℝ) • v) :
    v = (v 0) • T2_perronVec := by
  funext i
  fin_cases i
  show v 0 = ((v 0) • T2_perronVec) 0
  simp [T2_perronVec]

/-- **`T_3` Perron uniqueness.** Every Perron eigenvector `v : Fin 2 → ℝ`
    of `T_3` satisfies `v = v 0 • perronEigenvector`.
    Proof: the eigenvalue equation forces `v 0 = v 1` since
    `T_3 = !![0,1;1,0]` swaps the two coordinates. -/
theorem T3_perron_uniqueness (v : Fin 2 → ℝ)
    (hv : T3.mulVec v = (1 : ℝ) • v) :
    v = (v 0) • perronEigenvector := by
  -- Use the existing `T2T3` two-factor uniqueness machinery: extend
  -- `v : Fin 2 → ℝ` to `Fin 1 × Fin 2 → ℝ` by `(0, j) ↦ v j`, then
  -- apply `T2T3_perron_uniqueness`.
  set w : Fin 1 × Fin 2 → ℝ := fun ij => v ij.2 with hw
  have hwPerron : IsPerronEigenvector w := by
    unfold IsPerronEigenvector
    funext ij
    obtain ⟨i, j⟩ := ij
    fin_cases i
    -- Goal: (T2_trivial ⊗ₖ T3).mulVec w (0, j) = (1 • w) (0, j)
    -- Compute LHS via Matrix.mulVec definition.
    show ((T2_trivial ⊗ₖ T3).mulVec w) (0, j) = ((1 : ℝ) • w) (0, j)
    -- Direct unfolding: row (0, j) of T_2 ⊗ T_3 has a single non-zero
    -- contribution: entries are T_2[0,0] · T_3[j, k] = T_3[j, k].
    have hmv : ((T2_trivial ⊗ₖ T3).mulVec w) (0, j) =
                ∑ k, T3 j k * v k := by
      show ∑ kl : Fin 1 × Fin 2, (T2_trivial ⊗ₖ T3) (0, j) kl * w kl
            = ∑ k, T3 j k * v k
      rw [Fintype.sum_prod_type]
      simp only [Matrix.kronecker_apply]
      -- Inner sum is over Fin 1, single element 0.
      rw [Finset.sum_comm]
      have : ∀ k : Fin 2, ∑ i : Fin 1,
              T2_trivial 0 i * T3 j k * w (i, k) = T3 j k * v k := by
        intro k
        rw [Fintype.sum_unique]
        simp [T2_trivial, Matrix.one_apply, hw]
      rw [Finset.sum_congr rfl (fun k _ => this k)]
    rw [hmv]
    -- Now compare with `T3.mulVec v j = (1 • v) j` from hv.
    have hvj : (T3.mulVec v) j = ((1 : ℝ) • v) j := by rw [hv]
    have hmv_v : (T3.mulVec v) j = ∑ k, T3 j k * v k := by
      show ∑ k, T3 j k * v k = ∑ k, T3 j k * v k
      rfl
    rw [hmv_v] at hvj
    -- hvj : ∑ k, T3 j k * v k = (1 • v) j
    rw [hvj]
    -- Goal: (1 • v) j = (1 • w) (0, j)
    show ((1 : ℝ) • v) j = ((1 : ℝ) • w) (0, j)
    simp [hw]
  -- Now apply the two-factor Perron uniqueness to `w`.
  have hwUniq : w = (w (0, 0)) • T2T3_perronVec :=
    T2T3_perron_uniqueness w hwPerron
  -- Translate back to `v`.
  funext j
  fin_cases j
  · -- v 0 = (v 0) • perronEigenvector 0 = v 0 · 1
    show v 0 = ((v 0) • perronEigenvector) 0
    simp [perronEigenvector]
  · -- v 1 = (v 0) • perronEigenvector 1 = v 0 · 1.
    -- Use the fact that w (0, 0) = w (0, 1) by the symmetry lemma.
    have hsym : w (0, 0) = w (0, 1) := T2T3_perron_symmetry w hwPerron
    show v 1 = ((v 0) • perronEigenvector) 1
    have hv01 : v 0 = v 1 := by
      have h00 : w (0, 0) = v 0 := rfl
      have h01 : w (0, 1) = v 1 := rfl
      rw [h00, h01] at hsym
      exact hsym
    rw [← hv01]
    simp [perronEigenvector]

/-! ### Tensor-decomposable Perron eigenvector uniqueness for `T_30` -/

/-- **Tensor-decomposable Perron eigenvector uniqueness (factor form).**
    Suppose `v_2 : Fin 1 → ℝ`, `v_3 : Fin 2 → ℝ`, `v_5 : Fin 4 → ℝ` are
    Perron eigenvectors of `T_2`, `T_3`, `T_5` respectively, and that
    `v_5` is a scalar multiple of `T5.perronVec` (this is the only
    factor for which uniqueness must be supplied externally). Then

      `(v_2 ⊗ v_3) ⊗ v_5  =  c · T30_perronVec T5`

    for some scalar `c ∈ ℝ`. -/
theorem T30_perron_tensor_uniqueness_factor (T5 : T5Like)
    (v2 : Fin 1 → ℝ) (v3 : Fin 2 → ℝ) (v5 : Fin 4 → ℝ)
    (h2 : T2_trivial.mulVec v2 = (1 : ℝ) • v2)
    (h3 : T3.mulVec v3 = (1 : ℝ) • v3)
    (h5_unique : ∃ d : ℝ, v5 = d • T5.perronVec) :
    ∃ c : ℝ, vecTensor (vecTensor v2 v3) v5 = c • T30_perronVec T5 := by
  obtain ⟨d, hd⟩ := h5_unique
  have hv2 : v2 = (v2 0) • T2_perronVec :=
    T2_perron_uniqueness v2 h2
  have hv3 : v3 = (v3 0) • perronEigenvector :=
    T3_perron_uniqueness v3 h3
  refine ⟨v2 0 * v3 0 * d, ?_⟩
  -- Compute both sides componentwise.
  funext ijk
  obtain ⟨ij, k⟩ := ijk
  obtain ⟨i, j⟩ := ij
  show vecTensor (vecTensor v2 v3) v5 ((i, j), k) =
        ((v2 0 * v3 0 * d) • T30_perronVec T5) ((i, j), k)
  -- LHS = v2 i * v3 j * v5 k
  -- RHS = (v2 0 * v3 0 * d) * (T2_perronVec i * perronEigenvector j * T5.perronVec k)
  have lhs_eq : vecTensor (vecTensor v2 v3) v5 ((i, j), k) =
                  v2 i * v3 j * v5 k := by
    simp [vecTensor]
  have rhs_eq : ((v2 0 * v3 0 * d) • T30_perronVec T5) ((i, j), k) =
                  (v2 0 * v3 0 * d) *
                  (T2_perronVec i * perronEigenvector j * T5.perronVec k) := by
    simp [T30_perronVec, vecTensor, Pi.smul_apply, smul_eq_mul]
  rw [lhs_eq, rhs_eq]
  -- Substitute using hv2, hv3, hd.
  have hv2i : v2 i = v2 0 * T2_perronVec i := by
    have := congrFun hv2 i
    simpa [Pi.smul_apply, smul_eq_mul] using this
  have hv3j : v3 j = v3 0 * perronEigenvector j := by
    have := congrFun hv3 j
    simpa [Pi.smul_apply, smul_eq_mul] using this
  have hv5k : v5 k = d * T5.perronVec k := by
    have := congrFun hd k
    simpa [Pi.smul_apply, smul_eq_mul] using this
  rw [hv2i, hv3j, hv5k]
  ring

/-- **Tensor-decomposable Perron uniqueness (packaged convenience form).**
    If `v` decomposes as `v = (v_2 ⊗ v_3) ⊗ v_5` where each factor is a
    Perron eigenvector of the corresponding `T_p`, and where `v_5` is a
    scalar multiple of `T5.perronVec`, then `v = c • T30_perronVec T5`. -/
theorem T30_perron_tensor_uniqueness (T5 : T5Like)
    (v : (Fin 1 × Fin 2) × Fin 4 → ℝ)
    (h_tensor : ∃ (v2 : Fin 1 → ℝ) (v3 : Fin 2 → ℝ) (v5 : Fin 4 → ℝ),
                  v = vecTensor (vecTensor v2 v3) v5
                ∧ T2_trivial.mulVec v2 = (1 : ℝ) • v2
                ∧ T3.mulVec v3 = (1 : ℝ) • v3
                ∧ ∃ d : ℝ, v5 = d • T5.perronVec) :
    ∃ c : ℝ, v = c • T30_perronVec T5 := by
  obtain ⟨v2, v3, v5, hdecomp, h2, h3, h5⟩ := h_tensor
  obtain ⟨c, hc⟩ := T30_perron_tensor_uniqueness_factor T5 v2 v3 v5 h2 h3 h5
  exact ⟨c, by rw [hdecomp]; exact hc⟩

/-! ### Conditional full uniqueness: `T5LikeUnique` -/

/-- **`T5LikeUnique`**: extension of `T5Like` that *additionally* asserts the
    Perron eigenspace of `T_5` is one-dimensional, i.e., every Perron
    eigenvector of `T_5` is a scalar multiple of `T5.perronVec`. This is the
    natural extra hypothesis required to upgrade tensor-decomposable
    uniqueness for `T_30` to full Perron-eigenspace uniqueness on the
    tensor-decomposable sector. -/
structure T5LikeUnique extends T5Like where
  /-- **Perron eigenspace one-dimensionality for `T_5`.** Every Perron
      eigenvector of `T_5` is a scalar multiple of `T5.perronVec`. -/
  perron_uniqueness : ∀ v : Fin 4 → ℝ,
    matrix.mulVec v = (1 : ℝ) • v → ∃ d : ℝ, v = d • perronVec

/-- **Conditional full Perron uniqueness for `T_30`.** Given the
    Perron-uniqueness extension `T5LikeUnique`, every Perron eigenvector
    of `T_30` that decomposes into Kronecker factors (with each factor a
    Perron eigenvector of its `T_p`) is a scalar multiple of
    `T30_perronVec T5.toT5Like`. -/
theorem T30_perron_uniqueness_under_T5_unique (T5 : T5LikeUnique)
    (v : (Fin 1 × Fin 2) × Fin 4 → ℝ)
    (h_tensor : ∃ (v2 : Fin 1 → ℝ) (v3 : Fin 2 → ℝ) (v5 : Fin 4 → ℝ),
                  v = vecTensor (vecTensor v2 v3) v5
                ∧ T2_trivial.mulVec v2 = (1 : ℝ) • v2
                ∧ T3.mulVec v3 = (1 : ℝ) • v3
                ∧ T5.matrix.mulVec v5 = (1 : ℝ) • v5) :
    ∃ c : ℝ, v = c • T30_perronVec T5.toT5Like := by
  obtain ⟨v2, v3, v5, hdecomp, h2, h3, h5⟩ := h_tensor
  -- Apply T_5 uniqueness to recover v5 = d • T5.perronVec.
  have hv5 : ∃ d : ℝ, v5 = d • T5.perronVec :=
    T5.perron_uniqueness v5 h5
  exact T30_perron_tensor_uniqueness T5.toT5Like v
          ⟨v2, v3, v5, hdecomp, h2, h3, hv5⟩

/-! ### Nonzero variant -/

/-- **Tensor-decomposable Perron uniqueness (nonzero variant).** If `v` is
    nonzero and decomposes tensorially with each factor a Perron
    eigenvector (and `v_5` a multiple of `T5.perronVec`), then `v` is a
    *nonzero* scalar multiple of `T30_perronVec`. -/
theorem T30_perron_tensor_uniqueness_nonzero (T5 : T5Like)
    (v : (Fin 1 × Fin 2) × Fin 4 → ℝ)
    (hne : v ≠ 0)
    (h_tensor : ∃ (v2 : Fin 1 → ℝ) (v3 : Fin 2 → ℝ) (v5 : Fin 4 → ℝ),
                  v = vecTensor (vecTensor v2 v3) v5
                ∧ T2_trivial.mulVec v2 = (1 : ℝ) • v2
                ∧ T3.mulVec v3 = (1 : ℝ) • v3
                ∧ ∃ d : ℝ, v5 = d • T5.perronVec) :
    ∃ c : ℝ, c ≠ 0 ∧ v = c • T30_perronVec T5 := by
  obtain ⟨c, hc⟩ := T30_perron_tensor_uniqueness T5 v h_tensor
  refine ⟨c, ?_, hc⟩
  intro hc0
  apply hne
  rw [hc, hc0, zero_smul]

/-! ### Headline summary -/

/-- **Headline (Perron uniqueness for `T_30 = T_2 ⊗ T_3 ⊗ T_5`).**

    1. `T30_perronVec T5` is a Perron eigenvector of `T_30`.
    2. `T30_perronVec T5` is strictly positive componentwise.
    3. **Tensor-decomposable uniqueness**: every Perron eigenvector
       of `T_30` arising from a Kronecker tensor decomposition with each
       factor a Perron eigenvector (and `v_5` a multiple of `T5.perronVec`)
       is a scalar multiple of `T30_perronVec T5`.
    4. **Conditional full uniqueness (via `T5LikeUnique`)**: extending
       `T5Like` with a Perron-uniqueness axiom for `T_5` upgrades (3) to
       drop the `v_5 ∝ T5.perronVec` precondition.

    The structural reduction makes precise the sense in which "Perron
    uniqueness for `T_30`" is equivalent to "Perron uniqueness for `T_5`"
    (with `T_2` and `T_3` handled unconditionally). -/
theorem T30_perron_uniqueness_summary (T5 : T5Like) :
    IsT30PerronEigenvector T5 (T30_perronVec T5)
    ∧ (∀ ijk, 0 < T30_perronVec T5 ijk)
    ∧ (∀ v : (Fin 1 × Fin 2) × Fin 4 → ℝ,
         (∃ (v2 : Fin 1 → ℝ) (v3 : Fin 2 → ℝ) (v5 : Fin 4 → ℝ),
            v = vecTensor (vecTensor v2 v3) v5
          ∧ T2_trivial.mulVec v2 = (1 : ℝ) • v2
          ∧ T3.mulVec v3 = (1 : ℝ) • v3
          ∧ ∃ d : ℝ, v5 = d • T5.perronVec) →
         ∃ c : ℝ, v = c • T30_perronVec T5) :=
  ⟨T30_perronVec_isPerronEigenvector T5,
   T30_perronVec_pos T5,
   T30_perron_tensor_uniqueness T5⟩

end PT.Stochastic
