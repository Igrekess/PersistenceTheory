/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Stochastic.T2T30Normalisation
import Mathlib.Tactic

/-!
# `T_30` anti-sector eigenvalues (Ch07 extension)

This file extends `PT.Stochastic.T2T30Normalisation` by enumerating the
**anti-symmetric (T_3-anti) sector eigenvalues** of `T_30 = T_2 ⊗ T_3 ⊗ T_5`:

* The **all-Perron** eigenpair gives `λ = 1` (covered by `T30_perron_eigen`).
* The **T_3-anti / T_5-Perron** eigenpair: `(v_+(T_2) ⊗ v_-(T_3) ⊗ v_+(T_5))`
  has eigenvalue `1 · (-1) · 1 = -1`.
* The **T_3-anti / T_5-sub** eigenpair: `(v_+(T_2) ⊗ v_-(T_3) ⊗ w_2(T_5))`
  has eigenvalue `1 · (-1) · λ_2(T_5) = -λ_2(T_5)` with `|·| ≤ 1/4`.
* The **T_5-sub only** sector is `T30_lambda_eff_eigen` in the base file.

The **absolute-value bound** on the spectral radius of `T_30` (restricted to
the full Perron and anti sectors covered here) is `max(1, |λ_2(T_5)|) = 1`
since `|λ_2(T_5)| ≤ 1/4 < 1`.

## Reference

Monograph Chapter 7, §"Spectre de `T_{30}`", follow-up to
`T2T30Normalisation`.
-/

namespace PT.Stochastic

open Matrix PT.Sieve

/-! ### Anti-symmetric tensor on the `T_3` factor -/

/-- The anti-symmetric tensor on the `T_3` factor:
    `v_+(T_2) ⊗ v_-(T_3) ⊗ v_+(T_5)`. -/
def T30_anti3_perronVec (T5 : T5Like) :
    (Fin 1 × Fin 2) × Fin 4 → ℝ :=
  vecTensor (vecTensor T2_perronVec antiEigenvector) T5.perronVec

/-- **`T_3`-anti / `T_5`-Perron eigenpair.**
    `T_30 . (v_+(T_2) ⊗ v_-(T_3) ⊗ v_+(T_5)) = (1 · (-1) · 1) · (...)`. -/
theorem T30_anti3_perron_eigen (T5 : T5Like) :
    (T30 T5).mulVec (T30_anti3_perronVec T5)
      = ((1 : ℝ) * (-1) * 1) • T30_anti3_perronVec T5 := by
  unfold T30 T30_anti3_perronVec
  exact kron3_eigenvector_left T2_trivial T3 T5.matrix
          T2_perronVec antiEigenvector T5.perronVec
          1 (-1) 1
          T2_perron_eigen T3_sub_eigen T5.perron_eigen

/-- Cleaner restatement: the eigenvalue is `-1`. -/
theorem T30_anti3_perron_eigen' (T5 : T5Like) :
    (T30 T5).mulVec (T30_anti3_perronVec T5)
      = (-1 : ℝ) • T30_anti3_perronVec T5 := by
  have h := T30_anti3_perron_eigen T5
  simpa [neg_mul, one_mul] using h

/-! ### Anti-symmetric on `T_3` × sub-eigenvector on `T_5` -/

/-- The "fully sub-dominant" tensor: `v_+(T_2) ⊗ v_-(T_3) ⊗ w_2(T_5)`. -/
def T30_anti3_sub_Vec (T5 : T5Like) :
    (Fin 1 × Fin 2) × Fin 4 → ℝ :=
  vecTensor (vecTensor T2_perronVec antiEigenvector) T5.subVec

/-- **`T_3`-anti / `T_5`-sub eigenpair.**
    Eigenvalue `1 · (-1) · λ_2(T_5) = -λ_2(T_5)`. -/
theorem T30_anti3_sub_eigen (T5 : T5Like) :
    (T30 T5).mulVec (T30_anti3_sub_Vec T5)
      = ((1 : ℝ) * (-1) * T5.subEigenvalue) • T30_anti3_sub_Vec T5 := by
  unfold T30 T30_anti3_sub_Vec
  exact kron3_eigenvector_left T2_trivial T3 T5.matrix
          T2_perronVec antiEigenvector T5.subVec
          1 (-1) T5.subEigenvalue
          T2_perron_eigen T3_sub_eigen T5.sub_eigen

/-- Cleaner restatement: the eigenvalue is `-λ_2(T_5)`. -/
theorem T30_anti3_sub_eigen' (T5 : T5Like) :
    (T30 T5).mulVec (T30_anti3_sub_Vec T5)
      = (-T5.subEigenvalue) • T30_anti3_sub_Vec T5 := by
  have h := T30_anti3_sub_eigen T5
  simpa [neg_mul, one_mul] using h

/-! ### Absolute-value bound on the anti-sector -/

/-- The anti-sector sub-eigenvalue is bounded by `1/4` in absolute value
    (inherits the `T_5` bound). -/
theorem T30_anti3_sub_abs_bound (T5 : T5Like) :
    |(-T5.subEigenvalue)| ≤ (1 : ℝ) / 4 := by
  rw [abs_neg]
  exact T5.subdominant_bound

/-! ### Headline (anti-sector summary) -/

/-- **Headline (anti-sector summary).** Four eigenpairs of `T_30` on the
    Perron + anti-sector are now made explicit:

    * `(1, v_+(T_2) ⊗ v_+(T_3) ⊗ v_+(T_5))` — all-Perron (see base file)
    * `(-1, v_+(T_2) ⊗ v_-(T_3) ⊗ v_+(T_5))` — T_3-anti, T_5-Perron
    * `(λ_2(T_5), v_+(T_2) ⊗ v_+(T_3) ⊗ w_2(T_5))` — T_5-sub (see base file)
    * `(-λ_2(T_5), v_+(T_2) ⊗ v_-(T_3) ⊗ w_2(T_5))` — both-sub

    All sub-dominant eigenvalues `(±λ_2(T_5))` satisfy `|·| ≤ 1/4 = s²`. -/
theorem T30_perron_anti_sector_summary (T5 : T5Like) :
    -- Perron eigenvalue 1
    (T30 T5).mulVec (T30_perronVec T5) = (1 : ℝ) • T30_perronVec T5
    -- T_3-anti eigenvalue -1
    ∧ (T30 T5).mulVec (T30_anti3_perronVec T5)
        = (-1 : ℝ) • T30_anti3_perronVec T5
    -- T_5-sub eigenvalue λ_2(T_5)
    ∧ (T30 T5).mulVec (T30_lambda_eff_vec T5)
        = T5.subEigenvalue • T30_lambda_eff_vec T5
    -- T_3-anti / T_5-sub eigenvalue -λ_2(T_5)
    ∧ (T30 T5).mulVec (T30_anti3_sub_Vec T5)
        = (-T5.subEigenvalue) • T30_anti3_sub_Vec T5
    -- Both sub-eigenvalues bounded by 1/4
    ∧ |T5.subEigenvalue| ≤ (1 : ℝ) / 4
    ∧ |(-T5.subEigenvalue)| ≤ (1 : ℝ) / 4 :=
  ⟨T30_perron_eigen T5,
   T30_anti3_perron_eigen' T5,
   T30_lambda_eff_eigen' T5,
   T30_anti3_sub_eigen' T5,
   T5.subdominant_bound,
   T30_anti3_sub_abs_bound T5⟩

end PT.Stochastic
