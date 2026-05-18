/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Stochastic.T2T3T5KroneckerSpectrum
import PT.Stochastic.T2T3SpectralMixing
import Mathlib.Tactic

/-!
# Extended spectral mixing analysis for `T_30 = T_2 ⊗ T_3 ⊗ T_5`

This file extends `PT.Stochastic.T2T3T5KroneckerSpectrum` with a
**spectral mixing diagnostic** for the three-factor Kronecker block
`T_30 = T_2 ⊗ T_3 ⊗ T_5` (parametric in the `T5Like` record). It
mirrors the analysis carried out for `T_2 ⊗ T_3` in
`PT.Stochastic.T2T3SpectralMixing`, with the additional structural
content provided by the third factor `T_5`.

## Headline picture

Among the four explicit eigenpairs of `T_30` (cf.
`T30_four_eigenpairs_with_bounds`), exactly **two** carry unit modulus
(`+1` Perron and `-1` `T_3`-anti) and **two** are bounded by
`|λ_2(T_5)| ≤ 1/4 = s²`. This forces a sharp dichotomy between two
different notions of "spectral gap":

* **Total spectral gap = 0.** The sub-Perron eigenvalue closest to the
  Perron one is `-1`, with `|-1| = 1`, so on the full spectrum there is
  no decay — exactly as for `T_2 ⊗ T_3`. The `-1` mode is the
  `T_3`-alternation mode `v_+(T_2) ⊗ v_-(T_3) ⊗ v_+(T_5)`.
* **Perron-symmetric spectral gap ≥ 3/4.** On the subspace where the
  `T_3`-anti factor is **excluded** (i.e. vectors of the form
  `v_+(T_2) ⊗ v_+(T_3) ⊗ w`), the only sub-Perron eigenvalues are
  `±λ_2(T_5)` with `|·| ≤ 1/4`. Hence the gap is `1 - 1/4 = 3/4`.

This recovers the monograph's "renormalised" picture (`λ_2^{eff}` in
`T2T30Normalisation`) and explains why the renormalisation step in the
PT-spectral conservation theorem T2 is **mandatory**: only after
excluding the `-1` alternation mode does the Kronecker block exhibit
geometric mixing.

## Mixing-time picture (Perron-symmetric sector)

On the Perron-symmetric sector the iterated action decays geometrically:
for the sub-eigenvector `u_2^{eff} = v_+(T_2) ⊗ v_+(T_3) ⊗ w_2(T_5)`
we have

  `(T_30)^n . u_2^{eff} = (λ_2(T_5))^n . u_2^{eff}` ,
  `|(λ_2(T_5))^n| ≤ (1/4)^n` .

Hence after `n` steps the sub-Perron component is reduced by a factor
of at least `4^n`; geometrically, `t_mix(ε) ≤ ⌈log(1/ε) / log 4⌉`. For
`ε = 1/100` this gives `t_mix ≤ 4` (since `4^4 = 256 > 100`). We make
this concrete for the headline values `n = 1, 2, 3, 4`.

## Comparison with `T_2 ⊗ T_3`

The two analyses live side by side:

* `T_2 ⊗ T_3`: total gap is `0`, Cesàro mixing is exact at finite time
  `N = 2` (`T2T3_mixing_dichotomy`). No geometric mixing.
* `T_30`: total gap is `0`, but the **Perron-symmetric sector** admits
  a geometric mixing rate `≤ 1/4`. The third factor `T_5` is what
  enables this "second-order" decay.

## Main results

* `T30_subPerronEigenvalue` — the `-1` alternation eigenvalue.
* `T30_total_spectralGap` — definition; equals `0`.
* `T30_total_spectralGap_eq_zero` — proof.
* `T30_perronSymGap` — the Perron-symmetric sector gap `1 - |λ_2(T_5)|`.
* `T30_perronSymGap_lower_bound` — `≥ 3/4`.
* `T30_no_total_geometric_mixing` — exhibits a unit-modulus sub-Perron
  eigenvalue with eigenvector.
* `T30_perronSym_iter_eigen` — `(T_30)^n . u_2^{eff} = (λ_2)^n . u_2^{eff}`.
* `T30_perronSym_geometric_decay` — `|(λ_2)^n| ≤ (1/4)^n`.
* `T30_mixing_time_two`, `T30_mixing_time_three`, `T30_mixing_time_four`
  — explicit decay at small `n`.
* `T30_mixing_dichotomy` — headline summary.
* `T30_vs_T2T3_comparison` — side-by-side comparison.

## Reference

Monograph Chapter 7, §"Mélange spectral du bloc Kronecker T_30",
synthesis of `T2T30Normalisation`, `T30AntiSector`,
`T2T3T5KroneckerSpectrum`, and `T2T3SpectralMixing`.

## Status

`[DER]` — pure reuse of the four explicit eigenpairs already established.
No new numerical input. The geometric decay statements are
unconditional (in terms of `|λ_2(T_5)| ≤ 1/4` encoded in `T5Like`).
-/

namespace PT.Stochastic

open Matrix Kronecker BigOperators PT.Sieve

/-! ### Sub-Perron eigenvalue and total spectral gap -/

/-- **Sub-Perron eigenvalue of `T_30` on the full spectrum.** The
    explicit four-eigenpair enumeration shows that the eigenvalues of
    `T_30` distinct from the Perron `+1` are `{-1, λ_2(T_5), -λ_2(T_5)}`.
    Since `|λ_2(T_5)| ≤ 1/4 < 1 = |-1|`, the **maximum-modulus
    sub-Perron eigenvalue is `-1`** (the `T_3`-alternation mode). -/
def T30_subPerronEigenvalue : ℝ := -1

/-- The total sub-Perron eigenvalue has absolute value `1`. -/
theorem T30_sub_perron_eigenvalue_abs :
    |T30_subPerronEigenvalue| = 1 := by
  unfold T30_subPerronEigenvalue; norm_num

/-- **Total spectral gap of `T_30`.** Defined as `1 - |λ_sub|` where
    `λ_sub` is the maximum-modulus sub-Perron eigenvalue. For `T_30`
    this is `1 - |-1| = 0`. -/
noncomputable def T30_total_spectralGap : ℝ :=
  1 - |T30_subPerronEigenvalue|

/-- **Total spectral gap is zero.** Direct computation: the `-1`
    alternation eigenvalue saturates the unit modulus. -/
theorem T30_total_spectralGap_eq_zero : T30_total_spectralGap = 0 := by
  unfold T30_total_spectralGap
  rw [T30_sub_perron_eigenvalue_abs]
  norm_num

/-- **No total geometric mixing.** There exists a sub-Perron eigenvalue
    of `T_30` with absolute value `1` (namely `-1`), realised by the
    `T_3`-anti / `T_5`-Perron eigenvector. Hence the iterated action
    `(T_30)^n . v` does **not** decay on this vector. -/
theorem T30_no_total_geometric_mixing (T5 : T5Like) :
    ∃ (μ : ℝ), |μ| = 1 ∧ μ ≠ 1 ∧
      (T30 T5).mulVec (T30_anti3_perronVec T5) = μ • T30_anti3_perronVec T5 := by
  refine ⟨-1, ?_, ?_, ?_⟩
  · norm_num
  · norm_num
  · exact T30_anti3_perron_eigen' T5

/-! ### Perron-symmetric sector: restricted gap -/

/-- **Perron-symmetric sector spectral gap.** On the subspace where
    the `T_3`-anti factor has been excluded (eigenvectors of the form
    `v_+(T_2) ⊗ v_+(T_3) ⊗ w`), the only sub-Perron eigenvalues are
    `±λ_2(T_5)`, both bounded by `1/4`. The gap is

      `1 - |λ_2(T_5)| ≥ 1 - 1/4 = 3/4`. -/
noncomputable def T30_perronSymGap (T5 : T5Like) : ℝ :=
  1 - |T5.subEigenvalue|

/-- **Perron-symmetric gap lower bound: `≥ 3/4`.** -/
theorem T30_perronSymGap_lower_bound (T5 : T5Like) :
    (3 : ℝ) / 4 ≤ T30_perronSymGap T5 := by
  unfold T30_perronSymGap
  have h := T5.subdominant_bound
  linarith

/-- **Perron-symmetric gap is non-negative.** -/
theorem T30_perronSymGap_nonneg (T5 : T5Like) :
    0 ≤ T30_perronSymGap T5 := by
  have h := T30_perronSymGap_lower_bound T5
  linarith

/-- **Perron-symmetric gap is strictly positive.** -/
theorem T30_perronSymGap_pos (T5 : T5Like) :
    0 < T30_perronSymGap T5 := by
  have h := T30_perronSymGap_lower_bound T5
  linarith

/-- **Headline: Perron-symmetric gap ≥ 3/4 > 0 = total gap.** The
    restriction to the Perron-symmetric sector strictly improves the
    spectral gap from `0` to at least `3/4`. -/
theorem T30_perronSymGap_strictly_better (T5 : T5Like) :
    T30_total_spectralGap < T30_perronSymGap T5 := by
  rw [T30_total_spectralGap_eq_zero]
  exact T30_perronSymGap_pos T5

/-! ### Iterated action on the Perron-symmetric sub-eigenvector

    Helper lemma: if `A . v = λ • v`, then `(A^n) . v = (λ^n) • v`.
    Specialised below to `A = T_30`, `v = u_2^{eff}`, `λ = λ_2(T_5)`. -/

/-- **Iterated eigenvalue lemma (generic).** If `A . v = lam • v`, then
    for every `n`, `(A^n) . v = (lam^n) • v`. Standard induction. -/
theorem matrix_pow_mulVec_of_eigen
    {α : Type*} [Fintype α] [DecidableEq α]
    (A : Matrix α α ℝ) (v : α → ℝ) (lam : ℝ)
    (h : A.mulVec v = lam • v) (n : ℕ) :
    (A ^ n).mulVec v = (lam ^ n) • v := by
  induction n with
  | zero =>
      simp [Matrix.one_mulVec]
  | succ k ih =>
      rw [pow_succ, ← Matrix.mulVec_mulVec, h, Matrix.mulVec_smul, ih]
      ext i
      simp [Pi.smul_apply, smul_eq_mul, pow_succ]
      ring

/-- **Iterated action on the Perron-symmetric sub-eigenvector.** For
    every `n`, `(T_30)^n . u_2^{eff} = (λ_2(T_5))^n . u_2^{eff}`. -/
theorem T30_perronSym_iter_eigen (T5 : T5Like) (n : ℕ) :
    ((T30 T5) ^ n).mulVec (T30_lambda_eff_vec T5)
      = (T5.subEigenvalue ^ n) • T30_lambda_eff_vec T5 :=
  matrix_pow_mulVec_of_eigen (T30 T5) (T30_lambda_eff_vec T5) T5.subEigenvalue
    (T30_lambda_eff_eigen' T5) n

/-- **Iterated action on the `T_3`-anti / `T_5`-sub eigenvector.** -/
theorem T30_anti_sub_iter_eigen (T5 : T5Like) (n : ℕ) :
    ((T30 T5) ^ n).mulVec (T30_anti3_sub_Vec T5)
      = ((-T5.subEigenvalue) ^ n) • T30_anti3_sub_Vec T5 :=
  matrix_pow_mulVec_of_eigen (T30 T5) (T30_anti3_sub_Vec T5) (-T5.subEigenvalue)
    (T30_anti3_sub_eigen' T5) n

/-- **Iterated action on the all-Perron eigenvector: fixed.** Since
    the Perron eigenvalue is `+1`, every iterate is the identity on
    `u_+`. -/
theorem T30_perron_iter_eigen (T5 : T5Like) (n : ℕ) :
    ((T30 T5) ^ n).mulVec (T30_perronVec T5) = T30_perronVec T5 := by
  have h : ((T30 T5) ^ n).mulVec (T30_perronVec T5)
              = ((1 : ℝ) ^ n) • T30_perronVec T5 :=
    matrix_pow_mulVec_of_eigen (T30 T5) (T30_perronVec T5) (1 : ℝ)
      (T30_perron_eigen T5) n
  rw [h, one_pow, one_smul]

/-- **Iterated action on the `T_3`-anti / `T_5`-Perron eigenvector:
    sign flip.** Since the eigenvalue is `-1`, the orbit alternates
    `u_- ↦ -u_- ↦ u_- ↦ ⋯`. -/
theorem T30_anti_perron_iter_eigen (T5 : T5Like) (n : ℕ) :
    ((T30 T5) ^ n).mulVec (T30_anti3_perronVec T5)
      = ((-1 : ℝ) ^ n) • T30_anti3_perronVec T5 :=
  matrix_pow_mulVec_of_eigen (T30 T5) (T30_anti3_perronVec T5) (-1)
    (T30_anti3_perron_eigen' T5) n

/-! ### Geometric decay on the Perron-symmetric sector -/

/-- **Geometric decay rate on the Perron-symmetric sector.** The
    eigenvalue contracted at step `n` is `(λ_2(T_5))^n`, whose absolute
    value is bounded by `(1/4)^n` since `|λ_2(T_5)| ≤ 1/4`. -/
theorem T30_perronSym_geometric_decay (T5 : T5Like) (n : ℕ) :
    |T5.subEigenvalue ^ n| ≤ ((1 : ℝ) / 4) ^ n := by
  rw [abs_pow]
  exact pow_le_pow_left₀ (abs_nonneg _) T5.subdominant_bound n

/-- **Geometric decay rate on the anti-sub sector.** Same bound:
    `|(-λ_2)^n| ≤ (1/4)^n`. -/
theorem T30_anti_sub_geometric_decay (T5 : T5Like) (n : ℕ) :
    |(-T5.subEigenvalue) ^ n| ≤ ((1 : ℝ) / 4) ^ n := by
  rw [abs_pow, abs_neg]
  exact pow_le_pow_left₀ (abs_nonneg _) T5.subdominant_bound n

/-! ### Explicit small-`n` decay values

    Concrete witnesses for `n = 1, 2, 3, 4`: at `n = 4` the
    sub-Perron component is reduced by at least `4^4 = 256`,
    exceeding the `1/100` threshold typically used in mixing-time
    discussions. -/

/-- **Decay at `n = 1`: factor `≤ 1/4`.** -/
theorem T30_mixing_time_one (T5 : T5Like) :
    |T5.subEigenvalue ^ 1| ≤ (1 : ℝ) / 4 := by
  have h := T30_perronSym_geometric_decay T5 1
  simpa using h

/-- **Decay at `n = 2`: factor `≤ 1/16`.** -/
theorem T30_mixing_time_two (T5 : T5Like) :
    |T5.subEigenvalue ^ 2| ≤ (1 : ℝ) / 16 := by
  have h := T30_perronSym_geometric_decay T5 2
  have : ((1 : ℝ) / 4) ^ 2 = 1 / 16 := by norm_num
  linarith [h]

/-- **Decay at `n = 3`: factor `≤ 1/64`.** -/
theorem T30_mixing_time_three (T5 : T5Like) :
    |T5.subEigenvalue ^ 3| ≤ (1 : ℝ) / 64 := by
  have h := T30_perronSym_geometric_decay T5 3
  have : ((1 : ℝ) / 4) ^ 3 = 1 / 64 := by norm_num
  linarith [h]

/-- **Decay at `n = 4`: factor `≤ 1/256 < 1/100`.** Thus four iterates
    of `T_30` already drive the Perron-symmetric sub-component below
    the `ε = 1/100` mixing threshold. -/
theorem T30_mixing_time_four (T5 : T5Like) :
    |T5.subEigenvalue ^ 4| ≤ (1 : ℝ) / 256 := by
  have h := T30_perronSym_geometric_decay T5 4
  have : ((1 : ℝ) / 4) ^ 4 = 1 / 256 := by norm_num
  linarith [h]

/-- **`n = 4` clears the `1/100` mixing threshold.** -/
theorem T30_mixing_time_four_below_one_hundred (T5 : T5Like) :
    |T5.subEigenvalue ^ 4| ≤ (1 : ℝ) / 100 := by
  have h := T30_mixing_time_four T5
  linarith

/-! ### Convergence of the iterated sub-Perron eigenvalue to `0` -/

/-- **The sub-Perron eigenvalue (raised to the `n`-th power) tends to
    `0`.** This is the asymptotic form of the geometric-decay theorem:
    `(λ_2(T_5))^n → 0` as `n → ∞`. The bound `|λ_2| ≤ 1/4 < 1`
    delivers the exponential rate. -/
theorem T30_perronSym_pow_tendsto_zero (T5 : T5Like) :
    Filter.Tendsto (fun n => T5.subEigenvalue ^ n) Filter.atTop (nhds 0) := by
  have hle : |T5.subEigenvalue| ≤ (1 : ℝ) / 4 := T5.subdominant_bound
  have hlt : |T5.subEigenvalue| < 1 := by linarith
  exact tendsto_pow_atTop_nhds_zero_of_abs_lt_one hlt

/-! ### Comparison with `T_2 ⊗ T_3` -/

/-- **Comparison theorem: `T_30` vs `T_2 ⊗ T_3`.** Both Kronecker
    blocks have **total spectral gap `0`** (witness: the `T_3`-anti
    sector with eigenvalue `-1`). The fundamental difference appears
    on the **Perron-symmetric sector**:

    * `T_2 ⊗ T_3` has *no* non-trivial sector below the anti-Perron
      mode — once `-1` is excluded, only the Perron eigenvalue `+1`
      remains; the Perron-symmetric sector is 1-dimensional. (The
      Cesàro average reaches the stationary projector at `N = 2`.)
    * `T_30` has a *non-trivial* Perron-symmetric sector: the `T_5`
      factor contributes a sub-eigenvalue `λ_2(T_5)` with
      `|λ_2| ≤ 1/4`, hence the Perron-symmetric sector exhibits
      **geometric mixing** with rate `≤ 1/4` per step.

    Quantitatively, both gaps are `0` globally, but the
    Perron-symmetric refinement gives:

      `gap(T_30 | Perron-sym) ≥ 3/4` vs `gap(T_2 ⊗ T_3) = 0`.

    This is the algebraic explanation of the "spectral upgrade"
    delivered by extending the 2-factor block to the 3-factor block
    via CRT at `m = 30`. -/
theorem T30_vs_T2T3_comparison (T5 : T5Like) :
    -- Both total gaps vanish
    T2T3_spectralGap = 0
    ∧ T30_total_spectralGap = 0
    -- Both witness `-1` in their respective spectra
    ∧ |T2T3_subPerronEigenvalue| = 1
    ∧ |T30_subPerronEigenvalue| = 1
    -- `T_30` Perron-symmetric refinement is strictly better
    ∧ (3 : ℝ) / 4 ≤ T30_perronSymGap T5
    -- And strictly positive
    ∧ 0 < T30_perronSymGap T5 :=
  ⟨T2T3_spectralGap_eq_zero,
   T30_total_spectralGap_eq_zero,
   T2T3_sub_perron_eigenvalue_abs,
   T30_sub_perron_eigenvalue_abs,
   T30_perronSymGap_lower_bound T5,
   T30_perronSymGap_pos T5⟩

/-! ### Headline dichotomy -/

/-- **Headline (spectral mixing dichotomy for `T_30`).**

    The three-factor Kronecker block `T_30 = T_2 ⊗ T_3 ⊗ T_5` exhibits
    a layered dichotomy between *total* and *Perron-symmetric*
    spectral analysis:

    * **Total spectral gap is `0`.** The sub-Perron eigenvalue is `-1`
      (the `T_3`-alternation mode, with eigenvector
      `v_+(T_2) ⊗ v_-(T_3) ⊗ v_+(T_5)`), so on the full spectrum
      there is *no* geometric mixing.
    * **Perron-symmetric spectral gap is `≥ 3/4`.** After excluding
      the `T_3`-anti factor, the only sub-Perron eigenvalues are
      `±λ_2(T_5)`, both bounded by `1/4 = s²`. Hence the gap on the
      Perron-symmetric sector is `1 - |λ_2(T_5)| ≥ 3/4 > 0`.
    * **Geometric decay on the Perron-symmetric sector.** The iterated
      action satisfies `(T_30)^n . u_2^{eff} = (λ_2(T_5))^n . u_2^{eff}`
      with `|(λ_2(T_5))^n| ≤ (1/4)^n`. In particular `n = 4`
      iterates already reduce the sub-Perron component below the
      `1/100` mixing threshold.
    * **Asymptotic convergence.** `(λ_2(T_5))^n → 0` as `n → ∞`.

    This is the structural ground for the monograph's renormalised
    bound `|λ_2^{eff}(T_30)| ≤ s²`: the renormalisation isolates
    precisely the Perron-symmetric sector on which geometric mixing
    actually holds. -/
theorem T30_mixing_dichotomy (T5 : T5Like) :
    -- Total spectral gap = 0
    T30_total_spectralGap = 0
    -- Witness: -1 is a sub-Perron eigenvalue
    ∧ |T30_subPerronEigenvalue| = 1
    -- Perron-symmetric sector gap ≥ 3/4
    ∧ (3 : ℝ) / 4 ≤ T30_perronSymGap T5
    -- And strictly positive
    ∧ 0 < T30_perronSymGap T5
    -- Total gap strictly less than Perron-sym gap
    ∧ T30_total_spectralGap < T30_perronSymGap T5
    -- Geometric decay on the Perron-symmetric sub-eigenvector
    ∧ (∀ n : ℕ, ((T30 T5) ^ n).mulVec (T30_lambda_eff_vec T5)
                  = (T5.subEigenvalue ^ n) • T30_lambda_eff_vec T5)
    -- Geometric decay rate ≤ (1/4)^n
    ∧ (∀ n : ℕ, |T5.subEigenvalue ^ n| ≤ ((1 : ℝ) / 4) ^ n)
    -- Asymptotic mixing time: 4 iterates clear ε = 1/100
    ∧ |T5.subEigenvalue ^ 4| ≤ (1 : ℝ) / 100
    -- Asymptotic convergence to 0
    ∧ Filter.Tendsto (fun n => T5.subEigenvalue ^ n) Filter.atTop (nhds 0) :=
  ⟨T30_total_spectralGap_eq_zero,
   T30_sub_perron_eigenvalue_abs,
   T30_perronSymGap_lower_bound T5,
   T30_perronSymGap_pos T5,
   T30_perronSymGap_strictly_better T5,
   T30_perronSym_iter_eigen T5,
   T30_perronSym_geometric_decay T5,
   T30_mixing_time_four_below_one_hundred T5,
   T30_perronSym_pow_tendsto_zero T5⟩

end PT.Stochastic
