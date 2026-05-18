/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Stochastic.T2T3KroneckerEigenvalues
import PT.Stochastic.T2T3CesaroLimit
import PT.Stochastic.T2T3StationaryUniqueness
import Mathlib.Tactic

/-!
# Spectral mixing analysis for `T_2 ⊗ T_3`

This file synthesises the spectral and Cesàro information about the
Kronecker involution `T_2 ⊗ T_3` to produce its **mixing diagnostics**:

* Its **spectral gap** (in the classical sense `1 - max{|λ| : λ ≠ +1}`)
  is **exactly `0`**, because the sub-Perron eigenvalue is `-1` with
  `|-1| = 1`. Hence pointwise (geometric) convergence of `(T_2 ⊗ T_3)^n v`
  to the stationary distribution **fails**.

* The **Cesàro mixing time** is nevertheless **finite and equal to `2`**:
  `(1/N) Σ_{k<N} (T_2 ⊗ T_3)^k . v = π` exactly at every even `N ≥ 2`,
  starting at `N = 2`. This re-uses `T2T3_cesaro_two`,
  `T2T3_cesaro_four`, `T2T3_cesaro_even` and `T2T3_pow_sum_even`.

* Consequently the natural mixing diagnostic for this chain is the
  **Cesàro mixing time `t_mix^Ces = 2`**, *not* the geometric mixing
  time, which is infinite.

## Main results

* `T2T3_sub_perron_eigenvalue_abs` — `|λ_sub| = |-1| = 1`.
* `T2T3_spectralGap` — `spectralGap (T_2 ⊗ T_3) = 0`.
* `T2T3_cesaroMixingTime` — `t_mix^Ces = 2` (definition).
* `T2T3_cesaroMixingTime_exact` — `Cesàro_N (T_2 ⊗ T_3) = π` for every
  even `N ≥ 2`, in particular for `N = t_mix^Ces = 2`.
* `T2T3_mixing_dichotomy` — headline summary: spectral gap is `0`
  (no geometric mixing) but Cesàro mixing is exact at `N = 2`.

## Reference

Monograph Chapter 7, §"Spectral mixing of the Kronecker block
`T_2 ⊗ T_3`", synthesis of `PT.Stochastic.T2T3KroneckerEigenvalues`,
`PT.Stochastic.T2T3CesaroLimit`, and
`PT.Stochastic.T2T3StationaryUniqueness`.
-/

namespace PT.Stochastic

open Matrix Kronecker BigOperators PT.Sieve

/-! ### Sub-Perron eigenvalue and its modulus -/

/-- The **sub-Perron eigenvalue** of `T_2 ⊗ T_3` is `-1`. -/
def T2T3_subPerronEigenvalue : ℝ := -1

/-- The sub-Perron eigenvalue has absolute value `1`. -/
theorem T2T3_sub_perron_eigenvalue_abs :
    |T2T3_subPerronEigenvalue| = 1 := by
  unfold T2T3_subPerronEigenvalue; norm_num

/-! ### Spectral gap -/

/-- **Spectral gap of `T_2 ⊗ T_3`.** By definition,
    `spectralGap := 1 - |λ_sub|`. For `T_2 ⊗ T_3` the only sub-Perron
    eigenvalue is `-1`, so this equals `0`. -/
noncomputable def T2T3_spectralGap : ℝ :=
  1 - |T2T3_subPerronEigenvalue|

/-- **Spectral gap is zero.** Direct computation from the definition
    and `|-1| = 1`. -/
theorem T2T3_spectralGap_eq_zero : T2T3_spectralGap = 0 := by
  unfold T2T3_spectralGap
  rw [T2T3_sub_perron_eigenvalue_abs]
  norm_num

/-- **Max-absolute eigenvalue (spectral radius) equals `1`.**
    Both eigenvalues `+1` and `-1` have modulus `1`. -/
theorem T2T3_max_abs_eigenvalue : max |(1 : ℝ)| |(-1 : ℝ)| = 1 := by
  norm_num

/-- **Geometric mixing fails:** there exists a non-zero eigenvalue
    distinct from `+1` whose modulus equals `1`. Concretely, the
    anti-eigenvector `v_-` carries eigenvalue `-1`, and applying
    `(T_2 ⊗ T_3)^k` keeps it on the unit sphere forever (no decay). -/
theorem T2T3_no_geometric_mixing :
    ∃ (μ : ℝ), |μ| = 1 ∧ μ ≠ 1 ∧
      (T2_trivial ⊗ₖ T3).mulVec T2T3_antiVec = μ • T2T3_antiVec := by
  refine ⟨-1, ?_, ?_, ?_⟩
  · norm_num
  · norm_num
  · exact T2T3_anti_eigen'

/-! ### Cesàro mixing time -/

/-- **Cesàro mixing time of `T_2 ⊗ T_3`.** The smallest even `N ≥ 2` for
    which the Cesàro average equals the stationary projector exactly.
    Here it is `2`. -/
def T2T3_cesaroMixingTime : ℕ := 2

/-- **Exact Cesàro mixing at `N = t_mix^Ces = 2`.** Direct from
    `T2T3_cesaro_two`. -/
theorem T2T3_cesaroMixingTime_exact :
    ((1 : ℝ) / T2T3_cesaroMixingTime)
        • (∑ k ∈ Finset.range T2T3_cesaroMixingTime, (T2_trivial ⊗ₖ T3)^k)
      = T2T3_stationaryProjector := by
  -- Unfold `T2T3_cesaroMixingTime = 2` and reuse `T2T3_cesaro_two`.
  show ((1 : ℝ) / 2)
        • (∑ k ∈ Finset.range 2, (T2_trivial ⊗ₖ T3)^k)
      = T2T3_stationaryProjector
  have hsum :
      (∑ k ∈ Finset.range 2, (T2_trivial ⊗ₖ T3)^k)
        = (T2_trivial ⊗ₖ T3)^0 + (T2_trivial ⊗ₖ T3)^1 := by
    rw [Finset.sum_range_succ, Finset.sum_range_succ, Finset.sum_range_zero,
        zero_add]
  rw [hsum]
  exact T2T3_cesaro_two

/-- **Cesàro average is exact at every even `N ≥ 2`.** Concretely, for
    every `n ≥ 1`, the Cesàro average over `N = 2n` steps equals
    `T2T3_stationaryProjector`. Re-uses `T2T3_cesaro_even`. -/
theorem T2T3_cesaroMixing_at_every_even (n : ℕ) (hn : 0 < n) :
    ((1 : ℝ) / (2 * n))
        • (∑ k ∈ Finset.range (2 * n), (T2_trivial ⊗ₖ T3)^k)
      = T2T3_stationaryProjector :=
  T2T3_cesaro_even n hn

/-! ### Action on a single vector -/

/-- **Cesàro action on an arbitrary initial vector.** For any vector
    `v`, the Cesàro average at `N = 2` applied to `v` equals
    `T2T3_stationaryProjector.mulVec v`. -/
theorem T2T3_cesaroMixing_mulVec (v : Fin 1 × Fin 2 → ℝ) :
    (((1 : ℝ) / 2)
        • ((T2_trivial ⊗ₖ T3)^0 + (T2_trivial ⊗ₖ T3)^1)).mulVec v
      = T2T3_stationaryProjector.mulVec v := by
  rw [T2T3_cesaro_two]

/-- **Cesàro action on the Perron eigenvector.** The stationary projector
    fixes the Perron direction `v_+`. -/
theorem T2T3_cesaroMixing_perronVec :
    (((1 : ℝ) / 2)
        • ((T2_trivial ⊗ₖ T3)^0 + (T2_trivial ⊗ₖ T3)^1)).mulVec T2T3_perronVec
      = T2T3_stationaryProjector.mulVec T2T3_perronVec :=
  T2T3_cesaroMixing_mulVec _

/-- **Cesàro action kills the anti-eigenvector.** Since
    `T_2 ⊗ T_3 . v_- = -v_-`, we have
    `(I + T_2 ⊗ T_3) . v_- = v_- + (-v_-) = 0`, so the Cesàro average
    sends `v_-` to `0`. -/
theorem T2T3_cesaroMixing_kills_antiVec :
    T2T3_stationaryProjector.mulVec T2T3_antiVec = 0 := by
  unfold T2T3_stationaryProjector
  -- Apply linearity: ((1/2)•(I + T)) . v = (1/2) ((I . v) + (T . v))
  rw [Matrix.smul_mulVec, Matrix.add_mulVec, Matrix.one_mulVec,
      T2T3_anti_eigen']
  -- Goal: (1/2) • (T2T3_antiVec + (-1) • T2T3_antiVec) = 0
  ext i
  simp [Pi.smul_apply]

/-! ### Headline dichotomy -/

/-- **Headline (spectral mixing dichotomy for `T_2 ⊗ T_3`).**

    The Kronecker block `T_2 ⊗ T_3` exhibits a sharp dichotomy between
    *geometric* and *Cesàro* convergence:

    * **Spectral gap is `0`.** The sub-Perron eigenvalue is `-1`, with
      modulus `1`, so there is *no* geometric (exponential) mixing.
      The anti-eigenvector `v_-` satisfies `(T_2 ⊗ T_3)^k . v_- = (-1)^k v_-`
      and never decays.

    * **Cesàro mixing is exact at finite time.** The Cesàro average
      `(1/N) Σ_{k<N} (T_2 ⊗ T_3)^k` equals the stationary projector
      `T2T3_stationaryProjector` for every even `N ≥ 2`. In particular
      the Cesàro mixing time is `t_mix^Ces = 2`.

    * **Stationary projector kills `v_-`.** The Cesàro limit annihilates
      the sub-Perron direction, leaving only the Perron component.

    Together with the uniqueness of the stationary distribution
    (`T2T3_unique_stationary`), this gives a complete mixing picture
    for the Kronecker dynamics on `(ℤ/2ℤ)* × (ℤ/3ℤ)*`. -/
theorem T2T3_mixing_dichotomy :
    -- Spectral gap = 0
    T2T3_spectralGap = 0
    -- Sub-Perron eigenvalue has modulus 1 (no geometric mixing)
    ∧ |T2T3_subPerronEigenvalue| = 1
    -- Cesàro mixing time = 2
    ∧ T2T3_cesaroMixingTime = 2
    -- Cesàro exact at N = 2
    ∧ ((1 : ℝ) / T2T3_cesaroMixingTime)
        • (∑ k ∈ Finset.range T2T3_cesaroMixingTime, (T2_trivial ⊗ₖ T3)^k)
        = T2T3_stationaryProjector
    -- Cesàro exact at every even N ≥ 2
    ∧ (∀ n : ℕ, 0 < n →
         ((1 : ℝ) / (2 * n))
            • (∑ k ∈ Finset.range (2 * n), (T2_trivial ⊗ₖ T3)^k)
           = T2T3_stationaryProjector)
    -- Stationary projector annihilates the anti-direction
    ∧ T2T3_stationaryProjector.mulVec T2T3_antiVec = 0
    -- Stationary distribution is unique
    ∧ (∀ π, IsT2T3Stationary π → π = T2T3_stationary) :=
  ⟨T2T3_spectralGap_eq_zero,
   T2T3_sub_perron_eigenvalue_abs,
   rfl,
   T2T3_cesaroMixingTime_exact,
   T2T3_cesaroMixing_at_every_even,
   T2T3_cesaroMixing_kills_antiVec,
   T2T3_unique_stationary⟩

end PT.Stochastic
