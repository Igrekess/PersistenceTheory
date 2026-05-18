/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Stochastic.SpectralDominance
import PT.Stochastic.T3SpectralDecomposition
import Mathlib.Tactic

/-!
# `T_3` Cesàro limit — Exact equality at even `N` (Ch07 extension)

For the involution `T_3` (`T_3² = I`), the **Cesàro average**

  `C_N := (1/N) · ∑_{k=0}^{N-1} T_3^k`

equals the stationary projector `Π = (I + T_3) / 2` **exactly** when `N`
is even, and approaches it as `1/N → 0` for odd `N`.

This file proves the exact equality for the smallest even `N`:

* `N = 2`: `C_2 = (1/2)(I + T_3) = Π` (immediate).
* For any `N = 2k`: `C_{2k} = (1/2)(I + T_3) = Π` (since the sum
  telescopes: `k · I + k · T_3 = k · (I + T_3)`, then divide by `2k`).

The result is the **finite-time Cesàro equality** for an involution —
no asymptotic limit needed.

## Reference

Monograph Chapter 7 §"Limite Cesàro à temps fini",
follow-up to `PT/Stochastic/SpectralDominance.lean`.
-/

namespace PT.Stochastic

open Matrix PT.Sieve

/-! ### Cesàro average for `N = 2` -/

/-- `T_3^0 + T_3^1 = I + T_3 = 2 · stationaryProjector`. -/
theorem T3_pow_sum_two :
    T3^0 + T3^1 = (1 : Matrix (Fin 2) (Fin 2) ℝ) + T3 := by
  simp [pow_zero, pow_one]

/-- The Cesàro average at `N = 2`:
    `(1/2) · (T_3^0 + T_3^1) = stationaryProjector`. -/
theorem T3_cesaro_two :
    ((1 : ℝ) / 2) • (T3^0 + T3^1) = stationaryProjector := by
  rw [T3_pow_sum_two, add_comm]
  exact (stationaryProjector_eq_average).symm

/-! ### Cesàro average for `N = 4` -/

/-- `T_3^0 + T_3^1 + T_3^2 + T_3^3 = 2 · I + 2 · T_3
                                     = 2 · (I + T_3)`. -/
theorem T3_pow_sum_four :
    T3^0 + T3^1 + T3^2 + T3^3
      = (2 : ℝ) • ((1 : Matrix (Fin 2) (Fin 2) ℝ) + T3) := by
  have h2 : T3^2 = 1 := by rw [sq]; exact T3_involution
  have h3 : T3^3 = T3 := by
    have : T3^3 = T3^2 * T3 := by rw [pow_succ]
    rw [this, h2, one_mul]
  rw [pow_zero, pow_one, h2, h3]
  ext i j
  fin_cases i <;> fin_cases j <;> simp [T3]
  all_goals norm_num

/-- The Cesàro average at `N = 4`:
    `(1/4) · (T_3^0 + ... + T_3^3) = (1/2)(I + T_3) = stationaryProjector`. -/
theorem T3_cesaro_four :
    ((1 : ℝ) / 4) • (T3^0 + T3^1 + T3^2 + T3^3) = stationaryProjector := by
  rw [T3_pow_sum_four]
  -- (1/4) · (2 · (I + T_3)) = (1/2) · (I + T_3) = stationaryProjector
  rw [smul_smul]
  have : (1 : ℝ) / 4 * 2 = 1 / 2 := by norm_num
  rw [this, add_comm]
  exact (stationaryProjector_eq_average).symm

/-! ### Generalised even-N Cesàro -/

/-- For any `n ≥ 1`, the Cesàro sum
    `T_3^0 + T_3^1 + ... + T_3^(2n-1)` equals `n · (I + T_3)`. -/
theorem T3_pow_sum_even (n : ℕ) :
    ∑ k ∈ Finset.range (2 * n), T3^k
      = (n : ℝ) • ((1 : Matrix (Fin 2) (Fin 2) ℝ) + T3) := by
  induction n with
  | zero => simp
  | succ k ih =>
    have hstep : 2 * (k + 1) = 2 * k + 2 := by ring
    rw [hstep, Finset.sum_range_succ, Finset.sum_range_succ, ih]
    have heven : T3^(2 * k) = (1 : Matrix (Fin 2) (Fin 2) ℝ) :=
      T3_pow_even k
    have hodd : T3^(2 * k + 1) = T3 := T3_pow_odd k
    rw [heven, hodd]
    push_cast
    -- (k : ℝ) • (I + T_3) + I + T_3 = (k + 1) • (I + T_3)
    rw [add_smul, one_smul]
    abel

/-! ### Headline -/

/-- **Headline (Cesàro convergence at finite time).**

    For the involution `T_3`:

    * The Cesàro average at `N = 2` equals `stationaryProjector` exactly:
      `(1/2)(I + T_3) = Π`.
    * The Cesàro average at `N = 4` equals `stationaryProjector` exactly.
    * In general, for `N = 2n` (even `N ≥ 2`), the sum
      `∑_{k=0}^{2n-1} T_3^k = n · (I + T_3)`, so the Cesàro average is
      `(1/2)(I + T_3) = Π`.

    The convergence to the stationary projector is therefore **exact at
    every even time step**, not just asymptotic. -/
theorem T3_cesaro_summary :
    -- N = 2 exact
    ((1 : ℝ) / 2) • (T3^0 + T3^1) = stationaryProjector
    -- N = 4 exact
    ∧ ((1 : ℝ) / 4) • (T3^0 + T3^1 + T3^2 + T3^3) = stationaryProjector
    -- General even-N formula
    ∧ (∀ n : ℕ, ∑ k ∈ Finset.range (2 * n), T3^k
        = (n : ℝ) • ((1 : Matrix (Fin 2) (Fin 2) ℝ) + T3)) :=
  ⟨T3_cesaro_two, T3_cesaro_four, T3_pow_sum_even⟩

end PT.Stochastic
