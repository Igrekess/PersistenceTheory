/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Stochastic.T3CesaroLimit
import PT.Stochastic.T2T30Normalisation
import PT.Stochastic.CRTFactorizationT2
import Mathlib.Tactic

/-!
# `T_2 ⊗ T_3` Cesàro limit — Finite-time exact equality on the 2D Kronecker block

This file extends `PT/Stochastic/T3CesaroLimit.lean` to the *Kronecker
product* `T_2 ⊗ T_3`, where `T_2 = (1)` is the trivial `1 × 1` block (see
`T2T30Normalisation`). Since `T_2² = I` (trivially) and `T_3² = I`
(`T3_involution`), the Kronecker product is again an involution
(`kronecker_involution`):

  `(T_2 ⊗ T_3)² = T_2² ⊗ T_3² = I ⊗ I = I.`

Hence the **Cesàro average** of `T_2 ⊗ T_3` equals the *stationary
projector* `Π_{T_2 ⊗ T_3} := (1/2)(I + T_2 ⊗ T_3)` **exactly** at every
even `N`, in direct analogy with `T3_cesaro_two`, `T3_cesaro_four`, and
`T3_pow_sum_even`.

## Main results

* `T2T3_involution` — the Kronecker product is an involution.
* `T2T3_pow_even`, `T2T3_pow_odd` — orbit alternates between `I` and
  `T_2 ⊗ T_3`.
* `T2T3_stationaryProjector := (1/2) • (I + T_2 ⊗ T_3)` — the Cesàro
  projector for the Kronecker dynamics.
* `T2T3_cesaro_two` — `(1/2) • ((T_2 ⊗ T_3)^0 + (T_2 ⊗ T_3)^1) = Π`.
* `T2T3_cesaro_four` — `(1/4) • ∑_{k<4} (T_2 ⊗ T_3)^k = Π`.
* `T2T3_pow_sum_even` — general even-`N` formula
  `∑_{k<2n} (T_2 ⊗ T_3)^k = n • (I + T_2 ⊗ T_3)`.

## Reference

Monograph Chapter 7 §"Limite Cesàro à temps fini" (Kronecker extension),
follow-up to `PT.Stochastic.T3CesaroLimit` and
`PT.Stochastic.T2T30Normalisation`.
-/

namespace PT.Stochastic

open Matrix Kronecker BigOperators PT.Sieve

/-! ### Involution of `T_2 ⊗ T_3` -/

/-- **`T_2` is trivially an involution.** Since `T_2 = (1)` is the
    `1 × 1` identity matrix, `T_2 · T_2 = 1 · 1 = 1 = I`. -/
theorem T2_trivial_involution :
    T2_trivial * T2_trivial = (1 : Matrix (Fin 1) (Fin 1) ℝ) := by
  simp [T2_trivial]

/-- **`T_2 ⊗ T_3` is an involution.** Direct consequence of
    `kronecker_involution` applied to `T_2² = I` and `T_3² = I`. -/
theorem T2T3_involution :
    (T2_trivial ⊗ₖ T3) * (T2_trivial ⊗ₖ T3)
      = (1 : Matrix (Fin 1 × Fin 2) (Fin 1 × Fin 2) ℝ) :=
  kronecker_involution T2_trivial T3 T2_trivial_involution T3_involution

/-! ### Powers of the Kronecker involution -/

/-- **Even powers of `T_2 ⊗ T_3`.** For every `n : ℕ`,
    `(T_2 ⊗ T_3)^(2n) = I`. -/
theorem T2T3_pow_even (n : ℕ) :
    (T2_trivial ⊗ₖ T3) ^ (2 * n)
      = (1 : Matrix (Fin 1 × Fin 2) (Fin 1 × Fin 2) ℝ) := by
  induction n with
  | zero => simp
  | succ k ih =>
    have hstep : 2 * (k + 1) = 2 * k + 2 := by ring
    have hsq : (T2_trivial ⊗ₖ T3) ^ 2
                 = (1 : Matrix (Fin 1 × Fin 2) (Fin 1 × Fin 2) ℝ) := by
      rw [sq]; exact T2T3_involution
    rw [hstep, pow_add, ih, one_mul, hsq]

/-- **Odd powers of `T_2 ⊗ T_3`.** For every `n : ℕ`,
    `(T_2 ⊗ T_3)^(2n+1) = T_2 ⊗ T_3`. -/
theorem T2T3_pow_odd (n : ℕ) :
    (T2_trivial ⊗ₖ T3) ^ (2 * n + 1) = T2_trivial ⊗ₖ T3 := by
  rw [pow_add, T2T3_pow_even, pow_one, one_mul]

/-! ### Stationary projector of `T_2 ⊗ T_3` -/

/-- **Stationary projector of `T_2 ⊗ T_3`.**
    `Π_{T_2 ⊗ T_3} := (1/2) • (I + T_2 ⊗ T_3)`. Because `T_2 ⊗ T_3` is
    an involution, this is the Cesàro / Perron-Frobenius limit (also
    achieved at finite time, see below). -/
noncomputable def T2T3_stationaryProjector :
    Matrix (Fin 1 × Fin 2) (Fin 1 × Fin 2) ℝ :=
  ((1 : ℝ) / 2) • ((1 : Matrix (Fin 1 × Fin 2) (Fin 1 × Fin 2) ℝ)
                    + T2_trivial ⊗ₖ T3)

/-- Defining equation of `T2T3_stationaryProjector`. -/
theorem T2T3_stationaryProjector_def :
    T2T3_stationaryProjector
      = ((1 : ℝ) / 2) • ((1 : Matrix (Fin 1 × Fin 2) (Fin 1 × Fin 2) ℝ)
                          + T2_trivial ⊗ₖ T3) := rfl

/-! ### Cesàro average for `N = 2` -/

/-- `(T_2 ⊗ T_3)^0 + (T_2 ⊗ T_3)^1 = I + T_2 ⊗ T_3`. -/
theorem T2T3_pow_sum_two :
    (T2_trivial ⊗ₖ T3)^0 + (T2_trivial ⊗ₖ T3)^1
      = (1 : Matrix (Fin 1 × Fin 2) (Fin 1 × Fin 2) ℝ) + T2_trivial ⊗ₖ T3 := by
  simp [pow_zero, pow_one]

/-- **Cesàro average at `N = 2`.**
    `(1/2) • ((T_2 ⊗ T_3)^0 + (T_2 ⊗ T_3)^1) = T2T3_stationaryProjector`. -/
theorem T2T3_cesaro_two :
    ((1 : ℝ) / 2) • ((T2_trivial ⊗ₖ T3)^0 + (T2_trivial ⊗ₖ T3)^1)
      = T2T3_stationaryProjector := by
  rw [T2T3_pow_sum_two, T2T3_stationaryProjector_def]

/-! ### Cesàro average for `N = 4` -/

/-- `(T_2 ⊗ T_3)^0 + (T_2 ⊗ T_3)^1 + (T_2 ⊗ T_3)^2 + (T_2 ⊗ T_3)^3
    = 2 • (I + T_2 ⊗ T_3)`. -/
theorem T2T3_pow_sum_four :
    (T2_trivial ⊗ₖ T3)^0 + (T2_trivial ⊗ₖ T3)^1
      + (T2_trivial ⊗ₖ T3)^2 + (T2_trivial ⊗ₖ T3)^3
      = (2 : ℝ) • ((1 : Matrix (Fin 1 × Fin 2) (Fin 1 × Fin 2) ℝ)
                    + T2_trivial ⊗ₖ T3) := by
  have h2 : (T2_trivial ⊗ₖ T3)^2
              = (1 : Matrix (Fin 1 × Fin 2) (Fin 1 × Fin 2) ℝ) := by
    rw [sq]; exact T2T3_involution
  have h3 : (T2_trivial ⊗ₖ T3)^3 = T2_trivial ⊗ₖ T3 := by
    have hstep : (T2_trivial ⊗ₖ T3)^3 = (T2_trivial ⊗ₖ T3)^2 * (T2_trivial ⊗ₖ T3) := by
      rw [pow_succ]
    rw [hstep, h2, one_mul]
  rw [pow_zero, pow_one, h2, h3]
  -- LHS: I + T + I + T = 2 • (I + T)
  rw [two_smul]
  abel

/-- **Cesàro average at `N = 4`.**
    `(1/4) • ((T_2 ⊗ T_3)^0 + ... + (T_2 ⊗ T_3)^3) = T2T3_stationaryProjector`. -/
theorem T2T3_cesaro_four :
    ((1 : ℝ) / 4) • ((T2_trivial ⊗ₖ T3)^0 + (T2_trivial ⊗ₖ T3)^1
                      + (T2_trivial ⊗ₖ T3)^2 + (T2_trivial ⊗ₖ T3)^3)
      = T2T3_stationaryProjector := by
  rw [T2T3_pow_sum_four, T2T3_stationaryProjector_def, smul_smul]
  have hrw : (1 : ℝ) / 4 * 2 = 1 / 2 := by norm_num
  rw [hrw]

/-! ### Generalised even-`N` Cesàro -/

/-- **Even-`N` Cesàro sum.** For any `n : ℕ`,
    `∑_{k<2n} (T_2 ⊗ T_3)^k = n • (I + T_2 ⊗ T_3)`. -/
theorem T2T3_pow_sum_even (n : ℕ) :
    ∑ k ∈ Finset.range (2 * n), (T2_trivial ⊗ₖ T3)^k
      = (n : ℝ) • ((1 : Matrix (Fin 1 × Fin 2) (Fin 1 × Fin 2) ℝ)
                    + T2_trivial ⊗ₖ T3) := by
  induction n with
  | zero => simp
  | succ k ih =>
    have hstep : 2 * (k + 1) = 2 * k + 2 := by ring
    rw [hstep, Finset.sum_range_succ, Finset.sum_range_succ, ih]
    have heven : (T2_trivial ⊗ₖ T3)^(2 * k)
                   = (1 : Matrix (Fin 1 × Fin 2) (Fin 1 × Fin 2) ℝ) :=
      T2T3_pow_even k
    have hodd : (T2_trivial ⊗ₖ T3)^(2 * k + 1) = T2_trivial ⊗ₖ T3 :=
      T2T3_pow_odd k
    rw [heven, hodd]
    push_cast
    -- (k : ℝ) • (I + T) + I + T = (k + 1) • (I + T)
    rw [add_smul, one_smul]
    abel

/-- **Cesàro at any even `N = 2n` (with `n ≥ 1`) equals the projector.**
    Dividing by `2n` gives exactly `(1/2)(I + T_2 ⊗ T_3) = Π`. -/
theorem T2T3_cesaro_even (n : ℕ) (hn : 0 < n) :
    ((1 : ℝ) / (2 * n)) • (∑ k ∈ Finset.range (2 * n), (T2_trivial ⊗ₖ T3)^k)
      = T2T3_stationaryProjector := by
  rw [T2T3_pow_sum_even, T2T3_stationaryProjector_def, smul_smul]
  have hnpos : (n : ℝ) ≠ 0 := by exact_mod_cast hn.ne'
  have : (1 : ℝ) / (2 * n) * (n : ℝ) = 1 / 2 := by
    field_simp
  rw [this]

/-! ### Headline -/

/-- **Headline (Cesàro convergence at finite time, Kronecker extension).**

    For the Kronecker involution `T_2 ⊗ T_3`:

    * `(T_2 ⊗ T_3)² = I` — it is an involution.
    * The Cesàro average at `N = 2` equals `T2T3_stationaryProjector` exactly.
    * The Cesàro average at `N = 4` equals `T2T3_stationaryProjector` exactly.
    * For any `N = 2n` (even `N ≥ 2`), the sum
      `∑_{k<2n} (T_2 ⊗ T_3)^k = n • (I + T_2 ⊗ T_3)`, so the Cesàro
      average equals `(1/2)(I + T_2 ⊗ T_3) = Π`.

    The convergence to the stationary projector of the Kronecker dynamics
    is therefore **exact at every even time step**, in direct analogy with
    the single-factor case `T_3`. -/
theorem T2T3_cesaro_summary :
    -- Involution
    (T2_trivial ⊗ₖ T3) * (T2_trivial ⊗ₖ T3)
        = (1 : Matrix (Fin 1 × Fin 2) (Fin 1 × Fin 2) ℝ)
    -- N = 2 exact
    ∧ ((1 : ℝ) / 2) • ((T2_trivial ⊗ₖ T3)^0 + (T2_trivial ⊗ₖ T3)^1)
        = T2T3_stationaryProjector
    -- N = 4 exact
    ∧ ((1 : ℝ) / 4) • ((T2_trivial ⊗ₖ T3)^0 + (T2_trivial ⊗ₖ T3)^1
                        + (T2_trivial ⊗ₖ T3)^2 + (T2_trivial ⊗ₖ T3)^3)
        = T2T3_stationaryProjector
    -- General even-N formula
    ∧ (∀ n : ℕ, ∑ k ∈ Finset.range (2 * n), (T2_trivial ⊗ₖ T3)^k
        = (n : ℝ) • ((1 : Matrix (Fin 1 × Fin 2) (Fin 1 × Fin 2) ℝ)
                      + T2_trivial ⊗ₖ T3)) :=
  ⟨T2T3_involution, T2T3_cesaro_two, T2T3_cesaro_four, T2T3_pow_sum_even⟩

end PT.Stochastic
