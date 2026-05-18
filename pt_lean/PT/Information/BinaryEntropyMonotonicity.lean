/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Information.EntropyOfBinaryDistribution
import Mathlib.Analysis.SpecialFunctions.Log.NegMulLog
import Mathlib.Tactic

/-!
# Monotonicity envelope of binary entropy

The binary entropy `H_bin(p) := -p log p - (1-p) log (1-p)` is symmetric
about `p = 1/2`, where it attains its global maximum `log 2`
(`binEntropy_half`). The map is strictly increasing on `[0, 1/2]` and
strictly decreasing on `[1/2, 1]`. A complete formalisation of the
strict monotonicity requires the derivative `H'(p) = log((1-p)/p)`,
which is non-trivial in Lean.

This file delivers a **realistic** subset of monotonicity results that
follow algebraically from the strict concavity of `negMulLog`
(Mathlib: `Real.strictConcaveOn_negMulLog`) plus elementary log
inequalities. Concretely:

* `binEntropy_lt_log_two_of_ne_half` — for `p ∈ [0, 1]` with `p ≠ 1/2`,
  `H_bin(p) < log 2`. This strengthens `binEntropy_le_log_two` from
  the Bekenstein bound module.
* `binEntropy_third_lt_log_two` — `H_bin(1/3) < log 2` (specialisation).
* `binEntropy_third_pos` — `0 = H_bin(0) < H_bin(1/3)`,
  using the explicit value `log 3 - (2/3) log 2 > 0`.
* `log_three_gt_two_thirds_log_two` — the auxiliary inequality
  `(2/3) log 2 < log 3`, proved via `3 log 3 > 2 log 2` ↔ `27 > 4`.
* `binEntropy_two_thirds_eq_third` — symmetry corollary
  `H_bin(2/3) = H_bin(1/3)`.
* `binary_entropy_comparisons_headline` — packages four explicit
  comparisons at the canonical points `0, 1/3, 1/2, 2/3`.

The strict-monotonicity-on-half-intervals statement is **deliberately
left open**: it is not used downstream in the present scope.

## References

Monograph Chapter 4 §"Bornes tendues" (Bekenstein); see also the
parent module `EntropyOfBinaryDistribution`.
-/

namespace PT.Information

open Real

/-! ### Auxiliary log inequality: `(2/3) log 2 < log 3` -/

/-- **Auxiliary inequality.** `3 · log 3 > 2 · log 2`.

    Equivalent to `log 27 > log 4`, which follows from `27 > 4` and the
    strict monotonicity of `Real.log` on `ℝ>0`. -/
lemma three_log_three_gt_two_log_two :
    2 * Real.log 2 < 3 * Real.log 3 := by
  have h1 : Real.log 4 < Real.log 27 := by
    apply Real.log_lt_log (by norm_num : (0 : ℝ) < 4)
    norm_num
  -- `log 4 = 2 log 2` and `log 27 = 3 log 3`.
  have h4 : Real.log 4 = 2 * Real.log 2 := by
    have : (4 : ℝ) = 2 ^ 2 := by norm_num
    rw [this, Real.log_pow]
    ring
  have h27 : Real.log 27 = 3 * Real.log 3 := by
    have : (27 : ℝ) = 3 ^ 3 := by norm_num
    rw [this, Real.log_pow]
    ring
  rw [h4, h27] at h1
  exact h1

/-- **Headline inequality.** `(2/3) · log 2 < log 3`.

    Direct rescaling of `three_log_three_gt_two_log_two`. -/
lemma log_three_gt_two_thirds_log_two :
    (2 : ℝ) / 3 * Real.log 2 < Real.log 3 := by
  have h := three_log_three_gt_two_log_two
  linarith

/-! ### Strict bound: `H_bin(p) < log 2` for `p ≠ 1/2` -/

/-- **Strict Bekenstein bound on `Fin 2`.**

    For `p ∈ [0, 1]` with `p ≠ 1/2`, `H_bin(p) < log 2`.

    Proof: strict concavity of `negMulLog` on `Set.Ici 0` (Mathlib's
    `Real.strictConcaveOn_negMulLog`) applied to the pair `(p, 1 - p)`
    with weights `(1/2, 1/2)`, whose midpoint is `1/2`. Since `p ≠ 1/2`
    forces `p ≠ 1 - p`, the strict inequality applies. -/
theorem binEntropy_lt_log_two_of_ne_half
    {p : ℝ} (hp0 : 0 ≤ p) (hp1 : p ≤ 1) (hne : p ≠ 1 / 2) :
    binEntropy p < Real.log 2 := by
  rw [binEntropy_eq]
  have h1p_nn : (0 : ℝ) ≤ 1 - p := by linarith
  have hp_mem : p ∈ Set.Ici (0 : ℝ) := hp0
  have h1p_mem : (1 - p) ∈ Set.Ici (0 : ℝ) := h1p_nn
  -- `p ≠ 1 - p` since `p ≠ 1/2`.
  have hpne : p ≠ 1 - p := by
    intro h
    apply hne
    linarith
  have key := Real.strictConcaveOn_negMulLog.2 hp_mem h1p_mem hpne
    (show (0 : ℝ) < 1 / 2 by norm_num)
    (show (0 : ℝ) < 1 / 2 by norm_num)
    (show (1 : ℝ) / 2 + 1 / 2 = 1 by norm_num)
  simp only [smul_eq_mul] at key
  have hmid : (1 : ℝ) / 2 * p + (1 : ℝ) / 2 * (1 - p) = 1 / 2 := by ring
  rw [hmid] at key
  have hneg_half : Real.negMulLog ((1 : ℝ) / 2) = Real.log 2 / 2 := by
    unfold Real.negMulLog
    have hlog : Real.log ((1 : ℝ) / 2) = -Real.log 2 := by
      rw [one_div]; exact Real.log_inv 2
    rw [hlog]; ring
  rw [hneg_half] at key
  linarith

/-- **Specialisation at `p = 1/3`.** `H_bin(1/3) < log 2`. -/
theorem binEntropy_third_lt_log_two :
    binEntropy ((1 : ℝ) / 3) < Real.log 2 := by
  apply binEntropy_lt_log_two_of_ne_half
  · norm_num
  · norm_num
  · norm_num

/-- **Specialisation at `p = 2/3`.** `H_bin(2/3) < log 2`. -/
theorem binEntropy_two_thirds_lt_log_two :
    binEntropy ((2 : ℝ) / 3) < Real.log 2 := by
  apply binEntropy_lt_log_two_of_ne_half
  · norm_num
  · norm_num
  · norm_num

/-! ### Symmetry corollary -/

/-- **Symmetry corollary.** `H_bin(2/3) = H_bin(1/3)`. -/
theorem binEntropy_two_thirds_eq_third :
    binEntropy ((2 : ℝ) / 3) = binEntropy ((1 : ℝ) / 3) := by
  have hsym := binEntropy_symm ((1 : ℝ) / 3)
  -- `binEntropy (1/3) = binEntropy (1 - 1/3) = binEntropy (2/3)`.
  have h23 : (1 : ℝ) - 1 / 3 = 2 / 3 := by norm_num
  rw [h23] at hsym
  exact hsym.symm

/-! ### Strict positivity: `0 < H_bin(1/3)` -/

/-- **Strict positivity at `1/3`.** `H_bin(0) = 0 < H_bin(1/3)`.

    Uses the explicit value `H_bin(1/3) = log 3 - (2/3) log 2` and
    the auxiliary inequality `(2/3) log 2 < log 3`. -/
theorem binEntropy_third_pos :
    binEntropy (0 : ℝ) < binEntropy ((1 : ℝ) / 3) := by
  rw [binEntropy_zero, binEntropy_third]
  linarith [log_three_gt_two_thirds_log_two]

/-- **Strict positivity at `2/3`.** `H_bin(0) = 0 < H_bin(2/3)`. -/
theorem binEntropy_two_thirds_pos :
    binEntropy (0 : ℝ) < binEntropy ((2 : ℝ) / 3) := by
  rw [binEntropy_two_thirds_eq_third]
  exact binEntropy_third_pos

/-! ### Comparisons at `1` -/

/-- **Strict positivity at `1/3` versus `1`.** `H_bin(1) = 0 < H_bin(1/3)`. -/
theorem binEntropy_one_lt_third :
    binEntropy (1 : ℝ) < binEntropy ((1 : ℝ) / 3) := by
  rw [binEntropy_one]
  have := binEntropy_third_pos
  rw [binEntropy_zero] at this
  exact this

/-! ### Strict ordering at the midpoint -/

/-- **Strict ordering: `H_bin(1/3) < H_bin(1/2)`.** -/
theorem binEntropy_third_lt_half :
    binEntropy ((1 : ℝ) / 3) < binEntropy ((1 : ℝ) / 2) := by
  rw [binEntropy_half]
  exact binEntropy_third_lt_log_two

/-- **Strict ordering: `H_bin(2/3) < H_bin(1/2)`.** -/
theorem binEntropy_two_thirds_lt_half :
    binEntropy ((2 : ℝ) / 3) < binEntropy ((1 : ℝ) / 2) := by
  rw [binEntropy_two_thirds_eq_third]
  exact binEntropy_third_lt_half

/-! ### Headline -/

/-- **Headline (binary entropy comparisons at PT-canonical points).**

    Packaging the explicit order envelope at `p ∈ {0, 1/3, 1/2, 2/3, 1}`:

    * `H_bin(0) = H_bin(1) = 0`.
    * `0 < H_bin(1/3)` and `0 < H_bin(2/3)` (strict).
    * `H_bin(2/3) = H_bin(1/3)` (symmetry).
    * `H_bin(1/3) < H_bin(1/2) = log 2` (strict, Bekenstein-tight).
    * `H_bin(2/3) < H_bin(1/2) = log 2` (strict). -/
theorem binary_entropy_comparisons_headline :
    binEntropy (0 : ℝ) < binEntropy ((1 : ℝ) / 3)
    ∧ binEntropy ((1 : ℝ) / 3) < binEntropy ((1 : ℝ) / 2)
    ∧ binEntropy ((2 : ℝ) / 3) = binEntropy ((1 : ℝ) / 3)
    ∧ binEntropy ((2 : ℝ) / 3) < binEntropy ((1 : ℝ) / 2)
    ∧ binEntropy (1 : ℝ) < binEntropy ((1 : ℝ) / 3) :=
  ⟨binEntropy_third_pos, binEntropy_third_lt_half,
   binEntropy_two_thirds_eq_third, binEntropy_two_thirds_lt_half,
   binEntropy_one_lt_third⟩

end PT.Information
