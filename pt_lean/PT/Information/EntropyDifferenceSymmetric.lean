/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Information.EntropyOfBinaryDistribution
import PT.Information.BinaryEntropyMonotonicity
import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.Tactic

/-!
# Entropy differences at PT-canonical points

For the binary entropy `H_bin(p) := -p log p - (1-p) log (1-p)` we
collect explicit, closed-form values of **differences** at the
PT-canonical points `p ∈ {0, 1/3, 1/2, 2/3, 1}`. Three categories
appear:

* **Trivial difference (symmetry).** `H_bin(p) - H_bin(1-p) = 0`.
  This is the algebraic content of the symmetry
  `binEntropy_symm` packaged as a vanishing difference.
* **Non-trivial difference (Bekenstein gap).**
  `H_bin(1/2) - H_bin(1/3) = (5/3) log 2 - log 3`.
  This is the *strict* Bekenstein gap at the flat T1 prior. By
  `log_three_gt_two_thirds_log_two`, this is positive, so it
  coincides with `|H_bin(1/2) - H_bin(1/3)|`.
* **Cascade difference.**
  `(H_bin(1/3) - H_bin(0)) - (H_bin(1/2) - H_bin(1/3))
     = 2 log 3 - (7/3) log 2 > 0`.
  Equivalently, the *first step* `H(0) → H(1/3)` covers more
  entropy than the *second step* `H(1/3) → H(1/2)`. This is the
  concrete numerical signature of the strict concavity of
  `H_bin`.
* **Lipschitz-style envelope.** For `p, q ∈ [0, 1]`,
  `|H_bin(p) - H_bin(q)| ≤ log 2`. This is the elementary
  envelope deduced from the Bekenstein bound
  `binEntropy_le_log_two` and the non-negativity of `H_bin`
  (`binEntropy_zero ≤ H_bin(p)` is folklore but **not** proved
  here as an inequality; we use the symmetric two-sided form
  via `binEntropy_le_log_two` on both sides).

All proofs rely only on:

* `binEntropy_symm`, `binEntropy_zero`, `binEntropy_half`,
  `binEntropy_third`, `binEntropy_two_thirds_eq_third`
  (parent `EntropyOfBinaryDistribution`);
* `log_three_gt_two_thirds_log_two`,
  `three_log_three_gt_two_log_two` (parent
  `BinaryEntropyMonotonicity`);
* elementary algebra and `Real.log`-arithmetic.

## References

Monograph Chapter 4 §"Bornes tendues" (Bekenstein); see also
`EntropyOfBinaryDistribution` and `BinaryEntropyMonotonicity`.
-/

namespace PT.Information

open Real

/-! ### Trivial symmetric difference -/

/-- **Trivial difference (symmetry).** `H_bin(p) - H_bin(1 - p) = 0`. -/
@[simp] theorem binEntropy_sub_symm (p : ℝ) :
    binEntropy p - binEntropy (1 - p) = 0 := by
  rw [← binEntropy_symm]
  ring

/-- **Specialisation.** `H_bin(1/3) - H_bin(2/3) = 0`. -/
theorem binEntropy_third_sub_two_thirds :
    binEntropy ((1 : ℝ) / 3) - binEntropy ((2 : ℝ) / 3) = 0 := by
  rw [binEntropy_two_thirds_eq_third]; ring

/-- **Specialisation.** `H_bin(0) - H_bin(1) = 0`. -/
theorem binEntropy_zero_sub_one :
    binEntropy (0 : ℝ) - binEntropy (1 : ℝ) = 0 := by
  rw [binEntropy_zero, binEntropy_one]; ring

/-! ### Non-trivial difference: `H_bin(1/2) − H_bin(1/3)` -/

/-- **Bekenstein gap at `1/3`.**

    `H_bin(1/2) - H_bin(1/3) = (5/3) log 2 - log 3`.

    Direct algebraic consequence of `binEntropy_half` and
    `binEntropy_third`:
    `log 2 - (log 3 - (2/3) log 2) = (5/3) log 2 - log 3`. -/
theorem binEntropy_half_sub_third :
    binEntropy ((1 : ℝ) / 2) - binEntropy ((1 : ℝ) / 3)
      = (5 : ℝ) / 3 * Real.log 2 - Real.log 3 := by
  rw [binEntropy_half, binEntropy_third]
  ring

/-- **The Bekenstein gap at `1/3` is positive.**

    `0 < H_bin(1/2) - H_bin(1/3) = (5/3) log 2 - log 3`.

    Proof: from `log_three_gt_two_thirds_log_two`, i.e.
    `(2/3) log 2 < log 3`, combined with `0 < log 2`. We have
    `(5/3) log 2 - log 3 = log 2 + ((2/3) log 2 - log 3) + (log 3 - log 3)`
    ... cleaner: rewrite as `log 2 - (log 3 - (2/3) log 2) > 0`
    iff `log 3 - (2/3) log 2 < log 2`, which is
    `log 3 < (5/3) log 2`. We instead use the direct strict
    comparison `binEntropy_third_lt_half`. -/
theorem binEntropy_half_sub_third_pos :
    0 < binEntropy ((1 : ℝ) / 2) - binEntropy ((1 : ℝ) / 3) := by
  have h := binEntropy_third_lt_half
  linarith

/-- **Absolute Bekenstein gap at `1/3`.**

    `|H_bin(1/2) - H_bin(1/3)| = (5/3) log 2 - log 3`.

    Combines `binEntropy_half_sub_third` with the positivity
    `binEntropy_half_sub_third_pos`. -/
theorem abs_binEntropy_half_sub_third :
    |binEntropy ((1 : ℝ) / 2) - binEntropy ((1 : ℝ) / 3)|
      = (5 : ℝ) / 3 * Real.log 2 - Real.log 3 := by
  rw [abs_of_pos binEntropy_half_sub_third_pos]
  exact binEntropy_half_sub_third

/-! ### Non-trivial difference: `H_bin(1/3) − H_bin(0)` -/

/-- **First-step gap.** `H_bin(1/3) - H_bin(0) = log 3 - (2/3) log 2`. -/
theorem binEntropy_third_sub_zero :
    binEntropy ((1 : ℝ) / 3) - binEntropy (0 : ℝ)
      = Real.log 3 - (2 : ℝ) / 3 * Real.log 2 := by
  rw [binEntropy_third, binEntropy_zero]
  ring

/-- **The first-step gap is positive.**
    `0 < H_bin(1/3) - H_bin(0)`. -/
theorem binEntropy_third_sub_zero_pos :
    0 < binEntropy ((1 : ℝ) / 3) - binEntropy (0 : ℝ) := by
  have h := binEntropy_third_pos
  linarith

/-- **Absolute first-step gap.**
    `|H_bin(1/3) - H_bin(0)| = log 3 - (2/3) log 2`. -/
theorem abs_binEntropy_third_sub_zero :
    |binEntropy ((1 : ℝ) / 3) - binEntropy (0 : ℝ)|
      = Real.log 3 - (2 : ℝ) / 3 * Real.log 2 := by
  rw [abs_of_pos binEntropy_third_sub_zero_pos]
  exact binEntropy_third_sub_zero

/-! ### Cascade difference (concavity signature) -/

/-- **Auxiliary inequality.** `(7/3) log 2 < 2 log 3`.

    Equivalent to `log (2^(7/3)) < log (3^2) = log 9`, i.e.
    `2^7 < 9^3`, i.e. `128 < 729`. We prove the rational form
    `7 log 2 < 6 log 3`, i.e. `log 128 < log 729`. -/
lemma seven_log_two_lt_six_log_three :
    7 * Real.log 2 < 6 * Real.log 3 := by
  have h1 : Real.log 128 < Real.log 729 := by
    apply Real.log_lt_log (by norm_num : (0 : ℝ) < 128)
    norm_num
  have h128 : Real.log 128 = 7 * Real.log 2 := by
    have heq : (128 : ℝ) = 2 ^ 7 := by norm_num
    rw [heq, Real.log_pow]; push_cast; ring
  have h729 : Real.log 729 = 6 * Real.log 3 := by
    have heq : (729 : ℝ) = 3 ^ 6 := by norm_num
    rw [heq, Real.log_pow]; push_cast; ring
  rw [h128, h729] at h1
  exact h1

/-- **Headline form.** `(7/3) log 2 < 2 · log 3`. -/
lemma seven_thirds_log_two_lt_two_log_three :
    (7 : ℝ) / 3 * Real.log 2 < 2 * Real.log 3 := by
  have h := seven_log_two_lt_six_log_three
  linarith

/-- **Cascade difference.**

    `(H_bin(1/2) - H_bin(1/3)) - (H_bin(1/3) - H_bin(0))
        = (7/3) log 2 - 2 log 3`.

    The right-hand side is **negative**: the first step
    `H_bin(0) → H_bin(1/3)` carries strictly more entropy than
    the second step `H_bin(1/3) → H_bin(1/2)`. This is the
    quantitative signature of strict concavity of `H_bin` on
    `[0, 1/2]`. -/
theorem binEntropy_cascade_diff :
    (binEntropy ((1 : ℝ) / 2) - binEntropy ((1 : ℝ) / 3))
        - (binEntropy ((1 : ℝ) / 3) - binEntropy (0 : ℝ))
      = (7 : ℝ) / 3 * Real.log 2 - 2 * Real.log 3 := by
  rw [binEntropy_half, binEntropy_third, binEntropy_zero]
  ring

/-- **Cascade difference is negative.**

    `(H_bin(1/2) - H_bin(1/3)) - (H_bin(1/3) - H_bin(0)) < 0`,
    equivalently
    `|H_bin(1/2) - H_bin(1/3)| < |H_bin(1/3) - H_bin(0)|`. -/
theorem binEntropy_cascade_diff_neg :
    (binEntropy ((1 : ℝ) / 2) - binEntropy ((1 : ℝ) / 3))
        - (binEntropy ((1 : ℝ) / 3) - binEntropy (0 : ℝ)) < 0 := by
  rw [binEntropy_cascade_diff]
  linarith [seven_thirds_log_two_lt_two_log_three]

/-- **Cascade difference of absolute values is negative.**

    Since both gaps are positive (`binEntropy_half_sub_third_pos`,
    `binEntropy_third_sub_zero_pos`), the cascade difference of
    absolute values equals the raw cascade difference:
    `|H_bin(1/2) - H_bin(1/3)| - |H_bin(1/3) - H_bin(0)|
       = (7/3) log 2 - 2 log 3 < 0`. -/
theorem abs_binEntropy_cascade_diff :
    |binEntropy ((1 : ℝ) / 2) - binEntropy ((1 : ℝ) / 3)|
        - |binEntropy ((1 : ℝ) / 3) - binEntropy (0 : ℝ)|
      = (7 : ℝ) / 3 * Real.log 2 - 2 * Real.log 3 := by
  rw [abs_binEntropy_half_sub_third, abs_binEntropy_third_sub_zero]
  ring

/-! ### Lipschitz-like envelope on `[0, 1]` -/

/-- **Non-negativity of `binEntropy` on `[0, 1]`.**

    For `p ∈ [0, 1]`, `0 ≤ H_bin(p)`.

    Proof: `H_bin(p) = negMulLog p + negMulLog (1 - p)`, and each
    summand is non-negative on `[0, 1]` since `negMulLog x = -x log x`
    with `log x ≤ 0` for `x ∈ [0, 1]` (and `x ≥ 0`). -/
theorem binEntropy_nonneg {p : ℝ} (hp0 : 0 ≤ p) (hp1 : p ≤ 1) :
    0 ≤ binEntropy p := by
  rw [binEntropy_eq_log]
  -- show `0 ≤ -(p log p) - ((1-p) log (1-p))`
  have h1p0 : (0 : ℝ) ≤ 1 - p := by linarith
  have h1p1 : (1 - p) ≤ 1 := by linarith
  -- log p ≤ 0 if 0 ≤ p ≤ 1 (with the convention log 0 = 0).
  have hlog_p : Real.log p ≤ 0 := by
    rcases eq_or_lt_of_le hp0 with hp_eq | hp_pos
    · rw [← hp_eq]; simp
    · exact Real.log_nonpos hp0 hp1
  have hlog_1p : Real.log (1 - p) ≤ 0 := by
    rcases eq_or_lt_of_le h1p0 with h1p_eq | h1p_pos
    · rw [← h1p_eq]; simp
    · exact Real.log_nonpos h1p0 h1p1
  -- so p log p ≤ 0 and (1-p) log (1-p) ≤ 0.
  have hp_log_p : p * Real.log p ≤ 0 :=
    mul_nonpos_iff.mpr (Or.inl ⟨hp0, hlog_p⟩)
  have h1p_log_1p : (1 - p) * Real.log (1 - p) ≤ 0 :=
    mul_nonpos_iff.mpr (Or.inl ⟨h1p0, hlog_1p⟩)
  linarith

/-- **Lipschitz-style envelope.**

    For `p, q ∈ [0, 1]`,
    `|H_bin(p) - H_bin(q)| ≤ log 2`.

    Proof: both `H_bin(p)` and `H_bin(q)` lie in `[0, log 2]`
    (combining `binEntropy_nonneg` with `binEntropy_le_log_two`),
    hence their difference lies in `[-log 2, log 2]`. -/
theorem abs_binEntropy_sub_le_log_two
    {p q : ℝ} (hp0 : 0 ≤ p) (hp1 : p ≤ 1)
    (hq0 : 0 ≤ q) (hq1 : q ≤ 1) :
    |binEntropy p - binEntropy q| ≤ Real.log 2 := by
  have hpL : binEntropy p ≤ Real.log 2 := binEntropy_le_log_two hp0 hp1
  have hpN : 0 ≤ binEntropy p := binEntropy_nonneg hp0 hp1
  have hqL : binEntropy q ≤ Real.log 2 := binEntropy_le_log_two hq0 hq1
  have hqN : 0 ≤ binEntropy q := binEntropy_nonneg hq0 hq1
  rw [abs_sub_le_iff]
  refine ⟨?_, ?_⟩
  · linarith
  · linarith

/-- **Saturated envelope at the canonical extremes.**

    The Lipschitz-style envelope is **tight**: it is attained at
    the pair `(p, q) = (1/2, 0)`, where
    `|H_bin(1/2) - H_bin(0)| = log 2`. -/
theorem abs_binEntropy_half_sub_zero :
    |binEntropy ((1 : ℝ) / 2) - binEntropy (0 : ℝ)| = Real.log 2 := by
  rw [binEntropy_half, binEntropy_zero]
  have : (Real.log 2 - 0 : ℝ) = Real.log 2 := by ring
  rw [this, abs_of_nonneg]
  exact Real.log_nonneg (by norm_num)

/-! ### Headline -/

/-- **Headline (entropy differences at PT-canonical points).**

    The binary entropy `H_bin` has the following closed-form
    differences at the PT-canonical points `p ∈ {0, 1/3, 1/2}`:

    * **Symmetric (trivial) difference.** `H_bin(p) - H_bin(1-p) = 0`.
    * **Bekenstein gap at `1/3`.**
      `H_bin(1/2) - H_bin(1/3) = (5/3) log 2 - log 3 > 0`.
    * **First-step gap.**
      `H_bin(1/3) - H_bin(0) = log 3 - (2/3) log 2 > 0`.
    * **Cascade difference (concavity signature).**
      `(H_bin(1/2) - H_bin(1/3)) - (H_bin(1/3) - H_bin(0))
        = (7/3) log 2 - 2 log 3 < 0`,
      so the first step carries strictly more entropy than
      the second.
    * **Saturated envelope.**
      `|H_bin(1/2) - H_bin(0)| = log 2`. -/
theorem binary_entropy_differences_headline :
    (∀ p, binEntropy p - binEntropy (1 - p) = 0)
    ∧ binEntropy ((1 : ℝ) / 2) - binEntropy ((1 : ℝ) / 3)
        = (5 : ℝ) / 3 * Real.log 2 - Real.log 3
    ∧ binEntropy ((1 : ℝ) / 3) - binEntropy (0 : ℝ)
        = Real.log 3 - (2 : ℝ) / 3 * Real.log 2
    ∧ (binEntropy ((1 : ℝ) / 2) - binEntropy ((1 : ℝ) / 3))
          - (binEntropy ((1 : ℝ) / 3) - binEntropy (0 : ℝ))
        = (7 : ℝ) / 3 * Real.log 2 - 2 * Real.log 3
    ∧ |binEntropy ((1 : ℝ) / 2) - binEntropy (0 : ℝ)| = Real.log 2 :=
  ⟨binEntropy_sub_symm,
   binEntropy_half_sub_third,
   binEntropy_third_sub_zero,
   binEntropy_cascade_diff,
   abs_binEntropy_half_sub_zero⟩

end PT.Information
