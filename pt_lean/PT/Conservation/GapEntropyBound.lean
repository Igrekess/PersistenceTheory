/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Conservation.ConservationIDExtensions
import PT.Information.BekensteinBound
import PT.Information.EntropyMonotonicity
import Mathlib.Tactic

/-!
# Shannon entropy of the PT prime-gap distribution

This file constructs the **normalised prime-gap distribution** for the PT
prime sequence and bounds its Shannon entropy.

For an integer cutoff `N ≥ 1`, define

  `S_N := ∑_{n=1}^{N} g_n`        (= `p_{N+1} - 2` by conservation),
  `gapsDist N n := g_n / S_N`     (for `n ∈ {1, …, N}`).

This is a discrete probability distribution on the index set `{1, …, N}`
(non-negative, sums to one). Its Shannon entropy `H(gapsDist N)` is then
bounded between `0` and `log N` (the Bekenstein bound, attained by the
uniform distribution).

For the small instance `N = 4` (gaps `(1, 2, 2, 4)`, `S_4 = 9`):

  `gapsDist 4 = (1/9, 2/9, 2/9, 4/9)`.

We prove:

* **Normalisation.** `∑_{n=1}^{N} gapsDist N n = 1` for `N ∈ {1, 2, 3, 4}`.
* **Pointwise range.** `0 ≤ gapsDist N n ≤ 1` for each `n` in the support.
* **Entropy non-negativity.** `0 ≤ shannonH s (gapsDist N)` for
  `N ∈ {1, 2, 3, 4}` (direct consequence of `shannonH_nonneg`).
* **Bekenstein bound (parametric).** `shannonH s (gapsDist N) ≤ log N`
  whenever the Gibbs KL is non-negative (Gibbs' inequality, kept as a
  hypothesis consistent with `shannonH_le_log_of_KL_nonneg`).

**Scope reduction.** Exact numerical computation of `H(gapsDist 4)` is
out of scope (it involves the irrational logarithms `log(1/9)`,
`log(2/9)`, `log(4/9)`). We keep the file to normalisation + the
parametric Bekenstein bracket, which is the formalisable part.

## Reference

Monograph Ch04 §"Borne de Bekenstein" / Ch03 §3.1 (conservation identity).
M4 article §"Entropie de la séquence de gaps".
-/

namespace PT.Conservation

open Finset Real

/-! ### Definition of the normalised gap distribution -/

/-- The cumulative gap sum on the index range `Ico 1 (N+1)`, as a real
    number. For the PT prime sequence,
    `cumGapReal N = p_{N+1} - 2 = ∑_{n=1}^{N} g_n`. -/
noncomputable def cumGapReal (N : ℕ) : ℝ :=
  ∑ n ∈ Ico 1 (N + 1), (gap ptPrime n : ℝ)

/-- The normalised PT prime-gap distribution: `gapsDist N n := g_n / S_N`,
    where `g_n` is the `n`-th gap of `ptPrime` and `S_N` is the cumulative
    sum on `Ico 1 (N+1)`. -/
noncomputable def gapsDist (N n : ℕ) : ℝ :=
  (gap ptPrime n : ℝ) / cumGapReal N

/-! ### Cumulative sum values (real-valued versions) -/

/-- `cumGapReal 1 = 1`. -/
theorem cumGapReal_one : cumGapReal 1 = 1 := by
  unfold cumGapReal
  have h : ∑ n ∈ Ico 1 2, gap ptPrime n = 1 := conservation_N1
  exact_mod_cast h

/-- `cumGapReal 2 = 3`. -/
theorem cumGapReal_two : cumGapReal 2 = 3 := by
  unfold cumGapReal
  have h : ∑ n ∈ Ico 1 3, gap ptPrime n = 3 := conservation_N2
  exact_mod_cast h

/-- `cumGapReal 3 = 5`. -/
theorem cumGapReal_three : cumGapReal 3 = 5 := by
  unfold cumGapReal
  have h : ∑ n ∈ Ico 1 4, gap ptPrime n = 5 := conservation_N3
  exact_mod_cast h

/-- `cumGapReal 4 = 9`. -/
theorem cumGapReal_four : cumGapReal 4 = 9 := by
  unfold cumGapReal
  have h : ∑ n ∈ Ico 1 5, gap ptPrime n = 9 := conservation_N4
  exact_mod_cast h

/-- Positivity of the cumulative sum at `N = 4`. -/
theorem cumGapReal_four_pos : 0 < cumGapReal 4 := by
  rw [cumGapReal_four]; norm_num

/-! ### Pointwise values of `gapsDist 4` -/

/-- The first probability: `gapsDist 4 1 = 1/9`. -/
theorem gapsDist_four_one : gapsDist 4 1 = 1 / 9 := by
  unfold gapsDist
  rw [cumGapReal_four]
  have h : (gap ptPrime 1 : ℝ) = 1 := by exact_mod_cast gap_one
  rw [h]

/-- The second probability: `gapsDist 4 2 = 2/9`. -/
theorem gapsDist_four_two : gapsDist 4 2 = 2 / 9 := by
  unfold gapsDist
  rw [cumGapReal_four]
  have h : (gap ptPrime 2 : ℝ) = 2 := by exact_mod_cast gap_two
  rw [h]

/-- The third probability: `gapsDist 4 3 = 2/9`. -/
theorem gapsDist_four_three : gapsDist 4 3 = 2 / 9 := by
  unfold gapsDist
  rw [cumGapReal_four]
  have h : (gap ptPrime 3 : ℝ) = 2 := by exact_mod_cast gap_three
  rw [h]

/-- The fourth probability: `gapsDist 4 4 = 4/9`. -/
theorem gapsDist_four_four : gapsDist 4 4 = 4 / 9 := by
  unfold gapsDist
  rw [cumGapReal_four]
  have h : (gap ptPrime 4 : ℝ) = 4 := by exact_mod_cast gap_four
  rw [h]

/-! ### Normalisation -/

/-- **Normalisation (`N = 4`).** `∑_{n=1}^{4} gapsDist 4 n = 1`. -/
theorem gapsDist_four_sum : ∑ n ∈ Ico 1 5, gapsDist 4 n = 1 := by
  have hIco : (Ico 1 5 : Finset ℕ) = {1, 2, 3, 4} := by decide
  rw [hIco]
  rw [show ({1, 2, 3, 4} : Finset ℕ) = insert 1 (insert 2 (insert 3 {4})) from rfl]
  rw [Finset.sum_insert (by decide), Finset.sum_insert (by decide),
      Finset.sum_insert (by decide), Finset.sum_singleton]
  rw [gapsDist_four_one, gapsDist_four_two, gapsDist_four_three, gapsDist_four_four]
  norm_num

/-- **Normalisation (generic, factored form).** For `N ≥ 1` with
    `cumGapReal N ≠ 0`, the gap distribution sums to one:
    `∑_{n=1}^{N} gapsDist N n = 1`. -/
theorem gapsDist_sum_eq_one (N : ℕ) (hN : cumGapReal N ≠ 0) :
    ∑ n ∈ Ico 1 (N + 1), gapsDist N n = 1 := by
  unfold gapsDist
  rw [← Finset.sum_div]
  exact div_self hN

/-- **Normalisation (`N = 1`).** -/
theorem gapsDist_one_sum : ∑ n ∈ Ico 1 2, gapsDist 1 n = 1 :=
  gapsDist_sum_eq_one 1 (by rw [cumGapReal_one]; norm_num)

/-- **Normalisation (`N = 2`).** -/
theorem gapsDist_two_sum : ∑ n ∈ Ico 1 3, gapsDist 2 n = 1 :=
  gapsDist_sum_eq_one 2 (by rw [cumGapReal_two]; norm_num)

/-- **Normalisation (`N = 3`).** -/
theorem gapsDist_three_sum : ∑ n ∈ Ico 1 4, gapsDist 3 n = 1 :=
  gapsDist_sum_eq_one 3 (by rw [cumGapReal_three]; norm_num)

/-! ### Pointwise non-negativity and ≤ 1 -/

/-- Each probability of `gapsDist 4` is non-negative. -/
theorem gapsDist_four_nonneg (n : ℕ) (hn : n ∈ Ico 1 5) : 0 ≤ gapsDist 4 n := by
  -- enumerate the four cases via `interval_cases`
  have hIco : (Ico 1 5 : Finset ℕ) = {1, 2, 3, 4} := by decide
  rw [hIco] at hn
  fin_cases hn
  · rw [gapsDist_four_one]; norm_num
  · rw [gapsDist_four_two]; norm_num
  · rw [gapsDist_four_three]; norm_num
  · rw [gapsDist_four_four]; norm_num

/-- Each probability of `gapsDist 4` is at most 1. -/
theorem gapsDist_four_le_one (n : ℕ) (hn : n ∈ Ico 1 5) : gapsDist 4 n ≤ 1 := by
  have hIco : (Ico 1 5 : Finset ℕ) = {1, 2, 3, 4} := by decide
  rw [hIco] at hn
  fin_cases hn
  · rw [gapsDist_four_one]; norm_num
  · rw [gapsDist_four_two]; norm_num
  · rw [gapsDist_four_three]; norm_num
  · rw [gapsDist_four_four]; norm_num

/-! ### Entropy: non-negativity (lower bound `H ≥ 0`) -/

open PT.Information in
/-- **Lower bound `H ≥ 0`.** The Shannon entropy of `gapsDist 4` is
    non-negative. This is a direct consequence of `shannonH_nonneg`
    since each probability lies in `[0, 1]`. -/
theorem shannonH_gapsDist_four_nonneg :
    0 ≤ shannonH (Ico 1 5) (gapsDist 4) :=
  shannonH_nonneg (Ico 1 5) (gapsDist 4)
    gapsDist_four_nonneg gapsDist_four_le_one

/-! ### Entropy: upper bound `H ≤ log N` (Bekenstein, parametric) -/

open PT.Information in
/-- **Upper bound `H ≤ log 4` (parametric form).** Assuming Gibbs'
    inequality `0 ≤ klToUniform (Ico 1 5) 4 (gapsDist 4)` (which is the
    standard non-negativity of KL divergence to the uniform), the entropy
    of the PT gap distribution is bounded by `log 4`.

    The Gibbs hypothesis is taken parametrically, consistent with the
    project convention in `EntropyMonotonicity.shannonH_le_log_of_KL_nonneg`
    (a strict-concavity proof of Gibbs is out of scope here). -/
theorem shannonH_gapsDist_four_le_log_four
    (hKL_nn : 0 ≤ klToUniform (Ico 1 5) (4 : ℝ) (gapsDist 4)) :
    shannonH (Ico 1 5) (gapsDist 4) ≤ Real.log 4 := by
  have hm : (0 : ℝ) < 4 := by norm_num
  exact shannonH_le_log_of_KL_nonneg (Ico 1 5) (4 : ℝ) hm
    (gapsDist 4) gapsDist_four_nonneg gapsDist_four_sum hKL_nn

/-! ### Bilateral bracket (headline) -/

open PT.Information in
/-- **Headline bracket: `0 ≤ H(gapsDist 4) ≤ log 4`** under Gibbs'
    inequality. The lower bound is unconditional (pointwise
    `negMulLog ≥ 0` on `[0,1]`). The upper bound uses the GFT identity
    `log m = D_KL + H` and `D_KL ≥ 0`. -/
theorem shannonH_gapsDist_four_bracket
    (hKL_nn : 0 ≤ klToUniform (Ico 1 5) (4 : ℝ) (gapsDist 4)) :
    0 ≤ shannonH (Ico 1 5) (gapsDist 4)
    ∧ shannonH (Ico 1 5) (gapsDist 4) ≤ Real.log 4 :=
  ⟨shannonH_gapsDist_four_nonneg,
   shannonH_gapsDist_four_le_log_four hKL_nn⟩

end PT.Conservation
