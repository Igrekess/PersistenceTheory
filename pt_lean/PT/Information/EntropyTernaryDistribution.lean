/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Information.EntropyOfBinaryDistribution
import PT.Information.EntropyBoundsTight
import PT.Information.GFTSpecialMValues
import PT.Information.ShannonEntropyConcavity
import Mathlib.Analysis.SpecialFunctions.Log.NegMulLog
import Mathlib.Tactic

/-!
# Entropy of ternary distributions at PT-canonical points

For reals `p, q`, the **ternary distribution** `terDist p q : Fin 3 → ℝ`
is defined by

* `terDist p q 0 = p`,
* `terDist p q 1 = q`,
* `terDist p q 2 = 1 - p - q`.

Its Shannon entropy is
`H_ter(p, q) := negMulLog p + negMulLog q + negMulLog (1 - p - q)`.

The ternary case is the PT-canonical entropy on the **three active
primes** `{3, 5, 7}` — the next sector after the binary `s = 1/2`
case treated in `EntropyOfBinaryDistribution`.

We give:

* the definition `terDist` and the elementary lemma that it sums to `1`;
* the definition `terEntropy p q := shannonH univ (terDist p q)`;
* the equilibrium identity `terEntropy (1/3) (1/3) = log 3`;
* the **binary degeneracy**: `terEntropy p (1 - p) = binEntropy p`
  (the third coordinate is `0` and `negMulLog 0 = 0`);
* the endpoint value `terEntropy 1 0 = 0` (delta mass at coordinate `0`);
* the **Bekenstein bound** `terEntropy p q ≤ log 3` on the standard
  2-simplex (`p, q ≥ 0`, `p + q ≤ 1`);
* the **symmetry** `terEntropy p q = terEntropy q p` (swap of the first
  two arguments).

## Reference

Monograph Chapter 4 §"Bornes tendues" (Bekenstein on `Fin 3`),
companion to `EntropyOfBinaryDistribution`.
-/

namespace PT.Information

open Real Finset

/-! ### Definition of the ternary distribution -/

/-- The ternary distribution on `Fin 3` with parameters `p, q`. -/
noncomputable def terDist (p q : ℝ) : Fin 3 → ℝ
  | 0 => p
  | 1 => q
  | 2 => 1 - p - q

@[simp] lemma terDist_zero (p q : ℝ) : terDist p q 0 = p := rfl
@[simp] lemma terDist_one  (p q : ℝ) : terDist p q 1 = q := rfl
@[simp] lemma terDist_two  (p q : ℝ) : terDist p q 2 = 1 - p - q := rfl

/-- `terDist p q` sums to `1` over `Fin 3` (no hypothesis needed). -/
lemma terDist_sum_eq_one (p q : ℝ) :
    ∑ i, terDist p q i = 1 := by
  rw [Fin.sum_univ_three]
  simp [terDist]

/-- `terDist p q` is non-negative iff `(p, q)` lies in the standard
    2-simplex `{p ≥ 0, q ≥ 0, p + q ≤ 1}`. -/
lemma terDist_nonneg {p q : ℝ}
    (hp : 0 ≤ p) (hq : 0 ≤ q) (hpq : p + q ≤ 1) :
    ∀ i ∈ (Finset.univ : Finset (Fin 3)), 0 ≤ terDist p q i := by
  intro i _
  fin_cases i
  · simpa using hp
  · simpa using hq
  · simp; linarith

/-! ### Ternary entropy -/

/-- The Shannon entropy of the ternary distribution `terDist p q`. -/
noncomputable def terEntropy (p q : ℝ) : ℝ :=
  shannonH (Finset.univ : Finset (Fin 3)) (terDist p q)

/-- Pointwise expansion of `terEntropy` in terms of `negMulLog`. -/
lemma terEntropy_eq (p q : ℝ) :
    terEntropy p q =
      Real.negMulLog p + Real.negMulLog q + Real.negMulLog (1 - p - q) := by
  unfold terEntropy shannonH
  rw [Fin.sum_univ_three]
  simp [terDist]

/-- Pointwise expansion of `terEntropy` in terms of `Real.log`. -/
lemma terEntropy_eq_log (p q : ℝ) :
    terEntropy p q =
      -(p * Real.log p) - (q * Real.log q)
        - ((1 - p - q) * Real.log (1 - p - q)) := by
  rw [terEntropy_eq]
  unfold Real.negMulLog
  ring

/-! ### Symmetry under swap of the first two arguments -/

/-- **Symmetry (swap of the first two arguments).**
    `H_ter(p, q) = H_ter(q, p)`. -/
theorem terEntropy_swap (p q : ℝ) : terEntropy p q = terEntropy q p := by
  rw [terEntropy_eq, terEntropy_eq]
  have h : (1 : ℝ) - q - p = 1 - p - q := by ring
  rw [h]
  ring

/-! ### Binary degeneracy: `terEntropy p (1 - p) = binEntropy p` -/

/-- **Binary degeneracy.** When the third coordinate is `0` (i.e.
    `q = 1 - p`), the ternary entropy reduces to the binary entropy.
    Uses Mathlib's convention `negMulLog 0 = 0`. -/
theorem terEntropy_degenerate_binary (p : ℝ) :
    terEntropy p (1 - p) = binEntropy p := by
  rw [terEntropy_eq, binEntropy_eq]
  -- 1 - p - (1 - p) = 0, and negMulLog 0 = 0.
  have h : (1 : ℝ) - p - (1 - p) = 0 := by ring
  rw [h]
  simp [Real.negMulLog]

/-! ### Endpoint value: `terEntropy 1 0 = 0` -/

/-- `terEntropy 1 0 = 0` — the delta mass at coordinate `0` has zero
    entropy. -/
@[simp] theorem terEntropy_one_zero : terEntropy 1 0 = 0 := by
  rw [terEntropy_eq]
  -- 1 - 1 - 0 = 0
  have h : (1 : ℝ) - 1 - 0 = 0 := by ring
  rw [h]
  simp [Real.negMulLog]

/-- `terEntropy 0 1 = 0` — the delta mass at coordinate `1`. -/
@[simp] theorem terEntropy_zero_one : terEntropy 0 1 = 0 := by
  rw [terEntropy_eq]
  have h : (1 : ℝ) - 0 - 1 = 0 := by ring
  rw [h]
  simp [Real.negMulLog]

/-- `terEntropy 0 0 = 0` — the delta mass at coordinate `2`. -/
@[simp] theorem terEntropy_zero_zero : terEntropy 0 0 = 0 := by
  rw [terEntropy_eq]
  have h : (1 : ℝ) - 0 - 0 = 1 := by ring
  rw [h]
  simp [Real.negMulLog]

/-! ### Equilibrium: `terEntropy (1/3) (1/3) = log 3` -/

/-- **Equilibrium identity.** `H_ter(1/3, 1/3) = log 3`.

    This is the Bekenstein maximum at the uniform PT-active-prime
    distribution `(1/3, 1/3, 1/3)`. Follows from
    `shannonH_uniform_eq_log` on `Fin 3`. -/
theorem terEntropy_uniform :
    terEntropy ((1 : ℝ) / 3) ((1 : ℝ) / 3) = Real.log 3 := by
  -- At `p = q = 1/3`, `terDist p q = U_3` on `Fin 3`.
  have hext : terDist ((1 : ℝ) / 3) ((1 : ℝ) / 3)
      = (fun _ : Fin 3 => (1 : ℝ) / 3) := by
    funext i
    fin_cases i
    · rfl
    · rfl
    · simp [terDist]; ring
  unfold terEntropy
  rw [hext]
  apply shannonH_uniform_eq_log (Finset.univ : Finset (Fin 3)) 3
  · norm_num
  · simp

/-! ### Bekenstein bound: `terEntropy p q ≤ log 3` on the 2-simplex -/

/-- **Ternary Bekenstein bound.** For `(p, q)` in the standard
    2-simplex (`p, q ≥ 0`, `p + q ≤ 1`), `H_ter(p, q) ≤ log 3`.

    Derived from `shannonH_le_log_card_fin` on `Fin 3`. -/
theorem terEntropy_le_log_three
    {p q : ℝ} (hp : 0 ≤ p) (hq : 0 ≤ q) (hpq : p + q ≤ 1) :
    terEntropy p q ≤ Real.log 3 := by
  unfold terEntropy
  -- Apply the generic Jensen-Bekenstein bound on `Fin 3`.
  have h := shannonH_le_log_card_fin (n := 3) (by norm_num)
    (terDist p q) (terDist_nonneg hp hq hpq)
    (terDist_sum_eq_one p q)
  -- h : shannonH univ (terDist p q) ≤ Real.log ((3 : ℕ) : ℝ)
  have hcast : ((3 : ℕ) : ℝ) = 3 := by norm_num
  rw [hcast] at h
  exact h

/-! ### Headline -/

/-- **Headline (ternary entropy at PT-canonical points).**

    The ternary entropy
    `H_ter(p, q) := -p log p - q log q - (1-p-q) log (1-p-q)`
    on the standard 2-simplex `{p ≥ 0, q ≥ 0, p + q ≤ 1}` satisfies:

    * `H_ter` is symmetric under swap: `H_ter(p, q) = H_ter(q, p)`.
    * Delta corners vanish: `H_ter(1, 0) = H_ter(0, 1) = H_ter(0, 0) = 0`.
    * `H_ter(1/3, 1/3) = log 3` (Bekenstein maximum on `Fin 3`).
    * `H_ter(p, q) ≤ log 3` on the 2-simplex.
    * **Binary degeneracy**: `H_ter(p, 1 - p) = H_bin(p)` — when the
      third coordinate vanishes, ternary reduces to binary entropy. -/
theorem ternary_entropy_headline :
    terEntropy 1 0 = 0
    ∧ terEntropy 0 1 = 0
    ∧ terEntropy 0 0 = 0
    ∧ terEntropy ((1 : ℝ) / 3) ((1 : ℝ) / 3) = Real.log 3
    ∧ (∀ p q, terEntropy p q = terEntropy q p)
    ∧ (∀ p, terEntropy p (1 - p) = binEntropy p) :=
  ⟨terEntropy_one_zero, terEntropy_zero_one, terEntropy_zero_zero,
   terEntropy_uniform, terEntropy_swap, terEntropy_degenerate_binary⟩

end PT.Information
