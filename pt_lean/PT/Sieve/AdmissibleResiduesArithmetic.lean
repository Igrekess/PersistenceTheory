/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Sieve.Bimodality
import PT.Sieve.PrimorialCoprime
import Mathlib.Tactic

/-!
# Admissible residues — Arithmetic invariants (App N extension)

This file collects **arithmetic invariants** of the admissible residue set
`R = {1, 7, 11, 13, 17, 19, 23, 29}`:

* Sum: `1 + 7 + 11 + 13 + 17 + 19 + 23 + 29 = 120 = 4 · 30`.
  (The sum equals `4 · primorial 3 = (|R| / 2) · primorial 3`.)
* All elements are odd.
* All elements are coprime to `30 = primorial 3`.
* All elements satisfy `r ≡ ±1 mod 6` (since they are coprime to `6`).
* The set is **symmetric** about `15`: `r ↔ 30 - r` is an involution.

## Reference

Monograph Appendix N §"Invariants arithmétiques de R".
-/

namespace PT.Sieve

open Finset

/-! ### Sum and product -/

/-- The sum of admissible residues equals `120 = 4 · 30`. -/
theorem admissibleResidues_sum :
    admissibleResidues.sum id = 120 := by
  decide

/-- The sum factors as `4 · 30 = (|R| / 2) · primorial 3`. -/
theorem admissibleResidues_sum_factored :
    admissibleResidues.sum id = 4 * 30 := by
  decide

/-- The product of admissible residues. -/
theorem admissibleResidues_prod :
    admissibleResidues.prod id = 1 * 7 * 11 * 13 * 17 * 19 * 23 * 29 := by
  decide

/-! ### Parity: every admissible residue is odd -/

theorem admissible_all_odd (r : ℕ) (hr : r ∈ admissibleResidues) :
    r % 2 = 1 := by
  unfold admissibleResidues at hr
  fin_cases hr <;> decide

/-! ### Coprime to 6 (corollary of coprime to 30) -/

theorem admissible_coprime_6 (r : ℕ) (hr : r ∈ admissibleResidues) :
    Nat.Coprime r 6 := by
  unfold admissibleResidues at hr
  fin_cases hr <;> decide

/-! ### Symmetry: `r ↔ 30 - r` -/

/-- The set `R` is **symmetric** about `15`: `1 + 29 = 30`, `7 + 23 = 30`,
    `11 + 19 = 30`, `13 + 17 = 30`. -/
theorem admissible_symmetric_pairs :
    (1 + 29 : ℕ) = 30
    ∧ (7 + 23 : ℕ) = 30
    ∧ (11 + 19 : ℕ) = 30
    ∧ (13 + 17 : ℕ) = 30 := by
  refine ⟨?_, ?_, ?_, ?_⟩ <;> decide

/-- The involution `r ↦ 30 - r` maps `R` to itself. -/
theorem admissible_involution_30_sub :
    (30 - 1 : ℕ) ∈ admissibleResidues
    ∧ (30 - 7 : ℕ) ∈ admissibleResidues
    ∧ (30 - 11 : ℕ) ∈ admissibleResidues
    ∧ (30 - 13 : ℕ) ∈ admissibleResidues
    ∧ (30 - 17 : ℕ) ∈ admissibleResidues
    ∧ (30 - 19 : ℕ) ∈ admissibleResidues
    ∧ (30 - 23 : ℕ) ∈ admissibleResidues
    ∧ (30 - 29 : ℕ) ∈ admissibleResidues := by
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩ <;> decide

/-! ### Mean equals 15 = μ* -/

/-- The arithmetic mean of `R` equals `120 / 8 = 15 = μ*`. -/
theorem admissibleResidues_mean :
    admissibleResidues.sum id / admissibleResidues.card = 15 := by
  rw [admissibleResidues_sum]; decide

/-! ### Headline -/

/-- **Headline (arithmetic invariants).** The 8 admissible residues
    `R = {1, 7, 11, 13, 17, 19, 23, 29}` satisfy:

    * Sum = `120 = 4 · 30 = (|R|/2) · primorial 3`.
    * Mean = `15 = μ*`.
    * All elements are odd.
    * All elements are coprime to `30` (hence coprime to `6` as well).
    * The set is symmetric about `15` via the involution `r ↦ 30 - r`. -/
theorem admissibleResidues_arithmetic_summary :
    admissibleResidues.sum id = 120
    ∧ admissibleResidues.sum id = 4 * 30
    ∧ admissibleResidues.sum id / admissibleResidues.card = 15
    ∧ (∀ r ∈ admissibleResidues, r % 2 = 1)
    ∧ (∀ r ∈ admissibleResidues, Nat.Coprime r 6)
    ∧ (1 + 29 : ℕ) = 30
    ∧ (7 + 23 : ℕ) = 30
    ∧ (11 + 19 : ℕ) = 30
    ∧ (13 + 17 : ℕ) = 30 :=
  ⟨admissibleResidues_sum,
   admissibleResidues_sum_factored,
   admissibleResidues_mean,
   admissible_all_odd,
   admissible_coprime_6,
   admissible_symmetric_pairs.1,
   admissible_symmetric_pairs.2.1,
   admissible_symmetric_pairs.2.2.1,
   admissible_symmetric_pairs.2.2.2⟩

end PT.Sieve
