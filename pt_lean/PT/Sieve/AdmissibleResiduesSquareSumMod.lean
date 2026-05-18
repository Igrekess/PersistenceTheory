/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Sieve.Bimodality
import PT.Sieve.AdmissibleResiduesArithmetic
import PT.Sieve.CoprimeAdmissibilityProduct
import Mathlib.Tactic

/-!
# Admissible residues — Sum of squares modulo various `q` (App N extension)

This file collects the **sum of squares invariants** of the admissible
residue set `R = {1, 7, 11, 13, 17, 19, 23, 29} = (ℤ/30ℤ)*`:

$$\sum_{r \in R} r^2 \;=\; 1 + 49 + 121 + 169 + 289 + 361 + 529 + 841
                   \;=\; 2360.$$

Companion to:
* `AdmissibleResiduesArithmetic.lean` — additive `∑ r = 120`.
* `CoprimeAdmissibilityProduct.lean` — multiplicative `∏ r`, squares mod 30.

## Main results

* `admissibleResidues_sumSq`: `∑ r² = 2360`.
* `admissibleResidues_sumSq_mod_*`: closed-form residue of `2360` modulo
  every divisor of `60 = 2² · 3 · 5` and modulo `8` (the orthogonal
  factor coming from `(ℤ/8ℤ)*`).
* `admissibleResidues_sumSq_variance_index`: the "index of dispersion"
  identity `∑ r² − (∑ r)² / |R| = 2360 − 1800 = 560`.
* `admissibleResidues_sumSq_summary`: headline table of residues `2360 mod q`.

## Patterns (analogues of the Wilson collapses)

The residue `2360 mod q` exhibits the following "collapse" pattern across
the divisors of `60`:

| `q`  |  `2360 mod q` |  Note                                  |
|------|---------------|----------------------------------------|
|  2   |  0            |  all `r` odd ⇒ all `r²` odd, sum of 8 odds is even (in fact ≡ 0 mod 8). |
|  3   |  2            |  Legendre count: 8 residues, all coprime to 3, squares ∈ {1}. |
|  5   |  0            |  squares mod 5 hit `{1, 4}` evenly (4+4), sum `≡ 4·1+4·4 = 20 ≡ 0`. |
|  6   |  2            |  CRT(2,3): (0, 2) ↦ 2. |
|  8   |  0            |  all `r²` odd ≡ 1 mod 8 (units of ℤ/8), 8·1 = 8 ≡ 0 mod 8. |
|  10  |  0            |  CRT(2,5): (0, 0). |
|  15  |  5            |  CRT(3,5): (2, 0) ↦ 5. |
|  30  |  20           |  CRT(2,3,5): (0,2,0) ↦ 20. |
|  60  |  20           |  CRT(4,3,5) with 2360 = 4·590. |

The **`q = 8`** collapse `Σ r² ≡ 0 mod 8` is the squares-mod-8 analogue
of the Wilson collapse: every `r ∈ R` is odd, hence `r² ≡ 1 mod 8`,
and `|R| = 8` cancels the residue.

## Reference

Monograph Appendix N §"Invariants arithmétiques de R" (square-sum part).
-/

namespace PT.Sieve

open Finset

/-! ### The sum of squares: exact value -/

/-- The sum of squares of the admissible residues equals `2360`. -/
theorem admissibleResidues_sumSq :
    admissibleResidues.sum (fun r => r * r) = 2360 := by
  decide

/-- The sum of squares, term-by-term decomposition (for documentation). -/
theorem admissibleResidues_sumSq_terms :
    admissibleResidues.sum (fun r => r * r)
      = 1 * 1 + 7 * 7 + 11 * 11 + 13 * 13
        + 17 * 17 + 19 * 19 + 23 * 23 + 29 * 29 := by
  decide

/-- The sum of squares using `r ^ 2`. -/
theorem admissibleResidues_sumSq_pow :
    admissibleResidues.sum (fun r => r ^ 2) = 2360 := by
  decide

/-! ### Residues of `2360` modulo divisors of `60` (and `8`) -/

/-- `2360 ≡ 0 mod 2` (even). -/
theorem admissibleResidues_sumSq_mod_2 :
    admissibleResidues.sum (fun r => r * r) % 2 = 0 := by
  decide

/-- `2360 ≡ 2 mod 3`. -/
theorem admissibleResidues_sumSq_mod_3 :
    admissibleResidues.sum (fun r => r * r) % 3 = 2 := by
  decide

/-- `2360 ≡ 0 mod 5`. -/
theorem admissibleResidues_sumSq_mod_5 :
    admissibleResidues.sum (fun r => r * r) % 5 = 0 := by
  decide

/-- `2360 ≡ 2 mod 6`. -/
theorem admissibleResidues_sumSq_mod_6 :
    admissibleResidues.sum (fun r => r * r) % 6 = 2 := by
  decide

/-- **Wilson-like collapse mod 8.** `2360 ≡ 0 mod 8`: every admissible
    `r` is odd, so `r² ≡ 1 mod 8`; summing 8 copies of `1` gives
    `8 ≡ 0 mod 8`. -/
theorem admissibleResidues_sumSq_mod_8 :
    admissibleResidues.sum (fun r => r * r) % 8 = 0 := by
  decide

/-- `2360 ≡ 0 mod 10`. -/
theorem admissibleResidues_sumSq_mod_10 :
    admissibleResidues.sum (fun r => r * r) % 10 = 0 := by
  decide

/-- `2360 ≡ 5 mod 15`. -/
theorem admissibleResidues_sumSq_mod_15 :
    admissibleResidues.sum (fun r => r * r) % 15 = 5 := by
  decide

/-- `2360 ≡ 20 mod 30`. -/
theorem admissibleResidues_sumSq_mod_30 :
    admissibleResidues.sum (fun r => r * r) % 30 = 20 := by
  decide

/-- `2360 ≡ 20 mod 60`. Note `2360 = 39 · 60 + 20`. -/
theorem admissibleResidues_sumSq_mod_60 :
    admissibleResidues.sum (fun r => r * r) % 60 = 20 := by
  decide

/-! ### Variance index: `∑ r² − (∑ r)² / |R|` -/

/-- **Variance-of-index identity.** With `∑ r = 120` and `|R| = 8`:
    `∑ r² − (∑ r)² / |R| = 2360 − 14400 / 8 = 2360 − 1800 = 560`.

    This `560` is the (non-normalised) numerator of the population
    variance of `R`: `Var(R) = 560 / 8 = 70`. -/
theorem admissibleResidues_sumSq_variance_index :
    admissibleResidues.sum (fun r => r * r)
      - (admissibleResidues.sum id) ^ 2 / admissibleResidues.card
      = 560 := by
  decide

/-- Explicit form of the variance numerator: `560 = 2³ · 70 = 2³ · 2 · 5 · 7`. -/
theorem admissibleResidues_variance_factored :
    (560 : ℕ) = 8 * 70 ∧ (560 : ℕ) = 2 * 2 * 2 * 2 * 5 * 7 := by
  refine ⟨?_, ?_⟩ <;> decide

/-- The (integer) population variance of `R` equals `70`. -/
theorem admissibleResidues_population_variance :
    (admissibleResidues.sum (fun r => r * r)
      - (admissibleResidues.sum id) ^ 2 / admissibleResidues.card)
      / admissibleResidues.card = 70 := by
  decide

/-! ### Headline -/

/-- **Headline (sum-of-squares invariants).** The admissible residue set
    `R = {1, 7, 11, 13, 17, 19, 23, 29}` satisfies:

    * `∑ r² = 2360`.
    * Modular table:
      `2360 ≡ 0 mod 2`, `≡ 2 mod 3`, `≡ 0 mod 5`,
      `≡ 2 mod 6`, `≡ 0 mod 8`, `≡ 0 mod 10`,
      `≡ 5 mod 15`, `≡ 20 mod 30`, `≡ 20 mod 60`.
    * Variance-of-index: `∑ r² − (∑ r)² / |R| = 560`,
      population variance `Var(R) = 70`.

    The collapses `mod 2`, `mod 5`, `mod 8`, `mod 10` echo the
    Wilson-style cancellation of `(ℤ/30ℤ)*` (every `r` is a unit, and
    `|R| = 8` divides the residue in each case). -/
theorem admissibleResidues_sumSq_summary :
    admissibleResidues.sum (fun r => r * r) = 2360
    ∧ admissibleResidues.sum (fun r => r * r) % 2 = 0
    ∧ admissibleResidues.sum (fun r => r * r) % 3 = 2
    ∧ admissibleResidues.sum (fun r => r * r) % 5 = 0
    ∧ admissibleResidues.sum (fun r => r * r) % 6 = 2
    ∧ admissibleResidues.sum (fun r => r * r) % 8 = 0
    ∧ admissibleResidues.sum (fun r => r * r) % 10 = 0
    ∧ admissibleResidues.sum (fun r => r * r) % 15 = 5
    ∧ admissibleResidues.sum (fun r => r * r) % 30 = 20
    ∧ admissibleResidues.sum (fun r => r * r) % 60 = 20
    ∧ admissibleResidues.sum (fun r => r * r)
        - (admissibleResidues.sum id) ^ 2 / admissibleResidues.card
        = 560 :=
  ⟨admissibleResidues_sumSq,
   admissibleResidues_sumSq_mod_2,
   admissibleResidues_sumSq_mod_3,
   admissibleResidues_sumSq_mod_5,
   admissibleResidues_sumSq_mod_6,
   admissibleResidues_sumSq_mod_8,
   admissibleResidues_sumSq_mod_10,
   admissibleResidues_sumSq_mod_15,
   admissibleResidues_sumSq_mod_30,
   admissibleResidues_sumSq_mod_60,
   admissibleResidues_sumSq_variance_index⟩

end PT.Sieve
