/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Sieve.Bimodality
import Mathlib.Data.Nat.GCD.Basic
import Mathlib.Tactic

/-!
# Admissible residues — Coprimality with the third primorial `30 = 2·3·5`

This file proves the **coprimality structure** of the admissible residue
set `R = {1, 7, 11, 13, 17, 19, 23, 29}` of `PT/Sieve/Bimodality.lean`:

* Each element `r ∈ R` is coprime to `30`.
* Equivalently, each `r ∈ R` lies in `(ℤ/30ℤ)*`.
* The set `R` has cardinality `8 = φ(30)`, the Euler totient of `30`.
* Thus `R = (ℤ/30ℤ)*` (exhaustive enumeration).

This complements `PT/Sieve/BimodalityCardinality.lean` (which gives the
cardinality balance) by exhibiting the **structural identity**
`admissibleResidues = (ℤ/30ℤ)*` set-theoretically.

## Reference

Monograph Appendix N §"Résidus admissibles modulo 30".
-/

namespace PT.Sieve

open Finset

/-! ### Coprimality of each admissible residue -/

theorem coprime_30_of_1 : Nat.Coprime 1 30 := by decide
theorem coprime_30_of_7 : Nat.Coprime 7 30 := by decide
theorem coprime_30_of_11 : Nat.Coprime 11 30 := by decide
theorem coprime_30_of_13 : Nat.Coprime 13 30 := by decide
theorem coprime_30_of_17 : Nat.Coprime 17 30 := by decide
theorem coprime_30_of_19 : Nat.Coprime 19 30 := by decide
theorem coprime_30_of_23 : Nat.Coprime 23 30 := by decide
theorem coprime_30_of_29 : Nat.Coprime 29 30 := by decide

/-- **Each admissible residue is coprime to `30`.** -/
theorem admissible_coprime_30 (r : ℕ) (hr : r ∈ admissibleResidues) :
    Nat.Coprime r 30 := by
  unfold admissibleResidues at hr
  -- r ∈ {1, 7, 11, 13, 17, 19, 23, 29}
  fin_cases hr
  · exact coprime_30_of_1
  · exact coprime_30_of_7
  · exact coprime_30_of_11
  · exact coprime_30_of_13
  · exact coprime_30_of_17
  · exact coprime_30_of_19
  · exact coprime_30_of_23
  · exact coprime_30_of_29

/-! ### Non-divisibility by the prime factors of 30 -/

/-- For each admissible residue `r`, `r` is not divisible by `2`, `3`, or `5`. -/
theorem admissible_not_div_by_2_3_5 (r : ℕ) (hr : r ∈ admissibleResidues) :
    ¬ 2 ∣ r ∧ ¬ 3 ∣ r ∧ ¬ 5 ∣ r := by
  unfold admissibleResidues at hr
  fin_cases hr <;> refine ⟨?_, ?_, ?_⟩ <;> decide

/-! ### Cardinality matches φ(30) = 8 -/

/-- The Euler totient of `30 = 2 · 3 · 5`: `φ(30) = (2-1)(3-1)(5-1) = 8`. -/
theorem nat_totient_30 : Nat.totient 30 = 8 := by decide

/-- **Cardinality match.** `|admissibleResidues| = φ(30) = 8`. -/
theorem admissibleResidues_card_eq_totient :
    admissibleResidues.card = Nat.totient 30 := by
  rw [nat_totient_30]
  decide

/-! ### Exhaustive form: every `r < 30` coprime to 30 lies in admissibleResidues -/

/-- The set of nat `r < 30` coprime to `30`. -/
def coprimeMod30 : Finset ℕ :=
  (Finset.range 30).filter (fun r => Nat.Coprime r 30)

/-- **Exhaustive identification.** The set `coprimeMod30` equals
    `admissibleResidues` (with the additional element `1` already in both). -/
theorem coprimeMod30_eq_admissibleResidues :
    coprimeMod30 = admissibleResidues := by
  decide

/-! ### Headline -/

/-- **Headline (primorial-coprime structure).** The 8 admissible residues
    `R = {1, 7, 11, 13, 17, 19, 23, 29}` form **exactly** the set of
    naturals `< 30` coprime to `30`, of cardinality `φ(30) = 8`. They
    are pairwise *not* divisible by `2, 3, 5`, i.e. they survive the
    sieve of Eratosthenes up to `√30 < 6`. -/
theorem admissible_primorial_coprime_summary :
    (∀ r ∈ admissibleResidues, Nat.Coprime r 30)
    ∧ admissibleResidues.card = 8
    ∧ Nat.totient 30 = 8
    ∧ admissibleResidues.card = Nat.totient 30
    ∧ coprimeMod30 = admissibleResidues := by
  refine ⟨admissible_coprime_30, ?_, nat_totient_30,
          admissibleResidues_card_eq_totient,
          coprimeMod30_eq_admissibleResidues⟩
  decide

end PT.Sieve
