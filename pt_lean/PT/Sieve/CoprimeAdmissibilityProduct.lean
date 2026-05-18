/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Sieve.Bimodality
import PT.Sieve.AdmissibleResiduesArithmetic
import PT.Sieve.PrimorialCoprime
import PT.Sieve.TotientCascadeIdentities
import Mathlib.Tactic

/-!
# Admissible residues — Multiplicative invariants (App N extension)

This file collects **multiplicative invariants** of the admissible residue
set `R = {1, 7, 11, 13, 17, 19, 23, 29} = (ℤ/30ℤ)*`.

Companion to `AdmissibleResiduesArithmetic.lean` (additive invariants).

## Main results

* `admissibleResidues_prod_value`: `∏ r ∈ R, r = 215656441`.
* `admissibleResidues_prod_mod_30`: `(∏ r ∈ R, r) ≡ 1 mod 30`.
  This is the Wilson-like identity: the product of the units of a finite
  abelian group is the product of its involutions (here equal to `1`).
* `inverse_pair_*`: explicit inverse pairs in `(ℤ/30ℤ)*`:
  `7·13 ≡ 1`, `17·23 ≡ 1` (cross pairs);
  `1·1 ≡ 1`, `11·11 ≡ 1`, `19·19 ≡ 1`, `29·29 ≡ 1` (auto-inverses).
* `selfInverses_mod_30`: exactly 4 elements of `R` are their own inverses
  mod 30, namely `{1, 11, 19, 29}` (these are the involutions of
  `(ℤ/30ℤ)* ≅ (ℤ/2ℤ)³`).
* `square_mod_30_dichotomy`: for every `r ∈ R`, `r² mod 30 ∈ {1, 19}`;
  the image of squaring is `{1, 19}` (a subgroup of order 2).

## Reference

Monograph Appendix N §"Invariants arithmétiques de R" (multiplicative part).
-/

namespace PT.Sieve

open Finset

/-! ### Product of admissible residues -/

/-- The product of admissible residues is exactly `215656441`. -/
theorem admissibleResidues_prod_value :
    admissibleResidues.prod id = 215656441 := by
  decide

/-- The product of admissible residues modulo `30` equals `1`.
    (Wilson-like identity for `(ℤ/30ℤ)*`.) -/
theorem admissibleResidues_prod_mod_30 :
    admissibleResidues.prod id % 30 = 1 := by
  decide

/-- The product of admissible residues, computed in `ZMod 30`, equals `1`. -/
theorem admissibleResidues_prod_zmod_30 :
    ((admissibleResidues.prod id : ℕ) : ZMod 30) = 1 := by
  decide

/-! ### Inverse pairs in `(ℤ/30ℤ)*` -/

/-- `7 · 13 = 91 ≡ 1 mod 30`: the inverse of `7` in `(ℤ/30ℤ)*` is `13`. -/
theorem inverse_pair_7_13 : (7 * 13 : ℕ) % 30 = 1 := by decide

/-- `17 · 23 = 391 ≡ 1 mod 30`: the inverse of `17` in `(ℤ/30ℤ)*` is `23`. -/
theorem inverse_pair_17_23 : (17 * 23 : ℕ) % 30 = 1 := by decide

/-- `1 · 1 ≡ 1 mod 30`: `1` is its own inverse. -/
theorem inverse_self_1 : (1 * 1 : ℕ) % 30 = 1 := by decide

/-- `11² = 121 ≡ 1 mod 30`: `11` is its own inverse. -/
theorem inverse_self_11 : (11 * 11 : ℕ) % 30 = 1 := by decide

/-- `19² = 361 ≡ 1 mod 30`: `19` is its own inverse. -/
theorem inverse_self_19 : (19 * 19 : ℕ) % 30 = 1 := by decide

/-- `29² = 841 ≡ 1 mod 30`: `29` is its own inverse. -/
theorem inverse_self_29 : (29 * 29 : ℕ) % 30 = 1 := by decide

/-! ### Existence of inverse in `R` for every `r ∈ R` -/

/-- Every admissible residue `r ∈ R` has a multiplicative inverse `s ∈ R`
    modulo `30`. -/
theorem admissible_has_inverse_in_R (r : ℕ) (hr : r ∈ admissibleResidues) :
    ∃ s ∈ admissibleResidues, (r * s) % 30 = 1 := by
  unfold admissibleResidues at hr
  fin_cases hr
  · exact ⟨1, by decide, by decide⟩
  · exact ⟨13, by decide, by decide⟩
  · exact ⟨11, by decide, by decide⟩
  · exact ⟨7, by decide, by decide⟩
  · exact ⟨23, by decide, by decide⟩
  · exact ⟨19, by decide, by decide⟩
  · exact ⟨17, by decide, by decide⟩
  · exact ⟨29, by decide, by decide⟩

/-! ### Self-inverses (involutions of `(ℤ/30ℤ)*`) -/

/-- The set of self-inverses of `R` modulo `30`. -/
def selfInverses : Finset ℕ := {1, 11, 19, 29}

/-- The 4 self-inverse residues form a subset of `R`. -/
theorem selfInverses_subset :
    selfInverses ⊆ admissibleResidues := by
  decide

/-- The cardinality of the self-inverse set is `4 = 2³ / 2 = φ(30) / 2`. -/
theorem selfInverses_card : selfInverses.card = 4 := by decide

/-- Characterisation: `r ∈ R` is self-inverse mod `30` iff `r ∈ {1, 11, 19, 29}`. -/
theorem mem_selfInverses_iff (r : ℕ) (hr : r ∈ admissibleResidues) :
    (r * r) % 30 = 1 ↔ r ∈ selfInverses := by
  unfold admissibleResidues at hr
  fin_cases hr <;> simp [selfInverses]

/-- Every element of `selfInverses` squares to `1` mod `30`. -/
theorem selfInverses_sq_one (r : ℕ) (hr : r ∈ selfInverses) :
    (r * r) % 30 = 1 := by
  unfold selfInverses at hr
  fin_cases hr <;> decide

/-! ### Squares in `(ℤ/30ℤ)*` -/

/-- For every `r ∈ R`, the square `r² mod 30` lies in `{1, 19}`. -/
theorem square_mod_30_dichotomy (r : ℕ) (hr : r ∈ admissibleResidues) :
    (r * r) % 30 = 1 ∨ (r * r) % 30 = 19 := by
  unfold admissibleResidues at hr
  fin_cases hr <;> decide

/-- Explicit square table for `R`:
    `1²≡1, 7²≡19, 11²≡1, 13²≡19, 17²≡19, 19²≡1, 23²≡19, 29²≡1` (all mod `30`). -/
theorem square_table_mod_30 :
    (1 * 1 : ℕ) % 30 = 1
    ∧ (7 * 7 : ℕ) % 30 = 19
    ∧ (11 * 11 : ℕ) % 30 = 1
    ∧ (13 * 13 : ℕ) % 30 = 19
    ∧ (17 * 17 : ℕ) % 30 = 19
    ∧ (19 * 19 : ℕ) % 30 = 1
    ∧ (23 * 23 : ℕ) % 30 = 19
    ∧ (29 * 29 : ℕ) % 30 = 1 := by
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩ <;> decide

/-- The image of the squaring map `r ↦ r² mod 30` restricted to `R` is `{1, 19}`. -/
theorem squares_image :
    (admissibleResidues.image (fun r => (r * r) % 30)) = {1, 19} := by
  decide

/-! ### Headline -/

/-- **Headline (multiplicative invariants).** The admissible residue set
    `R = {1, 7, 11, 13, 17, 19, 23, 29} = (ℤ/30ℤ)*` satisfies:

    * `∏ r ∈ R, r = 215656441 ≡ 1 mod 30` (Wilson-like for `(ℤ/30ℤ)*`).
    * Every `r ∈ R` has a multiplicative inverse `s ∈ R` mod `30`.
    * Exactly `4` elements are self-inverse: `{1, 11, 19, 29}` (involutions).
    * The set of squares `{r² mod 30 : r ∈ R}` is `{1, 19}` (order-2 subgroup). -/
theorem admissibleResidues_multiplicative_summary :
    admissibleResidues.prod id = 215656441
    ∧ admissibleResidues.prod id % 30 = 1
    ∧ selfInverses.card = 4
    ∧ selfInverses ⊆ admissibleResidues
    ∧ (admissibleResidues.image (fun r => (r * r) % 30)) = {1, 19} :=
  ⟨admissibleResidues_prod_value,
   admissibleResidues_prod_mod_30,
   selfInverses_card,
   selfInverses_subset,
   squares_image⟩

end PT.Sieve
