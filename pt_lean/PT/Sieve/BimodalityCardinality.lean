/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Sieve.LegendreLogParity
import PT.Sieve.BimodalityT1Projection
import Mathlib.Data.Finset.Card
import Mathlib.Tactic

/-!
# Bimodality — cardinality balance from `(ℤ/5ℤ)*` character orthogonality

This file extends `PT.Sieve.BimodalityT1Projection` by giving the
**cardinality balance** witness for the bimodality dichotomy as the direct
projection of `(ℤ/5ℤ)*`-character orthogonality:

* `|QR mod 5| = |NR mod 5| = 2` — the two character classes have equal size.
* `|lowResidues| = |highResidues| = 4` — both bimodality classes have equal
  size on the admissible-residue set.
* The cardinality balance `2 · |QR| = |R| / 2 = |lowResidues|` makes the
  factor of `2 = p₁` explicit in the bridge between `(ℤ/5ℤ)*` and the
  `(ℤ/30ℤ)*`-admissible classes.

These are pure finset cardinality identities, proven by `decide`.

## Reference

Appendix N of the monograph, follow-up to `\label{cor:bimodality_T1}`.
Closes audit row #42 (cardinality side).
-/

namespace PT.Sieve

open Finset

/-! ### Cardinality of QR/NR sets in `(ℤ/5ℤ)*` -/

/-- The QR-finset `{1, 4} ⊂ ZMod 5` has cardinality 2. -/
theorem qr_mod5_card : (({1, 4} : Finset (ZMod 5))).card = 2 := by decide

/-- The NR-finset `{2, 3} ⊂ ZMod 5` has cardinality 2. -/
theorem nr_mod5_card : (({2, 3} : Finset (ZMod 5))).card = 2 := by decide

/-- **Character orthogonality (cardinality form).** `|QR| = |NR| = 2`. -/
theorem qr_card_eq_nr_card :
    (({1, 4} : Finset (ZMod 5))).card = (({2, 3} : Finset (ZMod 5))).card := by
  decide

/-- Total number of nonzero residues mod 5 is `4 = |QR| + |NR|`. -/
theorem nonzero_mod5_card :
    (({1, 4} : Finset (ZMod 5))).card + (({2, 3} : Finset (ZMod 5))).card = 4 := by
  decide

/-! ### Admissible residue set splits 4 + 4 -/

/-- The admissible-residue set `R` has `|R| = 8`. -/
theorem admissibleResidues_card : admissibleResidues.card = 8 := by decide

/-- **Cardinality balance corollary.** `|R| = 2 · |lowResidues| = 2 · |highResidues|`. -/
theorem admissibleResidues_card_eq_two_lowcard :
    admissibleResidues.card = 2 * lowResidues.card := by decide

theorem admissibleResidues_card_eq_two_highcard :
    admissibleResidues.card = 2 * highResidues.card := by decide

/-! ### Explicit factor of `p₁ = 2` -/

/-- **Headline.** The admissible-residue set has `|R| = 2 · p₁ · |QR mod 5|`:
    factor `2` for the parity split (low / high), factor `p₁ = 2` for the
    bimodality amplitude, factor `|QR| = 2` for the character classes. -/
theorem admissibleResidues_card_factored :
    admissibleResidues.card = 2 * p1 * (({1, 4} : Finset (ZMod 5))).card := by
  decide

/-! ### Fibre interpretation: `lowResidues` and `highResidues` reduce to QR / NR -/

/-- Casting `lowResidues` modulo 5, each element lands in the QR set `{1, 4}`. -/
theorem lowResidues_image_in_QR (r : ℕ) (hr : r ∈ lowResidues) :
    to5 r ∈ ({1, 4} : Finset (ZMod 5)) := by
  revert r hr
  decide

/-- Casting `highResidues` modulo 5, each element lands in the NR set `{2, 3}`. -/
theorem highResidues_image_in_NR (r : ℕ) (hr : r ∈ highResidues) :
    to5 r ∈ ({2, 3} : Finset (ZMod 5)) := by
  revert r hr
  decide

/-! ### Bimodality completeness as a multiset partition -/

/-- The full `R` partitions into 8 admissible residues, split 4 + 4 between
    the two `δ̄`-values. -/
theorem bimodality_partition_card :
    lowResidues.card + highResidues.card = admissibleResidues.card := by
  decide

/-- **Disjoint union form.** `lowResidues ∪ highResidues = R` and they are
    disjoint, so `|R| = |low| + |high| = 4 + 4 = 8`. -/
theorem bimodality_partition_complete :
    lowResidues ∪ highResidues = admissibleResidues
    ∧ Disjoint lowResidues highResidues
    ∧ lowResidues.card + highResidues.card = 8 := by
  refine ⟨low_union_high_eq_admissible, low_disjoint_high, ?_⟩
  decide

end PT.Sieve
