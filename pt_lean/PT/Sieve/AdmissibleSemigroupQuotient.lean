/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Sieve.Bimodality
import PT.Sieve.PrimorialCoprime
import PT.Sieve.CoprimeAdmissibilityProduct
import PT.Sieve.AdmissibleSquareGroupStructure
import PT.Sieve.KleinFourEmbedding
import Mathlib.Tactic

/-!
# Admissible residues — Quotient `R / selfInverses ≅ ℤ/2ℤ` (App N)

This file formalises the **quotient group** `R / selfInverses`, where
`R = (ℤ/30ℤ)* = {1, 7, 11, 13, 17, 19, 23, 29}` (cardinality 8) and
`selfInverses = {1, 11, 19, 29}` is the Klein-four subgroup of involutions
(cardinality 4, index 2).

By Lagrange, `|R / selfInverses| = 8 / 4 = 2`, hence the quotient is
isomorphic to `ℤ/2ℤ`. The two cosets are:

* `selfInverses = {1, 11, 19, 29}` — the **identity** coset of the quotient.
* `R \ selfInverses = {7, 13, 17, 23}` — the **non-identity** coset.

## Strategy

We encode the quotient by a Boolean classifier `quotientClass : ℕ → Bool`
that returns `true` iff `r ∈ selfInverses`. The two values of `quotientClass`
on `R` then index the two cosets.

The multiplicative law on `R / selfInverses ≅ ℤ/2ℤ` translates, under this
encoding, to the rule

  `quotientClass (a · b mod 30) = (quotientClass a == quotientClass b)`

(`true` ↔ same coset, `false` ↔ different coset), which is exactly
`XNOR` on `Bool`. Equivalently, with the additive identification
`true ↦ 0, false ↦ 1`, the rule is addition mod 2.

## Main results

* `quotientClass` — the Bool classifier `r ∈ selfInverses ↦ true`.
* `quotientClass_true_set` / `quotientClass_false_set` — the two cosets
  explicitly: `{1, 11, 19, 29}` and `{7, 13, 17, 23}`.
* `quotientClass_true_card` / `quotientClass_false_card` — each coset has
  cardinality `4`.
* `quotientClass_partition` — the two cosets partition `R`.
* `quotientClass_mul` — multiplicative stability: for `a, b ∈ R`,
  `quotientClass ((a·b) % 30) = (quotientClass a == quotientClass b)`.
* `quotientClass_addMod2` — additive recasting via the map `Bool → ZMod 2`,
  `true ↦ 0, false ↦ 1`: the map is a group homomorphism `R → ZMod 2`.
* `quotient_card_eq_2` — the quotient `R / selfInverses` has 2 classes.
* `quotient_iso_Z2` — bundled isomorphism `R / selfInverses ≅ ℤ/2ℤ`.

## Reference

Monograph Appendix N §"Invariants arithmétiques de R" — subgroup lattice
`{1} < {1,19} < {1,11,19,29} < R`. The quotient `R / selfInverses ≅ ℤ/2`
is the top-level factor in this chain.
-/

namespace PT.Sieve

open Finset

/-! ### The Boolean classifier `quotientClass` -/

/-- The Boolean classifier of the quotient `R / selfInverses`:
    `quotientClass r = true ↔ r ∈ selfInverses`.

    On `R`, this distinguishes the two cosets of `selfInverses`:
    * `true` ↔ `{1, 11, 19, 29}` (identity coset);
    * `false` ↔ `{7, 13, 17, 23}` (non-identity coset). -/
def quotientClass (r : ℕ) : Bool :=
  decide (r ∈ selfInverses)

/-- Explicit values of `quotientClass` on `selfInverses`. -/
theorem quotientClass_on_selfInverses :
    quotientClass 1 = true
    ∧ quotientClass 11 = true
    ∧ quotientClass 19 = true
    ∧ quotientClass 29 = true := by
  refine ⟨?_, ?_, ?_, ?_⟩ <;> decide

/-- Explicit values of `quotientClass` on `R \ selfInverses = {7, 13, 17, 23}`. -/
theorem quotientClass_on_complement :
    quotientClass 7 = false
    ∧ quotientClass 13 = false
    ∧ quotientClass 17 = false
    ∧ quotientClass 23 = false := by
  refine ⟨?_, ?_, ?_, ?_⟩ <;> decide

/-! ### The two cosets, explicit -/

/-- The **identity coset** of the quotient: the set of `r ∈ R` with
    `quotientClass r = true`. Equals `selfInverses = {1, 11, 19, 29}`. -/
def quotientClass_true_set : Finset ℕ := {1, 11, 19, 29}

/-- The **non-identity coset** of the quotient: the set of `r ∈ R` with
    `quotientClass r = false`. Equals `R \ selfInverses = {7, 13, 17, 23}`. -/
def quotientClass_false_set : Finset ℕ := {7, 13, 17, 23}

/-- The identity coset coincides with `selfInverses`. -/
theorem quotientClass_true_set_eq_selfInverses :
    quotientClass_true_set = selfInverses := by decide

/-- The non-identity coset is the complement of `selfInverses` in `R`. -/
theorem quotientClass_false_set_eq_complement :
    quotientClass_false_set = admissibleResidues \ selfInverses := by decide

/-- Filter form: the identity coset is exactly the subset of `R` where
    `quotientClass` evaluates to `true`. -/
theorem quotientClass_true_set_filter :
    admissibleResidues.filter (fun r => quotientClass r = true)
      = quotientClass_true_set := by decide

/-- Filter form: the non-identity coset is exactly the subset of `R` where
    `quotientClass` evaluates to `false`. -/
theorem quotientClass_false_set_filter :
    admissibleResidues.filter (fun r => quotientClass r = false)
      = quotientClass_false_set := by decide

/-! ### Cardinalities -/

/-- The identity coset has cardinality `4`. -/
theorem quotientClass_true_card : quotientClass_true_set.card = 4 := by decide

/-- The non-identity coset has cardinality `4`. -/
theorem quotientClass_false_card : quotientClass_false_set.card = 4 := by decide

/-- The two cosets are disjoint. -/
theorem quotientClass_cosets_disjoint :
    Disjoint quotientClass_true_set quotientClass_false_set := by decide

/-- **Partition of `R`.** The two cosets partition `R`:
    `R = {1, 11, 19, 29} ⊔ {7, 13, 17, 23}`. -/
theorem quotientClass_partition :
    quotientClass_true_set ∪ quotientClass_false_set = admissibleResidues := by
  decide

/-! ### Multiplicative stability — verification table -/

/-- `7 · 13 = 91 ≡ 1 (mod 30) ∈ selfInverses` (non-id · non-id = id). -/
theorem quotientClass_mul_7_13 :
    quotientClass ((7 * 13 : ℕ) % 30)
      = (quotientClass 7 == quotientClass 13) := by decide

/-- `7 · 17 = 119 ≡ 29 (mod 30) ∈ selfInverses` (non-id · non-id = id). -/
theorem quotientClass_mul_7_17 :
    quotientClass ((7 * 17 : ℕ) % 30)
      = (quotientClass 7 == quotientClass 17) := by decide

/-- `7 · 23 = 161 ≡ 11 (mod 30) ∈ selfInverses` (non-id · non-id = id). -/
theorem quotientClass_mul_7_23 :
    quotientClass ((7 * 23 : ℕ) % 30)
      = (quotientClass 7 == quotientClass 23) := by decide

/-- `13 · 17 = 221 ≡ 11 (mod 30) ∈ selfInverses` (non-id · non-id = id). -/
theorem quotientClass_mul_13_17 :
    quotientClass ((13 * 17 : ℕ) % 30)
      = (quotientClass 13 == quotientClass 17) := by decide

/-- `13 · 23 = 299 ≡ 29 (mod 30) ∈ selfInverses` (non-id · non-id = id). -/
theorem quotientClass_mul_13_23 :
    quotientClass ((13 * 23 : ℕ) % 30)
      = (quotientClass 13 == quotientClass 23) := by decide

/-- `17 · 23 = 391 ≡ 1 (mod 30) ∈ selfInverses` (non-id · non-id = id). -/
theorem quotientClass_mul_17_23 :
    quotientClass ((17 * 23 : ℕ) % 30)
      = (quotientClass 17 == quotientClass 23) := by decide

/-- `11 · 7 = 77 ≡ 17 (mod 30) ∈ R \ selfInverses` (id · non-id = non-id). -/
theorem quotientClass_mul_11_7 :
    quotientClass ((11 * 7 : ℕ) % 30)
      = (quotientClass 11 == quotientClass 7) := by decide

/-- `19 · 13 = 247 ≡ 7 (mod 30) ∈ R \ selfInverses` (id · non-id = non-id). -/
theorem quotientClass_mul_19_13 :
    quotientClass ((19 * 13 : ℕ) % 30)
      = (quotientClass 19 == quotientClass 13) := by decide

/-- **Multiplicative stability of `quotientClass`.** For all `a, b ∈ R`,

      `quotientClass ((a · b) mod 30) = (quotientClass a == quotientClass b)`

    in `Bool`. This is exactly the rule of `ℤ/2ℤ` under the identification
    `true ↦ 0, false ↦ 1` (then `==` becomes addition mod 2).

    Equivalently: the class of `a · b` is the identity coset iff `a` and `b`
    sit in the *same* coset (both in `selfInverses`, or both in its
    complement). -/
theorem quotientClass_mul (a b : ℕ)
    (ha : a ∈ admissibleResidues) (hb : b ∈ admissibleResidues) :
    quotientClass ((a * b) % 30)
      = (quotientClass a == quotientClass b) := by
  unfold admissibleResidues at ha hb
  fin_cases ha <;> fin_cases hb <;> decide

/-! ### Additive recasting via `Bool → ZMod 2` -/

/-- Additive encoding of the quotient: `true ↦ 0, false ↦ 1`.

    With this convention, the multiplication on `R / selfInverses` becomes
    addition in `ℤ/2ℤ` (the identity coset goes to `0`, the non-identity
    coset goes to `1`). -/
def quotientCode (r : ℕ) : ZMod 2 :=
  if quotientClass r then 0 else 1

/-- Explicit values of `quotientCode` on the identity coset (`→ 0`). -/
theorem quotientCode_on_selfInverses :
    quotientCode 1 = 0
    ∧ quotientCode 11 = 0
    ∧ quotientCode 19 = 0
    ∧ quotientCode 29 = 0 := by
  refine ⟨?_, ?_, ?_, ?_⟩ <;> decide

/-- Explicit values of `quotientCode` on the non-identity coset (`→ 1`). -/
theorem quotientCode_on_complement :
    quotientCode 7 = 1
    ∧ quotientCode 13 = 1
    ∧ quotientCode 17 = 1
    ∧ quotientCode 23 = 1 := by
  refine ⟨?_, ?_, ?_, ?_⟩ <;> decide

/-- **`quotientCode` is a group homomorphism `R → ℤ/2ℤ`.** For all
    `a, b ∈ R`:

      `quotientCode ((a · b) mod 30) = quotientCode a + quotientCode b`. -/
theorem quotientCode_addMod2 (a b : ℕ)
    (ha : a ∈ admissibleResidues) (hb : b ∈ admissibleResidues) :
    quotientCode ((a * b) % 30) = quotientCode a + quotientCode b := by
  unfold admissibleResidues at ha hb
  fin_cases ha <;> fin_cases hb <;> decide

/-- The kernel of `quotientCode` on `R` is exactly `selfInverses`. -/
theorem quotientCode_kernel :
    admissibleResidues.filter (fun r => quotientCode r = 0) = selfInverses := by
  decide

/-! ### Cardinality of the quotient -/

/-- **The quotient has 2 classes.** The image of `quotientCode` on `R` has
    cardinality `2` (it is exactly `{0, 1} = ZMod 2`). -/
theorem quotient_card_eq_2 :
    (admissibleResidues.image quotientCode).card = 2 := by decide

/-- The image of `quotientCode` on `R` is the full `ZMod 2`. -/
theorem quotientCode_image :
    admissibleResidues.image quotientCode = ({0, 1} : Finset (ZMod 2)) := by
  decide

/-- Cardinality identity: `|R| = 2 · |selfInverses| = 2 · 4 = 8`, matching
    the index of `selfInverses` in `R`. -/
theorem admissible_card_eq_two_times_selfInverses :
    admissibleResidues.card = 2 * selfInverses.card := by decide

/-- **Lagrange identity for the quotient.** `|R| / |selfInverses| = 2`,
    equal to the cardinality of the image of `quotientCode`. -/
theorem quotient_card_lagrange :
    admissibleResidues.card / selfInverses.card
      = (admissibleResidues.image quotientCode).card := by decide

/-! ### Bundled isomorphism `R / selfInverses ≅ ℤ/2ℤ` -/

/-- **Bundled statement: `R / selfInverses ≅ ℤ/2ℤ`.**

    The encoding `quotientCode : R → ℤ/2ℤ` is:
    1. A **group homomorphism** for the multiplication mod 30 (its image
       satisfies the additive rule of `ZMod 2`).
    2. Has **kernel** exactly `selfInverses`.
    3. Has **image** exactly `{0, 1} = ZMod 2` (cardinality `2`).
    4. The **partition** of `R` into its two fibres matches the two cosets
       of `selfInverses`.
    5. **Cardinality match**: `|R / selfInverses| = 2 = |ZMod 2|`. -/
theorem quotient_iso_Z2 :
    -- (1) homomorphism
    (∀ a b, a ∈ admissibleResidues → b ∈ admissibleResidues →
        quotientCode ((a * b) % 30) = quotientCode a + quotientCode b)
    -- (2) kernel = selfInverses
    ∧ admissibleResidues.filter (fun r => quotientCode r = 0) = selfInverses
    -- (3) image = {0, 1}
    ∧ admissibleResidues.image quotientCode = ({0, 1} : Finset (ZMod 2))
    -- (4) coset partition
    ∧ quotientClass_true_set ∪ quotientClass_false_set = admissibleResidues
    -- (5) cardinality match
    ∧ admissibleResidues.card / selfInverses.card = 2
    ∧ Fintype.card (ZMod 2) = 2 :=
  ⟨quotientCode_addMod2,
   quotientCode_kernel,
   quotientCode_image,
   quotientClass_partition,
   selfInverses_quotient_card,
   by decide⟩

/-! ### Headline -/

/-- **Headline (quotient `R / selfInverses`).** The quotient of
    `R = (ℤ/30ℤ)*` by its Klein-four subgroup of involutions
    `selfInverses = {1, 11, 19, 29}` is isomorphic to `ℤ/2ℤ`.

    Explicitly, the Boolean classifier `quotientClass r = decide (r ∈ selfInverses)`
    partitions `R` into the two cosets

      * `selfInverses = {1, 11, 19, 29}` (`true`, identity in the quotient),
      * `R \ selfInverses = {7, 13, 17, 23}` (`false`, non-identity).

    The composition law follows the `XNOR` rule on `Bool`
    (`quotientClass (a·b) = (quotientClass a == quotientClass b)`),
    equivalently addition mod 2 under `true ↦ 0, false ↦ 1`. -/
theorem quotient_R_selfInverses_iso_Z2 :
    -- two cosets, of size 4 each, partitioning R
    quotientClass_true_set = selfInverses
    ∧ quotientClass_false_set = admissibleResidues \ selfInverses
    ∧ quotientClass_true_set.card = 4
    ∧ quotientClass_false_set.card = 4
    ∧ quotientClass_true_set ∪ quotientClass_false_set = admissibleResidues
    -- multiplicative stability (the quotient is a group of order 2)
    ∧ (∀ a b, a ∈ admissibleResidues → b ∈ admissibleResidues →
        quotientClass ((a * b) % 30)
          = (quotientClass a == quotientClass b))
    -- index is 2
    ∧ admissibleResidues.card / selfInverses.card = 2 :=
  ⟨quotientClass_true_set_eq_selfInverses,
   quotientClass_false_set_eq_complement,
   quotientClass_true_card,
   quotientClass_false_card,
   quotientClass_partition,
   quotientClass_mul,
   selfInverses_quotient_card⟩

end PT.Sieve
