/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Sieve.Bimodality
import PT.Sieve.BimodalityT1Projection
import PT.Sieve.BimodalityCardinality
import PT.Sieve.PrimorialCoprime
import PT.Sieve.AdmissibleResiduesArithmetic
import PT.Sieve.CoprimeAdmissibilityProduct
import PT.Sieve.AdmissibleSquareGroupStructure
import PT.Sieve.KleinFourEmbedding
import PT.Sieve.AdmissibleSemigroupQuotient
import PT.Sieve.TotientCascadeIdentities
import Mathlib.Tactic

/-!
# Admissible residues `R` — master characterisation (App N synthesis)

This file is a **pure aggregator** that bundles, in one master statement,
the seven independent characterisations of the admissible residue set

  `R := {1, 7, 11, 13, 17, 19, 23, 29}`

modulo `30`, each previously established in its own module of
`PT/Sieve/`. No new calculation is performed: every clause below cites
a theorem that already lives in the corpus.

## The seven characterisations

1. **Cardinality + coprimality** (`PrimorialCoprime`):
   `|R| = 8` and every `r ∈ R` is coprime to `30`. Equivalently,
   `R = coprimeMod30 = (ℤ/30ℤ)*` set-theoretically.

2. **Totient** (`PrimorialCoprime` + `TotientCascadeIdentities`):
   `|R| = φ(30) = 8`, matching the Euler totient of the third primorial.

3. **Parity** (`AdmissibleResiduesArithmetic`):
   every `r ∈ R` is odd (`r % 2 = 1`).

4. **Sum** (`AdmissibleResiduesArithmetic`):
   `Σ r ∈ R, r = 120 = 4 · 30 = (|R|/2) · primorial 3`.

5. **Product** (`CoprimeAdmissibilityProduct`):
   `(∏ r ∈ R, r) ≡ 1 (mod 30)` — Wilson-like identity for `(ℤ/30ℤ)*`.

6. **Group structure** (`AdmissibleSquareGroupStructure`,
   `KleinFourEmbedding`, `AdmissibleSemigroupQuotient`):
   `R` is the group `(ℤ/30ℤ)*` under multiplication mod 30. It contains
   the Klein-four subgroup of involutions `selfInverses = {1, 11, 19, 29}`
   (index 2), itself containing the squares subgroup
   `squaresMod30 = {1, 19}` (index 4). The quotient
   `R / selfInverses ≅ ℤ/2ℤ`.

7. **Bimodality** (`Bimodality` + `BimodalityCardinality`):
   `R` splits 4 + 4 according to `δ̄(r) ∈ {7, 11}` (equivalently,
   according to the Legendre symbol `(r/5) ∈ {+1, -1}`).

## Master theorem

`admissibleR_master_characterisation` collects all seven simultaneously
in a single big conjunction. It is the canonical reference statement to
cite when an external module needs "everything we know about `R`".

## Reference

Monograph Appendix N (synthesis chapter). Closes the App N family of
modules `PT/Sieve/Admissible*` and `PT/Sieve/Bimodality*`.
-/

namespace PT.Sieve

open Finset

/-! ### Characterisation 1 — cardinality + coprimality -/

/-- **Characterisation 1.** `R` has cardinality 8 and every element is
    coprime to `30`; moreover `R` is **exactly** the set of naturals
    `< 30` coprime to `30`. -/
theorem charac_card_coprime :
    admissibleResidues.card = 8
    ∧ (∀ r ∈ admissibleResidues, Nat.Coprime r 30)
    ∧ coprimeMod30 = admissibleResidues :=
  ⟨admissibleResidues_card,
   admissible_coprime_30,
   coprimeMod30_eq_admissibleResidues⟩

/-- **Uniqueness corollary of Characterisation 1.** If `S ⊆ Finset.range 30`
    has cardinality 8 and every element of `S` is coprime to `30`, then
    `S = admissibleResidues`. (Proof: by `coprimeMod30_eq_admissibleResidues`,
    `S ⊆ coprimeMod30 = admissibleResidues`; both have cardinality 8.) -/
theorem admissibleResidues_unique_card8_coprime
    (S : Finset ℕ) (hS_sub : S ⊆ Finset.range 30) (hS_card : S.card = 8)
    (hS_coprime : ∀ r ∈ S, Nat.Coprime r 30) :
    S = admissibleResidues := by
  have hS_in : S ⊆ coprimeMod30 := by
    intro r hr
    unfold coprimeMod30
    rw [Finset.mem_filter]
    exact ⟨hS_sub hr, hS_coprime r hr⟩
  have hC : coprimeMod30 = admissibleResidues := coprimeMod30_eq_admissibleResidues
  have hS_in' : S ⊆ admissibleResidues := hC ▸ hS_in
  exact Finset.eq_of_subset_of_card_le hS_in'
    (by rw [admissibleResidues_card, hS_card])

/-! ### Characterisation 2 — totient -/

/-- **Characterisation 2.** `|R| = φ(30) = 8`. -/
theorem charac_totient :
    admissibleResidues.card = Nat.totient 30
    ∧ Nat.totient 30 = 8 :=
  ⟨admissibleResidues_card_eq_totient, nat_totient_30⟩

/-! ### Characterisation 3 — parity -/

/-- **Characterisation 3.** Every `r ∈ R` is odd. -/
theorem charac_parity :
    ∀ r ∈ admissibleResidues, r % 2 = 1 :=
  admissible_all_odd

/-! ### Characterisation 4 — sum -/

/-- **Characterisation 4.** `Σ R = 120 = 4 · 30 = (|R|/2) · 30`. -/
theorem charac_sum :
    admissibleResidues.sum id = 120
    ∧ admissibleResidues.sum id = 4 * 30 :=
  ⟨admissibleResidues_sum, admissibleResidues_sum_factored⟩

/-! ### Characterisation 5 — product -/

/-- **Characterisation 5.** `(∏ R) ≡ 1 (mod 30)` (Wilson-like). -/
theorem charac_product :
    admissibleResidues.prod id % 30 = 1
    ∧ admissibleResidues.prod id = 215656441 :=
  ⟨admissibleResidues_prod_mod_30, admissibleResidues_prod_value⟩

/-! ### Characterisation 6 — group structure -/

/-- **Characterisation 6.** `R` is a multiplicative group modulo 30:

    * every `r ∈ R` has a multiplicative inverse `s ∈ R` mod 30;
    * the involutions form a Klein-four subgroup
      `selfInverses = {1, 11, 19, 29}` of index 2;
    * the squares form a subgroup `squaresMod30 = {1, 19}` of index 4;
    * the quotient `R / selfInverses ≅ ℤ/2ℤ`. -/
theorem charac_group_structure :
    -- existence of inverses inside R
    (∀ r ∈ admissibleResidues, ∃ s ∈ admissibleResidues, (r * s) % 30 = 1)
    -- self-inverse subgroup (Klein-four)
    ∧ selfInverses.card = 4
    ∧ selfInverses ⊆ admissibleResidues
    ∧ (∀ a b, a ∈ selfInverses → b ∈ selfInverses → (a * b) % 30 ∈ selfInverses)
    -- squares subgroup (order 2)
    ∧ squaresMod30.card = 2
    ∧ squaresMod30 ⊆ admissibleResidues
    ∧ (∀ a b, a ∈ squaresMod30 → b ∈ squaresMod30 →
        (a * b) % 30 ∈ squaresMod30)
    -- subgroup chain
    ∧ squaresMod30 ⊆ selfInverses
    -- indices
    ∧ admissibleResidues.card / selfInverses.card = 2
    ∧ admissibleResidues.card / squaresMod30.card = 4
    -- quotient R / selfInverses ≅ ℤ/2ℤ (additive code is a homomorphism)
    ∧ (∀ a b, a ∈ admissibleResidues → b ∈ admissibleResidues →
        quotientCode ((a * b) % 30) = quotientCode a + quotientCode b)
    ∧ admissibleResidues.image quotientCode = ({0, 1} : Finset (ZMod 2)) :=
  ⟨admissible_has_inverse_in_R,
   selfInverses_card,
   selfInverses_subset,
   selfInverses_mul_closed,
   squaresMod30_card,
   squaresMod30_subset_admissible,
   squaresMod30_mul_closed,
   squaresMod30_subset_selfInverses,
   selfInverses_quotient_card,
   quotient_card,
   quotientCode_addMod2,
   quotientCode_image⟩

/-! ### Characterisation 7 — bimodality -/

/-- **Characterisation 7.** `R` splits 4 + 4 according to `δ̄(r) ∈ {7, 11}`
    (equivalently, the Legendre symbol `(r/5)`). -/
theorem charac_bimodality :
    -- The two `δ̄`-classes have cardinality 4 each
    lowResidues.card = 4
    ∧ highResidues.card = 4
    -- They are disjoint and union to R
    ∧ Disjoint lowResidues highResidues
    ∧ lowResidues ∪ highResidues = admissibleResidues
    -- And carry the explicit `δ̄ · 6` values 42 / 66
    ∧ deltaBarTimes6 (1  : ZMod 30) = 42
    ∧ deltaBarTimes6 (11 : ZMod 30) = 42
    ∧ deltaBarTimes6 (19 : ZMod 30) = 42
    ∧ deltaBarTimes6 (29 : ZMod 30) = 42
    ∧ deltaBarTimes6 (7  : ZMod 30) = 66
    ∧ deltaBarTimes6 (13 : ZMod 30) = 66
    ∧ deltaBarTimes6 (17 : ZMod 30) = 66
    ∧ deltaBarTimes6 (23 : ZMod 30) = 66 :=
  ⟨low_card,
   high_card,
   low_disjoint_high,
   low_union_high_eq_admissible,
   deltaBar_1, deltaBar_11, deltaBar_19, deltaBar_29,
   deltaBar_7, deltaBar_13, deltaBar_17, deltaBar_23⟩

/-! ### Master theorem -/

/-- **Master characterisation of `R = admissibleResidues`.**

    The admissible residue set `R = {1, 7, 11, 13, 17, 19, 23, 29}` modulo
    `30` is simultaneously characterised by **seven** independent properties:

    1. (cardinality + coprimality) `|R| = 8`; every `r ∈ R` is coprime to
       `30`; in fact `R` is **exactly** the set of `r < 30` coprime to `30`.
    2. (totient) `|R| = φ(30) = 8`.
    3. (parity) every `r ∈ R` is odd.
    4. (sum) `Σ R = 120 = 4 · 30`.
    5. (product) `(∏ R) ≡ 1 (mod 30)`.
    6. (group) `R` is a multiplicative group mod 30, with Klein-four
       subgroup of involutions `{1, 11, 19, 29}` of index 2, squares
       subgroup `{1, 19}` of index 4, and quotient
       `R / selfInverses ≅ ℤ/2ℤ`.
    7. (bimodality) `R` partitions 4 + 4 into `lowResidues = {1, 11, 19, 29}`
       (with `δ̄·6 = 42`) and `highResidues = {7, 13, 17, 23}`
       (with `δ̄·6 = 66`).

    Each clause is a literal restatement of a theorem proven elsewhere in
    `PT/Sieve/`. The bundling is the **canonical reference** for downstream
    modules that need any combination of these facts. -/
theorem admissibleR_master_characterisation :
    -- (1) cardinality + coprimality
    (admissibleResidues.card = 8
      ∧ (∀ r ∈ admissibleResidues, Nat.Coprime r 30)
      ∧ coprimeMod30 = admissibleResidues)
    -- (2) totient
    ∧ (admissibleResidues.card = Nat.totient 30
      ∧ Nat.totient 30 = 8)
    -- (3) parity
    ∧ (∀ r ∈ admissibleResidues, r % 2 = 1)
    -- (4) sum
    ∧ (admissibleResidues.sum id = 120
      ∧ admissibleResidues.sum id = 4 * 30)
    -- (5) product
    ∧ (admissibleResidues.prod id % 30 = 1
      ∧ admissibleResidues.prod id = 215656441)
    -- (6) group structure
    ∧ ((∀ r ∈ admissibleResidues, ∃ s ∈ admissibleResidues, (r * s) % 30 = 1)
      ∧ selfInverses.card = 4
      ∧ selfInverses ⊆ admissibleResidues
      ∧ (∀ a b, a ∈ selfInverses → b ∈ selfInverses →
          (a * b) % 30 ∈ selfInverses)
      ∧ squaresMod30.card = 2
      ∧ squaresMod30 ⊆ admissibleResidues
      ∧ (∀ a b, a ∈ squaresMod30 → b ∈ squaresMod30 →
          (a * b) % 30 ∈ squaresMod30)
      ∧ squaresMod30 ⊆ selfInverses
      ∧ admissibleResidues.card / selfInverses.card = 2
      ∧ admissibleResidues.card / squaresMod30.card = 4
      ∧ (∀ a b, a ∈ admissibleResidues → b ∈ admissibleResidues →
          quotientCode ((a * b) % 30) = quotientCode a + quotientCode b)
      ∧ admissibleResidues.image quotientCode = ({0, 1} : Finset (ZMod 2)))
    -- (7) bimodality
    ∧ (lowResidues.card = 4
      ∧ highResidues.card = 4
      ∧ Disjoint lowResidues highResidues
      ∧ lowResidues ∪ highResidues = admissibleResidues
      ∧ deltaBarTimes6 (1  : ZMod 30) = 42
      ∧ deltaBarTimes6 (11 : ZMod 30) = 42
      ∧ deltaBarTimes6 (19 : ZMod 30) = 42
      ∧ deltaBarTimes6 (29 : ZMod 30) = 42
      ∧ deltaBarTimes6 (7  : ZMod 30) = 66
      ∧ deltaBarTimes6 (13 : ZMod 30) = 66
      ∧ deltaBarTimes6 (17 : ZMod 30) = 66
      ∧ deltaBarTimes6 (23 : ZMod 30) = 66) :=
  ⟨charac_card_coprime,
   charac_totient,
   charac_parity,
   charac_sum,
   charac_product,
   charac_group_structure,
   charac_bimodality⟩

/-! ### Equivalence-style consequence: any of (1) alone pins `R` down -/

/-- **Sharpened uniqueness.** Among subsets of `{0, 1, …, 29}` of
    cardinality 8 whose elements are all coprime to 30, there is **exactly
    one**: `R = admissibleResidues`. So characterisation 1 alone (without
    any of 2–7) already determines `R`. -/
theorem admissibleR_unique_from_charac1
    (S : Finset ℕ) (hS_sub : S ⊆ Finset.range 30) (hS_card : S.card = 8)
    (hS_coprime : ∀ r ∈ S, Nat.Coprime r 30) :
    S = admissibleResidues :=
  admissibleResidues_unique_card8_coprime S hS_sub hS_card hS_coprime

end PT.Sieve
