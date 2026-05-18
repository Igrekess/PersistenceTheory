/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Sieve.Bimodality
import PT.Sieve.PrimorialCoprime
import PT.Sieve.AdmissibleSquareGroupStructure
import Mathlib.Tactic

/-!
# Klein-four embedding `selfInverses ↪ (ℤ/2ℤ)²` (App N)

This file makes the Klein-four structure of the involutions subgroup
`selfInverses = {1, 11, 19, 29} ⊂ R = (ℤ/30ℤ)*` **explicit and arithmetic**,
via the Chinese Remainder decomposition `30 = 2 · 3 · 5`.

## Strategic comment on the group structure

The textbook decomposition is

$$
(ℤ/30ℤ)^* \;\cong\; (ℤ/2ℤ)^* \times (ℤ/3ℤ)^* \times (ℤ/5ℤ)^*
\;\cong\; \{1\} \times ℤ/2 \times ℤ/4.
$$

So `R` is **not** elementary abelian `(ℤ/2)³` — the `(ℤ/5ℤ)^*` factor is
cyclic of order 4. The (folklore) `(ℤ/2)³` claim that appears elsewhere in
the corpus is a **cardinality coincidence** (`φ(30) = 8 = 2³`), not a
structural statement.

What is genuinely elementary abelian of rank 2 is the **2-torsion subgroup**
`selfInverses = {r ∈ R : r² ≡ 1 (mod 30)} = {1, 11, 19, 29}`. This is the
direct product of the 2-torsion of each cyclic factor:

  * the `ℤ/2` factor (mod 3) contributes its unique non-trivial involution
    (the class of `-1 ≡ 2 mod 3`);
  * the `ℤ/4` factor (mod 5) contributes its unique non-trivial involution
    (the class of `-1 ≡ 4 mod 5`, namely the unique element of order 2 in
    a cyclic group of order 4).

Hence `selfInverses ≅ (ℤ/2ℤ)²`, the Klein-four group.

## Main results

* `crtTriple` — the CRT projection `r ↦ (r % 2, r % 3, r % 5)`.
* `crtTriple_inj_on_admissible` — `crtTriple` is injective on `R` (the
  arithmetic content of CRT for `30 = 2 · 3 · 5`).
* `klein4Code` — the explicit `(ℤ/2ℤ) × (ℤ/2ℤ)` code of an element of
  `selfInverses`, defined by `1 ↦ (0,0), 11 ↦ (1,0), 19 ↦ (0,1), 29 ↦ (1,1)`.
* `klein4Code_mod3_mod5` — semantic interpretation: the code is
  `(sgn(r mod 3), sgn(r mod 5))` where `sgn(1) = 0`, `sgn(-1) = 1`.
* `selfInverses_mul_table_FULL` — complete multiplication table for
  `selfInverses` (16 entries).
* `klein4Code_mul` — `klein4Code (a · b mod 30) = klein4Code a + klein4Code b`
  for `a, b ∈ selfInverses` (group homomorphism).
* `klein4Code_inj` — `klein4Code` is injective on `selfInverses`.
* `klein4Code_surj` — `klein4Code` reaches all `4` points of `ZMod 2 × ZMod 2`.
* `selfInverses_klein_four` — bundled statement: `selfInverses` is closed,
  has cardinality 4, every non-identity element has order 2, and the product
  of two distinct non-identity elements equals the third non-identity.

## Reference

Monograph Appendix N §"Invariants arithmétiques de R", and the M3+ note on
the 2-torsion of `(ℤ/p₁p₂p₃)*` for the PT primorial `30 = p₁p₂p₃`.
-/

namespace PT.Sieve

/-! ### CRT projection -/

/-- The CRT projection `r ↦ (r mod 2, r mod 3, r mod 5)`.
    Encodes an element of `ℤ/30ℤ` by its triple of components under
    `30 = 2 · 3 · 5`. -/
def crtTriple (r : ℕ) : ℕ × ℕ × ℕ := (r % 2, r % 3, r % 5)

/-- The CRT projection of each element of `R` is explicit. -/
theorem crtTriple_values :
    crtTriple 1 = (1, 1, 1)
    ∧ crtTriple 7 = (1, 1, 2)
    ∧ crtTriple 11 = (1, 2, 1)
    ∧ crtTriple 13 = (1, 1, 3)
    ∧ crtTriple 17 = (1, 2, 2)
    ∧ crtTriple 19 = (1, 1, 4)
    ∧ crtTriple 23 = (1, 2, 3)
    ∧ crtTriple 29 = (1, 2, 4) := by
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩ <;> decide

/-- The CRT projection is **injective on `R`** (arithmetic content of the
    Chinese remainder theorem for `30 = 2 · 3 · 5`: an element of `R` is
    determined by its residues mod 2, mod 3, mod 5). -/
theorem crtTriple_inj_on_admissible :
    ∀ a ∈ admissibleResidues, ∀ b ∈ admissibleResidues,
      crtTriple a = crtTriple b → a = b := by
  intro a ha b hb hab
  unfold admissibleResidues at ha hb
  fin_cases ha <;> fin_cases hb <;> simp [crtTriple] at hab <;> rfl

/-! ### CRT projection on `selfInverses` -/

/-- All elements of `selfInverses` have residue 1 modulo 2 (they are odd).
    Their interesting CRT data lives in the `(mod 3) × (mod 5)` factor. -/
theorem selfInverses_mod2 :
    ∀ r ∈ selfInverses, r % 2 = 1 := by
  intro r hr
  unfold selfInverses at hr
  fin_cases hr <;> decide

/-- The mod-3 / mod-5 projection of each element of `selfInverses` is
    one of the four sign combinations:

    | r  | r mod 3 | r mod 5 |
    |----|---------|---------|
    | 1  | 1       | 1       |
    | 11 | 2 (=-1) | 1       |
    | 19 | 1       | 4 (=-1) |
    | 29 | 2 (=-1) | 4 (=-1) | -/
theorem selfInverses_mod3_mod5 :
    (1 % 3, 1 % 5) = (1, 1)
    ∧ (11 % 3, 11 % 5) = (2, 1)
    ∧ (19 % 3, 19 % 5) = (1, 4)
    ∧ (29 % 3, 29 % 5) = (2, 4) := by
  refine ⟨?_, ?_, ?_, ?_⟩ <;> decide

/-! ### Klein-four code -/

/-- The explicit `(ℤ/2ℤ) × (ℤ/2ℤ)` encoding of `selfInverses`:

    `1 ↦ (0, 0)`, `11 ↦ (1, 0)`, `19 ↦ (0, 1)`, `29 ↦ (1, 1)`.

    Anything outside `selfInverses` is mapped to `(0, 0)` (the value is
    irrelevant on the complement; only the restriction to `selfInverses`
    matters for the embedding statements below). -/
def klein4Code (r : ℕ) : ZMod 2 × ZMod 2 :=
  if r = 11 then (1, 0)
  else if r = 19 then (0, 1)
  else if r = 29 then (1, 1)
  else (0, 0)

/-- Explicit values of `klein4Code` on `selfInverses`. -/
theorem klein4Code_values :
    klein4Code 1 = (0, 0)
    ∧ klein4Code 11 = (1, 0)
    ∧ klein4Code 19 = (0, 1)
    ∧ klein4Code 29 = (1, 1) := by
  refine ⟨?_, ?_, ?_, ?_⟩ <;> decide

/-- **Semantic reading of `klein4Code`** : the first coordinate detects
    whether `r ≡ -1 (mod 3)` (i.e. `r % 3 = 2`), the second detects whether
    `r ≡ -1 (mod 5)` (i.e. `r % 5 = 4`). This is the structural meaning of
    the Klein-four embedding, factoring `selfInverses` through

      `2-torsion of (ℤ/3)* × 2-torsion of (ℤ/5)* = ℤ/2 × ℤ/2`. -/
theorem klein4Code_mod3_mod5 :
    klein4Code 1  = ((if 1  % 3 = 2 then 1 else 0), (if 1  % 5 = 4 then 1 else 0))
    ∧ klein4Code 11 = ((if 11 % 3 = 2 then 1 else 0), (if 11 % 5 = 4 then 1 else 0))
    ∧ klein4Code 19 = ((if 19 % 3 = 2 then 1 else 0), (if 19 % 5 = 4 then 1 else 0))
    ∧ klein4Code 29 = ((if 29 % 3 = 2 then 1 else 0), (if 29 % 5 = 4 then 1 else 0)) := by
  refine ⟨?_, ?_, ?_, ?_⟩ <;> decide

/-! ### Full multiplication table of `selfInverses` -/

/-- The complete 4 × 4 multiplication table of `selfInverses` modulo 30
    (16 entries — already provided pointwise in
    `AdmissibleSquareGroupStructure.lean`, bundled here for the
    Klein-four embedding statements). -/
theorem selfInverses_mul_table_FULL :
    (1 * 1 : ℕ) % 30 = 1 ∧ (1 * 11 : ℕ) % 30 = 11
    ∧ (1 * 19 : ℕ) % 30 = 19 ∧ (1 * 29 : ℕ) % 30 = 29
    ∧ (11 * 1 : ℕ) % 30 = 11 ∧ (11 * 11 : ℕ) % 30 = 1
    ∧ (11 * 19 : ℕ) % 30 = 29 ∧ (11 * 29 : ℕ) % 30 = 19
    ∧ (19 * 1 : ℕ) % 30 = 19 ∧ (19 * 11 : ℕ) % 30 = 29
    ∧ (19 * 19 : ℕ) % 30 = 1 ∧ (19 * 29 : ℕ) % 30 = 11
    ∧ (29 * 1 : ℕ) % 30 = 29 ∧ (29 * 11 : ℕ) % 30 = 19
    ∧ (29 * 19 : ℕ) % 30 = 11 ∧ (29 * 29 : ℕ) % 30 = 1 := by
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_,
          ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩ <;> decide

/-! ### Klein-four homomorphism: `klein4Code` preserves multiplication -/

/-- **Klein-four embedding (homomorphism property).** For all
    `a, b ∈ selfInverses`,

      `klein4Code (a · b mod 30) = klein4Code a + klein4Code b`

    in `(ℤ/2ℤ) × (ℤ/2ℤ)`. This is the explicit isomorphism with the
    Klein-four group. -/
theorem klein4Code_mul (a b : ℕ)
    (ha : a ∈ selfInverses) (hb : b ∈ selfInverses) :
    klein4Code ((a * b) % 30) = klein4Code a + klein4Code b := by
  unfold selfInverses at ha hb
  fin_cases ha <;> fin_cases hb <;> decide

/-- **Identity is `(0, 0)`.** -/
theorem klein4Code_one : klein4Code 1 = (0, 0) := by decide

/-- **Every non-identity element of `selfInverses` has order 2** in the code:
    `klein4Code r + klein4Code r = (0, 0)` for `r ∈ {11, 19, 29}` (in fact,
    for any `r` since `ZMod 2 + ZMod 2 = 0` for all elements). -/
theorem klein4Code_order_two (r : ℕ) (hr : r ∈ selfInverses) :
    klein4Code r + klein4Code r = (0, 0) := by
  unfold selfInverses at hr
  fin_cases hr <;> decide

/-- **Klein-four product rule.** The product of two distinct non-identity
    involutions equals the third non-identity:

      `11 · 19 ≡ 29`, `11 · 29 ≡ 19`, `19 · 29 ≡ 11` (mod 30). -/
theorem klein4_non_identity_product :
    (11 * 19 : ℕ) % 30 = 29
    ∧ (11 * 29 : ℕ) % 30 = 19
    ∧ (19 * 29 : ℕ) % 30 = 11 := by
  refine ⟨?_, ?_, ?_⟩ <;> decide

/-! ### Injectivity, surjectivity, bijection -/

/-- `klein4Code` is **injective on `selfInverses`**. -/
theorem klein4Code_inj :
    ∀ a ∈ selfInverses, ∀ b ∈ selfInverses,
      klein4Code a = klein4Code b → a = b := by
  intro a ha b hb hab
  unfold selfInverses at ha hb
  fin_cases ha <;> fin_cases hb <;> simp [klein4Code] at hab <;> rfl

/-- `klein4Code` **hits every point** of `(ℤ/2ℤ) × (ℤ/2ℤ)`. -/
theorem klein4Code_surj :
    ∀ x : ZMod 2 × ZMod 2, ∃ r ∈ selfInverses, klein4Code r = x := by
  intro x
  fin_cases x
  · exact ⟨1, by decide, by decide⟩
  · exact ⟨19, by decide, by decide⟩
  · exact ⟨11, by decide, by decide⟩
  · exact ⟨29, by decide, by decide⟩

/-- **Cardinality match.** Both `selfInverses` and `ZMod 2 × ZMod 2` have
    cardinality `4`. -/
theorem klein4Code_card_match :
    selfInverses.card = 4
    ∧ Fintype.card (ZMod 2 × ZMod 2) = 4 := by
  refine ⟨?_, ?_⟩
  · decide
  · decide

/-! ### Bundled Klein-four statement -/

/-- **`selfInverses` is the Klein-four group `(ℤ/2ℤ)²`.** The bundled
    statement combines:

    1. **Closure** under multiplication mod 30.
    2. **Cardinality** `|selfInverses| = 4`.
    3. **2-torsion** : every element squares to `1 mod 30`.
    4. **Klein-four product rule**: any two distinct non-identity elements
       multiply to the third non-identity.
    5. **Homomorphism** `klein4Code : selfInverses → (ℤ/2ℤ)²` preserves
       multiplication.
    6. **Bijection** : `klein4Code` is injective with full image.

    Together, (5) + (6) constitute the explicit group isomorphism

      `selfInverses ≅ (ℤ/2ℤ) × (ℤ/2ℤ)`. -/
theorem selfInverses_klein_four :
    -- (1) closure
    (∀ a b, a ∈ selfInverses → b ∈ selfInverses → (a * b) % 30 ∈ selfInverses)
    -- (2) cardinality
    ∧ selfInverses.card = 4
    -- (3) 2-torsion
    ∧ (∀ r ∈ selfInverses, (r * r) % 30 = 1)
    -- (4) Klein-four product rule
    ∧ ((11 * 19 : ℕ) % 30 = 29
      ∧ (11 * 29 : ℕ) % 30 = 19
      ∧ (19 * 29 : ℕ) % 30 = 11)
    -- (5) homomorphism
    ∧ (∀ a b, a ∈ selfInverses → b ∈ selfInverses →
        klein4Code ((a * b) % 30) = klein4Code a + klein4Code b)
    -- (6) bijection: injective + surjective onto ZMod 2 × ZMod 2
    ∧ (∀ a ∈ selfInverses, ∀ b ∈ selfInverses,
        klein4Code a = klein4Code b → a = b)
    ∧ (∀ x : ZMod 2 × ZMod 2, ∃ r ∈ selfInverses, klein4Code r = x) :=
  ⟨selfInverses_mul_closed,
   selfInverses_card,
   selfInverses_sq_one,
   klein4_non_identity_product,
   klein4Code_mul,
   klein4Code_inj,
   klein4Code_surj⟩

/-! ### Consequence: quotient `R / selfInverses` has index 2 -/

/-- **Index 2 of `selfInverses` in `R`.** The cardinality quotient
    `|R| / |selfInverses| = 8 / 4 = 2`, consistent with the abstract
    structure

      `R ≅ ℤ/2 × ℤ/4`,
      `selfInverses ≅ 2-torsion of R = ℤ/2 × ℤ/2`,
      `R / selfInverses ≅ ℤ/2` (the order-4 cyclic factor mod its 2-torsion). -/
theorem selfInverses_quotient_card :
    admissibleResidues.card / selfInverses.card = 2 := by decide

/-- **Cardinality identity matching the abstract structure.**
    `8 = 2 · 4 = 2 · |selfInverses|`, with the factor `2` being the index
    `[R : selfInverses]`. -/
theorem admissible_eq_two_times_selfInverses_card :
    admissibleResidues.card = 2 * selfInverses.card := by decide

/-! ### Headline -/

/-- **Headline (Klein-four embedding).** The involutions subgroup
    `selfInverses = {1, 11, 19, 29}` of `R = (ℤ/30ℤ)*` is, via the CRT
    decomposition `30 = 2 · 3 · 5`, **explicitly isomorphic to the
    Klein-four group `(ℤ/2ℤ) × (ℤ/2ℤ)`**.

    The isomorphism is the map `klein4Code` defined by
    `1 ↦ (0,0), 11 ↦ (1,0), 19 ↦ (0,1), 29 ↦ (1,1)`, which factors through
    the sign of `r mod 3` (first coordinate) and the sign of `r mod 5`
    (second coordinate). It preserves multiplication
    (`klein4Code_mul`), is injective (`klein4Code_inj`) and surjective
    (`klein4Code_surj`) onto a 4-element target. -/
theorem klein_four_embedding_headline :
    -- The map is a group isomorphism in the algebraic sense:
    (∀ a b, a ∈ selfInverses → b ∈ selfInverses →
        klein4Code ((a * b) % 30) = klein4Code a + klein4Code b)
    ∧ (∀ a ∈ selfInverses, ∀ b ∈ selfInverses,
        klein4Code a = klein4Code b → a = b)
    ∧ (∀ x : ZMod 2 × ZMod 2, ∃ r ∈ selfInverses, klein4Code r = x)
    -- Cardinality match on both sides:
    ∧ selfInverses.card = 4
    ∧ Fintype.card (ZMod 2 × ZMod 2) = 4
    -- Index in R is 2 (consistent with R ≅ ℤ/2 × ℤ/4):
    ∧ admissibleResidues.card / selfInverses.card = 2 :=
  ⟨klein4Code_mul, klein4Code_inj, klein4Code_surj,
   selfInverses_card,
   klein4Code_card_match.2,
   selfInverses_quotient_card⟩

end PT.Sieve
