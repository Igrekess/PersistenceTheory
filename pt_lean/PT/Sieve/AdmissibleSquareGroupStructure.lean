/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Sieve.Bimodality
import PT.Sieve.PrimorialCoprime
import PT.Sieve.CoprimeAdmissibilityProduct
import Mathlib.Tactic

/-!
# Admissible residues — Subgroup structure of squares modulo 30 (App N)

This file formalises the **subgroup structure** of the squares of
`R = (ℤ/30ℤ)* = {1, 7, 11, 13, 17, 19, 23, 29}` modulo `30`.

Companion to `CoprimeAdmissibilityProduct.lean`, which proves that the
image of the squaring map `r ↦ r² mod 30` on `R` is exactly `{1, 19}`
(`squares_image`).

We show here that:

## Main results

* `squaresMod30 = {1, 19}`: the squares of `R` mod `30`.
* `squaresMod30_subset_admissible`: `{1, 19} ⊆ R`.
* `squaresMod30_card`: `|{1, 19}| = 2`.
* `squaresMod30_mul_closed`: `{1, 19}` is closed under multiplication mod
  `30`. Concretely: `19 · 19 = 361 ≡ 1 mod 30`, the only non-trivial case.
* `squaresMod30_subgroup`: full multiplication table, witnessing the
  subgroup property.
* `quotient_card`: `|R| / |{1, 19}| = 8 / 2 = 4`. The index is `4`.
* `cosets_explicit`: the 4 cosets of `{1, 19}` in `R` are
  `{1, 19}`, `{7, 13}`, `{11, 29}`, `{17, 23}`.
* `selfInverses_mul_closed`: the involutions `{1, 11, 19, 29}` are also a
  subgroup (Klein-four), of index `2`.

The point of these statements is to make explicit the lattice of
characteristic subgroups of `R ≅ (ℤ/2ℤ)³` (an elementary abelian
2-group with three involutions and seven proper non-trivial subgroups):
* the squares subgroup `{1, 19}` (order 2, index 4);
* the involutions subgroup `{1, 11, 19, 29}` (Klein-four, order 4, index 2).

## Reference

Monograph Appendix N §"Invariants arithmétiques de R" (subgroup lattice).
-/

namespace PT.Sieve

open Finset

/-! ### The squares set `{1, 19}` -/

/-- The set of squares modulo `30` of admissible residues:
    `{r² mod 30 : r ∈ R} = {1, 19}`.

    See `CoprimeAdmissibilityProduct.squares_image` for the proof that this
    enumeration agrees with the image of `r ↦ r² mod 30` on `R`. -/
def squaresMod30 : Finset ℕ := {1, 19}

/-- The squares set has cardinality `2`. -/
theorem squaresMod30_card : squaresMod30.card = 2 := by decide

/-- The squares set is a subset of `R`: `{1, 19} ⊆ {1, 7, 11, 13, 17, 19, 23, 29}`. -/
theorem squaresMod30_subset_admissible :
    squaresMod30 ⊆ admissibleResidues := by decide

/-- The squares set is a strict (proper) subset of `R`. -/
theorem squaresMod30_ssubset_admissible :
    squaresMod30 ⊂ admissibleResidues := by decide

/-- Consistency with `squares_image`: `squaresMod30` is exactly the image of
    `r ↦ r² mod 30` on `R`. -/
theorem squaresMod30_eq_image :
    squaresMod30 = admissibleResidues.image (fun r => (r * r) % 30) := by
  rw [squares_image]; rfl

/-! ### Multiplicative closure of `squaresMod30` -/

/-- `1 · 1 ≡ 1 (mod 30)`. -/
theorem squaresMod30_mul_1_1 : (1 * 1 : ℕ) % 30 = 1 := by decide

/-- `1 · 19 ≡ 19 (mod 30)`. -/
theorem squaresMod30_mul_1_19 : (1 * 19 : ℕ) % 30 = 19 := by decide

/-- `19 · 1 ≡ 19 (mod 30)`. -/
theorem squaresMod30_mul_19_1 : (19 * 1 : ℕ) % 30 = 19 := by decide

/-- `19 · 19 = 361 = 12·30 + 1 ≡ 1 (mod 30)`: the only non-trivial closure. -/
theorem squaresMod30_mul_19_19 : (19 * 19 : ℕ) % 30 = 1 := by decide

/-- **Multiplicative closure of `squaresMod30`.**
    For all `a, b ∈ {1, 19}`, the product `a · b mod 30` again lies in `{1, 19}`. -/
theorem squaresMod30_mul_closed (a b : ℕ)
    (ha : a ∈ squaresMod30) (hb : b ∈ squaresMod30) :
    (a * b) % 30 ∈ squaresMod30 := by
  unfold squaresMod30 at ha hb ⊢
  fin_cases ha <;> fin_cases hb <;> decide

/-- The full multiplication table for `squaresMod30 = {1, 19}` mod `30`,
    witnessing the subgroup property:

    | ·  | 1  | 19 |
    |----|----|----|
    | 1  | 1  | 19 |
    | 19 | 19 | 1  | -/
theorem squaresMod30_subgroup :
    (1 * 1 : ℕ) % 30 = 1
    ∧ (1 * 19 : ℕ) % 30 = 19
    ∧ (19 * 1 : ℕ) % 30 = 19
    ∧ (19 * 19 : ℕ) % 30 = 1 :=
  ⟨squaresMod30_mul_1_1, squaresMod30_mul_1_19,
   squaresMod30_mul_19_1, squaresMod30_mul_19_19⟩

/-- The squares set contains the identity `1`. -/
theorem one_mem_squaresMod30 : (1 : ℕ) ∈ squaresMod30 := by decide

/-- Every element of `squaresMod30` is its own inverse mod `30`. -/
theorem squaresMod30_self_inverse (a : ℕ) (ha : a ∈ squaresMod30) :
    (a * a) % 30 = 1 := by
  unfold squaresMod30 at ha
  fin_cases ha <;> decide

/-- `squaresMod30 ⊆ selfInverses`: every square is an involution. -/
theorem squaresMod30_subset_selfInverses :
    squaresMod30 ⊆ selfInverses := by decide

/-! ### Index 4: quotient cardinality `|R| / |{1,19}| = 4` -/

/-- The cardinality of `R` is `8`. -/
theorem admissibleResidues_card_eq_8 : admissibleResidues.card = 8 := by decide

/-- **Quotient cardinality.** `|R| / |squaresMod30| = 8 / 2 = 4`.
    This is the index of `squaresMod30` in `R`. -/
theorem quotient_card :
    admissibleResidues.card / squaresMod30.card = 4 := by decide

/-! ### The four cosets of `squaresMod30` in `R` -/

/-- Coset `1 · {1, 19} = {1, 19}` (the trivial coset). -/
def coset₁ : Finset ℕ := {1, 19}

/-- Coset `7 · {1, 19} mod 30 = {7, 13}` (since `7·19 = 133 ≡ 13`). -/
def coset₇ : Finset ℕ := {7, 13}

/-- Coset `11 · {1, 19} mod 30 = {11, 29}` (since `11·19 = 209 ≡ 29`). -/
def coset₁₁ : Finset ℕ := {11, 29}

/-- Coset `17 · {1, 19} mod 30 = {17, 23}` (since `17·19 = 323 ≡ 23`). -/
def coset₁₇ : Finset ℕ := {17, 23}

/-- Coset `7 · {1, 19} mod 30 = {7, 13}` — explicit verification. -/
theorem coset₇_eq :
    ({1, 19} : Finset ℕ).image (fun s => (7 * s) % 30) = coset₇ := by decide

/-- Coset `11 · {1, 19} mod 30 = {11, 29}` — explicit verification. -/
theorem coset₁₁_eq :
    ({1, 19} : Finset ℕ).image (fun s => (11 * s) % 30) = coset₁₁ := by decide

/-- Coset `17 · {1, 19} mod 30 = {17, 23}` — explicit verification. -/
theorem coset₁₇_eq :
    ({1, 19} : Finset ℕ).image (fun s => (17 * s) % 30) = coset₁₇ := by decide

/-- Coset `1 · {1, 19} mod 30 = {1, 19}` — explicit verification. -/
theorem coset₁_eq :
    ({1, 19} : Finset ℕ).image (fun s => (1 * s) % 30) = coset₁ := by decide

/-- Each coset has cardinality `2`. -/
theorem cosets_card :
    coset₁.card = 2 ∧ coset₇.card = 2
    ∧ coset₁₁.card = 2 ∧ coset₁₇.card = 2 := by
  refine ⟨?_, ?_, ?_, ?_⟩ <;> decide

/-- The four cosets are pairwise disjoint. -/
theorem cosets_disjoint :
    Disjoint coset₁ coset₇
    ∧ Disjoint coset₁ coset₁₁
    ∧ Disjoint coset₁ coset₁₇
    ∧ Disjoint coset₇ coset₁₁
    ∧ Disjoint coset₇ coset₁₇
    ∧ Disjoint coset₁₁ coset₁₇ := by
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_⟩ <;> decide

/-- **Coset decomposition.** The four cosets of `squaresMod30` in `R`
    partition `R`:
    `R = {1,19} ⊔ {7,13} ⊔ {11,29} ⊔ {17,23}`. -/
theorem cosets_partition :
    coset₁ ∪ coset₇ ∪ coset₁₁ ∪ coset₁₇ = admissibleResidues := by decide

/-- **Index 4 — full statement.** The four cosets are explicit, pairwise
    disjoint, of cardinality `2`, and partition `R`. Hence
    `[R : squaresMod30] = 4`. -/
theorem index_4 :
    coset₁.card = 2 ∧ coset₇.card = 2 ∧ coset₁₁.card = 2 ∧ coset₁₇.card = 2
    ∧ Disjoint coset₁ coset₇
    ∧ Disjoint coset₁ coset₁₁
    ∧ Disjoint coset₁ coset₁₇
    ∧ Disjoint coset₇ coset₁₁
    ∧ Disjoint coset₇ coset₁₇
    ∧ Disjoint coset₁₁ coset₁₇
    ∧ coset₁ ∪ coset₇ ∪ coset₁₁ ∪ coset₁₇ = admissibleResidues := by
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩ <;> decide

/-! ### Link with the structure `(ℤ/2ℤ)³` -/

/-- `R` has order `2³ = 8 = φ(30)`, consistent with the isomorphism
    `R ≅ (ℤ/2ℤ)³` (an elementary abelian 2-group of rank 3).

    The subgroup `squaresMod30 = {1, 19}` has order `2`, so it has index
    `4 = 2²`, matching the corank `2` of a 1-dimensional 𝔽₂-subspace in
    a 3-dimensional 𝔽₂-vector space. -/
theorem klein_eight_structure :
    admissibleResidues.card = 2 ^ 3
    ∧ squaresMod30.card = 2
    ∧ admissibleResidues.card / squaresMod30.card = 2 ^ 2 := by
  refine ⟨?_, ?_, ?_⟩ <;> decide

/-! ### Alternative subgroup: the involutions `selfInverses = {1, 11, 19, 29}` -/

/-- `1 · 1 ≡ 1`. -/
theorem selfInverses_mul_1_1 : (1 * 1 : ℕ) % 30 = 1 := by decide
/-- `1 · 11 ≡ 11`. -/
theorem selfInverses_mul_1_11 : (1 * 11 : ℕ) % 30 = 11 := by decide
/-- `1 · 19 ≡ 19`. -/
theorem selfInverses_mul_1_19 : (1 * 19 : ℕ) % 30 = 19 := by decide
/-- `1 · 29 ≡ 29`. -/
theorem selfInverses_mul_1_29 : (1 * 29 : ℕ) % 30 = 29 := by decide
/-- `11 · 11 = 121 ≡ 1`. -/
theorem selfInverses_mul_11_11 : (11 * 11 : ℕ) % 30 = 1 := by decide
/-- `11 · 19 = 209 ≡ 29`. -/
theorem selfInverses_mul_11_19 : (11 * 19 : ℕ) % 30 = 29 := by decide
/-- `11 · 29 = 319 ≡ 19`. -/
theorem selfInverses_mul_11_29 : (11 * 29 : ℕ) % 30 = 19 := by decide
/-- `19 · 19 = 361 ≡ 1`. -/
theorem selfInverses_mul_19_19 : (19 * 19 : ℕ) % 30 = 1 := by decide
/-- `19 · 29 = 551 ≡ 11`. -/
theorem selfInverses_mul_19_29 : (19 * 29 : ℕ) % 30 = 11 := by decide
/-- `29 · 29 = 841 ≡ 1`. -/
theorem selfInverses_mul_29_29 : (29 * 29 : ℕ) % 30 = 1 := by decide

/-- **Multiplicative closure of `selfInverses`.**
    For all `a, b ∈ {1, 11, 19, 29}`, the product `a · b mod 30` again lies
    in `selfInverses`. This is the Klein-four subgroup of `R`. -/
theorem selfInverses_mul_closed (a b : ℕ)
    (ha : a ∈ selfInverses) (hb : b ∈ selfInverses) :
    (a * b) % 30 ∈ selfInverses := by
  unfold selfInverses at ha hb ⊢
  fin_cases ha <;> fin_cases hb <;> decide

/-- **Index 2.** `|R| / |selfInverses| = 8 / 4 = 2`. -/
theorem selfInverses_index_2 :
    admissibleResidues.card / selfInverses.card = 2 := by decide

/-- The squares subgroup is contained in the involutions subgroup
    (a chain `{1} < {1,19} < {1,11,19,29} < R` in the subgroup lattice). -/
theorem squaresMod30_subgroup_chain :
    squaresMod30 ⊆ selfInverses
    ∧ selfInverses ⊆ admissibleResidues := by
  exact ⟨squaresMod30_subset_selfInverses, selfInverses_subset⟩

/-! ### Headline -/

/-- **Headline (subgroup structure).** The squares set
    `squaresMod30 = {1, 19}` of `R = (ℤ/30ℤ)*`:

    * has cardinality `2` and is a subset of `R`;
    * is closed under multiplication mod `30` (subgroup of order `2`);
    * has index `4` in `R`, witnessed by 4 explicit cosets of size 2;
    * fits the chain `{1, 19} ⊂ {1, 11, 19, 29} ⊂ R`, matching the
      filtration of `R ≅ (ℤ/2ℤ)³` by 𝔽₂-subspaces. -/
theorem squaresMod30_structure_summary :
    squaresMod30.card = 2
    ∧ squaresMod30 ⊆ admissibleResidues
    ∧ (∀ a b, a ∈ squaresMod30 → b ∈ squaresMod30 →
        (a * b) % 30 ∈ squaresMod30)
    ∧ admissibleResidues.card / squaresMod30.card = 4
    ∧ coset₁ ∪ coset₇ ∪ coset₁₁ ∪ coset₁₇ = admissibleResidues
    ∧ squaresMod30 ⊆ selfInverses :=
  ⟨squaresMod30_card,
   squaresMod30_subset_admissible,
   squaresMod30_mul_closed,
   quotient_card,
   cosets_partition,
   squaresMod30_subset_selfInverses⟩

end PT.Sieve
