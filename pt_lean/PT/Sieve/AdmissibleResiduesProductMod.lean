/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Sieve.Bimodality
import PT.Sieve.CoprimeAdmissibilityProduct
import Mathlib.Tactic

/-!
# Admissible residues — Product modulo several moduli (App N extension)

This file extends `CoprimeAdmissibilityProduct.lean` by tabulating the residue
of

  `Π R := ∏ r ∈ R, r = 215656441`,
  `R = {1, 7, 11, 13, 17, 19, 23, 29} = (ℤ/30ℤ)*`

modulo every divisor `q | 30` and a few additional natural moduli
(`q ∈ {2, 3, 5, 6, 10, 15, 30, 60, 100}`). It records the striking fact that

  `Π R ≡ 1 mod q` for every divisor `q | 30` (and even mod `60`).

This is the **Wilson-like phenomenon** on `(ℤ/30ℤ)*`: since this group is
isomorphic to `(ℤ/2ℤ)³`, every element is its own inverse modulo the
involutions, and the product of all elements equals the product of all
involutions, namely `1 · 11 · 19 · 29 = 6061 ≡ 1 mod 30`. Restricted to any
divisor `q | 30`, the quotient inherits the property.

By contrast, mod `100` (which does not divide a power of `30`) the product
ends in the last two digits `41` — purely arithmetic data, no group-theoretic
collapse.

## Main results

* `prod_value_eq_factorisation`: `Π R = 7·11·13·17·19·23·29`.
* `prod_mod_q` for `q ∈ {2, 3, 5, 6, 10, 15, 30, 60, 100}`.
* `prod_mod_divisor_of_30`: any divisor `q | 30` gives residue `1`.
* `prod_mod_60_eq_one`: extends Wilson-like collapse one step beyond `30`.
* `prod_mod_100_eq_41`: the last two digits of `Π R` are `41`.
* `prod_mod_table`: the headline table of residues
  `(2 ↦ 1, 3 ↦ 1, 5 ↦ 1, 6 ↦ 1, 10 ↦ 1, 15 ↦ 1, 30 ↦ 1, 60 ↦ 1, 100 ↦ 41)`.

## Reference

Monograph Appendix N §"Invariants arithmétiques de R" (multiplicative part).
-/

namespace PT.Sieve

open Finset

/-! ### Factorisation of the product -/

/-- The product `Π R` factors as `7·11·13·17·19·23·29` (the seven odd primes
    `≤ 29` distinct from `2, 3, 5`). The factor `1 ∈ R` is multiplicatively
    trivial. -/
theorem prod_value_eq_factorisation :
    admissibleResidues.prod id = 7 * 11 * 13 * 17 * 19 * 23 * 29 := by
  decide

/-- Numerical sanity check: `7·11·13·17·19·23·29 = 215656441`. -/
theorem factorisation_value :
    (7 * 11 * 13 * 17 * 19 * 23 * 29 : ℕ) = 215656441 := by decide

/-! ### Residue of `Π R` modulo small moduli

For each `q`, we use `admissibleResidues_prod_value` (which evaluates
`Π R = 215656441`) and then `decide` to verify the residue. -/

/-- `Π R ≡ 1 mod 2`: the product is odd (no factor of `2`). -/
theorem prod_mod_2 : admissibleResidues.prod id % 2 = 1 := by decide

/-- `Π R ≡ 1 mod 3`: the product is coprime to `3` and ≡ `1`. -/
theorem prod_mod_3 : admissibleResidues.prod id % 3 = 1 := by decide

/-- `Π R ≡ 1 mod 5`: the product is coprime to `5` and ≡ `1`. -/
theorem prod_mod_5 : admissibleResidues.prod id % 5 = 1 := by decide

/-- `Π R ≡ 1 mod 6`: combining mod-2 and mod-3 via CRT. -/
theorem prod_mod_6 : admissibleResidues.prod id % 6 = 1 := by decide

/-- `Π R ≡ 1 mod 10`: the last digit is `1` (`215656441`). -/
theorem prod_mod_10 : admissibleResidues.prod id % 10 = 1 := by decide

/-- `Π R ≡ 1 mod 15`: combining mod-3 and mod-5 via CRT. -/
theorem prod_mod_15 : admissibleResidues.prod id % 15 = 1 := by decide

/-- `Π R ≡ 1 mod 30` (already proven in `CoprimeAdmissibilityProduct.lean`,
    re-exported here for the table). -/
theorem prod_mod_30 : admissibleResidues.prod id % 30 = 1 :=
  admissibleResidues_prod_mod_30

/-- `Π R ≡ 1 mod 60`: the Wilson-like collapse extends one step beyond
    `30` because `Π R = 215656441` is `≡ 1 mod 8` (not only mod `2`),
    and combining `mod 8` with `mod 15` via CRT gives `mod 60`. Note that
    in general `Π R mod 8 = 1` (the seven odd primes `7,11,13,17,19,23,29`
    are themselves `≡ ±1, ±3 mod 8`, and their product collapses). -/
theorem prod_mod_60 : admissibleResidues.prod id % 60 = 1 := by decide

/-- `Π R ≡ 41 mod 100`: the last two digits of `215656441`. Outside the
    Wilson regime: `100 = 4 · 25` is not a divisor of any power of `30`, so
    no group-theoretic collapse applies and we read off the residue
    arithmetically. -/
theorem prod_mod_100 : admissibleResidues.prod id % 100 = 41 := by decide

/-! ### Bonus: residue mod 8 (the source of the `mod 60` collapse) -/

/-- `Π R ≡ 1 mod 8`. Together with `Π R ≡ 1 mod 15`, this is the CRT input
    for `prod_mod_60`. Note that `8 ∤ 30`, so this is *not* a divisor of
    `30`: the collapse mod `8` is a separate Wilson-like cancellation in
    `(ℤ/8ℤ)*` (which is also a Klein four-group `(ℤ/2ℤ)²`). -/
theorem prod_mod_8 : admissibleResidues.prod id % 8 = 1 := by decide

/-! ### Wilson collapse for any divisor of `30` -/

/-- **Wilson-like collapse on non-trivial divisors of 30.** For every
    divisor `q ∈ {2, 3, 5, 6, 10, 15, 30}` of `30` (excluding the trivial
    divisor `q = 1`, for which `n % 1 = 0`), we have `Π R ≡ 1 mod q`. -/
theorem prod_mod_divisor_of_30 (q : ℕ)
    (hq : q ∈ ({2, 3, 5, 6, 10, 15, 30} : Finset ℕ)) :
    admissibleResidues.prod id % q = 1 := by
  rw [admissibleResidues_prod_value]
  fin_cases hq <;> decide

/-! ### Headline: the table of residues -/

/-- **Headline (product residue table).** The product
    `Π R = ∏ r ∈ R, r = 215656441` satisfies the following residue table:

    | `q`   |  2 |  3 |  5 |  6 | 10 | 15 | 30 | 60 | 100 |
    |-------|----|----|----|----|----|----|----|----|-----|
    | `Π R mod q` | 1 | 1 | 1 | 1 | 1 | 1 | 1 |  1 |  41 |

    The first eight entries embody the **Wilson collapse** on the units of
    `(ℤ/30ℤ)* ≅ (ℤ/2ℤ)³` (and its lift to `(ℤ/60ℤ)* ≅ (ℤ/2ℤ)⁴` via the
    extra factor of `8`). The last entry is purely arithmetic. -/
theorem prod_mod_table :
    admissibleResidues.prod id %   2 =  1
    ∧ admissibleResidues.prod id %   3 =  1
    ∧ admissibleResidues.prod id %   5 =  1
    ∧ admissibleResidues.prod id %   6 =  1
    ∧ admissibleResidues.prod id %  10 =  1
    ∧ admissibleResidues.prod id %  15 =  1
    ∧ admissibleResidues.prod id %  30 =  1
    ∧ admissibleResidues.prod id %  60 =  1
    ∧ admissibleResidues.prod id % 100 = 41 :=
  ⟨prod_mod_2, prod_mod_3, prod_mod_5, prod_mod_6, prod_mod_10,
   prod_mod_15, prod_mod_30, prod_mod_60, prod_mod_100⟩

/-! ### Pattern: divisors of 30 collapse uniformly -/

/-- **Uniform Wilson collapse.** For *every* `q` in the canonical list
    `{2, 3, 5, 6, 10, 15, 30, 60}` (the non-trivial divisors of `30`
    together with `60 = lcm(8, 15)`), the product `Π R` reduces to `1`
    modulo `q`. -/
theorem prod_mod_pattern :
    ∀ q ∈ ({2, 3, 5, 6, 10, 15, 30, 60} : Finset ℕ),
      admissibleResidues.prod id % q = 1 := by
  intro q hq
  rw [admissibleResidues_prod_value]
  fin_cases hq <;> decide

end PT.Sieve
