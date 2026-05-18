/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.FixedPoint.T7MuStar
import PT.FixedPoint.T7GlobalUniqueness
import PT.FixedPoint.DimensionProtection
import PT.FixedPoint.CrystallisationBinary
import PT.Sieve.N3aFactorisation
import PT.Sieve.TotientCascadeIdentities
import Mathlib.Tactic

/-!
# FixedPoint Master Theorem — `μ⋆ = 15` and its full arithmetic ecosystem

This file is a **pure synthesis** master theorem for the canonical PT
distinguished integer `μ⋆ = 15`. It does *not* introduce any new
computation; instead it aggregates the headline results that have been
proved separately in:

* `PT/FixedPoint/T7MuStar.lean` — `muStar = 15`, `T7_muStar_isFixed`,
  `persActiveAt15 = {3,5,7}`, `T7_unique_fixedpoint_small`.
* `PT/FixedPoint/T7GlobalUniqueness.lean` — `T7_muStar_unique_global`,
  `T7_fixedpoint_iff`, `T7_all_fixed_points`.
* `PT/FixedPoint/DimensionProtection.lean` — `persActive_card_eq_three`,
  `persActive_card_le_three`, `dimension_protection`,
  `Fpers_eq_sum_persActive`.
* `PT/FixedPoint/CrystallisationBinary.lean` — cosmogonic identities
  `17 - 2 = 15`, `2 + 3 + 5 + 7 = 17`, `3 + 5 + 7 = 15`,
  `crystallisation_binary`.
* `PT/Sieve/N3aFactorisation.lean` — `nat_factorise_muStar : 15 = 3 · 5`,
  `muStar_primeFactors`, `active_primes_sum_eq_muStar`.
* `PT/Sieve/TotientCascadeIdentities.lean` — `Nat.totient 15 = 8`,
  `totient_mu_star_eq_admissible_card`.

The result is a **single master bundle** characterising `μ⋆ = 15` along
its five canonical axes — **fixed point**, **factorisation**,
**dimension**, **totient resonance**, and **cosmogonic crystallisation** —
plus a *uniqueness* statement and a *characterisation* iff.

This file is the natural top-level reference for any downstream module
that needs to invoke "the master `μ⋆` theorem" without having to recall
which sub-file proves which sub-statement.

## Reference

Monograph chapters `ch07_convergence.tex` (T7 fixed-point chain) and
`ch08_fixed_point.tex` (crystallisation, dimension protection), Appendix N
(totient cascade); audit rows #3 (factorisation), #34 (global uniqueness),
#35 (dimension protection).
-/

namespace PT.FixedPoint

open Nat PT.Sieve

/-! ### Section 1 — Master fixed-point statements

The persistence-active map `F_pers` has `μ⋆ = 15` as its unique fixed
point on the non-degenerate range `μ ≥ 8`, and *exactly* two fixed
points on `μ ≥ 2`: the degenerate one-prime cascade `μ = 3` and the
PT canonical `μ⋆ = 15`. -/

/-- **Master fixed point — existence.** `μ⋆ = 15` is a fixed point of
`F_pers`: `F_pers 15 = 15`. -/
theorem master_fixedpoint_existence : Fpers muStar = muStar :=
  T7_muStar_isFixed

/-- **Master fixed point — global uniqueness on `μ ≥ 8`.** Any integer
`μ ≥ 8` satisfying `F_pers μ = μ` must equal `15`. -/
theorem master_fixedpoint_unique (μ : ℕ) (h : 8 ≤ μ) (hfix : Fpers μ = μ) :
    μ = 15 :=
  T7_muStar_unique_global μ h hfix

/-- **Master fixed point — iff characterisation on `μ ≥ 8`.** -/
theorem master_fixedpoint_iff (μ : ℕ) (h : 8 ≤ μ) :
    Fpers μ = μ ↔ μ = 15 :=
  T7_fixedpoint_iff μ h

/-- **Master fixed point — full classification on `μ ≥ 2`.** On the
whole non-trivial range, the fixed points of `F_pers` are *exactly*
`{3, 15}`: the degenerate one-prime cascade and `μ⋆`. -/
theorem master_fixedpoint_full_classification (μ : ℕ) (h : 2 ≤ μ)
    (hfix : Fpers μ = μ) : μ = 3 ∨ μ = 15 :=
  T7_all_fixed_points μ h hfix

/-! ### Section 2 — Master factorisation statements

`μ⋆ = 15` admits two canonical decompositions: a *multiplicative* one
(`15 = 3 · 5`, the product of the two smallest active primes) and an
*additive* one (`15 = 3 + 5 + 7`, the sum of the three active primes).
These two decompositions are *distinct* — `3 · 5 ≠ 3 + 5 + 7` would be
mis-stated; the point is rather that the **same** integer is reached by
two different arithmetic procedures over the same prime set, which is a
PT signature. -/

/-- **Master factorisation — multiplicative.** `μ⋆ = 15 = 3 · 5`
(product of the two smallest active primes). -/
theorem master_factorisation_multiplicative : muStar = 3 * 5 :=
  nat_factorise_muStar

/-- **Master factorisation — additive.** `μ⋆ = 15 = 3 + 5 + 7`
(sum of all three active primes). -/
theorem master_factorisation_additive : muStar = 3 + 5 + 7 :=
  muStar_eq_three_plus_five_plus_seven

/-- **Master factorisation — prime factors.** The set of prime factors of
`μ⋆ = 15` is exactly `{3, 5}` (the two smallest active primes). The third
active prime `7` enters only additively. -/
theorem master_factorisation_primeFactors :
    muStar.primeFactors = ({3, 5} : Finset ℕ) := by
  rw [muStar_eq]; exact muStar_primeFactors

/-- **Master factorisation — composite.** `μ⋆ = 15` is *not* prime. -/
theorem master_factorisation_not_prime : ¬ Nat.Prime muStar := by
  rw [muStar_eq]; exact muStar_not_prime

/-- **Master factorisation — combined identity.** The multiplicative and
additive decompositions both reach `μ⋆ = 15` from the active prime
set `{3, 5, 7}`. -/
theorem master_factorisation_dual :
    muStar = 3 * 5 ∧ muStar = 3 + 5 + 7 :=
  ⟨master_factorisation_multiplicative, master_factorisation_additive⟩

/-! ### Section 3 — Master dimension statements

The persistence-active set has cardinality **exactly 3** on the
non-degenerate domain `μ ≥ 7`, and is **bounded by 3** everywhere. This
is the *dimension protection* property: the sieve dynamics does not
change the count of active primes. -/

/-- **Master dimension — cardinality at `μ⋆`.** `|persActive 15| = 3`. -/
theorem master_dimension_card_at_muStar :
    (persActive muStar).card = 3 := by
  rw [muStar_eq]; exact persActive_card_eq_three 15 (by decide)

/-- **Master dimension — constancy on `μ ≥ 7`.** For every `μ ≥ 7`,
`|persActive μ| = 3`. -/
theorem master_dimension_constancy (μ : ℕ) (h : 7 ≤ μ) :
    (persActive μ).card = 3 :=
  persActive_card_eq_three μ h

/-- **Master dimension — universal cap.** For every `μ`,
`|persActive μ| ≤ 3`. The cascade dimension never exceeds three. -/
theorem master_dimension_cap (μ : ℕ) : (persActive μ).card ≤ 3 :=
  persActive_card_le_three μ

/-- **Master dimension — explicit set at `μ⋆`.** `persActive 15 = {3,5,7}`. -/
theorem master_dimension_set_at_muStar :
    persActive muStar = ({3, 5, 7} : Finset ℕ) := by
  rw [muStar_eq]; exact persActive_eq_three_primes 15 (by decide)

/-- **Master dimension — sum identity.** `F_pers μ = ∑ p ∈ persActive μ, p`
for every `μ`. -/
theorem master_dimension_sum_identity (μ : ℕ) :
    Fpers μ = (persActive μ).sum id :=
  Fpers_eq_sum_persActive μ

/-- **Master dimension — protection bundle at the fixed point.** At
`μ⋆ = 15`, the persistence-active set is `{3, 5, 7}`, its cardinality is
`3`, and the sieve dynamics evaluates to `15`. -/
theorem master_dimension_protection_at_muStar :
    persActive muStar = ({3, 5, 7} : Finset ℕ)
    ∧ (persActive muStar).card = 3
    ∧ Fpers muStar = muStar := by
  refine ⟨master_dimension_set_at_muStar,
          master_dimension_card_at_muStar, ?_⟩
  rw [muStar_eq]
  exact Fpers_eq_fifteen_of_ge_seven 15 (by decide)

/-! ### Section 4 — Master totient resonance

The Euler totient of `μ⋆` is `8`, which coincides with the cardinality
of the admissible-residue set modulo `30 = 2 · μ⋆`. This is the
**arithmetic resonance** at the heart of the PT cascade. -/

/-- **Master totient — value at `μ⋆`.** `φ(μ⋆) = φ(15) = 8`. -/
theorem master_totient_value : Nat.totient muStar = 8 := by
  rw [muStar_eq]; exact nat_totient_15

/-- **Master totient — split over prime factors.**
`φ(μ⋆) = φ(3) · φ(5) = 2 · 4 = 8`. -/
theorem master_totient_split :
    Nat.totient muStar = Nat.totient 3 * Nat.totient 5 := by
  rw [muStar_eq]; exact totient_15_eq_factor_product

/-- **Master totient — resonance with admissible residues.**
`φ(μ⋆) = |admissibleResidues| = 8`. -/
theorem master_totient_admissible_resonance :
    Nat.totient muStar = admissibleResidues.card := by
  rw [muStar_eq]; exact totient_mu_star_eq_admissible_card

/-- **Master totient — doubling stability.** `φ(μ⋆) = φ(2 · μ⋆) = φ(30) = 8`. -/
theorem master_totient_doubling :
    Nat.totient muStar = Nat.totient 30 := by
  rw [muStar_eq]; exact totient_15_eq_totient_30

/-! ### Section 5 — Master cosmogonic / crystallisation statements

The **binary crystallisation** transition `μ = 17 → μ⋆ = 15` is the
PT-canonical mechanism that produces `μ⋆` from the raw four-prime sum
by removing the binary infrastructure prime `p = 2`. -/

/-- **Master crystallisation — raw fixed point.** The raw four-prime
sum is `2 + 3 + 5 + 7 = 17`. -/
theorem master_crystallisation_raw : (2 : ℕ) + 3 + 5 + 7 = 17 :=
  cosmogonic_seventeen

/-- **Master crystallisation — reduced fixed point = `μ⋆`.**
`3 + 5 + 7 = μ⋆ = 15`. -/
theorem master_crystallisation_reduced : (3 : ℕ) + 5 + 7 = muStar :=
  cosmogonic_fifteen

/-- **Master crystallisation — subtractive identity.**
`17 - p₁ = 17 - 2 = μ⋆`. -/
theorem master_crystallisation_subtract : (17 : ℕ) - 2 = muStar :=
  crystallisation_subtract

/-- **Master crystallisation — gap is `p₁`.** `17 - μ⋆ = p₁ = 2`. -/
theorem master_crystallisation_gap : (17 : ℕ) - muStar = 2 :=
  crystallisation_gap_is_p1

/-! ### Section 6 — The master ecosystem theorem

A single statement that bundles the canonical facts about `μ⋆ = 15`
into one theorem. Anyone who imports `FixedPointMasterTheorem` and
opens `PT.FixedPoint` can refer to this as "the master `μ⋆` theorem". -/

/-- **Master ecosystem — `μ⋆ = 15`.** The complete arithmetic
characterisation of `μ⋆`:

1. **Identity.** `μ⋆ = 15`.
2. **Fixed point.** `F_pers μ⋆ = μ⋆`.
3. **Multiplicative factorisation.** `μ⋆ = 3 · 5`.
4. **Additive factorisation.** `μ⋆ = 3 + 5 + 7`.
5. **Prime factors.** `μ⋆.primeFactors = {3, 5}`.
6. **Active set.** `persActive μ⋆ = {3, 5, 7}`.
7. **Dimension.** `|persActive μ⋆| = 3`.
8. **Dimension cap.** `|persActive μ| ≤ 3` for all `μ`.
9. **Totient.** `φ(μ⋆) = 8 = |admissibleResidues|`.
10. **Crystallisation.** `μ⋆ = 17 - 2` (binary infrastructure removed). -/
theorem fixed_point_master_ecosystem :
    muStar = 15
    ∧ Fpers muStar = muStar
    ∧ muStar = 3 * 5
    ∧ muStar = 3 + 5 + 7
    ∧ muStar.primeFactors = ({3, 5} : Finset ℕ)
    ∧ persActive muStar = ({3, 5, 7} : Finset ℕ)
    ∧ (persActive muStar).card = 3
    ∧ (∀ μ, (persActive μ).card ≤ 3)
    ∧ Nat.totient muStar = 8
    ∧ Nat.totient muStar = admissibleResidues.card
    ∧ muStar = 17 - 2 := by
  refine ⟨muStar_eq, master_fixedpoint_existence,
          master_factorisation_multiplicative,
          master_factorisation_additive,
          master_factorisation_primeFactors,
          master_dimension_set_at_muStar,
          master_dimension_card_at_muStar,
          master_dimension_cap,
          master_totient_value,
          master_totient_admissible_resonance,
          crystallisation_subtract.symm⟩

/-- **Master ecosystem — uniqueness companion.** On the non-degenerate
range `μ ≥ 8`, `μ⋆ = 15` is the unique fixed point of `F_pers`; on the
full range `μ ≥ 2` the only other fixed point is the degenerate
`μ = 3`. -/
theorem fixed_point_master_uniqueness :
    (∀ μ, 8 ≤ μ → (Fpers μ = μ ↔ μ = 15))
    ∧ (∀ μ, 2 ≤ μ → Fpers μ = μ → μ = 3 ∨ μ = 15) :=
  ⟨T7_fixedpoint_iff, T7_all_fixed_points⟩

/-- **Master summary — final aggregated headline.** Bundles the
ecosystem theorem with the uniqueness companion. This is the single
top-level theorem of the file: the complete arithmetic specification of
the PT canonical integer `μ⋆ = 15`. -/
theorem fixed_point_master_summary :
    -- Identity & fixed point
    muStar = 15
    ∧ Fpers muStar = muStar
    -- Factorisation (multiplicative ∧ additive)
    ∧ muStar = 3 * 5
    ∧ muStar = 3 + 5 + 7
    -- Dimension protection
    ∧ persActive muStar = ({3, 5, 7} : Finset ℕ)
    ∧ (persActive muStar).card = 3
    ∧ (∀ μ, (persActive μ).card ≤ 3)
    -- Totient resonance
    ∧ Nat.totient muStar = 8
    ∧ Nat.totient muStar = admissibleResidues.card
    -- Crystallisation
    ∧ muStar = 17 - 2
    -- Uniqueness
    ∧ (∀ μ, 8 ≤ μ → (Fpers μ = μ ↔ μ = 15))
    ∧ (∀ μ, 2 ≤ μ → Fpers μ = μ → μ = 3 ∨ μ = 15) := by
  refine ⟨muStar_eq, master_fixedpoint_existence,
          master_factorisation_multiplicative,
          master_factorisation_additive,
          master_dimension_set_at_muStar,
          master_dimension_card_at_muStar,
          master_dimension_cap,
          master_totient_value,
          master_totient_admissible_resonance,
          crystallisation_subtract.symm,
          T7_fixedpoint_iff, T7_all_fixed_points⟩

end PT.FixedPoint
