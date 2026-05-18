/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Sieve.PrimitiveRootCascade
import PT.Sieve.PrimorialCoprime
import Mathlib.Data.Nat.Totient
import Mathlib.Tactic

/-!
# Euler-totient identities for the PT cascade

This file collects the **closed-form values** of `Nat.totient` on the
arithmetic skeleton of the Persistence Theory:

* the PT-active primes `{3, 5, 7}` and their immediate neighbours
  `{2, 11, 13}` (boundary of the active cascade);
* the **third primorial** `30 = 2·3·5`;
* the **active-prime primorial** `105 = 3·5·7`;
* the **canonical dimension** `μ⋆ = 15 = 3·5`.

It complements:

* `PT/Sieve/PrimorialCoprime.lean` — the identity `φ(30) = 8 = |R|`;
* `PT/Sieve/PrimitiveRootCascade.lean` — the `phiPT` table on `{3, 5, 7}`
  and the strict cascade `φ(3) < φ(5) < φ(7)`.

All values are computed by `decide`, which is sound for `Nat.totient` on
small explicit arguments. The multiplicativity identities
`φ(3·5) = φ(3)·φ(5)`, `φ(3·5·7) = φ(3)·φ(5)·φ(7)`, etc., are stated as
plain numerical equalities and discharged by `decide` rather than by
invoking `Nat.totient_mul` (which would require explicit coprimality
hypotheses). This keeps the file self-contained.

## Reference

Monograph Appendix N, sections on the third primorial and the
dimensional cascade `μ⋆ = 15`.
-/

namespace PT.Sieve

/-! ### Explicit values of `Nat.totient` on the PT skeleton

These pin down `Nat.totient` on every integer that appears in the PT
cascade: the active primes `3, 5, 7`, their neighbours `2, 11, 13`,
the dimension `μ⋆ = 15`, the third primorial `30`, and the
active-prime primorial `105`. -/

theorem nat_totient_2   : Nat.totient 2   = 1   := by decide
theorem nat_totient_3   : Nat.totient 3   = 2   := by decide
theorem nat_totient_5   : Nat.totient 5   = 4   := by decide
theorem nat_totient_7   : Nat.totient 7   = 6   := by decide
theorem nat_totient_11  : Nat.totient 11  = 10  := by decide
theorem nat_totient_13  : Nat.totient 13  = 12  := by decide
theorem nat_totient_15  : Nat.totient 15  = 8   := by decide
theorem nat_totient_105 : Nat.totient 105 = 48  := by decide

/-! Re-export of `nat_totient_30 = 8` from `PrimorialCoprime` as a
local convenience alias. -/

theorem nat_totient_30' : Nat.totient 30 = 8 := nat_totient_30

/-! ### Multiplicativity on coprime factors of the PT skeleton

For each multiplicative split of an integer in the PT skeleton into
pairwise coprime parts, we record the corresponding factorisation of
the totient. We discharge each by `decide` (both sides are concrete
small naturals). -/

/-- `φ(3·5) = φ(3) · φ(5) = 2 · 4 = 8`. -/
theorem totient_mul_3_5 :
    Nat.totient (3 * 5) = Nat.totient 3 * Nat.totient 5 := by decide

/-- `φ(2·3) = φ(2) · φ(3) = 1 · 2 = 2`. -/
theorem totient_mul_2_3 :
    Nat.totient (2 * 3) = Nat.totient 2 * Nat.totient 3 := by decide

/-- `φ(2·15) = φ(2) · φ(15) = 1 · 8 = 8 = φ(30)`. -/
theorem totient_mul_2_15 :
    Nat.totient (2 * 15) = Nat.totient 2 * Nat.totient 15 := by decide

/-- `φ(3·5·7) = φ(3) · φ(5) · φ(7) = 2 · 4 · 6 = 48`. -/
theorem totient_mul_3_5_7 :
    Nat.totient (3 * 5 * 7) = Nat.totient 3 * Nat.totient 5 * Nat.totient 7 := by
  decide

/-- `φ(2·3·5) = φ(2) · φ(3) · φ(5) = 1 · 2 · 4 = 8`. -/
theorem totient_mul_2_3_5 :
    Nat.totient (2 * 3 * 5) = Nat.totient 2 * Nat.totient 3 * Nat.totient 5 := by
  decide

/-! ### `μ⋆ = 15` and the third primorial `30`

These pinpoint the two key identifications:

* `μ⋆ = 15 = 3 · 5` and `φ(μ⋆) = φ(3) · φ(5) = 2 · 4 = 8`;
* `30 = 2 · 15` and `φ(30) = φ(2) · φ(15) = 1 · 8 = 8`.

In particular `φ(μ⋆) = φ(30) = 8`: the doubling `15 → 30` does **not**
change the totient (because `φ(2) = 1`). -/

/-- **`μ⋆ = 15` from `3 · 5`.** -/
theorem totient_15_eq_factor_product :
    Nat.totient 15 = Nat.totient 3 * Nat.totient 5 := by decide

/-- **The PT canonical-dimension identity.**
    `φ(μ⋆) = φ(15) = φ(3) · φ(5) = 2 · 4 = 8`. -/
theorem totient_mu_star : Nat.totient 15 = 2 * 4 := by decide

/-- **Doubling identity.** `φ(15) = φ(30)` and both equal `8`. -/
theorem totient_15_eq_totient_30 :
    Nat.totient 15 = Nat.totient 30 := by decide

/-- **Ratio identity.** `φ(30) / φ(2) = φ(15)`. (Stated as
    `φ(30) = φ(2) · φ(15)`, with `φ(2) = 1`.) -/
theorem totient_30_eq_totient_2_mul_15 :
    Nat.totient 30 = Nat.totient 2 * Nat.totient 15 := by decide

/-- **Active-prime primorial identity.**
    `φ(105) = φ(3) · φ(5) · φ(7) = 48`. -/
theorem totient_105_eq_active_prime_product :
    Nat.totient 105 = Nat.totient 3 * Nat.totient 5 * Nat.totient 7 := by
  decide

/-! ### Strict cascade on the PT-active primes and their neighbours

`φ` is strictly increasing on `2 < 3 < 5 < 7 < 11 < 13` because each
move to the next prime adds `≥ 2` to the totient. -/

/-- **Strict totient cascade on the prime skeleton.**
    `φ(2) < φ(3) < φ(5) < φ(7) < φ(11) < φ(13)`, i.e.
    `1 < 2 < 4 < 6 < 10 < 12`. -/
theorem totient_cascade_strict :
    Nat.totient 2  < Nat.totient 3
    ∧ Nat.totient 3  < Nat.totient 5
    ∧ Nat.totient 5  < Nat.totient 7
    ∧ Nat.totient 7  < Nat.totient 11
    ∧ Nat.totient 11 < Nat.totient 13 := by
  refine ⟨?_, ?_, ?_, ?_, ?_⟩ <;> decide

/-- **Strict totient cascade, numerical form.** -/
theorem totient_cascade_values :
    Nat.totient 2  = 1
    ∧ Nat.totient 3  = 2
    ∧ Nat.totient 5  = 4
    ∧ Nat.totient 7  = 6
    ∧ Nat.totient 11 = 10
    ∧ Nat.totient 13 = 12 :=
  ⟨nat_totient_2, nat_totient_3, nat_totient_5, nat_totient_7,
   nat_totient_11, nat_totient_13⟩

/-! ### Bridge to `phiPT` (PrimitiveRootCascade)

The bespoke partial function `phiPT : ℕ → ℕ` from
`PrimitiveRootCascade.lean` agrees with `Nat.totient` on each
PT-active prime. -/

/-- **Compatibility of `phiPT` with `Nat.totient` on active primes.** -/
theorem phiPT_eq_nat_totient_active :
    phiPT 3 = Nat.totient 3
    ∧ phiPT 5 = Nat.totient 5
    ∧ phiPT 7 = Nat.totient 7 := by
  refine ⟨?_, ?_, ?_⟩ <;> decide

/-! ### μ⋆-residue identity

The headline identity: the canonical dimension `μ⋆ = 15` satisfies
`φ(μ⋆) = 8`, which is exactly the cardinality `|R|` of the admissible
residues modulo `30` (themselves equal to `(ℤ/30ℤ)*`). So the
**number of admissible residues** is determined twice over by the
arithmetic: once as `|R| = 8` (combinatorial), once as
`φ(15) = φ(μ⋆) = 8` (number-theoretic). -/

/-- **μ⋆-residue identity (numerical form).**
    `φ(μ⋆) = φ(15) = 8 = |R| = |admissibleResidues|`. -/
theorem totient_mu_star_eq_admissible_card :
    Nat.totient 15 = admissibleResidues.card := by
  rw [nat_totient_15]
  decide

/-- **μ⋆ vs third primorial identity (chain).**
    `φ(μ⋆) = φ(15) = φ(30) = |admissibleResidues| = 8`. -/
theorem totient_mu_star_chain :
    Nat.totient 15 = Nat.totient 30
    ∧ Nat.totient 30 = admissibleResidues.card
    ∧ admissibleResidues.card = 8 := by
  refine ⟨totient_15_eq_totient_30, ?_, ?_⟩
  · exact admissibleResidues_card_eq_totient.symm
  · decide

/-! ### Headline -/

/-- **Headline (PT totient cascade).** All the explicit totient values
    needed by Persistence Theory:

    * Active primes: `φ(3) = 2`, `φ(5) = 4`, `φ(7) = 6`.
    * Neighbours: `φ(2) = 1`, `φ(11) = 10`, `φ(13) = 12`.
    * Canonical dimension: `φ(μ⋆) = φ(15) = 8`.
    * Third primorial: `φ(30) = 8`.
    * Active-prime primorial: `φ(105) = 48`.
    * Multiplicativity: `φ(3·5·7) = φ(3)·φ(5)·φ(7) = 48`.
    * Bridge to admissibles: `φ(μ⋆) = |R| = 8`.
    * Strict cascade: `1 < 2 < 4 < 6 < 10 < 12`. -/
theorem pt_totient_cascade_summary :
    -- Active primes
    Nat.totient 3 = 2
    ∧ Nat.totient 5 = 4
    ∧ Nat.totient 7 = 6
    -- Neighbours
    ∧ Nat.totient 2  = 1
    ∧ Nat.totient 11 = 10
    ∧ Nat.totient 13 = 12
    -- μ⋆ and primorials
    ∧ Nat.totient 15  = 8
    ∧ Nat.totient 30  = 8
    ∧ Nat.totient 105 = 48
    -- Multiplicativity on coprime triple
    ∧ Nat.totient (3 * 5 * 7)
        = Nat.totient 3 * Nat.totient 5 * Nat.totient 7
    -- μ⋆-residue identity
    ∧ Nat.totient 15 = admissibleResidues.card
    -- Strict cascade
    ∧ Nat.totient 2  < Nat.totient 3
    ∧ Nat.totient 3  < Nat.totient 5
    ∧ Nat.totient 5  < Nat.totient 7
    ∧ Nat.totient 7  < Nat.totient 11
    ∧ Nat.totient 11 < Nat.totient 13 := by
  refine ⟨nat_totient_3, nat_totient_5, nat_totient_7,
          nat_totient_2, nat_totient_11, nat_totient_13,
          nat_totient_15, nat_totient_30', nat_totient_105,
          totient_105_eq_active_prime_product,
          totient_mu_star_eq_admissible_card,
          ?_, ?_, ?_, ?_, ?_⟩ <;> decide

end PT.Sieve
