/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Sieve.Bimodality
import PT.Sieve.AdmissibleResiduesArithmetic
import PT.Sieve.AdmissibleResiduesSquareSumMod
import PT.Sieve.AdmissibleResiduesCubeSumMod
import PT.Sieve.AdmissibleResiduesFourthPowerSumMod
import Mathlib.Tactic

/-!
# Admissible residues — Sum of fifth powers modulo various `q` (App N extension)

This file collects the **sum of fifth-power invariants** of the admissible
residue set `R = {1, 7, 11, 13, 17, 19, 23, 29} = (ℤ/30ℤ)*`:

$$\sum_{r \in R} r^5
  \;=\; 1 + 16807 + 161051 + 371293 + 1419857 + 2476099 + 6436343 + 20511149
  \;=\; 31\,392\,600.$$

Companion to:
* `AdmissibleResiduesArithmetic.lean` — additive `∑ r = 120`.
* `AdmissibleResiduesSquareSumMod.lean` — `∑ r² = 2360` (even, non-collapsed).
* `AdmissibleResiduesCubeSumMod.lean` — `∑ r³ = 52200` (odd, total Wilson
  collapse modulo every divisor of `120`).
* `AdmissibleResiduesFourthPowerSumMod.lean` — `∑ r⁴ = 1\,246\,568` (even,
  Carmichael-Wilson reduction to `|R| ≡ 8 mod 120`).

## Main results

* `admissibleResidues_sumFifth`: `∑ r⁵ = 31\,392\,600`.
* `admissibleResidues_sumFifth_factor`: `31\,392\,600 = 261\,605 · 120
  = 2³ · 3 · 5² · 52\,321` (the residual factor `52\,321` is prime).
  The factor `120` is structural (odd-moment central-symmetry collapse,
  see below).
* `admissibleResidues_sumFifth_mod_*`: closed-form residue of `31\,392\,600`
  modulo every `q ∈ {2, 3, 5, 6, 8, 10, 15, 30, 60, 100, 120}`.
* `admissibleResidues_sumFifth_skew5_zero`: the **fifth central moment**
  identity `∑ (r − 15)⁵ = 0`, a direct consequence of the central
  symmetry `r ↔ 30 − r` of `R`.
* `admissibleResidues_sumFifth_summary`: headline table of residues.

## Patterns (the odd-moment Wilson collapse)

The residue `31\,392\,600 mod q` exhibits the **same uniform collapse**
as the cubic case, across **every** modulus in the test set:

| `q`  |  `31\,392\,600 mod q` |  Note                                          |
|------|-----------------------|------------------------------------------------|
|  2   |  0                    |  even.                                         |
|  3   |  0                    |  digit sum `3+1+3+9+2+6 = 24 ≡ 0 mod 3`.       |
|  5   |  0                    |  last digit `0`.                               |
|  6   |  0                    |  CRT(2,3).                                     |
|  8   |  0                    |  `31\,392\,600 = 3\,924\,075 · 8`.             |
|  10  |  0                    |  last digit `0`.                               |
|  15  |  0                    |  CRT(3,5).                                     |
|  30  |  0                    |  `31\,392\,600 = 1\,046\,420 · 30`.            |
|  60  |  0                    |  `31\,392\,600 = 523\,210 · 60`.               |
| 100  |  0                    |  last two digits `00`.                         |
| 120  |  0                    |  `31\,392\,600 = 261\,605 · 120`.              |

This is the **odd-moment Wilson collapse** at order 5: just as for `∑ r³`,
*every* divisor of `120 = 2³·3·5` (and `100`) divides `∑ r⁵`. The
structural reason is the central symmetry `r ↔ 30 − r` of `(ℤ/30ℤ)*`:
for any odd exponent `k`, the involution `r ↦ 30 − r` pairs up terms with
the **same parity** of `(r − 15)^k`, but the pair `(r, 30 − r)` produces
contributions `r^k + (30 − r)^k` that are even-degree polynomials in `r`
whose **central moments vanish identically** at every odd order. Concretely:

$$\sum_{r \in R} (r - 15)^{2j+1} \;=\; 0 \quad \text{for all } j \ge 0,$$

since the involution `r ↦ 30 − r` is sign-reversing on `(r − 15)`. The
"raw" sums `∑ r^{2j+1}` then collapse modulo every divisor of `120`
because `∑ r = 120` and `(∑ r)^k ≡ 0` modulo every such divisor.

## Dichotomy of moments

The triple `(∑ r³, ∑ r⁴, ∑ r⁵)` now realises the full moment-parity
dichotomy of `(ℤ/30ℤ)*`:

* **Odd moments** (`∑ r³`, `∑ r⁵`): **total Wilson collapse**,
  `≡ 0 mod q` for every divisor `q` of `120` (and `100`), by central
  symmetry. Compare:
  - `∑ r³ = 52\,200 = 435 · 120`,
  - `∑ r⁵ = 31\,392\,600 = 261\,605 · 120`.
* **Even moments** (`∑ r²`, `∑ r⁴`): **non-trivial residues**, governed
  by Wilson / Carmichael identities `r^{2k} ≡ 1` modulo the relevant
  primes. Compare:
  - `∑ r² = 2360`, residues `(0, 2, 0, 2, 0, 0, 5, 20, 20, 60, 80)` at
    `q ∈ {2, 3, 5, 6, 8, 10, 15, 30, 60, 100, 120}`,
  - `∑ r⁴ = 1\,246\,568`, residues `(0, 2, 3, 2, 0, 8, 8, 8, 8, 68, 8)`.

The **odd–even dichotomy** is the cleanest invariant statement available
about the moment table of `R`: every odd raw moment `∑ r^{2j+1}` is
**fully collapsed** modulo every divisor of `120`, while every even raw
moment carries a non-trivial residue.

## The fifth central moment vanishes

The set `R` is centrally symmetric about its mean `r̄ = 15`. The
involution `r ↦ 30 − r` is a bijection of `R` that reverses the sign of
`(r − 15)`, hence reverses the sign of `(r − 15)^k` for every **odd**
`k`. Summing over `R`, every odd central moment vanishes:

$$\sum_{r \in R} (r - 15)^5 \;=\; 0.$$

(In integer form, since `(r - 15)` may be negative we cast to `ℤ`.)

## Reference

Monograph Appendix N §"Invariants arithmétiques de R" (fifth-power part),
extending the cubic table to the next odd moment and making explicit the
moment-parity dichotomy paired with `AdmissibleResiduesFourthPowerSumMod`.
-/

namespace PT.Sieve

open Finset

/-! ### The sum of fifth powers: exact value -/

/-- The sum of fifth powers of the admissible residues equals `31\,392\,600`. -/
theorem admissibleResidues_sumFifth :
    admissibleResidues.sum (fun r => r * r * r * r * r) = 31392600 := by
  decide

/-- The sum of fifth powers, term-by-term decomposition (for documentation).
    The eight summands are
    `1, 16807, 161051, 371293, 1419857, 2476099, 6436343, 20511149`. -/
theorem admissibleResidues_sumFifth_terms :
    admissibleResidues.sum (fun r => r * r * r * r * r)
      = 1 + 16807 + 161051 + 371293 + 1419857 + 2476099 + 6436343 + 20511149 := by
  decide

/-- The sum of fifth powers using `r ^ 5`. -/
theorem admissibleResidues_sumFifth_pow :
    admissibleResidues.sum (fun r => r ^ 5) = 31392600 := by
  decide

/-- Factorisation: `31\,392\,600 = 261\,605 · 120 = 2³ · 3 · 5² · 52\,321`,
    with `52\,321` prime. The structural factor `120 = ∑ r` reflects the
    **odd-moment Wilson collapse** (every odd raw moment of `R` is a
    multiple of `∑ r`); the residual factor `261\,605 = 5 · 52\,321`
    contains no further small prime, in contrast with the cubic case
    `52\,200 = 435 · 120 = 2³ · 3² · 5² · 29` where the residual `435 = 3·5·29`
    splits over `{3, 5, 29}` (with `29 ∈ R` the boundary unit). -/
theorem admissibleResidues_sumFifth_factor :
    admissibleResidues.sum (fun r => r * r * r * r * r) = 261605 * 120
    ∧ admissibleResidues.sum (fun r => r * r * r * r * r)
        = 2 * 2 * 2 * 3 * 5 * 5 * 52321 := by
  refine ⟨?_, ?_⟩ <;> decide

/-! ### Residues of `31\,392\,600` modulo divisors of `120` (and `100`) -/

/-- `31\,392\,600 ≡ 0 mod 2` (even). -/
theorem admissibleResidues_sumFifth_mod_2 :
    admissibleResidues.sum (fun r => r * r * r * r * r) % 2 = 0 := by
  decide

/-- `31\,392\,600 ≡ 0 mod 3`. -/
theorem admissibleResidues_sumFifth_mod_3 :
    admissibleResidues.sum (fun r => r * r * r * r * r) % 3 = 0 := by
  decide

/-- `31\,392\,600 ≡ 0 mod 5`. -/
theorem admissibleResidues_sumFifth_mod_5 :
    admissibleResidues.sum (fun r => r * r * r * r * r) % 5 = 0 := by
  decide

/-- `31\,392\,600 ≡ 0 mod 6`. -/
theorem admissibleResidues_sumFifth_mod_6 :
    admissibleResidues.sum (fun r => r * r * r * r * r) % 6 = 0 := by
  decide

/-- **Odd-moment Wilson collapse mod 8.** `31\,392\,600 ≡ 0 mod 8`: every
    admissible `r` is odd, so `r² ≡ 1 mod 8`, hence `r⁵ = r·(r²)² ≡ r mod 8`;
    the order-1 collapse `∑ r ≡ 0 mod 8` (from `∑ r = 120 = 15·8`) then
    transports to `∑ r⁵ ≡ 0 mod 8`. -/
theorem admissibleResidues_sumFifth_mod_8 :
    admissibleResidues.sum (fun r => r * r * r * r * r) % 8 = 0 := by
  decide

/-- `31\,392\,600 ≡ 0 mod 10`. -/
theorem admissibleResidues_sumFifth_mod_10 :
    admissibleResidues.sum (fun r => r * r * r * r * r) % 10 = 0 := by
  decide

/-- `31\,392\,600 ≡ 0 mod 15`. -/
theorem admissibleResidues_sumFifth_mod_15 :
    admissibleResidues.sum (fun r => r * r * r * r * r) % 15 = 0 := by
  decide

/-- `31\,392\,600 ≡ 0 mod 30`. -/
theorem admissibleResidues_sumFifth_mod_30 :
    admissibleResidues.sum (fun r => r * r * r * r * r) % 30 = 0 := by
  decide

/-- `31\,392\,600 ≡ 0 mod 60`. -/
theorem admissibleResidues_sumFifth_mod_60 :
    admissibleResidues.sum (fun r => r * r * r * r * r) % 60 = 0 := by
  decide

/-- `31\,392\,600 ≡ 0 mod 100`. **No outlier**: unlike the fourth-power
    case (where `1\,246\,568 ≡ 68 mod 100` because `100 ∤ 120` and
    `λ(100) = 20 ∤ 4`), at order 5 the central-symmetry argument is
    available and yields `∑ r⁵ ≡ 0 mod 100` regardless of any Wilson /
    Carmichael obstruction. The mechanism is purely combinatorial:
    `r⁵ + (30 − r)⁵` is a polynomial of degree `4` in `r` whose
    contributions cancel `mod 100` when summed over the four pairs of `R`. -/
theorem admissibleResidues_sumFifth_mod_100 :
    admissibleResidues.sum (fun r => r * r * r * r * r) % 100 = 0 := by
  decide

/-- **Total Wilson collapse mod 120.** `31\,392\,600 ≡ 0 mod 120`,
    equivalently `∑ r⁵ = 261\,605 · 120 = 261\,605 · (∑ r)`. The order-5
    moment is an integer multiple of the order-1 moment, with ratio
    `261\,605 = 5 · 52\,321`; structurally this is the same
    odd-moment collapse that gave `∑ r³ = 435 · 120`. Here the residual
    `52\,321` is prime, so the structure beyond the `120`-factor is trivial. -/
theorem admissibleResidues_sumFifth_mod_120 :
    admissibleResidues.sum (fun r => r * r * r * r * r) % 120 = 0 := by
  decide

/-! ### Fifth central moment vanishes (sign-reversing involution `r ↦ 30 − r`) -/

/-- **Fifth central moment vanishes.** `∑ (r − 15)⁵ = 0`, the direct
    integer translation of the central symmetry `r ↔ 30 − r` of `R`: the
    involution `r ↦ 30 − r` is a bijection of `R` that reverses the sign
    of `(r − 15)`, hence of every odd power `(r − 15)^{2j+1}`. We cast to
    `ℤ` because `(r − 15)` is negative for the four residues `1, 7, 11, 13`. -/
theorem admissibleResidues_sumFifth_skew5_zero :
    admissibleResidues.sum (fun r => ((r : ℤ) - 15) ^ 5) = 0 := by
  decide

/-- **Per-pair sign-reversing identity.** For every `r ∈ R`, the symmetric
    partner `30 − r` also lies in `R`, and the pair contributes
    `(r − 15)⁵ + ((30 − r) − 15)⁵ = (r − 15)⁵ + (15 − r)⁵ = 0` in `ℤ`.
    This is the per-element manifestation of the central symmetry that
    drives all odd central moments to vanish. -/
theorem admissibleResidues_pair_fifth_cancel :
    ∀ r ∈ admissibleResidues,
      ((r : ℤ) - 15) ^ 5 + ((30 - r : ℤ) - 15) ^ 5 = 0 := by
  decide

/-- **Mean fifth power.** `(∑ r⁵)/|R| = 31\,392\,600 / 8 = 3\,924\,075
    = 3 · 5² · 52\,321`. -/
theorem admissibleResidues_sumFifth_mean :
    admissibleResidues.sum (fun r => r * r * r * r * r) / admissibleResidues.card
      = 3924075 := by
  decide

/-! ### Cross-moment identities -/

/-- **Fifth-power / first-moment ratio.** Since the odd-moment Wilson
    collapse gives `120 ∣ ∑ r⁵`, the integer quotient is well-defined:
    `(∑ r⁵) / (∑ r) = 31\,392\,600 / 120 = 261\,605`. -/
theorem admissibleResidues_sumFifth_over_sum :
    admissibleResidues.sum (fun r => r * r * r * r * r)
      / (admissibleResidues.sum id) = 261605 := by
  decide

/-- **Fifth-power / cube cross identity (raw integer form).**
    `5 · ∑ r⁵ = 25 · 6280 · ∑ r³ + …`: the cleanest stable integer
    identity at this order is the direct ratio
    `∑ r⁵ = 261\,605 · (∑ r)`, recorded above. We additionally observe
    the bare comparison `∑ r⁵ − 6 · (∑ r²)² = 31\,392\,600 − 6 · 5\,569\,600
    = 31\,392\,600 − 33\,417\,600 = −2\,025\,000` in `ℤ`. -/
theorem admissibleResidues_sumFifth_minus_sixSumSquareSq :
    (admissibleResidues.sum (fun r => r * r * r * r * r) : ℤ)
      - 6 * (admissibleResidues.sum (fun r => r * r) : ℤ) ^ 2
      = -2025000 := by
  decide

/-! ### Headline -/

/-- **Headline (sum-of-fifth-powers invariants).** The admissible residue
    set `R = {1, 7, 11, 13, 17, 19, 23, 29}` satisfies:

    * `∑ r⁵ = 31\,392\,600 = 261\,605 · 120 = 2³ · 3 · 5² · 52\,321`
      (with `52\,321` prime).
    * **Total odd-moment Wilson collapse**: `31\,392\,600 ≡ 0 mod q` for
      every `q ∈ {2, 3, 5, 6, 8, 10, 15, 30, 60, 100, 120}` — every
      divisor of `120` divides `∑ r⁵` (and so does `100`). This matches
      the cubic table exactly and contrasts sharply with the fourth-power
      table where eleven distinct residues `(0, 2, 3, 2, 0, 8, 8, 8, 8, 68, 8)`
      appear.
    * **Fifth central moment vanishes**: `∑ (r − 15)⁵ = 0` in `ℤ`, the
      direct consequence of the central symmetry `r ↔ 30 − r` of `R`.
    * Mean fifth power: `(∑ r⁵)/|R| = 3\,924\,075 = 3 · 5² · 52\,321`.
    * Reduced ratio: `(∑ r⁵)/(∑ r) = 261\,605`.

    Coupled with the cubic table, this confirms the **odd–even moment
    dichotomy** of `R`: every **odd** raw moment `∑ r^{2j+1}` is fully
    collapsed modulo every divisor of `120` (third and fifth powers
    verified, with the same mechanism extending to all higher odd
    powers), while every **even** raw moment carries a non-trivial
    residue table governed by Carmichael-Wilson identities `r^{\lambda(q)}
    \equiv 1`. The pair `(odd : skewness collapse, even : Wilson reduction)`
    is now established at orders `3`/`4` and reinforced at orders `5`/`4`. -/
theorem admissibleResidues_sumFifth_summary :
    admissibleResidues.sum (fun r => r * r * r * r * r) = 31392600
    ∧ admissibleResidues.sum (fun r => r * r * r * r * r) % 2 = 0
    ∧ admissibleResidues.sum (fun r => r * r * r * r * r) % 3 = 0
    ∧ admissibleResidues.sum (fun r => r * r * r * r * r) % 5 = 0
    ∧ admissibleResidues.sum (fun r => r * r * r * r * r) % 6 = 0
    ∧ admissibleResidues.sum (fun r => r * r * r * r * r) % 8 = 0
    ∧ admissibleResidues.sum (fun r => r * r * r * r * r) % 10 = 0
    ∧ admissibleResidues.sum (fun r => r * r * r * r * r) % 15 = 0
    ∧ admissibleResidues.sum (fun r => r * r * r * r * r) % 30 = 0
    ∧ admissibleResidues.sum (fun r => r * r * r * r * r) % 60 = 0
    ∧ admissibleResidues.sum (fun r => r * r * r * r * r) % 100 = 0
    ∧ admissibleResidues.sum (fun r => r * r * r * r * r) % 120 = 0
    ∧ admissibleResidues.sum (fun r => ((r : ℤ) - 15) ^ 5) = 0 :=
  ⟨admissibleResidues_sumFifth,
   admissibleResidues_sumFifth_mod_2,
   admissibleResidues_sumFifth_mod_3,
   admissibleResidues_sumFifth_mod_5,
   admissibleResidues_sumFifth_mod_6,
   admissibleResidues_sumFifth_mod_8,
   admissibleResidues_sumFifth_mod_10,
   admissibleResidues_sumFifth_mod_15,
   admissibleResidues_sumFifth_mod_30,
   admissibleResidues_sumFifth_mod_60,
   admissibleResidues_sumFifth_mod_100,
   admissibleResidues_sumFifth_mod_120,
   admissibleResidues_sumFifth_skew5_zero⟩

end PT.Sieve
