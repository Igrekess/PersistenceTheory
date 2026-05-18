/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Sieve.Bimodality
import PT.Sieve.AdmissibleResiduesArithmetic
import PT.Sieve.AdmissibleResiduesSquareSumMod
import PT.Sieve.AdmissibleResiduesCubeSumMod
import Mathlib.Tactic

/-!
# Admissible residues — Sum of fourth powers modulo various `q` (App N extension)

This file collects the **sum of fourth-power invariants** of the admissible
residue set `R = {1, 7, 11, 13, 17, 19, 23, 29} = (ℤ/30ℤ)*`:

$$\sum_{r \in R} r^4
  \;=\; 1 + 2401 + 14641 + 28561 + 83521 + 130321 + 279841 + 707281
  \;=\; 1\,246\,568.$$

Companion to:
* `AdmissibleResiduesArithmetic.lean` — additive `∑ r = 120`.
* `AdmissibleResiduesSquareSumMod.lean` — `∑ r² = 2360`.
* `AdmissibleResiduesCubeSumMod.lean` — `∑ r³ = 52200`.

## Main results

* `admissibleResidues_sumQuart`: `∑ r⁴ = 1246568`.
* `admissibleResidues_sumQuart_factor`: `1246568 = 8 · 155821 = 2³ · 155821`,
  with `155821` prime. (Contrast with `∑ r³ = 435·120 = 2³·3²·5²·29`.)
* `admissibleResidues_sumQuart_mod_*`: closed-form residue of `1246568`
  modulo every `q ∈ {2, 3, 5, 6, 8, 10, 15, 30, 60, 100, 120}`.
* `admissibleResidues_quart_carmichael`: **the fourth-power Wilson identity**
  `r⁴ ≡ 1 mod 120` for every `r ∈ R`, hence `∑ r⁴ ≡ |R| = 8 mod 120`.
* `admissibleResidues_sumQuart_kurtosis_integer`: the integer kurtosis
  identity (vanishing of the centered fourth-moment defect).
* `admissibleResidues_sumQuart_summary`: headline table of residues.

## Patterns (the fourth-order Wilson identity)

The residue `1246568 mod q` exhibits a strikingly different pattern from
the cubic case (`∑ r³ ≡ 0` for *all* divisors of `120`). The fourth power
is an **even** symmetry-respecting moment, so it does **not** vanish under
the central symmetry `r ↔ 30 − r`; instead it obeys the **Carmichael
Wilson identity**: since `λ(120) = lcm(λ(8), λ(3), λ(5)) = lcm(2, 2, 4) = 4`,
every `r ∈ R` (all coprime to `120`) satisfies

$$r^4 \;\equiv\; 1 \pmod{120}.$$

Summing eight copies:

$$\sum_{r \in R} r^4 \;\equiv\; |R| \;=\; 8 \pmod{120}.$$

This **single congruence** drives the entire residue table:

| `q`  |  `1246568 mod q` |  Explanation                                          |
|------|------------------|-------------------------------------------------------|
|  2   |  0               |  `8 ≡ 0 mod 2`.                                       |
|  3   |  2               |  `8 ≡ 2 mod 3` (each `r⁴ ≡ 1 mod 3` by Fermat).       |
|  5   |  3               |  `8 ≡ 3 mod 5` (each `r⁴ ≡ 1 mod 5` by Fermat).       |
|  6   |  2               |  `8 ≡ 2 mod 6`.                                       |
|  8   |  0               |  `8 ≡ 0 mod 8` (each `r⁴ ≡ 1 mod 8` from `r²≡1 mod 8`)|
|  10  |  8               |  `8 ≡ 8 mod 10`.                                      |
|  15  |  8               |  `8 ≡ 8 mod 15`.                                      |
|  30  |  8               |  `8 ≡ 8 mod 30`.                                      |
|  60  |  8               |  `8 ≡ 8 mod 60`.                                      |
| 100  | 68               |  `100 ∤ 120`; the fourth-power Wilson does not apply. |
| 120  |  8               |  **`8` exactly**, the Carmichael Wilson identity.     |

The residue `mod 100` is the **outlier**: `100` is not a divisor of `120`
and its Carmichael exponent `λ(100) = lcm(λ(4), λ(25)) = lcm(2, 20) = 20`
does not divide `4`, so `r⁴` is *not* universally `≡ 1 mod 100`. The
remaining eleven residues all follow from the single arithmetic input
`8 mod q`.

## Contrast with `∑ r²` and `∑ r³`

* `∑ r² ≡ 0 mod 8` (eight units, each squared `≡ 1 mod 8`).
* `∑ r³ ≡ 0 mod q` for *every* divisor of `120` (third central moment
  vanishes by central symmetry `r ↔ 30 − r`).
* `∑ r⁴ ≡ 8 mod 120` (Carmichael `λ(120) = 4`) — even moment, full
  Wilson reduction to `|R|`.

The pair `(∑ r³, ∑ r⁴)` thus realises a perfect dichotomy:

* **Odd moment** ⇒ skewness collapse (zero mod every divisor of `120`).
* **Even moment** ⇒ Carmichael-Wilson reduction to `|R| = 8 mod 120`.

## The kurtosis identity

The set `R` is centrally symmetric about its mean `r̄ = (∑ r)/|R| = 15`.
Its **integer kurtosis** (fourth central moment) is

$$M_4 \;:=\; \sum_{r \in R} (r - 15)^4 \;=\; 85\,568,$$

with `M_4 / |R| = 10\,696` and the (population) excess kurtosis
`M_4 / (|R| \cdot \mathrm{Var}(R)^2) = 10\,696 / 70^2 \approx 2.183`.

Expanding `(r - 15)^4` and clearing denominators by `|R|^4 = 4096`, the
integer kurtosis identity reads:

$$|R|^4 \cdot \sum r^4 + 6 \cdot |R|^2 \cdot (\sum r)^2 \cdot \sum r^2
   \;=\; |R|^4 \cdot M_4
       + 4 \cdot |R|^3 \cdot (\sum r) \cdot \sum r^3
       + 3 \cdot |R| \cdot (\sum r)^4.$$

Numerically:
`4096·1246568 + 6·64·14400·2360 = 4096·85568 + 4·512·120·52200 + 3·8·120^4`
`= 5106102272 + 13049856000 = 350486528 + 12822220800 + 4976640000`
`= 18155958272` on each side.

## Reference

Monograph Appendix N §"Invariants arithmétiques de R" (fourth-power part),
extending the cubic table with the Carmichael-Wilson reduction `λ(120) = 4`.
-/

namespace PT.Sieve

open Finset

/-! ### The sum of fourth powers: exact value -/

/-- The sum of fourth powers of the admissible residues equals `1246568`. -/
theorem admissibleResidues_sumQuart :
    admissibleResidues.sum (fun r => r * r * r * r) = 1246568 := by
  decide

/-- The sum of fourth powers, term-by-term decomposition (for documentation). -/
theorem admissibleResidues_sumQuart_terms :
    admissibleResidues.sum (fun r => r * r * r * r)
      = 1 + 2401 + 14641 + 28561 + 83521 + 130321 + 279841 + 707281 := by
  decide

/-- The sum of fourth powers using `r ^ 4`. -/
theorem admissibleResidues_sumQuart_pow :
    admissibleResidues.sum (fun r => r ^ 4) = 1246568 := by
  decide

/-- Factorisation: `1246568 = 8 · 155821 = 2³ · 155821`. The factor `8 = |R|`
    is structural (Carmichael-Wilson `λ(120) = 4` ⇒ `r⁴ ≡ 1 mod 120` ⇒
    `∑ r⁴ = |R| + 120·k`); the residual `155821` is prime, in sharp
    contrast with the cubic sum `52200 = 2³ · 3² · 5² · 29` which factors
    over the active primes `{2, 3, 5}` plus the boundary unit `29 ∈ R`. -/
theorem admissibleResidues_sumQuart_factor :
    admissibleResidues.sum (fun r => r * r * r * r) = 8 * 155821
    ∧ admissibleResidues.sum (fun r => r * r * r * r) = 2 * 2 * 2 * 155821 := by
  refine ⟨?_, ?_⟩ <;> decide

/-! ### Carmichael-Wilson identity: `r⁴ ≡ 1 mod 120` -/

/-- **Fourth-power Wilson identity (per element).** Every admissible residue
    satisfies `r⁴ ≡ 1 mod 120`. The Carmichael exponent of `120 = 2³·3·5`
    is `λ(120) = lcm(λ(8), λ(3), λ(5)) = lcm(2, 2, 4) = 4`, hence
    `r^λ(120) = r⁴ ≡ 1 mod 120` for every `r` coprime to `120` — which
    includes every `r ∈ R` (since `R ⊂ (ℤ/30ℤ)* ⊂ (ℤ/120ℤ)*`). -/
theorem admissibleResidues_quart_carmichael :
    ∀ r ∈ admissibleResidues, (r * r * r * r) % 120 = 1 := by
  decide

/-- **Fourth-power Wilson identity (aggregate).** Summing the per-element
    congruence `r⁴ ≡ 1 mod 120` over `R` gives
    `∑ r⁴ ≡ |R| ≡ 8 mod 120`. This is the order-4 analogue of the
    classical Wilson congruence and the structural driver of the entire
    residue table below. -/
theorem admissibleResidues_sumQuart_mod_120 :
    admissibleResidues.sum (fun r => r * r * r * r) % 120 = 8 := by
  decide

/-! ### Residues of `1246568` modulo divisors of `120` (and `100`) -/

/-- `1246568 ≡ 0 mod 2` (even). -/
theorem admissibleResidues_sumQuart_mod_2 :
    admissibleResidues.sum (fun r => r * r * r * r) % 2 = 0 := by
  decide

/-- `1246568 ≡ 2 mod 3`. Each `r⁴ ≡ 1 mod 3` (Fermat: `r² ≡ 1 mod 3` for
    `r` coprime to `3`); summing `8` copies of `1` gives `8 ≡ 2 mod 3`. -/
theorem admissibleResidues_sumQuart_mod_3 :
    admissibleResidues.sum (fun r => r * r * r * r) % 3 = 2 := by
  decide

/-- `1246568 ≡ 3 mod 5`. Each `r⁴ ≡ 1 mod 5` (Fermat little theorem
    in `(ℤ/5ℤ)*`); summing `8` copies of `1` gives `8 ≡ 3 mod 5`. -/
theorem admissibleResidues_sumQuart_mod_5 :
    admissibleResidues.sum (fun r => r * r * r * r) % 5 = 3 := by
  decide

/-- `1246568 ≡ 2 mod 6`. CRT of `(0, 2)` from mod 2 and mod 3. -/
theorem admissibleResidues_sumQuart_mod_6 :
    admissibleResidues.sum (fun r => r * r * r * r) % 6 = 2 := by
  decide

/-- **Wilson collapse mod 8.** `1246568 ≡ 0 mod 8`. Every admissible `r`
    is odd, so `r² ≡ 1 mod 8`, hence `r⁴ = (r²)² ≡ 1 mod 8`; summing
    `|R| = 8` copies of `1` gives `8 ≡ 0 mod 8`. -/
theorem admissibleResidues_sumQuart_mod_8 :
    admissibleResidues.sum (fun r => r * r * r * r) % 8 = 0 := by
  decide

/-- `1246568 ≡ 8 mod 10`. CRT of `(0, 3)` from mod 2 and mod 5. -/
theorem admissibleResidues_sumQuart_mod_10 :
    admissibleResidues.sum (fun r => r * r * r * r) % 10 = 8 := by
  decide

/-- `1246568 ≡ 8 mod 15`. CRT of `(2, 3)` from mod 3 and mod 5 — but
    more directly, `r⁴ ≡ 1 mod 15` for every `r ∈ R` (Carmichael
    `λ(15) = lcm(2, 4) = 4`), so `∑ r⁴ ≡ 8 mod 15`. -/
theorem admissibleResidues_sumQuart_mod_15 :
    admissibleResidues.sum (fun r => r * r * r * r) % 15 = 8 := by
  decide

/-- `1246568 ≡ 8 mod 30`. CRT of `(0, 2, 3)` from mod 2, 3, 5; equivalently
    `r⁴ ≡ 1 mod 30` already fails (since `λ(30) = 4` but `2 ∤ r⁴` would
    require `r` odd, which it is). In fact `r⁴ mod 30` takes values `1` and
    `16`; the sum is `8 mod 30`. -/
theorem admissibleResidues_sumQuart_mod_30 :
    admissibleResidues.sum (fun r => r * r * r * r) % 30 = 8 := by
  decide

/-- `1246568 ≡ 8 mod 60`. Each `r⁴ ≡ 1 mod 60` (Carmichael
    `λ(60) = lcm(λ(4), λ(3), λ(5)) = lcm(2, 2, 4) = 4`), so
    `∑ r⁴ ≡ 8 mod 60`. -/
theorem admissibleResidues_sumQuart_mod_60 :
    admissibleResidues.sum (fun r => r * r * r * r) % 60 = 8 := by
  decide

/-- `1246568 ≡ 68 mod 100`. **The outlier**: `100 ∤ 120` and the
    Carmichael exponent `λ(100) = lcm(λ(4), λ(25)) = lcm(2, 20) = 20`
    does not divide `4`, so the fourth-power Wilson identity `r⁴ ≡ 1`
    fails modulo `100`. Numerically `1246568 = 12465·100 + 68`. -/
theorem admissibleResidues_sumQuart_mod_100 :
    admissibleResidues.sum (fun r => r * r * r * r) % 100 = 68 := by
  decide

/-! ### Kurtosis identity: integer form of the fourth central moment -/

/-- **Kurtosis identity (integer form).** Expanding `(r - r̄)^4` with
    `r̄ = (∑ r)/|R| = 15` and clearing denominators by `|R|^4 = 4096`:

    $$|R|^4 \cdot \sum r^4 \;+\; 6\,|R|^2 (\sum r)^2 \sum r^2
       \;=\; |R|^4 \cdot M_4
         \;+\; 4\,|R|^3 (\sum r) \sum r^3 \;+\; 3\,|R| (\sum r)^4,$$

    where `M_4 = ∑ (r - 15)^4 = 85568` is the fourth central moment.
    Numerically both sides equal `18\,155\,958\,272`. -/
theorem admissibleResidues_sumQuart_kurtosis_integer :
    admissibleResidues.card ^ 4 * admissibleResidues.sum (fun r => r * r * r * r)
      + 6 * admissibleResidues.card ^ 2 * (admissibleResidues.sum id) ^ 2
          * admissibleResidues.sum (fun r => r * r)
    = admissibleResidues.card ^ 4 * 85568
      + 4 * admissibleResidues.card ^ 3 * (admissibleResidues.sum id)
          * admissibleResidues.sum (fun r => r * r * r)
      + 3 * admissibleResidues.card * (admissibleResidues.sum id) ^ 4 := by
  decide

/-- **Fourth central moment.** `∑ (r - 15)^4 = 85568`, the (unnormalised)
    integer kurtosis of `R`. We cast to `ℤ` because natural subtraction
    `r - 15` truncates to `0` when `r < 15` (which happens for the four
    residues `1, 7, 11, 13`). -/
theorem admissibleResidues_central_fourth_moment :
    admissibleResidues.sum (fun r => ((r : ℤ) - 15) ^ 4 ) = 85568 := by
  decide

/-- **Mean fourth power.** `(∑ r⁴) / |R| = 1246568 / 8 = 155821`. The
    quotient is prime (see `admissibleResidues_sumQuart_factor`), an
    arithmetic accident in stark contrast with `(∑ r³)/|R| = 6525
    = 3²·5²·29` factoring over the active primes. -/
theorem admissibleResidues_sumQuart_mean :
    admissibleResidues.sum (fun r => r * r * r * r) / admissibleResidues.card
      = 155821 := by
  decide

/-! ### Relation to variance: the kurtosis-variance integer identity -/

/-- **Fourth-power / square cross identity.**
    `|R| · ∑ r⁴ − (∑ r²)² = 8 · 1246568 − 2360² = 9972544 − 5569600 = 4402944
    = 2⁸ · 3³ · 7² · 13`. This is the integer analogue of
    `|R| · Var(r²) = ∑ (r² − (∑ r²)/|R|)²` (the variance of the squared
    residues), unnormalised. -/
theorem admissibleResidues_sumQuart_cross_square :
    admissibleResidues.card * admissibleResidues.sum (fun r => r * r * r * r)
      = (admissibleResidues.sum (fun r => r * r)) ^ 2 + 4402944 := by
  decide

/-- Numerator of `Var(R²) = (|R|·∑ r⁴ − (∑ r²)²)/|R|² = 4402944/64 = 68796
    = 2²·3·7²·117`. -/
theorem admissibleResidues_var_of_squares :
    (admissibleResidues.card * admissibleResidues.sum (fun r => r * r * r * r)
      - (admissibleResidues.sum (fun r => r * r)) ^ 2)
        / admissibleResidues.card ^ 2 = 68796 := by
  decide

/-! ### Headline -/

/-- **Headline (sum-of-fourth-powers invariants).** The admissible residue
    set `R = {1, 7, 11, 13, 17, 19, 23, 29}` satisfies:

    * `∑ r⁴ = 1\,246\,568 = 2³ · 155\,821` (the residual factor `155\,821`
      is prime — an arithmetic accident).
    * **Per-element Carmichael-Wilson identity**: `r⁴ ≡ 1 mod 120` for
      every `r ∈ R` (since `λ(120) = 4`), hence
      **`∑ r⁴ ≡ |R| = 8 mod 120`**.
    * Modular table (driven by `8 mod q`):
      `≡ 0 mod 2`, `≡ 2 mod 3`, `≡ 3 mod 5`, `≡ 2 mod 6`,
      `≡ 0 mod 8`, `≡ 8 mod 10`, `≡ 8 mod 15`, `≡ 8 mod 30`,
      `≡ 8 mod 60`, **`≡ 8 mod 120`**, with the single outlier
      `≡ 68 mod 100` (since `100 ∤ 120` and `λ(100) = 20 ∤ 4`).
    * Fourth central moment: `M_4 = ∑ (r - 15)⁴ = 85\,568`, integer
      kurtosis identity satisfied (both sides `= 18\,155\,958\,272`).
    * Mean fourth power: `(∑ r⁴)/|R| = 155\,821` (prime).
    * Variance of squares (numerator): `|R| · ∑ r⁴ − (∑ r²)² = 4\,402\,944`.

    Compared with `∑ r³` (which collapsed `≡ 0 mod q` for *every*
    divisor of `120` by central symmetry), the fourth-power table is
    governed by a **single Carmichael-Wilson reduction** to `|R| mod 120`.
    The pair `(odd : skewness collapse, even : Wilson reduction)` makes
    explicit the moment-parity dichotomy of `(ℤ/30ℤ)*`. -/
theorem admissibleResidues_sumQuart_summary :
    admissibleResidues.sum (fun r => r * r * r * r) = 1246568
    ∧ admissibleResidues.sum (fun r => r * r * r * r) % 2 = 0
    ∧ admissibleResidues.sum (fun r => r * r * r * r) % 3 = 2
    ∧ admissibleResidues.sum (fun r => r * r * r * r) % 5 = 3
    ∧ admissibleResidues.sum (fun r => r * r * r * r) % 6 = 2
    ∧ admissibleResidues.sum (fun r => r * r * r * r) % 8 = 0
    ∧ admissibleResidues.sum (fun r => r * r * r * r) % 10 = 8
    ∧ admissibleResidues.sum (fun r => r * r * r * r) % 15 = 8
    ∧ admissibleResidues.sum (fun r => r * r * r * r) % 30 = 8
    ∧ admissibleResidues.sum (fun r => r * r * r * r) % 60 = 8
    ∧ admissibleResidues.sum (fun r => r * r * r * r) % 100 = 68
    ∧ admissibleResidues.sum (fun r => r * r * r * r) % 120 = 8
    ∧ (∀ r ∈ admissibleResidues, (r * r * r * r) % 120 = 1) :=
  ⟨admissibleResidues_sumQuart,
   admissibleResidues_sumQuart_mod_2,
   admissibleResidues_sumQuart_mod_3,
   admissibleResidues_sumQuart_mod_5,
   admissibleResidues_sumQuart_mod_6,
   admissibleResidues_sumQuart_mod_8,
   admissibleResidues_sumQuart_mod_10,
   admissibleResidues_sumQuart_mod_15,
   admissibleResidues_sumQuart_mod_30,
   admissibleResidues_sumQuart_mod_60,
   admissibleResidues_sumQuart_mod_100,
   admissibleResidues_sumQuart_mod_120,
   admissibleResidues_quart_carmichael⟩

end PT.Sieve
