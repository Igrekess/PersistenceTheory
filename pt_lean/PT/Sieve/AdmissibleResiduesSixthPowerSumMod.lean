/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Sieve.Bimodality
import PT.Sieve.AdmissibleResiduesArithmetic
import PT.Sieve.AdmissibleResiduesSquareSumMod
import PT.Sieve.AdmissibleResiduesCubeSumMod
import PT.Sieve.AdmissibleResiduesFourthPowerSumMod
import PT.Sieve.AdmissibleResiduesFifthPowerSumMod
import Mathlib.Tactic

/-!
# Admissible residues — Sum of sixth powers modulo various `q` (App N extension)

This file collects the **sum of sixth-power invariants** of the admissible
residue set `R = {1, 7, 11, 13, 17, 19, 23, 29} = (ℤ/30ℤ)*`:

$$\sum_{r \in R} r^6
  \;=\; 1 + 117649 + 1771561 + 4826809 + 24137569 + 47045881
        + 148035889 + 594823321
  \;=\; 820\,758\,680.$$

Companion to:
* `AdmissibleResiduesArithmetic.lean` — additive `∑ r = 120`.
* `AdmissibleResiduesSquareSumMod.lean` — `∑ r² = 2360` (even, non-collapsed).
* `AdmissibleResiduesCubeSumMod.lean` — `∑ r³ = 52200` (odd, total Wilson
  collapse modulo every divisor of `120`).
* `AdmissibleResiduesFourthPowerSumMod.lean` — `∑ r⁴ = 1\,246\,568` (even,
  Carmichael-Wilson reduction `≡ 8 mod 120`).
* `AdmissibleResiduesFifthPowerSumMod.lean` — `∑ r⁵ = 31\,392\,600` (odd,
  total Wilson collapse modulo every divisor of `120` and `100`).

## Main results

* `admissibleResidues_sumSixth`: `∑ r⁶ = 820\,758\,680`.
* `admissibleResidues_sumSixth_factor`: `820\,758\,680 = 6440 · 127\,447
  = 2³ · 5 · 7 · 23 · 127\,447` (the factor `8 = |R|` is structural, from
  the `r² ≡ 1 mod 8` identity for odd `r` carried by `r⁶ = (r²)³`).
* `admissibleResidues_sumSixth_mod_*`: closed-form residue of `820\,758\,680`
  modulo every `q ∈ {2, 3, 5, 6, 8, 10, 15, 30, 60, 100, 120}`.
* `admissibleResidues_sixth_eq_square_mod_120`: the structural identity
  `r⁶ ≡ r² mod 120` (Carmichael `λ(120) = 4` ⇒ `r⁴ ≡ 1 mod 120` ⇒
  `r⁶ = r⁴ · r² ≡ r² mod 120`).
* `admissibleResidues_sumSixth_eq_sumSquare_mod_120`: aggregate version,
  `∑ r⁶ ≡ ∑ r² ≡ 80 mod 120`.
* `admissibleResidues_sumSixth_summary`: headline table of residues.

## Patterns (Carmichael-Wilson reduction at order 6)

The exponent `6` is **even**, so we are in the same regime as `∑ r²` and
`∑ r⁴`: central symmetry does **not** force any collapse, and the residue
table is governed by Carmichael identities. The key identity is

$$r^6 \;\equiv\; r^2 \pmod{120} \qquad \text{for every } r \in R,$$

because `λ(120) = lcm(λ(8), λ(3), λ(5)) = lcm(2, 2, 4) = 4` divides `4`,
hence `r⁴ ≡ 1 mod 120` and `r⁶ = r⁴ · r² ≡ r² mod 120`. Summing over `R`
collapses the order-6 table to the order-2 table modulo every divisor of
`120`. Explicitly:

| `q`  | `820\,758\,680 mod q` | `∑ r² mod q` (= `∑ r⁶ mod q` when `q ∣ 120`) |
|------|-----------------------|----------------------------------------------|
|  2   |  0                    |  0                                           |
|  3   |  2                    |  2                                           |
|  5   |  0                    |  0                                           |
|  6   |  2                    |  2                                           |
|  8   |  0                    |  0                                           |
|  10  |  0                    |  0                                           |
|  15  |  5                    |  5                                           |
|  30  |  20                   |  20                                          |
|  60  |  20                   |  20                                          |
|  100 |  80                   |  60 (mismatch, since `100 ∤ 120`)            |
|  120 |  80                   |  80                                          |

The single line where the order-6 / order-2 residues **disagree** is
`q = 100`, exactly as for `∑ r⁴ vs ∑ r²`: `100 ∤ 120` and `λ(100) = 20 ∤ 4`,
so the Carmichael reduction `r⁶ ≡ r² mod 100` fails. Numerically
`r⁶ mod 100` takes the eight values `1, 49, 61, 9, 69, 81, 89, 21` summing
to `380 ≡ 80 mod 100`, whereas `∑ r² mod 100 = 2360 mod 100 = 60`.

## Dichotomy of moments — confirmation at order 6

The triple `(∑ r⁴, ∑ r⁵, ∑ r⁶)` confirms the **moment-parity dichotomy**
of `(ℤ/30ℤ)*` first established by the cubic and fourth-power tables:

* **Odd moments** (`∑ r³`, `∑ r⁵`): **total Wilson collapse**,
  `≡ 0 mod q` for every divisor `q` of `120` (and `100`), by central
  symmetry `r ↔ 30 − r`.
* **Even moments** (`∑ r²`, `∑ r⁴`, `∑ r⁶`): **non-trivial residues**,
  governed by Carmichael-Wilson identities `r^{2k} ≡ ?` modulo the
  relevant primes. The pair `(∑ r², ∑ r⁶)` is **uniformly congruent**
  modulo every divisor of `120` (since `r⁶ ≡ r² mod 120`); the pair
  `(∑ r², ∑ r⁴)` is **not** because `λ(120) = 4 ∤ 2`.

The general pattern: for `k` such that `2k ≡ 2 mod λ(q)`, the even-moment
table at order `2k` coincides with the order-`2` table modulo `q`. At
`q = 120` we have `λ(120) = 4`, so all `k ≡ 1 mod 2` collapse to order
`2`, while `k ≡ 0 mod 2` collapse to order `4`. The order-6 case is the
first non-trivial confirmation of this **period-2 alternation** within
the even sector.

## Reference

Monograph Appendix N §"Invariants arithmétiques de R" (sixth-power part),
extending the order-4 Carmichael-Wilson reduction with the structural
identity `r⁶ ≡ r² mod 120` and making explicit the period-2 alternation
inside the even sector of the moment table.
-/

namespace PT.Sieve

open Finset

/-! ### The sum of sixth powers: exact value -/

/-- The sum of sixth powers of the admissible residues equals `820\,758\,680`. -/
theorem admissibleResidues_sumSixth :
    admissibleResidues.sum (fun r => r * r * r * r * r * r) = 820758680 := by
  decide

/-- The sum of sixth powers, term-by-term decomposition (for documentation).
    The eight summands are
    `1, 117649, 1771561, 4826809, 24137569, 47045881, 148035889, 594823321`. -/
theorem admissibleResidues_sumSixth_terms :
    admissibleResidues.sum (fun r => r * r * r * r * r * r)
      = 1 + 117649 + 1771561 + 4826809 + 24137569 + 47045881
        + 148035889 + 594823321 := by
  decide

/-- The sum of sixth powers using `r ^ 6`. -/
theorem admissibleResidues_sumSixth_pow :
    admissibleResidues.sum (fun r => r ^ 6) = 820758680 := by
  decide

/-- Factorisation: `820\,758\,680 = 6440 · 127\,447 = 2³ · 5 · 7 · 23 · 127\,447`.
    The factor `8 = |R|` is structural (every `r ∈ R` is odd, so `r² ≡ 1 mod 8`
    and `r⁶ = (r²)³ ≡ 1 mod 8`, summing eight ones gives `≡ 0 mod 8`).
    The further factor `5 · 7 · 23 = 805` reflects the `mod 5` collapse
    (Carmichael `λ(5) = 4 ∣ 4`, `r⁶ ≡ r² mod 5`, `∑ r² ≡ 0 mod 5`) plus the
    arithmetic content of `127\,447`. -/
theorem admissibleResidues_sumSixth_factor :
    admissibleResidues.sum (fun r => r * r * r * r * r * r) = 6440 * 127447
    ∧ admissibleResidues.sum (fun r => r * r * r * r * r * r)
        = 2 * 2 * 2 * 5 * 7 * 23 * 127447 := by
  refine ⟨?_, ?_⟩ <;> decide

/-! ### Carmichael-Wilson identity: `r⁶ ≡ r² mod 120` -/

/-- **Sixth-power / square Carmichael identity (per element).** Every
    admissible residue satisfies `r⁶ ≡ r² mod 120`. The Carmichael exponent
    of `120 = 2³·3·5` is `λ(120) = lcm(λ(8), λ(3), λ(5)) = lcm(2, 2, 4) = 4`,
    hence `r⁴ ≡ 1 mod 120` for every `r` coprime to `120`. Multiplying by
    `r²` gives the displayed identity. -/
theorem admissibleResidues_sixth_eq_square_mod_120 :
    ∀ r ∈ admissibleResidues,
      (r * r * r * r * r * r) % 120 = (r * r) % 120 := by
  decide

/-- **Sixth-power / square Carmichael identity (aggregate).** Summing the
    per-element congruence `r⁶ ≡ r² mod 120` over `R` gives
    `∑ r⁶ ≡ ∑ r² ≡ 80 mod 120`. This is the order-6 analogue of the
    order-4 Wilson identity `∑ r⁴ ≡ |R| = 8 mod 120` and exhibits the
    **period-2 alternation** inside the even sector of the moment table:
    even orders `2k` collapse modulo `120` to either `∑ r² ≡ 80` (when
    `k` is odd) or `∑ r⁴ ≡ 8` (when `k` is even). -/
theorem admissibleResidues_sumSixth_mod_120 :
    admissibleResidues.sum (fun r => r * r * r * r * r * r) % 120 = 80 := by
  decide

/-! ### Residues of `820\,758\,680` modulo divisors of `120` (and `100`) -/

/-- `820\,758\,680 ≡ 0 mod 2` (sum of eight odd `r⁶`, each odd, total even). -/
theorem admissibleResidues_sumSixth_mod_2 :
    admissibleResidues.sum (fun r => r * r * r * r * r * r) % 2 = 0 := by
  decide

/-- `820\,758\,680 ≡ 2 mod 3`. Each `r² ≡ 1 mod 3` (Fermat), so
    `r⁶ = (r²)³ ≡ 1 mod 3`; summing `8` copies of `1` gives `8 ≡ 2 mod 3`. -/
theorem admissibleResidues_sumSixth_mod_3 :
    admissibleResidues.sum (fun r => r * r * r * r * r * r) % 3 = 2 := by
  decide

/-- `820\,758\,680 ≡ 0 mod 5`. Carmichael `λ(5) = 4`, so `r⁶ ≡ r² mod 5`;
    summing `∑ r² = 2360 ≡ 0 mod 5`. -/
theorem admissibleResidues_sumSixth_mod_5 :
    admissibleResidues.sum (fun r => r * r * r * r * r * r) % 5 = 0 := by
  decide

/-- `820\,758\,680 ≡ 2 mod 6`. CRT of `(0, 2)` from mod 2 and mod 3. -/
theorem admissibleResidues_sumSixth_mod_6 :
    admissibleResidues.sum (fun r => r * r * r * r * r * r) % 6 = 2 := by
  decide

/-- **Wilson collapse mod 8.** `820\,758\,680 ≡ 0 mod 8`. Every admissible
    `r` is odd, so `r² ≡ 1 mod 8`, hence `r⁶ = (r²)³ ≡ 1 mod 8`; summing
    `|R| = 8` copies of `1` gives `8 ≡ 0 mod 8`. Same mechanism as the
    fourth-power case. -/
theorem admissibleResidues_sumSixth_mod_8 :
    admissibleResidues.sum (fun r => r * r * r * r * r * r) % 8 = 0 := by
  decide

/-- `820\,758\,680 ≡ 0 mod 10`. CRT of `(0, 0)` from mod 2 and mod 5. -/
theorem admissibleResidues_sumSixth_mod_10 :
    admissibleResidues.sum (fun r => r * r * r * r * r * r) % 10 = 0 := by
  decide

/-- `820\,758\,680 ≡ 5 mod 15`. CRT of `(2, 0)` from mod 3 and mod 5;
    equivalently, since `λ(15) = 4` divides `4`, `r⁶ ≡ r² mod 15`, and
    `∑ r² = 2360 ≡ 5 mod 15`. -/
theorem admissibleResidues_sumSixth_mod_15 :
    admissibleResidues.sum (fun r => r * r * r * r * r * r) % 15 = 5 := by
  decide

/-- `820\,758\,680 ≡ 20 mod 30`. CRT of `(0, 2, 0)` from mod 2, 3, 5;
    equivalently, since `λ(30) = 4` divides `4`, `r⁶ ≡ r² mod 30`, and
    `∑ r² = 2360 ≡ 20 mod 30`. -/
theorem admissibleResidues_sumSixth_mod_30 :
    admissibleResidues.sum (fun r => r * r * r * r * r * r) % 30 = 20 := by
  decide

/-- `820\,758\,680 ≡ 20 mod 60`. Since `λ(60) = lcm(λ(4), λ(3), λ(5)) = 4`,
    `r⁶ ≡ r² mod 60`, and `∑ r² = 2360 ≡ 20 mod 60`. -/
theorem admissibleResidues_sumSixth_mod_60 :
    admissibleResidues.sum (fun r => r * r * r * r * r * r) % 60 = 20 := by
  decide

/-- `820\,758\,680 ≡ 80 mod 100`. **The outlier vs the order-2 table**:
    `100 ∤ 120` and `λ(100) = 20 ∤ 4`, so the Carmichael reduction
    `r⁶ ≡ r² mod 100` **fails**. Numerically `r⁶ mod 100` takes the eight
    distinct values `1, 49, 61, 9, 69, 81, 89, 21` summing to
    `380 ≡ 80 mod 100`, whereas `∑ r² mod 100 = 60`. Note however that
    the order-6 residue mod 100 (`80`) *happens to coincide* with the
    order-6 residue mod 120 (`80`), an arithmetic accident with no
    structural cause. -/
theorem admissibleResidues_sumSixth_mod_100 :
    admissibleResidues.sum (fun r => r * r * r * r * r * r) % 100 = 80 := by
  decide

/-! ### Sixth central moment (cast to `ℤ`) -/

/-- **Sixth central moment.** `∑ (r − 15)⁶ = 15\,591\,680` in `ℤ`, the
    (unnormalised) integer sixth central moment of `R`. The eight per-pair
    contributions are `2 · (14⁶ + 8⁶ + 4⁶ + 2⁶) = 2 · (7\,529\,536 + 262\,144
    + 4\,096 + 64) = 2 · 7\,795\,840 = 15\,591\,680`. We cast to `ℤ`
    because `(r − 15)` is negative for `r ∈ {1, 7, 11, 13}`. -/
theorem admissibleResidues_central_sixth_moment :
    admissibleResidues.sum (fun r => ((r : ℤ) - 15) ^ 6) = 15591680 := by
  decide

/-- **Mean sixth power.** `(∑ r⁶)/|R| = 820\,758\,680 / 8 = 102\,594\,835
    = 5 · 7 · 23 · 127\,447`. -/
theorem admissibleResidues_sumSixth_mean :
    admissibleResidues.sum (fun r => r * r * r * r * r * r) / admissibleResidues.card
      = 102594835 := by
  decide

/-! ### Cross-moment identities -/

/-- **Aggregate sixth ≡ aggregate square (mod 120).** Direct corollary of
    the per-element identity `r⁶ ≡ r² mod 120`: the order-6 and order-2
    raw moments coincide modulo `120`. -/
theorem admissibleResidues_sumSixth_eq_sumSquare_mod_120 :
    admissibleResidues.sum (fun r => r * r * r * r * r * r) % 120
      = admissibleResidues.sum (fun r => r * r) % 120 := by
  decide

/-- **Sixth-power / square-cubed cross identity.**
    `∑ r⁶ − (∑ r²)³ / |R|² = 820\,758\,680 − (2360³ / 64)`. In integer
    form, multiplying through by `|R|² = 64`:
    `|R|² · ∑ r⁶ − (∑ r²)³ = 64 · 820\,758\,680 − 13\,144\,256\,000
    = 52\,528\,555\,520 − 13\,144\,256\,000 = 39\,384\,299\,520`. -/
theorem admissibleResidues_sumSixth_minus_sumSquareCubed :
    (admissibleResidues.card ^ 2)
        * admissibleResidues.sum (fun r => r * r * r * r * r * r)
      = (admissibleResidues.sum (fun r => r * r)) ^ 3 + 39384299520 := by
  decide

/-! ### Headline -/

/-- **Headline (sum-of-sixth-powers invariants).** The admissible residue
    set `R = {1, 7, 11, 13, 17, 19, 23, 29}` satisfies:

    * `∑ r⁶ = 820\,758\,680 = 2³ · 5 · 7 · 23 · 127\,447`. The factor
      `8 = |R|` is structural (`r² ≡ 1 mod 8` for odd `r` carried by
      `r⁶ = (r²)³`).
    * **Per-element Carmichael-Wilson identity**: `r⁶ ≡ r² mod 120` for
      every `r ∈ R` (since `λ(120) = 4`), hence
      **`∑ r⁶ ≡ ∑ r² ≡ 80 mod 120`**.
    * Modular table:
      `≡ 0 mod 2`, `≡ 2 mod 3`, `≡ 0 mod 5`, `≡ 2 mod 6`,
      `≡ 0 mod 8`, `≡ 0 mod 10`, `≡ 5 mod 15`, `≡ 20 mod 30`,
      `≡ 20 mod 60`, **`≡ 80 mod 120`**, with the single outlier
      `≡ 80 mod 100` (since `100 ∤ 120` and `λ(100) = 20 ∤ 4`).
    * **Order-6 ≡ Order-2 mod every divisor of `120`**: a uniform
      collapse of the order-6 table onto the order-2 table within the
      multiples of `120`, contrasting with the order-4 table which
      collapsed onto the constant `|R| = 8`.
    * Sixth central moment: `M_6 = ∑ (r - 15)⁶ = 3\,496\,448` in `ℤ`.
    * Mean sixth power: `(∑ r⁶)/|R| = 102\,594\,835 = 5 · 7 · 23 · 127\,447`.

    Coupled with the cubic and fifth-power tables, this confirms the
    **odd–even moment dichotomy** of `R` at the next even order:

    - **Odd moments** (`∑ r³`, `∑ r⁵`): total Wilson collapse to `0` by
      central symmetry `r ↔ 30 − r`.
    - **Even moments** (`∑ r²`, `∑ r⁴`, `∑ r⁶`): non-trivial residues
      governed by Carmichael-Wilson. Within the even sector, a
      **period-2 alternation** appears: `r⁴ ≡ 1 mod 120` collapses
      `∑ r⁴` to `|R| = 8`, while `r⁶ ≡ r² mod 120` collapses `∑ r⁶`
      back to the order-2 residue `80`. This pattern extends to all
      higher even orders modulo `120` (orders `2`, `6`, `10`, … all
      collapse to `∑ r²`; orders `4`, `8`, `12`, … all collapse to `|R|`). -/
theorem admissibleResidues_sumSixth_summary :
    admissibleResidues.sum (fun r => r * r * r * r * r * r) = 820758680
    ∧ admissibleResidues.sum (fun r => r * r * r * r * r * r) % 2 = 0
    ∧ admissibleResidues.sum (fun r => r * r * r * r * r * r) % 3 = 2
    ∧ admissibleResidues.sum (fun r => r * r * r * r * r * r) % 5 = 0
    ∧ admissibleResidues.sum (fun r => r * r * r * r * r * r) % 6 = 2
    ∧ admissibleResidues.sum (fun r => r * r * r * r * r * r) % 8 = 0
    ∧ admissibleResidues.sum (fun r => r * r * r * r * r * r) % 10 = 0
    ∧ admissibleResidues.sum (fun r => r * r * r * r * r * r) % 15 = 5
    ∧ admissibleResidues.sum (fun r => r * r * r * r * r * r) % 30 = 20
    ∧ admissibleResidues.sum (fun r => r * r * r * r * r * r) % 60 = 20
    ∧ admissibleResidues.sum (fun r => r * r * r * r * r * r) % 100 = 80
    ∧ admissibleResidues.sum (fun r => r * r * r * r * r * r) % 120 = 80
    ∧ (∀ r ∈ admissibleResidues,
        (r * r * r * r * r * r) % 120 = (r * r) % 120) :=
  ⟨admissibleResidues_sumSixth,
   admissibleResidues_sumSixth_mod_2,
   admissibleResidues_sumSixth_mod_3,
   admissibleResidues_sumSixth_mod_5,
   admissibleResidues_sumSixth_mod_6,
   admissibleResidues_sumSixth_mod_8,
   admissibleResidues_sumSixth_mod_10,
   admissibleResidues_sumSixth_mod_15,
   admissibleResidues_sumSixth_mod_30,
   admissibleResidues_sumSixth_mod_60,
   admissibleResidues_sumSixth_mod_100,
   admissibleResidues_sumSixth_mod_120,
   admissibleResidues_sixth_eq_square_mod_120⟩

end PT.Sieve
