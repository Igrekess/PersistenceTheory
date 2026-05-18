/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Sieve.Bimodality
import PT.Sieve.AdmissibleResiduesArithmetic
import PT.Sieve.AdmissibleResiduesSquareSumMod
import Mathlib.Tactic

/-!
# Admissible residues — Sum of cubes modulo various `q` (App N extension)

This file collects the **sum of cubes invariants** of the admissible
residue set `R = {1, 7, 11, 13, 17, 19, 23, 29} = (ℤ/30ℤ)*`:

$$\sum_{r \in R} r^3 \;=\; 1 + 343 + 1331 + 2197 + 4913 + 6859 + 12167 + 24389
                   \;=\; 52200.$$

Companion to:
* `AdmissibleResiduesArithmetic.lean` — additive `∑ r = 120`.
* `AdmissibleResiduesSquareSumMod.lean` — `∑ r² = 2360` and its residues.

## Main results

* `admissibleResidues_sumCube`: `∑ r³ = 52200`.
* `admissibleResidues_sumCube_factor`: `52200 = 435 · 120 = 2³ · 3² · 5² · 29`.
* `admissibleResidues_sumCube_mod_*`: closed-form residue of `52200` modulo
  every divisor of `120 = 2³ · 3 · 5` and the auxiliary moduli `100`.
* `admissibleResidues_sumCube_skew_zero`: the **third central moment** identity
  `∑ r³ − 3·(∑ r²)·(∑ r)/|R| + 2·(∑ r)³/|R|² = 0`
  (consequence of the central symmetry `R = 30 − R`).
* `admissibleResidues_sumCube_summary`: headline table of residues.

## Patterns (Wilson collapse at order 3)

The residue `52200 mod q` exhibits the following **uniform collapse**
across **every** divisor of `120`:

| `q`  |  `52200 mod q`  |  Note                                                |
|------|-----------------|------------------------------------------------------|
|  2   |  0              |  pair (52200 = 26100·2).                             |
|  3   |  0              |  digit sum 5+2+2+0+0 = 9 ≡ 0 mod 3.                  |
|  5   |  0              |  last digit 0.                                       |
|  6   |  0              |  CRT(2,3).                                           |
|  8   |  0              |  52200 = 6525·8.                                     |
|  10  |  0              |  last digit 0.                                       |
|  15  |  0              |  CRT(3,5).                                           |
|  30  |  0              |  52200 = 1740·30.                                    |
|  60  |  0              |  52200 = 870·60.                                     |
| 100  |  0              |  52200 = 522·100.                                    |
| 120  |  0              |  52200 = 435·120.                                    |

This is the **cubic Wilson collapse**: at order 3, *every* divisor `q` of
`120 = 2³·3·5` divides `∑ r³`. The structural reason is the central
symmetry `r ↔ 30 − r` of `(ℤ/30ℤ)*` combined with the affine identities
`∑ r = 120` and `∑ r² = 2360`. The collapse is *stronger* than for
`∑ r²` (whose residues are non-zero modulo `3, 6, 15`), reflecting the
fact that the odd cubic moment of a centrally symmetric set vanishes.

## The skewness identity

The set `R` is centrally symmetric about its mean `(∑ r)/|R| = 15`:
the involution `r ↦ 30 − r` is a bijection of `R`. The third central
moment therefore vanishes:

$$\sum_{r \in R} (r - \bar r)^3 \;=\; 0,$$

equivalently in raw moments:

$$\sum r^3 \;-\; 3\,\bar r\,\sum r^2 \;+\; 2\,|R|\,\bar r^{\,3} \;=\; 0,$$

which after clearing denominators by `|R|² = 64` becomes:

$$|R|^2 \sum r^3 \;-\; 3\,(\sum r)(\sum r^2)\,|R| \;+\; 2\,(\sum r)^3 \;=\; 0.$$

Numerically: `64·52200 − 3·120·2360·8 + 2·120³ = 3340800 − 6796800 + 3456000 = 0`.

## Reference

Monograph Appendix N §"Invariants arithmétiques de R" (cube-sum part).
-/

namespace PT.Sieve

open Finset

/-! ### The sum of cubes: exact value -/

/-- The sum of cubes of the admissible residues equals `52200`. -/
theorem admissibleResidues_sumCube :
    admissibleResidues.sum (fun r => r * r * r) = 52200 := by
  decide

/-- The sum of cubes, term-by-term decomposition (for documentation). -/
theorem admissibleResidues_sumCube_terms :
    admissibleResidues.sum (fun r => r * r * r)
      = 1 + 343 + 1331 + 2197 + 4913 + 6859 + 12167 + 24389 := by
  decide

/-- The sum of cubes using `r ^ 3`. -/
theorem admissibleResidues_sumCube_pow :
    admissibleResidues.sum (fun r => r ^ 3) = 52200 := by
  decide

/-- Prime factorisation: `52200 = 2³ · 3² · 5² · 29`. Note the appearance
    of `29 ∈ R` as the only prime factor `> 5` — a numerical coincidence
    reflecting `52200 = 1800 · 29 = (∑ r)² / 8 · 29`. -/
theorem admissibleResidues_sumCube_factor :
    admissibleResidues.sum (fun r => r * r * r) = 435 * 120
    ∧ admissibleResidues.sum (fun r => r * r * r) = 2 * 2 * 2 * 3 * 3 * 5 * 5 * 29 := by
  refine ⟨?_, ?_⟩ <;> decide

/-! ### Residues of `52200` modulo divisors of `120` (and `100`) -/

/-- `52200 ≡ 0 mod 2` (even). -/
theorem admissibleResidues_sumCube_mod_2 :
    admissibleResidues.sum (fun r => r * r * r) % 2 = 0 := by
  decide

/-- `52200 ≡ 0 mod 3`. -/
theorem admissibleResidues_sumCube_mod_3 :
    admissibleResidues.sum (fun r => r * r * r) % 3 = 0 := by
  decide

/-- `52200 ≡ 0 mod 5`. -/
theorem admissibleResidues_sumCube_mod_5 :
    admissibleResidues.sum (fun r => r * r * r) % 5 = 0 := by
  decide

/-- `52200 ≡ 0 mod 6`. -/
theorem admissibleResidues_sumCube_mod_6 :
    admissibleResidues.sum (fun r => r * r * r) % 6 = 0 := by
  decide

/-- **Cubic Wilson collapse mod 8.** `52200 ≡ 0 mod 8`: every admissible
    `r` is odd, so `r³ ≡ r mod 8`; the order-1 collapse `∑ r ≡ 0 mod 8`
    (from `∑ r = 120 = 15·8`) then transports to `∑ r³ ≡ 0 mod 8`. -/
theorem admissibleResidues_sumCube_mod_8 :
    admissibleResidues.sum (fun r => r * r * r) % 8 = 0 := by
  decide

/-- `52200 ≡ 0 mod 10`. -/
theorem admissibleResidues_sumCube_mod_10 :
    admissibleResidues.sum (fun r => r * r * r) % 10 = 0 := by
  decide

/-- `52200 ≡ 0 mod 15`. Contrast with `∑ r² ≡ 5 mod 15`: the cubic moment
    collapses where the quadratic one does not, by the central symmetry. -/
theorem admissibleResidues_sumCube_mod_15 :
    admissibleResidues.sum (fun r => r * r * r) % 15 = 0 := by
  decide

/-- `52200 ≡ 0 mod 30`. -/
theorem admissibleResidues_sumCube_mod_30 :
    admissibleResidues.sum (fun r => r * r * r) % 30 = 0 := by
  decide

/-- `52200 ≡ 0 mod 60`. -/
theorem admissibleResidues_sumCube_mod_60 :
    admissibleResidues.sum (fun r => r * r * r) % 60 = 0 := by
  decide

/-- `52200 ≡ 0 mod 100`. -/
theorem admissibleResidues_sumCube_mod_100 :
    admissibleResidues.sum (fun r => r * r * r) % 100 = 0 := by
  decide

/-- **Total Wilson collapse mod 120.** `52200 ≡ 0 mod 120`, equivalently
    `∑ r³ = 435 · 120 = 435 · (∑ r)`. The order-3 moment is an integer
    multiple of the order-1 moment, with ratio `435 = 3 · 5 · 29`. -/
theorem admissibleResidues_sumCube_mod_120 :
    admissibleResidues.sum (fun r => r * r * r) % 120 = 0 := by
  decide

/-! ### Skewness identity: third central moment vanishes -/

/-- **Skewness identity (integer form).** Multiplying through by `|R|² = 64`
    to keep everything in `ℕ`:

    $$|R|^2 \cdot \sum r^3 \;+\; 2 (\sum r)^3
       \;=\; 3 \cdot |R| \cdot (\sum r) (\sum r^2).$$

    Numerically: `64 · 52200 + 2 · 120³ = 3 · 8 · 120 · 2360`,
    i.e. `3340800 + 3456000 = 6796800`. This expresses the vanishing of
    the third central moment `∑(r - r̄)³ = 0`, a direct consequence of
    the central symmetry `r ↔ 30 - r` of `R`. -/
theorem admissibleResidues_sumCube_skew_zero :
    admissibleResidues.card ^ 2 * admissibleResidues.sum (fun r => r * r * r)
      + 2 * (admissibleResidues.sum id) ^ 3
      = 3 * admissibleResidues.card * (admissibleResidues.sum id)
          * (admissibleResidues.sum (fun r => r * r)) := by
  decide

/-- **Mean cube.** The arithmetic mean of `r³ over R` equals `6525 = 3² · 5² · 29`.
    Compare with the cube of the arithmetic mean `15³ = 3375`: the difference
    `6525 − 3375 = 3150 = 3·1050` is `(3·Var(R))·r̄ + Skew(R)`, with the
    skewness component vanishing by the symmetry identity above. -/
theorem admissibleResidues_sumCube_mean :
    admissibleResidues.sum (fun r => r * r * r) / admissibleResidues.card = 6525 := by
  decide

/-! ### Headline -/

/-- **Headline (sum-of-cubes invariants).** The admissible residue set
    `R = {1, 7, 11, 13, 17, 19, 23, 29}` satisfies:

    * `∑ r³ = 52200 = 435 · 120 = 2³ · 3² · 5² · 29`.
    * **Total Wilson collapse**: `52200 ≡ 0 mod q` for every
      `q ∈ {2, 3, 5, 6, 8, 10, 15, 30, 60, 100, 120}` — every divisor
      of `120` divides `∑ r³` (and so does `100`).
    * **Third central moment vanishes**: by the central symmetry
      `r ↔ 30 − r`, the integer identity
      `|R|²·∑ r³ + 2·(∑ r)³ = 3·|R|·(∑ r)·(∑ r²)` holds, i.e.
      `Σ(r − r̄)³ = 0` and the population skewness of `R` is zero.
    * Mean cube: `(∑ r³)/|R| = 6525`.

    Compared with the quadratic table (`AdmissibleResiduesSquareSumMod`),
    the cubic table is "fully collapsed": where `∑ r² mod q` left non-zero
    residues `{2, 2, 5, 20, 20}` at `q ∈ {3, 6, 15, 30, 60}`, the cubic
    moment is identically `0` modulo every divisor of `120`. This is the
    Wilson collapse upgraded by the antisymmetry of `r ↦ r³` about the
    symmetry centre `r = 15`. -/
theorem admissibleResidues_sumCube_summary :
    admissibleResidues.sum (fun r => r * r * r) = 52200
    ∧ admissibleResidues.sum (fun r => r * r * r) % 2 = 0
    ∧ admissibleResidues.sum (fun r => r * r * r) % 3 = 0
    ∧ admissibleResidues.sum (fun r => r * r * r) % 5 = 0
    ∧ admissibleResidues.sum (fun r => r * r * r) % 6 = 0
    ∧ admissibleResidues.sum (fun r => r * r * r) % 8 = 0
    ∧ admissibleResidues.sum (fun r => r * r * r) % 10 = 0
    ∧ admissibleResidues.sum (fun r => r * r * r) % 15 = 0
    ∧ admissibleResidues.sum (fun r => r * r * r) % 30 = 0
    ∧ admissibleResidues.sum (fun r => r * r * r) % 60 = 0
    ∧ admissibleResidues.sum (fun r => r * r * r) % 100 = 0
    ∧ admissibleResidues.sum (fun r => r * r * r) % 120 = 0
    ∧ admissibleResidues.card ^ 2 * admissibleResidues.sum (fun r => r * r * r)
        + 2 * (admissibleResidues.sum id) ^ 3
        = 3 * admissibleResidues.card * (admissibleResidues.sum id)
            * (admissibleResidues.sum (fun r => r * r)) :=
  ⟨admissibleResidues_sumCube,
   admissibleResidues_sumCube_mod_2,
   admissibleResidues_sumCube_mod_3,
   admissibleResidues_sumCube_mod_5,
   admissibleResidues_sumCube_mod_6,
   admissibleResidues_sumCube_mod_8,
   admissibleResidues_sumCube_mod_10,
   admissibleResidues_sumCube_mod_15,
   admissibleResidues_sumCube_mod_30,
   admissibleResidues_sumCube_mod_60,
   admissibleResidues_sumCube_mod_100,
   admissibleResidues_sumCube_mod_120,
   admissibleResidues_sumCube_skew_zero⟩

end PT.Sieve
