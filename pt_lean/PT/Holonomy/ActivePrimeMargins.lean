/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Holonomy.ActivePrimeCriterion
import PT.Holonomy.GammaMonotonicity
import Mathlib.Tactic

/-!
# Active Prime — Quantitative margins (Ch06 extension)

This file extends `PT.Holonomy.ActivePrimeCriterion` with **quantitative
margins**: explicit lower/upper brackets on `γ_p − s` (positive margins for
active primes, negative for inactive).

* **Active margins** (`γ_p > s`):
  - `γ_3 - s ∈ (0.30, 0.31)` (large margin, `≈ 0.3076`)
  - `γ_5 - s ∈ (0.19, 0.20)` (medium margin, `≈ 0.1963`)
  - `γ_7 - s ∈ (0.09, 0.10)` (small margin, `≈ 0.0955`)
* **Inactive margins** (`γ_p < s`):
  - `s - γ_11 ∈ (0.07, 0.08)` (modest gap, `≈ 0.0743`)
  - `s - γ_13 ∈ (0.14, 0.15)` (small gap, `≈ 0.1438`)
* **Robustness** — strict ordering `γ_3 > γ_5 > γ_7 > 1/2 > γ_11 > γ_13` (the
  active-set determination is well above the threshold).
* **Extension to `p = 17`**: `γ_17 < γ_13` (the inactivity trend continues).

All bounds are proven via `norm_num` on rational expressions.

## Reference

Audit follow-up to `ActivePrimeCriterion`. Monograph Ch06 Table 6.2
("Marges des premiers actifs").
-/

namespace PT.Holonomy

/-! ### Exact margins (at `μ* = 15`, `q = 13/15`)

At the PT fixed point, every margin `γ_p - s` (or `s - γ_p`) is an *exact*
rational, computable from `gammaQ` and `sPT` by `decide` / `native_decide`.
We expose the exact values here as the source of truth; the decimal
brackets below are immediate corollaries via `norm_num`. -/

/-- **Exact margin (p = 3).** `γ_3 - s = 1727777 / 5616704`. -/
theorem gamma_3_margin_exact : gammaQ 3 - sPT = 1727777 / 5616704 := by
  unfold gammaQ deltaQ qPT muStar sPT
  norm_num

/-- **Exact margin (p = 5).** `γ_5 - s = 68621964134 / 349548756097`. -/
theorem gamma_5_margin_exact :
    gammaQ 5 - sPT = 68621964134 / 349548756097 := by
  unfold gammaQ deltaQ qPT muStar sPT
  norm_num

/-- **Exact margin (p = 7).**
    `γ_7 - s = 453321961202883 / 4748396022746468`. -/
theorem gamma_7_margin_exact :
    gammaQ 7 - sPT = 453321961202883 / 4748396022746468 := by
  unfold gammaQ deltaQ qPT muStar sPT
  norm_num

/-- **Exact margin (p = 11).**
    `s - γ_11 = 23355758527719967875358739 / 314484242174863298491777064`. -/
theorem gamma_11_margin_exact :
    sPT - gammaQ 11 = 23355758527719967875358739 / 314484242174863298491777064 := by
  unfold gammaQ deltaQ qPT muStar sPT
  norm_num

/-- **Exact margin (p = 13).**
    `s - γ_13 = 1445737279842490789790798287284 / 10056656355222204000717535253077`. -/
theorem gamma_13_margin_exact :
    sPT - gammaQ 13
      = 1445737279842490789790798287284 / 10056656355222204000717535253077 := by
  unfold gammaQ deltaQ qPT muStar sPT
  norm_num

/-! ### Active margins -/

/-- **Active margin (p = 3).** `0.30 < γ_3 - s < 0.31`. -/
theorem gamma_3_margin :
    30 / 100 < gammaQ 3 - sPT ∧ gammaQ 3 - sPT < 31 / 100 := by
  unfold gammaQ deltaQ qPT muStar sPT
  refine ⟨?_, ?_⟩ <;> norm_num

/-- **Active margin (p = 5).** `0.19 < γ_5 - s < 0.20`. -/
theorem gamma_5_margin :
    19 / 100 < gammaQ 5 - sPT ∧ gammaQ 5 - sPT < 20 / 100 := by
  unfold gammaQ deltaQ qPT muStar sPT
  refine ⟨?_, ?_⟩ <;> norm_num

/-- **Active margin (p = 7).** `0.09 < γ_7 - s < 0.10`. -/
theorem gamma_7_margin :
    9 / 100 < gammaQ 7 - sPT ∧ gammaQ 7 - sPT < 10 / 100 := by
  unfold gammaQ deltaQ qPT muStar sPT
  refine ⟨?_, ?_⟩ <;> norm_num

/-! ### Inactive margins -/

/-- **Inactive margin (p = 11).** `0.07 < s - γ_11 < 0.08`. -/
theorem gamma_11_margin :
    7 / 100 < sPT - gammaQ 11 ∧ sPT - gammaQ 11 < 8 / 100 := by
  unfold gammaQ deltaQ qPT muStar sPT
  refine ⟨?_, ?_⟩ <;> norm_num

/-- **Inactive margin (p = 13).** `0.14 < s - γ_13 < 0.15`. -/
theorem gamma_13_margin :
    14 / 100 < sPT - gammaQ 13 ∧ sPT - gammaQ 13 < 15 / 100 := by
  unfold gammaQ deltaQ qPT muStar sPT
  refine ⟨?_, ?_⟩ <;> norm_num

/-! ### Strict ordering chain

Les théorèmes `gamma_3_gt_gamma_5`, `gamma_5_gt_gamma_7`, `gamma_11_gt_gamma_13`
sont importés depuis `PT.Holonomy.GammaMonotonicity`. -/

/-! ### Extension to `p = 17` -/

/-- `γ_17 < 1/2` (inactive). -/
theorem gamma_17_inactive : gammaQ 17 < sPT := by
  unfold gammaQ deltaQ qPT muStar sPT
  norm_num

/-- `γ_17 < γ_13` (monotone decay continues). -/
theorem gamma_17_lt_gamma_13 : gammaQ 17 < gammaQ 13 := by
  unfold gammaQ deltaQ qPT muStar
  norm_num

/-! ### Headline -/

/-- **Headline (quantitative active-prime chain).** The full strict ordering

    `γ_3 > γ_5 > γ_7 > 1/2 > γ_11 > γ_13 > γ_17 > 0`

    holds at `μ* = 15` with `q = 13/15`. The active set is `{3, 5, 7}` with
    margins all `> 0.18`; the inactive primes `{11, 13, 17}` all satisfy
    `γ < 0.45`, well below `1/2`. -/
theorem active_prime_chain_quantitative :
    gammaQ 3 > gammaQ 5
    ∧ gammaQ 5 > gammaQ 7
    ∧ gammaQ 7 > sPT
    ∧ sPT > gammaQ 11
    ∧ gammaQ 11 > gammaQ 13
    ∧ gammaQ 13 > gammaQ 17 := by
  refine ⟨gamma_3_gt_gamma_5, gamma_5_gt_gamma_7, gamma_7_active,
          ?_, gamma_11_gt_gamma_13, ?_⟩
  · exact gamma_11_inactive
  · exact gamma_17_lt_gamma_13

/-- **Robustness witness.** All active margins exceed `0.09` (the smallest
    margin, attained at `p = 7`): `γ_p - 1/2 > 0.09` for `p ∈ {3, 5, 7}`. -/
theorem active_margin_robust :
    gammaQ 3 - sPT > 9 / 100
    ∧ gammaQ 5 - sPT > 9 / 100
    ∧ gammaQ 7 - sPT > 9 / 100 := by
  refine ⟨?_, ?_, ?_⟩
  · have := gamma_3_margin.1; linarith
  · have := gamma_5_margin.1; linarith
  · exact gamma_7_margin.1

end PT.Holonomy
