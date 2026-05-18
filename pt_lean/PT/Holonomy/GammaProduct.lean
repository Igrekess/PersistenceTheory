/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Holonomy.ActivePrimeCriterion
import PT.Holonomy.ActivePrimeMargins
import PT.Holonomy.GammaTablesExtended
import PT.Holonomy.CouplingReconstruction
import Mathlib.Tactic

/-!
# Product of anomalous dimensions `∏ γ_p` over the active cascade `{3, 5, 7}`

**Statement (paper-level, Ch06 §"Cascade arithmétique" and Ch09 supplement).**
At the PT fixed point `μ* = 15`, `q_+ = 13/15`, the anomalous-dimension
factor at one full loop of the cascade is the **product**

$$\Gamma_{\rm active} \;:=\; \prod_{p \in \{3, 5, 7\}} \gamma_p
   \;=\; \gamma_3 \, \gamma_5 \, \gamma_7.$$

This is a *different* product from the bare coupling
`α_bare = ∏ sin²θ_p` (`PT.Holonomy.CouplingReconstruction.alphaBareQ`).
Both share the active set `{3, 5, 7}` but evaluate distinct functions of
`δ_p`:

* `γ_p   = 4 q^{p-1}(1 - δ_p) / (μ* δ_p (2 - δ_p))`
* `sin²θ_p = δ_p (2 - δ_p)`

so `γ_p · sin²θ_p = 4 q^{p-1}(1 - δ_p) / μ*`. In particular,
`Γ_active` and `α_bare` are **never equal** at `q_+ = 13/15`; their ratio
is an exact rational, given here.

This module records:

1. The **definition** `gammaProductActive` (product over `{3, 5, 7}`).
2. **Positivity** `0 < gammaProductActive`, deduced from each
   `γ_p > 1/2 > 0` (`ActivePrimeCriterion.gamma_*_active`).
3. **Tight bracket** `0.334 < Γ_active < 0.335` matching the numerical
   value `Γ_active ≈ 0.3348`.
4. **Upper bound** `gammaProductActive < 1`, a one-line consequence of
   the individual `γ_p < 1` cap (which we also prove here).
5. **Comparison to `α_bare`**: `gammaProductActive ≠ alphaBareQ`, plus
   the exact ratio bracket `45.6 < Γ_active / α_bare < 45.7`.
6. **Extension to the 5-prime tail**: `gammaProductExtended = γ_3 γ_5
   γ_7 γ_11 γ_13`, with bracket `0.050 < · < 0.051`.
7. **Strict-decreasing partial products** `P_1 > P_2 > P_3 > P_4 > P_5`,
   a direct consequence of `γ_p < 1` for every odd prime in the table.

All proofs are pure rational arithmetic — `unfold` + `norm_num`.

## Reference

* Monograph Chapter 6, §"Cascade arithmétique", `\label{thm:gamma_product}`.
* Monograph Chapter 9, §"Reconstruction des couplages" (companion to
  `α_bare`).
* `PT.Holonomy.GammaTablesExtended` for the underlying γ-table.
-/

namespace PT.Holonomy

/-! ### Definition: product over the active cascade `{3, 5, 7}` -/

/-- The PT anomalous-dimension product over the active primes `{3, 5, 7}`
    at `q_+ = 13/15`, `μ* = 15`:

    `gammaProductActive := γ_3 · γ_5 · γ_7`. -/
def gammaProductActive : ℚ := gammaQ 3 * gammaQ 5 * gammaQ 7

/-- Definitional unfolding: the active product is literally the triple
    product. -/
theorem gammaProductActive_eq_prod :
    gammaProductActive = gammaQ 3 * gammaQ 5 * gammaQ 7 := rfl

/-! ### Positivity from `γ_p > 1/2 > 0` -/

/-- `γ_3 > 0` (consequence of `γ_3 > 1/2`). -/
theorem gammaQ_3_pos : 0 < gammaQ 3 :=
  lt_trans (by norm_num : (0 : ℚ) < 1 / 2) gamma_3_active

/-- `γ_5 > 0` (consequence of `γ_5 > 1/2`). -/
theorem gammaQ_5_pos : 0 < gammaQ 5 :=
  lt_trans (by norm_num : (0 : ℚ) < 1 / 2) gamma_5_active

/-- `γ_7 > 0` (consequence of `γ_7 > 1/2`). -/
theorem gammaQ_7_pos : 0 < gammaQ 7 :=
  lt_trans (by norm_num : (0 : ℚ) < 1 / 2) gamma_7_active

/-- `γ_11 > 0` (inactive but still positive). -/
theorem gammaQ_11_pos : 0 < gammaQ 11 := by
  unfold gammaQ deltaQ qPT muStar
  norm_num

/-- `γ_13 > 0` (inactive but still positive). -/
theorem gammaQ_13_pos : 0 < gammaQ 13 := by
  unfold gammaQ deltaQ qPT muStar
  norm_num

/-- **Positivity of the active product.** `0 < γ_3 · γ_5 · γ_7`. -/
theorem gammaProductActive_pos : 0 < gammaProductActive := by
  unfold gammaProductActive
  exact mul_pos (mul_pos gammaQ_3_pos gammaQ_5_pos) gammaQ_7_pos

/-! ### Individual upper bounds `γ_p < 1` for `p ≥ 3` -/

/-- `γ_3 < 1`. -/
theorem gammaQ_3_lt_one : gammaQ 3 < 1 := by
  unfold gammaQ deltaQ qPT muStar
  norm_num

/-- `γ_5 < 1`. -/
theorem gammaQ_5_lt_one : gammaQ 5 < 1 := by
  unfold gammaQ deltaQ qPT muStar
  norm_num

/-- `γ_7 < 1`. -/
theorem gammaQ_7_lt_one : gammaQ 7 < 1 := by
  unfold gammaQ deltaQ qPT muStar
  norm_num

/-- `γ_11 < 1`. -/
theorem gammaQ_11_lt_one : gammaQ 11 < 1 := by
  unfold gammaQ deltaQ qPT muStar
  norm_num

/-- `γ_13 < 1`. -/
theorem gammaQ_13_lt_one : gammaQ 13 < 1 := by
  unfold gammaQ deltaQ qPT muStar
  norm_num

/-! ### Tight rational bracket on the active product -/

/-- **Tight bracket for `Γ_active`.**
    `0.334 < γ_3 · γ_5 · γ_7 < 0.335`. Numerical value `≈ 0.3348`. -/
theorem gammaProductActive_bracket :
    334 / 1000 < gammaProductActive ∧ gammaProductActive < 335 / 1000 := by
  unfold gammaProductActive gammaQ deltaQ qPT muStar
  refine ⟨?_, ?_⟩ <;> norm_num

/-- **Loose upper bound.** `Γ_active < 1`.

    Stated via the bracket; a direct combinatorial proof from `γ_p < 1`
    would also work but the bracket is sharper and zero-cost here. -/
theorem gammaProductActive_lt_one : gammaProductActive < 1 :=
  lt_trans gammaProductActive_bracket.2 (by norm_num : (335 / 1000 : ℚ) < 1)

/-- **Lower bound separating from zero.** `Γ_active > 1/4`.

    Useful because it shows `Γ_active` is a *macroscopic* fraction, not
    a perturbative residue: the active cascade contributes O(1) holonomy. -/
theorem gammaProductActive_gt_quarter : gammaProductActive > 1 / 4 :=
  lt_trans (by norm_num : (1 / 4 : ℚ) < 334 / 1000) gammaProductActive_bracket.1

/-! ### Comparison with the bare coupling `α_bare` -/

/-- **`gammaProductActive ≠ alphaBareQ`.**

    The two products share the active set `{3, 5, 7}` but evaluate
    different functions of `δ_p`. The bracket on `α_bare` is
    `7.3 × 10⁻³ < α_bare < 7.4 × 10⁻³`
    (`CouplingReconstruction.alphaBareQ_bracket`), well below
    `0.334 < Γ_active`, so the two cannot coincide. -/
theorem gammaProductActive_ne_alphaBare :
    gammaProductActive ≠ alphaBareQ := by
  intro h
  have h1 : (334 : ℚ) / 1000 < gammaProductActive :=
    gammaProductActive_bracket.1
  have h2 : alphaBareQ < 74 / 10000 := alphaBareQ_bracket.2
  rw [h] at h1
  -- 334/1000 < α_bare < 74/10000, contradiction
  have : (334 : ℚ) / 1000 < 74 / 10000 := lt_trans h1 h2
  norm_num at this

/-- **Exact ratio bracket.** `45.6 < Γ_active / α_bare < 45.7`.

    The ratio is the rational
    `γ_3 γ_5 γ_7 / (sin²θ_3 sin²θ_5 sin²θ_7)`, which simplifies on each
    factor to `(4 q_+^{p-1} (1 - δ_p)) / (μ* δ_p (2 - δ_p)^2)` — but
    numerically we only need the bracket. -/
theorem gammaProductActive_over_alphaBare_bracket :
    456 / 10 < gammaProductActive / alphaBareQ
    ∧ gammaProductActive / alphaBareQ < 457 / 10 := by
  unfold gammaProductActive alphaBareQ sinSqQ gammaQ deltaQ qPT muStar
  refine ⟨?_, ?_⟩ <;> norm_num

/-! ### Extension to the 5-prime tail `{3, 5, 7, 11, 13}` -/

/-- The 5-prime extended product `γ_3 γ_5 γ_7 γ_11 γ_13`. -/
def gammaProductExtended : ℚ :=
  gammaQ 3 * gammaQ 5 * gammaQ 7 * gammaQ 11 * gammaQ 13

/-- The extended product factors through the active product. -/
theorem gammaProductExtended_eq :
    gammaProductExtended = gammaProductActive * gammaQ 11 * gammaQ 13 := by
  unfold gammaProductExtended gammaProductActive
  ring

/-- **Positivity of the extended product.** -/
theorem gammaProductExtended_pos : 0 < gammaProductExtended := by
  unfold gammaProductExtended
  exact mul_pos (mul_pos (mul_pos
    (mul_pos gammaQ_3_pos gammaQ_5_pos) gammaQ_7_pos) gammaQ_11_pos) gammaQ_13_pos

/-- **Tight bracket for the extended product.**
    `0.050 < γ_3 γ_5 γ_7 γ_11 γ_13 < 0.051`. Numerical value `≈ 0.05078`. -/
theorem gammaProductExtended_bracket :
    50 / 1000 < gammaProductExtended ∧ gammaProductExtended < 51 / 1000 := by
  unfold gammaProductExtended gammaQ deltaQ qPT muStar
  refine ⟨?_, ?_⟩ <;> norm_num

/-- **Extended < active.** Multiplying by two more sub-unit factors `γ_11`
    and `γ_13` strictly decreases the product. -/
theorem gammaProductExtended_lt_active :
    gammaProductExtended < gammaProductActive := by
  rw [gammaProductExtended_eq]
  -- gammaProductActive * γ_11 * γ_13 < gammaProductActive
  -- since 0 < gammaProductActive, γ_11 < 1, γ_13 < 1
  have hpos : 0 < gammaProductActive := gammaProductActive_pos
  have h11lt1 : gammaQ 11 < 1 := gammaQ_11_lt_one
  have h13lt1 : gammaQ 13 < 1 := gammaQ_13_lt_one
  have h11pos : 0 < gammaQ 11 := gammaQ_11_pos
  -- gammaProductActive * gammaQ 11 < gammaProductActive
  have step1 : gammaProductActive * gammaQ 11 < gammaProductActive := by
    have := mul_lt_mul_of_pos_left h11lt1 hpos
    simpa using this
  have step1_pos : 0 < gammaProductActive * gammaQ 11 :=
    mul_pos hpos h11pos
  -- (gammaProductActive * γ_11) * γ_13 < gammaProductActive * γ_11
  have step2 : (gammaProductActive * gammaQ 11) * gammaQ 13
                < gammaProductActive * gammaQ 11 := by
    have := mul_lt_mul_of_pos_left h13lt1 step1_pos
    simpa using this
  exact lt_trans step2 step1

/-! ### Strict-decreasing chain of partial products -/

/-- Partial product up to prime `3`. -/
def partialP1 : ℚ := gammaQ 3
/-- Partial product up to prime `5`. -/
def partialP2 : ℚ := gammaQ 3 * gammaQ 5
/-- Partial product up to prime `7` — equals `gammaProductActive`. -/
def partialP3 : ℚ := gammaQ 3 * gammaQ 5 * gammaQ 7
/-- Partial product up to prime `11`. -/
def partialP4 : ℚ := gammaQ 3 * gammaQ 5 * gammaQ 7 * gammaQ 11
/-- Partial product up to prime `13` — equals `gammaProductExtended`. -/
def partialP5 : ℚ :=
  gammaQ 3 * gammaQ 5 * gammaQ 7 * gammaQ 11 * gammaQ 13

/-- `partialP3 = gammaProductActive`. -/
theorem partialP3_eq_active : partialP3 = gammaProductActive := rfl

/-- `partialP5 = gammaProductExtended`. -/
theorem partialP5_eq_extended : partialP5 = gammaProductExtended := rfl

/-- Strict decrease `P_1 > P_2` (because `γ_5 < 1`, `P_1 > 0`). -/
theorem partialP1_gt_P2 : partialP1 > partialP2 := by
  unfold partialP1 partialP2
  have h : gammaQ 3 * gammaQ 5 < gammaQ 3 * 1 :=
    mul_lt_mul_of_pos_left gammaQ_5_lt_one gammaQ_3_pos
  simpa using h

/-- Strict decrease `P_2 > P_3`. -/
theorem partialP2_gt_P3 : partialP2 > partialP3 := by
  unfold partialP2 partialP3
  have hpos : 0 < gammaQ 3 * gammaQ 5 := mul_pos gammaQ_3_pos gammaQ_5_pos
  have h : gammaQ 3 * gammaQ 5 * gammaQ 7 < gammaQ 3 * gammaQ 5 * 1 :=
    mul_lt_mul_of_pos_left gammaQ_7_lt_one hpos
  simpa using h

/-- Strict decrease `P_3 > P_4`. -/
theorem partialP3_gt_P4 : partialP3 > partialP4 := by
  unfold partialP3 partialP4
  have hpos : 0 < gammaQ 3 * gammaQ 5 * gammaQ 7 :=
    mul_pos (mul_pos gammaQ_3_pos gammaQ_5_pos) gammaQ_7_pos
  have h : gammaQ 3 * gammaQ 5 * gammaQ 7 * gammaQ 11
            < gammaQ 3 * gammaQ 5 * gammaQ 7 * 1 :=
    mul_lt_mul_of_pos_left gammaQ_11_lt_one hpos
  simpa using h

/-- Strict decrease `P_4 > P_5`. -/
theorem partialP4_gt_P5 : partialP4 > partialP5 := by
  unfold partialP4 partialP5
  have hpos : 0 < gammaQ 3 * gammaQ 5 * gammaQ 7 * gammaQ 11 :=
    mul_pos (mul_pos (mul_pos gammaQ_3_pos gammaQ_5_pos) gammaQ_7_pos)
      gammaQ_11_pos
  have h : gammaQ 3 * gammaQ 5 * gammaQ 7 * gammaQ 11 * gammaQ 13
            < gammaQ 3 * gammaQ 5 * gammaQ 7 * gammaQ 11 * 1 :=
    mul_lt_mul_of_pos_left gammaQ_13_lt_one hpos
  simpa using h

/-- **Full strict-decreasing chain** of partial products
    `γ_3 > γ_3 γ_5 > γ_3 γ_5 γ_7 > γ_3 γ_5 γ_7 γ_11 > γ_3 γ_5 γ_7 γ_11 γ_13`. -/
theorem partial_products_chain :
    partialP1 > partialP2
    ∧ partialP2 > partialP3
    ∧ partialP3 > partialP4
    ∧ partialP4 > partialP5 :=
  ⟨partialP1_gt_P2, partialP2_gt_P3, partialP3_gt_P4, partialP4_gt_P5⟩

/-! ### Headline -/

/-- **Headline (γ-product cascade).** At `q_+ = 13/15`, `μ* = 15`:

    1. The active product satisfies `0 < Γ_active < 1` and lies in the
       tight bracket `(0.334, 0.335)`.
    2. `Γ_active ≠ α_bare` — the two products are distinct, with ratio
       `Γ_active / α_bare ∈ (45.6, 45.7)`.
    3. The extended product `Γ_ext = Γ_active · γ_11 · γ_13` lies in
       `(0.050, 0.051)` and is strictly below `Γ_active`.
    4. Partial products are strictly decreasing across the table
       `P_1 > P_2 > P_3 > P_4 > P_5` (each `γ_p < 1`). -/
theorem gamma_product_summary :
    0 < gammaProductActive
    ∧ gammaProductActive < 1
    ∧ 334 / 1000 < gammaProductActive ∧ gammaProductActive < 335 / 1000
    ∧ gammaProductActive ≠ alphaBareQ
    ∧ 456 / 10 < gammaProductActive / alphaBareQ
    ∧ gammaProductActive / alphaBareQ < 457 / 10
    ∧ 50 / 1000 < gammaProductExtended
    ∧ gammaProductExtended < 51 / 1000
    ∧ gammaProductExtended < gammaProductActive
    ∧ partialP1 > partialP2
    ∧ partialP2 > partialP3
    ∧ partialP3 > partialP4
    ∧ partialP4 > partialP5 :=
  ⟨gammaProductActive_pos,
   gammaProductActive_lt_one,
   gammaProductActive_bracket.1,
   gammaProductActive_bracket.2,
   gammaProductActive_ne_alphaBare,
   gammaProductActive_over_alphaBare_bracket.1,
   gammaProductActive_over_alphaBare_bracket.2,
   gammaProductExtended_bracket.1,
   gammaProductExtended_bracket.2,
   gammaProductExtended_lt_active,
   partialP1_gt_P2,
   partialP2_gt_P3,
   partialP3_gt_P4,
   partialP4_gt_P5⟩

end PT.Holonomy
