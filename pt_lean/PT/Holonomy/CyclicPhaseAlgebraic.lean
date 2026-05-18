/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Holonomy.ActivePrimeCriterion
import PT.Holonomy.CouplingReconstruction
import PT.Holonomy.CyclicPhaseIdentity
import PT.Holonomy.CyclicPhaseInversion
import Mathlib.Tactic
import Mathlib.Data.Rat.Defs

/-!
# Cyclic Phase — Algebraic content of Route 3 (Fisher), audit row #25

**Statement (paper-level, Ch06b).** Route 3 of the cyclic-phase theorem
interprets `sin²θ_p` as a component of the information-geometric Fisher
curvature on the probability simplex. The full Riemannian-geometric
derivation requires `Mathlib.Geometry.Manifold.SmoothManifoldsWithCorners`
maturity and is therefore HARD (audit row #25, deferred). However, the
*algebraic substrate* of Route 3 — the closed-form expression of `sin²θ_p`
directly in terms of the gap fraction `δ_p(q) = (1 - q^p)/p` and its
monotone dependence on the branch parameter `q` — is fully tractable.

This file records:

1. **Direct algebraic form.** `sin²θ_p = δ_p(q)(2 - δ_p(q))` rephrased as
   a polynomial in `q^p` only (no `cos θ` intermediary):
   `sin²θ_p = (1 - q^{2p}) / p² · p / 1 - …` — concretely
   `δ(2 - δ) = 1 - (1 - δ)²` and `1 - δ_p = 1 - (1-q^p)/p = (p - 1 + q^p)/p`,
   giving `sin²θ_p = 1 - ((p - 1 + q^p)/p)²`.
2. **Strict monotonicity in `q` (Route 3 algebraic content).** On the
   PT-relevant interval `q ∈ [0, 1]`, the squared sine `sin²θ_p` is a
   *strictly decreasing* function of `q` (since `δ_p` decreases in `q^p`
   and `δ ↦ δ(2 - δ)` is strictly increasing on `[0, 1]` — `CyclicPhaseInversion`).
3. **Exact rational values.** For `q = qPT = 13/15` and `p ∈ {3, 5, 7}`,
   `δ_p(qPT)` and `sin²θ_p(qPT)` are exact rationals; we compute them and
   show they agree with the brackets of `CouplingReconstruction.sinSqQ`.

Note: Route 3 algebraic monotonicity here is in the **branch parameter `q`**,
parallel to (and independent from) `CyclicPhaseInversion`'s monotonicity in
`δ` and `CouplingReconstruction`'s monotonicity in `p` (across active primes).

## Reference

* Monograph Chapter 6b, §"Route 3: Fisher". Audit row #25 algebraic content.
* `app_e_constants.tex`: numerical values `δ_3, δ_5, δ_7`.
-/

namespace PT.Holonomy

/-! ### Algebraic form of `sin²θ_p` directly in `q^p`

We give a closed-form polynomial expression of `sin²θ_p` in the variable
`q^p`, without any intermediate `cos θ`. -/

/-- **Direct algebraic identity (no `cos`).** For any real `q` and any
    positive natural `p`, with `δ_p = (1 - q^p)/p`,

    `δ_p (2 - δ_p) = 1 - ((p - 1 + q^p)/p)²`.

    Polynomial identity in `q^p` and `p`. -/
theorem sin_sq_algebraic_form (q : ℝ) (p : ℕ) (hp : 0 < (p : ℝ)) :
    let δ := (1 - q^p) / (p : ℝ)
    δ * (2 - δ) = 1 - ((((p : ℝ) - 1 + q^p) / (p : ℝ))) ^ 2 := by
  intro δ
  have hp_ne : (p : ℝ) ≠ 0 := ne_of_gt hp
  -- δ(2-δ) = 1 - (1-δ)², and 1 - δ = 1 - (1-q^p)/p = (p - 1 + q^p)/p
  have h1 : 1 - δ = ((p : ℝ) - 1 + q^p) / (p : ℝ) := by
    show 1 - (1 - q^p) / (p : ℝ) = ((p : ℝ) - 1 + q^p) / (p : ℝ)
    field_simp
    ring
  have h2 : δ * (2 - δ) = 1 - (1 - δ)^2 := by ring
  rw [h2, h1]

/-- **Direct algebraic identity (rational version, `ℚ`).** Same as
    `sin_sq_algebraic_form` but in exact rationals. -/
theorem sin_sq_algebraic_form_rat (q : ℚ) (p : ℕ) (hp : 0 < (p : ℚ)) :
    let δ := (1 - q^p) / (p : ℚ)
    δ * (2 - δ) = 1 - ((((p : ℚ) - 1 + q^p) / (p : ℚ))) ^ 2 := by
  intro δ
  have hp_ne : (p : ℚ) ≠ 0 := ne_of_gt hp
  have h1 : 1 - δ = ((p : ℚ) - 1 + q^p) / (p : ℚ) := by
    show 1 - (1 - q^p) / (p : ℚ) = ((p : ℚ) - 1 + q^p) / (p : ℚ)
    field_simp
    ring
  have h2 : δ * (2 - δ) = 1 - (1 - δ)^2 := by ring
  rw [h2, h1]

/-! ### Strict monotonicity in `q` (Route 3 algebraic content)

For fixed `p ≥ 1`, `δ_p(q) = (1 - q^p)/p` is strictly decreasing on
`q ∈ [0, 1]`. Composing with the strictly increasing `δ ↦ δ(2 - δ)` on
`[0, 1]` yields the strict-decrease statement for `sin²θ_p` in `q`. -/

/-- For `p ≥ 1` and `q, q' ∈ [0, 1]` with `q < q'`, we have
    `q^p < q'^p`. -/
theorem pow_strictMono_on_unit (p : ℕ) (hp : 1 ≤ p) :
    ∀ q q' : ℝ, 0 ≤ q → q < q' → q' ≤ 1 → q^p < q'^p := by
  intros q q' h0q hqq' _
  exact pow_lt_pow_left₀ hqq' h0q (by omega : p ≠ 0)

/-- For `p ≥ 1`, the gap fraction `δ_p(q) = (1 - q^p)/p` is strictly
    decreasing in `q` on `[0, 1]`. -/
theorem delta_strictAnti_q (p : ℕ) (hp : 1 ≤ p) :
    ∀ q q' : ℝ, 0 ≤ q → q < q' → q' ≤ 1 →
      (1 - q'^p) / (p : ℝ) < (1 - q^p) / (p : ℝ) := by
  intros q q' h0q hqq' hq'1
  have hp_pos : (0 : ℝ) < (p : ℝ) := by exact_mod_cast hp
  have hpow_lt : q^p < q'^p := pow_strictMono_on_unit p hp q q' h0q hqq' hq'1
  -- (1 - q'^p)/p < (1 - q^p)/p ⇔ 1 - q'^p < 1 - q^p ⇔ q^p < q'^p
  exact (div_lt_div_iff_of_pos_right hp_pos).mpr (by linarith)

/-- For `p ≥ 1` and `q ∈ [0, 1]`, the gap fraction `δ_p(q)` lies in `[0, 1/p]`,
    hence in particular in `[0, 1]`. -/
theorem delta_mem_unit (p : ℕ) (hp : 1 ≤ p) :
    ∀ q : ℝ, 0 ≤ q → q ≤ 1 →
      0 ≤ (1 - q^p) / (p : ℝ) ∧ (1 - q^p) / (p : ℝ) ≤ 1 := by
  intros q h0q hq1
  have hp_pos : (0 : ℝ) < (p : ℝ) := by exact_mod_cast hp
  have hpow_le : q^p ≤ 1 := pow_le_one₀ h0q hq1
  have hpow_nn : 0 ≤ q^p := pow_nonneg h0q p
  refine ⟨?_, ?_⟩
  · exact div_nonneg (by linarith) (le_of_lt hp_pos)
  · -- (1 - q^p)/p ≤ 1 ⇔ 1 - q^p ≤ p, which follows since 1 - q^p ≤ 1 ≤ p.
    have h1 : 1 - q^p ≤ 1 := by linarith
    have hp_ge_one : (1 : ℝ) ≤ (p : ℝ) := by exact_mod_cast hp
    have : 1 - q^p ≤ (p : ℝ) := le_trans h1 hp_ge_one
    exact (div_le_one hp_pos).mpr this

/-- **Strict monotonicity in `q` (Route 3 algebraic statement).**
    For `p ≥ 1` and `q, q' ∈ [0, 1]` with `q < q'`,

    `sin²θ_p(q') < sin²θ_p(q)`,

    i.e. the squared sine is a strictly *decreasing* function of the branch
    parameter `q`. Proof: `δ_p` is strictly decreasing in `q` on `[0, 1]`
    (`delta_strictAnti_q`), both `δ_p(q), δ_p(q')` lie in `[0, 1]`
    (`delta_mem_unit`), and the map `δ ↦ δ(2 - δ)` is strictly increasing on
    `[0, 1]` (`delta_two_minus_delta_strictMono`). -/
theorem sin_sq_strictAnti_q (p : ℕ) (hp : 1 ≤ p) :
    ∀ q q' : ℝ, 0 ≤ q → q < q' → q' ≤ 1 →
      let δ' := (1 - q'^p) / (p : ℝ)
      let δ  := (1 - q^p)  / (p : ℝ)
      δ' * (2 - δ') < δ * (2 - δ) := by
  intros q q' h0q hqq' hq'1 δ' δ
  have h0q' : 0 ≤ q' := le_of_lt (lt_of_le_of_lt h0q hqq')
  have hq_le : q ≤ 1 := le_of_lt (lt_of_lt_of_le hqq' hq'1)
  obtain ⟨h0δ, h1δ⟩ := delta_mem_unit p hp q h0q hq_le
  obtain ⟨h0δ', h1δ'⟩ := delta_mem_unit p hp q' h0q' hq'1
  have hδlt : δ' < δ := delta_strictAnti_q p hp q q' h0q hqq' hq'1
  exact delta_two_minus_delta_strictMono δ' δ h0δ' h1δ' h0δ h1δ hδlt

/-! ### Exact rational values at `q = qPT = 13/15`

We give the explicit closed-form rationals for `δ_p(qPT)` and `sin²θ_p(qPT)`
for `p ∈ {3, 5, 7}`. These are the values that underlie the brackets
proved in `CouplingReconstruction.sinSq_{3,5,7}_bracket`. -/

/-- Exact rational value of `δ_3(qPT)` in power-form representation:
    `δ_3 = (1 - (13/15)^3)/3 = (15^3 - 13^3)/(3 · 15^3) = 1178/(3·3375) = 1178/10125`.

    Renamed from `deltaQ_3_eq` (2026-05-17) to avoid namespace collision
    with the decimal-form `deltaQ_3_eq` in `CyclicPhaseTable.lean`. -/
theorem deltaQ_3_eq_pow_form : deltaQ 3 = 1178 / 10125 := by
  unfold deltaQ qPT
  norm_num

/-- Exact rational value of `δ_5(qPT)` in power-form representation
    `(15^5 - 13^5)/(5 · 15^5)`. Decimal form `388082/3796875` in
    `CyclicPhaseTable.lean`. -/
theorem deltaQ_5_eq_pow_form : deltaQ 5 = (15^5 - 13^5) / (5 * 15^5) := by
  unfold deltaQ qPT
  norm_num

/-- Exact rational value of `δ_7(qPT)` in power-form representation
    `(15^7 - 13^7)/(7 · 15^7)`. Decimal form `108110858/1196015625` in
    `CyclicPhaseTable.lean`. -/
theorem deltaQ_7_eq_pow_form : deltaQ 7 = (15^7 - 13^7) / (7 * 15^7) := by
  unfold deltaQ qPT
  norm_num

/-- Numerical bracket for `δ_3(qPT)`: `0.116 < δ_3 < 0.117`
    (published value: `δ_3 ≈ 0.1163`). -/
theorem deltaQ_3_bracket : 116 / 1000 < deltaQ 3 ∧ deltaQ 3 < 117 / 1000 := by
  unfold deltaQ qPT
  refine ⟨?_, ?_⟩ <;> norm_num

/-- Numerical bracket for `δ_5(qPT)`: `0.102 < δ_5 < 0.103`. -/
theorem deltaQ_5_bracket : 102 / 1000 < deltaQ 5 ∧ deltaQ 5 < 103 / 1000 := by
  unfold deltaQ qPT
  refine ⟨?_, ?_⟩ <;> norm_num

/-- Numerical bracket for `δ_7(qPT)`: `0.090 < δ_7 < 0.091`. -/
theorem deltaQ_7_bracket : 90 / 1000 < deltaQ 7 ∧ deltaQ 7 < 91 / 1000 := by
  unfold deltaQ qPT
  refine ⟨?_, ?_⟩ <;> norm_num

/-! ### Consistency with `CouplingReconstruction.sinSqQ`

`sinSqQ p` (from `CouplingReconstruction`) is exactly `δ_p(2 - δ_p)`, hence
the algebraic form derived above coincides with it for `p ∈ {3, 5, 7}`. -/

/-- `sin²θ_p` (in the `sinSqQ` form) equals the direct algebraic form
    `1 - ((p - 1 + q^p)/p)²` for any positive prime exponent `p`. -/
theorem sinSqQ_eq_algebraic (p : ℕ) (hp : 0 < (p : ℚ)) :
    sinSqQ p = 1 - ((((p : ℚ) - 1 + qPT^p) / (p : ℚ))) ^ 2 := by
  unfold sinSqQ
  exact sin_sq_algebraic_form_rat qPT p hp

/-- **Closed form for `sin²θ_3(qPT)`.** Exact rational value. -/
theorem sinSqQ_3_eq_algebraic :
    sinSqQ 3 = 1 - ((((3 : ℚ) - 1 + qPT^3) / (3 : ℚ))) ^ 2 := by
  exact sinSqQ_eq_algebraic 3 (by norm_num)

/-- **Closed form for `sin²θ_5(qPT)`.** -/
theorem sinSqQ_5_eq_algebraic :
    sinSqQ 5 = 1 - ((((5 : ℚ) - 1 + qPT^5) / (5 : ℚ))) ^ 2 := by
  exact sinSqQ_eq_algebraic 5 (by norm_num)

/-- **Closed form for `sin²θ_7(qPT)`.** -/
theorem sinSqQ_7_eq_algebraic :
    sinSqQ 7 = 1 - ((((7 : ℚ) - 1 + qPT^7) / (7 : ℚ))) ^ 2 := by
  exact sinSqQ_eq_algebraic 7 (by norm_num)

/-! ### Headline: Route 3 algebraic content -/

/-- **Headline (Route 3 algebraic content).** The cyclic-phase squared
    sine admits a direct closed-form expression as a rational function of
    `q^p` (no `cos θ` intermediary), and on the unit interval `q ∈ [0, 1]`
    it is a strictly monotone (decreasing) function of `q` for every
    `p ≥ 1`. The Riemannian-geometric (Fisher) interpretation of this
    monotonicity — namely that `sin²θ_p` is the squared geodesic distance
    component on the simplex `Δ_p` — is the HARD content deferred to
    audit row #25. -/
theorem cyclic_phase_route3_algebraic_summary :
    -- direct algebraic form (real)
    (∀ q : ℝ, ∀ p : ℕ, 0 < (p : ℝ) →
        let δ := (1 - q^p) / (p : ℝ)
        δ * (2 - δ) = 1 - ((((p : ℝ) - 1 + q^p) / (p : ℝ))) ^ 2)
    -- strict antitone in q
    ∧ (∀ p : ℕ, 1 ≤ p → ∀ q q' : ℝ, 0 ≤ q → q < q' → q' ≤ 1 →
        let δ' := (1 - q'^p) / (p : ℝ)
        let δ  := (1 - q^p)  / (p : ℝ)
        δ' * (2 - δ') < δ * (2 - δ)) := by
  refine ⟨sin_sq_algebraic_form, sin_sq_strictAnti_q⟩

end PT.Holonomy
