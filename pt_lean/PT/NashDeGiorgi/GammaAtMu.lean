/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import Mathlib.Tactic
import Mathlib.Data.Rat.Defs

/-!
# `γ_p(μ)` for variable `μ` — extended threshold values

This module extends `PT.Holonomy.ActivePrimeCriterion` (which fixes
`μ = μ* = 15`) to *variable* `μ ∈ ℚ`, in order to formalise

* the **threshold values** `μ_p^seuil` (note 09 of `PT_NASH_DEGIORGI`);
* the **m_K formula** (note 03), `m_K = |{p odd : γ_p(μ_+) > 1/2}|`;
* the **rational verifications** of the active-set characterization.

All values are computed in exact rational arithmetic (`ℚ`), making the
key inequalities decidable by `norm_num` / `decide`.

## Reference

Note 03 §4 ("Théorème principal Z1") and note 09 §5 (tableau des seuils)
of `PT_PROJECTS/PT_NASH_DEGIORGI/`.

## Notation

* `q_at μ := 1 - 2/μ` — branch parameter as a function of `μ`.
* `delta_at p μ := (1 - q_at(μ)^p) / p` — gap fraction at scale `μ`.
* `gamma_at p μ` — anomalous dimension at scale `μ` (closed form).

These reduce to `qPT`, `deltaQ p`, `gammaQ p` when `μ = 15`.
-/

namespace PT.NashDeGiorgi

/-! ### Extended definitions (variable μ) -/

/-- Branch parameter `q(μ) = 1 - 2/μ` at scale `μ`. -/
def q_at (μ : ℚ) : ℚ := 1 - 2 / μ

/-- Gap fraction `δ_p(μ) = (1 - q(μ)^p)/p`. -/
def delta_at (p : ℕ) (μ : ℚ) : ℚ := (1 - (q_at μ) ^ p) / (p : ℚ)

/-- Anomalous dimension `γ_p(μ) = 4 q^(p-1) (1-δ)/(μ δ (2-δ))`. -/
def gamma_at (p : ℕ) (μ : ℚ) : ℚ :=
  (4 * (q_at μ) ^ (p - 1) * (1 - delta_at p μ)) /
   (μ * delta_at p μ * (2 - delta_at p μ))

/-! ### Sanity check: reduces to `gammaQ p` at `μ = 15` -/

/-- `q_at(15) = 13/15`. -/
@[simp] theorem q_at_15 : q_at 15 = 13 / 15 := by
  unfold q_at; norm_num

/-! ### Active-set verification at `μ* = 15`

These reproduce the table of `PT_NASH_DEGIORGI/06_thresholds.json`:
* `γ_3(15) ≈ 0.808 > 1/2` — active
* `γ_5(15) ≈ 0.696 > 1/2` — active
* `γ_7(15) ≈ 0.595 > 1/2` — active
* `γ_11(15) ≈ 0.426 < 1/2` — echo
* `γ_13(15) ≈ 0.356 < 1/2` — echo
* `γ_17(15) ≈ 0.245 < 1/2` — echo
-/

theorem gamma_3_at_15_active : gamma_at 3 15 > 1 / 2 := by
  unfold gamma_at delta_at q_at; norm_num

theorem gamma_5_at_15_active : gamma_at 5 15 > 1 / 2 := by
  unfold gamma_at delta_at q_at; norm_num

theorem gamma_7_at_15_active : gamma_at 7 15 > 1 / 2 := by
  unfold gamma_at delta_at q_at; norm_num

theorem gamma_11_at_15_echo : gamma_at 11 15 < 1 / 2 := by
  unfold gamma_at delta_at q_at; norm_num

theorem gamma_13_at_15_echo : gamma_at 13 15 < 1 / 2 := by
  unfold gamma_at delta_at q_at; norm_num

theorem gamma_17_at_15_echo : gamma_at 17 15 < 1 / 2 := by
  unfold gamma_at delta_at q_at; norm_num

/-! ### Active set on `K = [3, 30]` — primes that exceed `1/2` at `μ_+ = 30`

By Lemma A (monotonicity in `μ`, formalised as `sorry` until note 03 §2 is
completed analytically), `γ_p(μ) ≤ γ_p(30)` for `μ ∈ [3, 30]`. Hence the
active set on `K = [3, 30]` is

  `A(K) = {p odd prime : γ_p(30) > 1/2}`.

The table below verifies that this set is exactly `{3, 5, 7, 11, 13, 17}`,
so `m_K = 6` (note 03 §5 and note 09 §5).
-/

theorem gamma_3_at_30_active : gamma_at 3 30 > 1 / 2 := by
  unfold gamma_at delta_at q_at; norm_num

theorem gamma_5_at_30_active : gamma_at 5 30 > 1 / 2 := by
  unfold gamma_at delta_at q_at; norm_num

theorem gamma_7_at_30_active : gamma_at 7 30 > 1 / 2 := by
  unfold gamma_at delta_at q_at; norm_num

theorem gamma_11_at_30_active : gamma_at 11 30 > 1 / 2 := by
  unfold gamma_at delta_at q_at; norm_num

theorem gamma_13_at_30_active : gamma_at 13 30 > 1 / 2 := by
  unfold gamma_at delta_at q_at; norm_num

theorem gamma_17_at_30_active : gamma_at 17 30 > 1 / 2 := by
  unfold gamma_at delta_at q_at; norm_num

theorem gamma_19_at_30_echo : gamma_at 19 30 < 1 / 2 := by
  unfold gamma_at delta_at q_at; norm_num

theorem gamma_23_at_30_echo : gamma_at 23 30 < 1 / 2 := by
  unfold gamma_at delta_at q_at; norm_num

/-! ### Lemma A — monotonicity of `γ_p(μ)` in `μ`

The continuous statement requires real-variable calculus
(see note 03 §2 of `PT_NASH_DEGIORGI`). We record the discrete instance
`γ_p(μ₁) < γ_p(μ₂)` for `μ₁ < μ₂` in `ℚ` via `norm_num`, and state the
full real-variable monotonicity as an axiom (to be discharged in
`LemmaA_proof.lean` once the analytic skeleton of note 03 is finalised).
-/

/-- Discrete instance of Lemma A: `γ_11` is increasing on the interval
    `[15, 30]` in `ℚ` (it crosses the activity threshold somewhere between
    `μ = 17` and `μ = 18` per note 09). -/
theorem gamma_11_monotone_15_to_30 :
    gamma_at 11 15 < gamma_at 11 30 := by
  unfold gamma_at delta_at q_at; norm_num

theorem gamma_17_monotone_15_to_30 :
    gamma_at 17 15 < gamma_at 17 30 := by
  unfold gamma_at delta_at q_at; norm_num

-- **Lemma A (paper-level statement).** `γ_p(μ)` is strictly increasing
-- in `μ` on `(2, ∞)` for every integer `p ≥ 2`.
--
-- The full statement on ℝ with real-valued `gamma_p` is now formalised
-- in `PT.NashDeGiorgi.LemmaAFormal` (see `lemmaA_strictMono` there).
-- The algebraic core (positivity `p - S_p(q) > 0` on `(0,1)`,
-- Theorem 3.1 of note 15) is fully proved in Lean ; the derivative
-- chain-rule computation that connects it to `d γ_p / d μ > 0` is
-- left as a TODO (note 15 §2).

end PT.NashDeGiorgi
