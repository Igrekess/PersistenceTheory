/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.EML.QSheffer
import PT.EML.EMLIdentities
import Mathlib.Analysis.SpecialFunctions.Exp
import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.Topology.ContinuousOn
import Mathlib.Tactic

/-!
# EML — Algebraic properties (App O extension)

This file collects **structural algebraic properties** of the EML operator
`eml(x, y) := exp(x) - log(y)` (cf. `PT.EML.QSheffer`), complementing the
identities of `EMLIdentities` and the depth-3 constructions of `EMLDepth3`:

* **Non-associativity (concrete witness)** — `eml(eml(1,1), 1) ≠ eml(1, eml(1,1))`.
  The left-hand side equals `exp(exp 1)` (very large); the right-hand side
  equals `exp 1 - 1` (small). A clean structural negative result.
* **No left/right identity** — for every candidate `e ∈ ℝ`, we exhibit a
  witness `y` for which `eml(e, y) ≠ y` and `eml(y, e) ≠ y`. This formalises
  that the EML algebra `(ℝ, eml)` has neither a left nor a right neutral
  element.
* **Iterated `eml(·, 1)`** — the right-unit specialisation gives an
  exponentiation tower: `eml(eml(x, 1), 1) = exp(exp x)`.
* **Joint continuity** on the half-plane `y > 0` — extends the per-variable
  continuity lemmas of `EMLIdentities`.
* **First-argument bijection** — for any fixed `y`, the map `x ↦ eml(x, y)`
  is a bijection from `ℝ` onto its range `(-log y, +∞)` (a translated copy
  of `Set.Ioi 0`, the image of `exp`).

## Reference

Appendix O of the monograph (`app:eml`), structural follow-up to
Propositions O.1–O.3. Odrzywołek, *EML: a continuous Sheffer primitive*,
2026, §2 (algebraic properties).
-/

namespace PT.EML

open Real

/-! ### Non-associativity (concrete witness) -/

/-- `eml(eml(1, 1), 1) = exp(exp 1)`. The two right-unit collapses
    compose into an exponentiation tower. -/
theorem eml_eml_one_one_one_left :
    eml (eml 1 1) 1 = Real.exp (Real.exp 1) := by
  rw [eml_arg2_one, eml_arg2_one]

/-- `eml(1, eml(1, 1)) = exp 1 - 1`. Inner right-unit collapse, then the
    outer EML exponentiates `1` and subtracts `log(exp 1) = 1`. -/
theorem eml_one_eml_one_one_right :
    eml 1 (eml 1 1) = Real.exp 1 - 1 := by
  rw [eml_arg2_one]
  unfold eml
  rw [Real.log_exp]

/-- **Non-associativity (concrete witness at `a = b = c = 1`).**
    `eml(eml(1, 1), 1) = exp(exp 1)` is large (`> e > 2`), while
    `eml(1, eml(1, 1)) = exp 1 - 1 < 2`. -/
theorem eml_not_associative :
    eml (eml 1 1) 1 ≠ eml 1 (eml 1 1) := by
  rw [eml_eml_one_one_one_left, eml_one_eml_one_one_right]
  -- Goal: exp (exp 1) ≠ exp 1 - 1
  intro h
  -- We have exp 1 > 1, hence exp(exp 1) > exp 1, and so > exp 1 - 1.
  have h_exp1_pos : (0 : ℝ) < Real.exp 1 := Real.exp_pos 1
  have h_one_lt_exp1 : (1 : ℝ) < Real.exp 1 := by
    have : Real.exp 0 < Real.exp 1 :=
      Real.exp_lt_exp.mpr (by norm_num : (0 : ℝ) < 1)
    rw [Real.exp_zero] at this
    exact this
  have h_exp1_lt : Real.exp 1 < Real.exp (Real.exp 1) :=
    Real.exp_lt_exp.mpr h_one_lt_exp1
  -- Combine: exp(exp 1) > exp 1 > exp 1 - 1
  linarith

/-! ### No left or right identity -/

/-- **No right identity.** For every candidate `e ∈ ℝ`, there exists `y`
    with `eml(y, e) ≠ y`. Witness: `y = 0`. Then `eml(0, e) = 1 - log e`,
    which equals `0` only if `log e = 1`, i.e. `e = exp 1`. But even then
    we have a different obstruction: take `y = 1`, so `eml(1, e) = exp 1 - log e`.
    For this to equal `1`, we need `log e = exp 1 - 1 ≠ 0 = log 1`, i.e.
    `e ≠ 1`. Combining `y = 0` (fails unless `log e = 1`) with `y = 1`
    (fails unless `log e = exp 1 - 1`) covers all `e` since the two
    required values of `log e` differ. -/
theorem eml_no_right_identity (e : ℝ) : ∃ y : ℝ, eml y e ≠ y := by
  -- Try y = 0 first. eml(0, e) = 1 - log e ≠ 0 iff log e ≠ 1.
  by_cases hloge : Real.log e = 1
  · -- Then try y = 1. eml(1, e) = exp 1 - log e = exp 1 - 1 ≠ 1
    -- since exp 1 ≠ 2 (which we prove by exp 1 < 3 and a direct bound,
    -- or more simply: exp 1 - 1 = 1 ↔ exp 1 = 2, but exp 1 > 2 by
    -- the strict-convexity bound exp x > 1 + x at x = 1: actually
    -- Real.add_one_lt_exp gives exp 1 > 2, so exp 1 - 1 > 1.
    refine ⟨1, ?_⟩
    unfold eml
    rw [hloge]
    -- Goal: exp 1 - 1 ≠ 1
    intro h
    have hexp1 : Real.exp 1 = 2 := by linarith
    -- But exp 1 > 1 + 1 = 2 strictly (Real.add_one_lt_exp at x = 1, x ≠ 0)
    have h_strict : (1 : ℝ) + 1 < Real.exp 1 :=
      Real.add_one_lt_exp (by norm_num : (1 : ℝ) ≠ 0)
    linarith
  · refine ⟨0, ?_⟩
    unfold eml
    rw [Real.exp_zero]
    -- Goal: 1 - log e ≠ 0
    intro h
    apply hloge
    linarith

/-- **No left identity.** For every candidate `e ∈ ℝ`, there exists `y`
    with `eml(e, y) ≠ y`. Witness: choose `y = 0`. Then
    `eml(e, 0) = exp e - log 0 = exp e - 0 = exp e ≠ 0` since `exp e > 0`. -/
theorem eml_no_left_identity (e : ℝ) : ∃ y : ℝ, eml e y ≠ y := by
  refine ⟨0, ?_⟩
  unfold eml
  rw [Real.log_zero, sub_zero]
  -- Goal: exp e ≠ 0
  exact ne_of_gt (Real.exp_pos e)

/-! ### Iterated `eml(·, 1)`: the exponential tower -/

/-- **Iterated right-unit collapse.** `eml(eml(x, 1), 1) = exp(exp x)`. -/
theorem eml_eml_arg2_one_arg2_one_general (x : ℝ) :
    eml (eml x 1) 1 = Real.exp (Real.exp x) := by
  rw [eml_arg2_one, eml_arg2_one]

/-- **Triple iteration.** `eml(eml(eml(x, 1), 1), 1) = exp(exp(exp x))`. -/
theorem eml_eml_eml_arg2_one (x : ℝ) :
    eml (eml (eml x 1) 1) 1 = Real.exp (Real.exp (Real.exp x)) := by
  rw [eml_arg2_one, eml_arg2_one, eml_arg2_one]

/-! ### Joint continuity on the open half-plane `y > 0` -/

/-- **Joint continuity of `eml` on the open half-plane `{(x, y) : y > 0}`.**
    Combines `Real.continuous_exp` with `Real.continuousAt_log` on the
    nonzero condition. -/
theorem continuousOn_eml_jointly :
    ContinuousOn (fun p : ℝ × ℝ => eml p.1 p.2) {p : ℝ × ℝ | 0 < p.2} := by
  unfold eml
  refine ContinuousOn.sub ?_ ?_
  · -- exp ∘ p.1 is continuous everywhere, in particular on the set.
    exact (Real.continuous_exp.comp continuous_fst).continuousOn
  · -- log ∘ p.2 is continuous at every point where p.2 ≠ 0.
    intro p hp
    have hp2 : p.2 ≠ 0 := ne_of_gt hp
    have h_log : ContinuousAt Real.log p.2 := Real.continuousAt_log hp2
    have h_snd : ContinuousAt (fun q : ℝ × ℝ => q.2) p := continuous_snd.continuousAt
    exact (h_log.comp h_snd).continuousWithinAt

/-! ### First-argument bijection onto a translated half-line -/

/-- **First-argument injectivity.** For any fixed `y`, the map
    `x ↦ eml(x, y)` is injective (it is strictly monotone, hence injective). -/
theorem eml_left_injective (y : ℝ) :
    Function.Injective (fun x => eml x y) :=
  (eml_strictMono_left y).injective

/-- **Range of `eml(·, y)`.** The image of `x ↦ eml(x, y) = exp x - log y`
    is exactly `Set.Ioi (-Real.log y)`: the translated copy of `(0, +∞)`
    obtained by shifting the range of `exp` down by `log y`. -/
theorem range_eml_left (y : ℝ) :
    Set.range (fun x => eml x y) = Set.Ioi (-(Real.log y)) := by
  ext z
  constructor
  · rintro ⟨x, hx⟩
    -- z = exp x - log y > 0 - log y = -log y
    simp only [Set.mem_Ioi]
    have h_pos : (0 : ℝ) < Real.exp x := Real.exp_pos x
    have hx' : Real.exp x - Real.log y = z := hx
    linarith
  · intro hz
    -- z > -log y, i.e. z + log y > 0; pick x = log (z + log y)
    simp only [Set.mem_Ioi] at hz
    have h_pos : (0 : ℝ) < z + Real.log y := by linarith
    refine ⟨Real.log (z + Real.log y), ?_⟩
    show eml (Real.log (z + Real.log y)) y = z
    unfold eml
    rw [Real.exp_log h_pos]
    ring

/-- **Bijection onto range.** For any fixed `y`, the map `x ↦ eml(x, y)`
    is a bijection from `ℝ` onto `Set.Ioi (-(Real.log y))`. -/
theorem eml_left_bijOn (y : ℝ) :
    Set.BijOn (fun x => eml x y) Set.univ (Set.Ioi (-(Real.log y))) := by
  refine ⟨?_, ?_, ?_⟩
  · -- MapsTo
    intro x _
    have h : eml x y ∈ Set.range (fun x => eml x y) := ⟨x, rfl⟩
    rw [range_eml_left] at h
    exact h
  · -- InjOn
    exact (eml_left_injective y).injOn
  · -- SurjOn
    intro z hz
    have hrange : z ∈ Set.range (fun x => eml x y) := by
      rw [range_eml_left]; exact hz
    obtain ⟨x, hx⟩ := hrange
    exact ⟨x, Set.mem_univ x, hx⟩

/-! ### Headline -/

/-- **Headline.** Algebraic structure of `(ℝ, eml)`:

    * **Non-associative**: `eml(eml(1,1), 1) ≠ eml(1, eml(1,1))`
      (concrete witness, `exp(exp 1)` vs `exp 1 - 1`).
    * **No two-sided identity**: for every `e`, there exist witnesses
      `y_L, y_R` with `eml(e, y_R) ≠ y_R` and `eml(y_L, e) ≠ y_L`.
    * **Right-unit tower**: iterating `eml(·, 1)` builds an exponential
      tower `exp(exp(... exp x))`.
    * **Jointly continuous** on `{y > 0}`.
    * **Bijection in the first variable** onto the translated half-line
      `(-log y, +∞)`.

    Together with the depth hierarchy of `EMLDepth3`, this paints the EML
    algebra as a non-associative, identity-free, but topologically and
    order-theoretically well-behaved Sheffer primitive — exactly the
    profile required by App O.
-/
theorem eml_algebra_headline :
    -- Non-associativity at (1, 1, 1)
    eml (eml 1 1) 1 ≠ eml 1 (eml 1 1)
    -- No right identity
    ∧ (∀ e : ℝ, ∃ y : ℝ, eml y e ≠ y)
    -- No left identity
    ∧ (∀ e : ℝ, ∃ y : ℝ, eml e y ≠ y)
    -- Right-unit iteration produces an exponential tower
    ∧ (∀ x : ℝ, eml (eml x 1) 1 = Real.exp (Real.exp x))
    -- First-argument injectivity
    ∧ (∀ y : ℝ, Function.Injective (fun x => eml x y)) := by
  refine ⟨eml_not_associative, eml_no_right_identity, eml_no_left_identity,
          eml_eml_arg2_one_arg2_one_general, eml_left_injective⟩

end PT.EML
