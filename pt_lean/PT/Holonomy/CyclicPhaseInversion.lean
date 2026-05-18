/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Holonomy.CyclicPhaseIdentity
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Basic
import Mathlib.Tactic

/-!
# Cyclic Phase — Inversion and monotonicity (Route 2 proxy, Ch06 #24)

This file extends `PT.Holonomy.CyclicPhaseIdentity` with:

* **Inversion (Jordan-like)** — given a target `s = sin²θ`, recover the gap
  fraction `δ` ∈ `[0, 1]` via `δ = 1 - √(1 - s)`.
* **Strict monotonicity** — the map `δ ↦ δ(2 - δ)` is strictly increasing
  on `[0, 1]`. Hence `sin²θ_p` is a strictly increasing function of `δ_p`
  on the PT-relevant interval `δ_p ∈ (0, 1)`.
* **Range** — `δ(2 - δ) ∈ [0, 1]` for `δ ∈ [0, 1]`, with the endpoints
  attained exactly at `δ = 0` and `δ = 1`.

These are the "Route 2 proxy" of the cyclic-phase theorem: the audit-row
#24 ("Cyclic phase Route 2 (spectral)" — Fourier on `ℤ/pℤ`) requires
character analysis on prime cyclic groups, which is beyond Mathlib's
present prime-group Fourier infrastructure. The **algebraic monotonicity /
inversion content** of Route 2 — i.e., the bijection
`δ ↔ sin²θ` on `[0, 1]` and the strict monotone correspondence
`(δ_p < δ_p') ⇔ (sin²θ_p < sin²θ_p')` — is captured here.

## Reference

Monograph Chapter 6, §6.3 ("Inversion et monotonie"). Audit row #24
algebraic content.
-/

namespace PT.Holonomy

open Real

/-! ### Algebraic auxiliaries -/

/-- Range bound: for `δ ∈ [0, 1]`, `δ(2 - δ) ∈ [0, 1]`. -/
theorem delta_two_minus_delta_le_one (δ : ℝ) (_h₀ : 0 ≤ δ) (h₁ : δ ≤ 1) :
    δ * (2 - δ) ≤ 1 := by
  nlinarith

/-- `δ(2 - δ) = 1 - (1 - δ)²`. -/
theorem delta_two_minus_delta_eq_one_sub_sq (δ : ℝ) :
    δ * (2 - δ) = 1 - (1 - δ)^2 := by ring

/-- Non-negativity of `δ(2 - δ)` on `[0, 2]`. -/
theorem delta_two_minus_delta_nonneg (δ : ℝ) (h₀ : 0 ≤ δ) (h₁ : δ ≤ 2) :
    0 ≤ δ * (2 - δ) := by
  have h2 : 0 ≤ 2 - δ := by linarith
  exact mul_nonneg h₀ h2

/-! ### Strict monotonicity on `[0, 1]` -/

/-- **Strict monotonicity (key lemma).** On `[0, 1]`, the map `δ ↦ δ(2 - δ)`
    is strictly increasing.

    Proof: `δ'(2-δ') - δ(2-δ) = 2(δ' - δ) - (δ'² - δ²) = (δ' - δ)(2 - δ - δ')`,
    and `2 - δ - δ' > 0` since `δ + δ' < 2` when both are in `[0, 1)`. -/
theorem delta_two_minus_delta_strictMono :
    ∀ δ δ' : ℝ, 0 ≤ δ → δ ≤ 1 → 0 ≤ δ' → δ' ≤ 1 → δ < δ' →
    δ * (2 - δ) < δ' * (2 - δ') := by
  intros δ δ' h0 h1 h0' h1' hlt
  nlinarith

/-- **Monotone (non-strict) version.** On `[0, 1]`, `δ ≤ δ' ⇒ δ(2-δ) ≤ δ'(2-δ')`. -/
theorem delta_two_minus_delta_monotone :
    ∀ δ δ' : ℝ, 0 ≤ δ → δ ≤ 1 → 0 ≤ δ' → δ' ≤ 1 → δ ≤ δ' →
    δ * (2 - δ) ≤ δ' * (2 - δ') := by
  intros δ δ' h0 h1 h0' h1' hle
  rcases eq_or_lt_of_le hle with heq | hlt
  · subst heq; rfl
  · exact le_of_lt (delta_two_minus_delta_strictMono δ δ' h0 h1 h0' h1' hlt)

/-! ### Inversion: from `sin²θ` back to `δ` -/

/-- **Inversion identity.** Given a target `s ∈ [0, 1]`, the value
    `δ := 1 - √(1 - s)` lies in `[0, 1]` and satisfies `δ(2 - δ) = s`. -/
theorem delta_two_minus_delta_inversion (s : ℝ) (h0 : 0 ≤ s) (h1 : s ≤ 1) :
    let δ := 1 - Real.sqrt (1 - s)
    0 ≤ δ ∧ δ ≤ 1 ∧ δ * (2 - δ) = s := by
  set δ := 1 - Real.sqrt (1 - s) with hδdef
  have h1ms : 0 ≤ 1 - s := by linarith
  have hsqrt_nn : 0 ≤ Real.sqrt (1 - s) := Real.sqrt_nonneg _
  have hsqrt_le_one : Real.sqrt (1 - s) ≤ 1 := by
    have h_le : Real.sqrt (1 - s) ≤ Real.sqrt 1 := Real.sqrt_le_sqrt (by linarith)
    rwa [Real.sqrt_one] at h_le
  refine ⟨?_, ?_, ?_⟩
  · -- δ ≥ 0: 1 - sqrt ≥ 0 ⇔ sqrt ≤ 1
    linarith
  · -- δ ≤ 1: 1 - sqrt ≤ 1 ⇔ sqrt ≥ 0
    linarith
  · -- δ(2 - δ) = s
    -- 2 - δ = 1 + sqrt(1-s)
    have h2sub : 2 - δ = 1 + Real.sqrt (1 - s) := by
      show 2 - (1 - Real.sqrt (1 - s)) = 1 + Real.sqrt (1 - s); ring
    rw [hδdef, h2sub]
    -- (1 - sqrt)(1 + sqrt) = 1 - sqrt² = 1 - (1 - s) = s
    have hsq : Real.sqrt (1 - s) ^ 2 = 1 - s := Real.sq_sqrt h1ms
    nlinarith [hsq]

/-! ### Surjectivity of `sin²θ_p` on `[0, 1]` via cyclic phase -/

/-- **Surjectivity (cyclic phase form).** Every `s ∈ [0, 1]` is attained by
    `sin²θ` for some `θ` whose cosine has the PT form `cos θ = 1 - δ` with
    `δ ∈ [0, 1]`. -/
theorem cyclic_phase_surjective (s : ℝ) (h0 : 0 ≤ s) (h1 : s ≤ 1) :
    ∃ δ : ℝ, 0 ≤ δ ∧ δ ≤ 1 ∧ δ * (2 - δ) = s := by
  obtain ⟨h0δ, h1δ, hδ⟩ := delta_two_minus_delta_inversion s h0 h1
  exact ⟨1 - Real.sqrt (1 - s), h0δ, h1δ, hδ⟩

/-! ### Headline (Route 2 algebraic content) -/

/-- **Headline.** The map `δ ↦ δ(2 - δ)` is a strictly monotone bijection
    `[0, 1] ↔ [0, 1]`, with inverse `s ↦ 1 - √(1 - s)`. Combined with
    `CyclicPhaseIdentity.sin_sq_of_cos_eq_one_sub`, this gives the cyclic-phase
    isomorphism `δ_p ↔ sin²θ_p` on the PT-relevant interval. -/
theorem cyclic_phase_bijection_summary :
    -- monotone in `δ`
    (∀ δ δ' : ℝ, 0 ≤ δ → δ ≤ 1 → 0 ≤ δ' → δ' ≤ 1 → δ < δ' →
        δ * (2 - δ) < δ' * (2 - δ'))
    -- inversion via sqrt
    ∧ (∀ s : ℝ, 0 ≤ s → s ≤ 1 →
        let δ := 1 - Real.sqrt (1 - s)
        0 ≤ δ ∧ δ ≤ 1 ∧ δ * (2 - δ) = s) :=
  ⟨delta_two_minus_delta_strictMono, delta_two_minus_delta_inversion⟩

end PT.Holonomy
