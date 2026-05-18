/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Analysis.W7SpiralIdentity
import Mathlib.MeasureTheory.Integral.IntervalIntegral.FundThmCalculus
import Mathlib.Analysis.Calculus.Deriv.MeanValue
import Mathlib.Topology.Basic

/-!
# W7-1 — Reverse direction of the spiral identity

This module formalises the **reverse direction** of Theorem 6.3 of the
cuspidal Berry–Keating preprint (`PT_PROJECTS/PT_CH/paper_phase1/preprint1_cusp_fr.md`,
§6.6.4):

  *For `k ≥ 1`, `J σ k = J σ 0` implies `σ² = π(k+1)`.*

The case `k = 0` is degenerate: `J σ 0 = J σ 0` holds trivially for any
`σ`, so the identity gives no constraint; the conclusion `σ² = π · 1 = π`
is false in general. Hence the reverse direction is meaningful only for
`k ≥ 1`, and the statement is parameterised accordingly.

## Proof outline

1. Apply the centring substitution `v = x - σ²` (lemma `J_eq_centered`).
2. Use the parity of `gaussKer σ` to rewrite
   `J σ 0 = ∫ v in (σ² - 2π)..σ², gaussKer σ v` (`J_zero_reflected`).
3. Define the **window integral**
   `gaussWindow σ a := ∫ v in a..(a + 2π), gaussKer σ v`.
   Then `J σ k = gaussWindow σ (2πk - σ²)` and (reflected)
   `J σ 0 = gaussWindow σ (σ² - 2π)`.
4. Show `gaussWindow σ` is **strictly unimodal** with unique maximum at
   `a = -π`, via FTC + sign of derivative on each half-line + MVT-based
   `strictMonoOn_of_hasDerivWithinAt_pos` / `strictAntiOn_of_hasDerivWithinAt_neg`.
5. Combined with the reflection symmetry `F(σ, -2π - a) = F(σ, a)`, the
   level-set is `{a, -2π - a}`, so equality forces either `a₁ = a₂`
   (giving `σ² = π(k+1)`) or `a₁ + a₂ = -2π` (which collapses to `k = 0`).
-/

namespace PT.Analysis.W7

open Real intervalIntegral MeasureTheory Topology Filter

/-! ### Auxiliary: continuity & integrability of the centred Gaussian -/

/-- `gaussKer σ` is continuous on `ℝ` for any fixed `σ`. -/
lemma continuous_gaussKer (σ : ℝ) : Continuous (gaussKer σ) := by
  unfold gaussKer
  exact Real.continuous_exp.comp (by fun_prop)

/-- `gaussKer σ` is integrable on every bounded interval. -/
lemma intervalIntegrable_gaussKer (σ a b : ℝ) :
    IntervalIntegrable (gaussKer σ) MeasureTheory.volume a b :=
  (continuous_gaussKer σ).intervalIntegrable a b

/-! ### Reflection identity for `J σ 0` -/

/-- `J σ 0` written as an integral over the *reflected* window
`[σ² − 2π, σ²]`. This is the key reformulation that makes both turn
integrals comparable as windows of width `2π`. -/
lemma J_zero_reflected (σ : ℝ) :
    J σ 0 = ∫ v in (σ ^ 2 - 2 * π)..(σ ^ 2), gaussKer σ v := by
  rw [J_eq_centered σ 0]
  simp only [Nat.cast_zero, mul_zero, zero_sub, zero_add, mul_one]
  -- Goal: ∫ v in (-σ²)..(2π - σ²), gaussKer σ v
  --     = ∫ v in (σ² - 2π)..σ², gaussKer σ v
  have hpar : ∀ v : ℝ, gaussKer σ v = gaussKer σ (-v) :=
    fun v => (gaussKer_neg σ v).symm
  rw [show (∫ v in (-σ ^ 2)..(2 * π - σ ^ 2), gaussKer σ v)
        = (∫ v in (-σ ^ 2)..(2 * π - σ ^ 2), gaussKer σ (-v))
       from intervalIntegral.integral_congr (fun x _hx => hpar x)]
  rw [intervalIntegral.integral_comp_neg (gaussKer σ)]
  congr 1 <;> ring

/-! ### The window integral and its derivative -/

/-- The **window integral**
`gaussWindow σ a := ∫ v in a..(a + 2π), gaussKer σ v`.

This is the natural common form of `J σ k` (after centring) and
`J σ 0` (after centring + reflection). Both have width `2π`. -/
noncomputable def gaussWindow (σ a : ℝ) : ℝ :=
  ∫ v in a..(a + 2 * π), gaussKer σ v

/-- `J σ k` written as a window. -/
lemma J_eq_window (σ : ℝ) (k : ℕ) :
    J σ k = gaussWindow σ (2 * π * (k : ℝ) - σ ^ 2) := by
  rw [J_eq_centered σ k]
  unfold gaussWindow
  congr 1
  ring

/-- `J σ 0` written as a window with the reflected base point. -/
lemma J_zero_eq_window (σ : ℝ) :
    J σ 0 = gaussWindow σ (σ ^ 2 - 2 * π) := by
  rw [J_zero_reflected σ]
  unfold gaussWindow
  congr 1
  ring

/-! ### Sign of `gaussKer σ (a+2π) − gaussKer σ a` -/

/-- Sign of `gaussKer σ (a+2π) − gaussKer σ a`: strictly **positive** iff `a < -π`. -/
lemma gaussKer_diff_pos_iff (σ : ℝ) (hσ : 0 < σ) (a : ℝ) :
    0 < gaussKer σ (a + 2 * π) - gaussKer σ a ↔ a < -π := by
  unfold gaussKer
  rw [sub_pos, Real.exp_lt_exp]
  have hσsq : 0 < 4 * σ ^ 2 := by positivity
  rw [div_lt_div_iff_of_pos_right hσsq, neg_lt_neg_iff]
  constructor
  · intro h; nlinarith [Real.pi_pos]
  · intro h; nlinarith [Real.pi_pos]

/-- Sign of `gaussKer σ (a+2π) − gaussKer σ a`: strictly **negative** iff `-π < a`. -/
lemma gaussKer_diff_neg_iff (σ : ℝ) (hσ : 0 < σ) (a : ℝ) :
    gaussKer σ (a + 2 * π) - gaussKer σ a < 0 ↔ -π < a := by
  unfold gaussKer
  rw [sub_neg, Real.exp_lt_exp]
  have hσsq : 0 < 4 * σ ^ 2 := by positivity
  rw [div_lt_div_iff_of_pos_right hσsq, neg_lt_neg_iff]
  constructor
  · intro h; nlinarith [Real.pi_pos]
  · intro h; nlinarith [Real.pi_pos]

/-! ### Derivative of the window integral via FTC -/

/-- **FTC for the window integral.** `gaussWindow σ` has derivative
`gaussKer σ (a + 2π) − gaussKer σ a` at every point `a`. -/
lemma hasDerivAt_gaussWindow (σ a : ℝ) :
    HasDerivAt (gaussWindow σ)
      (gaussKer σ (a + 2 * π) - gaussKer σ a) a := by
  set g := gaussKer σ with hg_def
  -- Define the primitive G(u) = ∫ 0..u, g v dv.
  set G : ℝ → ℝ := fun u => ∫ v in (0 : ℝ)..u, g v with hG_def
  have hg_cont : Continuous g := continuous_gaussKer σ
  -- G has derivative g(u) at every point u.
  have hG_deriv : ∀ u : ℝ, HasDerivAt G (g u) u := by
    intro u
    have hint : IntervalIntegrable g MeasureTheory.volume 0 u :=
      (hg_cont.intervalIntegrable 0 u)
    have hcont_at : ContinuousAt g u := hg_cont.continuousAt
    have hsmeas : StronglyMeasurableAtFilter g (𝓝 u) MeasureTheory.volume :=
      hg_cont.aestronglyMeasurable.stronglyMeasurableAtFilter
    exact intervalIntegral.integral_hasDerivAt_right hint hsmeas hcont_at
  -- gaussWindow σ a = G(a + 2π) - G a.
  have hF_eq : ∀ a : ℝ, gaussWindow σ a = G (a + 2 * π) - G a := by
    intro a
    show (∫ v in a..(a + 2 * π), g v) = G (a + 2 * π) - G a
    have h1 : IntervalIntegrable g MeasureTheory.volume 0 a :=
      hg_cont.intervalIntegrable 0 a
    have h2 : IntervalIntegrable g MeasureTheory.volume a (a + 2 * π) :=
      hg_cont.intervalIntegrable a (a + 2 * π)
    have hadj := intervalIntegral.integral_add_adjacent_intervals h1 h2
    -- hadj : (∫ 0..a, g) + (∫ a..(a+2π), g) = ∫ 0..(a+2π), g
    show (∫ v in a..(a + 2 * π), g v) = (∫ v in (0:ℝ)..(a + 2 * π), g v) - (∫ v in (0:ℝ)..a, g v)
    linarith
  -- Derivative of u ↦ G(u + 2π) at a.
  have hcomp : HasDerivAt (fun a : ℝ => G (a + 2 * π)) (g (a + 2 * π)) a := by
    have h1 : HasDerivAt (fun a : ℝ => a + 2 * π) (1 : ℝ) a :=
      (hasDerivAt_id a).add_const _
    have h2 : HasDerivAt G (g (a + 2 * π)) (a + 2 * π) := hG_deriv _
    have hc := h2.comp a h1
    simpa using hc
  have hbase : HasDerivAt G (g a) a := hG_deriv a
  have hsub : HasDerivAt (fun a => G (a + 2 * π) - G a)
      (g (a + 2 * π) - g a) a := hcomp.sub hbase
  -- Transport via the equality.
  have heq : gaussWindow σ = fun a => G (a + 2 * π) - G a := funext hF_eq
  rw [heq]
  exact hsub

/-! ### Strict monotonicity on each half-line -/

/-- `gaussWindow σ` is **strictly increasing** on `(-∞, -π]`. -/
lemma gaussWindow_strictMonoOn_Iic (σ : ℝ) (hσ : 0 < σ) :
    StrictMonoOn (gaussWindow σ) (Set.Iic (-π)) := by
  apply strictMonoOn_of_hasDerivWithinAt_pos (convex_Iic _)
  · -- ContinuousOn on Iic (-π).
    intro x _
    exact (hasDerivAt_gaussWindow σ x).continuousAt.continuousWithinAt
  · -- HasDerivWithinAt on the interior.
    intro x _
    have h := hasDerivAt_gaussWindow σ x
    exact h.hasDerivWithinAt
  · -- Derivative is strictly positive on interior(Iic (-π)) = Iio (-π).
    intro x hx
    rw [interior_Iic] at hx
    exact (gaussKer_diff_pos_iff σ hσ x).mpr hx

/-- `gaussWindow σ` is **strictly decreasing** on `[-π, +∞)`. -/
lemma gaussWindow_strictAntiOn_Ici (σ : ℝ) (hσ : 0 < σ) :
    StrictAntiOn (gaussWindow σ) (Set.Ici (-π)) := by
  apply strictAntiOn_of_hasDerivWithinAt_neg (convex_Ici _)
  · intro x _
    exact (hasDerivAt_gaussWindow σ x).continuousAt.continuousWithinAt
  · intro x _
    exact (hasDerivAt_gaussWindow σ x).hasDerivWithinAt
  · intro x hx
    rw [interior_Ici] at hx
    exact (gaussKer_diff_neg_iff σ hσ x).mpr hx

/-! ### Reflection symmetry: `F(σ, −2π − a) = F(σ, a)` -/

/-- The window integral is symmetric about `a = -π`:
`gaussWindow σ (-2π - a) = gaussWindow σ a`. -/
lemma gaussWindow_reflect (σ a : ℝ) :
    gaussWindow σ (-2 * π - a) = gaussWindow σ a := by
  unfold gaussWindow
  have hpar : ∀ v : ℝ, gaussKer σ v = gaussKer σ (-v) :=
    fun v => (gaussKer_neg σ v).symm
  -- LHS bounds: lo = -2π - a, hi = lo + 2π = -a.
  rw [show (∫ v in (-2 * π - a)..(-2 * π - a + 2 * π), gaussKer σ v)
        = (∫ v in (-2 * π - a)..(-2 * π - a + 2 * π), gaussKer σ (-v))
       from intervalIntegral.integral_congr (fun x _hx => hpar x)]
  rw [intervalIntegral.integral_comp_neg (gaussKer σ)]
  -- Bounds:  -(-2π - a + 2π) = a   and   -(-2π - a) = 2π + a = a + 2π.
  congr 1 <;> ring

/-! ### Main level-set characterisation -/

/-- **Level-set lemma.** For `σ > 0` and any `a₁, a₂`,
`gaussWindow σ a₁ = gaussWindow σ a₂` iff `a₁ = a₂` or `a₁ + a₂ = -2π`. -/
lemma gaussWindow_eq_iff (σ : ℝ) (hσ : 0 < σ) (a₁ a₂ : ℝ) :
    gaussWindow σ a₁ = gaussWindow σ a₂ ↔ a₁ = a₂ ∨ a₁ + a₂ = -2 * π := by
  refine ⟨fun h => ?_, ?_⟩
  · -- (⇒) Forward: distinguish cases on a₁ ≤ -π vs > -π and similarly for a₂.
    -- WLOG via symmetry under (a₁ ↔ a₂); handled by case analysis directly.
    rcases lt_trichotomy a₁ a₂ with hlt | heq | hgt
    · -- a₁ < a₂.
      right
      -- Use reflection: F(a₁) = F(-2π - a₁) = F(a₂). Then if both -2π - a₁ and a₂
      -- lie in a region of strict monotonicity, injectivity gives -2π - a₁ = a₂.
      by_cases hb1 : a₁ ≤ -π
      · -- a₁ ≤ -π, so a₁' := -2π - a₁ ≥ -π.
        have ha1' : -π ≤ -2 * π - a₁ := by linarith
        have hreflect : gaussWindow σ (-2 * π - a₁) = gaussWindow σ a₁ :=
          gaussWindow_reflect σ a₁
        have hFeq : gaussWindow σ (-2 * π - a₁) = gaussWindow σ a₂ := by
          rw [hreflect]; exact h
        by_cases hb2 : -π ≤ a₂
        · -- both -2π - a₁ and a₂ in Ici (-π) → injectivity ⇒ equality.
          have hInj : Set.InjOn (gaussWindow σ) (Set.Ici (-π)) :=
            (gaussWindow_strictAntiOn_Ici σ hσ).injOn
          have := hInj ha1' hb2 hFeq
          linarith
        · -- a₂ < -π. Then a₁ < a₂ < -π, use strict mono on Iic (-π).
          simp only [not_le] at hb2
          have hInj : Set.InjOn (gaussWindow σ) (Set.Iic (-π)) :=
            (gaussWindow_strictMonoOn_Iic σ hσ).injOn
          have ha1_in : a₁ ∈ Set.Iic (-π) := hb1
          have ha2_in : a₂ ∈ Set.Iic (-π) := le_of_lt hb2
          have := hInj ha1_in ha2_in h
          exact absurd this (ne_of_lt hlt)
      · -- a₁ > -π.
        simp only [not_le] at hb1
        have ha1_in : a₁ ∈ Set.Ici (-π) := le_of_lt hb1
        by_cases hb2 : -π ≤ a₂
        · -- both in Ici (-π): injectivity ⇒ a₁ = a₂. Contradicts a₁ < a₂.
          have hInj : Set.InjOn (gaussWindow σ) (Set.Ici (-π)) :=
            (gaussWindow_strictAntiOn_Ici σ hσ).injOn
          exact absurd (hInj ha1_in hb2 h) (ne_of_lt hlt)
        · simp only [not_le] at hb2
          -- a₁ > -π and a₂ < -π, but a₁ < a₂ ⇒ a₁ < -π. Contradiction.
          linarith
    · -- a₁ = a₂.
      left; exact heq
    · -- a₂ < a₁. Symmetric: invoke the previous branch with roles swapped.
      right
      have hsym : gaussWindow σ a₂ = gaussWindow σ a₁ := h.symm
      by_cases hb2 : a₂ ≤ -π
      · have ha2' : -π ≤ -2 * π - a₂ := by linarith
        have hreflect : gaussWindow σ (-2 * π - a₂) = gaussWindow σ a₂ :=
          gaussWindow_reflect σ a₂
        have hFeq : gaussWindow σ (-2 * π - a₂) = gaussWindow σ a₁ := by
          rw [hreflect]; exact hsym
        by_cases hb1 : -π ≤ a₁
        · have hInj : Set.InjOn (gaussWindow σ) (Set.Ici (-π)) :=
            (gaussWindow_strictAntiOn_Ici σ hσ).injOn
          have := hInj ha2' hb1 hFeq
          linarith
        · simp only [not_le] at hb1
          have hInj : Set.InjOn (gaussWindow σ) (Set.Iic (-π)) :=
            (gaussWindow_strictMonoOn_Iic σ hσ).injOn
          have ha2_in : a₂ ∈ Set.Iic (-π) := hb2
          have ha1_in : a₁ ∈ Set.Iic (-π) := le_of_lt hb1
          exact absurd (hInj ha2_in ha1_in hsym) (ne_of_lt hgt)
      · simp only [not_le] at hb2
        have ha2_in : a₂ ∈ Set.Ici (-π) := le_of_lt hb2
        by_cases hb1 : -π ≤ a₁
        · have hInj : Set.InjOn (gaussWindow σ) (Set.Ici (-π)) :=
            (gaussWindow_strictAntiOn_Ici σ hσ).injOn
          exact absurd (hInj ha2_in hb1 hsym) (ne_of_lt hgt)
        · simp only [not_le] at hb1
          linarith
  · -- (⇐) Reverse: trivial from reflection symmetry.
    rintro (rfl | hsum)
    · rfl
    · -- a₁ + a₂ = -2π ⇒ a₂ = -2π - a₁
      have h2 : a₂ = -2 * π - a₁ := by linarith
      rw [h2, gaussWindow_reflect]

/-! ### Main theorem (reverse direction) -/

/-- **Theorem W7-1 (reverse direction).** For `k ≥ 1` and `σ > 0`, the
spiral identity `J σ k = J σ 0` forces `σ² = π(k+1)`.

For `k = 0` the hypothesis `J σ 0 = J σ 0` is trivial and the conclusion
`σ² = π` is false in general; hence the hypothesis `1 ≤ k`. -/
theorem W7_1_spiral_identity_reverse (k : ℕ) (hk : 1 ≤ k) (σ : ℝ) (hσ : 0 < σ)
    (h_eq : J σ k = J σ 0) :
    σ ^ 2 = π * ((k : ℝ) + 1) := by
  -- Both sides as windows.
  rw [J_eq_window σ k, J_zero_eq_window σ] at h_eq
  -- Level-set lemma.
  have hcases := (gaussWindow_eq_iff σ hσ
    (2 * π * (k : ℝ) - σ ^ 2) (σ ^ 2 - 2 * π)).mp h_eq
  rcases hcases with heq | hsum
  · -- Case (a): 2πk - σ² = σ² - 2π ⇒ σ² = π(k+1).
    linarith
  · -- Case (b): (2πk - σ²) + (σ² - 2π) = -2π ⇒ 2πk = 0 ⇒ k = 0.
    have h0 : 2 * π * (k : ℝ) = 0 := by linarith
    have hπpos : (0 : ℝ) < π := Real.pi_pos
    have h2π_pos : (0 : ℝ) < 2 * π := by linarith
    have hk_zero : (k : ℝ) = 0 := by
      have := mul_eq_zero.mp h0
      rcases this with h | h
      · exfalso; linarith
      · exact h
    have hk0 : k = 0 := by exact_mod_cast hk_zero
    omega

end PT.Analysis.W7
