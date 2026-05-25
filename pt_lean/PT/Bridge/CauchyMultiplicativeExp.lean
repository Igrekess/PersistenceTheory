/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.Analysis.SpecialFunctions.Exp
import Mathlib.Topology.Instances.RealVectorSpace
import Mathlib.Tactic

/-!
# Multiplicative Cauchy functional equation : `f(x+y) = f(x)·f(y)` ⇒ `f = exp(c·x)`

This module formalises the **mathematical core** of the K4 closure
chain (Higgs self-coupling `λ_H = 1/8`) in PT.

## Statement

If `f : ℝ → ℝ` is **continuous**, **strictly positive**, and satisfies the
**additive-to-multiplicative Cauchy equation**
`f(x + y) = f(x) · f(y)` with `f(0) = 1`, then
`f(x) = exp(c · x)` for a unique constant `c = log(f(1))`.

This is the classical theorem of Aczél (*Lectures on Functional
Equations and Their Applications* §2.1.2). It is the **last classical
mathematical step** in the chain that closes the PT Higgs self-coupling
`λ_H = 1/8` (K4).

## PT context

The relevance of this functional equation in PT is the following:

* CRT factorisation of the PT Dirac operator
  `D_PT² = ⊕_p D_p²` (proven PT theorem, `PT/Stochastic/CRTFactorizationT2.lean`)
* + Shore–Johnson G1 axiom (system independence, used throughout
  `PT/Information/`)
* ⟹ the spectral cutoff `f` of the Chamseddine–Connes spectral action
  must satisfy `f(u₁+u₂)·f(0) = f(u₁)·f(u₂)` to preserve CRT factorisation
  (otherwise it introduces cross-channel correlations that contradict G1).
* By the present theorem `cauchy_mult_exp`, this forces
  `f(u) = exp(c·u)` uniquely.
* Specialised by (P1) decay (`c < 0`), (P3) normalisation `f(0) = 1`,
  (P4) mean = `N_b` (cardinality of the `q_+/q_-` bifurcation, equal to
  `2`), one concludes `f_PT(u) = exp(-u/N_b) = exp(-u/2)`.
* This is the cutoff function that, plugged into the algebraic closure
  formula `(f₀/f₂)·(v/Λ)² = 21/100`, yields `λ_H = 1/8` (K4).

## Architecture of the proof

The proof here is a 50-line wrapper around the existing PT_LEAN
infrastructure:

1. Define `h(x) := log(f(x))`, well-defined and continuous since `f > 0`.
2. Show `h` is additive: `h(x+y) = h(x) + h(y)` (follows from
   `f(x+y) = f(x)·f(y)` and `Real.log_mul`).
3. Package `h` as `AddMonoidHom ℝ →+ ℝ`.
4. Apply Mathlib's `map_real_smul` (continuous additive map between
   real vector spaces is `ℝ`-linear) to deduce `h(x) = x · h(1)`.
5. Conclude `f(x) = exp(h(x)) = exp(h(1) · x)`.

This mirrors the pattern of `PT.Information.cauchy_log_equation`
(`PT/Information/T6cChencov.lean` lines 512-552), which proved the
**dual** statement (multiplicative-input Cauchy
`g(xy) = g(x) + g(y)` ⇒ `g = k · log`).

## References

* Aczél, *Lectures on Functional Equations and Their Applications*
  (Academic Press, 1966), Chapter 2, §2.1.2.
* PT monograph: appendix Y §Y.13.6 (`appendices_fr/app_y_higgs_zeta_duality.tex`),
  `thm:E8_formal_uniqueness` and `thm:E8_classification`.
-/

namespace PT.Bridge

open Real

/-! ## The main theorem -/

/-- **Multiplicative Cauchy functional equation on ℝ.**
    If `f : ℝ → ℝ` is continuous, strictly positive, and satisfies
    `f(x+y) = f(x) · f(y)` with `f(0) = 1`, then there exists `c ∈ ℝ` such
    that `f(x) = exp(c · x)` for every `x`. The constant `c` equals
    `log(f(1))`.

    **Proof strategy.** Substitute `h(x) := log(f(x))`. The hypotheses imply
    `h` is continuous, additive (`h(x+y) = h(x) + h(y)`) and satisfies
    `h(0) = 0`. The Mathlib lemma `map_real_smul` (continuous additive
    endomorphism of `ℝ` is `ℝ`-linear) then yields `h(x) = x · h(1)`. The
    inverse `f = exp ∘ h` finishes the proof. -/
theorem cauchy_mult_exp
    (f : ℝ → ℝ) (hf_cont : Continuous f)
    (hf_pos : ∀ x, 0 < f x)
    (hf_mult : ∀ x y, f (x + y) = f x * f y) :
    ∃ c : ℝ, ∀ x, f x = Real.exp (c * x) := by
  -- Define `h(x) := log(f(x))`.
  set h : ℝ → ℝ := fun x => Real.log (f x) with h_def
  -- `f(0) = 1` follows from `f(0+0) = f(0)·f(0)` and positivity.
  have hf_zero : f 0 = 1 := by
    have key : f 0 = f 0 * f 0 := by
      have := hf_mult 0 0
      simpa using this
    -- `a = a*a` with `a > 0` ⟹ `a = 1`.
    have hf0_pos : 0 < f 0 := hf_pos 0
    have hf0_ne : f 0 ≠ 0 := ne_of_gt hf0_pos
    -- From `f 0 = f 0 * f 0`, divide both sides by `f 0`.
    have : 1 = f 0 := by
      have := key
      field_simp at this
      linarith
    linarith
  -- `h` is continuous: composition of continuous `log` (on positives) and `f`.
  have h_cont : Continuous h := by
    refine Continuous.log hf_cont ?_
    intro x
    exact ne_of_gt (hf_pos x)
  -- `h(0) = log(f(0)) = log 1 = 0`.
  have h_zero : h 0 = 0 := by
    show Real.log (f 0) = 0
    rw [hf_zero, Real.log_one]
  -- `h(x+y) = log(f(x)·f(y)) = log(f x) + log(f y) = h x + h y`.
  have h_add : ∀ x y : ℝ, h (x + y) = h x + h y := by
    intro x y
    show Real.log (f (x + y)) = Real.log (f x) + Real.log (f y)
    rw [hf_mult x y, Real.log_mul (ne_of_gt (hf_pos x)) (ne_of_gt (hf_pos y))]
  -- Package `h` as an `AddMonoidHom`.
  let hHom : ℝ →+ ℝ :=
    { toFun := h, map_zero' := h_zero, map_add' := h_add }
  have hHom_cont : Continuous hHom := h_cont
  -- Continuous additive endomorphism of `ℝ` is `ℝ`-linear.
  have h_linear : ∀ x : ℝ, h x = x * h 1 := by
    intro x
    have := map_real_smul hHom hHom_cont x (1 : ℝ)
    simp only [smul_eq_mul, mul_one] at this
    exact this
  -- Conclude with `c := h 1 = log(f 1)`.
  refine ⟨h 1, fun x => ?_⟩
  -- `f x = exp(log(f x)) = exp(h x) = exp(x · h 1) = exp(h 1 · x)`.
  have step1 : f x = Real.exp (h x) := (Real.exp_log (hf_pos x)).symm
  rw [step1, h_linear x, mul_comm]

/-! ## Uniqueness corollary -/

/-- **Uniqueness of the constant `c`.** If `f(x) = exp(c₁·x) = exp(c₂·x)`
    for every `x`, then `c₁ = c₂`. Follows by evaluating at `x = 1` and
    taking `log`. -/
theorem cauchy_mult_exp_unique
    (c₁ c₂ : ℝ)
    (h_eq : ∀ x : ℝ, Real.exp (c₁ * x) = Real.exp (c₂ * x)) :
    c₁ = c₂ := by
  have key : Real.exp c₁ = Real.exp c₂ := by
    have := h_eq 1
    simpa using this
  exact Real.exp_injective key

/-! ## Specialisation to decaying cutoffs (PT case) -/

/-- **PT specialisation.** Under the additional assumption that `f` is
    **decreasing on `[0, ∞)`** (the PT regulator condition (P1)), the
    constant `c` is strictly negative, so `f(x) = exp(-x/τ)` for some
    `τ > 0` (with `τ = -1/c`). -/
theorem cauchy_mult_exp_decay
    (f : ℝ → ℝ) (hf_cont : Continuous f)
    (hf_pos : ∀ x, 0 < f x)
    (hf_mult : ∀ x y, f (x + y) = f x * f y)
    (hf_decr : f 1 < f 0) :
    ∃ τ : ℝ, 0 < τ ∧ ∀ x, f x = Real.exp (-x / τ) := by
  obtain ⟨c, hc⟩ := cauchy_mult_exp f hf_cont hf_pos hf_mult
  -- `c < 0`: from `exp(c·1) = f(1) < f(0) = exp(0) = 1`.
  have hc_neg : c < 0 := by
    have hf1 : f 1 = Real.exp c := by simpa using hc 1
    have hf0 : f 0 = 1 := by
      have key : f 0 = Real.exp (c * 0) := hc 0
      simp at key
      exact key
    rw [hf0, hf1] at hf_decr
    -- `exp c < 1 = exp 0` ⟹ `c < 0`.
    have := Real.exp_lt_exp.mp (by simpa using hf_decr)
    exact this
  have h_pos : (0 : ℝ) < -1 / c := by
    rw [neg_div, neg_pos]
    exact div_neg_of_pos_of_neg one_pos hc_neg
  refine ⟨-1/c, h_pos, fun x => ?_⟩
  rw [hc x]
  congr 1
  field_simp

/-! ## PT-specific cutoff: `f(u) = exp(-u/N_b)`

The PT spectral cutoff arises as the case `c = -1/N_b` of
`cauchy_mult_exp_decay`, where `N_b = 2` is the cardinality of the
`q_+/q_-` bifurcation. This wraps the abstract theorem into the PT-named
identifier used in the monograph annexe Y. -/

/-- The number of branches of the PT bifurcation `q_+ / q_-`. By
    convention this equals `2`: there are exactly two thermal/structural
    branches whose difference `Δq = q_- - q_+` defines the Higgs field
    via the off-diagonal Dirac `D_F = Δq · σ_x`. -/
noncomputable def N_b : ℝ := 2

@[simp]
theorem N_b_eq : N_b = 2 := rfl

theorem N_b_pos : 0 < N_b := by
  unfold N_b; norm_num

/-- **PT spectral cutoff (paper-level).** The unique cutoff function
    `f : ℝ → ℝ` that is continuous, strictly positive, multiplicative
    (`f(x+y) = f(x)·f(y)`), decreasing, and has *mean* `N_b` (equivalently
    `f'(0) / f(0) = -1/N_b` for the corresponding probability density) is

    $$ f_{\mathrm{PT}}(u) = \exp(-u/N_b) = \exp(-u/2). $$

    The "mean = `N_b`" condition is the structural PT input (P4): it
    reflects the cardinality of the bifurcation, the unique natural
    scale of the finite spectral triple `ST_F = (ℂ², ℂ², m·σ_x, σ_x)`. -/
noncomputable def cutoffPT : ℝ → ℝ := fun u => Real.exp (-u / N_b)

theorem cutoffPT_pos : ∀ u, 0 < cutoffPT u := by
  intro u; unfold cutoffPT; exact Real.exp_pos _

theorem cutoffPT_mult : ∀ u v, cutoffPT (u + v) = cutoffPT u * cutoffPT v := by
  intro u v
  unfold cutoffPT
  rw [← Real.exp_add]
  congr 1
  ring

theorem cutoffPT_zero : cutoffPT 0 = 1 := by
  unfold cutoffPT
  simp [Real.exp_zero]

theorem cutoffPT_continuous : Continuous cutoffPT := by
  unfold cutoffPT
  fun_prop

/-! ## Moments `f₀, f₂` of the PT cutoff

The moments `f₀ := f(0)` and `f₂ := ∫₀^∞ u · f(u) du` of the PT cutoff
play a key role in the algebraic closure of K4. For
`f_PT(u) = exp(-u/N_b)`, these are:

* `f₀ = 1`
* `f₂ = N_b²`
* hence `f₀ / f₂ = 1/N_b² = 1/4`

which is precisely the prefactor required by the E6c closure
`(f₀/f₂)·(v/Λ)² = 21/100 = (p₁·p₃)/(N_b²·p₂²)`.

The lemma `cutoffPT_f0_div_f2` records the algebraic identity at the
finite level. The integral statement is left to a future formalisation
that imports `Mathlib.MeasureTheory.Integral`. -/

/-- **Boundary moment of the PT cutoff.** `f₀ := f_PT(0) = 1`. -/
theorem cutoffPT_f0 : cutoffPT 0 = 1 := cutoffPT_zero

/-- **Algebraic ratio `f₀ / f₂` for the PT cutoff.** The classical
    moments of `exp(-u/τ)` over `[0, ∞)` are `f₀ = 1`, `f₂ = τ²`. For
    `τ = N_b`, we have `f₀ / f₂ = 1 / N_b² = 1/4`. The integral computation
    is standard (`∫₀^∞ u · exp(-u/τ) du = τ²`) and is recorded at the
    algebraic level here. -/
theorem cutoffPT_f0_div_f2 :
    (1 : ℝ) / (N_b * N_b) = 1 / 4 := by
  unfold N_b; norm_num

end PT.Bridge
