/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import Mathlib.MeasureTheory.Integral.IntervalIntegral.Basic
import Mathlib.Analysis.SpecialFunctions.Exp
import Mathlib.Analysis.SpecialFunctions.Pow.Real

/-!
# W7-1 — The Spiral Identity of the Weil Explicit Formula

This module formalises Theorem 6.3 of the cuspidal Berry-Keating preprint
(`PT_PROJECTS/PT_CH/paper_phase1/preprint1_cusp_fr.md`, §6.6.3), known as
the **spiral identity** of Persistence Theory.

## Statement

Define the **shifted Gaussian** centred at `σ²` and width `σ` by
`gaussKer σ y := exp(-y²/(4σ²))`, and the integral over a turn of the
log-polar spiral of Archimedes by

  `J σ k := ∫ x in (2π k)..(2π (k+1)), gaussKer σ (x - σ²)`.

This is, up to an `e^{σ²/4}` prefactor and the discrete-vs-continuum
correction, the contribution of the `k`-th turn of the log-polar spiral
to the prime sum `T_G` in the Weil explicit formula (Iwaniec-Kowalski
2004 §5.5).

**Theorem W7-1 (continuum form, this file).** *Under the substitution
`u = log p`, dp/log p ↦ e^u du/u, and the PNT density approximation, the
contributions of distinct turns of the spiral satisfy*

  `σ² = π(k+1) → J σ k = J σ 0.`

In particular, the cascade `σ_crit^(k) = √(π(k+1))` defines a sequence
of universal seuils intrinsic to the period `2π` of the spiral.

## What is formalised here

- **`gaussKer`**: the centred Gaussian kernel.
- **`gaussKer_neg`**: parity (`gaussKer σ (-y) = gaussKer σ y`).
- **`J`**: the spiral-turn integral.
- **`sigmaCrit`**: the cascade `√(π(k+1))`.
- **`sigmaCrit_sq`**: `(σ_crit k)² = π(k+1)`.
- **`W7_1_spiral_identity`**: the **forward direction** of Theorem 6.3
  — `σ² = π(k+1)` implies `J σ k = J σ 0`.

The reverse direction requires a monotonicity argument and is left as a
future-work TODO; it is not on the PT critical path.

## Triple numerical validation

The continuum identity is empirically validated to <1 % at three orders:

| `k` | `σ_crit^(k) = √(π(k+1))` | observed | relative error | `#primes` in turn `k` |
|---|---|---|---|---|
| 1 | √(2π) ≈ 2.5066 | 2.4909 | −0.626 % | 24 877 |
| 2 | √(3π) ≈ 3.0700 | 3.0630 | −0.228 % | ≈ 8.6 × 10⁶ |
| 3 | √(4π) ≈ 3.5449 | 3.5419 | −0.085 % | ≈ 3.4 × 10⁹ |

(Validation k=3 via `primesieve` CLI on turn 3 ⊂ [e^{6π}, e^{8π}],
83 chunks of range 10⁹, 12.4 min CPU.)

## References

- Y. Senez, *L'opérateur de Berry-Keating sur un cusp hyperbolique*
  (`PT_PROJECTS/PT_CH/paper_phase1/preprint1_cusp_fr.md`), §6.6.
- Y. Senez, *Modèle spectral canonique des zéros de ζ via la dynamique
  de persistance PT* (`PT_New_Math_Consolidation_FR/article_PT_zeta_FR.md`),
  §7.7.
- Session log `PT_RH_WEIL_SPIRAL/analysis/HP_WEIL_POSITIVITY_SPIRAL.md`.
- Weil, A., 1952, *Sur les "formules explicites" de la théorie des
  nombres premiers*, Comm. Sém. Math. Univ. Lund.
- Iwaniec, H., Kowalski, E., 2004, *Analytic Number Theory*, AMS Coll. 53.
-/

namespace PT.Analysis.W7

open Real intervalIntegral MeasureTheory

/-! ### Definitions -/

/-- The centred Gaussian kernel `g(y) = exp(-y² / (4 σ²))`. Even in `y`. -/
noncomputable def gaussKer (σ y : ℝ) : ℝ := Real.exp (-y ^ 2 / (4 * σ ^ 2))

/-- The shifted Gaussian `gauss(σ, x) = gaussKer(σ, x − σ²)`. This is the
integrand in `J σ k` and corresponds to the natural integrand of the
Weil prime sum after substitution `u = log p` and completion of the
square (cf. preprint cusp §6.6.3, formula (6.12)). -/
noncomputable def gaussShifted (σ x : ℝ) : ℝ := gaussKer σ (x - σ ^ 2)

/-- The spiral-turn integral
`J σ k = ∫ x in (2π k)..(2π (k+1)), gaussShifted σ x`. -/
noncomputable def J (σ : ℝ) (k : ℕ) : ℝ :=
  ∫ x in (2 * π * (k : ℝ))..(2 * π * ((k : ℝ) + 1)), gaussShifted σ x

/-- The cascade of critical scales: `σ_crit^(k) = √(π (k+1))`. -/
noncomputable def sigmaCrit (k : ℕ) : ℝ := Real.sqrt (π * ((k : ℝ) + 1))

/-! ### Elementary lemmas -/

/-- Parity of the centred Gaussian kernel: `gaussKer σ (-y) = gaussKer σ y`. -/
lemma gaussKer_neg (σ y : ℝ) : gaussKer σ (-y) = gaussKer σ y := by
  unfold gaussKer
  congr 1
  ring

/-- `(σ_crit k)² = π(k+1)`. -/
lemma sigmaCrit_sq (k : ℕ) : (sigmaCrit k) ^ 2 = π * ((k : ℝ) + 1) := by
  unfold sigmaCrit
  rw [sq, Real.mul_self_sqrt]
  exact mul_nonneg Real.pi_nonneg (by positivity)

/-- `σ_crit k ≥ 0`. -/
lemma sigmaCrit_nonneg (k : ℕ) : 0 ≤ sigmaCrit k := Real.sqrt_nonneg _

/-- `σ_crit k > 0`. -/
lemma sigmaCrit_pos (k : ℕ) : 0 < sigmaCrit k := by
  unfold sigmaCrit
  apply Real.sqrt_pos.mpr
  have : (0 : ℝ) < π := Real.pi_pos
  positivity

/-! ### Centering substitution -/

/-- **Centering substitution.** Under `v = x − σ²`,
`J σ k = ∫ v in (2π k − σ²)..(2π(k+1) − σ²), gaussKer σ v`. -/
lemma J_eq_centered (σ : ℝ) (k : ℕ) :
    J σ k =
    ∫ v in (2 * π * (k : ℝ) - σ ^ 2)..(2 * π * ((k : ℝ) + 1) - σ ^ 2), gaussKer σ v := by
  unfold J gaussShifted
  exact intervalIntegral.integral_comp_sub_right (gaussKer σ) (σ ^ 2)

/-! ### Theorem W7-1 (forward direction) -/

/-- **Theorem W7-1 (forward).** If `σ² = π(k+1)`, then `J σ k = J σ 0`.

*Proof sketch.* After centring (lemma `J_eq_centered`), both integrals
become integrals of the even function `gaussKer σ` over intervals which
are symmetric about `0` when `σ² = π(k+1)`:

* `J σ k` becomes `∫ v in π(k−1)..π(k+1), gaussKer σ v`,
* `J σ 0` becomes `∫ v in −π(k+1)..−π(k−1), gaussKer σ v`.

By parity of `gaussKer σ` (lemma `gaussKer_neg`) and
`intervalIntegral.integral_comp_neg`, the two integrals are equal. □ -/
theorem W7_1_spiral_identity (k : ℕ) (σ : ℝ) (_hσ : 0 < σ)
    (hkσ : σ ^ 2 = π * ((k : ℝ) + 1)) :
    J σ k = J σ 0 := by
  -- Step 1: centre both integrals via the substitution v = x - σ².
  rw [J_eq_centered σ k, J_eq_centered σ 0, hkσ]
  -- Goal:
  --   ∫ v in (2π·k - π(k+1))..(2π(k+1) - π(k+1)), gaussKer σ v
  --   = ∫ v in (2π·↑0 - π(k+1))..(2π(↑0+1) - π(k+1)), gaussKer σ v
  -- Normalize the Nat.cast 0 in the RHS bounds.
  simp only [Nat.cast_zero, mul_zero, zero_add, mul_one] at *
  -- Goal now:
  --   ∫ v in (2π·k - π(k+1))..(2π(k+1) - π(k+1)), gaussKer σ v
  --   = ∫ v in (0 - π(k+1))..(2π - π(k+1)), gaussKer σ v
  -- Step 2: Transform the RHS via parity then `integral_comp_neg`.
  symm
  -- Replace the integrand on the LHS using parity.
  have hpar : ∀ v : ℝ, gaussKer σ v = gaussKer σ (-v) :=
    fun v => (gaussKer_neg σ v).symm
  rw [show (∫ v in (0 - π * ((k : ℝ) + 1))..(2 * π - π * ((k : ℝ) + 1)), gaussKer σ v)
        = (∫ v in (0 - π * ((k : ℝ) + 1))..(2 * π - π * ((k : ℝ) + 1)), gaussKer σ (-v))
       from intervalIntegral.integral_congr (fun x _hx => hpar x)]
  -- Apply integral_comp_neg: ∫ v in a..b, f(-v) = ∫ v in -b..-a, f v
  rw [intervalIntegral.integral_comp_neg (gaussKer σ)]
  -- Bounds match by ring algebra:
  --   -(2π - π(k+1)) = π(k+1) - 2π = π·k - π = 2π·k - π(k+1)
  --   -(0 - π(k+1)) = π(k+1) = 2π(k+1) - π(k+1)
  congr 1 <;> ring

/-- **Corollary.** At the cascade critical scale `σ = σ_crit k = √(π(k+1))`,
the spiral identity holds: `J (σ_crit k) k = J (σ_crit k) 0`. -/
theorem J_at_sigmaCrit_eq_J_zero (k : ℕ) :
    J (sigmaCrit k) k = J (sigmaCrit k) 0 := by
  apply W7_1_spiral_identity
  · exact sigmaCrit_pos k
  · exact sigmaCrit_sq k

/-! ### Future work

The reverse direction of the spiral identity (`J σ k = J σ 0 ⟹ σ² = π(k+1)`)
requires showing strict monotonicity of `σ ↦ J σ k / J σ 0` for fixed
`k ≥ 1`. This is true but needs an explicit derivative argument and is
left as a TODO. It is not on the PT critical path: the forward direction
suffices for all downstream applications, in particular for identifying
the cascade `σ_crit^(k)` as a sequence of natural seuils.

The discrete-to-continuum estimate `|σ_crit_obs^(k) − √(π(k+1))| =
O(1/log(#primes))` is not formalised here; it is an empirical statement
relying on quantitative versions of the Prime Number Theorem (Vinogradov-
Korobov bounds), which are not yet in Mathlib in the required form. -/

end PT.Analysis.W7
