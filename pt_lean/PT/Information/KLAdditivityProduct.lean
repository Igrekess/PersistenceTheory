/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Information.GFTIdentity
import PT.Information.BekensteinExtensions
import PT.Information.GFTSpecialisations
import Mathlib.Tactic

/-!
# KL divergence — Additivity at the delta corner (Ch04 extension)

For the standard KL divergence to the uniform, the **tensor product**
behavior is a key compatibility check:

  `D_KL(δ_(r₀, s₀) ‖ U_{m·n}) = log(m · n) = log m + log n
                              = D_KL(δ_{r₀} ‖ U_m) + D_KL(δ_{s₀} ‖ U_n).`

This file proves this **additivity at the delta corner** explicitly,
restricted to the cases where one or both distributions are deltas (where
the equality reduces to the additive identity `log(mn) = log m + log n`).

The full additivity for arbitrary `P ⊗ Q` requires a non-trivial product
construction on `α × β`-distributions, which we leave to a future module.

## Reference

Monograph Chapter 4 §"Additivité KL", `\label{prop:KL_additivity}`.
-/

namespace PT.Information

open Real Finset

/-! ### Logarithm additivity (the algebraic core) -/

/-- **Additivity of `log` over multiplication.** `log(m · n) = log m + log n`
    for `m, n > 0`. This is the algebraic identity behind KL additivity. -/
theorem log_mul_add (m n : ℝ) (hm : 0 < m) (hn : 0 < n) :
    Real.log (m * n) = Real.log m + Real.log n :=
  Real.log_mul (ne_of_gt hm) (ne_of_gt hn)

/-! ### Delta-corner KL additivity -/

/-- **Delta-corner KL additivity.**
    `D_KL(δ_{r₀} ‖ U_m) + D_KL(δ_{s₀} ‖ U_n) = log(m · n)`. -/
theorem klDelta_add_klDelta_eq_log_mul
    {α β : Type*} [DecidableEq α] [DecidableEq β]
    (s : Finset α) (t : Finset β) (m n : ℝ)
    (hm : 0 < m) (hn : 0 < n)
    (r₀ : α) (hr : r₀ ∈ s) (s₀ : β) (hs : s₀ ∈ t) :
    klToUniform s m (deltaDist r₀) + klToUniform t n (deltaDist s₀)
      = Real.log (m * n) := by
  rw [klToUniform_delta s m hm r₀ hr]
  rw [klToUniform_delta t n hn s₀ hs]
  rw [log_mul_add m n hm hn]

/-! ### Specialised: equal-modulus case `m = n` -/

/-- **Equal-modulus delta additivity.** When `m = n`,
    `D_KL(δ ‖ U_m) + D_KL(δ ‖ U_m) = 2 · log m = log(m²)`. -/
theorem klDelta_double_eq_log_sq
    {α : Type*} [DecidableEq α]
    (s : Finset α) (m : ℝ) (hm : 0 < m)
    (r₀ : α) (hr : r₀ ∈ s) :
    2 * klToUniform s m (deltaDist r₀) = Real.log (m * m) := by
  rw [klToUniform_delta s m hm r₀ hr]
  rw [Real.log_mul (ne_of_gt hm) (ne_of_gt hm)]
  ring

/-! ### Specialised at PT moduli `m, n ∈ {2, 8, 30}` -/

/-- At `(m, n) = (2, 8)`: `D_KL(δ ‖ U_2) + D_KL(δ ‖ U_8) = log 16`. -/
theorem klDelta_2_plus_klDelta_8 :
    klToUniform (Finset.univ : Finset (Fin 2)) 2 (deltaDist 0)
      + klToUniform (Finset.univ : Finset (Fin 8)) 8 (deltaDist 0)
        = Real.log 16 := by
  have h := klDelta_add_klDelta_eq_log_mul
    (Finset.univ : Finset (Fin 2)) (Finset.univ : Finset (Fin 8))
    2 8 (by norm_num) (by norm_num) 0 (by decide) 0 (by decide)
  convert h using 2
  norm_num

/-- At `(m, n) = (8, 30)`: `D_KL(δ ‖ U_8) + D_KL(δ ‖ U_30) = log 240`. -/
theorem klDelta_8_plus_klDelta_30 :
    klToUniform (Finset.univ : Finset (Fin 8)) 8 (deltaDist 0)
      + klToUniform (Finset.univ : Finset (Fin 30)) 30 (deltaDist 0)
        = Real.log 240 := by
  have h := klDelta_add_klDelta_eq_log_mul
    (Finset.univ : Finset (Fin 8)) (Finset.univ : Finset (Fin 30))
    8 30 (by norm_num) (by norm_num) 0 (by decide) 0 (by decide)
  convert h using 2
  norm_num

/-! ### Headline -/

/-- **Headline (KL additivity at delta).** The KL divergence to the uniform
    is **additive** over independent moduli at the delta corner:

      `D_KL(δ_{r₀} ‖ U_m) + D_KL(δ_{s₀} ‖ U_n) = log(m · n)`.

    Specialised:

    * `(m, n) = (2, 8)`: `log 2 + log 8 = log 16`.
    * `(m, n) = (8, 30)`: `log 8 + log 30 = log 240`.
    * Equal-modulus: `2 · D_KL(δ ‖ U_m) = log(m²)`. -/
theorem KL_additivity_summary :
    -- General delta additivity
    (∀ {α β : Type*} [DecidableEq α] [DecidableEq β]
        (s : Finset α) (t : Finset β) (m n : ℝ),
      0 < m → 0 < n → ∀ (r₀ : α), r₀ ∈ s → ∀ (s₀ : β), s₀ ∈ t →
        klToUniform s m (deltaDist r₀) + klToUniform t n (deltaDist s₀)
          = Real.log (m * n))
    -- (2, 8) instance
    ∧ klToUniform (Finset.univ : Finset (Fin 2)) 2 (deltaDist 0)
       + klToUniform (Finset.univ : Finset (Fin 8)) 8 (deltaDist 0)
        = Real.log 16
    -- (8, 30) instance
    ∧ klToUniform (Finset.univ : Finset (Fin 8)) 8 (deltaDist 0)
       + klToUniform (Finset.univ : Finset (Fin 30)) 30 (deltaDist 0)
        = Real.log 240 :=
  ⟨@klDelta_add_klDelta_eq_log_mul,
   klDelta_2_plus_klDelta_8, klDelta_8_plus_klDelta_30⟩

end PT.Information
