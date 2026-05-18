/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Information.GFTIdentity

/-!
# Bekenstein Bound — `D_KL(P ‖ U_m) ≤ log m`

**Statement (paper-level, Ch04 §"Borne de Bekenstein").**
For any probability distribution `P = (p_r)` on a finite set of cardinality
`m`,

$$D_{KL}(P \,\Vert\, U_m) \;\le\; \log m,$$

with equality iff `P` is a delta (point-mass) on a single class.

**Proof.** From the GFT identity `log m = D_KL + H` and the
non-negativity of Shannon entropy `H(P) ≥ 0` (each summand
`-p_r log p_r ≥ 0` on `[0,1]`), one obtains `D_KL = log m - H ≤ log m`.

## References

Monograph chapter 4, §"Borne de Bekenstein", `\label{thm:bekenstein}`.
M4 article (PT_ARTICLES/PT_MATHEMATICS/M4), `\label{thm:bekenstein}`.
-/

namespace PT.Information

open Real Finset

/-- The Shannon entropy of a probability distribution `p : α → ℝ` supported
    in `[0, 1]` (per-bin) is non-negative. -/
theorem shannonH_nonneg {α : Type*} (s : Finset α) (p : α → ℝ)
    (h0 : ∀ r ∈ s, 0 ≤ p r) (h1 : ∀ r ∈ s, p r ≤ 1) :
    0 ≤ shannonH s p := by
  unfold shannonH
  apply Finset.sum_nonneg
  intro r hr
  exact Real.negMulLog_nonneg (h0 r hr) (h1 r hr)

/-- **Bekenstein bound (GFT corollary).**
    For any probability distribution `p : α → ℝ` on a finite set `s` of
    cardinality `m`, `D_KL(P ‖ U_m) ≤ log m`. -/
theorem bekenstein_bound {α : Type*} (s : Finset α) (m : ℝ) (hm : 0 < m)
    (p : α → ℝ)
    (hp_nonneg : ∀ r ∈ s, 0 ≤ p r)
    (hp_le_one : ∀ r ∈ s, p r ≤ 1)
    (hp_sum : ∑ r ∈ s, p r = 1) :
    klToUniform s m p ≤ Real.log m := by
  have hid := GFT_identity s m hm p hp_nonneg hp_sum
  -- hid : klToUniform + shannonH = log m
  have hH : 0 ≤ shannonH s p := shannonH_nonneg s p hp_nonneg hp_le_one
  linarith

/-- **Bekenstein bound for a natural-number modulus `m ≥ 1`.** -/
theorem bekenstein_bound_nat {α : Type*} (s : Finset α) (m : ℕ) (hm : 1 ≤ m)
    (p : α → ℝ)
    (hp_nonneg : ∀ r ∈ s, 0 ≤ p r)
    (hp_le_one : ∀ r ∈ s, p r ≤ 1)
    (hp_sum : ∑ r ∈ s, p r = 1) :
    klToUniform s (m : ℝ) p ≤ Real.log (m : ℝ) :=
  bekenstein_bound s (m : ℝ)
    (by exact_mod_cast Nat.lt_of_lt_of_le Nat.zero_lt_one hm)
    p hp_nonneg hp_le_one hp_sum

end PT.Information
