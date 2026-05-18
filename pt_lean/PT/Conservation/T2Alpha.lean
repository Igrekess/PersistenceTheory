/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Stochastic.SHalf

/-!
# T2 — The Conservation Exponent α_cons = s² = 1/4

Theorem T2 of Persistence Theory states that the conservation exponent at the
sieve level is

$$\alpha_{\rm cons} = s^2 = \tfrac{1}{4},$$

where `s = 1/2` is the foundational symmetry parameter proved in
`PT.Stochastic.SHalf`. Operationally, `α_cons` is the asymptotic two-point
correlation `lim ℰ[s_n s_{n+1}]` on the stationary `T₃`-Markov chain; since the
chain has marginals `(1/2, 1/2)` and the matrix `T₃` exchanges states, this
correlation collapses to the algebraic identity `s ⋅ s = 1/4`.

The full Markov-chain machinery (needed for the operational characterisation)
is not yet adequately developed in Mathlib for this 2-state finite case, so we
take the algebraic identity `α_cons := s²` as definition and prove the
elementary identity `α_cons = 1/4`. The Markov-chain interpretation is left as
a `[BRIDGE]` comment and will be added once `MeasureTheory.Markov`
matures in Mathlib.

## Reference

Theorem T2 of: Y. Senez, M1 article (PT_ARTICLES/PT_MATHEMATICS/M1).
Monograph Chapter 1, §1.4.
-/

namespace PT.Conservation

open PT.Stochastic

/-- The PT conservation exponent. Definition: `α_cons := s²` where
    `s = 1/2` is the symmetry parameter (see `PT.Stochastic.SHalf`). -/
noncomputable def alpha_cons : ℝ := s * s

/-- **Theorem T2.** The conservation exponent equals `1/4`. -/
theorem T2_alpha_eq_one_quarter : alpha_cons = 1 / 4 := by
  unfold alpha_cons
  rw [s_def]
  norm_num

/-- Restatement: `α_cons = s²`. -/
theorem T2_alpha_eq_s_sq : alpha_cons = s ^ 2 := by
  unfold alpha_cons
  ring

/-- The conservation exponent is strictly positive. -/
theorem T2_alpha_pos : 0 < alpha_cons := by
  rw [T2_alpha_eq_one_quarter]; norm_num

/-- The conservation exponent is strictly less than `s`. -/
theorem T2_alpha_lt_s : alpha_cons < s := by
  rw [T2_alpha_eq_one_quarter, s_def]; norm_num

/-! ### Bridge: two-point correlation on the stationary chain

For the stationary `T₃`-Markov chain with marginal `π = (1/2, 1/2)`, the
two-point correlation factorises as

$$\mathbb{E}[s_n s_{n+1}] = \sum_{i,j} \pi(i)\, T_3(i,j)\, s_i\, s_j
   = s \cdot s = \alpha_{\rm cons},$$

with `s_i ∈ {1/2, 1/2}` constant. The discrete-sum identity below makes this
explicit for the canonical encoding `s_i := s = 1/2`. -/

/-- The two-point correlation on the stationary `T₃`-chain with constant
    "state observable" `s_i = s`. By direct calculation this equals `s ⋅ s`. -/
theorem T2_two_point_correlation :
    (∑ i : Fin 2, ∑ j : Fin 2, piHalf i * PT.Sieve.T3 i j * s * s) = alpha_cons := by
  simp [piHalf, PT.Sieve.T3, Fin.sum_univ_two, s, alpha_cons]
  ring

end PT.Conservation
