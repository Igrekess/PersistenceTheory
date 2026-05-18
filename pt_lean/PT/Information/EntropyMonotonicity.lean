/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Information.GFTIdentity
import PT.Information.GFTSpecialisations
import PT.Information.BekensteinBound
import PT.Information.BekensteinExtensions
import PT.Information.EntropyBoundsTight
import PT.Information.MutualInformationBasic
import Mathlib.Analysis.SpecialFunctions.Log.NegMulLog
import Mathlib.Tactic

/-!
# Monotonicity properties of Shannon entropy

This file consolidates **monotonic / order** properties of the Shannon
entropy `shannonH s p`, re-exposing facts already established in the
`PT.Information.*` modules and adding a small number of corollaries
that follow algebraically from them.

The picture we assemble is the **complete order envelope**
`0 ≤ shannonH s p ≤ log m` on probability distributions over a finite
set `s` of cardinality `m`, with the two endpoints attained by deltas
(lower) and the uniform distribution (upper).

## Headline results

* `shannonH_nonneg_re` — `H ≥ 0` for any probability distribution
  supported in `[0, 1]` (re-exposes `BekensteinBound.shannonH_nonneg`).
* `shannonH_delta_eq_zero` — `H(δ_{r₀}) = 0` (re-exposes
  `GFTSpecialisations.shannonH_delta`).
* `shannonH_uniform_eq_log_re` — `H(U_m) = log m` (re-exposes
  `EntropyBoundsTight.shannonH_uniform_eq_log`).
* `shannonH_joint_eq_add` — `H(p ⊗ q) = H(p) + H(q)` (re-exposes
  `MutualInformationBasic.shannonH_joint`).
* `shannonH_joint_uniform_add_log` — `H(p ⊗ U_m) = H(p) + log m`,
  the **direct corollary** of additivity and the uniform-entropy formula.
* `shannonH_le_log_of_KL_nonneg` — `H ≤ log m` whenever the KL term is
  non-negative (Bekenstein restated as an entropy upper bound).
* `shannonH_strict_between` — strict bilateral bound `0 < H < log m`
  given `H > 0` and `D_KL > 0` (parametric form).
* `shannonH_monotonicity_summary` — single headline bundle.

## Scope (reduction)

We deliberately stay within the *algebraic* corollaries of
`GFT_identity` and `shannonH_joint`. The hard backward characterisations

* `H = 0 ↔ p = δ_{r₀}` and
* `H = log m ↔ p = U_m`

require Gibbs' inequality / strict-concavity arguments that are not
in scope for this file. We provide:

* the **forward** directions (`p = δ → H = 0`, `p = U_m → H = log m`)
  unconditionally, and
* the **strict bilateral bound** in *parametric* form, taking
  `0 < D_KL` as a hypothesis (which a separate Gibbs proof would
  discharge for `p ≠ U_m`).

This is consistent with the surrounding files
(`BekensteinExtensions.bekenstein_strict_of_entropy_pos`,
`EntropyBoundsTight.entropy_pos_of_strict_bekenstein`), which also
state strictness in parametric form.
-/

namespace PT.Information

open Real Finset

/-! ### 1. Non-negativity `H ≥ 0` (re-exposed) -/

/-- **`shannonH ≥ 0` (re-exposure of `BekensteinBound.shannonH_nonneg`).**
    For any distribution `p : α → ℝ` with values in `[0, 1]`,
    `shannonH s p ≥ 0`. -/
theorem shannonH_nonneg_re {α : Type*} (s : Finset α) (p : α → ℝ)
    (h0 : ∀ r ∈ s, 0 ≤ p r) (h1 : ∀ r ∈ s, p r ≤ 1) :
    0 ≤ shannonH s p :=
  shannonH_nonneg s p h0 h1

/-! ### 2. Delta-distribution lower endpoint -/

/-- **`H(δ_{r₀}) = 0` (re-exposure of `GFTSpecialisations.shannonH_delta`).** -/
theorem shannonH_delta_eq_zero {α : Type*} [DecidableEq α]
    (s : Finset α) (r₀ : α) (hr : r₀ ∈ s) :
    shannonH s (deltaDist r₀) = 0 :=
  shannonH_delta s r₀ hr

/-! ### 3. Additivity of entropy under independence (re-exposed) -/

/-- **`H(p ⊗ q) = H(p) + H(q)` (re-exposure of
    `MutualInformationBasic.shannonH_joint`).** -/
theorem shannonH_joint_eq_add
    {α β : Type*} (s : Finset α) (t : Finset β)
    (p : α → ℝ) (q : β → ℝ)
    (hp_sum : ∑ a ∈ s, p a = 1) (hq_sum : ∑ b ∈ t, q b = 1) :
    shannonH (s ×ˢ t) (joint p q) = shannonH s p + shannonH t q :=
  shannonH_joint s t p q hp_sum hq_sum

/-! ### 4. Adjoining an independent uniform shifts entropy by `log m` -/

/-- **Adjoining an independent uniform component adds `log m`.**
    `H(p ⊗ U_m) = H(p) + log m` whenever `∑ p = 1`, `0 < m`, and the
    cardinality of `t` equals `m`.

    This is the direct corollary of `shannonH_joint` (additivity) and
    `shannonH_uniform_eq_log` (`H(U_m) = log m`). -/
theorem shannonH_joint_uniform_add_log
    {α β : Type*} (s : Finset α) (t : Finset β)
    (p : α → ℝ) (hp_sum : ∑ a ∈ s, p a = 1)
    (m : ℝ) (hm : 0 < m) (hcard_t : (t.card : ℝ) = m) :
    shannonH (s ×ˢ t) (joint p (fun _ : β => (1 : ℝ) / m))
      = shannonH s p + Real.log m := by
  -- The uniform distribution on `t` sums to `1`.
  have hq_sum : ∑ _b ∈ t, (1 : ℝ) / m = 1 := by
    rw [Finset.sum_const, nsmul_eq_mul, hcard_t, mul_one_div]
    exact div_self (ne_of_gt hm)
  -- Apply additivity, then specialise the uniform's entropy to `log m`.
  rw [shannonH_joint s t p (fun _ : β => (1 : ℝ) / m) hp_sum hq_sum]
  rw [shannonH_uniform_eq_log t m hm hcard_t]

/-! ### 5. Uniform attains the upper endpoint `log m` -/

/-- **`H(U_m) = log m` (re-exposure of
    `EntropyBoundsTight.shannonH_uniform_eq_log`).** -/
theorem shannonH_uniform_eq_log_re {α : Type*} (s : Finset α) (m : ℝ)
    (hm : 0 < m) (hcard : (s.card : ℝ) = m) :
    shannonH s (fun _ : α => (1 : ℝ) / m) = Real.log m :=
  shannonH_uniform_eq_log s m hm hcard

/-! ### 6. Upper bound `H ≤ log m` (entropy form of Bekenstein) -/

/-- **`shannonH ≤ log m` — entropy form of the Bekenstein bound.**

    From the GFT identity `klToUniform + shannonH = log m`, the
    non-negativity of `klToUniform` gives the upper bound
    `shannonH ≤ log m`.

    *Note.* The hypothesis `0 ≤ klToUniform s m p` is Gibbs' inequality,
    which holds for all probability distributions on a finite set; a
    proof of it is out of scope for this file (it requires
    strict-concavity / Jensen). We therefore state the bound in
    parametric form, consistent with the existing
    `BekensteinExtensions.bekenstein_strict_of_entropy_pos`. -/
theorem shannonH_le_log_of_KL_nonneg
    {α : Type*} (s : Finset α) (m : ℝ) (hm : 0 < m)
    (p : α → ℝ)
    (hp_nonneg : ∀ r ∈ s, 0 ≤ p r)
    (hp_sum : ∑ r ∈ s, p r = 1)
    (hKL_nn : 0 ≤ klToUniform s m p) :
    shannonH s p ≤ Real.log m := by
  -- shannonH = log m - klToUniform ≤ log m
  have hid := GFT_identity s m hm p hp_nonneg hp_sum
  linarith

/-! ### 7. Strict bilateral bound (parametric) -/

/-- **Strict bilateral bound on `shannonH`** (parametric form):
    if `H > 0` and `D_KL > 0`, then `0 < H < log m`.

    The strictness on the right is obtained from the GFT identity
    `shannonH = log m - klToUniform` and `klToUniform > 0`. -/
theorem shannonH_strict_between
    {α : Type*} (s : Finset α) (m : ℝ) (hm : 0 < m)
    (p : α → ℝ)
    (hp_nonneg : ∀ r ∈ s, 0 ≤ p r)
    (hp_sum : ∑ r ∈ s, p r = 1)
    (hH_pos : 0 < shannonH s p)
    (hKL_pos : 0 < klToUniform s m p) :
    0 < shannonH s p ∧ shannonH s p < Real.log m := by
  refine ⟨hH_pos, ?_⟩
  have hid := GFT_identity s m hm p hp_nonneg hp_sum
  linarith

/-! ### 8. Headline monotonicity summary -/

/-- **Headline (Shannon entropy monotonicity envelope).**

    For any probability distribution `p : α → ℝ` on `s : Finset α` of
    cardinality `m`, with `0 ≤ p ≤ 1` and `∑ p = 1`:

    * **Lower endpoint reached by deltas.** `H(δ_{r₀}) = 0`, and
      `H(p) ≥ 0` for every probability distribution.
    * **Upper endpoint reached by uniform.** `H(U_m) = log m`.
    * **Additivity under independence.** `H(p ⊗ q) = H(p) + H(q)`.
    * **Adjoining a uniform shifts entropy by `log m`.**
      `H(p ⊗ U_m) = H(p) + log m`. -/
theorem shannonH_monotonicity_summary
    {α β : Type*} [DecidableEq α]
    (s : Finset α) (t : Finset β)
    (m : ℝ) (hm : 0 < m) (hcard_s : (s.card : ℝ) = m)
    (hcard_t : (t.card : ℝ) = m)
    (p : α → ℝ)
    (hp_nonneg : ∀ r ∈ s, 0 ≤ p r)
    (hp_le_one : ∀ r ∈ s, p r ≤ 1)
    (hp_sum : ∑ r ∈ s, p r = 1)
    (r₀ : α) (hr : r₀ ∈ s)
    (q : β → ℝ) (hq_sum : ∑ b ∈ t, q b = 1) :
    -- 1. Non-negativity
    0 ≤ shannonH s p
    -- 2. Delta endpoint
    ∧ shannonH s (deltaDist r₀) = 0
    -- 3. Uniform endpoint
    ∧ shannonH s (fun _ : α => (1 : ℝ) / m) = Real.log m
    -- 4. Additivity under independence
    ∧ shannonH (s ×ˢ t) (joint p q) = shannonH s p + shannonH t q
    -- 5. Adjoining a uniform shifts by `log m`
    ∧ shannonH (s ×ˢ t) (joint p (fun _ : β => (1 : ℝ) / m))
        = shannonH s p + Real.log m := by
  refine ⟨?_, ?_, ?_, ?_, ?_⟩
  · exact shannonH_nonneg_re s p hp_nonneg hp_le_one
  · exact shannonH_delta_eq_zero s r₀ hr
  · exact shannonH_uniform_eq_log_re s m hm hcard_s
  · exact shannonH_joint_eq_add s t p q hp_sum hq_sum
  · exact shannonH_joint_uniform_add_log s t p hp_sum m hm hcard_t

end PT.Information
