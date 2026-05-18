/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.FixedPoint.T7MuStar
import PT.FixedPoint.T7GlobalUniqueness
import Mathlib.Data.Finset.Card
import Mathlib.Tactic

/-!
# Dimension protection (audit row #35)

**Statement (paper-level).**
The cardinality `|P_active|` of the **persistence-active prime set** is
**invariant** under the PT sieve dynamics on the non-degenerate domain
`μ ≥ 7`:

$$|P_{\text{active}}(\mu)| \;=\; 3 \qquad\forall\,\mu\ge 7.$$

This is the **dimension-protection** statement of monograph Ch08 (audit row
#35): the dimension count of the active set is stable under the sieve
dynamics. At the fixed point `μ* = 15` the cardinality is `|{3,5,7}| = 3`,
and `F_pers` is constant at `15` for all `μ ≥ 7` (cf.
`PT.FixedPoint.T7GlobalUniqueness.Fpers_eq_fifteen_of_ge_seven`), so the
active set never gains or loses elements.

In this file we:

* Define the (combinatorial) persistence-active set at `μ` as
  `{3} · 1_{μ ≥ 3} ∪ {5} · 1_{μ ≥ 5} ∪ {7} · 1_{μ ≥ 7}`.
* Prove its cardinality is `3` for every `μ ≥ 7`.
* Prove a "dimension cap" at three: `|P_active(μ)| ≤ 3` for every `μ`.
* Express `F_pers` as the sum over `P_active`.

## Reference

Audit row #35 ("Dimension protection"), monograph Ch08 §"Protection
dimensionnelle". Companion to `T7MuStar` and `T7GlobalUniqueness`.
-/

namespace PT.FixedPoint

open Finset

/-! ### The persistence-active set as an explicit finset -/

/-- The persistence-active prime set at integer `μ`, in finset form.
    Element `p ∈ {3, 5, 7}` is included iff `p ≤ μ`. -/
def persActive (μ : ℕ) : Finset ℕ :=
  (if 3 ≤ μ then ({3} : Finset ℕ) else ∅) ∪
  (if 5 ≤ μ then ({5} : Finset ℕ) else ∅) ∪
  (if 7 ≤ μ then ({7} : Finset ℕ) else ∅)

/-! ### Small-`μ` snapshots -/

theorem persActive_lt3 (μ : ℕ) (h : μ < 3) : persActive μ = ∅ := by
  unfold persActive
  have h3 : ¬ 3 ≤ μ := by omega
  have h5 : ¬ 5 ≤ μ := by omega
  have h7 : ¬ 7 ≤ μ := by omega
  simp [h3, h5, h7]

theorem persActive_4 : persActive 4 = ({3} : Finset ℕ) := by
  unfold persActive; decide

theorem persActive_5 : persActive 5 = ({3, 5} : Finset ℕ) := by
  unfold persActive; decide

theorem persActive_6 : persActive 6 = ({3, 5} : Finset ℕ) := by
  unfold persActive; decide

theorem persActive_7 : persActive 7 = ({3, 5, 7} : Finset ℕ) := by
  unfold persActive; decide

theorem persActive_15 : persActive 15 = ({3, 5, 7} : Finset ℕ) := by
  unfold persActive; decide

/-! ### Dimension constancy on `μ ≥ 7` -/

/-- **Dimension protection (constancy).** For every `μ ≥ 7`,
    `persActive μ = {3, 5, 7}`. -/
theorem persActive_eq_three_primes (μ : ℕ) (h : 7 ≤ μ) :
    persActive μ = ({3, 5, 7} : Finset ℕ) := by
  unfold persActive
  have h3 : 3 ≤ μ := by omega
  have h5 : 5 ≤ μ := by omega
  rw [if_pos h3, if_pos h5, if_pos h]
  ext n
  simp

/-- **Dimension protection (cardinality).** For every `μ ≥ 7`,
    `|persActive μ| = 3`. -/
theorem persActive_card_eq_three (μ : ℕ) (h : 7 ≤ μ) :
    (persActive μ).card = 3 := by
  rw [persActive_eq_three_primes μ h]
  decide

/-! ### Dimension cap: `|persActive| ≤ 3` always -/

/-- **Dimension cap.** For every `μ`, `|persActive μ| ≤ 3`. The combinatorial
    active set never exceeds three primes, regardless of `μ`. -/
theorem persActive_card_le_three (μ : ℕ) : (persActive μ).card ≤ 3 := by
  unfold persActive
  split_ifs <;> decide

/-! ### Sum over persistence-active = `F_pers` -/

/-- **Sum identity.** `F_pers μ = ∑ p ∈ persActive μ, p` for every `μ`.

    This is the constructive translation: the persistence-active map is the
    sum over the persistence-active finset. -/
theorem Fpers_eq_sum_persActive (μ : ℕ) :
    Fpers μ = (persActive μ).sum id := by
  unfold Fpers persActive
  split_ifs with h3 h5 h7
  all_goals (simp; try omega)

/-! ### Headline dimension-protection statement -/

/-- **Dimension protection (headline, audit row #35).**

    On the non-degenerate domain `μ ≥ 7`, the persistence-active set is
    constant: `persActive μ = {3, 5, 7}`, cardinality `3`. Equivalently, the
    sieve dynamics `F_pers` is constant at `15` and the active-set dimension
    is protected at exactly `3` (the cascade dimension predicted by PT). -/
theorem dimension_protection (μ : ℕ) (h : 7 ≤ μ) :
    persActive μ = ({3, 5, 7} : Finset ℕ)
    ∧ (persActive μ).card = 3
    ∧ Fpers μ = 15 :=
  ⟨persActive_eq_three_primes μ h,
   persActive_card_eq_three μ h,
   Fpers_eq_fifteen_of_ge_seven μ h⟩

/-! ### Connection to `T7MuStar.persActiveAt15` -/

/-- The finset `persActive 15` matches `persActiveAt15` from `T7MuStar`. -/
theorem persActive_15_eq_persActiveAt15 :
    persActive 15 = persActiveAt15 := by
  rw [persActive_15]
  rfl

end PT.FixedPoint
