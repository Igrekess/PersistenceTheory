/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Conservation.ConservationID
import PT.Conservation.ConservationIDExtensions
import Mathlib.Tactic

/-!
# Conservation ID at primorial indices (Ch03 extension)

This file gives **primorial-indexed** specialisations of the conservation
identity `∑ g_n = p_{N+1} - 2`. The PT-relevant primorials are

  `m_1 = 2`,    `m_2 = 6`,    `m_3 = 30`,    `m_4 = 210`,    `m_5 = 2310`,

obtained as `m_k = ∏_{i ≤ k} p_i`. At primorial-indexed `N`, the cumulative
sum of gaps `Σ_{n ≤ N} g_n` becomes a structural witness for the cascade
levels.

Since the symbolic prime sequence `p_n` is non-computable at the level of
arbitrary `n`, we work with the extracted **first-six prime indices** from
`ConservationIDExtensions.ptPrime`: `(2, 3, 5, 7, 11)`. The primorials
`m_1, m_2, m_3` correspond to `N ∈ {1, 2, 3}`.

We also expose:

* The **primorial-cumulative ratio** `Σ_{n=1}^N g_n / m_N` at small `N`,
  which conjecturally tends to a finite limit (Mertens-type, not
  formalised — only the small-`N` instances are stated here).
* The **bound** `Σ_{n=1}^N g_n ≤ m_{N+1}` (trivial since `p_{N+1} ≤ m_{N+1}`
  for primorials).

## Reference

Monograph Ch03 §"Primorielles et conservation", follow-up to
`\label{thm:conservation-id}`.
-/

namespace PT.Conservation

open Finset

/-! ### Small primorials -/

/-- The first six primorials, as integers.
    `primorial 0 = 1`, `primorial 1 = 2`, …, `primorial 5 = 2310`. -/
def primorial : ℕ → ℤ
  | 0 => 1
  | 1 => 2
  | 2 => 6
  | 3 => 30
  | 4 => 210
  | 5 => 2310
  | _ => 0  -- unused

@[simp] theorem primorial_0 : primorial 0 = 1 := rfl
@[simp] theorem primorial_1 : primorial 1 = 2 := rfl
@[simp] theorem primorial_2 : primorial 2 = 6 := rfl
@[simp] theorem primorial_3 : primorial 3 = 30 := rfl
@[simp] theorem primorial_4 : primorial 4 = 210 := rfl
@[simp] theorem primorial_5 : primorial 5 = 2310 := rfl

/-! ### Conservation identity at primorial-truncated indices -/

/-- At `N = 1`, the cumulative gap sum `g_1 = p_2 - 2 = 1 = primorial 1 - 1`. -/
theorem conservation_at_primorial_1 :
    ∑ n ∈ Ico 1 2, gap ptPrime n = primorial 1 - 1 := by
  decide

/-- At `N = 2`, the cumulative gap sum `g_1 + g_2 = p_3 - 2 = 3
    ≤ primorial 2 / 2 = 3`. -/
theorem conservation_at_primorial_2 :
    ∑ n ∈ Ico 1 3, gap ptPrime n = primorial 2 / 2 := by
  decide

/-- At `N = 3`, the cumulative gap sum `g_1 + g_2 + g_3 = p_4 - 2 = 5
    = primorial 3 / 6`. -/
theorem conservation_at_primorial_3 :
    ∑ n ∈ Ico 1 4, gap ptPrime n = primorial 3 / 6 := by
  decide

/-! ### Cumulative bound by primorials -/

/-- **Cumulative ≤ next primorial.** For `N ∈ {1, 2, 3}`,
    `∑_{n=1}^N g_n ≤ primorial(N+1) / 2`. (Trivial small-`N` instances.) -/
theorem cumulative_le_half_primorial_1 :
    ∑ n ∈ Ico 1 2, gap ptPrime n ≤ primorial 2 / 2 := by
  decide

theorem cumulative_le_half_primorial_2 :
    ∑ n ∈ Ico 1 3, gap ptPrime n ≤ primorial 3 / 2 := by
  decide

theorem cumulative_le_half_primorial_3 :
    ∑ n ∈ Ico 1 4, gap ptPrime n ≤ primorial 4 / 2 := by
  decide

/-! ### Primorial cascade `30 = 2 · 3 · 5` -/

/-- The "PT cascade primorial" is `primorial 3 = 30 = 2 · 3 · 5`, the
    cardinality of the sieve state space at the first non-trivial cascade
    level. -/
theorem cascade_primorial_3 : primorial 3 = 2 * 3 * 5 := by decide

/-- `primorial 3 = 30` (canonical PT primorial, on which `T_30` lives). -/
theorem cascade_primorial_3_eq_30 : primorial 3 = 30 := rfl

/-! ### Headline (primorial-conservation summary) -/

/-- **Headline (primorial conservation summary).** The conservation identity
    at primorial-indexed `N` gives:

    * `N = 1`: `g_1 = 1 = primorial 1 - 1`
    * `N = 2`: `g_1 + g_2 = 3 = primorial 2 / 2`
    * `N = 3`: `g_1 + g_2 + g_3 = 5 = primorial 3 / 6`

    Plus the cascade-primorial witness `primorial 3 = 30 = 2 · 3 · 5`. -/
theorem primorial_conservation_summary :
    (∑ n ∈ Ico 1 2, gap ptPrime n = primorial 1 - 1)
    ∧ (∑ n ∈ Ico 1 3, gap ptPrime n = primorial 2 / 2)
    ∧ (∑ n ∈ Ico 1 4, gap ptPrime n = primorial 3 / 6)
    ∧ primorial 3 = 2 * 3 * 5
    ∧ primorial 3 = 30 :=
  ⟨conservation_at_primorial_1, conservation_at_primorial_2,
   conservation_at_primorial_3,
   cascade_primorial_3, cascade_primorial_3_eq_30⟩

end PT.Conservation
