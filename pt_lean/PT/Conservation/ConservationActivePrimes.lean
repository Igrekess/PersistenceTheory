/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Conservation.ConservationID
import PT.Conservation.ConservationIDExtensions
import Mathlib.Algebra.BigOperators.Intervals
import Mathlib.Tactic

/-!
# Conservation Identity at Active Primes — `{3, 5, 7}` and `μ* = 15`

This file specialises the conservation identity
`∑_{n=1}^{N} g_n = p_{N+1} - 2`
(see `PT.Conservation.ConservationID`) to the **PT active prime indices**
`{3, 5, 7}` — the canonical PT cascade.

Two arithmetic quantities coexist on the first few primes:

* the **cumulative gap sum** `Σ_{n=1}^{N} g_n = p_{N+1} - 2` (conservation
  identity, telescoping),
* the **sum of active primes** `3 + 5 + 7 = 15 = μ*` (T7).

This module exposes the bridge between the two pictures (gap-cumulative
vs active-prime-sum) and proves that they are *distinct* arithmetic
witnesses (`9 ≠ 15` at `N = 4`).

All identities are concrete finite arithmetic on the first five `ptPrime`
values `(2, 3, 5, 7, 11)`; everything is discharged by `decide` or
`unfold gap; decide`.

## Reference

* Monograph Ch03 §3.1 (`\label{thm:conservation-id}`).
* `PT/FixedPoint/T7MuStar.lean` — `muStar = 15`, active set `{3, 5, 7}`.
* `PT/Holonomy/ActivePrimeCriterion.lean` — `IsActive`.
-/

namespace PT.Conservation

open Finset

/-! ### Active gaps — gaps that *bridge* an active prime

In PT the active prime indices are `{3, 5, 7}`. In the `ptPrime`
sequence `(p_1, p_2, p_3, p_4, p_5) = (2, 3, 5, 7, 11)`, the active
primes occupy positions `2`, `3`, `4`. We single out three "active
gaps" — the gaps that *land on* an active prime:

* `gap ptPrime 1 = p_2 - p_1 = 3 - 2 = 1` (lands on `3`),
* `gap ptPrime 2 = p_3 - p_2 = 5 - 3 = 2` (lands on `5`),
* `gap ptPrime 3 = p_4 - p_3 = 7 - 5 = 2` (lands on `7`).

Their sum is `1 + 2 + 2 = 5 = p_4 - 2`, recovering `conservation_N3`.
-/

/-- Sum of the three "active gaps" — the gaps whose right endpoint is an
    active prime `{3, 5, 7}` of the PT cascade. -/
def activeGapsSum : ℤ :=
  gap ptPrime 1 + gap ptPrime 2 + gap ptPrime 3

/-- The active-gaps sum equals `5 = p_4 - 2`. -/
theorem activeGapsSum_eq_five : activeGapsSum = 5 := by
  unfold activeGapsSum gap; decide

/-- The active-gaps sum recovers the `N = 3` conservation cumulative. -/
theorem activeGapsSum_eq_cumulative_N3 :
    activeGapsSum = ∑ n ∈ Ico 1 4, gap ptPrime n := by
  unfold activeGapsSum
  decide

/-- The active-gaps sum equals `p_4 - p_1 = 7 - 2 = 5` (a direct
    application of the telescoping identity). -/
theorem activeGapsSum_eq_telescope :
    activeGapsSum = ptPrime 4 - ptPrime 1 := by
  unfold activeGapsSum gap; decide

/-! ### Sum of the active primes themselves: `3 + 5 + 7 = 15 = μ*`

This is *not* a gap sum — it is the sum of the **values** of the active
primes. Numerically it equals `μ*`, the PT fixed-point characteristic
length (`T7MuStar.lean`). It is recovered here as a direct arithmetic
identity on `ptPrime 2, ptPrime 3, ptPrime 4`.
-/

/-- Sum of the active primes `p_2 + p_3 + p_4 = 3 + 5 + 7`. -/
def activePrimesSum : ℤ :=
  ptPrime 2 + ptPrime 3 + ptPrime 4

/-- The active-primes sum equals `15`, i.e. `μ*` of the PT fixed point. -/
theorem activePrimesSum_eq_muStar : activePrimesSum = 15 := by
  unfold activePrimesSum; decide

/-- Restated: `3 + 5 + 7 = 15`. -/
theorem three_plus_five_plus_seven : (3 : ℤ) + 5 + 7 = 15 := by decide

/-! ### Gap-cumulative vs active-prime-sum — they differ

The two arithmetic witnesses count *different things*:

* `cumulative N` counts the **total displacement** from `p_1 = 2` up to
  `p_{N+1}`.
* `activePrimesSum` counts the **sum of values** of the active primes.

At `N = 3` (the natural stopping point of the active cascade), the
cumulative is `5 = p_4 - 2`, while the active-prime sum is `15`. The
difference is `2 · p_1 + ?` — explicitly the value `15 - 5 = 10`. We
record this as a structural identity below.
-/

/-- The active-prime sum equals the cumulative gap sum on `Ico 1 4`,
    *plus* `2 · p_1 + (p_3 - p_1)` — i.e. `activePrimesSum
    = cumulative_N3 + 10`. -/
theorem activePrimesSum_eq_cumulative_plus_ten :
    activePrimesSum = (∑ n ∈ Ico 1 4, gap ptPrime n) + 10 := by
  unfold activePrimesSum
  decide

/-- **Distinctness witness.** The cumulative gap sum on `Ico 1 4`
    (i.e. `p_4 - 2 = 5`) is strictly less than the active-prime sum
    `3 + 5 + 7 = 15`. -/
theorem cumulative_lt_activePrimesSum :
    (∑ n ∈ Ico 1 4, gap ptPrime n) < activePrimesSum := by
  unfold activePrimesSum
  decide

/-- **Distinctness witness, `N = 4`.** The cumulative gap sum on
    `Ico 1 5` (i.e. `p_5 - 2 = 9`) is *still* strictly less than the
    active-prime sum `15`. -/
theorem cumulative_N4_lt_activePrimesSum :
    (∑ n ∈ Ico 1 5, gap ptPrime n) < activePrimesSum := by
  unfold activePrimesSum
  decide

/-- **Concrete numerical separation.** `9 ≠ 15`: the cumulative gap sum
    at `N = 4` (i.e. `9`) and the active-prime sum `μ* = 15` are
    distinct integers. -/
theorem cumulative_N4_ne_muStar :
    (∑ n ∈ Ico 1 5, gap ptPrime n) ≠ 15 := by
  decide

/-! ### Bridge identities: gaps expressed via `ptPrime`

These are immediate by definition of `gap`; we record them as named
lemmas so downstream files can rewrite freely. -/

/-- `gap ptPrime 1 = ptPrime 2 - ptPrime 1`. -/
theorem gap_one_eq_diff : gap ptPrime 1 = ptPrime 2 - ptPrime 1 := by
  unfold gap; decide

/-- `gap ptPrime 2 = ptPrime 3 - ptPrime 2`. -/
theorem gap_two_eq_diff : gap ptPrime 2 = ptPrime 3 - ptPrime 2 := by
  unfold gap; decide

/-- `gap ptPrime 3 = ptPrime 4 - ptPrime 3`. -/
theorem gap_three_eq_diff : gap ptPrime 3 = ptPrime 4 - ptPrime 3 := by
  unfold gap; decide

/-- `gap ptPrime 4 = ptPrime 5 - ptPrime 4`. -/
theorem gap_four_eq_diff : gap ptPrime 4 = ptPrime 5 - ptPrime 4 := by
  unfold gap; decide

/-! ### Headline summary -/

/-- **Headline (active-prime conservation summary).** Putting everything
    together:

    * `activeGapsSum = 5`              (cumulative on `{1, 2, 3}` of gaps)
    * `activeGapsSum = p_4 - p_1`      (telescoping)
    * `activePrimesSum = 15 = μ*`      (sum of the active prime values)
    * cumulative at `N = 4` is `9`, distinct from `μ* = 15`. -/
theorem active_primes_conservation_summary :
    activeGapsSum = 5
    ∧ activeGapsSum = ptPrime 4 - ptPrime 1
    ∧ activePrimesSum = 15
    ∧ (∑ n ∈ Ico 1 5, gap ptPrime n) = 9
    ∧ (∑ n ∈ Ico 1 5, gap ptPrime n) ≠ activePrimesSum :=
  ⟨activeGapsSum_eq_five,
   activeGapsSum_eq_telescope,
   activePrimesSum_eq_muStar,
   conservation_N4,
   by unfold activePrimesSum; decide⟩

end PT.Conservation
