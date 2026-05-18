/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Conservation.ConservationID
import PT.Conservation.ConservationIDExtensions
import PT.Conservation.ConservationIDPrimorial
import PT.Conservation.ConservationActivePrimes
import PT.Conservation.CumulativeBoundsExtended
import PT.Conservation.PrimeGapMoments
import PT.Conservation.GapBoundedBelow
import PT.Conservation.GapMaxBound
import Mathlib.Tactic

/-!
# Gap-sum master identity — synthesis aggregator

This file is a **pure synthesis aggregator**. It does not introduce any
new computational content: every constituent identity is already proved
in one of the imported modules. The point of this module is to give a
single canonical reference theorem that bundles all the gap-sum
identities and bounds proved across the `PT.Conservation` namespace.

## What is bundled

The master synthesis collects:

1. **Generic telescoping** (master form, from `ConservationID`):
   `∀ p N, ∑_{n=1}^N gap p n = p (N+1) - p 1`.

2. **PT specialisation** (with `ptPrime 1 = 2`):
   `∀ N, ∑_{n=1}^N gap ptPrime n = ptPrime (N+1) - 2`.

3. **Four explicit cumulative values** for `N ∈ {1, 2, 3, 4}`
   (from `ConservationIDExtensions`).

4. **Three primorial-aligned values** for `N ∈ {1, 2, 3}`
   (from `ConservationIDPrimorial`).

5. **Active primes sum** `3 + 5 + 7 = 15 = μ*` and active-gaps sum
   (from `ConservationActivePrimes`).

6. **Six extended cumulative values** for `N ∈ {5, …, 10}`
   (from `CumulativeBoundsExtended`).

7. **Individual gap bounds**: `g_n ≥ 1` for `n ∈ {1, …, 4}`
   (from `GapBoundedBelow`) and `g_n ≤ 6` for `n ∈ {1, …, 10}`
   (from `GapMaxBound`).

8. **Higher moments at `N = 4`**: square sum `= 25`, cube sum `= 81`
   (from `PrimeGapMoments`).

## Unified provenance

Items (1)–(6) all follow from the single source `sum_gap_telescope`
applied to the appropriate sequence (`ptPrime` for items 2–5/8;
`ptPrimeExt` for item 6). Item (7) follows from finite case analysis
on the values `(g_1, …, g_10) = (1, 2, 2, 4, 2, 4, 2, 4, 6, 2)`.

## Reference

* Monograph Chapter 3, §3.1, `\label{thm:conservation-id}`.
* M1 article (PT_ARTICLES/PT_MATHEMATICS/M1), Conservation identity.
* This file is the canonical entry point for any external reference
  to "the PT gap-sum identities".
-/

namespace PT.Conservation

open Finset

/-! ### 1. Master telescoping (generic form, re-export) -/

/-- **Master telescoping identity (generic form).** This is a verbatim
    re-export of `sum_gap_telescope`, named to signal that it is the
    canonical entry point for the master synthesis. -/
theorem master_telescope_identity (p : ℕ → ℤ) (N : ℕ) :
    ∑ n ∈ Ico 1 (N + 1), gap p n = p (N + 1) - p 1 :=
  sum_gap_telescope p N

/-! ### 2. PT specialisation -/

/-- **Master PT telescoping.** Since `ptPrime 1 = 2`, the generic
    telescope specialises to `∑ g_n = ptPrime (N+1) - 2`. -/
theorem master_pt_telescope (N : ℕ) :
    ∑ n ∈ Ico 1 (N + 1), gap ptPrime n = ptPrime (N + 1) - 2 := by
  have h := master_telescope_identity ptPrime N
  rw [ptPrime_one] at h
  exact h

/-- **Master PT-extended telescoping.** Same identity for `ptPrimeExt`,
    valid for `N ∈ {0, …, 10}` where `ptPrimeExt` is concretely
    defined. -/
theorem master_ptExt_telescope (N : ℕ) :
    ∑ n ∈ Ico 1 (N + 1), gap ptPrimeExt n = ptPrimeExt (N + 1) - 2 := by
  have h := master_telescope_identity ptPrimeExt N
  rw [ptPrimeExt_one] at h
  exact h

/-! ### 3. Cross-check at `N = 4` via the master form -/

/-- `conservation_N4` follows directly from `master_pt_telescope`. -/
theorem conservation_N4_from_master :
    ∑ n ∈ Ico 1 5, gap ptPrime n = 9 := by
  have h := master_pt_telescope 4
  simp [ptPrime_five] at h
  exact h

/-- `cumulativeExt_N10` follows directly from `master_ptExt_telescope`. -/
theorem cumulativeExt_N10_from_master :
    ∑ n ∈ Ico 1 11, gap ptPrimeExt n = 29 := by
  have h := master_ptExt_telescope 10
  simp [ptPrimeExt_eleven] at h
  exact h

/-! ### 4. Headline master synthesis -/

/-- **Gap sum master summary.** The single canonical theorem bundling
    every gap-sum identity proved in `PT.Conservation`:

    * **(A) Generic telescoping**: for any integer-valued sequence `p`,
      `∑_{n=1}^N gap p n = p (N+1) - p 1`.
    * **(B) PT specialisation**: with `ptPrime 1 = 2`,
      `∑_{n=1}^N gap ptPrime n = ptPrime (N+1) - 2`.
    * **(C) Four explicit values** `N ∈ {1, 2, 3, 4}`:
      `(1, 3, 5, 9)`.
    * **(D) Three primorial-aligned values** `N ∈ {1, 2, 3}`:
      `(primorial 1 - 1, primorial 2 / 2, primorial 3 / 6)
        = (1, 3, 5)`.
    * **(E) Active primes sum**: `activePrimesSum = 15 = μ*`, and
      `activeGapsSum = ptPrime 4 - ptPrime 1 = 5`.
    * **(F) Six extended values** `N ∈ {5, …, 10}`:
      `(11, 15, 17, 21, 27, 29)`.
    * **(G) Individual gap bounds**: `g_n ≥ 1` for `n ∈ {1, …, 4}`
      and `g_n ≤ 6` for `n ∈ {1, …, 10}`.
    * **(H) Higher moments at `N = 4`**: `sqSumGaps 4 = 25`,
      `cubeSumGaps 4 = 81`. -/
theorem gap_sum_master_summary :
    -- (A) generic telescoping
    (∀ (p : ℕ → ℤ) (N : ℕ),
        ∑ n ∈ Ico 1 (N + 1), gap p n = p (N + 1) - p 1)
    -- (B) PT specialisation
    ∧ (∀ N : ℕ, ∑ n ∈ Ico 1 (N + 1), gap ptPrime n = ptPrime (N + 1) - 2)
    -- (C) four explicit cumulative values
    ∧ (∑ n ∈ Ico 1 2, gap ptPrime n = 1)
    ∧ (∑ n ∈ Ico 1 3, gap ptPrime n = 3)
    ∧ (∑ n ∈ Ico 1 4, gap ptPrime n = 5)
    ∧ (∑ n ∈ Ico 1 5, gap ptPrime n = 9)
    -- (D) three primorial-aligned values
    ∧ (∑ n ∈ Ico 1 2, gap ptPrime n = primorial 1 - 1)
    ∧ (∑ n ∈ Ico 1 3, gap ptPrime n = primorial 2 / 2)
    ∧ (∑ n ∈ Ico 1 4, gap ptPrime n = primorial 3 / 6)
    -- (E) active primes
    ∧ (activePrimesSum = 15)
    ∧ (activeGapsSum = ptPrime 4 - ptPrime 1)
    -- (F) six extended cumulative values
    ∧ (∑ n ∈ Ico 1 6,  gap ptPrimeExt n = 11)
    ∧ (∑ n ∈ Ico 1 7,  gap ptPrimeExt n = 15)
    ∧ (∑ n ∈ Ico 1 8,  gap ptPrimeExt n = 17)
    ∧ (∑ n ∈ Ico 1 9,  gap ptPrimeExt n = 21)
    ∧ (∑ n ∈ Ico 1 10, gap ptPrimeExt n = 27)
    ∧ (∑ n ∈ Ico 1 11, gap ptPrimeExt n = 29)
    -- (G) individual gap bounds
    ∧ (gap ptPrime 1 ≥ 1 ∧ gap ptPrime 2 ≥ 1
        ∧ gap ptPrime 3 ≥ 1 ∧ gap ptPrime 4 ≥ 1)
    ∧ (gap ptPrimeExt 1 ≤ 6 ∧ gap ptPrimeExt 2 ≤ 6
        ∧ gap ptPrimeExt 3 ≤ 6 ∧ gap ptPrimeExt 4 ≤ 6
        ∧ gap ptPrimeExt 5 ≤ 6 ∧ gap ptPrimeExt 6 ≤ 6
        ∧ gap ptPrimeExt 7 ≤ 6 ∧ gap ptPrimeExt 8 ≤ 6
        ∧ gap ptPrimeExt 9 ≤ 6 ∧ gap ptPrimeExt 10 ≤ 6)
    -- (H) higher moments at N = 4
    ∧ (sqSumGaps 4 = 25)
    ∧ (cubeSumGaps 4 = 81) :=
  ⟨sum_gap_telescope,
   master_pt_telescope,
   conservation_N1, conservation_N2, conservation_N3, conservation_N4,
   conservation_at_primorial_1,
   conservation_at_primorial_2,
   conservation_at_primorial_3,
   activePrimesSum_eq_muStar,
   activeGapsSum_eq_telescope,
   cumulativeExt_N5, cumulativeExt_N6, cumulativeExt_N7,
   cumulativeExt_N8, cumulativeExt_N9, cumulativeExt_N10,
   gap_pos_small_n,
   gap_max_le_six,
   sqSumGaps_four, cubeSumGaps_four⟩

/-! ### 5. Unified provenance witness

The following theorem makes explicit that items (B)–(F) of the master
summary are all instances of a single underlying identity: the generic
`sum_gap_telescope`. Concretely, every cumulative value listed in the
master summary equals `p (N+1) - p 1` for the appropriate sequence
`p ∈ {ptPrime, ptPrimeExt}` and the appropriate `N`. -/

/-- **Unified provenance witness.** Each of the explicit cumulative
    values bundled in `gap_sum_master_summary` is an instance of the
    generic telescoping identity `master_telescope_identity`, applied
    to either `ptPrime` (for `N ∈ {1, 2, 3, 4}`) or `ptPrimeExt` (for
    `N ∈ {5, …, 10}`). -/
theorem gap_sum_master_provenance :
    -- ptPrime instances
    (∑ n ∈ Ico 1 2, gap ptPrime n = ptPrime 2 - ptPrime 1)
    ∧ (∑ n ∈ Ico 1 3, gap ptPrime n = ptPrime 3 - ptPrime 1)
    ∧ (∑ n ∈ Ico 1 4, gap ptPrime n = ptPrime 4 - ptPrime 1)
    ∧ (∑ n ∈ Ico 1 5, gap ptPrime n = ptPrime 5 - ptPrime 1)
    -- ptPrimeExt instances
    ∧ (∑ n ∈ Ico 1 6,  gap ptPrimeExt n = ptPrimeExt 6  - ptPrimeExt 1)
    ∧ (∑ n ∈ Ico 1 7,  gap ptPrimeExt n = ptPrimeExt 7  - ptPrimeExt 1)
    ∧ (∑ n ∈ Ico 1 8,  gap ptPrimeExt n = ptPrimeExt 8  - ptPrimeExt 1)
    ∧ (∑ n ∈ Ico 1 9,  gap ptPrimeExt n = ptPrimeExt 9  - ptPrimeExt 1)
    ∧ (∑ n ∈ Ico 1 10, gap ptPrimeExt n = ptPrimeExt 10 - ptPrimeExt 1)
    ∧ (∑ n ∈ Ico 1 11, gap ptPrimeExt n = ptPrimeExt 11 - ptPrimeExt 1) :=
  ⟨by simpa using master_telescope_identity ptPrime 1,
   by simpa using master_telescope_identity ptPrime 2,
   by simpa using master_telescope_identity ptPrime 3,
   by simpa using master_telescope_identity ptPrime 4,
   by simpa using master_telescope_identity ptPrimeExt 5,
   by simpa using master_telescope_identity ptPrimeExt 6,
   by simpa using master_telescope_identity ptPrimeExt 7,
   by simpa using master_telescope_identity ptPrimeExt 8,
   by simpa using master_telescope_identity ptPrimeExt 9,
   by simpa using master_telescope_identity ptPrimeExt 10⟩

end PT.Conservation
