/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Conservation.ConservationID
import PT.Conservation.ConservationIDExtensions
import PT.Conservation.ConservationIDPrimorial
import PT.Conservation.ConservationActivePrimes
import PT.Conservation.CumulativeBoundsExtended
import PT.Conservation.GapBoundedBelow
import PT.Conservation.GapMaxBound
import PT.Conservation.GapParityDecomposition
import PT.Conservation.GapSumMasterIdentity
import PT.Conservation.GapDistributionMomentsSummary
import PT.Sieve.N3aFactorisation
import Mathlib.Tactic

/-!
# PT prime structural theorem — master synthesis aggregator

This file is a **pure synthesis aggregator** for the PT prime sequence
(`ptPrime` / `ptPrimeExt`). It does not introduce any new computational
content: every constituent identity is already proved in one of the
imported modules. The point of this module is to gather, in five
canonical bundles, the complete structural characterisation of the PT
prime sequence and its gap dynamics.

## What is bundled

The five master bundles are:

1. **`master_pt_prime_seq`** — the explicit values of `ptPrime n` for
   `n ∈ {1, …, 5}` and `ptPrimeExt n` for `n ∈ {1, …, 11}`, together
   with the coherence statement that the two sequences agree on the
   common range `{1, …, 5}`.

2. **`master_primorial_structure`** — the explicit values of the small
   primorials `primorial k` for `k ∈ {0, …, 5}`, together with the
   factorisation `primorial 3 = 2 · 3 · 5 = 30` and the extended values
   `primorialExt k` for `k ∈ {6, …, 11}`.

3. **`master_gap_properties`** — the ten gap values
   `(g_1, …, g_10) = (1, 2, 2, 4, 2, 4, 2, 4, 6, 2)`, the lower bound
   `g_n ≥ 1` on `{1, …, 4}`, the upper bound `g_n ≤ 6` on `{1, …, 10}`,
   the parity statement (`g_1` is the unique odd gap among the first
   four), and the strict maximum at `n = 9`.

4. **`master_conservation_identities`** — the generic telescoping
   identity `∀ p N, ∑ gap p n = p (N+1) - p 1`, the PT specialisation
   `∑ gap ptPrime n = ptPrime (N+1) - 2`, the ten explicit cumulative
   values `S_N` for `N ∈ {1, …, 10}`, and the active-prime sum identity
   `3 + 5 + 7 = 15 = μ*`.

5. **`master_arithmetic_invariants`** — the distributional moments
   `(μ_N, σ²_N, m_{3,N}, m_{4,N})` for `N ∈ {5, …, 10}` together with
   the anchor pair `(μ_4, σ²_4) = (3, 10/9)` on `ptPrime`.

6. **`pt_prime_structural_summary`** — the single citeable theorem
   bundling all five master bundles.

## Unified provenance

All bundles are pure re-exports / conjunctive packagings. No new proof
is performed beyond `refine ⟨…⟩` / direct reuse of the proved bundles.
The single underlying analytic ingredient is `sum_gap_telescope`
(`ConservationID`); everything else is finite case analysis.

## Reference

* Monograph Chapter 3, §3.1, `\label{thm:conservation-id}`.
* M1 article (PT_ARTICLES/PT_MATHEMATICS/M1), Conservation identity.
* This file is the canonical Tier 1 entry point for any external
  reference to "the PT prime structural characterisation".
-/

namespace PT.Conservation

open Finset

/-! ### 1. Master prime sequence -/

/-- **Master PT prime sequence bundle.** The PT prime sequence is fully
    specified on its concrete range by:

    * `ptPrime n` for `n ∈ {1, …, 5}` taking the values `(2, 3, 5, 7, 11)`;
    * `ptPrimeExt n` for `n ∈ {1, …, 11}` taking the values
      `(2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31)`;
    * coherence: the two sequences agree on the common range `{1, …, 5}`. -/
theorem master_pt_prime_seq :
    -- ptPrime values (n = 1..5)
    (ptPrime 1 = 2)
    ∧ (ptPrime 2 = 3)
    ∧ (ptPrime 3 = 5)
    ∧ (ptPrime 4 = 7)
    ∧ (ptPrime 5 = 11)
    -- ptPrimeExt values (n = 1..11)
    ∧ (ptPrimeExt 1  = 2)
    ∧ (ptPrimeExt 2  = 3)
    ∧ (ptPrimeExt 3  = 5)
    ∧ (ptPrimeExt 4  = 7)
    ∧ (ptPrimeExt 5  = 11)
    ∧ (ptPrimeExt 6  = 13)
    ∧ (ptPrimeExt 7  = 17)
    ∧ (ptPrimeExt 8  = 19)
    ∧ (ptPrimeExt 9  = 23)
    ∧ (ptPrimeExt 10 = 29)
    ∧ (ptPrimeExt 11 = 31)
    -- Coherence on the common range
    ∧ (ptPrimeExt 1 = ptPrime 1
        ∧ ptPrimeExt 2 = ptPrime 2
        ∧ ptPrimeExt 3 = ptPrime 3
        ∧ ptPrimeExt 4 = ptPrime 4
        ∧ ptPrimeExt 5 = ptPrime 5) :=
  ⟨ptPrime_one, ptPrime_two, ptPrime_three, ptPrime_four, ptPrime_five,
   ptPrimeExt_one,  ptPrimeExt_two,    ptPrimeExt_three,
   ptPrimeExt_four, ptPrimeExt_five,   ptPrimeExt_six,
   ptPrimeExt_seven, ptPrimeExt_eight, ptPrimeExt_nine,
   ptPrimeExt_ten,   ptPrimeExt_eleven,
   ptPrimeExt_eq_ptPrime_small⟩

/-! ### 2. Master primorial structure -/

/-- **Master primorial structure bundle.** The PT primorials are fully
    specified on the concrete range by:

    * Small primorials `primorial k` for `k ∈ {0, …, 5}`, taking values
      `(1, 2, 6, 30, 210, 2310)`;
    * Extended primorials `primorialExt k` for `k ∈ {6, …, 11}`, taking
      values `(30030, 510510, 9699690, 223092870, 6469693230,
      200560490130)`;
    * The structural identity `primorial 3 = 2 · 3 · 5 = 30`
      (the cascade primorial). -/
theorem master_primorial_structure :
    -- small primorials (k = 0..5)
    (primorial 0 = 1)
    ∧ (primorial 1 = 2)
    ∧ (primorial 2 = 6)
    ∧ (primorial 3 = 30)
    ∧ (primorial 4 = 210)
    ∧ (primorial 5 = 2310)
    -- extended primorials (k = 6..11)
    ∧ (primorialExt 6  = 30030)
    ∧ (primorialExt 7  = 510510)
    ∧ (primorialExt 8  = 9699690)
    ∧ (primorialExt 9  = 223092870)
    ∧ (primorialExt 10 = 6469693230)
    ∧ (primorialExt 11 = 200560490130)
    -- cascade-primorial factorisation
    ∧ (primorial 3 = 2 * 3 * 5)
    -- natural-number factorisations (squarefree primorials)
    ∧ ((6 : ℕ) = 2 * 3)
    ∧ ((30 : ℕ) = 2 * 3 * 5)
    ∧ ((210 : ℕ) = 2 * 3 * 5 * 7)
    ∧ ((2310 : ℕ) = 2 * 3 * 5 * 7 * 11) :=
  ⟨primorial_0, primorial_1, primorial_2, primorial_3,
   primorial_4, primorial_5,
   primorialExt_6, primorialExt_7, primorialExt_8,
   primorialExt_9, primorialExt_10, primorialExt_11,
   cascade_primorial_3,
   PT.Sieve.nat_factorise_primorial_2,
   PT.Sieve.nat_factorise_primorial_3,
   PT.Sieve.nat_factorise_primorial_4,
   PT.Sieve.nat_factorise_primorial_5⟩

/-! ### 3. Master gap properties -/

/-- **Master gap properties bundle.** The PT gap sequence
    `g_n = ptPrime(Ext)(n+1) - ptPrime(Ext)(n)` is fully characterised on
    its concrete range by:

    * Exact values `(g_1, …, g_10) = (1, 2, 2, 4, 2, 4, 2, 4, 6, 2)`.
    * Lower bound `g_n ≥ 1` for `n ∈ {1, …, 4}` (strict positivity).
    * Upper bound `g_n ≤ 6` for `n ∈ {1, …, 10}` (global maximum).
    * Parity: `g_1 = 1` is the unique odd gap among `(g_1, g_2, g_3, g_4)`;
      `g_2, g_3, g_4` are all even.
    * Strict maximum: `g_9 = 6` is strictly larger than every other
      `g_k` for `k ∈ {1, …, 8} ∪ {10}`. -/
theorem master_gap_properties :
    -- (i) exact values on n ∈ {1, …, 10} (via ptPrimeExt)
    (gap ptPrimeExt 1  = 1
      ∧ gap ptPrimeExt 2  = 2
      ∧ gap ptPrimeExt 3  = 2
      ∧ gap ptPrimeExt 4  = 4
      ∧ gap ptPrimeExt 5  = 2
      ∧ gap ptPrimeExt 6  = 4
      ∧ gap ptPrimeExt 7  = 2
      ∧ gap ptPrimeExt 8  = 4
      ∧ gap ptPrimeExt 9  = 6
      ∧ gap ptPrimeExt 10 = 2)
    -- (ii) lower bound g_n ≥ 1 on n ∈ {1, …, 4} (via ptPrime)
    ∧ (gap ptPrime 1 ≥ 1 ∧ gap ptPrime 2 ≥ 1
        ∧ gap ptPrime 3 ≥ 1 ∧ gap ptPrime 4 ≥ 1)
    -- (iii) upper bound g_n ≤ 6 on n ∈ {1, …, 10}
    ∧ (gap ptPrimeExt 1 ≤ 6 ∧ gap ptPrimeExt 2 ≤ 6
        ∧ gap ptPrimeExt 3 ≤ 6 ∧ gap ptPrimeExt 4 ≤ 6
        ∧ gap ptPrimeExt 5 ≤ 6 ∧ gap ptPrimeExt 6 ≤ 6
        ∧ gap ptPrimeExt 7 ≤ 6 ∧ gap ptPrimeExt 8 ≤ 6
        ∧ gap ptPrimeExt 9 ≤ 6 ∧ gap ptPrimeExt 10 ≤ 6)
    -- (iv) parity: only g_1 is odd among (g_1, …, g_4)
    ∧ (gap ptPrime 1 % 2 = 1
        ∧ gap ptPrime 2 % 2 = 0
        ∧ gap ptPrime 3 % 2 = 0
        ∧ gap ptPrime 4 % 2 = 0)
    -- (v) strict max at n = 9
    ∧ (gap ptPrimeExt 9 = 6
        ∧ gap ptPrimeExt 9 > gap ptPrimeExt 1
        ∧ gap ptPrimeExt 9 > gap ptPrimeExt 2
        ∧ gap ptPrimeExt 9 > gap ptPrimeExt 3
        ∧ gap ptPrimeExt 9 > gap ptPrimeExt 4
        ∧ gap ptPrimeExt 9 > gap ptPrimeExt 5
        ∧ gap ptPrimeExt 9 > gap ptPrimeExt 6
        ∧ gap ptPrimeExt 9 > gap ptPrimeExt 7
        ∧ gap ptPrimeExt 9 > gap ptPrimeExt 8
        ∧ gap ptPrimeExt 9 > gap ptPrimeExt 10) :=
  ⟨gap_values_tuple,
   gap_pos_small_n,
   gap_max_le_six,
   only_g1_is_odd,
   gap_strict_max_at_nine⟩

/-! ### 4. Master conservation identities -/

/-- **Master conservation identities bundle.** The full conservation
    content of `PT.Conservation` is bundled here:

    * **(A) Generic telescoping**: `∀ p N, ∑ gap p n = p (N+1) - p 1`.
    * **(B) PT specialisation**: `∑ gap ptPrime n = ptPrime (N+1) - 2`.
    * **(C) PT-Ext specialisation**: same identity for `ptPrimeExt`.
    * **(D) Ten explicit cumulative values** for `N ∈ {1, …, 10}`,
      taking values `(1, 3, 5, 9, 11, 15, 17, 21, 27, 29)`.
    * **(E) Active-prime sum** `3 + 5 + 7 = 15 = μ*`.
    * **(F) Active-gaps sum** equals `ptPrime 4 - ptPrime 1 = 5`. -/
theorem master_conservation_identities :
    -- (A) Generic telescoping
    (∀ (p : ℕ → ℤ) (N : ℕ),
        ∑ n ∈ Ico 1 (N + 1), gap p n = p (N + 1) - p 1)
    -- (B) PT specialisation
    ∧ (∀ N : ℕ, ∑ n ∈ Ico 1 (N + 1), gap ptPrime n = ptPrime (N + 1) - 2)
    -- (C) PT-Ext specialisation
    ∧ (∀ N : ℕ, ∑ n ∈ Ico 1 (N + 1), gap ptPrimeExt n = ptPrimeExt (N + 1) - 2)
    -- (D) ten explicit cumulative values
    ∧ (∑ n ∈ Ico 1 2,  gap ptPrime n = 1)
    ∧ (∑ n ∈ Ico 1 3,  gap ptPrime n = 3)
    ∧ (∑ n ∈ Ico 1 4,  gap ptPrime n = 5)
    ∧ (∑ n ∈ Ico 1 5,  gap ptPrime n = 9)
    ∧ (∑ n ∈ Ico 1 6,  gap ptPrimeExt n = 11)
    ∧ (∑ n ∈ Ico 1 7,  gap ptPrimeExt n = 15)
    ∧ (∑ n ∈ Ico 1 8,  gap ptPrimeExt n = 17)
    ∧ (∑ n ∈ Ico 1 9,  gap ptPrimeExt n = 21)
    ∧ (∑ n ∈ Ico 1 10, gap ptPrimeExt n = 27)
    ∧ (∑ n ∈ Ico 1 11, gap ptPrimeExt n = 29)
    -- (E) active-prime sum = 15
    ∧ (activePrimesSum = 15)
    -- (F) active-gaps sum = ptPrime 4 - ptPrime 1
    ∧ (activeGapsSum = ptPrime 4 - ptPrime 1) :=
  ⟨sum_gap_telescope,
   master_pt_telescope,
   master_ptExt_telescope,
   conservation_N1, conservation_N2, conservation_N3, conservation_N4,
   cumulativeExt_N5,  cumulativeExt_N6,  cumulativeExt_N7,
   cumulativeExt_N8,  cumulativeExt_N9,  cumulativeExt_N10,
   activePrimesSum_eq_muStar,
   activeGapsSum_eq_telescope⟩

/-! ### 5. Master arithmetic invariants -/

/-- **Master arithmetic invariants bundle.** The first four distributional
    moments of the gap-weighted index distribution `n ↦ g_n / S_N` are
    bundled here:

    * **Anchor (`N = 4`, on `ptPrime`)**: `(μ_4, σ²_4) = (3, 10/9)`,
      contrasted with the classical value-based variance `19/16`.
    * **Extended table (`N ∈ {5, …, 10}`, on `ptPrimeExt`)**: the six
      four-tuples `(μ_N, σ²_N, m_{3,N}, m_{4,N})`. -/
theorem master_arithmetic_invariants :
    -- anchor at N = 4 on ptPrime
    (gapsDistMeanQ 4 = 3)
    ∧ (gapsDistVarianceQ 4 = 10 / 9)
    ∧ (varianceGapsTo 4 = 19 / 16)
    ∧ (gapsDistVarianceQ 4 < varianceGapsTo 4)
    -- six extended four-tuples on ptPrimeExt
    ∧ (gapsDistMomentsExt_tuple 5  =
        ( 37 / 11, 182 / 121, -1038 / 1331, 70754 / 14641 ))
    ∧ (gapsDistMomentsExt_tuple 6  =
        ( 61 / 15, 554 / 225, -4138 / 3375, 208034 / 16875 ))
    ∧ (gapsDistMomentsExt_tuple 7  =
        ( 75 / 17, 886 / 289, -6522 / 4913, 1604890 / 83521 ))
    ∧ (gapsDistMomentsExt_tuple 8  =
        ( 107 / 21, 1970 / 441, -16238 / 9261, 2540354 / 64827 ))
    ∧ (gapsDistMomentsExt_tuple 9  =
        ( 161 / 27, 4454 / 729, -92342 / 19683, 12445406 / 177147 ))
    ∧ (gapsDistMomentsExt_tuple 10 =
        ( 181 / 29, 5664 / 841, -133584 / 24389, 61313616 / 707281 )) :=
  gapsDistMoments_full_summary

/-! ### 6. PT prime structural summary — the master aggregator -/

/-- **PT prime structural summary.** The single citeable theorem bundling
    the five master synthesis bundles:

    1. `master_pt_prime_seq` — the prime values themselves
       (`ptPrime n` for `n ∈ {1, …, 5}`, `ptPrimeExt n` for
       `n ∈ {1, …, 11}`, with coherence on the overlap).
    2. `master_primorial_structure` — the primorials and their
       factorisation, including the cascade primorial `30 = 2 · 3 · 5`.
    3. `master_gap_properties` — the gap values, bounds (`1 ≤ g_n ≤ 6`),
       parity, and strict maximum at `n = 9`.
    4. `master_conservation_identities` — generic telescoping, PT
       specialisation, ten explicit cumulative values, and the
       active-prime sum `15 = μ*`.
    5. `master_arithmetic_invariants` — the first four distributional
       moments for `N ∈ {4, …, 10}`.

    This is the Tier-1 entry point for any external reference to "the
    PT prime structural theorem". No new proof beyond the conjunction
    of the five bundles. -/
theorem pt_prime_structural_summary :
    -- Bundle 1: prime sequence
    (ptPrime 1 = 2 ∧ ptPrime 2 = 3 ∧ ptPrime 3 = 5
      ∧ ptPrime 4 = 7 ∧ ptPrime 5 = 11
      ∧ ptPrimeExt 1 = 2 ∧ ptPrimeExt 2 = 3 ∧ ptPrimeExt 3 = 5
      ∧ ptPrimeExt 4 = 7 ∧ ptPrimeExt 5 = 11
      ∧ ptPrimeExt 6 = 13 ∧ ptPrimeExt 7 = 17 ∧ ptPrimeExt 8 = 19
      ∧ ptPrimeExt 9 = 23 ∧ ptPrimeExt 10 = 29 ∧ ptPrimeExt 11 = 31
      ∧ (ptPrimeExt 1 = ptPrime 1 ∧ ptPrimeExt 2 = ptPrime 2
          ∧ ptPrimeExt 3 = ptPrime 3 ∧ ptPrimeExt 4 = ptPrime 4
          ∧ ptPrimeExt 5 = ptPrime 5))
    -- Bundle 2: primorial structure (cascade-primorial factorisation)
    ∧ (primorial 3 = 30 ∧ primorial 3 = 2 * 3 * 5
        ∧ (30 : ℕ) = 2 * 3 * 5)
    -- Bundle 3: gap properties (key witnesses — full values, max at 9)
    ∧ ((gap ptPrimeExt 1 = 1 ∧ gap ptPrimeExt 2 = 2
          ∧ gap ptPrimeExt 3 = 2 ∧ gap ptPrimeExt 4 = 4
          ∧ gap ptPrimeExt 5 = 2 ∧ gap ptPrimeExt 6 = 4
          ∧ gap ptPrimeExt 7 = 2 ∧ gap ptPrimeExt 8 = 4
          ∧ gap ptPrimeExt 9 = 6 ∧ gap ptPrimeExt 10 = 2)
        ∧ gap ptPrime 1 % 2 = 1
        ∧ gap ptPrimeExt 9 = 6)
    -- Bundle 4: conservation identities (generic + PT + active-prime sum)
    ∧ ((∀ (p : ℕ → ℤ) (N : ℕ),
          ∑ n ∈ Ico 1 (N + 1), gap p n = p (N + 1) - p 1)
        ∧ (∀ N : ℕ, ∑ n ∈ Ico 1 (N + 1), gap ptPrime n = ptPrime (N + 1) - 2)
        ∧ activePrimesSum = 15)
    -- Bundle 5: arithmetic invariants (anchor at N = 4)
    ∧ (gapsDistMeanQ 4 = 3
        ∧ gapsDistVarianceQ 4 = 10 / 9) :=
  ⟨⟨ptPrime_one, ptPrime_two, ptPrime_three, ptPrime_four, ptPrime_five,
    ptPrimeExt_one, ptPrimeExt_two, ptPrimeExt_three,
    ptPrimeExt_four, ptPrimeExt_five, ptPrimeExt_six,
    ptPrimeExt_seven, ptPrimeExt_eight, ptPrimeExt_nine,
    ptPrimeExt_ten, ptPrimeExt_eleven,
    ptPrimeExt_eq_ptPrime_small⟩,
   ⟨cascade_primorial_3_eq_30, cascade_primorial_3,
    PT.Sieve.nat_factorise_primorial_3⟩,
   ⟨gap_values_tuple, gap_one_odd, gapExt_nine⟩,
   ⟨sum_gap_telescope, master_pt_telescope, activePrimesSum_eq_muStar⟩,
   ⟨gapsDistMeanQ_four, gapsDistVarianceQ_four⟩⟩

end PT.Conservation
