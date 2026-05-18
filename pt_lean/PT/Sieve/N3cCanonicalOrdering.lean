/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import Mathlib.Data.Nat.Prime.Nth
import Mathlib.Data.Nat.Prime.Infinite
import Mathlib.Data.Nat.PrimeFin
import Mathlib.Data.Nat.Nth
import Mathlib.Tactic

/-!
# N3c — Canonical Ordering of Primes: `2 < 3 < 5 < 7 < ⋯`

**Statement (paper-level, Ch02 §"Canonical ordering").**
Among all orderings of `Nat.Prime`, the *canonical* one — the standard
order `<` on `ℕ` restricted to the primes — is the **unique** ordering
preserving the **DPI monotonicity property**: every divisibility-prime
inference (DPI), namely the fact that "the n-th prime grows" with `n`,
is monotone in the standard order and only in that order.

Operationalised: there is a *unique* strictly monotone enumeration
`f : ℕ → ℕ` of the (infinite) set of primes, namely
`f = Nat.nth Nat.Prime` (a/k/a "the n-th prime"); see Mathlib's
`Nat.nth` and `Nat.Prime.Nth`. The first values are exactly the PT
prime alphabet:

* `nth Prime 0 = 2`
* `nth Prime 1 = 3`
* `nth Prime 2 = 5` ← first active prime in PT
* `nth Prime 3 = 7` ← last active prime in PT
* `nth Prime 4 = 11` ← first inactive prime (γ_p < s = 1/2)

This file formalises:

1. **Strict monotonicity** `nth Prime` is strictly monotone (since the
   prime set is infinite).
2. **Uniqueness**: any strictly monotone bijection `g : ℕ → {p : ℕ | p.Prime}`
   (equivalently, any strictly monotone enumeration of primes) is equal
   to `nth Prime` componentwise.
3. **Concrete values** for the first five primes (used elsewhere as the
   PT prime alphabet `{2, 3, 5, 7, 11}`).
4. **Range characterisation**: `nth Prime n` is itself a prime, and every
   prime arises this way.

## Reference

Monograph Chapter 2, §"Canonical ordering", `\label{thm:N3c}`.
Audit row #5 (MEDIUM, ~3 sessions): "DPI = multiplicative divergence
monotonicity (Mathlib import)."

## Strategy

We piggyback on Mathlib's `Nat.nth p` (the `n`-th value where `p` holds),
which is strictly monotone whenever the predicate is satisfied infinitely
often. The PT canonical ordering of primes is then *definitionally*
`Nat.nth Nat.Prime`, and the uniqueness is the standard uniqueness of
order embeddings between two well-ordered infinite subsets of `ℕ`.

## Mathlib lemmas used

* `Nat.infinite_setOf_prime` — the prime set is infinite.
* `Nat.nth_strictMono` — `nth` is strictly monotone for infinite predicates.
* `Nat.nth_mem_of_infinite` — `nth p n` satisfies `p`.
* `Nat.range_nth_of_infinite` — `nth p` is surjective onto `{x | p x}`.
* `Nat.nth_prime_*_eq_*` — concrete small values.
-/

namespace PT.Sieve

open Nat

/-! ### Canonical PT prime enumeration -/

/-- The PT canonical prime enumeration `pₙ`: the `n`-th prime in the
    standard order `2, 3, 5, 7, 11, …`. This is just `Nat.nth Nat.Prime`. -/
noncomputable def ptPrime (n : ℕ) : ℕ := Nat.nth Nat.Prime n

/-- `pt_prime 0 = 2`. -/
@[simp] theorem ptPrime_zero : ptPrime 0 = 2 := by
  unfold ptPrime; exact_mod_cast Nat.nth_prime_zero_eq_two

/-- `pt_prime 1 = 3`. -/
@[simp] theorem ptPrime_one : ptPrime 1 = 3 := by
  unfold ptPrime; exact Nat.nth_prime_one_eq_three

/-- `pt_prime 2 = 5` (first PT active prime). -/
@[simp] theorem ptPrime_two : ptPrime 2 = 5 := by
  unfold ptPrime; exact Nat.nth_prime_two_eq_five

/-- `pt_prime 3 = 7` (last PT active prime). -/
@[simp] theorem ptPrime_three : ptPrime 3 = 7 := by
  unfold ptPrime; exact Nat.nth_prime_three_eq_seven

/-- `pt_prime 4 = 11` (first PT inactive prime, `γ_{11} < 1/2`). -/
@[simp] theorem ptPrime_four : ptPrime 4 = 11 := by
  unfold ptPrime; exact Nat.nth_prime_four_eq_eleven

/-! ### Strict monotonicity of the canonical ordering -/

/-- **DPI monotonicity (existence).** The PT canonical prime enumeration
    is **strictly monotone**: `n < m → pₙ < pₘ`. This is the positive
    statement of N3c: the standard order on `ℕ` restricted to primes is
    strictly increasing. -/
theorem ptPrime_strictMono : StrictMono ptPrime := by
  unfold ptPrime
  exact Nat.nth_strictMono Nat.infinite_setOf_prime

/-- Explicit form: `2 = p₀ < p₁ = 3`. -/
theorem ptPrime_zero_lt_one : ptPrime 0 < ptPrime 1 :=
  ptPrime_strictMono (by norm_num : (0 : ℕ) < 1)

/-- `3 = p₁ < p₂ = 5`. -/
theorem ptPrime_one_lt_two : ptPrime 1 < ptPrime 2 :=
  ptPrime_strictMono (by norm_num : (1 : ℕ) < 2)

/-- `5 = p₂ < p₃ = 7` (gap between PT active primes). -/
theorem ptPrime_two_lt_three : ptPrime 2 < ptPrime 3 :=
  ptPrime_strictMono (by norm_num : (2 : ℕ) < 3)

/-- `7 = p₃ < p₄ = 11` (active/inactive boundary). -/
theorem ptPrime_three_lt_four : ptPrime 3 < ptPrime 4 :=
  ptPrime_strictMono (by norm_num : (3 : ℕ) < 4)

/-- **Canonical cascade `2 < 3 < 5 < 7 < 11`**, the first five values of
    the PT prime alphabet, in strict order. -/
theorem ptPrime_cascade :
    ptPrime 0 < ptPrime 1 ∧ ptPrime 1 < ptPrime 2
      ∧ ptPrime 2 < ptPrime 3 ∧ ptPrime 3 < ptPrime 4 :=
  ⟨ptPrime_zero_lt_one, ptPrime_one_lt_two,
   ptPrime_two_lt_three, ptPrime_three_lt_four⟩

/-! ### Membership: each `pₙ` is actually prime -/

/-- Every `pₙ` is a prime number. -/
theorem ptPrime_isPrime (n : ℕ) : Nat.Prime (ptPrime n) := by
  unfold ptPrime
  exact Nat.nth_mem_of_infinite Nat.infinite_setOf_prime n

/-! ### Surjectivity: every prime is `pₙ` for some `n` -/

/-- **Range characterisation.** Every prime is `pₙ` for exactly one `n`:
    the canonical enumeration is *bijective* onto the prime set. -/
theorem ptPrime_surjective_on_primes (q : ℕ) (hq : Nat.Prime q) :
    ∃ n : ℕ, ptPrime n = q := by
  -- `Nat.range_nth_of_infinite`: `Set.range (Nat.nth p) = {x | p x}`
  -- when `{x | p x}` is infinite.
  have hrange : Set.range (Nat.nth Nat.Prime) = { x | Nat.Prime x } :=
    Nat.range_nth_of_infinite Nat.infinite_setOf_prime
  have : q ∈ Set.range (Nat.nth Nat.Prime) := by
    rw [hrange]; exact hq
  obtain ⟨n, hn⟩ := this
  exact ⟨n, hn⟩

/-! ### Uniqueness of the canonical ordering

The PT claim "the order `2 < 3 < 5 < ⋯` is the unique DPI-monotone
ordering of primes" is operationalised as: any strictly monotone function
`f : ℕ → ℕ` whose range is exactly the prime set must equal `ptPrime`. -/

/-- **Uniqueness of the canonical ordering (N3c headline).**
    Any strictly monotone function `f : ℕ → ℕ` whose range is exactly the
    set of primes coincides pointwise with the PT canonical enumeration
    `ptPrime = Nat.nth Nat.Prime`.

    Intuition: two strictly monotone bijections of `ℕ` onto the same
    well-ordered infinite subset of `ℕ` must agree (by induction on `n`:
    both must send `0` to the least prime, then `1` to the next, etc.).
    Mathlib's `Nat.nth_eq_orderIsoOfNat` makes this precise, but we
    derive the simpler pointwise statement directly. -/
theorem ptPrime_unique_strictMono_enum
    (f : ℕ → ℕ)
    (hmono : StrictMono f)
    (hrange : Set.range f = { x | Nat.Prime x }) :
    ∀ n : ℕ, f n = ptPrime n := by
  -- Strategy: both `f` and `ptPrime` are strictly monotone with the
  -- same range (the primes). By induction on `n`, the value `f n` is
  -- forced to be the (n+1)-th smallest prime.
  intro n
  induction n using Nat.strong_induction_on with
  | _ n ih =>
    -- Step 1: show `f n ∈ {primes}`.
    have hfn_prime : Nat.Prime (f n) := by
      have hmem : f n ∈ Set.range f := ⟨n, rfl⟩
      rw [hrange] at hmem
      exact hmem
    -- Step 2: show `ptPrime n` is prime.
    have hptn_prime : Nat.Prime (ptPrime n) := ptPrime_isPrime n
    -- Step 3: show all `f k` for `k < n` equal `ptPrime k` (induction).
    have hbelow : ∀ k < n, f k = ptPrime k := fun k hk => ih k hk
    -- Step 4: use the order-theoretic characterisation of `nth p n`:
    -- it is the smallest prime strictly greater than `nth p k` for all `k < n`.
    -- We rely on `Nat.nth_set_eq_setOf` / `Nat.nth_set` (`Nat.nth p n`
    -- is the minimum of the set of `q` with `p q` and `q ≠ nth p k`
    -- for all `k < n`).
    -- For brevity, we use antisymmetry via two inequalities.
    have hle1 : ptPrime n ≤ f n := by
      -- `ptPrime n = nth Prime n` is the minimum of primes not in
      -- `{ptPrime k | k < n} = {f k | k < n}`. We show `f n` is such
      -- a prime: it's prime (hfn_prime) and `f n > f k = ptPrime k`
      -- for `k < n` (strict monotonicity of `f`).
      -- Therefore `f n ≥ nth Prime n = ptPrime n`.
      by_contra hlt
      push Not at hlt
      -- f n < ptPrime n. By surjectivity, f n = ptPrime m for some m.
      obtain ⟨m, hm⟩ := ptPrime_surjective_on_primes (f n) hfn_prime
      -- Then ptPrime m < ptPrime n, hence m < n by strict monotonicity
      -- (contrapositive).
      have hm_lt_n : m < n := by
        have : ptPrime m < ptPrime n := by rw [hm]; exact hlt
        by_contra hge
        push Not at hge
        have : ptPrime n ≤ ptPrime m := ptPrime_strictMono.le_iff_le.mpr hge
        omega
      -- But then f m = ptPrime m = f n, contradicting injectivity of f.
      have hfm : f m = ptPrime m := ih m hm_lt_n
      have hfm_eq_fn : f m = f n := by rw [hfm, hm]
      have : m = n := hmono.injective hfm_eq_fn
      omega
    have hle2 : f n ≤ ptPrime n := by
      by_contra hlt
      push Not at hlt
      -- ptPrime n < f n. By surjectivity (of f), ptPrime n = f k for some k.
      have hptn_in_range : ptPrime n ∈ Set.range f := by
        rw [hrange]; exact hptn_prime
      obtain ⟨k, hk⟩ := hptn_in_range
      have hfk_lt_fn : f k < f n := by rw [hk]; exact hlt
      have hk_lt_n : k < n := hmono.lt_iff_lt.mp hfk_lt_fn
      -- f k = ptPrime k by IH, so ptPrime k = ptPrime n.
      have hfk_eq : f k = ptPrime k := ih k hk_lt_n
      have : ptPrime k = ptPrime n := by rw [← hfk_eq]; exact hk
      have : k = n := ptPrime_strictMono.injective this
      omega
    exact (le_antisymm hle2 hle1)

/-! ### N3c headline -/

/-- **N3c — Canonical ordering of primes.** Aggregated statement:

    1. The standard order on `ℕ` restricted to primes is strictly monotone
       (existence of canonical ordering).
    2. Every prime is enumerated exactly once (`ptPrime` is bijective onto
       the prime set).
    3. The first five values are `2, 3, 5, 7, 11`.
    4. Any strictly monotone enumeration of the primes coincides with
       `ptPrime` (uniqueness).

    This is the "DPI monotonicity preservation" claim of N3c, in its
    operational form. -/
theorem N3c_canonical_ordering :
    StrictMono ptPrime
    ∧ (∀ q : ℕ, Nat.Prime q → ∃ n : ℕ, ptPrime n = q)
    ∧ ptPrime 0 = 2
    ∧ ptPrime 1 = 3
    ∧ ptPrime 2 = 5
    ∧ ptPrime 3 = 7
    ∧ ptPrime 4 = 11
    ∧ (∀ f : ℕ → ℕ, StrictMono f → Set.range f = { x | Nat.Prime x }
        → ∀ n, f n = ptPrime n) :=
  ⟨ptPrime_strictMono,
   ptPrime_surjective_on_primes,
   ptPrime_zero, ptPrime_one, ptPrime_two, ptPrime_three, ptPrime_four,
   ptPrime_unique_strictMono_enum⟩

end PT.Sieve
