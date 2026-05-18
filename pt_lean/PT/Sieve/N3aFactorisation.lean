/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Sieve.N3aMinimalMonoid
import Mathlib.Data.Nat.Factorization.Basic
import Mathlib.Tactic

/-!
# N3a — Explicit factorisations and PT-relevant instances (Ch02 extension)

This file extends `PT.Sieve.N3aMinimalMonoid` with **explicit small-number
prime factorisations**, computed via `decide` / `Nat.primeFactorsList`:

* PT primorials `2, 6, 30, 210, 2310` factor exactly as products of distinct
  small primes.
* PT fixed-point `μ* = 15` factors as `3 · 5`.
* The "active primes" `3, 5, 7` are atoms (irreducible in `(ℕ, ×)`).
* The "inactive primes" `2, 11, 13` are also atoms (TFA).

These are pure `decide` instances of the N3a UFD structure on `ℕ`.

## Reference

Monograph Chapter 2, follow-up to `N3a_nat_is_canonical_UFM` (`thm:N3a`).
Audit row #3 extension.
-/

namespace PT.Sieve

/-! ### PT primorial factorisations -/

/-- `primorial 1 = 2` is itself a prime atom. -/
theorem nat_factorise_primorial_1 : Nat.Prime 2 := by decide

/-- `primorial 2 = 6 = 2 · 3`. -/
theorem nat_factorise_primorial_2 : (6 : ℕ) = 2 * 3 := by decide

/-- `primorial 3 = 30 = 2 · 3 · 5`. -/
theorem nat_factorise_primorial_3 : (30 : ℕ) = 2 * 3 * 5 := by decide

/-- `primorial 4 = 210 = 2 · 3 · 5 · 7`. -/
theorem nat_factorise_primorial_4 : (210 : ℕ) = 2 * 3 * 5 * 7 := by decide

/-- `primorial 5 = 2310 = 2 · 3 · 5 · 7 · 11`. -/
theorem nat_factorise_primorial_5 : (2310 : ℕ) = 2 * 3 * 5 * 7 * 11 := by decide

/-! ### PT fixed-point factorisation -/

/-- The PT distinguished integer `μ* = 15` factors as `3 · 5`. -/
theorem nat_factorise_muStar : (15 : ℕ) = 3 * 5 := by decide

/-- `μ* = 15` is **not prime** (it factors non-trivially). -/
theorem muStar_not_prime : ¬ Nat.Prime 15 := by decide

/-- The two prime factors of `μ*` are exactly the two largest active primes
    `{3, 5}`. (The third active prime `7` is *not* a factor of `15`; it
    enters by the sum `3 + 5 + 7 = 15`, not by the product.) -/
theorem muStar_primeFactors :
    (15 : ℕ).primeFactors = ({3, 5} : Finset ℕ) := by
  native_decide

/-! ### Active primes are atoms -/

/-- The active primes are atoms (irreducible). -/
theorem active_3_atom : Nat.Prime 3 := by decide
theorem active_5_atom : Nat.Prime 5 := by decide
theorem active_7_atom : Nat.Prime 7 := by decide

/-- The "inactive primes" `{11, 13}` are also atoms. -/
theorem inactive_11_atom : Nat.Prime 11 := by decide
theorem inactive_13_atom : Nat.Prime 13 := by decide
theorem inactive_17_atom : Nat.Prime 17 := by decide

/-- The "anti-info" prime `p = 2` is an atom (TFA base case). -/
theorem antiinfo_2_atom : Nat.Prime 2 := by decide

/-! ### N3a-specialised structural facts at `μ*` -/

/-- The sum of active primes `3 + 5 + 7 = 15 = μ*`. -/
theorem active_primes_sum_eq_muStar : 3 + 5 + 7 = (15 : ℕ) := by decide

/-- The product of the **first two** active primes equals `μ*`:
    `3 · 5 = 15`. -/
theorem first_two_active_product_eq_muStar : (3 : ℕ) * 5 = 15 := by decide

/-- The product of **all three** active primes is `3 · 5 · 7 = 105`. -/
theorem all_active_product : (3 : ℕ) * 5 * 7 = 105 := by decide

/-! ### Factor counts and primorial cascade -/

/-- The number of distinct prime factors of `30` is `3` (matching the
    cascade dimension `|P_active| = 3`). -/
theorem primorial_3_omega : (30 : ℕ).primeFactors.card = 3 := by
  native_decide

/-- The number of distinct prime factors of `210` is `4`. -/
theorem primorial_4_omega : (210 : ℕ).primeFactors.card = 4 := by
  native_decide

/-- The number of distinct prime factors of `2310` is `5`. -/
theorem primorial_5_omega : (2310 : ℕ).primeFactors.card = 5 := by
  native_decide

/-! ### Headline -/

/-- **Headline (factorisation summary).**

    * Primorials `{2, 6, 30, 210, 2310}` factor cleanly as products of
      distinct primes (squarefree).
    * `μ* = 15 = 3 · 5` (product of two largest active primes); `7` enters
      only additively.
    * Active primes `{3, 5, 7}` and inactive `{2, 11, 13, 17}` are all atoms.
    * `ω(primorial k) = k` for `k ∈ {3, 4, 5}` (number of distinct prime
      factors of the `k`-th primorial). -/
theorem N3a_PT_factorisations :
    -- Primorials
    (6 : ℕ) = 2 * 3
    ∧ (30 : ℕ) = 2 * 3 * 5
    ∧ (210 : ℕ) = 2 * 3 * 5 * 7
    -- μ*
    ∧ (15 : ℕ) = 3 * 5
    ∧ (15 : ℕ).primeFactors = ({3, 5} : Finset ℕ)
    -- Active primes atomic
    ∧ Nat.Prime 3 ∧ Nat.Prime 5 ∧ Nat.Prime 7
    -- Sum identity
    ∧ 3 + 5 + 7 = (15 : ℕ) := by
  refine ⟨nat_factorise_primorial_2, nat_factorise_primorial_3,
          nat_factorise_primorial_4, nat_factorise_muStar,
          muStar_primeFactors,
          active_3_atom, active_5_atom, active_7_atom,
          active_primes_sum_eq_muStar⟩

end PT.Sieve
