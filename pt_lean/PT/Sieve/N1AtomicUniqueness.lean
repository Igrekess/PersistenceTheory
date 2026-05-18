/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import Mathlib.Data.Nat.Prime.Basic
import Mathlib.Tactic

/-!
# N1 — Atomic uniqueness of primes in `(ℕ_{≥2}, ×)`

**Statement (paper-level, Ch02 §"N1").**
The primes `ℙ ⊆ ℕ_{≥2}` are exactly the **atoms** of the multiplicative
monoid `(ℕ_{≥1}, ×)`: an element `n ≥ 2` is *atomic* (no proper
factorisation `n = a · b` with `a, b ≥ 2`) iff `n` is prime.

Proof: classical, via the elementary characterisation of primes
("a prime has no divisor strictly between `1` and `n`"). Reduces to
`Nat.prime_def_lt` in Mathlib.

## Reference

Monograph chapter 2, §"N1 : Unicité algébrique des nombres premiers",
`\label{thm:N1}`.
-/

namespace PT.Sieve

/-- A natural number `n ≥ 2` is **atomic** in `(ℕ, ×)` if it has no proper
    multiplicative factorisation `n = a · b` with `a, b ≥ 2`. -/
def IsAtomic (n : ℕ) : Prop :=
  2 ≤ n ∧ ∀ a b : ℕ, 2 ≤ a → 2 ≤ b → n ≠ a * b

/-- **N1 — Atomic ⟺ Prime.**
    For `n ≥ 2`, `n` is atomic iff `n` is prime in the usual sense. -/
theorem isAtomic_iff_prime {n : ℕ} (hn : 2 ≤ n) :
    IsAtomic n ↔ n.Prime := by
  constructor
  · -- Atomic → Prime
    rintro ⟨_, hatom⟩
    rw [Nat.prime_def_lt]
    refine ⟨hn, ?_⟩
    intro m hmn hmd
    -- Suppose m ∣ n with m < n. Show m = 1.
    rcases Nat.lt_or_ge m 2 with hm_lt | hm_ge
    · -- m < 2 means m = 0 or m = 1
      interval_cases m
      · -- m = 0: 0 ∣ n implies n = 0, contradicting n ≥ 2
        exact absurd (Nat.zero_dvd.mp hmd) (by omega)
      · rfl
    · -- m ≥ 2 and m ∣ n means n = m * (n/m); we get a non-trivial factorisation
      obtain ⟨k, hk⟩ := hmd
      -- hk : n = m * k
      -- m < n and m ≥ 2 imply k ≥ 2 (since k = 1 gives m = n)
      have hk2 : 2 ≤ k := by
        rcases Nat.lt_or_ge k 2 with hk_lt | hk_ge
        · interval_cases k
          · -- k = 0 forces n = 0
            simp at hk; omega
          · -- k = 1 forces n = m, but m < n
            simp at hk; omega
        · exact hk_ge
      exact absurd hk (hatom m k hm_ge hk2)
  · -- Prime → Atomic
    intro hp
    refine ⟨hn, ?_⟩
    intro a b ha hb hne
    -- hne : n = a * b with a, b ≥ 2; contradicts primality.
    rw [Nat.prime_def_lt] at hp
    -- a ∣ n and 1 < a < n
    have ha_dvd : a ∣ n := ⟨b, hne⟩
    have ha_lt : a < n := by
      have hb_pos : 1 < b := hb
      have ha_pos : 0 < a := by omega
      have hmul : a * 1 < a * b :=
        (Nat.mul_lt_mul_left ha_pos).mpr hb_pos
      rw [Nat.mul_one] at hmul
      rw [hne]; exact hmul
    have : a = 1 := hp.2 a ha_lt ha_dvd
    omega

/-- **N1 — Primes are atoms.** Every prime `p` is atomic. -/
theorem prime_isAtomic {p : ℕ} (hp : p.Prime) : IsAtomic p :=
  (isAtomic_iff_prime hp.two_le).mpr hp

/-- **N1 — Atoms are primes.** Every atom `n ≥ 2` is prime. -/
theorem isAtomic_prime {n : ℕ} (h : IsAtomic n) : n.Prime :=
  (isAtomic_iff_prime h.1).mp h

end PT.Sieve
