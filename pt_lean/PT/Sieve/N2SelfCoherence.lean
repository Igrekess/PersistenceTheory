/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Sieve.N1AtomicUniqueness

/-!
# N2 — Self-coherence: `S(ℙ) = ℙ`

**Statement (paper-level, Ch02 §"N2").**
Let `S(G)` denote the multiplicative sieve on `G ⊆ ℕ_{≥2}`: an element
`n ≥ 2` belongs to `S(G)` iff *no* element `g ∈ G` with `1 < g < n`
divides `n`. Then `ℙ` (the set of primes) is **self-coherent**:
`S(ℙ) = ℙ`. Moreover, `ℙ` is the *unique* self-coherent subset of
`ℕ_{≥2}` (case analysis on `G ⊊ ℙ`, `G ⊋ ℙ`).

This file formalises the *existence* part `S(ℙ) = ℙ` cleanly. The
uniqueness statement (case analysis on `G ⊊ ℙ`, `G ⊋ ℙ`, `G ∦ ℙ`)
is **not yet formalised here**: it is treated informally in the
monograph (ch02 §N2, three-case argument) and tagged
`\leanProvedScope{existence ; unicité informelle}` accordingly.
Closing it would yield a `S_unique_fixed_point` companion lemma.

## Reference

Monograph chapter 2, §"N2 : Auto-cohérence (unique point fixe)",
`\label{thm:N2}`.
-/

namespace PT.Sieve

/-- The multiplicative sieve operator: `n ∈ S G` iff `n ≥ 2` and no
    `g ∈ G` with `1 < g < n` divides `n`. -/
def sieveStep (G : Set ℕ) (n : ℕ) : Prop :=
  2 ≤ n ∧ ∀ g ∈ G, 1 < g → g < n → ¬ (g ∣ n)

/-- **N2 — Existence: `S(ℙ) = ℙ`.**
    A natural `n ≥ 2` belongs to `S(ℙ)` iff `n` is prime. -/
theorem sieveStep_primes_iff (n : ℕ) :
    sieveStep {p | Nat.Prime p} n ↔ n.Prime := by
  constructor
  · -- direction →: n is in S(ℙ) ⇒ n is prime
    rintro ⟨hn2, hnd⟩
    -- If n were composite, write n = a * b with 2 ≤ a, 2 ≤ b
    rw [Nat.prime_def_lt]
    refine ⟨hn2, ?_⟩
    intro m hmn hmd
    rcases Nat.lt_or_ge m 2 with hm_lt | hm_ge
    · interval_cases m
      · exact absurd (Nat.zero_dvd.mp hmd) (by omega)
      · rfl
    · -- m ≥ 2 and m ∣ n with m < n; pick a prime factor p of m
      -- p ≤ m < n, p prime, p ∣ m ∣ n, gives p ∈ ℙ with 1 < p < n and p ∣ n
      obtain ⟨p, hp_prime, hp_dvd⟩ := Nat.exists_prime_and_dvd (show m ≠ 1 by omega)
      have hp_le_m : p ≤ m := Nat.le_of_dvd (by omega) hp_dvd
      have hp_lt_n : p < n := lt_of_le_of_lt hp_le_m hmn
      have hp_dvd_n : p ∣ n := dvd_trans hp_dvd hmd
      have hp_in : p ∈ {q | Nat.Prime q} := hp_prime
      have hp_gt_one : 1 < p := hp_prime.one_lt
      exact absurd hp_dvd_n (hnd p hp_in hp_gt_one hp_lt_n)
  · -- direction ←: n is prime ⇒ n is in S(ℙ)
    intro hp
    refine ⟨hp.two_le, ?_⟩
    intro g _hg_prime hg_gt_one hg_lt_n hg_dvd
    -- g ∣ n with 1 < g < n contradicts primality (Nat.prime_def_lt second part)
    rw [Nat.prime_def_lt] at hp
    have := hp.2 g hg_lt_n hg_dvd
    omega

/-- **N2 — Self-coherence as set-equality.** The primes are exactly the
    fixed-point set of the sieve operator. -/
theorem S_primes_eq_primes :
    {n : ℕ | sieveStep {p | Nat.Prime p} n} = {p | Nat.Prime p} := by
  ext n
  exact sieveStep_primes_iff n

end PT.Sieve
