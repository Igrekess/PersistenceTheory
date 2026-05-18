/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Sieve.N2SelfCoherence

/-!
# N2 — Uniqueness: `ℙ` is the unique self-coherent subset of `ℕ_{≥2}`

This file completes `PT.Sieve.N2SelfCoherence` (existence, `S(ℙ) = ℙ`)
with the **uniqueness** half: if `G ⊆ {n | 2 ≤ n}` and `sieveStep G`
equals `G` as a set, then `G = {n | n.Prime}`.

The proof mirrors the three-case argument of the monograph ch02 §N2:

* **Step 1 (`primes_subset_of_self_coherent`).** Any prime `p` lies in
  `S(G) = G`: by primality, no proper divisor of `p` lies in
  `Set.Ici 2`, so the divisor-clause of `sieveStep` is vacuous.

* **Step 2 (`self_coherent_subset_primes`).** If `c ∈ G` were composite,
  take a prime factor `p ∣ c` with `p < c`. By Step 1, `p ∈ G`, so
  `sieveStep G c` fails — contradicting `c ∈ G = {n | sieveStep G n}`.

* **Conclusion (`S_unique_fixed_point`).** Antisymmetry of `⊆`.

## Reference

Monograph chapter 2, §"N2 : Auto-cohérence (unique point fixe)",
`\label{thm:N2}`, three-case proof.
-/

namespace PT.Sieve

/-- **Step 1.** Every prime lies in any self-coherent set `G ⊆ ℕ_{≥2}`. -/
lemma primes_subset_of_self_coherent (G : Set ℕ)
    (h_fixed : {n | sieveStep G n} = G) :
    {n | n.Prime} ⊆ G := by
  intro p hp
  -- `hp : p.Prime`
  have hp_prime : p.Prime := hp
  -- Show `sieveStep G p`, then transport along `h_fixed`.
  have h_step : sieveStep G p := by
    refine ⟨hp_prime.two_le, ?_⟩
    intro g _hg_in hg_gt_one hg_lt_p hg_dvd
    -- A divisor of a prime is `1` or the prime itself; both contradict
    -- `1 < g < p`.
    rcases (Nat.dvd_prime hp_prime).mp hg_dvd with h1 | hgp
    · omega
    · omega
  -- `sieveStep G p` means `p ∈ {n | sieveStep G n} = G`.
  have : p ∈ ({n | sieveStep G n} : Set ℕ) := h_step
  simpa [h_fixed] using this

/-- **Step 2.** A self-coherent set `G ⊆ ℕ_{≥2}` consists of primes only. -/
lemma self_coherent_subset_primes (G : Set ℕ)
    (hG : G ⊆ {n | 2 ≤ n})
    (h_fixed : {n | sieveStep G n} = G) :
    G ⊆ {n | n.Prime} := by
  intro c hc
  -- `c ∈ G`, hence `2 ≤ c` and `sieveStep G c`.
  have hc_ge_two : 2 ≤ c := hG hc
  have hc_step : sieveStep G c := by
    have : c ∈ ({n | sieveStep G n} : Set ℕ) := by
      simpa [h_fixed] using hc
    exact this
  -- Suppose for contradiction `c` is not prime.
  by_contra h_not_prime
  simp only [Set.mem_setOf_eq] at h_not_prime
  -- Since `2 ≤ c` and `c` is not prime, `c` admits a proper divisor `m`
  -- with `1 < m < c`.
  rw [Nat.prime_def_lt] at h_not_prime
  push Not at h_not_prime
  obtain ⟨m, hm_lt_c, hm_dvd, hm_ne_one⟩ := h_not_prime hc_ge_two
  -- `m ≥ 2` (since `m ∣ c`, `c ≥ 2`, so `m ≠ 0`; combined with `m ≠ 1`).
  have hm_ne_zero : m ≠ 0 := by
    intro h
    subst h
    have : c = 0 := Nat.zero_dvd.mp hm_dvd
    omega
  have hm_ge_two : 2 ≤ m := by
    rcases Nat.lt_or_ge m 2 with h | h
    · interval_cases m
      · exact absurd rfl hm_ne_zero
      · exact absurd rfl hm_ne_one
    · exact h
  -- Pick a prime factor `p ∣ m`.
  obtain ⟨p, hp_prime, hp_dvd_m⟩ :=
    Nat.exists_prime_and_dvd (show m ≠ 1 from hm_ne_one)
  have hp_le_m : p ≤ m := Nat.le_of_dvd (by omega) hp_dvd_m
  have hp_lt_c : p < c := lt_of_le_of_lt hp_le_m hm_lt_c
  have hp_dvd_c : p ∣ c := dvd_trans hp_dvd_m hm_dvd
  -- By Step 1, `p ∈ G`.
  have hp_in_G : p ∈ G :=
    primes_subset_of_self_coherent G h_fixed hp_prime
  -- But then `sieveStep G c` forbids this divisor.
  have hp_gt_one : 1 < p := hp_prime.one_lt
  exact hc_step.2 p hp_in_G hp_gt_one hp_lt_c hp_dvd_c

/-- **N2 — Uniqueness.** The set of primes is the *unique* self-coherent
    (`sieveStep G = G`) subset of `ℕ_{≥2}`. -/
theorem S_unique_fixed_point (G : Set ℕ)
    (hG_subset : G ⊆ {n | 2 ≤ n})
    (h_fixed : {n | sieveStep G n} = G) :
    G = {n | n.Prime} := by
  apply Set.Subset.antisymm
  · exact self_coherent_subset_primes G hG_subset h_fixed
  · exact primes_subset_of_self_coherent G h_fixed

end PT.Sieve
