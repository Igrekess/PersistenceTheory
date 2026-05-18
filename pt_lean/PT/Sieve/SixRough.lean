/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import Mathlib.Tactic

/-!
# 6-Rough Integers — Definition and Residue Characterisation

A natural number is **6-rough** if it is positive and coprime to 6 (i.e. coprime
to both 2 and 3). The 6-rough integers are the candidate primes after the
classical sieve removes multiples of 2 and 3.

This file collects the basic definitions and residue lemmas used by the
sieve-level transfer matrix `T₃` of Persistence Theory. The main theorem T1
("forbidden self-transitions") is proved in `PT.Sieve.T1ForbiddenTransitions`.

## Main definitions

* `PT.Sieve.SixRough` : the 6-rough predicate on `ℕ`
* `PT.Sieve.NextSixRough` : `b` is the smallest 6-rough greater than `a`

## Main results

* `SixRough.mod6` : 6-rough ⟹ residue mod 6 ∈ {1, 5}
* `SixRough.mod3` : 6-rough ⟹ residue mod 3 ∈ {1, 2}
* `nextSixRough_of_mod6_one` : from `1 mod 6`, next 6-rough is at distance 4
* `nextSixRough_of_mod6_five` : from `5 mod 6`, next 6-rough is at distance 2

## Reference

Section 2 of: Y. Senez, *Forbidden transitions and informational structure of
prime residues under modular projection* (PT_ARTICLES/PT_MATHEMATICS/M1).
-/

namespace PT.Sieve

/-- A natural number is **6-rough** if it is positive and coprime to 2 and 3.
    Equivalently, `n % 6 ∈ {1, 5}`. -/
def SixRough (n : ℕ) : Prop :=
  0 < n ∧ n % 2 ≠ 0 ∧ n % 3 ≠ 0

instance : DecidablePred SixRough :=
  fun n => inferInstanceAs (Decidable (0 < n ∧ n % 2 ≠ 0 ∧ n % 3 ≠ 0))

/-- `b` is the smallest 6-rough integer strictly greater than `a`. -/
def NextSixRough (a b : ℕ) : Prop :=
  SixRough a ∧ SixRough b ∧ a < b ∧ ∀ c, a < c → c < b → ¬ SixRough c

/-! ### Residue characterisation -/

/-- A 6-rough number has residue 1 or 5 modulo 6. -/
theorem SixRough.mod6 {n : ℕ} (h : SixRough n) : n % 6 = 1 ∨ n % 6 = 5 := by
  obtain ⟨_, _, _⟩ := h
  omega

/-- A 6-rough number has residue 1 or 2 modulo 3. -/
theorem SixRough.mod3 {n : ℕ} (h : SixRough n) : n % 3 = 1 ∨ n % 3 = 2 := by
  obtain ⟨_, _, _⟩ := h
  omega

/-! ### Gap structure -/

/-- From residue `1 mod 6`, the next 6-rough is at distance 4.
    Intermediate values `a+1` (even), `a+2` (divisible by 3), `a+3` (even) are
    eliminated by the coprimality conditions. -/
theorem nextSixRough_of_mod6_one {a : ℕ} (ha : SixRough a) (h6 : a % 6 = 1) :
    NextSixRough a (a + 4) := by
  refine ⟨ha, ⟨by omega, by omega, by omega⟩, by omega, fun c hac hcb => ?_⟩
  rintro ⟨_, hc2, hc3⟩
  have : c = a + 1 ∨ c = a + 2 ∨ c = a + 3 := by omega
  rcases this with rfl | rfl | rfl <;> omega

/-- From residue `5 mod 6`, the next 6-rough is at distance 2.
    The sole intermediate value `a+1` is divisible by 6. -/
theorem nextSixRough_of_mod6_five {a : ℕ} (ha : SixRough a) (h6 : a % 6 = 5) :
    NextSixRough a (a + 2) := by
  refine ⟨ha, ⟨by omega, by omega, by omega⟩, by omega, fun c hac hcb => ?_⟩
  rintro ⟨_, hc2, _⟩
  have : c = a + 1 := by omega
  subst this; omega

/-! ### Uniqueness of the successor -/

/-- If `a % 6 = 1`, the next 6-rough after `a` is exactly `a + 4`. -/
theorem next_eq_add_four {a b : ℕ} (h : NextSixRough a b) (h6 : a % 6 = 1) :
    b = a + 4 := by
  obtain ⟨ha, hb, hab, hcons⟩ := h
  by_contra hne
  rcases show b < a + 4 ∨ a + 4 < b by omega with hlt | hgt
  · obtain ⟨_, _, _⟩ := hb; omega
  · exact absurd (nextSixRough_of_mod6_one ha h6).2.1 (hcons _ (by omega) hgt)

/-- If `a % 6 = 5`, the next 6-rough after `a` is exactly `a + 2`. -/
theorem next_eq_add_two {a b : ℕ} (h : NextSixRough a b) (h6 : a % 6 = 5) :
    b = a + 2 := by
  obtain ⟨ha, hb, hab, hcons⟩ := h
  by_contra hne
  rcases show b < a + 2 ∨ a + 2 < b by omega with hlt | hgt
  · obtain ⟨_, _, _⟩ := hb; omega
  · exact absurd (nextSixRough_of_mod6_five ha h6).2.1 (hcons _ (by omega) hgt)

/-- The only gaps between consecutive 6-rough integers are 2 and 4. -/
theorem sixRough_gap (a b : ℕ) (h : NextSixRough a b) :
    b - a = 2 ∨ b - a = 4 := by
  rcases h.1.mod6 with h6 | h6
  · right; have := next_eq_add_four h h6; omega
  · left;  have := next_eq_add_two h h6; omega

end PT.Sieve
