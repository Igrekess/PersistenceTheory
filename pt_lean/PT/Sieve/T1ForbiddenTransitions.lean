/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Sieve.SixRough

/-!
# T1 — Forbidden Self-Transitions in the 6-Rough Sequence

Consecutive 6-rough integers (positive integers coprime to 6) always switch
their residue modulo 3. The sieve-level transfer matrix on `{1, 2} mod 3` is
therefore the anti-diagonal permutation matrix

$$T_3 = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}.$$

Its eigenvalues are `{+1, −1}` and its unique stationary distribution is
`π = (1/2, 1/2)`, yielding the symmetry parameter **s = 1/2**.

The matrix-level statements are formalised in `PT.Sieve.T3Antidiagonal`,
`PT.Stochastic.SHalf` and `PT.Conservation.T2Alpha`. This file contains the
purely number-theoretic statement at the level of `ℕ`.

## Reference

Theorem T1 of: Y. Senez, *Forbidden transitions and informational structure of
prime residues under modular projection* (PT_ARTICLES/PT_MATHEMATICS/M1).
-/

namespace PT.Sieve.T1

open PT.Sieve

/-! ### Theorem T1 -/

/-- **Theorem T1 (Forbidden self-transitions).**
    Consecutive 6-rough integers always have different residues modulo 3.
    The diagonal of the sieve-level transfer matrix `T₃` is zero. -/
theorem T1_forbidden_self_transition (a b : ℕ) (h : NextSixRough a b) :
    a % 3 ≠ b % 3 := by
  rcases h.1.mod6 with h6 | h6
  · -- a ≡ 1 (mod 6)  ⟹  b = a + 4  ⟹  residues 1 → 2
    have := next_eq_add_four h h6; omega
  · -- a ≡ 5 (mod 6)  ⟹  b = a + 2  ⟹  residues 2 → 1
    have := next_eq_add_two h h6; omega

/-! ### Corollaries -/

/-- The mod-3 residue alternates deterministically: `1 → 2` or `2 → 1`. -/
theorem T1_residue_switch (a b : ℕ) (h : NextSixRough a b) :
    (a % 3 = 1 ∧ b % 3 = 2) ∨ (a % 3 = 2 ∧ b % 3 = 1) := by
  rcases h.1.mod6 with h6 | h6
  · left;  exact ⟨by omega, by have := next_eq_add_four h h6; omega⟩
  · right; exact ⟨by omega, by have := next_eq_add_two h h6; omega⟩

/-- The transfer matrix is the involution `σ(1) = 2, σ(2) = 1`,
    i.e. `T₃ = antidiag(1, 1)`. -/
theorem T1_antidiag (a b : ℕ) (h : NextSixRough a b) :
    (a % 3 = 1 → b % 3 = 2) ∧ (a % 3 = 2 → b % 3 = 1) := by
  rcases T1_residue_switch a b h with ⟨_, h2⟩ | ⟨_, h2⟩ <;>
    exact ⟨fun _ => by omega, fun _ => by omega⟩

/-! ### Concrete witnesses (verified by the kernel) -/

/-- 5 → 7: gap 2, residues 2 → 1. -/
example : NextSixRough 5 7 := by
  refine ⟨⟨by omega, by omega, by omega⟩,
          ⟨by omega, by omega, by omega⟩, by omega, fun c _ _ => ?_⟩
  rintro ⟨_, _, _⟩; omega

/-- 7 → 11: gap 4, residues 1 → 2. -/
example : NextSixRough 7 11 := by
  refine ⟨⟨by omega, by omega, by omega⟩,
          ⟨by omega, by omega, by omega⟩, by omega, fun c _ _ => ?_⟩
  rintro ⟨_, _, _⟩; omega

/-- 11 → 13: gap 2, residues 2 → 1. -/
example : NextSixRough 11 13 := by
  refine ⟨⟨by omega, by omega, by omega⟩,
          ⟨by omega, by omega, by omega⟩, by omega, fun c _ _ => ?_⟩
  rintro ⟨_, _, _⟩; omega

/-- 13 → 17: gap 4, residues 1 → 2. -/
example : NextSixRough 13 17 := by
  refine ⟨⟨by omega, by omega, by omega⟩,
          ⟨by omega, by omega, by omega⟩, by omega, fun c _ _ => ?_⟩
  rintro ⟨_, _, _⟩; omega

-- Residue switches verified numerically:
example : 5  % 3 ≠ 7  % 3 := by decide
example : 7  % 3 ≠ 11 % 3 := by decide
example : 11 % 3 ≠ 13 % 3 := by decide
example : 13 % 3 ≠ 17 % 3 := by decide
example : 17 % 3 ≠ 19 % 3 := by decide

end PT.Sieve.T1
