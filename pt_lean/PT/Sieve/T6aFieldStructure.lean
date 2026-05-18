/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import Mathlib.Algebra.Field.ZMod
import Mathlib.RingTheory.Ideal.Lattice
import Mathlib.Tactic

/-!
# T6a — Field-structure uniqueness: `R_p = {0}` in `ℤ/pℤ`

**Statement (paper-level, Ch02 §"Théorème T6a").**
For every prime `p`, the only proper ideal of `ℤ/pℤ` is `{0}`. Consequently,
the Eratosthenes elimination set `R_p` in `ℤ/pℤ` is `{0} = {[0]}`.

**Proof.** `ZMod p` (`p` prime) is a field; in any (division) ring all
ideals are either `⊥` or `⊤`; the *proper* ideal is `⊥`, i.e. `{0}`.

This is one of the three independent proofs of the T6-family uniqueness
(`T6a` field, `T6b` 4-axiom, `T6c` Shore-Johnson-Čencov). Here we discharge
`T6a` via Mathlib's `Ideal.eq_bot_or_top` for division (semi)rings.

## Reference

Monograph chapter 2, §"Théorème T6a", `\label{thm:T6_restricted}`.
M2 article (PT_ARTICLES/PT_MATHEMATICS/M2), section "T6a".
-/

namespace PT.Sieve

open ZMod

variable (p : ℕ) [Fact p.Prime]

/-- **T6a — Field-structure dichotomy on `ZMod p`.**
    Every ideal of `ℤ/pℤ` is either `⊥` or `⊤`. -/
theorem T6a_ideal_dichotomy (I : Ideal (ZMod p)) : I = ⊥ ∨ I = ⊤ :=
  Ideal.eq_bot_or_top I

/-- **T6a — Field structure (uniqueness of the proper ideal).**
    The only proper ideal of `ℤ/pℤ` is `⊥`. -/
theorem T6a_unique_proper_ideal (I : Ideal (ZMod p)) (hne : I ≠ ⊤) : I = ⊥ :=
  (T6a_ideal_dichotomy p I).resolve_right hne

/-- **T6a — Bottom ideal = `{0}` characterisation.**
    The bottom ideal of `ℤ/pℤ` consists exactly of the zero element. -/
theorem T6a_bot_eq_zero (x : ZMod p) : x ∈ (⊥ : Ideal (ZMod p)) ↔ x = 0 :=
  Ideal.mem_bot

/-- **T6a (final form).** For any ideal `I ≠ ⊤` of `ℤ/pℤ`,
    `I = {0}` as a subset of `ZMod p`. -/
theorem T6a_field_structure (I : Ideal (ZMod p)) (hne : I ≠ ⊤) :
    ∀ x : ZMod p, x ∈ I ↔ x = 0 := by
  have h := T6a_unique_proper_ideal p I hne
  intro x
  rw [h]
  exact T6a_bot_eq_zero p x

end PT.Sieve
