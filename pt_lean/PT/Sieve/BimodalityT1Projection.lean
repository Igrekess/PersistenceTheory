/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Sieve.Bimodality
import PT.Sieve.LegendreLogParity
import PT.Sieve.SixRough
import Mathlib.Tactic

/-!
# Bimodality as T1 projection (Appendix N, #42)

**Statement (paper-level, App N).**
The bimodality dichotomy `δ̄(r) ∈ {7, 11}` on admissible residues `r` modulo
`30` is the projection along the `mod 5`-fibre of the underlying T1 bifurcation
factor `2 = p₁` (the first prime forbidden as a "transition" in the
antidiagonal `T₃`). Concretely:

* The amplitude `Δ = 4` of `δ̄` factors as `Δ = 2 · p₁ = 2 · 2`.
* The fibre splits the 8 admissible residues `R ⊂ (ℤ/30ℤ)*` into the two
  classes of the Legendre symbol mod `5`: 4 residues with `(r/5) = +1`
  (giving `δ̄ = 7`) and 4 with `(r/5) = -1` (giving `δ̄ = 11`).

This file formalises the **bridge** between the bimodality combinatorial
core (`PT/Sieve/Bimodality.lean`) and the Legendre = log-parity identity
(`PT/Sieve/LegendreLogParity.lean`), and quantifies the amplitude in terms
of the T1 bifurcation prime `p₁ = 2`.

## Reference

Appendix N of the monograph, `\label{cor:bimodality_T1}`. Item #42 of
`AUDIT_FORMALISABLE.md`.
-/

namespace PT.Sieve

/-! ### The T1 bifurcation prime -/

/-- The T1 bifurcation prime `p₁ = 2`. This is the smallest prime, and the
    unique prime "transition" outlawed by the antidiagonal `T₃` of T1
    (cf. `PT.Sieve.T1ForbiddenTransitions`). -/
def p1 : ℕ := 2

/-- `p₁ = 2`. -/
@[simp] theorem p1_eq_two : p1 = 2 := rfl

/-- `p₁` is the smallest prime. -/
theorem p1_prime : Nat.Prime p1 := by decide

/-! ### Amplitude factorisation -/

/-- The bimodality amplitude `Δ̄` in units of `δ̄` (after dividing the
    `deltaBarTimes6`-difference by `6`): `(66 - 42) / 6 = 4`. -/
def amplitude : ℕ := (66 - 42) / 6

/-- `amplitude = 4`. -/
@[simp] theorem amplitude_eq_four : amplitude = 4 := by decide

/-- **Corollary (App N #42) — Bimodality amplitude = `2 · p₁`.**
    The bimodality amplitude `Δ̄ = 4` factors as `2 · p₁ = 2 · 2`, where
    `p₁ = 2` is the T1 bifurcation prime. -/
theorem bimodality_amplitude_eq_two_p1 : amplitude = 2 * p1 := by
  decide

/-! ### QR/NR partition of admissible residues -/

/-- The four admissible residues with `δ̄ = 7` (i.e. `deltaBarTimes6 = 42`). -/
def lowResidues : Finset ℕ := {1, 11, 19, 29}

/-- The four admissible residues with `δ̄ = 11` (i.e. `deltaBarTimes6 = 66`). -/
def highResidues : Finset ℕ := {7, 13, 17, 23}

/-- `low ⊔ high = admissibleResidues`, and they are disjoint. -/
theorem low_union_high_eq_admissible :
    lowResidues ∪ highResidues = admissibleResidues := by
  decide

theorem low_disjoint_high : Disjoint lowResidues highResidues := by
  decide

/-! ### Legendre projection: the dichotomy is `(r/5)` -/

/-- Cast a `ℕ`-residue (modulo 30) to `ZMod 5` by reducing modulo `5`. -/
def to5 (r : ℕ) : ZMod 5 := (r : ZMod 5)

/-- The four low-residues reduce modulo `5` to quadratic residues
    `{1, 4} ⊂ (ℤ/5ℤ)*`. -/
theorem low_residues_are_QR :
    ∀ r ∈ lowResidues, legendre5 (to5 r) = 1 := by
  decide

/-- The four high-residues reduce modulo `5` to non-residues
    `{2, 3} ⊂ (ℤ/5ℤ)*`. -/
theorem high_residues_are_NR :
    ∀ r ∈ highResidues, legendre5 (to5 r) = -1 := by
  decide

/-! ### Headline: bimodality as T1 projection -/

/-- **Theorem (App N #42) — Bimodality dichotomy = Legendre projection mod 5.**

    The 8 admissible residues `r ∈ R ⊂ (ℤ/30ℤ)*` partition into two classes
    of size `4` each, witnessed by `δ̄(r)`:

    * `δ̄(r) = 7` (i.e. `deltaBarTimes6 r = 42`) iff `r` is a QR mod 5
      (Legendre `(r/5) = +1`);
    * `δ̄(r) = 11` (i.e. `deltaBarTimes6 r = 66`) iff `r` is a NR mod 5
      (Legendre `(r/5) = -1`).

    Combined with `bimodality_amplitude_eq_two_p1`, the bimodality theorem
    is the projection of the T1 bifurcation factor `2 = p₁` along the
    `mod 5`-fibre. -/
theorem bimodality_eq_legendre_projection :
    (∀ r ∈ lowResidues, deltaBarTimes6 (r : ZMod 30) = 42
                        ∧ legendre5 (to5 r) = 1) ∧
    (∀ r ∈ highResidues, deltaBarTimes6 (r : ZMod 30) = 66
                         ∧ legendre5 (to5 r) = -1) := by
  refine ⟨?_, ?_⟩
  · intro r hr; refine ⟨?_, ?_⟩
    · revert r hr; decide
    · revert r hr; decide
  · intro r hr; refine ⟨?_, ?_⟩
    · revert r hr; decide
    · revert r hr; decide

/-! ### Cardinality balance -/

/-- The two `δ̄`-classes have equal cardinality `|low| = |high| = 4`. -/
theorem low_card : lowResidues.card = 4 := by decide
theorem high_card : highResidues.card = 4 := by decide

/-- The cardinality balance `|low| = |high|` is the projection of the
    `(ℤ/5ℤ)*`-equidistribution `|QR| = |NR| = 2` (factor of `2` accounts for
    the `(ℤ/3ℤ)*`-fibre, which is trivial since the admissible class is fixed
    modulo `3` once `mod 30`). -/
theorem bimodality_class_balance :
    lowResidues.card = highResidues.card := by decide

end PT.Sieve
