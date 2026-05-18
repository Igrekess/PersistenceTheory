/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Sieve.T6aFieldStructure

/-!
# T6b — Axioms C1–C4 force `E_p = {0}` (Eratosthenes from ring congruence)

**Statement (paper-level, Ch02 §"Théorème T6b").**
Let `E_p ⊆ ℤ/pℤ` be an elimination procedure satisfying the four axioms

* **C1 (ring congruence):** `E_p` is an ideal of `ℤ/pℤ`;
* **C2 (absorption and non-triviality):** `(p : ℤ/pℤ) = [0] ∈ E_p`
  (so `E_p ≠ ∅`) and `E_p ≠ ℤ/pℤ` (the procedure does not delete
  everything);
* **C3 (irreducibility):** `E_p` is not further decomposable;
* **C4 (static completeness):** the surviving residues are closed under
  multiplication.

Then `E_p = {[0]}` is the **unique** procedure. (This is the M2 article's
formal kernel; the proof here keeps the structural backbone.)

We formalise the **C1 ∧ C2 ⇒ E_p = {[0]}** kernel, where the key step is
the field dichotomy already proved in `T6aFieldStructure`. C3 and C4 are
automatic on this single-element ideal (and serve to exclude composite
moduli or non-multiplicative survivor sets — those structural claims are
recorded as remarks).

The `4/4 independence` (each of C1, C2, C3, C4 is necessary) is documented
in the paper but not formalised here (each counter-model would be its own
small theorem; the kernel uniqueness from C1+C2 is the main result).

## Reference

Monograph chapter 2, §"Théorème T6b", `\label{thm:T6_universal}`.
M2 article (PT_ARTICLES/PT_MATHEMATICS/M2), section "T6b".
-/

namespace PT.Sieve

variable (p : ℕ) [Fact p.Prime]

/-! ### Axioms C1 and C2 as structural predicates on ideals -/

/-- **C1 (ring congruence).** Vacuous: in our formalisation we model `E_p`
    *as* an ideal of `ℤ/pℤ`, so C1 is built into the type. -/
abbrev C1 (_E : Ideal (ZMod p)) : Prop := True

/-- **C2 (absorption and non-triviality).** The zero residue is in `E`,
    and `E` is *not* the entire ring `ℤ/pℤ`. -/
def C2 (E : Ideal (ZMod p)) : Prop :=
  ((0 : ZMod p) ∈ E) ∧ (E ≠ ⊤)

/-! ### T6b kernel: C1 ∧ C2 ⇒ E_p = {[0]} -/

/-- **T6b (kernel form).** Any ideal `E` of `ℤ/pℤ` that contains `[0]` and
    is *not* the full ring is equal to `⊥` (i.e. `{[0]}`).

    Proof: by the field dichotomy `T6a_ideal_dichotomy`,
    `E ∈ {⊥, ⊤}`. C2 (`E ≠ ⊤`) selects `E = ⊥`. The absorption part of
    C2 (`[0] ∈ E`) is automatic from `0 ∈ ⊥`, hence consistent. -/
theorem T6b_unique (E : Ideal (ZMod p)) (_h1 : C1 p E) (h2 : C2 p E) :
    E = ⊥ :=
  T6a_unique_proper_ideal p E h2.2

/-- **T6b (final form).** The unique elimination procedure on `ℤ/pℤ`
    satisfying C1 ∧ C2 is the singleton `{[0]}`, characterised as
    `x ∈ E ↔ x = 0`. -/
theorem T6b_iff_zero (E : Ideal (ZMod p)) (_h1 : C1 p E) (h2 : C2 p E)
    (x : ZMod p) : x ∈ E ↔ x = 0 := by
  rw [T6b_unique p E _h1 h2]
  exact T6a_bot_eq_zero p x

/-! ### Independence of the four axioms (informal record)

The paper records *four* independent counter-models showing that each of
C1, C2, C3, C4 is necessary. Each counter-model is a small structural
theorem about a specific deviation from the canonical sieve. We do *not*
formalise the four counter-models here — they each require an
independent setup, and their formalisation is left as future work.

* `¬C1`: a non-symmetric removal set (distinguishing the residue `1` from
  its negative) violates bilateral symmetry, hence does not yield an ideal.
* `¬C2`: the "complementary" set `E_p = 1 + pℤ` (remove residue `1` instead
  of `0`) is a non-trivial coset, *not* an ideal containing zero.
* `¬C3`: the union `E = pℤ ∪ {4}` decomposes into a multiplicative part
  and an additive correction, violating irreducibility.
* `¬C4`: a survivor set not closed under multiplication.
-/

end PT.Sieve
