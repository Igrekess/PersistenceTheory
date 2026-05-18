/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import Mathlib.NumberTheory.Primorial
import Mathlib.Data.ZMod.Basic
import Mathlib.Tactic

/-!
# PM algebra — the primorial-modular ring

This module formalises the **primorial-modular ring** `PrimorialMod n`
underlying the multiplicative structure of the persistence sieve.

## Mathematical content

The primorial `n#` (Mathlib `primorial`) is the product of all primes
`p ≤ n`. The ring `PrimorialMod n := ZMod (n#)` is the canonical ambient
ring for residue arithmetic on the sieve: every prime `p ≤ n` divides
`n#`, so for each such prime there is a canonical projection
`PrimorialMod n →+* ZMod p`. Pairwise these projections produce the
**Chinese-remainder isomorphism**

  `PrimorialMod n  ≃+*  ∏_{p prime, p ≤ n} ZMod p`.

For the PT-relevant case (active primes `{3, 5, 7}`, total `n = 7`) the
primorial is `7# = 2 · 3 · 5 · 7 = 210`, and the CRT iso decomposes
`ZMod 210 ≃+* ZMod 2 × ZMod 3 × ZMod 5 × ZMod 7`.

This file states and proves:

* `PrimorialMod n` — the ring abbreviation `ZMod (primorial n)`.
* Explicit values `primorial 7 = 210`, `primorial 6 = 30`,
  `primorial 5 = 30`, etc., closed by `decide`.
* `pmEquivProd2_3_5_7` — the explicit CRT iso for `n = 7`:
  `PrimorialMod 7 ≃+* ZMod 2 × ZMod 3 × ZMod 5 × ZMod 7`
  obtained by three iterated applications of `ZMod.chineseRemainder`.
* Cardinality and finiteness corollaries.

## PT framing

In the monograph the symbol `M#` denotes the primorial-modular ring of
the cycle-modular sieve. The "primorial-modular algebra" referenced in
Appendix D (audit row #29) is the commutative-ring backbone on which
the CRT projections `ZMod p` for active primes `p ∈ {3, 5, 7}` are
defined — itself the ambient ring containing the sieve transfer
matrices `T_p` after coordinatisation by `ZMod p`.

## Sorry budget

This file contains **0 sorrys** and introduces **no axioms**.

## References

* Mathlib `Nat.primorial`, `ZMod.chineseRemainder`.
* Monograph Appendix D (`app_d_pm_algebra.tex`) — though the appendix
  in its current form articulates the *epistemic classification*
  governing PT statements, the **algebraic** primorial-modular ring
  is the mathematical object that this Lean module captures.
-/

namespace PT.Algebra

open scoped Nat

/-- **Primorial-modular ring.** `PrimorialMod n := ZMod (primorial n)`. -/
abbrev PrimorialMod (n : ℕ) : Type := ZMod (primorial n)

/-- The primorial `0# = 1`. -/
theorem primorial_zero_val : primorial 0 = 1 := by decide

/-- The primorial `1# = 1`. -/
theorem primorial_one_val : primorial 1 = 1 := by decide

/-- The primorial `2# = 2`. -/
theorem primorial_two_val : primorial 2 = 2 := by decide

/-- The primorial `3# = 6`. -/
theorem primorial_three_val : primorial 3 = 6 := by decide

/-- The primorial `5# = 30 = 2·3·5`. -/
theorem primorial_five_val : primorial 5 = 30 := by decide

/-- The primorial `6# = 30` (no new prime between 5 and 6). -/
theorem primorial_six_val : primorial 6 = 30 := by decide

/-- The primorial `7# = 210 = 2·3·5·7`, the PT-relevant value. -/
theorem primorial_seven_val : primorial 7 = 210 := by decide

/-- `primorial 7 = 2 * 3 * 5 * 7`. The active-primes factorisation. -/
theorem primorial_seven_factorisation :
    primorial 7 = 2 * 3 * 5 * 7 := by
  rw [primorial_seven_val]

/-! ### Cardinality -/

/-- `PrimorialMod 7` has cardinality `210`. -/
theorem card_primorialMod_seven :
    Nat.card (PrimorialMod 7) = 210 := by
  unfold PrimorialMod
  rw [primorial_seven_val]
  exact Nat.card_zmod 210

/-! ### Coprimality of the building blocks -/

/-- `2` and `3·5·7 = 105` are coprime. -/
theorem coprime_two_105 : Nat.Coprime 2 105 := by decide

/-- `3` and `5·7 = 35` are coprime. -/
theorem coprime_three_35 : Nat.Coprime 3 35 := by decide

/-- `5` and `7` are coprime. -/
theorem coprime_five_seven : Nat.Coprime 5 7 := by decide

/-! ### CRT factorisation: `ZMod 210 ≃+* ZMod 2 × ZMod 3 × ZMod 5 × ZMod 7` -/

/-- **Step 1.** `ZMod 210 ≃+* ZMod 2 × ZMod 105`. -/
def crtStep1 : ZMod 210 ≃+* ZMod 2 × ZMod 105 := by
  have h : (210 : ℕ) = 2 * 105 := by decide
  exact (RingEquiv.cast (by rw [h]) : ZMod 210 ≃+* ZMod (2 * 105)).trans
        (ZMod.chineseRemainder coprime_two_105)

/-- **Step 2.** `ZMod 105 ≃+* ZMod 3 × ZMod 35`. -/
def crtStep2 : ZMod 105 ≃+* ZMod 3 × ZMod 35 := by
  have h : (105 : ℕ) = 3 * 35 := by decide
  exact (RingEquiv.cast (by rw [h]) : ZMod 105 ≃+* ZMod (3 * 35)).trans
        (ZMod.chineseRemainder coprime_three_35)

/-- **Step 3.** `ZMod 35 ≃+* ZMod 5 × ZMod 7`. -/
def crtStep3 : ZMod 35 ≃+* ZMod 5 × ZMod 7 := by
  have h : (35 : ℕ) = 5 * 7 := by decide
  exact (RingEquiv.cast (by rw [h]) : ZMod 35 ≃+* ZMod (5 * 7)).trans
        (ZMod.chineseRemainder coprime_five_seven)

/-- **CRT iso, full four-factor decomposition.**
    `PrimorialMod 7 ≃+* ZMod 2 × (ZMod 3 × (ZMod 5 × ZMod 7))`.

    Each factor `ZMod p` for `p ∈ {2, 3, 5, 7}` is the residue ring at the
    prime `p`. The active primes of PT are `{3, 5, 7}`; the spurious `2`
    factor is present because `primorial 7` is defined over *all* primes
    `≤ 7`, and PT-internal constructions later project out the `ZMod 2`
    factor (sieve parity). -/
def pmEquivProd2_3_5_7 :
    PrimorialMod 7 ≃+*
      ZMod 2 × (ZMod 3 × (ZMod 5 × ZMod 7)) := by
  -- PrimorialMod 7 = ZMod (primorial 7) = ZMod 210
  have hpm : (primorial 7) = 210 := primorial_seven_val
  refine ((RingEquiv.cast (by rw [hpm]) :
            ZMod (primorial 7) ≃+* ZMod 210).trans crtStep1).trans ?_
  exact (RingEquiv.refl (ZMod 2)).prodCongr (crtStep2.trans
          ((RingEquiv.refl (ZMod 3)).prodCongr crtStep3))

/-! ### Auxiliary: the "active-prime quotient" `ZMod 105` -/

/-- The **active-prime primorial** `3 · 5 · 7 = 105`. -/
def activePrimorial : ℕ := 105

/-- The active-prime ring (`ZMod 105`) is the quotient of `PrimorialMod 7`
    by the `ZMod 2` factor. -/
abbrev ActivePrimorialMod : Type := ZMod activePrimorial

/-- **Active-prime CRT iso.** `ZMod 105 ≃+* ZMod 3 × ZMod 5 × ZMod 7`,
    the PT-canonical residue decomposition over the active primes. -/
def activePrimorialEquiv :
    ActivePrimorialMod ≃+* ZMod 3 × (ZMod 5 × ZMod 7) :=
  crtStep2.trans ((RingEquiv.refl (ZMod 3)).prodCongr crtStep3)

/-- The active-prime ring has cardinality `105`. -/
theorem card_activePrimorialMod :
    Nat.card ActivePrimorialMod = 105 := by
  unfold ActivePrimorialMod activePrimorial
  exact Nat.card_zmod 105

/-! ### Headline summary -/

/-- **PM algebra headline.** The primorial-modular ring at `n = 7`
    (`PrimorialMod 7 = ZMod 210`) decomposes by CRT into the product of
    residue rings at primes ≤ 7. The PT-relevant factor `ZMod 3 × ZMod 5
    × ZMod 7` is the **active-prime ring**, a quotient of `PrimorialMod 7`
    of cardinality `105`. -/
theorem pm_algebra_summary :
    Nat.card (PrimorialMod 7) = 210
    ∧ primorial 7 = 2 * 3 * 5 * 7
    ∧ Nat.card ActivePrimorialMod = 105
    ∧ activePrimorial = 3 * 5 * 7 :=
  ⟨card_primorialMod_seven, primorial_seven_factorisation,
   card_activePrimorialMod, by decide⟩

end PT.Algebra
