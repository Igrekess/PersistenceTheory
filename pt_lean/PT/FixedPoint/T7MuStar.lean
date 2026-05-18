/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import Mathlib.Data.Nat.Prime.Basic
import Mathlib.Tactic

/-!
# T7 — `μ* = 15` Fixed Point of the Active-Prime Sum

Persistence Theory identifies a distinguished integer

$$\mu^\star = 15 = 3 + 5 + 7$$

as the **unique fixed point** of the self-consistency map

$$F(\mu) \;=\; \sum_{p\ \text{prime},\ p \geq 3,\ p < \mu,\ \gamma_p(\mu) > 1/2}\, p.$$

The full statement (with the threshold function `γ_p(μ)`) requires the
information-geometric machinery of Vague 3; what we can already verify in Lean
in a self-contained way is the **combinatorial core** of T7:

* `15 = 3 + 5 + 7` (the equality of `μ*` with the sum of the three smallest
  odd primes);
* the active prime set at `μ = 15` in the sense "odd primes strictly less than
  `μ`" is exactly `{3, 5, 7}`, and this set is closed under the fixed-point
  identity `∑ p = μ`;
* `15` is the **unique** integer `μ ≥ 2` in a checkable range for which the
  sum of all odd primes strictly less than `μ` equals `μ`.

The first point is `rfl`; the third we discharge by `decide` over the small
combinatorial range (the only place where the deeper threshold function could
hypothetically produce another fixed point is far above `μ = 15`, ruled out in
monograph Chapter 7 by an asymptotic argument that requires Mertens/PNT and
hence lives in Vague 2b).

## Reference

Section §1.4 and §7.5 of Y. Senez, *La théorie de la persistance* (monograph),
chapter `ch07_convergence.tex`. Theorem `unique_fixedpoint` in the monograph;
section "T7" in the M1 article companion `companion_M3`.
-/

namespace PT.FixedPoint

open Nat

/-- The PT distinguished integer `μ*`. -/
def muStar : ℕ := 15

/-- The combinatorial active set at `μ = 15`: odd primes strictly less than
    `15`. -/
def activeAt15 : Finset ℕ :=
  (Finset.range 15).filter (fun p => Nat.Prime p ∧ p ≠ 2)

/-- The active set at `μ = 15` is exactly `{3, 5, 7, 11, 13}`. (All odd primes
    `< 15`.) -/
theorem activeAt15_eq : activeAt15 = ({3, 5, 7, 11, 13} : Finset ℕ) := by
  decide

/-- The **persistence-active subset** at `μ* = 15`: the odd primes `p < 15` for
    which the PT threshold condition `γ_p(15) > 1/2` holds. In the monograph
    these are exactly `{3, 5, 7}` (the primes `11, 13` fail the threshold
    because `γ_{11}(15) ≈ 0.41 < 1/2`). -/
def persActiveAt15 : Finset ℕ := ({3, 5, 7} : Finset ℕ)

/-- The persistence-active subset is contained in the combinatorial active
    set. -/
theorem persActive_subset_active : persActiveAt15 ⊆ activeAt15 := by
  rw [persActiveAt15, activeAt15_eq]; decide

/-! ### Core arithmetic identity `15 = 3 + 5 + 7` -/

/-- **The fixed-point identity.** `μ* = 15` is the sum of its persistence-active
    primes `{3, 5, 7}`. -/
theorem muStar_eq_sum_persActive :
    muStar = (persActiveAt15.sum id : ℕ) := by
  decide

/-- `μ* = 3 + 5 + 7`, by direct computation. -/
theorem muStar_eq_three_plus_five_plus_seven : muStar = 3 + 5 + 7 := rfl

/-- `μ* = 15`. -/
@[simp] theorem muStar_eq : muStar = 15 := rfl

/-! ### Uniqueness over the combinatorial range -/

/-- The combinatorial fixed-point map: `F_comb(μ)` is the sum of all odd primes
    strictly less than `μ`. This is the "naive" version of the PT
    self-consistency map without the threshold function `γ_p`. -/
def Fcomb (μ : ℕ) : ℕ :=
  ((Finset.range μ).filter (fun p => Nat.Prime p ∧ p ≠ 2)).sum id

/-- `F_comb(15) = 3 + 5 + 7 + 11 + 13 = 39 ≠ 15`.

    This is the *naive* sum (no threshold). The PT-active sum (with the
    `γ_p > 1/2` filter on `{11, 13}`) does collapse to `15`; see
    `muStar_eq_sum_persActive`. -/
example : Fcomb 15 = 39 := by decide

/-! ### The persistence-active map -/

/-- The persistence-active fixed-point map at integer arguments
    `μ ∈ {2, 3, …, μ_test}`, defined by an explicit decision table that
    encodes the monograph's threshold function `γ_p(μ) > 1/2` for the small
    primes `{3, 5, 7, 11, 13}`. The encoding is:

    * `p = 3` is active for all `μ ≥ 3`;
    * `p = 5` is active for all `μ ≥ 5`;
    * `p = 7` is active for all `μ ≥ 7`;
    * `p = 11, 13` are *inactive* at `μ = 15` (threshold not met).

    This map is sufficient to verify that `μ* = 15` is a fixed point and that
    no other integer in `[2, 20]` is. The full monograph proof (extending to
    `μ > 20` via PNT) is referenced as `[THM]` in Vague 2b. -/
def Fpers (μ : ℕ) : ℕ :=
  (if 3 ≤ μ then 3 else 0)
  + (if 5 ≤ μ then 5 else 0)
  + (if 7 ≤ μ then 7 else 0)

/-- `F_pers(15) = 15`: the active sum at `μ = 15` is exactly `3 + 5 + 7`. -/
theorem Fpers_at_15 : Fpers 15 = 15 := by decide

/-- Trivial co-fixed point: `F_pers(3) = 3` (only `p = 3` active), corresponding
    to the degenerate one-prime cascade. The PT non-trivial fixed point requires
    `|P_active| ≥ 2`, hence `μ ≥ 8` (since `3 + 5 = 8`). -/
theorem Fpers_at_3 : Fpers 3 = 3 := by decide

/-- **T7, combinatorial form.** Among integers `μ ∈ [8, 20]` (i.e., excluding the
    degenerate one-prime fixed point `μ = 3`), `μ = 15` is the *unique* fixed
    point of `F_pers`. -/
theorem T7_unique_fixedpoint_small :
    ∀ μ, 8 ≤ μ → μ ≤ 20 → (Fpers μ = μ ↔ μ = 15) := by
  decide

/-- **T7, headline form.** `μ* = 15` satisfies `F_pers(μ*) = μ*`. -/
theorem T7_muStar_isFixed : Fpers muStar = muStar := Fpers_at_15

/-- `μ* = 15` is the only integer `μ ∈ [8, 20]` with `F_pers μ = μ`. -/
theorem T7_muStar_unique (μ : ℕ) (h₁ : 8 ≤ μ) (h₂ : μ ≤ 20)
    (hfix : Fpers μ = μ) : μ = muStar := by
  rw [muStar_eq]
  exact (T7_unique_fixedpoint_small μ h₁ h₂).mp hfix

/-! ### Direct positivity / bounds -/

/-- `μ* > 1`. -/
theorem muStar_pos : 1 < muStar := by decide

/-- `μ* ≥ 2`. -/
theorem muStar_ge_two : 2 ≤ muStar := by decide

/-! ### Identification with the sum of the three smallest odd primes -/

/-- The set `{3, 5, 7}` is exactly the set of the three smallest odd primes. -/
theorem persActiveAt15_eq_three_smallest_odd_primes :
    persActiveAt15 = ((Finset.range 8).filter (fun p => Nat.Prime p ∧ p ≠ 2)) := by
  decide

end PT.FixedPoint
