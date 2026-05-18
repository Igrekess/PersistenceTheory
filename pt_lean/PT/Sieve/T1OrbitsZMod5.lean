/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Sieve.LegendreLogParity
import Mathlib.Tactic

/-!
# T1 — Orbits of the doubling map on `(ℤ/5ℤ)*`

The sieve transfer associated to the prime `p = 5` acts on the four non-zero
residues `(ℤ/5ℤ)^* = {1, 2, 3, 4}` by multiplication by the primitive root
`g = 2` (see `PT.Sieve.LegendreLogParity.two_primitive_root_mod5`). The
resulting permutation is a single 4-cycle:

  `1 ↦ 2 ↦ 4 ↦ 3 ↦ 1`.

This file formalises the **orbit structure** of the doubling map on
`(ℤ/5ℤ)^*`. It is the cyclic-4 analogue of
`PT.Sieve.T1AntidiagOrbits` (which covers the cyclic-2 case of `(ℤ/3ℤ)^*`).

* `sigma5 : Fin 4 → Fin 4` — the underlying 4-cycle, indexed so that
  `Fin 4 ≃ (1, 2, 4, 3) ∈ (ℤ/5ℤ)^*` in cycle order.
* `sigma5^[4] = id` — the permutation has order `4`.
* `sigma5 i ≠ i` for every `i` — no fixed point.
* `sigma5^[1], sigma5^[2], sigma5^[3]` are all distinct from `id`.
* Bridge to `ZMod 5`: the iterates of `2` traverse `{2, 4, 3, 1}`.
* **Dimensional comparison with `T1AntidiagOrbits`.** The cyclic-2 case
  has order `2`; the cyclic-4 case has order `4` (in particular,
  `sigma5^[2] ≠ id`).

## Reference

Theorem T1 (sieve transfer orbit structure) generalised from `p = 3`
to `p = 5`, Chapter 7 of the monograph, and
`PT.Sieve.LegendreLogParity.two_primitive_root_mod5`.
-/

namespace PT.Sieve

/-! ### The induced permutation `sigma5` on `Fin 4` -/

/-- The permutation of `Fin 4` induced by multiplication by `2` on
    `(ℤ/5ℤ)^*`, written in cycle order. With the indexing
    `0 ↔ 1, 1 ↔ 2, 2 ↔ 4, 3 ↔ 3`, the action `x ↦ 2x mod 5` becomes
    `0 ↦ 1 ↦ 2 ↦ 3 ↦ 0` on `Fin 4`. -/
def sigma5 : Fin 4 → Fin 4
  | 0 => 1
  | 1 => 2
  | 2 => 3
  | 3 => 0

@[simp] lemma sigma5_zero  : sigma5 0 = 1 := rfl
@[simp] lemma sigma5_one   : sigma5 1 = 2 := rfl
@[simp] lemma sigma5_two   : sigma5 2 = 3 := rfl
@[simp] lemma sigma5_three : sigma5 3 = 0 := rfl

/-! ### Iterates of `sigma5` -/

/-- Pointwise: `sigma5^[2] i` cycles `0 → 2, 1 → 3, 2 → 0, 3 → 1`. -/
@[simp] theorem sigma5_iterate_two (i : Fin 4) :
    sigma5^[2] i = sigma5 (sigma5 i) := by
  simp [Function.iterate_succ, Function.comp_apply]

/-- Pointwise: `sigma5^[3] i = sigma5 (sigma5 (sigma5 i))`. -/
theorem sigma5_iterate_three (i : Fin 4) :
    sigma5^[3] i = sigma5 (sigma5 (sigma5 i)) := by
  simp [Function.iterate_succ, Function.comp_apply]

/-- Pointwise: `sigma5^[4] i = i`. -/
@[simp] theorem sigma5_iterate_four_apply (i : Fin 4) :
    sigma5^[4] i = i := by
  fin_cases i <;> decide

/-- **Order 4.** `sigma5^[4] = id`. -/
theorem sigma5_iterate_four : sigma5^[4] = id := by
  funext i
  exact sigma5_iterate_four_apply i

/-! ### No fixed point -/

/-- `sigma5` has no fixed point: `sigma5 i ≠ i` for all `i`. -/
theorem sigma5_no_fixed_point (i : Fin 4) : sigma5 i ≠ i := by
  fin_cases i <;> decide

/-! ### Intermediate iterates are non-trivial -/

/-- The first iterate is non-trivial: `sigma5 ≠ id`. -/
theorem sigma5_iterate_one_ne_id : sigma5 ≠ id := by
  intro h
  have h0 : sigma5 0 = (id : Fin 4 → Fin 4) 0 := by rw [h]
  simp at h0

/-- The second iterate is non-trivial: `sigma5^[2] ≠ id`. -/
theorem sigma5_iterate_two_ne_id : sigma5^[2] ≠ id := by
  intro h
  have h0 : sigma5^[2] 0 = (id : Fin 4 → Fin 4) 0 := by rw [h]
  simp at h0

/-- The third iterate is non-trivial: `sigma5^[3] ≠ id`. -/
theorem sigma5_iterate_three_ne_id : sigma5^[3] ≠ id := by
  intro h
  have h0 : sigma5^[3] 0 = (id : Fin 4 → Fin 4) 0 := by rw [h]
  simp [sigma5_iterate_three] at h0

/-! ### Bijectivity and inverse structure -/

/-- `sigma5` is a bijection of `Fin 4`. -/
theorem sigma5_bijective : Function.Bijective sigma5 := by
  refine ⟨?_, ?_⟩
  · intro i j h
    fin_cases i <;> fin_cases j <;> simp_all
  · intro j
    fin_cases j
    · exact ⟨3, rfl⟩
    · exact ⟨0, rfl⟩
    · exact ⟨1, rfl⟩
    · exact ⟨2, rfl⟩

/-- `sigma5^[3]` is the inverse of `sigma5`. -/
theorem sigma5_inv_eq_iterate_three (i : Fin 4) :
    sigma5 (sigma5^[3] i) = i := by
  have h : sigma5^[4] i = i := sigma5_iterate_four_apply i
  have hstep : sigma5^[4] i = sigma5 (sigma5^[3] i) := by
    show sigma5^[3 + 1] i = sigma5 (sigma5^[3] i)
    rw [Function.iterate_succ_apply']
  rw [← hstep]; exact h

/-! ### Bridge with `ZMod 5`: powers of `2` -/

/-- **Bridge.** The iterates of multiplication by `2` on `(ℤ/5ℤ)^*` starting
    at `1` traverse `{2, 4, 3, 1}`, matching the four powers
    `2^1, 2^2, 2^3, 2^4`. -/
theorem pow_two_mod5_one : (2 : ZMod 5)^1 * 1 = 2 := by decide
theorem pow_two_mod5_two : (2 : ZMod 5)^2 * 1 = 4 := by decide
theorem pow_two_mod5_three : (2 : ZMod 5)^3 * 1 = 3 := by decide
theorem pow_two_mod5_four : (2 : ZMod 5)^4 * 1 = 1 := by decide

/-- Embedding `Fin 4 ↪ (ℤ/5ℤ)^*` matching the cycle order
    `0 ↦ 1, 1 ↦ 2, 2 ↦ 4, 3 ↦ 3`. -/
def toZMod5 : Fin 4 → ZMod 5
  | 0 => 1
  | 1 => 2
  | 2 => 4
  | 3 => 3

@[simp] lemma toZMod5_zero  : toZMod5 0 = 1 := rfl
@[simp] lemma toZMod5_one   : toZMod5 1 = 2 := rfl
@[simp] lemma toZMod5_two   : toZMod5 2 = 4 := rfl
@[simp] lemma toZMod5_three : toZMod5 3 = 3 := rfl

/-- **Conjugation identity.** The combinatorial cycle `sigma5` corresponds, via
    the embedding `toZMod5`, to multiplication by `2` on `(ℤ/5ℤ)^*`. -/
theorem toZMod5_sigma5 (i : Fin 4) :
    toZMod5 (sigma5 i) = 2 * toZMod5 i := by
  fin_cases i <;> decide

/-- The embedding never hits `0`. -/
theorem toZMod5_ne_zero (i : Fin 4) : toZMod5 i ≠ 0 := by
  fin_cases i <;> decide

/-! ### Dimensional comparison with the cyclic-2 case `T1AntidiagOrbits` -/

/-- **Dimensional comparison.** Whereas the cyclic-2 permutation `sigma` of
    `Fin 2` (cf. `PT.Sieve.T1AntidiagOrbits.sigma`) satisfies `sigma ∘ sigma
    = id`, the cyclic-4 permutation `sigma5` of `Fin 4` satisfies
    `sigma5^[4] = id` but `sigma5^[2] ≠ id`. The "dimension" of the orbit
    structure (the order of the underlying cycle) is therefore strictly
    larger here. -/
theorem sigma5_order_strictly_greater_than_two :
    sigma5^[2] ≠ id ∧ sigma5^[4] = id :=
  ⟨sigma5_iterate_two_ne_id, sigma5_iterate_four⟩

/-- **Minimal period is exactly 4.** None of the iterates `sigma5^[k]` for
    `k ∈ {1, 2, 3}` is the identity, while `sigma5^[4] = id`. -/
theorem sigma5_minimal_period_four :
    sigma5 ≠ id
    ∧ sigma5^[2] ≠ id
    ∧ sigma5^[3] ≠ id
    ∧ sigma5^[4] = id :=
  ⟨sigma5_iterate_one_ne_id,
   sigma5_iterate_two_ne_id,
   sigma5_iterate_three_ne_id,
   sigma5_iterate_four⟩

/-! ### Headline: orbit structure of the doubling map mod 5 -/

/-- **Headline.** Orbit structure of the doubling map on `(ℤ/5ℤ)^*`,
    encoded as the permutation `sigma5` on `Fin 4`:

    * `sigma5` is a single 4-cycle `(0 1 2 3)` with no fixed point.
    * The minimal period is exactly `4`: `sigma5^[k] ≠ id` for `k ∈ {1, 2, 3}`
      and `sigma5^[4] = id`.
    * Bridge to `ZMod 5`: `toZMod5 ∘ sigma5 = (2 *) ∘ toZMod5`, i.e. `sigma5`
      conjugates to the multiplication-by-`2` map on `(ℤ/5ℤ)^*`.
    * **Dimensional comparison.** Strictly stronger than the cyclic-2
      `T1AntidiagOrbits.sigma`: order `4` vs order `2`. -/
theorem sigma5_orbit_structure :
    -- Cycle properties:
    sigma5^[4] = id
    ∧ (∀ i, sigma5 i ≠ i)
    ∧ sigma5 0 = 1 ∧ sigma5 1 = 2 ∧ sigma5 2 = 3 ∧ sigma5 3 = 0
    -- Minimal period 4:
    ∧ sigma5 ≠ id ∧ sigma5^[2] ≠ id ∧ sigma5^[3] ≠ id
    -- Bridge to ZMod 5:
    ∧ (∀ i, toZMod5 (sigma5 i) = 2 * toZMod5 i)
    -- Bijectivity:
    ∧ Function.Bijective sigma5 :=
  ⟨sigma5_iterate_four,
   sigma5_no_fixed_point,
   rfl, rfl, rfl, rfl,
   sigma5_iterate_one_ne_id,
   sigma5_iterate_two_ne_id,
   sigma5_iterate_three_ne_id,
   toZMod5_sigma5,
   sigma5_bijective⟩

end PT.Sieve
