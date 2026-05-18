/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Sieve.T1OrbitsZMod5
import PT.Sieve.T1AntidiagOrbits
import Mathlib.Tactic

/-!
# T1 — Orbits of the multiplication-by-3 map on `(ℤ/7ℤ)*`

The sieve transfer associated to the prime `p = 7` acts on the six non-zero
residues `(ℤ/7ℤ)^* = {1, 2, 3, 4, 5, 6}` by multiplication by the primitive
root `g = 3`. Indeed,

  `3^1 = 3, 3^2 = 2, 3^3 = 6, 3^4 = 4, 3^5 = 5, 3^6 = 1`  (mod 7),

so the resulting permutation is a single 6-cycle:

  `1 ↦ 3 ↦ 2 ↦ 6 ↦ 4 ↦ 5 ↦ 1`.

This file formalises the **orbit structure** of the multiplication-by-3
map on `(ℤ/7ℤ)^*`. It is the cyclic-6 analogue of
`PT.Sieve.T1OrbitsZMod5` (cyclic-4) and `PT.Sieve.T1AntidiagOrbits`
(cyclic-2).

* `sigma7 : Fin 6 → Fin 6` — the underlying 6-cycle, indexed so that
  `Fin 6 ≃ (1, 3, 2, 6, 4, 5) ∈ (ℤ/7ℤ)^*` in cycle order.
* `sigma7^[6] = id` — the permutation has order `6`.
* `sigma7 i ≠ i` for every `i` — no fixed point.
* `sigma7^[k] ≠ id` for `k ∈ {1, 2, 3, 4, 5}` — minimal period exactly `6`.
* Bridge to `ZMod 7`: the iterates of `3` traverse `{3, 2, 6, 4, 5, 1}`.
* **Dimensional comparison.** Order `6 > 4 > 2` going from `p = 7` to
  `p = 5` to `p = 3`.

## Reference

Theorem T1 (sieve transfer orbit structure) generalised from `p ∈ {3, 5}`
to `p = 7`, Chapter 7 of the monograph.
-/

namespace PT.Sieve

/-! ### The induced permutation `sigma7` on `Fin 6` -/

/-- The permutation of `Fin 6` induced by multiplication by `3` on
    `(ℤ/7ℤ)^*`, written in cycle order. With the indexing
    `0 ↔ 1, 1 ↔ 3, 2 ↔ 2, 3 ↔ 6, 4 ↔ 4, 5 ↔ 5`, the action `x ↦ 3x mod 7`
    becomes `0 ↦ 1 ↦ 2 ↦ 3 ↦ 4 ↦ 5 ↦ 0` on `Fin 6`. -/
def sigma7 : Fin 6 → Fin 6
  | 0 => 1
  | 1 => 2
  | 2 => 3
  | 3 => 4
  | 4 => 5
  | 5 => 0

@[simp] lemma sigma7_zero  : sigma7 0 = 1 := rfl
@[simp] lemma sigma7_one   : sigma7 1 = 2 := rfl
@[simp] lemma sigma7_two   : sigma7 2 = 3 := rfl
@[simp] lemma sigma7_three : sigma7 3 = 4 := rfl
@[simp] lemma sigma7_four  : sigma7 4 = 5 := rfl
@[simp] lemma sigma7_five  : sigma7 5 = 0 := rfl

/-! ### Iterates of `sigma7` -/

/-- Pointwise: `sigma7^[2] i = sigma7 (sigma7 i)`. -/
@[simp] theorem sigma7_iterate_two (i : Fin 6) :
    sigma7^[2] i = sigma7 (sigma7 i) := by
  simp [Function.iterate_succ, Function.comp_apply]

/-- Pointwise: `sigma7^[3] i = sigma7 (sigma7 (sigma7 i))`. -/
theorem sigma7_iterate_three (i : Fin 6) :
    sigma7^[3] i = sigma7 (sigma7 (sigma7 i)) := by
  simp [Function.iterate_succ, Function.comp_apply]

/-- Pointwise: `sigma7^[4] i = sigma7 (sigma7 (sigma7 (sigma7 i)))`. -/
theorem sigma7_iterate_four (i : Fin 6) :
    sigma7^[4] i = sigma7 (sigma7 (sigma7 (sigma7 i))) := by
  simp [Function.iterate_succ, Function.comp_apply]

/-- Pointwise: `sigma7^[5] i = sigma7 (sigma7 (sigma7 (sigma7 (sigma7 i))))`. -/
theorem sigma7_iterate_five (i : Fin 6) :
    sigma7^[5] i = sigma7 (sigma7 (sigma7 (sigma7 (sigma7 i)))) := by
  simp [Function.iterate_succ, Function.comp_apply]

/-- Pointwise: `sigma7^[6] i = i`. -/
@[simp] theorem sigma7_iterate_six_apply (i : Fin 6) :
    sigma7^[6] i = i := by
  fin_cases i <;> decide

/-- **Order 6.** `sigma7^[6] = id`. -/
theorem sigma7_iterate_six : sigma7^[6] = id := by
  funext i
  exact sigma7_iterate_six_apply i

/-! ### No fixed point -/

/-- `sigma7` has no fixed point: `sigma7 i ≠ i` for all `i`. -/
theorem sigma7_no_fixed_point (i : Fin 6) : sigma7 i ≠ i := by
  fin_cases i <;> decide

/-! ### Intermediate iterates are non-trivial (minimal period 6) -/

/-- The first iterate is non-trivial: `sigma7 ≠ id`. -/
theorem sigma7_iterate_one_ne_id : sigma7 ≠ id := by
  intro h
  have h0 : sigma7 0 = (id : Fin 6 → Fin 6) 0 := by rw [h]
  simp at h0

/-- The second iterate is non-trivial: `sigma7^[2] ≠ id`. -/
theorem sigma7_iterate_two_ne_id : sigma7^[2] ≠ id := by
  intro h
  have h0 : sigma7^[2] 0 = (id : Fin 6 → Fin 6) 0 := by rw [h]
  simp at h0

/-- The third iterate is non-trivial: `sigma7^[3] ≠ id`. -/
theorem sigma7_iterate_three_ne_id : sigma7^[3] ≠ id := by
  intro h
  have h0 : sigma7^[3] 0 = (id : Fin 6 → Fin 6) 0 := by rw [h]
  simp [sigma7_iterate_three] at h0

/-- The fourth iterate is non-trivial: `sigma7^[4] ≠ id`. -/
theorem sigma7_iterate_four_ne_id : sigma7^[4] ≠ id := by
  intro h
  have h0 : sigma7^[4] 0 = (id : Fin 6 → Fin 6) 0 := by rw [h]
  simp [sigma7_iterate_four] at h0

/-- The fifth iterate is non-trivial: `sigma7^[5] ≠ id`. -/
theorem sigma7_iterate_five_ne_id : sigma7^[5] ≠ id := by
  intro h
  have h0 : sigma7^[5] 0 = (id : Fin 6 → Fin 6) 0 := by rw [h]
  simp [sigma7_iterate_five] at h0

/-! ### Bijectivity and inverse structure -/

/-- `sigma7` is a bijection of `Fin 6`. -/
theorem sigma7_bijective : Function.Bijective sigma7 := by
  refine ⟨?_, ?_⟩
  · intro i j h
    fin_cases i <;> fin_cases j <;> simp_all
  · intro j
    fin_cases j
    · exact ⟨5, rfl⟩
    · exact ⟨0, rfl⟩
    · exact ⟨1, rfl⟩
    · exact ⟨2, rfl⟩
    · exact ⟨3, rfl⟩
    · exact ⟨4, rfl⟩

/-- `sigma7^[5]` is the inverse of `sigma7`. -/
theorem sigma7_inv_eq_iterate_five (i : Fin 6) :
    sigma7 (sigma7^[5] i) = i := by
  have h : sigma7^[6] i = i := sigma7_iterate_six_apply i
  have hstep : sigma7^[6] i = sigma7 (sigma7^[5] i) := by
    show sigma7^[5 + 1] i = sigma7 (sigma7^[5] i)
    rw [Function.iterate_succ_apply']
  rw [← hstep]; exact h

/-! ### Bridge with `ZMod 7`: powers of `3` -/

/-- **Bridge.** The iterates of multiplication by `3` on `(ℤ/7ℤ)^*` starting
    at `1` traverse `{3, 2, 6, 4, 5, 1}`, matching the six powers
    `3^1, 3^2, 3^3, 3^4, 3^5, 3^6`. -/
theorem pow_three_mod7_one   : (3 : ZMod 7)^1 * 1 = 3 := by decide
theorem pow_three_mod7_two   : (3 : ZMod 7)^2 * 1 = 2 := by decide
theorem pow_three_mod7_three : (3 : ZMod 7)^3 * 1 = 6 := by decide
theorem pow_three_mod7_four  : (3 : ZMod 7)^4 * 1 = 4 := by decide
theorem pow_three_mod7_five  : (3 : ZMod 7)^5 * 1 = 5 := by decide
theorem pow_three_mod7_six   : (3 : ZMod 7)^6 * 1 = 1 := by decide

/-- **3 is a primitive root mod 7.** `(3 : ZMod 7)^6 = 1`. -/
theorem three_pow_six_mod7 : (3 : ZMod 7)^6 = 1 := by decide

/-- **3 is a primitive root mod 7 (minimality).** None of the lower powers
    `(3 : ZMod 7)^k` for `k ∈ {1, 2, 3, 4, 5}` equals `1`. -/
theorem three_primitive_root_mod7 :
    (3 : ZMod 7)^1 ≠ 1
    ∧ (3 : ZMod 7)^2 ≠ 1
    ∧ (3 : ZMod 7)^3 ≠ 1
    ∧ (3 : ZMod 7)^4 ≠ 1
    ∧ (3 : ZMod 7)^5 ≠ 1
    ∧ (3 : ZMod 7)^6 = 1 := by
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_⟩ <;> decide

/-- Embedding `Fin 6 ↪ (ℤ/7ℤ)^*` matching the cycle order
    `0 ↦ 1, 1 ↦ 3, 2 ↦ 2, 3 ↦ 6, 4 ↦ 4, 5 ↦ 5`. -/
def toZMod7 : Fin 6 → ZMod 7
  | 0 => 1
  | 1 => 3
  | 2 => 2
  | 3 => 6
  | 4 => 4
  | 5 => 5

@[simp] lemma toZMod7_zero  : toZMod7 0 = 1 := rfl
@[simp] lemma toZMod7_one   : toZMod7 1 = 3 := rfl
@[simp] lemma toZMod7_two   : toZMod7 2 = 2 := rfl
@[simp] lemma toZMod7_three : toZMod7 3 = 6 := rfl
@[simp] lemma toZMod7_four  : toZMod7 4 = 4 := rfl
@[simp] lemma toZMod7_five  : toZMod7 5 = 5 := rfl

/-- **Conjugation identity.** The combinatorial cycle `sigma7` corresponds, via
    the embedding `toZMod7`, to multiplication by `3` on `(ℤ/7ℤ)^*`. -/
theorem toZMod7_sigma7 (i : Fin 6) :
    toZMod7 (sigma7 i) = 3 * toZMod7 i := by
  fin_cases i <;> decide

/-- The embedding never hits `0`. -/
theorem toZMod7_ne_zero (i : Fin 6) : toZMod7 i ≠ 0 := by
  fin_cases i <;> decide

/-! ### Dimensional comparison with the cyclic-2 and cyclic-4 cases -/

/-- **Dimensional comparison vs the cyclic-4 case.** Whereas the cyclic-4
    permutation `sigma5` of `Fin 4` satisfies `sigma5^[4] = id`, the cyclic-6
    permutation `sigma7` of `Fin 6` satisfies `sigma7^[6] = id` but already
    `sigma7^[4] ≠ id`. The "dimension" of the orbit structure (the order of
    the underlying cycle) is therefore strictly larger here. -/
theorem sigma7_order_strictly_greater_than_four :
    sigma7^[4] ≠ id ∧ sigma7^[6] = id :=
  ⟨sigma7_iterate_four_ne_id, sigma7_iterate_six⟩

/-- **Dimensional comparison vs the cyclic-2 case.** A fortiori,
    `sigma7^[2] ≠ id`, so the cyclic-6 case is also strictly stronger than
    the cyclic-2 antidiagonal case. -/
theorem sigma7_order_strictly_greater_than_two :
    sigma7^[2] ≠ id ∧ sigma7^[6] = id :=
  ⟨sigma7_iterate_two_ne_id, sigma7_iterate_six⟩

/-- **Minimal period is exactly 6.** None of the iterates `sigma7^[k]` for
    `k ∈ {1, 2, 3, 4, 5}` is the identity, while `sigma7^[6] = id`. -/
theorem sigma7_minimal_period_six :
    sigma7 ≠ id
    ∧ sigma7^[2] ≠ id
    ∧ sigma7^[3] ≠ id
    ∧ sigma7^[4] ≠ id
    ∧ sigma7^[5] ≠ id
    ∧ sigma7^[6] = id :=
  ⟨sigma7_iterate_one_ne_id,
   sigma7_iterate_two_ne_id,
   sigma7_iterate_three_ne_id,
   sigma7_iterate_four_ne_id,
   sigma7_iterate_five_ne_id,
   sigma7_iterate_six⟩

/-- **Cascade comparison.** Order grows strictly through `p = 3, 5, 7`:
    `2 < 4 < 6`. This file extends the cyclic structure already recorded for
    `(ℤ/3ℤ)^*` (cf. `PT.Sieve.T1AntidiagOrbits`, order `2`) and
    `(ℤ/5ℤ)^*` (cf. `PT.Sieve.T1OrbitsZMod5`, order `4`). -/
theorem sigma7_dimensional_cascade :
    sigma^[2] = id
    ∧ sigma5^[4] = id ∧ sigma5^[2] ≠ id
    ∧ sigma7^[6] = id ∧ sigma7^[4] ≠ id := by
  refine ⟨?_, sigma5_iterate_four, sigma5_iterate_two_ne_id,
          sigma7_iterate_six, sigma7_iterate_four_ne_id⟩
  funext i
  fin_cases i <;> rfl

/-! ### Headline: orbit structure of the multiplication-by-3 map mod 7 -/

/-- **Headline.** Orbit structure of the multiplication-by-3 map on
    `(ℤ/7ℤ)^*`, encoded as the permutation `sigma7` on `Fin 6`:

    * `sigma7` is a single 6-cycle `(0 1 2 3 4 5)` with no fixed point.
    * The minimal period is exactly `6`: `sigma7^[k] ≠ id` for
      `k ∈ {1, 2, 3, 4, 5}` and `sigma7^[6] = id`.
    * Bridge to `ZMod 7`: `toZMod7 ∘ sigma7 = (3 *) ∘ toZMod7`, i.e. `sigma7`
      conjugates to the multiplication-by-`3` map on `(ℤ/7ℤ)^*`.
    * `3` is a primitive root mod `7`: `3^6 = 1` and `3^k ≠ 1` for
      `k ∈ {1, ..., 5}`.
    * **Dimensional cascade.** Strictly stronger than `T1OrbitsZMod5`
      (order `4`) and `T1AntidiagOrbits` (order `2`). -/
theorem sigma7_orbit_structure :
    -- Cycle properties:
    sigma7^[6] = id
    ∧ (∀ i, sigma7 i ≠ i)
    ∧ sigma7 0 = 1 ∧ sigma7 1 = 2 ∧ sigma7 2 = 3
    ∧ sigma7 3 = 4 ∧ sigma7 4 = 5 ∧ sigma7 5 = 0
    -- Minimal period 6:
    ∧ sigma7 ≠ id ∧ sigma7^[2] ≠ id ∧ sigma7^[3] ≠ id
    ∧ sigma7^[4] ≠ id ∧ sigma7^[5] ≠ id
    -- Bridge to ZMod 7:
    ∧ (∀ i, toZMod7 (sigma7 i) = 3 * toZMod7 i)
    -- Bijectivity:
    ∧ Function.Bijective sigma7 :=
  ⟨sigma7_iterate_six,
   sigma7_no_fixed_point,
   rfl, rfl, rfl, rfl, rfl, rfl,
   sigma7_iterate_one_ne_id,
   sigma7_iterate_two_ne_id,
   sigma7_iterate_three_ne_id,
   sigma7_iterate_four_ne_id,
   sigma7_iterate_five_ne_id,
   toZMod7_sigma7,
   sigma7_bijective⟩

end PT.Sieve
