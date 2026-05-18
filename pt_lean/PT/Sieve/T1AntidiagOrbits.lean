/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Sieve.T3Antidiagonal
import PT.Stochastic.SpectralDominance
import Mathlib.Tactic

/-!
# T1 — Orbits of the antidiagonal transfer `T_3`

The sieve-level transfer matrix `T_3 = !![0,1;1,0]` (see
`PT.Sieve.T3Antidiagonal`) is an involution: `T_3^2 = I`. Hence the action
of `T_3` on the two non-zero residues `(ℤ/3ℤ)^* = {1, 2} ≃ Fin 2` is a
permutation of order dividing `2`, i.e. either the identity or the unique
transposition. Since `T_3` is antidiagonal, the induced permutation is the
non-trivial transposition `σ = (0 1)`.

This file formalises the **orbit structure** of `T_3` acting on `Fin 2`:

* `σ : Fin 2 → Fin 2` — the underlying permutation, `0 ↦ 1`, `1 ↦ 0`.
* `σ ∘ σ = id` — the permutation is an involution.
* `σ` is a single 2-cycle: it has no fixed point and exchanges `0` and `1`.
* `orbit_T3 i j` — reachability predicate: `j` is reachable from `i` under
  iteration of `σ`.
* The orbit of any `i` is the whole space `{0, 1}`.
* **Iteration parity.** After an even number of steps we are back at `i`;
  after an odd number we are at `σ i`.

These facts are the discrete (combinatorial) counterpart of the spectral
decomposition of `T_3` recorded in `PT.Stochastic.SpectralDominance`:
the matrix-level statement `T_3^2 = I` is reflected at the level of
permutations by `σ ∘ σ = id`, and the two eigenmodes `v_±` correspond to
the symmetric/antisymmetric combinations of the two orbit points.

## Reference

Theorem T1 (forbidden transitions and orbit structure of the sieve transfer),
Chapter 7 of the monograph, and `PT_ARTICLES/PT_MATHEMATICS/M1`.
-/

namespace PT.Sieve

open Matrix

/-! ### The induced permutation `σ` on `Fin 2` -/

/-- The permutation of `Fin 2` induced by the antidiagonal transfer `T_3`:
    `σ 0 = 1`, `σ 1 = 0`. -/
def sigma : Fin 2 → Fin 2
  | 0 => 1
  | 1 => 0

@[simp] lemma sigma_zero : sigma 0 = 1 := rfl
@[simp] lemma sigma_one  : sigma 1 = 0 := rfl

/-- `σ` is an involution: `σ ∘ σ = id`. -/
theorem sigma_involution : sigma ∘ sigma = id := by
  funext i
  fin_cases i <;> rfl

/-- Pointwise involution: `σ (σ i) = i`. -/
@[simp] theorem sigma_sigma (i : Fin 2) : sigma (sigma i) = i := by
  fin_cases i <;> rfl

/-- `σ` is a bijection of `Fin 2` (it is its own inverse). -/
theorem sigma_bijective : Function.Bijective sigma := by
  refine ⟨?_, ?_⟩
  · intro i j h
    have : sigma (sigma i) = sigma (sigma j) := by rw [h]
    simpa using this
  · intro j
    exact ⟨sigma j, by simp⟩

/-! ### `σ` is a single 2-cycle (no fixed point) -/

/-- `σ` has no fixed point: `σ i ≠ i` for all `i`. -/
theorem sigma_no_fixed_point (i : Fin 2) : sigma i ≠ i := by
  fin_cases i <;> decide

/-- `σ` is a 2-cycle: it is non-trivial (`σ ≠ id`) and squares to the
    identity. This is the cycle decomposition `σ = (0 1)`. -/
theorem sigma_is_two_cycle :
    sigma ≠ id ∧ sigma ∘ sigma = id := by
  refine ⟨?_, sigma_involution⟩
  intro h
  have h0 : sigma 0 = (id : Fin 2 → Fin 2) 0 := by rw [h]
  simp at h0

/-- **Cycle structure.** `σ` exchanges `0` and `1`: explicitly, `σ 0 = 1`
    and `σ 1 = 0`, which is the standard transposition `(0 1)`. -/
theorem sigma_cycle_decomposition :
    sigma 0 = 1 ∧ sigma 1 = 0 :=
  ⟨rfl, rfl⟩

/-! ### Action of `T_3` on the standard basis -/

/-- The standard basis vector `δ_i : Fin 2 → ℝ` defined by `δ_i j = 1` if
    `j = i` and `0` otherwise. -/
def delta (i : Fin 2) : Fin 2 → ℝ := fun j => if j = i then 1 else 0

@[simp] lemma delta_self (i : Fin 2) : delta i i = 1 := by
  simp [delta]

lemma delta_zero_zero : delta 0 0 = 1 := by simp
lemma delta_zero_one  : delta 0 1 = 0 := by simp [delta]
lemma delta_one_zero  : delta 1 0 = 0 := by simp [delta]
lemma delta_one_one   : delta 1 1 = 1 := by simp

/-- **Matrix–permutation correspondence.** Acting by `T_3` on the indicator
    vector `δ_i` produces the indicator vector `δ_{σ i}`. -/
theorem T3_mulVec_delta (i : Fin 2) :
    T3.mulVec (delta i) = delta (sigma i) := by
  ext j
  fin_cases i <;> fin_cases j <;>
    simp [T3, delta, sigma, Matrix.mulVec, dotProduct]

/-! ### Iterated action and parity -/

/-- Iterating `σ` an even number of times gives the identity. -/
theorem sigma_iterate_even (i : Fin 2) (n : ℕ) :
    sigma^[2 * n] i = i := by
  induction n with
  | zero => simp
  | succ k ih =>
    have hstep : 2 * (k + 1) = 2 * k + 2 := by ring
    rw [hstep, Function.iterate_add_apply]
    simp [Function.comp_apply, ih, sigma_sigma]

/-- Iterating `σ` an odd number of times gives `σ i`. -/
theorem sigma_iterate_odd (i : Fin 2) (n : ℕ) :
    sigma^[2 * n + 1] i = sigma i := by
  rw [Function.iterate_add_apply, sigma_iterate_even]
  rfl

/-! ### Reachability and orbits -/

/-- **Reachability predicate.** `orbit_T3 i j` is `true` iff `j` is reachable
    from `i` by a (possibly empty) finite sequence of `σ`-iterations.

    Since `σ` is an involution, only `j = i` (0 steps) and `j = σ i`
    (1 step) are reachable; this is captured combinatorially. -/
def orbit_T3 (i j : Fin 2) : Bool :=
  decide (j = i ∨ j = sigma i)

@[simp] lemma orbit_T3_refl (i : Fin 2) : orbit_T3 i i = true := by
  simp [orbit_T3]

lemma orbit_T3_sigma (i : Fin 2) : orbit_T3 i (sigma i) = true := by
  simp [orbit_T3]

/-- The orbit of `0` under `σ` is the whole space: both `0` and `1`
    are reachable. -/
theorem orbit_from_zero :
    orbit_T3 0 0 = true ∧ orbit_T3 0 1 = true := by
  refine ⟨by simp, ?_⟩
  have : sigma 0 = 1 := rfl
  simp [orbit_T3, this]

/-- The orbit of `1` under `σ` is the whole space: both `0` and `1`
    are reachable. -/
theorem orbit_from_one :
    orbit_T3 1 0 = true ∧ orbit_T3 1 1 = true := by
  refine ⟨?_, by simp⟩
  have : sigma 1 = 0 := rfl
  simp [orbit_T3, this]

/-- **Orbits cover everything.** From any starting point `i ∈ Fin 2`, the
    orbit of `σ` reaches every element of `Fin 2`. -/
theorem orbit_T3_complete (i j : Fin 2) : orbit_T3 i j = true := by
  fin_cases i <;> fin_cases j <;> decide

/-- **Orbit symmetry.** Reachability is symmetric: if `j` is reachable from
    `i`, then `i` is reachable from `j`. -/
theorem orbit_T3_symm (i j : Fin 2) :
    orbit_T3 i j = orbit_T3 j i := by
  fin_cases i <;> fin_cases j <;> decide

/-! ### Iteration parity in terms of `T_3^n`-action on `δ_i` -/

/-- **Even-step return.** After `2n` iterations of `σ`, every point returns
    to itself. -/
theorem reach_even (i : Fin 2) (n : ℕ) :
    sigma^[2 * n] i = i :=
  sigma_iterate_even i n

/-- **Odd-step displacement.** After `2n+1` iterations of `σ`, every point
    is sent to its `σ`-image. -/
theorem reach_odd (i : Fin 2) (n : ℕ) :
    sigma^[2 * n + 1] i = sigma i :=
  sigma_iterate_odd i n

/-- **General parity dichotomy.** For every `n`, `σ^[n] i` is either `i`
    (when `n` is even) or `σ i` (when `n` is odd). -/
theorem sigma_iterate_parity (i : Fin 2) (n : ℕ) :
    sigma^[n] i = i ∨ sigma^[n] i = sigma i := by
  rcases Nat.even_or_odd n with h | h
  · left
    obtain ⟨k, hk⟩ := h
    have hk' : n = 2 * k := by omega
    rw [hk', sigma_iterate_even]
  · right
    obtain ⟨k, hk⟩ := h
    have hk' : n = 2 * k + 1 := by omega
    rw [hk', sigma_iterate_odd]

/-! ### Headline: orbit structure of the antidiagonal transfer -/

/-- **Headline.** The orbit structure of the sieve transfer `T_3` on the
    two non-zero residues `Fin 2`:

    * The induced permutation `σ` is a 2-cycle `(0 1)` with no fixed
      point.
    * `σ` is its own inverse: `σ ∘ σ = id`.
    * The orbit of every starting point is the entire space `{0, 1}`.
    * Iteration is determined by parity: `σ^[2n] i = i` and
      `σ^[2n+1] i = σ i`.
    * Matrix-level correspondence: `T_3 · δ_i = δ_{σ i}`. -/
theorem T3_orbit_structure :
    -- Permutation properties:
    sigma ∘ sigma = id
    ∧ (∀ i, sigma i ≠ i)
    ∧ sigma 0 = 1 ∧ sigma 1 = 0
    -- Full reachability:
    ∧ (∀ i j, orbit_T3 i j = true)
    -- Parity dichotomy:
    ∧ (∀ i n, sigma^[2 * n] i = i)
    ∧ (∀ i n, sigma^[2 * n + 1] i = sigma i)
    -- Matrix-level correspondence:
    ∧ (∀ i, T3.mulVec (delta i) = delta (sigma i)) :=
  ⟨sigma_involution, sigma_no_fixed_point, rfl, rfl,
   orbit_T3_complete, sigma_iterate_even, sigma_iterate_odd,
   T3_mulVec_delta⟩

end PT.Sieve
