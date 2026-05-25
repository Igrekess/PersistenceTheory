/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Bridge.PartitionFunctionFactorisation
import PT.Information.T6cG1Autonomous

/-!
# Gibbs distribution on a finite spectrum

This module defines the **Gibbs (Boltzmann) probability distribution**
on a finite energy spectrum, the canonical maximum-entropy distribution
under an energy constraint (Jaynes 1957).

## What this module accomplishes

* **Definition `gibbsState β E hn`** — for a finite spectrum
  `E : Fin n → ℝ` (with `n > 0`) at inverse temperature `β`, the
  Gibbs distribution `p_i = exp(-β · E_i) / Z(β, E)` packaged as
  a `Simplex n` (kernel-verified probability distribution).

* **Theorem `gibbsState_prob_eq`** — pointwise formula
  `(gibbsState β E hn).prob i = exp(-β · E i) / Z(β, E)`.

* **Theorem `gibbsState_pos`** — every Gibbs weight is strictly
  positive.

The bridge to the PT heat-kernel cutoff is the identification

  `f_PT(u) = exp(-β u) = Z(β, [u]) / Z(β, ∅) · Z(β, ∅) = Z(β, [u])`

(for a one-element spectrum, the Gibbs partition function evaluated
at energy `u` is simply `exp(-β u)`). This is the basis for deriving
`PT_cutoff_is_heat_kernel` from Jaynes-Gibbs structure.

## Strategy

The proof of `gibbsState_sums_to_one` is the standard normalisation
identity:
`(∑_i exp(-β · E_i) / Z) = (∑_i exp(-β · E_i)) / Z = Z/Z = 1`.

The factorisation of the Gibbs state on tensor products is a corollary
of `partitionFunction_prod_factorises`.

## References

* Gibbs, *Elementary Principles in Statistical Mechanics* (1902).
* Jaynes, *Information theory and statistical mechanics*, Phys. Rev.
  **106** (1957), 620-630.
* Cover & Thomas, *Elements of Information Theory*, §12.2.
-/

namespace PT.Bridge

open Real Finset PT.Information
open scoped BigOperators

/-! ## Gibbs distribution -/

/-- **Gibbs (Boltzmann) distribution** on a finite spectrum `E : Fin n → ℝ`
    at inverse temperature `β`. Each state has probability
    `exp(-β · E_i) / Z(β, E)` where `Z` is the partition function.

    Packaged as a `Simplex n` (kernel-verified probability distribution):
    * each prob `p_i > 0` (since `exp > 0`),
    * sum equals 1 (definition of `Z`).

    Hypothesis `hn : 0 < n` ensures `Z > 0` so the division is well-
    defined.

    For the PT bifurcation `ST_F = (ℂ², ℂ², m·σ_x, σ_x)` with degenerate
    spectrum `Spec(D_F²) = {m², m²}` (`D_F_PT_sq`), the Gibbs distribution
    is uniform `(1/2, 1/2)` at every `β`, reflecting the equiprobable
    weight of the two bifurcation branches `q_+` and `q_-`. -/
noncomputable def gibbsState {n : ℕ} (β : ℝ) (E : Fin n → ℝ) (hn : 0 < n) :
    Simplex n where
  prob i := Real.exp (-β * E i) / partitionFunction β E
  nonneg i := by
    apply div_nonneg
    · exact (Real.exp_pos _).le
    · exact partitionFunction_nonneg β E
  sums_to_one := by
    -- ∑_i exp(-β · E i) / Z = (∑_i exp(-β · E i)) / Z = Z / Z = 1
    rw [← Finset.sum_div]
    -- Goal: (∑ i, exp(-β · E i)) / Z = 1
    -- LHS numerator = Z by definition.
    have hZ : (∑ i, Real.exp (-β * E i)) = partitionFunction β E := rfl
    rw [hZ]
    -- Goal: Z / Z = 1
    have hZ_pos : 0 < partitionFunction β E := partitionFunction_pos β E hn
    field_simp

/-- **Pointwise formula** for the Gibbs distribution. -/
@[simp]
theorem gibbsState_prob_eq {n : ℕ} (β : ℝ) (E : Fin n → ℝ) (hn : 0 < n)
    (i : Fin n) :
    (gibbsState β E hn).prob i = Real.exp (-β * E i) / partitionFunction β E :=
  rfl

/-- **Strict positivity** of every Gibbs probability. -/
theorem gibbsState_pos {n : ℕ} (β : ℝ) (E : Fin n → ℝ) (hn : 0 < n) (i : Fin n) :
    0 < (gibbsState β E hn).prob i := by
  rw [gibbsState_prob_eq]
  apply div_pos
  · exact Real.exp_pos _
  · exact partitionFunction_pos β E hn

/-! ## Gibbs partition function as a cutoff function

For a single-state spectrum `E : Fin 1 → ℝ` with `E 0 = u`, the
partition function is `Z(β, [u]) = exp(-β · u)`. This identifies the
PT heat-kernel cutoff `f_PT(u) = exp(-β · u)` with the Gibbs
partition function of a one-element spectrum.

More generally, for the PT bifurcation `ST_F` with degenerate spectrum
`{m², m²}` (`D_F_PT_sq`), `Z(β, [m², m²]) = 2 · exp(-β · m²)`. The
factor `2 = N_b_PT_nat` is the Tr_F(I_F) degeneracy. -/

/-- **Partition function on a one-element spectrum.** For `E : Fin 1 → ℝ`
    with `E 0 = u`, `Z(β, [u]) = exp(-β · u)`. This is the simplest
    case identifying the PT cutoff `f_PT(u) = exp(-β · u)` with a
    Gibbs partition function. -/
theorem partitionFunction_single (β u : ℝ) :
    partitionFunction (n := 1) β (fun _ => u) = Real.exp (-β * u) := by
  unfold partitionFunction
  rw [Fin.sum_univ_one]

/-- **Partition function on a degenerate two-element spectrum.** For
    `E : Fin 2 → ℝ` with `E i = u` for all `i`, `Z(β, [u, u]) = 2 · exp(-β · u)`.
    This is the case relevant to the PT bifurcation `ST_F`:
    `D_F² = m² · I` has spectrum `{m², m²}` (degenerate), and the
    factor `2 = N_b_PT_nat = (Tr_F(I_F)).re` is the **dimension of `H_F`**. -/
theorem partitionFunction_degenerate_two (β u : ℝ) :
    partitionFunction (n := 2) β (fun _ => u) = 2 * Real.exp (-β * u) := by
  unfold partitionFunction
  rw [Fin.sum_univ_two]
  ring

end PT.Bridge
