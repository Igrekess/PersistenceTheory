/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import Mathlib.Analysis.SpecialFunctions.Exp
import Mathlib.Analysis.Complex.Exponential
import Mathlib.Algebra.BigOperators.Pi
import Mathlib.Data.Fintype.BigOperators

/-!
# Multiplicativity of the partition function on tensor products

This module provides the algebraic core that justifies the multilinear
factorisation of the Chamseddine-Connes spectral cutoff (already
proven in `PT/Bridge/HeatKernelPostulate.lean` from the postulate).

## What this module accomplishes

* **Definition `partitionFunction β E`** — the canonical Gibbs partition
  function `Z(β, E) = ∑_i exp(-β · E i)` for a finite-spectrum energy
  function `E : Fin n → ℝ` at inverse temperature `β`.

* **Theorem `partitionFunction_pos`** — `Z > 0` since each summand is
  strictly positive.

* **Theorem `partitionFunction_factorises`** — kernel-verified
  **multiplicativity of `Z` on tensor products**:

      `Z(β, E_1 ⊕ E_2) = Z(β, E_1) · Z(β, E_2)`

  where `(E_1 ⊕ E_2)(i,j) := E_1(i) + E_2(j)` is the additive energy on
  the tensor product `Fin m × Fin n`.

  This is the **classical statistical-mechanics theorem** that
  underlies the SJ G1 factorisation of the spectral cutoff: independent
  subsystems have multiplicative partition functions.

## Epistemic significance

The multilinear factorisation of the PT spectral cutoff (the axiom
`SJG1_spectral_cutoff_factorises` of `ShoreJohnsonG1Spectral.lean`,
derived as theorem from heat kernel in `HeatKernelPostulate.lean`) has
its **classical statistical-mechanics origin** in the multiplicativity
of the Gibbs partition function on tensor products of independent
subsystems. This module records that fact in Lean.

The Lean-formalised multiplicativity is `partitionFunction_factorises`,
which together with the heat-kernel identification
`PT_cutoff_is_heat_kernel` (postulated, NCG standard) yields the
factorisation as a theorem.

## Strategy

The proof is a direct double-sum manipulation using `Fintype.sum_prod_type`
+ `Real.exp_sum` (or the simpler `Real.exp_add` distributed via
distributivity of `*` over the inner sum):

```
Z(β, E_1 ⊕ E_2)
  = ∑_{(i,j)} exp(-β · (E_1 i + E_2 j))
  = ∑_{(i,j)} exp(-β · E_1 i) · exp(-β · E_2 j)        -- exp(a+b) = exp(a)·exp(b)
  = (∑_i exp(-β · E_1 i)) · (∑_j exp(-β · E_2 j))      -- distributivity
  = Z(β, E_1) · Z(β, E_2).
```

## References

* Gibbs, *Elementary Principles in Statistical Mechanics* (1902),
  §I-II (partition function, Boltzmann distribution).
* Jaynes, *Information theory and statistical mechanics*, Phys. Rev.
  **106** (1957), 620-630.
* Cover & Thomas, *Elements of Information Theory*, §12.2 (maximum
  entropy distributions; the discrete Gibbs theorem).
* Connes & Marcolli, *Noncommutative Geometry, Quantum Fields and
  Motives*, AMS Colloquium Publications **55** (2008), §1.18 (KMS
  states on spectral algebras).
-/

namespace PT.Bridge

open Real Finset
open scoped BigOperators

/-! ## The partition function -/

/-- **Gibbs partition function** at inverse temperature `β` for a finite
    spectrum `E : Fin n → ℝ`:

      `Z(β, E) := ∑_{i ∈ Fin n} exp(-β · E i)`.

    This is the canonical normalising constant of the Gibbs distribution
    `p_i = exp(-β · E i) / Z(β, E)`. It is the unique normalised maximum-
    entropy distribution on `Fin n` with energy constraint
    `⟨E⟩ = -d log Z / dβ` (Jaynes 1957, Cover-Thomas §12.2). -/
noncomputable def partitionFunction {n : ℕ} (β : ℝ) (E : Fin n → ℝ) : ℝ :=
  ∑ i, Real.exp (-β * E i)

/-! ## Basic properties -/

/-- **Strict positivity of `Z`.** Each summand `exp(-β · E i) > 0`, hence
    the sum is strictly positive for `n ≥ 1`.

    For `n = 0` the sum is empty and `Z = 0`. The PT applications use
    `n = N_b_PT_nat = 2 > 0`. -/
theorem partitionFunction_pos {n : ℕ} (β : ℝ) (E : Fin n → ℝ) (hn : 0 < n) :
    0 < partitionFunction β E := by
  unfold partitionFunction
  -- Sum of strictly positive terms over a non-empty finite index.
  have h_nonempty : (Finset.univ : Finset (Fin n)).Nonempty := by
    rw [Finset.univ_nonempty_iff]
    exact ⟨⟨0, hn⟩⟩
  apply Finset.sum_pos
  · intro i _
    exact Real.exp_pos _
  · exact h_nonempty

/-- **Non-negativity of `Z`.** Always true (trivially for `n = 0` since
    `Z = 0`). -/
theorem partitionFunction_nonneg {n : ℕ} (β : ℝ) (E : Fin n → ℝ) :
    0 ≤ partitionFunction β E := by
  unfold partitionFunction
  apply Finset.sum_nonneg
  intro i _
  exact (Real.exp_pos _).le

/-! ## Additive energy on tensor products -/

/-- **Additive tensor-product energy.** For `E_1 : Fin m → ℝ` and
    `E_2 : Fin n → ℝ`, the energy on the tensor product
    `Fin m × Fin n` is the sum:

      `(E_1 ⊕ E_2)(i, j) := E_1(i) + E_2(j)`.

    This encodes the standard physical postulate that energies of
    independent subsystems add. -/
def addTensor {m n : ℕ} (E_1 : Fin m → ℝ) (E_2 : Fin n → ℝ) :
    Fin m × Fin n → ℝ :=
  fun ij => E_1 ij.1 + E_2 ij.2

/-! ## Multiplicativity theorem (Gibbs partition function on tensors)

This is the central theorem of this module: the partition function of
the additive tensor-product energy factorises as the product of partition
functions. It is the kernel-verified analogue of the
`SJG1_spectral_cutoff_factorises_from_heat_kernel` theorem in
`HeatKernelPostulate.lean`, but stated directly at the level of partition
functions (without reference to a cutoff function `f`).

This is the **classical statistical-mechanics origin** of SJ G1 spectral
factorisation: independent subsystems have multiplicative partition
functions. -/

/-- **Partition function multiplicativity on Fin m × Fin n.**

    For any inverse temperature `β` and any spectra
    `E_1 : Fin m → ℝ`, `E_2 : Fin n → ℝ`:

        `∑_{(i,j)} exp(-β · (E_1 i + E_2 j))
           = (∑_i exp(-β · E_1 i)) · (∑_j exp(-β · E_2 j))`.

    **Proof.** Direct double-sum manipulation:
    `exp(-β(a+b)) = exp(-βa) · exp(-βb)` (distributivity of `exp`),
    then `Fintype.sum_prod_type` + `Finset.sum_mul_sum`. -/
theorem partitionFunction_prod_factorises
    {m n : ℕ} (β : ℝ) (E_1 : Fin m → ℝ) (E_2 : Fin n → ℝ) :
    (∑ ij : Fin m × Fin n, Real.exp (-β * (E_1 ij.1 + E_2 ij.2)))
      = (∑ i, Real.exp (-β * E_1 i)) * (∑ j, Real.exp (-β * E_2 j)) := by
  -- Step 1: convert LHS to a nested sum (Fintype.sum_prod_type).
  -- Step 2: convert RHS to a nested sum (Finset.sum_mul_sum / Fintype.sum_mul_sum).
  -- Step 3: term-by-term equality via mul_add + Real.exp_add.
  rw [Fintype.sum_prod_type, Fintype.sum_mul_sum]
  refine Finset.sum_congr rfl fun i _ => ?_
  refine Finset.sum_congr rfl fun j _ => ?_
  rw [mul_add, Real.exp_add]

/-- **Partition function multiplicativity (named form).** Using the
    `addTensor` notation for the additive tensor-product energy. -/
theorem partitionFunction_addTensor
    {m n : ℕ} (β : ℝ) (E_1 : Fin m → ℝ) (E_2 : Fin n → ℝ) :
    (∑ ij : Fin m × Fin n, Real.exp (-β * (addTensor E_1 E_2 ij)))
      = partitionFunction β E_1 * partitionFunction β E_2 := by
  unfold partitionFunction addTensor
  exact partitionFunction_prod_factorises β E_1 E_2

/-! ## Multilinear factorisation (k-fold partition function)

The multilinear factorisation `Z(∑ E_i) = ∏ Z(E_i)` for k independent
subsystems is the iterated form of `partitionFunction_prod_factorises`.

For the PT application, where the cutoff is the partition function of
the heat kernel and the spectrum factorises via CRT
(`PT/Stochastic/CRTFactorizationT2.lean`), this gives the full
multilinear factorisation `f(∑u) = ∏ f(u)` that was previously the
axiom `SJG1_spectral_cutoff_factorises`. -/

/-- **Multilinear factorisation of `exp(-β · _)` on sums.** Specialisation
    of `Real.exp_sum` with the factor `-β`. This is the abstract form
    that, applied to the partition function `Z(β, E)` and identifying
    `f(u) = exp(-β u)`, gives the multilinear cutoff factorisation
    `f(∑u_i) = ∏ f(u_i)`. -/
theorem exp_neg_mul_sum_factorises (β : ℝ) {k : ℕ} (u : Fin k → ℝ) :
    Real.exp (-β * (∑ i, u i)) = ∏ i, Real.exp (-β * u i) := by
  rw [Finset.mul_sum]
  exact Real.exp_sum _ _

end PT.Bridge
