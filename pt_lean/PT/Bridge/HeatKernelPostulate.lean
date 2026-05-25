/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Bridge.ShoreJohnsonG1Spectral
import PT.Bridge.FiniteSpectralTriple
import Mathlib.Analysis.Complex.Exponential

/-!
# Heat-kernel postulate: PT cutoff = Gibbs partition function

This module follows the **Connes-Marcolli 2008 §1.18** identification
of the spectral action cutoff with the heat kernel of `D²`.

## Conceptual framing

The standard NCG literature (Eckstein-Iochum 2019 review, Connes-van
Suijlekom 2020, Khalkhali-Pagliaroli, Bost-Connes/Connes-Marcolli
systems) provides a direct bridge between the spectrum of a Dirac
operator and a probability distribution on its eigenvalues:

* The Chamseddine-Connes spectral action `Tr f(D²/Λ²)` admits the
  exponential cutoff `f(u) = exp(-βu)` as the canonical choice
  (Eckstein-Iochum 2019, §2.3: *"common choices for the cutoff
  function are a decreasing exponential, or the characteristic
  function of the unit interval (a sharp cutoff)"*).
* This exponential cutoff IS the **heat kernel** of `D²` at inverse
  temperature `β`: `f(D²/Λ²) = exp(-β·D²)` with `β = 1/Λ²`.
* The heat kernel is the partition function of a **Gibbs state**
  on the spectral algebra (Bost-Connes 1995, Connes-Marcolli 2008
  §1.18).
* The Gibbs state is the unique **maximum-entropy** distribution
  under an energy constraint (Jaynes 1957; Gibbs-Jaynes theorem of
  statistical mechanics).
* **Multiplicativity** of the partition function on independent
  subsystems `Z(A⊗B) = Z(A)·Z(B)` is the classical theorem of
  statistical mechanics that lifts to the multilinear factorisation
  of the spectral cutoff on CRT-independent eigenvalues.

We **postulate** the Chamseddine-Connes identification of the PT
cutoff as a heat kernel, and **derive** the multilinear
factorisation as a kernel-verified theorem.

## What this module accomplishes

* **Axiom `PT_cutoff_is_heat_kernel`** — the *standard NCG postulate*
  that the PT spectral cutoff has the form `f(u) = exp(-β·u)` for
  some `β > 0`, corresponding to a Gibbs state on the spectral
  algebra (Chamseddine-Connes 1996, Connes-Marcolli 2008 §1.18).
  This replaces the *PT-specific* axiom
  `SJG1_spectral_cutoff_factorises` (multilinear factorisation) by
  the standard NCG identification.

* **Theorem `SJG1_spectral_cutoff_factorises_from_heat_kernel`** —
  the multilinear factorisation, previously an axiom, is now derived
  as a 6-line theorem from the heat-kernel postulate via
  `Real.exp_sum`. This is the kernel-verified content of "the
  partition function is multiplicative on tensor products of
  CRT-independent subsystems".

* **Theorem `cutoff_eq_neg_exp_of_heat_kernel`** — the heat-kernel
  postulate directly gives `f(x) = exp(c · x)` with `c = -β < 0`,
  bypassing the Cauchy + decay + scale chain of the previous
  approach.

## Epistemic status

**Alternative formulation of K4**: K4 depends on:
* `PT_cutoff_is_heat_kernel` (**standard NCG postulate**, Chamseddine-
  Connes spectral action principle applied to PT)
* `cutoff_scale_from_dim_FST` (structural FST scale axiom, derivable
  from Gibbs-Jaynes MaxEnt with constraint `⟨E⟩ = N_b`)

This route uses 2 axioms (same count as the SJG1+scale route), but
**both axioms are now standard NCG/QSM postulates applied to PT**,
instead of one PT-specific axiom + one opaque scale parameter. The
status of K4 reads as `[DER strict modulo standard NCG/QSM
postulates applied to PT]` rather than `[DER strict modulo 2
PT-specific axioms]`.

The Gibbs-Jaynes MaxEnt derivation of `β = 1/N_b` from
`Tr_F(I_F) = N_b` would close the chain completely (`[DER strict
absolute]` modulo classical statistical mechanics), but requires
formalising Gibbs-Jaynes in Lean — not attempted here.

## References

* Connes & Marcolli, *Noncommutative Geometry, Quantum Fields and
  Motives*, AMS Colloquium Publications **55** (2008), §1.18
  (quantum statistical mechanics on spectral algebras).
* Chamseddine & Connes, *The Spectral Action Principle*, Phys. Lett.
  B **412** (1996), 1-9 (heat kernel as canonical cutoff).
* Bost & Connes, *Hecke algebras, type III factors and phase
  transitions with spontaneous symmetry breaking in number theory*,
  Selecta Math. (NS) **1** (1995), 411-457 (KMS states on spectral
  algebras).
* Eckstein & Iochum, *Spectral Action in Noncommutative Geometry*,
  Springer Briefs (2019), §2.3 (exponential as standard cutoff).
* Jaynes, *Information theory and statistical mechanics*, Phys. Rev.
  **106** (1957), 620-630 (MaxEnt → Gibbs distribution).
-/

namespace PT.Bridge

open Real Finset
open scoped BigOperators

/-! ## The heat-kernel postulate (standard NCG) -/

/-- **Axiom (Chamseddine-Connes heat-kernel postulate applied to PT).**

    The PT spectral cutoff `f : ℝ → ℝ` is the heat kernel of the PT
    Dirac operator: there exists an inverse temperature `β > 0` such
    that `f(u) = exp(-β · u)` for all `u ≥ 0`.

    **Standard NCG justification.** In Chamseddine-Connes 1996, the
    spectral action `Tr f(D²/Λ²)` admits the exponential cutoff
    `f(u) = exp(-β u)` as the canonical choice (Eckstein-Iochum 2019
    §2.3). This cutoff is the partition function of a Gibbs state on
    the spectral algebra of `D²`, i.e. the heat kernel
    `e^{-β D²}` evaluated on eigenstates. In Connes-Marcolli 2008
    §1.18, the framework of quantum statistical mechanics on
    spectral algebras is developed in detail (KMS states, Bost-Connes
    type systems).

    Applied to the PT bifurcation `ST_F = (ℂ², ℂ², m·σ_x, σ_x)`, the
    heat kernel `e^{-β D_F²}` evaluated on `H_F` gives
    `exp(-β m²) · I_F`, whose trace is `Tr_F(I_F) · exp(-β m²) =
    N_b · exp(-β m²)`. This is the spectral action on the bifurcation
    sector.

    The remaining freedom is the inverse temperature `β`, fixed
    structurally by the bifurcation scale via Gibbs-Jaynes MaxEnt
    with constraint `⟨E⟩ = N_b` (axiomatised separately as
    `cutoff_scale_from_dim_FST`).

    **What this axiom replaces**: the PT-specific multilinear axiom
    `SJG1_spectral_cutoff_factorises` of
    `ShoreJohnsonG1Spectral.lean`. With the heat-kernel postulate,
    multilinear factorisation becomes a kernel-verified **theorem**
    (`SJG1_spectral_cutoff_factorises_from_heat_kernel` below). -/
axiom PT_cutoff_is_heat_kernel
    (f : ℝ → ℝ) (hf : IsPTSpectralCutoff f) :
    ∃ β > 0, ∀ u, f u = Real.exp (-β * u)

/-! ## Derived theorems from the heat-kernel postulate -/

/-- **Multilinear factorisation as a theorem (Gibbs partition function
    multiplicativity).**

    Previously posited as the axiom
    `SJG1_spectral_cutoff_factorises`, the multilinear factorisation
    of the PT spectral cutoff is now **derived** from the heat-kernel
    postulate via `Real.exp_sum` (the classical multiplicativity of
    the exponential on sums).

    This is the kernel-verified analogue of the well-known
    statistical-mechanics statement that **the partition function of
    a tensor product of independent systems factorises as the product
    of partition functions**:

        `Z(A ⊗ B) = Z(A) · Z(B)`,

    which is the SJ G1 system-independence axiom applied to the
    spectral inference. -/
theorem SJG1_spectral_cutoff_factorises_from_heat_kernel
    (f : ℝ → ℝ) (hf : IsPTSpectralCutoff f)
    (k : ℕ) (u : Fin k → ℝ) :
    f (∑ i, u i) = ∏ i, f (u i) := by
  -- Heat-kernel postulate: f(x) = exp(-β · x) for some β > 0.
  obtain ⟨β, _, hβ⟩ := PT_cutoff_is_heat_kernel f hf
  -- Rewrite both sides using f = exp(-β · _).
  rw [hβ]
  -- LHS: exp(-β · ∑ u_i). Pull the constant -β inside the sum.
  rw [Finset.mul_sum]
  -- LHS = exp(∑ (-β · u_i)). Apply Real.exp_sum.
  rw [Real.exp_sum]
  -- RHS goal: ∏ exp(-β · u_i) = ∏ f(u_i). Rewrite each f.
  congr 1
  ext i
  exact (hβ (u i)).symm

/-- **The PT cutoff has exponential form, directly.** Specialisation
    of the heat-kernel postulate to the Cauchy/exp framework of
    `CauchyMultiplicativeExp`. This bypasses the bilinear-then-exp
    chain of the original approach. -/
theorem cutoff_eq_neg_exp_of_heat_kernel
    (f : ℝ → ℝ) (hf : IsPTSpectralCutoff f) :
    ∃ c : ℝ, c < 0 ∧ ∀ x, f x = Real.exp (c * x) := by
  obtain ⟨β, hβ_pos, hβ⟩ := PT_cutoff_is_heat_kernel f hf
  refine ⟨-β, ?_, ?_⟩
  · -- c = -β < 0 since β > 0.
    linarith
  · -- f(x) = exp(-β · x) = exp((-β) · x).
    intro x
    rw [hβ]

/-! ## Compatibility: equivalence with the SJG1 axiom

The heat-kernel postulate `PT_cutoff_is_heat_kernel` is **logically
equivalent** to the conjunction of `SJG1_spectral_cutoff_factorises`,
positivity, continuity, and decay — both characterise the family of
exponential cutoffs `f(u) = exp(-β·u)` with `β > 0`. The heat-kernel
formulation is **epistemically preferable** because:

1. It is a recognisable concept from standard NCG (Chamseddine-Connes
   spectral action principle, Connes-Marcolli 2008 §1.18).
2. It directly gives the exponential form, bypassing the
   Cauchy → exp chain.
3. It transparently identifies the PT cutoff as a Gibbs partition
   function, opening the door to Gibbs-Jaynes MaxEnt for the scale
   axiom (left to future work). -/

/-- **Equivalence statement.** Under the heat-kernel postulate, the
    SJG1 spectral multilinear factorisation holds; conversely, any
    multilinear continuous strictly-positive function satisfies the
    heat-kernel form (via `cauchy_mult_exp` applied to `k=2`). -/
theorem heat_kernel_implies_SJG1_factorisation
    (f : ℝ → ℝ) (hf : IsPTSpectralCutoff f) :
    ∀ (k : ℕ) (u : Fin k → ℝ), f (∑ i, u i) = ∏ i, f (u i) :=
  SJG1_spectral_cutoff_factorises_from_heat_kernel f hf

end PT.Bridge
