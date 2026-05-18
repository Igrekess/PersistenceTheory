/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Sieve.T3Antidiagonal
import PT.Stochastic.SHalf
import PT.FixedPoint.T7MuStar
import PT.FixedPoint.T7GlobalUniqueness
import PT.Holonomy.ActivePrimeCriterion
import PT.Holonomy.CyclicPhaseIdentity
import PT.Holonomy.CouplingReconstruction
import PT.Holonomy.CouplingReconstructionBounds

/-!
# The PT Derivation Chain — Cross-Axes Synthesis Theorem

This file does **no new computation**. Its purpose is to *expose* the full
Persistence Theory derivation as a **single conjunctive theorem** that chains
the eight foundational steps already proven in their respective modules:

```
  Step 1.  T_3 = !![0,1;1,0]                              (T1 sieve antidiagonal)
            │
            │  (unique trace-0 doubly-stochastic 2×2)
            ▼
  Step 2.  T_3 has a unique stationary distribution π = (1/2, 1/2)
            │
            │  (forced by πᵀ T_3 = πᵀ)
            ▼
  Step 3.  s := π(0) = π(1) = 1/2                         (foundational axiom-as-theorem)
            │
            │  (anomalous-dimension fixed-point map F_pers)
            ▼
  Step 4.  μ⋆ = 15 is the unique non-degenerate fixed point of F_pers on μ ≥ 8
            │
            │  (γ_p > 1/2 criterion at q_+ = 13/15)
            ▼
  Step 5.  persActive = {3, 5, 7}                         (γ_p > s for p ∈ {3,5,7};
            │                                               γ_p < s for p ∈ {11,13})
            │
            │  (Pythagorean identity, cos θ_p = 1 - δ_p)
            ▼
  Step 6.  sin² θ_p = δ_p (2 - δ_p)                       (cyclic-phase identity)
            │
            │  (product over active primes)
            ▼
  Step 7.  α_bare = ∏_{p ∈ {3,5,7}} sin² θ_p              (coupling reconstruction)
            │
            │  (exact rational evaluation at q_+ = 13/15)
            ▼
  Step 8.  α_bare ∈ (7335/10⁶, 7340/10⁶)                  (matches published 7.338 × 10⁻³;
                                                            α_EM via Koide-running, outside
                                                            this kernel)
```

The headline `PT_full_derivation_chain` reexports all eight steps as a single
proposition, providing a self-contained *certificate* that the PT cascade
"`T_3 antidiagonal` ⇒ ... ⇒ `α_bare` bracket" is fully formalised in Lean
with no `sorry` and no axiom beyond Mathlib classical foundations.

Two intermediate synthesis theorems (`PT_sieve_to_s_chain` and
`PT_muStar_to_alphaBare_chain`) bundle adjacent steps for readability and for
downstream re-use.

## Reference

* Monograph Chapters 1, 6, 7, 9 — the canonical PT derivation arc.
* `AUDIT_FORMALISABLE.md` — rows for T1, SHalf, T7, ActivePrimeCriterion,
  CyclicPhaseIdentity, CouplingReconstruction, CouplingReconstructionBounds.
-/

namespace PT.Bridge

open Matrix PT.Sieve PT.Stochastic PT.FixedPoint PT.Holonomy

/-! ## Intermediate synthesis: sieve → symmetry -/

/-- **Sieve-to-symmetry chain (Steps 1 → 2 → 3).**
    Reexports the chain "`T₃` is the antidiagonal permutation
    ⇒ `T₃` has a unique stationary distribution ⇒ that distribution is
    `(1/2, 1/2)` ⇒ `s = 1/2`" as a single conjunction. -/
theorem PT_sieve_to_s_chain :
    -- (1) T_3 is the antidiagonal permutation, square = identity:
    (T3 = !![0,1;1,0] ∧ T3 * T3 = (1 : Matrix (Fin 2) (Fin 2) ℝ))
    -- (2) Uniqueness of the stationary distribution:
    ∧ (∀ π : Fin 2 → ℝ, IsStationary π → π = piHalf)
    -- (3) The foundational identity `s = 1/2`:
    ∧ (PT.Stochastic.s = 1 / 2)
    -- (3') Every stationary distribution takes the value `s`:
    ∧ (∀ π : Fin 2 → ℝ, IsStationary π → ∀ i : Fin 2, π i = PT.Stochastic.s) := by
  refine ⟨⟨rfl, T3_sq_eq_one⟩, ?_, s_eq_one_half, ?_⟩
  · intro π h; exact T3_unique_stationary π h
  · intro π h i; exact s_is_T3_stationary_value π h i

/-! ## Intermediate synthesis: fixed point → coupling -/

/-- **Fixed-point-to-coupling chain (Steps 4 → 5 → 6 → 7 → 8).**
    Reexports the chain "`μ⋆ = 15` is the unique fixed point of `F_pers`
    ⇒ `persActive = {3, 5, 7}` (γ_p > 1/2 for p ∈ {3,5,7}, γ_p < 1/2 for
    p ∈ {11,13}) ⇒ `sin² θ_p = δ_p (2 - δ_p)` is the cyclic-phase
    identity ⇒ `α_bare = ∏ sin² θ_p` ⇒ `α_bare ∈ (7335/10⁶, 7340/10⁶)`"
    as a single conjunction. -/
theorem PT_muStar_to_alphaBare_chain :
    -- (4) μ⋆ = 15 = 3 + 5 + 7, and is the unique fixed point of F_pers on μ ≥ 8:
    (PT.FixedPoint.muStar = 15
       ∧ PT.FixedPoint.muStar = 3 + 5 + 7
       ∧ Fpers PT.FixedPoint.muStar = PT.FixedPoint.muStar
       ∧ ∀ μ : ℕ, 8 ≤ μ → Fpers μ = μ → μ = 15)
    -- (5) Active prime criterion: persActive = {3, 5, 7}
    ∧ (IsActive 3 ∧ IsActive 5 ∧ IsActive 7
       ∧ ¬ IsActive 11 ∧ ¬ IsActive 13)
    -- (6) Cyclic-phase identity (Pythagorean, scalar form)
    ∧ (∀ θ δ : ℝ, Real.cos θ = 1 - δ → Real.sin θ ^ 2 = δ * (2 - δ))
    -- (7) Coupling reconstruction: α_bare factorises into the product
    ∧ (alphaBareQ = sinSqQ 3 * sinSqQ 5 * sinSqQ 7)
    -- (7') Real-valued coupling reconstruction
    ∧ (∀ θ3 θ5 θ7 : ℝ,
        Real.cos θ3 = 1 - (deltaQ 3 : ℝ) →
        Real.cos θ5 = 1 - (deltaQ 5 : ℝ) →
        Real.cos θ7 = 1 - (deltaQ 7 : ℝ) →
        Real.sin θ3 ^ 2 * Real.sin θ5 ^ 2 * Real.sin θ7 ^ 2 = (alphaBareQ : ℝ))
    -- (8) Tight numerical bracket on the bare coupling
    ∧ (7335 / 1000000 < alphaBareQ ∧ alphaBareQ < 7340 / 1000000) := by
  refine ⟨?_, active_primes_3_5_7, ?_, alphaBareQ_eq_prod, ?_,
          alphaBareQ_bracket_tight⟩
  · exact ⟨muStar_eq, muStar_eq_three_plus_five_plus_seven,
           T7_muStar_isFixed, T7_muStar_unique_global⟩
  · intro θ δ h; exact sin_sq_of_cos_eq_one_sub θ δ h
  · intro θ3 θ5 θ7 h3 h5 h7
    exact coupling_reconstruction_real θ3 θ5 θ7 h3 h5 h7

/-! ## The full chain -/

/-- **`PT_full_derivation_chain` — the master synthesis theorem.**

    This is the single conjunction certifying the entire PT cascade

    ```
    T_3 antidiagonal (T1)
      → T_3 unique stationary = (1/2, 1/2) (SHalf)
      → s = 1/2 (SHalf)
      → μ⋆ = 15 (T7)
      → persActive = {3, 5, 7} (T7MuStar + ActivePrimeCriterion)
      → γ_p > 1/2 for active primes (ActivePrimeCriterion)
      → sin² θ_p = δ_p (2 - δ_p) (CyclicPhaseIdentity)
      → α_bare = ∏ sin² θ_p (CouplingReconstruction)
      → α_bare ∈ (7335/10⁶, 7340/10⁶) (CouplingReconstructionBounds)
    ```

    All eight links are previously formalised in the PT library; this theorem
    merely composes them in one statement. No new derivation; no `sorry`. -/
theorem PT_full_derivation_chain :
    -- Step 1: T_3 is the antidiagonal sieve transfer matrix
    (T3 = !![0,1;1,0] ∧ T3.trace = 0 ∧ T3.det = -1)
    -- Step 2: T_3 admits a unique stationary distribution = (1/2, 1/2)
    ∧ (IsStationary piHalf ∧ ∀ π, IsStationary π → π = piHalf)
    -- Step 3: s = 1/2 is forced
    ∧ (PT.Stochastic.s = 1 / 2)
    -- Step 4: μ⋆ = 15, global uniqueness on μ ≥ 8
    ∧ (PT.FixedPoint.muStar = 15
       ∧ Fpers PT.FixedPoint.muStar = PT.FixedPoint.muStar
       ∧ ∀ μ : ℕ, 8 ≤ μ → Fpers μ = μ → μ = 15)
    -- Step 5: persActive = {3, 5, 7}
    ∧ (persActiveAt15 = ({3, 5, 7} : Finset ℕ)
       ∧ PT.FixedPoint.muStar = (persActiveAt15.sum id : ℕ))
    -- Step 6: γ_p > s = 1/2 for active primes, γ_p < s for p ∈ {11,13}
    ∧ (gammaQ 3 > sPT ∧ gammaQ 5 > sPT ∧ gammaQ 7 > sPT
       ∧ gammaQ 11 < sPT ∧ gammaQ 13 < sPT)
    -- Step 7: Cyclic-phase identity sin² θ_p = δ_p (2 - δ_p)
    ∧ (∀ θ δ : ℝ, Real.cos θ = 1 - δ → Real.sin θ ^ 2 = δ * (2 - δ))
    -- Step 8: Coupling reconstruction α_bare = ∏_{p ∈ {3,5,7}} sin² θ_p
    ∧ (alphaBareQ = sinSqQ 3 * sinSqQ 5 * sinSqQ 7
       ∧ ∀ θ3 θ5 θ7 : ℝ,
           Real.cos θ3 = 1 - (deltaQ 3 : ℝ) →
           Real.cos θ5 = 1 - (deltaQ 5 : ℝ) →
           Real.cos θ7 = 1 - (deltaQ 7 : ℝ) →
           Real.sin θ3 ^ 2 * Real.sin θ5 ^ 2 * Real.sin θ7 ^ 2
             = (alphaBareQ : ℝ))
    -- Step 9: Tight numerical bracket α_bare ∈ (7335/10⁶, 7340/10⁶)
    ∧ (0 < alphaBareQ
       ∧ 7335 / 1000000 < alphaBareQ ∧ alphaBareQ < 7340 / 1000000) := by
  refine ⟨?step1, ?step2, s_eq_one_half, ?step4, ?step5, ?step6, ?step7,
          ?step8, ?step9⟩
  · exact ⟨rfl, T3_trace_zero, T3_det_neg_one⟩
  · exact ⟨piHalf_isStationary, fun π h => T3_unique_stationary π h⟩
  · exact ⟨muStar_eq, T7_muStar_isFixed, T7_muStar_unique_global⟩
  · exact ⟨rfl, muStar_eq_sum_persActive⟩
  · exact ⟨gamma_3_active, gamma_5_active, gamma_7_active,
           gamma_11_inactive, gamma_13_inactive⟩
  · intro θ δ h; exact sin_sq_of_cos_eq_one_sub θ δ h
  · refine ⟨alphaBareQ_eq_prod, ?_⟩
    intro θ3 θ5 θ7 h3 h5 h7
    exact coupling_reconstruction_real θ3 θ5 θ7 h3 h5 h7
  · exact ⟨alphaBareQ_pos,
           alphaBareQ_bracket_tight.1, alphaBareQ_bracket_tight.2⟩

/-! ## Narrative headline -/

/-- **Narrative headline.** A short, human-readable certificate that the
    full PT derivation chain holds. Equivalent to (but more concise than)
    `PT_full_derivation_chain`: it gives the "load-bearing" identities of
    each step. -/
theorem PT_cascade_headline :
    -- (s = 1/2) ∧ (μ⋆ = 3 + 5 + 7) ∧ (active primes are exactly {3,5,7})
    -- ∧ (α_bare = ∏ sin² θ_p ∈ tight bracket)
    PT.Stochastic.s = 1 / 2
    ∧ PT.FixedPoint.muStar = 3 + 5 + 7
    ∧ (IsActive 3 ∧ IsActive 5 ∧ IsActive 7
       ∧ ¬ IsActive 11 ∧ ¬ IsActive 13)
    ∧ alphaBareQ = sinSqQ 3 * sinSqQ 5 * sinSqQ 7
    ∧ 7335 / 1000000 < alphaBareQ ∧ alphaBareQ < 7340 / 1000000 := by
  refine ⟨s_eq_one_half, muStar_eq_three_plus_five_plus_seven,
          active_primes_3_5_7, alphaBareQ_eq_prod,
          alphaBareQ_bracket_tight.1, alphaBareQ_bracket_tight.2⟩

end PT.Bridge
