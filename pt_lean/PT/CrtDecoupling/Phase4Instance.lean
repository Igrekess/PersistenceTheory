/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.CrtDecoupling.Phase4
import Mathlib.Tactic

/-!
# Phase 4 — Existence of a populated `SieveErgodicSetup`

The Phase 4 file `PT.CrtDecoupling.Phase4` defines the
`SieveErgodicSetup` structure that packages all the ingredients of
the sieve dynamical system: state space, path measure, shift,
ergodicity. The conclusion of Phase 4
(`sieve_empirical_invariance`) holds for *any* populated instance of
this structure.

This file *closes the loop* by exhibiting a populated instance —
demonstrating that `SieveErgodicSetup` is non-empty and that the
final theorem `sieve_empirical_invariance` is *applicable*.

## The trivial instance

The instance we construct uses the trivial state space `α = Unit`.
On this state space:

* the path space `pathSpace Unit = ℕ → Unit` is a singleton type
  (`Subsingleton (ℕ → Unit)`), because every function to `Unit` is
  equal to the constant function;
* every map `f : (ℕ → Unit) → (ℕ → Unit)` is the identity on this
  singleton — in particular the shift is the identity;
* the only probability measure on `(ℕ → Unit)` is the Dirac measure
  at the unique element, and the identity map is trivially ergodic
  for this measure.

This is the **smallest** instance of `SieveErgodicSetup`: it shows
that the structure is non-empty, hence the Phase 4 closure theorem
`sieve_empirical_invariance` is not vacuous. In particular, the
chain of Phases 1 → 2 → 2.5 → 3 → 4 is **logically consistent**
and *applicable*.

## What this does and doesn't establish

This instance establishes that:

* `SieveErgodicSetup` is non-empty (existence).
* The Phase 4 closure theorem `sieve_empirical_invariance` produces a
  non-vacuous statement when applied to it.

It does *not* establish:

* The empirical content of Theorem 6.1 for the *actual* prime-sieve
  Markov measure on the path space of prime gap residues. That
  requires a different, non-trivial instance of `SieveErgodicSetup`
  built from
  - `Mathlib.Probability.Kernel.IonescuTulcea` (for the Markov measure
    construction from the transfer-matrix kernel `T_p`),
  - `PT.Stochastic.T30FullSpectralAnalysis.master_T30_*` (for the
    primitivity of `T_p`),
  - the classical theorem "primitive Markov chain ⇒ ergodic", which
    is not yet formalised in Mathlib at the level of finite state
    spaces (this is the single remaining technical gap).

The non-trivial instance is the natural continuation of this file
and is the single remaining concrete task to close the
formalisation entirely. The structural skeleton (this file) shows
that the work plan is well-formed.

## Main results

* `trivialSetup` — a populated `SieveErgodicSetup` with `α = Unit`.
* `trivialSetup_empirical_invariance` — applying
  `sieve_empirical_invariance` to `trivialSetup`: any shift-invariant
  null-measurable observable on the trivial path space is a.e.
  constant under the Dirac measure.
-/

namespace PT.CrtDecoupling.Phase4

open MeasureTheory Function

/-! ## The path space `ℕ → Unit` is a singleton -/

/-- The trivial path space `ℕ → Unit` has at most one element. -/
instance subsingleton_path_unit : Subsingleton (pathSpace Unit) :=
  ⟨fun _ _ => funext fun _ => rfl⟩

/-- The shift on `pathSpace Unit` is the identity (in fact, every
    self-map is the identity on a Subsingleton). -/
theorem shift_unit_eq_id : (shift : pathSpace Unit → pathSpace Unit) = id := by
  funext ω; exact Subsingleton.elim _ _

/-! ## The Dirac measure is shift-invariant and ergodic on `pathSpace Unit` -/

/-- The canonical Dirac probability measure on `pathSpace Unit`,
    concentrated at the unique element. -/
noncomputable def diracUnitPath : Measure (pathSpace Unit) :=
  Measure.dirac (fun _ => PUnit.unit)

instance diracUnitPath_isProbabilityMeasure :
    IsProbabilityMeasure diracUnitPath := by
  unfold diracUnitPath
  exact Measure.dirac.isProbabilityMeasure

/-- The shift on `pathSpace Unit` preserves the Dirac measure.
    Immediate from `shift_unit_eq_id` and `MeasurePreserving.id`. -/
theorem shift_measurePreserving_unit :
    MeasurePreserving (shift : pathSpace Unit → pathSpace Unit) diracUnitPath diracUnitPath := by
  rw [shift_unit_eq_id]
  exact MeasurePreserving.id _

/-- The shift on `pathSpace Unit` is pre-ergodic for the Dirac measure.
    On a singleton, every measurable set is either ∅ or the whole
    space (modulo the Dirac measure), and both are trivially
    invariant. -/
theorem shift_preErgodic_unit :
    PreErgodic (shift : pathSpace Unit → pathSpace Unit) diracUnitPath where
  aeconst_set {s} _ _ := by
    -- On a Subsingleton, every set is either empty or the universe.
    rw [Filter.eventuallyConst_set']
    rcases Set.eq_empty_or_nonempty s with hs | hs
    · -- s = ∅
      left; rw [hs]
    · -- s nonempty: since Subsingleton, s = univ
      right
      have h_univ : s = Set.univ := by
        ext ω
        simp only [Set.mem_univ, iff_true]
        obtain ⟨ω', hω'⟩ := hs
        exact Subsingleton.elim ω' ω ▸ hω'
      rw [h_univ]

/-- The shift on `pathSpace Unit` is ergodic for the Dirac measure. -/
theorem shift_ergodic_unit :
    Ergodic (shift : pathSpace Unit → pathSpace Unit) diracUnitPath :=
  { shift_measurePreserving_unit with
    aeconst_set := shift_preErgodic_unit.aeconst_set }

/-! ## The trivial `SieveErgodicSetup` instance -/

/-- **The trivial populated instance of `SieveErgodicSetup`.**

    State space `α = Unit`, path measure = Dirac at the unique
    element of `pathSpace Unit`, shift ergodicity = the trivial
    ergodicity of the identity on a singleton.

    This is the smallest populated instance; it demonstrates that
    `SieveErgodicSetup` is non-empty and that the Phase 4 closure
    theorem `sieve_empirical_invariance` is applicable to *some*
    actual structure. -/
noncomputable def trivialSetup : SieveErgodicSetup where
  α := Unit
  π_path := diracUnitPath
  π_path_prob := diracUnitPath_isProbabilityMeasure
  hErgodic := shift_ergodic_unit

/-! ## Applying the Phase 4 closure theorem -/

/-- **The chain closes: Phase 4 closure theorem applied to the trivial
    instance.** For any parameter-indexed family of shift-invariant
    null-measurable observables on the trivial path space, each
    observable is a.e. constant under the Dirac measure.

    This is the concrete consequence of populating `SieveErgodicSetup`
    with `trivialSetup`. The same theorem applies to any other
    populated instance — in particular to the (still-to-be-built)
    instance for the actual prime-sieve Markov measure. -/
theorem trivialSetup_empirical_invariance
    {γ : Type*} [MeasurableSpace γ] [Nonempty γ]
    [MeasurableSpace.CountablyGenerated γ] [MeasurableSingletonClass γ]
    (G : ℝ → pathSpace Unit → γ)
    (hGm : ∀ μ_param, NullMeasurable (G μ_param) diracUnitPath)
    (hG  : ∀ μ_param,
      (G μ_param) ∘ shift =ᵐ[diracUnitPath] (G μ_param)) :
    ∀ μ_param, ∃ c : γ,
      (G μ_param) =ᵐ[diracUnitPath] Function.const _ c :=
  sieve_empirical_invariance trivialSetup G hGm hG

end PT.CrtDecoupling.Phase4
