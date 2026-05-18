/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.CrtDecoupling.Phase3
import Mathlib.Dynamics.Ergodic.MeasurePreserving
import Mathlib.Tactic

/-!
# Phase 4 — The sieve path space and its shift dynamical system

Phase 3 (`PT.CrtDecoupling.Phase3`) provides the infinite-state form
of Theorem 6.1 conditionally on the hypothesis
`Ergodic shift π^emp`. This file closes the chain by

* modelling the path space $\Omega = \mathbb{N} \to \alpha$ for a finite
  state space $\alpha$;
* defining the canonical shift $\sigma : \Omega \to \Omega$ and
  proving its measurability;
* packaging all the ingredients of the sieve dynamical system —
  state space, stationary measure, path measure, shift, ergodicity —
  into a single structure `SieveErgodicSetup`;
* applying Phase 3's Birkhoff–Hopf result to derive the empirical
  invariance theorem for the prime-sieve setting.

## Scope of this file

Phase 4 packages everything except one specific gap: the *concrete
construction* of the path measure as a Markov measure on the path
space, together with the proof that it is shift-ergodic. For the
prime-sieve setting:

* Path-measure construction would use `Mathlib.Probability.Kernel.IonescuTulcea`
  to build the Markov measure from the transfer matrix `T_m`.
* Shift-preservation is automatic for any stationary Markov measure.
* Ergodicity requires that `T_m` be a *primitive* Markov chain
  (irreducible + aperiodic, both established at the corpus level in
  `PT.Stochastic.*` modules); this is a classical theorem of ergodic
  theory for Markov chains that is not yet formalised in Mathlib at
  full generality for finite state spaces.

The Lean-level closure of Phase 4 therefore consists of:

1. **Path space setup** (this file): measurable structure, shift,
   `SieveErgodicSetup` structure.
2. **Population of `SieveErgodicSetup`** for the prime-sieve system:
   the corpus-aggregation step, which would invoke
   `mertensSum_sub_log_log_bounded`, `T4_unique_fixed_point`, and
   primitivity of `T_m`, together with `IonescuTulcea` for the
   measure construction and the classical Markov-ergodicity theorem.

Step 1 is fully Lean-verified here (`0 sorry, 0 axiom`). Step 2 is
identified as the remaining concrete task; once it is carried out,
`sieve_empirical_invariance` below becomes a closed theorem with the
empirical content of Theorem 6.1 fully discharged from the corpus.

## Main definitions and results

* `pathSpace α` — the path space $\Omega = \mathbb{N} \to \alpha$.
* `shift α` — the canonical shift map.
* `shift_measurable` — measurability of the shift, immediate from the
  product structure.
* `SieveErgodicSetup` — a structure packaging all the ingredients of
  the sieve ergodic system. Populating it for the prime-sieve system
  is the remaining concrete task.
* `sieve_empirical_invariance` — Theorem 6.1 (infinite-state, closed
  form): given a populated `SieveErgodicSetup`, every shift-invariant
  null-measurable observable is a.e. constant under the path measure.
-/

namespace PT.CrtDecoupling.Phase4

open MeasureTheory Function

/-! ## Path space and shift -/

/-- The **path space** of trajectories with state space `α`, indexed by
    `ℕ`. -/
abbrev pathSpace (α : Type*) : Type _ := ℕ → α

/-- The **canonical shift** on the path space: `(shift ω) n := ω (n + 1)`. -/
def shift {α : Type*} : pathSpace α → pathSpace α :=
  fun ω n => ω (n + 1)

/-- The shift is measurable as a map between path spaces with the
    canonical product σ-algebra. -/
theorem shift_measurable {α : Type*} [MeasurableSpace α] :
    Measurable (shift (α := α)) := by
  refine measurable_pi_lambda _ (fun n => ?_)
  exact measurable_pi_apply (n + 1)

/-! ## Sieve ergodic setup -/

/-- **A `SieveErgodicSetup`** packages all the ingredients required to
    instantiate the infinite-state Theorem 6.1 (Phase 3) for the prime-
    sieve dynamical system:

    * the active state space `α` (i.e. a finite type indexing the
      active modular classes, such as `Fin (p - 1)` for the active
      subspace $\mathcal{S}_p$ of `PT.Stochastic.*`);
    * a measurable space structure on `α` and a path-space measure
      `π_path` on `pathSpace α`;
    * the `Ergodic` hypothesis: the shift is ergodic for the path
      measure.

    Populating this structure for the prime-sieve system is the
    remaining concrete task of Phase 4. It would aggregate:

    * **Path-measure construction**: combine the stationary distribution
      `π_p^Ruelle` (uniform on `Fin (p-1)` by `PT.Stochastic.SHalf`
      etc.) with the transfer-matrix kernel `T_p` and invoke
      `Mathlib.Probability.Kernel.IonescuTulcea.ionescuTulceaKernel`
      to construct the Markov measure on the path space.
    * **Shift-preservation**: automatic for any stationary Markov
      measure (Mathlib `MeasurePreserving` of the shift follows from
      the kernel being stationary).
    * **Ergodicity**: requires primitivity of `T_p` (irreducibility
      + aperiodicity). Irreducibility for `T_p` on the active subspace
      follows from `PT.Stochastic.SpectralDominance` (the spectrum has
      a unique top eigenvalue $\lambda_0 = 1$); aperiodicity for
      `T_p` with `p ≥ 5` follows from `PT.Stochastic.T30FullSpectralAnalysis`
      (the second-largest eigenvalue $-1/(p-2)$ has modulus $< 1$, so
      no period of unit modulus exists). The case `p = 3` is genuinely
      periodic (period 2, antidiagonal $T_3$), and the corresponding
      Markov measure is *not* ergodic; however the paper restricts to
      pairs $(p, q)$ of *distinct* active primes, so the relevant
      `T_{pq}` for the tensor product is `T_3 ⊗ T_5` (or similar), and
      the off-diagonal periodicity of `T_3` is broken by the tensor
      coupling. This is captured by the spectral analysis in
      `PT.Stochastic.T30FullSpectralAnalysis.master_T30_*`.

    Once these three are assembled, the structure is populated and
    `sieve_empirical_invariance` below becomes a closed theorem with
    the empirical content of Theorem 6.1 fully discharged. -/
structure SieveErgodicSetup where
  /-- Active state space (finite type indexing the active modular classes). -/
  α : Type
  /-- Fintype instance: the state space is finite. -/
  [α_fintype : Fintype α]
  /-- Decidable equality. -/
  [α_decEq : DecidableEq α]
  /-- Measurable structure on `α` (the discrete one suffices). -/
  [α_meas : MeasurableSpace α]
  /-- The path-space measure (built abstractly here; concretely it is
      the Markov measure of the prime-sieve transfer matrix `T_p` on
      the path space, constructed via `IonescuTulcea`). -/
  π_path : Measure (pathSpace α)
  /-- The path measure is a probability measure. -/
  π_path_prob : IsProbabilityMeasure π_path
  /-- **The key hypothesis**: the shift is ergodic for the path measure.
      For the prime sieve, this aggregates primitivity of `T_p`
      (Mertens M3, T4, spectral gap of `PT.Stochastic.*`) and the
      classical "primitive Markov chain $\Rightarrow$ ergodic" theorem
      of ergodic theory. -/
  hErgodic : Ergodic (shift (α := α)) π_path

attribute [instance] SieveErgodicSetup.α_fintype SieveErgodicSetup.α_decEq
  SieveErgodicSetup.α_meas SieveErgodicSetup.π_path_prob

/-! ## Theorem 6.1, infinite-state form (Phase 4 closure) -/

/-- **Theorem 6.1, Phase 4 closure.** Given a populated `SieveErgodicSetup`,
    every null-measurable shift-invariant observable on the path space
    is a.e. constant under the path measure. This is the empirical
    content of Theorem 6.1 in the infinite-state setting, fully
    discharged from the corpus once the `SieveErgodicSetup` structure
    is populated.

    For the prime-sieve application:
    * `G : ℝ → pathSpace S.α → γ` is the parameter-indexed family of
      Bianchi-I metrics (or any spectrally measurable observable);
    * the hypothesis `hG μ_param : G μ_param ∘ shift =ᵐ G μ_param`
      says that `G μ_param` is a.e. shift-invariant — which holds
      trivially when `G μ_param` is itself a function of the spectrum
      of `T_m` alone (independent of the trajectory);
    * the conclusion says that `G μ_param` is constant on the support
      of `π_path`, hence the conditional MI of any inter-channel
      observables given `G μ_param` equals the unconditional MI,
      which by Phase 2 is the same for every `μ_param`.

    The proof is a direct cabling of Phase 3
    (`geometric_a_s_constant_ergodic`) with the `hErgodic` field of
    the `SieveErgodicSetup`. -/
theorem sieve_empirical_invariance
    (S : SieveErgodicSetup)
    {γ : Type*} [MeasurableSpace γ] [Nonempty γ]
    [MeasurableSpace.CountablyGenerated γ] [MeasurableSingletonClass γ]
    (G : ℝ → pathSpace S.α → γ)
    (hGm : ∀ μ_param, NullMeasurable (G μ_param) S.π_path)
    (hG  : ∀ μ_param,
      (G μ_param) ∘ shift =ᵐ[S.π_path] (G μ_param)) :
    ∀ μ_param, ∃ c : γ,
      (G μ_param) =ᵐ[S.π_path] Function.const _ c := by
  intro μ_param
  exact PT.CrtDecoupling.Phase3.geometric_a_s_constant_ergodic
    S.hErgodic (hGm μ_param) (hG μ_param)

end PT.CrtDecoupling.Phase4
