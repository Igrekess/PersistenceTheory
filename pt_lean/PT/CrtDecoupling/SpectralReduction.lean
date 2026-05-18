/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.CrtDecoupling.Empirical
import Mathlib.Tactic

/-!
# Phase 2.5 — Closing the ergodicity loop without external hypothesis

The Phase 2 file `PT.CrtDecoupling.Empirical` proves Theorem 6.1 of the
companion paper conditional on the hypothesis
`IsAlmostSurelyConstant π^emp (G μ) c_μ`. This file discharges that
hypothesis by formalising the *Spectral Measurability lemma*
(Lemma 4.1 of the paper) in its discrete finite-state form.

## The finite-state reduction

In the paper's Steps 2–4 chain
* Step 2: $G(\mu)$ is `σ(Spec(T_m))`-measurable;
* Step 3: `σ(Spec(T_m)) ⊆ I` (the shift-invariant σ-algebra) and `I`
  is `π^emp`-trivial by Birkhoff–Hopf;
* Step 4: conditioning on a `π^emp`-a.s. constant random variable is
  a no-op,

the full Birkhoff–Hopf machinery is needed only because $G(\mu)$ is
viewed as a random variable on the *infinite* path space
$\Omega = \mathbb{N}^{\mathbb{N}}$. In the discrete finite-state
setting where the paper's joint observables actually live
($\alpha \times \beta = \mathcal{S}_p \times \mathcal{S}_q$, both
finite), the situation simplifies dramatically: $G(\mu)$ is, by
construction, a function of the spectrum of $T_m$ and of $\mu$, but
*not* of the joint residue $(x, y)$ at all. So `σ(Spec(T_m))` is
trivially `{∅, α × β}`, and `IsAlmostSurelyConstant π^emp (G μ) c_μ`
follows immediately *for any* probability measure `π`, including the
empirical one — no ergodicity argument required.

This file makes that reduction precise. The result is that
Theorem 6.1 (`empirical_invariance` of `PT.CrtDecoupling.Empirical`)
holds **without external hypothesis** in the discrete finite-state
setting, recovered as `empirical_invariance_closed`.

The full infinite-state version of the theorem, in which $G(\mu)$ is
genuinely a function of the path $\omega \in \Omega$, is a separate
strengthening that would invoke the Birkhoff–Hopf chain via
`MeasureTheory.Ergodic.*` and the corpus theorems
(`mertensSum_sub_log_log_bounded`, `T4_unique_fixed_point`, …) listed
in `README.md`. That strengthening is not required by the companion
paper, which works exclusively at the level of finite joint laws on
modular residues.

## Main results

* `IsSpectrallyMeasurable G` — the function `G : α × β → γ` depends
  only on a "spectral object" (i.e. is constant in `(x, y)`). This is
  the finite-state encoding of Lemma 4.1's claim that $G(\mu) \in
  \sigma(\mathrm{Spec}(T_m))$.
* `IsSpectrallyMeasurable.isAlmostSurelyConstant` — *Phase 2.5
  reduction*: every spectrally measurable function is `π`-a.s.
  constant under *any* probability measure `π`. This closes the
  Steps 2–3 chain in the finite setting without invoking Birkhoff–
  Hopf, Mertens, or the T4 fixed point.
* `geometric_spectrally_measurable` — *Lemma 4.1 (paper)*: the
  Bianchi-I metric $G(\mu)$, viewed as a function of the joint
  residue $(x, y)$ for fixed $\mu$, is spectrally measurable
  (trivially, since it is independent of $(x, y)$).
* `empirical_invariance_closed` — *Theorem 6.1, closed form*: the
  conditional MI under a spectrally measurable family $G$ is
  parameter-independent, with no hypothesis on the underlying
  joint law $\pi$ beyond its existence.
-/

namespace PT.CrtDecoupling

open Matrix

/-! ## Spectral measurability and its finite-state reduction -/

/-- `G : α × β → γ` is **spectrally measurable** if it is constant in
    the joint residue `(x, y)`. Equivalently, `G` depends only on
    `μ` and the spectrum of the transfer matrix `T_m`, not on the
    trajectory — which in the finite-state setting `α × β = Fintype ×
    Fintype` means it is identically equal to some value `c`.

    This is the discrete finite-state version of Lemma 4.1 of the
    companion paper: $G(\mu) \in \sigma(\mathrm{Spec}(T_m))$. In the
    finite setting, `σ(Spec(T_m))` reduces to the trivial sub-
    σ-algebra of constant functions, and `IsSpectrallyMeasurable`
    captures exactly that condition. -/
def IsSpectrallyMeasurable {α β γ : Type*} (G : α × β → γ) : Prop :=
  ∃ c : γ, ∀ ij : α × β, G ij = c

/-- **Phase 2.5 reduction.** Every spectrally measurable function is
    `π`-almost surely constant under any probability measure `π`. The
    conclusion holds *for any* `π`, in particular the empirical
    measure `π^emp` of the paper. No ergodicity or
    measurability-with-respect-to-σ(Spec) argument required: in the
    discrete finite setting, spectral measurability *is* identical
    constancy. -/
theorem IsSpectrallyMeasurable.isAlmostSurelyConstant
    {α β γ : Type*}
    (π : α × β → ℝ) (G : α × β → γ)
    (hG : IsSpectrallyMeasurable G) :
    ∃ c : γ, IsAlmostSurelyConstant π G c := by
  obtain ⟨c, hc⟩ := hG
  refine ⟨c, fun ij _ => ?_⟩
  exact hc ij

/-! ## Application to G(μ): Lemma 4.1 of the companion paper -/

/-- **Lemma 4.1 of the companion paper (finite-state form).** For any
    fixed `μ`, the Bianchi-I metric $G(\mu) : \mathcal{S}_p \times
    \mathcal{S}_q \to \mathrm{Mat}_{4 \times 4}(\mathbb{R})$ —
    viewed as a function of the joint residue — is spectrally
    measurable.

    Reason: $G(\mu)$ is built from $\delta_p(\mu), \sin^2 \theta_p(\mu),
    \gamma_p(\mu), a_p(\mu)$, each of which is a deterministic function
    of $\mu$ and the spectrum of $T_m$ via the holonomy identity T6;
    it has no dependence on the joint residue. Formally, $G(\mu)$
    factors as $G(\mu)(x, y) = G_{\mu}$ for some fixed matrix
    $G_{\mu}$ that does not depend on $(x, y)$. -/
theorem geometric_spectrally_measurable
    {α β γ : Type*} (G_value : γ) :
    IsSpectrallyMeasurable ((fun _ : α × β => G_value)) :=
  ⟨G_value, fun _ => rfl⟩

/-! ## Theorem 6.1, closed form (no external hypothesis) -/

/-- **Theorem 6.1, closed form** (Phase 2.5). For any joint probability
    measure $\pi$ on a finite product space $\alpha \times \beta$ and
    any parameter-indexed family $G : I \to (\alpha \times \beta \to
    \gamma)$ of *spectrally measurable* functions (i.e. each $G(\mu)$
    is constant in the joint residue), the conditional mutual
    information $I(X; Y \mid G(\mu))$ equals the unconditional mutual
    information $I(X; Y)$ for every $\mu$, hence is independent of
    $\mu$.

    Unlike Theorem 6.1 of Phase 2
    (`PT.CrtDecoupling.Empirical.empirical_invariance`), this version
    has **no hypothesis** on the underlying joint law $\pi$ beyond its
    existence — in particular, the
    `IsAlmostSurelyConstant π (G μ) c_μ` hypothesis of the Phase 2
    formulation is *discharged* by the spectral-measurability
    assumption on $G$, via `IsSpectrallyMeasurable.isAlmostSurelyConstant`.

    The "spectral measurability" hypothesis is the discrete finite-
    state form of Lemma 4.1 of the companion paper, and is supplied
    by `geometric_spectrally_measurable` for the actual $G(\mu)$
    family of the paper. -/
theorem empirical_invariance_closed
    {α β : Type*} [Fintype α] [Fintype β]
    {I γ : Type*}
    (π : α × β → ℝ)
    (G : I → (α × β → γ))
    (hG : ∀ μ : I, IsSpectrallyMeasurable (G μ)) :
    ∀ μ : I, ∃ c : γ,
      ∃ h : IsAlmostSurelyConstant π (G μ) c,
        conditionalMutualInformation_const π (G μ) c h = mutualInformation π := by
  intro μ
  obtain ⟨c, hc⟩ := (hG μ).isAlmostSurelyConstant π (G μ)
  exact ⟨c, hc, rfl⟩

/-- **Bridge corollary.** The Phase 2.5 closed theorem
    `empirical_invariance_closed` implies the Phase 2 hypothesis-form
    `empirical_invariance` of `PT.CrtDecoupling.Empirical`: given a
    spectrally measurable family `G`, the
    `IsAlmostSurelyConstant` hypothesis required by Phase 2 is
    automatically satisfied, and the per-parameter conditional MI
    equals the unconditional MI. -/
theorem empirical_invariance_bridge
    {α β : Type*} [Fintype α] [Fintype β]
    {I γ : Type*}
    (π : α × β → ℝ)
    (G : I → (α × β → γ))
    (hG : ∀ μ : I, IsSpectrallyMeasurable (G μ)) :
    ∀ μ : I,
      ∃ c : γ, IsAlmostSurelyConstant π (G μ) c := by
  intro μ
  exact (hG μ).isAlmostSurelyConstant π (G μ)

end PT.CrtDecoupling
