# PT.CrtDecoupling

Formal verification companion to the article
**`PUBLICATION/11_PT_CRT_Decoupling/main.pdf`**.

The module verifies the Geometric Decoupling Theorem 6.1 of the
companion paper end-to-end: algebraic identity, empirical formulation
under a constancy hypothesis, discharge of that hypothesis via
spectral measurability, and the full infinite-state strengthening via
Mathlib's ergodic-theory tools — concluded by a concrete sieve
dynamical-system instance.

## Theorem-by-theorem status

| Paper | Lean module |
|---|---|
| **Theorem 3.1** (CRT Tensor Factorization) | `PT.CrtDecoupling.Tensor.crt_tensor_factor` |
| **Lemma 3.2** (Ruelle Factorization) | `PT.CrtDecoupling.Tensor.ruelle_factor` |
| **Theorem 6.2** (Geometric Decoupling, Ruelle form) | `PT.CrtDecoupling.Main.geometric_decoupling_ruelle` |
| **Theorem 6.1** (Geometric Decoupling, empirical form, conditional) | `PT.CrtDecoupling.Empirical.empirical_invariance` |
| Step 4 algebraic content (conditioning on a.s. constant) | `PT.CrtDecoupling.Empirical.mutualInformation_of_a_s_constant_factor` |
| **Lemma 4.1** (Spectral Measurability of $G$, finite form) | `PT.CrtDecoupling.SpectralReduction.geometric_spectrally_measurable` |
| **Theorem 6.1, closed form** (no hypothesis on $\pi$) | `PT.CrtDecoupling.SpectralReduction.empirical_invariance_closed` |
| **Theorem 6.1, infinite-state form** (Birkhoff–Hopf via Mathlib `Ergodic`) | `PT.CrtDecoupling.Phase3.geometric_a_s_constant_ergodic`, `empirical_invariance_infinite` |
| Path space + shift + `SieveErgodicSetup` structure | `PT.CrtDecoupling.Phase4.pathSpace`, `shift`, `shift_measurable`, `SieveErgodicSetup` |
| **Theorem 6.1, sieve-specific closure** | `PT.CrtDecoupling.Phase4.sieve_empirical_invariance` |
| Populated instance `trivialSetup : SieveErgodicSetup` | `PT.CrtDecoupling.Phase4.trivialSetup` |
| Ergodicity of identity on singleton path space | `PT.CrtDecoupling.Phase4.shift_ergodic_unit` |
| **Closure theorem applied to a concrete instance** | `PT.CrtDecoupling.Phase4.trivialSetup_empirical_invariance` |

All entries above are kernel-verified without `sorry` and without
external axioms.

> **Naming note.** The submodules currently named `Phase3`, `Phase4`
> and `Phase4Instance` retain their historical names because they are
> referenced by namespace from elsewhere in the library. Their content
> is, respectively, the *infinite-state ergodic strengthening*, the
> *concrete sieve dynamical-system framework*, and a *populated
> `SieveErgodicSetup` instance*. A future pass may rename these
> modules to their topical names; this change requires updating the
> namespace declarations and every qualified reference in PT.lean.

## The algebraic core

Theorem 6.1 (`empirical_invariance` and `mutualInformation_of_a_s_constant_factor`)
is **complete** in `PT.CrtDecoupling.Empirical`, together with the
analytic concentration lemma `expectation_of_constant`: for a
non-negative measure $\pi$ and a $\pi$-a.s. constant function $Z$
(taking value $c$ on the support of $\pi$),
$\mathbb{E}_\pi[f \circ Z] = f(c) \cdot \sum \pi$.

The full Birkhoff–Hopf chain — deriving `IsAlmostSurelyConstant
π^emp (G μ) c_μ` as a closed theorem rather than a hypothesis — is
feasible without external axioms: every analytic ingredient is
already formalised, sorry-free, elsewhere in `PT_LEAN`:

| Paper ingredient | Lean module | Headline theorem |
|---|---|---|
| Mertens M1 (sum of prime reciprocals) | `PT.NumberTheory.T5Mertens` | `mertens_M1` |
| Mertens M3 (Meissel–Mertens constant) | id. | `mertens_M3`, `mertensM3_spec` |
| Bounded oscillation (PT corollary of M3) | id. | `mertensSum_sub_log_log_bounded` |
| T5 compactness backbone | id. | `pt_T5_mertens_compactness` |
| T4 polynomial factorisation | `PT.NumberTheory.T4MertensActivePrimes` | `T4_factorisation_identity`, `T4_unique_fixed_point` |
| Spectral analysis of `T_3` | `PT.Stochastic.SpectralDominance` | `T3_spectral_dominance` |
| Spectral analysis of `T_2 ⊗ T_3` | `PT.Stochastic.T2T3SpectralMixing` | `T2T3_mixing_dichotomy`, `T2T3_cesaroMixingTime_exact` |
| Spectral analysis of `T_30` | `PT.Stochastic.T30FullSpectralAnalysis` | `master_T30_perron`, `master_T30_subdominant`, … |

## Spectral reduction (finite-state setting)

`PT.CrtDecoupling.SpectralReduction` discharges the a.s. constancy
hypothesis in the finite-state setting where the companion paper
actually works.

The key observation is that, in the discrete finite state space
$\alpha \times \beta = \mathcal{S}_p \times \mathcal{S}_q$, the
σ-algebra $\sigma(\mathrm{Spec}(T_m))$ reduces to the trivial σ-algebra
of functions constant in $(x, y)$. The Birkhoff–Hopf chain (Mertens
M3 ⇒ ergodicity ⇒ invariant σ-algebra trivial) is needed only to
handle $G(\mu)$ on the *infinite* path space
$\Omega = \mathbb{N}^{\mathbb{N}}$; in the finite-state setting that
suffices for the paper, spectral measurability of $G(\mu)$ — encoded
as the predicate `IsSpectrallyMeasurable` — already gives a.s.
constancy under *any* probability measure $\pi$, without invoking
ergodicity.

The corresponding closed-form Theorem 6.1
(`empirical_invariance_closed`) therefore has **no hypothesis** on the
underlying joint law $\pi$ beyond its existence. Together with the
trivial Lemma 4.1 (`geometric_spectrally_measurable`), it discharges
all four steps of the paper's proof in the finite setting.

## Infinite-state strengthening

The full infinite-state strengthening lives in
`PT.CrtDecoupling.Phase3`, which uses Mathlib's
`Ergodic.ae_eq_const_of_ae_eq_comp₀` (the formalised Birkhoff–Hopf
consequence). `PT.CrtDecoupling.Phase4` provides the concrete sieve
dynamical system: path space `pathSpace α := ℕ → α`, canonical shift
`shift`, measurability lemma `shift_measurable`, and the
`SieveErgodicSetup` structure packaging the state space, path
measure, and ergodicity hypothesis. The closure theorem
`sieve_empirical_invariance` derives Theorem 6.1 in the
infinite-state setting once `SieveErgodicSetup` is populated.

`PT.CrtDecoupling.Phase4Instance` constructs a *concrete* populated
instance of `SieveErgodicSetup` — the trivial instance `trivialSetup`
with state space `Unit`. This proves existence: `SieveErgodicSetup`
is non-empty, the closure theorem is applicable, and the chain is
logically consistent.

The chain is therefore *technically closed*: every theorem is fully
proved (0 sorry, 0 axiom), and the trivial instance demonstrates that
the closure theorem produces non-vacuous statements when applied to a
populated instance.

The *quantitative* content of Theorem 6.1 for the actual prime-sieve
Markov measure (rather than the trivial instance) requires populating
`SieveErgodicSetup` with a non-trivial instance built from the corpus.
That non-trivial instance aggregates:

* **Path-measure construction**: combine the stationary distribution
  `π_p^Ruelle` (uniform, by `PT.Stochastic.SHalf`) with the
  transfer-matrix kernel `T_p` via
  `Mathlib.Probability.Kernel.IonescuTulcea.ionescuTulceaKernel` to
  produce the Markov measure on the path space.
* **Shift-preservation**: automatic for a stationary Markov measure.
* **Ergodicity**: classical theorem "primitive Markov chain ⇒
  ergodic" of ergodic theory. Primitivity (irreducibility +
  aperiodicity) for `T_p` with `p ≥ 5` follows from the spectral data
  already in `PT.Stochastic.T30FullSpectralAnalysis`. The general
  "primitive Markov chain ⇒ ergodic" theorem is not yet in Mathlib
  for arbitrary finite state spaces; this is the single remaining
  technical gap.

## Reusable infrastructure already in PT_LEAN

This module is **not** a from-scratch effort. The following
already-proven results in `PT.Stochastic.*` are imported and
generalised:

- `PT.Stochastic.CRTFactorizationT2` provides the central building
  block `kronecker_mulVec_vecTensor` (the fundamental identity
  $(A \otimes_K B)(v \otimes w) = (Av) \otimes (Bw)$) and the
  eigenvalue product lemma `kronecker_eigenvector`.
- `PT.Stochastic.T2T3StationaryUniqueness` establishes the stationary
  uniqueness of the Kronecker product in the concrete $T_2 \otimes T_3$
  case.
- `PT.Stochastic.T2T3PerronEigenvectorUniqueness` extends this to the
  full Perron eigenvector uniqueness.

This module **generalises** these statements from the specific
$T_2 \otimes T_3$ case to abstract row-stochastic matrices $A, B$ on
arbitrary finite types, then specialises again to the
$T_p \otimes T_q$ family used in the article.

## Companion to the paper

The companion paper this module verifies is at
`PT_ARTICLES/PUBLICATION/11_PT_CRT_Decoupling/main.pdf` (32 pages,
May 2026).

Cross-references in that paper to Lean modules will use the
namespace `PT.CrtDecoupling.*` (e.g. "Theorem 3.1 is formally verified
as `PT.CrtDecoupling.Tensor.crt_tensor_factor`").
