# PT Lean Formalisation

A Lean 4 + Mathlib package formalising the foundational theorems of
**Persistence Theory** (PT).

Twenty-two foundational theorems on the critical path

```
T1 → T3 → s = 1/2 → T2 → L0 → T7 → pt_T5_mertens_compactness → W7-1
```

are *kernel-verified without `sorry`*. Around 160 additional secondary
modules cover the surrounding ecosystem (conservation identities, GFT,
Bekenstein bound, cyclic phase, active-prime criterion, EML, CRT
decoupling, …). Across the whole tree, **186 of 187 modules are entirely
sorry-free** — a single module (`G3FisherUniqueness`) keeps one
documented `sorry` for a classical uniqueness result that the monograph
already marks as an external dependency (`\leanExternal`).

## Status

Last verified build: **2026-05-18**, Lean `v4.30.0-rc2`,
Mathlib `e12bcd0ff19fedc29dbd4dbbab360ea3a0f47a0a`.

- Foundational count (critical path T1 → T7 → W7-1, both directions):
  **30 named theorems** formally proved without `sorry`, including the
  extensions added on 2026-05-18: N2 uniqueness, T6b under the full C1-C4
  axiom system, L0 strict uniqueness, Bekenstein equality case, T2
  spectral bridge via Kronecker factorisation, W7-1 reverse direction,
  and active-prime monotonicity for all primes `p ≥ 11`.
- Full ecosystem: **187 modules / ~2 850 declarations** across the `PT/`
  sub-directories. **One** documented `sorry` remains, in
  `PT/Information/G3FisherUniqueness.lean` (Fisher metric uniqueness up
  to scale — a classical result, marked `\leanExternal` in the
  monograph).
- Build status: `lake build` completes successfully (3343 jobs on the
  newly-extended tree).

## Module catalog by topic

The library is organised by mathematical topic, mirroring the directory
layout under `PT/`. Each entry lists the theorem statement and the file
that proves it.

### Sieve (`PT/Sieve/`)

The discrete arithmetic core: the 6-rough predicate, the forbidden
mod-3 transition, and the antidiagonal transfer matrix.

| Theorem | File | Statement |
|---|---|---|
| **T1** | `PT/Sieve/T1ForbiddenTransitions.lean` | Consecutive 6-rough integers always switch residue mod 3. |
| **T3** | `PT/Sieve/T3Antidiagonal.lean` | The sieve-level transfer matrix is `!![0,1;1,0]`: doubly stochastic, trace 0, det −1, eigenvalues `{+1,−1}`, the **unique** trace-0 doubly-stochastic 2×2 matrix. |
| **T6a** | `PT/Sieve/T6aFieldStructure.lean` | Field structure of the cyclic phase carrier. |
| **T6b** | `PT/Sieve/T6bAxiomsC1C4.lean` | Axioms C1–C4 of the sieve operator. |
| **Bimodality** | `PT/Sieve/Bimodality.lean` | Bimodal structure of the sieve walk under T1 projection. |
| **N1** | `PT/Sieve/N1AtomicUniqueness.lean` | Atomic uniqueness of the prime cascade. |
| **N2** | `PT/Sieve/N2SelfCoherence.lean` | Self-coherence of the sieve. |

### Stochastic (`PT/Stochastic/`)

| Theorem | File | Statement |
|---|---|---|
| **`s = 1/2`** | `PT/Stochastic/SHalf.lean` | `T₃` admits the **unique** stationary distribution `(1/2, 1/2)`; the PT symmetry parameter `s = 1/2` is forced, not axiomatised. |

### Conservation (`PT/Conservation/`)

| Theorem | File | Statement |
|---|---|---|
| **T2** | `PT/Conservation/T2Alpha.lean` | The conservation exponent satisfies `α_cons = s² = 1/4`. |
| **Conservation ID** | `PT/Conservation/ConservationID.lean` | `∑ g_n = p_{N+1} − 2`. |

### Information (`PT/Information/`)

| Theorem | File | Statement |
|---|---|---|
| **L0** (algebraic kernel) | `PT/Information/L0MaxEntropy.lean` | Geometric distribution `(1−q) q^{k−1}` normalises to 1; algebraic identity `mean = μ ⇔ q = 1 − 1/μ`; PT specialisation `q_+ = 13/15` at `μ = 15/2`. |
| **L0** (headline) | same file, `L0_geometric_maximises_entropy` | Among all distributions of given mean, entropy is maximised by the geometric one. Proved via Gibbs inequality `1 − 1/y ≤ log y`. |
| **GFT identity** | `PT/Information/GFTIdentity.lean` | `log₂ m = D_KL + H`. |
| **Bekenstein bound** | `PT/Information/BekensteinBound.lean` | `D_KL ≤ log₂ m`. |
| **T6c (Shore–Johnson–Čencov)** | `PT/Information/T6cChencov.lean` | Shannon entropy + the four SJČ axioms (Shannon basics, max-uniform, additivity, Cauchy logarithmic equation, binary characterisation). Sorry-free. |

### Holonomy (`PT/Holonomy/`)

| Theorem | File | Statement |
|---|---|---|
| **Cyclic phase** | `PT/Holonomy/CyclicPhaseIdentity.lean` | `sin² θ_p = δ_p (2 − δ_p)`. |
| **Active-prime criterion** | `PT/Holonomy/ActivePrimeCriterion.lean` | `γ_p > 1/2 ⇔ p ∈ {3, 5, 7}`. |

### Fixed point (`PT/FixedPoint/`)

| Theorem | File | Statement |
|---|---|---|
| **T7** | `PT/FixedPoint/T7MuStar.lean` | `μ* = 15 = 3 + 5 + 7`; the unique fixed point of the persistence-active sum `F_pers` on `μ ∈ [8, 20]`. |

### Number theory (`PT/NumberTheory/`)

| Theorem | File | Statement |
|---|---|---|
| `mertensSum`, `mertensLog`, `mertensM3` (defs) | `PT/NumberTheory/T5Mertens.lean` | Standard definitions + boundary cases + non-negativity + `mertensM3_spec`. |
| **`pt_T5_mertens_compactness`** | same file | PT-needed: `‖mertensSum x − log log x‖` bounded. Alias of `mertensSum_sub_log_log_bounded`. |
| `mertens_M1`, `mertens_M3_exists`, `mertens_M3` | same file | Full Mertens asymptotic. **Sorry-free.** |

**PT pipeline note**: the article M1 (`m3_convergence_fr.tex`) only
uses Mertens for the *compactness* step `α_k = 1/2 − O(1/ln p_k)`.
The relevant theorem for the PT chain is `pt_T5_mertens_compactness`.
The full asymptotic is now also formalised end-to-end.

### EML (`PT/EML/`)

| Theorem | File | Statement |
|---|---|---|
| **`q±` Sheffer asymmetry** | `PT/EML/QSheffer.lean` | Two-branch Sheffer-class identity. |

### Analysis (`PT/Analysis/`)

| Theorem | File | Statement |
|---|---|---|
| **W7-1 (spiral identity, forward)** | `PT/Analysis/W7SpiralIdentity.lean` | Spiral identity of the Weil prime sum: under the PNT continuum density, the integral over the `k`-th turn of the log-polar Archimedean spiral equals the integral over turn 0 iff `σ² = π(k+1)`. Yields the cascade `σ_crit^(k) = √(π(k+1))`. |
| **W7-1 (reverse, `k ≥ 1`)** | `PT/Analysis/W7SpiralIdentityReverse.lean` | The converse: `J σ k = J σ 0 ⟹ σ² = π(k+1)` for `k ≥ 1`, proved via strict monotonicity (mean-value theorem on the Gaussian window). The `k = 0` case is degenerate and excluded explicitly. |

The forward proof (~10 tactic lines) uses substitution `v = x − σ²`,
parity of the centred Gaussian, an `intervalIntegral.integral_congr`
swap, and reflection `v ↦ −v`. The reverse direction (333 lines) uses
`gaussWindow σ a := ∫_a^{a+2π} gaussKer σ v`, FTC for the derivative,
and `strictMonoOn` / `strictAntiOn` on the two half-lines around `a = -π`.

Source paper: `PT_PROJECTS/PT_CH/paper_phase1/preprint1_cusp_fr.md`
§6.6 (Theorem 6.3). Numerical validation (out of Lean): k = 1, 2, 3 at
−0.626 %, −0.228 %, −0.085 % relative error respectively, with
precision growing in `k` (via the `primesieve` CLI for k = 3 over
3.4 × 10⁹ primes).

### CRT decoupling (`PT/CrtDecoupling/`)

A self-contained subsystem with its own [`PT/CrtDecoupling/README.md`](PT/CrtDecoupling/README.md).
Covers the geometric decoupling theorem 6.1 of the PT–CRT preprint
(`PT_ARTICLES/PUBLICATION/11_PT_CRT_Decoupling/`), with `0 sorry, 0 axiom`
through the abstract algebraic stage, the empirical stage, the spectral
reduction, the infinite-state ergodic strengthening, and a concrete
sieve dynamical-system instance.

### Other secondary modules

`PT/Algebra/`, `PT/Anomaly/`, `PT/Appendix/`, `PT/Bridge/`, `PT/CH/`,
`PT/Conservation/` (extensions and corollaries beyond the T2/ConservationID
core), … carry the 171 secondary modules (across 2 701 declarations total) that elaborate the
ecosystem but are not on the critical path.

## Build instructions

Requires `elan` (the Lean toolchain manager). The Lean version is pinned
in `lean-toolchain`.

```sh
cd PT_LEAN
lake exe cache get      # download Mathlib oleans (~1 GB, ~5 min)
lake build              # compile the PT library (~5 min cold, ~10 s cached)
```

To check a specific module:

```sh
lake build PT.Sieve.T1ForbiddenTransitions
lake build PT.Sieve.T3Antidiagonal
lake build PT.Stochastic.SHalf
lake build PT.Conservation.T2Alpha
lake build PT.Information.L0MaxEntropy
lake build PT.FixedPoint.T7MuStar
lake build PT.Analysis.W7SpiralIdentity
lake build PT           # umbrella import
```

## Roadmap

Across 187 modules, a single documented `sorry` remains. PT's core
derivation chain is unaffected.

| Open goal | File | Note |
|---|---|---|
| Fisher metric uniqueness up to scale (classical) | `PT/Information/G3FisherUniqueness.lean` | The monograph already marks this result `\leanExternal` (a classical theorem). Closing it would internalise the only remaining external dependency. |

The W7-1 reverse direction, previously listed as open, was formalised
on 2026-05-18 in `PT/Analysis/W7SpiralIdentityReverse.lean` (strict
monotonicity via MVT). Other gaps that the monograph still flags but
that are not on the critical path: N1 part 2 (uniqueness of the
constructive Eratosthenes procedure — meta-mathematical), and the
asymptotic equidistribution step of the bimodality theorem (likely
conditional on Hardy–Littlewood, marked `[COND]`-style if formalised).

Closing G3 does not affect the PT derivation chain itself.

### Downstream physics and chemistry

The PT physics and chemistry modules (`PT_PHYSICS`, `PT_CHEMISTRY`)
live *downstream* of these foundational results. Their numerical
derivations live in Python (`PUBLIC/PT_PHYSICS`, `PUBLIC/PT_CHEMISTRY`);
the Lean package captures only the foundational algebra and arithmetic.
Anything that ultimately reduces to `α_EM = ∏ sin² θ_p`,
`sin² θ_W = γ_7² / Σ γ_p²`, etc. is in principle formalisable once L0,
T5, and T7 are in place.

## What Lean caught during the initial bring-up

Six concrete issues, none of which were visible to human review:

1. `PT.lean`: a `/-! ... -/` module docstring appearing **before** the
   `import` block — illegal in Lean 4. The hand-written file had imports
   at line 22, not at the top.
2. `T3Antidiagonal.lean`: import path `Mathlib.Data.Matrix.Notation` does
   not exist in current Mathlib; the file is
   `Mathlib.LinearAlgebra.Matrix.Notation`.
3. `T3Antidiagonal.lean`: `Matrix.dotProduct` is not in the `Matrix`
   namespace; the correct symbol is the top-level `dotProduct` (defined
   in `Mathlib.Data.Matrix.Mul`, which must be imported explicitly).
4. `T3_det_neg_one`: the tactic `rw [Matrix.det_fin_two]; ring` does not
   close the goal because the rewritten
   `M 0 0 * M 1 1 − M 0 1 * M 1 0` doesn't automatically reduce concrete
   index look-ups; `simp [T3, Matrix.det_fin_two]` works.
5. `SHalf.lean`: `def piHalf : Fin 2 → ℝ := ![1/2, 1/2]` must be marked
   `noncomputable` because `1/2 : ℝ` depends on
   `Real.instDivInvMonoid`, which is noncomputable.
6. Two transient race-condition build failures when downstream modules
   tried to read `.olean` files that were still being written by a
   parallel build. These are flaky Lake errors, fixed by re-running
   `lake build`.

Items 1–5 are exactly the class of "small" mistakes that paper proofs
hide. The dependency on `Mathlib.Data.Matrix.Mul` (item 3) in particular
would not have been noticed without Lean's strict import resolution.

## References

- Y. Senez, *Forbidden transitions and informational structure of prime
  residues under modular projection* — `PT_ARTICLES/PT_MATHEMATICS/M1/`.
  Theorems T1, T3, `s = 1/2`, and T2 are the headline results of this
  paper; this Lean package is the formal verification.
- Y. Senez, *La théorie de la persistance* (monograph) —
  `PT_MONOGRAPHY/`, Chapter 1 §§1.3–1.4, Appendix Z for the W7-1 proof
  in the explicit-formula context.
- Y. Senez, *Berry–Keating operator at the cusp $\Sigma_{\rm pers}$* —
  `PT_PROJECTS/PT_CH/paper_phase1/preprint1_cusp_fr.md`, §6.6 for the
  spiral identity theorem.
- Mathlib (`master` @ `e12bcd0`) —
  https://github.com/leanprover-community/mathlib4

## License

Apache 2.0, matching Mathlib's convention.
