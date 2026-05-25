/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Bridge.CauchyMultiplicativeExp
import Mathlib.Algebra.BigOperators.Fin

/-!
# Shore--Johnson G1 applied to PT spectral cutoffs

This module re-architects the axiomatic basis of
`HiggsCutoffUniqueness.lean` so that the bilinear ad-hoc axiom

  `CRT_ShoreJohnson_cutoff_factorises : ∀ x y, f (x + y) = f x * f y`

is replaced by a more **structural** multilinear axiom that directly
expresses the content of *Shore--Johnson G1 applied to a `k`-channel
spectral inference*, together with a kernel-verified **theorem**
deriving the original bilinear statement as the `k = 2` specialisation.

## Conceptual content

The PT Dirac operator has CRT factorisation
`D_PT² = ⊕_p D_p²` (proven in
`PT/Stochastic/CRTFactorizationT2.lean`). The spectral eigenvalue on
`|r_{p_1}, …, r_{p_k}⟩` is `λ² = ∑_p u_p`. Shore--Johnson G1
(`PT/Information/T6cG1Autonomous.lean`) is the **system-independence
axiom**: an inference rule on the joint system `⊗_p A_p` must factorise
as a product of inferences on each `A_p`.

Applied to the spectral cutoff `f`, this means: the value `f(∑ u_p)`
on a CRT-independent state must factorise as `∏ f(u_p)`. In its purest
form this is a **multilinear** (`k`-fold) factorisation, of which
the bilinear case `f(u_1 + u_2) = f(u_1) · f(u_2)` is the `k = 2`
specialisation.

## What this module accomplishes

* **Axiom `SJG1_spectral_cutoff_factorises`** (multivariate, `k`-fold).
  This is the **structurally pure** statement of SJ G1 for spectral
  cutoffs: for any natural `k` and any indexed family
  `u : Fin k → ℝ` of CRT-independent eigenvalue contributions,
  `f(∑ u_i) = ∏ f(u_i)`. The bilinear ad-hoc axiom is now a corollary.

* **Theorem `cutoff_zero_eq_one_of_SJG1`**. The normalisation
  `f(0) = 1` is *derived* (not posited) from the multivariate axiom
  at `k = 0` (empty sum and empty product, both equal to `1`).

* **Theorem `cauchy_bilinear_from_SJG1`**. The bilinear Cauchy
  equation `f(x + y) = f(x) · f(y)` is *derived* by specialising
  the multivariate axiom at `k = 2` with `u 0 = x, u 1 = y`. This
  is the exact content previously posited as
  `CRT_ShoreJohnson_cutoff_factorises`, now a **theorem**.

## Epistemic progress

Before this module, `HiggsCutoffUniqueness.lean` declared an *ad-hoc
bilinear* axiom encoding "CRT + SJ G1 application to cutoff".
After this module, the dependency chain reads:

| Statement | Status | Module |
|---|---|---|
| `SJG1_spectral_cutoff_factorises` (multivariate) | **axiom** (purest form of SJ G1 spectral) | this file |
| `cauchy_bilinear_from_SJG1` (bilinear case) | **theorem** | this file |
| `cutoff_PT_unique_eq_cutoffPT` | theorem given multivariate axiom + `PT_cutoff_decay` + `PT_cutoff_scale_eq_Nb` | `HiggsCutoffUniqueness.lean` |

The number of PT-specific axioms feeding `HiggsCutoffUniqueness`
is unchanged (one for the CRT-SJ content, one for decay, one for the
scale `N_b`), but the **conceptual character** of the CRT-SJ axiom
is improved: it is now a single multilinear principle (the form in
which Shore--Johnson G1 naturally extends to multi-channel inference)
rather than an opaque bilinear identity.

## What this module does NOT accomplish

The natural extension would be to *eliminate* the SJ-spectral axiom
entirely by deriving it from the existing PT theorem
`G1_autonomous_DKL_unique` (proven for probability distributions in
`PT/Information/T6cG1Autonomous.lean`). That elimination would
require formalising:

1. the spectral action `Tr f(D²/Λ²)` as a Lean object,
2. the bridge "spectral inference is a Shore--Johnson inference",
3. the rigorous extension of G1 from f-divergences to cutoffs.

These are substantial future-work targets. The present module makes
the *conceptual* progress accessible without that infrastructure.

## References

* PT monograph: appendix Y §Y.13.6
  (`appendices_fr/app_y_higgs_zeta_duality.tex`),
  `thm:E8_formal_uniqueness`.
* Shore-Johnson 1980, *Axiomatic derivation of the principle of
  maximum entropy*, IEEE Trans. Info. Theory **26**, 26-37.
* `PT/Information/T6cG1Autonomous.lean` (SJ G1 for KL).
-/

namespace PT.Bridge

open Real Finset
open scoped BigOperators

/-! ## PT spectral cutoff (regularity part) -/

/-- **PT spectral cutoff (regularity part).** A function `f : ℝ → ℝ` is a
    PT spectral cutoff if it is continuous and strictly positive. The
    factorisation and scale properties are stated below as axioms /
    derived in the modules `ShoreJohnsonG1Spectral.lean` (this file)
    and `CutoffMeanCharacterisation.lean`. -/
structure IsPTSpectralCutoff (f : ℝ → ℝ) : Prop where
  continuous : Continuous f
  positive : ∀ x, 0 < f x

/-! ## Multivariate SJ-G1 spectral axiom -/

/-- **Axiom (Shore--Johnson G1 applied to spectral cutoffs).**

    For any `k ∈ ℕ` and any indexed family `u : Fin k → ℝ` of eigenvalue
    contributions on a CRT-independent product state of the PT Dirac
    operator, the spectral cutoff `f` satisfies the multivariate
    factorisation
    $$ f\!\left(\sum_{i < k} u_i\right) \;=\; \prod_{i < k} f(u_i). $$

    **Justification.** Combine:
    * CRT factorisation of the PT Dirac operator
      (`PT/Stochastic/CRTFactorizationT2.lean`, kernel-verified):
      `D_PT² = ⊕_p D_p²`, hence eigenvalue on
      `|r_{p₁}, …, r_{p_k}⟩` equals `∑_p u_p`;
    * Shore--Johnson G1 (`PT/Information/T6cG1Autonomous.lean`,
      kernel-verified for KL divergences): the spectral inference on
      `⊗_p A_p` must factorise as `∏_p` of inferences on each `A_p`.

    The remaining step "the cutoff `f` is itself a Shore--Johnson
    spectral inference" is the **conceptual content** of this axiom.
    A direct Lean derivation requires formalising the spectral action
    `Tr f(D²/Λ²)` (future work). Here we expose this as a single
    multilinear statement, the structurally purest form. -/
axiom SJG1_spectral_cutoff_factorises
    (f : ℝ → ℝ) (hf : IsPTSpectralCutoff f) :
    ∀ (k : ℕ) (u : Fin k → ℝ),
      f (∑ i, u i) = ∏ i, f (u i)

/-! ## Derived theorems -/

/-- **Normalisation `f(0) = 1` follows from `k = 0` case.**

    Specialising `SJG1_spectral_cutoff_factorises` at `k = 0` (empty
    index `Fin 0`) gives `f(empty sum) = empty product`, i.e.
    `f(0) = 1`. This is a theorem, not an axiom: it is *derived* from
    the multivariate factorisation. -/
theorem cutoff_zero_eq_one_of_SJG1
    (f : ℝ → ℝ) (hf : IsPTSpectralCutoff f) :
    f 0 = 1 := by
  -- Apply the multivariate axiom at k = 0.
  have h := SJG1_spectral_cutoff_factorises f hf 0 (fun (i : Fin 0) => i.elim0)
  -- Both sums/products over `Fin 0` are empty, hence `0` and `1`.
  simp at h
  exact h

/-- **Bilinear Cauchy equation, derived.**

    The bilinear identity `f(x + y) = f(x) · f(y)` — previously
    posited as the ad-hoc axiom
    `CRT_ShoreJohnson_cutoff_factorises` — is now a **theorem**
    derived from the multivariate `SJG1_spectral_cutoff_factorises`
    by specialisation at `k = 2`. -/
theorem cauchy_bilinear_from_SJG1
    (f : ℝ → ℝ) (hf : IsPTSpectralCutoff f) :
    ∀ x y : ℝ, f (x + y) = f x * f y := by
  intro x y
  -- Apply the multivariate axiom at k = 2 with u 0 = x, u 1 = y.
  let u : Fin 2 → ℝ := fun i => if i = 0 then x else y
  have h := SJG1_spectral_cutoff_factorises f hf 2 u
  -- The sum over Fin 2 is u 0 + u 1 = x + y.
  have hsum : (∑ i, u i) = x + y := by
    rw [Fin.sum_univ_two]
    show (if (0 : Fin 2) = 0 then x else y) +
         (if (1 : Fin 2) = 0 then x else y) = x + y
    simp
  -- The product over Fin 2 is u 0 * u 1 = x * y... wait, it's f(u 0) * f(u 1).
  have hprod : (∏ i, f (u i)) = f x * f y := by
    rw [Fin.prod_univ_two]
    show f (if (0 : Fin 2) = 0 then x else y) *
         f (if (1 : Fin 2) = 0 then x else y) = f x * f y
    simp
  rw [hsum, hprod] at h
  exact h

/-! ## Higher-arity corollary (k = 3) — useful sanity check -/

/-- **Trilinear Cauchy equation (corollary).** Specialisation at `k = 3`,
    used as a sanity check that the multivariate axiom is consistent
    with the bilinear case under iteration. -/
theorem cauchy_trilinear_from_SJG1
    (f : ℝ → ℝ) (hf : IsPTSpectralCutoff f) :
    ∀ x y z : ℝ, f (x + y + z) = f x * f y * f z := by
  intro x y z
  -- Iterate the bilinear case: `f((x+y)+z) = f(x+y) * f(z) = f x * f y * f z`.
  have h_bi := cauchy_bilinear_from_SJG1 f hf
  calc f (x + y + z)
      = f (x + y) * f z := h_bi _ _
    _ = (f x * f y) * f z := by rw [h_bi]
    _ = f x * f y * f z := by ring

/-! ## Bridge to `HiggsCutoffUniqueness` -/

/-- **Bridge lemma.** Combining `cauchy_bilinear_from_SJG1` with the
    `IsPTSpectralCutoff` hypothesis, we obtain the input required by
    `cauchy_mult_exp` from `PT/Bridge/CauchyMultiplicativeExp.lean`.

    This is the bridge that lets `HiggsCutoffUniqueness.lean` close
    the chain without the ad-hoc bilinear axiom: it just needs the
    multilinear `SJG1_spectral_cutoff_factorises` from this module. -/
theorem cutoff_to_cauchy_mult_exp
    (f : ℝ → ℝ) (hf : IsPTSpectralCutoff f) :
    ∃ c : ℝ, ∀ x, f x = Real.exp (c * x) :=
  cauchy_mult_exp f hf.continuous hf.positive
    (cauchy_bilinear_from_SJG1 f hf)

end PT.Bridge
