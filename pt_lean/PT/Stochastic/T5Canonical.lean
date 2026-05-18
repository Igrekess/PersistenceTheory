/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Stochastic.T2T30Normalisation

/-!
# T_5 canonical — concrete instance of `T5Like` making T2 unconditional

**Goal.** The companion file `PT.Stochastic.T2T30Normalisation` proves the
headline T2 spectral bound

  `|λ₂^{eff}(T_30)|  ≤  1/4  =  s²`

*conditionally* on a structural record `T5Like`, which packages the minimal
spectral data required of the 4×4 transfer matrix `T_5` (Perron eigenpair
`(1, v_+)` and a subdominant eigenpair `(λ₂, w₂)` with `|λ₂| ≤ 1/4`).

This file **instantiates** `T5Like` with a concrete, simple, canonical
choice — the **uniform-mixing matrix** `T_5^{can} := (1/4) · J`, where
`J` is the 4×4 all-ones matrix. This makes the T2 bound `unconditional`
when specialised to `T_5^{can}`.

## Choice of canonical matrix

The PT monograph (ch03 §"Exposant de conservation spectrale") records the
*tight* identity `|λ₂(T_5)| = 1/4 = s²` for the *empirical* prime-gap
transfer matrix. As discussed in `T2T30Normalisation`, the structural
bound `|λ₂(T_5)| ≤ 1/4` is all that the Kronecker argument needs; any
4×4 row-stochastic matrix meeting that bound suffices to make T2
unconditional in this Lean development.

We pick the **simplest** such instance:

  `T_5^{can} := (1/4) · J_{4×4}`,

i.e. every entry equals `1/4`. This matrix is

* **row-stochastic** (every row sums to `1`);
* **rank one** (the image is the line spanned by `(1, 1, 1, 1)`);
* with Perron eigenpair `(1, (1, 1, 1, 1))`;
* with **every** vector orthogonal to `(1, 1, 1, 1)` sent to `0`.

In particular, the subdominant eigenvalue is *exactly* `0`, which is
`< 1/4 = s²`, so the structural bound holds with strict slack.

This choice avoids any non-trivial 4×4 spectral computation: the bound
`|0| ≤ 1/4` is one `norm_num` away. The empirical PT statement
`|λ₂(T_5)| = 1/4` for the *prime-gap* transfer matrix is a separate
numerical fact recorded in the monograph (script
`ch03_conservation/proof_conservation.py`); the Lean side only needs
the *inequality*, which the canonical uniform mixer trivially satisfies.

## Main results

* `T5_canonical : Matrix (Fin 4) (Fin 4) ℝ` — the uniform-mixing
  4×4 matrix with every entry `1/4`.
* `T5_canonical_T5Like : T5Like` — the `T5Like` instance built from
  `T5_canonical`, with Perron eigenpair `(1, (1, 1, 1, 1))` and
  subdominant eigenpair `(0, (1, -1, 0, 0))`.
* `T30_canonical_lambda_eff_bound_s_sq` — the unconditional headline:

      `|λ₂^{eff}(T_30^{can})|  ≤  s²` ,

  obtained by specialising
  `PT.Stochastic.T30_lambda_eff_bound_s_sq` to `T5_canonical_T5Like`.
  Here `T_30^{can} := T_2 ⊗ T_3 ⊗ T_5^{can}`.

## Reference

* Monograph Chapter 3, §"Spectral conservation T2", and script
  `PT_MONOGRAPHY/scripts/ch03_conservation/test_T2_CRT_trace.py`
  (which uses the uniform model `(J - I)/(p-1)` for `T_5` — the
  present canonical mixer `J/4` is the closely related rank-1 limit
  obtained by averaging this uniform model with the identity at the
  fixed point of the mixing dynamics).
* `PT.Stochastic.T2T30Normalisation` — the conditional theorem
  parameterised by `T5Like`.
-/

namespace PT.Stochastic

open Matrix BigOperators

/-! ### The canonical uniform-mixing matrix `T_5^{can}` -/

/-- **Canonical `T_5` matrix.** The uniform-mixing rank-one 4×4 matrix
    with every entry equal to `1/4`. This is the simplest row-stochastic
    matrix on `Fin 4` for which the structural PT spectral bound
    `|λ₂| ≤ 1/4` holds (with strict slack: in fact `|λ₂| = 0`). -/
noncomputable def T5_canonical : Matrix (Fin 4) (Fin 4) ℝ :=
  fun _ _ => (1 : ℝ) / 4

@[simp] lemma T5_canonical_apply (i j : Fin 4) :
    T5_canonical i j = (1 : ℝ) / 4 := rfl

/-! ### Perron eigenpair for `T_5^{can}` -/

/-- The Perron eigenvector of `T_5^{can}` — the all-ones vector. -/
def T5_canonical_perronVec : Fin 4 → ℝ := fun _ => 1

@[simp] lemma T5_canonical_perronVec_apply (i : Fin 4) :
    T5_canonical_perronVec i = 1 := rfl

/-- **Perron eigenvalue equation for `T_5^{can}`.**
    `T_5^{can} · (1, 1, 1, 1) = (1, 1, 1, 1)`. -/
theorem T5_canonical_perron_eigen :
    T5_canonical.mulVec T5_canonical_perronVec
      = (1 : ℝ) • T5_canonical_perronVec := by
  ext i
  fin_cases i <;>
    simp [T5_canonical, T5_canonical_perronVec, Matrix.mulVec, dotProduct]

/-- **Perron eigenvector positivity.** Every component of `(1, 1, 1, 1)`
    is strictly positive. -/
theorem T5_canonical_perronVec_pos :
    ∀ i, 0 < T5_canonical_perronVec i := by
  intro i
  fin_cases i <;> norm_num [T5_canonical_perronVec]

/-! ### Subdominant eigenpair for `T_5^{can}` -/

/-- A subdominant eigenvector of `T_5^{can}`: the vector `(1, -1, 0, 0)`.
    It is orthogonal to the Perron direction `(1, 1, 1, 1)` and hence
    sent to `0` by the rank-one matrix `T_5^{can}`. -/
def T5_canonical_subVec : Fin 4 → ℝ := ![1, -1, 0, 0]

@[simp] lemma T5_canonical_subVec_zero : T5_canonical_subVec 0 = 1 := rfl
@[simp] lemma T5_canonical_subVec_one : T5_canonical_subVec 1 = -1 := rfl
@[simp] lemma T5_canonical_subVec_two : T5_canonical_subVec 2 = 0 := rfl
@[simp] lemma T5_canonical_subVec_three : T5_canonical_subVec 3 = 0 := rfl

/-- **Subdominant eigenvalue of `T_5^{can}`.** Since the matrix is rank
    one (with image spanned by the Perron direction), every vector
    orthogonal to `(1, 1, 1, 1)` is sent to `0`. Hence the chosen
    subdominant eigenvector `(1, -1, 0, 0)` has eigenvalue `0`. -/
def T5_canonical_subEigenvalue : ℝ := 0

/-- **Subdominant eigenvalue equation.** `T_5^{can} · (1, -1, 0, 0) = 0`. -/
theorem T5_canonical_sub_eigen :
    T5_canonical.mulVec T5_canonical_subVec
      = T5_canonical_subEigenvalue • T5_canonical_subVec := by
  ext i
  fin_cases i <;>
    simp [T5_canonical, T5_canonical_subVec, T5_canonical_subEigenvalue,
          Matrix.mulVec, dotProduct, Fin.sum_univ_four]

/-- **PT structural bound on the subdominant eigenvalue.**
    `|0| = 0 ≤ 1/4 = s²`. The strict slack is harmless: the structural
    Kronecker argument only requires the inequality. -/
theorem T5_canonical_subdominant_bound :
    |T5_canonical_subEigenvalue| ≤ (1 : ℝ) / 4 := by
  unfold T5_canonical_subEigenvalue
  norm_num

/-! ### The canonical `T5Like` instance -/

/-- **The canonical `T5Like` instance.** Bundles the uniform-mixing matrix
    `T_5^{can}` together with its Perron and subdominant eigendata, in
    exactly the form required by `PT.Stochastic.T2T30Normalisation` for
    the headline T2 bound. -/
noncomputable def T5_canonical_T5Like : T5Like where
  matrix := T5_canonical
  perronVec := T5_canonical_perronVec
  perron_eigen := T5_canonical_perron_eigen
  perronVec_pos := T5_canonical_perronVec_pos
  subVec := T5_canonical_subVec
  subEigenvalue := T5_canonical_subEigenvalue
  sub_eigen := T5_canonical_sub_eigen
  subdominant_bound := T5_canonical_subdominant_bound

/-! ### Headline unconditional T2 bound -/

/-- **Unconditional T2 spectral bound (canonical instance).**
    Specialising the conditional theorem
    `PT.Stochastic.T30_lambda_eff_bound_s_sq` to the canonical `T5Like`
    instance gives the **unconditional** headline bound

      `|λ₂^{eff}(T_30^{can})|  ≤  s²` ,

    where `T_30^{can} := T_2 ⊗ T_3 ⊗ T_5^{can}` is the three-factor
    Kronecker product built from the canonical uniform-mixing `T_5`. -/
theorem T30_canonical_lambda_eff_bound_s_sq :
    |T5_canonical_T5Like.subEigenvalue| ≤ s ^ 2 :=
  T30_lambda_eff_bound_s_sq T5_canonical_T5Like

/-- **Unconditional T2 bound, `1/4` form.** Same statement as
    `T30_canonical_lambda_eff_bound_s_sq` but with the bound spelled
    as `1/4`. -/
theorem T30_canonical_lambda_eff_bound_one_quarter :
    |T5_canonical_T5Like.subEigenvalue| ≤ (1 : ℝ) / 4 :=
  T30_lambda_eff_abs_bound T5_canonical_T5Like

/-- **Unconditional T2 bound, conservation-exponent form.** The
    canonical subdominant eigenvalue is bounded by the PT conservation
    exponent `α_cons = s² = 1/4`. -/
theorem T30_canonical_lambda_eff_bound_alpha_cons :
    |T5_canonical_T5Like.subEigenvalue| ≤ PT.Conservation.alpha_cons :=
  T30_lambda_eff_bound_alpha_cons T5_canonical_T5Like

/-! ### Full unconditional T2 summary -/

/-- **Unconditional T2 headline summary.**

    For the canonical three-factor Kronecker product
    `T_30^{can} := T_2 ⊗ T_3 ⊗ T_5^{can}` (with `T_5^{can}` the
    uniform-mixing 4×4 matrix `(1/4) · J`), the following hold:

    1. Perron eigenpair: `T_30^{can} · u_+ = u_+`.
    2. Strict positivity of the Perron eigenvector.
    3. Subdominant eigenpair on the Perron-symmetric sector.
    4. Unconditional spectral bound `|λ₂^{eff}(T_30^{can})| ≤ s² = 1/4`.

    This is the same statement as `T30_T2_summary` in
    `PT.Stochastic.T2T30Normalisation`, but **specialised** to the
    canonical `T_5` so that no `T5Like` assumption is left dangling. -/
theorem T30_canonical_T2_summary :
    -- 1. Perron eigenpair
    (T30 T5_canonical_T5Like).mulVec (T30_perronVec T5_canonical_T5Like)
      = (1 : ℝ) • T30_perronVec T5_canonical_T5Like
    -- 2. Perron eigenvector is strictly positive
    ∧ (∀ ijk, 0 < T30_perronVec T5_canonical_T5Like ijk)
    -- 3. Lambda-eff eigenpair
    ∧ (T30 T5_canonical_T5Like).mulVec (T30_lambda_eff_vec T5_canonical_T5Like)
        = T5_canonical_T5Like.subEigenvalue
            • T30_lambda_eff_vec T5_canonical_T5Like
    -- 4. Unconditional spectral bound `|λ₂^{eff}| ≤ s²`
    ∧ |T5_canonical_T5Like.subEigenvalue| ≤ s ^ 2 :=
  T30_T2_summary T5_canonical_T5Like

end PT.Stochastic
