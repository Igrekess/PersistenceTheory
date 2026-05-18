/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Stochastic.T2T30Normalisation

/-!
# T_5 canonical (tight) — concrete instance of `T5Like` saturating `|λ₂| = 1/4`

**Goal.** The companion file `PT.Stochastic.T5Canonical` instantiates `T5Like`
with the rank-one uniform mixer `(1/4) · J`, whose subdominant eigenvalue is
exactly `0`. That choice trivially satisfies the structural bound
`|λ₂| ≤ 1/4 = s²` but does so with *strict slack*.

The PT monograph (ch03 §"Exposant de conservation spectrale") records the
**tight** identity `|λ₂(T_5)| = 1/4 = s²` for the empirical prime-gap
transfer matrix. This file delivers a Lean witness that **saturates** that
bound on the nose, with all eigendata verified by direct calculation.

## Tight choice of canonical matrix

We pick

  `T_5^{tight} := (3/16) · J + (1/4) · I`,

i.e. the row-stochastic 4×4 matrix with diagonal entries `7/16` and
off-diagonal entries `3/16`. Each row sums to `7/16 + 3·(3/16) = 16/16 = 1`,
so the matrix is row-stochastic. Its spectrum is `{1, 1/4, 1/4, 1/4}`:

* Perron eigenpair `(1, (1,1,1,1))` — each row sums to `1`.
* Subdominant eigenpair `(1/4, (1,-1,0,0))` — and the same eigenvalue is
  shared by `(1,0,-1,0)` and `(1,0,0,-1)`, since `(3/16) · J + (1/4) · I`
  on the zero-sum hyperplane reduces to `(1/4) · I`.

In particular, `|λ_2| = 1/4` **exactly**: the PT spectral bound is saturated.

## Main results

* `T5_tight : Matrix (Fin 4) (Fin 4) ℝ` — the tight canonical matrix.
* `T5_tight_T5Like : T5Like` — the `T5Like` record built from `T5_tight`
  with Perron eigenpair `(1, (1,1,1,1))` and subdominant eigenpair
  `(1/4, (1,-1,0,0))`.
* `T5_tight_subEigenvalue_eq` — the witnessing identity
  `T5_tight_T5Like.subEigenvalue = 1/4`.
* `T30_tight_lambda_eff_bound_s_sq` — the unconditional headline
  `|λ_2^{eff}(T_30^{tight})| ≤ s²`, **saturated**: the bound holds with
  equality (modulo the implicit `s² = 1/4`).

## Comparison with `T5_canonical`

| Instance        | matrix                       | `subEigenvalue` | bound `≤ 1/4` |
|-----------------|------------------------------|------------------|---------------|
| `T5_canonical`  | `(1/4) · J` (rank-one mixer) | `0`              | strict slack  |
| `T5_tight`      | `(3/16) · J + (1/4) · I`     | `1/4`            | saturated     |

Both make the headline T2 bound unconditional; `T5_tight` is the matching
witness for the **tight** monograph identity `|λ_2(T_5)| = 1/4`.

## Reference

* Monograph Chapter 3, §"Spectral conservation T2".
* `PT.Stochastic.T2T30Normalisation` — the conditional theorem
  parameterised by `T5Like`.
* `PT.Stochastic.T5Canonical` — the rank-one slack-witness companion.
-/

namespace PT.Stochastic

open Matrix BigOperators

/-! ### The tight canonical matrix `T_5^{tight}` -/

/-- **Tight canonical `T_5` matrix.** The row-stochastic 4×4 matrix with
    diagonal `7/16` and off-diagonal `3/16`. Equivalently
    `(3/16) · J + (1/4) · I`, where `J` is the all-ones matrix. Its
    spectrum is `{1, 1/4, 1/4, 1/4}` — saturating the PT structural bound
    `|λ_2| ≤ 1/4 = s²`. -/
noncomputable def T5_tight : Matrix (Fin 4) (Fin 4) ℝ := !![
  7/16, 3/16, 3/16, 3/16;
  3/16, 7/16, 3/16, 3/16;
  3/16, 3/16, 7/16, 3/16;
  3/16, 3/16, 3/16, 7/16]

/-! ### Perron eigenpair for `T_5^{tight}` -/

/-- Perron eigenvector of `T_5^{tight}` — the all-ones vector. -/
def T5_tight_perronVec : Fin 4 → ℝ := fun _ => 1

@[simp] lemma T5_tight_perronVec_apply (i : Fin 4) :
    T5_tight_perronVec i = 1 := rfl

/-- **Perron eigenvalue equation for `T_5^{tight}`.** Each row of
    `T_5^{tight}` sums to `1`, so the all-ones vector is fixed:
    `T_5^{tight} · (1,1,1,1) = (1,1,1,1)`. -/
theorem T5_tight_perron_eigen :
    T5_tight.mulVec T5_tight_perronVec = (1 : ℝ) • T5_tight_perronVec := by
  ext i
  fin_cases i <;>
    simp [T5_tight, T5_tight_perronVec, Matrix.mulVec, dotProduct,
          Fin.sum_univ_four] <;>
    norm_num

/-- **Perron eigenvector positivity.** Every component of `(1,1,1,1)`
    is strictly positive. -/
theorem T5_tight_perronVec_pos :
    ∀ i, 0 < T5_tight_perronVec i := by
  intro i
  fin_cases i <;> norm_num [T5_tight_perronVec]

/-! ### Subdominant eigenpair for `T_5^{tight}` -/

/-- A subdominant eigenvector of `T_5^{tight}`: the vector `(1, -1, 0, 0)`.
    Since `T_5^{tight} - (1/4) · I = (3/16) · J` has every zero-sum vector
    in its kernel, this is an eigenvector with eigenvalue `1/4`. -/
def T5_tight_subVec : Fin 4 → ℝ := ![1, -1, 0, 0]

@[simp] lemma T5_tight_subVec_zero : T5_tight_subVec 0 = 1 := rfl
@[simp] lemma T5_tight_subVec_one : T5_tight_subVec 1 = -1 := rfl
@[simp] lemma T5_tight_subVec_two : T5_tight_subVec 2 = 0 := rfl
@[simp] lemma T5_tight_subVec_three : T5_tight_subVec 3 = 0 := rfl

/-- **Subdominant eigenvalue of `T_5^{tight}`.** Exactly `1/4`, saturating
    the PT structural bound `|λ_2(T_5)| ≤ 1/4 = s²`. -/
noncomputable def T5_tight_subEigenvalue : ℝ := 1 / 4

/-- **Subdominant eigenvalue equation.** `T_5^{tight} · (1,-1,0,0) =
    (1/4) · (1,-1,0,0)`. Row-by-row:

    * row 0: `7/16 - 3/16 + 0 + 0 = 4/16 = 1/4`
    * row 1: `3/16 - 7/16 + 0 + 0 = -4/16 = -1/4`
    * row 2: `3/16 - 3/16 + 0 + 0 = 0`
    * row 3: `3/16 - 3/16 + 0 + 0 = 0`
-/
theorem T5_tight_sub_eigen :
    T5_tight.mulVec T5_tight_subVec
      = T5_tight_subEigenvalue • T5_tight_subVec := by
  ext i
  fin_cases i <;>
    simp [T5_tight, T5_tight_subVec, T5_tight_subEigenvalue,
          Matrix.mulVec, dotProduct, Fin.sum_univ_four] <;>
    norm_num

/-- **PT structural bound on the subdominant eigenvalue — saturated.**
    `|1/4| = 1/4 ≤ 1/4 = s²`. The bound holds with equality. -/
theorem T5_tight_subdominant_bound :
    |T5_tight_subEigenvalue| ≤ (1 : ℝ) / 4 := by
  unfold T5_tight_subEigenvalue
  norm_num

/-! ### The tight canonical `T5Like` instance -/

/-- **The tight canonical `T5Like` instance.** Bundles `T_5^{tight}`
    together with its Perron eigenpair `(1, (1,1,1,1))` and its
    subdominant eigenpair `(1/4, (1,-1,0,0))`, saturating the PT
    structural bound `|λ_2| ≤ 1/4`. -/
noncomputable def T5_tight_T5Like : T5Like where
  matrix := T5_tight
  perronVec := T5_tight_perronVec
  perron_eigen := T5_tight_perron_eigen
  perronVec_pos := T5_tight_perronVec_pos
  subVec := T5_tight_subVec
  subEigenvalue := T5_tight_subEigenvalue
  sub_eigen := T5_tight_sub_eigen
  subdominant_bound := T5_tight_subdominant_bound

/-- **Subdominant eigenvalue of `T5_tight_T5Like` is exactly `1/4`.** This
    is the witness that the structural bound `|λ_2| ≤ 1/4` is *saturated*
    by this canonical instance, in contrast with `T5_canonical_T5Like`
    where the subdominant eigenvalue is `0`. -/
@[simp] theorem T5_tight_subEigenvalue_eq :
    T5_tight_T5Like.subEigenvalue = (1 : ℝ) / 4 := rfl

/-! ### Headline unconditional T2 bound (tight) -/

/-- **Unconditional T2 spectral bound (tight canonical instance).**
    Specialising the conditional theorem
    `PT.Stochastic.T30_lambda_eff_bound_s_sq` to the tight `T5Like`
    instance gives the **unconditional** headline bound

      `|λ_2^{eff}(T_30^{tight})|  ≤  s²` ,

    where `T_30^{tight} := T_2 ⊗ T_3 ⊗ T_5^{tight}` is the three-factor
    Kronecker product built from the tight canonical `T_5`. This bound
    is *saturated*: by `T5_tight_subEigenvalue_eq`, the subdominant
    eigenvalue equals `1/4 = s²` exactly. -/
theorem T30_tight_lambda_eff_bound_s_sq :
    |T5_tight_T5Like.subEigenvalue| ≤ s ^ 2 :=
  T30_lambda_eff_bound_s_sq T5_tight_T5Like

/-- **Unconditional T2 bound, `1/4` form.** Same statement as
    `T30_tight_lambda_eff_bound_s_sq` but with the bound spelled as
    `1/4`. -/
theorem T30_tight_lambda_eff_bound_one_quarter :
    |T5_tight_T5Like.subEigenvalue| ≤ (1 : ℝ) / 4 :=
  T30_lambda_eff_abs_bound T5_tight_T5Like

/-- **Unconditional T2 bound, conservation-exponent form.** The
    tight subdominant eigenvalue is bounded by the PT conservation
    exponent `α_cons = s² = 1/4`. -/
theorem T30_tight_lambda_eff_bound_alpha_cons :
    |T5_tight_T5Like.subEigenvalue| ≤ PT.Conservation.alpha_cons :=
  T30_lambda_eff_bound_alpha_cons T5_tight_T5Like

/-- **Saturating identity.** The absolute value of the subdominant
    eigenvalue of the tight canonical instance is exactly `1/4`, matching
    the monograph identity `|λ_2(T_5)| = 1/4 = s²` on the nose. -/
@[simp] theorem T30_tight_lambda_eff_abs_eq_one_quarter :
    |T5_tight_T5Like.subEigenvalue| = (1 : ℝ) / 4 := by
  simp [T5_tight_T5Like, T5_tight_subEigenvalue]

/-! ### Full unconditional T2 summary (tight) -/

/-- **Unconditional T2 headline summary (tight canonical instance).**

    For the canonical three-factor Kronecker product
    `T_30^{tight} := T_2 ⊗ T_3 ⊗ T_5^{tight}` (with `T_5^{tight}` the
    tight 4×4 matrix `(3/16) · J + (1/4) · I`), the following hold:

    1. Perron eigenpair: `T_30^{tight} · u_+ = u_+`.
    2. Strict positivity of the Perron eigenvector.
    3. Subdominant eigenpair on the Perron-symmetric sector, with
       eigenvalue `1/4` saturating the PT structural bound.
    4. Unconditional spectral bound `|λ_2^{eff}(T_30^{tight})| ≤ s² = 1/4`.

    This is the same statement as `T30_T2_summary` in
    `PT.Stochastic.T2T30Normalisation`, but **specialised** to the
    tight `T_5` so that no `T5Like` assumption is left dangling and the
    structural bound is saturated. -/
theorem T30_tight_T2_summary :
    -- 1. Perron eigenpair
    (T30 T5_tight_T5Like).mulVec (T30_perronVec T5_tight_T5Like)
      = (1 : ℝ) • T30_perronVec T5_tight_T5Like
    -- 2. Perron eigenvector is strictly positive
    ∧ (∀ ijk, 0 < T30_perronVec T5_tight_T5Like ijk)
    -- 3. Lambda-eff eigenpair
    ∧ (T30 T5_tight_T5Like).mulVec (T30_lambda_eff_vec T5_tight_T5Like)
        = T5_tight_T5Like.subEigenvalue
            • T30_lambda_eff_vec T5_tight_T5Like
    -- 4. Unconditional spectral bound `|λ_2^{eff}| ≤ s²`
    ∧ |T5_tight_T5Like.subEigenvalue| ≤ s ^ 2 :=
  T30_T2_summary T5_tight_T5Like

end PT.Stochastic
