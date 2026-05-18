/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Conservation.T2Alpha
import PT.Stochastic.T5CanonicalTight
import PT.Stochastic.T2T3T5KroneckerSpectrum

/-!
# T2 — Spectral bridge: `|λ_2(T_30)| = s² = 1/4 = α_cons`

This file closes **gap C.1** of `LEAN_MONOGRAPHIE_GAPS.md` by exhibiting
explicitly the bridge between

* the *algebraic* identity `α_cons = s · s = 1/4` already proved in
  `PT.Conservation.T2Alpha` (a trivial scalar fact), and
* the *spectral* identity

      `|λ_2(T_30)|  =  s²  =  1/4`

  obtained from the Chinese-Remainder Kronecker factorisation
  `T_30 ≅ T_2 ⊗ T_3 ⊗ T_5` (monograph Chapter 3, Remark `rem:T2_trace`).

The two-factor Kronecker route is fully formalised upstream in
`PT.Stochastic.CRTFactorizationT2`, the three-factor extension and the
PT subdominant bound `|λ_2(T_5)| ≤ 1/4` in
`PT.Stochastic.T2T30Normalisation`, the *tight* witness saturating the
bound in `PT.Stochastic.T5CanonicalTight`, and the four-eigenpair
partition (two of modulus `1`, two of modulus `≤ 1/4`) in
`PT.Stochastic.T2T3T5KroneckerSpectrum`. This file pulls those pieces
together into a **single named bridge theorem** that exhibits

  `|λ_2^{eff}(T_30^{tight})|  =  α_cons  =  s²  =  1/4`

without any free `T5Like` parameter and without any `sorry`.

## Why this is the bridge

The pre-existing `T2_alpha_eq_one_quarter` only proves the scalar identity
`(1/2) · (1/2) = 1/4`. Per se this carries no spectral content: it is a
numerical equation between rationals. The bridge theorem here promotes
that scalar identity into a *spectral* statement about the CRT-factorised
transfer matrix `T_30 = T_2 ⊗ T_3 ⊗ T_5^{tight}`:

* `T_2` contributes the trivial Perron eigenvalue `+1` (since
  `(ℤ/2ℤ)* = {1}`);
* `T_3 = !![0,1;1,0]` contributes the involution spectrum `{+1, -1}`;
* `T_5^{tight} = (3/16) · J + (1/4) · I` contributes the spectrum
  `{+1, 1/4, 1/4, 1/4}` (rank-one Perron + thrice-degenerate
  subdominant `1/4`).

Combining via the Kronecker product, the spectrum of `T_30^{tight}` is
the set of all products `λ(T_2) · λ(T_3) · λ(T_5^{tight})`, and on the
Perron-symmetric sector of `T_3` (the one relevant for the *convergence*
of `T_30^n`) the second-largest eigenvalue is `1 · 1 · 1/4 = 1/4 = s²`.
This is the *spectral* meaning of `α_cons`.

## Main results

* `T2_spectral_bridge` — the headline:
  `|λ_2^{eff}(T_30^{tight})| = α_cons = s² = 1/4`, with explicit witness.
* `T2_spectral_bridge_eigenpair` — the witnessing eigenpair:
  `T_30^{tight} · u_2 = (s²) · u_2`, where `u_2 = v_+(T_2) ⊗ v_+(T_3) ⊗
  (1, -1, 0, 0)`.
* `T2_spectral_bridge_full` — packaged summary: Perron eigenpair,
  subdominant eigenpair, equality `|λ_2| = s²`, equality `s² = α_cons`,
  in one statement.

## Reference

* Monograph Chapter 3, §"Exposant de conservation spectrale" and
  Remark `rem:T2_trace`.
* `LEAN_MONOGRAPHIE_GAPS.md`, gap C.1.

## Status

`[THM]` — kernel-verified, no `sorry`, no free parameters. The
numerical input `|λ_2(T_5)| = 1/4` is supplied unconditionally by the
tight canonical instance `T5_tight_T5Like` of
`PT.Stochastic.T5CanonicalTight`.
-/

namespace PT.Conservation

open PT.Stochastic

/-! ### The spectral bridge

The bridge has three layers:

  layer 1 (algebraic): `α_cons := s · s`, so `α_cons = 1/4`.
  layer 2 (CRT structure): `T_30 = T_2 ⊗ T_3 ⊗ T_5` via Kronecker.
  layer 3 (spectral): `|λ_2^{eff}(T_30)| = |λ_2(T_5)| = 1/4`.

The first layer is the file `T2Alpha.lean`; layers 2 and 3 are the
files `T2T30Normalisation.lean` and `T5CanonicalTight.lean`. This
section ties them. -/

/-- **Tight canonical `T_30`.** Notation for the three-factor Kronecker
    product built from the tight canonical `T_5`. Its spectrum on the
    Perron-symmetric sector of `T_3` is `{1, 1/4}`, saturating the PT
    structural bound. -/
noncomputable abbrev T30_tight :
    Matrix ((Fin 1 × Fin 2) × Fin 4) ((Fin 1 × Fin 2) × Fin 4) ℝ :=
  T30 T5_tight_T5Like

/-- **The convergence-controlling eigenvector of `T_30^{tight}`.**
    `u_2 := v_+(T_2) ⊗ v_+(T_3) ⊗ (1, -1, 0, 0)`. -/
noncomputable abbrev T30_tight_lambda_eff_vec :
    (Fin 1 × Fin 2) × Fin 4 → ℝ :=
  T30_lambda_eff_vec T5_tight_T5Like

/-- **Spectral bridge — eigenpair witness.** The vector
    `u_2 = v_+(T_2) ⊗ v_+(T_3) ⊗ (1,-1,0,0)` is an eigenvector of the
    tight three-factor Kronecker product `T_30^{tight}` with eigenvalue
    exactly `s² = 1/4`:

    `T_30^{tight} · u_2  =  s² · u_2`.

    This is the *spectral* incarnation of the scalar identity
    `α_cons = s · s = 1/4` from `T2Alpha`. -/
theorem T2_spectral_bridge_eigenpair :
    T30_tight.mulVec T30_tight_lambda_eff_vec
      = (PT.Stochastic.s ^ 2) • T30_tight_lambda_eff_vec := by
  have h := T30_lambda_eff_eigen' T5_tight_T5Like
  -- `h : T_30 · u_2 = T5_tight_subEigenvalue • u_2`
  -- `T5_tight_subEigenvalue = 1/4 = s²`
  have hs : T5_tight_T5Like.subEigenvalue = PT.Stochastic.s ^ 2 := by
    simp [T5_tight_subEigenvalue_eq, PT.Stochastic.s_def]
    norm_num
  rw [hs] at h
  exact h

/-- **Spectral bridge — absolute value identity.** The convergence-
    controlling eigenvalue of the tight canonical `T_30` has absolute
    value exactly `s² = 1/4`. This **saturates** the structural Kronecker
    bound `|λ_2^{eff}(T_30)| ≤ s²` (proved unconditionally in
    `T30_lambda_eff_bound_s_sq`) on the explicit tight witness. -/
theorem T2_spectral_bridge_abs :
    |T5_tight_T5Like.subEigenvalue| = PT.Stochastic.s ^ 2 := by
  rw [T30_tight_lambda_eff_abs_eq_one_quarter,
      ← PT.Stochastic.one_quarter_eq_s_sq]

/-- **Spectral bridge — alignment with `α_cons`.** The convergence-
    controlling eigenvalue of the tight canonical `T_30` has absolute
    value exactly equal to the PT conservation exponent
    `α_cons = s² = 1/4` of `T2Alpha`:

    `|λ_2^{eff}(T_30^{tight})|  =  α_cons`. -/
theorem T2_spectral_bridge_alpha_cons :
    |T5_tight_T5Like.subEigenvalue| = alpha_cons := by
  rw [T2_spectral_bridge_abs, ← T2_alpha_eq_s_sq]

/-- **Headline of gap C.1 — the spectral bridge.**

    On the Perron-symmetric sector of `T_3`, the convergence-controlling
    eigenvalue `λ_2^{eff}` of the CRT-factorised transfer matrix
    `T_30^{tight} = T_2 ⊗ T_3 ⊗ T_5^{tight}` has absolute value exactly
    equal to the **conservation exponent** of T2:

    `|λ_2^{eff}(T_30^{tight})|  =  α_cons  =  s²  =  1/4`.

    The chain of identities reads, left to right:

    * the absolute value equals `α_cons` — by the bridge above;
    * `α_cons = s²` — definition of `α_cons` (`T2_alpha_eq_s_sq`);
    * `s² = 1/4` — `s = 1/2` (`SHalf.s_def`).

    Compared with the pre-existing `T2_alpha_eq_one_quarter`, which is
    just the scalar identity `(1/2)·(1/2) = 1/4`, this theorem now ties
    `α_cons` to the **spectral** content of `T_30` via the Kronecker
    factorisation — closing gap C.1 of `LEAN_MONOGRAPHIE_GAPS.md`. -/
theorem T2_spectral_bridge :
    |T5_tight_T5Like.subEigenvalue| = alpha_cons
    ∧ alpha_cons = PT.Stochastic.s ^ 2
    ∧ PT.Stochastic.s ^ 2 = (1 : ℝ) / 4 :=
  ⟨T2_spectral_bridge_alpha_cons,
   T2_alpha_eq_s_sq,
   by rw [PT.Stochastic.s_def]; norm_num⟩

/-! ### Bundled summary -/

/-- **Spectral bridge — bundled summary.**

    For the CRT-factorised three-factor Kronecker product
    `T_30^{tight} = T_2 ⊗ T_3 ⊗ T_5^{tight}` we package, in a single
    statement:

    1. The Perron eigenpair `T_30^{tight} · u_+ = u_+`, with the Perron
       eigenvector strictly positive.
    2. The convergence-controlling eigenpair
       `T_30^{tight} · u_2 = s² · u_2`,
       with eigenvalue exactly `s²` (saturating the structural bound).
    3. The headline equality `|λ_2^{eff}(T_30^{tight})| = s² = α_cons
       = 1/4`.

    This is the bridge from the scalar identity `α_cons = s · s = 1/4`
    of `T2Alpha` to the *spectral* identity for `T_30`. -/
theorem T2_spectral_bridge_full :
    -- 1. Perron eigenpair (positive)
    (T30_tight.mulVec (T30_perronVec T5_tight_T5Like)
        = (1 : ℝ) • T30_perronVec T5_tight_T5Like)
    ∧ (∀ ijk, 0 < T30_perronVec T5_tight_T5Like ijk)
    -- 2. Subdominant eigenpair with eigenvalue s²
    ∧ T30_tight.mulVec T30_tight_lambda_eff_vec
        = (PT.Stochastic.s ^ 2) • T30_tight_lambda_eff_vec
    -- 3. Headline equalities |λ_2| = s² = α_cons = 1/4
    ∧ |T5_tight_T5Like.subEigenvalue| = PT.Stochastic.s ^ 2
    ∧ PT.Stochastic.s ^ 2 = alpha_cons
    ∧ alpha_cons = (1 : ℝ) / 4 :=
  ⟨T30_perron_eigen T5_tight_T5Like,
   T30_perronVec_pos T5_tight_T5Like,
   T2_spectral_bridge_eigenpair,
   T2_spectral_bridge_abs,
   T2_alpha_eq_s_sq.symm,
   T2_alpha_eq_one_quarter⟩

end PT.Conservation
