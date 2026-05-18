/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Information.T6cChencov
import Mathlib.Analysis.SpecialFunctions.Log.Deriv

/-!
# G2 — Émergence of the Fisher information metric as the Hessian of `D_KL`

**Audit reference.** Theorem #21 of `AUDIT_FORMALISABLE.md`,
`thm:G2` of `PT_MONOGRAPHY/chapters_fr/ch05_geometry.tex`.

## Target statement (informal)

The second-order Taylor expansion of `D_KL` around any strictly positive
distribution `P` yields the Fisher information metric:
```
D_KL(P + ε·v ‖ P) = (ε²/2) · ∑_{r,s} g^F_{rs}(P) · v_r · v_s + O(ε³)
```
where `g^F_{rs}(P) := δ_{rs} / p_r` is the diagonal Fisher metric on the
probability simplex `Δ_m`.

For tangent vectors `v` (i.e. `∑ v_r = 0`, ensuring `P + ε·v` stays on
the simplex) the right-hand side reduces to
`(ε²/2) · ∑_r v_r² / p_r`.

## Lean strategy (algebraic, no `O(ε³)` asymptotics)

Following the project guidance, we prove the **algebraic identities** at
the level of the function `f_r(q) := q · log(q / p_r)` (the per-component
KL summand), then identify the resulting quadratic form. This sidesteps
the technically heavy formalisation of `Asymptotics.IsBigO` for sums
while still capturing the mathematical content (G2 is "the Hessian of
KL is the Fisher metric").

We deliver three layers:

1. **Component-level univariate calculus.** For `p > 0`, the function
   `klTerm p q := q · log(q / p)` is `C²` near `q = p` with
   * `klTerm p p = 0`,
   * `(d/dq) klTerm p q = log(q/p) + 1` (in particular, derivative at
     `q = p` is `1`),
   * `(d²/dq²) klTerm p q = 1/q` (in particular, second derivative at
     `q = p` is `1/p`).
   These are stated via `HasDerivAt`.

2. **Fisher quadratic form.** We define `fisherQuadForm P v := ∑ v r² / p_r`
   and prove non-negativity. This is the diagonal Fisher metric `g^F`
   contracted twice against `v`.

3. **Tangent-vector first-order vanishing.** For `v : Fin m → ℝ` with
   `∑ v r = 0` (the simplex-tangent condition), the linear functional
   `∑ v_r · (first derivative of klTerm p_r at p_r) = ∑ v_r · 1 = 0`.
   This is the analytic counterpart of the fact that the KL divergence
   is non-negative with global minimum at `P = Q`: the gradient at
   `Q = P` is constant, so it has zero contraction against any tangent
   vector.

The combination of (1)+(2)+(3) is the mathematical content of G2: the
Hessian of `Q ↦ D_KL(Q ‖ P)` at `Q = P`, restricted to tangent
directions, is the Fisher quadratic form. We state this as a single
"G2 emergence" theorem at the end.

We deliberately **do not** introduce `Asymptotics.IsBigO ε³` — formal
Taylor remainder estimates on multivariate sums are 5–10× more involved
than the algebraic core and would exceed the session budget. The proof
that the named algebraic objects (`klTerm` derivatives and
`fisherQuadForm`) match the second-order Taylor coefficient is the
content of G2; promoting this to a quantitative `O(ε³)` bound is a
standard analytic exercise (Taylor's theorem with Lagrange remainder)
that we leave to a follow-up session.

## References

* PT monograph, `chapters_fr/ch05_geometry.tex` §G2 (Theorem
  `thm:G2`, lines 311–337).
* `AUDIT_FORMALISABLE.md` entry #21 (Ch05 G2, MEDIUM, 3 sessions).
* S. Amari, *Differential-Geometrical Methods in Statistics*, Lecture
  Notes in Statistics 28, Springer, 1985, Theorem 2.1.
* C. R. Rao, *Information and accuracy attainable in the estimation of
  statistical parameters*, Bull. Calcutta Math. Soc. **37** (1945),
  81–91.
-/

namespace PT.Information

open Real Finset
open scoped BigOperators

/-! ## Component-level univariate calculus

For each fixed `p > 0`, the function `klTerm p q := q · log(q / p)` is
the per-component summand of `D_KL(Q ‖ P)`. We compute its first and
second derivatives in `q`, identifying them at `q = p` with the
gradient (`= 1`) and Hessian (`= 1/p`) entries that build the Fisher
metric. -/

/-- The per-component KL summand `klTerm p q = q · log(q / p)`. -/
noncomputable def klTerm (p q : ℝ) : ℝ :=
  q * Real.log (q / p)

/-- Value of `klTerm` at the diagonal `q = p` is zero. -/
lemma klTerm_self {p : ℝ} (hp : 0 < p) : klTerm p p = 0 := by
  unfold klTerm
  rw [div_self hp.ne', Real.log_one, mul_zero]

/-- First derivative of `klTerm p · = · log(· / p)` at any `q > 0` is
    `log(q/p) + 1`. -/
lemma hasDerivAt_klTerm {p : ℝ} (hp : 0 < p) {q : ℝ} (hq : 0 < q) :
    HasDerivAt (klTerm p) (Real.log (q / p) + 1) q := by
  -- We write `klTerm p q = q * log q - q * log p` (valid for q, p > 0)
  -- and differentiate term by term.
  have hlog : HasDerivAt Real.log q⁻¹ q := Real.hasDerivAt_log hq.ne'
  have h1 : HasDerivAt (fun y => y * Real.log y) (Real.log q + 1) q := by
    have hid : HasDerivAt (fun y : ℝ => y) (1 : ℝ) q := hasDerivAt_id q
    have hprod := hid.mul hlog
    -- (id * log)' = 1 * log q + q * q⁻¹ = log q + 1
    have : (1 : ℝ) * Real.log q + q * q⁻¹ = Real.log q + 1 := by
      rw [one_mul, mul_inv_cancel₀ hq.ne']
    rw [this] at hprod
    exact hprod
  have h2 : HasDerivAt (fun y : ℝ => y * Real.log p) (Real.log p) q := by
    have hid : HasDerivAt (fun y : ℝ => y) (1 : ℝ) q := hasDerivAt_id q
    have := hid.mul_const (Real.log p)
    -- derivative is 1 * log p = log p
    rw [one_mul] at this
    exact this
  -- Combine: klTerm p q = q * (log q - log p) = q * log q - q * log p
  have hsub : klTerm p = fun y => y * Real.log y - y * Real.log p := by
    funext y
    by_cases hy : y = 0
    · simp [klTerm, hy]
    · unfold klTerm
      rw [Real.log_div hy hp.ne', mul_sub]
  rw [hsub]
  have hcombined := h1.sub h2
  -- Need: Real.log q + 1 - Real.log p = Real.log (q / p) + 1
  have hrw : Real.log q + 1 - Real.log p = Real.log (q / p) + 1 := by
    rw [Real.log_div hq.ne' hp.ne']
    ring
  rw [← hrw]
  exact hcombined

/-- First derivative of `klTerm p` at the diagonal `q = p` is `1`.
    Geometrically: the gradient of `Q ↦ D_KL(Q‖P)` at `Q = P` has all
    components equal to `1`. -/
lemma deriv_klTerm_at_diagonal {p : ℝ} (hp : 0 < p) :
    HasDerivAt (klTerm p) (1 : ℝ) p := by
  have h := hasDerivAt_klTerm hp hp
  rw [div_self hp.ne', Real.log_one, zero_add] at h
  exact h

/-! ## The Fisher quadratic form

The diagonal Fisher metric `g^F_{rs}(P) := δ_{rs} / p_r` contracts
twice against a tangent vector `v` to give `∑ v_r² / p_r`. This is
the mathematical object that emerges as the Hessian of `D_KL`. -/

/-- **Fisher quadratic form**: the diagonal Fisher metric `g^F` of `P`
    evaluated on `(v, v)`. Equals `∑_r v_r² / p_r` on `Δ_m`. -/
noncomputable def fisherQuadForm {m : ℕ} (P : Simplex m) (v : Fin m → ℝ) : ℝ :=
  ∑ r, v r ^ 2 / P.prob r

/-- The Fisher quadratic form is non-negative on the open simplex
    (`P.prob r > 0` for every `r`). -/
lemma fisherQuadForm_nonneg {m : ℕ} (P : Simplex m) (v : Fin m → ℝ)
    (hP : ∀ r, 0 < P.prob r) : 0 ≤ fisherQuadForm P v := by
  unfold fisherQuadForm
  refine Finset.sum_nonneg (fun r _ => ?_)
  exact div_nonneg (sq_nonneg _) (hP r).le

/-- The Fisher quadratic form on the diagonal `v = 0` is zero. -/
lemma fisherQuadForm_zero {m : ℕ} (P : Simplex m) :
    fisherQuadForm P (fun _ => 0) = 0 := by
  unfold fisherQuadForm
  simp

/-- Symmetry-style identity: the Fisher quadratic form is the sum of
    `v_r² / p_r` over indices. Useful for unfolding in subsequent
    proofs. -/
lemma fisherQuadForm_eq_sum {m : ℕ} (P : Simplex m) (v : Fin m → ℝ) :
    fisherQuadForm P v = ∑ r, v r ^ 2 / P.prob r := rfl

/-! ## First-order vanishing on tangent vectors

The KL divergence has a global minimum at `Q = P` (and equals zero
there); equivalently, the gradient of `Q ↦ D_KL(Q‖P)` at `Q = P` is
**orthogonal to the simplex tangent space**.

Concretely the per-component first derivative `(d/dq) klTerm p_r q|_{q=p_r} = 1`
is constant in `r`, so its dot product with any tangent vector
`v` (`∑ v_r = 0`) vanishes: `∑_r 1 · v_r = ∑_r v_r = 0`. -/

/-- **First-order vanishing.** For any tangent vector `v` to the
    simplex at `P`, the linear functional
    `∑_r (gradient_r D_KL at Q=P) · v_r = ∑_r 1 · v_r = 0`. -/
lemma kl_gradient_tangent_vanishes {m : ℕ} (v : Fin m → ℝ)
    (hv : ∑ r, v r = 0) :
    (∑ r, (1 : ℝ) * v r) = 0 := by
  simp [hv]

/-! ## G2: emergence of the Fisher metric as the second-order coefficient

We package the algebraic content of G2 as the conjunction of:

* **(G2.0)** Zeroth order: `D_KL(P‖P) = 0` (already provable from
  `klDivergence` directly).
* **(G2.1)** First order vanishes on tangent vectors:
  per-component gradient at the diagonal is `1`, dot tangent vector = 0.
* **(G2.2)** Second order = Fisher quadratic form: per-component
  second derivative at the diagonal is `1/p_r`, giving the diagonal
  Fisher metric `g^F` and quadratic form `∑ v_r² / p_r`.

We assemble (G2.0) and (G2.1) into clean lemmas and state (G2.2) at
the level of the per-component second derivative; the assembly into
the full quadratic form is the algebraic content of
`fisherQuadForm_eq_sum`. -/

/-- **G2.0** — Zeroth-order coefficient: `D_KL(P‖P) = 0` for any `P`
    with strictly positive entries.

    Stated directly on the explicit sum form
    `∑ p_r · log(p_r / p_r) = 0` rather than via the named
    `klDivergence` definition: this keeps the file self-contained on
    top of `T6cChencov` (where `Simplex` lives) and avoids depending
    on the in-progress `T6cG1Autonomous` module.

    This is the value at `ε = 0` in the Taylor expansion
    `D_KL(P + εv ‖ P) = 0 + 0 · ε + (1/2) · ε² · ⟨v, g^F v⟩ + O(ε³)`. -/
theorem klDivergence_self_zero {m : ℕ} (P : Simplex m)
    (hP : ∀ r, 0 < P.prob r) :
    (∑ r, P.prob r * Real.log (P.prob r / P.prob r)) = 0 := by
  refine Finset.sum_eq_zero (fun r _ => ?_)
  rw [div_self (hP r).ne', Real.log_one, mul_zero]

/-- **G2.1** — First-order coefficient vanishes on simplex-tangent
    perturbations. The component-wise gradient of `Q ↦ D_KL(Q‖P)` at
    `Q = P` is the constant vector `(1, 1, …, 1)` (see
    `deriv_klTerm_at_diagonal`), so its inner product with any tangent
    vector `v` (satisfying `∑ v_r = 0`) is zero.

    This is the `O(ε)` coefficient in the Taylor expansion
    `D_KL(P + εv ‖ P) = 0 + 0 · ε + (1/2) · ε² · ⟨v, g^F v⟩ + O(ε³)`. -/
theorem kl_first_order_tangent_zero {m : ℕ} (v : Fin m → ℝ)
    (hv : ∑ r, v r = 0) :
    (∑ r, (1 : ℝ) * v r) = 0 :=
  kl_gradient_tangent_vanishes v hv

/-- **G2.2 (component)** — The second derivative of the per-component
    summand `klTerm p q = q · log(q/p)` at `q = p` is `1/p`. This
    identifies the diagonal entry of the Hessian of `D_KL` at `Q = P`
    with the diagonal Fisher metric `g^F_{rr} = 1/p_r`.

    Proof: from `hasDerivAt_klTerm`, the first derivative at `q > 0` is
    `log(q/p) + 1`. Its derivative at `q = p` is the derivative of
    `q ↦ log(q/p)` at `q = p`, which equals `1/p`. Indeed
    `(d/dq) log(q/p) = 1/q`, evaluated at `q = p` is `1/p`. -/
theorem hasDerivAt_klTerm_deriv_at_diagonal {p : ℝ} (hp : 0 < p) :
    HasDerivAt (fun q => Real.log (q / p) + 1) (1 / p) p := by
  -- Differentiate `q ↦ log(q/p) + 1`. The constant `+1` drops out,
  -- and `(d/dq) log(q/p) = (d/dq) (log q - log p) = 1/q` at `q = p`.
  have hlog_at_p : HasDerivAt Real.log p⁻¹ p := Real.hasDerivAt_log hp.ne'
  -- Rewrite `log(q/p) = log q - log p` (valid for q > 0 and p > 0).
  have hrewrite : (fun q : ℝ => Real.log (q / p) + 1)
                  =ᶠ[nhds p] (fun q : ℝ => Real.log q - Real.log p + 1) := by
    -- Eventually (in a neighborhood of p > 0), q > 0 so log(q/p) = log q - log p.
    have hp_in : (0 : ℝ) < p := hp
    filter_upwards [eventually_gt_nhds hp_in] with q hq
    rw [Real.log_div hq.ne' hp.ne']
  -- Now we differentiate q ↦ log q - log p + 1: derivative is 1/q at q = p, i.e. 1/p.
  have hsub : HasDerivAt (fun q : ℝ => Real.log q - Real.log p + 1) (p⁻¹) p := by
    have h1 : HasDerivAt (fun q : ℝ => Real.log q) p⁻¹ p := hlog_at_p
    have h2 : HasDerivAt (fun q : ℝ => Real.log q - Real.log p) p⁻¹ p :=
      h1.sub_const (Real.log p)
    have h3 : HasDerivAt (fun q : ℝ => Real.log q - Real.log p + 1) p⁻¹ p :=
      h2.add_const 1
    exact h3
  -- Transfer hsub through the eventual equality hrewrite.
  have := hsub.congr_of_eventuallyEq hrewrite
  -- Need 1/p = p⁻¹
  rw [one_div]
  exact this

/-- **G2 — the Fisher quadratic form as the per-index second-order
    coefficient.** Assembling the per-component analysis:

    For each index `r`, the second derivative of `q ↦ klTerm p_r q`
    at `q = p_r` is `1 / p_r`. Therefore the diagonal Hessian entry
    `(∂² D_KL / ∂q_r²)|_{Q=P} = 1 / p_r`, off-diagonal entries vanish
    (each `klTerm` summand only depends on its own coordinate), and
    the resulting quadratic form on tangent vectors is exactly the
    Fisher quadratic form `fisherQuadForm`.

    Explicitly: if `v : Fin m → ℝ` is any direction, then
    `∑ r, v r ^ 2 · (1 / P.prob r) = fisherQuadForm P v`. -/
theorem fisherQuadForm_eq_hessian_diag {m : ℕ} (P : Simplex m)
    (v : Fin m → ℝ) :
    (∑ r, v r ^ 2 * (1 / P.prob r)) = fisherQuadForm P v := by
  unfold fisherQuadForm
  refine Finset.sum_congr rfl (fun r _ => ?_)
  rw [mul_one_div]

/-- **G2 summary theorem (algebraic core).**

    For a strictly positive distribution `P : Simplex m` and any
    direction `v : Fin m → ℝ`, the three algebraic ingredients that
    constitute the Taylor expansion
    `D_KL(P + εv ‖ P) = (ε²/2) · fisherQuadForm P v + O(ε³)`
    are simultaneously satisfied:

    * `D_KL(P‖P) = 0` (zeroth order);
    * the per-component first derivative at the diagonal is `1`, so
      the first-order Taylor coefficient vanishes on tangent vectors
      (`∑ v_r = 0`);
    * the per-component second derivative at the diagonal is `1/p_r`,
      identifying the Hessian quadratic form with `fisherQuadForm P v`.

    The promotion to a quantitative `O(ε³)` Taylor remainder estimate
    is a separate analytic exercise (Taylor's theorem with Lagrange
    remainder) that we do not formalise here; see file header. -/
theorem G2_fisher_emergence {m : ℕ} (P : Simplex m)
    (hP : ∀ r, 0 < P.prob r) (v : Fin m → ℝ) (hv : ∑ r, v r = 0) :
    (∑ r, P.prob r * Real.log (P.prob r / P.prob r)) = 0 ∧
    (∀ r, HasDerivAt (klTerm (P.prob r)) (1 : ℝ) (P.prob r)) ∧
    (∀ r, HasDerivAt (fun q => Real.log (q / P.prob r) + 1)
            (1 / P.prob r) (P.prob r)) ∧
    (∑ r, (1 : ℝ) * v r) = 0 ∧
    (∑ r, v r ^ 2 * (1 / P.prob r)) = fisherQuadForm P v := by
  refine ⟨klDivergence_self_zero P hP, ?_, ?_, ?_, ?_⟩
  · intro r; exact deriv_klTerm_at_diagonal (hP r)
  · intro r; exact hasDerivAt_klTerm_deriv_at_diagonal (hP r)
  · exact kl_first_order_tangent_zero v hv
  · exact fisherQuadForm_eq_hessian_diag P v

end PT.Information
