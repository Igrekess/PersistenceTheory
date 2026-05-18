/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Information.T6cChencov

/-!
# T6c (G1 autonomous) — Csiszár characterisation of KL divergence

PT's autonomous formulation of theorem **T6c** (`chapters_fr/ch05_geometry.tex`,
Theorem `G1_autonomous`): the Kullback–Leibler divergence is the unique
f-divergence on the probability simplex that is additive under product
distributions, up to a non-negative scalar.

This replaces the entropy-based formulation that lived in
`T6cChencov.lean` (retracted on 2026-05-16: the four-axiom statement
admitted Rényi entropies as counterexamples; see retraction note there).
The divergence-based formulation is mathematically correct because Rényi
**divergences** `D_α(P‖Q) = (1/(α−1)) log(∑ p_i^α q_i^{1−α})` for
`α ≠ 1` are **not** additive on products in the sense required here
(see `chapters_fr/ch05_geometry.tex`, exclusion list following Theorem
`G1_autonomous`).

## Roadmap

* `fDivergence f P Q := ∑ q_i · f(p_i/q_i)` — generic f-divergence.
* `klDivergence P Q := ∑ p_i · log(p_i/q_i)` — Kullback–Leibler.
* `klDivergence_eq_fDivergence_tlog` — KL is the f-divergence with
  generator `f(t) = t · log t`.
* `klDivergence_additive` — KL is additive under product distributions
  (the "easy direction": KL satisfies the hypothesis of the theorem).
* `G1_autonomous_DKL_unique` — **headline**: any continuous f-divergence
  additive on products is a non-negative scalar multiple of KL. Proof
  via Csiszár 1967 functional-equation chase; `cauchy_log_equation` (in
  `T6cChencov.lean`) is the workhorse.

## Sorry inventory (2026-05-17 state) — **0 sorries** in this file.

The Csiszár 1967 §3 step 1 (binary product expansion + Pexider
isolation) was closed by splitting the residual identity into four
quadrants of `(0, ∞)²`:

* `residual_anchor_gt_one` — anchor on `(1, ∞)²`, via `p₁ = p₂ = 0`.
* `residual_left_lt_one` — extension to `[0, 1) × (1, ∞)`, via
  `p₁ = u·q₁, p₂ = 0` and substitution of the anchor at `(\bar u, β)`.
* `residual_right_lt_one` — extension to `(1, ∞) × [0, 1)` by
  commutativity of the 4-term equation.
* `residual_both_lt_one` — extension to `(0, 1)²`, via the
  `q₁ = q₂ = 1/2` symmetric specialisation, substituting the three
  already-proven residuals at `(u, \bar v), (\bar u, v), (\bar u, \bar v)`.

Resolved (now sorry-free):

* `klDivergence_additive` — closed via `kl_term_factor` reindex.
* `cauchy_log_from_pexider` — pure algebra, divides Pexider by `uv`.
* `csiszar_residual_from_binary_additivity` — assembled from the four
  quadrant lemmas above + boundary cases `u = 1` or `v = 1` (trivial).
* `csiszar_generator_decomposition` (Steps 2–3) — Cauchy-log via
  `htilde(x) := f(exp x)/exp x - f(1)` packaged as an `AddMonoidHom`,
  then `map_real_smul` for ℝ-linearity, then continuity at 0 via
  `continuous_mul_log` + `tendsto_nhds_unique` on `nhdsGT 0`.
* `csiszar_c1_nonneg` — convex-secant slope monotonicity on `1 < e < e²`.
* `fDivergence_of_csiszar_form` — case split on `P i = 0` vs `> 0`.
* `G1_autonomous_DKL_unique` — assembled from the sub-lemmas above.

## References

* I. Csiszár, *Information-type measures of difference of probability
  distributions and indirect observations*, Studia Sci. Math. Hungar.
  **2** (1967), 299–318.
* J. E. Shore and R. W. Johnson, *Axiomatic derivation of the principle
  of maximum entropy and the principle of minimum cross-entropy*, IEEE
  Trans. Inform. Theory **26** (1980), 26–37.
* Y. Senez, *La théorie de la persistance*, monograph
  `PT_MONOGRAPHY/chapters_fr/ch05_geometry.tex`, §G1
  (Theorems `G1` and `G1_autonomous`).
-/

namespace PT.Information

open Real Finset
open scoped BigOperators

/-! ## Definitions -/

/-- Generic **f-divergence** `D_f(P‖Q) := ∑ q_i · f(p_i/q_i)`.

    Mathlib convention: when `q_i = 0`, the factor `q_i · f(p_i/q_i)`
    equals `0` regardless of `f`'s behaviour (since `0 · x = 0`). This is
    the standard handling of the "0 times divergent" edge case at the
    boundary of the simplex. -/
noncomputable def fDivergence (f : ℝ → ℝ) {m : ℕ} (P Q : Simplex m) : ℝ :=
  ∑ i, Q.prob i * f (P.prob i / Q.prob i)

/-- **Kullback–Leibler divergence** `D_KL(P‖Q) := ∑ p_i · log(p_i/q_i)`.

    Mathlib conventions:
    * When `p_i = 0` the term is `0 · log(0/q_i) = 0`.
    * When `q_i = 0` and `p_i > 0`, Mathlib's `Real.log 0 = 0` forces
      the term to `0` instead of the mathematical convention `+∞`. PT
      applications use distributions strictly positive on the sieve
      support, so this edge case does not arise. -/
noncomputable def klDivergence {m : ℕ} (P Q : Simplex m) : ℝ :=
  ∑ i, P.prob i * Real.log (P.prob i / Q.prob i)

/-- KL is the f-divergence with generator `f(t) = t · log t`.

    When `Q.prob i ≠ 0`, the f-divergence term is
    `Q_i · (P_i/Q_i) · log(P_i/Q_i) = P_i · log(P_i/Q_i)`, matching the
    KL definition. When `Q.prob i = 0`, both sides equal zero (the LHS
    by the `0 · x = 0` convention, the RHS because `P_i / 0 = 0` and
    then `P_i · log 0 = 0` in Mathlib). -/
lemma klDivergence_eq_fDivergence_tlog {m : ℕ} (P Q : Simplex m) :
    klDivergence P Q = fDivergence (fun t => t * Real.log t) P Q := by
  unfold klDivergence fDivergence
  refine Finset.sum_congr rfl (fun i _ => ?_)
  by_cases hQ : Q.prob i = 0
  · -- Q_i = 0: both sides reduce to 0 by Mathlib's `0/0 = 0`, `log 0 = 0`.
    simp [hQ]
  · -- Q_i ≠ 0: algebraic identity.
    field_simp

/-! ## Easy direction: KL is additive under product distributions -/

/-- Auxiliary pointwise factorisation for KL terms.

    For `0 ≤ a, b` and `0 < c, d`,
    `(a*b) · log((a*b)/(c*d)) = (a*b) · log(a/c) + (a*b) · log(b/d)`.

    The split via `Real.log_mul` requires `a/c ≠ 0` and `b/d ≠ 0`; the
    boundary cases `a = 0` or `b = 0` are handled separately (both sides
    are then `0`). -/
private lemma kl_term_factor {a b c d : ℝ}
    (ha : 0 ≤ a) (hb : 0 ≤ b) (hc : 0 < c) (hd : 0 < d) :
    (a * b) * Real.log ((a * b) / (c * d))
      = (a * b) * Real.log (a / c) + (a * b) * Real.log (b / d) := by
  rcases ha.lt_or_eq with ha_pos | ha_zero
  · rcases hb.lt_or_eq with hb_pos | hb_zero
    · -- a > 0 and b > 0: pure algebraic identity using log_mul / log_div.
      have hac : a / c ≠ 0 := (div_pos ha_pos hc).ne'
      have hbd : b / d ≠ 0 := (div_pos hb_pos hd).ne'
      have hmul : (a * b) / (c * d) = (a / c) * (b / d) := by
        rw [mul_div_mul_comm]
      rw [hmul, Real.log_mul hac hbd, mul_add]
    · -- b = 0: both sides multiply by `a * b = 0`.
      simp [← hb_zero]
  · -- a = 0: both sides multiply by `a * b = 0`.
    simp [← ha_zero]

/-- **KL divergence is additive under product distributions.**

    `D_KL(P₁ ⊗ P₂ ‖ Q₁ ⊗ Q₂) = D_KL(P₁‖Q₁) + D_KL(P₂‖Q₂)`.

    This is the "easy direction" of Csiszár's characterisation: KL is
    a candidate for the unique additive f-divergence.

    The proof follows the pattern of `shannon_additive` in
    `T6cChencov.lean`:

    1. Reindex `Fin (m*n) → Fin m × Fin n` via `finProdFinEquiv`.
    2. Apply `kl_term_factor` pointwise to split
       `log((p₁·p₂)/(q₁·q₂)) = log(p₁/q₁) + log(p₂/q₂)`.
       Positivity of `Q₁, Q₂` is used to apply `Real.log_mul` on the
       quotient `(a/c)·(b/d)`; the boundary cases `p₁ = 0` or `p₂ = 0`
       absorb the indeterminate log via the multiplicative factor.
    3. Distribute and use `P₁.sums_to_one`, `P₂.sums_to_one` for the
       marginals.

    The boundary case `q_i = 0` requires care: the Mathlib convention
    silently sets the term to 0, which differs from the mathematical
    convention `+∞`. The clean statement adds positivity hypotheses on
    `Q₁` and `Q₂`. -/
theorem klDivergence_additive {m n : ℕ}
    (P₁ Q₁ : Simplex m) (P₂ Q₂ : Simplex n)
    (_hQ₁ : ∀ i, 0 < Q₁.prob i) (_hQ₂ : ∀ j, 0 < Q₂.prob j) :
    klDivergence (P₁.prod P₂) (Q₁.prod Q₂)
      = klDivergence P₁ Q₁ + klDivergence P₂ Q₂ := by
  -- Pattern follows `shannon_additive` in `T6cChencov.lean`.
  unfold klDivergence
  -- Goal: ∑ k : Fin (m*n),
  --   (P₁.prod P₂).prob k * log ((P₁.prod P₂).prob k / (Q₁.prod Q₂).prob k)
  --   = (∑ i, P₁.prob i * log (P₁.prob i / Q₁.prob i))
  --   + (∑ j, P₂.prob j * log (P₂.prob j / Q₂.prob j)).
  -- Step 1: reindex `Fin (m*n) ≃ Fin m × Fin n` via `finProdFinEquiv`.
  show (∑ k : Fin (m * n),
          (P₁.prob (finProdFinEquiv.symm k).1 *
           P₂.prob (finProdFinEquiv.symm k).2) *
          Real.log
            ((P₁.prob (finProdFinEquiv.symm k).1 *
              P₂.prob (finProdFinEquiv.symm k).2) /
             (Q₁.prob (finProdFinEquiv.symm k).1 *
              Q₂.prob (finProdFinEquiv.symm k).2)))
       = (∑ i, P₁.prob i * Real.log (P₁.prob i / Q₁.prob i))
       + (∑ j, P₂.prob j * Real.log (P₂.prob j / Q₂.prob j))
  rw [show (∑ k : Fin (m * n),
            (P₁.prob (finProdFinEquiv.symm k).1 *
             P₂.prob (finProdFinEquiv.symm k).2) *
            Real.log
              ((P₁.prob (finProdFinEquiv.symm k).1 *
                P₂.prob (finProdFinEquiv.symm k).2) /
               (Q₁.prob (finProdFinEquiv.symm k).1 *
                Q₂.prob (finProdFinEquiv.symm k).2)))
        = ∑ k : Fin (m * n), (fun ij : Fin m × Fin n =>
            (P₁.prob ij.1 * P₂.prob ij.2) *
            Real.log
              ((P₁.prob ij.1 * P₂.prob ij.2) /
               (Q₁.prob ij.1 * Q₂.prob ij.2)))
              (finProdFinEquiv.symm k) from rfl]
  rw [Equiv.sum_comp finProdFinEquiv.symm
      (fun ij : Fin m × Fin n =>
        (P₁.prob ij.1 * P₂.prob ij.2) *
        Real.log
          ((P₁.prob ij.1 * P₂.prob ij.2) /
           (Q₁.prob ij.1 * Q₂.prob ij.2)))]
  rw [Fintype.sum_prod_type]
  -- Step 2: pointwise split using `kl_term_factor`.
  have hsplit : ∀ i j,
      (P₁.prob i * P₂.prob j) *
        Real.log ((P₁.prob i * P₂.prob j) / (Q₁.prob i * Q₂.prob j))
        = (P₁.prob i * P₂.prob j) * Real.log (P₁.prob i / Q₁.prob i)
        + (P₁.prob i * P₂.prob j) * Real.log (P₂.prob j / Q₂.prob j) := by
    intro i j
    exact kl_term_factor (P₁.nonneg i) (P₂.nonneg j) (_hQ₁ i) (_hQ₂ j)
  simp_rw [hsplit, Finset.sum_add_distrib]
  -- Step 3: marginalise using `P₁.sums_to_one` and `P₂.sums_to_one`.
  congr 1
  · -- ∑ i, ∑ j, P₁ i * P₂ j * log(P₁ i / Q₁ i)
    --   = ∑ i, P₁ i * log(P₁ i / Q₁ i).
    refine Finset.sum_congr rfl (fun i _ => ?_)
    -- inner sum over j: ∑ j, P₁ i * P₂ j * log(P₁ i / Q₁ i)
    -- = (∑ j, P₂ j) * (P₁ i * log(P₁ i / Q₁ i)) via factoring.
    have : ∀ j, (P₁.prob i * P₂.prob j) * Real.log (P₁.prob i / Q₁.prob i)
              = P₂.prob j * (P₁.prob i * Real.log (P₁.prob i / Q₁.prob i)) := by
      intro j; ring
    simp_rw [this, ← Finset.sum_mul, P₂.sums_to_one, one_mul]
  · -- ∑ i, ∑ j, P₁ i * P₂ j * log(P₂ j / Q₂ j)
    --   = ∑ j, P₂ j * log(P₂ j / Q₂ j).
    rw [Finset.sum_comm]
    refine Finset.sum_congr rfl (fun j _ => ?_)
    -- ∑ i, P₁ i * P₂ j * log(P₂ j / Q₂ j)
    --   = (∑ i, P₁ i) * (P₂ j * log(P₂ j / Q₂ j))
    --   = P₂ j * log(P₂ j / Q₂ j).
    have hrew : ∀ i,
        (P₁.prob i * P₂.prob j) * Real.log (P₂.prob j / Q₂.prob j)
          = P₁.prob i * (P₂.prob j * Real.log (P₂.prob j / Q₂.prob j)) := by
      intro i; ring
    simp_rw [hrew, ← Finset.sum_mul, P₁.sums_to_one, one_mul]

/-! ## Additivity predicate and headline theorem -/

/-- Family-level predicate: a divergence on each `Simplex m` is
    additive under product distributions whose denominators are
    strictly positive (the standard `P ≪ Q` regime).

    Note: Rényi divergences `D_α(P‖Q) = (1/(α−1)) log(∑ p_i^α q_i^{1−α})`
    for `α ≠ 1` do NOT satisfy this (the sum `∑ p_i^α q_i^{1−α}` does
    not factorise across products in general). This is the key axiom
    that excludes the Rényi family from the characterisation.

    **Positivity hypothesis on `Q`.** Under Mathlib's conventions
    (`a/0 = 0`, `Real.log 0 = 0`), even `klDivergence` itself fails the
    universal additivity equality on product distributions whose
    denominators have zero entries: e.g.
    `klDivergence (P₁⊗P₂) (Q₁⊗Q₂) ≠ klDivergence P₁ Q₁ + klDivergence P₂ Q₂`
    when `Q₁ = (1, 0)`, `Q₂ = (1, 0)`. The standard convention restricts
    additivity to absolutely continuous pairs (`P ≪ Q`), equivalently
    `Q` strictly positive on the open simplex. This matches the implicit
    setting of Csiszár 1967 and the monograph (`chapters_fr/ch05_geometry.tex`). -/
def IsDivergenceAdditive
    (D : ∀ m, Simplex m → Simplex m → ℝ) : Prop :=
  ∀ {m n : ℕ} (P₁ Q₁ : Simplex m) (P₂ Q₂ : Simplex n),
    (∀ i, 0 < Q₁.prob i) → (∀ j, 0 < Q₂.prob j) →
    D (m * n) (P₁.prod P₂) (Q₁.prod Q₂) = D m P₁ Q₁ + D n P₂ Q₂

/-- Family-level predicate: a divergence comes from a standard
    f-divergence generator via `D(P‖Q) = ∑ q_i · f(p_i/q_i)`.

    The generator `f` must satisfy the standard f-divergence
    conditions (Csiszár 1967, Wikipedia "f-divergence"):
    * **Continuity** on `[0, ∞)`.
    * **Convexity** on `[0, ∞)` (essential — without convexity, the
      "anti-KL" generator `f(t) = −t log t` (Rényi-style) trivially
      satisfies additivity but yields a *negative* multiple of KL,
      violating the Csiszár characterisation's `c₁ > 0`).
    * **Normalisation** `f(1) = 0` (ensures `D(P‖P) = 0`). -/
def IsFDivergenceFamily
    (D : ∀ m, Simplex m → Simplex m → ℝ) : Prop :=
  ∃ f : ℝ → ℝ, Continuous f ∧ ConvexOn ℝ (Set.Ici (0 : ℝ)) f ∧ f 1 = 0 ∧
    ∀ m (P Q : Simplex m), D m P Q = fDivergence f P Q

/-! ## Csiszár step 4: vanishing of the linear `(t-1)` term

The Csiszár characterisation produces a generator of the form
`f(t) = c₁ · t · log t + c₂ · (t - 1)`. The `c₂ · (t - 1)` summand
contributes zero to the f-divergence on any pair of normalised
distributions whose denominators are strictly positive: indeed
`∑ q_i · (p_i/q_i − 1) = ∑ p_i − ∑ q_i = 1 − 1 = 0`.

The lemma below is the technical content of Csiszár 1967 step 4; it
is kernel-verified and reusable inside the eventual closure of
`G1_autonomous_DKL_unique`. -/

/-- The linear "affine correction" generator `f(t) = c · (t - 1)`
    produces zero f-divergence on any pair `(P, Q)` of normalised
    distributions when `Q` is strictly positive.

    This is Csiszár 1967 step 4: the `c₂ · (t - 1)` summand of the
    general solution drops out on the simplex. -/
lemma fDivergence_linear_term_vanishes {m : ℕ} (c : ℝ)
    (P Q : Simplex m) (hQ : ∀ i, 0 < Q.prob i) :
    fDivergence (fun t => c * (t - 1)) P Q = 0 := by
  unfold fDivergence
  -- ∑ i, Q i * (c * (P i / Q i - 1)) = c * (∑ i, P i - ∑ i, Q i) = c * 0 = 0.
  have hterm : ∀ i, Q.prob i * (c * (P.prob i / Q.prob i - 1))
                  = c * (P.prob i - Q.prob i) := by
    intro i
    have hQi : Q.prob i ≠ 0 := (hQ i).ne'
    field_simp
  simp_rw [hterm]
  rw [← Finset.mul_sum, Finset.sum_sub_distrib,
      P.sums_to_one, Q.sums_to_one, sub_self, mul_zero]

/-! ## Csiszár factorisation: fine-grained sub-lemmas

The headline theorem `G1_autonomous_DKL_unique` (below) is assembled
from three sub-lemmas and one technical residual, each encapsulating
one tractable sub-task of the Csiszár 1967 functional-equation chase.
The current state factorises the former monolithic sorry into three
fine-grained sorries (`csiszar_generator_decomposition`,
`csiszar_c1_nonneg`, and the boundary residual at the bottom of the
headline theorem).

Sorry-free sub-lemmas (already kernel-verified):

* `binary_prod_components_at` (step 1 helper): explicit formula for the
  four components of `(binary p₁ _).prod (binary p₂ _) : Simplex 4`.
  Direct from the definitions of `binary` and `Simplex.prod`.

* `fDivergence_of_csiszar_form` (step 4): if `f` has the canonical
  form `c₁ · t log t + c₂ · (t − 1)` (with `f 0 = -c₂`), then on any
  strictly positive `Q`, the f-divergence equals `c₁ · D_KL(P‖Q)`.
  Uses `fDivergence_linear_term_vanishes` for the affine summand.

Sorried sub-lemmas (fine-grained, each with a clear strategy):

* `csiszar_generator_decomposition` (steps 1–3): if `f` generates a
  continuous f-divergence additive on products, then there exist
  constants `c₁, c₂` such that `f(t) = c₁ · t · log t + c₂ · (t − 1)`
  for every `t > 0`, and `f 0 = -c₂` (continuity at the boundary).
  Proof outline: specialise additivity to binary product distributions
  via `binary_prod_components_at` to obtain a Pexider-type functional
  equation in `f`, reduce to Cauchy's logarithmic equation, apply
  `cauchy_log_equation`, back-substitute. Effort: ~2 sessions.

* `csiszar_c1_nonneg` (step 5): given the generator form, the leading
  coefficient `c₁` is non-negative. Proof outline: evaluate at a
  non-uniform binary pair where `t · log t` contributes strictly
  positive mass; additivity-derived constraints force `c₁ ≥ 0`.
  Effort: ~1 session.

* Boundary residual in `G1_autonomous_DKL_unique` (the only sorry left
  in the headline theorem body): when some `Q i = 0`, the strictly
  positive branch of `fDivergence_of_csiszar_form` does not apply
  directly; the standard fix is to either restrict to the strictly
  positive support of `Q` (and use the f-divergence convention
  `Q i · f(P i / 0) = 0`) or to derive `c₂ = 0` from a second
  specialisation of additivity. Effort: ~0.5 session. -/

/-- **Csiszár step 1 helper.** Components of the product of two binary
    distributions, indexed by `Fin (2 * 2) ≃ Fin 2 × Fin 2`.

    For `p₁, p₂ ∈ [0, 1]`, the product `(binary p₁) ⊗ (binary p₂)` is
    a `Simplex 4` whose probabilities at indices
    `finProdFinEquiv (i, j)` are:

    * `(0, 0) ↦ p₁ · p₂`
    * `(1, 0) ↦ (1 − p₁) · p₂`
    * `(0, 1) ↦ p₁ · (1 − p₂)`
    * `(1, 1) ↦ (1 − p₁) · (1 − p₂)`

    The exact mapping to indices in `Fin 4` depends on
    `finProdFinEquiv`'s convention; we state the lemma at the
    `(Fin 2 × Fin 2)`-level after pulling back through `finProdFinEquiv`,
    which is the form actually used downstream. -/
private lemma binary_prod_components_at
    {p₁ p₂ : ℝ}
    (hp1₀ : 0 ≤ p₁) (hp1₁ : p₁ ≤ 1)
    (hp2₀ : 0 ≤ p₂) (hp2₁ : p₂ ≤ 1)
    (ij : Fin 2 × Fin 2) :
    ((binary p₁ hp1₀ hp1₁).prod (binary p₂ hp2₀ hp2₁)).prob
        (finProdFinEquiv ij)
      = (binary p₁ hp1₀ hp1₁).prob ij.1 *
        (binary p₂ hp2₀ hp2₁).prob ij.2 := by
  -- `Simplex.prod` is defined via `(finProdFinEquiv.symm k).1` and
  -- `(finProdFinEquiv.symm k).2`. The lemma is essentially `rfl` after
  -- collapsing `finProdFinEquiv.symm (finProdFinEquiv ij) = ij`.
  show (binary p₁ hp1₀ hp1₁).prob (finProdFinEquiv.symm (finProdFinEquiv ij)).1 *
       (binary p₂ hp2₀ hp2₁).prob (finProdFinEquiv.symm (finProdFinEquiv ij)).2
       = _
  simp

/-! ### Csiszár step 1 — gauge-corrected Pexider equation

If `f` generates a continuous f-divergence additive on product
distributions, then there exists a real constant `c` such that the
*gauge-corrected* Pexider equation

  `f(uv) - u·f(v) - v·f(u) + uv·f(1) = -c·(u-1)·(v-1)`

holds for every `u, v > 0`. The constant `c` is the coefficient
of the `(t-1)` gauge in the eventual decomposition
`f(t) = c₁·t·log t + c·(t-1)`.

**Why the gauge term.** The f-divergence is invariant under
`f → f + c·(t-1)` (the linear summand vanishes on normalised
distributions: `∑ q·c·(p/q - 1) = c·(∑p - ∑q) = c·0 = 0`). This is
the standard f-divergence gauge ambiguity (Wikipedia "f-divergence":
*"D_{f(t)} = D_{f(t)+c·(t-1)}; this freedom is not only convenient,
but actually necessary"*). Therefore any functional equation derived
purely from additivity inherits this gauge freedom and **cannot**
have the pure Pexider form `... = 0`; the residual `-c(u-1)(v-1)`
captures the gauge.

**Verification.** For `f(t) = c₁·t·log t + c·(t-1)`, direct algebra:
`f(uv) - u·f(v) - v·f(u) + uv·f(1)`
`= c₁·uv·log(uv) + c·(uv-1) - u·[c₁·v·log v + c·(v-1)]`
`  - v·[c₁·u·log u + c·(u-1)] + uv·0`
`= c·(uv-1) - c·u·(v-1) - c·v·(u-1)`
`= c·(uv - 1 - uv + u - uv + v) = c·(u+v-uv-1) = -c·(u-1)·(v-1)`. ✓

Proof outline (Csiszár 1967, §3, step 1):

1. Specialise `_h_additive` to
   `(binary p₁) ⊗ (binary p₂)  ‖  (binary q₁) ⊗ (binary q₂)`
   with `q₁, q₂ ∈ (0, 1)` strictly.
2. Use `binary_prod_components_at` to expand the four product
   components; this yields a 4-term equation on `f`.
3. Specialising to `p₁ = q₁` collapses to `f(1) = 0`
   (already a hypothesis of the calling `csiszar_generator_decomposition`).
4. Free variation in `(p₁, q₁, p₂, q₂)` then isolates the
   gauge-corrected Pexider equation, where the gauge `c` is
   determined as a specific functional of `f` evaluated on a
   calibration pair.

Effort: ~1 session (the four-term expansion + algebra). -/

/-- **Helper.** Closed form for `fDivergence f` on binary distributions with
    strictly positive denominator. -/
private lemma fDivergence_binary
    (f : ℝ → ℝ) (p q : ℝ)
    (hp0 : 0 ≤ p) (hp1 : p ≤ 1) (hq0 : 0 ≤ q) (hq1 : q ≤ 1)
    (_hq_pos : 0 < q) (_hq_lt : q < 1) :
    fDivergence f (binary p hp0 hp1) (binary q hq0 hq1)
      = q * f (p / q) + (1 - q) * f ((1 - p) / (1 - q)) := by
  unfold fDivergence
  show ∑ i : Fin 2, (binary q hq0 hq1).prob i *
        f ((binary p hp0 hp1).prob i / (binary q hq0 hq1).prob i)
      = q * f (p / q) + (1 - q) * f ((1 - p) / (1 - q))
  rw [Fin.sum_univ_two]
  show (if (0 : Fin 2) = 0 then q else 1 - q) *
        f ((if (0 : Fin 2) = 0 then p else 1 - p)
          / (if (0 : Fin 2) = 0 then q else 1 - q))
      + (if (1 : Fin 2) = 0 then q else 1 - q) *
        f ((if (1 : Fin 2) = 0 then p else 1 - p)
          / (if (1 : Fin 2) = 0 then q else 1 - q))
      = q * f (p / q) + (1 - q) * f ((1 - p) / (1 - q))
  simp

/-- **Helper.** Closed form for `fDivergence f` on the product of two binary
    distributions with strictly positive denominators. The four product terms
    correspond to the four entries of `(binary p₁ ⊗ binary p₂)` against
    `(binary q₁ ⊗ binary q₂)`. -/
private lemma fDivergence_binary_prod
    (f : ℝ → ℝ) (p₁ q₁ p₂ q₂ : ℝ)
    (hp1₀ : 0 ≤ p₁) (hp1₁ : p₁ ≤ 1)
    (hq1₀ : 0 ≤ q₁) (hq1₁ : q₁ ≤ 1)
    (hp2₀ : 0 ≤ p₂) (hp2₁ : p₂ ≤ 1)
    (hq2₀ : 0 ≤ q₂) (hq2₁ : q₂ ≤ 1)
    (_hq1_pos : 0 < q₁) (_hq1_lt : q₁ < 1)
    (_hq2_pos : 0 < q₂) (_hq2_lt : q₂ < 1) :
    fDivergence f
        ((binary p₁ hp1₀ hp1₁).prod (binary p₂ hp2₀ hp2₁))
        ((binary q₁ hq1₀ hq1₁).prod (binary q₂ hq2₀ hq2₁))
      = q₁ * q₂ * f ((p₁ * p₂) / (q₁ * q₂))
        + q₁ * (1 - q₂) * f ((p₁ * (1 - p₂)) / (q₁ * (1 - q₂)))
        + (1 - q₁) * q₂ * f (((1 - p₁) * p₂) / ((1 - q₁) * q₂))
        + (1 - q₁) * (1 - q₂) *
            f (((1 - p₁) * (1 - p₂)) / ((1 - q₁) * (1 - q₂))) := by
  unfold fDivergence
  -- Reindex Fin (2*2) → Fin 2 × Fin 2 via finProdFinEquiv.
  show (∑ k : Fin (2 * 2),
          ((binary q₁ hq1₀ hq1₁).prob (finProdFinEquiv.symm k).1 *
           (binary q₂ hq2₀ hq2₁).prob (finProdFinEquiv.symm k).2) *
          f (((binary p₁ hp1₀ hp1₁).prob (finProdFinEquiv.symm k).1 *
              (binary p₂ hp2₀ hp2₁).prob (finProdFinEquiv.symm k).2) /
             ((binary q₁ hq1₀ hq1₁).prob (finProdFinEquiv.symm k).1 *
              (binary q₂ hq2₀ hq2₁).prob (finProdFinEquiv.symm k).2)))
       = _
  rw [show (∑ k : Fin (2 * 2),
        ((binary q₁ hq1₀ hq1₁).prob (finProdFinEquiv.symm k).1 *
         (binary q₂ hq2₀ hq2₁).prob (finProdFinEquiv.symm k).2) *
        f (((binary p₁ hp1₀ hp1₁).prob (finProdFinEquiv.symm k).1 *
            (binary p₂ hp2₀ hp2₁).prob (finProdFinEquiv.symm k).2) /
           ((binary q₁ hq1₀ hq1₁).prob (finProdFinEquiv.symm k).1 *
            (binary q₂ hq2₀ hq2₁).prob (finProdFinEquiv.symm k).2)))
      = ∑ k : Fin (2 * 2), (fun ij : Fin 2 × Fin 2 =>
          ((binary q₁ hq1₀ hq1₁).prob ij.1 *
           (binary q₂ hq2₀ hq2₁).prob ij.2) *
          f (((binary p₁ hp1₀ hp1₁).prob ij.1 *
              (binary p₂ hp2₀ hp2₁).prob ij.2) /
             ((binary q₁ hq1₀ hq1₁).prob ij.1 *
              (binary q₂ hq2₀ hq2₁).prob ij.2)))
            (finProdFinEquiv.symm k) from rfl]
  rw [Equiv.sum_comp finProdFinEquiv.symm
      (fun ij : Fin 2 × Fin 2 =>
        ((binary q₁ hq1₀ hq1₁).prob ij.1 *
         (binary q₂ hq2₀ hq2₁).prob ij.2) *
        f (((binary p₁ hp1₀ hp1₁).prob ij.1 *
            (binary p₂ hp2₀ hp2₁).prob ij.2) /
           ((binary q₁ hq1₀ hq1₁).prob ij.1 *
            (binary q₂ hq2₀ hq2₁).prob ij.2)))]
  rw [Fintype.sum_prod_type]
  -- Expand the four (Fin 2 × Fin 2) cases explicitly.
  rw [Fin.sum_univ_two, Fin.sum_univ_two, Fin.sum_univ_two]
  -- After unfolding `binary` (`prob 0 = p`, `prob 1 = 1-p`).
  show q₁ * q₂ * f (p₁ * p₂ / (q₁ * q₂))
       + q₁ * (1 - q₂) * f (p₁ * (1 - p₂) / (q₁ * (1 - q₂)))
       + ((1 - q₁) * q₂ * f ((1 - p₁) * p₂ / ((1 - q₁) * q₂))
       + (1 - q₁) * (1 - q₂) * f ((1 - p₁) * (1 - p₂) / ((1 - q₁) * (1 - q₂))))
     = _
  ring

/-- **Helper.** The strict-positivity hypotheses on `binary q _ _` reduce to
    `0 < q` and `q < 1`. -/
private lemma binary_pos_iff {q : ℝ} (hq0 : 0 ≤ q) (hq1 : q ≤ 1)
    (hq_pos : 0 < q) (hq_lt : q < 1) :
    ∀ i, 0 < (binary q hq0 hq1).prob i := by
  intro i
  show 0 < (if i = 0 then q else 1 - q)
  by_cases h : i = 0
  · simp [h, hq_pos]
  · simp [h]; linarith

/-- **Csiszár step 1 — fundamental 4-term equation from binary additivity.**

    Specialising `IsDivergenceAdditive` to the product of two binary
    distributions yields, for every `(p_i, q_i) ∈ [0,1] × (0,1)` with
    `i = 1, 2`:

    `q₁q₂ f(p₁p₂/(q₁q₂)) + q₁(1-q₂) f(p₁(1-p₂)/(q₁(1-q₂)))`
    `+ (1-q₁)q₂ f((1-p₁)p₂/((1-q₁)q₂)) + (1-q₁)(1-q₂) f((1-p₁)(1-p₂)/((1-q₁)(1-q₂)))`
    `= q₁ f(p₁/q₁) + (1-q₁) f((1-p₁)/(1-q₁)) + q₂ f(p₂/q₂) + (1-q₂) f((1-p₂)/(1-q₂))`. -/
private lemma binary_additivity_4_term
    (f : ℝ → ℝ)
    (h_additive : IsDivergenceAdditive
        (fun _ P Q => fDivergence f P Q))
    (p₁ q₁ p₂ q₂ : ℝ)
    (hp1₀ : 0 ≤ p₁) (hp1₁ : p₁ ≤ 1)
    (hq1₀ : 0 ≤ q₁) (hq1₁ : q₁ ≤ 1)
    (hp2₀ : 0 ≤ p₂) (hp2₁ : p₂ ≤ 1)
    (hq2₀ : 0 ≤ q₂) (hq2₁ : q₂ ≤ 1)
    (hq1_pos : 0 < q₁) (hq1_lt : q₁ < 1)
    (hq2_pos : 0 < q₂) (hq2_lt : q₂ < 1) :
    q₁ * q₂ * f ((p₁ * p₂) / (q₁ * q₂))
      + q₁ * (1 - q₂) * f ((p₁ * (1 - p₂)) / (q₁ * (1 - q₂)))
      + (1 - q₁) * q₂ * f (((1 - p₁) * p₂) / ((1 - q₁) * q₂))
      + (1 - q₁) * (1 - q₂) *
          f (((1 - p₁) * (1 - p₂)) / ((1 - q₁) * (1 - q₂)))
      = q₁ * f (p₁ / q₁) + (1 - q₁) * f ((1 - p₁) / (1 - q₁))
        + (q₂ * f (p₂ / q₂) + (1 - q₂) * f ((1 - p₂) / (1 - q₂))) := by
  -- Apply additivity at the appropriate binary distributions.
  have hQ₁ := binary_pos_iff hq1₀ hq1₁ hq1_pos hq1_lt
  have hQ₂ := binary_pos_iff hq2₀ hq2₁ hq2_pos hq2_lt
  have hadd := h_additive
      (binary p₁ hp1₀ hp1₁) (binary q₁ hq1₀ hq1₁)
      (binary p₂ hp2₀ hp2₁) (binary q₂ hq2₀ hq2₁) hQ₁ hQ₂
  -- Beta-reduce the lambda wrapper.
  simp only at hadd
  -- Rewrite LHS via `fDivergence_binary_prod`.
  rw [fDivergence_binary_prod f p₁ q₁ p₂ q₂ hp1₀ hp1₁ hq1₀ hq1₁ hp2₀ hp2₁ hq2₀ hq2₁
        hq1_pos hq1_lt hq2_pos hq2_lt] at hadd
  -- Rewrite RHS via `fDivergence_binary` (twice).
  rw [fDivergence_binary f p₁ q₁ hp1₀ hp1₁ hq1₀ hq1₁ hq1_pos hq1_lt,
      fDivergence_binary f p₂ q₂ hp2₀ hp2₁ hq2₀ hq2₁ hq2_pos hq2_lt] at hadd
  exact hadd

/-- **Csiszár step 1 — `f(1) = 0` is automatic.**

    The additivity hypothesis alone forces `f(1) = 0` (no normalisation
    hypothesis required). Specialise the 4-term equation at
    `p_i = q_i = 1/2`: each ratio becomes `1`, so the LHS collapses to
    `f(1)·(q₁q₂ + q₁(1-q₂) + (1-q₁)q₂ + (1-q₁)(1-q₂)) = f(1)` while the
    RHS becomes `2 f(1)`, forcing `f(1) = 0`. -/
private lemma f_one_zero_from_additivity
    (f : ℝ → ℝ)
    (h_additive : IsDivergenceAdditive
        (fun _ P Q => fDivergence f P Q)) :
    f 1 = 0 := by
  -- Use the 4-term equation at p₁ = q₁ = 1/2, p₂ = q₂ = 1/2.
  have h := binary_additivity_4_term f h_additive
      (1/2 : ℝ) (1/2 : ℝ) (1/2 : ℝ) (1/2 : ℝ)
      (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (by norm_num) (by norm_num) (by norm_num) (by norm_num)
  -- Every ratio is 1, so each term is `(coeff) * f 1`. After collecting,
  -- LHS = f 1 and RHS = 2 f 1, hence f 1 = 0.
  have hkey : f 1 = 2 * f 1 := by
    have h1 : (1/2 : ℝ) * (1/2 : ℝ) / ((1/2 : ℝ) * (1/2 : ℝ)) = 1 := by norm_num
    have h2 : (1/2 : ℝ) * (1 - 1/2) / ((1/2 : ℝ) * (1 - 1/2)) = 1 := by norm_num
    have h3 : (1 - (1/2 : ℝ)) * (1/2 : ℝ) / ((1 - 1/2) * (1/2 : ℝ)) = 1 := by norm_num
    have h4 : (1 - (1/2 : ℝ)) * (1 - 1/2) / ((1 - 1/2) * (1 - 1/2)) = 1 := by norm_num
    have h5 : (1/2 : ℝ) / (1/2 : ℝ) = 1 := by norm_num
    have h6 : (1 - (1/2 : ℝ)) / (1 - (1/2 : ℝ)) = 1 := by norm_num
    rw [h1, h2, h3, h4, h5, h6] at h
    linarith
  linarith

/-- **Csiszár step 1 — anchor on `(1, ∞)²`.**

    Specialising the 4-term binary additivity equation at `p₁ = p₂ = 0`
    yields, with `α := 1/(1-q₁)`, `β := 1/(1-q₂)` (both `> 1`):

      `f(αβ) - α·f(β) - β·f(α) = (α-1)(β-1)·f(0)`.

    Equivalently, setting `c := -f(0)`, and using `f(1) = 0`:

      `f(αβ) - α·f(β) - β·f(α) + αβ·f(1) = -c·(α-1)·(β-1)`. -/
private lemma residual_anchor_gt_one
    (f : ℝ → ℝ)
    (h_additive : IsDivergenceAdditive
        (fun _ P Q => fDivergence f P Q))
    (α β : ℝ) (hα : 1 < α) (hβ : 1 < β) :
    f (α * β) - α * f β - β * f α + α * β * f 1 = f 0 * (α - 1) * (β - 1) := by
  -- Define q₁ = 1 - 1/α, q₂ = 1 - 1/β; then α = 1/(1-q₁), β = 1/(1-q₂).
  have hα_pos : 0 < α := lt_trans one_pos hα
  have hβ_pos : 0 < β := lt_trans one_pos hβ
  have hα_ne : α ≠ 0 := hα_pos.ne'
  have hβ_ne : β ≠ 0 := hβ_pos.ne'
  set q₁ : ℝ := 1 - 1/α with hq₁_def
  set q₂ : ℝ := 1 - 1/β with hq₂_def
  have hq₁_pos : 0 < q₁ := by
    show 0 < 1 - 1/α
    rw [sub_pos, div_lt_one hα_pos]; exact hα
  have hq₁_lt : q₁ < 1 := by
    show 1 - 1/α < 1
    have : 0 < 1/α := by positivity
    linarith
  have hq₂_pos : 0 < q₂ := by
    show 0 < 1 - 1/β
    rw [sub_pos, div_lt_one hβ_pos]; exact hβ
  have hq₂_lt : q₂ < 1 := by
    show 1 - 1/β < 1
    have : 0 < 1/β := by positivity
    linarith
  have h_1mq₁_α : (1 - q₁) * α = 1 := by
    show (1 - (1 - 1/α)) * α = 1
    have : (1 - (1 - 1/α)) = 1/α := by ring
    rw [this]; field_simp
  have h_1mq₂_β : (1 - q₂) * β = 1 := by
    show (1 - (1 - 1/β)) * β = 1
    have : (1 - (1 - 1/β)) = 1/β := by ring
    rw [this]; field_simp
  -- Apply the 4-term equation at p₁ = p₂ = 0.
  have h4 := binary_additivity_4_term f h_additive
      (0 : ℝ) q₁ (0 : ℝ) q₂
      (le_refl 0) (by linarith) hq₁_pos.le hq₁_lt.le
      (le_refl 0) (by linarith) hq₂_pos.le hq₂_lt.le
      hq₁_pos hq₁_lt hq₂_pos hq₂_lt
  -- Simplify ratios.
  have h00 : (0 : ℝ) * 0 / (q₁ * q₂) = 0 := by norm_num
  have h01 : (0 : ℝ) * (1 - 0) / (q₁ * (1 - q₂)) = 0 := by norm_num
  have h10 : (1 - 0) * (0 : ℝ) / ((1 - q₁) * q₂) = 0 := by norm_num
  have h_1mq₁_pos : 0 < 1 - q₁ := by linarith
  have h_1mq₂_pos : 0 < 1 - q₂ := by linarith
  have h_1mq₁_ne : 1 - q₁ ≠ 0 := h_1mq₁_pos.ne'
  have h_1mq₂_ne : 1 - q₂ ≠ 0 := h_1mq₂_pos.ne'
  have h11 : (1 - 0) * (1 - 0) / ((1 - q₁) * (1 - q₂)) = α * β := by
    have h1mq₁ : 1 - q₁ = 1/α := by show 1 - (1 - 1/α) = 1/α; ring
    have h1mq₂ : 1 - q₂ = 1/β := by show 1 - (1 - 1/β) = 1/β; ring
    rw [h1mq₁, h1mq₂]
    field_simp
    ring
  have h_p1_div : (0 : ℝ) / q₁ = 0 := zero_div _
  have h_p2_div : (0 : ℝ) / q₂ = 0 := zero_div _
  have h_1mp1 : (1 - 0) / (1 - q₁) = α := by
    have h1mq₁ : 1 - q₁ = 1/α := by show 1 - (1 - 1/α) = 1/α; ring
    rw [h1mq₁]; rw [sub_zero]; field_simp
  have h_1mp2 : (1 - 0) / (1 - q₂) = β := by
    have h1mq₂ : 1 - q₂ = 1/β := by show 1 - (1 - 1/β) = 1/β; ring
    rw [h1mq₂]; rw [sub_zero]; field_simp

  rw [h00, h01, h10, h11, h_p1_div, h_p2_div, h_1mp1, h_1mp2] at h4
  -- h4 now reads:
  -- q₁ q₂ f 0 + q₁ (1 - q₂) f 0 + (1 - q₁) q₂ f 0 + (1 - q₁)(1 - q₂) f(α β)
  --   = q₁ f 0 + (1 - q₁) f α + (q₂ f 0 + (1 - q₂) f β)
  -- Goal: f(αβ) - α f β - β f α + αβ f 1 = f 0 (α - 1)(β - 1).
  have hf1 : f 1 = 0 := f_one_zero_from_additivity f h_additive
  -- Express q₁ via α and q₂ via β: q₁ = 1 - 1/α and q₂ = 1 - 1/β.
  -- Substitute these into h4 to obtain a polynomial identity in α, β.
  have hq1_alt : q₁ = 1 - 1/α := rfl
  have hq2_alt : q₂ = 1 - 1/β := rfl
  -- Multiply both sides of h4 by α*β to clear denominators.
  have h4_scaled : (q₁ * q₂ * f 0 + q₁ * (1 - q₂) * f 0
                  + (1 - q₁) * q₂ * f 0 + (1 - q₁) * (1 - q₂) * f (α * β))
                * (α * β)
              = (q₁ * f 0 + (1 - q₁) * f α + (q₂ * f 0 + (1 - q₂) * f β)) * (α * β) := by
    rw [h4]
  -- Substitute q₁ = (α-1)/α, q₂ = (β-1)/β so that (1-q₁) = 1/α etc.
  rw [hq1_alt, hq2_alt] at h4_scaled
  have h4_poly : (α - 1) * (β - 1) * f 0 + (α - 1) * f 0 + (β - 1) * f 0 + f (α * β)
              = (α - 1) * β * f 0 + β * f α + α * (β - 1) * f 0 + α * f β := by
    have := h4_scaled
    have hαβ_pos : (0 : ℝ) < α * β := mul_pos hα_pos hβ_pos
    have hαβ_ne : α * β ≠ 0 := hαβ_pos.ne'
    field_simp at this
    linarith
  linear_combination h4_poly + α * β * hf1

/-- **Csiszár step 2 — extension to `[0, 1) × (1, ∞)`.**

    For `u ∈ [0, 1)` and `β > 1`, the residual identity
    `f(uβ) - u·f(β) - β·f(u) + uβ·f(1) = f(0)·(u-1)·(β-1)` holds.

    Proof: apply the 4-term equation at `p₁ = u·q₁`, `p₂ = 0`,
    with `q₁ := 1/2` (any value in `(0, 1)` works), `q₂ := 1 - 1/β`.
    Then `a_0 = u`, `a_1 = (1 - u q₁)/(1 - q₁) =: \bar u > 1`,
    `b_0 = 0`, `b_1 = β`. Substituting the anchor identity at
    `(\bar u, β)` and simplifying (using the identity
    `β(1-q₂) = 1`, which makes the `f(\bar u)` terms cancel)
    yields the residual at `(u, β)`. -/
private lemma residual_left_lt_one
    (f : ℝ → ℝ)
    (h_additive : IsDivergenceAdditive
        (fun _ P Q => fDivergence f P Q))
    (u β : ℝ) (hu_pos : 0 < u) (hu_lt_one : u < 1) (hβ : 1 < β) :
    f (u * β) - u * f β - β * f u + u * β * f 1 = f 0 * (u - 1) * (β - 1) := by
  -- Set up: q₁ := 1/2, q₂ := 1 - 1/β, p₁ := u·q₁ = u/2, p₂ := 0.
  set q₁ : ℝ := 1/2 with hq₁_def
  have hq₁_pos : 0 < q₁ := by show (0 : ℝ) < 1/2; norm_num
  have hq₁_lt : q₁ < 1 := by show (1 : ℝ)/2 < 1; norm_num
  have h_1mq₁_pos : 0 < 1 - q₁ := by show 0 < 1 - (1:ℝ)/2; norm_num
  have h_1mq₁_ne : (1 - q₁) ≠ 0 := h_1mq₁_pos.ne'
  have hq₁_ne : q₁ ≠ 0 := hq₁_pos.ne'
  have hβ_pos : 0 < β := lt_trans one_pos hβ
  have hβ_ne : β ≠ 0 := hβ_pos.ne'
  set q₂ : ℝ := 1 - 1/β with hq₂_def
  have hq₂_pos : 0 < q₂ := by
    show 0 < 1 - 1/β
    rw [sub_pos, div_lt_one hβ_pos]; exact hβ
  have hq₂_lt : q₂ < 1 := by
    show 1 - 1/β < 1
    have : 0 < 1/β := by positivity
    linarith
  have h_1mq₂_eq : 1 - q₂ = 1/β := by show 1 - (1 - 1/β) = 1/β; ring
  have h_1mq₂_β : (1 - q₂) * β = 1 := by rw [h_1mq₂_eq]; field_simp
  -- p₁ := u·q₁.
  set p₁ : ℝ := u * q₁ with hp₁_def
  have hp₁_nonneg : 0 ≤ p₁ := mul_nonneg hu_pos.le hq₁_pos.le
  have hp₁_le_one : p₁ ≤ 1 := by
    show u * q₁ ≤ 1
    have : u * q₁ ≤ 1 * q₁ := by
      exact mul_le_mul_of_nonneg_right hu_lt_one.le hq₁_pos.le
    linarith [this]
  -- Apply the 4-term equation.
  have h4 := binary_additivity_4_term f h_additive
      p₁ q₁ (0 : ℝ) q₂
      hp₁_nonneg hp₁_le_one hq₁_pos.le hq₁_lt.le
      (le_refl 0) (by linarith) hq₂_pos.le hq₂_lt.le
      hq₁_pos hq₁_lt hq₂_pos hq₂_lt
  -- Simplify ratios.
  -- p₁ * 0 / (q₁ * q₂) = 0
  have h_pp_div : p₁ * 0 / (q₁ * q₂) = 0 := by norm_num
  -- p₁ * (1 - 0) / (q₁ * (1 - q₂)) = u * β
  have h_pq_div : p₁ * (1 - 0) / (q₁ * (1 - q₂)) = u * β := by
    rw [sub_zero, mul_one, h_1mq₂_eq, hp₁_def]
    -- u * q₁ / (q₁ * (1/β)) = u * β
    rw [div_eq_iff (by positivity : q₁ * (1/β) ≠ 0)]
    field_simp
  -- (1 - p₁) * 0 / ((1 - q₁) * q₂) = 0
  have h_1pp_div : (1 - p₁) * 0 / ((1 - q₁) * q₂) = 0 := by norm_num
  -- ū := (1 - p₁) / (1 - q₁)
  set ū : ℝ := (1 - p₁) / (1 - q₁) with hū_def
  have hū_gt_one : 1 < ū := by
    rw [hū_def, hp₁_def, hq₁_def]
    rw [lt_div_iff₀ (by norm_num : (0 : ℝ) < 1 - 1/2)]
    have : u * (1/2 : ℝ) < 1 * (1/2 : ℝ) := by
      exact mul_lt_mul_of_pos_right hu_lt_one (by norm_num)
    linarith [this]
  have hū_pos : 0 < ū := lt_trans one_pos hū_gt_one
  -- (1 - p₁) * (1 - 0) / ((1 - q₁) * (1 - q₂)) = ū * β
  have h_1p1q_div : (1 - p₁) * (1 - 0) / ((1 - q₁) * (1 - q₂)) = ū * β := by
    rw [sub_zero, mul_one, h_1mq₂_eq]
    rw [div_eq_iff (by positivity : (1 - q₁) * (1/β) ≠ 0)]
    rw [hū_def]
    field_simp
  -- p₁/q₁ = u
  have h_p1q1_div : p₁ / q₁ = u := by
    rw [hp₁_def]; field_simp
  -- (1 - p₁)/(1 - q₁) = ū (by definition)
  have h_1p1q1_div : (1 - p₁) / (1 - q₁) = ū := rfl
  -- 0/q₂ = 0
  have h_0q2 : (0 : ℝ) / q₂ = 0 := zero_div _
  -- (1-0)/(1-q₂) = β
  have h_1mz_1mq₂ : (1 - 0) / (1 - q₂) = β := by
    rw [h_1mq₂_eq, sub_zero]; field_simp
  rw [h_pp_div, h_pq_div, h_1pp_div, h_1p1q_div, h_p1q1_div,
      h_0q2, h_1mz_1mq₂] at h4
  -- h4 now reads:
  -- q₁ q₂ f 0 + q₁ (1 - q₂) f(uβ) + (1 - q₁) q₂ f 0 + (1 - q₁)(1 - q₂) f(ūβ)
  --  = q₁ f u + (1 - q₁) f ū + (q₂ f 0 + (1 - q₂) f β)
  -- Anchor: f(ūβ) = ū f β + β f ū - f(1)·ū·β + f 0 (ū-1)(β-1)
  --  i.e. f(ūβ) - ū f β - β f ū + ū β f 1 = f 0 (ū-1)(β-1).
  have hanchor := residual_anchor_gt_one f h_additive ū β hū_gt_one hβ
  have hf1 : f 1 = 0 := f_one_zero_from_additivity f h_additive
  -- Use hanchor to substitute f(ūβ); then derive the target from h4.
  -- Multiply h4 by β = 1/(1-q₂); we'll use h_1mq₂_β and the constraint
  -- (1-q₁)ū = 1 - u q₁ (which holds since q₁ u + (1-q₁) ū = 1, equivalent
  -- to (1-q₁) ū = 1 - q₁ p₁/q₁... wait p₁ = u q₁, so 1 - q₁ u, so
  -- ū = (1 - u q₁)/(1 - q₁) and (1-q₁) ū = 1 - u q₁).
  have h_constraint : q₁ * u + (1 - q₁) * ū = 1 := by
    rw [hū_def, hp₁_def]
    field_simp
    ring
  -- Strategy: scale h4 by β, substitute hanchor (multiplied by appropriate
  -- coefficient), and use linear_combination.
  -- After scaling h4 by β and using h_1mq₂_β = 1 (i.e., (1-q₂)β = 1):
  -- β q₁ q₂ f 0 + q₁ f(uβ) + β(1-q₁) q₂ f 0 + (1-q₁) f(ūβ)
  --   = β q₁ f u + β(1-q₁) f ū + (β q₂ f 0 + f β)
  -- Substitute f(ūβ) = ū f β + β f ū + f 0 (ū-1)(β-1) (from hanchor, hf1):
  -- β q₁ q₂ f 0 + q₁ f(uβ) + β(1-q₁) q₂ f 0 + (1-q₁)·[ū f β + β f ū + f 0 (ū-1)(β-1)]
  --   = β q₁ f u + β(1-q₁) f ū + β q₂ f 0 + f β
  -- (1-q₁) β f ū on LHS cancels with β(1-q₁) f ū on RHS.
  -- (1-q₁) ū f β = (1 - q₁ u) f β (using h_constraint).
  -- So: β q₁ q₂ f 0 + q₁ f(uβ) + β(1-q₁) q₂ f 0 + (1 - q₁ u) f β + (1-q₁)(ū-1)(β-1) f 0
  --   = β q₁ f u + β q₂ f 0 + f β
  -- => q₁ f(uβ) - q₁ u f β - q₁ β f u
  --    = β q₂ f 0 - β q₁ q₂ f 0 - β(1-q₁) q₂ f 0 - (1-q₁)(ū-1)(β-1) f 0
  --    = β q₂ f 0 (1 - q₁ - (1-q₁)) f 0 - (1-q₁)(ū-1)(β-1) f 0
  --    = 0 f 0 - (1-q₁)(ū-1)(β-1) f 0
  -- Wait β q₂ f 0 (1 - q₁ - (1-q₁)) = 0. So
  -- q₁ [f(uβ) - u f β - β f u] = -(1-q₁)(ū-1)(β-1) f 0
  -- Now (1-q₁)(ū-1) = (1-q₁)ū - (1-q₁) = (1 - q₁ u) - (1 - q₁) = q₁ - q₁ u = q₁(1-u)
  -- So q₁ [f(uβ) - u f β - β f u] = -q₁(1-u)(β-1) f 0
  -- => f(uβ) - u f β - β f u = -(1-u)(β-1) f 0 = (u-1)(β-1) f 0. ✓
  linear_combination
    (2 * β) * h4 + (-(1 : ℝ)) * hanchor + (2 * β) * hf1
    + (-f (u * β) - f (ū * β) + 2 * f β) * h_1mq₂_β
    + (-2 * f β + 2 * β * f 1 - 2 * (β - 1) * f 0) * h_constraint

/-- **Csiszár step 2′ — extension to `(1, ∞) × [0, 1)`.**

    Symmetric counterpart of `residual_left_lt_one`. Proof: same strategy
    with the roles of `(p₁, q₁)` and `(p₂, q₂)` swapped. -/
private lemma residual_right_lt_one
    (f : ℝ → ℝ)
    (h_additive : IsDivergenceAdditive
        (fun _ P Q => fDivergence f P Q))
    (α v : ℝ) (hα : 1 < α) (hv_pos : 0 < v) (hv_lt_one : v < 1) :
    f (α * v) - α * f v - v * f α + α * v * f 1 = f 0 * (α - 1) * (v - 1) := by
  -- Use commutativity: swap the roles to reduce to `residual_left_lt_one`.
  have h := residual_left_lt_one f h_additive v α hv_pos hv_lt_one hα
  rw [show α * v = v * α from mul_comm α v]
  linear_combination h

/-- **Csiszár step 3 — extension to `(0, 1) × (0, 1)`.**

    Apply the 4-term equation at `p₁ = u/2`, `p₂ = v/2` (i.e.,
    `q₁ = q₂ = 1/2`). Then `ū = 2 - u, vb = 2 - v` are both `> 1`.
    Substituting the three already-proven residuals (anchor at
    `(ū, vb)`, left at `(u, vb)`, right at `(ū, v)`) and using
    `f(1) = 0`, the target residual at `(u, v)` follows by direct
    algebra (the cross terms collapse via `u + ū = 2`, `v + vb = 2`). -/
private lemma residual_both_lt_one
    (f : ℝ → ℝ)
    (h_additive : IsDivergenceAdditive
        (fun _ P Q => fDivergence f P Q))
    (u v : ℝ) (hu_pos : 0 < u) (hu_lt_one : u < 1)
              (hv_pos : 0 < v) (hv_lt_one : v < 1) :
    f (u * v) - u * f v - v * f u + u * v * f 1 = f 0 * (u - 1) * (v - 1) := by
  -- Set q = 1/2.
  set q : ℝ := 1/2 with hq_def
  have hq_pos : 0 < q := by show (0 : ℝ) < 1/2; norm_num
  have hq_lt : q < 1 := by show (1 : ℝ)/2 < 1; norm_num
  have h_1mq_pos : 0 < 1 - q := by show 0 < 1 - (1:ℝ)/2; norm_num
  -- p₁ := u/2, p₂ := v/2.
  set p₁ : ℝ := u * q with hp₁_def
  have hp₁_nonneg : 0 ≤ p₁ := mul_nonneg hu_pos.le hq_pos.le
  have hp₁_le_one : p₁ ≤ 1 := by
    show u * q ≤ 1
    have : u * q ≤ 1 * q := mul_le_mul_of_nonneg_right hu_lt_one.le hq_pos.le
    linarith
  set p₂ : ℝ := v * q with hp₂_def
  have hp₂_nonneg : 0 ≤ p₂ := mul_nonneg hv_pos.le hq_pos.le
  have hp₂_le_one : p₂ ≤ 1 := by
    show v * q ≤ 1
    have : v * q ≤ 1 * q := mul_le_mul_of_nonneg_right hv_lt_one.le hq_pos.le
    linarith
  have h4 := binary_additivity_4_term f h_additive
      p₁ q p₂ q
      hp₁_nonneg hp₁_le_one hq_pos.le hq_lt.le
      hp₂_nonneg hp₂_le_one hq_pos.le hq_lt.le
      hq_pos hq_lt hq_pos hq_lt
  -- ū := 2 - u, vb := 2 - v (since q = 1/2).
  set ū : ℝ := (1 - p₁) / (1 - q) with hū_def
  set vb : ℝ := (1 - p₂) / (1 - q) with hvb_def
  have h_ū_eq : ū = 2 - u := by
    show (1 - u * q) / (1 - q) = 2 - u
    rw [hq_def]
    field_simp; ring
  have h_vb_eq : vb = 2 - v := by
    show (1 - v * q) / (1 - q) = 2 - v
    rw [hq_def]
    field_simp; ring
  have hū_gt_one : 1 < ū := by rw [h_ū_eq]; linarith
  have hvb_gt_one : 1 < vb := by rw [h_vb_eq]; linarith
  -- Simplify the four ratios.
  have h_a0b0 : p₁ * p₂ / (q * q) = u * v := by
    rw [hp₁_def, hp₂_def]; field_simp
  have h_a0b1 : p₁ * (1 - p₂) / (q * (1 - q)) = u * vb := by
    show p₁ * (1 - p₂) / (q * (1 - q)) = u * ((1 - p₂) / (1 - q))
    rw [hp₁_def]
    field_simp
  have h_a1b0 : (1 - p₁) * p₂ / ((1 - q) * q) = ū * v := by
    show (1 - p₁) * p₂ / ((1 - q) * q) = ((1 - p₁) / (1 - q)) * v
    rw [hp₂_def]
    field_simp
  have h_a1b1 : (1 - p₁) * (1 - p₂) / ((1 - q) * (1 - q)) = ū * vb := by
    show (1 - p₁) * (1 - p₂) / ((1 - q) * (1 - q))
         = ((1 - p₁) / (1 - q)) * ((1 - p₂) / (1 - q))
    field_simp
  have h_p1q : p₁ / q = u := by rw [hp₁_def]; field_simp
  have h_p2q : p₂ / q = v := by rw [hp₂_def]; field_simp
  -- (1 - p₁)/(1 - q) = ū and (1 - p₂)/(1 - q) = vb by definition.
  rw [h_a0b0, h_a0b1, h_a1b0, h_a1b1, h_p1q, h_p2q] at h4
  -- h4 now reads:
  -- q² f(uv) + q(1-q) f(uvb) + (1-q) q f(ūv) + (1-q)² f(ūvb)
  --   = q f u + (1-q) f ū + q f v + (1-q) f vb
  -- The three residual substitutions:
  have h_anchor := residual_anchor_gt_one f h_additive ū vb hū_gt_one hvb_gt_one
  have h_left := residual_left_lt_one f h_additive u vb hu_pos hu_lt_one hvb_gt_one
  have h_right := residual_right_lt_one f h_additive ū v hū_gt_one hv_pos hv_lt_one
  have hf1 : f 1 = 0 := f_one_zero_from_additivity f h_additive
  -- The q = 1/2 substitutions give u + ū = 2 and v + vb = 2.
  have h_uū : u + ū = 2 := by rw [h_ū_eq]; ring
  have h_vvb : v + vb = 2 := by rw [h_vb_eq]; ring
  -- Scale h4 by 4 (= 1/q²) and substitute the three residuals.
  -- After algebra (see commentary above): goal follows from
  -- 4 h4 + (-2) h_left + (-2) h_right + (-1) h_anchor
  -- plus appropriate multiples of h_uū and h_vvb and hf1.
  -- Let G = goal. We computed:
  -- 4 h4 ≡ (4·q²) f(uv) + (4q(1-q)) f(uvb) + (4(1-q)q) f(ūv) + (4(1-q)²) f(ūvb) - ...
  --       = f(uv) + f(uvb) + f(ūv) + f(ūvb) - 2 f u - 2 f ū - 2 f v - 2 f vb
  -- substitute the three residuals on f(uvb), f(ūv), f(ūvb) to eliminate
  -- those f-terms.
  linear_combination
    4 * h4 - h_left - h_right - h_anchor
    - ((f v + f vb) + (v + vb - 2) * f 0) * h_uū
    - (f u + f ū) * h_vvb
    + (u + ū) * (v + vb) * hf1

/-- **Csiszár step 1 (residual extraction).**

    The gauge-corrected Pexider equation follows from the binary
    4-term equation. With `c := -f(0)`, the identity

      `f(uv) - u·f(v) - v·f(u) + uv·f(1) = -c·(u-1)·(v-1)`

    holds for every `u, v > 0`. This is established by:

    1. **Anchor** (`residual_anchor_gt_one`): the identity holds for
       `(u, v) ∈ (1, ∞)²`, by specialising the 4-term equation at
       `p₁ = p₂ = 0` and substituting `α = 1/(1-q₁), β = 1/(1-q₂)`.

    2. **Extension to `[0, 1) × (1, ∞)`**: specialise the 4-term equation
       at `p₁ = u·q₁` (free, `u < 1` so that `\bar u > 1`) and `p₂ = 0`;
       substitute the anchor identity at `(\bar u, β)` and simplify. The
       `f(\bar u)` terms cancel exactly via the identity
       `β·(1-q₂) - 1 = 0`, leaving the residual at `(u, β)`.

    3. **Extension to `(1, ∞) × [0, 1)`**: symmetric to step 2.

    4. **Extension to `[0, 1) × [0, 1)`**: chain the previous extensions.

    The implementation currently completes step 1 (anchor) but leaves
    the iterated chain (steps 2–4) as `sorry` for time reasons. The
    `csiszar_generator_decomposition` caller only invokes
    `h_residual u v` for arbitrary `u, v > 0`, so closing the chain
    completes the proof. -/
private lemma csiszar_residual_from_binary_additivity
    (f : ℝ → ℝ)
    (h_additive : IsDivergenceAdditive
        (fun _ P Q => fDivergence f P Q)) :
    ∃ c : ℝ, ∀ u v : ℝ, 0 < u → 0 < v →
      f (u * v) - u * f v - v * f u + u * v * f 1 = -c * (u - 1) * (v - 1) := by
  -- The candidate gauge is c := -f(0). With this choice, the target
  -- becomes f(uv) - u f v - v f u + uv f 1 = f(0)·(u-1)·(v-1), which is
  -- exactly the statement of `residual_anchor_gt_one` on (1, ∞)² and
  -- propagates to all (0, ∞)² via the iterated chain described in the
  -- doc above.
  refine ⟨-f 0, ?_⟩
  intro u v hu hv
  -- Convert -(-f 0) = f 0.
  have hgauge : (-(-f 0)) = f 0 := by ring
  -- Split by quadrants in (0, ∞)² based on positions relative to 1.
  rcases lt_trichotomy u 1 with hu1 | hu1 | hu1
  · -- u < 1
    rcases lt_trichotomy v 1 with hv1 | hv1 | hv1
    · -- u < 1, v < 1: corner case.
      have := residual_both_lt_one f h_additive u v hu hu1 hv hv1
      linear_combination this
    · -- u < 1, v = 1: trivial residual.
      subst hv1
      simp only [mul_one]
      ring
    · -- u < 1, v > 1.
      have := residual_left_lt_one f h_additive u v hu hu1 hv1
      linear_combination this
  · -- u = 1: trivial.
    subst hu1
    simp only [one_mul]
    ring
  · -- u > 1
    rcases lt_trichotomy v 1 with hv1 | hv1 | hv1
    · -- u > 1, v < 1.
      have := residual_right_lt_one f h_additive u v hu1 hv hv1
      linear_combination this
    · -- u > 1, v = 1: trivial.
      subst hv1
      simp only [mul_one]
      ring
    · -- u > 1, v > 1: anchor.
      have := residual_anchor_gt_one f h_additive u v hu1 hv1
      linear_combination this

/-- **Csiszár step 2 — Cauchy-log reduction.**

    From the Pexider equation `f(uv) = u f(v) + v f(u) - uv f(1)` for
    `u, v > 0`, the function

    `H : ℝ → ℝ`,  `H(t) := if 0 < t then f t / t - f 1 else 0`

    is continuous on `(0, ∞)` (because `f` is continuous and `t > 0`)
    and satisfies the multiplicative Cauchy equation
    `H(uv) = H(u) + H(v)` for every `u, v > 0`.

    Indeed, dividing the Pexider equation by `uv`:
    `f(uv)/(uv) = f(v)/v + f(u)/u - f(1)`, so
    `H(uv) = f(uv)/(uv) - f(1) = (f(u)/u - f(1)) + (f(v)/v - f(1)) = H(u) + H(v)`. -/
private lemma cauchy_log_from_pexider
    (f : ℝ → ℝ)
    (h_pexider : ∀ u v : ℝ, 0 < u → 0 < v →
        f (u * v) - u * f v - v * f u + u * v * f 1 = 0) :
    ∀ u v : ℝ, 0 < u → 0 < v →
      (f (u * v) / (u * v) - f 1) = (f u / u - f 1) + (f v / v - f 1) := by
  intro u v hu hv
  have huv_pos : 0 < u * v := mul_pos hu hv
  have huv_ne : u * v ≠ 0 := huv_pos.ne'
  have hu_ne : u ≠ 0 := hu.ne'
  have hv_ne : v ≠ 0 := hv.ne'
  have hpex := h_pexider u v hu hv
  -- From hpex: f(uv) = u f v + v f u - uv f 1.
  have hfuv : f (u * v) = u * f v + v * f u - u * v * f 1 := by linarith
  -- Divide by uv.
  have key : f (u * v) / (u * v) = f v / v + f u / u - f 1 := by
    rw [hfuv]
    field_simp
  linarith

/-- **Csiszár steps 1–3 (functional equation, Cauchy-log, back-substitution).**

    Given a generator `f` continuous on `ℝ` with `f 1 = 0` such that the
    associated f-divergence is additive under product distributions,
    there exist real constants `c₁, c₂` such that

    `f(t) = c₁ · t · log t + c₂ · (t − 1)` for every `t > 0`,
    and `f 0 = -c₂`.

    Under the normalisation `f 1 = 0`, the back-substitution yields
    `c₂ = 0` (so the affine term vanishes) and `f(0) = 0 = -c₂` by
    continuity of `f` at the boundary `t = 0` combined with the
    classical limit `lim_{t→0⁺} t · log t = 0` (`continuous_mul_log`).

    The proof has four movements:

    * **Step 0 (gauge extraction).** `csiszar_residual_from_binary_additivity`
      extracts the gauge constant `c` (the c₂ coefficient) and the
      gauge-corrected Pexider residual
      `f(uv) − u·f(v) − v·f(u) + uv·f(1) = -c·(u-1)·(v-1)`.

    * **Step 1 (gauge subtraction).** Define `g(t) := f(t) − c·(t-1)`.
      Then `g(1) = 0` (from `f(1) = 0`), `g` is continuous, and `g`
      satisfies the *pure* Pexider equation (no gauge term), because
      the residual on `f` and the gauge term on `g(uv) − u·g(v) − v·g(u)`
      cancel exactly (algebraic identity `(uv-1) − u(v-1) − v(u-1) =
      -(u-1)(v-1)`).

    * **Step 2 (Cauchy-log on g).** `cauchy_log_from_pexider g hg_pexider`
      gives the multiplicative Cauchy equation for `H(t) := g(t)/t`.
      Package as `AddMonoidHom ℝ →+ ℝ` via `t ↦ g(exp t)/exp t`, apply
      `map_real_smul` for ℝ-linearity, conclude `g(t) = k · t · log t`.

    * **Step 3 (back-substitution + boundary).** With `c₁ := k` and
      `c₂ := c`, the form `f(t) = c₁·t·log t + c·(t-1)` holds on
      `t > 0`. Continuity of `f` at `0` plus the form's limit
      `lim_{t→0⁺} (k·t·log t + c·(t-1)) = 0 − c = -c` forces
      `f(0) = -c = -c₂`. -/
private lemma csiszar_generator_decomposition
    (f : ℝ → ℝ) (hf_cont : Continuous f)
    (hf_one : f 1 = 0)
    (h_additive : IsDivergenceAdditive
        (fun _ P Q => fDivergence f P Q)) :
    ∃ c₁ c₂ : ℝ,
      (∀ t : ℝ, 0 < t → f t = c₁ * (t * Real.log t) + c₂ * (t - 1)) ∧
      f 0 = -c₂ := by
  -- Step 0: extract the gauge constant `c` from the residual Pexider.
  obtain ⟨c, h_residual⟩ := csiszar_residual_from_binary_additivity f h_additive
  -- Step 1: define `g(t) := f(t) - c·(t-1)`. Then `g(1) = 0` and `g`
  -- satisfies the *pure* Pexider equation, by direct algebra cancelling
  -- the gauge residual.
  set g : ℝ → ℝ := fun t => f t - c * (t - 1) with hg_def
  have hg_cont : Continuous g := by
    show Continuous (fun t => f t - c * (t - 1))
    exact hf_cont.sub (continuous_const.mul (continuous_id.sub continuous_const))
  have hg_one : g 1 = 0 := by
    show f 1 - c * (1 - 1) = 0
    rw [hf_one]; ring
  have hg_pexider : ∀ u v : ℝ, 0 < u → 0 < v →
      g (u * v) - u * g v - v * g u + u * v * g 1 = 0 := by
    intro u v hu hv
    have hres := h_residual u v hu hv
    -- hres: f(uv) - u f v - v f u + uv f 1 = -c (u-1)(v-1)
    -- After expanding g, the goal becomes a polynomial identity in
    -- f(uv), f(v), f(u), f(1), u, v, c that follows from `hres` via
    -- `ring` once we treat f's values as opaque constants.
    show f (u * v) - c * (u * v - 1) - u * (f v - c * (v - 1))
         - v * (f u - c * (u - 1)) + u * v * (f 1 - c * (1 - 1)) = 0
    linear_combination hres
  -- Step 2: apply the existing cauchy-log chain to g (not f).
  let htilde : ℝ → ℝ := fun x => g (Real.exp x) / Real.exp x - g 1
  have htilde_cont : Continuous htilde := by
    show Continuous (fun x => g (Real.exp x) / Real.exp x - g 1)
    refine Continuous.sub ?_ continuous_const
    refine Continuous.div ?_ Real.continuous_exp ?_
    · exact hg_cont.comp Real.continuous_exp
    · intro x; exact (Real.exp_pos x).ne'
  have htilde_add : ∀ s t : ℝ, htilde (s + t) = htilde s + htilde t := by
    intro s t
    show g (Real.exp (s + t)) / Real.exp (s + t) - g 1
        = (g (Real.exp s) / Real.exp s - g 1) + (g (Real.exp t) / Real.exp t - g 1)
    have hes : 0 < Real.exp s := Real.exp_pos s
    have het : 0 < Real.exp t := Real.exp_pos t
    have hcl := cauchy_log_from_pexider g hg_pexider (Real.exp s) (Real.exp t) hes het
    rw [Real.exp_add]
    linarith
  have htilde_zero : htilde 0 = 0 := by
    show g (Real.exp 0) / Real.exp 0 - g 1 = 0
    rw [Real.exp_zero, div_one, hg_one]
    ring
  let hHom : ℝ →+ ℝ :=
    { toFun := htilde, map_zero' := htilde_zero, map_add' := htilde_add }
  have hHom_cont : Continuous hHom := htilde_cont
  have htilde_linear : ∀ c' : ℝ, htilde c' = c' * htilde 1 := by
    intro c'
    have := map_real_smul hHom hHom_cont c' (1 : ℝ)
    simp only [smul_eq_mul, mul_one] at this
    exact this
  set k : ℝ := htilde 1 with hk_def
  have hg_on_pos : ∀ t : ℝ, 0 < t → g t = k * (t * Real.log t) := by
    intro t ht
    have hexp_log : Real.exp (Real.log t) = t := Real.exp_log ht
    have h1 : htilde (Real.log t) = g t / t - g 1 := by
      show g (Real.exp (Real.log t)) / Real.exp (Real.log t) - g 1
          = g t / t - g 1
      rw [hexp_log]
    have h2 : htilde (Real.log t) = Real.log t * k := htilde_linear (Real.log t)
    rw [hg_one, sub_zero] at h1
    have ht_ne : t ≠ 0 := ht.ne'
    have hgt_div : g t / t = k * Real.log t := by
      rw [← h1, h2, mul_comm]
    have hgt_eq : g t = k * Real.log t * t := by
      field_simp at hgt_div
      linarith
    rw [hgt_eq]; ring
  -- Step 3: back-substitute g → f and handle the boundary.
  refine ⟨k, c, ?_, ?_⟩
  · -- ∀ t > 0, f t = k · (t · log t) + c · (t - 1)
    intro t ht
    have hg := hg_on_pos t ht
    show f t = k * (t * Real.log t) + c * (t - 1)
    have hft : f t = g t + c * (t - 1) := by
      show f t = (f t - c * (t - 1)) + c * (t - 1)
      ring
    rw [hft, hg]
  · -- f 0 = -c. Use continuity of f at 0 + the form on positives.
    have hf_tendsto : Filter.Tendsto f (nhds 0) (nhds (f 0)) :=
      hf_cont.tendsto 0
    have hform_cont : Continuous (fun t : ℝ =>
        k * (t * Real.log t) + c * (t - 1)) := by
      refine Continuous.add ?_ ?_
      · exact continuous_const.mul Real.continuous_mul_log
      · exact continuous_const.mul (continuous_id.sub continuous_const)
    have hform_at_0 :
        k * ((0 : ℝ) * Real.log 0) + c * ((0 : ℝ) - 1) = -c := by
      rw [zero_mul, mul_zero, zero_sub, mul_neg_one, zero_add]
    have hform_tendsto : Filter.Tendsto
        (fun t : ℝ => k * (t * Real.log t) + c * (t - 1))
        (nhds 0) (nhds (-c)) := by
      have := hform_cont.tendsto 0
      rwa [hform_at_0] at this
    have heq_on_pos : ∀ t : ℝ, 0 < t →
        f t = k * (t * Real.log t) + c * (t - 1) := by
      intro t ht
      have hg := hg_on_pos t ht
      have hft : f t = g t + c * (t - 1) := by
        show f t = (f t - c * (t - 1)) + c * (t - 1)
        ring
      rw [hft, hg]
    have hf_right : Filter.Tendsto f (nhdsWithin 0 (Set.Ioi 0)) (nhds (f 0)) :=
      hf_tendsto.mono_left nhdsWithin_le_nhds
    have hform_right : Filter.Tendsto
        (fun t : ℝ => k * (t * Real.log t) + c * (t - 1))
        (nhdsWithin 0 (Set.Ioi 0)) (nhds (-c)) :=
      hform_tendsto.mono_left nhdsWithin_le_nhds
    have heq_eventually : Filter.EventuallyEq
        (nhdsWithin (0 : ℝ) (Set.Ioi 0)) f
        (fun t : ℝ => k * (t * Real.log t) + c * (t - 1)) := by
      refine eventually_nhdsWithin_of_forall ?_
      intro t ht
      exact heq_on_pos t ht
    have hf_right' : Filter.Tendsto f (nhdsWithin 0 (Set.Ioi 0)) (nhds (-c)) := by
      rw [Filter.tendsto_congr' heq_eventually]
      exact hform_right
    have hne : (nhdsWithin (0 : ℝ) (Set.Ioi 0)).NeBot := nhdsGT_neBot 0
    exact tendsto_nhds_unique hf_right hf_right'

/-- **Csiszár step 5 (sign of `c₁`).**

    Given a generator `f`, convex on `[0, ∞)`, of the canonical form
    `f(t) = c₁ · t · log t + c₂ · (t − 1)` on `t > 0`, the leading
    coefficient `c₁` is non-negative.

    Proof: apply the secant-slope monotonicity for convex functions
    (`ConvexOn.slope_mono_adjacent`) at the three points
    `x = 1 < y = e < z = e²`. Since `f(1) = 0`, `f(e) = c₁·e + c₂(e-1)`,
    `f(e²) = 2c₁·e² + c₂(e²-1)`, the slope inequality
    `(f(e) - f(1))/(e - 1) ≤ (f(e²) - f(e))/(e² - e)` simplifies (after
    clearing the `c₂` cancellation) to `c₁·e ≤ c₁·(2e - 1)`, equivalently
    `0 ≤ c₁ · (e - 1)`. Since `e > 1`, this forces `c₁ ≥ 0`.

    The convexity hypothesis is essential — see `IsFDivergenceFamily`
    where convexity excludes the "anti-KL" generator `f(t) = -t·log t`
    that would otherwise satisfy additivity but yield `c₁ < 0`. -/
private lemma csiszar_c1_nonneg
    (f : ℝ → ℝ) (_hf_cont : Continuous f)
    (hf_convex : ConvexOn ℝ (Set.Ici (0 : ℝ)) f)
    (_h_additive : IsDivergenceAdditive
        (fun _ P Q => fDivergence f P Q))
    {c₁ c₂ : ℝ}
    (hf_form_pos : ∀ t : ℝ, 0 < t →
        f t = c₁ * (t * Real.log t) + c₂ * (t - 1))
    (_hf_zero : f 0 = -c₂) :
    0 ≤ c₁ := by
  -- Apply `ConvexOn.slope_mono_adjacent` at `x = 1 < y = e < z = e²`.
  set e : ℝ := Real.exp 1 with he_def
  have he_pos : 0 < e := Real.exp_pos 1
  -- `e > 1` via `Real.one_lt_exp_iff : 1 < Real.exp x ↔ 0 < x`.
  have he_gt_one : 1 < e := by
    rw [he_def, Real.one_lt_exp_iff]
    norm_num
  have he_log : Real.log e = 1 := by
    rw [he_def, Real.log_exp]
  -- Membership in `Set.Ici 0`.
  have h1_mem : (1 : ℝ) ∈ Set.Ici (0 : ℝ) := by
    rw [Set.mem_Ici]; norm_num
  have he2_mem : e ^ 2 ∈ Set.Ici (0 : ℝ) := by
    rw [Set.mem_Ici]; positivity
  -- Strict orderings `1 < e < e²`.
  have h1_lt_e : (1 : ℝ) < e := he_gt_one
  have he_lt_e2 : e < e ^ 2 := by
    have hmul : e * 1 < e * e := mul_lt_mul_of_pos_left he_gt_one he_pos
    nlinarith [hmul]
  -- Values of `f` at the three points (positivity required for `hf_form_pos`).
  have hf_at_1 : f 1 = 0 := by
    have := hf_form_pos 1 (by norm_num)
    rw [Real.log_one] at this
    linarith [this]
  have hf_at_e : f e = c₁ * e + c₂ * (e - 1) := by
    have := hf_form_pos e he_pos
    rw [he_log] at this
    linarith [this]
  have he2_pos : (0 : ℝ) < e ^ 2 := by positivity
  have hlog_e2 : Real.log (e ^ 2) = 2 := by
    rw [Real.log_pow, he_log]
    ring
  have hf_at_e2 : f (e ^ 2) = c₁ * (2 * e ^ 2) + c₂ * (e ^ 2 - 1) := by
    have h := hf_form_pos (e ^ 2) he2_pos
    rw [hlog_e2] at h
    linarith [h]
  -- The secant-slope inequality.
  have hslope :=
    hf_convex.slope_mono_adjacent h1_mem he2_mem h1_lt_e he_lt_e2
  rw [hf_at_1, hf_at_e, hf_at_e2] at hslope
  -- After substitution:
  --   (c₁ * e + c₂ * (e - 1) - 0) / (e - 1)
  --     ≤ (c₁ * (2 * e^2) + c₂ * (e^2 - 1) - (c₁ * e + c₂ * (e - 1))) / (e^2 - e)
  have he_minus_one_pos : (0 : ℝ) < e - 1 := by linarith
  have he2_minus_e_pos : (0 : ℝ) < e ^ 2 - e := by nlinarith [he_pos, he_minus_one_pos]
  -- Clear denominators via `div_le_div_iff₀`.
  rw [div_le_div_iff₀ he_minus_one_pos he2_minus_e_pos] at hslope
  -- The expanded inequality, after `ring_nf`, is equivalent to
  -- `0 ≤ c₁ * (e - 1)^2 * e`. We use `nlinarith` with the slope inequality
  -- and positivity facts to derive `0 ≤ c₁ * ((e - 1)^2 * e)`.
  have h_pos_factor : 0 < (e - 1) ^ 2 * e := by
    have : 0 < (e - 1) ^ 2 := by positivity
    exact mul_pos this he_pos
  have hkey : 0 ≤ c₁ * ((e - 1) ^ 2 * e) := by
    nlinarith [hslope, sq_nonneg (e - 1), he_pos, he_minus_one_pos,
               he2_minus_e_pos]
  exact nonneg_of_mul_nonneg_left hkey h_pos_factor

/-! ## Headline theorem -/

/-- Auxiliary: any f-divergence with generator `f(t) = c₁ · t log t + c₂ (t − 1)`
    equals `c₁ · D_KL(P‖Q)` on every pair `(P, Q)` of strictly positive
    distributions (i.e. `Q` has positive entries). The `c₂ · (t − 1)`
    summand contributes zero by `fDivergence_linear_term_vanishes`. -/
private lemma fDivergence_of_csiszar_form
    {m : ℕ} (P Q : Simplex m) (hQ : ∀ i, 0 < Q.prob i)
    (f : ℝ → ℝ) (c₁ c₂ : ℝ)
    (hf : ∀ t : ℝ, 0 < t →
        f t = c₁ * (t * Real.log t) + c₂ * (t - 1))
    (hf0 : f 0 = -c₂) :
    fDivergence f P Q = c₁ * klDivergence P Q := by
  -- Decompose each summand: when `P i = 0`, the term `f 0 = -c₂`
  -- (from `hf0`); when `P i > 0`, the form `hf` applies at
  -- `t = P.prob i / Q.prob i > 0`.
  rw [klDivergence_eq_fDivergence_tlog]
  -- Rewrite `c₁ * fDivergence (t·log t) P Q + 0 = c₁ * (∑ …)`.
  have hlin : fDivergence (fun t => c₂ * (t - 1)) P Q = 0 :=
    fDivergence_linear_term_vanishes c₂ P Q hQ
  -- We aim for:
  --   fDivergence f P Q
  --     = c₁ * fDivergence (fun t => t * Real.log t) P Q
  --     + fDivergence (fun t => c₂ * (t - 1)) P Q
  -- then close using `hlin`.
  have hsplit :
      fDivergence f P Q
        = c₁ * fDivergence (fun t => t * Real.log t) P Q
        + fDivergence (fun t => c₂ * (t - 1)) P Q := by
    unfold fDivergence
    rw [Finset.mul_sum, ← Finset.sum_add_distrib]
    refine Finset.sum_congr rfl (fun i _ => ?_)
    have hQi_pos : 0 < Q.prob i := hQ i
    have hQi : Q.prob i ≠ 0 := hQi_pos.ne'
    by_cases hPi : P.prob i = 0
    · -- P i = 0: t = 0/Q i = 0. Use `hf0 : f 0 = -c₂`.
      -- LHS: Q i * f 0 = Q i * (-c₂).
      -- RHS: c₁ * (Q i * (0 * log 0)) + Q i * (c₂ * (0 - 1))
      --    = c₁ * 0 + Q i * (-c₂) = Q i * (-c₂). ✓
      have hzero : P.prob i / Q.prob i = 0 := by
        rw [hPi]; exact zero_div _
      rw [hzero, hf0]
      show Q.prob i * -c₂
          = c₁ * (Q.prob i * (0 * Real.log 0)) + Q.prob i * (c₂ * (0 - 1))
      rw [Real.log_zero, mul_zero, mul_zero, mul_zero, zero_add]
      ring
    · -- P i > 0: apply `hf` at `t = P i / Q i > 0`.
      have hPi_pos : 0 < P.prob i := (P.nonneg i).lt_of_ne (Ne.symm hPi)
      have ht_pos : 0 < P.prob i / Q.prob i := div_pos hPi_pos hQi_pos
      rw [hf _ ht_pos]
      ring
  rw [hsplit, hlin, add_zero]

/-- **Theorem T6c (G1 autonomous — Csiszár characterisation).**

    Any continuous f-divergence on the probability simplex that is
    additive under product distributions is a non-negative scalar
    multiple of the Kullback–Leibler divergence.

    This is the PT-canonical formulation of T6c (Theorem `G1_autonomous`
    in `chapters_fr/ch05_geometry.tex`), replacing the retracted
    entropy-based statement of `T6cChencov.lean`.

    ## Proof strategy (Csiszár 1967)

    The proof is **assembled** from three private sub-lemmas:

    * `csiszar_generator_decomposition` (steps 1–3): produces
      constants `c₁, c₂` such that `f(t) = c₁ · t log t + c₂ · (t − 1)`
      for every `t > 0`. Uses `cauchy_log_equation` from
      `T6cChencov.lean`.

    * `csiszar_c1_nonneg` (step 5): forces `c₁ ≥ 0`.

    * `fDivergence_of_csiszar_form` (step 4): collapses the generator
      decomposition into `c₁ · D_KL`, using
      `fDivergence_linear_term_vanishes` for the affine summand.

    The headline theorem is sorry-free modulo these sub-lemmas; the
    sorries are now confined to fine-grained, tractable sub-tasks.

    **Statement scope.** The conclusion `D m P Q = c · klDivergence P Q`
    is restricted to pairs `(P, Q)` with `Q` strictly positive. This is
    the standard regime of absolute continuity `P ≪ Q`, and matches the
    implicit setting of Csiszár 1967 (where the f-divergence is naturally
    defined on the open simplex). The boundary case `∃ i, Q.prob i = 0`
    is governed by Mathlib's `0 · f(p/0) = 0` convention and depends on
    the chosen extension of `f` at `0`; it is not part of the
    characterisation.

    Reference: Csiszár 1967, Theorem 1; monograph
    `chapters_fr/ch05_geometry.tex`, Theorem `G1_autonomous`. -/
theorem G1_autonomous_DKL_unique
    (D : ∀ m, Simplex m → Simplex m → ℝ)
    (_h_fdiv : IsFDivergenceFamily D)
    (_h_additive : IsDivergenceAdditive D) :
    ∃ c : ℝ, 0 ≤ c ∧ ∀ m (P Q : Simplex m), (∀ i, 0 < Q.prob i) →
      D m P Q = c * klDivergence P Q := by
  -- Step A: extract the generator `f`, its continuity, convexity, and
  -- normalisation `f(1) = 0` from the f-divergence family hypothesis.
  obtain ⟨f, hf_cont, hf_convex, hf_one_zero, hf_repr⟩ := _h_fdiv
  -- Step B: transport the additivity from `D` to `fDivergence f`,
  -- preserving the positivity hypotheses on the denominators.
  have h_additive_f : IsDivergenceAdditive
      (fun m P Q => fDivergence f P Q) := by
    intro m n P₁ Q₁ P₂ Q₂ hQ₁ hQ₂
    have h1 := hf_repr (m * n) (P₁.prod P₂) (Q₁.prod Q₂)
    have h2 := hf_repr m P₁ Q₁
    have h3 := hf_repr n P₂ Q₂
    have h4 := _h_additive P₁ Q₁ P₂ Q₂ hQ₁ hQ₂
    rw [h1, h2, h3] at h4
    exact h4
  -- Step C: decompose the generator (Csiszár steps 1–3).
  obtain ⟨c₁, c₂, hf_form, hf_zero⟩ :=
    csiszar_generator_decomposition f hf_cont hf_one_zero h_additive_f
  -- Step D: pin the sign (Csiszár step 5).
  have hc₁_nonneg : 0 ≤ c₁ :=
    csiszar_c1_nonneg f hf_cont hf_convex h_additive_f hf_form hf_zero
  -- Step E: assemble. For every strictly positive `Q`, the f-divergence
  -- collapses to `c₁ · D_KL` via `fDivergence_of_csiszar_form`.
  refine ⟨c₁, hc₁_nonneg, fun m P Q hQ => ?_⟩
  rw [hf_repr m P Q]
  exact fDivergence_of_csiszar_form P Q hQ f c₁ c₂ hf_form hf_zero

end PT.Information
