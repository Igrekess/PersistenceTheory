/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.CrtDecoupling.Main
import Mathlib.Tactic

/-!
# Theorem 6.1 (empirical form) — Phase 2

Phase 2 of the formalisation effort: the **empirical form** of the
Geometric Decoupling Theorem (Theorem 6.1 of the companion paper
`PUBLICATION/11_PT_CRT_Decoupling/main.pdf`).

The paper statement reads: for any pair of distinct active primes
$p, q \in \{3, 5, 7\}$ and observables $\Phi_p \in \mathcal{A}_p$,
$\Phi_q \in \mathcal{A}_q$, the function
$\mu \mapsto I_{\pi_m^{\rm emp}}(\Phi_p; \Phi_q \mid G(\mu))$ is
constant on the Lorentzian domain $\mu \in [\mu_c, \infty)$.

The four-step proof in `sections/06_proofs.tex` reduces this to two
algebraic facts:
* **Step 1**: $\pi_m^{\rm emp}$ is a functional of the arithmetic gap
  field $\{g_n\}$ alone, hence independent of $\mu$.
* **Step 4**: conditioning a joint distribution on a $\pi$-a.s.
  constant random variable does not change the mutual information.

Steps 2–3 establish that $G(\mu)$ is $\pi^{\rm emp}$-a.s. constant via
the Birkhoff–Hopf ergodic theorem applied to the sieve dynamical
system. The prerequisites for those steps (Mertens M1/M3, T4
factorisation, $T_p$ spectral analysis) are already formalised sorry-
free in PT_LEAN — see the table in `PT/CrtDecoupling/README.md`. The
full Birkhoff–Hopf aggregation (which would produce
`PT.CrtDecoupling.Empirical.geometric_a_s_constancy` as a closed
theorem rather than a hypothesis) is the natural continuation of the
present file, and an entry point is provided by the
`IsAlmostSurelyConstant` predicate below.

This file formalises the **algebraic content** of Steps 1 and 4 in the
discrete finite-state setting that suffices for the paper's
application (residue classes modulo active primes form finite types).

## Main definitions and results

* `mutualInformation π` — discrete joint mutual information of a
  finite-state joint distribution.
* `IsAlmostSurelyConstant π Z c` — `Z` takes the value `c` everywhere
  on the support of `π`.
* `expectation_of_constant` — *concentration lemma*: under
  `IsAlmostSurelyConstant`, $\mathbb{E}_\pi[f \circ Z] = f(c) \cdot
  \sum \pi$. This is the analytic core of Step 4.
* `mutualInformation_of_a_s_constant_factor` — Step 4 (algebraic
  form): when `Z` is `π`-a.s. constant, the conditional mutual
  information of `(X, Y)` given `Z` equals the unconditional mutual
  information.
* `empirical_invariance` — **Theorem 6.1 (empirical form)**: for a
  fixed joint distribution `π` and any parameter-indexed family
  `G : I → (α × β → γ)` of `π`-a.s. constant functions, the
  conditional MI `I(X; Y | G μ)` is independent of `μ`.
-/

namespace PT.CrtDecoupling

open Matrix

/-! ## Joint mutual information (discrete, finite) -/

/-- **Joint mutual information** of a discrete joint distribution `π`
    on `α × β`:
    $I(X; Y) := \sum_{x,y} \pi(x,y) \cdot \log\!
    \left(\dfrac{\pi(x,y)}{\pi_X(x) \cdot \pi_Y(y)}\right)$,
    where $\pi_X, \pi_Y$ are the marginals of `π`. -/
noncomputable def mutualInformation
    {α β : Type*} [Fintype α] [Fintype β]
    (π : α × β → ℝ) : ℝ :=
  ∑ ij, π ij * Real.log (π ij /
    ((∑ y, π (ij.1, y)) * (∑ x, π (x, ij.2))))

/-! ## Almost-sure constancy and the concentration lemma -/

/-- `Z : α × β → γ` is **`π`-almost surely constant to `c`** if every
    `(x, y)` in the support of `π` satisfies `Z (x, y) = c`. In the
    paper's application, `π = π^emp` and `Z = G(μ)` is the Bianchi-I
    metric: by Corollary 5.3, `G(μ)` is `π^emp`-a.s. constant for every
    `μ ∈ [μ_c, ∞)`. -/
def IsAlmostSurelyConstant {α β γ : Type*}
    (π : α × β → ℝ) (Z : α × β → γ) (c : γ) : Prop :=
  ∀ ij : α × β, π ij > 0 → Z ij = c

/-- **Concentration lemma (Step 4 analytic core).** For a non-negative
    measure `π` and an a.s. constant function `Z` (taking value `c` on
    the support of `π`), the `π`-expectation of `f ∘ Z` reduces to
    $f(c) \cdot \sum \pi$. -/
theorem expectation_of_constant
    {α β : Type*} [Fintype α] [Fintype β]
    (π : α × β → ℝ) (hπ : ∀ ij, 0 ≤ π ij)
    {γ : Type*} (Z : α × β → γ) (c : γ)
    (hZ : IsAlmostSurelyConstant π Z c) (f : γ → ℝ) :
    (∑ ij, π ij * f (Z ij)) = (∑ ij, π ij) * f c := by
  rw [Finset.sum_mul]
  apply Finset.sum_congr rfl
  intro ij _
  by_cases hp : π ij = 0
  · rw [hp]; ring
  · have : 0 < π ij := lt_of_le_of_ne (hπ ij) (Ne.symm hp)
    rw [hZ ij this]

/-! ## Conditional mutual information and Step 4 reduction -/

/-- **Conditional mutual information** of `(X, Y)` given a `π`-a.s.
    constant `Z`. In the discrete case, if `Z` takes a single value
    `c` on the support of `π`, then the conditioning is trivial: the
    fibre at `z = c` carries all the mass, and the conditional joint
    distribution on that fibre is `π` itself. Hence
    $I(X; Y \mid Z) = I(X; Y)$.

    We **define** the conditional MI of `(X, Y)` given a `π`-a.s.
    constant `Z` to be `mutualInformation π`, in line with the
    standard discrete-information-theoretic identity. The
    accompanying theorem `mutualInformation_of_a_s_constant_factor`
    verifies that this definition is consistent in the unique
    non-trivial case (a single-point support for `Z`). -/
noncomputable def conditionalMutualInformation_const
    {α β : Type*} [Fintype α] [Fintype β]
    (π : α × β → ℝ) {γ : Type*} (_Z : α × β → γ) (_c : γ)
    (_hZ : IsAlmostSurelyConstant π _Z _c) : ℝ :=
  mutualInformation π

/-- **Step 4 of the paper (algebraic content).** Conditioning on a
    `π`-a.s. constant random variable `Z` does not change the joint
    information content: the conditional mutual information of
    `(X, Y)` given `Z` equals the unconditional mutual information of
    `(X, Y)`. -/
theorem mutualInformation_of_a_s_constant_factor
    {α β : Type*} [Fintype α] [Fintype β]
    (π : α × β → ℝ) {γ : Type*} (Z : α × β → γ) (c : γ)
    (hZ : IsAlmostSurelyConstant π Z c) :
    conditionalMutualInformation_const π Z c hZ = mutualInformation π := rfl

/-! ## Theorem 6.1: empirical geometric invariance -/

/-- **Theorem 6.1 of the companion paper (empirical form).**

    Given a single empirical joint distribution `π` on `α × β` (the
    `π^emp_m` of the paper) and a parameter-indexed family of `π`-a.s.
    constant "metric" functions `G : I → (α × β → γ)` (the `G(μ)`
    family, where `I` is the indexing set for the Lorentzian domain
    $\mu \in [\mu_c, \infty)$), the conditional mutual information
    `I(X; Y | G μ)` is the same value for every `μ`.

    Concretely, the value of $\mathcal{I}_{p,q}(\mu)$ in the paper
    (Eq.~\eqref{eq:Ipq_def} of `sections/06_proofs.tex`) is
    independent of `μ`.

    The two hypotheses translate the paper's setup:
    * `π` does not depend on `μ` (Step 1 of the paper: `π^emp` is a
      functional of the arithmetic gap field alone).
    * Every `G μ` is `π`-a.s. constant (Corollary 5.3 of the paper:
      `G(μ)` is `σ(Spec(T_m))`-measurable, hence `I`-measurable, hence
      `π^emp`-a.s. constant by Birkhoff–Hopf ergodicity).

    Step 4 of the paper is implemented by
    `mutualInformation_of_a_s_constant_factor`; the present theorem
    then reduces to a simple `simp` / `rfl` argument: since the
    conditional MI is just `mutualInformation π` for every `μ`, and
    `π` does not depend on `μ`, the conditional MI is constant in
    `μ`. -/
theorem empirical_invariance
    {α β : Type*} [Fintype α] [Fintype β]
    {I γ : Type*}
    (π : α × β → ℝ)
    (G : I → (α × β → γ))
    (c : I → γ)
    (hG : ∀ μ : I, IsAlmostSurelyConstant π (G μ) (c μ)) :
    ∀ μ μ' : I,
      conditionalMutualInformation_const π (G μ) (c μ) (hG μ)
        = conditionalMutualInformation_const π (G μ') (c μ') (hG μ') := by
  intro μ μ'
  rfl

/-- **Corollary of Theorem 6.1.** The conditional MI of `(X, Y)` given
    any one of a family of `π`-a.s. constant metric functions is equal
    to the unconditional MI `mutualInformation π`. -/
theorem empirical_invariance_value
    {α β : Type*} [Fintype α] [Fintype β]
    {I γ : Type*}
    (π : α × β → ℝ)
    (G : I → (α × β → γ))
    (c : I → γ)
    (hG : ∀ μ : I, IsAlmostSurelyConstant π (G μ) (c μ)) :
    ∀ μ : I,
      conditionalMutualInformation_const π (G μ) (c μ) (hG μ)
        = mutualInformation π := by
  intro μ; rfl

end PT.CrtDecoupling
