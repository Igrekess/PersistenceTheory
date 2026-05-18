/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import Mathlib.Analysis.SpecialFunctions.Exp
import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.Tactic

/-!
# Bifurcation vertex-edge — two canonical PT parameters

**Statement (paper-level, Ch03 §"Bifurcation sommet/arête").**
Two canonical branch parameters emerge from the conservation laws at the
persistence fixed point `μ* = 15`:

* **Vertex parameter** `q_+ = 1 - 2/μ*` (with the L0 mean convention
  `μ_+ = μ*/2 = 15/2`, equivalently `q_+ = 1 - 1/μ_+ = 13/15`).
  Characterised by the **Lagrange / max-entropy** condition: among
  distributions of given mean `μ_+`, the geometric distribution
  `p_k = (1-q_+) q_+^{k-1}` maximises Shannon entropy.

* **Edge parameter** `q_- = exp(-1/μ*) = exp(-1/15)`.
  Characterised by the **Boltzmann / Gibbs self-consistency** condition:
  the thermal weight `q_- = exp(-β)` with inverse temperature `β = 1/μ*`,
  derived from the soliton constraint on the (continuous, edge-labelled)
  cascade.

**Bifurcation theorem.** `q_+ ≠ q_-` at the PT fixed point `μ* = 15`. The
two parameters are genuinely distinct: PT exhibits a vertex/edge bifurcation
with no degenerate limit where the two collapse onto one another.

This module formalises the algebraic core:

1. The **two equations** that each parameter satisfies (the *vertex
   equation* `q_+ μ_+ = μ_+ - 1` and the *edge equation* `log q_- = -1/μ*`).
2. The **separation theorem** `q_+ ≠ q_-` at `μ* = 15`.
3. The **quantitative gap** `q_- - q_+ ≥ 1/15` via the Taylor lower bound
   `exp(x) ≥ 1 + x`.

## Reference

Monograph Chapter 3 §"Bifurcation vertex/edge", `\label{prop:bifurcation}`.
Audit row #15 (EASY, 1 session): "Two routes: Lagrange, self-consistent
Boltzmann."
Cross-refs: `PT.Information.L0MaxEntropy` (max-entropy route);
`PT.EML.QSheffer` (EML depth comparison route).

## Strategy

Pure algebra on `ℝ`. The vertex equation is `(1 - q_+) μ_+ = 1`, an identity
checked by `norm_num` at the PT specialisation. The edge equation is
`log q_- = -1/15`, immediate from `Real.log_exp`. The separation uses the
strict Taylor bound `Real.add_one_le_exp` applied at `x = -1/15`:
`q_- = exp(-1/15) ≥ 1 + (-1/15) = 14/15 > 13/15 = q_+`.
-/

namespace PT.Conservation

open Real

/-! ### Definition of the two PT parameters at the fixed point `μ* = 15` -/

/-- The PT fixed point. -/
def muStarBif : ℝ := 15

/-- The L0 mean convention (half of `μ*`, governing the jump-distribution lattice). -/
noncomputable def muPlusBif : ℝ := muStarBif / 2  -- = 15/2

/-- The vertex parameter `q_+`. -/
noncomputable def qPlusBif : ℝ := 1 - 1 / muPlusBif  -- = 1 - 2/15 = 13/15

/-- The edge parameter `q_-`. -/
noncomputable def qMinusBif : ℝ := Real.exp (-(1 / muStarBif))  -- = exp(-1/15)

@[simp] lemma muStarBif_eq : muStarBif = 15 := rfl

@[simp] lemma muPlusBif_eq : muPlusBif = 15 / 2 := by
  unfold muPlusBif muStarBif; norm_num

/-! ### Explicit value of the vertex parameter `q_+` -/

/-- The vertex parameter `q_+ = 13/15`. -/
theorem qPlus_value : qPlusBif = 13 / 15 := by
  unfold qPlusBif muPlusBif muStarBif
  norm_num

/-- `q_+` lies in the open unit interval. -/
theorem qPlus_lt_one : qPlusBif < 1 := by
  rw [qPlus_value]; norm_num

/-- `q_+` is non-negative. -/
theorem qPlus_nonneg : 0 ≤ qPlusBif := by
  rw [qPlus_value]; norm_num

/-! ### The vertex equation: `(1 - q_+) μ_+ = 1` (Lagrange / max-entropy) -/

/-- **Vertex equation.** `q_+` satisfies the max-entropy self-consistency
    relation `(1 - q_+) μ_+ = 1`, equivalent to `q_+ = 1 - 1/μ_+`.

    This is the Lagrange / KKT condition for the geometric distribution
    `p_k = (1-q) q^{k-1}` to maximise Shannon entropy on `ℕ⁺` under a
    prescribed-mean constraint `∑ k p_k = μ_+`. -/
theorem qPlus_vertex_equation : (1 - qPlusBif) * muPlusBif = 1 := by
  unfold qPlusBif muPlusBif muStarBif
  norm_num

/-! ### The edge equation: `log q_- = -1/μ*` (Boltzmann / Gibbs) -/

/-- The edge parameter `q_-` is strictly positive. -/
theorem qMinus_pos : 0 < qMinusBif := by
  unfold qMinusBif
  exact Real.exp_pos _

/-- **Edge equation.** `q_-` satisfies the Gibbs self-consistency relation
    `log q_- = -1/μ*`, equivalent to `q_- = exp(-1/μ*)`.

    This is the Boltzmann / soliton condition for the continuous (edge-
    labelled) cascade: each edge carries a thermal weight `exp(-β)` with
    inverse temperature `β = 1/μ*`. -/
theorem qMinus_edge_equation : Real.log qMinusBif = -(1 / muStarBif) := by
  unfold qMinusBif
  exact Real.log_exp _

/-! ### Lower bound on `q_-` via the Taylor inequality `exp(x) ≥ 1 + x` -/

/-- `q_- ≥ 14/15`: the linear Taylor lower bound at `x = -1/15`. -/
theorem qMinus_ge_14_15 : (14 : ℝ) / 15 ≤ qMinusBif := by
  unfold qMinusBif muStarBif
  have h := Real.add_one_le_exp (-((1 : ℝ) / 15))
  -- h : 1 + (-1/15) ≤ exp(-1/15)
  have habs : 1 + -((1 : ℝ) / 15) = 14 / 15 := by norm_num
  linarith [h]

/-! ### Separation theorem: `q_- ≠ q_+` at `μ* = 15` -/

/-- **Bifurcation separation (gap form).** The gap `q_- - q_+` is at least
    `1/15` at `μ* = 15`. (Equivalently `q_- ≥ 14/15 > 13/15 = q_+`.) -/
theorem bifurcation_gap : qPlusBif + 1 / 15 ≤ qMinusBif := by
  rw [qPlus_value]
  have h := qMinus_ge_14_15
  have heq : (13 : ℝ) / 15 + 1 / 15 = 14 / 15 := by norm_num
  linarith [heq]

/-- **Bifurcation separation (strict).** `q_+ < q_-` at the PT fixed point. -/
theorem qPlus_lt_qMinus : qPlusBif < qMinusBif := by
  have h := bifurcation_gap
  linarith

/-- **Bifurcation theorem — `q_+ ≠ q_-`.** The vertex and edge parameters are
    genuinely distinct at `μ* = 15`; no degenerate limit collapses them. -/
theorem bifurcation_nondegenerate : qPlusBif ≠ qMinusBif :=
  ne_of_lt qPlus_lt_qMinus

/-! ### Headline: the two canonical parameters and their equations -/

/-- **Bifurcation headline.** Aggregated statement: at `μ* = 15`, PT carries
    two canonical branch parameters

    * `q_+ = 13/15` (vertex) satisfying the Lagrange / max-entropy
      equation `(1 - q_+) μ_+ = 1`;
    * `q_- = exp(-1/15)` (edge) satisfying the Gibbs / Boltzmann equation
      `log q_- = -1/μ*`;

    and the two are strictly separated (`q_+ < q_-`).

    Equivalently, the conservation flow on the persistence sieve undergoes a
    **vertex/edge bifurcation** at `μ*`: the two attractors do not coalesce.
    -/
theorem bifurcation_vertex_edge :
    qPlusBif = 13 / 15
    ∧ (1 - qPlusBif) * muPlusBif = 1
    ∧ Real.log qMinusBif = -(1 / muStarBif)
    ∧ qPlusBif < qMinusBif :=
  ⟨qPlus_value, qPlus_vertex_equation, qMinus_edge_equation, qPlus_lt_qMinus⟩

end PT.Conservation
