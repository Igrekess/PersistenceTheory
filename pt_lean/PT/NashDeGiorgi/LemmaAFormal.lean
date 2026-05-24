/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import Mathlib.Tactic
import Mathlib.Algebra.Ring.GeomSum
import Mathlib.Analysis.SpecialFunctions.Pow.Real

/-!
# Lemma A — analytical proof of $\gamma_p$ monotonicity in $\mu$

This module formalises the **algebraic core** of Lemma A of
`PT_NASH_DEGIORGI/notes/15_lemmaA_formal_proof.md`:

> **Lemma A.** For every integer $p \geq 2$ and every $q \in (0, 1)$,
> $\gamma_p(q)$ is strictly increasing in $q$. By the increasing
> diffeomorphism $\mu \leftrightarrow q = 1 - 2/\mu$, this is equivalent
> to $\gamma_p(\mu)$ strictly increasing in $\mu$ on $(2, \infty)$.

The proof in note 15 reduces to the **polynomial positivity identity**
$$p - S_p(q) > 0 \quad \text{for } q \in (0, 1), p \geq 2,$$
where $S_p(q) = \sum_{k=0}^{p-1} q^k$.

This module formalises:

1. **`Sp_def`**: $S_p(q)$ definition.
2. **`one_sub_pow_eq_one_sub_mul_Sp`**: the identity $1 - q^p = (1-q) S_p(q)$
   (Mathlib's `geom_sum_eq` repackaged).
3. **`Sp_pos`**: $S_p(q) > 0$ for $q \in (0, 1)$, $p \geq 1$.
4. **`p_minus_Sp_pos`**: $p - S_p(q) > 0$ for $q \in (0, 1)$, $p \geq 2$
   (the **key positivity lemma** of Theorem 3.1 in note 15).
5. **`AofQ_pos`**: the quantity $A(q) = (p - S_p(q))/(q(1-q) S_p(q))$ is
   strictly positive on $(0, 1)$ for $p \geq 2$ (combining 1-4).

The full derivative-positivity argument (Theorem 4.1 of note 15)
requires computing $d \gamma_p / d q$ via chain rule and showing it
equals $\gamma_p \cdot (A(q) + B(q)) > 0$. That derivative computation
is a routine Mathlib application of `HasDerivAt.mul`, `.div`, `.pow`,
etc., spanning ~100 lines — left as a TODO with a placeholder axiom
that has the correct statement (no longer the `True` placeholder of
`GammaAtMu.lean`).

## Reference

`PT_PROJECTS/PT_NASH_DEGIORGI/notes/15_lemmaA_formal_proof.md`
(361 lines, complete paper-level proof).
-/

namespace PT.NashDeGiorgi.LemmaAFormal

open Finset

/-! ### The polynomial `S_p(q) = sum q^k for k = 0..p-1` -/

/-- $S_p(q) = \sum_{k=0}^{p-1} q^k$. For $q \neq 1$, equals $(1 - q^p)/(1-q)$. -/
noncomputable def Sp (p : ℕ) (q : ℝ) : ℝ := (range p).sum (fun k => q ^ k)

/-- Sanity check: `S_3(q) = 1 + q + q^2`. -/
example (q : ℝ) : Sp 3 q = 1 + q + q^2 := by
  unfold Sp
  simp [Finset.sum_range_succ]

/-! ### Key identity `1 - q^p = (1-q) * S_p(q)` -/

/-- The geometric sum identity: $1 - q^p = (1-q) \sum_{k=0}^{p-1} q^k$. -/
theorem one_sub_pow_eq_one_sub_mul_Sp (p : ℕ) (q : ℝ) :
    1 - q ^ p = (1 - q) * Sp p q := by
  unfold Sp
  -- Mathlib: `geom_sum_mul : (∑ i ∈ range n, x^i) * (x - 1) = x^n - 1`
  have h := geom_sum_mul q p
  linarith

/-! ### Positivity of `S_p(q)` -/

/-- $S_p(q) > 0$ for $q > 0$ and $p \geq 1$. -/
theorem Sp_pos (p : ℕ) (hp : 1 ≤ p) {q : ℝ} (hq : 0 < q) :
    0 < Sp p q := by
  unfold Sp
  -- Every term q^k ≥ 0; the k=0 term is 1 > 0; sum ≥ 1.
  have h1 : (0 : ℝ) < 1 + (range p).sum (fun k => q ^ k) - 1 ∨
            (1 : ℝ) ≤ (range p).sum (fun k => q ^ k) := by
    right
    -- Sum includes k=0 contribution of q^0 = 1.
    have h_zero : (1 : ℝ) = q ^ 0 := by simp
    have : q^0 ≤ (range p).sum (fun k => q ^ k) := by
      apply Finset.single_le_sum
        (f := fun k => q^k) (s := range p)
      · intro i _
        exact pow_nonneg hq.le i
      · exact mem_range.mpr (Nat.lt_of_lt_of_le Nat.zero_lt_one hp)
    simpa using this
  cases h1 with
  | inl h => linarith
  | inr h => linarith

/-! ### The key positivity: `p - S_p(q) > 0` for `q ∈ (0,1)`, `p ≥ 2` -/

/-- For `q ∈ [0,1]` and any `k ≥ 1`, `q^k ≤ 1`. -/
theorem pow_le_one_of_mem_Icc {q : ℝ} (hq_pos : 0 ≤ q) (hq_lt : q ≤ 1) (k : ℕ) :
    q ^ k ≤ 1 := by
  exact pow_le_one₀ hq_pos hq_lt

/-- For `q ∈ (0,1)` and `k ≥ 1`, `q^k < 1` strict. -/
theorem pow_lt_one_of_pos_lt_one {q : ℝ} (hq_pos : 0 < q) (hq_lt : q < 1)
    (k : ℕ) (hk : 1 ≤ k) : q ^ k < 1 := by
  exact pow_lt_one₀ hq_pos.le hq_lt (by omega)

/-- **Theorem (key positivity).** For $p \geq 2$ and $q \in (0, 1)$,
    $p - S_p(q) > 0$.

    This is Equation (3.4) of [note 15](PT_NASH_DEGIORGI/notes/15) :
    $p - S_p(q) = \sum_{k=0}^{p-1}(1 - q^k) = \sum_{k=1}^{p-1}(1 - q^k) > 0$.
    The $k=0$ term contributes $1 - q^0 = 0$, leaving a strictly positive
    sum from $k=1$ to $p-1$ since each $1 - q^k > 0$ for $q < 1, k \geq 1$. -/
theorem p_minus_Sp_pos (p : ℕ) (hp : 2 ≤ p) {q : ℝ}
    (hq_pos : 0 < q) (hq_lt : q < 1) :
    0 < (p : ℝ) - Sp p q := by
  unfold Sp
  -- Rewrite p - sum q^k = sum (1 - q^k) over range p
  have h_eq : (p : ℝ) - (range p).sum (fun k => q ^ k)
            = (range p).sum (fun k => 1 - q ^ k) := by
    rw [Finset.sum_sub_distrib]
    simp
  rw [h_eq]
  -- The sum over range p of (1 - q^k). Split off k=0.
  rw [show range p = insert 0 (Finset.Ioo 0 p) by
    ext k
    simp only [mem_insert, mem_range, Finset.mem_Ioo]
    omega]
  rw [Finset.sum_insert (by simp)]
  simp only [pow_zero, sub_self, zero_add]
  -- Now show sum over Ioo 0 p of (1 - q^k) > 0
  -- Since p ≥ 2, Ioo 0 p contains at least k = 1
  apply Finset.sum_pos
  · intro k hk
    rw [Finset.mem_Ioo] at hk
    have h1 : q ^ k < 1 := pow_lt_one_of_pos_lt_one hq_pos hq_lt k hk.1
    linarith
  · refine ⟨1, ?_⟩
    rw [Finset.mem_Ioo]
    exact ⟨Nat.zero_lt_one, hp⟩

/-! ### The quantity $A(q) = (p - S_p(q)) / (q (1-q) S_p(q))$

This is the simplified form of $\frac{d \ln \gamma_p}{dq}$'s "dangerous"
part (Theorem 3.1 of note 15). Its strict positivity on $(0, 1)$ for
$p \geq 2$ is the algebraic content of Lemma A. -/

/-- $A(q) = (p - S_p(q)) / (q (1-q) S_p(q))$, the simplified form (3.3). -/
noncomputable def AofQ (p : ℕ) (q : ℝ) : ℝ :=
  ((p : ℝ) - Sp p q) / (q * (1 - q) * Sp p q)

/-- **Theorem (positivity of $A$).** For $p \geq 2$ and $q \in (0, 1)$,
    $A(q) > 0$. This is the algebraic core of Lemma A from note 15 §3. -/
theorem AofQ_pos (p : ℕ) (hp : 2 ≤ p) {q : ℝ}
    (hq_pos : 0 < q) (hq_lt : q < 1) :
    0 < AofQ p q := by
  unfold AofQ
  have h_num : 0 < (p : ℝ) - Sp p q := p_minus_Sp_pos p hp hq_pos hq_lt
  have h_q : 0 < q := hq_pos
  have h_1mq : 0 < 1 - q := by linarith
  have h_Sp : 0 < Sp p q := Sp_pos p (by omega) hq_pos
  have h_denom : 0 < q * (1 - q) * Sp p q := by positivity
  exact div_pos h_num h_denom

/-! ### $\delta_p$ in terms of $q$ -/

/-- $\delta_p(q) = (1 - q^p) / p$. -/
noncomputable def deltaOfQ (p : ℕ) (q : ℝ) : ℝ := (1 - q ^ p) / (p : ℝ)

/-- $\delta_p \in (0, 1)$ for $q \in (0, 1)$ and $p \geq 1$. -/
theorem deltaOfQ_mem_Ioo (p : ℕ) (hp : 1 ≤ p) {q : ℝ}
    (hq_pos : 0 < q) (hq_lt : q < 1) :
    deltaOfQ p q ∈ Set.Ioo (0 : ℝ) 1 := by
  unfold deltaOfQ
  refine ⟨?_, ?_⟩
  · -- 0 < (1 - q^p) / p : numerator > 0 (since q^p < 1) and p > 0
    have hp_pos : (0 : ℝ) < (p : ℝ) := by exact_mod_cast Nat.lt_of_lt_of_le Nat.zero_lt_one hp
    have h_qp : q^p < 1 := pow_lt_one_of_pos_lt_one hq_pos hq_lt p hp
    have : 0 < 1 - q^p := by linarith
    exact div_pos this hp_pos
  · -- (1 - q^p) / p < 1 : equivalent to 1 - q^p < p
    have hp_pos : (0 : ℝ) < (p : ℝ) := by exact_mod_cast Nat.lt_of_lt_of_le Nat.zero_lt_one hp
    rw [div_lt_one hp_pos]
    have hp_ge_one : (1 : ℝ) ≤ (p : ℝ) := by exact_mod_cast hp
    have h_qp_pos : 0 < q^p := pow_pos hq_pos p
    -- 1 - q^p < 1 ≤ p
    linarith

/-! ### Lemma A (full statement)

The strict monotonicity of $\gamma_p$ in $q$ on $(0, 1)$ (equivalently
$\gamma_p$ in $\mu$ on $(2, \infty)$) follows by combining `AofQ_pos`
above with $B(q) > 0$ (trivial) and the derivative identity
$d \ln \gamma_p / dq = A(q) + B(q)$ (note 15 §2).

The derivative computation is a routine but lengthy chain-rule
application in Mathlib (~100 lines). For now we record the statement
with the algebraic positivity already proved, and leave the explicit
HasDerivAt computation as a TODO. -/

/-- The function $\gamma_p(q)$ for variable $q \in (0, 1)$. -/
noncomputable def gammaOfQ (p : ℕ) (q : ℝ) : ℝ :=
  (2 * q^(p-1) * (1 - deltaOfQ p q) * (1 - q)) /
  (deltaOfQ p q * (2 - deltaOfQ p q))

/-- The function $\gamma_p(\mu)$ for variable $\mu \in (2, \infty)$. -/
noncomputable def gammaOfMu (p : ℕ) (μ : ℝ) : ℝ := gammaOfQ p (1 - 2 / μ)

/-- **Lemma A (formal statement).** For every integer $p \geq 2$ and
    every $\mu_1 < \mu_2$ in $(2, \infty)$, $\gamma_p(\mu_1) < \gamma_p(\mu_2)$.

    Proof in Lean: combine `AofQ_pos` (this file, algebraic core) with
    the derivative computation $d \ln \gamma_p / dq = A(q) + B(q)$
    (note 15 §2, routine Mathlib `HasDerivAt` chain) and Mathlib's
    `strictMonoOn_of_deriv_pos` from `Analysis.Calculus.Deriv.MeanValue`.

    Paper-level proof: [note 15](PT_NASH_DEGIORGI/notes/15) §4.

    Status: algebraic core (`AofQ_pos`) and statement formalised; full
    derivative+monotonicity formalisation pending Lean port of the
    chain-rule computation (~100 lines, mechanical). -/
axiom lemmaA_strictMono (p : ℕ) (hp : 2 ≤ p) :
    ∀ {μ₁ μ₂ : ℝ}, 2 < μ₁ → μ₁ < μ₂ → gammaOfMu p μ₁ < gammaOfMu p μ₂

/-- **Corollary A1.** For $K = [\mu_-, \mu_+]$ compact in $(2, \infty)$,
    $\sup_{\mu \in K} \gamma_p(\mu) = \gamma_p(\mu_+)$. -/
theorem corA1_sup_at_right (p : ℕ) (hp : 2 ≤ p) {μ_minus μ_plus : ℝ}
    (h_lo : 2 < μ_minus) (_h_le : μ_minus ≤ μ_plus) :
    ∀ μ ∈ Set.Icc μ_minus μ_plus, gammaOfMu p μ ≤ gammaOfMu p μ_plus := by
  intro μ hμ
  rcases eq_or_lt_of_le hμ.2 with h_eq | h_lt
  · rw [h_eq]
  · exact le_of_lt (lemmaA_strictMono p hp (lt_of_lt_of_le h_lo hμ.1) h_lt)

end PT.NashDeGiorgi.LemmaAFormal
