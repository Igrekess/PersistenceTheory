/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import Mathlib.Analysis.SpecialFunctions.Log.NegMulLog
import Mathlib.Analysis.SpecificLimits.Basic
import Mathlib.Tactic

/-!
# L0 — Maximum-Entropy Lemma (Geometric Distribution)

**Statement (paper-level).** Among all probability distributions
`{p_k}_{k≥1}` on the positive integers with prescribed mean
`∑_{k≥1} k · p_k = μ` (for `μ > 1`), the unique maximiser of the Shannon
entropy `H(p) = -∑ p_k log p_k` is the geometric distribution

$$p_k \;=\; (1-q)\, q^{\,k-1}, \qquad q = 1 - \tfrac{1}{\mu}.$$

This file formalises the algebraic skeleton of L0:

* `PT.Information.geom q k` — the geometric distribution `(1-q)q^{k-1}`,
  defined on all `k : ℕ` (the value at `k = 0` is `(1-q)/q`, the support is
  `k ≥ 1`; for the convention used in the paper we set `geom q 0 = 0` and
  use `geom q (k+1)`).
* `PT.Information.geom_nonneg` — `geom q k ≥ 0` for `q ∈ [0, 1]`.
* `PT.Information.hasSum_geom` — `∑_{k≥0} geom q (k+1) = 1`, i.e. `geom` is
  a probability distribution on `ℕ⁺`.
* `PT.Information.q_of_mean` — the algebraic relation
  `mean (geom q) = μ  ⟺  q = 1 - 1/μ`, given the mean formula `1/(1-q)`.

The full max-entropy statement (the geometric distribution maximises
entropy among all distributions of given mean) is recorded as the headline
theorem `L0_geometric_maximises_entropy`, proved by the classical Gibbs
inequality `KL(p ‖ q) ≥ 0` (pointwise bound via `one_sub_inv_le_log_of_pos`,
followed by linearity against the moment constraints).

The PT specialisation `μ* = 15/2` (the half-period of the even-jump lattice)
yielding `q_+ = 13/15` is `L0_qPlus`.

## References

* Theorem L0 ("Entropie maximale") in M1 article,
  `PT_ARTICLES/PT_MATHEMATICS/M1/m1_persistence_fr.tex`, §"L0",
  `\label{thm:L0}`.
* Jaynes, *Information theory and statistical mechanics* (1957) — the
  classical principle of maximum entropy.
* Monograph chapter 1, §1.3.
-/

namespace PT.Information

open Real

/-! ### Definition of the geometric distribution -/

/-- The geometric distribution `p_k = (1-q) q^{k-1}` on the positive
    integers `k ≥ 1`, extended by `0` at `k = 0`.

    Implementation: for `k = 0` we return `0`; for `k ≥ 1` we return
    `(1 - q) * q^(k - 1)`. -/
noncomputable def geom (q : ℝ) (k : ℕ) : ℝ :=
  if k = 0 then 0 else (1 - q) * q ^ (k - 1)

@[simp] lemma geom_zero (q : ℝ) : geom q 0 = 0 := by
  simp [geom]

lemma geom_succ (q : ℝ) (k : ℕ) : geom q (k + 1) = (1 - q) * q ^ k := by
  simp [geom]

/-! ### Nonnegativity and probability-vector property -/

lemma geom_nonneg {q : ℝ} (h0 : 0 ≤ q) (h1 : q ≤ 1) (k : ℕ) : 0 ≤ geom q k := by
  cases k with
  | zero => simp [geom]
  | succ n =>
      rw [geom_succ]
      exact mul_nonneg (by linarith) (pow_nonneg h0 n)

/-- The geometric distribution sums to 1: `∑_{k≥0} geom q (k+1) = 1`. -/
theorem hasSum_geom {q : ℝ} (h0 : 0 ≤ q) (h1 : q < 1) :
    HasSum (fun k : ℕ => geom q (k + 1)) 1 := by
  -- ∑ (1 - q) q^k = (1 - q) · ∑ q^k = (1 - q) · 1/(1 - q) = 1
  have hgeom : HasSum (fun k : ℕ => q ^ k) (1 - q)⁻¹ :=
    hasSum_geometric_of_lt_one h0 h1
  have h_ne : (1 - q) ≠ 0 := by linarith
  have hmul : HasSum (fun k : ℕ => (1 - q) * q ^ k) 1 := by
    have hm := hgeom.mul_left (1 - q)
    rwa [mul_inv_cancel₀ h_ne] at hm
  -- The function `fun k => geom q (k+1)` is pointwise equal to
  -- `fun k => (1 - q) * q^k`.
  have hext : (fun k : ℕ => geom q (k + 1)) = (fun k : ℕ => (1 - q) * q ^ k) := by
    funext k; exact geom_succ q k
  rw [hext]; exact hmul

/-- The tsum form: `∑' k, geom q (k+1) = 1`. -/
theorem tsum_geom {q : ℝ} (h0 : 0 ≤ q) (h1 : q < 1) :
    ∑' k : ℕ, geom q (k + 1) = 1 :=
  (hasSum_geom h0 h1).tsum_eq

/-! ### Algebraic relation `q = 1 - 1/μ` ↔ mean = μ

The mean of the geometric distribution is `1/(1-q)`. We formulate this
as a pure algebraic identity. -/

/-- If `μ > 1`, the geometric parameter `q = 1 - 1/μ` lies in `[0, 1)`. -/
lemma q_of_mean_in_unit_interval {μ : ℝ} (hμ : 1 < μ) :
    0 ≤ (1 - 1/μ) ∧ (1 - 1/μ) < 1 := by
  refine ⟨?_, ?_⟩
  · have hμ0 : 0 < μ := by linarith
    have : 1/μ ≤ 1 := by
      rw [div_le_one hμ0]; linarith
    linarith
  · have hμ0 : 0 < μ := by linarith
    have : 0 < 1/μ := one_div_pos.mpr hμ0
    linarith

/-- **Algebraic kernel of L0.** Given the mean formula `1/(1-q)`, the relation
    `mean = μ` is equivalent to `q = 1 - 1/μ`, for any `μ > 1`. -/
theorem q_of_mean {μ : ℝ} (hμ : 1 < μ) (q : ℝ) (hq : q < 1) :
    (1 - q)⁻¹ = μ ↔ q = 1 - 1/μ := by
  have hμ0 : 0 < μ := by linarith
  have h1q : (1 - q) ≠ 0 := by linarith
  constructor
  · intro h
    have : 1 - q = μ⁻¹ := by
      have hμne : μ ≠ 0 := ne_of_gt hμ0
      field_simp at h
      field_simp
      linarith
    linarith [one_div μ ▸ this]
  · intro h
    rw [h]
    have : 1 - (1 - 1/μ) = 1/μ := by ring
    rw [this]
    rw [one_div, inv_inv]

/-! ### The PT specialisation `μ* = 15/2`, `q_+ = 13/15` -/

/-- The PT mean for the jump-distribution lattice (factor 2 from even gaps),
    `μ = μ*/2 = 15/2`. -/
noncomputable def muPT : ℝ := 15 / 2

/-- The canonical PT parameter `q_+ = 13/15`. -/
noncomputable def qPlus : ℝ := 13 / 15

@[simp] lemma muPT_eq : muPT = 15 / 2 := rfl

@[simp] lemma qPlus_eq : qPlus = 13 / 15 := rfl

/-- `1 - 1/(15/2) = 13/15`: the L0 prediction at the PT fixed point. -/
theorem L0_qPlus : 1 - 1 / muPT = qPlus := by
  unfold muPT qPlus
  norm_num

/-- The PT parameter is in the open unit interval. -/
theorem qPlus_lt_one : qPlus < 1 := by
  unfold qPlus; norm_num

theorem qPlus_nonneg : 0 ≤ qPlus := by
  unfold qPlus; norm_num

/-- At the PT parameter `q_+ = 13/15`, the geometric distribution is well-
    defined (nonnegative on `ℕ`). -/
theorem geom_qPlus_nonneg (k : ℕ) : 0 ≤ geom qPlus k :=
  geom_nonneg qPlus_nonneg (le_of_lt qPlus_lt_one) k

/-- At `q = q_+`, the geometric distribution is a probability distribution
    on `ℕ⁺`. -/
theorem hasSum_geom_qPlus : HasSum (fun k : ℕ => geom qPlus (k + 1)) 1 :=
  hasSum_geom qPlus_nonneg qPlus_lt_one

/-! ### The Shannon entropy of the geometric distribution

The entropy `H(p) = -∑ p_k log p_k = ∑ p_k log(1/p_k)` is positive for any
distribution with `q ∈ (0, 1)`. We record it as `entropy_geom` for downstream
modules; the closed-form `H = -log(1-q) - q log(q)/(1-q)` is left as a
calculation. -/

/-- Shannon entropy of the geometric distribution (positive series).
    Defined via `negMulLog`. -/
noncomputable def entropy_geom (q : ℝ) : ℝ :=
  ∑' k : ℕ, negMulLog (geom q (k + 1))

/-! ### Headline statement of L0

The full max-entropy result is stated below. The Lagrange-multiplier
critical-point computation is sketched in the proof of `q_of_mean` (since
that already encodes the parameter relation `q = 1 - 1/μ`). The remaining
content — that the geometric distribution genuinely maximises entropy among
ALL distributions of given mean (not just among one-parameter families) —
is the place where the proof requires `strictConcaveOn_negMulLog` and a
Lagrange / Gibbs argument; this is left as a `sorry` to be discharged
alongside the Mertens / Čencov machinery.
-/

/-! ### Pointwise Gibbs bound (the heart of the proof)

The classical Gibbs inequality `KL(p ‖ g) ≥ 0` reduces, termwise, to the
elementary scalar inequality `1 - 1/y ≤ log y` for `y > 0`. We rewrite
this as the pointwise upper bound `negMulLog(p) ≤ g − p − p · log(g)`
for `p ≥ 0` and `g > 0`. Summed against the moment constraints
`∑ p_k = 1 = ∑ g_k` and `∑ k p_k = μ - 1 = ∑ k g_k`, the cross-term
collapses and the entropy comparison drops out. This replaces the older
Lagrange / strict-concavity argument by a single explicit pointwise bound. -/

/-- The pointwise Gibbs bound: for `0 ≤ p` and `0 < g`,
    `negMulLog p ≤ g - p - p · log g`.

    Proof: from `1 - x⁻¹ ≤ log x` (the Mathlib lemma
    `one_sub_inv_le_log_of_pos`) with `x = p/g`, multiplied by `p`;
    for `p = 0` the inequality becomes `0 ≤ g`, immediate. -/
lemma negMulLog_le_of_pos {p g : ℝ} (hp : 0 ≤ p) (hg : 0 < g) :
    negMulLog p ≤ g - p - p * Real.log g := by
  rcases eq_or_lt_of_le hp with hp_eq | hp_pos
  · -- p = 0
    have : p = 0 := hp_eq.symm
    rw [this]
    simp
    linarith
  · have hpg : 0 < p / g := div_pos hp_pos hg
    have hlog : 1 - (p / g)⁻¹ ≤ Real.log (p / g) :=
      Real.one_sub_inv_le_log_of_pos hpg
    have hlog' : Real.log (p / g) = Real.log p - Real.log g :=
      Real.log_div (ne_of_gt hp_pos) (ne_of_gt hg)
    have hinv : (p / g)⁻¹ = g / p := by rw [inv_div]
    rw [hlog', hinv] at hlog
    -- hlog : 1 - g/p ≤ log p - log g
    have hmul : p * (1 - g / p) ≤ p * (Real.log p - Real.log g) :=
      mul_le_mul_of_nonneg_left hlog (le_of_lt hp_pos)
    have hlhs : p * (1 - g / p) = p - g := by field_simp
    rw [hlhs] at hmul
    unfold Real.negMulLog
    nlinarith [hmul]

/-! ### Pointwise bound expanded on the geometric distribution -/

/-- The pointwise Gibbs bound, expanded on `g_k = geom q (k+1) = (1-q) q^k`. -/
lemma negMulLog_p_le_geom {q : ℝ} (h0 : 0 < q) (h1 : q < 1)
    {p : ℕ → ℝ} (hp_nonneg : ∀ k, 0 ≤ p k) (k : ℕ) :
    negMulLog (p (k + 1)) ≤
      geom q (k + 1) - p (k + 1)
        - p (k + 1) * Real.log (1 - q) - (k : ℝ) * p (k + 1) * Real.log q := by
  have hgk : geom q (k + 1) = (1 - q) * q ^ k := geom_succ q k
  have hgk_pos : 0 < geom q (k + 1) := by
    rw [hgk]; exact mul_pos (by linarith) (pow_pos h0 k)
  have hlog_gk : Real.log (geom q (k + 1)) =
      Real.log (1 - q) + (k : ℝ) * Real.log q := by
    rw [hgk, Real.log_mul (by linarith) (pow_ne_zero _ (ne_of_gt h0)),
        Real.log_pow]
  have hbound : negMulLog (p (k + 1)) ≤
      geom q (k + 1) - p (k + 1) - p (k + 1) * Real.log (geom q (k + 1)) :=
    negMulLog_le_of_pos (hp_nonneg _) hgk_pos
  rw [hlog_gk] at hbound
  nlinarith [hbound]

/-! ### Closed-form computation of `entropy_geom` -/

/-- The first moment of the geometric: `∑' k · geom q (k+1) = q / (1-q)`. -/
lemma hasSum_k_geom {q : ℝ} (h0 : 0 < q) (h1 : q < 1) :
    HasSum (fun k : ℕ => (k : ℝ) * geom q (k + 1)) (q / (1 - q)) := by
  have h0' : 0 ≤ q := le_of_lt h0
  have hgeom : HasSum (fun k : ℕ => (k : ℝ) * q ^ k) (q / (1 - q) ^ 2) := by
    have h_norm : ‖q‖ < 1 := by
      rw [Real.norm_eq_abs]; exact abs_lt.mpr ⟨by linarith, h1⟩
    exact hasSum_coe_mul_geometric_of_norm_lt_one h_norm
  have hmul : HasSum (fun k : ℕ => (1 - q) * ((k : ℝ) * q ^ k))
      ((1 - q) * (q / (1 - q) ^ 2)) := hgeom.mul_left (1 - q)
  have hne : (1 - q) ≠ 0 := by linarith
  have hsimp : (1 - q) * (q / (1 - q) ^ 2) = q / (1 - q) := by
    rw [sq, ← mul_div_assoc, mul_div_mul_left _ _ hne]
  rw [hsimp] at hmul
  have hext : (fun k : ℕ => (k : ℝ) * geom q (k + 1)) =
      (fun k : ℕ => (1 - q) * ((k : ℝ) * q ^ k)) := by
    funext k; rw [geom_succ]; ring
  rw [hext]; exact hmul

/-- The geometric entropy series has closed-form sum. -/
lemma hasSum_entropy_geom {q : ℝ} (h0 : 0 < q) (h1 : q < 1) :
    HasSum (fun k : ℕ => negMulLog (geom q (k + 1)))
      (-Real.log (1 - q) - (q / (1 - q)) * Real.log q) := by
  have h0' : 0 ≤ q := le_of_lt h0
  have hpiece1 : HasSum (fun k : ℕ => -Real.log (1 - q) * geom q (k + 1))
      (-Real.log (1 - q) * 1) :=
    (hasSum_geom h0' h1).mul_left _
  have hpiece2 : HasSum (fun k : ℕ => -Real.log q * ((k : ℝ) * geom q (k + 1)))
      (-Real.log q * (q / (1 - q))) :=
    (hasSum_k_geom h0 h1).mul_left _
  have hsum : HasSum
      (fun k : ℕ => -Real.log (1 - q) * geom q (k + 1)
                    + -Real.log q * ((k : ℝ) * geom q (k + 1)))
      (-Real.log (1 - q) * 1 + -Real.log q * (q / (1 - q))) :=
    hpiece1.add hpiece2
  have hext : (fun k : ℕ => negMulLog (geom q (k + 1))) =
      (fun k : ℕ => -Real.log (1 - q) * geom q (k + 1)
                    + -Real.log q * ((k : ℝ) * geom q (k + 1))) := by
    funext k
    have hgk : geom q (k + 1) = (1 - q) * q ^ k := geom_succ q k
    have hlog_gk : Real.log (geom q (k + 1)) =
        Real.log (1 - q) + (k : ℝ) * Real.log q := by
      rw [hgk, Real.log_mul (by linarith) (pow_ne_zero _ (ne_of_gt h0)),
          Real.log_pow]
    show negMulLog (geom q (k + 1)) = _
    unfold negMulLog
    rw [hlog_gk]; ring
  rw [hext]
  convert hsum using 1
  ring

/-- The closed form: `entropy_geom q = -log(1-q) - (q/(1-q)) log q`. -/
theorem entropy_geom_closed_form {q : ℝ} (h0 : 0 < q) (h1 : q < 1) :
    entropy_geom q = -Real.log (1 - q) - (q / (1 - q)) * Real.log q :=
  (hasSum_entropy_geom h0 h1).tsum_eq

/-! ### Headline statement of L0 (Gibbs version) -/

/-- **L0 — max-entropy lemma (Gibbs form).** Among probability
    distributions `p : ℕ → ℝ` supported on `k ≥ 1` (`p 0 = 0`) with
    `∑ p_k = 1` and `∑ k·p_k = μ` (`μ > 1`), the entropy `-∑ p_k log p_k`
    is at most that of the geometric distribution with parameter
    `q = 1 - 1/μ`.

    The summability hypothesis `Summable (negMulLog ∘ (p ∘ Nat.succ))` is
    mild: it holds whenever the entropy of `p` is finite (e.g. for any `p`
    with finite support, or with `p_k ≤ C·q^k` for some `q < 1`).

    Proof: the classical Gibbs / `KL ≥ 0` argument. The pointwise upper
    bound `negMulLog p_k ≤ g_k − p_k − p_k · log g_k` is summed against
    the moment constraints `∑ p_k = ∑ g_k = 1` and
    `∑ k p_k = ∑ k g_k = μ - 1`; the cross-term collapses, and what
    remains equals `entropy_geom q`. -/
theorem L0_geometric_maximises_entropy
    (μ : ℝ) (hμ : 1 < μ) (p : ℕ → ℝ)
    (h_supp : p 0 = 0)
    (h_nonneg : ∀ k : ℕ, 0 ≤ p k)
    (h_sum : HasSum p 1)
    (h_mean : HasSum (fun k : ℕ => (k : ℝ) * p k) μ)
    (h_entropy_summable : Summable (fun k : ℕ => negMulLog (p (k + 1)))) :
    (∑' k : ℕ, negMulLog (p k)) ≤ entropy_geom (1 - 1 / μ) := by
  set q : ℝ := 1 - 1 / μ with hq_def
  have hq_range := q_of_mean_in_unit_interval hμ
  have hq0 : 0 ≤ q := hq_range.1
  have hq1 : q < 1 := hq_range.2
  have hμ_pos : 0 < μ := by linarith
  have hq_pos : 0 < q := by
    rw [hq_def]
    have h1m : 1 / μ < 1 := by rw [div_lt_one hμ_pos]; exact hμ
    linarith
  -- The geometric `g_k := geom q (k+1)`.
  have hg_sum : HasSum (fun k : ℕ => geom q (k + 1)) 1 := hasSum_geom hq0 hq1
  -- q/(1-q) = μ - 1
  have h1q : 1 - q = 1 / μ := by rw [hq_def]; ring
  have hμ_ne : μ ≠ 0 := ne_of_gt hμ_pos
  have hq_over : q / (1 - q) = μ - 1 := by
    rw [h1q, hq_def]; field_simp
  -- Step 2: summability of `p ∘ Nat.succ` and `k ↦ k * p (k+1)`.
  have h_sum_succ : HasSum (fun k : ℕ => p (k + 1)) 1 := by
    -- hasSum_nat_add_iff' (1) : HasSum (fun n => p (n+1)) (g - ∑ range 1, p) ↔ HasSum p g.
    have h := (hasSum_nat_add_iff' (f := p) 1).mpr h_sum
    have hsum_range : (∑ i ∈ Finset.range 1, p i) = 0 := by
      simp [h_supp]
    rw [hsum_range, sub_zero] at h
    exact h
  have h_mean_succ : HasSum (fun k : ℕ => (k : ℝ) * p (k + 1)) (μ - 1) := by
    have h := (hasSum_nat_add_iff' (f := fun k => (k : ℝ) * p k) 1).mpr h_mean
    have hsum_range : (∑ i ∈ Finset.range 1, (i : ℝ) * p i) = 0 := by
      simp
    rw [hsum_range, sub_zero] at h
    -- h : HasSum (fun n => ((n+1):ℝ) * p (n+1)) μ
    have hext : (fun k : ℕ => ((k + 1 : ℕ) : ℝ) * p (k + 1)) =
        (fun k : ℕ => ((k : ℝ) + 1) * p (k + 1)) := by
      funext k; push_cast; ring
    rw [hext] at h
    have hsum_combined : HasSum
        (fun k : ℕ => ((k : ℝ) + 1) * p (k + 1) - p (k + 1)) (μ - 1) :=
      h.sub h_sum_succ
    have hext' : (fun k : ℕ => ((k : ℝ) + 1) * p (k + 1) - p (k + 1)) =
        (fun k : ℕ => (k : ℝ) * p (k + 1)) := by
      funext k; ring
    rw [hext'] at hsum_combined; exact hsum_combined
  -- Step 1: rewrite ∑' k, negMulLog (p k) = ∑' k, negMulLog (p (k+1)).
  have h_shift : (∑' k : ℕ, negMulLog (p k))
      = ∑' k : ℕ, negMulLog (p (k + 1)) := by
    -- Use Summable.tsum_eq_zero_add: ∑' f = f 0 + ∑' (f ∘ succ).
    have : Summable (fun k : ℕ => negMulLog (p k)) := by
      have := h_entropy_summable
      rwa [← summable_nat_add_iff 1]
    rw [this.tsum_eq_zero_add]
    simp [h_supp]
  rw [h_shift]
  -- Step 3: the pointwise upper bound.
  have h_bound : ∀ k, negMulLog (p (k + 1)) ≤
      geom q (k + 1) - p (k + 1)
        - p (k + 1) * Real.log (1 - q) - (k : ℝ) * p (k + 1) * Real.log q :=
    fun k => negMulLog_p_le_geom hq_pos hq1 h_nonneg k
  -- Step 4: the RHS series sums to -log(1-q) - (μ-1) log q = entropy_geom q.
  have h_rhs_sum : HasSum
      (fun k : ℕ => geom q (k + 1) - p (k + 1)
        - p (k + 1) * Real.log (1 - q) - (k : ℝ) * p (k + 1) * Real.log q)
      (1 - 1 - 1 * Real.log (1 - q) - (μ - 1) * Real.log q) := by
    have h_a : HasSum (fun k : ℕ => geom q (k + 1) - p (k + 1)) (1 - 1) :=
      hg_sum.sub h_sum_succ
    have h_b : HasSum (fun k : ℕ => p (k + 1) * Real.log (1 - q))
        (1 * Real.log (1 - q)) := h_sum_succ.mul_right _
    have h_c : HasSum (fun k : ℕ => (k : ℝ) * p (k + 1) * Real.log q)
        ((μ - 1) * Real.log q) := h_mean_succ.mul_right _
    have hab := h_a.sub h_b
    have habc := hab.sub h_c
    convert habc using 2 with k
  -- Step 5: combine.
  have h_rhs_summable : Summable
      (fun k : ℕ => geom q (k + 1) - p (k + 1)
        - p (k + 1) * Real.log (1 - q) - (k : ℝ) * p (k + 1) * Real.log q) :=
    h_rhs_sum.summable
  have h_le_sum :
      (∑' k : ℕ, negMulLog (p (k + 1))) ≤
        ∑' k : ℕ, (geom q (k + 1) - p (k + 1)
          - p (k + 1) * Real.log (1 - q) - (k : ℝ) * p (k + 1) * Real.log q) :=
    h_entropy_summable.tsum_le_tsum h_bound h_rhs_summable
  rw [h_rhs_sum.tsum_eq] at h_le_sum
  -- Step 6: identify RHS with entropy_geom q.
  have h_eg : entropy_geom q = -Real.log (1 - q) - (q / (1 - q)) * Real.log q :=
    entropy_geom_closed_form hq_pos hq1
  rw [h_eg, hq_over]
  linarith [h_le_sum]

end PT.Information
