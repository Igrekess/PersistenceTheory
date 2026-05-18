/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Information.L0MaxEntropy

/-!
# L0 — Uniqueness of the maximum-entropy distribution

This file complements `PT.Information.L0MaxEntropy` (which proves the Gibbs
upper bound `H(p) ≤ entropy_geom q` for `q = 1 - 1/μ`) by proving the
**uniqueness** statement: any distribution `p : ℕ → ℝ` supported on `k ≥ 1`
with `∑ p_k = 1`, `∑ k·p_k = μ`, and entropy **equal** to
`entropy_geom (1 - 1/μ)` must coincide pointwise with the geometric
distribution `geom (1 - 1/μ)`.

## Proof sketch

The pointwise Gibbs bound
`negMulLog (p_k) ≤ g_k − p_k − p_k log g_k`
(see `negMulLog_le_of_pos` in `L0MaxEntropy`) is **strict** unless
`p_k = g_k`. Saturating the Gibbs sum therefore forces termwise equality.
Concretely:

1. Define the (nonneg) gap `d k := (g_k − p_k − p_k log g_k) − negMulLog (p_k)`.
2. The Gibbs bound says `d k ≥ 0`.
3. The summed bound and the entropy hypothesis force `∑ d k = 0`.
4. By `Summable.tsum_pos` (in contrapositive), `d k = 0` for every `k`.
5. The strict Gibbs lemma `negMulLog_eq_iff` yields `p (k+1) = geom q (k+1)`.

The headline result `L0_geom_is_unique_maximiser` is the pointwise equality
`p k = geom q k` for all `k`, with `q = 1 - 1/μ`.

## References

* `PT.Information.L0MaxEntropy` (Gibbs upper bound).
* Monograph §1.3 (max-entropy lemma, uniqueness clause).
-/

namespace PT.Information

open Real

/-! ### Strict Gibbs lemma

The pointwise bound `negMulLog p ≤ g − p − p log g` (for `p ≥ 0`, `g > 0`)
becomes an equality iff `p = g`. -/

/-- **Strict Gibbs lemma.** For `p ≥ 0` and `g > 0`,
    `negMulLog p = g − p − p · log g` iff `p = g`. -/
lemma negMulLog_eq_iff_of_pos {p g : ℝ} (hp : 0 ≤ p) (hg : 0 < g) :
    negMulLog p = g - p - p * Real.log g ↔ p = g := by
  refine ⟨fun heq => ?_, fun hpg => ?_⟩
  · -- Forward: saturation forces p = g.
    rcases eq_or_lt_of_le hp with hp_eq | hp_pos
    · -- p = 0 case: bound becomes 0 = g, contradiction with g > 0.
      have hp0 : p = 0 := hp_eq.symm
      rw [hp0] at heq
      simp at heq
      -- heq : 0 = g - 0 - 0 * log g, which simplifies to 0 = g
      linarith
    · -- p > 0: strict Gibbs `log x < x - 1` for x ≠ 1 (x > 0), applied at x = g/p.
      by_contra hne
      have hpg_inv_pos : 0 < g / p := div_pos hg hp_pos
      have hpg_inv_ne_one : g / p ≠ 1 := by
        intro h
        apply hne
        field_simp at h
        linarith
      have hstrict : Real.log (g / p) < g / p - 1 :=
        Real.log_lt_sub_one_of_pos hpg_inv_pos hpg_inv_ne_one
      have hlog_div : Real.log (g / p) = Real.log g - Real.log p :=
        Real.log_div (ne_of_gt hg) (ne_of_gt hp_pos)
      rw [hlog_div] at hstrict
      -- hstrict : log g - log p < g / p - 1
      -- Multiply both sides by p > 0.
      have hmul : p * (Real.log g - Real.log p) < p * (g / p - 1) :=
        mul_lt_mul_of_pos_left hstrict hp_pos
      have hrhs : p * (g / p - 1) = g - p := by field_simp
      -- From `heq` (Gibbs equality): negMulLog p = g - p - p log g,
      -- i.e. -p · log p = g - p - p · log g, hence
      --      p · (log g - log p) = -(g - p) - p · log g + p · log g + ...
      -- Cleaner: heq gives  p * log p - p * log g = p - g  (a single linear identity).
      have hkey : p * Real.log p - p * Real.log g = p - g := by
        have hnml : Real.negMulLog p = -(p * Real.log p) := by
          unfold Real.negMulLog; ring
        rw [hnml] at heq
        linarith
      -- So p (log g - log p) = -(p log p - p log g) = -(p - g) = g - p.
      have hlhs : p * (Real.log g - Real.log p) = g - p := by
        have : p * (Real.log g - Real.log p)
            = -(p * Real.log p - p * Real.log g) := by ring
        rw [this, hkey]; ring
      rw [hlhs, hrhs] at hmul
      -- Now hmul : g - p < g - p, contradiction.
      exact lt_irrefl _ hmul
  · -- Backward: p = g gives negMulLog g = g - g - g * log g = -g * log g.
    subst hpg
    unfold Real.negMulLog
    ring

/-! ### Headline uniqueness theorem -/

/-- **L0 — Uniqueness of the max-entropy distribution.**
    Among probability distributions `p : ℕ → ℝ` supported on `k ≥ 1`
    (`p 0 = 0`) with `∑ p_k = 1` and `∑ k·p_k = μ` (`μ > 1`), saturation
    of the Gibbs bound — i.e. `H(p) = entropy_geom (1 − 1/μ)` — forces
    `p k = geom (1 − 1/μ) k` for every `k`.

    In particular `p_0 = 0` and, for `k ≥ 1`,
    `p_k = (1 − q) · q^(k−1)` with `q = 1 − 1/μ`. -/
theorem L0_geom_is_unique_maximiser
    (μ : ℝ) (hμ : 1 < μ) (p : ℕ → ℝ)
    (h_supp : p 0 = 0)
    (h_nonneg : ∀ k : ℕ, 0 ≤ p k)
    (h_sum : HasSum p 1)
    (h_mean : HasSum (fun k : ℕ => (k : ℝ) * p k) μ)
    (h_entropy_summable : Summable (fun k : ℕ => Real.negMulLog (p (k + 1))))
    (h_eq : (∑' k : ℕ, Real.negMulLog (p k)) = entropy_geom (1 - 1 / μ)) :
    ∀ k : ℕ, p k = geom (1 - 1 / μ) k := by
  set q : ℝ := 1 - 1 / μ with hq_def
  have hq_range := q_of_mean_in_unit_interval hμ
  have hq0 : 0 ≤ q := hq_range.1
  have hq1 : q < 1 := hq_range.2
  have hμ_pos : 0 < μ := by linarith
  have hq_pos : 0 < q := by
    rw [hq_def]
    have h1m : 1 / μ < 1 := by rw [div_lt_one hμ_pos]; exact hμ
    linarith
  -- Geometric distribution g_k = geom q (k+1) = (1-q) q^k.
  have hg_sum : HasSum (fun k : ℕ => geom q (k + 1)) 1 := hasSum_geom hq0 hq1
  have h1q : 1 - q = 1 / μ := by rw [hq_def]; ring
  have hμ_ne : μ ≠ 0 := ne_of_gt hμ_pos
  have hq_over : q / (1 - q) = μ - 1 := by
    rw [h1q, hq_def]; field_simp
  -- Shift the sums to k ≥ 1.
  have h_sum_succ : HasSum (fun k : ℕ => p (k + 1)) 1 := by
    have h := (hasSum_nat_add_iff' (f := p) 1).mpr h_sum
    have hsum_range : (∑ i ∈ Finset.range 1, p i) = 0 := by simp [h_supp]
    rw [hsum_range, sub_zero] at h
    exact h
  have h_mean_succ : HasSum (fun k : ℕ => (k : ℝ) * p (k + 1)) (μ - 1) := by
    have h := (hasSum_nat_add_iff' (f := fun k => (k : ℝ) * p k) 1).mpr h_mean
    have hsum_range : (∑ i ∈ Finset.range 1, (i : ℝ) * p i) = 0 := by simp
    rw [hsum_range, sub_zero] at h
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
  -- Shift the entropy sum.
  have h_shift : (∑' k : ℕ, Real.negMulLog (p k))
      = ∑' k : ℕ, Real.negMulLog (p (k + 1)) := by
    have h_summable_full : Summable (fun k : ℕ => Real.negMulLog (p k)) := by
      have := h_entropy_summable
      rwa [← summable_nat_add_iff 1]
    rw [h_summable_full.tsum_eq_zero_add]
    simp [h_supp]
  -- Pointwise Gibbs bound on geom.
  have h_bound : ∀ k, Real.negMulLog (p (k + 1)) ≤
      geom q (k + 1) - p (k + 1)
        - p (k + 1) * Real.log (1 - q) - (k : ℝ) * p (k + 1) * Real.log q :=
    fun k => negMulLog_p_le_geom hq_pos hq1 h_nonneg k
  -- The RHS sums to entropy_geom q.
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
  -- Identify the RHS scalar with entropy_geom q.
  have h_eg : entropy_geom q = -Real.log (1 - q) - (q / (1 - q)) * Real.log q :=
    entropy_geom_closed_form hq_pos hq1
  have h_rhs_value :
      1 - 1 - 1 * Real.log (1 - q) - (μ - 1) * Real.log q = entropy_geom q := by
    rw [h_eg, hq_over]; ring
  -- Define the (nonneg) gap.
  let d : ℕ → ℝ := fun k =>
    (geom q (k + 1) - p (k + 1)
      - p (k + 1) * Real.log (1 - q) - (k : ℝ) * p (k + 1) * Real.log q)
    - Real.negMulLog (p (k + 1))
  have hd_nonneg : ∀ k, 0 ≤ d k := fun k => sub_nonneg.mpr (h_bound k)
  have hd_summable : Summable d := by
    exact h_rhs_sum.summable.sub h_entropy_summable
  -- Sum of d k is 0.
  have hd_sum_zero : ∑' k, d k = 0 := by
    have hd_hasSum : HasSum d (entropy_geom q - ∑' k, Real.negMulLog (p (k + 1))) := by
      have hh : HasSum (fun k : ℕ => Real.negMulLog (p (k + 1)))
          (∑' k, Real.negMulLog (p (k + 1))) := h_entropy_summable.hasSum
      have h_rhs_val : HasSum
        (fun k : ℕ => geom q (k + 1) - p (k + 1)
          - p (k + 1) * Real.log (1 - q) - (k : ℝ) * p (k + 1) * Real.log q)
        (entropy_geom q) := by
        rw [← h_rhs_value]; exact h_rhs_sum
      exact h_rhs_val.sub hh
    have h_eq2 : ∑' k, Real.negMulLog (p (k + 1)) = entropy_geom q := by
      rw [← h_shift]; exact h_eq
    rw [h_eq2, sub_self] at hd_hasSum
    exact hd_hasSum.tsum_eq
  -- Each d k = 0 (contrapositive of Summable.tsum_pos).
  have hd_zero : ∀ k, d k = 0 := by
    intro k
    by_contra hne
    have hdk_pos : 0 < d k :=
      lt_of_le_of_ne (hd_nonneg k) (Ne.symm hne)
    have : 0 < ∑' j, d j := hd_summable.tsum_pos hd_nonneg k hdk_pos
    linarith
  -- Each d k = 0 means equality in the Gibbs bound at k+1.
  have h_eq_pw : ∀ k, Real.negMulLog (p (k + 1)) =
      geom q (k + 1) - p (k + 1)
        - p (k + 1) * Real.log (1 - q) - (k : ℝ) * p (k + 1) * Real.log q := by
    intro k
    have := hd_zero k
    show Real.negMulLog (p (k + 1)) = _
    linarith [this, hd_nonneg k]
  -- Convert termwise Gibbs-equality to p (k+1) = geom q (k+1).
  -- The Gibbs bound came from `negMulLog p ≤ g - p - p log g` with
  -- `g = geom q (k+1) = (1-q) q^k` and `log g = log(1-q) + k log q`.
  have h_pk_eq : ∀ k, p (k + 1) = geom q (k + 1) := by
    intro k
    have hgk : geom q (k + 1) = (1 - q) * q ^ k := geom_succ q k
    have hgk_pos : 0 < geom q (k + 1) := by
      rw [hgk]; exact mul_pos (by linarith) (pow_pos hq_pos k)
    have hlog_gk : Real.log (geom q (k + 1)) =
        Real.log (1 - q) + (k : ℝ) * Real.log q := by
      rw [hgk, Real.log_mul (by linarith) (pow_ne_zero _ (ne_of_gt hq_pos)),
          Real.log_pow]
    -- Restate the termwise equality in the form needed by the strict Gibbs.
    have h_rewrite : Real.negMulLog (p (k + 1)) =
        geom q (k + 1) - p (k + 1)
          - p (k + 1) * Real.log (geom q (k + 1)) := by
      rw [hlog_gk]
      have := h_eq_pw k
      linarith
    exact (negMulLog_eq_iff_of_pos (h_nonneg (k + 1)) hgk_pos).mp h_rewrite
  -- Final assembly.
  intro k
  cases k with
  | zero => rw [h_supp, geom_zero]
  | succ n => exact h_pk_eq n

end PT.Information
