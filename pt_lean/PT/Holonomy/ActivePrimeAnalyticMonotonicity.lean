/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Holonomy.ActivePrimeCriterion
import PT.Holonomy.ActivePrimeMonotonicity
import Mathlib.Tactic

/-!
# Active Prime Analytic Monotonicity — uniform inactivity for every `p ≥ 11`

**Statement (paper-level, Ch06 §"Critère de premier actif", D.1 général).**
For *every* natural `p ≥ 11`, the closed-form anomalous dimension

$$\gamma_p = \frac{4\,q^{p-1}(1 - \delta_p)}{\mu^* \,\delta_p\,(2 - \delta_p)},
   \qquad q = \tfrac{13}{15}, \qquad \delta_p = \tfrac{1 - q^p}{p}$$

lies strictly below the activity threshold `s = 1/2`. This closes the analytic
side of the active-prime criterion (the wave-1 module
`ActivePrimeMonotonicity` enumerated only `11 ≤ p ≤ 251` via `norm_num`).

## Strategy

The proof factors through an effective upper bound on `γ_p`:

1. **Coarse bounds on `δ_p`** (for `p ≥ 11`):
   * `δ_p ≤ 1/p`  (since `q^p ≥ 0`);
   * `δ_p ≥ (1 - q^{11})/p`  (since `q^p ≤ q^{11}` when `p ≥ 11`, `0 < q < 1`);
   * `1 - δ_p ≤ 1`  (since `δ_p ≥ 0`);
   * `2 - δ_p ≥ 21/11`  (since `δ_p ≤ 1/p ≤ 1/11`).

2. **Effective bound**:
   `γ_p ≤ (44 / (315 · (1 - q^{11}))) · p · q^{p-1} =: K · p · q^{p-1}`.

3. **Monotonicity of `p ↦ p · q^{p-1}` for `p ≥ 11`**:
   The successor ratio is `q · (p+1)/p`, which is `< 1` iff `q(p+1) < p`,
   i.e. `13(p+1) < 15 p`, i.e. `p ≥ 7`. So the sequence `p · q^{p-1}` is
   strictly decreasing for `p ≥ 7`.

4. **Numerical anchor at `p = 11`**:
   `K · 11 · q^{10} < 1/2` — an exact rational inequality decided by
   `norm_num`.

Steps (3)+(4) yield `K · p · q^{p-1} ≤ K · 11 · q^{10} < 1/2` for all `p ≥ 11`.

## Reference

Monograph Chapter 6, §"Critère de premier actif", gap D.1 (général).
`LEAN_MONOGRAPHIE_GAPS.md`.
-/

namespace PT.Holonomy

/-! ### Step A — Positivity and monotonicity of `q^p` -/

lemma qPT_pos : (0 : ℚ) < qPT := by
  unfold qPT; norm_num

lemma qPT_lt_one : qPT < 1 := by
  unfold qPT; norm_num

lemma qPT_nonneg : (0 : ℚ) ≤ qPT := le_of_lt qPT_pos

lemma qPT_le_one : qPT ≤ 1 := le_of_lt qPT_lt_one

lemma qPT_pow_pos (p : ℕ) : (0 : ℚ) < qPT ^ p :=
  pow_pos qPT_pos p

lemma qPT_pow_nonneg (p : ℕ) : (0 : ℚ) ≤ qPT ^ p :=
  le_of_lt (qPT_pow_pos p)

lemma qPT_pow_le_one (p : ℕ) : qPT ^ p ≤ 1 :=
  pow_le_one₀ qPT_nonneg qPT_le_one

/-- `q^p` is antitone in `p` (with `0 < q < 1`). -/
lemma qPT_pow_antitone {p q : ℕ} (h : p ≤ q) : qPT ^ q ≤ qPT ^ p :=
  pow_le_pow_of_le_one qPT_nonneg qPT_le_one h

/-! ### Step B — Bounds on `δ_p` -/

/-- `δ_p > 0` for `p ≥ 1`. -/
lemma deltaQ_pos {p : ℕ} (hp : 1 ≤ p) : 0 < deltaQ p := by
  unfold deltaQ
  have hp' : (0 : ℚ) < (p : ℚ) := by exact_mod_cast Nat.lt_of_lt_of_le Nat.zero_lt_one hp
  have hnum : (0 : ℚ) < 1 - qPT ^ p := by
    have := qPT_pow_le_one p
    have hlt : qPT ^ p < 1 := by
      have : qPT ^ p ≤ qPT ^ 1 := qPT_pow_antitone hp
      simp at this
      calc qPT ^ p ≤ qPT := this
        _ < 1 := qPT_lt_one
    linarith
  exact div_pos hnum hp'

/-- Upper bound: `δ_p ≤ 1/p` for `p ≥ 1`. -/
lemma deltaQ_upper_bound {p : ℕ} (hp : 1 ≤ p) : deltaQ p ≤ 1 / (p : ℚ) := by
  unfold deltaQ
  have hp' : (0 : ℚ) < (p : ℚ) := by exact_mod_cast Nat.lt_of_lt_of_le Nat.zero_lt_one hp
  rw [div_le_div_iff₀ hp' hp']
  have := qPT_pow_nonneg p
  nlinarith

/-- Lower bound: `δ_p ≥ (1 - q^{11})/p` for `p ≥ 11`. -/
lemma deltaQ_lower_bound {p : ℕ} (hp : 11 ≤ p) :
    (1 - qPT ^ 11) / (p : ℚ) ≤ deltaQ p := by
  unfold deltaQ
  have hp' : (0 : ℚ) < (p : ℚ) := by
    exact_mod_cast Nat.lt_of_lt_of_le (by norm_num : 0 < 11) hp
  rw [div_le_div_iff₀ hp' hp']
  have hmono : qPT ^ p ≤ qPT ^ 11 := qPT_pow_antitone hp
  nlinarith

/-- `1 - δ_p ≤ 1` for `p ≥ 1`. -/
lemma one_sub_deltaQ_le_one {p : ℕ} (hp : 1 ≤ p) : 1 - deltaQ p ≤ 1 := by
  have := deltaQ_pos hp
  linarith

/-- `2 - δ_p ≥ 21/11` for `p ≥ 11`. -/
lemma two_sub_deltaQ_lower {p : ℕ} (hp : 11 ≤ p) : (21 : ℚ) / 11 ≤ 2 - deltaQ p := by
  have hp1 : 1 ≤ p := le_trans (by norm_num) hp
  have hub := deltaQ_upper_bound hp1
  have hp' : (0 : ℚ) < (p : ℚ) := by
    exact_mod_cast Nat.lt_of_lt_of_le (by norm_num : 0 < 11) hp
  -- 1/p ≤ 1/11 when p ≥ 11
  have h11 : (1 : ℚ) / (p : ℚ) ≤ 1 / 11 := by
    rw [div_le_div_iff₀ hp' (by norm_num : (0 : ℚ) < 11)]
    have : (11 : ℚ) ≤ (p : ℚ) := by exact_mod_cast hp
    linarith
  have : deltaQ p ≤ 1 / 11 := le_trans hub h11
  linarith

/-- `2 - δ_p > 0` for `p ≥ 11`. -/
lemma two_sub_deltaQ_pos {p : ℕ} (hp : 11 ≤ p) : 0 < 2 - deltaQ p := by
  have := two_sub_deltaQ_lower hp
  linarith

/-! ### Step C — Effective bound `γ_p ≤ K · p · q^{p-1}` -/

/-- The effective constant `K = 44 / (315 · (1 - q^{11}))`. -/
def Kbound : ℚ := 44 / (315 * (1 - qPT ^ 11))

lemma Kbound_pos : 0 < Kbound := by
  unfold Kbound
  apply div_pos (by norm_num)
  have : qPT ^ 11 < 1 := by
    calc qPT ^ 11 ≤ qPT := by
          have := qPT_pow_antitone (show 1 ≤ 11 by norm_num)
          simpa using this
      _ < 1 := qPT_lt_one
  have h11pos : (0 : ℚ) < 1 - qPT ^ 11 := by linarith
  have : (0 : ℚ) < 315 := by norm_num
  positivity

/-- Effective upper bound on `γ_p` for `p ≥ 11`. -/
lemma gammaQ_le_Kbound_pq {p : ℕ} (hp : 11 ≤ p) :
    gammaQ p ≤ Kbound * (p : ℚ) * qPT ^ (p - 1) := by
  have hp1 : 1 ≤ p := le_trans (by norm_num) hp
  have hp' : (0 : ℚ) < (p : ℚ) := by
    exact_mod_cast Nat.lt_of_lt_of_le (by norm_num : 0 < 11) hp
  have hdp : 0 < deltaQ p := deltaQ_pos hp1
  have hdp_ub : deltaQ p ≤ 1 / (p : ℚ) := deltaQ_upper_bound hp1
  have hdp_lb : (1 - qPT ^ 11) / (p : ℚ) ≤ deltaQ p := deltaQ_lower_bound hp
  have h1md : 1 - deltaQ p ≤ 1 := one_sub_deltaQ_le_one hp1
  have h2md : (21 : ℚ) / 11 ≤ 2 - deltaQ p := two_sub_deltaQ_lower hp
  have h2md_pos : 0 < 2 - deltaQ p := two_sub_deltaQ_pos hp
  have h11pos : (0 : ℚ) < 1 - qPT ^ 11 := by
    have : qPT ^ 11 < 1 := by
      calc qPT ^ 11 ≤ qPT := by
            have := qPT_pow_antitone (show 1 ≤ 11 by norm_num)
            simpa using this
        _ < 1 := qPT_lt_one
    linarith
  have hqp : 0 < qPT ^ (p - 1) := qPT_pow_pos _
  have hqp_nn : 0 ≤ qPT ^ (p - 1) := le_of_lt hqp
  have h1md_nn : 0 ≤ 1 - deltaQ p := by
    have := deltaQ_upper_bound hp1
    have : deltaQ p ≤ 1 := by
      calc deltaQ p ≤ 1 / (p : ℚ) := this
        _ ≤ 1 := by
              rw [div_le_one hp']
              exact_mod_cast hp1
    linarith
  -- We want: γ_p = N / D ≤ K * p * q^(p-1), where N, D > 0
  -- N = 4 * q^(p-1) * (1 - δ_p), D = 15 * δ_p * (2 - δ_p)
  -- Rewriting: γ_p * D = N, and N ≤ 4 * q^(p-1) (since 1 - δ_p ≤ 1)
  -- So γ_p ≤ 4 * q^(p-1) / D ≤ 4 * q^(p-1) / (15 * δ_p * (2 - δ_p))
  -- ≤ 4 * q^(p-1) / (15 * (1-q^11)/p * 21/11)
  -- = 4 * 11 * p * q^(p-1) / (15 * 21 * (1 - q^11))
  -- = (44 / (315 * (1 - q^11))) * p * q^(p-1) = K * p * q^(p-1)
  unfold gammaQ muStar Kbound
  rw [div_le_iff₀ (by positivity)]
  -- Goal: 4 * q^(p-1) * (1 - δ_p) ≤ (44 / (315*(1-q^11))) * p * q^(p-1) * (15 * δ_p * (2 - δ_p))
  -- Use 1 - δ_p ≤ 1 on LHS, then bound from below the RHS factors
  have hDpos : 0 < 15 * deltaQ p * (2 - deltaQ p) := by positivity
  have h15dpos : (0 : ℚ) < 15 := by norm_num
  -- Step 1: LHS ≤ 4 * q^(p-1)
  have step1 : 4 * qPT ^ (p - 1) * (1 - deltaQ p) ≤ 4 * qPT ^ (p - 1) := by
    have h4qp_nn : 0 ≤ 4 * qPT ^ (p - 1) := by positivity
    nlinarith [h1md, h1md_nn]
  -- Step 2: 4 * q^(p-1) ≤ K * p * q^(p-1) * (15 * δ_p * (2 - δ_p))
  -- i.e. 4 ≤ K * p * 15 * δ_p * (2 - δ_p)  (after dividing by q^(p-1))
  -- K * p * δ_p ≥ K * p * (1-q^11)/p = K * (1-q^11) = 44/315
  -- K * p * δ_p * (2 - δ_p) ≥ (44/315) * (21/11) = (44*21)/(315*11) = 924/3465 = 4/15
  -- So K * p * δ_p * (2 - δ_p) * 15 ≥ 4. ✓
  have step2 : 4 * qPT ^ (p - 1) ≤
      44 / (315 * (1 - qPT ^ 11)) * (p : ℚ) * qPT ^ (p - 1) * (15 * deltaQ p * (2 - deltaQ p)) := by
    -- Factor q^(p-1) on both sides (both nonneg)
    have hgoal :
        4 ≤ 44 / (315 * (1 - qPT ^ 11)) * (p : ℚ) * (15 * deltaQ p * (2 - deltaQ p)) := by
      -- Use δ_p ≥ (1-q^11)/p:
      -- p * δ_p ≥ 1 - q^11
      have hpd : (1 - qPT ^ 11) ≤ (p : ℚ) * deltaQ p := by
        have := hdp_lb
        rw [div_le_iff₀ hp'] at this
        linarith
      -- Then (p * δ_p) * (2 - δ_p) ≥ (1 - q^11) * 21/11
      have hpd2 : (1 - qPT ^ 11) * (21 / 11) ≤ (p : ℚ) * deltaQ p * (2 - deltaQ p) := by
        have h2md_nn : (0 : ℚ) ≤ 2 - deltaQ p := le_of_lt h2md_pos
        have hpdnn : 0 ≤ (p : ℚ) * deltaQ p := by positivity
        have h1q11nn : 0 ≤ 1 - qPT ^ 11 := le_of_lt h11pos
        nlinarith [hpd, h2md, h2md_nn]
      -- Multiply by 15: 15 * (p * δ_p) * (2 - δ_p) ≥ 15 * (1 - q^11) * 21/11
      -- = 315 * (1-q^11) / 11 * (21/315) ... let me just nlinarith
      -- Then K * 15 * p * δ_p * (2-δ_p) ≥ (44/(315*(1-q^11))) * 15 * (1-q^11) * 21/11
      --                                  = 44*15*21/(315*11) = 4
      have hKpos : 0 < 44 / (315 * (1 - qPT ^ 11)) := by
        apply div_pos (by norm_num)
        have : (0 : ℚ) < 315 := by norm_num
        positivity
      -- Multiply both sides by K (positive):
      have hmul : 44 / (315 * (1 - qPT ^ 11)) * ((1 - qPT ^ 11) * (21 / 11)) ≤
                  44 / (315 * (1 - qPT ^ 11)) * ((p : ℚ) * deltaQ p * (2 - deltaQ p)) := by
        exact mul_le_mul_of_nonneg_left hpd2 (le_of_lt hKpos)
      -- LHS simplifies to 4/15? Let's check: 44/(315*(1-q^11)) * (1-q^11) * 21/11
      -- = 44 * 21 / (315 * 11) = 924 / 3465 = 4/15
      have hsimp : 44 / (315 * (1 - qPT ^ 11)) * ((1 - qPT ^ 11) * (21 / 11)) = 4 / 15 := by
        have hne : (1 - qPT ^ 11) ≠ 0 := ne_of_gt h11pos
        field_simp
        ring
      rw [hsimp] at hmul
      -- Now goal: 4 ≤ K * p * (15 * δ_p * (2 - δ_p))
      -- We have: 4/15 ≤ K * (p * δ_p * (2 - δ_p))
      -- So 4 ≤ K * 15 * p * δ_p * (2 - δ_p) = K * p * (15 * δ_p * (2 - δ_p))
      linarith [mul_le_mul_of_nonneg_left hmul (by norm_num : (0:ℚ) ≤ 15)]
    -- Multiply hgoal by q^(p-1) ≥ 0
    have := mul_le_mul_of_nonneg_right hgoal hqp_nn
    -- Rearrange to match goal
    nlinarith [this, hqp_nn, hgoal]
  linarith

/-! ### Step D — Monotonicity of `p ↦ p · q^{p-1}` for `p ≥ 7` -/

/-- Inductive monotonicity step: for `p ≥ 7`, the next term is no larger. -/
lemma succ_pq_le {p : ℕ} (hp : 7 ≤ p) :
    ((p + 1 : ℕ) : ℚ) * qPT ^ ((p + 1) - 1) ≤ (p : ℚ) * qPT ^ (p - 1) := by
  -- (p+1) - 1 = p, so RHS factor is q^p = q^(p-1) * q
  have hp1 : 1 ≤ p := le_trans (by norm_num) hp
  have hpsub : p - 1 + 1 = p := Nat.sub_add_cancel hp1
  have hsucc : (p + 1) - 1 = p := by omega
  rw [hsucc]
  -- Goal: (p+1) * q^p ≤ p * q^(p-1)
  -- i.e. (p+1) * q^(p-1) * q ≤ p * q^(p-1)
  -- i.e. q * (p+1) ≤ p   (dividing by q^(p-1) > 0)
  -- i.e. (13/15) * (p+1) ≤ p
  -- i.e. 13 * (p+1) ≤ 15 p
  -- i.e. 13 ≤ 2 p
  -- i.e. p ≥ 7 (since p natural, 2p ≥ 14 ≥ 13)
  have hq_eq : qPT ^ p = qPT ^ (p - 1) * qPT := by
    conv_lhs => rw [← hpsub, pow_succ]
  rw [hq_eq]
  have hqp_nn : 0 ≤ qPT ^ (p - 1) := qPT_pow_nonneg _
  have hkey : qPT * ((p : ℚ) + 1) ≤ (p : ℚ) := by
    unfold qPT
    have hpQ : (7 : ℚ) ≤ (p : ℚ) := by exact_mod_cast hp
    -- (13/15)(p+1) ≤ p ↔ 13(p+1) ≤ 15 p ↔ 13 ≤ 2 p ↔ p ≥ 7 (rational)
    nlinarith
  -- ((p+1):ℕ → ℚ) = (p:ℚ) + 1
  push_cast
  nlinarith [hqp_nn, hkey, qPT_nonneg]

/-- For `p ≥ 11`, `p · q^{p-1} ≤ 11 · q^{10}`. -/
lemma pq_le_eleven_q10 {p : ℕ} (hp : 11 ≤ p) :
    (p : ℚ) * qPT ^ (p - 1) ≤ (11 : ℚ) * qPT ^ 10 := by
  -- Induct on p - 11
  induction p, hp using Nat.le_induction with
  | base =>
      norm_num
  | succ n hn ih =>
      -- n ≥ 11 ≥ 7, so succ_pq_le applies
      have h7 : 7 ≤ n := le_trans (by norm_num) hn
      have step := succ_pq_le h7
      exact le_trans step ih

/-! ### Step E — Numerical anchor: `K · 11 · q^{10} < 1/2` -/

lemma Kbound_anchor : Kbound * 11 * qPT ^ 10 < sPT := by
  unfold Kbound qPT sPT
  norm_num

/-! ### Headline theorem -/

/-- **General analytic monotonicity (gap D.1 général).** For every natural
    `p ≥ 11`, the closed-form anomalous dimension `γ_p` lies strictly below
    the activity threshold `s = 1/2`. Combined with `gamma_3_active`,
    `gamma_5_active`, `gamma_7_active`, this shows that the active set is
    exactly `{3, 5, 7}` — uniformly in `p`, no enumeration needed. -/
theorem gammaQ_lt_half_of_ge_eleven (p : ℕ) (hp : 11 ≤ p) : gammaQ p < sPT := by
  -- γ_p ≤ K · p · q^{p-1} ≤ K · 11 · q^{10} < 1/2.
  have h1 : gammaQ p ≤ Kbound * (p : ℚ) * qPT ^ (p - 1) := gammaQ_le_Kbound_pq hp
  have h2 : (p : ℚ) * qPT ^ (p - 1) ≤ (11 : ℚ) * qPT ^ 10 := pq_le_eleven_q10 hp
  have hKnn : 0 ≤ Kbound := le_of_lt Kbound_pos
  have h3 : Kbound * ((p : ℚ) * qPT ^ (p - 1)) ≤ Kbound * ((11 : ℚ) * qPT ^ 10) :=
    mul_le_mul_of_nonneg_left h2 hKnn
  have h3' : Kbound * (p : ℚ) * qPT ^ (p - 1) ≤ Kbound * 11 * qPT ^ 10 := by
    have ea : Kbound * (p : ℚ) * qPT ^ (p - 1) = Kbound * ((p : ℚ) * qPT ^ (p - 1)) := by ring
    have eb : Kbound * 11 * qPT ^ 10 = Kbound * ((11 : ℚ) * qPT ^ 10) := by ring
    rw [ea, eb]; exact h3
  have h4 : Kbound * 11 * qPT ^ 10 < sPT := Kbound_anchor
  linarith

/-- **Final form.** The active set at `μ* = 15` consists of exactly the odd
    primes `{3, 5, 7}` — for every other natural `p ≥ 11`, `γ_p < 1/2`. -/
theorem active_primes_complete (p : ℕ) (hp : 11 ≤ p) : ¬ IsActive p := by
  intro h
  have := gammaQ_lt_half_of_ge_eleven p hp
  exact absurd this (not_lt.mpr (le_of_lt h))

end PT.Holonomy
