/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import Mathlib.NumberTheory.PrimeCounting
import Mathlib.NumberTheory.SumPrimeReciprocals
import Mathlib.NumberTheory.Chebyshev
import Mathlib.NumberTheory.AbelSummation
import Mathlib.NumberTheory.ArithmeticFunction.VonMangoldt
import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.Analysis.SpecialFunctions.Stirling
import Mathlib.Analysis.Asymptotics.Defs
import Mathlib.Tactic

/-!
# T5 — Mertens' Third Theorem (Meissel–Mertens Constant)

This module scaffolds the formalisation of **Mertens' Third Theorem (M3)**, which
gives the second-order asymptotic for the partial sums of prime reciprocals:

$$\sum_{p \le x,\ p\ \text{prime}} \frac{1}{p}
  \;=\; \log\log x \,+\, M_3 \,+\, O\!\left(\frac{1}{\log x}\right),$$

where the constant `M₃ ≈ 0.26149721…` is the **Meissel–Mertens constant**.

## Role of T5 in Persistence Theory

In the PT corpus, "T5" is the *self-consistency theorem* of the M3 article
(`PT_ARTICLES/PT_MATHEMATICS/M3/m3_convergence.tex`): it proves that the sieve
admits a unique stable fixed point at `μ* = 15`. The analytic backbone of T5
is *the classical Mertens theorem*, used as an **external, non-PT-specific
ingredient** for one purpose only — to ensure **compactness** of the orbit of
transition matrices `T(k)` in the compact set of `3 × 3` stochastic matrices,
which forces continuous functionals of `T` to converge.

The PT article is explicit (`m3_convergence_fr.tex`, Remark `mertens_role`):

> Le théorème de Mertens (classique, indépendant de PT) est utilisé ici
> *uniquement* pour la compacité […]. Il fournit l'asymptotique de fond
> `α_k = 1/2 - O(1/ln p_k)`. C'est une conséquence du théorème des nombres
> premiers et tient pour *tout* crible multiplicatif.

So the Lean-level dependency from T5 (M3) onto Mertens M3 is the **bounded
oscillation** statement: `mertensSum x - log log x` is *bounded* for `x ≥ 2`.
The full asymptotic with the explicit constant `M₃` is the headline; the
weaker bounded-oscillation corollary is what PT actually needs.

## Sorry budget (scaffolding pass)

This file is **sorry-free** for all classical Mertens results (M1, M3-existence,
M3-quantitative, compactness alias).

| Name                              | Strategy                       | Status     |
|-----------------------------------|--------------------------------|------------|
| `log_factorial_eq_sum_vonMangoldt_mul_floor` | Induction via `Nat.succ_div` + `vonMangoldt_sum` | **closed** (foundation for M1) |
| `stirling_log_factorial_effective` | `Stirling.log_stirlingSeq_formula` + antitone + bounded constant | **closed** (Sub-lemma 1 for M1) |
| `mertensLog_eq_mertensLog_floor`  | Direct from `mertensLog` def    | **closed** (Sub-lemma 4 prep) |
| `abs_log_sub_log_floor_le`        | `Nat.floor_le` + `Nat.lt_floor_add_one` + `log_le_log`     | **closed** (Sub-lemma 4 prep) |
| `mertens_M1_of_nat`               | Real-to-integer reduction via floor | **closed** (Sub-lemma 4) |
| `summable_vonMangoldt_prime_power_tail` | Re-derivation of Mathlib's `summable_F''` (private) via `Real.log_le_rpow_div` + geometric + p-series + `prodNatEquiv` | **closed** (Sub-lemma 2, 2026-05-17) |
| `mertensLog_floor_eq_log_floor_add_bounded` | Legendre + Stirling + ψ-bound + prime-power tail, decomposition `⌊N/d⌋ = N/d − {N/d}` | **closed** (Sub-lemma 3, 2026-05-17) |
| `mertens_M1`                      | `mertens_M1_of_nat ∘ mertensLog_floor_eq_log_floor_add_bounded` | **closed** (2026-05-17) |
| `mertensSum_eq_mertensLog_div_log_add_integral` | Abel pivot for M3, `f(t) = (log t)⁻¹` | **closed** (Sub-lemma 5, 2026-05-17) |
| `integral_one_div_t_log_t`        | Antiderivative `∫ 1/(t log t) = log log t` via FTC | **closed** (Sub-lemma 6, 2026-05-17) |
| `integral_one_div_t_log_sq`       | Antiderivative `∫ 1/(t log²t) = -1/log t` via FTC | **closed** (Sub-lemma 7, 2026-05-17) |
| `mertensLog_monotone`             | `Finset.sum_le_sum_of_subset_of_nonneg` on `primesBelow` | **closed** (helper for M3, 2026-05-17) |
| `mertensLog_div_intervalIntegrable` | `Monotone.intervalIntegrable` + `mul_continuousOn` | **closed** (helper for M3, 2026-05-17) |
| `mertensSum_sub_log_log_oscillation` | Abel pivot + M1 + closed-form integrals, `|F(x) - F(y)| ≤ 3C/log y` | **closed** (Cauchy step for M3, 2026-05-17) |
| `mertens_M3_exists`               | Cauchy criterion via oscillation bound + `cauchySeq_tendsto_of_complete` | **closed** (2026-05-17) |
| `mertens_M3`                      | Quantitative refinement: oscillation between `x` and `z → ∞` | **closed** (2026-05-17) |
| `mertensSum_sub_log_log_bounded`  | Boundedness from `mertens_M3_exists` + monotonicity on `[2, X]` | **closed** |
| `pt_T5_mertens_compactness`       | Direct alias of `_bounded`     | **closed** (alias) |

Total: T5Mertens **entirely sorry-free**.

## Main definitions

* `PT.NumberTheory.mertensSum x` — partial sum `∑ 1/p` over primes `p ≤ x`.
* `PT.NumberTheory.mertensLog x` — auxiliary partial sum `∑ (log p)/p`.
* `PT.NumberTheory.mertensM3` — the Meissel–Mertens constant, defined as the
  classical limit (existence stubbed by `mertens_M3_exists`).

## Main theorems

* `mertens_M1` — Mertens' First Theorem: `∑_{p ≤ x} (log p)/p = log x + O(1)`.
* `mertens_M3_exists` *[closed]* — the limit defining `M₃` exists.
* `mertens_M3` *[closed]* — headline asymptotic with rate `O(1/log x)`.
* `mertensSum_sub_log_log_bounded` *[closed]* — the PT-needed compactness
  corollary (proved from `mertens_M3_exists` + finite-range monotonicity).
* `pt_T5_mertens_compactness` *[closed]* — alias used by the PT corpus,
  refers to `mertensSum_sub_log_log_bounded`.

## References

* F. Mertens, *Ein Beitrag zur analytischen Zahlentheorie*, J. Reine Angew. Math.
  78 (1874), 46–62.
* H. Davenport, *Multiplicative Number Theory*, 3rd ed., GTM 74 (2000), §6.
* G. Tenenbaum, *Introduction to Analytic and Probabilistic Number Theory*,
  3rd ed., Graduate Studies in Mathematics 163, §I.1.4–I.1.6.
* Y. Senez, *La théorie de la persistance*, monograph chapter `ch07_convergence`;
  Article M3 (`PT_ARTICLES/PT_MATHEMATICS/M3/m3_convergence.tex`), Remark
  `rem:mertens_role`.
-/

namespace PT.NumberTheory

open Real Nat Finset
open scoped BigOperators

/-! ### Partial sums -/

/-- **`mertensSum x`**: the partial sum `∑_{p ≤ x, prime} 1/p`, indexed by the
    natural-number cutoff `⌊x⌋₊`. We use `Nat.primesBelow` (primes strictly less
    than `n`) on `n = ⌊x⌋₊ + 1`, which equals the set of primes `≤ ⌊x⌋₊`.

    For `x < 0` we get `⌊x⌋₊ = 0`, and the sum is empty. For `x ∈ [0, 2)` the
    sum is still empty (no primes ≤ 1). -/
noncomputable def mertensSum (x : ℝ) : ℝ :=
  ∑ p ∈ Nat.primesBelow (⌊x⌋₊ + 1), (1 : ℝ) / p

/-- **`mertensLog x`**: the auxiliary partial sum `∑_{p ≤ x, prime} (log p)/p`,
    which appears in Mertens' First Theorem (`mertens_M1` below). -/
noncomputable def mertensLog (x : ℝ) : ℝ :=
  ∑ p ∈ Nat.primesBelow (⌊x⌋₊ + 1), Real.log p / p

/-! ### Trivial boundary cases (no `sorry`) -/

/-- For `x < 2`, the natural-number floor `⌊x⌋₊` is at most `1`. -/
private lemma floor_le_one_of_lt_two {x : ℝ} (hx : x < 2) : ⌊x⌋₊ ≤ 1 := by
  by_cases hx0 : 0 ≤ x
  · -- 0 ≤ x < 2 ⇒ ⌊x⌋₊ < 2 (via Nat.floor_lt) ⇒ ⌊x⌋₊ ≤ 1
    have : ⌊x⌋₊ < 2 := by
      rw [Nat.floor_lt hx0]
      exact_mod_cast hx
    omega
  · -- x < 0 ⇒ ⌊x⌋₊ = 0
    have hx0' : x < 0 := lt_of_not_ge hx0
    have : ⌊x⌋₊ = 0 := Nat.floor_eq_zero.mpr (by linarith)
    omega

/-- `primesBelow n = ∅` whenever `n ≤ 2`. -/
private lemma primesBelow_eq_empty_of_le_two {n : ℕ} (hn : n ≤ 2) :
    Nat.primesBelow n = ∅ := by
  rw [Nat.primesBelow_eq_filter_range]
  apply Finset.filter_eq_empty_iff.mpr
  intro p hp
  have hpr : p < n := Finset.mem_range.mp hp
  have hp2 : p < 2 := lt_of_lt_of_le hpr hn
  intro hprime
  exact absurd (hprime.two_le) (by omega)

/-- The Mertens sum on the empty range `x < 2` is `0`: there are no primes
    `≤ x` when `x < 2`. -/
theorem mertensSum_lt_two {x : ℝ} (hx : x < 2) : mertensSum x = 0 := by
  unfold mertensSum
  have hx_floor : ⌊x⌋₊ ≤ 1 := floor_le_one_of_lt_two hx
  rw [primesBelow_eq_empty_of_le_two (by omega : ⌊x⌋₊ + 1 ≤ 2), Finset.sum_empty]

/-- Likewise, the log-weighted sum vanishes for `x < 2`. -/
theorem mertensLog_lt_two {x : ℝ} (hx : x < 2) : mertensLog x = 0 := by
  unfold mertensLog
  have hx_floor : ⌊x⌋₊ ≤ 1 := floor_le_one_of_lt_two hx
  rw [primesBelow_eq_empty_of_le_two (by omega : ⌊x⌋₊ + 1 ≤ 2), Finset.sum_empty]

/-- `mertensSum` is non-negative for every `x` (each summand `1/p` is non-negative
    when `p` is a positive prime; trivially `0` below `x = 2`). -/
theorem mertensSum_nonneg (x : ℝ) : 0 ≤ mertensSum x := by
  unfold mertensSum
  apply Finset.sum_nonneg
  intro p hp
  have hp_pos : 0 < p := by
    rw [Nat.primesBelow_eq_filter_range, Finset.mem_filter] at hp
    exact hp.2.pos
  positivity

/-- `mertensLog` is non-negative for every `x` (each summand `log p / p ≥ 0`
    since primes are `≥ 2`, so `log p ≥ log 2 > 0`). -/
theorem mertensLog_nonneg (x : ℝ) : 0 ≤ mertensLog x := by
  unfold mertensLog
  apply Finset.sum_nonneg
  intro p hp
  have hp_prime : Nat.Prime p := by
    rw [Nat.primesBelow_eq_filter_range, Finset.mem_filter] at hp
    exact hp.2
  have hp_ge_two : (2 : ℝ) ≤ p := by exact_mod_cast hp_prime.two_le
  have hp_pos : (0 : ℝ) < p := by linarith
  have hlog : 0 ≤ Real.log p := Real.log_nonneg (by linarith)
  positivity

/-! ### Building blocks toward M1: Legendre's log-factorial identity

The classical Mertens 1874 proof of M1 (predating PNT) rests on combining
Stirling's approximation `log n! = n log n − n + O(log n)` with the
**Legendre summatory identity**

$$\log N! \;=\; \sum_{d=1}^{N} \Lambda(d) \,\lfloor N/d\rfloor,$$

which is the summatory form of `Λ * 1 = log` (Mathlib's
`ArithmeticFunction.vonMangoldt_sum`). Mathlib has the pointwise
identity at each `n` but not the summatory version over `Icc 1 N`. The
lemma below provides it by induction on `N`, using `Nat.succ_div` to
control the increment `⌊(N+1)/d⌋ − ⌊N/d⌋ = [d ∣ N+1]`. -/

open scoped ArithmeticFunction in
/-- **Legendre's log-factorial identity (summatory form).** For every
    natural `N`,

    `log N! = ∑_{d ∈ Icc 1 N} Λ(d) · ⌊N/d⌋`.

    Equivalently, `(Λ * 1)` summed over `[1, N]` equals `log N!`. This is
    the analytic identity behind Mertens' First Theorem and is the
    summatory dual of `ArithmeticFunction.vonMangoldt_sum`. -/
lemma log_factorial_eq_sum_vonMangoldt_mul_floor (N : ℕ) :
    Real.log ((N.factorial : ℕ) : ℝ)
      = ∑ d ∈ Finset.Icc 1 N,
          ArithmeticFunction.vonMangoldt d * ((N / d : ℕ) : ℝ) := by
  induction N with
  | zero => simp
  | succ n ih =>
    -- `(n+1)! = (n+1) · n!`, so `log (n+1)! = log (n+1) + log n!`.
    have h_split : Real.log (((n + 1).factorial : ℕ) : ℝ)
        = Real.log ((n + 1 : ℕ) : ℝ) + Real.log ((n.factorial : ℕ) : ℝ) := by
      have h1 : (((n + 1).factorial : ℕ) : ℝ) = ((n + 1 : ℕ) : ℝ) * ((n.factorial : ℕ) : ℝ) := by
        push_cast [Nat.factorial_succ]; ring
      rw [h1, Real.log_mul (by exact_mod_cast (Nat.succ_ne_zero n))
            (by exact_mod_cast n.factorial_ne_zero)]
    rw [h_split, ih]
    -- Goal:
    --   log (n+1) + ∑_{d ∈ Icc 1 n} Λ(d) · (n/d) = ∑_{d ∈ Icc 1 (n+1)} Λ(d) · ((n+1)/d).
    -- Strategy: split off `d = n+1` on the RHS, expand `(n+1)/d` via `Nat.succ_div`.
    -- First rewrite Icc 1 (n+1) as Icc 1 n ∪ {n+1}.
    have h_icc_split : Finset.Icc 1 (n + 1) = insert (n + 1) (Finset.Icc 1 n) := by
      ext k
      simp only [Finset.mem_insert, Finset.mem_Icc]
      omega
    rw [h_icc_split, Finset.sum_insert (by simp : n + 1 ∉ Finset.Icc 1 n)]
    -- RHS becomes Λ(n+1) · ((n+1)/(n+1)) + ∑_{d ∈ Icc 1 n} Λ(d) · ((n+1)/d).
    -- Simplify (n+1)/(n+1) = 1.
    rw [show ((n + 1) / (n + 1) : ℕ) = 1 by exact Nat.div_self (Nat.succ_pos n)]
    -- Use Nat.succ_div : (n+1)/d = n/d + (if d ∣ n+1 then 1 else 0).
    have h_inner : ∀ d ∈ Finset.Icc 1 n,
        ArithmeticFunction.vonMangoldt d * (((n + 1) / d : ℕ) : ℝ)
          = ArithmeticFunction.vonMangoldt d * ((n / d : ℕ) : ℝ)
            + ArithmeticFunction.vonMangoldt d
                * (if d ∣ n + 1 then (1 : ℝ) else 0) := by
      intro d _
      rw [Nat.succ_div]
      by_cases hdvd : d ∣ n + 1
      · simp [hdvd]; ring
      · simp [hdvd]
    rw [Finset.sum_congr rfl h_inner, Finset.sum_add_distrib]
    -- Now: Λ(n+1) · 1 + (∑_{d ∈ Icc 1 n} Λ(d) · ((n+1)/d)-from-split)
    --     = log(n+1) + ∑_{d ∈ Icc 1 n} Λ(d) · (n/d).
    -- The "if-1-else-0" piece selects divisors d of n+1 in Icc 1 n; adding Λ(n+1) gives
    -- the sum over (n+1).divisors, which equals log(n+1) by `vonMangoldt_sum`.
    have h_div_sum : ArithmeticFunction.vonMangoldt (n + 1)
        + ∑ d ∈ Finset.Icc 1 n,
            ArithmeticFunction.vonMangoldt d * (if d ∣ n + 1 then (1 : ℝ) else 0)
        = Real.log ((n + 1 : ℕ) : ℝ) := by
      -- Drop the indicator: sum over divisors of n+1 lying in Icc 1 n.
      have h_indicator : ∀ d ∈ Finset.Icc 1 n,
          ArithmeticFunction.vonMangoldt d * (if d ∣ n + 1 then (1 : ℝ) else 0)
            = if d ∣ n + 1 then ArithmeticFunction.vonMangoldt d else 0 := by
        intro d _; split_ifs <;> ring
      rw [Finset.sum_congr rfl h_indicator]
      -- Convert filter form, then use Nat.divisors.
      rw [← Finset.sum_filter]
      -- {d ∈ Icc 1 n | d ∣ n+1} ∪ {n+1} = (n+1).divisors
      have h_div_eq : (n + 1).divisors
          = insert (n + 1) ((Finset.Icc 1 n).filter (fun d => d ∣ n + 1)) := by
        ext d
        rw [Nat.mem_divisors]
        simp only [Finset.mem_insert, Finset.mem_filter, Finset.mem_Icc]
        constructor
        · rintro ⟨hdvd, _⟩
          have hd_pos : 1 ≤ d := by
            rcases Nat.eq_zero_or_pos d with rfl | h
            · exact absurd (Nat.eq_zero_of_zero_dvd hdvd) (Nat.succ_ne_zero n)
            · exact h
          have hd_le : d ≤ n + 1 := Nat.le_of_dvd (Nat.succ_pos n) hdvd
          by_cases h : d = n + 1
          · left; exact h
          · right; exact ⟨⟨hd_pos, by omega⟩, hdvd⟩
        · rintro (rfl | ⟨⟨_, _⟩, hdvd⟩)
          · exact ⟨dvd_refl _, Nat.succ_ne_zero _⟩
          · exact ⟨hdvd, Nat.succ_ne_zero _⟩
      have h_not_mem : n + 1 ∉ (Finset.Icc 1 n).filter (fun d => d ∣ n + 1) := by
        intro hmem
        rw [Finset.mem_filter, Finset.mem_Icc] at hmem
        omega
      rw [show ArithmeticFunction.vonMangoldt (n + 1)
            + ∑ d ∈ (Finset.Icc 1 n).filter (fun d => d ∣ n + 1),
                ArithmeticFunction.vonMangoldt d
          = ∑ d ∈ insert (n + 1) ((Finset.Icc 1 n).filter (fun d => d ∣ n + 1)),
              ArithmeticFunction.vonMangoldt d
          from (Finset.sum_insert h_not_mem).symm]
      rw [← h_div_eq]
      exact ArithmeticFunction.vonMangoldt_sum
    -- Finish by linear combination.
    have h_one_mul : ArithmeticFunction.vonMangoldt (n + 1) * (1 : ℝ)
        = ArithmeticFunction.vonMangoldt (n + 1) := mul_one _
    linarith [h_div_sum]

/-! ### Abel summation pivot: `mertensLog` as `θ(x)/x` plus an integral

The classical Mertens 1874 proof of M1 starts by **partial summation** rewriting
`∑_{p ≤ x} (log p)/p` in terms of the Chebyshev function `θ(x) := ∑_{p ≤ x} log p`:

$$\sum_{p \le x} \frac{\log p}{p}
   \;=\; \frac{\theta(x)}{x} \,+\, \int_2^x \frac{\theta(t)}{t^2}\,dt.$$

This is the analogue, for `f(t) = 1/t`, of Mathlib's
`Chebyshev.primeCounting_eq_theta_div_log_add_integral` (which uses
`f(t) = 1/log t`). Once available, combining this identity with the elementary
Chebyshev bound `θ(t) = O(t)` yields M1: the boundary term is `O(1)` and the
integral is `O(log x) - O(1)`, but, more precisely for M1, the identity together
with `θ(t) = t + o(t)` gives `mertensLog x = log x + O(1)`.
-/

open Asymptotics Filter MeasureTheory in
/-- **Abel-summation pivot for Mertens M1.**

For every `x ≥ 2`,
$$\sum_{p \le x} \frac{\log p}{p}
   \;=\; \frac{\theta(x)}{x} \,+\, \int_2^x \frac{\theta(t)}{t^2}\,dt.$$

This is the Mertens analogue of Mathlib's
`Chebyshev.primeCounting_eq_theta_div_log_add_integral` (applied to
`f(t) = 1/t` instead of `f(t) = 1/\log t`). It is the standard analytic
backbone of **Mertens' First Theorem**: combined with Chebyshev's bound
`θ(t) ≤ (\log 4) · t`, the right-hand side is `log x + O(1)`. -/
theorem mertensLog_eq_theta_div_x_add_integral {x : ℝ} (hx : 2 ≤ x) :
    mertensLog x = Chebyshev.theta x / x
      + ∫ t in (2 : ℝ)..x, Chebyshev.theta t / t ^ 2 := by
  -- Rewrite `mertensLog x` in the form to which Abel summation applies.
  unfold mertensLog
  rw [Nat.primesBelow_eq_filter_range, Nat.range_succ_eq_Icc_zero, Finset.sum_filter]
  -- The Abel-summation "sequence" `a(n) = log n · [n prime]`.
  let a : ℕ → ℝ := Set.indicator (setOf Nat.Prime) (fun n ↦ Real.log n)
  -- Step 1: rewrite the conditional `log p / p` as `(1/n) * a n`.
  trans ∑ n ∈ Finset.Icc 0 ⌊x⌋₊, ((n : ℝ))⁻¹ * a n
  · refine Finset.sum_congr rfl fun n _hn ↦ ?_
    by_cases h : Nat.Prime n
    · simp [a, h, Set.indicator_of_mem, div_eq_mul_inv, mul_comm]
    · simp [a, h, Set.indicator_of_notMem]
  -- Step 2: apply Abel summation (`sum_mul_eq_sub_integral_mul₁`).
  rw [sum_mul_eq_sub_integral_mul₁ a (f := fun t ↦ t⁻¹)
        (by simp [a]) (by simp [a, Nat.not_prime_one]),
      ← intervalIntegral.integral_of_le hx]
  · -- Step 3: simplify the boundary term and the integrand.
    -- The derivative `(d/dt) t⁻¹ = -t⁻²`.
    have hderiv : ∀ u ∈ Set.uIcc (2 : ℝ) x, deriv (fun t : ℝ ↦ t⁻¹) u = -(u ^ 2)⁻¹ := by
      intro u _; simp [deriv_inv']
    have int_deriv (g : ℝ → ℝ) :
        ∫ u in (2 : ℝ)..x, deriv (fun t : ℝ ↦ t⁻¹) u * g u
          = ∫ u in (2 : ℝ)..x, g u * -(u ^ 2)⁻¹ :=
      intervalIntegral.integral_congr fun u hu ↦ by
        rw [hderiv u hu]; ring
    rw [int_deriv]
    -- The boundary term: `f x · ∑ a = (1/x) · θ(x) = θ(x)/x`.
    have hθ_sum : ∑ k ∈ Finset.Icc 0 ⌊x⌋₊, a k = Chebyshev.theta x := by
      rw [Chebyshev.theta_eq_sum_Icc, Finset.sum_filter]
      refine Finset.sum_congr rfl fun n _ ↦ ?_
      by_cases h : Nat.Prime n
      · simp [a, h, Set.indicator_of_mem]
      · simp [a, h, Set.indicator_of_notMem]
    -- The integrand: `∑ a over [0..⌊t⌋] = θ(t)` for every `t`.
    have hθ_partial : ∀ t : ℝ,
        ∑ k ∈ Finset.Icc 0 ⌊t⌋₊, a k = Chebyshev.theta t := by
      intro t
      rw [Chebyshev.theta_eq_sum_Icc, Finset.sum_filter]
      refine Finset.sum_congr rfl fun n _ ↦ ?_
      by_cases h : Nat.Prime n
      · simp [a, h, Set.indicator_of_mem]
      · simp [a, h, Set.indicator_of_notMem]
    rw [hθ_sum]
    -- Now: `(1/x) * θ(x) - ∫ θ(t) * -(t^2)⁻¹  = θ(x)/x + ∫ θ(t)/t^2`.
    have hint_eq :
        ∫ u in (2 : ℝ)..x, Chebyshev.theta u * -(u ^ 2)⁻¹
          = -∫ u in (2 : ℝ)..x, Chebyshev.theta u / u ^ 2 := by
      rw [← intervalIntegral.integral_neg]
      refine intervalIntegral.integral_congr fun u _ ↦ ?_
      simp [div_eq_mul_inv]
    -- Push the inner-sum-equals-θ rewriting into the integrand:
    have hint_eq' :
        ∫ u in (2 : ℝ)..x, (∑ k ∈ Finset.Icc 0 ⌊u⌋₊, a k) * -(u ^ 2)⁻¹
          = ∫ u in (2 : ℝ)..x, Chebyshev.theta u * -(u ^ 2)⁻¹ :=
      intervalIntegral.integral_congr fun u _ ↦ by rw [hθ_partial u]
    rw [hint_eq', hint_eq]
    -- Final algebra: x⁻¹ * θ(x) - (-∫…) = θ(x)/x + ∫…
    rw [sub_neg_eq_add]
    congr 1
    rw [mul_comm, div_eq_mul_inv]
  · -- Differentiability of `t ↦ t⁻¹` on `[2, x]` (avoiding `0`).
    intro z hz
    have hz_pos : (0 : ℝ) < z := by
      have : (2 : ℝ) ≤ z := hz.1
      linarith
    have hzne : z ≠ 0 := ne_of_gt hz_pos
    exact differentiableAt_inv hzne
  · -- Integrability of `deriv (·⁻¹) = -(·^2)⁻¹` on `[2, x]`.
    -- Since `deriv_inv'` gives an unconditional pointwise equality of functions,
    -- we can rewrite the integrand and reduce to continuity of `-(z^2)⁻¹`.
    have hderiv_fn : (deriv fun t : ℝ ↦ t⁻¹) = fun z : ℝ ↦ -(z ^ 2)⁻¹ := deriv_inv'
    rw [hderiv_fn]
    refine ContinuousOn.integrableOn_Icc ?_
    intro z hz
    have hz_pos : (0 : ℝ) < z := by
      have : (2 : ℝ) ≤ z := hz.1
      linarith
    have hzne : z ≠ 0 := ne_of_gt hz_pos
    have hz2 : z ^ 2 ≠ 0 := pow_ne_zero 2 hzne
    exact ContinuousAt.continuousWithinAt (by fun_prop)

/-! ### Mertens' First Theorem (M1)

**Mertens' First Theorem.** The log-weighted partial sum of prime reciprocals
satisfies

$$\sum_{p \le x} \frac{\log p}{p} \;=\; \log x \,+\, O(1).$$

Equivalently, there exists a constant `C` such that
`|mertensLog x - log x| ≤ C` for all `x ≥ 2`.

This is the analytic input to Mertens' Third Theorem.

## Status (audit 2026-05-16)

The Abel-summation pivot `mertensLog_eq_theta_div_x_add_integral` (lines ~322–403,
`sorry`-free) gives the identity
`mertensLog x = θ(x)/x + ∫₂ˣ θ(t)/t² dt`.

**Route A (Abel + Chebyshev) is blocked by current Mathlib.**
Mathlib's elementary Chebyshev bounds (`Chebyshev.theta_le_log4_mul_x`,
`Chebyshev.theta_ge'`) sandwich `θ(t)` between `c·t` and `(log 4)·t` with
`c < 1 < log 4`. Plugging into the pivot gives
`mertensLog x = c'·log x + O(1)` for **some** `c' ∈ [log 2 / 2, log 4]` —
the coefficient is not pinned to `1`. To get coefficient exactly `1` requires
either:

* **PNT-level estimate** `θ(t) = t + O(t/log t)`, currently **absent** from
  Mathlib (no `Chebyshev.theta_sub_self_isBigO`, no
  `tendsto_theta_div_id_atTop`); or
* **Selberg-style elementary refinement** (Erdős–Selberg 1949), not in Mathlib.

**Route B (Legendre + Stirling, Mertens' original 1874 proof)** is viable but
unimplemented. Building blocks available in Mathlib + this file:

1. `log_factorial_eq_sum_vonMangoldt_mul_floor` (this file, **closed**):
   `log N! = ∑_{d ≤ N} Λ(d) · ⌊N/d⌋`.
2. `Stirling.le_log_factorial_stirling` (Mathlib, effective lower bound):
   `n·log n - n + log n / 2 + log(2π)/2 ≤ log n!`.
3. `Stirling.log_stirlingSeq'_antitone` + `log_stirlingSeq_bounded_by_constant`
   give an effective two-sided bound: `log(stirlingSeq n)` is bounded by an
   explicit constant for all `n ≥ 1`, hence
   `|log n! - (n·log n - n + (1/2)·log(2n))| ≤ C_S` (explicit `C_S`).
4. `Chebyshev.psi_le_const_mul_self` (Mathlib): `ψ(n) ≤ C·n`.

The proof then runs:

  (a) Use Stirling to get `log n! = n·log n - n + O(log n)`.
  (b) Split the Legendre sum: `∑ Λ(d) ⌊n/d⌋ = n·∑ Λ(d)/d - ∑ Λ(d)·{n/d}`,
      where the fractional-part sum is `≤ ψ(n) ≤ C·n`.
  (c) Bound the prime-power tail:
      `∑_{d ≤ n} Λ(d)/d - mertensLog n = ∑_{p^k ≤ n, k ≥ 2} (log p)/p^k`,
      which is `O(1)` (bounded by `∑_p (log p)/(p(p-1))`, convergent).
  (d) Divide by `n`, take `n = ⌊x⌋`, control `log x - log ⌊x⌋ ≤ log 2`.

Estimated effort: ~2 focused sessions (Stirling effective, fractional-part
identity, prime-power tail bound, real-to-integer reduction).

#### Sub-lemma 1: Stirling's approximation in effective form

We extract an effective version of Stirling's approximation from Mathlib's
`Stirling.log_stirlingSeq_formula`, `log_stirlingSeq_bounded_by_constant`, and
`log_stirlingSeq'_antitone`. The bound `|log n! − (n log n − n)| ≤ C · (1 + log n)`
is what Sub-lemma 3 (the Legendre+Stirling floor identity) consumes.
-/

/-- **Effective Stirling bound (PT form).** There exists a constant `C ≥ 0`
such that for every `n ≥ 1`,

  `|log n! − (n · log n − n)| ≤ C · (1 + log n)`.

This is the variant of Stirling's approximation needed for Mertens M1 via the
Legendre + Stirling route: the error term is `O(log n)`, which becomes `o(n)`
after dividing by `n` in the Mertens estimate. -/
lemma stirling_log_factorial_effective :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ n : ℕ, 1 ≤ n →
      |Real.log ((n.factorial : ℕ) : ℝ) - ((n : ℝ) * Real.log n - n)| ≤ C * (1 + Real.log n) := by
  -- Set the bounds for log(stirlingSeq m) for m ≥ 1.
  -- Lower: 1 - 1/12 - log 2 / 2 ≤ log(stirlingSeq (k+1)) (k ≥ 0).
  -- Upper: log(stirlingSeq (k+1)) ≤ log(stirlingSeq 1) = 1 - log 2 / 2 (antitone).
  set Clo : ℝ := 1 - 12⁻¹ - Real.log 2 / 2 with hClodef
  set Chi : ℝ := 1 - Real.log 2 / 2 with hChidef
  have hlog2_pos : 0 < Real.log 2 := Real.log_pos one_lt_two
  -- M := max |Clo| |Chi| bounds |log(stirlingSeq (k+1))|.
  set M : ℝ := |Clo| + |Chi| + 1 with hMdef
  have hM_nonneg : 0 ≤ M := by positivity
  -- The Stirling constant from the formula: |log n! - (n log n - n + log(2n)/2)| ≤ M.
  -- Bound from `log_stirlingSeq_formula`: log(stirlingSeq n) = log n! - log(2n)/2 - n*log(n/e).
  -- For n ≥ 1: n*log(n/e) = n*log n - n, so
  --   log n! = log(stirlingSeq n) + log(2n)/2 + n*log n - n.
  -- and |log(stirlingSeq n)| ≤ M (using bounds above).
  -- The full error |log n! - (n log n - n)| ≤ |log(stirlingSeq n)| + log(2n)/2 ≤ M + log n/2 + log 2/2.
  -- Take C := M + log 2 / 2 + 1/2.
  refine ⟨M + Real.log 2 / 2 + 1, ?_, ?_⟩
  · positivity
  intro n hn
  -- Case n ≥ 1: n = k+1 for some k ≥ 0.
  obtain ⟨k, rfl⟩ : ∃ k, n = k + 1 := ⟨n - 1, by omega⟩
  -- log_stirlingSeq_formula n: log(stirlingSeq n) = log n! - log(2n)/2 - n*log(n/e)
  have hform := Stirling.log_stirlingSeq_formula (k + 1)
  -- Compute n*log(n/e) = n*log n - n (since log(n/e) = log n - log e = log n - 1).
  have hkpos : (0 : ℝ) < (k + 1 : ℕ) := by exact_mod_cast Nat.succ_pos k
  have hkne : ((k + 1 : ℕ) : ℝ) ≠ 0 := ne_of_gt hkpos
  have hlog_div : Real.log (((k + 1 : ℕ) : ℝ) / Real.exp 1)
      = Real.log ((k + 1 : ℕ) : ℝ) - 1 := by
    rw [Real.log_div hkne (Real.exp_ne_zero 1), Real.log_exp]
  have h_n_log_div :
      ((k + 1 : ℕ) : ℝ) * Real.log (((k + 1 : ℕ) : ℝ) / Real.exp 1)
        = ((k + 1 : ℕ) : ℝ) * Real.log ((k + 1 : ℕ) : ℝ) - ((k + 1 : ℕ) : ℝ) := by
    rw [hlog_div]; ring
  -- log(stirlingSeq (k+1)) = log n! - log(2n)/2 - (n*log n - n).
  rw [h_n_log_div] at hform
  -- Bounds: Clo ≤ log(stirlingSeq (k+1)) ≤ Chi.
  have hlo : Clo ≤ Real.log (Stirling.stirlingSeq (k + 1)) := by
    have := Stirling.log_stirlingSeq_bounded_by_constant k
    -- Mathlib: 1 - 12⁻¹ - log 2 / 2 ≤ log(stirlingSeq (n+1))
    simpa [hClodef] using this
  have hhi : Real.log (Stirling.stirlingSeq (k + 1)) ≤ Chi := by
    -- antitone: log(stirlingSeq (m+1)) ≤ log(stirlingSeq 1) for m ≥ 0.
    have hanti := Stirling.log_stirlingSeq'_antitone (Nat.zero_le k)
    -- log_stirlingSeq'_antitone applies to `log ∘ stirlingSeq ∘ succ`,
    -- so hanti has type log(stirlingSeq (k+1)) ≤ log(stirlingSeq 1).
    -- compute Real.log (stirlingSeq 1) = 1 - log 2 / 2.
    have hs1 : Real.log (Stirling.stirlingSeq 1) = 1 - Real.log 2 / 2 := by
      rw [Stirling.stirlingSeq_one, Real.log_div (Real.exp_ne_zero 1)
        (by positivity : (Real.sqrt 2) ≠ 0), Real.log_exp, Real.log_sqrt (by norm_num)]
    -- hanti : Real.log (stirlingSeq (k+1)) ≤ Real.log (stirlingSeq 1)
    have hanti' : Real.log (Stirling.stirlingSeq (k + 1))
        ≤ Real.log (Stirling.stirlingSeq 1) := hanti
    rw [hs1] at hanti'
    exact hanti'.trans_eq hChidef.symm
  -- Derive |log(stirlingSeq (k+1))| ≤ M.
  have hClo_le : |Clo| ≤ M := by
    have : 0 ≤ |Chi| := abs_nonneg _
    have : (0 : ℝ) ≤ 1 := by norm_num
    show |Clo| ≤ |Clo| + |Chi| + 1
    linarith [abs_nonneg Chi]
  have hChi_le : |Chi| ≤ M := by
    show |Chi| ≤ |Clo| + |Chi| + 1
    linarith [abs_nonneg Clo]
  have habs_stir : |Real.log (Stirling.stirlingSeq (k + 1))| ≤ M := by
    rw [abs_le]
    refine ⟨?_, ?_⟩
    · have h1 : -|Clo| ≤ Clo := neg_abs_le _
      have h2 : -M ≤ -|Clo| := by linarith
      linarith
    · calc Real.log (Stirling.stirlingSeq (k + 1))
          ≤ Chi := hhi
        _ ≤ |Chi| := le_abs_self _
        _ ≤ M := hChi_le
  -- From hform: log(stirlingSeq (k+1)) = log n! - log(2n)/2 - n*log n + n
  -- So log n! - (n*log n - n) = log(stirlingSeq (k+1)) + log(2n)/2.
  have heq :
      Real.log ((((k + 1 : ℕ).factorial : ℕ) : ℝ))
        - (((k + 1 : ℕ) : ℝ) * Real.log ((k + 1 : ℕ) : ℝ) - ((k + 1 : ℕ) : ℝ))
      = Real.log (Stirling.stirlingSeq (k + 1))
        + (1 : ℝ) / 2 * Real.log (2 * ((k + 1 : ℕ) : ℝ)) := by
    linarith [hform]
  rw [heq]
  -- |log(stirlingSeq (k+1)) + log(2n)/2| ≤ M + log(2n)/2.
  have h_one_le_2kp1 : (1 : ℝ) ≤ 2 * ((k + 1 : ℕ) : ℝ) := by
    have : (1 : ℝ) ≤ ((k + 1 : ℕ) : ℝ) := by exact_mod_cast Nat.succ_le_succ (Nat.zero_le k)
    linarith
  have h_log_2n_nonneg : 0 ≤ Real.log (2 * ((k + 1 : ℕ) : ℝ)) :=
    Real.log_nonneg h_one_le_2kp1
  have h_log_2n_split : Real.log (2 * ((k + 1 : ℕ) : ℝ))
      = Real.log 2 + Real.log ((k + 1 : ℕ) : ℝ) :=
    Real.log_mul (by norm_num) hkne
  have h_log_n_nonneg : 0 ≤ Real.log ((k + 1 : ℕ) : ℝ) :=
    Real.log_nonneg (by exact_mod_cast Nat.succ_le_succ (Nat.zero_le k))
  have htri : |Real.log (Stirling.stirlingSeq (k + 1))
      + (1 : ℝ) / 2 * Real.log (2 * ((k + 1 : ℕ) : ℝ))|
      ≤ M + (1 : ℝ) / 2 * Real.log (2 * ((k + 1 : ℕ) : ℝ)) := by
    have habsadd := abs_add_le (Real.log (Stirling.stirlingSeq (k + 1)))
      ((1 : ℝ) / 2 * Real.log (2 * ((k + 1 : ℕ) : ℝ)))
    have habs_half : |(1 : ℝ) / 2 * Real.log (2 * ((k + 1 : ℕ) : ℝ))|
        = (1 : ℝ) / 2 * Real.log (2 * ((k + 1 : ℕ) : ℝ)) := by
      rw [abs_of_nonneg]
      positivity
    rw [habs_half] at habsadd
    linarith
  calc |Real.log (Stirling.stirlingSeq (k + 1))
        + (1 : ℝ) / 2 * Real.log (2 * ((k + 1 : ℕ) : ℝ))|
      ≤ M + (1 : ℝ) / 2 * Real.log (2 * ((k + 1 : ℕ) : ℝ)) := htri
    _ = M + (Real.log 2 + Real.log ((k + 1 : ℕ) : ℝ)) / 2 := by
        rw [h_log_2n_split]; ring
    _ ≤ (M + Real.log 2 / 2 + 1) * (1 + Real.log ((k + 1 : ℕ) : ℝ)) := by
        have hL2 : 0 ≤ Real.log 2 := le_of_lt hlog2_pos
        nlinarith [h_log_n_nonneg, hL2, hM_nonneg]

/-! #### Sub-lemma 4: real-to-integer reduction

Reduces the real form of Mertens M1 to its integer form via the identity
`mertensLog x = mertensLog ⌊x⌋₊` (immediate from the definition) and the
elementary `|log x − log ⌊x⌋₊| ≤ log 2` for `x ≥ 2`. -/

/-- The `mertensLog` partial sum depends on `x` only through `⌊x⌋₊`. -/
lemma mertensLog_eq_mertensLog_floor (x : ℝ) :
    mertensLog x = mertensLog ((⌊x⌋₊ : ℕ) : ℝ) := by
  unfold mertensLog
  congr 1
  rw [Nat.floor_natCast]

/-- For `x ≥ 2`, `log x − log ⌊x⌋₊` is bounded by `log 2`. -/
lemma abs_log_sub_log_floor_le {x : ℝ} (hx : 2 ≤ x) :
    |Real.log x - Real.log ((⌊x⌋₊ : ℕ) : ℝ)| ≤ Real.log 2 := by
  have h_floor_pos : (1 : ℝ) ≤ (⌊x⌋₊ : ℕ) := by
    have : (2 : ℕ) ≤ ⌊x⌋₊ := by
      rw [Nat.le_floor_iff (by linarith : (0 : ℝ) ≤ x)]
      exact_mod_cast hx
    exact_mod_cast (by omega : 1 ≤ ⌊x⌋₊)
  have h_floor_le_x : ((⌊x⌋₊ : ℕ) : ℝ) ≤ x := Nat.floor_le (by linarith)
  have h_x_lt_succ : x < ((⌊x⌋₊ + 1 : ℕ) : ℝ) := by
    have := Nat.lt_floor_add_one x
    push_cast
    push_cast at this
    linarith
  -- log ⌊x⌋₊ ≤ log x ≤ log(⌊x⌋₊ + 1) ≤ log(2 · ⌊x⌋₊) = log 2 + log ⌊x⌋₊.
  have h_floor_pos_r : (0 : ℝ) < ((⌊x⌋₊ : ℕ) : ℝ) := by linarith
  have hlog_ge : Real.log ((⌊x⌋₊ : ℕ) : ℝ) ≤ Real.log x :=
    Real.log_le_log h_floor_pos_r h_floor_le_x
  have hsucc_le : ((⌊x⌋₊ + 1 : ℕ) : ℝ) ≤ 2 * ((⌊x⌋₊ : ℕ) : ℝ) := by
    push_cast; linarith
  have hlog_le : Real.log x ≤ Real.log (2 * ((⌊x⌋₊ : ℕ) : ℝ)) := by
    apply Real.log_le_log (by linarith : (0 : ℝ) < x)
    linarith
  have hlog_2n_split : Real.log (2 * ((⌊x⌋₊ : ℕ) : ℝ))
      = Real.log 2 + Real.log ((⌊x⌋₊ : ℕ) : ℝ) :=
    Real.log_mul (by norm_num) (ne_of_gt h_floor_pos_r)
  rw [abs_le]
  refine ⟨?_, ?_⟩
  · -- −log 2 ≤ log x − log ⌊x⌋₊, equivalently log ⌊x⌋₊ − log 2 ≤ log x.
    have hlog2_nonneg : 0 ≤ Real.log 2 := le_of_lt (Real.log_pos one_lt_two)
    linarith
  · -- log x − log ⌊x⌋₊ ≤ log 2.
    rw [hlog_2n_split] at hlog_le
    linarith

/-- **Real-to-integer reduction.** Given a uniform bound on `|mertensLog N − log N|`
for natural numbers `N ≥ 2`, the same bound (shifted by `log 2`) holds for all
real `x ≥ 2`. -/
lemma mertens_M1_of_nat
    (h : ∃ C : ℝ, ∀ N : ℕ, 2 ≤ N → |mertensLog ((N : ℕ) : ℝ) - Real.log ((N : ℕ) : ℝ)| ≤ C) :
    ∃ C : ℝ, ∀ x : ℝ, 2 ≤ x → |mertensLog x - Real.log x| ≤ C := by
  obtain ⟨C, hC⟩ := h
  refine ⟨C + Real.log 2, ?_⟩
  intro x hx
  have h_floor_ge : (2 : ℕ) ≤ ⌊x⌋₊ := by
    rw [Nat.le_floor_iff (by linarith : (0 : ℝ) ≤ x)]
    exact_mod_cast hx
  have hbnd := hC ⌊x⌋₊ h_floor_ge
  have hreduce := mertensLog_eq_mertensLog_floor x
  have habs := abs_log_sub_log_floor_le hx
  -- |mertensLog x - log x| = |mertensLog ⌊x⌋₊ - log x|
  --                       ≤ |mertensLog ⌊x⌋₊ - log ⌊x⌋₊| + |log ⌊x⌋₊ - log x|
  --                       ≤ C + log 2.
  calc |mertensLog x - Real.log x|
      = |mertensLog ((⌊x⌋₊ : ℕ) : ℝ) - Real.log x| := by rw [hreduce]
    _ = |mertensLog ((⌊x⌋₊ : ℕ) : ℝ) - Real.log ((⌊x⌋₊ : ℕ) : ℝ)
          + (Real.log ((⌊x⌋₊ : ℕ) : ℝ) - Real.log x)| := by ring_nf
    _ ≤ |mertensLog ((⌊x⌋₊ : ℕ) : ℝ) - Real.log ((⌊x⌋₊ : ℕ) : ℝ)|
          + |Real.log ((⌊x⌋₊ : ℕ) : ℝ) - Real.log x| :=
        abs_add_le _ _
    _ ≤ C + Real.log 2 := by
        have h1 : |Real.log ((⌊x⌋₊ : ℕ) : ℝ) - Real.log x|
            = |Real.log x - Real.log ((⌊x⌋₊ : ℕ) : ℝ)| := abs_sub_comm _ _
        rw [h1]
        linarith

/-! ### Sub-lemma 2: prime-power tail bounded

The prime-power tail (the contribution of `Λ(p^k) = log p` for `k ≥ 2`)
to `∑ Λ(n)/n` is **absolutely summable**. This is the second analytic
input to Mertens M1 via the Legendre + Stirling route.

The proof re-derives the unnamed (private) Mathlib machinery from
`Mathlib.NumberTheory.LSeries.PrimesInAP` (lemmas `F''_le` and `summable_F''`,
lines 178–220 of that file), repackaged here as public lemmas usable from
this file. The bound is uniform in `n`:

  `(Λ(p^{k+2}) / p^{k+2}) = log p / p^{k+2} ≤ 2 · p^{-(k + 3/2)}`,

and the double sum over `(p, k) ∈ Primes × ℕ` is dominated by the product
of a geometric series in `k` and the `p`-series `∑ 1/p^{3/2}`.
-/

open ArithmeticFunction in
/-- Pointwise bound: `Λ(p^{k+2}) / p^{k+2} ≤ 2 · p^{-(k + 3/2)}` for every prime `p`
    and every `k : ℕ`. Repackaging of Mathlib's private `F''_le`. -/
private lemma vonMangoldt_prime_pow_div_le (p : Nat.Primes) (k : ℕ) :
    Real.log (p.val : ℝ) * ((p.val : ℝ)⁻¹) ^ (k + 2)
      ≤ 2 * ((p.val : ℝ)⁻¹) ^ (k + 3 / 2 : ℝ) := by
  calc Real.log (p.val : ℝ) * ((p.val : ℝ)⁻¹) ^ (k + 2)
      ≤ (p.val : ℝ) ^ (1 / 2 : ℝ) / (1 / 2) * ((p.val : ℝ)⁻¹) ^ (k + 2) :=
        mul_le_mul_of_nonneg_right (Real.log_le_rpow_div p.val.cast_nonneg one_half_pos)
          (pow_nonneg (inv_nonneg_of_nonneg (Nat.cast_nonneg ↑p)) (k + 2))
    _ = 2 * ((p.val : ℝ)⁻¹) ^ (-1 / 2 : ℝ) * ((p.val : ℝ)⁻¹) ^ (k + 2) := by
        simp only [← div_mul, div_one, mul_comm, neg_div, Real.inv_rpow p.val.cast_nonneg,
          ← Real.rpow_neg p.val.cast_nonneg, neg_neg]
    _ = 2 * ((p.val : ℝ)⁻¹) ^ (k + 3 / 2 : ℝ) := by
        rw [mul_assoc, ← Real.rpow_natCast (((p.val : ℝ))⁻¹) (k + 2),
          ← Real.rpow_add <| by have := p.prop.pos; positivity, Nat.cast_add, Nat.cast_two,
          add_comm, add_assoc]
        norm_num

set_option maxHeartbeats 400000 in
/-- Summability of `(p, k) ↦ Λ(p^{k+2}) / p^{k+2}` over `Nat.Primes × ℕ`.
    Repackaging of Mathlib's private `summable_F''`. -/
private lemma summable_prime_pow_tail_prod :
    Summable (fun pk : Nat.Primes × ℕ =>
      Real.log (pk.1.val : ℝ) * ((pk.1.val : ℝ)⁻¹) ^ (pk.2 + 2)) := by
  have hp₀ (p : Nat.Primes) : 0 < (p.val : ℝ)⁻¹ :=
    inv_pos_of_pos (Nat.cast_pos.mpr p.prop.pos)
  have hp₁ (p : Nat.Primes) : (p.val : ℝ)⁻¹ < 1 :=
    (inv_lt_one₀ <| mod_cast p.prop.pos).mpr <| Nat.one_lt_cast.mpr <| p.prop.one_lt
  -- Bound the function by the rpow form which factors as a geometric × p-series.
  suffices Summable fun (pk : Nat.Primes × ℕ) ↦ ((pk.1.val : ℝ)⁻¹) ^ (pk.2 + 3 / 2 : ℝ) by
    refine (Summable.mul_left 2 this).of_nonneg_of_le (fun pk ↦ ?_)
      (fun pk ↦ vonMangoldt_prime_pow_div_le pk.1 pk.2)
    have hpos : (0 : ℝ) ≤ ((pk.1.val : ℝ))⁻¹ := le_of_lt (hp₀ pk.1)
    have h1 : (0 : ℝ) ≤ Real.log (pk.1.val : ℝ) :=
      Real.log_nonneg (by exact_mod_cast pk.1.prop.one_lt.le)
    exact mul_nonneg h1 (pow_nonneg hpos _)
  conv => enter [1, pk]; rw [Real.rpow_add <| hp₀ pk.1, Real.rpow_natCast]
  refine (summable_prod_of_nonneg (fun _ ↦ by positivity)).mpr ⟨(fun p ↦ ?_), ?_⟩
  · dsimp only
    exact Summable.mul_right _ <| summable_geometric_of_lt_one (hp₀ p).le (hp₁ p)
  · dsimp only
    conv => enter [1, p]; rw [tsum_mul_right, tsum_geometric_of_lt_one (hp₀ p).le (hp₁ p)]
    -- Now: ∑' p, (p⁻¹)^(3/2) * (1 - p⁻¹)⁻¹ ≤ 2 · ∑' p, p^(-3/2) — convergent p-series.
    -- We use the Subtype.val injection from Nat.Primes to ℕ to lift summable_nat_rpow.
    have h_summable_primes : Summable (fun p : Nat.Primes ↦ 2 * ((p.val : ℝ)⁻¹) ^ (3 / 2 : ℝ)) := by
      have hbase : Summable (fun n : ℕ ↦ (n : ℝ) ^ (-(3 / 2 : ℝ))) :=
        Real.summable_nat_rpow.mpr (by norm_num : -(3 / 2 : ℝ) < -1)
      have hbase2 : Summable (fun n : ℕ ↦ 2 * (n : ℝ) ^ (-(3 / 2 : ℝ))) :=
        hbase.mul_left 2
      have hinj : Function.Injective (fun p : Nat.Primes => (p.val : ℕ)) :=
        fun p q h => Subtype.ext h
      have h_comp : Summable (fun p : Nat.Primes ↦ 2 * ((p.val : ℕ) : ℝ) ^ (-(3 / 2 : ℝ))) :=
        hbase2.comp_injective hinj
      refine h_comp.congr ?_
      intro p
      rw [Real.inv_rpow p.val.cast_nonneg, Real.rpow_neg p.val.cast_nonneg]
    refine h_summable_primes.of_nonneg_of_le (fun p ↦ ?_) (fun p ↦ ?_)
    · positivity [sub_pos.mpr (hp₁ p)]
    · -- Goal: (p⁻¹)^(3/2) * (1 - p⁻¹)⁻¹ ≤ 2 * (p⁻¹)^(3/2)
      -- i.e., (1 - p⁻¹)⁻¹ ≤ 2, since (p⁻¹)^(3/2) ≥ 0.
      have h_factor_nn : 0 ≤ ((p.val : ℝ)⁻¹) ^ (3 / 2 : ℝ) :=
        Real.rpow_nonneg (le_of_lt (hp₀ p)) _
      have h_one_sub_pos : 0 < (1 - (p.val : ℝ)⁻¹) := sub_pos.mpr (hp₁ p)
      have h_inv_le_two : (1 - (p.val : ℝ)⁻¹)⁻¹ ≤ 2 := by
        rw [inv_le_comm₀ h_one_sub_pos zero_lt_two, le_sub_comm,
          show (1 : ℝ) - 2⁻¹ = 2⁻¹ by norm_num,
          inv_le_inv₀ (mod_cast p.prop.pos) zero_lt_two]
        exact Nat.ofNat_le_cast.mpr p.prop.two_le
      -- The current goal shape is `(1 - p⁻¹)⁻¹ * (p⁻¹)^(3/2) ≤ 2 * (p⁻¹)^(3/2)`.
      exact mul_le_mul_of_nonneg_right h_inv_le_two h_factor_nn

/-- Auxiliary function `F₀` from Mathlib's `PrimesInAP` private machinery,
    repackaged here. It is `0` on primes and on non-prime-powers, and equals
    `Λ(n) / n` on prime powers `p^k` with `k ≥ 2`. -/
private noncomputable def primePowerTailFn (n : ℕ) : ℝ :=
  (if n.Prime then 0 else ArithmeticFunction.vonMangoldt n) / (n : ℝ)

open ArithmeticFunction in
private lemma primePowerTailFn_nonneg (n : ℕ) : 0 ≤ primePowerTailFn n := by
  unfold primePowerTailFn
  split_ifs with h
  · simp
  · positivity [vonMangoldt_nonneg (n := n)]

set_option maxHeartbeats 400000 in
/-- **Summability of the prime-power tail.** The function `n ↦ Λ(n)/n`
    restricted to non-primes is summable. Equivalently, the contribution of
    prime powers `p^k` with `k ≥ 2` to `∑ Λ(n)/n` is finite. -/
private lemma summable_primePowerTailFn : Summable primePowerTailFn := by
  -- Strategy: factor through the equivalence
  --     Nat.Primes × ℕ ≃ {n // IsPrimePow n}
  -- and the injective shift `(p, j) ↦ (p, j + 1)` (which excludes `(p, 0)`,
  -- corresponding to `p^1 = p` — a prime, where `primePowerTailFn = 0`).
  have hF0_on_prime (p : Nat.Primes) : primePowerTailFn p.val = 0 := by
    simp only [primePowerTailFn, p.prop, ↓reduceIte, zero_div]
  have hF0_off_pp : ∀ n : ℕ, ¬ IsPrimePow n → primePowerTailFn n = 0 := by
    intro n hn
    have hΛ : ArithmeticFunction.vonMangoldt n = 0 :=
      ArithmeticFunction.vonMangoldt_eq_zero_iff.mpr hn
    have hnp : ¬ n.Prime := fun hp => hn hp.isPrimePow
    simp [primePowerTailFn, hΛ, hnp]
  -- Step 1: Show summability on the subtype `{n // IsPrimePow n}`, then transport.
  -- The function vanishes outside the IsPrimePow set, so this is enough.
  suffices h_sub : Summable (primePowerTailFn ∘ (Subtype.val : {n : ℕ // IsPrimePow n} → ℕ)) by
    have h_ind : Summable (({n : ℕ | IsPrimePow n}).indicator primePowerTailFn) :=
      (summable_subtype_iff_indicator (f := primePowerTailFn)
        (s := {n : ℕ | IsPrimePow n})).mp h_sub
    refine h_ind.congr ?_
    intro n
    -- Goal: indicator primePowerTailFn n = primePowerTailFn n.
    exact Set.indicator_apply_eq_self.mpr (fun hn => hF0_off_pp n hn)
  -- Step 2: Use the inverse equivalence `prodNatEquiv` to land on `Nat.Primes × ℕ`.
  rw [← Nat.Primes.prodNatEquiv.summable_iff]
  -- Goal: Summable (fun pk => primePowerTailFn (Nat.Primes.prodNatEquiv pk).val)
  -- Use coe_prodNatEquiv_apply : (prodNatEquiv (p, k) : ℕ) = p^(k+1).
  set g : Nat.Primes × ℕ → ℝ :=
    fun pk => primePowerTailFn ((pk.1.val : ℕ) ^ (pk.2 + 1)) with hg_def
  have h_goal_eq : ((primePowerTailFn ∘ Subtype.val) ∘ Nat.Primes.prodNatEquiv) = g := by
    funext pk
    obtain ⟨p, k⟩ := pk
    simp only [Function.comp_apply, g, Nat.Primes.coe_prodNatEquiv_apply]
  rw [h_goal_eq]
  -- Goal: Summable g
  -- The shift `Prod.map id (· + 1)` is injective, its image is `{pk | pk.2 ≥ 1}`,
  -- complement is `{pk | pk.2 = 0}`, where g vanishes (since p^1 = p is prime).
  have h_inj : Function.Injective
      ((Prod.map _root_.id (· + 1)) : Nat.Primes × ℕ → Nat.Primes × ℕ) :=
    Function.Injective.prodMap (fun ⦃_ _⦄ a ↦ a) (fun _ _ h => by omega)
  have h_zero_outside : ∀ pk ∉ Set.range
      ((Prod.map _root_.id (· + 1)) : Nat.Primes × ℕ → Nat.Primes × ℕ),
      g pk = 0 := by
    intro pk hpk
    have hpk2 : pk.2 = 0 := by
      by_contra hne
      apply hpk
      refine ⟨(pk.1, pk.2 - 1), ?_⟩
      simp only [Prod.map_apply, id_eq]
      ext
      · rfl
      · dsimp; omega
    simp only [g, hpk2, zero_add, pow_one, hF0_on_prime]
  rw [← Function.Injective.summable_iff h_inj h_zero_outside]
  -- Goal: Summable (g ∘ Prod.map id (· + 1))
  -- This equals our summable function from summable_prime_pow_tail_prod.
  refine summable_prime_pow_tail_prod.congr ?_
  intro pj
  obtain ⟨p, k⟩ := pj
  -- Goal: log p * (p⁻¹)^(k+2) = (g ∘ Prod.map id (·+1)) (p, k)
  --                          = g (p, k+1) = primePowerTailFn (p^(k+2)).
  show Real.log (p.val : ℝ) * ((p.val : ℝ)⁻¹) ^ (k + 2)
    = primePowerTailFn ((p.val : ℕ) ^ (k + 1 + 1))
  -- Unfold primePowerTailFn at p^(k+2) and reduce.
  have h_pow_eq : k + 1 + 1 = k + 2 := by omega
  rw [h_pow_eq]
  have h_ne_prime : ¬ ((p.val : ℕ) ^ (k + 2)).Prime :=
    Nat.Prime.not_prime_pow (by omega : 2 ≤ k + 2)
  have hΛ : ArithmeticFunction.vonMangoldt ((p.val : ℕ) ^ (k + 2))
      = Real.log (p.val : ℝ) := by
    rw [ArithmeticFunction.vonMangoldt_apply_pow (by omega : k + 2 ≠ 0),
      ArithmeticFunction.vonMangoldt_apply_prime p.prop]
  show Real.log (p.val : ℝ) * ((p.val : ℝ)⁻¹) ^ (k + 2)
    = (if ((p.val : ℕ) ^ (k + 2)).Prime then 0
        else ArithmeticFunction.vonMangoldt ((p.val : ℕ) ^ (k + 2))) / (((p.val : ℕ) ^ (k + 2) : ℕ) : ℝ)
  rw [if_neg h_ne_prime, hΛ]
  -- Goal: log p * (p⁻¹)^(k+2) = log p / (p^(k+2) : ℝ)
  push_cast
  rw [div_eq_mul_inv, inv_pow]

set_option maxHeartbeats 400000 in
/-- **Sub-lemma 2 (public form):** the prime-power tail of `∑ Λ(n)/n` is absolutely
    summable. Concretely, `Λ(n)/n − [n prime] · log n / n` is summable. -/
theorem summable_vonMangoldt_prime_power_tail :
    Summable (fun n : ℕ =>
      (ArithmeticFunction.vonMangoldt n
        - (if n.Prime then Real.log n else 0)) / (n : ℝ)) := by
  -- Pointwise, this expression equals `primePowerTailFn n`:
  --   if n.Prime: (log n − log n)/n = 0 = (if n.Prime then 0 else Λ n) / n.
  --   else:      (Λ n − 0)/n = Λ n / n = (if n.Prime then 0 else Λ n) / n.
  refine summable_primePowerTailFn.congr ?_
  intro n
  unfold primePowerTailFn
  by_cases h : n.Prime
  · simp only [h, ↓reduceIte, ArithmeticFunction.vonMangoldt_apply_prime h, sub_self,
      zero_div]
  · simp only [h, ↓reduceIte, sub_zero]

/-! #### Sub-lemma 3: Legendre + Stirling floor identity (integer form)

Combines `log_factorial_eq_sum_vonMangoldt_mul_floor` (Legendre),
`stirling_log_factorial_effective`, `Chebyshev.psi_le_const_mul_self`, and
`summable_vonMangoldt_prime_power_tail` to give a uniform `O(1)` bound on
`|mertensLog N − log N|` for all natural `N ≥ 2`.

The proof uses the decomposition `(N / d : ℕ) = (N : ℝ)/d − ((N % d : ℕ)/d)`
together with `∑ Λ(d) · (N % d / d) ≤ ψ(N)` and the prime-power tail bound. -/

/-- For every `N ≥ 1` and `d ≥ 1`, casting Mathlib's natural-number division
`(N / d : ℕ)` to `ℝ` equals `N/d − (N % d)/d`. -/
private lemma cast_nat_div_eq_real_div_sub_mod (N d : ℕ) (hd : 1 ≤ d) :
    (((N / d : ℕ) : ℕ) : ℝ) = (N : ℝ) / d - ((N % d : ℕ) : ℝ) / d := by
  have hd_pos : (0 : ℝ) < d := by exact_mod_cast hd
  have hd_ne : (d : ℝ) ≠ 0 := ne_of_gt hd_pos
  have h_eq : (N : ℝ) = (d : ℝ) * ((N / d : ℕ) : ℝ) + ((N % d : ℕ) : ℝ) := by
    have := Nat.div_add_mod N d
    have hcast : ((d * (N / d) + N % d : ℕ) : ℝ) = (N : ℝ) := by exact_mod_cast this
    push_cast at hcast
    linarith
  field_simp
  linarith

set_option maxHeartbeats 800000 in
/-- **Sub-lemma 3 (integer form).** There exists a constant `C ≥ 0` such that
    for every natural `N ≥ 2`,

    `|mertensLog N − log N| ≤ C`.

    This is the heart of Mertens' First Theorem: it follows from the
    Legendre identity `log N! = ∑ Λ(d) ⌊N/d⌋`, Stirling's effective form,
    Chebyshev's `ψ`-bound, and the absolute summability of the prime-power
    tail of `∑ Λ(n)/n`. -/
lemma mertensLog_floor_eq_log_floor_add_bounded :
    ∃ C : ℝ, ∀ N : ℕ, 2 ≤ N → |mertensLog ((N : ℕ) : ℝ) - Real.log ((N : ℕ) : ℝ)| ≤ C := by
  -- Set up the constants from the closed sub-lemmas.
  obtain ⟨C_S, hC_S_nn, hStirling⟩ := stirling_log_factorial_effective
  set T : ℝ := ∑' n : ℕ,
      (ArithmeticFunction.vonMangoldt n
        - (if n.Prime then Real.log n else 0)) / (n : ℝ) with hT_def
  have hT_nn : 0 ≤ T := by
    refine tsum_nonneg (fun n => ?_)
    -- The term `(Λ(n) − [n.Prime] log n) / n` equals `primePowerTailFn n`, which is ≥ 0.
    have hpp_nn : 0 ≤ primePowerTailFn n := primePowerTailFn_nonneg n
    rcases Nat.eq_zero_or_pos n with hn0 | hn_pos
    · subst hn0; simp
    · -- For n ≥ 1, split on whether n is prime.
      by_cases hp : n.Prime
      · simp [hp, ArithmeticFunction.vonMangoldt_apply_prime hp]
      · simp [hp]
        unfold primePowerTailFn at hpp_nn
        simp [hp] at hpp_nn
        exact hpp_nn
  set C_psi : ℝ := Real.log 4 + 4 with hC_psi_def
  have hC_psi_nn : 0 ≤ C_psi := by
    have : 0 ≤ Real.log 4 := Real.log_nonneg (by norm_num)
    linarith
  -- Final constant.
  refine ⟨C_S + C_psi + T + 1, ?_⟩
  intro N hN
  -- Notation: we work with N : ℕ, with N ≥ 2.
  have hN1 : 1 ≤ N := by omega
  have hN_pos : 0 < N := by omega
  have hN_R_pos : (0 : ℝ) < (N : ℝ) := by exact_mod_cast hN_pos
  have hN_R_ge_2 : (2 : ℝ) ≤ (N : ℝ) := by exact_mod_cast hN
  have hN_R_ne : (N : ℝ) ≠ 0 := ne_of_gt hN_R_pos
  have h_log_N_nn : 0 ≤ Real.log (N : ℝ) :=
    Real.log_nonneg (by linarith)
  -- Step 1: Legendre identity at N.
  have hLegendre := log_factorial_eq_sum_vonMangoldt_mul_floor N
  -- Step 2: Decompose Λ(d) * (N/d : ℕ) using cast_nat_div_eq_real_div_sub_mod.
  have h_decomp : ∀ d ∈ Finset.Icc 1 N,
      ArithmeticFunction.vonMangoldt d * (((N / d : ℕ) : ℕ) : ℝ)
        = (N : ℝ) * (ArithmeticFunction.vonMangoldt d / d)
          - ArithmeticFunction.vonMangoldt d * (((N % d : ℕ) : ℝ) / d) := by
    intro d hd
    rw [Finset.mem_Icc] at hd
    have hd_pos : 0 < d := hd.1
    have hd_R_pos : (0 : ℝ) < (d : ℝ) := by exact_mod_cast hd_pos
    have hd_R_ne : (d : ℝ) ≠ 0 := ne_of_gt hd_R_pos
    rw [cast_nat_div_eq_real_div_sub_mod N d hd.1]
    field_simp
  rw [Finset.sum_congr rfl h_decomp] at hLegendre
  rw [Finset.sum_sub_distrib, ← Finset.mul_sum] at hLegendre
  -- After this:
  -- log N! = N * ∑ Λ(d)/d - ∑ Λ(d) * ((N % d)/d).
  set S_main : ℝ := ∑ d ∈ Finset.Icc 1 N,
    ArithmeticFunction.vonMangoldt d / (d : ℝ) with hS_main_def
  set S_resid : ℝ := ∑ d ∈ Finset.Icc 1 N,
    ArithmeticFunction.vonMangoldt d * (((N % d : ℕ) : ℝ) / d) with hS_resid_def
  -- So hLegendre : log N! = N * S_main - S_resid.
  -- Step 3: Bound S_resid ≤ ψ(N) ≤ C_psi * N.
  have h_resid_nn_each : ∀ d ∈ Finset.Icc 1 N,
      (0 : ℝ) ≤ ArithmeticFunction.vonMangoldt d * (((N % d : ℕ) : ℝ) / d) := by
    intro d hd
    rw [Finset.mem_Icc] at hd
    have hd_pos : 0 < d := hd.1
    have hd_R_pos : (0 : ℝ) < (d : ℝ) := by exact_mod_cast hd_pos
    have hΛ_nn : 0 ≤ ArithmeticFunction.vonMangoldt d :=
      ArithmeticFunction.vonMangoldt_nonneg
    have hmod_nn : 0 ≤ ((N % d : ℕ) : ℝ) := by exact_mod_cast Nat.zero_le _
    positivity
  have h_resid_each_le : ∀ d ∈ Finset.Icc 1 N,
      ArithmeticFunction.vonMangoldt d * (((N % d : ℕ) : ℝ) / d)
        ≤ ArithmeticFunction.vonMangoldt d := by
    intro d hd
    rw [Finset.mem_Icc] at hd
    have hd_pos : 0 < d := hd.1
    have hd_R_pos : (0 : ℝ) < (d : ℝ) := by exact_mod_cast hd_pos
    have hΛ_nn : 0 ≤ ArithmeticFunction.vonMangoldt d :=
      ArithmeticFunction.vonMangoldt_nonneg
    have hmod_lt : N % d < d := Nat.mod_lt N hd_pos
    have hmod_le : ((N % d : ℕ) : ℝ) ≤ (d : ℝ) := by exact_mod_cast hmod_lt.le
    have h_frac_le_one : ((N % d : ℕ) : ℝ) / d ≤ 1 := by
      rw [div_le_one hd_R_pos]; exact hmod_le
    have h_frac_nn : 0 ≤ ((N % d : ℕ) : ℝ) / d := by
      apply div_nonneg
      · exact_mod_cast Nat.zero_le _
      · exact le_of_lt hd_R_pos
    calc ArithmeticFunction.vonMangoldt d * (((N % d : ℕ) : ℝ) / d)
        ≤ ArithmeticFunction.vonMangoldt d * 1 :=
          mul_le_mul_of_nonneg_left h_frac_le_one hΛ_nn
      _ = ArithmeticFunction.vonMangoldt d := by ring
  -- Bound: S_resid ≤ ψ(N) = ∑_{n ∈ Icc 0 ⌊N⌋₊} Λ n. Here ⌊N⌋₊ = N.
  have h_S_resid_nn : 0 ≤ S_resid :=
    Finset.sum_nonneg h_resid_nn_each
  have h_psi_eq : Chebyshev.psi (N : ℝ) = ∑ d ∈ Finset.Icc 0 N, ArithmeticFunction.vonMangoldt d := by
    rw [Chebyshev.psi_eq_sum_Icc, Nat.floor_natCast]
  have h_psi_le : Chebyshev.psi (N : ℝ) ≤ C_psi * (N : ℝ) :=
    Chebyshev.psi_le_const_mul_self (le_of_lt hN_R_pos)
  -- S_resid ≤ ∑_{d ∈ Icc 1 N} Λ(d) ≤ ∑_{d ∈ Icc 0 N} Λ(d) = ψ(N).
  have h_S_resid_le_psi : S_resid ≤ Chebyshev.psi (N : ℝ) := by
    rw [h_psi_eq]
    have h_sub : ∑ d ∈ Finset.Icc 1 N, ArithmeticFunction.vonMangoldt d
        ≤ ∑ d ∈ Finset.Icc 0 N, ArithmeticFunction.vonMangoldt d := by
      apply Finset.sum_le_sum_of_subset_of_nonneg
      · intro d hd
        rw [Finset.mem_Icc] at *
        omega
      · intros
        exact ArithmeticFunction.vonMangoldt_nonneg
    calc S_resid
        ≤ ∑ d ∈ Finset.Icc 1 N, ArithmeticFunction.vonMangoldt d :=
          Finset.sum_le_sum h_resid_each_le
      _ ≤ ∑ d ∈ Finset.Icc 0 N, ArithmeticFunction.vonMangoldt d := h_sub
  have h_S_resid_le : S_resid ≤ C_psi * (N : ℝ) :=
    le_trans h_S_resid_le_psi h_psi_le
  -- Step 4: Split S_main = mertensLog N + S_tail.
  -- mertensLog N = ∑ p ∈ primesBelow (N+1), log p / p.
  -- Express this as a sum over Icc 1 N of [d.Prime] * log d / d.
  have h_mertensLog_as_sum :
      mertensLog ((N : ℕ) : ℝ)
        = ∑ d ∈ Finset.Icc 1 N,
            (if d.Prime then Real.log (d : ℝ) else 0) / (d : ℝ) := by
    unfold mertensLog
    rw [Nat.floor_natCast]
    -- Goal: ∑ p ∈ primesBelow (N+1), log p / p
    --     = ∑ d ∈ Icc 1 N, [d.Prime] * log d / d.
    rw [Nat.primesBelow_eq_filter_range]
    -- LHS: ∑ p ∈ Finset.range (N+1) with p.Prime, log p / p.
    -- Convert RHS using sum_filter and adjust the index set.
    rw [show (∑ d ∈ Finset.Icc 1 N,
              (if d.Prime then Real.log (d : ℝ) else 0) / (d : ℝ))
          = ∑ d ∈ Finset.Icc 1 N with d.Prime, Real.log (d : ℝ) / (d : ℝ) from ?_]
    · -- LHS = ∑ p ∈ (range (N+1)).filter Prime, log p / p
      -- RHS = ∑ d ∈ (Icc 1 N).filter Prime, log d / d
      -- range (N+1) = {0,...,N}; Icc 1 N = {1,...,N}. They differ only at d=0,
      -- which is not prime, so the filters agree.
      apply Finset.sum_congr ?_ (fun _ _ => rfl)
      ext d
      simp only [Finset.mem_filter, Finset.mem_range, Finset.mem_Icc]
      constructor
      · rintro ⟨hd_lt, hd_prime⟩
        refine ⟨⟨hd_prime.one_lt.le, by omega⟩, hd_prime⟩
      · rintro ⟨⟨hd_ge, hd_le⟩, hd_prime⟩
        exact ⟨by omega, hd_prime⟩
    · -- Equality of the two RHS forms.
      rw [Finset.sum_filter]
      apply Finset.sum_congr rfl
      intro d _
      by_cases h : d.Prime
      · simp [h]
      · simp [h]
  -- Define S_tail = S_main - mertensLog N. This is the prime-power tail (finite).
  set S_tail : ℝ := ∑ d ∈ Finset.Icc 1 N,
      (ArithmeticFunction.vonMangoldt d
        - (if d.Prime then Real.log (d : ℝ) else 0)) / (d : ℝ) with hS_tail_def
  have h_S_main_decomp : S_main = mertensLog ((N : ℕ) : ℝ) + S_tail := by
    rw [hS_main_def, h_mertensLog_as_sum, hS_tail_def, ← Finset.sum_add_distrib]
    apply Finset.sum_congr rfl
    intro d hd
    rw [Finset.mem_Icc] at hd
    have hd_pos : 0 < d := hd.1
    have hd_R_pos : (0 : ℝ) < (d : ℝ) := by exact_mod_cast hd_pos
    have hd_R_ne : (d : ℝ) ≠ 0 := ne_of_gt hd_R_pos
    field_simp
    ring
  -- S_tail is bounded by T (the full tsum).
  have h_S_tail_nn : 0 ≤ S_tail := by
    apply Finset.sum_nonneg
    intro d hd
    rw [Finset.mem_Icc] at hd
    have hd_pos : 0 < d := hd.1
    have hd_R_pos : (0 : ℝ) < (d : ℝ) := by exact_mod_cast hd_pos
    -- (Λ(d) − [d.Prime] log d)/d = primePowerTailFn d ≥ 0.
    have hpp_nn : 0 ≤ primePowerTailFn d := primePowerTailFn_nonneg d
    by_cases hp : d.Prime
    · simp [hp, ArithmeticFunction.vonMangoldt_apply_prime hp]
    · simp [hp]
      unfold primePowerTailFn at hpp_nn
      simp [hp] at hpp_nn
      exact hpp_nn
  have h_S_tail_le_T : S_tail ≤ T := by
    rw [hT_def]
    -- The summand equals (Λ(n) - [n.Prime] log n) / n; the function is summable
    -- by summable_vonMangoldt_prime_power_tail. The finite sum is ≤ tsum.
    have h_sum_eq : S_tail = ∑ d ∈ Finset.Icc 1 N,
        (fun n : ℕ => (ArithmeticFunction.vonMangoldt n
          - (if n.Prime then Real.log n else 0)) / (n : ℝ)) d := by
      rfl
    rw [h_sum_eq]
    -- Use sum_le_tsum with the summable function and the non-negativity (outside Icc 1 N).
    apply Summable.sum_le_tsum (f := fun n : ℕ =>
      (ArithmeticFunction.vonMangoldt n
        - (if n.Prime then Real.log n else 0)) / (n : ℝ))
    · intro n _
      -- The pointwise nonneg follows from primePowerTailFn_nonneg.
      have hpp_nn : 0 ≤ primePowerTailFn n := primePowerTailFn_nonneg n
      rcases Nat.eq_zero_or_pos n with hn0 | hn_pos
      · subst hn0; simp
      · by_cases hp : n.Prime
        · simp [hp, ArithmeticFunction.vonMangoldt_apply_prime hp]
        · simp [hp]
          unfold primePowerTailFn at hpp_nn
          simp [hp] at hpp_nn
          exact hpp_nn
    · exact summable_vonMangoldt_prime_power_tail
  -- Step 5: Substitute the decomposition into the Legendre identity.
  -- log N! = N * (mertensLog N + S_tail) - S_resid
  --        = N * mertensLog N + N * S_tail - S_resid.
  have hLegendre' : Real.log ((N.factorial : ℕ) : ℝ)
      = (N : ℝ) * mertensLog ((N : ℕ) : ℝ) + (N : ℝ) * S_tail - S_resid := by
    rw [hLegendre, h_S_main_decomp]; ring
  -- Step 6: Use Stirling to relate log N! to N log N - N.
  have hStir := hStirling N hN1
  -- |log N! - (N log N - N)| ≤ C_S (1 + log N).
  -- So: log N! = N log N - N + δ, with |δ| ≤ C_S (1 + log N).
  set δ : ℝ := Real.log ((N.factorial : ℕ) : ℝ) - ((N : ℝ) * Real.log N - N) with hδ_def
  have hδ_abs : |δ| ≤ C_S * (1 + Real.log N) := hStir
  -- From hLegendre' and δ:
  -- N log N - N + δ = N * mertensLog N + N * S_tail - S_resid.
  -- Rearrange: N * mertensLog N = N log N - N - N * S_tail + S_resid + δ.
  -- Divide by N:
  -- mertensLog N = log N - 1 - S_tail + S_resid/N + δ/N.
  -- So: mertensLog N - log N = -1 - S_tail + S_resid/N + δ/N.
  -- |mertensLog N - log N| ≤ 1 + S_tail + S_resid/N + |δ|/N
  --                       ≤ 1 + T + C_psi + C_S * (1 + log N)/N
  --                       ≤ 1 + T + C_psi + C_S       (since (1+log N)/N ≤ 1 for N ≥ 1).
  have h_eqδ : Real.log ((N.factorial : ℕ) : ℝ) = (N : ℝ) * Real.log N - N + δ := by
    rw [hδ_def]; ring
  rw [h_eqδ] at hLegendre'
  -- (N) log N - N + δ = N * mertensLog N + N * S_tail - S_resid
  -- ⟹ N * mertensLog N - N * log N = -N + δ + S_resid - N * S_tail
  have h_key : (N : ℝ) * mertensLog ((N : ℕ) : ℝ) - (N : ℝ) * Real.log N
      = -(N : ℝ) + δ + S_resid - (N : ℝ) * S_tail := by linarith
  have h_div : mertensLog ((N : ℕ) : ℝ) - Real.log N
      = (- (N : ℝ) + δ + S_resid - (N : ℝ) * S_tail) / N := by
    have : (N : ℝ) * (mertensLog ((N : ℕ) : ℝ) - Real.log N)
        = -(N : ℝ) + δ + S_resid - (N : ℝ) * S_tail := by
      rw [mul_sub]; linarith
    field_simp at this ⊢
    linarith
  -- Bound the RHS.
  rw [h_div]
  rw [abs_div, abs_of_pos hN_R_pos]
  rw [div_le_iff₀ hN_R_pos]
  -- Goal: |-(N) + δ + S_resid - N * S_tail| ≤ (C_S + C_psi + T + 1) * N.
  have h_log_inv : (1 + Real.log (N : ℝ)) ≤ (N : ℝ) := by
    -- Real.add_one_le_iff: Real.log x ≤ x - 1, so 1 + log x ≤ x for x > 0.
    have h := Real.log_le_sub_one_of_pos hN_R_pos
    linarith
  have h_logN_div : C_S * (1 + Real.log (N : ℝ)) ≤ C_S * (N : ℝ) :=
    mul_le_mul_of_nonneg_left h_log_inv hC_S_nn
  have hδ_le_CSN : |δ| ≤ C_S * (N : ℝ) := le_trans hδ_abs h_logN_div
  -- Use triangle inequalities.
  -- |-(N) + δ + S_resid - N * S_tail| ≤ N + |δ| + |S_resid| + N * |S_tail|.
  have h_abs : |-(N : ℝ) + δ + S_resid - (N : ℝ) * S_tail|
      ≤ (N : ℝ) + |δ| + S_resid + (N : ℝ) * S_tail := by
    have h1 : |-(N : ℝ) + δ + S_resid - (N : ℝ) * S_tail|
        ≤ |-(N : ℝ) + δ + S_resid| + |(N : ℝ) * S_tail| := by
      have := abs_sub (-(N : ℝ) + δ + S_resid) ((N : ℝ) * S_tail)
      exact this
    have h2 : |-(N : ℝ) + δ + S_resid| ≤ |-(N : ℝ) + δ| + |S_resid| :=
      abs_add_le _ _
    have h3 : |-(N : ℝ) + δ| ≤ |-(N : ℝ)| + |δ| :=
      abs_add_le _ _
    have h4 : |-(N : ℝ)| = (N : ℝ) := by rw [abs_neg, abs_of_pos hN_R_pos]
    have h5 : |S_resid| = S_resid := abs_of_nonneg h_S_resid_nn
    have h6 : |(N : ℝ) * S_tail| = (N : ℝ) * S_tail := by
      rw [abs_mul, abs_of_pos hN_R_pos, abs_of_nonneg h_S_tail_nn]
    linarith
  -- Now bound the RHS pieces.
  have h_S_tail_le_T' : (N : ℝ) * S_tail ≤ (N : ℝ) * T :=
    mul_le_mul_of_nonneg_left h_S_tail_le_T (le_of_lt hN_R_pos)
  -- N + C_S * N + C_psi * N + N * T = (1 + C_S + C_psi + T) * N
  -- = (C_S + C_psi + T + 1) * N.
  have h_combine : (N : ℝ) + |δ| + S_resid + (N : ℝ) * S_tail
      ≤ (N : ℝ) + C_S * (N : ℝ) + C_psi * (N : ℝ) + (N : ℝ) * T := by
    linarith
  calc |-(N : ℝ) + δ + S_resid - (N : ℝ) * S_tail|
      ≤ (N : ℝ) + |δ| + S_resid + (N : ℝ) * S_tail := h_abs
    _ ≤ (N : ℝ) + C_S * (N : ℝ) + C_psi * (N : ℝ) + (N : ℝ) * T := h_combine
    _ = (C_S + C_psi + T + 1) * (N : ℝ) := by ring

/-- **Mertens' First Theorem (M1).** The log-weighted partial sum of prime
    reciprocals is `log x + O(1)`. -/
theorem mertens_M1 :
    ∃ C : ℝ, ∀ x : ℝ, 2 ≤ x → |mertensLog x - Real.log x| ≤ C :=
  mertens_M1_of_nat mertensLog_floor_eq_log_floor_add_bounded

/-! ### Mertens' Third Theorem (M3) -/

/-! #### Abel summation pivot for `mertensSum` (foundation for M3)

The classical proof of M3 starts from a second Abel-summation step, this time
expressing `mertensSum x = ∑_{p ≤ x} 1/p` in terms of `mertensLog`:

$$\sum_{p \le x} \frac{1}{p}
   \;=\; \frac{\mathcal{L}(x)}{\log x} \,+\, \int_2^x \frac{\mathcal{L}(t)}{t \log^2 t}\,dt,
\quad \text{where } \mathcal{L}(t) = \sum_{p \le t} \frac{\log p}{p}.$$

This is the analog of `mertensLog_eq_theta_div_x_add_integral`, but with
`f(t) = (\log t)^{-1}` (whose derivative is `-1/(t \log^2 t)` by
`Real.deriv_inv_log`) instead of `f(t) = t^{-1}`.

Combined with `mertens_M1` (`mertensLog x = log x + O(1)`), this identity
yields `mertensSum x = log log x + M_3 + O(1/\log x)`, the headline M3 theorem.
-/

open Asymptotics Filter MeasureTheory in
/-- **Abel-summation pivot for Mertens M3.**

For every `x ≥ 2`,
$$\sum_{p \le x} \frac{1}{p}
   \;=\; \frac{\mathcal{L}(x)}{\log x} \,+\, \int_2^x \frac{\mathcal{L}(t)}{t \log^2 t}\,dt,$$
where `\mathcal{L}(t) = mertensLog t = ∑_{p ≤ t} (log p)/p`.

This is the Mertens analogue for M3 of Mathlib's
`Chebyshev.primeCounting_eq_theta_div_log_add_integral` (applied to the
sequence `(log p / p) · [p prime]` instead of `log p · [p prime]`). It is the
analytic backbone of **Mertens' Third Theorem**. -/
theorem mertensSum_eq_mertensLog_div_log_add_integral {x : ℝ} (hx : 2 ≤ x) :
    mertensSum x = mertensLog x / Real.log x
      + ∫ t in (2 : ℝ)..x, mertensLog t / (t * (Real.log t) ^ 2) := by
  -- Rewrite `mertensSum x` in the form to which Abel summation applies.
  unfold mertensSum
  rw [Nat.primesBelow_eq_filter_range, Nat.range_succ_eq_Icc_zero, Finset.sum_filter]
  -- The Abel-summation "sequence" `b(n) = (log n / n) · [n prime]`.
  let b : ℕ → ℝ := Set.indicator (setOf Nat.Prime) (fun n ↦ Real.log n / n)
  -- Step 1: rewrite the conditional `1/p` as `(log n)⁻¹ * b n`.
  trans ∑ n ∈ Finset.Icc 0 ⌊x⌋₊, (Real.log n)⁻¹ * b n
  · refine Finset.sum_congr rfl fun n _hn ↦ ?_
    by_cases h : Nat.Prime n
    · have hn_pos : (1 : ℝ) < (n : ℝ) := by exact_mod_cast h.one_lt
      have hlog_ne : Real.log n ≠ 0 := by
        have hlog_pos : 0 < Real.log n := Real.log_pos hn_pos
        exact ne_of_gt hlog_pos
      have hn_ne : (n : ℝ) ≠ 0 := by
        have : 0 < (n : ℝ) := by exact_mod_cast h.pos
        exact ne_of_gt this
      simp only [b, Set.indicator_of_mem, Set.mem_setOf_eq, h, if_true]
      field_simp
    · simp [b, h, Set.indicator_of_notMem]
  -- Step 2: apply Abel summation (`sum_mul_eq_sub_integral_mul₁`).
  rw [sum_mul_eq_sub_integral_mul₁ b (f := fun t ↦ (Real.log t)⁻¹)
        (by simp [b]) (by simp [b, Nat.not_prime_one]),
      ← intervalIntegral.integral_of_le hx]
  · -- Step 3: simplify the boundary term and the integrand.
    -- The derivative `(d/dt) (log t)⁻¹ = -t⁻¹ / (log t)^2` (by `Real.deriv_inv_log`).
    have hderiv : ∀ u ∈ Set.uIcc (2 : ℝ) x,
        deriv (fun t : ℝ ↦ (Real.log t)⁻¹) u = -u⁻¹ / (Real.log u) ^ 2 := by
      intro u _; exact Real.deriv_inv_log
    have int_deriv (g : ℝ → ℝ) :
        ∫ u in (2 : ℝ)..x, deriv (fun t : ℝ ↦ (Real.log t)⁻¹) u * g u
          = ∫ u in (2 : ℝ)..x, g u * (-u⁻¹ / (Real.log u) ^ 2) :=
      intervalIntegral.integral_congr fun u hu ↦ by
        rw [hderiv u hu]; ring
    rw [int_deriv]
    -- The boundary term: `f x · ∑ b = (1/log x) · mertensLog x = mertensLog x / log x`.
    have hL_sum : ∑ k ∈ Finset.Icc 0 ⌊x⌋₊, b k = mertensLog x := by
      unfold mertensLog
      rw [Nat.primesBelow_eq_filter_range, Nat.range_succ_eq_Icc_zero, Finset.sum_filter]
      refine Finset.sum_congr rfl fun n _ ↦ ?_
      by_cases h : Nat.Prime n
      · simp [b, h, Set.indicator_of_mem]
      · simp [b, h, Set.indicator_of_notMem]
    -- The integrand: `∑ b over [0..⌊t⌋] = mertensLog t` for every `t`.
    have hL_partial : ∀ t : ℝ,
        ∑ k ∈ Finset.Icc 0 ⌊t⌋₊, b k = mertensLog t := by
      intro t
      unfold mertensLog
      rw [Nat.primesBelow_eq_filter_range, Nat.range_succ_eq_Icc_zero, Finset.sum_filter]
      refine Finset.sum_congr rfl fun n _ ↦ ?_
      by_cases h : Nat.Prime n
      · simp [b, h, Set.indicator_of_mem]
      · simp [b, h, Set.indicator_of_notMem]
    rw [hL_sum]
    -- Now: (log x)⁻¹ * mertensLog x - ∫ mertensLog t * (-t⁻¹ / (log t)^2)
    --     = mertensLog x / log x + ∫ mertensLog t / (t * (log t)^2).
    have hint_eq :
        ∫ u in (2 : ℝ)..x, mertensLog u * (-u⁻¹ / (Real.log u) ^ 2)
          = -∫ u in (2 : ℝ)..x, mertensLog u / (u * (Real.log u) ^ 2) := by
      rw [← intervalIntegral.integral_neg]
      refine intervalIntegral.integral_congr fun u _ ↦ ?_
      simp only [neg_div]
      field_simp
    -- Push the inner-sum-equals-mertensLog rewriting into the integrand:
    have hint_eq' :
        ∫ u in (2 : ℝ)..x, (∑ k ∈ Finset.Icc 0 ⌊u⌋₊, b k) * (-u⁻¹ / (Real.log u) ^ 2)
          = ∫ u in (2 : ℝ)..x, mertensLog u * (-u⁻¹ / (Real.log u) ^ 2) :=
      intervalIntegral.integral_congr fun u _ ↦ by rw [hL_partial u]
    rw [hint_eq', hint_eq]
    -- Final algebra: (log x)⁻¹ * mertensLog x - (-∫…) = mertensLog x / log x + ∫…
    rw [sub_neg_eq_add]
    congr 1
    rw [mul_comm, div_eq_mul_inv]
  · -- Differentiability of `t ↦ (log t)⁻¹` on `[2, x]` (avoiding zeros of `log`).
    intro z hz
    have hz_ge_2 : (2 : ℝ) ≤ z := hz.1
    have hz_pos : (0 : ℝ) < z := by linarith
    have hz_ne : z ≠ 0 := ne_of_gt hz_pos
    have hz_gt_1 : (1 : ℝ) < z := by linarith
    have hlog_pos : 0 < Real.log z := Real.log_pos hz_gt_1
    have hlog_ne : Real.log z ≠ 0 := ne_of_gt hlog_pos
    exact (Real.differentiableAt_log hz_ne).inv hlog_ne
  · -- Integrability of `deriv ((log ·)⁻¹) = -(·)⁻¹ / (log ·)^2` on `[2, x]`.
    -- We rewrite the integrand via the pointwise equality `Real.deriv_inv_log`,
    -- then reduce to continuity on the compact interval.
    refine ContinuousOn.integrableOn_Icc ?_
    intro z hz
    have hz_ge_2 : (2 : ℝ) ≤ z := hz.1
    have hz_pos : (0 : ℝ) < z := by linarith
    have hz_ne : z ≠ 0 := ne_of_gt hz_pos
    have hz_gt_1 : (1 : ℝ) < z := by linarith
    have hlog_pos : 0 < Real.log z := Real.log_pos hz_gt_1
    have hlog_ne : Real.log z ≠ 0 := ne_of_gt hlog_pos
    have hlog2_ne : Real.log z ^ 2 ≠ 0 := pow_ne_zero 2 hlog_ne
    -- Continuity of `deriv (fun t ↦ (log t)⁻¹)` at z = continuity of `-z⁻¹ / (log z)^2` at z.
    refine ContinuousAt.continuousWithinAt ?_
    have h_eq : (deriv fun t : ℝ ↦ (Real.log t)⁻¹)
        = fun y : ℝ ↦ -y⁻¹ / (Real.log y) ^ 2 := by
      funext y; exact Real.deriv_inv_log
    rw [h_eq]
    fun_prop

/-! #### Antiderivative `∫ 1/(t log t) dt = log log t`

This evaluates the principal-term integral that arises when combining the
Abel pivot `mertensSum = mertensLog/log + ∫ mertensLog/(t log² t)` with
`mertens_M1` (`mertensLog t = log t + O(1)`). Since the leading-order
contribution to the integrand is `(log t)/(t log² t) = 1/(t log t)`, its
integral on `[2, x]` is exactly `log log x - log log 2`. -/

open Asymptotics Filter MeasureTheory in
/-- **Antiderivative `∫ 1/(t log t) dt = log log t` on `[2, x]`.**

For every `x ≥ 2`,
$$\int_2^x \frac{1}{t \log t}\,dt \;=\; \log\log x \,-\, \log\log 2.$$

The Lean proof uses the FTC (`integral_deriv_eq_sub'`) together with the chain
rule `(d/dt)(\log\log t) = 1/(t \log t)` (for `t > 1`). -/
theorem integral_one_div_t_log_t {x : ℝ} (hx : 2 ≤ x) :
    ∫ t in (2 : ℝ)..x, 1 / (t * Real.log t)
      = Real.log (Real.log x) - Real.log (Real.log 2) := by
  -- The antiderivative is `F(t) = log (log t)`.
  -- We have `F'(t) = (1/log t) · (1/t) = 1/(t · log t)` for `t > 1`.
  set F : ℝ → ℝ := fun t ↦ Real.log (Real.log t) with hF
  have hderiv : ∀ t ∈ Set.uIcc (2 : ℝ) x, HasDerivAt F (1 / (t * Real.log t)) t := by
    intro t ht
    rw [Set.uIcc_of_le hx] at ht
    have ht_pos : 0 < t := by linarith [ht.1]
    have ht_ne : t ≠ 0 := ne_of_gt ht_pos
    have ht_gt_1 : 1 < t := by linarith [ht.1]
    have hlog_pos : 0 < Real.log t := Real.log_pos ht_gt_1
    have hlog_ne : Real.log t ≠ 0 := ne_of_gt hlog_pos
    -- chain rule
    have h1 : HasDerivAt Real.log t⁻¹ t := Real.hasDerivAt_log ht_ne
    have h2 : HasDerivAt (fun y ↦ Real.log (Real.log y)) (t⁻¹ / Real.log t) t :=
      h1.log hlog_ne
    convert h2 using 1
    field_simp
  -- Apply FTC.
  have h_cont : ContinuousOn (fun t : ℝ ↦ 1 / (t * Real.log t)) (Set.uIcc 2 x) := by
    rw [Set.uIcc_of_le hx]
    intro t ht
    have ht_pos : 0 < t := by linarith [ht.1]
    have ht_ne : t ≠ 0 := ne_of_gt ht_pos
    have ht_gt_1 : 1 < t := by linarith [ht.1]
    have hlog_pos : 0 < Real.log t := Real.log_pos ht_gt_1
    have hlog_ne : Real.log t ≠ 0 := ne_of_gt hlog_pos
    have hprod_ne : t * Real.log t ≠ 0 := mul_ne_zero ht_ne hlog_ne
    refine ContinuousAt.continuousWithinAt ?_
    fun_prop
  have := intervalIntegral.integral_eq_sub_of_hasDerivAt hderiv h_cont.intervalIntegrable
  simpa [F] using this

/-- **Antiderivative `∫ 1/(t log² t) dt = -1/log t` on `[2, x]`.**

For every `x ≥ 2`,
$$\int_2^x \frac{1}{t \log^2 t}\,dt \;=\; \frac{1}{\log 2} \,-\, \frac{1}{\log x}.$$

This is the **error-term integral** for Mertens M3: combining with the
`mertens_M1` bound `|mertensLog t - log t| ≤ C` gives an absolutely convergent
correction, and hence the existence of the limit `M_3`. -/
theorem integral_one_div_t_log_sq {x : ℝ} (hx : 2 ≤ x) :
    ∫ t in (2 : ℝ)..x, 1 / (t * (Real.log t) ^ 2)
      = 1 / Real.log 2 - 1 / Real.log x := by
  -- The antiderivative is `G(t) = -1/log t = -(log t)⁻¹`.
  -- We have `G'(t) = -(d/dt)(log t)⁻¹ = -(-(t⁻¹)/(log t)^2) = 1/(t·(log t)^2)`.
  set G : ℝ → ℝ := fun t ↦ -(Real.log t)⁻¹ with hG
  have hderiv : ∀ t ∈ Set.uIcc (2 : ℝ) x, HasDerivAt G (1 / (t * (Real.log t) ^ 2)) t := by
    intro t ht
    rw [Set.uIcc_of_le hx] at ht
    have ht_pos : 0 < t := by linarith [ht.1]
    have ht_ne : t ≠ 0 := ne_of_gt ht_pos
    have ht_gt_1 : 1 < t := by linarith [ht.1]
    have hlog_pos : 0 < Real.log t := Real.log_pos ht_gt_1
    have hlog_ne : Real.log t ≠ 0 := ne_of_gt hlog_pos
    have hlog2_ne : Real.log t ^ 2 ≠ 0 := pow_ne_zero 2 hlog_ne
    -- d/dt (log t)⁻¹ = -t⁻¹ / (log t)^2 (Real.deriv_inv_log)
    -- So d/dt (-(log t)⁻¹) = t⁻¹ / (log t)^2 = 1/(t · (log t)^2)
    have h_inv_log_diff : DifferentiableAt ℝ (fun y ↦ (Real.log y)⁻¹) t :=
      (Real.differentiableAt_log ht_ne).inv hlog_ne
    have h_deriv_inv_log : deriv (fun y ↦ (Real.log y)⁻¹) t = -t⁻¹ / (Real.log t) ^ 2 :=
      Real.deriv_inv_log
    have h1 : HasDerivAt (fun y ↦ (Real.log y)⁻¹) (-t⁻¹ / (Real.log t) ^ 2) t := by
      rw [← h_deriv_inv_log]; exact h_inv_log_diff.hasDerivAt
    have h2 : HasDerivAt G (-(-t⁻¹ / (Real.log t) ^ 2)) t := h1.neg
    convert h2 using 1
    field_simp
  -- Apply FTC.
  have h_cont : ContinuousOn (fun t : ℝ ↦ 1 / (t * (Real.log t) ^ 2)) (Set.uIcc 2 x) := by
    rw [Set.uIcc_of_le hx]
    intro t ht
    have ht_pos : 0 < t := by linarith [ht.1]
    have ht_ne : t ≠ 0 := ne_of_gt ht_pos
    have ht_gt_1 : 1 < t := by linarith [ht.1]
    have hlog_pos : 0 < Real.log t := Real.log_pos ht_gt_1
    have hlog_ne : Real.log t ≠ 0 := ne_of_gt hlog_pos
    have hlog2_ne : Real.log t ^ 2 ≠ 0 := pow_ne_zero 2 hlog_ne
    have hprod_ne : t * (Real.log t) ^ 2 ≠ 0 := mul_ne_zero ht_ne hlog2_ne
    refine ContinuousAt.continuousWithinAt ?_
    fun_prop
  have := intervalIntegral.integral_eq_sub_of_hasDerivAt hderiv h_cont.intervalIntegrable
  -- `G x - G 2 = -(log x)⁻¹ - (-(log 2)⁻¹) = 1/log 2 - 1/log x`.
  simp only [G] at this
  rw [this]
  ring

/-! #### Cauchy step: `F(x) := mertensSum x - log log x` has uniformly small
oscillation on `[Y, ∞)` for large `Y`.

The proof of `mertens_M3_exists` follows the standard 19th-century recipe:
combine the M3-Abel-pivot `mertensSum_eq_mertensLog_div_log_add_integral` with
M1 (`|mertensLog t - log t| ≤ C`) and the closed-form antiderivatives
`integral_one_div_t_log_t`, `integral_one_div_t_log_sq` to bound
`|F(x) - F(y)|` by `2·C_{M1} / \log y` for `2 ≤ y ≤ x`. The pivot tells us:

```
F(x) - F(y)
  = mertensLog x / log x − mertensLog y / log y
    + ∫_y^x mertensLog(t) / (t · log²t)  dt
    − ∫_y^x 1/(t · log t)  dt.
```

Using `mertensLog t = log t + δ(t)` with `|δ(t)| ≤ C_{M1}` from M1, this
collapses to:

```
F(x) - F(y) = δ(x)/log x − δ(y)/log y + ∫_y^x δ(t) / (t · log²t) dt,
```

and the closed-form `∫_y^x 1/(t log² t) = 1/log y − 1/log x` then gives
`|F(x) - F(y)| ≤ 2 C_{M1} / log y`. -/

-- The full proof of `mertens_M3_exists` is given below, after the
-- supporting helper lemmas `mertensLog_monotone`,
-- `mertensLog_div_intervalIntegrable`, and `mertensSum_sub_log_log_oscillation`.

/-- Monotonicity of `mertensLog` on `ℝ` (as a function of the cutoff). -/
private lemma mertensLog_monotone : Monotone mertensLog := by
  intro a b hab
  unfold mertensLog
  apply Finset.sum_le_sum_of_subset_of_nonneg
  · intro p hp
    rw [Nat.primesBelow_eq_filter_range, Finset.mem_filter] at hp ⊢
    refine ⟨?_, hp.2⟩
    have hfl : ⌊a⌋₊ ≤ ⌊b⌋₊ := Nat.floor_le_floor hab
    exact Finset.mem_range.mpr (lt_of_lt_of_le (Finset.mem_range.mp hp.1) (by omega))
  · intro p hp _
    have hp_prime : Nat.Prime p := by
      rw [Nat.primesBelow_eq_filter_range, Finset.mem_filter] at hp
      exact hp.2
    have hp_pos : 0 < (p : ℝ) := by exact_mod_cast hp_prime.pos
    have hp_ge2 : 2 ≤ (p : ℝ) := by exact_mod_cast hp_prime.two_le
    have hlogp_nn : 0 ≤ Real.log p :=
      Real.log_nonneg (by linarith)
    positivity

/-- Integrability of `mertensLog t / (t · log²t)` on `[2, x]`. -/
private lemma mertensLog_div_intervalIntegrable {x : ℝ} (hx : 2 ≤ x) :
    IntervalIntegrable (fun t => mertensLog t / (t * (Real.log t) ^ 2))
      MeasureTheory.volume 2 x := by
  -- Rewrite as mertensLog t * (1/(t · log²t)).
  have hrw : (fun t => mertensLog t / (t * (Real.log t) ^ 2))
           = (fun t => mertensLog t * (1 / (t * (Real.log t) ^ 2))) := by
    funext t; ring
  rw [hrw]
  refine IntervalIntegrable.mul_continuousOn ?_ ?_
  · exact mertensLog_monotone.intervalIntegrable
  · rw [Set.uIcc_of_le hx]
    intro z hz
    have hz_ge_2 : (2 : ℝ) ≤ z := hz.1
    have hz_pos : (0 : ℝ) < z := by linarith
    have hz_ne : z ≠ 0 := ne_of_gt hz_pos
    have hz_gt_1 : (1 : ℝ) < z := by linarith
    have hlog_pos : 0 < Real.log z := Real.log_pos hz_gt_1
    have hlog_ne : Real.log z ≠ 0 := ne_of_gt hlog_pos
    have hlog2_ne : Real.log z ^ 2 ≠ 0 := pow_ne_zero 2 hlog_ne
    have hprod_ne : z * (Real.log z) ^ 2 ≠ 0 := mul_ne_zero hz_ne hlog2_ne
    refine ContinuousAt.continuousWithinAt ?_
    fun_prop

/-- **Oscillation bound for `F(x) := mertensSum x - log log x`.**

For `2 ≤ y ≤ x`, the difference `F(x) - F(y)` is bounded in absolute value
by `3 C / log y` where `C ≥ 0` is the M1 constant. -/
private lemma mertensSum_sub_log_log_oscillation
    {y x : ℝ} (hy : 2 ≤ y) (hyx : y ≤ x) :
    |(mertensSum x - Real.log (Real.log x)) -
      (mertensSum y - Real.log (Real.log y))| ≤
      3 * (|Classical.choose mertens_M1| + 1) / Real.log y := by
  set C₀ : ℝ := Classical.choose mertens_M1 with hC₀def
  have hC₀spec : ∀ x : ℝ, 2 ≤ x → |mertensLog x - Real.log x| ≤ C₀ :=
    Classical.choose_spec mertens_M1
  -- Use C = |C₀| + 1 to ensure C > 0, which simplifies positivity reasoning.
  set C : ℝ := |C₀| + 1 with hCdef
  have hC_pos : 0 < C := by positivity
  have hC_ge : ∀ t : ℝ, 2 ≤ t → |mertensLog t - Real.log t| ≤ C := by
    intro t ht
    have h1 : |mertensLog t - Real.log t| ≤ C₀ := hC₀spec t ht
    have h2 : C₀ ≤ |C₀| := le_abs_self _
    have h3 : |C₀| ≤ C := by simp [hCdef]
    linarith
  have hx : 2 ≤ x := le_trans hy hyx
  have hy_pos : 0 < y := by linarith
  have hx_pos : 0 < x := by linarith
  have hy_gt_1 : 1 < y := by linarith
  have hx_gt_1 : 1 < x := by linarith
  have hlogy_pos : 0 < Real.log y := Real.log_pos hy_gt_1
  have hlogx_pos : 0 < Real.log x := Real.log_pos hx_gt_1
  have hlogy_ne : Real.log y ≠ 0 := ne_of_gt hlogy_pos
  have hlogx_ne : Real.log x ≠ 0 := ne_of_gt hlogx_pos
  have hlogx_ge_logy : Real.log y ≤ Real.log x := Real.log_le_log hy_pos hyx
  -- Apply the Abel pivot at x and at y.
  have piv_x : mertensSum x = mertensLog x / Real.log x
      + ∫ t in (2 : ℝ)..x, mertensLog t / (t * (Real.log t) ^ 2) :=
    mertensSum_eq_mertensLog_div_log_add_integral hx
  have piv_y : mertensSum y = mertensLog y / Real.log y
      + ∫ t in (2 : ℝ)..y, mertensLog t / (t * (Real.log t) ^ 2) :=
    mertensSum_eq_mertensLog_div_log_add_integral hy
  -- The difference of the integrals = ∫_y^x mertensLog/(t · log²t).
  have hint_y : IntervalIntegrable (fun t => mertensLog t / (t * (Real.log t) ^ 2))
      MeasureTheory.volume 2 y := mertensLog_div_intervalIntegrable hy
  have hint_x : IntervalIntegrable (fun t => mertensLog t / (t * (Real.log t) ^ 2))
      MeasureTheory.volume 2 x := mertensLog_div_intervalIntegrable hx
  have hdiff_int : (∫ t in (2 : ℝ)..x, mertensLog t / (t * (Real.log t) ^ 2))
      - (∫ t in (2 : ℝ)..y, mertensLog t / (t * (Real.log t) ^ 2))
      = ∫ t in y..x, mertensLog t / (t * (Real.log t) ^ 2) :=
    intervalIntegral.integral_interval_sub_left hint_x hint_y
  -- mertensSum x - mertensSum y = (mertensLog x / log x - mertensLog y / log y) + ∫_y^x.
  have hsubpiv : mertensSum x - mertensSum y =
      (mertensLog x / Real.log x - mertensLog y / Real.log y)
      + ∫ t in y..x, mertensLog t / (t * (Real.log t) ^ 2) := by
    rw [piv_x, piv_y]; linarith
  -- log log x - log log y = ∫_y^x 1/(t · log t).
  have hint_log_2y : IntervalIntegrable (fun t => 1 / (t * Real.log t))
      MeasureTheory.volume 2 y := by
    refine ContinuousOn.intervalIntegrable ?_
    rw [Set.uIcc_of_le hy]
    intro z hz
    have hz_ge_2 : (2 : ℝ) ≤ z := hz.1
    have hz_pos : (0 : ℝ) < z := by linarith
    have hz_ne : z ≠ 0 := ne_of_gt hz_pos
    have hz_gt_1 : (1 : ℝ) < z := by linarith
    have hlog_pos : 0 < Real.log z := Real.log_pos hz_gt_1
    have hlog_ne : Real.log z ≠ 0 := ne_of_gt hlog_pos
    have hprod_ne : z * Real.log z ≠ 0 := mul_ne_zero hz_ne hlog_ne
    refine ContinuousAt.continuousWithinAt ?_
    fun_prop
  have hint_log_2x : IntervalIntegrable (fun t => 1 / (t * Real.log t))
      MeasureTheory.volume 2 x := by
    refine ContinuousOn.intervalIntegrable ?_
    rw [Set.uIcc_of_le hx]
    intro z hz
    have hz_ge_2 : (2 : ℝ) ≤ z := hz.1
    have hz_pos : (0 : ℝ) < z := by linarith
    have hz_ne : z ≠ 0 := ne_of_gt hz_pos
    have hz_gt_1 : (1 : ℝ) < z := by linarith
    have hlog_pos : 0 < Real.log z := Real.log_pos hz_gt_1
    have hlog_ne : Real.log z ≠ 0 := ne_of_gt hlog_pos
    have hprod_ne : z * Real.log z ≠ 0 := mul_ne_zero hz_ne hlog_ne
    refine ContinuousAt.continuousWithinAt ?_
    fun_prop
  have hint_log_split : ∫ t in y..x, 1 / (t * Real.log t)
      = Real.log (Real.log x) - Real.log (Real.log y) := by
    have h2x : ∫ t in (2 : ℝ)..x, 1 / (t * Real.log t)
        = Real.log (Real.log x) - Real.log (Real.log 2) :=
      integral_one_div_t_log_t hx
    have h2y : ∫ t in (2 : ℝ)..y, 1 / (t * Real.log t)
        = Real.log (Real.log y) - Real.log (Real.log 2) :=
      integral_one_div_t_log_t hy
    have hsplit := intervalIntegral.integral_interval_sub_left hint_log_2x hint_log_2y
    linarith
  -- Key algebraic identity: F(x) - F(y) = δ(x)/log x - δ(y)/log y + ∫_y^x δ(t)/(t·log²t)
  -- where δ(t) = mertensLog t - log t. We don't simplify symbolically — we just bound.
  -- We use:
  --   mertensLog x / log x = (log x + (mertensLog x - log x))/log x = 1 + δ(x)/log x.
  -- So mertensLog x / log x - mertensLog y / log y = δ(x)/log x - δ(y)/log y + (1 - 1)
  --                                                = δ(x)/log x - δ(y)/log y.
  have hsimpx : mertensLog x / Real.log x = 1 + (mertensLog x - Real.log x) / Real.log x := by
    rw [sub_div]; rw [div_self hlogx_ne]; ring
  have hsimpy : mertensLog y / Real.log y = 1 + (mertensLog y - Real.log y) / Real.log y := by
    rw [sub_div]; rw [div_self hlogy_ne]; ring
  -- ∫_y^x mertensLog/(t·log²t) = ∫_y^x log(t)/(t·log²t) + ∫_y^x δ(t)/(t·log²t)
  --                            = ∫_y^x 1/(t·log t) + ∫_y^x δ(t)/(t·log²t).
  -- We do this split via integrability of each part.
  have hint_log_yx : IntervalIntegrable (fun t => Real.log t / (t * (Real.log t) ^ 2))
      MeasureTheory.volume y x := by
    refine ContinuousOn.intervalIntegrable ?_
    rw [Set.uIcc_of_le hyx]
    intro z hz
    have hz_ge_y : y ≤ z := hz.1
    have hz_ge_2 : (2 : ℝ) ≤ z := le_trans hy hz_ge_y
    have hz_pos : (0 : ℝ) < z := by linarith
    have hz_ne : z ≠ 0 := ne_of_gt hz_pos
    have hz_gt_1 : (1 : ℝ) < z := by linarith
    have hlog_pos : 0 < Real.log z := Real.log_pos hz_gt_1
    have hlog_ne : Real.log z ≠ 0 := ne_of_gt hlog_pos
    have hlog2_ne : Real.log z ^ 2 ≠ 0 := pow_ne_zero 2 hlog_ne
    have hprod_ne : z * (Real.log z) ^ 2 ≠ 0 := mul_ne_zero hz_ne hlog2_ne
    refine ContinuousAt.continuousWithinAt ?_
    fun_prop
  have hint_mertensLog_yx : IntervalIntegrable
      (fun t => mertensLog t / (t * (Real.log t) ^ 2))
      MeasureTheory.volume y x := by
    -- restrict the integrability on [2, x] to [y, x].
    have := hint_x
    have hyx_in_2x : Set.uIcc y x ⊆ Set.uIcc 2 x := by
      rw [Set.uIcc_of_le hyx, Set.uIcc_of_le hx]
      intro t ht
      exact ⟨le_trans hy ht.1, ht.2⟩
    exact this.mono_set hyx_in_2x
  have hint_delta_yx : IntervalIntegrable
      (fun t => (mertensLog t - Real.log t) / (t * (Real.log t) ^ 2))
      MeasureTheory.volume y x := by
    have heq : (fun t => (mertensLog t - Real.log t) / (t * (Real.log t) ^ 2))
        = (fun t => mertensLog t / (t * (Real.log t) ^ 2)
                    - Real.log t / (t * (Real.log t) ^ 2)) := by
      funext t; ring
    rw [heq]
    exact hint_mertensLog_yx.sub hint_log_yx
  -- 1/(t · log t) = log(t)/(t·log²t).
  have hint_rec : (fun t : ℝ => Real.log t / (t * (Real.log t) ^ 2))
                  =ᵐ[MeasureTheory.volume.restrict (Set.Ioc y x)]
                  (fun t : ℝ => 1 / (t * Real.log t)) := by
    refine (MeasureTheory.ae_restrict_iff' measurableSet_Ioc).mpr ?_
    refine Filter.Eventually.of_forall (fun t ht => ?_)
    -- For t ∈ Ioc y x, we have y < t ≤ x, hence 2 ≤ t, log t ≠ 0.
    have ht_ge_y : y < t := ht.1
    have ht_ge_2 : (2 : ℝ) ≤ t := by linarith
    have ht_pos : 0 < t := by linarith
    have ht_ne : t ≠ 0 := ne_of_gt ht_pos
    have ht_gt_1 : 1 < t := by linarith
    have hlog_pos : 0 < Real.log t := Real.log_pos ht_gt_1
    have hlog_ne : Real.log t ≠ 0 := ne_of_gt hlog_pos
    have hlog2_ne : Real.log t ^ 2 ≠ 0 := pow_ne_zero 2 hlog_ne
    have hprod_ne : t * (Real.log t) ^ 2 ≠ 0 := mul_ne_zero ht_ne hlog2_ne
    have hprod1_ne : t * Real.log t ≠ 0 := mul_ne_zero ht_ne hlog_ne
    show Real.log t / (t * (Real.log t) ^ 2) = 1 / (t * Real.log t)
    rw [pow_two]; field_simp
  -- ∫_y^x log(t)/(t · log²t) = ∫_y^x 1/(t·log t) = log log x - log log y.
  have hint_log_eq : ∫ t in y..x, Real.log t / (t * (Real.log t) ^ 2)
                  = ∫ t in y..x, 1 / (t * Real.log t) := by
    rw [intervalIntegral.integral_of_le hyx, intervalIntegral.integral_of_le hyx]
    exact MeasureTheory.integral_congr_ae hint_rec
  -- Combine: ∫_y^x mertensLog/(t log²t) = (log log x - log log y) + ∫_y^x δ(t)/(t·log²t).
  have hint_decomp : ∫ t in y..x, mertensLog t / (t * (Real.log t) ^ 2)
      = (Real.log (Real.log x) - Real.log (Real.log y))
      + ∫ t in y..x, (mertensLog t - Real.log t) / (t * (Real.log t) ^ 2) := by
    have heq : (fun t : ℝ => mertensLog t / (t * (Real.log t) ^ 2))
             = (fun t : ℝ => Real.log t / (t * (Real.log t) ^ 2)
                + (mertensLog t - Real.log t) / (t * (Real.log t) ^ 2)) := by
      funext t; ring
    rw [show (fun t : ℝ => mertensLog t / (t * (Real.log t) ^ 2)) =
        (fun t : ℝ => Real.log t / (t * (Real.log t) ^ 2)
         + (mertensLog t - Real.log t) / (t * (Real.log t) ^ 2)) from heq]
    rw [intervalIntegral.integral_add hint_log_yx hint_delta_yx]
    rw [hint_log_eq, hint_log_split]
  -- Now combine all: F(x) - F(y) = (δ(x)/log x - δ(y)/log y) + ∫_y^x δ(t)/(t·log²t).
  have hF_diff :
      (mertensSum x - Real.log (Real.log x)) -
        (mertensSum y - Real.log (Real.log y))
      = (mertensLog x - Real.log x) / Real.log x
        - (mertensLog y - Real.log y) / Real.log y
        + ∫ t in y..x, (mertensLog t - Real.log t) / (t * (Real.log t) ^ 2) := by
    have : mertensSum x - mertensSum y =
        ((mertensLog x - Real.log x) / Real.log x + 1)
        - ((mertensLog y - Real.log y) / Real.log y + 1)
        + ((Real.log (Real.log x) - Real.log (Real.log y))
           + ∫ t in y..x, (mertensLog t - Real.log t) / (t * (Real.log t) ^ 2)) := by
      rw [hsubpiv]
      have hxsim : mertensLog x / Real.log x =
          (mertensLog x - Real.log x) / Real.log x + 1 := by
        rw [hsimpx]; ring
      have hysim : mertensLog y / Real.log y =
          (mertensLog y - Real.log y) / Real.log y + 1 := by
        rw [hsimpy]; ring
      rw [hxsim, hysim, hint_decomp]
    linarith
  rw [hF_diff]
  -- Bound each of the three terms.
  -- Term 1: |(mertensLog x - log x)/log x| ≤ C/log x ≤ C/log y.
  have hδx_abs : |mertensLog x - Real.log x| ≤ C := hC_ge x hx
  have hδy_abs : |mertensLog y - Real.log y| ≤ C := hC_ge y hy
  have hT1 : |(mertensLog x - Real.log x) / Real.log x| ≤ C / Real.log y := by
    rw [abs_div]
    have hlogx_abs : |Real.log x| = Real.log x := abs_of_pos hlogx_pos
    rw [hlogx_abs]
    have hbnd1 : |mertensLog x - Real.log x| / Real.log x ≤ C / Real.log x := by
      apply div_le_div_of_nonneg_right hδx_abs (le_of_lt hlogx_pos)
    have hbnd2 : C / Real.log x ≤ C / Real.log y := by
      apply div_le_div_of_nonneg_left (le_of_lt hC_pos) hlogy_pos hlogx_ge_logy
    linarith
  have hT2 : |(mertensLog y - Real.log y) / Real.log y| ≤ C / Real.log y := by
    rw [abs_div]
    have hlogy_abs : |Real.log y| = Real.log y := abs_of_pos hlogy_pos
    rw [hlogy_abs]
    apply div_le_div_of_nonneg_right hδy_abs (le_of_lt hlogy_pos)
  -- Term 3: |∫_y^x δ(t)/(t·log²t)| ≤ C·(1/log y - 1/log x) ≤ C/log y.
  have hT3 : |∫ t in y..x, (mertensLog t - Real.log t) / (t * (Real.log t) ^ 2)|
              ≤ C / Real.log y := by
    -- Bound the integrand: |δ(t)/(t·log²t)| ≤ C/(t·log²t) for t ∈ [y, x].
    have habs_bound : ∀ t ∈ Set.uIoc y x,
        ‖(mertensLog t - Real.log t) / (t * (Real.log t) ^ 2)‖
          ≤ C * (1 / (t * (Real.log t) ^ 2)) := by
      intro t ht
      rw [Set.uIoc_of_le hyx] at ht
      have ht_ge_y : y < t := ht.1
      have ht_ge_2 : (2 : ℝ) ≤ t := by linarith
      have ht_pos : 0 < t := by linarith
      have ht_ne : t ≠ 0 := ne_of_gt ht_pos
      have ht_gt_1 : 1 < t := by linarith
      have hlog_pos : 0 < Real.log t := Real.log_pos ht_gt_1
      have hlog_ne : Real.log t ≠ 0 := ne_of_gt hlog_pos
      have hlog2_pos : 0 < (Real.log t) ^ 2 := by positivity
      have hprod_pos : 0 < t * (Real.log t) ^ 2 := by positivity
      have hδt : |mertensLog t - Real.log t| ≤ C := hC_ge t ht_ge_2
      rw [Real.norm_eq_abs, abs_div]
      have h1 : |t * (Real.log t) ^ 2| = t * (Real.log t) ^ 2 := abs_of_pos hprod_pos
      rw [h1]
      rw [mul_one_div]
      exact (div_le_div_iff_of_pos_right hprod_pos).mpr hδt
    have hbound_int : IntervalIntegrable (fun t => C * (1 / (t * (Real.log t) ^ 2)))
        MeasureTheory.volume y x := by
      refine ContinuousOn.intervalIntegrable ?_
      rw [Set.uIcc_of_le hyx]
      intro z hz
      have hz_ge_y : y ≤ z := hz.1
      have hz_ge_2 : (2 : ℝ) ≤ z := le_trans hy hz_ge_y
      have hz_pos : (0 : ℝ) < z := by linarith
      have hz_ne : z ≠ 0 := ne_of_gt hz_pos
      have hz_gt_1 : (1 : ℝ) < z := by linarith
      have hlog_pos : 0 < Real.log z := Real.log_pos hz_gt_1
      have hlog_ne : Real.log z ≠ 0 := ne_of_gt hlog_pos
      have hlog2_ne : Real.log z ^ 2 ≠ 0 := pow_ne_zero 2 hlog_ne
      have hprod_ne : z * (Real.log z) ^ 2 ≠ 0 := mul_ne_zero hz_ne hlog2_ne
      refine ContinuousAt.continuousWithinAt ?_
      fun_prop
    have hineq : |∫ t in y..x, (mertensLog t - Real.log t) / (t * (Real.log t) ^ 2)|
        ≤ ∫ t in y..x, C * (1 / (t * (Real.log t) ^ 2)) := by
      have habs_bound' : ∀ᵐ (t : ℝ) ∂MeasureTheory.volume,
          t ∈ Set.Ioc y x →
            ‖(mertensLog t - Real.log t) / (t * (Real.log t) ^ 2)‖
              ≤ C * (1 / (t * (Real.log t) ^ 2)) := by
        refine Filter.Eventually.of_forall ?_
        intro t ht
        have : t ∈ Set.uIoc y x := by
          rw [Set.uIoc_of_le hyx]; exact ht
        exact habs_bound t this
      have hint_abs : ‖∫ t in y..x, (mertensLog t - Real.log t) / (t * (Real.log t) ^ 2)‖
          ≤ ∫ t in y..x, C * (1 / (t * (Real.log t) ^ 2)) :=
        intervalIntegral.norm_integral_le_of_norm_le hyx habs_bound' hbound_int
      have h_norm_eq : ‖∫ t in y..x, (mertensLog t - Real.log t) / (t * (Real.log t) ^ 2)‖
          = |∫ t in y..x, (mertensLog t - Real.log t) / (t * (Real.log t) ^ 2)| := by
        rfl
      linarith [hint_abs, h_norm_eq.symm.le, h_norm_eq.le]
    -- Compute ∫_y^x C · (1/(t·log²t)) = C · (1/log y - 1/log x).
    have h_int_sq_yx : ∫ t in y..x, 1 / (t * (Real.log t) ^ 2)
        = 1 / Real.log y - 1 / Real.log x := by
      have h1 : ∫ t in (2 : ℝ)..x, 1 / (t * (Real.log t) ^ 2) =
          1 / Real.log 2 - 1 / Real.log x := integral_one_div_t_log_sq hx
      have h2 : ∫ t in (2 : ℝ)..y, 1 / (t * (Real.log t) ^ 2) =
          1 / Real.log 2 - 1 / Real.log y := integral_one_div_t_log_sq hy
      have hint_sq_2y : IntervalIntegrable (fun t => 1 / (t * (Real.log t) ^ 2))
          MeasureTheory.volume 2 y := by
        refine ContinuousOn.intervalIntegrable ?_
        rw [Set.uIcc_of_le hy]
        intro z hz
        have hz_ge_2 : (2 : ℝ) ≤ z := hz.1
        have hz_pos : (0 : ℝ) < z := by linarith
        have hz_ne : z ≠ 0 := ne_of_gt hz_pos
        have hz_gt_1 : (1 : ℝ) < z := by linarith
        have hlog_pos : 0 < Real.log z := Real.log_pos hz_gt_1
        have hlog_ne : Real.log z ≠ 0 := ne_of_gt hlog_pos
        have hlog2_ne : Real.log z ^ 2 ≠ 0 := pow_ne_zero 2 hlog_ne
        have hprod_ne : z * (Real.log z) ^ 2 ≠ 0 := mul_ne_zero hz_ne hlog2_ne
        refine ContinuousAt.continuousWithinAt ?_
        fun_prop
      have hint_sq_2x : IntervalIntegrable (fun t => 1 / (t * (Real.log t) ^ 2))
          MeasureTheory.volume 2 x := by
        refine ContinuousOn.intervalIntegrable ?_
        rw [Set.uIcc_of_le hx]
        intro z hz
        have hz_ge_2 : (2 : ℝ) ≤ z := hz.1
        have hz_pos : (0 : ℝ) < z := by linarith
        have hz_ne : z ≠ 0 := ne_of_gt hz_pos
        have hz_gt_1 : (1 : ℝ) < z := by linarith
        have hlog_pos : 0 < Real.log z := Real.log_pos hz_gt_1
        have hlog_ne : Real.log z ≠ 0 := ne_of_gt hlog_pos
        have hlog2_ne : Real.log z ^ 2 ≠ 0 := pow_ne_zero 2 hlog_ne
        have hprod_ne : z * (Real.log z) ^ 2 ≠ 0 := mul_ne_zero hz_ne hlog2_ne
        refine ContinuousAt.continuousWithinAt ?_
        fun_prop
      have hsplit := intervalIntegral.integral_interval_sub_left hint_sq_2x hint_sq_2y
      linarith
    have h_int_eq : ∫ t in y..x, C * (1 / (t * (Real.log t) ^ 2))
        = C * (1 / Real.log y - 1 / Real.log x) := by
      rw [intervalIntegral.integral_const_mul, h_int_sq_yx]
    have hineq' : |∫ t in y..x, (mertensLog t - Real.log t) / (t * (Real.log t) ^ 2)|
        ≤ C * (1 / Real.log y - 1 / Real.log x) := by
      rw [h_int_eq] at hineq
      exact hineq
    have hbnd_final : C * (1 / Real.log y - 1 / Real.log x) ≤ C / Real.log y := by
      have h1 : 0 ≤ 1 / Real.log x := by positivity
      have ha : C * (1 / Real.log y - 1 / Real.log x) ≤ C * (1 / Real.log y) := by
        apply mul_le_mul_of_nonneg_left _ (le_of_lt hC_pos)
        linarith
      have hb : C * (1 / Real.log y) = C / Real.log y := by ring
      linarith
    linarith
  -- Triangle inequality on the three terms.
  have htri : |(mertensLog x - Real.log x) / Real.log x
              - (mertensLog y - Real.log y) / Real.log y
              + ∫ t in y..x, (mertensLog t - Real.log t) / (t * (Real.log t) ^ 2)|
      ≤ |(mertensLog x - Real.log x) / Real.log x|
        + |(mertensLog y - Real.log y) / Real.log y|
        + |∫ t in y..x, (mertensLog t - Real.log t) / (t * (Real.log t) ^ 2)| := by
    set A := (mertensLog x - Real.log x) / Real.log x
    set B := (mertensLog y - Real.log y) / Real.log y
    set I := ∫ t in y..x, (mertensLog t - Real.log t) / (t * (Real.log t) ^ 2)
    have h1 : |A - B + I| ≤ |A - B| + |I| := abs_add_le _ _
    have h2 : |A - B| ≤ |A| + |B| := by
      have := abs_add_le A (-B)
      simpa [sub_eq_add_neg, abs_neg] using this
    linarith
  have hcombined : |(mertensLog x - Real.log x) / Real.log x
          - (mertensLog y - Real.log y) / Real.log y
          + ∫ t in y..x, (mertensLog t - Real.log t) / (t * (Real.log t) ^ 2)|
      ≤ 3 * C / Real.log y := by
    have hsum : C / Real.log y + C / Real.log y + C / Real.log y = 3 * C / Real.log y := by
      ring
    linarith
  -- Convert 3 * C / log y to 3 * (|C₀| + 1) / log y.
  have : 3 * C / Real.log y = 3 * (|C₀| + 1) / Real.log y := by rw [hCdef]
  linarith

/--
**Mertens' Third Theorem — existence of the Meissel–Mertens limit.**

The sequence `x ↦ mertensSum x - log log x` converges as `x → ∞`. The limit is
the **Meissel–Mertens constant** `M₃ ≈ 0.26149721…`.

The proof is the standard Cauchy-criterion argument: combining the M3 Abel
pivot `mertensSum_eq_mertensLog_div_log_add_integral`, the M1 bound, and the
closed-form antiderivatives `integral_one_div_t_log_t`, `integral_one_div_t_log_sq`,
one shows that `|F(x) - F(y)| ≤ K / log y` for `2 ≤ y ≤ x`, where
`F(x) := mertensSum x - log log x`. Completeness of `ℝ` then yields the limit. -/
theorem mertens_M3_exists :
    ∃ M : ℝ, ∀ ε > (0 : ℝ), ∃ X : ℝ, ∀ x ≥ X,
      |mertensSum x - Real.log (Real.log x) - M| < ε := by
  -- The oscillation bound: |F(x) - F(y)| ≤ K / log y for 2 ≤ y ≤ x.
  set K : ℝ := 3 * (|Classical.choose mertens_M1| + 1) with hKdef
  have hK_pos : 0 < K := by
    have : 0 < |Classical.choose mertens_M1| + 1 := by positivity
    positivity
  -- Construct a sequence u : ℕ → ℝ that is Cauchy.
  set u : ℕ → ℝ := fun n => mertensSum (n + 2 : ℝ) - Real.log (Real.log (n + 2 : ℝ)) with hudef
  -- u is Cauchy: for m ≤ n, |u n - u m| ≤ K / log (m + 2).
  have hu_cauchy : CauchySeq u := by
    rw [Metric.cauchySeq_iff]
    intro ε hε
    -- Choose N large so that K / log (N + 2) < ε.
    -- Since log is unbounded, K / log (n + 2) → 0.
    have hlog_tend : Filter.Tendsto (fun n : ℕ => K / Real.log ((n : ℝ) + 2))
        Filter.atTop (nhds 0) := by
      have h_cast : Filter.Tendsto (fun n : ℕ => ((n : ℝ) + 2))
          Filter.atTop Filter.atTop :=
        Filter.tendsto_atTop_add_const_right _ 2 tendsto_natCast_atTop_atTop
      have h1 : Filter.Tendsto (fun n : ℕ => Real.log ((n : ℝ) + 2))
          Filter.atTop Filter.atTop := Real.tendsto_log_atTop.comp h_cast
      exact Filter.Tendsto.div_atTop (f := fun _ : ℕ => K)
        (tendsto_const_nhds (x := K) (f := Filter.atTop)) h1
    obtain ⟨N, hN⟩ := (tendsto_atTop_nhds.mp hlog_tend) (Set.Ioo (-ε) ε)
      (by constructor <;> [linarith; linarith]) isOpen_Ioo
    refine ⟨N, ?_⟩
    -- A helper: for m ≤ n, the bound holds in the form |u n - u m|.
    have helper : ∀ {m n : ℕ}, N ≤ m → m ≤ n →
        |(mertensSum ((n : ℝ) + 2) - Real.log (Real.log ((n : ℝ) + 2)))
          - (mertensSum ((m : ℝ) + 2) - Real.log (Real.log ((m : ℝ) + 2)))| < ε := by
      intro m n hmN hmn
      have hy : (2 : ℝ) ≤ ((m : ℝ) + 2) := by
        have : (0 : ℝ) ≤ m := by exact_mod_cast (Nat.zero_le m)
        linarith
      have hyx : ((m : ℝ) + 2) ≤ ((n : ℝ) + 2) := by
        have : (m : ℝ) ≤ (n : ℝ) := by exact_mod_cast hmn
        linarith
      have hosc := mertensSum_sub_log_log_oscillation
        (y := ((m : ℝ) + 2)) (x := ((n : ℝ) + 2)) hy hyx
      have hbnd :
          |(mertensSum ((n : ℝ) + 2) - Real.log (Real.log ((n : ℝ) + 2)))
            - (mertensSum ((m : ℝ) + 2) - Real.log (Real.log ((m : ℝ) + 2)))| ≤
          K / Real.log ((m : ℝ) + 2) := by
        have := hosc
        rw [hKdef]
        exact this
      have hN_app : K / Real.log ((m : ℝ) + 2) ∈ Set.Ioo (-ε) ε := hN m hmN
      have hKlt : K / Real.log ((m : ℝ) + 2) < ε := hN_app.2
      linarith
    intro m hm n hn
    have hu_eq_m : u m = mertensSum ((m : ℝ) + 2) - Real.log (Real.log ((m : ℝ) + 2)) := rfl
    have hu_eq_n : u n = mertensSum ((n : ℝ) + 2) - Real.log (Real.log ((n : ℝ) + 2)) := rfl
    by_cases hmn : m ≤ n
    · rw [Real.dist_eq, hu_eq_m, hu_eq_n, abs_sub_comm]
      exact helper hm hmn
    · have hmn' : n ≤ m := Nat.le_of_lt (Nat.lt_of_not_le hmn)
      rw [Real.dist_eq, hu_eq_m, hu_eq_n]
      exact helper hn hmn'
  -- Extract the limit.
  obtain ⟨M, hM⟩ := cauchySeq_tendsto_of_complete hu_cauchy
  refine ⟨M, ?_⟩
  intro ε hε
  -- Use the same idea: pick N₁ for u → M (ε/2), N₂ for oscillation bound (ε/2).
  have hε2 : 0 < ε / 2 := by linarith
  -- u n is within ε/2 of M eventually.
  rw [Metric.tendsto_atTop] at hM
  obtain ⟨N₁, hN₁⟩ := hM (ε / 2) hε2
  -- K / log (N + 2) < ε / 2 for N ≥ some N₂.
  have hlog_tend : Filter.Tendsto (fun n : ℕ => K / Real.log (n + 2 : ℝ))
      Filter.atTop (nhds 0) := by
    have h1 : Filter.Tendsto (fun n : ℕ => Real.log ((n : ℝ) + 2))
        Filter.atTop Filter.atTop := by
      have h_cast : Filter.Tendsto (fun n : ℕ => ((n : ℝ) + 2))
          Filter.atTop Filter.atTop :=
        Filter.tendsto_atTop_add_const_right _ 2 tendsto_natCast_atTop_atTop
      exact Real.tendsto_log_atTop.comp h_cast
    have h2 := Filter.Tendsto.div_atTop (f := fun _ : ℕ => K)
      (tendsto_const_nhds (x := K) (f := Filter.atTop)) h1
    exact h2
  obtain ⟨N₂, hN₂⟩ := (tendsto_atTop_nhds.mp hlog_tend) (Set.Ioo (-(ε/2)) (ε/2))
    (by constructor <;> [linarith; linarith]) isOpen_Ioo
  set N := max N₁ N₂ with hNdef
  refine ⟨(N + 2 : ℝ), ?_⟩
  intro x hx
  -- For x ≥ N + 2: write x = (N + 2) + (x - (N + 2)). Use oscillation between N+2 and x.
  have hy : (2 : ℝ) ≤ (N + 2 : ℝ) := by
    have : (0 : ℝ) ≤ N := by exact_mod_cast (Nat.zero_le N)
    linarith
  have hyx : ((N : ℝ) + 2) ≤ x := hx
  have hosc := mertensSum_sub_log_log_oscillation
    (y := ((N : ℝ) + 2)) (x := x) hy hyx
  -- |F(x) - F(N+2)| ≤ K / log (N+2).
  -- |u N - M| < ε / 2.
  have hN_use : N₁ ≤ N := le_max_left _ _
  have hN_to_M : dist (u N) M < ε / 2 := hN₁ N hN_use
  have hN_log : K / Real.log (N + 2 : ℝ) < ε / 2 := by
    have hN_use : N₂ ≤ N := le_max_right _ _
    have hN₂_app : K / Real.log (N + 2 : ℝ) ∈ Set.Ioo (-(ε/2)) (ε/2) := hN₂ N hN_use
    exact hN₂_app.2
  -- Triangular.
  have hu_eq : u N = mertensSum (N + 2 : ℝ) - Real.log (Real.log (N + 2 : ℝ)) := rfl
  rw [Real.dist_eq] at hN_to_M
  have hbnd : |(mertensSum x - Real.log (Real.log x))
      - (mertensSum (N + 2 : ℝ) - Real.log (Real.log (N + 2 : ℝ)))| ≤
      K / Real.log (N + 2 : ℝ) := by
    have := hosc
    rw [hKdef]
    exact this
  have := abs_add_le
    ((mertensSum x - Real.log (Real.log x)) -
      (mertensSum (N + 2 : ℝ) - Real.log (Real.log (N + 2 : ℝ))))
    ((mertensSum (N + 2 : ℝ) - Real.log (Real.log (N + 2 : ℝ))) - M)
  have hsimp :
    ((mertensSum x - Real.log (Real.log x)) -
      (mertensSum (N + 2 : ℝ) - Real.log (Real.log (N + 2 : ℝ)))) +
      ((mertensSum (N + 2 : ℝ) - Real.log (Real.log (N + 2 : ℝ))) - M) =
      mertensSum x - Real.log (Real.log x) - M := by ring
  have hfinal :
      |mertensSum x - Real.log (Real.log x) - M| ≤
      |(mertensSum x - Real.log (Real.log x)) -
        (mertensSum (N + 2 : ℝ) - Real.log (Real.log (N + 2 : ℝ)))| +
      |(mertensSum (N + 2 : ℝ) - Real.log (Real.log (N + 2 : ℝ))) - M| := by
    rw [← hsimp]
    exact this
  rw [hu_eq] at hN_to_M
  linarith

/-- **`mertensM3`**: the Meissel–Mertens constant `M₃ ≈ 0.26149721…`, defined as
    the limit guaranteed by `mertens_M3_exists`. -/
noncomputable def mertensM3 : ℝ := Classical.choose mertens_M3_exists

/-- Defining property of `mertensM3` (immediate from `Classical.choose_spec`).
    This lemma is `sorry`-free *given* `mertens_M3_exists`. -/
theorem mertensM3_spec :
    ∀ ε > (0 : ℝ), ∃ X : ℝ, ∀ x ≥ X,
      |mertensSum x - Real.log (Real.log x) - mertensM3| < ε :=
  Classical.choose_spec mertens_M3_exists

/--
**Headline Mertens M3 with explicit `O(1/log x)` rate.**

$$\sum_{p \le x} \frac{1}{p} \;=\; \log\log x \,+\, M_3 \,+\, O\!\left(\frac{1}{\log x}\right).$$

TODO (estimated effort: ~1 session after `mertens_M3_exists`):
* This is a refinement of `mertens_M3_exists`: replace the qualitative
  Cauchy-criterion form with the quantitative big-O bound.
* The rate comes from the error term in `mertens_M1` (`O(1)` becomes
  `O(1/log x)` after the Abel integration by parts because of the
  `1/(t log² t)` weight).
* Mathlib idiom: use `Asymptotics.IsBigO` with filter `Filter.atTop`:
  `(fun x => mertensSum x - Real.log (Real.log x) - mertensM3)
     =O[atTop] (fun x => 1 / Real.log x)`.
-/
theorem mertens_M3 :
    ∃ C : ℝ, ∀ x : ℝ, 2 ≤ x →
      |mertensSum x - Real.log (Real.log x) - mertensM3| ≤ C / Real.log x := by
  -- The oscillation bound applied with y := x and arbitrary z ≥ x gives:
  --   |F(x) - F(z)| ≤ K / log x.
  -- Taking z → ∞, F(z) → mertensM3, hence |F(x) - mertensM3| ≤ K / log x.
  set K : ℝ := 3 * (|Classical.choose mertens_M1| + 1) with hKdef
  have hK_pos : 0 < K := by
    have : 0 < |Classical.choose mertens_M1| + 1 := by positivity
    positivity
  refine ⟨K, ?_⟩
  intro x hx
  -- For x ≥ 2: log x > 0.
  have hx_pos : 0 < x := by linarith
  have hx_gt_1 : 1 < x := by linarith
  have hlogx_pos : 0 < Real.log x := Real.log_pos hx_gt_1
  have hlogx_ne : Real.log x ≠ 0 := ne_of_gt hlogx_pos
  -- Get the limit spec.
  have hspec := mertensM3_spec
  -- For any ε > 0, eventually |F(z) - mertensM3| < ε.
  -- Combining with oscillation: |F(x) - F(z)| ≤ K/log x, hence
  -- |F(x) - mertensM3| ≤ K/log x + ε. Letting ε → 0, conclude.
  -- We prove the bound by contradiction.
  by_contra hcontra
  have h : K / Real.log x < |mertensSum x - Real.log (Real.log x) - mertensM3| :=
    not_le.mp hcontra
  -- h : K / log x < |F(x) - mertensM3|.
  set δ : ℝ := |mertensSum x - Real.log (Real.log x) - mertensM3| - K / Real.log x with hδdef
  have hδ_pos : 0 < δ := by rw [hδdef]; linarith
  obtain ⟨X, hX⟩ := hspec (δ / 2) (by linarith)
  -- Pick z ≥ max x X.
  set z : ℝ := max x X + 1 with hzdef
  have hz_ge_x : x ≤ z := by
    rw [hzdef]; linarith [le_max_left x X]
  have hz_ge_X : X ≤ z := by
    rw [hzdef]; linarith [le_max_right x X]
  have hosc := mertensSum_sub_log_log_oscillation
    (y := x) (x := z) hx hz_ge_x
  have hbnd : |(mertensSum z - Real.log (Real.log z))
              - (mertensSum x - Real.log (Real.log x))| ≤
              K / Real.log x := by
    have := hosc
    rw [hKdef]
    exact this
  have hX_app := hX z hz_ge_X
  -- |F(z) - mertensM3| < δ/2.
  -- Therefore: |F(x) - mertensM3| ≤ |F(x) - F(z)| + |F(z) - mertensM3|
  --                              ≤ K/log x + δ/2.
  -- But |F(x) - mertensM3| = K/log x + δ.
  -- So δ ≤ δ/2, contradiction (since δ > 0).
  have htri := abs_add_le
    ((mertensSum x - Real.log (Real.log x)) - (mertensSum z - Real.log (Real.log z)))
    ((mertensSum z - Real.log (Real.log z)) - mertensM3)
  have hsimp :
      ((mertensSum x - Real.log (Real.log x)) - (mertensSum z - Real.log (Real.log z))) +
      ((mertensSum z - Real.log (Real.log z)) - mertensM3) =
      mertensSum x - Real.log (Real.log x) - mertensM3 := by ring
  have hT1 : |((mertensSum x - Real.log (Real.log x)) - (mertensSum z - Real.log (Real.log z)))| ≤
             K / Real.log x := by
    rw [abs_sub_comm]; exact hbnd
  have hT2 : |((mertensSum z - Real.log (Real.log z)) - mertensM3)| < δ / 2 := hX_app
  have : |mertensSum x - Real.log (Real.log x) - mertensM3| ≤
         K / Real.log x + δ / 2 := by
    have hsim2 :
      |mertensSum x - Real.log (Real.log x) - mertensM3| ≤
      |((mertensSum x - Real.log (Real.log x)) - (mertensSum z - Real.log (Real.log z)))|
        + |((mertensSum z - Real.log (Real.log z)) - mertensM3)| := by
      rw [← hsimp]; exact htri
    linarith
  -- Contradiction.
  linarith [hδ_pos]

/-! ### PT-specific corollary: compactness for T5 (M3 article) -/

/--
**The PT-needed compactness corollary.**

In the M3 article (Remark `rem:mertens_role`), the classical Mertens theorem is
used *exclusively* to deduce that the empirical sequence `α_k → 1/2` and the
transition matrices `T(k)` remain in a compact set. The Lean-level statement
of this dependency is that `mertensSum x - log log x` is **bounded** as
`x → ∞`. This is strictly weaker than `mertens_M3` (no rate needed) and is
all that the PT pipeline consumes.

TODO (estimated effort: trivial after `mertens_M3`):
* Direct consequence of `mertens_M3`: `|f(x)| ≤ |M₃| + C / log x` on `x ≥ 2`,
  and the RHS is bounded for `x ≥ e` (where `log x ≥ 1`).
* The cheap proof is `|f(x)| ≤ |M₃| + C` for `x ≥ e` and a finite case-split
  on `x ∈ [2, e]`.
-/
theorem mertensSum_sub_log_log_bounded :
    ∃ K : ℝ, ∀ x : ℝ, 2 ≤ x →
      |mertensSum x - Real.log (Real.log x)| ≤ K := by
  -- Strategy: from `mertens_M3_exists`, the function is bounded near ∞.
  -- For x in the finite interval [2, X], bound `mertensSum x` by `mertensSum X'`
  -- (monotone in `⌊x⌋₊`) and `|log log x|` by `max |log log 2| |log log X'|`.
  obtain ⟨M, hM⟩ := mertens_M3_exists
  obtain ⟨X, hX⟩ := hM 1 one_pos
  -- Define a finite cutoff X' ≥ max X 2.
  set X' : ℝ := max X 2 with hX'def
  have hX'_ge_X : X ≤ X' := le_max_left _ _
  have hX'_ge_2 : (2 : ℝ) ≤ X' := le_max_right _ _
  -- Bound on the tail x ≥ X': |f(x)| ≤ |M| + 1.
  have htail : ∀ x ≥ X', |mertensSum x - Real.log (Real.log x)| ≤ |M| + 1 := by
    intro x hx
    have hx_ge_X : X ≤ x := le_trans hX'_ge_X hx
    have hbnd := hX x hx_ge_X
    have htriangle : |mertensSum x - Real.log (Real.log x)|
        ≤ |mertensSum x - Real.log (Real.log x) - M| + |M| := by
      have := abs_add_le (mertensSum x - Real.log (Real.log x) - M) M
      simpa [sub_add_cancel] using this
    linarith
  -- Bound on the head x ∈ [2, X']: mertensSum is monotone in ⌊x⌋₊,
  -- bounded by the value at the natural number ⌊X'⌋₊ + 1.
  -- And |log log x| ≤ max |log log 2| |log log X'| because log∘log is monotone
  -- on [2, ∞) (since log 2 > 0 means log is monotone there).
  set N : ℕ := ⌊X'⌋₊ + 1 with hNdef
  -- Tail term bound:  ∑_{p ∈ primesBelow N} 1/p, a fixed finite number.
  set Sbound : ℝ := ∑ p ∈ Nat.primesBelow N, (1 : ℝ) / p with hSbound
  -- |log log x| bound on [2, X']: log log is monotone on [2, ∞) (where log ≥ log 2 > 0).
  have hlog2_pos : 0 < Real.log 2 := Real.log_pos one_lt_two
  set Lbound : ℝ := max |Real.log (Real.log 2)| |Real.log (Real.log X')| with hLbound
  refine ⟨max (Sbound + Lbound) (|M| + 1), ?_⟩
  intro x hx
  by_cases hxX' : X' ≤ x
  · -- Tail case: apply htail.
    exact le_trans (htail x hxX') (le_max_right _ _)
  · -- Head case: x ∈ [2, X').
    have hx_lt_X' : x < X' := lt_of_not_ge hxX'
    -- Bound mertensSum x ≤ Sbound: ⌊x⌋₊ + 1 ≤ N.
    have hfloor_le : ⌊x⌋₊ + 1 ≤ N := by
      have h1 : ⌊x⌋₊ ≤ ⌊X'⌋₊ := Nat.floor_le_floor (le_of_lt hx_lt_X')
      omega
    have hSum_le : mertensSum x ≤ Sbound := by
      unfold mertensSum
      apply Finset.sum_le_sum_of_subset_of_nonneg
      · -- primesBelow is monotone in n.
        intro p hp
        rw [Nat.primesBelow_eq_filter_range, Finset.mem_filter] at hp ⊢
        refine ⟨?_, hp.2⟩
        exact Finset.mem_range.mpr (lt_of_lt_of_le (Finset.mem_range.mp hp.1) hfloor_le)
      · intro p hp _
        have hp_prime : Nat.Prime p := by
          rw [Nat.primesBelow_eq_filter_range, Finset.mem_filter] at hp
          exact hp.2
        have hp_pos : 0 < (p : ℝ) := by exact_mod_cast hp_prime.pos
        positivity
    have hSum_nonneg : 0 ≤ mertensSum x := mertensSum_nonneg x
    -- Bound |log log x|: log log is monotone where log ≥ 1; but here log might be < 1.
    -- We just bound by Lbound = max |log log 2| |log log X'|.
    have hlog_x_pos : 0 < Real.log x := Real.log_pos (lt_of_lt_of_le one_lt_two hx)
    have hlog_2_le : Real.log 2 ≤ Real.log x := Real.log_le_log (by norm_num) hx
    have hlog_x_le : Real.log x ≤ Real.log X' :=
      Real.log_le_log (by linarith) (le_of_lt hx_lt_X')
    -- So log log 2 ≤ log log x ≤ log log X'.
    have hloglog_x_lb : Real.log (Real.log 2) ≤ Real.log (Real.log x) :=
      Real.log_le_log hlog2_pos hlog_2_le
    have hloglog_x_ub : Real.log (Real.log x) ≤ Real.log (Real.log X') :=
      Real.log_le_log hlog_x_pos hlog_x_le
    have hloglog_abs : |Real.log (Real.log x)| ≤ Lbound := by
      rw [abs_le]
      constructor
      · have h1 : -|Real.log (Real.log 2)| ≤ Real.log (Real.log 2) := neg_abs_le _
        have h2 : -Lbound ≤ -|Real.log (Real.log 2)| := by
          rw [neg_le_neg_iff]; exact le_max_left _ _
        linarith
      · have h1 : Real.log (Real.log X') ≤ |Real.log (Real.log X')| := le_abs_self _
        have h2 : |Real.log (Real.log X')| ≤ Lbound := le_max_right _ _
        linarith
    -- Combine: |mertensSum x - log log x| ≤ |mertensSum x| + |log log x| ≤ Sbound + Lbound.
    have : |mertensSum x - Real.log (Real.log x)|
        ≤ |mertensSum x| + |Real.log (Real.log x)| := abs_sub _ _
    have habs_sum : |mertensSum x| ≤ Sbound := by
      rw [abs_of_nonneg hSum_nonneg]; exact hSum_le
    have : |mertensSum x - Real.log (Real.log x)| ≤ Sbound + Lbound := by
      linarith
    exact le_trans this (le_max_left _ _)

/--
**Alias used by the PT M3 article.** `T5` in the M3 article refers to the
self-consistency theorem `μ* = 15`; its analytic prerequisite is precisely
`mertensSum_sub_log_log_bounded`. We expose it here under the PT name so that
downstream PT modules can `import PT.NumberTheory.T5Mertens` and `apply
pt_T5_mertens_compactness` without thinking about Mertens.

TODO (estimated effort: trivial — direct alias).
-/
theorem pt_T5_mertens_compactness :
    ∃ K : ℝ, ∀ x : ℝ, 2 ≤ x →
      |mertensSum x - Real.log (Real.log x)| ≤ K :=
  mertensSum_sub_log_log_bounded

end PT.NumberTheory
