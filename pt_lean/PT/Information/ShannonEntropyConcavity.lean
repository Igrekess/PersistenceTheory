/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Information.GFTIdentity
import PT.Information.EntropyBoundsTight
import PT.Information.EntropyOfBinaryDistribution
import Mathlib.Analysis.SpecialFunctions.Log.NegMulLog
import Mathlib.Analysis.Convex.Jensen
import Mathlib.Tactic

/-!
# Concavity of Shannon entropy

The Shannon entropy `H(p) = Σ_r negMulLog (p_r)` is **concave** in the
probability vector `p`. This file delivers a self-contained, Mathlib-
faithful chain of concavity results in the PT namespace:

1. **Pointwise concavity of `negMulLog`** (re-export of
   `Real.concaveOn_negMulLog` and its strict variant).
2. **Two-distribution concavity of `shannonH`**: for non-negative
   distributions `p, q : α → ℝ` and weights `a, b ≥ 0` with `a + b = 1`,
   `a · shannonH s p + b · shannonH s q ≤ shannonH s (a • p + b • q)`.
   Proved by termwise concavity of `negMulLog` and summation.
3. **Bekenstein via Jensen** (corner reformulation): on a finite set of
   cardinality `m > 0`, for any probability vector `p`,
   `shannonH s p ≤ log m`. This restates `EntropyBoundsTight` via the
   concavity-style argument and the GFT identity.
4. **Binary concavity**: `binEntropy (a · p + b · q) ≥ a · binEntropy p +
   b · binEntropy q` for `a, b ≥ 0` with `a + b = 1` and `p, q ∈ [0, 1]`.
5. **Mid-point specialisation**: `binEntropy ((1/3 + 1/2)/2) =
   binEntropy (5/12) ≥ (1/2) · (binEntropy (1/3) + binEntropy (1/2))`.

The general Jensen sum `H(Σ_i λ_i p_i) ≥ Σ_i λ_i H(p_i)` over arbitrary
convex combinations is a downstream corollary (left for a dedicated
module); here we cover the binary mixture, which is the bedrock case
used throughout the PT corpus (T6 averaging, Kähler-Fisher mixtures,
binary entropy at PT-canonical points).

## References

Mathlib: `Real.concaveOn_negMulLog`, `Real.strictConcaveOn_negMulLog`.
Cover-Thomas §2.7 (concavity of entropy); Amari-Nagaoka §3.5
(Fisher-Shannon hessian). Monograph Ch04 §"Bornes tendues".
-/

namespace PT.Information

open Real Finset

/-! ### (1) Pointwise concavity of `negMulLog` (re-export) -/

/-- **Concavity of `negMulLog` (re-export).** The map `x ↦ -x log x` is
    concave on `[0, ∞)`. -/
theorem concaveOn_negMulLog_PT :
    ConcaveOn ℝ (Set.Ici (0 : ℝ)) Real.negMulLog :=
  Real.concaveOn_negMulLog

/-- **Strict concavity of `negMulLog` (re-export).** The map
    `x ↦ -x log x` is strictly concave on `[0, ∞)`. -/
theorem strictConcaveOn_negMulLog_PT :
    StrictConcaveOn ℝ (Set.Ici (0 : ℝ)) Real.negMulLog :=
  Real.strictConcaveOn_negMulLog

/-- **Pointwise Jensen at two points for `negMulLog`.**
    For `x, y ≥ 0` and weights `a, b ≥ 0` with `a + b = 1`,
    `a · negMulLog x + b · negMulLog y ≤ negMulLog (a · x + b · y)`. -/
lemma negMulLog_two_point_jensen
    {x y a b : ℝ} (hx : 0 ≤ x) (hy : 0 ≤ y)
    (ha : 0 ≤ a) (hb : 0 ≤ b) (hab : a + b = 1) :
    a * Real.negMulLog x + b * Real.negMulLog y
      ≤ Real.negMulLog (a * x + b * y) := by
  have hx_mem : x ∈ Set.Ici (0 : ℝ) := hx
  have hy_mem : y ∈ Set.Ici (0 : ℝ) := hy
  have key := Real.concaveOn_negMulLog.2 hx_mem hy_mem ha hb hab
  simpa [smul_eq_mul] using key

/-! ### (2) Two-distribution concavity of `shannonH` -/

/-- **Two-distribution concavity of `shannonH`.**

    For non-negative distributions `p, q : α → ℝ` on a finite set `s`
    and weights `a, b ≥ 0` with `a + b = 1`,
    `a · shannonH s p + b · shannonH s q ≤ shannonH s (a • p + b • q)`.

    Proof: pointwise Jensen on `negMulLog` (Mathlib's
    `Real.concaveOn_negMulLog`), then summed over `s`. -/
theorem shannonH_concave_two
    {α : Type*} (s : Finset α) (p q : α → ℝ)
    (hp : ∀ r ∈ s, 0 ≤ p r) (hq : ∀ r ∈ s, 0 ≤ q r)
    {a b : ℝ} (ha : 0 ≤ a) (hb : 0 ≤ b) (hab : a + b = 1) :
    a * shannonH s p + b * shannonH s q
      ≤ shannonH s (fun r => a * p r + b * q r) := by
  unfold shannonH
  -- Reduce LHS to a single sum.
  rw [Finset.mul_sum, Finset.mul_sum, ← Finset.sum_add_distrib]
  -- Termwise Jensen.
  apply Finset.sum_le_sum
  intro r hr
  exact negMulLog_two_point_jensen (hp r hr) (hq r hr) ha hb hab

/-- **Midpoint specialisation of `shannonH_concave_two`.**
    `(1/2) · (shannonH s p + shannonH s q) ≤ shannonH s ((1/2) • (p + q))`. -/
theorem shannonH_concave_midpoint
    {α : Type*} (s : Finset α) (p q : α → ℝ)
    (hp : ∀ r ∈ s, 0 ≤ p r) (hq : ∀ r ∈ s, 0 ≤ q r) :
    (1 : ℝ) / 2 * shannonH s p + (1 : ℝ) / 2 * shannonH s q
      ≤ shannonH s (fun r => (1 : ℝ) / 2 * p r + (1 : ℝ) / 2 * q r) := by
  apply shannonH_concave_two s p q hp hq
  · norm_num
  · norm_num
  · norm_num

/-! ### (3) Bekenstein via Jensen / GFT (corner reformulation) -/

/-- **Bekenstein bound (Jensen corner).** For a finite set `s` of
    cardinality `m > 0` and a probability vector `p : α → ℝ`
    (non-negative, sums to 1), `shannonH s p ≤ log m`.

    Reformulation via the GFT identity: `klToUniform s m p ≥ 0` follows
    from Gibbs (here packaged by a single concavity step). We give the
    cleanest derivation: combine `GFT_identity` with the fact that the
    uniform distribution `U_m` maximises entropy among `s`-distributions,
    where the maximum is attained at `H(U_m) = log m`
    (`shannonH_uniform_eq_log`). The Jensen step is encapsulated in
    the inequality
    `negMulLog ((Σ p_r) / m) ≥ (1/m) · Σ negMulLog (p_r)`,
    which at `Σ p_r = 1` gives `negMulLog (1/m) ≥ (1/m) · H(p)`,
    i.e. `(log m)/m ≥ H(p)/m`, hence `H(p) ≤ log m`. -/
theorem shannonH_le_log_card_of_jensen
    {α : Type*} (s : Finset α) (m : ℝ) (hm : 0 < m)
    (hcard : (s.card : ℝ) = m)
    (p : α → ℝ) (hp_nonneg : ∀ r ∈ s, 0 ≤ p r)
    (hp_sum : ∑ r ∈ s, p r = 1) :
    shannonH s p ≤ Real.log m := by
  -- Use the GFT identity: `klToUniform + shannonH = log m`.
  have hid := GFT_identity s m hm p hp_nonneg hp_sum
  -- Suffices: `0 ≤ klToUniform s m p`.
  -- We prove the equivalent form: `shannonH s p ≤ shannonH s U_m`,
  -- where `U_m = fun _ => 1/m`. Combine with `shannonH_uniform_eq_log`.
  -- The KL non-negativity is established by Jensen applied to negMulLog,
  -- using cardinality weighting `1/m` over the support `s`.
  -- We bypass writing it explicitly and use the uniform-maximum form.
  have hm_ne : m ≠ 0 := ne_of_gt hm
  -- The uniform distribution `U_m` (constant `1/m`) on `s`.
  -- By concavity of `negMulLog` applied via the finite Jensen sum
  -- with equal weights `1/m`, Σ_r (1/m) negMulLog (p_r) ≤ negMulLog((1/m) Σ p_r).
  -- Multiplying by m and using Σ p_r = 1 gives `H(p) ≤ m · negMulLog(1/m) = log m`.
  -- We instantiate Mathlib's `ConcaveOn.le_inner_le_iff`-style on Finset using
  -- `Real.concaveOn_negMulLog.smul_le_sum`-pattern.
  -- Concretely:
  have key : ∑ r ∈ s, (1 / m : ℝ) * Real.negMulLog (p r)
              ≤ Real.negMulLog (∑ r ∈ s, (1 / m : ℝ) * p r) := by
    -- ConcaveOn.le_sum-style: but Mathlib bundles this as `ConcaveOn.inner_le_iff`
    -- or `ConcaveOn.smul_le_sum`. We use `ConcaveOn.le_sum` (or its dual).
    -- The exact lemma: `ConcaveOn.le_sum`: ∑ w_r * f(x_r) ≤ f(∑ w_r * x_r)
    -- under non-negative weights summing to 1.
    have hw_nn : ∀ r ∈ s, (0 : ℝ) ≤ (1 / m : ℝ) := by
      intro r _; positivity
    have hw_sum : ∑ _r ∈ s, (1 / m : ℝ) = 1 := by
      rw [Finset.sum_const, nsmul_eq_mul, hcard]
      field_simp
    have hx_mem : ∀ r ∈ s, p r ∈ Set.Ici (0 : ℝ) := fun r hr => hp_nonneg r hr
    have jensen := Real.concaveOn_negMulLog.le_map_sum
      (t := s) (w := fun _ => (1 / m : ℝ)) (p := p) hw_nn hw_sum hx_mem
    simpa [smul_eq_mul] using jensen
  -- Now: ∑ p_r / m = 1/m, so the RHS = negMulLog (1/m) = (log m)/m.
  have hRHS_arg : (∑ r ∈ s, (1 / m : ℝ) * p r) = 1 / m := by
    rw [← Finset.mul_sum, hp_sum]; ring
  rw [hRHS_arg] at key
  -- LHS: (1/m) · shannonH s p.
  have hLHS : ∑ r ∈ s, (1 / m : ℝ) * Real.negMulLog (p r)
              = (1 / m) * shannonH s p := by
    unfold shannonH
    rw [Finset.mul_sum]
  rw [hLHS] at key
  -- Compute negMulLog (1/m) = (log m)/m.
  have hneg : Real.negMulLog ((1 : ℝ) / m) = Real.log m / m := by
    unfold Real.negMulLog
    have hlog : Real.log (1 / m) = -Real.log m := by
      rw [one_div]; exact Real.log_inv m
    rw [hlog]; field_simp
  rw [hneg] at key
  -- Conclude: H ≤ log m via multiplication by m > 0.
  have hm_pos : (0 : ℝ) < m := hm
  have hinv_pos : (0 : ℝ) < 1 / m := by positivity
  -- key: (1/m) * H ≤ (log m) / m. Multiply both sides by m.
  have : shannonH s p ≤ Real.log m := by
    have hmul := mul_le_mul_of_nonneg_left key (le_of_lt hm_pos)
    -- m * ((1/m) * H) = H and m * ((log m)/m) = log m.
    have hL : m * ((1 / m) * shannonH s p) = shannonH s p := by
      field_simp
    have hR : m * (Real.log m / m) = Real.log m := by
      field_simp
    rw [hL, hR] at hmul
    exact hmul
  exact this

/-- **Bekenstein bound (Jensen corner) for `Finset.univ` on `Fin n`.** -/
theorem shannonH_le_log_card_fin
    {n : ℕ} (hn : 0 < n)
    (p : Fin n → ℝ) (hp_nonneg : ∀ r ∈ (Finset.univ : Finset (Fin n)), 0 ≤ p r)
    (hp_sum : ∑ r, p r = 1) :
    shannonH (Finset.univ : Finset (Fin n)) p ≤ Real.log (n : ℝ) := by
  apply shannonH_le_log_card_of_jensen (s := (Finset.univ : Finset (Fin n)))
    (m := (n : ℝ))
  · exact_mod_cast hn
  · simp
  · exact hp_nonneg
  · exact hp_sum

/-! ### (4) Binary concavity -/

/-- **Binary concavity (two-point Jensen for `binEntropy`).**

    For `p, q ∈ [0, 1]` and weights `a, b ≥ 0` with `a + b = 1`,
    `a · binEntropy p + b · binEntropy q ≤ binEntropy (a · p + b · q)`.

    Direct consequence of `negMulLog`'s concavity applied independently
    to the two coordinates `(p, q)` and `(1-p, 1-q)`. -/
theorem binEntropy_concave_two
    {p q : ℝ} (hp0 : 0 ≤ p) (hp1 : p ≤ 1) (hq0 : 0 ≤ q) (hq1 : q ≤ 1)
    {a b : ℝ} (ha : 0 ≤ a) (hb : 0 ≤ b) (hab : a + b = 1) :
    a * binEntropy p + b * binEntropy q ≤ binEntropy (a * p + b * q) := by
  -- Expand binEntropy via `negMulLog`.
  rw [binEntropy_eq, binEntropy_eq, binEntropy_eq]
  have h1p_nn : (0 : ℝ) ≤ 1 - p := by linarith
  have h1q_nn : (0 : ℝ) ≤ 1 - q := by linarith
  -- Concavity for the `p`-coordinate.
  have hP := negMulLog_two_point_jensen hp0 hq0 ha hb hab
  -- Concavity for the `(1-p)`-coordinate.
  have hQ := negMulLog_two_point_jensen h1p_nn h1q_nn ha hb hab
  -- The mixture argument in the `(1-p)` slot.
  have hmix : a * (1 - p) + b * (1 - q) = 1 - (a * p + b * q) := by
    have : a + b = 1 := hab
    linarith
  rw [hmix] at hQ
  -- Combine the two pointwise Jensens.
  linarith

/-- **Binary midpoint concavity.** For `p, q ∈ [0, 1]`,
    `(1/2) · (binEntropy p + binEntropy q) ≤ binEntropy ((p + q)/2)`. -/
theorem binEntropy_concave_midpoint
    {p q : ℝ} (hp0 : 0 ≤ p) (hp1 : p ≤ 1) (hq0 : 0 ≤ q) (hq1 : q ≤ 1) :
    (1 : ℝ) / 2 * binEntropy p + (1 : ℝ) / 2 * binEntropy q
      ≤ binEntropy ((p + q) / 2) := by
  have h := binEntropy_concave_two hp0 hp1 hq0 hq1
    (show (0 : ℝ) ≤ 1 / 2 by norm_num)
    (show (0 : ℝ) ≤ 1 / 2 by norm_num)
    (show (1 : ℝ) / 2 + 1 / 2 = 1 by norm_num)
  have harg : (1 : ℝ) / 2 * p + (1 : ℝ) / 2 * q = (p + q) / 2 := by ring
  rw [harg] at h
  exact h

/-! ### (5) Specialised mid-point: `1/3` and `1/2` -/

/-- **Mid-point specialisation.** The mixture of `(1/3, 1/2)` lies at
    `5/12`, and binary concavity gives
    `binEntropy (5/12) ≥ (1/2) · (binEntropy (1/3) + binEntropy (1/2))`. -/
theorem binEntropy_midpoint_third_half :
    (1 : ℝ) / 2 * binEntropy ((1 : ℝ) / 3) + (1 : ℝ) / 2 * binEntropy ((1 : ℝ) / 2)
      ≤ binEntropy ((5 : ℝ) / 12) := by
  have h := binEntropy_concave_midpoint
    (p := (1 : ℝ) / 3) (q := (1 : ℝ) / 2)
    (by norm_num) (by norm_num) (by norm_num) (by norm_num)
  have harg : ((1 : ℝ) / 3 + 1 / 2) / 2 = 5 / 12 := by norm_num
  rw [harg] at h
  exact h

/-- **Mid-point arithmetic.** `(1/3 + 1/2)/2 = 5/12`. -/
lemma midpoint_third_half : ((1 : ℝ) / 3 + 1 / 2) / 2 = 5 / 12 := by norm_num

/-! ### Headline -/

/-- **Headline (concavity of Shannon / binary entropy).**

    Bundling the principal concavity facts proved in this file:

    * `negMulLog` is concave on `[0, ∞)`.
    * `shannonH` is concave in the two-distribution sense.
    * Bekenstein bound `shannonH ≤ log m` follows from Jensen on `negMulLog`.
    * `binEntropy` is concave in the two-distribution sense.
    * Specialisation at `(1/3, 1/2)`: the average is bounded by
      `binEntropy (5/12)`. -/
theorem shannon_concavity_headline :
    -- (1) pointwise concavity of `negMulLog`
    ConcaveOn ℝ (Set.Ici (0 : ℝ)) Real.negMulLog
    -- (2) two-distribution concavity of `shannonH`
    ∧ (∀ (s : Finset (Fin 2)) (p q : Fin 2 → ℝ),
          (∀ r ∈ s, 0 ≤ p r) → (∀ r ∈ s, 0 ≤ q r) →
          ∀ (a b : ℝ), 0 ≤ a → 0 ≤ b → a + b = 1 →
          a * shannonH s p + b * shannonH s q
            ≤ shannonH s (fun r => a * p r + b * q r))
    -- (3) binary midpoint concavity
    ∧ (∀ {p q : ℝ}, 0 ≤ p → p ≤ 1 → 0 ≤ q → q ≤ 1 →
          (1 : ℝ) / 2 * binEntropy p + (1 : ℝ) / 2 * binEntropy q
            ≤ binEntropy ((p + q) / 2))
    -- (4) specialised mid-point inequality
    ∧ ((1 : ℝ) / 2 * binEntropy ((1 : ℝ) / 3)
         + (1 : ℝ) / 2 * binEntropy ((1 : ℝ) / 2)
        ≤ binEntropy ((5 : ℝ) / 12)) :=
  ⟨concaveOn_negMulLog_PT,
   fun s p q hp hq _a _b ha hb hab => shannonH_concave_two s p q hp hq ha hb hab,
   fun hp0 hp1 hq0 hq1 => binEntropy_concave_midpoint hp0 hp1 hq0 hq1,
   binEntropy_midpoint_third_half⟩

end PT.Information
