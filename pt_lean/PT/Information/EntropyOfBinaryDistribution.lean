/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Information.GFTIdentity
import PT.Information.BekensteinBound
import PT.Information.GFTSpecialMValues
import PT.Information.EntropyBoundsTight
import Mathlib.Analysis.SpecialFunctions.Log.NegMulLog
import Mathlib.Tactic

/-!
# Entropy of binary distributions at PT-canonical points

For a real `p`, the **binary distribution** `binDist p : Fin 2 → ℝ` is
defined by

* `binDist p 0 = p`,
* `binDist p 1 = 1 - p`.

Its Shannon entropy (in natural log) is
`H_bin(p) := -p log p - (1-p) log (1-p)`, equivalently
`H_bin(p) = negMulLog p + negMulLog (1-p)`.

We give:

* the definition `binDist` and the elementary lemma that it sums to `1`;
* the definition `binEntropy p := shannonH univ (binDist p)`;
* the equilibrium identity `binEntropy (1/2) = log 2`;
* the symmetry `binEntropy p = binEntropy (1 - p)`;
* the endpoint values `binEntropy 0 = 0` and `binEntropy 1 = 0`
  (under the Mathlib convention `Real.log 0 = 0`);
* the Bekenstein bound `binEntropy p ≤ log 2` on `p ∈ [0, 1]`;
* the explicit value at the PT-canonical point `p = 1/3`:
  `binEntropy (1/3) = log 3 - (2/3) log 2`;
* the explicit (algebraic) value at the PT-canonical point `p = 13/15`:
  `binEntropy (13/15) = -(13/15) log (13/15) - (2/15) log (2/15)`.

The three values `p ∈ {1/2, 1/3, 13/15}` correspond, respectively, to
the PT symmetry parameter `s = 1/2` (`PT.Stochastic.SHalf`), to the
flat T1 prior on the `Fin 3`-quotient (uniform on the active primes),
and to the equilibrium weight at `(p, 1-p) = (13/15, 2/15)` of the
T1 ratio on the active-prime sector.

## References

Monograph Chapter 4 §"Bornes tendues" (Bekenstein), and the auxiliary
files `GFTIdentity`, `GFTSpecialisations`, `EntropyBoundsTight`.
-/

namespace PT.Information

open Real Finset

/-! ### Definition of the binary distribution -/

/-- The binary distribution on `Fin 2` with parameter `p`. -/
noncomputable def binDist (p : ℝ) : Fin 2 → ℝ
  | 0 => p
  | 1 => 1 - p

@[simp] lemma binDist_zero (p : ℝ) : binDist p 0 = p := rfl
@[simp] lemma binDist_one  (p : ℝ) : binDist p 1 = 1 - p := rfl

/-- `binDist p` sums to `1` over `Fin 2` (no hypothesis needed). -/
lemma binDist_sum_eq_one (p : ℝ) :
    ∑ i, binDist p i = 1 := by
  rw [Fin.sum_univ_two]
  simp [binDist]

/-- `binDist p` is non-negative iff `0 ≤ p ≤ 1`. -/
lemma binDist_nonneg {p : ℝ} (hp0 : 0 ≤ p) (hp1 : p ≤ 1) :
    ∀ i ∈ (Finset.univ : Finset (Fin 2)), 0 ≤ binDist p i := by
  intro i _
  fin_cases i
  · simpa using hp0
  · simp; linarith

/-- `binDist p` is bounded above by `1` on `[0, 1]`. -/
lemma binDist_le_one {p : ℝ} (hp0 : 0 ≤ p) (hp1 : p ≤ 1) :
    ∀ i ∈ (Finset.univ : Finset (Fin 2)), binDist p i ≤ 1 := by
  intro i _
  fin_cases i
  · simpa using hp1
  · simp; linarith

/-! ### Binary entropy -/

/-- The Shannon entropy of the binary distribution `binDist p`. -/
noncomputable def binEntropy (p : ℝ) : ℝ :=
  shannonH (Finset.univ : Finset (Fin 2)) (binDist p)

/-- Pointwise expansion of `binEntropy` in terms of `negMulLog`. -/
lemma binEntropy_eq (p : ℝ) :
    binEntropy p = Real.negMulLog p + Real.negMulLog (1 - p) := by
  unfold binEntropy shannonH
  rw [Fin.sum_univ_two]
  simp [binDist]

/-- Pointwise expansion of `binEntropy` in terms of `Real.log`. -/
lemma binEntropy_eq_log (p : ℝ) :
    binEntropy p = -(p * Real.log p) - ((1 - p) * Real.log (1 - p)) := by
  rw [binEntropy_eq]
  unfold Real.negMulLog
  ring

/-! ### Symmetry -/

/-- **Symmetry of binary entropy.** `H_bin(p) = H_bin(1 - p)`. -/
theorem binEntropy_symm (p : ℝ) : binEntropy p = binEntropy (1 - p) := by
  rw [binEntropy_eq, binEntropy_eq]
  have h : (1 : ℝ) - (1 - p) = p := by ring
  rw [h, add_comm]

/-! ### Endpoint values: `binEntropy 0 = 0`, `binEntropy 1 = 0` -/

/-- `binEntropy 0 = 0` (under the convention `0 · log 0 = 0` and
    `Real.log 0 = 0` of Mathlib). -/
@[simp] theorem binEntropy_zero : binEntropy 0 = 0 := by
  rw [binEntropy_eq]
  simp [Real.negMulLog]

/-- `binEntropy 1 = 0`. -/
@[simp] theorem binEntropy_one : binEntropy 1 = 0 := by
  rw [binEntropy_eq]
  simp [Real.negMulLog]

/-! ### Equilibrium: `binEntropy (1/2) = log 2` -/

/-- **Equilibrium identity.** `H_bin(1/2) = log 2`.

    This is the Bekenstein maximum at the PT-symmetry point `s = 1/2`.
    It follows from `shannonH_uniform_eq_log` applied to `Fin 2`. -/
theorem binEntropy_half : binEntropy ((1 : ℝ) / 2) = Real.log 2 := by
  -- At `p = 1/2`, `binDist p = U_2` (the uniform distribution on `Fin 2`).
  have hext : binDist ((1 : ℝ) / 2) = (fun _ : Fin 2 => (1 : ℝ) / 2) := by
    funext i
    fin_cases i
    · rfl
    · simp [binDist]; ring
  unfold binEntropy
  rw [hext]
  apply shannonH_uniform_eq_log (Finset.univ : Finset (Fin 2)) 2
  · norm_num
  · simp

/-! ### Bekenstein bound: `binEntropy p ≤ log 2` on `[0, 1]` -/

/-- **Binary Bekenstein bound.** For `p ∈ [0, 1]`, `H_bin(p) ≤ log 2`.

    Proved by combining the GFT identity at `m = 2` with the
    non-negativity of KL divergence on `[0, 1]`. -/
theorem binEntropy_le_log_two {p : ℝ} (hp0 : 0 ≤ p) (hp1 : p ≤ 1) :
    binEntropy p ≤ Real.log 2 := by
  -- Apply the GFT identity at `m = 2`:
  --   klToUniform + binEntropy = log 2,
  -- and the fact that klToUniform ≥ 0 for distributions on [0, 1]
  -- (KL is non-negative; equivalently, `binEntropy ≤ log m`).
  -- We derive `binEntropy ≤ log 2` via the dual of Bekenstein:
  -- since `klToUniform ≥ 0` (negative-of-Bekenstein lower bound),
  -- `binEntropy = log 2 - klToUniform ≤ log 2`.
  have hid : klToUniform (Finset.univ : Finset (Fin 2)) 2 (binDist p)
              + shannonH (Finset.univ : Finset (Fin 2)) (binDist p)
              = Real.log 2 :=
    GFT_identity (Finset.univ : Finset (Fin 2)) 2 (by norm_num) (binDist p)
      (binDist_nonneg hp0 hp1) (binDist_sum_eq_one p)
  -- Lower-bound on `klToUniform`: use Bekenstein-symmetric reasoning.
  -- We show `0 ≤ klToUniform` directly via expansion.
  -- klToUniform = p * log(2p) + (1-p) * log(2(1-p)).
  -- Each term: p log (2p) ≥ -p log e ≥ ... actually simpler:
  -- We use the alternative route: H ≤ log m from `shannonH_nonneg` of
  -- the *complementary* allocation. But the cleanest path is to write
  -- H_bin = log 2 - klToUniform and bound klToUniform ≥ 0 by the
  -- Gibbs inequality on a 2-point sample.
  --
  -- Concretely: klToUniform = sum_i p_i log(2 p_i).
  -- For p_i ∈ [0,1] with sum 1, this is the KL divergence to U_2,
  -- which is ≥ 0 by Jensen / Gibbs. We give a direct proof using
  -- the convexity of `x ↦ x log x` via `negMulLog_le_log_of` style
  -- arguments — but in Mathlib the cleanest tool is to bound H by
  -- the uniform's entropy, which is `shannonH_uniform_eq_log = log 2`.
  -- That direction requires Jensen.  Here we take the most elementary
  -- route: combine `hid` with the fact that each summand of
  -- `klToUniform` is ≥ -1/e (and tighter), reducing to the
  -- specialised bound `negMulLog_le_one_div_e` chain — which is awkward.
  --
  -- Cleanest available route: invoke `Real.inner_le_nnorm_mul_nnorm`-style
  -- arguments? No: simplest is to reduce to the strict known result
  -- `negMulLog_le_log` (Mathlib): for `x ∈ [0,1]`, `-x log x ≤ log 2 / 2`
  -- ... that route is not in Mathlib in a directly applicable form.
  --
  -- Fallback: expand `H_bin` and reduce to the scalar
  -- `negMulLog p + negMulLog (1-p) ≤ log 2` via the explicit
  -- maximum at `p = 1/2`. We prove this via the AM-GM-style bound:
  -- using concavity of `negMulLog` (which is concave on [0,1])
  -- and Jensen at two points with weights (1/2, 1/2) on `(p, 1-p)`:
  -- (1/2) negMulLog p + (1/2) negMulLog (1-p) ≤ negMulLog ((p+(1-p))/2)
  --   = negMulLog (1/2) = (log 2)/2.
  -- Hence negMulLog p + negMulLog (1-p) ≤ log 2.
  clear hid
  rw [binEntropy_eq]
  -- Reduce to: `negMulLog p + negMulLog (1-p) ≤ log 2`.
  -- We use concavity of `negMulLog` (Mathlib: `concaveOn_negMulLog`)
  -- with weights `(1/2, 1/2)` on `(p, 1-p)`.
  have h1p_nn : (0 : ℝ) ≤ 1 - p := by linarith
  have hp_mem : p ∈ Set.Ici (0 : ℝ) := hp0
  have h1p_mem : (1 - p) ∈ Set.Ici (0 : ℝ) := h1p_nn
  have key := Real.concaveOn_negMulLog.2 hp_mem h1p_mem
    (show (0 : ℝ) ≤ 1 / 2 by norm_num)
    (show (0 : ℝ) ≤ 1 / 2 by norm_num)
    (show (1 : ℝ) / 2 + 1 / 2 = 1 by norm_num)
  -- key : (1/2) • negMulLog p + (1/2) • negMulLog (1-p)
  --        ≤ negMulLog ((1/2) • p + (1/2) • (1-p))
  simp only [smul_eq_mul] at key
  -- The midpoint is `1/2`.
  have hmid : (1 : ℝ) / 2 * p + (1 : ℝ) / 2 * (1 - p) = 1 / 2 := by ring
  rw [hmid] at key
  -- `negMulLog (1/2) = (log 2) / 2`.
  have hneg_half : Real.negMulLog ((1 : ℝ) / 2) = Real.log 2 / 2 := by
    unfold Real.negMulLog
    have hlog : Real.log ((1 : ℝ) / 2) = -Real.log 2 := by
      rw [one_div]; exact Real.log_inv 2
    rw [hlog]; ring
  rw [hneg_half] at key
  linarith

/-! ### Explicit value at `p = 1/3` -/

/-- **Explicit value at `p = 1/3`.**

    `H_bin(1/3) = log 3 - (2/3) log 2`.

    This is the entropy of the flat T1 prior on the three active primes
    seen as a `(1/3, 2/3)` split. -/
theorem binEntropy_third :
    binEntropy ((1 : ℝ) / 3) = Real.log 3 - ((2 : ℝ) / 3) * Real.log 2 := by
  rw [binEntropy_eq_log]
  -- 1 - 1/3 = 2/3
  have h1m : (1 : ℝ) - 1/3 = 2/3 := by ring
  rw [h1m]
  -- log(1/3) = -log 3
  have hlog3 : Real.log ((1 : ℝ) / 3) = -Real.log 3 := by
    rw [one_div]; exact Real.log_inv 3
  -- log(2/3) = log 2 - log 3
  have hlog23 : Real.log ((2 : ℝ) / 3) = Real.log 2 - Real.log 3 := by
    rw [Real.log_div (by norm_num) (by norm_num)]
  rw [hlog3, hlog23]
  ring

/-! ### Explicit value at `p = 13/15` -/

/-- **Explicit (algebraic) value at `p = 13/15`.**

    `H_bin(13/15) = -(13/15) log (13/15) - (2/15) log (2/15)`.

    Kept in algebraic form: the constituents `log 13`, `log 15`, `log 2`
    are not collapsed since `13` is prime and `15 = 3 · 5` shares no
    factor with `13`. -/
theorem binEntropy_thirteen_fifteenths :
    binEntropy ((13 : ℝ) / 15)
      = -((13 : ℝ) / 15) * Real.log ((13 : ℝ) / 15)
        - ((2 : ℝ) / 15) * Real.log ((2 : ℝ) / 15) := by
  rw [binEntropy_eq_log]
  -- 1 - 13/15 = 2/15
  have h1m : (1 : ℝ) - 13/15 = 2/15 := by norm_num
  rw [h1m]
  ring

/-! ### Headline -/

/-- **Headline (binary entropy at PT-canonical points).**

    The binary entropy `H_bin(p) := -p log p - (1-p) log (1-p)` satisfies:

    * `H_bin` is symmetric: `H_bin(p) = H_bin(1 - p)`.
    * `H_bin(0) = H_bin(1) = 0`.
    * `H_bin(1/2) = log 2` (Bekenstein maximum on `Fin 2`).
    * On `[0, 1]`, `H_bin(p) ≤ log 2`.
    * `H_bin(1/3) = log 3 - (2/3) log 2` (flat T1 prior split).
    * `H_bin(13/15)` is given by the algebraic form above (T1 equilibrium
      weight on the active-prime sector). -/
theorem binary_entropy_headline :
    binEntropy 0 = 0
    ∧ binEntropy 1 = 0
    ∧ binEntropy ((1 : ℝ) / 2) = Real.log 2
    ∧ (∀ p, binEntropy p = binEntropy (1 - p))
    ∧ binEntropy ((1 : ℝ) / 3) = Real.log 3 - ((2 : ℝ) / 3) * Real.log 2 :=
  ⟨binEntropy_zero, binEntropy_one, binEntropy_half,
   binEntropy_symm, binEntropy_third⟩

end PT.Information
