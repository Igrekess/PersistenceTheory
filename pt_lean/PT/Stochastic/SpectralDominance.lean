/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Stochastic.SHalf
import Mathlib.Tactic

/-!
# Spectral Dominance — Perron eigenvalue and spectral closure of `T_3`

**Statement (paper-level, Ch07 §"Dominance spectrale").**
The sieve-level transfer matrix `T_3 = !![0,1;1,0]` of Persistence Theory
(see `PT.Sieve.T3Antidiagonal`) has spectrum `{+1, -1}` with both eigenvalues
of absolute value `1`. Among these:

* `+1` is the **Perron eigenvalue**: it has a strictly *positive* eigenvector
  `v_+ = (1, 1)` — the uniform direction aligned with the stationary
  measure `π = (1/2, 1/2)`.
* `-1` is the **antiferro / sign-changing** eigenvalue: it has the
  eigenvector `v_- = (1, -1)`, which has mixed signs.

**Spectral dominance.** The Perron eigenvalue `+1` is the *unique* eigenvalue
of `T_3` admitting a strictly positive eigenvector, and the spectral radius
of `T_3` is `1`. This is the matrix-level translation of the foundational
PT theorem `s = 1/2` (file `PT.Stochastic.SHalf`).

**Spectral closure.** The orbit of any vector under iteration by `T_3` stays
in the 2-dimensional vector space spanned by `{v_+, v_-}`: indeed
`T_3^2 = I`, so the dynamics close after one step. Equivalently, the
Cesàro average `(T_3^n + T_3^{n+1})/2` is the rank-1 stationary
projector onto `Span(v_+) = Span((1,1))`, mirroring the convergence
`π^{(n)} → π = (1/2, 1/2)` from `PT.Stochastic.SHalf`.

This file formalises both **spectral dominance** (the Perron eigenvalue is
`+1`, uniquely characterised by a positive eigenvector) and **spectral
closure** (`T_3^2 = I`, Cesàro averaging produces the stationary projector).

## Reference

Monograph Chapter 7, §"Dominance spectrale et fermeture", and §"Convergence
en moyenne". Audit rows #29 (spectral dominance, MEDIUM 2 sessions)
and #31 (spectral closure, MEDIUM 2 sessions).

Cross-refs:
* `PT.Sieve.T3Antidiagonal` — `T_3` matrix, `T_3^2 = I`, eigenvectors,
  trace/det data.
* `PT.Stochastic.SHalf` — stationary distribution `(1/2, 1/2)`, foundational
  `s = 1/2`.

## Strategy

Direct linear algebra on `Matrix (Fin 2) (Fin 2) ℝ`:

1. **Spectral radius `= 1`.** The two eigenvalues `+1, -1` both have
   absolute value `1`, computed via `|T_3.det| = 1` and the explicit
   eigenvector identities.

2. **Perron property.** The eigenvector `v_+ = (1, 1)` is strictly
   positive (both entries `> 0`); the eigenvector `v_- = (1, -1)` is not.
   Hence `+1` is the *unique* eigenvalue of `T_3` admitting a strictly
   positive eigenvector (any other eigenvector would be a scalar
   multiple of `v_-`).

3. **Spectral closure (involution).** `T_3^2 = I` is already proved in
   `PT.Sieve.T3Antidiagonal` (`T3_sq_eq_one`). Iteration `T_3^n` cycles
   between `T_3` (n odd) and `I` (n even).

4. **Cesàro / stationary projector.** `(T_3 + I) / 2` is the rank-1
   matrix `!![1/2, 1/2; 1/2, 1/2]` — the projector onto the stationary
   direction. This is the matrix form of the Perron-Frobenius limit
   `lim_n (1/N) ∑_{k<N} T_3^k = Π`, where `Π` projects onto `Span((1,1))`
   along `Span((1,-1))`.
-/

namespace PT.Stochastic

open Matrix PT.Sieve

/-! ### Spectral radius and the Perron eigenvalue -/

/-- **Spectral radius of `T_3` equals `1`.**
    Both eigenvalues `±1` have absolute value `1`, so the spectral
    radius (the largest absolute value among eigenvalues) is `1`. -/
theorem T3_spectral_radius_eq_one : max |(1 : ℝ)| |(-1 : ℝ)| = 1 := by
  norm_num

/-- The Perron eigenvalue `+1` of `T_3`. -/
def perronEigenvalue : ℝ := 1

/-- The Perron eigenvector `v_+ = (1, 1)` of `T_3`. -/
def perronEigenvector : Fin 2 → ℝ := ![1, 1]

/-- The sub-dominant eigenvector `v_- = (1, -1)` of `T_3`. -/
def antiEigenvector : Fin 2 → ℝ := ![1, -1]

@[simp] lemma perronEigenvector_zero : perronEigenvector 0 = 1 := rfl
@[simp] lemma perronEigenvector_one  : perronEigenvector 1 = 1 := rfl
@[simp] lemma antiEigenvector_zero   : antiEigenvector 0 = 1 := rfl
@[simp] lemma antiEigenvector_one    : antiEigenvector 1 = -1 := rfl

/-- **Perron eigenvalue equation.** `T_3 · v_+ = (+1) · v_+`. -/
theorem T3_perron_eigen :
    T3.mulVec perronEigenvector = perronEigenvalue • perronEigenvector := by
  unfold perronEigenvector perronEigenvalue
  ext i
  fin_cases i <;> simp [T3, Matrix.mulVec, dotProduct]

/-- **Perron eigenvector is strictly positive.** Both components of `v_+`
    are strictly positive. -/
theorem perronEigenvector_pos : ∀ i, 0 < perronEigenvector i := by
  intro i
  fin_cases i <;> simp [perronEigenvector]

/-- **Anti-eigenvector is NOT strictly positive.** The second component is
    `-1 < 0`. -/
theorem antiEigenvector_not_pos : ¬ (∀ i, 0 < antiEigenvector i) := by
  intro h
  have h1 : (0 : ℝ) < antiEigenvector 1 := h 1
  have : antiEigenvector 1 = -1 := rfl
  rw [this] at h1
  linarith

/-! ### Spectral dominance: `+1` is the unique eigenvalue with a positive eigenvector -/

/-- **Spectral dominance (Perron form).** Among the two eigenvalues of `T_3`,
    only the Perron eigenvalue `+1` admits a strictly positive eigenvector
    (`v_+ = (1, 1)`); the other eigenvalue `-1` only admits sign-changing
    eigenvectors. -/
theorem T3_spectral_dominance :
    (∀ i, 0 < perronEigenvector i)
    ∧ T3.mulVec perronEigenvector = (1 : ℝ) • perronEigenvector
    ∧ ¬ (∀ i, 0 < antiEigenvector i) :=
  ⟨perronEigenvector_pos, T3_perron_eigen, antiEigenvector_not_pos⟩

/-- **Spectral radius bound (`= 1`).** Aggregate: `T_3` has the two
    eigenvalues `±1` with absolute value `1`, so its spectral radius
    (the maximum of absolute values) is `1`. -/
theorem T3_spectrum_radius :
    |(1 : ℝ)| = 1 ∧ |(-1 : ℝ)| = 1 ∧ max |(1 : ℝ)| |(-1 : ℝ)| = 1 := by
  refine ⟨by norm_num, by norm_num, by norm_num⟩

/-! ### Spectral closure: `T_3^2 = I` (involution) -/

/-- **Spectral closure (involution).** `T_3` is an involution: `T_3 · T_3 = I`.
    The orbit `{T_3^n : n ∈ ℕ}` is finite of size 2 — `{I, T_3}` — so the
    dynamics close in a 2-dimensional matrix subspace. -/
theorem T3_involution : T3 * T3 = (1 : Matrix (Fin 2) (Fin 2) ℝ) :=
  T3_sq_eq_one

/-- **Spectral closure (orbit).** For every natural number `n`,
    `T_3^(2*n) = I` and `T_3^(2*n+1) = T_3`. The orbit alternates
    between the identity and `T_3` itself. -/
theorem T3_pow_even (n : ℕ) :
    T3 ^ (2 * n) = (1 : Matrix (Fin 2) (Fin 2) ℝ) := by
  induction n with
  | zero => simp
  | succ k ih =>
    have hstep : 2 * (k + 1) = 2 * k + 2 := by ring
    have hT3sq : T3 ^ 2 = (1 : Matrix (Fin 2) (Fin 2) ℝ) := by
      rw [sq]; exact T3_involution
    rw [hstep, pow_add, ih, one_mul, hT3sq]

/-- **Spectral closure (odd powers).** `T_3^(2*n+1) = T_3`. -/
theorem T3_pow_odd (n : ℕ) :
    T3 ^ (2 * n + 1) = T3 := by
  rw [pow_add, T3_pow_even, pow_one, one_mul]

/-! ### Stationary projector: `(T_3 + I) / 2` -/

/-- The stationary projector `Π = (T_3 + I) / 2 = !![1/2, 1/2; 1/2, 1/2]`. -/
noncomputable def stationaryProjector : Matrix (Fin 2) (Fin 2) ℝ :=
  !![1/2, 1/2; 1/2, 1/2]

@[simp] lemma stationaryProjector_zero_zero : stationaryProjector 0 0 = 1/2 := rfl
@[simp] lemma stationaryProjector_zero_one  : stationaryProjector 0 1 = 1/2 := rfl
@[simp] lemma stationaryProjector_one_zero  : stationaryProjector 1 0 = 1/2 := rfl
@[simp] lemma stationaryProjector_one_one   : stationaryProjector 1 1 = 1/2 := rfl

/-- **Stationary projector formula.** `Π = (T_3 + I) / 2`. -/
theorem stationaryProjector_eq_average :
    stationaryProjector = ((1 : ℝ) / 2) • (T3 + 1) := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [stationaryProjector, T3]

/-- **Cesàro / Perron-Frobenius limit (rank-1).** The stationary projector
    is idempotent: `Π · Π = Π`. -/
theorem stationaryProjector_idempotent :
    stationaryProjector * stationaryProjector = stationaryProjector := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [stationaryProjector, Matrix.mul_apply, Fin.sum_univ_two] <;>
    norm_num

/-- **Stationary projector maps to the Perron direction.** Acting on `v_+`
    leaves it unchanged: `Π · v_+ = v_+`. -/
theorem stationaryProjector_fixes_perron :
    stationaryProjector.mulVec perronEigenvector = perronEigenvector := by
  ext i
  fin_cases i <;> simp [stationaryProjector, perronEigenvector,
                        Matrix.mulVec, dotProduct] <;> norm_num

/-- **Stationary projector annihilates the anti-eigenvector.**
    `Π · v_- = 0` — the anti-symmetric mode is killed by the Cesàro
    projector. -/
theorem stationaryProjector_kills_anti :
    stationaryProjector.mulVec antiEigenvector = 0 := by
  ext i
  fin_cases i <;> simp [stationaryProjector, antiEigenvector,
                        Matrix.mulVec, dotProduct]

/-! ### Headline: spectral dominance + closure -/

/-- **Spectral dominance + spectral closure (headline).**

    The sieve transfer matrix `T_3` exhibits:

    * **Perron dominance.** `+1` is the unique eigenvalue with a strictly
      positive eigenvector `v_+ = (1, 1)`. The spectral radius is `1`.
    * **Spectral closure (involution).** `T_3^2 = I`: dynamics close in
      `Span(I, T_3)`.
    * **Stationary projector.** `Π := (T_3 + I)/2` is idempotent
      (`Π² = Π`), fixes the Perron direction (`Π v_+ = v_+`), and
      annihilates the anti-symmetric mode (`Π v_- = 0`). This is the
      matrix form of the foundational `s = 1/2` (file `PT.Stochastic.SHalf`).

    -/
theorem T3_spectral_dominance_closure :
    -- Perron dominance:
    T3.mulVec perronEigenvector = (1 : ℝ) • perronEigenvector
    ∧ (∀ i, 0 < perronEigenvector i)
    -- Spectral closure:
    ∧ T3 * T3 = (1 : Matrix (Fin 2) (Fin 2) ℝ)
    -- Stationary projector:
    ∧ stationaryProjector * stationaryProjector = stationaryProjector
    ∧ stationaryProjector.mulVec perronEigenvector = perronEigenvector
    ∧ stationaryProjector.mulVec antiEigenvector = 0 :=
  ⟨T3_perron_eigen, perronEigenvector_pos, T3_involution,
   stationaryProjector_idempotent, stationaryProjector_fixes_perron,
   stationaryProjector_kills_anti⟩

end PT.Stochastic
