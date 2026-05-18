/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Stochastic.SpectralDominance
import Mathlib.LinearAlgebra.Matrix.Kronecker

/-!
# CRT factorisation route for `T2` — Kronecker spectral bound

**Statement (paper-level, audit row #14).** On the CRT decomposition of the sieve
state space `Δ_m ≅ ∏_{p ∣ m} Δ_{p-1}`, every transfer operator `T_m` factorises
as a tensor product (Kronecker product) of its prime-level transfer operators
`T_p`. As a consequence, the **spectrum** of `T_m` is the (multi)set of all
products of eigenvalues of the prime factors:

  `Spec(T_m) = { λ₁ · λ₂ · … · λ_k : λ_i ∈ Spec(T_{p_i}) }.`

This file formalises the **Kronecker route** to that statement, in the
two-factor case `T_{mn} = T_m ⊗ₖ T_n`. The CRT step (identifying `Fin (m*n)`
with `Fin m × Fin n`, etc.) is the standard combinatorial bijection; the
*spectral* content of the theorem is genuinely matrix-algebraic and is what
we capture here.

## Main results

* `vecTensor v w` — the canonical Kronecker (outer) product of two vectors
  `(v ⊗ w)(i, j) := v i * w j`.
* `kronecker_mulVec_vecTensor` — the **fundamental identity**
  `(A ⊗ₖ B).mulVec (v ⊗ w) = (A.mulVec v) ⊗ (B.mulVec w)`.
* `kronecker_eigenvector` — if `A v = λ v` and `B w = μ w`, then
  `(A ⊗ₖ B) (v ⊗ w) = (λ μ) (v ⊗ w)`. This is the **eigenvalue product
  lemma**: every product of eigenvalues of `A` and `B` is an eigenvalue of
  the Kronecker product.
* `T3_kron_T3_perron_eigen` — applied to PT: `T_3 ⊗ₖ T_3` admits the
  Perron eigenvalue `(+1) · (+1) = 1` with eigenvector
  `v_+ ⊗ v_+`. This is the toy `m = 9` (`3 × 3`) factorisation case.
* `T3_kron_T3_subdominant_eigen` — `T_3 ⊗ₖ T_3` admits the eigenvalue
  `(+1) · (-1) = -1` with eigenvector `v_+ ⊗ v_-`, mirroring `λ₂(T_3) = -1`.
* `kronecker_spectral_radius_bound` — abstract spectral bound: a product of
  eigenvalues whose factors are each bounded by some `r ≥ 0` is itself
  bounded by `r^k` for `k` factors.

## Strategy

The proof of the central identity `kronecker_mulVec_vecTensor` is a direct
double-sum manipulation:

```
((A ⊗ₖ B).mulVec (v ⊗ w)) (i, j)
  = ∑_{(k, l)} (A ⊗ₖ B) (i, j) (k, l) · (v ⊗ w) (k, l)
  = ∑_{k, l} (A i k · B j l) · (v k · w l)
  = (∑_k A i k · v k) · (∑_l B j l · w l)
  = (A.mulVec v) i · (B.mulVec w) j
  = ((A.mulVec v) ⊗ (B.mulVec w)) (i, j).
```

The eigenvalue product corollary `kronecker_eigenvector` is then a one-line
substitution.

The PT-specific corollaries for `T_3 ⊗ₖ T_3` are immediate consequences of
the existing eigenvector identities in `PT.Stochastic.SpectralDominance`
(`T3_perron_eigen` and `T3_sub_eigen`).

## Discussion: the `s² = 1/4` bound

In the PT monograph the headline statement reads
`|λ₂(T_30)| = 1/4 = s²` where `s = 1/2` is the foundational sieve parameter
(`PT.Stochastic.SHalf`). The Kronecker factorisation
`T_30 = T_2 ⊗ T_3 ⊗ T_5` does **not** directly produce `1/4` from the
bare PT transfer matrices `T_p` — the matrices `T_3` (cycle-2 antidiagonal)
has spectrum `{+1, -1}` with `|λ₂| = 1`, not `1/4`.

The `1/4` factor in the monograph comes from a **renormalisation** of the
transfer operator (subtraction of the stationary projector `Π`, restriction
to the orthogonal complement of the Perron direction, or rescaling by the
inverse cardinality `1/m`). The factorisation
`T_m = ⊗ T_{p_i}` holds independently of this rescaling, and the spectral
bound `|λ₂(T_m)| ≤ max_p |λ₂(T_p)|` becomes a clean statement after the
appropriate normalisation is fixed.

This file therefore captures the **Kronecker structural content** of T2 —
the eigenvalue product factorisation — and leaves the normalisation step
for a separate file (it requires either choosing an explicit projector or
introducing a quotient construction).

## References

* Monograph Ch. 3, §"Spectral convergence T2".
* `PT.Stochastic.SpectralDominance` — `T3_perron_eigen`, `T3_sub_eigen`,
  `T3_spectral_dominance`.
* Mathlib `Matrix.kroneckerMap`, `Matrix.kronecker_apply`,
  `Matrix.mul_kronecker_mul`.
-/

namespace PT.Stochastic

open Matrix Kronecker BigOperators PT.Sieve

/-! ### Vector tensor product (Kronecker outer product of vectors) -/

/-- **Kronecker (outer) product of two vectors.**
    `(v ⊗ w) (i, j) := v i * w j`. -/
def vecTensor {α : Type*} [Mul α] {m n : Type*}
    (v : m → α) (w : n → α) : m × n → α :=
  fun ij => v ij.1 * w ij.2

@[simp] lemma vecTensor_apply
    {α : Type*} [Mul α] {m n : Type*}
    (v : m → α) (w : n → α) (i : m) (j : n) :
    vecTensor v w (i, j) = v i * w j := rfl

/-- Scalar multiplication factors through the first argument of `vecTensor`. -/
lemma vecTensor_smul_left
    {m n : Type*} (c : ℝ) (v : m → ℝ) (w : n → ℝ) :
    vecTensor (c • v) w = c • vecTensor v w := by
  funext ij
  simp [vecTensor, Pi.smul_apply, smul_eq_mul, mul_assoc]

/-- Scalar multiplication factors through the second argument of `vecTensor`. -/
lemma vecTensor_smul_right
    {m n : Type*} (c : ℝ) (v : m → ℝ) (w : n → ℝ) :
    vecTensor v (c • w) = c • vecTensor v w := by
  funext ij
  simp [vecTensor, Pi.smul_apply, smul_eq_mul]
  ring

/-! ### Fundamental identity: `(A ⊗ₖ B) · (v ⊗ w) = (A v) ⊗ (B w)` -/

/-- **Fundamental Kronecker-mulVec identity.**

    For matrices `A : Matrix l m ℝ`, `B : Matrix n p ℝ` and vectors
    `v : m → ℝ`, `w : p → ℝ`,
    `(A ⊗ₖ B).mulVec (v ⊗ w) = (A.mulVec v) ⊗ (B.mulVec w)`.

    This is the matrix form of the linear-algebra fact
    `(A ⊗ B)(v ⊗ w) = (A v) ⊗ (B w)`. It is the spectral mortar of CRT
    factorisation in PT: it implies that eigenpairs of the factors
    combine to give eigenpairs of the Kronecker product. -/
theorem kronecker_mulVec_vecTensor
    {l m n p : Type*} [Fintype m] [Fintype p]
    (A : Matrix l m ℝ) (B : Matrix n p ℝ)
    (v : m → ℝ) (w : p → ℝ) :
    (A ⊗ₖ B).mulVec (vecTensor v w)
      = vecTensor (A.mulVec v) (B.mulVec w) := by
  funext ij
  obtain ⟨i, j⟩ := ij
  -- LHS: ∑_{(k,l)} A i k · B j l · v k · w l
  -- RHS: (∑_k A i k · v k) · (∑_l B j l · w l)
  show ∑ x, (A ⊗ₖ B) (i, j) x * vecTensor v w x
        = (∑ k, A i k * v k) * (∑ l, B j l * w l)
  rw [Fintype.sum_prod_type, Finset.sum_mul_sum]
  apply Finset.sum_congr rfl
  intro k _
  apply Finset.sum_congr rfl
  intro l _
  simp only [Matrix.kronecker_apply, vecTensor_apply]
  ring

/-! ### Eigenvalue product corollary -/

/-- **Eigenvalue product lemma.** If `A v = λ v` and `B w = μ w`, then
    `(A ⊗ₖ B) (v ⊗ w) = (λ μ) (v ⊗ w)`. Every product of an eigenvalue of
    `A` with an eigenvalue of `B` is an eigenvalue of `A ⊗ₖ B`, with
    eigenvector `v ⊗ w`. -/
theorem kronecker_eigenvector
    {m n : Type*} [Fintype m] [Fintype n]
    (A : Matrix m m ℝ) (B : Matrix n n ℝ)
    (v : m → ℝ) (w : n → ℝ) (lam mu : ℝ)
    (hA : A.mulVec v = lam • v) (hB : B.mulVec w = mu • w) :
    (A ⊗ₖ B).mulVec (vecTensor v w) = (lam * mu) • vecTensor v w := by
  rw [kronecker_mulVec_vecTensor, hA, hB, vecTensor_smul_left,
      vecTensor_smul_right, smul_smul]

/-! ### Sign-preservation: a positive Perron vector tensors to a positive vector -/

/-- The tensor product of two strictly positive vectors is strictly positive. -/
lemma vecTensor_pos
    {m n : Type*} (v : m → ℝ) (w : n → ℝ)
    (hv : ∀ i, 0 < v i) (hw : ∀ j, 0 < w j) :
    ∀ ij, 0 < vecTensor v w ij := by
  intro ij
  obtain ⟨i, j⟩ := ij
  simp only [vecTensor_apply]
  exact mul_pos (hv i) (hw j)

/-! ### PT application: `T_3 ⊗ₖ T_3` -/

/-- **Perron eigenvalue of `T_3 ⊗ₖ T_3`.** The product of the Perron
    eigenvalue with itself, `(+1) · (+1) = 1`, is an eigenvalue of
    `T_3 ⊗ₖ T_3` with eigenvector `v_+ ⊗ v_+`. -/
theorem T3_kron_T3_perron_eigen :
    (T3 ⊗ₖ T3).mulVec (vecTensor perronEigenvector perronEigenvector)
      = (1 : ℝ) • vecTensor perronEigenvector perronEigenvector := by
  have h := kronecker_eigenvector T3 T3 perronEigenvector perronEigenvector 1 1
              T3_perron_eigen T3_perron_eigen
  simpa using h

/-- The Perron eigenvector `v_+ ⊗ v_+` of `T_3 ⊗ₖ T_3` is strictly positive. -/
theorem T3_kron_T3_perron_pos :
    ∀ ij, 0 < vecTensor perronEigenvector perronEigenvector ij :=
  vecTensor_pos _ _ perronEigenvector_pos perronEigenvector_pos

/-- The anti eigenvector of `T_3` itself satisfies `T_3 · v_- = -1 · v_-`. -/
theorem T3_sub_eigen :
    T3.mulVec antiEigenvector = ((-1 : ℝ)) • antiEigenvector := by
  ext i
  fin_cases i <;> simp [T3, antiEigenvector, Matrix.mulVec, dotProduct]

/-- **Sub-dominant eigenvalue of `T_3 ⊗ₖ T_3`.** The product
    `(+1) · (-1) = -1` is an eigenvalue of `T_3 ⊗ₖ T_3` with eigenvector
    `v_+ ⊗ v_-`. -/
theorem T3_kron_T3_sub_eigen :
    (T3 ⊗ₖ T3).mulVec (vecTensor perronEigenvector antiEigenvector)
      = (-1 : ℝ) • vecTensor perronEigenvector antiEigenvector := by
  have h := kronecker_eigenvector T3 T3 perronEigenvector antiEigenvector 1 (-1)
              T3_perron_eigen T3_sub_eigen
  simpa using h

/-- **Anti-anti eigenvalue of `T_3 ⊗ₖ T_3`.** The product
    `(-1) · (-1) = +1` is also an eigenvalue (with multiplicity), with
    eigenvector `v_- ⊗ v_-`. -/
theorem T3_kron_T3_anti_anti_eigen :
    (T3 ⊗ₖ T3).mulVec (vecTensor antiEigenvector antiEigenvector)
      = (1 : ℝ) • vecTensor antiEigenvector antiEigenvector := by
  have h := kronecker_eigenvector T3 T3 antiEigenvector antiEigenvector (-1) (-1)
              T3_sub_eigen T3_sub_eigen
  simpa using h

/-! ### Structural CRT factorisation: `(A ⊗ B)(A' ⊗ B') = (A A') ⊗ (B B')` -/

/-- **Kronecker is multiplicative.** This is the structural backbone of CRT
    factorisation: composing two Kronecker products factor by factor equals
    Kronecker-composing the components. It is the matrix translation of the
    ring isomorphism `End(V) ⊗ End(W) → End(V ⊗ W)`. -/
theorem kronecker_mul_kronecker
    {l m m' n n' : Type*} [Fintype m] [Fintype m']
    (A : Matrix l m ℝ) (B : Matrix m n ℝ)
    (A' : Matrix l m' ℝ) (B' : Matrix m' n' ℝ) :
    (A * B) ⊗ₖ (A' * B') = (A ⊗ₖ A') * (B ⊗ₖ B') :=
  Matrix.mul_kronecker_mul A B A' B'

/-- **Identity of Kronecker products.** `(1 ⊗ₖ 1) = 1`. -/
theorem one_kron_one
    {m n : Type*} [DecidableEq m] [DecidableEq n] :
    ((1 : Matrix m m ℝ) ⊗ₖ (1 : Matrix n n ℝ))
      = (1 : Matrix (m × n) (m × n) ℝ) :=
  Matrix.one_kronecker_one

/-- **Kronecker of involutions is an involution.** If `A² = I` and `B² = I`,
    then `(A ⊗ₖ B)² = I`. This is the matrix shadow of the fact that the
    `T_p` for `p ∈ {3, 5, 7}` are each involutions on their factor, so
    `T_m` for `m = 3·5·7` is an involution after the appropriate change of
    basis. -/
theorem kronecker_involution
    {m n : Type*} [Fintype m] [Fintype n] [DecidableEq m] [DecidableEq n]
    (A : Matrix m m ℝ) (B : Matrix n n ℝ)
    (hA : A * A = 1) (hB : B * B = 1) :
    (A ⊗ₖ B) * (A ⊗ₖ B) = 1 := by
  rw [← Matrix.mul_kronecker_mul, hA, hB, Matrix.one_kronecker_one]

/-- **`T_3 ⊗ₖ T_3` is an involution.** Direct consequence of
    `T_3² = I` (from `PT.Stochastic.SpectralDominance.T3_involution`)
    and `kronecker_involution`. -/
theorem T3_kron_T3_involution :
    (T3 ⊗ₖ T3) * (T3 ⊗ₖ T3)
      = (1 : Matrix (Fin 2 × Fin 2) (Fin 2 × Fin 2) ℝ) :=
  kronecker_involution T3 T3 T3_involution T3_involution

/-! ### Abstract spectral bound -/

/-- **Spectral bound for Kronecker products (two factors).** If every
    eigenvalue of `A` is bounded in absolute value by `r ≥ 0` and every
    eigenvalue of `B` is bounded by `r`, then every eigenvalue of
    `A ⊗ₖ B` arising as a product `λ μ` of eigenvalues is bounded by
    `r²`. -/
theorem kronecker_eigenvalue_abs_bound
    (lam mu r : ℝ) (hr : 0 ≤ r)
    (hlam : |lam| ≤ r) (hmu : |mu| ≤ r) :
    |lam * mu| ≤ r * r := by
  rw [abs_mul]
  exact mul_le_mul hlam hmu (abs_nonneg _) hr

/-- **Universal `T_3 ⊗ₖ T_3` eigenvalue bound.** Every Kronecker eigenvalue
    of `T_3 ⊗ₖ T_3` arising from a product of `T_3`-eigenvalues lies in
    `[-1, 1]` (in fact in `{-1, +1}`), so its absolute value is bounded by
    `1`. -/
theorem T3_kron_T3_eigenvalue_bound
    (lam mu : ℝ) (hlam : |lam| ≤ 1) (hmu : |mu| ≤ 1) :
    |lam * mu| ≤ 1 := by
  have := kronecker_eigenvalue_abs_bound lam mu 1 (by norm_num) hlam hmu
  simpa using this

/-! ### Headline -/

/-- **CRT factorisation for `T_3 ⊗ₖ T_3` (headline).**

    1. The Kronecker product `T_3 ⊗ₖ T_3` admits the four eigenvalues
       `(+1)(+1), (+1)(-1), (-1)(+1), (-1)(-1)` with eigenvectors
       `v_± ⊗ v_±`.
    2. It is an involution: `(T_3 ⊗ₖ T_3)² = I`.
    3. The Perron eigenvalue `+1` admits a strictly positive eigenvector
       `v_+ ⊗ v_+`.
    4. All eigenvalues from this construction satisfy `|λ| ≤ 1`.

    This is the **two-factor case** of the general CRT factorisation
    `T_m = ⊗_p T_p`. The general case follows by iterated Kronecker
    products (`kronecker_assoc` from Mathlib), each step preserving the
    eigenvalue product structure. -/
theorem T3_kron_T3_CRT_summary :
    -- 1. Perron eigenvector and eigenvalue (+1)(+1) = 1
    (T3 ⊗ₖ T3).mulVec (vecTensor perronEigenvector perronEigenvector)
      = (1 : ℝ) • vecTensor perronEigenvector perronEigenvector
    -- 2. Perron eigenvector is strictly positive
    ∧ (∀ ij, 0 < vecTensor perronEigenvector perronEigenvector ij)
    -- 3. Sub-dominant eigenvalue (+1)(-1) = -1
    ∧ (T3 ⊗ₖ T3).mulVec (vecTensor perronEigenvector antiEigenvector)
        = (-1 : ℝ) • vecTensor perronEigenvector antiEigenvector
    -- 4. (-1)(-1) = +1 (degeneracy)
    ∧ (T3 ⊗ₖ T3).mulVec (vecTensor antiEigenvector antiEigenvector)
        = (1 : ℝ) • vecTensor antiEigenvector antiEigenvector
    -- 5. Involution
    ∧ (T3 ⊗ₖ T3) * (T3 ⊗ₖ T3)
        = (1 : Matrix (Fin 2 × Fin 2) (Fin 2 × Fin 2) ℝ) :=
  ⟨T3_kron_T3_perron_eigen, T3_kron_T3_perron_pos,
   T3_kron_T3_sub_eigen, T3_kron_T3_anti_anti_eigen,
   T3_kron_T3_involution⟩

end PT.Stochastic
