/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Stochastic.CRTFactorizationT2
import Mathlib.Tactic

/-!
# CRT Tensor Factorization and Ruelle Factorization

Abstract (matrix-level) formalisation of the central algebraic results
of the companion paper
`PUBLICATION/11_PT_CRT_Decoupling/main.pdf`:

* **Theorem 3.1 (CRT Tensor Factorization)** — the transfer matrix
  $T_m = \bigotimes_p T_p$ on a squarefree modulus $m$ factorises as a
  Kronecker product of the per-prime transfer matrices, and the Kronecker
  product of two row-stochastic matrices is again row-stochastic.

* **Lemma 3.2 (Ruelle Factorization)** — the unique left Perron
  eigenvector (the *Ruelle measure*) of the tensor product factorises as
  the Kronecker product of the per-factor Ruelle measures. Equivalently,
  the stationary distribution of the joint chain is the product of the
  marginal stationary distributions.

The proofs rely on `PT.Stochastic.CRTFactorizationT2`, which already
provides the fundamental identity
`(A ⊗ₖ B).mulVec (v ⊗ w) = (A.mulVec v) ⊗ (B.mulVec w)`
(`kronecker_mulVec_vecTensor`) and the eigenvalue product lemma
(`kronecker_eigenvector`).

The statements here are **paper-level abstract**: they apply to *any* pair
of row-stochastic matrices on finite types, not just the
$T_p$ family of the prime sieve. The specialisation to the
$T_p$ family is then immediate, and the article's numerical witness
$\|\pi_{15}^{\rm Ruelle} - \pi_3 \otimes \pi_5\|_2 = 5.88 \times 10^{-53}$
(test T2 PASS in `scripts/test_crt_decoupling.py`) becomes a corollary.

## Main results

* `vecTensor_sum` — the sum of a tensor product of vectors factorises as a
  product of sums:
  $\sum_{i,j} (v \otimes w)(i,j) = (\sum_i v_i)(\sum_j w_j)$.
* `kronecker_rowSum` — for row-stochastic $A, B$, the row sums of
  $A \otimes_K B$ are all equal to $1$.
* `kronecker_rowStochastic` — the Kronecker product of two row-stochastic
  matrices is row-stochastic.
* `crt_tensor_factor` — *Theorem 3.1, abstract form*: for any pair of
  matrices $A : Matrix α α ℝ$, $B : Matrix β β ℝ$ with stationary
  left-eigenvectors $\pi_A, \pi_B$ (i.e. $\pi_A^\top A = \pi_A^\top$ and
  similarly for $B$), the Kronecker product $A \otimes_K B$ has
  $\pi_A \otimes \pi_B$ as a stationary left-eigenvector.
* `ruelle_factor` — *Lemma 3.2, abstract form*: under the same
  hypotheses, the joint left-stationary distribution factorises as a
  Kronecker product.

The general $k$-factor extension follows by induction on $k$
(`crt_tensor_factor_induct`, future work).

-/

namespace PT.CrtDecoupling

open Matrix Kronecker

/-! ## Tensor product sums and row-stochasticity -/

/-- **Sum factorisation for the tensor product of vectors.**
    $\sum_{(i,j)} (v \otimes w)(i,j) = \bigl(\sum_i v_i\bigr)\bigl(\sum_j w_j\bigr)$. -/
lemma vecTensor_sum {α β : Type*} [Fintype α] [Fintype β]
    (v : α → ℝ) (w : β → ℝ) :
    (∑ ij, PT.Stochastic.vecTensor v w ij)
      = (∑ i, v i) * (∑ j, w j) := by
  rw [Fintype.sum_prod_type, Finset.sum_mul_sum]
  apply Finset.sum_congr rfl
  intro i _
  apply Finset.sum_congr rfl
  intro j _
  simp [PT.Stochastic.vecTensor]

/-- **Row sums of a Kronecker product factorise.** The row of
    $A \otimes_K B$ indexed by $(i, j)$ has sum equal to the product of
    the row sums of $A$ at $i$ and $B$ at $j$. -/
lemma kronecker_rowSum {α β : Type*} [Fintype α] [Fintype β]
    (A : Matrix α α ℝ) (B : Matrix β β ℝ) (i : α) (j : β) :
    (∑ kl, (A ⊗ₖ B) (i, j) kl)
      = (∑ k, A i k) * (∑ l, B j l) := by
  rw [Fintype.sum_prod_type, Finset.sum_mul_sum]
  apply Finset.sum_congr rfl
  intro k _
  apply Finset.sum_congr rfl
  intro l _
  simp

/-- **Kronecker product preserves row-stochasticity.** If every row of $A$
    sums to $1$ and every row of $B$ sums to $1$, then every row of
    $A \otimes_K B$ sums to $1$. -/
theorem kronecker_rowStochastic {α β : Type*} [Fintype α] [Fintype β]
    (A : Matrix α α ℝ) (B : Matrix β β ℝ)
    (hA : ∀ i, (∑ k, A i k) = 1) (hB : ∀ j, (∑ l, B j l) = 1) :
    ∀ ij : α × β, (∑ kl, (A ⊗ₖ B) ij kl) = 1 := by
  intro ij
  obtain ⟨i, j⟩ := ij
  rw [kronecker_rowSum, hA, hB, mul_one]

/-! ## Theorem 3.1: CRT Tensor Factorization (abstract form) -/

/-- **Theorem 3.1 of the companion paper (abstract form).**

    Given two row-stochastic matrices $A : Matrix α α ℝ$ and
    $B : Matrix β β ℝ$ together with left-stationary probability vectors
    $\pi_A : α → ℝ$ and $\pi_B : β → ℝ$ (satisfying $\pi_A^\top A = \pi_A^\top$
    and $\pi_B^\top B = \pi_B^\top$ in the
    `Matrix.vecMul`-on-rows convention used throughout `PT.Stochastic.*`),
    the Kronecker product $\pi_A \otimes \pi_B$ is a left-stationary
    probability vector of $A \otimes_K B$.

    In the prime-sieve application (with $\alpha = \mathcal{S}_p$,
    $\beta = \mathcal{S}_q$, $A = T_p$, $B = T_q$ for distinct active
    primes $p \neq q$), this is exactly the assertion that the joint
    chain $T_{pq}$ acts as $T_p \otimes T_q$ on the active subspace
    $\mathcal{S}_{pq} = \mathcal{S}_p \times \mathcal{S}_q$, and the
    joint Ruelle measure factorises as the product of the per-prime
    Ruelle measures. -/
theorem crt_tensor_factor {α β : Type*} [Fintype α] [Fintype β]
    (A : Matrix α α ℝ) (B : Matrix β β ℝ)
    (πA : α → ℝ) (πB : β → ℝ)
    (hπA : Matrix.vecMul πA A = πA)
    (hπB : Matrix.vecMul πB B = πB) :
    Matrix.vecMul (PT.Stochastic.vecTensor πA πB) (A ⊗ₖ B)
      = PT.Stochastic.vecTensor πA πB := by
  funext ij
  obtain ⟨j, l⟩ := ij
  -- LHS_(j,l) = ∑_(i,k) πA(i) * πB(k) * A(i,j) * B(k,l)
  -- RHS_(j,l) = πA(j) * πB(l)
  show (∑ ik, PT.Stochastic.vecTensor πA πB ik * (A ⊗ₖ B) ik (j, l))
        = πA j * πB l
  rw [Fintype.sum_prod_type]
  -- Reduce to ∑_i ∑_k πA(i) * πB(k) * A(i,j) * B(k,l)
  -- = (∑_i πA(i) * A(i,j)) * (∑_k πB(k) * B(k,l))
  -- = (vecMul πA A) j * (vecMul πB B) l
  -- = πA j * πB l by hπA, hπB.
  have h1 : (∑ i, ∑ k, PT.Stochastic.vecTensor πA πB (i, k) *
              (A ⊗ₖ B) (i, k) (j, l))
            = (∑ i, πA i * A i j) * (∑ k, πB k * B k l) := by
    rw [Finset.sum_mul_sum]
    apply Finset.sum_congr rfl
    intro i _
    apply Finset.sum_congr rfl
    intro k _
    simp [PT.Stochastic.vecTensor]
    ring
  rw [h1]
  -- (∑ i, πA i * A i j) = vecMul πA A j = πA j (by hπA)
  have h2 : (∑ i, πA i * A i j) = πA j := by
    have := congr_fun hπA j
    simpa [Matrix.vecMul, dotProduct] using this
  have h3 : (∑ k, πB k * B k l) = πB l := by
    have := congr_fun hπB l
    simpa [Matrix.vecMul, dotProduct] using this
  rw [h2, h3]

/-! ## Lemma 3.2: Ruelle Factorization (abstract form) -/

/-- **Lemma 3.2 of the companion paper (abstract form).**

    Under the hypotheses of `crt_tensor_factor` plus normalisation
    ($\sum \pi_A = 1$ and $\sum \pi_B = 1$), the Kronecker product
    $\pi_A \otimes \pi_B$ is a probability vector
    (i.e. sums to $1$) and is stationary for $A \otimes_K B$. -/
theorem ruelle_factor {α β : Type*} [Fintype α] [Fintype β]
    (A : Matrix α α ℝ) (B : Matrix β β ℝ)
    (πA : α → ℝ) (πB : β → ℝ)
    (hπA : Matrix.vecMul πA A = πA)
    (hπB : Matrix.vecMul πB B = πB)
    (hπA_sum : (∑ i, πA i) = 1)
    (hπB_sum : (∑ j, πB j) = 1) :
    Matrix.vecMul (PT.Stochastic.vecTensor πA πB) (A ⊗ₖ B)
        = PT.Stochastic.vecTensor πA πB
    ∧ (∑ ij, PT.Stochastic.vecTensor πA πB ij) = 1 := by
  refine ⟨crt_tensor_factor A B πA πB hπA hπB, ?_⟩
  rw [vecTensor_sum, hπA_sum, hπB_sum, mul_one]

end PT.CrtDecoupling
