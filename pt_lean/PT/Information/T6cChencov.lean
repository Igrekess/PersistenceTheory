/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import Mathlib.Analysis.SpecialFunctions.Log.NegMulLog
import Mathlib.Analysis.Convex.StdSimplex
import Mathlib.Analysis.MeanInequalities
import Mathlib.Topology.Instances.RealVectorSpace
import Mathlib.Tactic

/-!
# T6c — Shannon entropy infrastructure and Cauchy's logarithmic equation

**Scope revision (2026-05-16).** This module formerly carried a headline
theorem `T6c_chencov_characterisation` stating that any
`H : ∀ n, Simplex n → ℝ` satisfying four axioms (continuity, symmetry,
maximum-at-uniform, product-additivity) is a non-negative scalar
multiple of Shannon entropy. **That statement is mathematically false**:
the Rényi family `H_α(p) := (1/(1-α)) · log(∑ pᵢ^α)` satisfies all
four axioms for every `α > 0`, `α ≠ 1`, but is not proportional to
Shannon entropy. (For `α = 2`, `H_2(p, 1-p) = -log(p² + (1-p)²)` is
not a scalar multiple of `binaryEntropy p`.)

The literature documents this: Shore & Johnson 1980 use **five**
axioms (SJ1–SJ5), including a *subset-independence* axiom (SJ4) that
excludes Rényi. The monograph `chapters_fr/ch05_geometry.tex` makes the
exclusion of Rényi explicit ("Divergences de Rényi `D_α` (α ≠ 1) :
violent SJ1 (incohérence des mises à jour séquentielles pour α ≠ 1)").

The correct PT formulation is on **f-divergences** `D(P‖Q)`, not on
entropies directly, and uses Csiszár 1967 instead of the full
Shore–Johnson chain. That formulation lives in
`PT/Information/T6cG1Autonomous.lean`.

This file retains the kernel-verified infrastructure that survives the
revision: the `Simplex` structure and basic lemmas, Shannon entropy
and its properties under permutation / maximum / product additivity
(the "easy direction" — Shannon *does* satisfy the four axioms),
binary entropy parametrisation, and **`cauchy_log_equation`** (the
multiplicative-to-logarithmic Cauchy equation), which remains a
reusable building block for the divergence-side proof in the companion
file.

## Historical note

This module scaffolds the formalisation of **Theorem T6c** of the PT
monograph: on the state space of the persistence sieve, the
Kullback–Leibler divergence and the Fisher information metric are the
*unique* canonical choices, hence Shannon entropy is the unique
admissible information functional.

In the M1 article (`PT_ARTICLES/PT_MATHEMATICS/M1/m1_persistence_fr.tex`,
§"Théorème T6c : géométrie de l'information canonique") and in the
monograph (`PT_MONOGRAPHY/chapters_fr/ch02_uniqueness.tex`, §"Théorème
T6c"), T6c is itself the conjunction of two sub-results:

* **G1 (Shore–Johnson 1980)** — `D_KL` is the unique f-divergence on the
  state space `Δ_m` satisfying the five axioms A1 (coherence),
  A2 (permutation invariance), A3 (system independence),
  A4 (scaling / homogeneity), A5 (subset independence).
* **G3 (Čencov / Chentsov 1982)** — the Fisher information metric is the
  unique Riemannian metric (up to positive multiplicative constant) on
  `Δ_n` that contracts under Markov maps (information monotonicity).

Both are explicitly framed in PT as **reductions to classical external
theorems**, not autonomous PT derivations: the PT content is the
verification that the sieve state space (with its CRT / Markov product
structure) satisfies the premises of Shore–Johnson and Čencov.

Equivalently (and this is the cleanest Lean statement) the Shore–Johnson
axioms applied to an entropy functional `H : Simplex n → ℝ` force
`H = k · shannon` for some positive constant `k`. That is the form we
scaffold below.

## Scope of this first session (~90 min)

This file lays down the **structure** for T6c. The hard uniqueness
direction (Shore–Johnson) is `sorry`'d — its formalisation is estimated
at 10–15 sessions and goes substantially beyond what is in current
Mathlib. What we do close cleanly here:

* The `Simplex n` structure and basic lemmas (membership in `stdSimplex`,
  nonnegativity, sum-to-one).
* `shannon : Simplex n → ℝ`, defined via `Real.negMulLog`.
* The four Shore–Johnson-style axioms as Lean predicates on functionals
  `H : ∀ n, Simplex n → ℝ`.
* Lemmas showing `shannon` **satisfies** each of the four axioms (the
  "easy direction": Shannon is a valid candidate). Permutation
  invariance is closed cleanly; the others are stubbed with documented
  TODOs for the next session (they require routine `Fintype` /
  `negMulLog_mul` calculations).
* The **headline theorem** `T6c_chencov_characterisation` — typechecked,
  body `sorry`.
* The 2-state warm-up `T6c_binary_characterisation` — also `sorry`
  (the binary case is conceptually the same Shore–Johnson chain, just
  carried out on `Fin 2`).

## Sorry budget

| Name                                  | Strategy                                                | Effort      |
|---------------------------------------|---------------------------------------------------------|-------------|
| `shannon_symmetric`                   | `Equiv.sum_comp` over `Finset.univ`                     | 0.5 session |
| `shannon_max_uniform`                 | Jensen on `negMulLog` (concave) + `∑ = 1`               | 1 session   |
| `shannon_additive_product`            | `negMulLog_mul` + `Fintype.sum_prod_type`               | 1 session   |
| `T6c_binary_characterisation`         | Shore–Johnson axiom chain in dimension 2                | 4 sessions  |
| `T6c_chencov_characterisation`        | Full Shore–Johnson; sieve `Δ_m` satisfies A1–A5         | 10 sessions |

Total realistic remaining effort: **~16 sessions** for the full theorem,
or **~6 sessions** if one only wants the binary case plus the three
"easy direction" lemmas (i.e. a complete proof that Shannon is a
*candidate*, with the *uniqueness* direction admitted as an import).

The binary case is what the PT G3 (Čencov) statement actually uses:
the sieve state space at the elementary level is `Δ_2 = {(p, 1-p)}`
(§T6c, M1 article line 1698).

## References

* J. E. Shore and R. W. Johnson, *Axiomatic derivation of the principle
  of maximum entropy and the principle of minimum cross-entropy*, IEEE
  Trans. Inform. Theory **26** (1980), 26–37.
* N. N. Čencov (Chentsov), *Statistical Decision Rules and Optimal
  Inference*, Translations of Mathematical Monographs **53**, AMS, 1982.
* Y. Senez, *La théorie de la persistance*, monograph
  `PT_MONOGRAPHY/chapters_fr/ch02_uniqueness.tex`, §T6c (lines 575–637).
* M1 article `PT_ARTICLES/PT_MATHEMATICS/M1/m1_persistence_fr.tex`,
  §"Théorème T6c" (lines 1665–1728).
-/

namespace PT.Information

open Real Finset
open scoped BigOperators

/-! ## The discrete probability simplex -/

/-- The probability simplex on `Fin n`: nonneg vectors summing to 1.

    We use a `structure` rather than the Mathlib `stdSimplex ℝ (Fin n)`
    subset directly, because the latter is a `Set (Fin n → ℝ)`, which is
    less ergonomic for stating uniqueness theorems quantified over the
    "functional `H : Simplex n → ℝ`". The bridging lemma
    `Simplex.mem_stdSimplex` connects the two views. -/
structure Simplex (n : ℕ) where
  /-- The underlying probability vector. -/
  prob : Fin n → ℝ
  /-- Each entry is nonnegative. -/
  nonneg : ∀ i, 0 ≤ prob i
  /-- The entries sum to one. -/
  sums_to_one : ∑ i, prob i = 1

namespace Simplex

variable {n : ℕ}

/-- A `Simplex n` lies in Mathlib's standard simplex. -/
lemma mem_stdSimplex (p : Simplex n) : p.prob ∈ stdSimplex ℝ (Fin n) :=
  ⟨p.nonneg, p.sums_to_one⟩

/-- Each component of a `Simplex n` is at most `1`. -/
lemma prob_le_one (p : Simplex n) (i : Fin n) : p.prob i ≤ 1 := by
  have h := mem_Icc_of_mem_stdSimplex p.mem_stdSimplex i
  exact h.2

/-- Each component of a `Simplex n` is in `[0, 1]`. -/
lemma prob_mem_Icc (p : Simplex n) (i : Fin n) : p.prob i ∈ Set.Icc (0 : ℝ) 1 :=
  mem_Icc_of_mem_stdSimplex p.mem_stdSimplex i

/-- Extensionality: two simplex distributions are equal iff their probability
    vectors agree pointwise. -/
@[ext]
lemma ext {p q : Simplex n} (h : ∀ i, p.prob i = q.prob i) : p = q := by
  cases p; cases q
  congr
  exact funext h

/-- The uniform distribution on `Fin n` for `n ≥ 1`. -/
noncomputable def uniform (n : ℕ) (hn : 0 < n) : Simplex n where
  prob _ := 1 / (n : ℝ)
  nonneg _ := by positivity
  sums_to_one := by
    rw [Finset.sum_const, card_univ, Fintype.card_fin, nsmul_eq_mul]
    have hn' : (n : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr hn.ne'
    field_simp

/-- Permute the entries of a probability vector along an `Equiv` of `Fin n`. -/
def permute (p : Simplex n) (σ : Equiv.Perm (Fin n)) : Simplex n where
  prob i := p.prob (σ i)
  nonneg i := p.nonneg _
  sums_to_one := by
    -- Σ_i p (σ i) = Σ_i p i
    rw [Equiv.sum_comp σ (fun i => p.prob i), p.sums_to_one]

/-- The product distribution on `Fin (m * n)`, identifying `Fin (m*n)` with
    `Fin m × Fin n` via `finProdFinEquiv`. -/
noncomputable def prod {m n : ℕ} (p : Simplex m) (q : Simplex n) :
    Simplex (m * n) where
  prob k := p.prob (finProdFinEquiv.symm k).1 * q.prob (finProdFinEquiv.symm k).2
  nonneg k := mul_nonneg (p.nonneg _) (q.nonneg _)
  sums_to_one := by
    -- Σ_k p((π₁ ∘ finProdFinEquiv.symm) k) · q((π₂ ∘ finProdFinEquiv.symm) k) = 1
    have hsum := Equiv.sum_comp finProdFinEquiv.symm
      (fun ij : Fin m × Fin n => p.prob ij.1 * q.prob ij.2)
    -- hsum : ∑ k, (fun ij => p.prob ij.1 * q.prob ij.2) (finProdFinEquiv.symm k)
    --        = ∑ ij, p.prob ij.1 * q.prob ij.2
    rw [show (∑ k, p.prob (finProdFinEquiv.symm k).1 *
              q.prob (finProdFinEquiv.symm k).2)
          = ∑ k, (fun ij : Fin m × Fin n => p.prob ij.1 * q.prob ij.2)
              (finProdFinEquiv.symm k) from rfl]
    rw [hsum, Fintype.sum_prod_type]
    simp only [← Finset.mul_sum, ← Finset.sum_mul]
    rw [p.sums_to_one, q.sums_to_one, mul_one]

/-- The product of two uniform distributions is the uniform distribution on the
    product space. Concretely, `(1/m)·(1/n) = 1/(m·n)` at every index. -/
lemma uniform_prod_uniform {m n : ℕ} (hm : 0 < m) (hn : 0 < n) :
    (uniform m hm).prod (uniform n hn) = uniform (m * n) (Nat.mul_pos hm hn) := by
  ext k
  show (1 / (m : ℝ)) * (1 / (n : ℝ)) = 1 / ((m * n : ℕ) : ℝ)
  push_cast
  rw [one_div, one_div, one_div, ← mul_inv]

end Simplex

/-! ## Shannon entropy -/

/-- **Shannon entropy** on the simplex, `H(p) = -∑ p_i log p_i = ∑ negMulLog p_i`.

    Uses Mathlib's `Real.negMulLog x = -x log x`, which is continuous on `[0, ∞)`
    with the convention `negMulLog 0 = 0`. -/
noncomputable def shannon {n : ℕ} (p : Simplex n) : ℝ :=
  ∑ i, Real.negMulLog (p.prob i)

/-- `shannon` is nonnegative on the simplex. -/
lemma shannon_nonneg {n : ℕ} (p : Simplex n) : 0 ≤ shannon p := by
  unfold shannon
  exact Finset.sum_nonneg fun i _ =>
    Real.negMulLog_nonneg (p.nonneg i) (p.prob_le_one i)

/-- Shannon entropy of any deterministic distribution (point mass) is zero.
    Here we record the trivial special case where `n = 1`: there is only one
    distribution and it must be the point mass at the unique index. -/
lemma shannon_singleton (p : Simplex 1) : shannon p = 0 := by
  unfold shannon
  -- The unique index is 0, and p.prob 0 = 1 by sums_to_one
  have h : p.prob 0 = 1 := by
    have hs := p.sums_to_one
    simpa using hs
  simp [h]

/-! ## The four Shore–Johnson axioms as Lean predicates

We follow the formulation in `m1_persistence_fr.tex` lines 1682–1687
(A1 coherence, A2 permutation invariance, A3 system independence,
A4 scaling, A5 subset independence). For the Lean statement we package
A1 (coherence ≈ continuity in `p`), A2 (symmetry), A3 (additivity under
products), and the maximum-at-uniform axiom (a consequence in the
Shore–Johnson framework, taken as an axiom by some authors).

We state these as predicates on a family `H : ∀ n, Simplex n → ℝ`,
since the additivity axiom relates entropies in different dimensions.
-/

/-- **A1 (coherence ≈ continuity).** For each `n`, the functional
    `H n : Simplex n → ℝ` is continuous as a function of the underlying
    probability vector.

    Rather than equip `Simplex n` with a `TopologicalSpace` instance (which
    would force a long boilerplate detour), we state continuity directly in
    `ε–δ` form on the underlying probability vector. This is the operationally
    useful form: it says that small perturbations of the probabilities in
    sup-norm produce small perturbations of `H`. -/
def Continuous_axiom (H : ∀ n, Simplex n → ℝ) : Prop :=
  ∀ n, ∀ ⦃p₀ : Simplex n⦄, ∀ ε > (0 : ℝ), ∃ δ > (0 : ℝ),
    ∀ p : Simplex n,
      (∀ i, |p.prob i - p₀.prob i| < δ) → |H n p - H n p₀| < ε

/-- **A2 (permutation invariance / symmetry).** For each `n` and each
    permutation `σ` of `Fin n`, `H n (p ∘ σ) = H n p`. -/
def Symmetric_axiom (H : ∀ n, Simplex n → ℝ) : Prop :=
  ∀ n (p : Simplex n) (σ : Equiv.Perm (Fin n)),
    H n (p.permute σ) = H n p

/-- **A4 (maximum at uniform).** For each `n ≥ 1`, the uniform distribution
    maximises `H n`. -/
def Maximal_at_uniform_axiom (H : ∀ n, Simplex n → ℝ) : Prop :=
  ∀ n (hn : 0 < n) (p : Simplex n),
    H n p ≤ H n (Simplex.uniform n hn)

/-- **A3 (additivity / system independence).** For independent systems
    of sizes `m` and `n`, the entropy of the product distribution equals
    the sum of the entropies. This is Shore–Johnson's "system
    independence" axiom in its strong (additivity) form. -/
def Additive_axiom (H : ∀ n, Simplex n → ℝ) : Prop :=
  ∀ {m n : ℕ} (p : Simplex m) (q : Simplex n),
    H (m * n) (p.prod q) = H m p + H n q

/-! ## Shannon satisfies the axioms (easy direction)

These four lemmas establish that Shannon entropy is a valid candidate.
The hard direction of Shore–Johnson — uniqueness — is the headline
theorem below.
-/

/-- **Shannon is symmetric (A2).** -/
theorem shannon_symmetric : Symmetric_axiom @shannon := by
  intro n p σ
  unfold shannon Simplex.permute
  -- Σ_i negMulLog (p (σ i)) = Σ_i negMulLog (p i)
  exact Equiv.sum_comp σ (fun i => Real.negMulLog (p.prob i))

/-- **Shannon attains its maximum at the uniform distribution (A4).**

    Standard application of Jensen's inequality to the concave function
    `negMulLog` on `Set.Ici 0`, summed against the uniform weights.

    TODO (~1 session):
    * Use `ConcaveOn.le_sum_of_mem` (Jensen) for
      `Real.concaveOn_negMulLog` on the convex set `Set.Ici 0`.
    * Apply with weights `wᵢ = 1/n` and points `xᵢ = p.prob i`.
    * Conclude `(1/n) · ∑ negMulLog (p i) ≤ negMulLog ((1/n) · ∑ p i)
      = negMulLog (1/n)`, then multiply by `n` and observe
      `∑_i negMulLog (1/n) = n · negMulLog (1/n) = H(uniform)`.
    Mathlib API needed: `ConcaveOn.inner_le_sum` or
    `ConcaveOn.le_map_sum` (cf. `Mathlib.Analysis.Convex.Jensen`). -/
theorem shannon_max_uniform : Maximal_at_uniform_axiom @shannon := by
  intro n hn p
  have hn_pos : (0 : ℝ) < n := Nat.cast_pos.mpr hn
  have hn_ne : (n : ℝ) ≠ 0 := hn_pos.ne'
  -- Sum of uniform weights equals 1.
  have hsum_w : (∑ _ : Fin n, (1 / (n : ℝ))) = 1 := by
    rw [Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul,
        mul_one_div, div_self hn_ne]
  -- Jensen's inequality applied to `negMulLog` (concave on `Set.Ici 0`)
  -- with uniform weights `w i = 1/n` and points `p.prob i`.
  have jensen := Real.concaveOn_negMulLog.le_map_sum
    (t := (Finset.univ : Finset (Fin n)))
    (w := fun _ => (1 / (n : ℝ)))
    (p := p.prob)
    (fun _ _ => by positivity)
    hsum_w
    (fun i _ => Set.mem_Ici.mpr (p.nonneg i))
  -- Convert `•` (smul) to `*` (mul) on ℝ.
  simp only [smul_eq_mul] at jensen
  -- Factor the constant `(1/n)` out of both inner sums (LHS and RHS).
  simp_rw [← Finset.mul_sum] at jensen
  -- Simplify the inner sum on the RHS via `∑ p = 1`.
  rw [p.sums_to_one, mul_one] at jensen
  -- `jensen` now reads: (1/n) * shannon p ≤ negMulLog (1/n).
  -- Multiply by `n > 0` to get `shannon p ≤ n · negMulLog (1/n) = shannon uniform`.
  calc shannon p
      = (n : ℝ) * ((1 / (n : ℝ)) * shannon p) := by
        rw [← mul_assoc, mul_one_div, div_self hn_ne, one_mul]
    _ ≤ (n : ℝ) * Real.negMulLog (1 / (n : ℝ)) :=
        mul_le_mul_of_nonneg_left jensen hn_pos.le
    _ = shannon (Simplex.uniform n hn) := by
        show (n : ℝ) * Real.negMulLog (1 / (n : ℝ))
           = ∑ _ : Fin n, Real.negMulLog (1 / (n : ℝ))
        rw [Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul]

/-- **Shannon is additive on product distributions (A3).**

    Direct calculation using `negMulLog_mul` and `Fintype.sum_prod_type`.

    TODO (~1 session):
    * Unfold `Simplex.prod`; the components are `p.prob i * q.prob j`
      after reindexing through `finProdFinEquiv`.
    * Apply `Real.negMulLog_mul` pointwise:
      `negMulLog (xy) = y · negMulLog x + x · negMulLog y`.
    * Reindex via `Equiv.sum_comp` and `Fintype.sum_prod_type` to obtain
      a sum-of-products that factorises against `∑ q = 1` and `∑ p = 1`.
    No new Mathlib lemma needed; pure algebraic rearrangement. -/
theorem shannon_additive : Additive_axiom @shannon := by
  intro m n p q
  -- Goal (after unfolding `shannon` and `Simplex.prod`): the LHS is
  -- `∑ k : Fin (m*n), negMulLog (p.prob (finProdFinEquiv.symm k).1
  --                            * q.prob (finProdFinEquiv.symm k).2)`.
  show (∑ k : Fin (m * n), Real.negMulLog
          (p.prob (finProdFinEquiv.symm k).1 *
           q.prob (finProdFinEquiv.symm k).2))
     = (∑ i, Real.negMulLog (p.prob i))
     + (∑ i, Real.negMulLog (q.prob i))
  -- Step 1: reindex the sum over `Fin (m*n)` to a sum over `Fin m × Fin n`.
  rw [show (∑ k : Fin (m * n), Real.negMulLog
            (p.prob (finProdFinEquiv.symm k).1 *
             q.prob (finProdFinEquiv.symm k).2))
        = ∑ k : Fin (m * n), (fun ij : Fin m × Fin n =>
            Real.negMulLog (p.prob ij.1 * q.prob ij.2))
              (finProdFinEquiv.symm k) from rfl]
  rw [Equiv.sum_comp finProdFinEquiv.symm
      (fun ij : Fin m × Fin n => Real.negMulLog (p.prob ij.1 * q.prob ij.2))]
  rw [Fintype.sum_prod_type]
  -- Step 2: expand `negMulLog (xy) = y · negMulLog x + x · negMulLog y` and
  -- distribute the sums.
  simp_rw [Real.negMulLog_mul, Finset.sum_add_distrib]
  -- Step 3: factor out the constant inner sums using `q.sums_to_one`
  -- and `p.sums_to_one`.
  congr 1
  · -- ∑ i, ∑ j, q.prob j * negMulLog (p.prob i) = ∑ i, negMulLog (p.prob i)
    refine Finset.sum_congr rfl (fun i _ => ?_)
    rw [← Finset.sum_mul, q.sums_to_one, one_mul]
  · -- ∑ i, ∑ j, p.prob i * negMulLog (q.prob j) = ∑ j, negMulLog (q.prob j)
    simp_rw [← Finset.mul_sum]
    rw [← Finset.sum_mul, p.sums_to_one, one_mul]

/-! ## Auxiliary structure: the uniform-entropy diagonal `f(n) = H n (uniform n)`

These lemmas package facts about the diagonal `f H n := H n (uniform n)`
for an arbitrary 4-axiom family `H`:

* `f H 1 = 0` (from additivity applied to `uniform 1 ⊗ uniform 1`);
* `f H (m·n) = f H m + f H n` (multiplicativity from additivity);
* `f H (n^k) = k · f H n` (power form by induction).

These are kernel-verified; they were originally produced as steps 1a–1b
of a (now-retracted) Shore–Johnson chain. They remain as reusable
auxiliaries (a Shannon-style family does satisfy them, and they may be
useful for the divergence-based proof). The retracted chain is not
re-introduced; see file-top revision note for context.
-/

namespace ShoreJohnson

variable {H : ∀ n, Simplex n → ℝ}

/-- The "uniform-entropy diagonal" of a Shore–Johnson family. -/
noncomputable def f (H : ∀ n, Simplex n → ℝ) (n : ℕ) (hn : 0 < n) : ℝ :=
  H n (Simplex.uniform n hn)

/-- `f H 1 = 0`: the entropy of the unique 1-state distribution is zero.
    Follows from additivity applied to `uniform 1 ⊗ uniform 1 = uniform 1`. -/
lemma f_one_zero (h_additive : Additive_axiom H) :
    f H 1 Nat.one_pos = 0 := by
  -- additivity at (m=1, n=1) gives
  --   H (1*1) (uniform 1 . prod uniform 1) = f H 1 + f H 1.
  have hadd := h_additive (Simplex.uniform 1 Nat.one_pos)
                          (Simplex.uniform 1 Nat.one_pos)
  -- Rewrite the LHS using `uniform_prod_uniform`.
  rw [Simplex.uniform_prod_uniform] at hadd
  -- `hadd : H (1*1) (uniform (1*1) _) = f H 1 + f H 1`.
  -- The index `1*1` reduces to `1` definitionally; align the two sides.
  show H 1 (Simplex.uniform 1 Nat.one_pos) = 0
  have hsimp : H (1 * 1) (Simplex.uniform (1 * 1) (Nat.mul_pos Nat.one_pos Nat.one_pos))
              = H 1 (Simplex.uniform 1 Nat.one_pos) := by
    congr 1
  rw [hsimp] at hadd
  linarith

/-- `f H` is multiplicative on `ℕ`: `f H (m·n) = f H m + f H n`. -/
lemma f_multiplicative (h_additive : Additive_axiom H)
    {m n : ℕ} (hm : 0 < m) (hn : 0 < n) :
    f H (m * n) (Nat.mul_pos hm hn) = f H m hm + f H n hn := by
  have hadd := h_additive (Simplex.uniform m hm) (Simplex.uniform n hn)
  rw [Simplex.uniform_prod_uniform] at hadd
  exact hadd

/-- **Power form on naturals.** `f H (n^k) = k · f H n`.

    Direct induction from `f_multiplicative`: `n^(k+1) = n · n^k`, so
    `f(n^(k+1)) = f(n) + f(n^k) = f(n) + k · f(n) = (k+1) · f(n)`.
    The base case `k = 0` uses `n^0 = 1` and `f_one_zero`. -/
lemma f_pow (h_additive : Additive_axiom H)
    {n : ℕ} (hn : 0 < n) :
    ∀ k : ℕ, f H (n ^ k) (pow_pos hn k) = (k : ℝ) * f H n hn
  | 0 => by
      -- n^0 = 1, so f(n^0) = f(1) = 0 = 0 · f(n)
      have h1 : f H (n ^ 0) (pow_pos hn 0) = f H 1 Nat.one_pos := by congr 1
      rw [h1, f_one_zero h_additive]
      simp
  | (k + 1) => by
      -- n^(k+1) = n * n^k. By f_multiplicative, f(n^(k+1)) = f(n) + f(n^k).
      have ih := f_pow h_additive hn k
      have hnk_pos : 0 < n ^ k := pow_pos hn k
      have hmul := f_multiplicative h_additive hn hnk_pos
      have key : f H (n ^ (k + 1)) (pow_pos hn (k + 1))
               = f H (n * n ^ k) (Nat.mul_pos hn hnk_pos) := by
        congr 1
        ring
      rw [key, hmul, ih]
      push_cast
      ring

end ShoreJohnson

/-! ## Headline theorem — RETRACTED

The former `T6c_chencov_characterisation` (Shannon entropy
characterisation from four axioms) was **retracted on 2026-05-16**
because the four-axiom statement is mathematically false: the Rényi
family `H_α(p) = (1/(1-α)) log(∑ p_i^α)` for any `α > 0, α ≠ 1`
satisfies all four axioms (continuity, symmetry, max-at-uniform,
product-additivity) but is not proportional to Shannon entropy.

The correct PT statement is on **f-divergences** rather than
entropies; see `PT/Information/T6cG1Autonomous.lean` for the
Csiszár-based formulation that the monograph
(`chapters_fr/ch05_geometry.tex` Theorem `G1_autonomous`) uses.
-/

/-! ## Cauchy's multiplicative-to-logarithmic functional equation

A continuous function `g : ℝ → ℝ` satisfying `g(xy) = g(x) + g(y)` for all
positive `x, y` is a scalar multiple of `Real.log`. This is the
multiplicative version of Cauchy's functional equation; it is the missing
classical ingredient for the Shore–Johnson / Aczél–Daróczy chain (see the
header comment on `T6c_binary_characterisation`).

The proof is short modulo `Mathlib.Topology.Instances.RealVectorSpace`:
substitute `t ↦ exp t` to turn the multiplicative equation into an additive
one, then invoke the Mathlib lemma that a continuous additive endomorphism
of `ℝ` is `ℝ`-linear (`map_real_smul`). -/

/-- **Cauchy's logarithmic functional equation.** If `g : ℝ → ℝ` is
    continuous and satisfies `g(xy) = g(x) + g(y)` for all `x, y > 0`,
    then there exists `k` such that `g x = k · log x` for every `x > 0`.

    The constant `k` equals `g(e)`. -/
theorem cauchy_log_equation
    (g : ℝ → ℝ) (hg_cont : Continuous g)
    (hg_eq : ∀ x y, 0 < x → 0 < y → g (x * y) = g x + g y) :
    ∃ k : ℝ, ∀ x, 0 < x → g x = k * Real.log x := by
  -- Substitute `t ↦ g (exp t)` to convert the multiplicative equation into
  -- Cauchy's additive equation on all of `ℝ`.
  set h : ℝ → ℝ := fun t => g (Real.exp t) with h_def
  have h_cont : Continuous h := hg_cont.comp Real.continuous_exp
  -- `g 1 = 0` from `g(1·1) = g(1) + g(1)`.
  have hg_one : g 1 = 0 := by
    have key := hg_eq 1 1 zero_lt_one zero_lt_one
    rw [mul_one] at key
    linarith
  -- `h 0 = g (exp 0) = g 1 = 0`.
  have h_zero : h 0 = 0 := by
    show g (Real.exp 0) = 0
    rw [Real.exp_zero]; exact hg_one
  -- `h` is additive: `h(s+t) = g(e^{s+t}) = g(e^s · e^t) = g(e^s) + g(e^t)`.
  have h_add : ∀ s t : ℝ, h (s + t) = h s + h t := by
    intro s t
    show g (Real.exp (s + t)) = g (Real.exp s) + g (Real.exp t)
    rw [Real.exp_add]
    exact hg_eq _ _ (Real.exp_pos s) (Real.exp_pos t)
  -- Package `h` as an `AddMonoidHom ℝ →+ ℝ`.
  let hHom : ℝ →+ ℝ :=
    { toFun := h, map_zero' := h_zero, map_add' := h_add }
  have hHom_cont : Continuous hHom := h_cont
  -- A continuous additive endomorphism of `ℝ` is `ℝ`-linear: `h c = c · h 1`.
  have h_linear : ∀ c : ℝ, h c = c * h 1 := by
    intro c
    have := map_real_smul hHom hHom_cont c (1 : ℝ)
    simp only [smul_eq_mul, mul_one] at this
    exact this
  -- Conclude with `k := h 1 = g e`.
  refine ⟨h 1, fun x hx => ?_⟩
  -- `g x = g (exp (log x)) = h (log x) = (log x) · h 1 = h 1 · log x`.
  have hg_x : g x = h (Real.log x) := by
    show g x = g (Real.exp (Real.log x))
    rw [Real.exp_log hx]
  rw [hg_x, h_linear (Real.log x), mul_comm]

/-! ## 2-state warm-up

The binary case (`n = 2`) is the form that PT G3 (Čencov uniqueness of
Fisher) actually invokes (`m1_persistence_fr.tex` line 1698:
"$\Delta_2 = \{(p, 1-p)\}$"). It is conceptually the same Shore–Johnson
chain but in dimension 2, where the simplex is a 1-parameter family. -/

/-- Parametrisation of `Simplex 2` by a single real number `p ∈ [0, 1]`. -/
noncomputable def binary (p : ℝ) (hp0 : 0 ≤ p) (hp1 : p ≤ 1) : Simplex 2 where
  prob i := if i = 0 then p else 1 - p
  nonneg i := by
    by_cases h : i = 0
    · simp [h, hp0]
    · simp [h]; linarith
  sums_to_one := by
    simp [Fin.sum_univ_two]

/-- Binary Shannon entropy `h₂(p) = -p log p - (1-p) log(1-p)`. -/
noncomputable def binaryEntropy (p : ℝ) : ℝ :=
  Real.negMulLog p + Real.negMulLog (1 - p)

/-- Shannon entropy of the binary distribution equals binary entropy. -/
lemma shannon_binary (p : ℝ) (hp0 : 0 ≤ p) (hp1 : p ≤ 1) :
    shannon (binary p hp0 hp1) = binaryEntropy p := by
  unfold shannon binary binaryEntropy
  simp [Fin.sum_univ_two]

/-! The former `T6c_binary_characterisation` (binary case as a
corollary of `T6c_chencov_characterisation`) was retracted alongside
the n-ary statement on 2026-05-16. The binary case inherits the same
issue: the Rényi family `H_α(p, 1-p) = (1/(1-α)) log(p^α + (1-p)^α)`
satisfies the four axioms but is not a scalar multiple of
`binaryEntropy`. The Csiszár-based replacement on divergences lives
in `PT/Information/T6cG1Autonomous.lean`. -/

end PT.Information
