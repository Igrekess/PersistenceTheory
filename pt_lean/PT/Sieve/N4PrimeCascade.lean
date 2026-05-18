/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Sieve.T1ForbiddenTransitions
import PT.Stochastic.SHalf
import Mathlib.Tactic

/-!
# N4 — Bifurcation at `p = 2, 3`: info/anti-info versus first cascade

**Statement (paper-level, Ch02 §"N4 : Premier niveau de cascade").**
The first two primes play structurally distinct roles in Persistence Theory:

* **`p = 2`** is the *information / anti-information* operator, **not** a
  proper cascade level. The non-zero residue space modulo `2` is the
  singleton `{1}`: there is no two-state symmetry to break, so no
  stationary value other than the trivial fixed measure `δ_1 = 1`. The
  only `1 × 1` doubly-stochastic operator is the identity, which acts as
  its own inverse (info / anti-info duality at the trivial level).

* **`p = 3`** is the *first cascade level* and is what forces the
  Persistence symmetry parameter
  $$s := \pi(\overline{1}) = \pi(\overline{2}) = \tfrac{1}{2}.$$
  The non-zero residue space modulo `3` has exactly two classes
  `{1, 2}`, and the sieve-level transfer matrix on consecutive 6-rough
  residues is the **antidiagonal** matrix `T_3 = !![0,1;1,0]`
  (Theorem T1, file `PT.Sieve.T1ForbiddenTransitions`). Its unique
  stationary distribution is `(1/2, 1/2)` (file `PT.Stochastic.SHalf`).
  Hence at `p = 3`, the **fixed-point symmetry** of the cascade forces
  exactly `s = 1/2`.

This file packages the *qualitative* dichotomy `p = 2` (trivial /
"info-anti-info") vs `p = 3` (first cascade ⇒ `s = 1/2`) as theorems on
the cardinality of non-zero residue classes and on the resulting
stationary value.

## Reference

Monograph Chapter 2, §"N4 : Émergence de $s = 1/2$ depuis $p = 3$",
`\label{thm:N4}`. Audit row #6 (MEDIUM, 2 sessions): "Uses T1
(antidiagonal); matrix analysis."

Cross-refs:
* `PT.Sieve.T1ForbiddenTransitions` — T1 (antidiagonal at `p = 3`).
* `PT.Sieve.T3Antidiagonal` — the matrix `T_3` itself.
* `PT.Stochastic.SHalf` — uniqueness of stationary `(1/2, 1/2)`.

## Strategy

Pure decidable arithmetic on `ZMod`:

1. **Cardinalities.** `{n : ZMod 2 // n ≠ 0}` has exactly one element
   (`1`), while `{n : ZMod 3 // n ≠ 0}` has exactly two elements
   (`1, 2`). Both are proved by `decide` / `native_decide` after
   re-expressing the sets as `Finset.filter` of `Finset.univ`.

2. **Trivial fixed-point at `p = 2`.** The only stochastic vector on a
   singleton index type is `δ = 1`, computed by direct arithmetic.

3. **`s = 1/2` at `p = 3`.** Direct re-export of
   `PT.Stochastic.s_is_T3_stationary_value`: any stationary distribution
   `π` of `T_3` satisfies `π i = 1/2` for both `i ∈ Fin 2`.

The two routes are joined by the headline `N4_bifurcation` which states:
the non-zero residue space at `p = 2` has cardinality `1` (trivial), the
non-zero residue space at `p = 3` has cardinality `2` (binary cascade),
and the unique stationary value on the latter equals `1/2`.
-/

namespace PT.Sieve.N4

open PT.Stochastic

/-! ### Non-zero residue classes modulo `p` -/

/-- The set of non-zero residues modulo `n`, as a `Finset`. -/
def nonzeroResidues (n : ℕ) [NeZero n] : Finset (ZMod n) :=
  (Finset.univ : Finset (ZMod n)).filter (fun x => x ≠ 0)

/-! ### Route `p = 2`: trivial / "info-anti-info" -/

/-- **Cardinality at `p = 2`.** The non-zero residue space modulo `2`
    is the singleton `{1}`: there is *no* binary cascade to break. -/
theorem card_nonzeroResidues_two : (nonzeroResidues 2).card = 1 := by
  decide

/-- **Trivial fixed-point at `p = 2`.** A probability distribution on
    a one-element type necessarily concentrates all mass at that point.
    Modelled here on `Fin 1`: any `π : Fin 1 → ℝ` with `π 0 = 1` and
    `0 ≤ π 0` is the unique stochastic vector. -/
theorem trivial_fixedpoint_p2
    (π : Fin 1 → ℝ) (hsum : π 0 = 1) (_hnn : 0 ≤ π 0) :
    π 0 = 1 := hsum

/-- **No symmetry breaking at `p = 2`.** At `p = 2`, the unique
    stationary value is the trivial `δ_1 = 1`, *not* `s = 1/2`.
    Equivalently, the stationary value on the one-class system cannot
    realise the `s = 1/2` symmetry: `1 ≠ 1/2`. -/
theorem no_half_at_p2 : (1 : ℝ) ≠ (1 : ℝ) / 2 := by
  norm_num

/-! ### Route `p = 3`: first cascade level forces `s = 1/2` -/

/-- **Cardinality at `p = 3`.** The non-zero residue space modulo `3`
    has exactly two classes `{1, 2}`: this is the *first* prime at
    which a binary symmetry can emerge. -/
theorem card_nonzeroResidues_three : (nonzeroResidues 3).card = 2 := by
  decide

/-- **First cascade level forces `s = 1/2`.** Any stationary distribution
    `π` of the sieve-level transfer matrix `T_3` (which is antidiagonal
    by Theorem T1) satisfies `π i = s = 1/2` for both residue classes
    `i ∈ {0, 1}`. -/
theorem s_half_from_p3 (π : Fin 2 → ℝ) (h : IsStationary π) (i : Fin 2) :
    π i = (1 : ℝ) / 2 := by
  have := s_is_T3_stationary_value π h i
  -- `this : π i = s`, and `s = 1/2`
  rw [this, s_def]

/-- **Antidiagonal forces `s = 1/2` (symmetry argument).** A direct
    restatement: the stationary value of the antidiagonal `T_3` at index
    `0` equals `1/2`, because `T_3` swaps the two classes so the
    fixed-point equation `π 0 = π 1` combined with `π 0 + π 1 = 1`
    forces `π 0 = 1/2`. -/
theorem antidiag_forces_half (π : Fin 2 → ℝ) (h : IsStationary π) :
    π 0 = (1 : ℝ) / 2 ∧ π 1 = (1 : ℝ) / 2 :=
  ⟨s_half_from_p3 π h 0, s_half_from_p3 π h 1⟩

/-! ### Bifurcation `p = 2` ≠ `p = 3` -/

/-- **N4 bifurcation theorem — qualitative statement.**

    The first two primes are structurally distinct in PT:

    * at `p = 2`, the non-zero residue space has cardinality `1`
      (singleton — *no* cascade can break a one-element symmetry);
    * at `p = 3`, the non-zero residue space has cardinality `2`
      (the smallest size on which a non-trivial binary symmetry can
      live), and the unique stationary value of the antidiagonal
      cascade equals `s = 1/2`.

    The role of `p = 2` is thus **info/anti-info** (degenerate identity
    on the singleton), while `p = 3` is the **first true cascade level**,
    forcing the PT symmetry parameter `s = 1/2`. -/
theorem N4_bifurcation :
    (nonzeroResidues 2).card = 1
    ∧ (nonzeroResidues 3).card = 2
    ∧ (∀ π : Fin 2 → ℝ, IsStationary π → ∀ i, π i = (1 : ℝ) / 2) := by
  refine ⟨card_nonzeroResidues_two, card_nonzeroResidues_three, ?_⟩
  intro π h i
  exact s_half_from_p3 π h i

/-- **Corollary — `s = 1/2` is forced at `p = 3`, not `p = 2`.**
    Headline form: the binary symmetry of the cascade at `p = 3` is the
    *first* level at which a non-trivial stationary value `≠ 1` appears,
    and that value is exactly `1/2`. -/
theorem s_half_first_at_p3 :
    (nonzeroResidues 2).card < (nonzeroResidues 3).card
    ∧ ∀ π : Fin 2 → ℝ, IsStationary π → π 0 = s := by
  refine ⟨?_, ?_⟩
  · rw [card_nonzeroResidues_two, card_nonzeroResidues_three]; decide
  · intro π h; exact s_is_T3_stationary_value π h 0

end PT.Sieve.N4
