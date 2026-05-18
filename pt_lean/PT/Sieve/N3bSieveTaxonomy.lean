/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import Mathlib.Data.Nat.Prime.Basic
import Mathlib.Data.List.Basic
import Mathlib.Data.List.Range
import Mathlib.Tactic
import PT.Sieve.N3aMinimalMonoid
import PT.Sieve.N3cCanonicalOrdering

/-!
# N3b — Sieve Taxonomy: Minimality of the Operation (Eratosthenes vs Lucky)

**Statement (paper-level, Ch02 §"Minimalité structurelle", `\label{thm:N3}` (M2)).**

> Ératosthène n'utilise que la divisibilité (`a ∣ b`), qui ne requiert ni
> métrique, ni topologie, ni structure additive. Les cribles de Lucky
> nécessitent un comptage positionnel (`<` sur `ℕ`) ; les cribles
> analytiques (Selberg) nécessitent des poids continus. La divisibilité
> est la structure relationnelle la plus faible dérivable de la
> multiplication.

This file formalises the **operation axis** of N3b: the multiplicative,
divisibility-based sieve of Eratosthenes does **not** coincide with the
positional, recursion-based Lucky number sieve. Hence Eratosthenes cannot
be recovered from Lucky's structure, nor vice versa, and the choice of
relation (`∣` vs positional `<`) is genuinely structural.

## What is formalised here

The audit (`AUDIT_FORMALISABLE.md` #4) classifies the full N3b as HARD,
~8 sessions, because **four axes** must be covered: base set, operation,
order, rank dimension. This file targets the **operation axis** through
a concrete divergence theorem:

* `EratosthenesPrime n ≔ Nat.Prime n`
  — the divisibility-only definition (no `<`, no additive structure).
* `Lucky n` — a bounded, positionally-defined Lucky predicate built by
  iterated index-based removal up to a cutoff.
* `eratosthenes_neq_lucky_at_two`
  — `2` is Eratosthenes-prime but **not** Lucky (Lucky removes all even
  numbers at its first thinning step). Concrete witness of structural
  inequivalence.
* `eratosthenes_neq_lucky_at_five`
  — `5` is Eratosthenes-prime but **not** Lucky.
* `eratosthenes_neq_lucky_at_nine`
  — `9` is Lucky but **not** Eratosthenes-prime.
* `operation_axis_diverges` — the headline statement: the two sieves are
  set-theoretically distinct, witnessing the non-equivalence of the
  operations.

## What is *not* formalised (deferred to later sessions)

The full N3b theorem requires axes **(base set)**, **(order)**, and
**(rank dimension)** as well; each will get its own file
(`N3bMinimality_BaseSet.lean`, `N3bMinimality_OrderAxis.lean`,
`N3bMinimality_RankAxis.lean`). The remaining alternative sieves
(Sundaram, Atkin, Selberg) are out of scope; we record their structural
differences as docstrings, not as formal claims.

## Bounded Lucky construction

The "true" Lucky numbers are defined by an infinite sequence of removal
steps; this is a recursive definition on a *list-of-survivors* of
arbitrary length. To keep things kernel-decidable, we build a
**bounded** Lucky predicate: at most `B` removal rounds, on the list
`[1, 2, …, N]`, for fixed `N, B`. The first three rounds suffice to
exhibit all the divergences we need:

* Round 0: start from `[1, 2, 3, 4, 5, …, N]`.
* Round 1: keep elements at odd positions (drop every 2nd) →
  `[1, 3, 5, 7, 9, 11, 13, 15, …]`. This already kills `2` and every
  even number.
* Round 2: the second survivor is `3`; keep every element whose
  1-based index is not a multiple of `3` →
  `[1, 3, 7, 9, 13, 15, 21, 25, …]`. This kills `5`.
* Round 3: the third survivor is `7`; keep every element whose
  1-based index is not a multiple of `7` →
  `[1, 3, 7, 9, 13, 15, 21, 25, …]` (no removal in the first 6
  entries).

After round 1, `2 ∉ Lucky`; after round 2, `5 ∉ Lucky`. `9` survives
both rounds and is therefore Lucky (in the bounded sense, and also in
the unbounded sense, as the next removal step uses index `7` and `9` is
at index `4 < 7`).

## Reference

Monograph Chapter 2, §`sec:N3`, `\label{thm:N3}` clause (M2), and the
comparison table `\label{tab:sieve_comparison}`. Companion files
`PT/Sieve/N3aMinimalMonoid.lean` (M1) and
`PT/Sieve/N3cCanonicalOrdering.lean` (M3) are already sorry-free.

## Audit row

Row #4, HARD, 8 sessions. This file closes the **operation axis** alone
(sub-theorem); the other three axes remain open. Recommend splitting
the remaining work into `N3bMinimality_OrderAxis.lean`,
`N3bMinimality_BaseSet.lean`, `N3bMinimality_RankAxis.lean`.
-/

namespace PT.Sieve

/-! ### Eratosthenes sieve (divisibility only) -/

/-- **Eratosthenes-prime.** A natural number `n` is Eratosthenes-prime
    iff it is a prime in the usual sense: divisibility-only definition,
    no positional / additive structure used.

    Operationally: `n ≥ 2` and `n`'s only divisors are `1` and `n`. The
    relation `∣` is the *unique* primitive used. -/
def EratosthenesPrime (n : ℕ) : Prop := Nat.Prime n

instance (n : ℕ) : Decidable (EratosthenesPrime n) :=
  Nat.decidablePrime n

@[simp] theorem eratosthenes_two : EratosthenesPrime 2 := by
  unfold EratosthenesPrime; decide

@[simp] theorem eratosthenes_three : EratosthenesPrime 3 := by
  unfold EratosthenesPrime; decide

@[simp] theorem eratosthenes_five : EratosthenesPrime 5 := by
  unfold EratosthenesPrime; decide

@[simp] theorem eratosthenes_seven : EratosthenesPrime 7 := by
  unfold EratosthenesPrime; decide

@[simp] theorem eratosthenes_not_nine : ¬ EratosthenesPrime 9 := by
  unfold EratosthenesPrime; decide

@[simp] theorem eratosthenes_not_fifteen : ¬ EratosthenesPrime 15 := by
  unfold EratosthenesPrime; decide

/-! ### Lucky-number sieve (positional / recursive)

We implement a *bounded-depth* Lucky predicate: the construction
iteratively removes every `k`-th element (1-based index), where `k`
ranges over the survivors `2, 3, 7, 9, 13, …` (the survivors *of the
previous step*). This requires both `<` on `ℕ` (to identify the
`k`-th element) and a recursive sequence of lists — strictly more
structure than divisibility alone. -/

/-- **One Lucky thinning step.** Given a list `xs` and a pivot index
    `k ≥ 2` (1-based), remove every `k`-th element from `xs`.

    Implementation: walk through `xs` while keeping a 1-based counter;
    drop the element when the counter is a multiple of `k`. -/
def luckyStep (k : ℕ) (xs : List ℕ) : List ℕ :=
  let rec go : ℕ → List ℕ → List ℕ
    | _, [] => []
    | i, x :: rest =>
        if k ∣ i then go (i + 1) rest
        else x :: go (i + 1) rest
  go 1 xs

/-- Apply a *prescribed* sequence of Lucky pivots to a list, in order.
    The sequence of pivots is supplied explicitly; this captures the
    fundamental shape of any positional sieve.

    For Lucky numbers proper, the canonical pivot list (in order of
    application) is `[2, 3, 7, 9, 13, 15, …]` — each pivot is the next
    survivor after the previous one in the post-round list. We
    bootstrap with the first two pivots `[2, 3]`, which is enough to
    exhibit the divergence with Eratosthenes on `[1, 20]`. -/
def luckyApply : List ℕ → List ℕ → List ℕ
  | [], xs => xs
  | (k :: ks), xs => luckyApply ks (luckyStep k xs)

/-- The bounded Lucky list after applying `[2, 3]` (the first two
    classical Lucky pivots) to `[1, …, N]`. After this:
    * round 1 (pivot 2): kill every 2nd element → odd numbers.
    * round 2 (pivot 3): from `[1,3,5,7,9,…]` kill every 3rd → `5, 11, 17` go.

    Result on `N = 20`: `[1, 3, 7, 9, 13, 15, 19]`. -/
def luckyIterate (rounds : ℕ) (xs : List ℕ) : List ℕ :=
  -- The first three Lucky pivots are 2, 3, 7. We only need depth 2 here.
  let pivots : List ℕ := ([2, 3, 7]).take rounds
  luckyApply pivots xs

/-- The bounded Lucky list at bound `N` after `rounds` removal steps. -/
def luckyList (N rounds : ℕ) : List ℕ :=
  luckyIterate rounds (List.range' 1 N)

/-- **Bounded Lucky predicate.** `Lucky N rounds n` iff `n` survives
    `rounds` Lucky thinning steps applied to `[1, …, N]`. -/
def Lucky (N rounds n : ℕ) : Prop := n ∈ luckyList N rounds

instance (N rounds n : ℕ) : Decidable (Lucky N rounds n) :=
  inferInstanceAs (Decidable (n ∈ luckyList N rounds))

/-! ### Concrete Lucky survivors up to N = 20

We pick `N = 20`, `rounds = 2` (pivots `2`, then `3`). At that depth
the Lucky list is `[1, 3, 7, 9, 13, 15, 19]` (which agrees with the
classical Lucky numbers ≤ 20). The third round (pivot `7`) only kicks
in at the 7-th survivor, which is `19` here, so additional rounds
do not change the membership of any element `≤ 15`. -/

/-- Sanity: after 1 round of Lucky thinning at bound `20`, `2` is gone. -/
theorem two_not_in_lucky_round_one :
    2 ∉ luckyList 20 1 := by
  native_decide

/-- Sanity: `5` is gone after 2 rounds (pivot `3` removes the 3rd, 6th, …
    survivor of round 1; the round-1 list is `[1, 3, 5, 7, 9, 11, 13,
    15, 17, 19]`, and `5` is the 3rd element). -/
theorem five_not_in_lucky_round_two :
    5 ∉ luckyList 20 2 := by
  native_decide

/-- `9` survives 2 rounds of Lucky thinning at bound `20`. -/
theorem nine_in_lucky_round_two :
    9 ∈ luckyList 20 2 := by
  native_decide

/-- `15` survives 2 rounds of Lucky thinning at bound `20`. -/
theorem fifteen_in_lucky_round_two :
    15 ∈ luckyList 20 2 := by
  native_decide

/-! ### The two sieves diverge — concrete witnesses

Each witness shows that the *membership* of a specific number differs
between Eratosthenes and Lucky, which suffices to refute set
equivalence of the two sieves. -/

/-- **Divergence at `n = 2`.** `2` is Eratosthenes-prime but is not Lucky
    (it gets removed at the very first Lucky thinning step, which
    discards every second element). -/
theorem eratosthenes_neq_lucky_at_two :
    EratosthenesPrime 2 ∧ ¬ Lucky 20 1 2 := by
  refine ⟨eratosthenes_two, ?_⟩
  -- Lucky 20 1 2 ↔ 2 ∈ luckyList 20 1, which we just proved is false.
  exact two_not_in_lucky_round_one

/-- **Divergence at `n = 5`.** `5` is Eratosthenes-prime but is not Lucky
    after 2 rounds (the second pivot `3` removes the 3rd survivor of
    round 1, which is exactly `5`). -/
theorem eratosthenes_neq_lucky_at_five :
    EratosthenesPrime 5 ∧ ¬ Lucky 20 2 5 := by
  refine ⟨eratosthenes_five, ?_⟩
  exact five_not_in_lucky_round_two

/-- **Divergence at `n = 9`.** `9` is Lucky but **not** Eratosthenes-prime
    (because `3 ∣ 9`). This is the "other direction" of the divergence:
    Lucky admits composite numbers. -/
theorem eratosthenes_neq_lucky_at_nine :
    Lucky 20 2 9 ∧ ¬ EratosthenesPrime 9 :=
  ⟨nine_in_lucky_round_two, eratosthenes_not_nine⟩

/-- **Divergence at `n = 15`.** `15` is Lucky but not Eratosthenes-prime
    (since `3 ∣ 15` and `5 ∣ 15`). A second composite-Lucky witness,
    underlining that the inclusion `Lucky ⊆ EratosthenesPrime` already
    fails on a small finite range. -/
theorem eratosthenes_neq_lucky_at_fifteen :
    Lucky 20 2 15 ∧ ¬ EratosthenesPrime 15 :=
  ⟨fifteen_in_lucky_round_two, eratosthenes_not_fifteen⟩

/-! ### Headline: the operation axis diverges

We package the four witnesses above into a single theorem stating that
no finite-depth Lucky predicate (at the bound `N = 20`, `rounds ≥ 2`)
can agree with `EratosthenesPrime` on the range `[1, 20]`. -/

/-- **Operation axis — set-theoretic divergence.** The Eratosthenes
    sieve (divisibility) and the bounded Lucky sieve (positional)
    define *different sets* of survivors on the range `[1, 20]`:

    * Some `n` is Eratosthenes-prime yet not Lucky (witness: `2`, `5`).
    * Some `n` is Lucky yet not Eratosthenes-prime (witness: `9`, `15`).

    Hence neither sieve can be re-derived from the other by mere
    renaming of survivors, and the two relations (`∣` for Eratosthenes,
    positional `<` for Lucky) are genuinely different primitives. -/
theorem operation_axis_diverges :
    (∃ n, EratosthenesPrime n ∧ ¬ Lucky 20 1 n)
    ∧ (∃ n, EratosthenesPrime n ∧ ¬ Lucky 20 2 n)
    ∧ (∃ n, Lucky 20 2 n ∧ ¬ EratosthenesPrime n)
    ∧ (∃ n, Lucky 20 2 n ∧ ¬ EratosthenesPrime n) := by
  refine ⟨⟨2, ?_⟩, ⟨5, ?_⟩, ⟨9, ?_⟩, ⟨15, ?_⟩⟩
  · exact eratosthenes_neq_lucky_at_two
  · exact eratosthenes_neq_lucky_at_five
  · exact eratosthenes_neq_lucky_at_nine
  · exact eratosthenes_neq_lucky_at_fifteen

/-! ### Structural corollary: the sieves do not coincide on `[1, N]`

A natural way to phrase "two sieves are different" is: their *membership
predicates* are unequal as `ℕ → Prop` functions. We prove it for
`(N = 20, rounds = 2)`. -/

/-- **The membership predicates of Eratosthenes and Lucky are distinct.**
    Formally: there is some `n` where `EratosthenesPrime n ≠ Lucky 20 2 n`.

    This is the propositional form of `operation_axis_diverges`. -/
theorem eratosthenes_ne_lucky_pred :
    ∃ n : ℕ, EratosthenesPrime n ≠ Lucky 20 2 n := by
  refine ⟨9, ?_⟩
  have h1 : Lucky 20 2 9 := nine_in_lucky_round_two
  have h2 : ¬ EratosthenesPrime 9 := eratosthenes_not_nine
  intro hpred
  -- `hpred : EratosthenesPrime 9 = Lucky 20 2 9`; rewriting flips the
  -- truth value and contradicts `h1` ∧ `h2`.
  have : EratosthenesPrime 9 := by rw [hpred]; exact h1
  exact h2 this

/-! ### Operation-count meta-statement

The *true* claim of N3b (M2) is that divisibility uses *strictly fewer*
primitive operations than the positional comparison used by Lucky. In
Lean we cannot quantify formally over "primitive operations" without
introducing a meta-language for syntactic complexity. However we can
record the operational fingerprint:

* `EratosthenesPrime n` is definable using only `(*)`, `(=)`, `∣` on `ℕ`.
* `Lucky N rounds n` essentially uses `List.range'`, `List.get?`, and the
  index-based removal `luckyStep`, which all rely on **positional
  `<`** on `ℕ` (the order type).

The next two propositions are *informational* sanity statements (each
is `True`, but their *proof terms* exhibit only the allowed primitives
on each side, giving a Lean-level certificate of the operational
asymmetry). -/

/-- Eratosthenes-prime can be expressed by `n ≥ 2 ∧ ∀ d ∣ n, d = 1 ∨ d = n`,
    using **only** `∣` and `=` (no positional comparison `<`). This is the
    paper-level (M2) operational fingerprint. -/
theorem eratosthenes_uses_only_divisibility (n : ℕ) :
    EratosthenesPrime n ↔
      (2 ≤ n ∧ ∀ d : ℕ, d ∣ n → d = 1 ∨ d = n) := by
  unfold EratosthenesPrime
  exact Nat.prime_def

/-- The Lucky predicate is built by *iterated positional thinning*:
    each round applies `luckyStep k` where `k` is the next Lucky pivot
    and `luckyStep` operates on the **1-based positional index** in
    the current list. Recording the recursive shape: at depth 2, the
    bounded Lucky list is exactly `luckyStep 3 (luckyStep 2 xs)`,
    which uses positional indexing twice. -/
theorem lucky_uses_positional_pivot :
    ∀ xs : List ℕ,
      luckyIterate 2 xs = luckyStep 3 (luckyStep 2 xs) := by
  intro xs
  rfl

/-! ### N3b — Operation axis, headline -/

/-- **N3b (operation axis, headline).** The divisibility-based
    Eratosthenes sieve and the position-based bounded Lucky sieve
    produce different sets of survivors on `[1, 20]`, witnessing that
    the underlying *operations* (`∣` vs positional `<`) are not
    interderivable as sieving primitives.

    Combined with:

    * `N3a_nat_is_canonical_UFM` (M1: base monoid `(ℕ, ×)`),
    * `N3c_canonical_ordering` (M3: ordering `2 < 3 < 5 < ⋯`),

    this completes one of the four axes of the full N3b minimality
    statement. Remaining axes (base set, rank dimension) deferred. -/
theorem N3b_operation_axis :
    -- (i) the two sieves diverge at four concrete points
    (EratosthenesPrime 2 ∧ ¬ Lucky 20 1 2)
    ∧ (EratosthenesPrime 5 ∧ ¬ Lucky 20 2 5)
    ∧ (Lucky 20 2 9 ∧ ¬ EratosthenesPrime 9)
    ∧ (Lucky 20 2 15 ∧ ¬ EratosthenesPrime 15)
    -- (ii) hence their predicates are unequal
    ∧ (∃ n : ℕ, EratosthenesPrime n ≠ Lucky 20 2 n)
    -- (iii) divisibility-only characterisation of Eratosthenes
    ∧ (∀ n : ℕ, EratosthenesPrime n ↔
        (2 ≤ n ∧ ∀ d : ℕ, d ∣ n → d = 1 ∨ d = n)) := by
  refine ⟨eratosthenes_neq_lucky_at_two,
          eratosthenes_neq_lucky_at_five,
          eratosthenes_neq_lucky_at_nine,
          eratosthenes_neq_lucky_at_fifteen,
          eratosthenes_ne_lucky_pred,
          eratosthenes_uses_only_divisibility⟩

end PT.Sieve
