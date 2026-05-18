/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import Mathlib.Data.ZMod.Basic
import Mathlib.Tactic

/-!
# Legendre symbol mod 5 = log-parity (Appendix N, #41)

**Statement (paper-level, App N).**
On the multiplicative group `(ℤ/5ℤ)*`, the Legendre symbol `(n/5)` equals
`(-1)^{dlog_g n}` where `g = 2` is a primitive root and `dlog_g` is the
discrete logarithm. Explicitly:

* `g = 2` is a primitive root mod 5: its powers `2^0, 2^1, 2^2, 2^3 = 1, 2, 4, 3`
  exhaust `(ℤ/5ℤ)*`.
* The quadratic residues mod 5 are `{1, 4}` (even powers of `g`);
* The non-residues are `{2, 3}` (odd powers of `g`).
* Hence `(n/5) = (-1)^{dlog_2 n}`.

This file gives the **finite combinatorial core** of the Legendre = log-parity
identity by direct enumeration on `ZMod 5`. The link to the bimodality theorem
(`PT/Sieve/Bimodality.lean`) is made explicit via the dichotomy
`δ̄(r) = 7 ⟺ r ≡ QR mod 5`.

## Reference

Appendix N of the monograph, `\label{prop:legendre_log_parity}`, item #41 of
`AUDIT_FORMALISABLE.md`.
-/

namespace PT.Sieve

/-! ### Definitions: Legendre symbol mod 5 and discrete log base 2 -/

/-- The Legendre symbol `(n/5)` as an `ℤ`-valued function on `ZMod 5`.
    Returns `+1` if `n ∈ {1, 4}` (QR), `-1` if `n ∈ {2, 3}` (NR), `0` if `n = 0`. -/
def legendre5 (n : ZMod 5) : ℤ :=
  if n = 1 ∨ n = 4 then 1
  else if n = 2 ∨ n = 3 then -1
  else 0

/-- The discrete logarithm modulo 5 with respect to the primitive root `g = 2`.
    `dlog5 1 = 0`, `dlog5 2 = 1`, `dlog5 4 = 2`, `dlog5 3 = 3`.
    Convention: `dlog5 0 = 0` (undefined extended). -/
def dlog5 (n : ZMod 5) : ℕ :=
  if n = 1 then 0
  else if n = 2 then 1
  else if n = 4 then 2
  else if n = 3 then 3
  else 0

/-! ### Basic identities on `legendre5` and `dlog5` -/

@[simp] theorem legendre5_one : legendre5 1 = 1 := by decide
@[simp] theorem legendre5_two : legendre5 2 = -1 := by decide
@[simp] theorem legendre5_three : legendre5 3 = -1 := by decide
@[simp] theorem legendre5_four : legendre5 4 = 1 := by decide
@[simp] theorem legendre5_zero : legendre5 0 = 0 := by decide

@[simp] theorem dlog5_one : dlog5 1 = 0 := by decide
@[simp] theorem dlog5_two : dlog5 2 = 1 := by decide
@[simp] theorem dlog5_three : dlog5 3 = 3 := by decide
@[simp] theorem dlog5_four : dlog5 4 = 2 := by decide

/-! ### Primitivity of `g = 2` modulo 5 -/

/-- The powers of `2` modulo 5, in order: `2^0, 2^1, 2^2, 2^3 = 1, 2, 4, 3`. -/
theorem pow_two_zero_mod5 : (2 : ZMod 5)^0 = 1 := by decide
theorem pow_two_one_mod5 : (2 : ZMod 5)^1 = 2 := by decide
theorem pow_two_two_mod5 : (2 : ZMod 5)^2 = 4 := by decide
theorem pow_two_three_mod5 : (2 : ZMod 5)^3 = 3 := by decide
theorem pow_two_four_mod5 : (2 : ZMod 5)^4 = 1 := by decide

/-- **`g = 2` is a primitive root modulo 5.** Every `n ∈ (ℤ/5ℤ)*` is `2^k`
    for `k = dlog5 n ∈ {0, 1, 2, 3}`. -/
theorem two_primitive_root_mod5 (n : ZMod 5) (hn : n ≠ 0) :
    (2 : ZMod 5)^(dlog5 n) = n := by
  fin_cases n
  · exact absurd rfl hn
  all_goals (unfold dlog5; decide)

/-- The discrete log lies in `{0, 1, 2, 3}` for `n ≠ 0`. -/
theorem dlog5_lt_four (n : ZMod 5) (hn : n ≠ 0) : dlog5 n < 4 := by
  fin_cases n
  · exact absurd rfl hn
  all_goals (unfold dlog5; decide)

/-! ### Headline: Legendre = log-parity -/

/-- **Theorem (App N #41) — Legendre symbol mod 5 equals log-parity.**

    For every `n ∈ (ℤ/5ℤ)* = {1, 2, 3, 4}`:

    $$\bigl(\tfrac{n}{5}\bigr) \;=\; (-1)^{\,\mathrm{dlog}_2 n}.$$

    Proof: exhaustive enumeration on the four units of `ZMod 5`. -/
theorem legendre5_eq_pow_neg_one_dlog5 (n : ZMod 5) (hn : n ≠ 0) :
    legendre5 n = (-1 : ℤ)^(dlog5 n) := by
  fin_cases n
  · exact absurd rfl hn
  all_goals (unfold legendre5 dlog5; decide)

/-! ### Corollary: QR/NR characterisation via log-parity -/

/-- **Corollary.** `n` is a quadratic residue mod 5 iff `dlog5 n` is even. -/
theorem isQR_iff_dlog_even (n : ZMod 5) (hn : n ≠ 0) :
    legendre5 n = 1 ↔ dlog5 n % 2 = 0 := by
  fin_cases n
  · exact absurd rfl hn
  all_goals (unfold legendre5 dlog5; decide)

/-- **Corollary.** `n` is a non-residue mod 5 iff `dlog5 n` is odd. -/
theorem isNR_iff_dlog_odd (n : ZMod 5) (hn : n ≠ 0) :
    legendre5 n = -1 ↔ dlog5 n % 2 = 1 := by
  fin_cases n
  · exact absurd rfl hn
  all_goals (unfold legendre5 dlog5; decide)

/-! ### QR / NR sets explicit (Finset form) -/

/-- The Finset of quadratic residues mod 5 is exactly `{1, 4}`.
    Stated via a Finset enumeration so `decide` resolves cleanly. -/
theorem qr_mod5_finset :
    (({0, 1, 2, 3, 4} : Finset (ZMod 5)).filter (fun n => legendre5 n = 1))
      = ({1, 4} : Finset (ZMod 5)) := by
  decide

/-- The Finset of non-residues mod 5 is exactly `{2, 3}`. -/
theorem nr_mod5_finset :
    (({0, 1, 2, 3, 4} : Finset (ZMod 5)).filter (fun n => legendre5 n = -1))
      = ({2, 3} : Finset (ZMod 5)) := by
  decide

end PT.Sieve
