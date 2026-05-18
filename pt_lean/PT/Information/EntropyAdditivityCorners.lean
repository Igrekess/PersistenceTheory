/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Information.GFTIdentity
import PT.Information.GFTSpecialisations
import PT.Information.BekensteinExtensions
import PT.Information.EntropyBoundsTight
import PT.Information.KLAdditivityProduct
import Mathlib.Tactic

/-!
# Entropy / GFT additivity at the two corners (Ch04 extension)

This file proves **additivity** of the Shannon entropy and of the GFT
budget `D_KL + H` at the two extremal corners of the partition:

* **Delta corner** (`P = δ_{r₀}`): all mass on one residue. `H = 0`,
  `D_KL = log m`.
* **Uniform corner** (`P = U_m = 1/m`): mass spread evenly. `H = log m`,
  `D_KL = 0`.

The algebraic backbone is the elementary identity
`log m + log n = log (m · n)` (already exported by `KLAdditivityProduct`
as `log_mul_add`).

Combined with the corner specialisations of `shannonH` and `klToUniform`,
this gives **additivity** of the GFT total `D_KL + H` over independent
moduli `(m, n)` at both corners:

  `(D_KL + H)(δ_m) + (D_KL + H)(δ_n) = log m + log n = log (m · n)`,
  `(D_KL + H)(U_m) + (D_KL + H)(U_n) = log m + log n = log (m · n)`.

## Numerical specialisations

The PT-canonical moduli `(m, n) ∈ {(2, 8), (8, 30)}` are also given as
concrete instances.

## Reference

Monograph Chapter 4 §"Additivité aux coins", follow-up to
`KLAdditivityProduct` and `EntropyBoundsTight`.
-/

namespace PT.Information

open Real Finset

/-! ### Shannon-entropy additivity at the delta corner -/

/-- **Delta-corner Shannon additivity.** Both entropies vanish, so the
    sum is `0 = 0 + 0`. -/
theorem shannonH_delta_add
    {α β : Type*} [DecidableEq α] [DecidableEq β]
    (s : Finset α) (t : Finset β)
    (r₀ : α) (hr : r₀ ∈ s) (s₀ : β) (hs : s₀ ∈ t) :
    shannonH s (deltaDist r₀) + shannonH t (deltaDist s₀) = 0 := by
  rw [shannonH_delta s r₀ hr, shannonH_delta t s₀ hs]
  ring

/-! ### Shannon-entropy additivity at the uniform corner -/

/-- **Uniform-corner Shannon additivity.**
    `H(U_m) + H(U_n) = log m + log n = log (m · n)`. -/
theorem shannonH_uniform_add
    {α β : Type*} (s : Finset α) (t : Finset β)
    (m n : ℝ) (hm : 0 < m) (hn : 0 < n)
    (hcard_s : (s.card : ℝ) = m) (hcard_t : (t.card : ℝ) = n) :
    shannonH s (fun _ => (1 : ℝ) / m) + shannonH t (fun _ => (1 : ℝ) / n)
      = Real.log (m * n) := by
  rw [shannonH_uniform_eq_log s m hm hcard_s,
      shannonH_uniform_eq_log t n hn hcard_t]
  exact (log_mul_add m n hm hn).symm

/-! ### GFT total additivity at the delta corner -/

/-- **Delta-corner GFT additivity.** Combining the KL delta-additivity
    with the vanishing of Shannon at deltas:

      `(D_KL(δ_m) + H(δ_m)) + (D_KL(δ_n) + H(δ_n))
         = log m + log n = log (m · n)`. -/
theorem GFT_delta_add
    {α β : Type*} [DecidableEq α] [DecidableEq β]
    (s : Finset α) (t : Finset β) (m n : ℝ)
    (hm : 0 < m) (hn : 0 < n)
    (r₀ : α) (hr : r₀ ∈ s) (s₀ : β) (hs : s₀ ∈ t) :
    (klToUniform s m (deltaDist r₀) + shannonH s (deltaDist r₀))
      + (klToUniform t n (deltaDist s₀) + shannonH t (deltaDist s₀))
      = Real.log (m * n) := by
  rw [GFT_at_delta s m hm r₀ hr, GFT_at_delta t n hn s₀ hs]
  exact (log_mul_add m n hm hn).symm

/-! ### GFT total additivity at the uniform corner -/

/-- **Uniform-corner GFT additivity.** With KL = 0 and `H = log` at the
    uniform corner:

      `(D_KL(U_m) + H(U_m)) + (D_KL(U_n) + H(U_n))
         = 0 + log m + 0 + log n = log (m · n)`. -/
theorem GFT_uniform_add
    {α β : Type*} (s : Finset α) (t : Finset β)
    (m n : ℝ) (hm : 0 < m) (hn : 0 < n)
    (hcard_s : (s.card : ℝ) = m) (hcard_t : (t.card : ℝ) = n) :
    (klToUniform s m (fun _ => (1 : ℝ) / m)
        + shannonH s (fun _ => (1 : ℝ) / m))
      + (klToUniform t n (fun _ => (1 : ℝ) / n)
        + shannonH t (fun _ => (1 : ℝ) / n))
      = Real.log (m * n) := by
  rw [GFT_at_uniform s m hm hcard_s, GFT_at_uniform t n hn hcard_t]
  exact (log_mul_add m n hm hn).symm

/-! ### Total budget conservation -/

/-- **Total budget conservation.** Both corners obey the *same* additive
    identity `log m + log n = log (m · n)`. The internal allocation
    between `D_KL` and `H` differs:

    * Delta corner: `(log m, 0) + (log n, 0) = (log(m·n), 0)`.
    * Uniform corner: `(0, log m) + (0, log n) = (0, log(m·n))`.

    But the *total* `D_KL + H` is invariant: `log (m · n)` in both
    cases. -/
theorem GFT_corner_total_conservation
    {α β : Type*} [DecidableEq α] [DecidableEq β]
    (s : Finset α) (t : Finset β)
    (m n : ℝ) (hm : 0 < m) (hn : 0 < n)
    (hcard_s : (s.card : ℝ) = m) (hcard_t : (t.card : ℝ) = n)
    (r₀ : α) (hr : r₀ ∈ s) (s₀ : β) (hs : s₀ ∈ t) :
    -- Delta total = Uniform total
    (klToUniform s m (deltaDist r₀) + shannonH s (deltaDist r₀))
      + (klToUniform t n (deltaDist s₀) + shannonH t (deltaDist s₀))
    = (klToUniform s m (fun _ => (1 : ℝ) / m)
        + shannonH s (fun _ => (1 : ℝ) / m))
      + (klToUniform t n (fun _ => (1 : ℝ) / n)
        + shannonH t (fun _ => (1 : ℝ) / n)) := by
  rw [GFT_delta_add s t m n hm hn r₀ hr s₀ hs,
      GFT_uniform_add s t m n hm hn hcard_s hcard_t]

/-! ### PT-canonical instances `(m, n) ∈ {(2, 8), (8, 30)}` -/

/-- **Shannon uniform additivity at `(2, 8)`:**
    `H(U_2) + H(U_8) = log 16`. -/
theorem shannonH_uniform_2_plus_8 :
    shannonH (Finset.univ : Finset (Fin 2)) (fun _ => (1 : ℝ) / 2)
      + shannonH (Finset.univ : Finset (Fin 8)) (fun _ => (1 : ℝ) / 8)
      = Real.log 16 := by
  have h := shannonH_uniform_add
    (Finset.univ : Finset (Fin 2)) (Finset.univ : Finset (Fin 8))
    2 8 (by norm_num) (by norm_num) (by simp) (by simp)
  convert h using 2
  norm_num

/-- **Shannon uniform additivity at `(8, 30)`:**
    `H(U_8) + H(U_30) = log 240`. -/
theorem shannonH_uniform_8_plus_30 :
    shannonH (Finset.univ : Finset (Fin 8)) (fun _ => (1 : ℝ) / 8)
      + shannonH (Finset.univ : Finset (Fin 30)) (fun _ => (1 : ℝ) / 30)
      = Real.log 240 := by
  have h := shannonH_uniform_add
    (Finset.univ : Finset (Fin 8)) (Finset.univ : Finset (Fin 30))
    8 30 (by norm_num) (by norm_num) (by simp) (by simp)
  convert h using 2
  norm_num

/-- **GFT delta additivity at `(2, 8)`:** `(log 2 + 0) + (log 8 + 0) = log 16`. -/
theorem GFT_delta_2_plus_8 :
    (klToUniform (Finset.univ : Finset (Fin 2)) 2 (deltaDist 0)
        + shannonH (Finset.univ : Finset (Fin 2)) (deltaDist 0))
      + (klToUniform (Finset.univ : Finset (Fin 8)) 8 (deltaDist 0)
        + shannonH (Finset.univ : Finset (Fin 8)) (deltaDist 0))
      = Real.log 16 := by
  have h := GFT_delta_add
    (Finset.univ : Finset (Fin 2)) (Finset.univ : Finset (Fin 8))
    2 8 (by norm_num) (by norm_num) 0 (by decide) 0 (by decide)
  convert h using 2
  norm_num

/-- **GFT delta additivity at `(8, 30)`:** `(log 8 + 0) + (log 30 + 0) = log 240`. -/
theorem GFT_delta_8_plus_30 :
    (klToUniform (Finset.univ : Finset (Fin 8)) 8 (deltaDist 0)
        + shannonH (Finset.univ : Finset (Fin 8)) (deltaDist 0))
      + (klToUniform (Finset.univ : Finset (Fin 30)) 30 (deltaDist 0)
        + shannonH (Finset.univ : Finset (Fin 30)) (deltaDist 0))
      = Real.log 240 := by
  have h := GFT_delta_add
    (Finset.univ : Finset (Fin 8)) (Finset.univ : Finset (Fin 30))
    8 30 (by norm_num) (by norm_num) 0 (by decide) 0 (by decide)
  convert h using 2
  norm_num

/-- **GFT uniform additivity at `(2, 8)`:** `(0 + log 2) + (0 + log 8) = log 16`. -/
theorem GFT_uniform_2_plus_8 :
    (klToUniform (Finset.univ : Finset (Fin 2)) 2 (fun _ => (1 : ℝ) / 2)
        + shannonH (Finset.univ : Finset (Fin 2)) (fun _ => (1 : ℝ) / 2))
      + (klToUniform (Finset.univ : Finset (Fin 8)) 8 (fun _ => (1 : ℝ) / 8)
        + shannonH (Finset.univ : Finset (Fin 8)) (fun _ => (1 : ℝ) / 8))
      = Real.log 16 := by
  have h := GFT_uniform_add
    (Finset.univ : Finset (Fin 2)) (Finset.univ : Finset (Fin 8))
    2 8 (by norm_num) (by norm_num) (by simp) (by simp)
  convert h using 2
  norm_num

/-- **GFT uniform additivity at `(8, 30)`:** `(0 + log 8) + (0 + log 30) = log 240`. -/
theorem GFT_uniform_8_plus_30 :
    (klToUniform (Finset.univ : Finset (Fin 8)) 8 (fun _ => (1 : ℝ) / 8)
        + shannonH (Finset.univ : Finset (Fin 8)) (fun _ => (1 : ℝ) / 8))
      + (klToUniform (Finset.univ : Finset (Fin 30)) 30 (fun _ => (1 : ℝ) / 30)
        + shannonH (Finset.univ : Finset (Fin 30)) (fun _ => (1 : ℝ) / 30))
      = Real.log 240 := by
  have h := GFT_uniform_add
    (Finset.univ : Finset (Fin 8)) (Finset.univ : Finset (Fin 30))
    8 30 (by norm_num) (by norm_num) (by simp) (by simp)
  convert h using 2
  norm_num

/-! ### Headline -/

/-- **Headline (entropy / GFT additivity at the two corners).**
    The GFT total `D_KL + H` is additive over independent moduli at both
    extremal corners, with the *same* algebraic backbone
    `log m + log n = log (m · n)`:

    * **Delta corner** (`(log m, 0)`): trivial Shannon additivity
      (`0 + 0`), KL carries the entire `log (m · n)`.
    * **Uniform corner** (`(0, log m)`): KL vanishes, entropy carries
      the entire `log (m · n)`.

    Specialised instances at PT-canonical moduli `(2, 8)` and `(8, 30)`
    give `log 16` and `log 240` respectively. -/
theorem entropy_additivity_corners_summary
    {α β : Type*} [DecidableEq α] [DecidableEq β]
    (s : Finset α) (t : Finset β)
    (m n : ℝ) (hm : 0 < m) (hn : 0 < n)
    (hcard_s : (s.card : ℝ) = m) (hcard_t : (t.card : ℝ) = n)
    (r₀ : α) (hr : r₀ ∈ s) (s₀ : β) (hs : s₀ ∈ t) :
    -- Delta corner: Shannon additivity (= 0)
    shannonH s (deltaDist r₀) + shannonH t (deltaDist s₀) = 0
    -- Uniform corner: Shannon additivity (= log (m·n))
    ∧ shannonH s (fun _ => (1 : ℝ) / m) + shannonH t (fun _ => (1 : ℝ) / n)
        = Real.log (m * n)
    -- Delta corner: GFT total additivity (= log (m·n))
    ∧ (klToUniform s m (deltaDist r₀) + shannonH s (deltaDist r₀))
        + (klToUniform t n (deltaDist s₀) + shannonH t (deltaDist s₀))
        = Real.log (m * n)
    -- Uniform corner: GFT total additivity (= log (m·n))
    ∧ (klToUniform s m (fun _ => (1 : ℝ) / m)
          + shannonH s (fun _ => (1 : ℝ) / m))
        + (klToUniform t n (fun _ => (1 : ℝ) / n)
          + shannonH t (fun _ => (1 : ℝ) / n))
        = Real.log (m * n)
    -- Conservation: delta total = uniform total
    ∧ (klToUniform s m (deltaDist r₀) + shannonH s (deltaDist r₀))
        + (klToUniform t n (deltaDist s₀) + shannonH t (deltaDist s₀))
      = (klToUniform s m (fun _ => (1 : ℝ) / m)
          + shannonH s (fun _ => (1 : ℝ) / m))
        + (klToUniform t n (fun _ => (1 : ℝ) / n)
          + shannonH t (fun _ => (1 : ℝ) / n)) :=
  ⟨shannonH_delta_add s t r₀ hr s₀ hs,
   shannonH_uniform_add s t m n hm hn hcard_s hcard_t,
   GFT_delta_add s t m n hm hn r₀ hr s₀ hs,
   GFT_uniform_add s t m n hm hn hcard_s hcard_t,
   GFT_corner_total_conservation s t m n hm hn hcard_s hcard_t r₀ hr s₀ hs⟩

end PT.Information
