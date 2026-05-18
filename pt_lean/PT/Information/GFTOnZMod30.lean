/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Information.GFTIdentity
import PT.Information.BekensteinBound
import PT.Information.GFTSpecialisations
import PT.Information.BekensteinExtensions
import PT.Sieve.Bimodality
import Mathlib.Tactic

/-!
# GFT identity on the 8 admissible residues of `(ℤ/30ℤ)*` (Ch04 #18 application)

This file applies the **Generalised Fundamental Theorem of
Persistence-information** (`PT.Information.GFTIdentity.GFT_identity`) to a
concrete arithmetic setup: the **uniform distribution on the 8 admissible
residues `R = {1, 7, 11, 13, 17, 19, 23, 29}` modulo `30`** (cf.
`PT.Sieve.Bimodality.admissibleResidues`).

We work with the index type `Fin 8` so that `decide` / `native_decide`
discharges the cardinality goals.

## Contents

* `U8` — the uniform distribution `Fin 8 → ℝ`, `U8 i = 1/8`.
* `delta1` — the Dirac distribution on the first residue, `delta1 i = 1`
  iff `i = 0`, else `0`.
* `U8_sum`, `U8_nonneg` — `U8` is a valid probability vector.
* `delta1_sum`, `delta1_nonneg` — `delta1` is a valid probability vector.
* `klToUniform_U8` — `D_KL(U8 ‖ U_8) = 0`.
* `shannonH_U8` — `H(U8) = log 8` (via the GFT identity).
* `klToUniform_delta1` — `D_KL(δ₁ ‖ U_8) = log 8`.
* `shannonH_delta1` — `H(δ₁) = 0`.
* `GFT_uniform_vs_delta` — headline comparison: the uniform corner attains
  entropy maximum and KL minimum, while the delta corner saturates the
  Bekenstein bound and has zero entropy.

## Reference

Monograph Ch04 §"Saturation Bekenstein" + Appendix N (bimodality theorem
identifies these 8 residues as the active classes of the persistence sieve).
-/

namespace PT.Information

open Real Finset

/-! ### The index type and the cardinality `m = 8` -/

/-- The cardinality of the admissible residue set `R` as a real number. -/
noncomputable def m8 : ℝ := 8

lemma m8_pos : 0 < m8 := by unfold m8; norm_num

lemma m8_ge_one : (1 : ℝ) ≤ m8 := by unfold m8; norm_num

/-! ### The uniform distribution `U8` on `Fin 8` -/

/-- The uniform distribution on the 8 admissible residues. -/
noncomputable def U8 : Fin 8 → ℝ := fun _ => (1 : ℝ) / 8

/-- `U8` is non-negative. -/
theorem U8_nonneg : ∀ i ∈ (Finset.univ : Finset (Fin 8)), 0 ≤ U8 i := by
  intro i _
  unfold U8
  norm_num

/-- `U8` is bounded by `1` per coordinate. -/
theorem U8_le_one : ∀ i ∈ (Finset.univ : Finset (Fin 8)), U8 i ≤ 1 := by
  intro i _
  unfold U8
  norm_num

/-- `U8` sums to `1` over `Fin 8`. -/
theorem U8_sum : ∑ i ∈ (Finset.univ : Finset (Fin 8)), U8 i = 1 := by
  unfold U8
  rw [Finset.sum_const]
  simp

/-! ### The Dirac distribution `delta1` on the first residue -/

/-- The Dirac (delta) distribution concentrated on the first residue
    (index `0 : Fin 8`). -/
noncomputable def delta1 : Fin 8 → ℝ := fun i => if i = 0 then 1 else 0

/-- `delta1` is non-negative. -/
theorem delta1_nonneg : ∀ i ∈ (Finset.univ : Finset (Fin 8)), 0 ≤ delta1 i := by
  intro i _
  unfold delta1
  split_ifs <;> [exact zero_le_one; exact le_refl 0]

/-- `delta1` is bounded by `1` per coordinate. -/
theorem delta1_le_one : ∀ i ∈ (Finset.univ : Finset (Fin 8)), delta1 i ≤ 1 := by
  intro i _
  unfold delta1
  split_ifs <;> [exact le_refl 1; exact zero_le_one]

/-- `delta1` sums to `1` over `Fin 8`. -/
theorem delta1_sum : ∑ i ∈ (Finset.univ : Finset (Fin 8)), delta1 i = 1 := by
  unfold delta1
  rw [Finset.sum_ite_eq' (Finset.univ : Finset (Fin 8)) (0 : Fin 8) (fun _ => (1 : ℝ))]
  simp

/-! ### KL divergence and Shannon entropy at the uniform corner -/

/-- **KL divergence at the uniform corner.** `D_KL(U8 ‖ U_8) = 0`. -/
theorem klToUniform_U8 :
    klToUniform (Finset.univ : Finset (Fin 8)) m8 U8 = 0 := by
  unfold U8 m8
  exact klToUniform_uniform (Finset.univ : Finset (Fin 8)) (8 : ℝ) (by norm_num)

/-- **Shannon entropy at the uniform corner.** `H(U8) = log 8`. -/
theorem shannonH_U8 :
    shannonH (Finset.univ : Finset (Fin 8)) U8 = Real.log m8 := by
  have hid := GFT_identity (Finset.univ : Finset (Fin 8)) m8 m8_pos U8 U8_nonneg U8_sum
  have hkl := klToUniform_U8
  linarith

/-! ### KL divergence and Shannon entropy at the delta corner -/

/-- The first index `0 : Fin 8` is in the universe. -/
private lemma zero_mem_univ_fin8 : (0 : Fin 8) ∈ (Finset.univ : Finset (Fin 8)) :=
  Finset.mem_univ _

/-- The local `delta1` coincides with the generic `deltaDist 0` on `Fin 8`. -/
private lemma delta1_eq_deltaDist : delta1 = deltaDist (0 : Fin 8) := by
  funext i
  unfold delta1 deltaDist
  rfl

/-- **KL divergence at the delta corner.** `D_KL(δ₁ ‖ U_8) = log 8`. -/
theorem klToUniform_delta1 :
    klToUniform (Finset.univ : Finset (Fin 8)) m8 delta1 = Real.log m8 := by
  rw [delta1_eq_deltaDist]
  exact klToUniform_delta (Finset.univ : Finset (Fin 8)) m8 m8_pos
    (0 : Fin 8) zero_mem_univ_fin8

/-- **Shannon entropy at the delta corner.** `H(δ₁) = 0`. -/
theorem shannonH_delta1 :
    shannonH (Finset.univ : Finset (Fin 8)) delta1 = 0 := by
  rw [delta1_eq_deltaDist]
  exact shannonH_delta (Finset.univ : Finset (Fin 8))
    (0 : Fin 8) zero_mem_univ_fin8

/-! ### GFT identity for both corners on `Fin 8` -/

/-- **GFT identity at the uniform corner.** `0 + log 8 = log 8`. -/
theorem GFT_at_U8 :
    klToUniform (Finset.univ : Finset (Fin 8)) m8 U8
      + shannonH (Finset.univ : Finset (Fin 8)) U8 = Real.log m8 :=
  GFT_identity (Finset.univ : Finset (Fin 8)) m8 m8_pos U8 U8_nonneg U8_sum

/-- **GFT identity at the delta corner.** `log 8 + 0 = log 8`. -/
theorem GFT_at_delta1 :
    klToUniform (Finset.univ : Finset (Fin 8)) m8 delta1
      + shannonH (Finset.univ : Finset (Fin 8)) delta1 = Real.log m8 :=
  GFT_identity (Finset.univ : Finset (Fin 8)) m8 m8_pos delta1 delta1_nonneg delta1_sum

/-! ### Headline: uniform vs delta on the 8 admissible residues -/

/-- **Headline.** On the 8 admissible residues `R ⊂ (ℤ/30ℤ)*`, the GFT
    identity `log 8 = D_KL + H` admits two interpretable corners:

    * **Uniform corner** (`P = U8`): `D_KL = 0`, `H = log 8` — entropy
      is maximal, KL to uniform is minimal.
    * **Delta corner** (`P = δ₁`): `D_KL = log 8`, `H = 0` — KL
      saturates the Bekenstein bound, entropy collapses to zero.

    Both corners satisfy the GFT partition identity exactly. -/
theorem GFT_uniform_vs_delta :
    -- Uniform: D_KL = 0
    klToUniform (Finset.univ : Finset (Fin 8)) m8 U8 = 0
    -- Uniform: H = log 8
    ∧ shannonH (Finset.univ : Finset (Fin 8)) U8 = Real.log m8
    -- Delta: D_KL = log 8 (Bekenstein saturated)
    ∧ klToUniform (Finset.univ : Finset (Fin 8)) m8 delta1 = Real.log m8
    -- Delta: H = 0
    ∧ shannonH (Finset.univ : Finset (Fin 8)) delta1 = 0
    -- Both partitions equal log 8
    ∧ klToUniform (Finset.univ : Finset (Fin 8)) m8 U8
        + shannonH (Finset.univ : Finset (Fin 8)) U8 = Real.log m8
    ∧ klToUniform (Finset.univ : Finset (Fin 8)) m8 delta1
        + shannonH (Finset.univ : Finset (Fin 8)) delta1 = Real.log m8 :=
  ⟨klToUniform_U8, shannonH_U8, klToUniform_delta1, shannonH_delta1,
   GFT_at_U8, GFT_at_delta1⟩

/-! ### Bekenstein bound on the 8 residues -/

/-- **Bekenstein bound on the 8 admissible residues.** Any probability
    vector `p : Fin 8 → ℝ` satisfies `D_KL(p ‖ U_8) ≤ log 8`. -/
theorem bekenstein_bound_Fin8 (p : Fin 8 → ℝ)
    (hp_nonneg : ∀ i ∈ (Finset.univ : Finset (Fin 8)), 0 ≤ p i)
    (hp_le_one : ∀ i ∈ (Finset.univ : Finset (Fin 8)), p i ≤ 1)
    (hp_sum : ∑ i ∈ (Finset.univ : Finset (Fin 8)), p i = 1) :
    klToUniform (Finset.univ : Finset (Fin 8)) m8 p ≤ Real.log m8 :=
  bekenstein_bound (Finset.univ : Finset (Fin 8)) m8 m8_pos p
    hp_nonneg hp_le_one hp_sum

/-- **Bekenstein bound is saturated at `delta1`.** -/
theorem bekenstein_saturated_delta1 :
    klToUniform (Finset.univ : Finset (Fin 8)) m8 delta1 = Real.log m8 :=
  klToUniform_delta1

/-- **Bekenstein bound is attained at `0` at the uniform `U8`.** -/
theorem bekenstein_attained_U8 :
    klToUniform (Finset.univ : Finset (Fin 8)) m8 U8 = 0 :=
  klToUniform_U8

end PT.Information
