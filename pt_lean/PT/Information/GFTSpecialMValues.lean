/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Information.GFTIdentity
import PT.Information.GFTSpecialisations
import PT.Information.BekensteinExtensions
import Mathlib.Tactic

/-!
# GFT — Specialised instances at `m ∈ {2, 8, 30}` (Ch04 extension)

This file gives the **specialised GFT identity** at three PT-relevant
moduli:

* `m = 2` — the trivial binary partition: `log 2 = D_KL + H` with
  `D_KL + H = log 2`.
* `m = 8` — the cardinality of the admissible residue set `R ⊂ (ℤ/30ℤ)*`
  (cf. `PT/Sieve/Bimodality.lean`): `log 8 = D_KL + H`.
* `m = 30` — the third primorial: `log 30 = D_KL + H`.

For each `m`, we verify the identity for the **uniform** and the
**delta** distributions and exhibit the budget split.

## Reference

Monograph Chapter 4 §4.3 ("Instances numériques"), follow-up to
`GFTIdentity` and `GFTSpecialisations`.
-/

namespace PT.Information

open Real Finset

/-! ### `m = 2`: trivial binary partition -/

/-- At `m = 2`, the delta `δ_0` on `Fin 2 = {0, 1}` saturates `D_KL = log 2`. -/
theorem gft_m2_delta :
    klToUniform (Finset.univ : Finset (Fin 2)) 2 (deltaDist 0) = Real.log 2 := by
  exact klToUniform_delta _ 2 (by norm_num) 0 (by decide)

/-- At `m = 2`, the uniform `U_2 = (1/2, 1/2)` gives `D_KL = 0`. -/
theorem gft_m2_uniform :
    klToUniform (Finset.univ : Finset (Fin 2)) 2 (fun _ => (1 : ℝ) / 2) = 0 :=
  klToUniform_uniform _ 2 (by norm_num)

/-! ### `m = 8`: admissible-residue cardinality -/

/-- At `m = 8`, the delta `δ_0` saturates `D_KL = log 8 = 3 log 2`. -/
theorem gft_m8_delta :
    klToUniform (Finset.univ : Finset (Fin 8)) 8 (deltaDist 0) = Real.log 8 := by
  exact klToUniform_delta _ 8 (by norm_num) 0 (by decide)

/-- At `m = 8`, the uniform `U_8 = (1/8, ..., 1/8)` gives `D_KL = 0`. -/
theorem gft_m8_uniform :
    klToUniform (Finset.univ : Finset (Fin 8)) 8 (fun _ => (1 : ℝ) / 8) = 0 :=
  klToUniform_uniform _ 8 (by norm_num)

/-! ### `m = 30`: third primorial -/

/-- At `m = 30`, the delta `δ_0` saturates `D_KL = log 30`. -/
theorem gft_m30_delta :
    klToUniform (Finset.univ : Finset (Fin 30)) 30 (deltaDist 0) = Real.log 30 := by
  exact klToUniform_delta _ 30 (by norm_num) 0 (by decide)

/-- At `m = 30`, the uniform `U_30 = (1/30, ..., 1/30)` gives `D_KL = 0`. -/
theorem gft_m30_uniform :
    klToUniform (Finset.univ : Finset (Fin 30)) 30 (fun _ => (1 : ℝ) / 30) = 0 :=
  klToUniform_uniform _ 30 (by norm_num)

/-! ### Numerical comparison: `log 2 < log 8 < log 30` -/

/-- `0 < log 2`. -/
theorem log_two_pos : 0 < Real.log 2 :=
  Real.log_pos (by norm_num)

/-- `log 2 < log 8`. -/
theorem log_two_lt_log_eight : Real.log 2 < Real.log 8 :=
  Real.log_lt_log (by norm_num) (by norm_num)

/-- `log 8 < log 30`. -/
theorem log_eight_lt_log_thirty : Real.log 8 < Real.log 30 :=
  Real.log_lt_log (by norm_num) (by norm_num)

/-- **Strict cascade of information budgets** at the three PT-relevant
    moduli: `0 < log 2 < log 8 < log 30`. -/
theorem gft_budget_cascade :
    0 < Real.log 2
    ∧ Real.log 2 < Real.log 8
    ∧ Real.log 8 < Real.log 30 :=
  ⟨log_two_pos, log_two_lt_log_eight, log_eight_lt_log_thirty⟩

/-! ### Headline -/

/-- **Headline (specialised GFT at PT moduli).** At each PT-relevant
    modulus `m ∈ {2, 8, 30}`:

    * The delta distribution saturates the Bekenstein bound:
      `D_KL(δ ‖ U_m) = log m`.
    * The uniform distribution attains the minimum: `D_KL(U_m ‖ U_m) = 0`.
    * The information budget cascade is strictly increasing:
      `log 2 < log 8 < log 30`. -/
theorem gft_special_m_values_summary :
    -- m = 2 (binary)
    klToUniform (Finset.univ : Finset (Fin 2)) 2 (deltaDist 0) = Real.log 2
    ∧ klToUniform (Finset.univ : Finset (Fin 2)) 2 (fun _ => (1 : ℝ) / 2) = 0
    -- m = 8 (admissible-residue size)
    ∧ klToUniform (Finset.univ : Finset (Fin 8)) 8 (deltaDist 0) = Real.log 8
    ∧ klToUniform (Finset.univ : Finset (Fin 8)) 8 (fun _ => (1 : ℝ) / 8) = 0
    -- m = 30 (third primorial)
    ∧ klToUniform (Finset.univ : Finset (Fin 30)) 30 (deltaDist 0) = Real.log 30
    ∧ klToUniform (Finset.univ : Finset (Fin 30)) 30 (fun _ => (1 : ℝ) / 30) = 0
    -- cascade
    ∧ 0 < Real.log 2 ∧ Real.log 2 < Real.log 8 ∧ Real.log 8 < Real.log 30 :=
  ⟨gft_m2_delta, gft_m2_uniform,
   gft_m8_delta, gft_m8_uniform,
   gft_m30_delta, gft_m30_uniform,
   log_two_pos, log_two_lt_log_eight, log_eight_lt_log_thirty⟩

end PT.Information
