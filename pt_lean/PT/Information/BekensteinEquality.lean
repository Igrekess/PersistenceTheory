/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Information.BekensteinBound

/-!
# Bekenstein Equality — `D_KL(P ‖ U_m) = log m ↔ P` is a Dirac mass

**Statement (paper-level, Ch04 §"Borne de Bekenstein").**
Equality in the Bekenstein bound `D_KL(P ‖ U_m) ≤ log m` holds if and only if
`P` is a point mass concentrated on a single class `r₀ ∈ s`.

**Proof.** By the GFT identity `log m = D_KL + H`, equality `D_KL = log m`
holds iff `H = 0`. Since `H(P) = ∑_r negMulLog(p_r)` and each summand is
non-negative on `[0,1]`, vanishing of the sum forces every term to vanish.
`negMulLog x = 0 ⟺ x = 0 ∨ x = 1` on `[0,1]`. Combined with `∑_r p_r = 1`,
exactly one `p_{r₀} = 1` and all others are `0`.

## References

Monograph chapter 4, §"Borne de Bekenstein", `\label{thm:bekenstein}`.
-/

namespace PT.Information

open Real Finset

/-- `negMulLog x = 0 ↔ x = 0 ∨ x = 1` for `x ∈ [0, 1]`. -/
lemma negMulLog_eq_zero_iff_of_mem_Icc {x : ℝ} (h0 : 0 ≤ x) (h1 : x ≤ 1) :
    Real.negMulLog x = 0 ↔ x = 0 ∨ x = 1 := by
  constructor
  · intro h
    -- negMulLog x = -x * log x = 0 ⇒ x = 0 ∨ log x = 0
    have hmul : x * Real.log x = 0 := by
      have : -(x * Real.log x) = 0 := by
        simpa [Real.negMulLog_eq_neg] using h
      linarith
    rcases mul_eq_zero.mp hmul with hx0 | hlog0
    · exact Or.inl hx0
    · -- log x = 0 ⇒ x ∈ {0, 1, -1}; x ∈ [0,1] forces x = 0 ∨ x = 1
      rcases (Real.log_eq_zero).mp hlog0 with hx | hx | hx
      · exact Or.inl hx
      · exact Or.inr hx
      · -- x = -1 contradicts 0 ≤ x
        exfalso; linarith
  · rintro (rfl | rfl)
    · simp
    · simp

/-- **Bekenstein equality (real-modulus form).**
    For any probability distribution `p : α → ℝ` on a finite set `s` and
    any `m > 0`, equality `D_KL(P ‖ U_m) = log m` holds iff `p` is a
    Dirac mass at some `r₀ ∈ s`. -/
theorem bekenstein_eq_iff {α : Type*} (s : Finset α) (m : ℝ) (hm : 0 < m)
    (p : α → ℝ)
    (hp_nonneg : ∀ r ∈ s, 0 ≤ p r)
    (hp_le_one : ∀ r ∈ s, p r ≤ 1)
    (hp_sum : ∑ r ∈ s, p r = 1) :
    klToUniform s m p = Real.log m ↔
      ∃ r₀ ∈ s, p r₀ = 1 ∧ ∀ r ∈ s, r ≠ r₀ → p r = 0 := by
  have hid := GFT_identity s m hm p hp_nonneg hp_sum
  -- hid : klToUniform s m p + shannonH s p = log m
  -- Step 1: D_KL = log m ↔ H = 0.
  have h_eq_iff_H : klToUniform s m p = Real.log m ↔ shannonH s p = 0 := by
    constructor
    · intro h; linarith
    · intro h; linarith
  rw [h_eq_iff_H]
  -- Step 2: H = 0 ↔ ∀ r ∈ s, negMulLog (p r) = 0.
  have h_sum_zero_iff :
      shannonH s p = 0 ↔ ∀ r ∈ s, Real.negMulLog (p r) = 0 := by
    unfold shannonH
    refine Finset.sum_eq_zero_iff_of_nonneg ?_
    intro r hr
    exact Real.negMulLog_nonneg (hp_nonneg r hr) (hp_le_one r hr)
  rw [h_sum_zero_iff]
  -- Step 3: each negMulLog (p r) = 0 ↔ p r ∈ {0, 1}.
  have h_each_iff :
      (∀ r ∈ s, Real.negMulLog (p r) = 0) ↔ (∀ r ∈ s, p r = 0 ∨ p r = 1) := by
    refine ⟨?_, ?_⟩
    · intro h r hr
      exact (negMulLog_eq_zero_iff_of_mem_Icc (hp_nonneg r hr) (hp_le_one r hr)).mp (h r hr)
    · intro h r hr
      exact (negMulLog_eq_zero_iff_of_mem_Icc (hp_nonneg r hr) (hp_le_one r hr)).mpr (h r hr)
  rw [h_each_iff]
  -- Step 4: combined with ∑ p_r = 1, exactly one p_{r₀} = 1, rest = 0.
  classical
  constructor
  · intro h01
    -- Find r₀ with p r₀ = 1.
    by_contra hne
    push Not at hne
    -- hne : ∀ r ∈ s, p r = 1 → ∃ r' ∈ s, r' ≠ r ∧ p r' ≠ 0
    -- Show ∀ r ∈ s, p r = 0, contradicting ∑ p_r = 1.
    have hall_zero : ∀ r ∈ s, p r = 0 := by
      intro r hr
      rcases h01 r hr with h0 | h1
      · exact h0
      · exfalso
        obtain ⟨r', hr', hrr', hpr'⟩ := hne r hr h1
        rcases h01 r' hr' with h0' | h1'
        · exact hpr' h0'
        · -- both p r = 1 and p r' = 1, with r ≠ r'; sum ≥ 2.
          have hsum_ge : (2 : ℝ) ≤ ∑ x ∈ s, p x := by
            let pair : Finset α := insert r ({r'} : Finset α)
            have hpair : pair ⊆ s := by
              intro x hx
              rcases Finset.mem_insert.mp hx with rfl | hx'
              · exact hr
              · rw [Finset.mem_singleton] at hx'; subst hx'; exact hr'
            have hdisj : r ∉ ({r'} : Finset α) := by
              rw [Finset.mem_singleton]; exact (Ne.symm hrr')
            have hsum_pair : ∑ x ∈ pair, p x = p r + p r' := by
              simp [pair, Finset.sum_insert hdisj]
            calc (2 : ℝ) = p r + p r' := by rw [h1, h1']; norm_num
              _ = ∑ x ∈ pair, p x := hsum_pair.symm
              _ ≤ ∑ x ∈ s, p x :=
                  Finset.sum_le_sum_of_subset_of_nonneg hpair
                    (fun x hxs _ => hp_nonneg x hxs)
          linarith [hp_sum]
    -- Now ∑ p_r = 0, contradicting ∑ p_r = 1.
    have : ∑ r ∈ s, p r = 0 := by
      apply Finset.sum_eq_zero
      intro r hr
      exact hall_zero r hr
    linarith [hp_sum]
  · rintro ⟨r₀, hr₀, hp₀, hrest⟩
    intro r hr
    by_cases hr_eq : r = r₀
    · subst hr_eq; exact Or.inr hp₀
    · exact Or.inl (hrest r hr hr_eq)

/-- **Bekenstein equality for a natural-number modulus `m ≥ 1`.** -/
theorem bekenstein_eq_iff_nat {α : Type*} (s : Finset α) (m : ℕ) (hm : 1 ≤ m)
    (p : α → ℝ)
    (hp_nonneg : ∀ r ∈ s, 0 ≤ p r)
    (hp_le_one : ∀ r ∈ s, p r ≤ 1)
    (hp_sum : ∑ r ∈ s, p r = 1) :
    klToUniform s (m : ℝ) p = Real.log (m : ℝ) ↔
      ∃ r₀ ∈ s, p r₀ = 1 ∧ ∀ r ∈ s, r ≠ r₀ → p r = 0 :=
  bekenstein_eq_iff s (m : ℝ)
    (by exact_mod_cast Nat.lt_of_lt_of_le Nat.zero_lt_one hm)
    p hp_nonneg hp_le_one hp_sum

end PT.Information
