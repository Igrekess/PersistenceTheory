/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Stochastic.T30AntiSector

/-!
# `T_30` Kronecker spectrum: explicit 4 eigenpairs and absolute-value
  partition

The state space `(Fin 1 × Fin 2) × Fin 4` of `T_30 = T_2 ⊗ T_3 ⊗ T_5` has
dimension `1 · 2 · 4 = 8`. Eigenvalues of `T_30` are products
`λ_{T_2} · λ_{T_3} · λ_{T_5}` of factor eigenvalues. The `T5Like` record
only exposes the Perron eigenvalue `+1` and one subdominant eigenvalue
`λ_2(T_5)`. Combined with the full `T_3` spectrum `{+1, -1}` and the
trivial `T_2` spectrum `{+1}`, this yields **four explicit eigenpairs**
of `T_30`:

| Sector | Eigenvalue | Vector |
|---|---|---|
| All-Perron | `+1` | `v_+(T_2) ⊗ v_+(T_3) ⊗ v_+(T_5)` |
| `T_3`-anti / `T_5`-Perron | `-1` | `v_+(T_2) ⊗ v_-(T_3) ⊗ v_+(T_5)` |
| `T_5`-sub | `λ_2(T_5)` | `v_+(T_2) ⊗ v_+(T_3) ⊗ w_2(T_5)` |
| `T_3`-anti / `T_5`-sub | `-λ_2(T_5)` | `v_+(T_2) ⊗ v_-(T_3) ⊗ w_2(T_5)` |

Each is already established in `T2T30Normalisation.lean` and
`T30AntiSector.lean`. This file adds value by:

1. **Absolute-value partition** (`T30_known_abs_values_partition`):
   among the four known eigenvalues, exactly **two have `|·| = 1`**
   (the Perron and anti-Perron sectors) and the **two remaining are
   bounded by `1/4 = s²`** (the sub-dominant pair).
2. **Spectral gap** (`T30_known_spectral_gap`): the four known
   eigenvalues split cleanly with gap `1 - 1/4 = 3/4` between
   the unit-modulus pair and the sub-dominant pair.
3. **Top-pair comparison** (`T30_top_pair_abs_eq`): the Perron and
   anti-Perron eigenvalues have **equal absolute value `1`** — the
   `T_3`-anti sector does **not** strictly dominate Perron.
4. **Sub-pair comparison** (`T30_sub_pair_abs_eq`): the two sub-dominant
   eigenvalues `±λ_2(T_5)` have **equal absolute value `|λ_2(T_5)|`**.
5. **Aggregator** (`T30_four_eigenpairs_with_bounds`): the four
   eigenpairs collected with their absolute-value certificates in one
   statement.

## Status

`[DER]` — pure reuse of `T30AntiSector`. No numerical input beyond
`T5Like.subdominant_bound` already encoded in the base record.
-/

namespace PT.Stochastic

open PT.Sieve

/-! ### Absolute-value bounds, one per eigenvalue -/

/-- **Perron eigenvalue absolute value.** The all-Perron eigenpair of
    `T_30` has eigenvalue exactly `+1`, hence `|·| = 1`. -/
theorem T30_perron_abs : |(1 : ℝ)| = 1 := abs_one

/-- **`T_3`-anti / `T_5`-Perron absolute value.** The anti-Perron
    eigenvalue is `-1`, hence `|·| = 1`. -/
theorem T30_anti3_perron_abs : |(-1 : ℝ)| = 1 := by
  rw [abs_neg]; exact abs_one

/-- **`T_5`-sub absolute value bound.** Restatement of
    `T30_lambda_eff_abs_bound`: `|λ_2(T_5)| ≤ 1/4`. -/
theorem T30_lambda_eff_abs_bound' (T5 : T5Like) :
    |T5.subEigenvalue| ≤ (1 : ℝ) / 4 :=
  T30_lambda_eff_abs_bound T5

/-! ### Pairwise comparisons -/

/-- **Top pair: equal absolute values.** The Perron and anti-Perron
    eigenvalues of `T_30` have **equal absolute value `1`**. Neither
    strictly dominates the other; together they form the unit-modulus
    spectrum on the explicit four-eigenpair subset.

    Conceptually: `T_3`'s `±1` involution lifts into `T_30` as a sign
    flip that preserves Perron magnitude. -/
theorem T30_top_pair_abs_eq : |(1 : ℝ)| = |(-1 : ℝ)| := by
  rw [abs_one, abs_neg, abs_one]

/-- **Sub pair: equal absolute values.** The two sub-dominant
    eigenvalues `λ_2(T_5)` and `-λ_2(T_5)` (`T_5`-sub alone and
    `T_3`-anti combined with `T_5`-sub) have **equal absolute value
    `|λ_2(T_5)|`**.

    Conceptually: the `T_3`-anti sector only contributes a sign; it
    does not modify the radial spectral magnitude inherited from
    `T_5`. -/
theorem T30_sub_pair_abs_eq (T5 : T5Like) :
    |T5.subEigenvalue| = |(-T5.subEigenvalue)| := by
  rw [abs_neg]

/-! ### Spectral gap on the explicit 4-eigenpair subset -/

/-- **Spectral gap on the known subset.** Among the four explicit
    eigenpairs of `T_30`, the two top eigenvalues `±1` strictly
    dominate the two sub-dominant eigenvalues `±λ_2(T_5)`:

    `|±λ_2(T_5)| ≤ 1/4 < 3/4 = 1 - 1/4 = |±1| - 1/4`,

    so the gap between `|±1| = 1` and the sub-dominant magnitude is
    at least `3/4`. This is the discrete-spectrum analogue of the
    Perron–Frobenius gap restricted to the explicit 4-vector subset
    enumerated here. -/
theorem T30_known_spectral_gap (T5 : T5Like) :
    |T5.subEigenvalue| ≤ (1 : ℝ) - 3 / 4 := by
  have h := T5.subdominant_bound
  linarith

/-- **Same gap, for the anti-sub eigenvalue.** -/
theorem T30_known_spectral_gap_anti (T5 : T5Like) :
    |(-T5.subEigenvalue)| ≤ (1 : ℝ) - 3 / 4 := by
  rw [abs_neg]; exact T30_known_spectral_gap T5

/-! ### Partition of absolute values: two `= 1`, two `≤ 1/4` -/

/-- **Absolute-value partition of the four known eigenvalues.**
    Among the four explicit eigenpairs of `T_30`:

    * Exactly **two have `|·| = 1`**: the Perron `+1` and the
      anti-Perron `-1`.
    * The **other two are bounded by `1/4`**: `±λ_2(T_5)`.

    This is the structural content of "Perron + anti-Perron sector
    versus subdominant sector" packaged as a single absolute-value
    inequality system. -/
theorem T30_known_abs_values_partition (T5 : T5Like) :
    -- Top sector: both `|·| = 1`
    |(1 : ℝ)| = 1
    ∧ |(-1 : ℝ)| = 1
    -- Sub sector: both `|·| ≤ 1/4`
    ∧ |T5.subEigenvalue| ≤ (1 : ℝ) / 4
    ∧ |(-T5.subEigenvalue)| ≤ (1 : ℝ) / 4 :=
  ⟨T30_perron_abs,
   T30_anti3_perron_abs,
   T5.subdominant_bound,
   T30_anti3_sub_abs_bound T5⟩

/-! ### Aggregator: eigenpairs + their absolute-value certificates -/

/-- **Aggregator: four eigenpairs with absolute-value certificates.**
    All four explicit eigenpairs of `T_30` exposed by the `T5Like`
    parameterisation, together with their absolute values (`= 1` or
    `≤ 1/4`), in a single statement.

    Value added over `T30_perron_anti_sector_summary`:
    * the four eigenvalues are paired with **explicit `|·|`-witnesses**
      (rather than only the inequality on `|λ_2|`);
    * the two top eigenvalues are exhibited as **both having `|·| = 1`**;
    * the two sub eigenvalues are exhibited as **both bounded by `1/4`**.

    This is the "four-eigenpair Kronecker spectrum certificate" used
    downstream to argue that the convergence-controlling magnitude of
    `T_30` on the explicit subset is exactly `1/4`. -/
theorem T30_four_eigenpairs_with_bounds (T5 : T5Like) :
    -- (1) Perron eigenpair with `|·| = 1`
    ((T30 T5).mulVec (T30_perronVec T5) = (1 : ℝ) • T30_perronVec T5
     ∧ |(1 : ℝ)| = 1)
    -- (2) Anti-Perron eigenpair with `|·| = 1`
    ∧ ((T30 T5).mulVec (T30_anti3_perronVec T5)
         = (-1 : ℝ) • T30_anti3_perronVec T5
       ∧ |(-1 : ℝ)| = 1)
    -- (3) T_5-sub eigenpair with `|·| ≤ 1/4`
    ∧ ((T30 T5).mulVec (T30_lambda_eff_vec T5)
         = T5.subEigenvalue • T30_lambda_eff_vec T5
       ∧ |T5.subEigenvalue| ≤ (1 : ℝ) / 4)
    -- (4) Anti-sub eigenpair with `|·| ≤ 1/4`
    ∧ ((T30 T5).mulVec (T30_anti3_sub_Vec T5)
         = (-T5.subEigenvalue) • T30_anti3_sub_Vec T5
       ∧ |(-T5.subEigenvalue)| ≤ (1 : ℝ) / 4) :=
  ⟨⟨T30_perron_eigen T5, T30_perron_abs⟩,
   ⟨T30_anti3_perron_eigen' T5, T30_anti3_perron_abs⟩,
   ⟨T30_lambda_eff_eigen' T5, T5.subdominant_bound⟩,
   ⟨T30_anti3_sub_eigen' T5, T30_anti3_sub_abs_bound T5⟩⟩

/-! ### Headline: the spectral radius on the explicit subset is `1` -/

/-- **Headline.** The spectral radius restricted to the four explicit
    eigenpairs of `T_30` is `max(|1|, |-1|, |λ_2|, |-λ_2|) = 1`.

    Quantitatively: every one of the four eigenvalues has absolute
    value bounded by `1`, because `|±1| = 1` and `|±λ_2| ≤ 1/4 ≤ 1`. -/
theorem T30_known_spectral_radius_le_one (T5 : T5Like) :
    |(1 : ℝ)| ≤ 1
    ∧ |(-1 : ℝ)| ≤ 1
    ∧ |T5.subEigenvalue| ≤ 1
    ∧ |(-T5.subEigenvalue)| ≤ 1 := by
  refine ⟨?_, ?_, ?_, ?_⟩
  · rw [abs_one]
  · rw [abs_neg, abs_one]
  · have h := T5.subdominant_bound; linarith
  · rw [abs_neg]; have h := T5.subdominant_bound; linarith

end PT.Stochastic
