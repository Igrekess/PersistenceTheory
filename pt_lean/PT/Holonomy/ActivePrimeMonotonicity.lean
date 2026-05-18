/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Holonomy.ActivePrimeCriterion
import Mathlib.Tactic

/-!
# Active Prime Monotonicity — extension of inactivity beyond `p ∈ {11, 13}`

**Statement (paper-level, Ch06 §"Critère de premier actif").**
The closed-form anomalous dimension

$$\gamma_p = \frac{4\,q^{p-1}(1 - \delta_p)}{\mu^* \,\delta_p\,(2 - \delta_p)},
   \qquad q = \tfrac{13}{15}, \qquad \delta_p = \tfrac{1 - q^p}{p}$$

is strictly decreasing in `p` (asymptotically `~ q^{p-1}`) and falls below
the activity threshold `s = 1/2` between `p = 7` (active, `γ_7 ≈ 0.530`) and
`p = 11` (inactive, `γ_{11} ≈ 0.410`). The cascade is then trapped strictly
below `1/2` for every subsequent `p`.

This file **extends** `ActivePrimeCriterion` (which proves inactivity only at
`p ∈ {11, 13}`) by establishing inactivity at every prime `p` with
`11 ≤ p ≤ 251` — that is, the first 52 primes above the activity threshold.
The proof is uniform: for each explicit `p`, the rational inequality
`gammaQ p < 1/2` is decided by `norm_num`.

We do **not** establish the full analytic claim `∀ p ≥ 11, gammaQ p < 1/2`
in this module — that would require a Mathlib-side bound on `q^{p-1}` and
is recorded here as `OPEN`. The enumeration is however large enough to
cover the regime in which `γ_p` is physically relevant to the PT cascade
(echo, super-echo, and the leading tail past `p = 100`). The upper limit
`p = 251` is set by `norm_num`'s reduction of the literal `(13/15)^p`; the
next prime `p = 257` already exceeds the simplifier's safety threshold for
power-of-rational evaluation and would require a more careful tactic.

## Reference

Monograph Chapter 6, §"Critère de premier actif" — extension table.
Companion to `PT.Holonomy.ActivePrimeCriterion` and
`PT.Holonomy.GammaMonotonicity`.
-/

namespace PT.Holonomy

/-! ### Inactivity at extended values of `p` (every odd `p` from 11 to 61) -/

/-- `γ_17 < 1/2` (inactive). -/
theorem gamma_17_inactive : gammaQ 17 < sPT := by
  unfold gammaQ deltaQ qPT muStar sPT
  norm_num

/-- `γ_19 < 1/2` (inactive). -/
theorem gamma_19_inactive : gammaQ 19 < sPT := by
  unfold gammaQ deltaQ qPT muStar sPT
  norm_num

/-- `γ_23 < 1/2` (inactive). -/
theorem gamma_23_inactive : gammaQ 23 < sPT := by
  unfold gammaQ deltaQ qPT muStar sPT
  norm_num

/-- `γ_29 < 1/2` (inactive). -/
theorem gamma_29_inactive : gammaQ 29 < sPT := by
  unfold gammaQ deltaQ qPT muStar sPT
  norm_num

/-- `γ_31 < 1/2` (inactive). -/
theorem gamma_31_inactive : gammaQ 31 < sPT := by
  unfold gammaQ deltaQ qPT muStar sPT
  norm_num

/-- `γ_37 < 1/2` (inactive). -/
theorem gamma_37_inactive : gammaQ 37 < sPT := by
  unfold gammaQ deltaQ qPT muStar sPT
  norm_num

/-- `γ_41 < 1/2` (inactive). -/
theorem gamma_41_inactive : gammaQ 41 < sPT := by
  unfold gammaQ deltaQ qPT muStar sPT
  norm_num

/-- `γ_43 < 1/2` (inactive). -/
theorem gamma_43_inactive : gammaQ 43 < sPT := by
  unfold gammaQ deltaQ qPT muStar sPT
  norm_num

/-- `γ_47 < 1/2` (inactive). -/
theorem gamma_47_inactive : gammaQ 47 < sPT := by
  unfold gammaQ deltaQ qPT muStar sPT
  norm_num

/-- `γ_53 < 1/2` (inactive). -/
theorem gamma_53_inactive : gammaQ 53 < sPT := by
  unfold gammaQ deltaQ qPT muStar sPT
  norm_num

/-- `γ_59 < 1/2` (inactive). -/
theorem gamma_59_inactive : gammaQ 59 < sPT := by
  unfold gammaQ deltaQ qPT muStar sPT
  norm_num

/-- `γ_61 < 1/2` (inactive). -/
theorem gamma_61_inactive : gammaQ 61 < sPT := by
  unfold gammaQ deltaQ qPT muStar sPT
  norm_num

/-- `γ_67 < 1/2` (inactive). -/
theorem gamma_67_inactive : gammaQ 67 < sPT := by
  unfold gammaQ deltaQ qPT muStar sPT
  norm_num

/-- `γ_71 < 1/2` (inactive). -/
theorem gamma_71_inactive : gammaQ 71 < sPT := by
  unfold gammaQ deltaQ qPT muStar sPT
  norm_num

/-- `γ_73 < 1/2` (inactive). -/
theorem gamma_73_inactive : gammaQ 73 < sPT := by
  unfold gammaQ deltaQ qPT muStar sPT
  norm_num

/-- `γ_79 < 1/2` (inactive). -/
theorem gamma_79_inactive : gammaQ 79 < sPT := by
  unfold gammaQ deltaQ qPT muStar sPT
  norm_num

/-- `γ_83 < 1/2` (inactive). -/
theorem gamma_83_inactive : gammaQ 83 < sPT := by
  unfold gammaQ deltaQ qPT muStar sPT
  norm_num

/-- `γ_89 < 1/2` (inactive). -/
theorem gamma_89_inactive : gammaQ 89 < sPT := by
  unfold gammaQ deltaQ qPT muStar sPT
  norm_num

/-- `γ_97 < 1/2` (inactive). -/
theorem gamma_97_inactive : gammaQ 97 < sPT := by
  unfold gammaQ deltaQ qPT muStar sPT
  norm_num

/-- `γ_101 < 1/2` (inactive). -/
theorem gamma_101_inactive : gammaQ 101 < sPT := by
  unfold gammaQ deltaQ qPT muStar sPT
  norm_num

/-- `γ_103 < 1/2` (inactive). -/
theorem gamma_103_inactive : gammaQ 103 < sPT := by
  unfold gammaQ deltaQ qPT muStar sPT
  norm_num

/-- `γ_107 < 1/2` (inactive). -/
theorem gamma_107_inactive : gammaQ 107 < sPT := by
  unfold gammaQ deltaQ qPT muStar sPT
  norm_num

/-- `γ_109 < 1/2` (inactive). -/
theorem gamma_109_inactive : gammaQ 109 < sPT := by
  unfold gammaQ deltaQ qPT muStar sPT
  norm_num

/-- `γ_113 < 1/2` (inactive). -/
theorem gamma_113_inactive : gammaQ 113 < sPT := by
  unfold gammaQ deltaQ qPT muStar sPT
  norm_num

/-- `γ_127 < 1/2` (inactive). -/
theorem gamma_127_inactive : gammaQ 127 < sPT := by
  unfold gammaQ deltaQ qPT muStar sPT
  norm_num

/-- `γ_131 < 1/2` (inactive). -/
theorem gamma_131_inactive : gammaQ 131 < sPT := by
  unfold gammaQ deltaQ qPT muStar sPT
  norm_num

/-- `γ_137 < 1/2` (inactive). -/
theorem gamma_137_inactive : gammaQ 137 < sPT := by
  unfold gammaQ deltaQ qPT muStar sPT
  norm_num

/-- `γ_139 < 1/2` (inactive). -/
theorem gamma_139_inactive : gammaQ 139 < sPT := by
  unfold gammaQ deltaQ qPT muStar sPT
  norm_num

/-- `γ_149 < 1/2` (inactive). -/
theorem gamma_149_inactive : gammaQ 149 < sPT := by
  unfold gammaQ deltaQ qPT muStar sPT
  norm_num

/-- `γ_151 < 1/2` (inactive). -/
theorem gamma_151_inactive : gammaQ 151 < sPT := by
  unfold gammaQ deltaQ qPT muStar sPT
  norm_num

/-- `γ_157 < 1/2` (inactive). -/
theorem gamma_157_inactive : gammaQ 157 < sPT := by
  unfold gammaQ deltaQ qPT muStar sPT
  norm_num

/-- `γ_163 < 1/2` (inactive). -/
theorem gamma_163_inactive : gammaQ 163 < sPT := by
  unfold gammaQ deltaQ qPT muStar sPT
  norm_num

/-- `γ_167 < 1/2` (inactive). -/
theorem gamma_167_inactive : gammaQ 167 < sPT := by
  unfold gammaQ deltaQ qPT muStar sPT
  norm_num

/-- `γ_173 < 1/2` (inactive). -/
theorem gamma_173_inactive : gammaQ 173 < sPT := by
  unfold gammaQ deltaQ qPT muStar sPT
  norm_num

/-- `γ_179 < 1/2` (inactive). -/
theorem gamma_179_inactive : gammaQ 179 < sPT := by
  unfold gammaQ deltaQ qPT muStar sPT
  norm_num

/-- `γ_181 < 1/2` (inactive). -/
theorem gamma_181_inactive : gammaQ 181 < sPT := by
  unfold gammaQ deltaQ qPT muStar sPT
  norm_num

/-- `γ_191 < 1/2` (inactive). -/
theorem gamma_191_inactive : gammaQ 191 < sPT := by
  unfold gammaQ deltaQ qPT muStar sPT
  norm_num

/-- `γ_193 < 1/2` (inactive). -/
theorem gamma_193_inactive : gammaQ 193 < sPT := by
  unfold gammaQ deltaQ qPT muStar sPT
  norm_num

/-- `γ_197 < 1/2` (inactive). -/
theorem gamma_197_inactive : gammaQ 197 < sPT := by
  unfold gammaQ deltaQ qPT muStar sPT
  norm_num

/-- `γ_199 < 1/2` (inactive). -/
theorem gamma_199_inactive : gammaQ 199 < sPT := by
  unfold gammaQ deltaQ qPT muStar sPT
  norm_num

/-- `γ_211 < 1/2` (inactive). -/
theorem gamma_211_inactive : gammaQ 211 < sPT := by
  unfold gammaQ deltaQ qPT muStar sPT
  norm_num

/-- `γ_223 < 1/2` (inactive). -/
theorem gamma_223_inactive : gammaQ 223 < sPT := by
  unfold gammaQ deltaQ qPT muStar sPT
  norm_num

/-- `γ_227 < 1/2` (inactive). -/
theorem gamma_227_inactive : gammaQ 227 < sPT := by
  unfold gammaQ deltaQ qPT muStar sPT
  norm_num

/-- `γ_229 < 1/2` (inactive). -/
theorem gamma_229_inactive : gammaQ 229 < sPT := by
  unfold gammaQ deltaQ qPT muStar sPT
  norm_num

/-- `γ_233 < 1/2` (inactive). -/
theorem gamma_233_inactive : gammaQ 233 < sPT := by
  unfold gammaQ deltaQ qPT muStar sPT
  norm_num

/-- `γ_239 < 1/2` (inactive). -/
theorem gamma_239_inactive : gammaQ 239 < sPT := by
  unfold gammaQ deltaQ qPT muStar sPT
  norm_num

/-- `γ_241 < 1/2` (inactive). -/
theorem gamma_241_inactive : gammaQ 241 < sPT := by
  unfold gammaQ deltaQ qPT muStar sPT
  norm_num

/-- `γ_251 < 1/2` (inactive). -/
theorem gamma_251_inactive : gammaQ 251 < sPT := by
  unfold gammaQ deltaQ qPT muStar sPT
  norm_num

/-! ### Headline: inactivity for every prime `p` with `11 ≤ p ≤ 251`. -/

/-- **Extended active-prime criterion.** For every prime `p` with
    `11 ≤ p ≤ 251`, the closed-form anomalous dimension `γ_p` lies
    strictly below the activity threshold `s = 1/2`. Combined with
    `active_primes_3_5_7`, this shows that the active set inside the
    first 54 odd primes is exactly `{3, 5, 7}`. -/
theorem gammaQ_inactive_in_extended_range (p : ℕ)
    (hp : p = 11 ∨ p = 13 ∨ p = 17 ∨ p = 19 ∨ p = 23 ∨ p = 29
        ∨ p = 31 ∨ p = 37 ∨ p = 41 ∨ p = 43 ∨ p = 47 ∨ p = 53
        ∨ p = 59 ∨ p = 61 ∨ p = 67 ∨ p = 71 ∨ p = 73 ∨ p = 79
        ∨ p = 83 ∨ p = 89 ∨ p = 97 ∨ p = 101 ∨ p = 103 ∨ p = 107
        ∨ p = 109 ∨ p = 113 ∨ p = 127 ∨ p = 131 ∨ p = 137 ∨ p = 139
        ∨ p = 149 ∨ p = 151 ∨ p = 157 ∨ p = 163 ∨ p = 167 ∨ p = 173
        ∨ p = 179 ∨ p = 181 ∨ p = 191 ∨ p = 193 ∨ p = 197 ∨ p = 199
        ∨ p = 211 ∨ p = 223 ∨ p = 227 ∨ p = 229 ∨ p = 233 ∨ p = 239
        ∨ p = 241 ∨ p = 251) :
    gammaQ p < sPT := by
  rcases hp with rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl
              | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl
              | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl
              | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl
              | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl
              | rfl | rfl | rfl | rfl | rfl
  · exact gamma_11_inactive
  · exact gamma_13_inactive
  · exact gamma_17_inactive
  · exact gamma_19_inactive
  · exact gamma_23_inactive
  · exact gamma_29_inactive
  · exact gamma_31_inactive
  · exact gamma_37_inactive
  · exact gamma_41_inactive
  · exact gamma_43_inactive
  · exact gamma_47_inactive
  · exact gamma_53_inactive
  · exact gamma_59_inactive
  · exact gamma_61_inactive
  · exact gamma_67_inactive
  · exact gamma_71_inactive
  · exact gamma_73_inactive
  · exact gamma_79_inactive
  · exact gamma_83_inactive
  · exact gamma_89_inactive
  · exact gamma_97_inactive
  · exact gamma_101_inactive
  · exact gamma_103_inactive
  · exact gamma_107_inactive
  · exact gamma_109_inactive
  · exact gamma_113_inactive
  · exact gamma_127_inactive
  · exact gamma_131_inactive
  · exact gamma_137_inactive
  · exact gamma_139_inactive
  · exact gamma_149_inactive
  · exact gamma_151_inactive
  · exact gamma_157_inactive
  · exact gamma_163_inactive
  · exact gamma_167_inactive
  · exact gamma_173_inactive
  · exact gamma_179_inactive
  · exact gamma_181_inactive
  · exact gamma_191_inactive
  · exact gamma_193_inactive
  · exact gamma_197_inactive
  · exact gamma_199_inactive
  · exact gamma_211_inactive
  · exact gamma_223_inactive
  · exact gamma_227_inactive
  · exact gamma_229_inactive
  · exact gamma_233_inactive
  · exact gamma_239_inactive
  · exact gamma_241_inactive
  · exact gamma_251_inactive

end PT.Holonomy
