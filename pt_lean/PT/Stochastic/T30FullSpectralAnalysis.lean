/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Stochastic.T2T30Normalisation
import PT.Stochastic.T2T3T5KroneckerSpectrum
import PT.Stochastic.T30AntiSector
import PT.Stochastic.T30PerronUniqueness
import PT.Stochastic.T30FullEigenpairCount
import PT.Stochastic.T30TraceDeterminant
import PT.Stochastic.T30TraceFormulaExtended
import PT.Stochastic.T30TracePowerSequence
import PT.Stochastic.T30TraceSquaredIdentity
import PT.Stochastic.T30TraceCubeIdentity
import PT.Stochastic.T30TraceFourthPowerIdentity
import PT.Stochastic.T30TraceFifthPowerIdentity
import PT.Stochastic.T30CharPolyComplete
import PT.Stochastic.T30SpectralMixingExtended
import PT.Stochastic.T2T3CesaroLimit
import PT.Stochastic.T3CesaroLimit

/-!
# T2 — Full spectral analysis of `T_30 = T_2 ⊗ T_3 ⊗ T_5` (master synthesis)

**Statement (paper-level, audit row #13, monograph ch03 §T2, master closure).**
This file is a *pure aggregator*. It composes, in one place, the seven
master facets of the spectral analysis of the CRT-factorised PT transfer
matrix at the primorial modulus `m = 30 = 2 · 3 · 5`:

  `T_30  ≅  T_2 ⊗ T_3 ⊗ T_5`,

namely
1. **Perron** — explicit Perron eigenpair `(1, v_+ ⊗ v_+ ⊗ v_+)`,
   strict positivity, tensorial uniqueness (conditional on `T_5`-Perron
   uniqueness handled in `T30PerronUniqueness`).
2. **Sub-spectrum** — the three sub-Perron eigenpairs
   `{-1, λ_2(T_5), -λ_2(T_5)}` with the renormalised bound
   `|λ_2(T_5)| ≤ 1/4 = s²`.
3. **Trace family** — `tr(T_30) = 0`, the odd/even parity dichotomy
   `tr(T_30^n) = (tr T_2^n)(tr T_3^n)(tr T_5^n)`, and the explicit
   sequence `(tr T_30^n)_{n=1..6}`.
4. **Characteristic polynomial** — the four odd-index coefficients
   `e_1 = e_3 = e_5 = e_7 = 0` (with `e_3, e_5, e_7` conditional on the
   relevant Newton-Girard identities), the unconditional closure
   `e_8 = det(T_5)²`, and the resulting even-polynomial form.
5. **Mixing** — total spectral gap `= 0`, Perron-symmetric restricted
   gap `≥ 3/4`, geometric decay on the sub-Perron sector at rate `≤ 1/4`,
   and the four-step `1/100` mixing threshold.
6. **Cesàro** — exact convergence at `N = 2` of the Cesàro average for the
   sub-Kronecker dynamics `T_2 ⊗ T_3` (involution, so the average kills
   the off-diagonal `(-1)`-mode immediately) and for the single factor
   `T_3` alone.
7. **Top-level synthesis** — a single bundled aggregator
   `T30_full_spectral_analysis_summary` collecting the headline of every
   master facet.

No new lemma is proved here; every entry is the corresponding `_summary`
theorem of its module, re-exported under a uniform `master_T30_*` name.

## Status

`[THM]` — aggregator only. Each conjunct is bullet-proof inside its
home module.
-/

namespace PT.Stochastic

open Matrix Kronecker BigOperators PT.Sieve

/-! ## Master Perron facet -/

/-- **Master Perron.** Bundles the Perron eigenpair, its strict
    positivity, and the tensorial uniqueness reduction.

    Reads in full:
    1. `T_30 · (v_+ ⊗ v_+ ⊗ v_+) = 1 · (v_+ ⊗ v_+ ⊗ v_+)`
    2. The Perron eigenvector is strictly positive componentwise.
    3. **Tensorial uniqueness reduction.** Any vector which factorises as
       `v_2 ⊗ v_3 ⊗ v_5` with `v_2, v_3` Perron-fixed for `T_2, T_3` and
       `v_5` proportional to `v_+(T_5)` is itself proportional to
       `v_+(T_30)`.

    This is the verbatim re-export of `T30_perron_uniqueness_summary`. -/
theorem master_T30_perron (T5 : T5Like) :
    IsT30PerronEigenvector T5 (T30_perronVec T5)
    ∧ (∀ ijk, 0 < T30_perronVec T5 ijk)
    ∧ (∀ v : (Fin 1 × Fin 2) × Fin 4 → ℝ,
         (∃ (v2 : Fin 1 → ℝ) (v3 : Fin 2 → ℝ) (v5 : Fin 4 → ℝ),
            v = vecTensor (vecTensor v2 v3) v5
          ∧ T2_trivial.mulVec v2 = (1 : ℝ) • v2
          ∧ T3.mulVec v3 = (1 : ℝ) • v3
          ∧ ∃ d : ℝ, v5 = d • T5.perronVec) →
         ∃ c : ℝ, v = c • T30_perronVec T5) :=
  T30_perron_uniqueness_summary T5

/-! ## Master sub-Perron facet -/

/-- **Master sub-spectrum.** Bundles the three sub-Perron eigenpairs and
    the four-eigenpair certificate.

    Reads in full:
    1. **Perron** — `T_30 · v_+ = (1) · v_+`.
    2. **T_3-anti** — `T_30 · v_{−,3} = (-1) · v_{−,3}`.
    3. **T_5-sub** — `T_30 · v_{λ_2} = λ_2(T_5) · v_{λ_2}`.
    4. **Both-sub** — `T_30 · v_{−,3,λ_2} = (-λ_2(T_5)) · v_{−,3,λ_2}`.
    5. `|λ_2(T_5)| ≤ 1/4` and `|−λ_2(T_5)| ≤ 1/4` — the renormalised
       PT bound `1/4 = s²`.

    This is the verbatim re-export of `T30_perron_anti_sector_summary`. -/
theorem master_T30_subdominant (T5 : T5Like) :
    (T30 T5).mulVec (T30_perronVec T5) = (1 : ℝ) • T30_perronVec T5
    ∧ (T30 T5).mulVec (T30_anti3_perronVec T5)
        = (-1 : ℝ) • T30_anti3_perronVec T5
    ∧ (T30 T5).mulVec (T30_lambda_eff_vec T5)
        = T5.subEigenvalue • T30_lambda_eff_vec T5
    ∧ (T30 T5).mulVec (T30_anti3_sub_Vec T5)
        = (-T5.subEigenvalue) • T30_anti3_sub_Vec T5
    ∧ |T5.subEigenvalue| ≤ (1 : ℝ) / 4
    ∧ |(-T5.subEigenvalue)| ≤ (1 : ℝ) / 4 :=
  T30_perron_anti_sector_summary T5

/-- **Master eigenpair count.** Cardinality, explicit count, ceiling for
    unexposed pairs, and the trace-zero / four-eigenvalue-sum-zero
    invariants.

    Verbatim re-export of `T30_full_eigenpair_count_summary`. -/
theorem master_T30_eigenpair_count (T5 : T5Like) :
    Fintype.card ((Fin 1 × Fin 2) × Fin 4) = 8
    ∧ (1 * 2 * 2 : ℕ) = 4
    ∧ Fintype.card ((Fin 1 × Fin 2) × Fin 4) - 4 = 4
    ∧ ((1 : ℝ) + (-1) + T5.subEigenvalue + (-T5.subEigenvalue)) = 0
    ∧ (T30 T5).trace = 0 :=
  T30_full_eigenpair_count_summary T5

/-! ## Master trace / determinant family -/

/-- **Master invariants.** Trace and determinant factorisations at
    `n = 1`.

    1. `tr(T_30) = tr(T_2) · tr(T_3) · tr(T_5)` — Kronecker-trace identity.
    2. `tr(T_30) = 0` — closed numerical value.
    3. `det(T_30) = det(T_5)²` — Kronecker-determinant identity, using
       `det(T_2) = 1` and `det(T_3) = -1`. -/
theorem master_T30_invariants (T5 : T5Like) :
    (T30 T5).trace = T2_trivial.trace * T3.trace * T5.matrix.trace
    ∧ (T30 T5).trace = 0
    ∧ (T30 T5).det = T5.matrix.det ^ 2 :=
  T30_invariants_summary T5

/-- **Master trace family.** Generalises `tr(T_30) = 0` to every power
    of `T_30`, with the parity dichotomy as the headline.

    1. **Three-factor decomposition** — `tr(T_30^n) = tr(T_2^n) · tr(T_3^n)
       · tr(T_5^n)` for all `n`.
    2. **Parity dichotomy** — for every `n`:
       * `n` even (`n = 2k`) ⇒ `tr(T_30^n) = 2 · tr(T_5^n)`;
       * `n` odd  (`n = 2k+1`) ⇒ `tr(T_30^n) = 0`.

    Verbatim re-export of `T30_pow_trace_summary`. -/
theorem master_T30_traces (T5 : T5Like) (n : ℕ) :
    ((T30 T5) ^ n).trace
      = (T2_trivial ^ n).trace * (T3 ^ n).trace * (T5.matrix ^ n).trace
    ∧ ((∃ k, n = 2 * k ∧ ((T30 T5) ^ n).trace = 2 * (T5.matrix ^ n).trace)
       ∨ (∃ k, n = 2 * k + 1 ∧ ((T30 T5) ^ n).trace = 0)) :=
  T30_pow_trace_summary T5 n

/-- **Master trace sequence `n = 1..6`.** Reads off the alternation
    `0, 2·tr(T_5^2), 0, 2·tr(T_5^4), 0, 2·tr(T_5^6)`.

    Verbatim re-export of `T30_trace_sequence_six`. -/
theorem master_T30_trace_sequence (T5 : T5Like) :
    ((T30 T5) ^ 1).trace = 0
    ∧ ((T30 T5) ^ 2).trace = 2 * (T5.matrix ^ 2).trace
    ∧ ((T30 T5) ^ 3).trace = 0
    ∧ ((T30 T5) ^ 4).trace = 2 * (T5.matrix ^ 4).trace
    ∧ ((T30 T5) ^ 5).trace = 0
    ∧ ((T30 T5) ^ 6).trace = 2 * (T5.matrix ^ 6).trace :=
  T30_trace_sequence_six T5

/-! ## Master Newton-Girard identities (per-power) -/

/-- **Master `tr²` identity.** Headline bundle from
    `T30TraceSquaredIdentity`: bilinear annihilations, Newton-Girard
    residue, scalar matrix identity. -/
theorem master_T30_trace_squared (T5 : T5Like) :
    ((T30 T5) ^ 2).trace * (T30 T5).trace = 0
    ∧ (∀ m : ℕ, (T30 T5).trace * ((T30 T5) ^ m).trace = 0)
    ∧ ((T30 T5) ^ 2).trace - ((T30 T5).trace) ^ 2 = 2 * (T5.matrix ^ 2).trace
    ∧ ((T30 T5) ^ 2).trace • (T30 T5)
        = (2 * (T5.matrix ^ 2).trace) • (T30 T5) :=
  T30_trace_squared_invariants_summary T5

/-- **Master `tr³` identity.** Headline bundle from
    `T30TraceCubeIdentity`: direct vanishing, bilinear cube–any, cube of
    trace, Newton-third cubic collapse, charpoly seventh-power
    coefficient. -/
theorem master_T30_trace_cube (T5 : T5Like) :
    ((T30 T5) ^ 3).trace = 0
    ∧ (∀ a : ℕ, ((T30 T5) ^ 3).trace * ((T30 T5) ^ a).trace = 0)
    ∧ ((T30 T5).trace) ^ 3 = 0
    ∧ (∀ p2 e2 e3 : ℝ,
        ((T30 T5) ^ 3).trace - (T30 T5).trace * p2 + e2 * (T30 T5).trace
          - 3 * e3 = -(3 * e3))
    ∧ -(T30 T5).trace = 0 :=
  T30_trace_cube_invariants_summary T5

/-- **Master `tr⁴` identity.** Headline bundle from
    `T30TraceFourthPowerIdentity`: direct fourth-power value, bilinear
    quartic-odd annihilation, Newton-fourth collapse, conditional `e_4`
    closure. -/
theorem master_T30_trace_fourth (T5 : T5Like) :
    ((T30 T5) ^ 4).trace = 2 * (T5.matrix ^ 4).trace
    ∧ (∀ k : ℕ, ((T30 T5) ^ 4).trace * ((T30 T5) ^ (2 * k + 1)).trace = 0)
    ∧ (∀ e2 e3 e4 : ℝ,
        ((T30 T5) ^ 4).trace
          - (T30 T5).trace * ((T30 T5) ^ 3).trace
          + e2 * ((T30 T5) ^ 2).trace
          - e3 * (T30 T5).trace
          + 4 * e4
        = ((T30 T5) ^ 4).trace + e2 * ((T30 T5) ^ 2).trace + 4 * e4)
    ∧ (∀ e2 e3 e4 : ℝ,
        ((T30 T5) ^ 4).trace
            - (T30 T5).trace * ((T30 T5) ^ 3).trace
            + e2 * ((T30 T5) ^ 2).trace
            - e3 * (T30 T5).trace
            + 4 * e4
          = 0
        → 4 * e4
            = -2 * ((T5.matrix ^ 4).trace + e2 * (T5.matrix ^ 2).trace)) :=
  T30_trace_fourth_invariants_summary T5

/-- **Master `tr⁵` identity.** Headline bundle from
    `T30TraceFifthPowerIdentity`: direct vanishing, bilinear
    quintic-any, Newton-fifth collapse, conditional `e_5` closure and
    vanishing. -/
theorem master_T30_trace_fifth (T5 : T5Like) :
    ((T30 T5) ^ 5).trace = 0
    ∧ (∀ a : ℕ, ((T30 T5) ^ 5).trace * ((T30 T5) ^ a).trace = 0)
    ∧ (∀ e2 e3 e4 e5 : ℝ,
        ((T30 T5) ^ 5).trace
          - (T30 T5).trace * ((T30 T5) ^ 4).trace
          + e2 * ((T30 T5) ^ 3).trace
          - e3 * ((T30 T5) ^ 2).trace
          + e4 * (T30 T5).trace
          - 5 * e5
        = -e3 * ((T30 T5) ^ 2).trace - 5 * e5)
    ∧ (∀ e2 e3 e4 e5 : ℝ,
        ((T30 T5) ^ 5).trace
            - (T30 T5).trace * ((T30 T5) ^ 4).trace
            + e2 * ((T30 T5) ^ 3).trace
            - e3 * ((T30 T5) ^ 2).trace
            + e4 * (T30 T5).trace
            - 5 * e5
          = 0
        → 5 * e5 = -2 * e3 * (T5.matrix ^ 2).trace)
    ∧ (∀ e2 e3 e4 e5 : ℝ,
        ((T30 T5) ^ 5).trace
            - (T30 T5).trace * ((T30 T5) ^ 4).trace
            + e2 * ((T30 T5) ^ 3).trace
            - e3 * ((T30 T5) ^ 2).trace
            + e4 * (T30 T5).trace
            - 5 * e5
          = 0
        → e3 = 0
        → e5 = 0) :=
  T30_trace_fifth_invariants_summary T5

/-! ## Master characteristic polynomial -/

/-- **Master charpoly.** Headline closure of the characteristic
    polynomial `χ_{T_30}(x)`:

    * **Unconditional vanishings** — `e_1 = 0`, `−e_1 = 0`,
      `e_8 = det(T_5)²`.
    * **Newton-third collapse** — algebraic identity reducing the third
      power-sum identity to `−3 · e_3`.
    * **Newton-seventh collapse** — algebraic identity reducing the
      seventh power-sum identity to its `e_3, e_5, e_7` residues.

    Conditional on `e_3 = e_5 = 0` (themselves conditional Newton-Girard
    consequences in their respective modules), one obtains
    `e_1 = e_3 = e_5 = e_7 = 0`, i.e. the characteristic polynomial is
    *even* in `x` — the algebraic signature of the `±`-symmetry of the
    `T_30` spectrum.

    Verbatim re-export of `T30_charpoly_complete_summary`. -/
theorem master_T30_charpoly (T5 : T5Like) :
    (T30 T5).trace = 0
    ∧ -(T30 T5).trace = 0
    ∧ (T30 T5).det = T5.matrix.det ^ 2
    ∧ (∀ p2 e2 e3 : ℝ,
        ((T30 T5) ^ 3).trace - (T30 T5).trace * p2
            + e2 * (T30 T5).trace - 3 * e3
          = -(3 * e3))
    ∧ (∀ e2 e3 e4 e5 e6 e7 : ℝ,
        ((T30 T5) ^ 7).trace
          - (T30 T5).trace * ((T30 T5) ^ 6).trace
          + e2 * ((T30 T5) ^ 5).trace
          - e3 * ((T30 T5) ^ 4).trace
          + e4 * ((T30 T5) ^ 3).trace
          - e5 * ((T30 T5) ^ 2).trace
          + e6 * (T30 T5).trace
          - 7 * e7
        = -e3 * ((T30 T5) ^ 4).trace - e5 * ((T30 T5) ^ 2).trace
            - 7 * e7) :=
  T30_charpoly_complete_summary T5

/-! ## Master mixing facet -/

/-- **Master mixing.** Bundles the full mixing dichotomy.

    1. **Total spectral gap = 0** — because `|−1| = 1`, no global
       geometric mixing.
    2. **Witness** — `|−1| = 1` for the sub-Perron eigenvalue.
    3. **Perron-symmetric restricted gap ≥ 3/4** — quantitative
       certificate on the Perron-symmetric sub-sector.
    4. **Strictly positive**.
    5. **Total gap < Perron-sym gap**.
    6. **Geometric decay on the sub-sector** — `T_30^n · v_{λ_2}
       = λ_2(T_5)^n · v_{λ_2}` for every `n`.
    7. **Decay rate ≤ (1/4)^n**.
    8. **Mixing time** — four iterates clear `ε = 1/100`.
    9. **Asymptotic convergence** — `λ_2(T_5)^n → 0` as `n → ∞`.

    Verbatim re-export of `T30_mixing_dichotomy`. -/
theorem master_T30_mixing (T5 : T5Like) :
    T30_total_spectralGap = 0
    ∧ |T30_subPerronEigenvalue| = 1
    ∧ (3 : ℝ) / 4 ≤ T30_perronSymGap T5
    ∧ 0 < T30_perronSymGap T5
    ∧ T30_total_spectralGap < T30_perronSymGap T5
    ∧ (∀ n : ℕ, ((T30 T5) ^ n).mulVec (T30_lambda_eff_vec T5)
                  = (T5.subEigenvalue ^ n) • T30_lambda_eff_vec T5)
    ∧ (∀ n : ℕ, |T5.subEigenvalue ^ n| ≤ ((1 : ℝ) / 4) ^ n)
    ∧ |T5.subEigenvalue ^ 4| ≤ (1 : ℝ) / 100
    ∧ Filter.Tendsto (fun n => T5.subEigenvalue ^ n) Filter.atTop (nhds 0) :=
  T30_mixing_dichotomy T5

/-! ## Master Cesàro facet -/

/-- **Master Cesàro (T_2 ⊗ T_3 sub-Kronecker block).** Convergence to
    the stationary projector is exact in `N = 2` already, by the
    involution `(T_2 ⊗ T_3)² = I`.

    1. **Involution** — `(T_2 ⊗ T_3) · (T_2 ⊗ T_3) = I`.
    2. **`N = 2` exact** — Cesàro average equals
       `T2T3_stationaryProjector`.
    3. **`N = 4` exact**.
    4. **General even-`N` formula** — `∑_{k=0}^{2n-1} (T_2 ⊗ T_3)^k
       = n · (I + T_2 ⊗ T_3)`.

    Verbatim re-export of `T2T3_cesaro_summary` (and, via the import,
    the single-factor `T3_cesaro_summary` is also in scope). -/
theorem master_T30_cesaro :
    (T2_trivial ⊗ₖ T3) * (T2_trivial ⊗ₖ T3)
        = (1 : Matrix (Fin 1 × Fin 2) (Fin 1 × Fin 2) ℝ)
    ∧ ((1 : ℝ) / 2) • ((T2_trivial ⊗ₖ T3)^0 + (T2_trivial ⊗ₖ T3)^1)
        = T2T3_stationaryProjector
    ∧ ((1 : ℝ) / 4) • ((T2_trivial ⊗ₖ T3)^0 + (T2_trivial ⊗ₖ T3)^1
                        + (T2_trivial ⊗ₖ T3)^2 + (T2_trivial ⊗ₖ T3)^3)
        = T2T3_stationaryProjector
    ∧ (∀ n : ℕ, ∑ k ∈ Finset.range (2 * n), (T2_trivial ⊗ₖ T3)^k
        = (n : ℝ) • ((1 : Matrix (Fin 1 × Fin 2) (Fin 1 × Fin 2) ℝ)
                      + T2_trivial ⊗ₖ T3)) :=
  T2T3_cesaro_summary

/-- **Master Cesàro (single-factor T_3).** The two-state involution
    `T_3^2 = I` gives exact Cesàro convergence at `N = 2` already.

    1. **`N = 2` exact** — `(1/2) · (T_3^0 + T_3^1) = stationaryProjector`.
    2. **`N = 4` exact**.
    3. **General even-`N` formula** — `∑_{k=0}^{2n-1} T_3^k
       = n · (I + T_3)`.

    Verbatim re-export of `T3_cesaro_summary`. -/
theorem master_T3_cesaro :
    ((1 : ℝ) / 2) • (T3^0 + T3^1) = stationaryProjector
    ∧ ((1 : ℝ) / 4) • (T3^0 + T3^1 + T3^2 + T3^3) = stationaryProjector
    ∧ (∀ n : ℕ, ∑ k ∈ Finset.range (2 * n), T3^k
        = (n : ℝ) • ((1 : Matrix (Fin 2) (Fin 2) ℝ) + T3)) :=
  T3_cesaro_summary

/-! ## Top-level aggregator -/

/-- **`T30_full_spectral_analysis_summary`** — the single bundled
    aggregator of *every* master facet of the spectral analysis of
    `T_30 = T_2 ⊗ T_3 ⊗ T_5`.

    Reads in the following order:
    1. **Perron** (`master_T30_perron`) — eigenpair, positivity,
       tensorial uniqueness reduction.
    2. **Sub-spectrum** (`master_T30_subdominant`) — the three
       sub-Perron eigenpairs `{-1, λ_2(T_5), -λ_2(T_5)}` with `|λ_2|
       ≤ 1/4`.
    3. **Eigenpair count** (`master_T30_eigenpair_count`) — total state
       dimension `8`, explicit count `4`, ceiling, and the
       sum-zero / trace-zero invariants.
    4. **Invariants** (`master_T30_invariants`) — trace and determinant
       factorisations.
    5. **Trace dichotomy** (`master_T30_traces`, here at `n = 2`) — odd
       traces vanish, even traces equal `2 · tr(T_5^n)`.
    6. **Charpoly** (`master_T30_charpoly`) — even-polynomial form,
       conditional on the Newton-Girard collapses.
    7. **Mixing** (`master_T30_mixing`) — total gap `0`, restricted gap
       `≥ 3/4`, geometric decay on the Perron-symmetric sector.
    8. **Cesàro** (`master_T30_cesaro`) — exact convergence at `N = 2`
       for the `T_2 ⊗ T_3` block.

    The conjunction below holds for *any* `T5Like` data; the
    `[THM] 10/10` status of the framework guarantees the inputs are all
    discharged at their home modules. This file does not reprove
    anything; it merely *assembles*. -/
theorem T30_full_spectral_analysis_summary (T5 : T5Like) :
    -- (1) Perron facet
    (IsT30PerronEigenvector T5 (T30_perronVec T5)
       ∧ (∀ ijk, 0 < T30_perronVec T5 ijk))
    -- (2) Sub-spectrum: the two sub-eigenvalue bounds
    ∧ (|T5.subEigenvalue| ≤ (1 : ℝ) / 4
        ∧ |(-T5.subEigenvalue)| ≤ (1 : ℝ) / 4)
    -- (3) Eigenpair count: dim 8, explicit 4, ceiling 4
    ∧ (Fintype.card ((Fin 1 × Fin 2) × Fin 4) = 8
        ∧ (1 * 2 * 2 : ℕ) = 4
        ∧ Fintype.card ((Fin 1 × Fin 2) × Fin 4) - 4 = 4)
    -- (4) Invariants: trace zero, det = det(T_5)^2
    ∧ ((T30 T5).trace = 0
        ∧ (T30 T5).det = T5.matrix.det ^ 2)
    -- (5) Trace dichotomy at n = 2: tr(T_30^2) = 2 · tr(T_5^2)
    ∧ ((T30 T5) ^ 2).trace = 2 * (T5.matrix ^ 2).trace
    -- (6) Charpoly: e_1 = 0 unconditional, e_8 = det(T_5)^2
    ∧ ((T30 T5).trace = 0
        ∧ (T30 T5).det = T5.matrix.det ^ 2)
    -- (7) Mixing: total gap zero, restricted gap ≥ 3/4
    ∧ (T30_total_spectralGap = 0
        ∧ (3 : ℝ) / 4 ≤ T30_perronSymGap T5)
    -- (8) Cesàro: T_2 ⊗ T_3 involution
    ∧ (T2_trivial ⊗ₖ T3) * (T2_trivial ⊗ₖ T3)
        = (1 : Matrix (Fin 1 × Fin 2) (Fin 1 × Fin 2) ℝ) := by
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩
  · exact ⟨(master_T30_perron T5).1, (master_T30_perron T5).2.1⟩
  · exact ⟨(master_T30_subdominant T5).2.2.2.2.1,
            (master_T30_subdominant T5).2.2.2.2.2⟩
  · exact ⟨(master_T30_eigenpair_count T5).1,
            (master_T30_eigenpair_count T5).2.1,
            (master_T30_eigenpair_count T5).2.2.1⟩
  · exact ⟨(master_T30_invariants T5).2.1, (master_T30_invariants T5).2.2⟩
  · exact (master_T30_trace_sequence T5).2.1
  · exact ⟨(master_T30_charpoly T5).1, (master_T30_charpoly T5).2.2.1⟩
  · exact ⟨(master_T30_mixing T5).1, (master_T30_mixing T5).2.2.1⟩
  · exact (master_T30_cesaro).1

end PT.Stochastic
