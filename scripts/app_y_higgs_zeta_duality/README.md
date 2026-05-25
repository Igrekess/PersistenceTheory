# PT non-commutative bifurcation — numerical scripts

These ten Python scripts produce the numerical results underlying the
Higgs self-coupling derivation `λ_H = 1/8` (conjecture K4 of the
monograph; see Appendix Y).

They are the source of the figures, tables, and exact rational
identities cited in the monograph's K4 closure section. Each script
is self-contained and reproducible; running them in numerical order
(23 through 32) recovers the full chain of intermediate results.

## Dependencies

Python 3.10+, with `sympy`, `mpmath`, `numpy`, and `scipy` (see
`PUBLIC/PersistenceTheory/requirements.txt` for exact versions).

## Script overview

| # | Script | Purpose |
|---|---|---|
| 23 | `23_K4_correction_search.py` | Search for structural corrections to the `R = 17/16` ratio (gap 0.04 %). |
| 24 | `24_E5b_NLO_correction_search.py` | Search for an NLO / one-loop / dressing correction explaining `ε = 3.95 × 10⁻⁴` between `λ_H(K4) = 2R/17` and `λ_H(PT) = 1/8`. |
| 25 | `25_E5a_hyperbolic_full.py` | Compute on the hyperbolic cuspidal geometry of Σ_pers (Z1, Z2 theorems). |
| 26 | `26_K4_cusp_identity.py` | Establish the cuspidal geometric identity `R_cusp = p₂²/(p₁·p₃) = 25/21`. |
| 27 | `27_E6a_seeley_dewitt_full.py` | Full Seeley–DeWitt expansion to `a_4` for the PT Dirac operator, with identification of the non-minimal Higgs–Ricci coupling `ξ_PT = 5/12`. |
| 28 | `28_E6c1_higgs_metric_setup.py` | Setup of the Donoghue–Connes Higgs metric for PT. |
| 29 | `29_E6c23_kappa_cusp_propagation.py` | Computation of the κ factor with cuspidal geometry; propagation of `R_cusp` into the Higgs normalisation. |
| 30 | `30_E6c4_uniqueness_and_synthesis.py` | Test of uniqueness for the canonical identification; final synthesis of K4 = `λ_H = 1/8`. |
| 31 | `31_E7_uniqueness_proof.py` | Three uniqueness arguments for the canonical identifications. |
| 32 | `32_E8_formal_uniqueness.py` | Formal uniqueness of `f(u) = exp(-u/N_b)` via CRT factorisation + Shore–Johnson G1 + Cauchy's functional equation. |

## Relation to the Lean formalisation

The mathematical content these scripts establish at the numerical
level is independently formalised (kernel-verified) in the Lean
library under `pt_lean/PT/Bridge/`:

* `CauchyMultiplicativeExp.lean` — formal Cauchy → exp (the
  mathematical core of script 32).
* `HiggsCutoffUniqueness.lean` — PT cutoff uniqueness (the
  synthesis of scripts 31, 32).
* `K4LambdaH.lean` — algebraic closure of `λ_H = 1/8` (the synthesis
  of scripts 30, 26, 27).
* `FiniteSpectralTriple.lean` — the finite spectral triple
  `ST_F = (ℂ², ℂ², m·σ_x, σ_x)` with kernel-verified
  `Tr_F(I_F) = 2 = N_b`.
* `ScaleFromFiniteSpectralTriple.lean` — structural identification of
  the cutoff scale with the dimension of `H_F`.
* `SpectralAction.lean` — scalar-Dirac spectral action.
* `HeatKernelPostulate.lean` — Connes–Marcolli identification of the
  PT cutoff with the heat kernel of `D²`.
* `PartitionFunctionFactorisation.lean`,
  `GibbsDistribution.lean`, `JaynesPrinciple.lean` —
  Gibbs–Jaynes infrastructure.

See the corresponding `.lean` modules for the formal statements and
the kernel-verified proofs.
