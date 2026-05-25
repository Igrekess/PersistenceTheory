# Numerical scripts for the Theory of Persistence

This directory contains the Python scripts that reproduce all
quantitative results cited in the Theory of Persistence monograph
(`TheTheoryOfPersistence.pdf` / `TheorieDeLaPersistance_FR.pdf`).

## Organisation

The directory is organised by **monograph chapter**, with each
chapter directory containing the scripts that produce the figures,
tables and numerical identities cited in that chapter.

### Chapters (`chXX_<topic>/`)

| Directory | Monograph chapter |
|---|---|
| `ch01_sieve/` | Ch. 1 — The Eratosthenes sieve and T1 |
| `ch02_uniqueness/` | Ch. 2 — Uniqueness of the sieve (N1–N4, T6) |
| `ch03_conservation/` | Ch. 3 — Conservation, T2, spectral content |
| `ch04_gft/` | Ch. 4 — Gap Fundamental Theorem (GFT) |
| `ch05_geometry/` | Ch. 5 — Fisher geometry, G1–G3 |
| `ch06_holonomy/` | Ch. 6 — Holonomy, `sin²θ_p`, T6 |
| `ch07_convergence/` | Ch. 7 — Spectral convergence T4 |
| `ch08_fixed_point/` | Ch. 8 — Fixed point `μ* = 15`, T5/T7 |
| `ch09_bridge/` | Ch. 9 — Window transport, bridge BA0–BA5 |
| `ch10_fine_structure/` | Ch. 10 — `α_EM` cascade, fine structure |
| `ch11_couplings/` | Ch. 11 — Coupling reconstruction, Kähler-Fisher |
| `ch12_quantum/` | Ch. 12 — Quantum reconstruction |
| `ch13_relativity/` | Ch. 13 — Relativity, cosmology, QG spin foam |
| `ch14_thermodynamics/` | Ch. 14 — Thermodynamics |
| `ch15_sm_observables/` | Ch. 15 — SM observables (43 derivations) |
| `ch16_perturbative/` | Ch. 16 — Perturbative QFT (NLO, NNLO) |
| `ch17_feynman/` | Ch. 17 — Feynman programme |
| `ch18_quarkonium/` | Ch. 18 — Quarkonium |
| `ch19_hadrons/` | Ch. 19 — Light hadrons |
| `ch20_collider/` | Ch. 20 — Collider observables |
| `ch20b_*` to `ch20g_*` | Ch. 20 sub-sections (Wilson coefficients, BSM taxonomy, DM, neutrinos, dark energy, Hubble tension, etc.) |
| `ch21_predictions/` | Ch. 21 — Predictions and falsifiability |
| `ch22_chemistry/` | Ch. 22 — Chemistry (atomic, molecular, condensed matter) |
| `ch23_audit/` | Ch. 23 — Audit (epistemic categories, test suite) |
| `ch25_BA0_closing/` | Ch. 25 — BA0 closing theorem T0 |
| `ch37b_RH_analysis/` | Ch. 37b — Riemann hypothesis analysis (RH cusp, Berry-Keating, Hilbert-Polya) |
| `ch_PM/` | PM algebra (Predictive Mathematics) |
| `ch_math_structures/` | Mathematical structures |

Several chapters have specialised sub-directories that group related
scripts (e.g. `ch13_relativity/cosmology/`,
`ch15_sm_observables/derivations/`, `ch37b_RH_analysis/`,
`ch_math_structures/projective_limit/`, etc.).

### Appendices (`app_<letter>_<topic>/`)

| Directory | Monograph appendix |
|---|---|
| `app_p_eml/` | Appendix P — Effective Mass Limit (EML) |
| `app_q_completeness/` | Appendix Q — Completeness arguments |
| `app_r_gauge/` | Appendix R — Gauge structure |
| `app_s_eml_closure/` | Appendix S — EML closure |
| `app_v_eynard_orantin/` | Appendix V — Eynard–Orantin table |
| `app_y_higgs_zeta_duality/` | Appendix Y — Higgs–zeta duality (K4 closure, `λ_H = 1/8`) |

### Infrastructure

| File / Directory | Purpose |
|---|---|
| `lib/` | Shared library (utility modules used across chapters) |
| `tests/` | Test runners (`verify_sota.py`, `verify_molecules.py`) |
| `reports/` | Computed output reports per chapter (CSV, JSON) |
| `pt_constants.py` | Shared physical constants |
| `conftest.py` | pytest configuration |
| `run_all.py` | Top-level test runner |
| `test_runner.py` | Auxiliary test runner |
| `requirements.txt` | Python dependencies |

## Dependencies

Python 3.10+ with `sympy`, `mpmath`, `numpy`, `scipy`. See
`requirements.txt` for exact versions.

## Running scripts

From the repository root:

```bash
cd PUBLIC/PersistenceTheory
python scripts/ch01_sieve/proof_T1_sieve.py
```

Or run the full test suite via:

```bash
python scripts/run_all.py
```

## Relation to the Lean library

Many chapter scripts have a corresponding kernel-verified Lean module
under `pt_lean/PT/`. For example:

* `scripts/ch01_sieve/proof_T1_sieve.py` ↔ `pt_lean/PT/Sieve/T1*.lean`
* `scripts/ch03_conservation/*` ↔ `pt_lean/PT/Stochastic/T2T3*.lean`
* `scripts/ch08_fixed_point/*` ↔ `pt_lean/PT/FixedPoint/T7*.lean`
* `scripts/app_y_higgs_zeta_duality/*` ↔ `pt_lean/PT/Bridge/{Cauchy,Higgs,K4,FiniteSpectralTriple,…}.lean`

See `pt_lean/README.md` for the Lean library overview.
