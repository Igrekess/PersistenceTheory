# The Theory of Persistence

This repository is the public companion to the monograph
[Senez, Y. (2026). *The Theory of Persistence: A Complete Monograph (2026)*](https://zenodo.org/records/19520809).

The theory is currently released as a preprint and has not been peer reviewed yet.

It gathers the monograph PDF together with the companion scripts used to
reproduce, audit, and verify every numerical claim developed in the text.

Persistence Theory starts from **zero fitted parameters** and derives:

- an arithmetic core built from sieve residues, forbidden transitions, conservation, GFT, and holonomy
- a closure layer with convergence, the fixed point `mu* = 15`, and predictive mathematics
- a systematic reconstruction of physical structure from arithmetic
- reconstructions and validations across quantum theory, relativity, thermodynamics, Standard Model observables, collider phenomenology, and chemistry

`s = 1/2` is **derived** as Theorem T1 (not assumed). `m_e` is a dimensional translation factor (like `c = 3e8 m/s`), not a parameter.

## From Theory To Derived Repositories

This repository presents the full theoretical framework. It also serves as the common source behind several more focused public repositories:

- [*PT-MATHEMATICS (PTM)*](https://github.com/Igrekess/PT_MATHEMATICS): the mathematical core, structures, and theorem-oriented explorations
- [*PT-PHYSICS (PTP)*](https://github.com/Igrekess/PT_PHYSICS): the physical reconstruction and validation program, including a zero-parameter calculator and dashboard for particle, nuclear, collider, and heavy-ion observables
- [*PT-CHEMISTRY (PTC)*](https://github.com/Igrekess/PT_CHEMISTRY): the chemistry branch, including a computational chemistry calculator and interactive app for atoms, molecules, reactions, spectroscopy, electrochemistry, and materials
- [*Simplex Color Space (SCS)*](https://github.com/Igrekess/SimplexColorSpace): a perceptual color theory and color space derived from Persistence Theory

<p align="center">
  <img src="_images/archimedes_spectral_PT_500k_600px.png" alt="500,000 primes on the Archimedean spiral coloured by Hardy-Littlewood density." />
</p>

<p align="center">
  500,000 primes on the Archimedean spiral, coloured by Hardy-Littlewood density C(f).<br>
  Every colour value is computed entirely from <code>s = 1/2</code> (Theorem T1), with zero free parameters.
</p>

## Companion Scripts

The `scripts/` directory contains **41 self-contained Python scripts** organized by monograph chapter, each proving or verifying one theorem or demonstration. Together they execute **2,410 checks with 0 failures** in under 60 seconds.

### Structure

```
scripts/
  pt_constants.py              Central constants (all derived from s=1/2)
  run_all.py                   Master test runner + report generator
  conftest.py                  pytest integration
  lib/
    pt_check.py                Test framework (Checker + progress bar + JSON reports)
    _primes.py                 Prime generation (primesieve + pure-Python fallback)

  Part I -- Arithmetic Core (6 scripts, 352 checks)
    ch01_sieve/                proof_T1_sieve.py
    ch02_uniqueness/           proof_T6_uniqueness.py
    ch03_conservation/         proof_conservation.py
    ch04_gft/                  proof_GFT.py
    ch05_geometry/             proof_G1_G5.py
    ch06_holonomy/             proof_holonomy.py

  Part II -- PT Closure (5 scripts, 611 checks)
    ch07_convergence/          proof_T4_convergence.py
    ch08_fixed_point/          proof_T5_fixed_point.py
    ch25_BA0_closing/          proof_T0_closing.py
    ch_PM/                     test_PM_algebra.py
    ch_math_structures/        test_math_tools.py, test_complex_mechanics.py

  Part III -- Bridge to Physics (3 scripts, 140 checks)
    ch09_bridge/               proof_bridge_axioms.py
    ch10_fine_structure/       proof_alpha_EM.py
    ch11_couplings/            proof_couplings.py

  Part IV -- Physical Reconstruction (4 scripts, 359 checks)
    ch12_quantum/              proof_quantum.py
    ch13_relativity/           proof_relativity.py
    ch14_thermodynamics/       proof_thermodynamics.py
    ch15_sm_observables/       test_43_observables.py

  Part V -- Validation & Dynamics (6 scripts, 358 checks)
    ch16_perturbative/         test_NLO_NNLO.py
    ch17_feynman/              test_feynman_programme.py
    ch18_quarkonium/           test_quarkonium.py
    ch19_hadrons/              test_hadrons.py
    ch20_collider/             test_collider.py
    ch21_predictions/          test_predictions.py

  Part VI -- Chemistry (14 scripts, 576 checks)
    ch22_chemistry/            proof_C1_periodic_table.py
                               proof_C2_shannon_cap.py
                               proof_C3_permutation.py
                               proof_C7_transfer_matrix.py
                               proof_C8_band_gap.py
                               proof_C9_tier_cascade.py
                               proof_C12_completude.py
                               proof_C13_three_centers.py
                               proof_C14_thermo_limit.py
                               test_IE_atomic.py
                               test_condensed_matter.py
                               test_electrochemistry.py
                               test_molecular.py
                               test_nuclear.py

  Part VII -- Audit & SOTA (2 scripts, 14 checks)
    ch23_audit/                test_audit.py
    verify_sota/               verify_sota.py (requires ptc package)
```

### Running the scripts

```bash
# Install dependencies
cd scripts
pip install -r requirements.txt

# Run the full suite (41 scripts, ~55s)
python run_all.py

# Run a single chapter
python run_all.py ch10

# Show the script tree
python run_all.py --tree

# Regenerate summary CSV from existing JSON reports
python run_all.py --summary
```

Each script can also be run standalone:

```bash
python ch10_fine_structure/proof_alpha_EM.py
```

### Output

Every script produces:
- A **structured JSON report** in `reports/<chapter>/<script>.json` with all check results and computed values
- Console output with `[PASS]`/`[FAIL]` for each check, progress bars for long computations, and a final `BILAN` summary

The master runner `run_all.py` additionally generates:
- `reports/summary.csv` -- one row per script (chapter, checks, pass/fail, duration)
- `reports/values.csv` -- all numerical results (label, computed value, expected, error%, unit)
- A complete readable report with per-chapter tables and numerical results

### pytest integration

```bash
cd scripts
pytest -v              # run all via subprocess isolation
pytest -k "ch15"       # run a single chapter
```

## Prediction Snapshot

The monograph organizes the framework's falsifiable output into 28 named predictions (`P1`--`P28`):

### Structural predictions

| ID | Prediction | Type | Main test / status |
|---|---|---|---|
| P1 | Neutrinos are Dirac | PRED | LEGEND 2030 |
| P2 | Normal hierarchy | EXPL | JUNO 2027 |
| P3 | `theta_QCD = 0` exactly | PRED | ADMX ongoing |
| P4 | `delta_CP^PMNS = 197.358 deg` | PRED | DUNE 2032 |

### Numerical predictions

| ID | Prediction | Type | Main test / status |
|---|---|---|---|
| P5 | `m_nu3 = 0.0505 eV` | PRED | KATRIN 2027 |
| P6 | `sin^2(theta_23) = 0.5733` | EXPL | DUNE 2032 |
| P7 | `sin^2(theta_13) = 0.02222` | EXPL | JUNO 2027 |
| P8 | `lambda_H = 0.1295` | PRED | HL-LHC 2030 |
| P9 | `H_0 = 67.41 km/s/Mpc` | EXPL | Euclid 2030 |
| P15 | `alpha_GW < 10^-3` | PRED | Einstein Telescope 2035 |
| P16 | `n_s = 0.964` | EXPL | CMB-S4 ~2030 |
| P17 | No BSM muon `g-2` anomaly | PRED | Fermilab/J-PARC ~2026 |
| P18 | `m_W = 80.3635 GeV` | EXPL | ATLAS/CMS precision updates |
| P20 | `w_0 = -1.003`, `w_a ~ 0` | PRED | DESI/Euclid 2028-2030 |

### Negative predictions

| ID | Prediction | Type | Main test / status |
|---|---|---|---|
| P10 | No axion (QCD) | NEG | ADMX/IAXO ~2030 |
| P11 | No SUSY below `100 TeV` | NEG | FCC-hh ~2040 |
| P12 | `N_gen = 3` exactly | EXPL | Colliders ongoing |
| P13 | No proton decay | NEG | Hyper-K ~2028 |
| P14 | No WIMPs | NEG | LZ/XENONnT ~2028 |
| P19 | No extra spatial dimensions | NEG | LHC Run 3 / future colliders |

`PRED` = genuine prediction, `EXPL` = explanatory post-diction, `NEG` = negative prediction.

## Benchmark Snapshot

| Branch | Scope | Benchmark |
|---|---|---|
| Monograph scripts | Theorems T0-T6, GFT, G1-G5, 43 SM observables, chemistry C1-C14 | **41 scripts, 2,410 checks, 0 failures, 55s** |
| `PTP` core physics | Standard Model, nuclear, exotic observables | 57 scored, mean error `0.30%`, 0 fitted parameters |
| `PTP` collider | LEP/SLC/LHC + PT shower | 72/72 within `5%`, 55/72 within `1%` |
| `PTC` chemistry | Atoms, molecules, spectroscopy, electrochemistry, materials | ~1200 observables; IE `0.056%`, atomization `2.17%` |
| `SCS` color | Perceptual color space | COMBVD `r = 0.893` vs `0.878` for CIEDE2000 |

## Key Results

| Observable | PT value | PDG / experiment | Error |
|---|---|---|---|
| `1/alpha_EM` | 137.035999 | 137.035999084 | 0.004 ppb |
| `sin^2(theta_W)` | 0.23119 | 0.23121 | 0.010% |
| `alpha_s(m_Z)` | 0.11806 | 0.11800 | 0.048% |
| `m_W` | 80.3635 GeV | 80.369 GeV | 0.007% |
| `m_Z` | 91.1878 GeV | 91.1876 GeV | 0.0002% |
| `G_F` | 1.16638e-5 | 1.16638e-5 | 0.0001% |
| 43 SM observables | -- | -- | mean 0.29%, median 0.05% |

Zero fitted parameters. Zero inputs. `s = 1/2` is derived (Theorem T1).

## Contents

- [Monograph PDF](TheTheoryOfPersistence.pdf)
- [`scripts/`](scripts): 41 companion scripts organized by chapter (see Appendix F)
- [`scripts/run_all.py`](scripts/run_all.py): master test runner and report generator
- [`scripts/pt_constants.py`](scripts/pt_constants.py): all derived constants from `s = 1/2`
- [`requirements.txt`](requirements.txt): top-level dependencies

## Related Repositories

- Mathematics: [*PT-MATHEMATICS (PTM)*](https://github.com/Igrekess/PT_MATHEMATICS)
- Physics: [*PT-PHYSICS (PTP)*](https://github.com/Igrekess/PT_PHYSICS)
- Chemistry: [*PT-CHEMISTRY (PTC)*](https://github.com/Igrekess/PT_CHEMISTRY)
- Color: [*Simplex Color Space (SCS)*](https://github.com/Igrekess/SimplexColorSpace)
- Monograph: [Senez, Y. (2026). *The Theory of Persistence*](https://zenodo.org/records/19520809)

## Notes

- This public repository is meant to be readable, navigable, and reproducible.
- The monograph remains the authoritative narrative source.
- The companion scripts reproduce every numerical value and every theorem verification claimed in the text.
- Each script produces a structured JSON report for automated cross-checking.
- Large language models (LLMs) were used as research and writing assistants to support parts of the drafting process and accelerate the development of companion scripts.
