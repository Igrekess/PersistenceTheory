# Persistence Theory

**Persistence Theory (PT) derives the Standard Model of particle physics from
a single arithmetic axiom (`s = 1/2`).** This repository contains the
monograph (FR + EN PDFs), the Lean 4 formalisation of the foundational
theorems, the Python scripts that reproduce every numerical result, and the
public website source.

---

## Status

**Preprint.** Not yet peer-reviewed. Single-author work.

The mathematics is stated with explicit epistemic tags throughout:
`[THM]` theorem, `[DER]` derivation, `[DER-PHYS]` physical derivation under
an interpretation bridge, `[VAL]` validation against data, `[PRED]`
falsifiable prediction, `[COND]` conditional result. Numerical agreement is
reported as evidence, not as proof that the underlying ontology is correct.

## Author's note

> I have worked hard for a year, aided by modern tools, to formalise a
> theory I have been thinking about ever since reading Mark A. Ludwig and
> seeing the Ulam spiral thirty years ago. I have tried to be as rigorous
> and scrupulous as possible while remaining pedagogical, within the limits
> of my modest means. It is with humility that I add my stone to the edifice
> of human understanding, without claiming it is perfectly carved. Above
> all, I hope that some of my avenues of thought may inspire others, and
> that my work will persist in one form or another.
>
> — Yan Senez

---

## What you can do here

| Goal | Where to go |
|---|---|
| **Read the monograph (English)** | [`TheTheoryOfPersistence.pdf`](TheTheoryOfPersistence.pdf) — 906 pages |
| **Lire la monographie (français)** | [`TheorieDeLaPersistance_FR.pdf`](TheorieDeLaPersistance_FR.pdf) — 939 pages |
| **Verify the formal proofs** | [`pt_lean/`](pt_lean) — Lean 4 + Mathlib, `lake build` |
| **Reproduce the numerical results** | [`scripts/`](scripts) — `python run_all.py` |
| **Read companion research notes** | [`research_notes/`](research_notes) — six standalone drafts |
| **Browse the public website source** | [`website/`](website) — Astro 5 + MDX |

---

## Key results at a glance

- **One input** (`s = 1/2`, forced by the antidiagonal transfer matrix of the
  Eratosthenes sieve), **zero adjusted parameters**, **43 Standard Model
  observables** derived.
- **Mean relative error 0.303 %**, median 0.059 %, across the full set of
  derived observables (couplings, masses, mixing matrices, electroweak and
  QCD quantities, cosmological constants).
- **10 main theorems kernel-verified in Lean** along the critical path
  `T1 → T3 → s=1/2 → T2 → L0 → T7 → W7-1`, plus the BA0–BA5 bridge layer
  and N1–N4 / G1–G3 supporting theorems.
- **K4 closure** (Higgs self-coupling `λ_H = 1/8`) reduced from `[CONJ]` to
  `[DER strict modulo standard NCG postulates]`, with the algebraic core
  fully formalised — see [`scripts/app_y_higgs_zeta_duality/`](scripts/app_y_higgs_zeta_duality)
  and [`pt_lean/PT/Bridge/`](pt_lean/PT/Bridge).
- **28 Tier-1 falsifiable predictions** with explicit experimental windows
  (HL-LHC ~2030, DUNE ~2032, Euclid ~2028, KATRIN, LEGEND, LZ/XENONnT,
  Einstein Telescope, Cosmic Explorer).

The canonical numerical chain is short enough to recite:

```text
s = 1/2  →  μ* = 15 = 3 + 5 + 7  →  active primes {3, 5, 7}
         →  γ_p, sin² θ_p (cyclic phase, holonomy)
         →  α_EM⁻¹ ≈ ∏ over active primes  ≈ 137.036
         →  sin² θ_W, m_W/m_Z, CKM/PMNS, m_H, ...
```

No constant in this chain is fitted to data. The fixed point `μ* = 15` is
the unique solution of the persistence-active sum on `μ ∈ [8, 20]`; the
active prime set `{3, 5, 7}` is forced by `γ_p > 1/2 ⇔ p ∈ {3, 5, 7}`.

---

## Repository structure

```text
.
├── TheTheoryOfPersistence.pdf       English monograph (905 pages)
├── TheorieDeLaPersistance_FR.pdf    French monograph (938 pages)
├── README.md                        This file
├── requirements.txt                 Top-level Python dependencies
├── pt_lean/                         Lean 4 formalisation (~3600 modules)
├── scripts/                         Python scripts, organised by chapter
├── research_notes/                  Six standalone working drafts
├── website/                         Astro 5 site source (FR + EN)
├── _images/                         Repository images
└── _archive_linguistique/           Historical linguistic notes (archive)
```

The `pt_lean/` tree is a self-contained Lean 4 package; `scripts/` is the
Python companion; `website/` and `research_notes/` are self-describing
through their own READMEs.

---

## Building the Lean library

The Lean 4 formalisation contains the 10 main theorems on the PT critical
path (`T1, T3, s=1/2, T2, L0, T6a/b/c, T7, W7-1`), the N1–N4 uniqueness
chain, the G1–G3 Fisher-geometry layer, the BA0–BA5 window-transport
bridge, and the newly added Bridge modules that close K4 (`λ_H = 1/8`).

Requires [`elan`](https://github.com/leanprover/elan) (the Lean toolchain
manager). The exact Lean version is pinned in `pt_lean/lean-toolchain`.

```sh
cd pt_lean
lake exe cache get      # download Mathlib oleans (~1 GB, ~5 min)
lake build              # compile the PT library
```

Expected outcome: `lake build` completes with **3599 jobs PASS**. A single
documented `sorry` remains in `PT/Information/G3FisherUniqueness.lean`
(Fisher metric uniqueness up to scale, a classical Čencov result the
monograph marks `\leanExternal`); every other module is sorry-free.

To check a specific theorem:

```sh
lake build PT.Sieve.T1ForbiddenTransitions
lake build PT.Stochastic.SHalf
lake build PT.FixedPoint.T7MuStar
lake build PT.Bridge.K4LambdaH
lake build PT.Analysis.W7SpiralIdentity
lake build PT                       # umbrella import
```

See [`pt_lean/README.md`](pt_lean/README.md) for the full per-module
status table.

---

## Running the scripts

The `scripts/` directory is organised by monograph chapter
(`ch01_sieve/`, `ch10_fine_structure/`, …, `ch22_chemistry/`,
`app_y_higgs_zeta_duality/`, …). Each chapter directory contains the
scripts that reproduce the figures, tables, and numerical identities
cited in that chapter.

```sh
python3 -m venv venv
source venv/bin/activate          # macOS / Linux
# venv\Scripts\activate           # Windows

pip install -r scripts/requirements.txt

cd scripts
python run_all.py                 # run the full chapter suite
python run_all.py --summary       # one-line status per chapter
python run_all.py ch10            # run a single chapter
python ch10_fine_structure/proof_alpha_EM.py
pytest -v                         # run the test suite
```

Reports land under `scripts/reports/` as CSV / JSON. Some quantum-gravity
and Kerr checks can ingest external LVK / GWTC posterior data when
`PT_LVK_REMNANTS_DIR` is set; without that dataset, the structural checks
still run and the empirical branch exits cleanly.

Python 3.10+, `sympy`, `mpmath`, `numpy`, `scipy`. See
[`scripts/README.md`](scripts/README.md) for the full chapter map and the
script-to-Lean correspondence.

---

## A note on "no fitted parameter"

The phrase **no fitted parameter** means: no continuously fitted parameter
in the canonical PT chain. It does **not** mean that there are no units,
no structural choices, no perturbative ingredients, or no domain
assumptions. Those are stated and audited inside the monograph (Ch. 23
audit, Annex E status ledger).

The stronger claim is the careful one:

> The same persistence principle appears to generate rigid mathematical
> structures whose physical and chemical readings reproduce many observed
> regularities, with explicit and inspectable status tags.

PT is interesting only if it remains rigid when the tests become harder.
This repository is built to make that rigidity inspectable.

---

## Predictions snapshot

A small selection of the 28 Tier-1 falsifiable predictions:

| ID | Claim | Type | Test window |
|---|---|---|---|
| P1 | neutrinos are Dirac | `[PRED]` | LEGEND / neutrinoless ββ decay |
| P3 | `θ_QCD = 0` exactly | `[PRED]` | axion / strong-CP searches |
| P4 | `δ_CP^{PMNS} = 197.358°` | `[PRED]` | DUNE (~2032) |
| P5 | `m_{ν3} ≈ 0.0505 eV` | `[PRED]` | KATRIN / cosmology |
| P10 | no QCD axion | `[NEG]` | ADMX / IAXO |
| P11 | no SUSY below 100 TeV | `[NEG]` | future colliders |
| P14 | no WIMPs | `[NEG]` | LZ / XENONnT |
| P15 | `α_GW < 10⁻³` | `[PRED]` | Einstein Telescope / Cosmic Explorer |
| P29 | Casimir relative shift `ΔP/P = 1.31 × 10⁻⁴` | `[PRED]` | precision Casimir |

Known-value agreements in the Standard Model, chemistry, and cosmology are
reported in the monograph as validations or explanations — they are
**not** automatically reclassified as new predictions.

---

## How to cite

Each main artefact has a permanent Zenodo DOI. Please cite the DOI that
matches what you actually use.

| Artefact | DOI |
|---|---|
| Monograph (FR + EN PDFs) | [`10.5281/zenodo.18726591`](https://doi.org/10.5281/zenodo.18726591) |
| Mathematics articles bundle | [`10.5281/zenodo.19443954`](https://doi.org/10.5281/zenodo.19443954) |
| SCS — Sieve Color Space | [`10.5281/zenodo.19458652`](https://doi.org/10.5281/zenodo.19458652) |

```bibtex
@book{senez2026persistence,
  author    = {Senez, Yan},
  title     = {The Theory of Persistence: From the Sieve to the Standard Model},
  year      = {2026},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.18726591},
  url       = {https://doi.org/10.5281/zenodo.18726591}
}

@misc{senez2026persistencemathematics,
  author    = {Senez, Yan},
  title     = {Persistence Theory --- Mathematics articles},
  year      = {2026},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.19443954},
  url       = {https://doi.org/10.5281/zenodo.19443954}
}

@software{senez2026persistencescs,
  author    = {Senez, Yan},
  title     = {SCS --- Sieve Color Space},
  year      = {2026},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.19458652},
  url       = {https://doi.org/10.5281/zenodo.19458652}
}
```

---

## Companion projects

Focused companion repositories with narrower scope:

- [PT-MATHEMATICS](https://github.com/Igrekess/PT_MATHEMATICS)
- [PT-PHYSICS](https://github.com/Igrekess/PT_PHYSICS)
- [PT-CHEMISTRY](https://github.com/Igrekess/PT_CHEMISTRY)
- [Simplex Color Space](https://github.com/Igrekess/SimplexColorSpace)

---

## Website

The Astro 5 site source lives in [`website/`](website). It is designed as
a public introduction to the theory: a "5-minute idea" essay, the
mathematical theorem chain, physics and chemistry pages, the predictions
ledger, and explicit status / limitations pages.

```sh
cd website
npm install
npm run dev
```

---

## How to contribute or report issues

This is single-author work and the corpus is large, so the most useful
contributions are:

- **Bug reports on the scripts** — open a [GitHub
  issue](https://github.com/Igrekess/PersistenceTheory/issues) with the
  command you ran, the platform, and the failing output.
- **Counter-examples** to derivations or predictions — please cite the
  monograph section and the exact identity you contest. Status tags
  (`[THM]` / `[DER]` / `[PRED]`) frame what counts as a refutation.
- **Lean improvements** — closing the remaining `G3FisherUniqueness`
  `sorry` or porting a result currently flagged `\leanExternal` is
  welcome; please open a PR against `pt_lean/`.
- **Editorial corrections** to the monograph — typos, broken
  cross-references, ambiguous wording. Use an issue with the chapter and
  section.

The repository is read-only-by-default; please coordinate via issues
before opening PRs that change physical claims.

---

## AI assistance

Persistence Theory is primarily the product of human research and
interpretation. Large language models — ChatGPT, Claude — were used as
tools for drafting, coding, checking scripts, explaining unfamiliar
mathematics, and accelerating editorial work. Their use is acknowledged
here as assistance, not as authorship. The author does not claim that
current LLMs understand the full theory; they manipulate, test, and
explain parts of it, but they do not replace the human responsibility for
the claims.

---

## License

- **Lean library** (`pt_lean/`): Apache 2.0, matching Mathlib's
  convention.
- **Python scripts** (`scripts/`): MIT.
- **Website source** (`website/`): see
  [`website/LICENSE`](website/LICENSE) for the code license and
  [`website/LICENSE-CONTENT`](website/LICENSE-CONTENT) for the prose
  license.
- **Monograph PDFs and research notes**: © Yan Senez, 2026. Released as
  preprints; please cite via the Zenodo DOIs above.

---

## Contact

Yan Senez — <yan.senez@gmail.com>

For technical questions on the scripts or the Lean formalisation, prefer
a GitHub issue so the discussion stays public and searchable.
