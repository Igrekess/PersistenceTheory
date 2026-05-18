# Persistence Theory

Public repository for the monograph, website, and companion scripts of
**Persistence Theory**.

The monograph is a preprint. It has not yet been peer reviewed. The repository
is therefore meant to make the work readable, reproducible, and criticisable:
the mathematics is stated with explicit status tags, the numerical checks are
scripted, and the physical interpretation is presented as a hypothesis exposed
to future tests.

- English monograph: [TheTheoryOfPersistence.pdf](TheTheoryOfPersistence.pdf)
  (967 pages)
- French monograph: [TheorieDeLaPersistance_FR.pdf](TheorieDeLaPersistance_FR.pdf)
  (995 pages)
- Lean 4 formalisation: [`pt_lean/`](pt_lean) — 22 foundational
  theorems on the critical path T1 → T7 → W7-1 kernel-verified,
  spread across 8 modules, plus 171 additional secondary modules (179 modules in total).
- Research notes: [`research_notes/`](research_notes) — six
  standalone working drafts cited by the monograph (Hilbert–Pólya
  map, $p$-adic PT–Ramanujan, Casimir P29 prediction, Berry–Keating
  cusp, $A_{\mathrm{PT}}$ algebra, spectral model of zeros of $\zeta$).
- Website source: [`website/`](website)
- Companion scripts: [`scripts/`](scripts)

## The Starting Point

Persistence Theory does not begin with the Standard Model, the periodic table,
or even the primes.

It begins with a simpler question:

> When a system passes through constraints, what disappears, and what persists?

The shortest sentence is:

> What persists is what remains under constraint.

Here, "remaining" does not mean passive inertia. It means what survives
admissible filters without losing its structural identity. A pebble on a beach
is a useful first image: its shape was not chosen by the sea; it is what the
sea could not remove. Fragile edges disappeared, unstable roughness was taken
away, and a readable form remains.

PT proposes that this question is fundamental. A structure is not only present
or absent. It can dissipate, close, turn into noise, become an echo, or survive
successive filters. What survives is what can become stable, measurable,
transmissible, and eventually physical.

## The Fundamental Principle

In its most compact mathematical form, PT starts from the identity

```text
log2(m) = D_KL(P || U_m) + H(P).
```

For a non-specialist, the idea is this:

- a system has a total budget of possible distinctions;
- part of that budget becomes recognizable structure;
- the rest remains dispersed as entropy.

`D_KL` measures the structured part: what persists. `H` measures the dispersed
part: what remains entropic. The sum is conserved. This is why the identity
plays the role of the **fundamental principle of persistence**.

The rest of the theory asks what happens when this principle is made rigid,
computable, and repeatedly tested.

## Why The Sieve Appears

The sieve of Eratosthenes is not assumed to be "the world". It is used as the
minimal mathematical laboratory where persistence under constraint can be made
exact.

The usual school reflex is to look at what the sieve removes. PT asks the
opposite question: what remains, and why does it remain?

In that setting, the theory studies survivors, gaps, residue classes, cyclic
phases, fixed points, and anomalous dimensions. The discrete objects are read
as stable traces of a continuous mechanics under constraint. The sieve is thus
not an ontology by itself; it is the simplest place where the principle becomes
calculable.

## What The Theory Then Tests

Once the persistence principle has a computable form, PT asks whether the same
logic appears in domains that look unrelated at first.

### Mathematics

PT studies the arithmetic of survivors:

- prime gaps and circular gap structure;
- CRT decomposition and cyclic phases;
- the forced spin symmetry `s = 1/2`;
- the fixed point `mu* = 15`;
- Fisher geometry, anomalous dimensions, and holonomy;
- theorem-level and script-backed checks for the internal mathematical chain.

The mathematical layer is the place where the strongest unconditional claims
belong. It is also where the theory defines its vocabulary of persistence,
entropy, echo, closure, and bridge.

### Physics

PT then asks whether the structures that persist in the mathematical layer can
reconstruct known physical objects:

- coupling constants;
- masses and mixing matrices;
- electroweak and QCD observables;
- time, Fisher geometry, relativity, and quantum gravity constraints;
- cosmological quantities and negative predictions.

The physical layer is deliberately status-tagged. Some statements are
theorems inside the PT framework. Others are derivations under stated
assumptions, bridge identifications, validations against data, or predictions.
This distinction is essential: numerical agreement is evidence, not proof that
the ontology is correct.

### Chemistry

PT also tests the same persistence logic in chemistry:

- the periodic table and period lengths;
- shell/channel structure;
- ionization energies;
- electron affinities;
- molecular and condensed-matter checks;
- nuclear shell and magic-number modules.

The goal is not to replace all of quantum chemistry with a single toy formula.
It is to ask whether the same constrained persistence structure explains why
major chemical regularities have the shapes they do.

## What Is In This Repository

```text
.
├── TheTheoryOfPersistence.pdf          English monograph
├── TheorieDeLaPersistance_FR.pdf       French monograph
├── README.md                           This file
├── requirements.txt                    Top-level Python dependencies
├── scripts/                            Companion verification scripts
├── pt_lean/                            Lean 4 formalisation
├── research_notes/                     Standalone working drafts
├── website/                            Astro website source
└── _images/                            Repository images
```

The `scripts/` directory contains the companion checks used to audit the
monograph. They are organised by chapter and write machine-readable reports
under `scripts/reports/`.

The `pt_lean/` directory is a Lean 4 + Mathlib package. The 22
foundational theorems on the PT critical path
(T1 → T3 → `s=1/2` → T2 → L0 → T7 → W7-1) are kernel-verified
without `sorry`; about 150 additional modules cover downstream
material. Build with `lake build` after installing Lean via
`elan` (the `lean-toolchain` file pins the exact version). See
[`pt_lean/README.md`](pt_lean/README.md) for the per-module status
table.

The `research_notes/` directory collects six standalone drafts that
the monograph cites as companion artifacts (notably for the
Hilbert–Pólya map of Part III and for the Casimir P29 prediction).
None are published yet; each is a working document, status-tagged
on its own. See [`research_notes/README.md`](research_notes/README.md).

The public repository currently includes the clean script tree for the current
monograph, including arithmetic, physics, chemistry, audit, and registered
QG/Kerr checks. Large external datasets, exploratory notebooks, and early
research scratch files are not vendored here.

## Reproduce The Checks

```bash
git clone https://github.com/Igrekess/PersistenceTheory.git
cd PersistenceTheory

python3 -m venv venv
source venv/bin/activate       # macOS / Linux
# venv\Scripts\activate        # Windows

pip install -r scripts/requirements.txt

cd scripts
python run_all.py --summary
```

Useful commands:

```bash
python run_all.py --tree
python run_all.py ch10
python ch10_fine_structure/proof_alpha_EM.py
pytest -v
```

Some quantum-gravity / Kerr tests can use external LVK/GWTC posterior data when
`PT_LVK_REMNANTS_DIR` is set. Without that dataset, the structural checks still
run and the empirical branch exits cleanly.

## Current Epistemic Posture

PT must be read with its status labels:

- `[THM]`: theorem inside the mathematical framework;
- `[IDENTITY]`: exact identity or conservation statement;
- `[DER]`: derivation under stated assumptions;
- `[DER-PHYS]`: physical derivation requiring an interpretation bridge;
- `[VAL]`: validation against known data;
- `[PRED]`: falsifiable prediction;
- `[COND]`: conditional result;
- `[META]`: audit or methodological statement.

The phrase "no fitted parameter" means **no continuously fitted parameter** in
the canonical chain. It does not mean that there are no units, no structural
choices, no perturbative ingredients, or no domain assumptions. Those choices
must be stated and audited.

The theory does not claim that physical reality has already been proved to be
the arithmetic sieve. The stronger claim is more careful:

> the same persistence principle appears to generate rigid mathematical
> structures whose physical and chemical readings reproduce many observed
> regularities.

That is a serious hypothesis, not a completed consensus theory.

## Prediction Snapshot

The monograph lists named predictions and negative predictions. Their role is
to expose PT to failure, not to protect it.

Examples:

| ID | Claim | Type | Test window |
|---|---|---|---|
| P1 | neutrinos are Dirac | PRED | LEGEND / neutrinoless beta decay |
| P3 | `theta_QCD = 0` exactly | PRED | axion / strong CP searches |
| P4 | `delta_CP^PMNS = 197.358 deg` | PRED | DUNE |
| P5 | `m_nu3 ~= 0.0505 eV` | PRED | KATRIN / cosmology |
| P10 | no QCD axion | NEG | ADMX / IAXO |
| P11 | no SUSY below 100 TeV | NEG | future colliders |
| P14 | no WIMPs | NEG | LZ / XENONnT |
| P15 | `alpha_GW < 10^-3` | PRED | Einstein Telescope / Cosmic Explorer |

Known-value agreements in the Standard Model, chemistry, and cosmology are
reported in the monograph as validations or explanations, not automatically as
new predictions.

## Citation

Each main artefact has a permanent Zenodo DOI. Please use the DOI corresponding
to what you actually cite.

| Artefact | DOI |
|---|---|
| Monograph (FR + EN PDFs) | [`10.5281/zenodo.18726591`](https://doi.org/10.5281/zenodo.18726591) |
| Mathematics articles | [`10.5281/zenodo.19443954`](https://doi.org/10.5281/zenodo.19443954) |
| SCS — Sieve Color Space | [`10.5281/zenodo.19458652`](https://doi.org/10.5281/zenodo.19458652) |

Recommended BibTeX entries:

```bibtex
@book{senez2026persistencemonograph,
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

## Companion Projects

Focused companion repositories:

- [PT-MATHEMATICS](https://github.com/Igrekess/PT_MATHEMATICS)
- [PT-PHYSICS](https://github.com/Igrekess/PT_PHYSICS)
- [PT-CHEMISTRY](https://github.com/Igrekess/PT_CHEMISTRY)
- [Simplex Color Space](https://github.com/Igrekess/SimplexColorSpace)

## Website

The website lives in [`website/`](website). It is designed as a public
introduction to the theory:

- "The idea in 5 minutes" for the principle of persistence;
- mathematical pages for the theorem chain;
- physics pages for time, relativity, quantum gravity, cosmology, and
  predictions;
- chemistry pages for the periodic table, ionization energies, and electron
  affinities;
- status and limitations pages to keep the epistemic posture explicit.

Local development:

```bash
cd website
npm install
npm run dev
```

## AI Assistance

The theory is primarily the product of human research and interpretation.
Large language models, including ChatGPT and Claude, were used as tools for
drafting, coding, checking scripts, explaining unfamiliar mathematics, and
accelerating editorial work. Their use is acknowledged here as assistance, not
as authorship.

The author does not claim that current LLMs understand the full theory. They
can help manipulate, test, and explain parts of it; they do not replace the
human responsibility for the claims.

## How To Read This Repository

Start with the idea, then follow the evidence:

1. Read the principle in the website essay:
   `website/src/content/essays/en/idea-5min.mdx`.
2. Read the monograph introduction and status ledger.
3. Check the mathematical core.
4. Inspect the physical and chemical derivations with their status tags.
5. Run the scripts and compare the generated reports.
6. Evaluate the predictions as future failure points.

PT is interesting only if it remains rigid when the tests become harder. This
repository is built to make that rigidity inspectable.
