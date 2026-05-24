/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/

import PT.NashDeGiorgi.GammaAtMu
import PT.NashDeGiorgi.Z1ActiveSet
import PT.NashDeGiorgi.LemmaAFormal
import PT.NashDeGiorgi.MainTheorems

/-!
# PT.NashDeGiorgi — umbrella module

Formalisation of the four Nash + De Giorgi structural theorems of PT
(see `PT_PROJECTS/PT_NASH_DEGIORGI/` for the source notes and scripts).

* `GammaAtMu` — extended definitions `gamma_at p μ` for variable `μ ∈ ℚ`,
  plus rational verifications of the active-set table.
* `Z1ActiveSet` — the m_K formula and the eventually-exact convergence
  (combinatorial core), Theorem Z1.
* `MainTheorems` — collected statements of Z1, Z2a, Z3, Z4, BK-PT
  consolidated. The discrete cores are theorems; the functional-analytic
  parts are recorded as axioms with cross-references to the proofs in
  the corresponding notes.

## Status

Discrete/combinatorial cores: `[THM]` in Lean (proved here).
Functional-analytic statements: `[axiomatised]` with sketches in
`PT_NASH_DEGIORGI/notes/{07,08,11,12}.md`.

The cornerstone Lemma A (`γ_p(μ)` monotone in `μ`) is axiomatised
pending the analytical proof of note 03 §2.
-/
