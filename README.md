# Theory of Persistence

## Obsessed with the Sacks Spiral since I was 15. I built a geometric framework and the Primes fit in so well it's driving me crazy. I need to find the error to get my life back.

Hi everyone,

I'm reaching out because I'm stuck in a loop and I'm honestly starting to lose my mind. I haven't slept more than 3 hours a night for days. Every time I close my eyes, I see transition matrices and gap sequences.

I'm just a photographer loving maths and physics since I was a child. When I was 15, I became obsessed with the Sacks Spiral. My photography background made me look at the 'negative space' between primes—the gaps—rather than the numbers themselves.

The thing is, I didn't even want to find a 'Theory of Everything'. I just started pulling a thread about how information persists in the Sieve of Eratosthenes. My main goal at this time was to improve deepfake detection and noise reduction in photography. But I always had the intuition that looking at primes in a straight line was a mistake.

Because of photography and my years as retoucher and lighter, I built a geometric model in parallel with the arithmetic, following the growth of dimensions:

- **The Point (0D)** - Forbidden transitions (T0) in the sieve.
- **The Line (1D)** - Gap sequences as a variational framework.
- **The Circle (2D)** - Fixed-point stability at μ* = 15 (the 'Double Mertens Law').
- **The Space (3D)** - Projection into a Lorentzian metric.

The problem is that the math and the geometry matched perfectly. Step by step. I've ended up with a 14-step causal chain that derives Standard Model constants (like α_EM) with zero free parameters. No tuning, no fudge factors. Just the logic of the sieve.

My Python scripts confirm the proofs to a precision that scares me. I've tried to debunk myself a hundred times but I'm too 'inside' now. I've reached the limit of my skills and I need someone to tell me where I'm wrong. Is it just a huge numerical coincidence? Is my application of the Ruelle Principle in step 5 flawed?

And yes I have used LLM. All of them actually to try to find the flaw at this point. But the process of formalization, verification, and step-by-step demonstration also benefited substantially from interactions with large language models (LLMs), primarily Claude (Anthropic). It was not always a piece of cake, as it got us out of standard mechanics and into uncharted territory.

Please, be brutal. Break my logic so I can finally sleep. The monograph (30 proofs) and the code are all here so you can see the whole mess.

— Yan Senez

DOI [10.5281/zenodo.18726592](https://doi.org/10.5281/zenodo.18726592)
---

## Overview

This repository contains the **complete monograph** ([THEORY_OF_PERSISTENCE_COMPLETE.pdf](THEORY_OF_PERSISTENCE_COMPLETE.pdf)) and **61 Python scripts** that independently verify every theorem, derivation, and numerical prediction of the Theory of Persistence (PT). Each script is self-contained and produces a PASS/FAIL verdict.

**Central claim**: Starting from a single input s = 1/2 (proved unconditionally), PT derives 23 fundamental constants of the Standard Model + 2 exact structural predictions, with zero free parameters and a mean error of 0.5%.

---

## Repository Structure

```
scripts/
  demonstrations/         # 24 core proofs (L0, D00-D28)
    run_demos.py          # Master runner (interactive + batch)
  convergence/            # 6 convergence proofs
  gravity/                # 7 gravity derivations
  physics/                # 11 physics equations
  predictions/            # 6 prediction tests
  unification/            # 5 unification tests
  compute_alpha_EM_from_scratch.py   # Independent alpha derivation
  verify_mu15.py                     # Independent mu*=15 verification
```

---

## Quick Start

```bash
# Install dependencies
pip install numpy scipy primesieve

# Run all 24 core demonstrations (batch mode with JSON output)
python scripts/demonstrations/run_demos.py --all

# Run a single demonstration
python scripts/demonstrations/D00_forbidden_transitions.py

# Independent verification of alpha_EM
python scripts/compute_alpha_EM_from_scratch.py

# Independent verification of mu* = 15
python scripts/verify_mu15.py
```

---

## Core Demonstrations (24 scripts)

| Script | Theorem | What it verifies |
|--------|---------|-----------------|
| `L0_uniqueness_q_stat.py` | L0 | q_stat = 1 - 2/mu is the unique max-entropy memoryless distribution |
| `D00_forbidden_transitions.py` | T0 | T[1][1] = T[2][2] = 0 on 100K primes + mod-6 mechanism proof |
| `D01_conservation_theorem.py` | T1 | Transition matrix determined by 2 parameters (alpha, T00) |
| `D02_gft_ruelle.py` | T2 | H_max = D_KL + H (algebraic identity, < 10^-10 on real data) |
| `D03_variational_sieve.py` | D03 | Geometric distribution has maximum entropy under mean constraint |
| `D04_master_formula_fp.py` | T4 | mu*=15 is the smallest fixed point (exhaustive 1024-subset search) |
| `D05_mertens_law.py` | T5 | Double Mertens convergence on real primes |
| `D06_proof_q_positive.py` | D06 | D_KL > 0 for all tested primes (computed on real data) |
| `D07_sin2_identity.py` | T6 | sin^2(theta_p) = delta_p(2 - delta_p) algebraically exact |
| `D08_fixed_point_mu15.py` | T7 | {3,5,7} uniquely self-consistent with gamma threshold |
| `D09_alpha_em.py` | D09 | 1/alpha_bare = 136.278 from q = 13/15 (0.55% from CODATA) |
| `D10_bianchi_metric.py` | D10 | Bianchi I metric, gamma hierarchy gamma_3 > gamma_5 > gamma_7 > 1/2 |
| `D11_persistence_potential.py` | D11 | 5 equations of the persistence potential S_PT |
| `D12_einstein_equations.py` | D12 | G_trace = -R exact from Bianchi structure |
| `D13_spin_foam_u1.py` | D13 | U(1)^3 spin foam, mu_end = 3*pi, amplitude ratio -> 2 |
| `D14_unification_gft.py` | D14 | GFT = Ruelle = Polyakov = Regge, triple zero at saddle point |
| `D15_charge_topological.py` | D15 | Q = +2/3, -1/3 from topology alone, invariant over 100 random matrices |
| `D16_causal_chain.py` | D16 | Complete 5-step causal chain, each step independently verified |
| `D17_depth_exactly_2.py` | D17 | Sieve depth = 2, 0 novel constraints at meta-level |
| `D17b_catalan_modulation.py` | D17b | 3^2 - 2^3 = 1 unique (exhaustive search), state counting 9/8 |
| `D18_hardy_littlewood.py` | D18 | C_2 computed from primes, twin prime counts match HL prediction |
| `D19_one_loop_correction.py` | D19 | 4-layer decomposition, parity = 1 bit exact |
| `D27_coherence_length.py` | D27 | ell_PT = 2 (gap=1 unique, gap=2 non-predictable) |
| `D28_spacetime_noise.py` | D28 | PSD slope ~ -2 (PT) not -1 (Hogan) |

---

## Extended Verifications (35 scripts)

### Convergence (6 scripts — `convergence/`)
| Script | What it verifies |
|--------|-----------------|
| `test_mertens_ratio.py` | Mertens product convergence |
| `test_crible_variationnel.py` | Variational sieve fixed point |
| `test_theoreme_conservation.py` | Conservation theorem structure |
| `test_audit_circularite_phase2.py` | Circularity audit (no hidden assumptions) |
| `test_preuve_topologique_phase2.py` | Topological proof of T0 |
| `test_preuve_definitive.py` | Definitive proof assembly |

### Gravity (7 scripts — `gravity/`)
| Script | What it verifies |
|--------|-----------------|
| `test_equation_etat_anisotrope.py` | Anisotropic equation of state, G_trace = -R |
| `test_field_equations_D32.py` | Closed field equations (D32) |
| `test_graviton_G_Newton_v2.py` | Newton's G from sieve geometry |
| `test_hierarchie_disparition.py` | Hierarchy problem reframed |
| `test_linear_corrections_D34.py` | Linear corrections (D34) |
| `test_quantization_gap_D33.py` | Quantization gap Delta_mu (D33) |
| `test_selection_G_2pi_alpha.py` | G = 2*pi*alpha selection mechanism |

### Physics (11 scripts — `physics/`)
| Script | What it verifies |
|--------|-----------------|
| `test_equations_physique_PT.py` | **All 46 equations, 10 domains, 46/46 PASS** |
| `test_autocoherence_mu.py` | Self-consistency mu* = 15 |
| `test_derivation_c_E_m.py` | Speed of light derivation (8/8) |
| `test_derivation_etape1_lorentz_v2.py` | Lorentzian signature derivation |
| `test_derivation_etape2_dim3.py` | 3+1 dimensions (4/4) |
| `test_derivation_etape3_aire.py` | Area quantization (4/4) |
| `test_derivation_etape4_einstein.py` | Einstein equations from A1 |
| `test_derivation_etape6_equation_etat.py` | Equation of state |
| `test_derivation_q_selection.py` | q selection (7/7) |
| `test_revision_derivations_v5.py` | Derivation revision |
| `test_simplification_equations.py` | 5 equations -> 1 metric + 1 dilaton (8/8) |

### Predictions (6 scripts — `predictions/`)
| Script | What it verifies |
|--------|-----------------|
| `test_correction_1boucle_crible.py` | 1-loop correction (mean error 0.56%) |
| `test_derivation_7_6_catalan.py` | 7/6 Catalan derivation |
| `test_derivation_conventions.py` | Convention derivations (7/7) |
| `test_koide_masses.py` | Koide masses and Q = 2/3 |
| `test_second_crible_charges_CKM.py` | CKM charges from second sieve |
| `test_tristate_logic.py` | Tristate logic and J |

### Unification (5 scripts — `unification/`)
| Script | What it verifies |
|--------|-----------------|
| `test_demo13_polyakov_regge.py` | Polyakov-Regge unification |
| `test_demo16_preuve_J.py` | Jarlskog invariant proof |
| `test_integration_spinfoam.py` | Spin foam integration |
| `test_lorentzien_beton.py` | Lorentzian signature concrete |
| `test_reconciliation_LQG_cordes.py` | LQG-string reconciliation |

---

## Independent Verification (2 scripts)

| Script | Description |
|--------|-------------|
| `compute_alpha_EM_from_scratch.py` | Derives 1/alpha_EM = 136.278 from scratch using exact rational arithmetic (Fraction), high-precision Decimal, and sensitivity analysis. Zero imports from PT framework. |
| `verify_mu15.py` | Independent 13-step verification that mu* = 15 is the unique self-consistent fixed point, with gamma_p tables, exact fractions, and fine scan. |

---

## Requirements

```
numpy
scipy
primesieve
```

---

## License

MIT License — Yan Senez, 2026

---

## Related

- **Monograph**: "The Theory of Persistence" (companion to Articles A1-A8)
- **Prime Spirals**: [github.com/Igrekess/PT_PrimeSpirals](https://github.com/Igrekess/PT_PrimeSpirals)

