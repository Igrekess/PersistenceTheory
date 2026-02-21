# Demonstrations — Théorie de la Persistance / Persistence Theory

---

## Prérequis / Prerequisites

- **Python 3.8+**
- Dépendances :
  ```bash
  pip install primesieve numpy scipy
  ```

---

## Cloner le dépôt / Clone the repo

```bash
git clone https://github.com/Igrekess/PersistenceTheory.git
cd PersistenceTheory
pip install primesieve numpy scipy
```

---

## Lancer les démonstrations / Run the demonstrations

### Menu interactif / Interactive menu

```bash
python scripts/demonstrations/run_demos.py
```

Affiche un menu numéroté. Saisissez le numéro d'une démonstration pour
l'exécuter, `q` pour quitter. La sortie est affichée en temps réel avec
le statut PASS/FAIL à la fin.

Displays a numbered menu. Enter a number to run the corresponding
demonstration, `q` to quit. Output is streamed in real time with
PASS/FAIL status shown at the end.

### Batch complet / Full batch

```bash
python scripts/demonstrations/run_demos.py --all
```

Exécute les 24 scripts séquentiellement, affiche un tableau récapitulatif
et écrit les résultats dans `results/demo_results_YYYYMMDD_HHMMSS.json`.

Runs all 24 scripts sequentially, prints a summary table, and writes
results to `results/demo_results_YYYYMMDD_HHMMSS.json`.

### Chemin JSON personnalisé / Custom JSON path

```bash
python scripts/demonstrations/run_demos.py --all --json results/run.json
```

---

## Liste des 24 scripts / List of 24 scripts

| # | ID    | Script                        | Titre FR                              | Title EN                            |
|---|-------|-------------------------------|---------------------------------------|-------------------------------------|
| 0 | L0    | L0_uniqueness_q_stat.py       | Unicité de q_stat                     | Uniqueness of q_stat                |
| 1 | D00   | D00_forbidden_transitions.py  | Transitions interdites mod 3          | Forbidden transitions mod 3         |
| 2 | D01   | D01_conservation_theorem.py   | Théorème de conservation              | Conservation theorem                |
| 3 | D02   | D02_gft_ruelle.py             | GFT et opérateur de Ruelle            | GFT and Ruelle operator             |
| 4 | D03   | D03_variational_sieve.py      | Sieve variationnel                    | Variational sieve                   |
| 5 | D04   | D04_master_formula_fp.py      | Formule maître et point fixe          | Master formula and fixed point      |
| 6 | D05   | D05_mertens_law.py            | Loi de Mertens                        | Mertens law                         |
| 7 | D06   | D06_proof_q_positive.py       | Preuve q positif                      | Proof q positive                    |
| 8 | D07   | D07_sin2_identity.py          | Identité sin²                         | sin² identity                       |
| 9 | D08   | D08_fixed_point_mu15.py       | Point fixe μ* = 15                    | Fixed point μ* = 15                 |
|10 | D09   | D09_alpha_em.py               | Constante de structure fine α_EM      | Fine structure constant α_EM        |
|11 | D10   | D10_bianchi_metric.py         | Métrique de Bianchi                   | Bianchi metric                      |
|12 | D11   | D11_persistence_potential.py  | Potentiel de persistance              | Persistence potential               |
|13 | D12   | D12_einstein_equations.py     | Équations d'Einstein                  | Einstein equations                  |
|14 | D13   | D13_spin_foam_u1.py           | Spin foam U(1)                        | Spin foam U(1)                      |
|15 | D14   | D14_unification_gft.py        | Unification GFT                       | GFT unification                     |
|16 | D15   | D15_charge_topological.py     | Charge topologique                    | Topological charge                  |
|17 | D16   | D16_causal_chain.py           | Chaîne causale                        | Causal chain                        |
|18 | D17   | D17_depth_exactly_2.py        | Profondeur exactement 2               | Depth exactly 2                     |
|19 | D17b  | D17b_catalan_modulation.py    | Modulation de Catalan                 | Catalan modulation                  |
|20 | D18   | D18_hardy_littlewood.py       | Hardy-Littlewood et persistance       | Hardy-Littlewood and persistence    |
|21 | D19   | D19_one_loop_correction.py    | Correction à une boucle              | One-loop correction                 |
|22 | D27   | D27_coherence_length.py       | Longueur de cohérence                 | Coherence length                    |
|23 | D28   | D28_spacetime_noise.py        | Bruit de fond spatio-temporel PT      | Spacetime background noise PT       |

---

## Structure des résultats / Results structure

Le fichier JSON généré par `--all` a la structure suivante :

```json
{
  "run_date": "2026-02-19T14:30:00",
  "python": "3.11.0",
  "total": 24,
  "passed": 22,
  "failed": 2,
  "results": [
    {
      "id": "D08",
      "title_fr": "Point fixe mu* = 15",
      "title_en": "Fixed point mu* = 15",
      "script": "D08_fixed_point_mu15.py",
      "status": "PASS",
      "return_code": 0,
      "duration_s": 0.87,
      "output": "...",
      "error": ""
    }
  ]
}
```

Les fichiers JSON sont écrits dans `results/` (créé automatiquement à la
racine du dépôt). Chaque statut est l'un de : `PASS`, `FAIL`, `ERROR`,
`TIMEOUT`, `SKIP`.

---

Yan Senez — Février 2026
