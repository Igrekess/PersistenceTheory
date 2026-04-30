# PROMPT — Site web PT, V2 (interactivité) + V0.5 (vulgarisation)

**À coller dans une nouvelle session Claude Code.**

---

## Contexte (à connaître avant tout)

Tu travailles sur le **site web de la Théorie de la Persistance** (PT) à
`/Volumes/PT-YS-0326/LA THEORIE DE LA PERSITANCE/PT-WEBSITE/`.

Le site est conçu **bilingue FR/EN** (FR par défaut, EN sous `/en/`),
**statique-first** (Astro 5 + Tailwind + MDX + KaTeX + Pagefind),
**hosting GitHub Pages** sur le domaine **`persistencetheory.org`**, mais
**pas encore déployé** — on travaille en local.

Le **PLAN.md** à la racine du dossier détaille les phases V0.1 → V3.0 et
les décisions de l'utilisateur (notamment : MIT pour le code, CC BY 4.0
pour le contenu, pas de mailing-list pour l'instant mais architecture
extensible).

### État actuel (à la fin de la session précédente)

| Phase | Contenu | Statut |
|-------|---------|--------|
| V0.1 | Squelette Astro, layout bilingue, /monograph viewer, /about | ✅ |
| V1.0 | 9 théorèmes (T0-T6 + GFT + L0), 43 observables filtrable, 8 articles, audit transparent, Pagefind bilingue | ✅ |
| V0.5 | Vulgarisation (5 essais, crible animé, FAQ, glossaire popover) | ⏳ À faire |
| V2.0 | Interactivité (calculatrices, Pyodide, bouton profondeur, mailing) | ⏳ À faire |

**Build actuel** : 42 pages, 2 279 mots indexés Pagefind, ~7 s de build.
9 routes critiques toutes 200 (smoke test).

### Architecture en place (à NE PAS recasser)

```
PT-WEBSITE/
├── PLAN.md                          ← spec complète
├── package.json                     ← Astro 5 + Tailwind + MDX + KaTeX + Pagefind
├── astro.config.mjs                 ← i18n FR (défaut) + EN sous /en/
├── tailwind.config.mjs              ← palette épistémique (thmblue, derorange…)
├── tsconfig.json
├── LICENSE / LICENSE-CONTENT        ← MIT + CC BY 4.0
├── .github/workflows/deploy.yml.disabled  ← deploy DÉSACTIVÉ (à laisser)
├── public/
│   ├── CNAME                        ← persistencetheory.org
│   └── favicon.svg
├── src/
│   ├── i18n/strings.ts              ← UI bilingue, getLangFromUrl, alternateLangPath
│   ├── content/
│   │   ├── config.ts                ← collections schemas (theorems, articles)
│   │   ├── theorems/                ← 18 .mdx (9 × FR + EN)
│   │   └── articles/                ← 16 .mdx (8 × FR + EN)
│   ├── data/observables.json        ← 43 observables
│   ├── components/
│   │   ├── Tag.astro                ← badge épistémique (default import !)
│   │   ├── ObservablesTable.astro   ← table filtrable client-side
│   │   └── Search.astro             ← Pagefind avec filtre lang
│   ├── layouts/BaseLayout.astro     ← header + nav + lang switch + Search + footer + hreflang
│   ├── pages/                       ← 21 pages FR
│   │   ├── index.astro              landing
│   │   ├── monograph.astro          PDF iframe
│   │   ├── about.astro
│   │   ├── audit.astro
│   │   ├── observables.astro        table 43
│   │   ├── scripts.astro            stub (V2 = Pyodide)
│   │   ├── theorems/index.astro     liste
│   │   ├── theorems/[id].astro      page dynamique
│   │   ├── articles/index.astro     liste
│   │   └── articles/[slug].astro    page dynamique
│   └── pages/en/                    ← 21 pages EN miroirs
└── src/styles/global.css            ← Tailwind base + classes pt-tag, pt-card
```

### Conventions strictes

1. **Bilingue obligatoire dès le départ** : tout ce que tu crées doit avoir une
   version FR (chemin `src/pages/foo.astro`) ET une version EN
   (`src/pages/en/foo.astro`).
2. **i18n via `getLangFromUrl(Astro.url)`** dans BaseLayout ; les chaînes UI
   sont dans `src/i18n/strings.ts` (ajoute-les là si tu introduis un nouveau
   label réutilisable).
3. **Composants Astro = default import** (`import Tag from '@components/Tag.astro'`).
   Pas d'import nommé. Si tu débogues une erreur "X is not exported", c'est ça.
4. **MDX YAML frontmatter en single-quoted** (`'...'`) si tu mets du LaTeX :
   évite les YAML escapes `\\` cassés. Voir `src/content/theorems/T0.fr.mdx`
   pour modèle.
5. **Maths inline dans `.astro` casse** parce que `{` est interprété par JSX.
   Solution : ou bien `<sub>`/`<sup>` HTML simples, ou bien convertis la page
   en `.mdx` (qui supporte LaTeX via remark-math + rehype-katex).
6. **Aliases TS** : `@/`, `@layouts/`, `@components/`, `@content/` (voir
   `tsconfig.json`).
7. **Drift labels = 0** dans la monographie sous-jacente : si tu réutilises
   du contenu de `chapters/` ou `chapters_fr/`, garde les `\label{...}`
   exactement.
8. **Ne déploie pas**, ne pousse pas. Le workflow GitHub Pages est désactivé
   (`.github/workflows/deploy.yml.disabled`). On fait tout en local.
9. **Commit local OK** dans `PT-WEBSITE/` (pas encore initialisé en git, à
   faire si nécessaire), mais pas de push.

### Commandes (à lancer SÉPARÉMENT)

```bash
cd "/Volumes/PT-YS-0326/LA THEORIE DE LA PERSITANCE/PT-WEBSITE"
npm run dev       # http://localhost:4321 (sans search)
npm run build     # build complet + Pagefind
npm run preview   # http://localhost:4321 avec search
```

⚠️ Ne jamais coller les 3 sur la même ligne avec `#` (les commentaires zsh
cassent l'enchaînement).

---

## V2 — Interactivité (priorité 1, suit V1)

### V2.1 — Calculatrices PT (1–2 j)

Trois calculatrices interactives en JavaScript natif (pas besoin de Pyodide).

**Localisation** : `/calculator` FR + `/en/calculator` EN.

**Contenu** :

1. **Calculatrice γ_p(μ)** : input μ ∈ [3, 100], slider ; output table
   γ_3, γ_5, γ_7, γ_11, γ_13 + indicateur "actif" si γ_p > 1/2 ;
   point fixe trouvé si Σ_actif = μ.
   Formule (déjà dans la monographie) :
   ```
   q = 1 - 2/μ
   δ_p = (1 - q^p) / p
   γ_p = 4 p q^(p-1) (1 - δ_p) / [μ (1 - q^p) (2 - δ_p)]
   ```
2. **Calculatrice sin²(θ_p) sur les deux branches** : input p, μ ;
   output sin²(θ_p, q_+) et sin²(θ_p, q_−), où q_+ = 1-2/μ et
   q_− = exp(-1/μ). Formule : sin²(θ_p) = δ_p (2 - δ_p).
3. **Calculatrice α_EM "bare → dressed"** : à μ = 15, calcule
   α_bare = ∏_{p ∈ {3,5,7}} sin²(θ_p, q_+) ≈ 1/136.28, puis applique
   le facteur d'habillage (R51, voir ch10_fine_structure.tex de la
   monographie : Δ_1 + Δ_2 + Δ_3 + ghost VP) → 1/137.036.

**Tech** : composant Astro avec `<script>` inline (pas besoin de framework).
Pour les sliders/charts : SVG natif ou `<input type="range">`. Si chart
plus riche (γ_p en fonction de μ), Plotly.js via CDN.

**Référentiel** : les formules canoniques sont dans
`PT_MONOGRAPHY/chapters/ch06_holonomy.tex` lignes 305-403 et
`PT_MONOGRAPHY/chapters/ch10_fine_structure.tex`.

**Validation** : à μ = 15 les calculatrices doivent reproduire les valeurs
de `chapters/ch06_holonomy.tex` ligne 401-403 :
```
γ_3 = 4536129 / 5616704     = 0.80761...
γ_5 = 486792684365 / 699097512194  = 0.69632...
γ_7 = 2827519972576117 / 4748396022746468 = 0.59547...
```

### V2.2 — Pyodide runner pour 5 scripts emblématiques (2–3 j)

**Localisation** : `/scripts` (déjà stub) → enrichir avec runner.

**Sélection des 5 scripts** parmi `PUBLIC/PersistenceTheory/scripts/` :
- `ch08_fixed_point/proof_fixed_point.py` — démontre μ* = 15
- `ch06_holonomy/proof_holonomy.py` — vérifie sin² = δ(2−δ)
- `ch10_fine_structure/proof_alpha_dressed.py` — α_EM = 1/137.036
- `ch_math_structures/test_math_tools.py` — outils mathématiques
- `ch09_bridge/proof_bridge_axioms.py` — pont arithmétique→physique
  (113 PASS, déjà refactoré en Phase 4 audit ChatGPT)

**Tech** :
- Pyodide via CDN `https://cdn.jsdelivr.net/pyodide/v0.27.0/full/pyodide.js`
- Lazy-load (~10 MB WASM) au clic sur "Lancer le script"
- Charge le contenu du script depuis `PUBLIC/PersistenceTheory/scripts/...`
  (à servir comme assets statiques OU à inliner)
- Affiche stdout/stderr et compte les PASS/FAIL en live

**UX** : bouton "Exécuter" (warning : 10 MB), output console temps réel,
résultat final en grand (X/X PASS).

**Caveat** : la lib `lib/pt_check.py` doit être chargée aussi. Solution :
inliner les 2 fichiers comme strings dans le composant Astro, exécuter
ensemble.

### V2.3 — Bouton profondeur Vulgarisé / Standard / Technique (1 j)

Sur les pages théorèmes et observables, ajouter un toggle 3 niveaux qui
cache/affiche des sections selon `data-depth="L1|L2|L3"`. Stocké dans
localStorage.

Pages à enrichir : `theorems/[id]`, `observables`, `articles/[slug]`,
plus tard les essais V0.5.

### V2.4 — Sous-projets publics intégrés (1 j)

Ajouter sur la home et `/articles` des cards pour :
- **PT_CHEMISTRY** (PTC simulator, MAE 2,28%) — `PUBLIC/PT_CHEMISTRY/`
- **PT_PHYSICS** (`PUBLIC/PT_PHYSICS/dashboard.html` existe déjà)
- **PT_MATHEMATICS** (M1-M5, `PUBLIC/PT_MATHEMATICS/`)
- **SieveColorSpace** (`PUBLIC/SieveColorSpace/`)

Soit en cards renvoyant vers GitHub, soit en intégrant les launchers Python
(Pyodide pour les démos courtes).

### V2.5 — Architecture mailing-list (≈30 min, prépare V3)

Composant `<Subscribe />` aujourd'hui inactif :
- Form HTML avec champ email + checkbox consentement RGPD
- Variable `MAILING_PROVIDER` dans `.env` (vide en V2)
- Si vide : redirige vers `mailto:yan.senez@protonmail.com?subject=PT%20site`
- Si défini : POST vers le provider choisi (Buttondown, Listmonk, etc.)

À placer dans le footer du `BaseLayout.astro`.

---

## V0.5 — Vulgarisation (priorité 2, après V2)

### V0.5.1 — 5 essais L1/L2 bilingues (3–5 j)

Localisation : `/essays/[slug]` FR + `/en/essays/[slug]` EN.

Collection MDX `essays` à créer (similaire à `articles`).

Liste des 5 :

1. **L'idée en 5 minutes** (slug `idea-5min`)
   Source : reformuler `SYNTHESE_PT_ACCESSIBLE.md`. Ton vulgarisé, viser
   un lecteur curieux sans physique avancée. ≈ 800 mots.
2. **Pourquoi 3 dimensions ?** (slug `why-3-dimensions`)
   Réponse PT : car {3,5,7} sont les trois premiers actifs au point fixe.
   Croisement physique : Bianchi I aux 3 directions actives.
3. **D'où vient α_EM = 1/137 ?** (slug `where-alpha-comes-from`)
   Le produit ∏ sin²(θ_p, q_+) sur {3,5,7} donne 1/136,28, puis l'habillage
   donne 1/137,036. Visuel : 3 facteurs successifs avec leurs valeurs.
4. **Qu'est-ce que la persistance ?** (slug `what-is-persistence`)
   D_KL en bits, structure persistante vs bruit, GFT comme conservation.
5. **Pourquoi 3 générations de fermions ?** (slug `why-3-generations`)
   N_gen = N_c = 3 = nombre de premiers actifs. Le crible "ne peut pas"
   en avoir plus.

**Conventions de prose** :
- Niveau L1 (vulgarisé) en premier 1/3 ; L2 (plus technique) après ; option
  L3 dépliable (lien vers le théorème exact).
- Ton humain, opinions assumées (rappel : skill humanizer disponible —
  `/Users/yan/.claude/skills/humanizer/SKILL.md`).
- Bilingue impératif.

### V0.5.2 — Visualisation crible mod-3 animée (1–2 j)

Sur la home (en dessous de l'accroche), ajouter une animation SVG
interactive : nombres entiers défilant, élimination par 2, par 3, par 5…
puis passage en mod 3 et apparition des "transitions interdites" (T1).

**Tech** : SVG + petit script vanilla. ≈ 200 lignes.

**Référentiel inspiration** : voir `Articles_Fondateurs/A1_transitions/`
si présent, ou `JOURNAL_RECHERCHE_PERSISTANCE.md`.

### V0.5.3 — Diagramme chaîne causale interactif (1 j)

Localisation : `/chain` FR + `/en/chain` EN.

SVG cliquable : T0 → L0 → T1 → … → 43 observables. Chaque maillon ouvre
un popover ou redirige vers le théorème correspondant.

**Source figure** : `chapters/ch_PM.tex` ou `frontmatter/status_ledger.tex`
(le tikzpicture de chaîne critique).

### V0.5.4 — Glossaire popover (1 j)

Survol d'un terme PT (μ*, q_+, q_−, sin² θ_p, etc.) dans le texte = popover
définition. Source : `frontmatter/glossary.tex` + `frontmatter_fr/glossary.tex`.

**Tech** : parser les entrées en JSON au build time, composant `<Term term="μ*">`
qui wrappe avec un `<abbr title="...">`.

### V0.5.5 — FAQ 28 entrées (2 j)

Localisation : `/faq` FR + `/en/faq` EN.

Sources :
- 5 erreurs de lecture de l'audit ChatGPT (déjà documentées dans
  `audit.astro`)
- Sections "Ce que la théorie ne dit pas" du skill `PT-Expert`
  (`/Users/yan/.claude/skills/PT-Expert/SKILL.md`)
- Questions FAQ classiques (numérologie, peer review, statut épistémique,
  comparaison Modèle Standard, etc.)

Format : `<details>` HTML natif avec accordéon Tailwind.

---

## Ordre d'attaque recommandé

**Ne pas tout faire d'un coup.** Découpe en mini-livrables livrables un par un :

1. V2.1 calculatrices (vrai gain UX rapide)
2. V0.5.1 essai « L'idée en 5 min » (vulgarisation = haut impact public)
3. V0.5.2 visu crible animé
4. V2.2 Pyodide (ambitieux, à isoler)
5. V0.5.3 chaîne causale + V0.5.4 glossaire popover
6. V2.3 bouton profondeur (refactor sur tout le site)
7. V0.5.5 FAQ
8. V2.4 sous-projets publics
9. V2.5 architecture mailing-list

À chaque palier : `npm run build` + `npm run preview` + smoke test des
nouvelles routes. Pas de commit ni push (le repo n'est pas encore en git ;
demander à l'utilisateur s'il veut initialiser).

---

## Fichiers clés à lire pour s'orienter

| Pour comprendre… | Lire… |
|------------------|-------|
| Le plan global et les décisions | `PT-WEBSITE/PLAN.md` |
| L'i18n et les chaînes UI | `PT-WEBSITE/src/i18n/strings.ts` |
| Le layout de base | `PT-WEBSITE/src/layouts/BaseLayout.astro` |
| La structure des théorèmes | `PT-WEBSITE/src/content/theorems/T0.fr.mdx` |
| La table observables | `PT-WEBSITE/src/components/ObservablesTable.astro` |
| La recherche Pagefind | `PT-WEBSITE/src/components/Search.astro` |
| Le contenu source PT | `PT_MONOGRAPHY/chapters/` (EN), `chapters_fr/` (FR), `Articles_Fondateurs/`, `SYNTHESE_PT_ACCESSIBLE.md` |
| Le skill PT-Expert | `/Users/yan/.claude/skills/PT-Expert/SKILL.md` (active automatiquement) |
| Le skill humanizer | `/Users/yan/.claude/skills/humanizer/SKILL.md` (à appeler explicitement pour relire la prose vulgarisée) |

---

## Premiers gestes en arrivant

```bash
cd "/Volumes/PT-YS-0326/LA THEORIE DE LA PERSITANCE/PT-WEBSITE"
ls                          # voir structure
cat PLAN.md | head -80      # rappeler les décisions
npm run build               # confirmer que la base existe et compile
npm run preview             # site live sur :4321
```

Puis demander à l'utilisateur lequel des chantiers V2 / V0.5 il veut
attaquer en premier.

---

**Statut au démarrage de la prochaine session** : V0.1 + V1.0 livrées en
local, drift 0, 42 pages, 9 routes 200, Pagefind bilingue OK, build 7 s.
Pas de déploiement, pas de git, pas de push. Tout en local.
