# Plan — Site web Théorie de la Persistance

**Créé le 2026-04-26. Validé 2026-04-26.**

## Décisions utilisateur

| # | Question | Décision |
|---|----------|----------|
| 1 | Domaine | **`persistencetheory.org`** |
| 2 | Hosting | **GitHub Pages** |
| 3 | Priorité | **V1.0 catalogue scientifique** d'abord (saute V0.5) |
| 4 | Branding | aucun pour l'instant |
| 5 | Visu V0.5 | crible animé (différé puisque V1.0 prioritaire) |
| 6 | Articles thématiques | **sélection** des plus solides (8–10) |
| 7 | Pyodide | OK |
| 8 | Mailing-list | différée à V2+, mais **architecture extensible** dès V0.1 |
| 9 | Peer review ouverte | non en V1 |
| 10 | Licences | **MIT** pour le code, **CC BY 4.0** pour le contenu et les PDFs |

---

## 1. Vision

Un seul site, deux objectifs simultanés :

1. **Vulgarisation** — donner à n'importe quel curieux un accès clair, esthétique et progressif à l'idée centrale (« 43 observables du Modèle Standard dérivés d'un unique input s = 1/2 »).
2. **Vérification scientifique** — donner à n'importe quel chercheur un accès direct aux preuves complètes, aux scripts de vérification, et aux articles thématiques sélectionnés.

Le site **ne remplace pas** la monographie : il sert de **portail** vers elle, en l'incarnant sous une forme web navigable, interactive et bilingue.

**Bilingue FR/EN** dès le départ (les deux PDFs 909 pages existent déjà).

---

## 2. Audiences et parcours

Trois couches concentriques, chacune satisfaisable indépendamment :

| Couche | Audience | Temps d'engagement | Sortie attendue |
|--------|----------|-------------------|-----------------|
| **L1 — Découverte** | grand public, journalistes, étudiants | 2–10 min | « Je comprends ce que prétend la théorie » |
| **L2 — Compréhension** | étudiants avancés, scientifiques curieux | 30 min – 2 h | « Je sais comment la chaîne logique fonctionne » |
| **L3 — Vérification** | chercheurs, sceptiques, peer reviewers | journées à semaines | « Je peux reproduire chaque calcul et juger sur pièces » |

**Priorité validée : V1.0 (couche L3)** d'abord. La couche L1/L2 (V0.5) sera ajoutée par-dessus la fondation scientifique.

---

## 3. Architecture informationnelle

```
/                          Landing (FR par défaut)
/en                        Landing (EN)
├── /monograph             L3 — viewer PDF + download EN/FR
├── /theorems              L3 — index 8 théorèmes fondateurs
│   ├── T0-BA0-closing
│   ├── L0-uniqueness
│   ├── T1-forbidden-transitions
│   ├── T2-spectral-conservation
│   ├── T3-antidiagonal
│   ├── T4-convergence
│   ├── T5-fixed-point
│   ├── T6-holonomy
│   └── GFT-fundamental
├── /observables           L3 — table 43 observables filtrable
│   └── /[id]              dérivation + lien script
├── /articles              L3 — index 8-10 articles thématiques sélectionnés
│   └── /[slug]            résumé + PDF embed + lien GitHub
├── /scripts               L3 — Pyodide runner + index scripts
├── /audit                 L3 — audit épistémique transparent
├── /downloads             tout télécharger
├── /contact               email
└── /about                 auteur, licences, statut peer review
```

Couche L1/L2 (V0.5+) à ajouter ultérieurement : `/idea`, `/chain`, `/calculator`, `/glossary`, `/faq`.

---

## 4. Pile technique

| Couche | Choix | Raison |
|--------|-------|--------|
| **Framework** | Astro 5 | statique-first, MDX, *islands*, i18n natif, GitHub Pages friendly |
| **Styles** | Tailwind CSS | bundle minuscule + design system rapide |
| **Contenu** | MDX | Markdown + composants interactifs |
| **Math** | KaTeX (rehype) | rendu serveur, léger |
| **Code** | Shiki | coloration syntaxique server-side |
| **Visualisations** | D3.js (SVG) + Plotly | écosystème mature |
| **Recherche** | Pagefind | index full-text statique, bilingue |
| **Python in browser** | Pyodide | exécution scripts PT navigateur |
| **PDF viewer** | iframe + PDF.js | embed sans dépendance |
| **i18n** | Astro 5 natif (`/`, `/en/`) | FR par défaut |
| **Hosting** | GitHub Pages | CI gratuit, certif auto, deploy auto |
| **Domaine** | `persistencetheory.org` (CNAME) | DNS à configurer |
| **CI** | GitHub Actions | build + deploy à chaque push `main` |

---

## 5. Inventaire du contenu

### Source de vérité (déjà existant)

| Source | Localisation | Usage site |
|--------|--------------|------------|
| Monographie EN (909 p.) | `PUBLIC/PersistenceTheory/TheTheoryOfPersistence.pdf` | viewer + download |
| Monographie FR (909 p.) | `PUBLIC/PersistenceTheory/TheorieDeLaPersistance_FR.pdf` | viewer + download |
| Sources LaTeX | `PT_MONOGRAPHY/chapters/` (EN) + `chapters_fr/` (FR) | extraction MDX par chapitre |
| Synthèse vulgarisée | `SYNTHESE_PT_ACCESSIBLE.md` | base L1/L2 (V0.5) |
| Scripts vérification | `PUBLIC/PersistenceTheory/scripts/` (59 scripts, 2610 checks) | Pyodide runner |
| Articles thématiques | `PT_ARTICLES/PT_HL`, `PT_ICP`, `PT_GEO`, ..., 30+ projets | sélection 8-10 + PDF embed |
| Sous-projets publics | `PUBLIC/PT_CHEMISTRY`, `PT_PHYSICS`, `PT_MATHEMATICS`, `SieveColorSpace` | vitrines + démos |
| Articles fondateurs A1–A8 | `Articles_Fondateurs/` | théorèmes |
| Glossaire | `frontmatter/glossary.tex` (EN) + `frontmatter_fr/glossary.tex` (FR) | parsé en JSON |
| Notation | `frontmatter/notation.tex` + `frontmatter_fr/notation.tex` | tableaux symboles |

### Articles thématiques candidats à sélection (8-10 / 30+)

À discuter ensemble. Premiers candidats par solidité scientifique :

- `PT_HL` — Hardy-Littlewood singular series (article standalone)
- `PT_GRH` — GRH from sieve mixing
- `PT_ICP` — ICP brainweb (validation médicale)
- `PT_GEO` — Geophysics (predictions sismiques validées 2026-04-20 Japon M7.4)
- `PT_PROTEIN_ALLOSTERY` — biologie (validation 2026-04-18)
- `PT_BRIDGE` — Lemmas E/F/G (cœur du pont arithmétique→physique)
- `PT_LAMBDA` — cosmological constant
- `PT_CHESS` — application crible aux échecs (vulgarisation)

---

## 6. Phases de livraison (ordre validé)

### V0.1 — Squelette navigable (en cours)

- Init Astro 5 + Tailwind + MDX + i18n FR/EN
- Landing minimale FR/EN
- Page `/monograph` viewer PDF + download
- Page `/about` (licences, contact)
- GitHub Action deploy → `persistencetheory.org`
- **Livrable** : site en ligne, lecture PDFs OK

### V1.0 — Catalogue scientifique (priorité 1)

- 8 pages théorèmes (T0–T7 + GFT) extraites des chapitres LaTeX
- Table 43 observables filtrable par secteur
- 8-10 articles thématiques sélectionnés avec embed PDF
- Index scripts + Pyodide runner pour 5 scripts emblématiques
- Page audit transparent
- Pagefind search bilingue
- **Livrable** : un chercheur peut consulter le cœur scientifique entièrement

### V0.5 — Vulgarisation (priorité 2)

- 5 essais L1/L2 rédigés bilingues
- Visualisation crible mod-3 (animation D3.js/SVG)
- Diagramme chaîne causale (SVG cliquable)
- Glossaire popover
- FAQ 28 entrées
- **Livrable** : un curieux comprend la théorie en lisant le site

### V2.0 — Interactivité

- Calculatrices PT (γ_p, sin² θ_p, α_EM)
- Bouton profondeur Vulgarisé/Standard/Technique
- Sous-projets publics intégrés (PT_CHEMISTRY launcher)
- **Mailing-list** opérationnelle (techno à choisir : Buttondown, Listmonk self-host, ou fichier JSON statique avec Cloudflare Worker)
- **Livrable** : site outil

### V3.0 — Polish

- Visualisations 3D
- Tutoriels vidéo
- Notebook Jupyter Lab embed

---

## 7. Licences

| Élément | Licence |
|---------|---------|
| Code site (Astro, JS, configs) | **MIT** |
| Contenu rédactionnel (essais, FAQ, glossaire web) | **CC BY 4.0** |
| Monographie PDF | **CC BY 4.0** |
| Articles thématiques PDFs | **CC BY 4.0** |
| Scripts Python (PUBLIC/PersistenceTheory/scripts) | **MIT** (déjà en place) |

Fichiers à créer : `LICENSE` (MIT pour code), `LICENSE-CONTENT` (CC BY 4.0 pour textes), `CITATION.cff` au repo racine.

---

## 8. Architecture extensible mailing-list (Q8)

Pour permettre l'ajout futur d'une mailing-list **sans refactor**, V0.1 pose ces fondations :

- Form abstraction : composant `<Subscribe />` aujourd'hui inactif (renvoie vers contact email), demain pluggable sur n'importe quel backend
- Schema Zod pour Email + consentement RGPD
- Variable d'environnement `MAILING_PROVIDER` (vide en V0.1, à fixer en V2)

---

## 9. Risques et points d'attention

- **Maintenance** : update monographie → site (extraction LaTeX → MDX automatisable mais à mettre en place)
- **Cohérence FR/EN** : drift FR↔EN à surveiller (CI compare labels)
- **Performance Pyodide** : ~10 MB WASM = chargement lent en mobile ; lazy-load + cache
- **SEO** : `hreflang` correct dès V0.1
- **Accessibilité** : prose mathématique = lecture difficile ; viser WCAG AA
- **Domain DNS** : configuration `persistencetheory.org` à faire (CNAME → `igrekess.github.io` ou A records)

---

## 10. État actuel

- ✅ PLAN.md rédigé et validé
- 🚧 V0.1 en cours d'init (Astro skeleton)
