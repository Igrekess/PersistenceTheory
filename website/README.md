# persistencetheory.org — site web

Site web de la **Théorie de la Persistance** (PT). Bilingue FR/EN. Statique-first.

| | |
|---|---|
| **Statut** | V0.1 (squelette navigable) |
| **Domaine** | persistencetheory.org |
| **Hosting** | GitHub Pages |
| **Stack** | Astro 5 + Tailwind + MDX + KaTeX |
| **Licences** | MIT (code) + CC BY 4.0 (contenu) |

## Démarrage local

Prérequis : Node.js ≥ 20.

```bash
npm install
npm run dev      # http://localhost:4321
npm run build    # produit ./dist
npm run preview  # serveur de production local
```

## Structure

```
PT-WEBSITE/
├── PLAN.md                 plan détaillé du site (V0.1 → V3.0)
├── astro.config.mjs        config Astro + i18n FR/EN
├── tailwind.config.mjs     palette épistémique (thmblue, derorange…)
├── src/
│   ├── i18n/strings.ts     dictionnaire UI bilingue
│   ├── layouts/            layouts Astro
│   ├── pages/              FR par défaut (/), EN sous /en/
│   └── styles/global.css   base Tailwind + classes pt-tag, pt-card
├── public/
│   ├── CNAME               persistencetheory.org
│   └── favicon.svg
└── .github/workflows/deploy.yml   CI → GitHub Pages
```

## Routing i18n

| URL FR | URL EN |
|--------|--------|
| `/` | `/en/` |
| `/monograph` | `/en/monograph` |
| `/theorems` | `/en/theorems` |
| `/observables` | `/en/observables` |
| `/articles` | `/en/articles` |
| `/scripts` | `/en/scripts` |
| `/about` | `/en/about` |

## Roadmap

Voir [PLAN.md](./PLAN.md) pour les phases V0.1 → V3.0.

## Licences

- **Code** (Astro, JS, TS, configs) : MIT — voir [LICENSE](./LICENSE)
- **Contenu** (textes, PDFs) : CC BY 4.0 — voir [LICENSE-CONTENT](./LICENSE-CONTENT)
