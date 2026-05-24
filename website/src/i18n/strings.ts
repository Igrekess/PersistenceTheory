/**
 * Localized UI strings.
 *
 * Page-level content lives in MDX files under src/content/.
 * This file is for navigation, buttons, labels, footer.
 */

export const languages = {
  fr: 'Français',
  en: 'English',
} as const;

export type Lang = keyof typeof languages;

export const defaultLang: Lang = 'fr';

export const ui = {
  fr: {
    // Nav
    'nav.monograph': 'Monographie',
    'nav.start': 'Commencer ici',
    'nav.math': 'Mathématique',
    'nav.math_atlas': 'Atlas mathématique',
    'nav.math_explorations': 'Explorations mathématiques',
    'nav.pt_toolbox': 'Atlas des outils PT',
    'nav.rh_programme': 'Programme RH',
    'nav.physics': 'Physique',
    'nav.chemistry': 'Chimie',
    'nav.theorems': 'Théorèmes',
    'nav.first_principle': 'Principe fondamental',
    'nav.prime_gaps': 'Gaps premiers',
    'nav.observables': '43 observables',
    'nav.fundamental_questions': 'Questions fondamentales',
    'nav.time': 'Temps',
    'nav.relativity': 'Relativité',
    'nav.quantum_gravity': 'Gravité quantique',
    'nav.nuclear': 'Nucléaire',
    'nav.cosmology': 'Cosmologie',
    'nav.big_bang': 'Cosmogonie',
    'nav.periodic_table': 'Tableau périodique',
    'nav.interactive_proof': 'Preuve interactive',
    'nav.ionization_energies': 'Énergies d’ionisation',
    'nav.electron_affinities': 'Affinités électroniques',
    'nav.benchmark_ptc': 'Benchmark PTC',
    'nav.results': 'Résultats',
    'nav.result_status': 'Statuts des résultats',
    'nav.hard_to_dismiss': 'Difficile à rejeter',
    'nav.limits': 'Limites',
    'nav.articles': 'Articles',
    'nav.essays': 'En quelques mots',
    'nav.calculator': 'Calculatrices',
    'nav.chain': 'Chaîne causale',
    'nav.scripts': 'Scripts',
    'nav.lean': 'Formalisation Lean',
    'nav.audit': 'Audit',
    'nav.faq': 'FAQ',
    'nav.glossary': 'Glossaire',
    'nav.about': 'À propos',
    'nav.contact': 'Contact',
    'nav.principles': '7 principes',
    'nav.confirmations': 'Confirmations externes',
    'nav.critique': 'Critique & débat',
    'nav.research': 'Recherches',
    // Site meta
    'site.title': 'Théorie de la Persistance',
    'site.tagline': 'Du crible au Modèle Standard.',
    'site.description':
      "Théorie de la persistance : sous contrainte, une part se disperse en entropie, une autre persiste comme structure. La forme rend la contrainte visible.",
    // Buttons / actions
    'action.read_monograph': 'Lire la monographie',
    'action.download_pdf': 'Télécharger le PDF',
    'action.contact': 'Contact',
    'action.switch_lang': 'English',
    // Footer
    'footer.code_license': 'Code sous licence MIT',
    'footer.content_license': 'Contenu sous licence CC BY 4.0',
    'footer.author': 'yan senez',
    'footer.last_update': 'Dernière mise à jour',
    'footer.preprint_notice':
      'Préimpression — non révisée par les pairs.',
    'footer.dev_notice':
      'Site en cours de développement et moteur PTC en évolution active : des erreurs, incohérences, pages incomplètes ou chiffres datés peuvent subsister.',
    'footer.report_issue': 'Signaler une erreur',
    'depth.label': 'Profondeur',
    'depth.L1': 'Simple',
    'depth.L2': 'Standard',
    'depth.L3': 'Technique',
    'nav.more': 'Plus',
    'nav.menu': 'Menu',
    'nav.critical_path': 'Carte du chemin critique',
  },
  en: {
    'nav.monograph': 'Monograph',
    'nav.start': 'Start here',
    'nav.math': 'Mathematics',
    'nav.math_atlas': 'Mathematical atlas',
    'nav.math_explorations': 'Mathematical explorations',
    'nav.pt_toolbox': 'PT toolbox atlas',
    'nav.rh_programme': 'RH programme',
    'nav.physics': 'Physics',
    'nav.chemistry': 'Chemistry',
    'nav.theorems': 'Theorems',
    'nav.first_principle': 'Fundamental principle',
    'nav.prime_gaps': 'Prime gaps',
    'nav.observables': '43 observables',
    'nav.fundamental_questions': 'Fundamental questions',
    'nav.time': 'Time',
    'nav.relativity': 'Relativity',
    'nav.quantum_gravity': 'Quantum gravity',
    'nav.nuclear': 'Nuclear',
    'nav.cosmology': 'Cosmology',
    'nav.big_bang': 'Cosmogony',
    'nav.periodic_table': 'Periodic table',
    'nav.interactive_proof': 'Interactive proof',
    'nav.ionization_energies': 'Ionization energies',
    'nav.electron_affinities': 'Electron affinities',
    'nav.benchmark_ptc': 'PTC benchmark',
    'nav.results': 'Results',
    'nav.result_status': 'Result status',
    'nav.hard_to_dismiss': 'Hard to dismiss',
    'nav.limits': 'Limits',
    'nav.articles': 'Articles',
    'nav.essays': 'In a few words',
    'nav.calculator': 'Calculators',
    'nav.chain': 'Causal chain',
    'nav.scripts': 'Scripts',
    'nav.lean': 'Lean formalisation',
    'nav.audit': 'Audit',
    'nav.faq': 'FAQ',
    'nav.glossary': 'Glossary',
    'nav.about': 'About',
    'nav.contact': 'Contact',
    'nav.principles': '7 principles',
    'nav.confirmations': 'External confirmations',
    'nav.critique': 'Critique & debate',
    'nav.research': 'Research',
    'site.title': 'The Theory of Persistence',
    'site.tagline': 'From the sieve to the Standard Model.',
    'site.description':
      'Persistence Theory: under constraint, one part disperses as entropy, another persists as structure. Form is the structure made visible by constraint.',
    'action.read_monograph': 'Read the monograph',
    'action.download_pdf': 'Download PDF',
    'action.contact': 'Contact',
    'action.switch_lang': 'Français',
    'footer.code_license': 'Code under MIT licence',
    'footer.content_license': 'Content under CC BY 4.0',
    'footer.author': 'yan senez',
    'footer.last_update': 'Last update',
    'footer.preprint_notice':
      'Preprint — not peer-reviewed.',
    'footer.dev_notice':
      'Website under active development and PTC engine evolving: errors, inconsistencies, incomplete pages, or stale figures may remain.',
    'footer.report_issue': 'Report an issue',
    'depth.label': 'Depth',
    'depth.L1': 'Plain',
    'depth.L2': 'Standard',
    'depth.L3': 'Technical',
    'nav.more': 'More',
    'nav.menu': 'Menu',
    'nav.critical_path': 'Critical-path map',
  },
} as const;

export type UIKey = keyof (typeof ui)['fr'];

export function useTranslations(lang: Lang) {
  return function t(key: UIKey): string {
    return ui[lang][key] || ui[defaultLang][key];
  };
}

/**
 * Detect language from URL pathname.
 * `/en/foo` -> 'en'; `/foo` -> 'fr' (default).
 */
export function getLangFromUrl(url: URL): Lang {
  const [, segment] = url.pathname.split('/');
  if (segment === 'en') return 'en';
  return 'fr';
}

/**
 * Build the path to the same page in the other language.
 * `/foo` (FR) <-> `/en/foo` (EN); `/` <-> `/en/`.
 */
export function alternateLangPath(url: URL, target: Lang): string {
  const segments = url.pathname.split('/').filter(Boolean);
  const isCurrentEn = segments[0] === 'en';
  const rest = isCurrentEn ? segments.slice(1) : segments;
  if (target === 'en') return `/en/${rest.join('/')}` || '/en/';
  return `/${rest.join('/')}` || '/';
}
