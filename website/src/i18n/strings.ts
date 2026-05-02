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
    'nav.results': 'Résultats',
    'nav.result_status': 'Statuts des résultats',
    'nav.hard_to_dismiss': 'Difficile à rejeter',
    'nav.limits': 'Limites',
    'nav.articles': 'Articles',
    'nav.essays': 'En quelques mots',
    'nav.calculator': 'Calculatrices',
    'nav.chain': 'Chaîne',
    'nav.scripts': 'Scripts',
    'nav.audit': 'Audit',
    'nav.faq': 'FAQ',
    'nav.glossary': 'Glossaire',
    'nav.about': 'À propos',
    'nav.contact': 'Contact',
    // Site meta
    'site.title': 'Théorie de la Persistance',
    'site.tagline': 'Du crible au Modèle Standard.',
    'site.description':
      "La Théorie de la Persistance dérive 43 observables du Modèle Standard depuis la symétrie dérivée s = 1/2, sans paramètre continu ajusté. Monographie, preuves complètes, scripts de vérification.",
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
      'Site en cours de développement : des erreurs, incohérences ou pages incomplètes peuvent subsister.',
    'footer.report_issue': 'Signaler une erreur',
    'depth.label': 'Profondeur',
    'depth.L1': 'Vulgarisé',
    'depth.L2': 'Standard',
    'depth.L3': 'Technique',
    'nav.more': 'Plus',
    'nav.menu': 'Menu',
  },
  en: {
    'nav.monograph': 'Monograph',
    'nav.start': 'Start here',
    'nav.math': 'Mathematics',
    'nav.math_atlas': 'Mathematical atlas',
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
    'nav.results': 'Results',
    'nav.result_status': 'Result status',
    'nav.hard_to_dismiss': 'Hard to dismiss',
    'nav.limits': 'Limits',
    'nav.articles': 'Articles',
    'nav.essays': 'In a few words',
    'nav.calculator': 'Calculators',
    'nav.chain': 'Chain',
    'nav.scripts': 'Scripts',
    'nav.audit': 'Audit',
    'nav.faq': 'FAQ',
    'nav.glossary': 'Glossary',
    'nav.about': 'About',
    'nav.contact': 'Contact',
    'site.title': 'The Theory of Persistence',
    'site.tagline': 'From the sieve to the Standard Model.',
    'site.description':
      'The Theory of Persistence derives 43 Standard Model observables from the derived symmetry s = 1/2, with no fitted continuous parameter. Monograph, full proofs, verification scripts.',
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
      'This website is under active development: errors, inconsistencies, or incomplete pages may remain.',
    'footer.report_issue': 'Report an issue',
    'depth.label': 'Depth',
    'depth.L1': 'Plain',
    'depth.L2': 'Standard',
    'depth.L3': 'Technical',
    'nav.more': 'More',
    'nav.menu': 'Menu',
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
