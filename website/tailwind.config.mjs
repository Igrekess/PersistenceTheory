/** @type {import('tailwindcss').Config} */
export default {
  content: ['./src/**/*.{astro,html,js,jsx,md,mdx,svelte,ts,tsx,vue}'],
  theme: {
    extend: {
      // Palette épistémique modernisée : plus froide, plus nette à l'écran.
      colors: {
        // Epistemic palette (status tags)
        thmblue:    '#14566b', // Theorem
        idgreen:    '#25785a', // Identity
        condorange: '#b7791f', // Conditional
        derorange:  '#c96f24', // Derived
        bridgepurple: '#6856a6', // Bridge
        valgray:    '#64748b', // Validation
        predred:    '#b8323a', // Prediction
        closuregold: '#a8871b', // Closure
        // Domain palette (eyebrow / wayfinding) — used on hero eyebrow only
        'domain-math':      '#4338ca', // indigo
        'domain-physics':   '#0d9488', // teal
        'domain-chemistry': '#b45309', // ambre
        'domain-bio':       '#9333ea', // rose-violet
        'domain-cosmo':     '#1e3a8a', // bleu-nuit
        'domain-geo':       '#78350f', // brun terre
        'domain-meta':      '#64748b', // slate-500 (gris-bleu clair)
      },
      fontFamily: {
        serif: ['ui-serif', 'Georgia', 'Cambria', 'Times New Roman', 'serif'],
        sans: ['InterVariable', 'Inter', 'ui-sans-serif', 'system-ui', '-apple-system', 'BlinkMacSystemFont', 'Segoe UI', 'sans-serif'],
        mono: ['SFMono-Regular', 'ui-monospace', 'Cascadia Code', 'Menlo', 'monospace'],
      },
      typography: ({ theme }) => ({
        DEFAULT: {
          css: {
            maxWidth: '70ch',
          },
        },
      }),
    },
  },
  plugins: [],
};
