/** @type {import('tailwindcss').Config} */
export default {
  content: ['./src/**/*.{astro,html,js,jsx,md,mdx,svelte,ts,tsx,vue}'],
  theme: {
    extend: {
      // Palette épistémique modernisée : plus froide, plus nette à l'écran.
      colors: {
        thmblue:    '#14566b', // Theorem
        idgreen:    '#25785a', // Identity
        condorange: '#b7791f', // Conditional
        derorange:  '#c96f24', // Derived
        bridgepurple: '#6856a6', // Bridge
        valgray:    '#64748b', // Validation
        predred:    '#b8323a', // Prediction
        closuregold: '#a8871b', // Closure
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
