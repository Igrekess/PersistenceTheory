import { defineCollection, z } from 'astro:content';
import { glob } from 'astro/loaders';

const epistemicTag = z.enum([
  'THM',
  'IDENTITY',
  'DER',
  'COND',
  'BRIDGE',
  'VAL',
  'PRED',
  'CLOSURE',
]);

const lang = z.enum(['fr', 'en']);

/* ---------- Theorems ---------------------------------------------------- */

const theoremSchema = z.object({
  id: z.string(), // T0, T1, ..., GFT, L0, N1, ...
  name: z.string(),
  short: z.string(),
  statement: z.string(), // 1-line statement (markdown)
  status: epistemicTag,
  monograph_chapter: z.string(), // ch label, e.g. 'chap:fixed_point'
  order: z.number().int(),
  tier: z.enum(['fundamental', 'secondary']).default('fundamental'),
  lang,
});

const theorems_fr = defineCollection({
  loader: glob({ pattern: '*.mdx', base: './src/content/theorems/fr' }),
  schema: theoremSchema,
});

const theorems_en = defineCollection({
  loader: glob({ pattern: '*.mdx', base: './src/content/theorems/en' }),
  schema: theoremSchema,
});

/* ---------- Articles --------------------------------------------------- */

const articles = defineCollection({
  loader: glob({
    pattern: '**/*.mdx',
    base: './src/content/articles',
    generateId: ({ data }) => `${data.lang}/${data.slug}`,
  }),
  schema: z.object({
    slug: z.string(),
    title: z.string(),
    abstract: z.string(),
    pdf_url: z.string().url().optional(),
    repo_url: z.string().url().optional(),
    domain: z.enum(['math', 'physics', 'chemistry', 'biology', 'cosmology', 'geophysics']),
    status: epistemicTag,
    pages: z.number().int().optional(),
    year: z.number().int().default(2026),
    lang,
    order: z.number().int().default(99),
  }),
});

/* ---------- Essays (V0.5 vulgarisation) -------------------------------- */

const essaySchema = z.object({
  slug: z.string(),
  title: z.string(),
  summary: z.string(),
  level: z.enum(['L1', 'L2']).default('L1'),
  reading_minutes: z.number().int().default(5),
  related_theorems: z.array(z.string()).default([]),
  lang,
  order: z.number().int().default(99),
});

const essays_fr = defineCollection({
  loader: glob({ pattern: '*.mdx', base: './src/content/essays/fr' }),
  schema: essaySchema,
});

const essays_en = defineCollection({
  loader: glob({ pattern: '*.mdx', base: './src/content/essays/en' }),
  schema: essaySchema,
});

export const collections = { theorems_fr, theorems_en, articles, essays_fr, essays_en };
