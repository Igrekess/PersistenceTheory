/**
 * GitHub repository coordinates.
 * Single source of truth for repo URL and raw asset URLs.
 * Update OWNER / NAME here if the repo moves.
 */

export const REPO_OWNER = 'Igrekess';
export const REPO_NAME = 'PersistenceTheory';
export const REPO_BRANCH = 'main';

export const REPO_URL = `https://github.com/${REPO_OWNER}/${REPO_NAME}`;
export const REPO_RAW = `https://raw.githubusercontent.com/${REPO_OWNER}/${REPO_NAME}/${REPO_BRANCH}`;
export const REPO_ISSUES = `${REPO_URL}/issues`;

/** Build a raw URL for a file at `path` (relative to repo root). */
export function repoRaw(path: string): string {
  const clean = path.startsWith('/') ? path.slice(1) : path;
  return `${REPO_RAW}/${clean}`;
}

/** Convenience: monograph PDFs (at repo root). */
export const MONOGRAPH_PDF_EN = repoRaw('TheTheoryOfPersistence.pdf');
export const MONOGRAPH_PDF_FR = repoRaw('TheorieDeLaPersistance_FR.pdf');
