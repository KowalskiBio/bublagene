// Properties of a translated protein: mass, charge, pI, hydropathy, and a plain
// motif scan. These are textbook estimates (Kyte–Doolittle hydropathy, the
// EMBOSS pKa set for pI), not predictions from structure.

import { proteinMass } from './seq.js';
import { BUILTIN_LIBRARY } from './feature-library.js';

export const KYTE_DOOLITTLE = {
  A: 1.8, R: -4.5, N: -3.5, D: -3.5, C: 2.5, Q: -3.5, E: -3.5, G: -0.4, H: -3.2, I: 4.5,
  L: 3.8, K: -3.9, M: 1.9, F: 2.8, P: -1.6, S: -0.8, T: -0.7, W: -0.9, Y: -1.3, V: 4.2,
};
const PKA = { nterm: 8.6, cterm: 3.6, K: 10.8, R: 12.5, H: 6.5, D: 3.9, E: 4.1, C: 8.5, Y: 10.1 };

// translation ends at the first stop, so nothing after it is part of the protein
const clean = (p) => String(p).toUpperCase().replace(/\*[\s\S]*$/, '').replace(/[^A-Z]/g, '');

export function chargeAt(protein, pH) {
  const p = clean(protein);
  if (!p) return 0;
  const pos = (pka) => 1 / (1 + 10 ** (pH - pka));
  const neg = (pka) => -1 / (1 + 10 ** (pka - pH));
  let q = pos(PKA.nterm) + neg(PKA.cterm);
  for (const ch of p) {
    if (ch === 'K' || ch === 'R' || ch === 'H') q += pos(PKA[ch]);
    else if (ch === 'D' || ch === 'E' || ch === 'C' || ch === 'Y') q += neg(PKA[ch]);
  }
  return q;
}

export function isoelectricPoint(protein) {
  let lo = 0;
  let hi = 14;
  for (let i = 0; i < 60; i++) {
    const mid = (lo + hi) / 2;
    if (chargeAt(protein, mid) > 0) lo = mid; else hi = mid;
  }
  return (lo + hi) / 2;
}

export function composition(protein) {
  const counts = {};
  for (const ch of clean(protein)) counts[ch] = (counts[ch] ?? 0) + 1;
  return counts;
}

export function proteinProperties(protein) {
  const p = clean(protein);
  const c = composition(p);
  const gravy = p.length ? [...p].reduce((s, ch) => s + (KYTE_DOOLITTLE[ch] ?? 0), 0) / p.length : 0;
  return {
    length: p.length,
    mass: p.length ? proteinMass(p) : 0,
    pI: p.length ? isoelectricPoint(p) : 0,
    chargeAt7: chargeAt(p, 7.4),
    gravy,
    // Pace et al.: assumes all Cys pairs are cystines (reduced: without the C term)
    extinction: 5500 * (c.W ?? 0) + 1490 * (c.Y ?? 0),
    extinctionOxidised: 5500 * (c.W ?? 0) + 1490 * (c.Y ?? 0) + 125 * Math.floor((c.C ?? 0) / 2),
    composition: c,
    acidic: (c.D ?? 0) + (c.E ?? 0),
    basic: (c.K ?? 0) + (c.R ?? 0),
  };
}

// Mean hydropathy in a sliding window, one value per window centre.
export function hydropathyProfile(protein, window = 9) {
  const p = clean(protein);
  const half = Math.floor(window / 2);
  const out = [];
  for (let i = 0; i + window <= p.length; i++) {
    let sum = 0;
    for (let k = 0; k < window; k++) sum += KYTE_DOOLITTLE[p[i + k]] ?? 0;
    out.push({ position: i + half, value: sum / window });
  }
  return out;
}

// Stretches of ≥19 residues whose mean hydropathy stays above 1.6: the classic
// rule of thumb for a membrane-spanning helix. A hint, not a prediction.
export function hydrophobicSegments(protein, { window = 19, threshold = 1.6 } = {}) {
  const prof = hydropathyProfile(protein, window);
  const half = Math.floor(window / 2);
  const segs = [];
  let run = null;
  for (const { position, value } of prof) {
    if (value >= threshold) {
      if (!run) run = { start: position - half, end: position + half + 1 };
      else run.end = position + half + 1;
    } else if (run) { segs.push(run); run = null; }
  }
  if (run) segs.push(run);
  // a window that averages high can still reach over a charged residue at either
  // end; the helix itself starts and stops at hydrophobic residues
  const p = clean(protein);
  return segs
    .map(({ start, end }) => {
      const kd = (i) => KYTE_DOOLITTLE[p[i]] ?? 0;
      const mean3 = (i) => (kd(i) + kd(i + 1) + kd(i + 2)) / 3;
      let a = start;
      let b = end;
      while (b - a > 3 && mean3(a) < 0.5) a++;
      while (b - a > 3 && mean3(b - 3) < 0.5) b--;
      while (a < b && kd(a) < 0) a++;
      while (b > a && kd(b - 1) < 0) b--;
      return { start: a, end: b };
    })
    .filter((x) => x.end - x.start >= 12);
}

// Motifs worth pointing out in a protein: the known tags and cleavage sites,
// N-glycosylation sequons (N-X-S/T, X ≠ P), and possible membrane helices.
export function proteinMotifs(protein) {
  const p = clean(protein);
  const out = [];
  for (const e of BUILTIN_LIBRARY.filter((x) => x.kind === 'protein')) {
    for (let at = p.indexOf(e.seq); at >= 0; at = p.indexOf(e.seq, at + 1)) out.push({ name: e.name, start: at, end: at + e.seq.length, kind: 'tag or site' });
  }
  for (let i = 0; i + 3 <= p.length; i++) {
    if (p[i] === 'N' && p[i + 1] !== 'P' && (p[i + 2] === 'S' || p[i + 2] === 'T')) out.push({ name: 'N-glycosylation sequon', start: i, end: i + 3, kind: 'motif' });
  }
  for (const s of hydrophobicSegments(p)) out.push({ name: 'Possible transmembrane helix', start: s.start, end: s.end, kind: 'hydrophobic' });
  // a short site inside a longer tag (enterokinase inside FLAG) is part of it
  const kept = out.filter((m) => !out.some((o) => o !== m && o.kind === m.kind && o.start <= m.start && o.end >= m.end && o.end - o.start > m.end - m.start));
  return kept.sort((a, b) => a.start - b.start);
}
