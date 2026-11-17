// Bounded RNA secondary-structure prediction.
//
// Model (explicitly scoped): canonical Watson–Crick pairs plus G·U, a Turner
// subset of stacking and hairpin-loop energies, simple interior/bulge loops and
// a multiloop approximation. Sequences longer than `maxLength` are refused
// rather than silently truncated. This is not a full Turner 2004 implementation
// and makes no accuracy claim beyond the simple structures listed in the tests.

import { hairpin as hairpinDG, dimer as duplexDG } from './thermo.js';

export const FOLD_MODEL = 'Turner subset (WC + GU, stacking, hairpin, simple loops)';

const PAIR = {
  GC: -3.1, CG: -2.4, AU: -1.7, UA: -1.4, GU: -0.8, UG: -0.5,
};

const STACK = {
  GCGC: -3.3, GCCG: -2.4, GCUA: -2.1, GCAU: -2.2,
  CGGC: -2.4, CGCG: -3.3, CGAU: -2.1, CGUA: -2.0,
  AUUA: -1.1, AUAU: -1.2, AUGC: -2.1, AUCG: -2.1,
  UAAU: -1.2, UAUA: -1.1, UAGC: -2.2, UACG: -2.0,
  GUGU: -0.9, GUUG: -0.5, GUAU: -1.4, GUUA: -1.3,
  UGGU: -0.5, UGUG: -0.9, UGAU: -1.4, UGUA: -1.3,
};

function pairEnergy(a, b) {
  return PAIR[`${a}${b}`] ?? Infinity;
}

function hairpinEnergy(size) {
  if (size < 3) return Infinity;
  const table = { 3: 5.4, 4: 5.6, 5: 5.7, 6: 5.4, 7: 6.0, 8: 5.5, 9: 6.4 };
  if (table[size] != null) return table[size];
  return 6.0 + 1.75 * Math.log(size / 9);
}

function loopEnergy(u1, u2) {
  const size = u1 + u2;
  if (size === 0) return 0;
  if (size === 1) return 1.7;
  if (size === 2) return 1.5;
  return 1.3 + 0.6 * Math.log(size) + 0.2 * Math.abs(u1 - u2);
}

export function fold(seq, { maxLength = 250, minHairpin = 3 } = {}) {
  const s = String(seq).toUpperCase().replace(/U/g, 'T').replace(/[^ACGT]/g, '');
  if (s.length > maxLength) throw new Error(`Folding is limited to ${maxLength} nt (got ${s.length})`);
  const n = s.length;
  const INF = Infinity;
  const V = Array.from({ length: n }, () => new Float64Array(n).fill(INF));
  const W = Array.from({ length: n }, () => new Float64Array(n).fill(0));
  const vChoice = Array.from({ length: n }, () => new Int16Array(n).fill(-1));
  const wChoice = Array.from({ length: n }, () => new Int16Array(n).fill(-1));

  const canPair = (i, j) => j - i - 1 >= minHairpin && Number.isFinite(pairEnergy(s[i], s[j]));

  for (let len = 1; len <= n; len++) {
    for (let i = 0; i + len - 1 < n; i++) {
      const j = i + len - 1;
      // V (paired i,j) must be known before W can use V[i][j].
      const pairable = canPair(i, j);
      if (pairable) {
      // V: hairpin, stack, or a simple interior/bulge loop
      let vBest = hairpinEnergy(j - i - 1);
      let vc = -3;
      if (i + 1 < j - 1 && Number.isFinite(V[i + 1][j - 1])) {
        const key = `${s[i]}${s[j]}${s[i + 1]}${s[j - 1]}`;
        const stack = STACK[key] ?? (pairEnergy(s[i], s[j]) + pairEnergy(s[i + 1], s[j - 1]));
        const cand = stack + V[i + 1][j - 1];
        if (cand < vBest) { vBest = cand; vc = -1; }
      }
      for (let k = i + 1; k < j; k++) {
        for (let l = k + 1; l < j; l++) {
          if (!canPair(k, l)) continue;
          const u1 = k - i - 1;
          const u2 = j - l - 1;
          if (u1 + u2 > 30) continue;
          const cand = loopEnergy(u1, u2) + V[k][l];
          if (cand < vBest) { vBest = cand; vc = k * 1000 + l; }
        }
      }
      V[i][j] = vBest;
      vChoice[i][j] = vc;
      }
      // W: leave i unpaired, or pair it with some k
      let best = i + 1 <= j ? W[i + 1][j] : 0;
      let choice = -1;
      for (let k = i + minHairpin + 1; k <= j; k++) {
        if (!canPair(i, k)) continue;
        const after = k + 1 <= j ? W[k + 1][j] : 0;
        const cand = V[i][k] + after;
        if (cand < best) { best = cand; choice = k; }
      }
      W[i][j] = best;
      wChoice[i][j] = choice;
    }
  }

  const pairs = [];
  const traceW = (i, j) => {
    if (i > j) return;
    const c = wChoice[i][j];
    if (c === -1) traceW(i + 1, j);
    else if (c >= 0) { traceV(i, c); traceW(c + 1, j); }
  };
  const traceV = (i, j) => {
    if (i >= j) return;
    pairs.push([i, j]);
    const c = vChoice[i][j];
    if (c === -1) traceV(i + 1, j - 1);
    else if (c >= 1000) {
      const k = Math.floor(c / 1000);
      const l = c % 1000;
      traceV(k, l);
    }
  };
  traceW(0, n - 1);

  const structure = pairsToDotBracket(n, pairs);
  return {
    sequence: s,
    structure,
    energy: Math.round(W[0]?.[n - 1] * 100) / 100,
    pairs,
    length: n,
    model: FOLD_MODEL,
    maxLength,
    limits: 'Hairpins, stacks and simple interior/bulge loops; no pseudoknots or full multiloop parameters.',
  };
}

export function pairsToDotBracket(n, pairs) {
  const chars = Array(n).fill('.');
  for (const [i, j] of pairs) { chars[i] = '('; chars[j] = ')'; }
  return chars.join('');
}

export function dotBracketToPairs(structure) {
  const stack = [];
  const pairs = [];
  for (let i = 0; i < structure.length; i++) {
    if (structure[i] === '(') stack.push(i);
    else if (structure[i] === ')') {
      const j = stack.pop();
      if (j == null) throw new Error('Unbalanced dot-bracket structure');
      pairs.push([j, i]);
    }
  }
  if (stack.length) throw new Error('Unbalanced dot-bracket structure');
  return pairs.sort((a, b) => a[0] - b[0]);
}

// Rows of sequence + structure for a text view.
export function structureView(seq, structure, { width = 60 } = {}) {
  const rows = [];
  for (let i = 0; i < seq.length; i += width) {
    rows.push({ start: i, sequence: seq.slice(i, i + width), structure: structure.slice(i, i + width) });
  }
  return rows;
}

export function structureToText(seq, structure, { width = 60 } = {}) {
  return structureView(seq, structure, { width })
    .map((row) => `${String(row.start + 1).padStart(6)} ${row.sequence}\n       ${row.structure}`)
    .join('\n');
}

// SVG arcs over the sequence, suitable for a panel or export.
export function structureToSvg(seq, structure, { width = 900, rowHeight = 120, perRow = 60 } = {}) {
  const n = seq.length;
  const cols = Math.min(perRow, Math.max(1, n));
  const rows = Math.ceil(n / cols);
  const charWidth = Math.max(6, (width - 40) / cols);
  const height = rows * rowHeight + 30;
  const parts = [`<svg xmlns="http://www.w3.org/2000/svg" width="${width}" height="${height}" viewBox="0 0 ${width} ${height}">`];
  parts.push(`<rect width="100%" height="100%" fill="white"/>`);
  for (const [a, b] of dotBracketToPairs(structure)) {
    const ra = Math.floor(a / cols);
    const ca = a % cols;
    const rb = Math.floor(b / cols);
    const cb = b % cols;
    if (ra !== rb) continue;
    const x1 = 20 + ca * charWidth;
    const x2 = 20 + cb * charWidth;
    const y = ra * rowHeight + 40;
    const r = Math.max(6, (x2 - x1) / 2);
    parts.push(`<path d="M ${x1.toFixed(1)} ${y} A ${r.toFixed(1)} ${r.toFixed(1)} 0 0 1 ${x2.toFixed(1)} ${y}" fill="none" stroke="#5b8fd6" stroke-width="1.2"/>`);
  }
  for (let i = 0; i < n; i++) {
    const r = Math.floor(i / cols);
    const c = i % cols;
    const x = 20 + c * charWidth;
    const y = r * rowHeight + 90;
    parts.push(`<text x="${x.toFixed(1)}" y="${y}" font-family="monospace" font-size="11" fill="#333">${seq[i]}</text>`);
  }
  parts.push('</svg>');
  return parts.join('\n');
}

export { hairpinDG, duplexDG };
