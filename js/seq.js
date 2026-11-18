// Sequence utilities: complements, translation, IUPAC matching, GC, melting temperature.
// Coordinates everywhere in BublaGene: 0-based, end-exclusive; on circular
// sequences a range may wrap (end <= start means it crosses the origin).

import { translateWithTable, DEFAULT_TABLE } from './codon-tables.js';
import { memoLatest } from './cache.js';

const COMP = {
  A: 'T', T: 'A', G: 'C', C: 'G', U: 'A', N: 'N',
  R: 'Y', Y: 'R', S: 'S', W: 'W', K: 'M', M: 'K',
  B: 'V', V: 'B', D: 'H', H: 'D',
};

export function complement(s) {
  let out = '';
  for (const ch of s) {
    const up = ch.toUpperCase();
    const c = COMP[up] ?? up;
    out += ch === up ? c : c.toLowerCase();
  }
  return out;
}

export function revcomp(s) {
  return complement(s).split('').reverse().join('');
}

// Read a sequence pasted from anywhere: FASTA header lines, numbers, spaces and
// line breaks are dropped, U becomes T, and letters that are not IUPAC bases are
// left out and counted, so the caller can say so.
export function readPastedSequence(text) {
  const body = String(text).split(/\r?\n/).filter((line) => !line.trimStart().startsWith('>')).join('');
  let seq = '';
  let ignored = 0;
  for (const ch of body.toUpperCase()) {
    if (!/[A-Z]/.test(ch)) continue; // digits, spaces, dashes and punctuation are layout
    const c = ch === 'U' ? 'T' : ch;
    if ('ACGTNRYSWKMBDHV'.includes(c)) seq += c; else ignored++;
  }
  return { seq, ignored };
}

export function cleanSequence(s) {
  return s.replace(/[^A-Za-z]/g, '').toUpperCase();
}

// ---- ranges -------------------------------------------------------------

export function mod(a, n) {
  return ((a % n) + n) % n;
}

export function rangeLength(r, L) {
  if (r.end > r.start) return r.end - r.start;
  if (r.end === r.start) return r.full ? L : 0;
  return L - r.start + r.end;
}

export function subseq(seq, start, end, circular) {
  const L = seq.length;
  if (!circular || end > start) return seq.slice(start, end);
  return seq.slice(start) + seq.slice(0, end);
}

// Does position p (0-based) fall inside range r on a sequence of length L?
export function inRange(p, r, L) {
  if (r.end > r.start) return p >= r.start && p < r.end;
  if (r.end === r.start) return !!r.full;
  return p >= r.start || p < r.end;
}

// Is range a completely inside range b (both possibly wrapping)?
export function rangeInside(a, b, L) {
  const lb = rangeLength(b, L);
  const offA = mod(a.start - b.start, L);
  return offA + rangeLength(a, L) <= lb;
}

// ---- composition ----------------------------------------------------------

export function gcContent(s) {
  if (!s.length) return 0;
  let gc = 0;
  for (const ch of s.toUpperCase()) if (ch === 'G' || ch === 'C' || ch === 'S') gc++;
  return (100 * gc) / s.length;
}

// Counts of each IUPAC symbol, plus a plain A/C/G/T tally (U counted as T).
export function nucleotideCount(s) {
  const counts = {};
  let gc = 0;
  let at = 0;
  for (const ch of s.toUpperCase()) {
    counts[ch] = (counts[ch] ?? 0) + 1;
    if (ch === 'G' || ch === 'C' || ch === 'S') gc++;
    if (ch === 'A' || ch === 'T' || ch === 'U' || ch === 'W') at++;
  }
  return { ...counts, length: s.length, gc, at, gcPercent: s.length ? (100 * gc) / s.length : 0 };
}

// Average molecular weight of double-stranded DNA (sodium salt), Da.
export function dnaMolecularWeight(s) {
  return Math.round(s.length * 617.96 + 36.04);
}

// ---- translation (default genetic code: NCBI table 1) -------------------

const BASES = 'TCAG';
const AA = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG';
const CODON = {};
for (let i = 0; i < 64; i++) {
  CODON[BASES[i >> 4] + BASES[(i >> 2) & 3] + BASES[i & 3]] = AA[i];
}

export const THREE_LETTER = {
  A: 'Ala', R: 'Arg', N: 'Asn', D: 'Asp', C: 'Cys', Q: 'Gln', E: 'Glu', G: 'Gly',
  H: 'His', I: 'Ile', L: 'Leu', K: 'Lys', M: 'Met', F: 'Phe', P: 'Pro', S: 'Ser',
  T: 'Thr', W: 'Trp', Y: 'Tyr', V: 'Val', '*': 'Stop', X: 'Xaa',
};

// translate(dna [, tableId | { table, initiator, stopSymbol }])
// The default table keeps the historical one-argument behaviour.
/** @param {string} dna @param {number|{table?:number,transl_table?:number,initiator?:boolean,start?:boolean,stopSymbol?:string}} [options] */
export function translate(dna, options = DEFAULT_TABLE) {
  if (typeof options === 'number') return translateWithTable(dna, options);
  const opts = options ?? {};
  const table = opts.table ?? opts.transl_table ?? DEFAULT_TABLE;
  if (opts.initiator || opts.start) return translateWithTable(dna, table, { initiator: true, stopSymbol: opts.stopSymbol });
  if (opts.stopSymbol && opts.stopSymbol !== '*') return translateWithTable(dna, table, { stopSymbol: opts.stopSymbol });
  if (Number(table) !== DEFAULT_TABLE) return translateWithTable(dna, table);
  const s = dna.toUpperCase().replace(/U/g, 'T');
  let p = '';
  for (let i = 0; i + 3 <= s.length; i += 3) p += CODON[s.slice(i, i + 3)] ?? 'X';
  return p;
}

// Six reading frames, in the order +1,+2,+3,-1,-2,-3.
export function sixFrameTranslation(dna, tableId = DEFAULT_TABLE) {
  const s = dna.toUpperCase().replace(/U/g, 'T');
  const rc = revcomp(s);
  const out = [];
  for (const strand of [1, -1]) {
    const source = strand === 1 ? s : rc;
    for (let frame = 0; frame < 3; frame++) {
      out.push({ strand, frame: frame + 1, protein: translate(source.slice(frame), tableId) });
    }
  }
  return out;
}

const AA_MASS = {
  A: 71.0788, R: 156.1875, N: 114.1038, D: 115.0886, C: 103.1388, E: 129.1155,
  Q: 128.1307, G: 57.0519, H: 137.1411, I: 113.1594, L: 113.1594, K: 128.1741,
  M: 131.1926, F: 147.1766, P: 97.1167, S: 87.0782, T: 101.1051, W: 186.2132,
  Y: 163.1760, V: 99.1326,
};

export function proteinMass(p) {
  let m = 18.01524;
  for (const a of p) if (AA_MASS[a]) m += AA_MASS[a];
  return m;
}

// ---- IUPAC ------------------------------------------------------------------

export const IUPAC = {
  A: 'A', C: 'C', G: 'G', T: 'T', U: 'T', R: 'AG', Y: 'CT', S: 'CG', W: 'AT',
  K: 'GT', M: 'AC', B: 'CGT', D: 'AGT', H: 'ACT', V: 'ACG', N: 'ACGT',
};

export function iupacRegexSource(pattern) {
  let src = '';
  for (const ch of pattern.toUpperCase()) {
    const set = IUPAC[ch];
    if (!set) throw new Error(`Not a nucleotide code: ${ch}`);
    src += set.length === 1 ? set : `[${set}]`;
  }
  return src;
}

// All occurrences (overlapping) of an IUPAC pattern on both strands.
// Returns [{start, end, strand}] with start/end on the top strand.
export function findPattern(seq, pattern, circular) {
  const n = pattern.length;
  const L = seq.length;
  if (!n || n > L) return [];
  const hay = circular ? seq + seq.slice(0, n - 1) : seq;
  const hits = [];
  const scan = (pat, strand) => {
    const re = new RegExp(`(?=(${iupacRegexSource(pat)}))`, 'g');
    let m;
    while ((m = re.exec(hay))) {
      if (m.index >= L) break;
      hits.push({ start: m.index, end: (m.index + n) % (circular ? L : Infinity), strand });
      re.lastIndex = m.index + 1;
    }
  };
  scan(pattern, 1);
  const rc = revcomp(pattern);
  if (rc !== pattern.toUpperCase()) scan(rc, -1);
  for (const h of hits) if (!circular) h.end = h.start + n;
  return hits.sort((a, b) => a.start - b.start);
}

// ---- melting temperature ------------------------------------------------------
// Nearest-neighbour model, SantaLucia (1998) unified parameters,
// 50 mM Na+, 250 nM oligo by default.

const NN = {
  AA: [-7.9, -22.2], TT: [-7.9, -22.2],
  AT: [-7.2, -20.4],
  TA: [-7.2, -21.3],
  CA: [-8.5, -22.7], TG: [-8.5, -22.7],
  GT: [-8.4, -22.4], AC: [-8.4, -22.4],
  CT: [-7.8, -21.0], AG: [-7.8, -21.0],
  GA: [-8.2, -22.2], TC: [-8.2, -22.2],
  CG: [-10.6, -27.2],
  GC: [-9.8, -24.4],
  GG: [-8.0, -19.9], CC: [-8.0, -19.9],
};

export function meltingTemp(s, { na = 0.05, oligo = 250e-9 } = {}) {
  const q = s.toUpperCase().replace(/[^ACGT]/g, '');
  if (q.length < 2) return 0;
  let dH = 0;
  let dS = 0;
  for (let i = 0; i < q.length - 1; i++) {
    const p = NN[q.slice(i, i + 2)];
    dH += p[0];
    dS += p[1];
  }
  for (const end of [q[0], q[q.length - 1]]) {
    if (end === 'A' || end === 'T') { dH += 2.2; dS += 6.9; }
  }
  dS += 0.368 * (q.length - 1) * Math.log(na);
  const R = 1.987;
  const selfComp = q === revcomp(q);
  if (selfComp) dS += -1.4;
  const ct = selfComp ? oligo : oligo / 4;
  return (dH * 1000) / (dS + R * Math.log(ct)) - 273.15;
}

// ---- ORFs ---------------------------------------------------------------------

// Open reading frames (ATG ... stop) of at least minCodons codons, both strands.
// For each stop codon only the longest ORF (first ATG) is reported.
function scanOrfs(seq, circular, minCodons = 100) {
  const L = seq.length;
  const out = [];
  const scan = (s, strand) => {
    const hay = circular ? s + s + s.slice(0, 2) : s;
    const limit = circular ? 2 * L : L;
    const seenStops = new Set();
    for (let frame = 0; frame < 3; frame++) {
      let atg = -1;
      for (let i = frame; i + 3 <= hay.length && i < limit; i += 3) {
        // compare character codes: slicing a codon per step made this the slowest scan on long sequences
        const a = hay.charCodeAt(i);
        const b = hay.charCodeAt(i + 1);
        const g = hay.charCodeAt(i + 2);
        if (atg < 0 && a === 65 && b === 84 && g === 71) atg = i;
        if (a === 84 && ((b === 65 && (g === 65 || g === 71)) || (b === 71 && g === 65))) {
          if (atg >= 0 && atg < L) {
            const len = i + 3 - atg;
            const stopKey = (i + 3) % L;
            if (len / 3 - 1 >= minCodons && len <= L && !seenStops.has(stopKey)) {
              seenStops.add(stopKey);
              out.push({ s: atg, len, strand });
            }
          }
          atg = -1;
        }
      }
    }
  };
  const up = seq.toUpperCase();
  scan(up, 1);
  scan(revcomp(up), -1);
  return out.map(({ s, len, strand }) => {
    let start = strand === 1 ? s : L - s - len;
    start = mod(start, L);
    let end = start + len;
    if (circular) end = mod(end, L);
    return { start, end, strand, codons: len / 3 - 1 };
  }).sort((a, b) => b.codons - a.codons);
}

/** Cached: six-frame scans are the slowest thing a redraw can ask for. */
export const findOrfs = memoLatest(scanOrfs, (circular, minCodons = 100) => `${circular ? 1 : 0}|${minCodons}`);
