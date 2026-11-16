// Restriction enzymes, site finding and digestion.
//
// Each enzyme: [name, site, top, bottom]. `top`/`bottom` are cut positions
// counted from the first base of the site on the top strand, e.g. EcoRI
// G^AATT_C -> top 1, bottom 5. Type IIS enzymes cut outside the site
// (BsaI GGTCTC(1/5) -> top 7, bottom 11).

const RAW = [
  ['AatII', 'GACGTC', 5, 1], ['Acc65I', 'GGTACC', 1, 5], ['AflII', 'CTTAAG', 1, 5],
  ['AgeI', 'ACCGGT', 1, 5], ['AleI', 'CACNNNNGTG', 5, 5], ['AluI', 'AGCT', 2, 2],
  ['AlwNI', 'CAGNNNCTG', 6, 3], ['ApaI', 'GGGCCC', 5, 1], ['AscI', 'GGCGCGCC', 2, 6],
  ['AvrII', 'CCTAGG', 1, 5], ['BamHI', 'GGATCC', 1, 5], ['BbsI', 'GAAGAC', 8, 12],
  ['BfaI', 'CTAG', 1, 3], ['BglI', 'GCCNNNNNGGC', 7, 4], ['BglII', 'AGATCT', 1, 5],
  ['BlpI', 'GCTNAGC', 2, 5], ['BsaI', 'GGTCTC', 7, 11], ['BsiWI', 'CGTACG', 1, 5],
  ['BsmBI', 'CGTCTC', 7, 11], ['BspEI', 'TCCGGA', 1, 5], ['BspHI', 'TCATGA', 1, 5],
  ['BsrGI', 'TGTACA', 1, 5], ['BssHII', 'GCGCGC', 1, 5], ['BstBI', 'TTCGAA', 2, 4],
  ['BstEII', 'GGTNACC', 1, 6], ['Bsu36I', 'CCTNAGG', 2, 5], ['ClaI', 'ATCGAT', 2, 4],
  ['DpnII', 'GATC', 0, 4], ['DraI', 'TTTAAA', 3, 3], ['DraIII', 'CACNNNGTG', 6, 3],
  ['EagI', 'CGGCCG', 1, 5], ['EcoNI', 'CCTNNNNNAGG', 5, 6], ['EcoRI', 'GAATTC', 1, 5],
  ['EcoRV', 'GATATC', 3, 3], ['FseI', 'GGCCGGCC', 6, 2], ['HaeIII', 'GGCC', 2, 2],
  ['HhaI', 'GCGC', 3, 1], ['HincII', 'GTYRAC', 3, 3], ['HindIII', 'AAGCTT', 1, 5],
  ['HpaI', 'GTTAAC', 3, 3], ['KasI', 'GGCGCC', 1, 5], ['KpnI', 'GGTACC', 5, 1],
  ['MfeI', 'CAATTG', 1, 5], ['MluI', 'ACGCGT', 1, 5], ['MscI', 'TGGCCA', 3, 3],
  ['MseI', 'TTAA', 1, 3], ['MspI', 'CCGG', 1, 3], ['NaeI', 'GCCGGC', 3, 3],
  ['NcoI', 'CCATGG', 1, 5], ['NdeI', 'CATATG', 2, 4], ['NgoMIV', 'GCCGGC', 1, 5],
  ['NheI', 'GCTAGC', 1, 5], ['NlaIII', 'CATG', 4, 0], ['NotI', 'GCGGCCGC', 2, 6],
  ['NruI', 'TCGCGA', 3, 3], ['NsiI', 'ATGCAT', 5, 1], ['PacI', 'TTAATTAA', 5, 3],
  ['PciI', 'ACATGT', 1, 5], ['PflMI', 'CCANNNNNTGG', 7, 4], ['PmeI', 'GTTTAAAC', 4, 4],
  ['PspOMI', 'GGGCCC', 1, 5], ['PstI', 'CTGCAG', 5, 1], ['PvuI', 'CGATCG', 4, 2],
  ['PvuII', 'CAGCTG', 3, 3], ['RsaI', 'GTAC', 2, 2], ['SacI', 'GAGCTC', 5, 1],
  ['SacII', 'CCGCGG', 4, 2], ['SalI', 'GTCGAC', 1, 5], ['SapI', 'GCTCTTC', 8, 11],
  ['SbfI', 'CCTGCAGG', 6, 2], ['ScaI', 'AGTACT', 3, 3], ['SfiI', 'GGCCNNNNNGGCC', 8, 5],
  ['SfoI', 'GGCGCC', 3, 3], ['SmaI', 'CCCGGG', 3, 3], ['SnaBI', 'TACGTA', 3, 3],
  ['SpeI', 'ACTAGT', 1, 5], ['SphI', 'GCATGC', 5, 1], ['SspI', 'AATATT', 3, 3],
  ['StuI', 'AGGCCT', 3, 3], ['SwaI', 'ATTTAAAT', 4, 4], ['TaqI', 'TCGA', 1, 3],
  ['XbaI', 'TCTAGA', 1, 5], ['XhoI', 'CTCGAG', 1, 5], ['XmaI', 'CCCGGG', 1, 5],
  ['XmnI', 'GAANNNNTTC', 5, 5], ['ZraI', 'GACGTC', 3, 3],
];

import { iupacRegexSource, revcomp, mod } from './seq.js';

export const ENZYMES = RAW.map(([name, site, top, bottom]) => ({ name, site, top, bottom }));
export const ENZYME_BY_NAME = Object.fromEntries(ENZYMES.map((e) => [e.name, e]));

// A sensible default set: common 6+ cutters used for cloning.
export const COMMON_SET = [
  'AatII', 'AgeI', 'ApaI', 'AscI', 'AvrII', 'BamHI', 'BglII', 'BsaI', 'BsmBI', 'BsrGI',
  'ClaI', 'EcoRI', 'EcoRV', 'FseI', 'HindIII', 'HpaI', 'KpnI', 'MfeI', 'MluI', 'NcoI',
  'NdeI', 'NheI', 'NotI', 'NsiI', 'PacI', 'PmeI', 'PstI', 'PvuI', 'PvuII', 'SacI',
  'SacII', 'SalI', 'SbfI', 'ScaI', 'SmaI', 'SpeI', 'SphI', 'StuI', 'SwaI', 'XbaI',
  'XhoI', 'XmaI', 'ZraI', 'BstBI', 'AflII', 'SnaBI', 'SfiI',
];

export function isPalindrome(e) {
  return revcomp(e.site) === e.site;
}

export function overhangLabel(e) {
  const d = e.bottom - e.top;
  if (d === 0) return 'blunt';
  return d > 0 ? `5′ overhang (${d})` : `3′ overhang (${-d})`;
}

// Find every cut site of `enzyme`. Returns [{enzyme, site: {start, end}, strand, top, bottom}]
// where top/bottom are absolute cut positions (between base top-1 and top).
export function findSites(seq, circular, enzyme) {
  const L = seq.length;
  const n = enzyme.site.length;
  if (n > L) return [];
  const hay = circular ? seq + seq.slice(0, n - 1) : seq;
  const out = [];
  const scan = (pattern, strand) => {
    const re = new RegExp(`(?=${iupacRegexSource(pattern)})`, 'g');
    let m;
    while ((m = re.exec(hay))) {
      const p = m.index;
      if (p >= L) break;
      let top;
      let bottom;
      if (strand === 1) { top = p + enzyme.top; bottom = p + enzyme.bottom; } else { top = p + n - enzyme.bottom; bottom = p + n - enzyme.top; }
      if (circular) { top = mod(top, L); bottom = mod(bottom, L); } else if (top < 0 || top > L || bottom < 0 || bottom > L) {
        re.lastIndex = p + 1;
        continue;
      }
      out.push({ enzyme, site: { start: p, end: circular ? mod(p + n, L) || L : p + n }, strand, top, bottom });
      re.lastIndex = p + 1;
    }
  };
  scan(enzyme.site, 1);
  if (!isPalindrome(enzyme)) scan(revcomp(enzyme.site), -1);
  return out.sort((a, b) => a.top - b.top);
}

export function findAllSites(seq, circular, enzymes = ENZYMES) {
  const map = new Map();
  for (const e of enzymes) map.set(e.name, findSites(seq, circular, e));
  return map;
}

// Overhang left by a cut: { type: '5'|'3'|'blunt', oh: top-strand bases }
export function cutEnd(seq, circular, cut) {
  const { top, bottom } = cut;
  if (top === bottom) return { type: 'blunt', oh: '' };
  const sub = (a, b) => (circular && b < a ? seq.slice(a) + seq.slice(0, b) : seq.slice(a, b));
  if (circular) {
    const L = seq.length;
    const d = mod(bottom - top, L);
    if (d < L / 2) return { type: '5', oh: sub(top, bottom) };
    return { type: '3', oh: sub(bottom, top) };
  }
  return top < bottom ? { type: '5', oh: seq.slice(top, bottom) } : { type: '3', oh: seq.slice(bottom, top) };
}

// Digest with a list of cuts. Fragments run from one top-strand nick to the next.
// fragment = { start, end, len, left: {cut, end}, right: {cut, end} }  (start/end are top cut positions)
export function digest(seq, circular, cuts) {
  const L = seq.length;
  const sorted = [...cuts].sort((a, b) => a.top - b.top);
  // de-duplicate identical cut positions (isoschizomers)
  const uniq = sorted.filter((c, i) => i === 0 || c.top !== sorted[i - 1].top);
  const natural = (pos) => ({ cut: { top: pos, bottom: pos, enzyme: null }, end: { type: 'natural', oh: '' } });
  const withEnd = (c) => ({ cut: c, end: cutEnd(seq, circular, c) });
  const frags = [];
  if (!uniq.length) {
    return [{ start: 0, end: circular ? 0 : L, len: L, left: natural(0), right: natural(L), uncut: true }];
  }
  if (circular) {
    for (let i = 0; i < uniq.length; i++) {
      const a = uniq[i];
      const b = uniq[(i + 1) % uniq.length];
      const len = uniq.length === 1 ? L : mod(b.top - a.top, L);
      frags.push({ start: a.top, end: b.top, len, left: withEnd(a), right: withEnd(b) });
    }
  } else {
    let prev = natural(0);
    let prevPos = 0;
    for (const c of uniq) {
      if (c.top === 0 || c.top === L) continue;
      frags.push({ start: prevPos, end: c.top, len: c.top - prevPos, left: prev, right: withEnd(c) });
      prev = withEnd(c);
      prevPos = c.top;
    }
    frags.push({ start: prevPos, end: L, len: L - prevPos, left: prev, right: natural(L) });
  }
  return frags;
}

// ---- gel simulation ------------------------------------------------------------

export const LADDERS = {
  '1 kb ladder': [10000, 8000, 6000, 5000, 4000, 3000, 2000, 1500, 1000, 500],
  '100 bp ladder': [1517, 1200, 1000, 900, 800, 700, 600, 500, 400, 300, 200, 100],
  '1 kb plus ladder': [10000, 8000, 6000, 5000, 4000, 3000, 2000, 1500, 1200, 1000, 900, 800, 700, 600, 500, 400, 300, 200, 100],
};

// Relative migration (0 = well, 1 = bottom) for an agarose gel of the given %.
export function migration(bp, agarose = 1) {
  // log-linear model fitted so that the useful range scales with gel percentage
  const hi = Math.log10(20000 / agarose ** 1.5);
  const lo = Math.log10(Math.max(50, 80 * agarose));
  const x = (hi - Math.log10(Math.max(bp, 10))) / (hi - lo);
  return Math.min(1.02, Math.max(0.02, x));
}
