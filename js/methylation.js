// Methylation sensitivity of restriction enzymes.
//
// Only the well-established cases are listed; an enzyme that is not here is
// not claimed to be insensitive, so check the supplier's chart for anything
// critical.
//
//  Dam  methylates the A of GATC (most E. coli strains).
//  Dcm  methylates the inner C of CCWGG (most E. coli strains).
//  CpG  methylation (mammalian DNA, or M.SssI in vitro) of the C in CG.
//
// An effect is 'blocked' when the enzyme's own site contains the methylated
// motif, and 'overlap' when it only matters if the motif overlaps the site.

const RULES = {
  BclI: [['dam', 'blocked']], MboI: [['dam', 'blocked']], DpnII: [['dam', 'blocked']],
  ClaI: [['dam', 'overlap'], ['cpg', 'blocked']], XbaI: [['dam', 'overlap']], TaqI: [['dam', 'overlap']],
  EcoRII: [['dcm', 'blocked']], StuI: [['dcm', 'overlap']],
  HpaII: [['cpg', 'blocked']], SmaI: [['cpg', 'blocked']], XmaI: [['cpg', 'blocked']], NotI: [['cpg', 'blocked']],
  AscI: [['cpg', 'blocked']], BssHII: [['cpg', 'blocked']], MluI: [['cpg', 'blocked']], NruI: [['cpg', 'blocked']],
  EagI: [['cpg', 'blocked']], SacII: [['cpg', 'blocked']], SalI: [['cpg', 'blocked']], BstBI: [['cpg', 'blocked']],
  HhaI: [['cpg', 'blocked']],
};

export const KNOWN_ENZYMES = Object.keys(RULES);

const MOTIF = { dam: /GATC/g, dcm: /CC[AT]GG/g, cpg: /CG/g };
export const METHYLATION_LABEL = { dam: 'Dam', dcm: 'Dcm', cpg: 'CpG' };
export const HOST_NOTE = {
  dam: 'GATC methylated in most E. coli strains (not in dam⁻ strains)',
  dcm: 'CCWGG methylated in most E. coli strains (not in dcm⁻ strains)',
  cpg: 'CG methylated in mammalian DNA, or by M.SssI',
};

// The kinds of methylation an enzyme is known to care about (no sequence needed).
export function methylationKinds(name) {
  return (RULES[name] ?? []).map(([kind]) => kind);
}

// Does the motif sit inside, or overlap, [start, end) of a (possibly circular) sequence?
function motifTouches(seq, circular, kind, start, end) {
  const L = seq.length;
  const reach = 4; // longest motif
  const from = start - reach;
  const hay = circular ? seq + seq.slice(0, reach) : seq;
  const re = new RegExp(MOTIF[kind].source, 'g');
  const lo = circular ? from : Math.max(0, from);
  const win = circular
    ? ((seq + seq).slice(((from % L) + L) % L, ((from % L) + L) % L + (end - start) + 2 * reach))
    : hay.slice(lo, end + reach);
  const base = circular ? from : lo;
  let m;
  while ((m = re.exec(win))) {
    const a = base + m.index;
    const b = a + m[0].length;
    if (a < end && b > start) return true;
    re.lastIndex = m.index + 1;
  }
  return false;
}

/**
 * What methylation does to one cut site of one enzyme.
 * @param {{seq:string,circular?:boolean}} doc
 * @param {string} enzymeName
 * @param {{start:number,end:number}} site recognition site, as in findSites().site
 * @returns {{kind:string,effect:'blocked'|'overlap'}[]} the effects that apply at this site
 */
export function methylationAt(doc, enzymeName, site) {
  const out = [];
  const end = site.end > site.start ? site.end : site.end + doc.seq.length;
  for (const [kind, effect] of RULES[enzymeName] ?? []) {
    if (effect === 'blocked') { out.push({ kind, effect }); continue; }
    if (motifTouches(doc.seq, !!doc.circular, kind, site.start, end)) out.push({ kind, effect });
  }
  return out;
}

// For a list of cut sites: which kinds of methylation actually matter here.
export function methylationSummary(doc, enzymeName, sites) {
  const kinds = new Map();
  for (const s of sites) for (const { kind, effect } of methylationAt(doc, enzymeName, s.site)) {
    const prev = kinds.get(kind);
    if (!prev || (prev === 'overlap' && effect === 'blocked')) kinds.set(kind, effect);
  }
  return [...kinds].map(([kind, effect]) => ({ kind, effect }));
}
