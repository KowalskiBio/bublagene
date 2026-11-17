// Things drawn on the map and sequence view that are not features of the
// document: open reading frames (hollow, dashed) and primer binding sites.
// Picking one selects its span.

import { findOrfs } from '../seq.js';
import { findBindingSites } from '../primers.js';
import { colorFor } from '../model.js';

export const ORF_PREFIX = 'orf:';
export const isOrf = (f) => !!f && f.orf === true;
export const isPrimer = (f) => !!f && f.primer === true;
export const isOverlay = (f) => isOrf(f) || isPrimer(f);

let cache = null;

// findOrfs scans six frames, so keep the last answer while the sequence is unchanged
export function orfsFor(doc, minCodons) {
  if (cache && cache.seq === doc.seq && cache.circular === doc.circular && cache.min === minCodons) return cache.orfs;
  const orfs = findOrfs(doc.seq, doc.circular, minCodons);
  cache = { seq: doc.seq, circular: doc.circular, min: minCodons, orfs };
  return orfs;
}

export function orfFeatures(orfs) {
  return orfs.map((o, i) => ({
    id: `${ORF_PREFIX}${i}`,
    name: `ORF ${o.codons} aa`,
    type: 'ORF',
    strand: o.strand,
    color: 'transparent',
    orf: true,
    orfInfo: o,
    segments: [{ start: o.start, end: o.end }],
    qualifiers: [],
    extra: {},
  }));
}

// One bar per place each primer anneals; a 5' tail that does not match is left
// out of the bar, which covers only the annealed bases.
export function primerOverlays(doc) {
  const out = [];
  for (const p of doc.primers ?? []) {
    findBindingSites(doc.seq, doc.circular, p.seq).forEach((site, k) => {
      out.push({
        id: `primer:${p.id}:${k}`,
        name: p.name,
        type: 'primer_bind',
        strand: site.strand,
        color: colorFor('primer_bind'),
        primer: true,
        primerId: p.id,
        segments: [{ start: site.start, end: site.end }],
        qualifiers: [],
        extra: {},
      });
    });
  }
  return out;
}
