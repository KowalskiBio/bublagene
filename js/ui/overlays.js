// Things drawn on the map and sequence view that are not features of the
// document: open reading frames (hollow, dashed) and primer binding sites.
// Picking one selects its span.

import { findOrfs } from '../seq.js';
import { findBindingSites } from '../primers.js';
import { colorFor } from '../model.js';
import { deferredScan } from './deferred.js';

export const ORF_PREFIX = 'orf:';
export const isOrf = (f) => !!f && f.orf === true;
export const isPrimer = (f) => !!f && f.primer === true;
export const isOverlay = (f) => isOrf(f) || isPrimer(f);

// findOrfs scans six frames. For a long sequence the last answer stays on screen
// while the reader types; see deferred.js.
export function orfsFor(doc, minCodons) {
  return deferredScan(`orfs|${doc.id}|${doc.circular ? 1 : 0}|${minCodons}`, doc.seq, () => findOrfs(doc.seq, doc.circular, minCodons));
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
