// @ts-nocheck — legacy DOM wiring; the engines it calls are type-checked.
// What the map and sequence view draw, worked out from the page state.

import { ENZYME_BY_NAME, findSites } from '../enzymes.js';
import { alignments, state } from './core.js';
import { orfFeatures, orfsFor, primerOverlays } from './overlays.js';

export function shownSites(doc) {
  const out = [];
  for (const name of state.shownEnzymes) {
    const e = ENZYME_BY_NAME[name];
    if (e) out.push(...findSites(doc.seq, doc.circular, e));
  }
  return out.sort((a, b) => a.top - b.top);
}

// ORFs and primer sites to draw on top of the document's own features
export const overlaysFor = (doc) => {
  if (!doc) return [];
  return [
    ...(state.showOrfs ? orfFeatures(orfsFor(doc, state.orfMin)) : []),
    ...(state.showPrimers ? primerOverlays(doc) : []),
  ];
};

export const readsFor = (doc) => (doc && state.showReads ? alignments.get(doc.id) ?? [] : []);
