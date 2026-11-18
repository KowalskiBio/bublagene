// @ts-nocheck — legacy DOM wiring; the engines it calls are type-checked.
// BublaGene: shared page state.

import { COMMON_SET } from '../enzymes.js';

export const state = {
  docs: [],
  active: 0,
  view: 'map',
  panel: 'Features',
  selection: null,
  caret: null,        // insertion point in the sequence view, when nothing is selected
  highlight: null,
  enzymeSet: COMMON_SET,
  enzymeEdit: false,   // the set editor is open in the Enzymes panel
  enzymeFilter: '',
  shownEnzymes: new Set(),
  digestWith: [],
  agarose: 1,
  ladder: '1 kb ladder',
  showTranslation: true,
  showSixFrame: false,
  showOrfs: false,
  showPrimers: true,
  orfMin: 100,
  showComplement: true,
  featureFilter: '',
  showReads: true,
  trace: 0,
  alignRow: 0,      // which aligned sequence the Alignment view shows the chromatogram of
  alignCw: 10,      // pixels per column in the Alignment view
  alignScroll: null,
  traceBase: null,
  selectedFeatures: new Set(),
  alignTarget: null,
  alignLocal: true,
  alignBand: 0,
  alignBusy: false,
  alignResult: null,
  crispr: null,
  fold: null,
  clipboard: null,
  mapStyle: 'auto',     // 'auto' follows the file; or force 'circular' / 'linear'
  gelLanes: [],         // saved lanes for the digest gel: { id, docId, label, enzymes }
  tm: { na: 50, mg: 0, oligo: 250 }, // melting-temperature conditions: mM, mM, nM
  primerPair: { a: null, b: null },
  plot: { track: 'gc', window: 100 },
  protein: { feature: null, structure: null, colour: 'structure', highlight: '' },
  structures: [],       // opened PDB files: { id, structure }
  contigs: { ids: [], minOverlap: 30, maxMismatch: 10, useTraces: true, result: null, error: null },
  msa: { ids: [], result: null, error: null },
  search: { query: '', scope: 'file', results: null, saved: null, error: null },
  library: { records: null, loading: false, name: 'My files', error: null },
  mutate: { mode: 'aa', feature: null, residue: 1, newAa: 'A', bases: '', method: 'overlap', targetTm: '', result: null, error: null },
  detect: null,        // last feature-detection run: { docId, hits, rate, useDocs }
  detectRate: 5,       // allowed mismatches, percent
  detectUseDocs: true,
  featureClip: null,   // annotations copied from the Features panel
  pxPerBase: 14,
  cloning: {
    mode: 'restriction',
    vector: { docId: null, enzymes: [], fragment: null },
    insert: { docId: null, enzymes: [], fragment: null },
    circular: true,
    orientation: 'auto',
    gibson: [],
    minOverlap: 20,
    infusionOverlap: 15,
    goldenGate: { enzyme: 'BsaI', docs: [] },
    recomb: { insert: null, acceptor: null, byproduct: null },
    topo: { insert: null, kind: 'TA', orientation: 'forward' },
  },
  restorable: null, // a saved session found on startup, until the user answers
};

// Undo history per document: doc id -> History.
export const undo = new Map();
// Aligned reads and chromatograms are decoded from a file's own bytes rather than
// stored on the document, so they are never written to the session or saved back.
export const alignments = new Map(); // doc id -> reads
export const traces = new Map();     // doc id -> chromatograms

export const el = (id) => document.getElementById(id);
export const current = () => state.docs[state.active] ?? null;
export const isDirty = (doc) => !!(undo.get(doc.id)?.canUndo?.() || doc.restoredDirty);

// Almost everything that changes state ends by redrawing the page, but the code
// that draws it also calls the code that changes state. Registering the renderer
// here, rather than importing it, keeps those two from importing each other.
let renderer = () => {};
export function render() { renderer(); }
export function setRenderer(fn) { renderer = fn; }
