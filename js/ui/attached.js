// @ts-nocheck — session and project plumbing.
// The sequencing reads and chromatograms a reader opened against a document live
// beside it, not inside it. A session or project that kept only the documents
// would bring back the reference and lose every read, so they are collected here
// and put back with it. Reads that come out of a .dna file's own segments are not
// kept: the file re-supplies them.

import { alignments, state, traces } from './core.js';

const rowForStorage = (row, { plain }) => ({
  ...row,
  insertions: row.insertions instanceof Map ? (plain ? [...row.insertions.entries()] : row.insertions) : (row.insertions ?? []),
});

/** What the documents currently open carry that is not part of them. */
export function collectAttached({ plain = false, withTraces = true, onlyAttached = true } = {}) {
  const ids = new Set(state.docs.map((d) => d.id));
  const rows = [];
  for (const [docId, list] of alignments) {
    if (!ids.has(docId) || !list?.length) continue;
    const keep = list.filter((r) => !onlyAttached || r.attached);
    if (keep.length) rows.push({ docId, reads: keep.map((r) => rowForStorage(r, { plain })) });
  }
  const found = [];
  if (withTraces) {
    for (const [docId, list] of traces) {
      if (!ids.has(docId) || !list?.length) continue;
      const keep = list.filter((t) => !onlyAttached || t.attached);
      if (keep.length) found.push({ docId, traces: keep });
    }
  }
  return { alignments: rows, traces: found };
}

/** Put saved reads and chromatograms back, marked as the reader's own. */
export function restoreAttached({ alignments: rows = [], traces: found = [] } = {}) {
  for (const { docId, reads } of rows) {
    const restored = (reads ?? []).map((r) => ({
      ...r, attached: true, insertions: r.insertions instanceof Map ? r.insertions : new Map(r.insertions ?? []),
    }));
    alignments.set(docId, [...(alignments.get(docId) ?? []).filter((r) => !r.attached), ...restored]
      .sort((a, b) => a.start - b.start || String(a.label).localeCompare(String(b.label))));
  }
  for (const { docId, traces: list } of found) {
    traces.set(docId, [...(traces.get(docId) ?? []).filter((t) => !t.attached), ...(list ?? []).map((t) => ({ ...t, attached: true }))]);
  }
}
