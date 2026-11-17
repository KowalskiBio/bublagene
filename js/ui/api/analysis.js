// @ts-nocheck — legacy DOM wiring; the engines it calls are type-checked.
// The panels' view of the app: display options, protein and structures, contigs, search, saved collections, files and export.

import { getFromCollection, listCollection, removeFromCollection, saveToCollection, searchCollection } from '../../collection.js';
import { assembleContigs } from '../../contig.js';
import { addHistory, cloneDoc, colorFor, featureSpan, newDocument, newId } from '../../model.js';
import { multipleAlign } from '../../msa.js';
import { projectToJSON } from '../../project.js';
import { searchAll } from '../../search.js';
import { saveDocument } from '../../storage.js';
import { current, el, render, state, traces } from '../core.js';
import { openDoc } from '../documents.js';
import { fmt, toast } from '../dom.js';
import { download, exportableSvg, pdfBlob, pngBlob, svgBlob } from '../export.js';
import { openProject, save, saveProject } from '../files.js';

export function analysisApi(api) {
  return {
    setMapStyle(v) { state.mapStyle = v; render(); },
    setPlot(patch) { state.plot = { ...state.plot, ...patch }; render(); },
    // ---- protein and structures ----
    setProtein(patch) { state.protein = { ...state.protein, ...patch }; render(); },
    pickFile() { el('fileInput').click(); },
    removeStructure(id) { state.structures = state.structures.filter((x) => x.id !== id); render(); },
    // ---- contigs and multiple alignment ----
    traceCount(doc) { return traces.get(doc.id)?.length ?? 0; },
    setContigs(patch) { state.contigs = { ...state.contigs, result: null, error: null, ...patch }; render(); },
    runContigs() {
      const c = state.contigs;
      const reads = [];
      for (const id of c.ids) {
        const d = state.docs.find((x) => x.id === id);
        if (!d) continue;
        const found = c.useTraces ? traces.get(d.id) : null;
        if (found?.length) found.forEach((t, i) => reads.push({ name: `${d.name}:${t.name || i + 1}`, seq: t.bases, qualities: t.confidence }));
        else reads.push({ name: d.name, seq: d.seq });
      }
      try {
        if (!reads.length) throw new Error('Pick the files to assemble first');
        const all = assembleContigs(reads, { minOverlap: c.minOverlap, maxMismatchRate: c.maxMismatch / 100 });
        const joined = all.filter((x) => x.reads.length > 1);
        const contigs = joined.length ? joined : all;
        const unplaced = joined.length ? all.filter((x) => x.reads.length === 1).flatMap((x) => x.reads) : [];
        state.contigs = { ...c, result: { contigs, unplaced, readCount: reads.length }, error: null };
      } catch (e) { state.contigs = { ...c, result: null, error: e.message }; }
      render();
    },
    openContig(i) {
      const ct = state.contigs.result?.contigs[i];
      if (!ct) return;
      const doc = newDocument({ name: `Contig ${i + 1}`, seq: ct.seq, circular: false });
      for (const r of ct.layout) {
        doc.features.push({
          id: newId(), name: r.name, type: 'misc_feature', strand: r.reverse ? -1 : 1, color: colorFor('misc_feature'),
          segments: [{ start: r.offset, end: Math.min(ct.length, r.offset + r.length) }], qualifiers: [], extra: {},
        });
      }
      addHistory(doc, 'assembled', `${ct.reads.length} reads, ${ct.length} bp, ${ct.depth}× deep`);
      openDoc(doc);
    },
    setMsa(patch) { state.msa = { ...state.msa, result: null, error: null, ...patch }; render(); },
    runMsa() {
      const items = state.msa.ids.map((id) => state.docs.find((d) => d.id === id)).filter(Boolean).map((d) => ({ name: d.name, seq: d.seq }));
      try { state.msa = { ...state.msa, result: multipleAlign(items), error: null }; } catch (e) { state.msa = { ...state.msa, result: null, error: e.message }; }
      render();
    },
    // ---- search ----
    setSearch(patch) { state.search = { ...state.search, ...patch }; render(); },
    runSearch(query) {
      const q = String(query ?? '').trim();
      const scope = state.search.scope;
      state.search = { ...state.search, query: q, results: null, saved: null, error: null };
      if (!q) { render(); return; }
      const docsToSearch = scope === 'file' ? [current()] : scope === 'open' ? state.docs : [];
      const results = docsToSearch.filter(Boolean).map((d) => ({
        docId: d.id, name: d.name,
        items: searchAll(d, q).map((hit) => {
          if (hit.kind === 'sequence') return { start: hit.start, end: hit.end, strand: hit.strand, label: 'sequence match' };
          const f = d.features.find((x) => x.id === hit.featureId);
          const span = f ? featureSpan(f) : { start: 0, end: 0 };
          return { start: span.start, end: span.end, strand: f?.strand ?? 0, label: `${hit.name}${hit.kind === 'qualifier' ? ` · ${hit.field}: ${hit.text}` : ''}` };
        }),
      })).filter((r) => r.items.length || scope === 'file');
      state.search.results = scope === 'saved' ? null : results;
      render();
      if (scope === 'saved') {
        searchCollection(q).then((saved) => { state.search = { ...state.search, saved }; render(); },
          (e) => { state.search = { ...state.search, error: `Could not search the saved files: ${e.message}` }; render(); });
      }
    },
    // ---- saved collections ----
    setLibraryName(name) { state.library = { ...state.library, name }; },
    refreshLibrary() {
      state.library = { ...state.library, loading: true, error: null };
      listCollection().then((records) => { state.library = { ...state.library, records, loading: false }; render(); },
        (e) => { state.library = { ...state.library, records: [], loading: false, error: `Saved files are not available here: ${e.message}` }; render(); });
    },
    async saveCurrentToCollection() {
      const doc = current();
      if (!doc) return;
      try { await saveToCollection(doc, state.library.name); toast(`Saved ${doc.name}`); api.refreshLibrary(); } catch (e) { toast(`Could not save: ${e.message}`, 'error'); }
    },
    async saveAllToCollection() {
      try { for (const d of state.docs) await saveToCollection(d, state.library.name); toast(`Saved ${state.docs.length} files`); api.refreshLibrary(); } catch (e) { toast(`Could not save: ${e.message}`, 'error'); }
    },
    async openSaved(id) {
      try {
        const doc = await getFromCollection(id);
        if (!doc) { toast('That file is no longer saved', 'error'); return; }
        const copy = cloneDoc(doc);
        copy.id = newId(); // the saved copy stays as it was
        openDoc(copy);
      } catch (e) { toast(`Could not open: ${e.message}`, 'error'); }
    },
    async deleteSaved(id) { try { await removeFromCollection(id); api.refreshLibrary(); } catch (e) { toast(`Could not remove: ${e.message}`, 'error'); } },
    async exportCollection() {
      try {
        const records = state.library.records ?? [];
        const docs = [];
        for (const r of records) { const d = await getFromCollection(r.id); if (d) docs.push(d); }
        const result = await saveDocument({ data: projectToJSON(docs, { active: 0, settings: {} }), name: `${state.library.name.replace(/[^\w.-]+/g, '_') || 'collection'}.json`, types: [{ description: 'BublaGene project', accept: { 'application/json': ['.json'] } }] });
        toast(result.method === 'download' ? 'Project downloaded' : 'Project saved');
      } catch (e) { toast(`Could not export: ${e.message}`, 'error'); }
    },
    // ---- export and print ----
    async exportMap(kind) {
      const doc = current();
      const svg = el('view').querySelector('.mapwrap svg');
      if (!doc || !svg) { toast('Switch to the map view to export it', 'error'); return; }
      const exp = exportableSvg(svg);
      const base = doc.name.replace(/[^\w.-]+/g, '_') || 'map';
      try {
        if (kind === 'svg') download(svgBlob(exp), `${base}.svg`);
        else if (kind === 'png') download(await pngBlob(exp), `${base}.png`);
        else if (kind === 'pdf') download(await pdfBlob(exp, { title: doc.name, subtitle: `${fmt(doc.seq.length)} bp · ${doc.circular ? 'circular' : 'linear'}` }), `${base}.pdf`);
        toast(`Exported ${base}.${kind}`);
      } catch (e) { toast(`Could not export: ${e.message}`, 'error'); }
    },
    print() { window.print(); },
    save(kind) { save(kind); },
    saveProject,
    openProject,
  };
}
