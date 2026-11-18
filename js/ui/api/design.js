// @ts-nocheck — legacy DOM wiring; the engines it calls are type-checked.
// The panels' view of the app: Cloning, feature detection, alignment, CRISPR and RNA folding.

import { alignPair } from '../../align.js';
import { offload } from '../../compute.js';
import { editPrimers, findGuides, guidesToFeatures, hdrEdit, nhejEdit } from '../../crispr.js';
import { BUILTIN_LIBRARY, detectFeatures, hitToFeature, libraryFromDocs } from '../../detect.js';
import { fold, structureToSvg } from '../../fold.js';
import { newId } from '../../model.js';
import { subseq } from '../../seq.js';
import { callVariants, variantsToFeatures } from '../../variants.js';
import { selectedSeq } from '../clipboard.js';
import { placeRead, rowFromPairwise } from '../../readalign.js';
import { current, render, state, traces } from '../core.js';
import { attachRead, edit, openDoc } from '../documents.js';
import { h, toast } from '../dom.js';
import { loadUserLibrary, storeUserLibrary } from '../local-libraries.js';

export function designApi(api) {
  return {
    clearCrispr() { state.crispr = null; render(); },
    setCloning(patch) { state.cloning = { ...state.cloning, ...patch }; render(); },
    userLibrary() { return loadUserLibrary(); },
    setDetectOptions(patch) { Object.assign(state, patch); render(); },
    async runDetect() {
      const doc = current();
      if (!doc || state.detectBusy) return;
      const others = state.docs.filter((d) => d.id !== doc.id);
      const library = [...BUILTIN_LIBRARY, ...loadUserLibrary(), ...(state.detectUseDocs ? libraryFromDocs(others) : [])];
      const { seq } = doc;
      state.detectBusy = true;
      render();
      let hits = null;
      let failure = null;
      try {
        // a long sequence takes a while; do it off the page's thread
        const { result } = await offload(
          ({ doc: d, library: lib, options }) => detectFeatures(d, lib, options),
          { kind: 'detect', doc: { seq: doc.seq, circular: doc.circular, features: doc.features }, library, options: { maxMismatchRate: state.detectRate / 100 } },
          { workerUrl: new URL('../../workers/compute-worker.js', import.meta.url), timeoutMs: 120000 },
        );
        if (result?.__error) throw new Error(result.__error);
        hits = result.map((h, i) => ({ ...h, key: i, checked: true }));
      } catch (e) { failure = e.message; }
      state.detectBusy = false;
      // the reader may have edited or closed the file while this ran; its answer no longer fits
      if (current()?.id === doc.id && current().seq === seq) {
        if (failure) toast(`Could not detect: ${failure}`, 'error');
        else state.detect = { docId: doc.id, hits, libraryCount: library.length };
      }
      render();
    },
    toggleDetected(key) {
      const h = state.detect?.hits.find((x) => x.key === key);
      if (h) { h.checked = !h.checked; render(); }
    },
    setDetectedAll(on) { for (const h of state.detect?.hits ?? []) h.checked = on; render(); },
    addDetected() {
      const doc = current();
      const picked = state.detect?.docId === doc?.id ? state.detect.hits.filter((h) => h.checked) : [];
      if (!picked.length) return;
      edit(`Add ${picked.length} detected feature${picked.length === 1 ? '' : 's'}`, (d) => { d.features.push(...picked.map(hitToFeature)); });
      toast(`Added ${picked.length} feature${picked.length === 1 ? '' : 's'}`);
      state.panel = 'Features'; // edit() has already cleared the detection results
      render();
    },
    saveSelectionToLibrary() {
      const seq = selectedSeq();
      if (seq.length < 12) { toast('Select at least 12 bases to save as an element', 'error'); return; }
      const name = window.prompt('Name this element:');
      if (!name?.trim()) return;
      const type = window.prompt('Feature type (promoter, terminator, rep_origin, misc_feature …):', 'misc_feature') || 'misc_feature';
      const list = loadUserLibrary().filter((e) => e.name !== name.trim());
      list.push({ name: name.trim(), type: type.trim(), kind: 'dna', seq });
      if (!storeUserLibrary(list)) { toast('This browser would not keep the element', 'error'); return; }
      toast(`Saved “${name.trim()}” to your library`);
      render();
    },
    removeFromLibrary(name) { storeUserLibrary(loadUserLibrary().filter((e) => e.name !== name)); render(); },
    setAlignTarget(id) { state.alignTarget = id; state.alignResult = null; render(); },
    setAlignLocal(v) { state.alignLocal = v; render(); },
    setAlignBand(v) { state.alignBand = v; render(); },
    async runAlignment() {
      const doc = current();
      const others = state.docs.filter((d) => d.id !== doc?.id);
      const target = others.find((d) => d.id === state.alignTarget) ?? others[0];
      if (!doc || !target) return;
      state.alignBusy = true;
      render();
      try {
        const options = { local: state.alignLocal, band: state.alignBand > 0 ? state.alignBand : null };
        const { result } = await offload(
          ({ a, b, options: opts }) => alignPair(a, b, opts),
          { kind: 'align', a: doc.seq, b: target.seq, options },
          { workerUrl: new URL('../../workers/compute-worker.js', import.meta.url) },
        );
        if (result?.__error) throw new Error(result.__error);
        state.alignResult = { docId: doc.id, targetId: target.id, result };
        // show it in the Alignment view of the reference: its map, this sequence aligned to it, and the chromatogram if there is one
        const own = (traces.get(doc.id) ?? []).find((t) => t.bases?.length === doc.seq.length) ?? null;
        const placed = own && doc.seq.length <= 3000 ? placeRead(target.seq, own, { label: doc.name }) : null;
        const row = placed ?? rowFromPairwise(result, doc.seq.length, doc.name);
        attachRead(target, row, own ? { ...own } : null);
      } catch (e) {
        state.alignResult = { docId: doc.id, targetId: target.id, error: e.message };
      }
      state.alignBusy = false;
      render();
    },
    alignmentVariants() {
      const doc = current();
      const result = state.alignResult;
      const target = state.docs.find((d) => d.id === result?.targetId);
      if (!doc || !target) return;
      const calls = callVariants(target.seq, doc.seq);
      edit('Add variant features', (d) => { d.features.push(...variantsToFeatures(d, calls.variants)); });
      toast(`${calls.variants.length} variant(s) added as features`);
    },
    runCrispr(opts) {
      const doc = current();
      if (!doc) return;
      try {
        const guides = findGuides(doc.seq, {
          enzyme: opts.enzyme, guideLength: opts.guideLength, minGC: opts.minGC, maxGC: opts.maxGC,
          circular: doc.circular, maxOffTargets: opts.offTargets ? 3 : null,
        });
        state.crispr = { docId: doc.id, ...opts, guides, error: null };
      } catch (e) {
        state.crispr = { docId: doc.id, ...opts, guides: [], error: e.message };
      }
      render();
    },
    addCrisprFeatures() {
      const doc = current();
      const guides = state.crispr?.docId === doc?.id ? state.crispr.guides : [];
      if (!guides.length) return;
      edit('Add CRISPR guide features', (d) => { d.features.push(...guidesToFeatures(d, guides)); });
    },
    simulateNhej(guide) {
      const doc = current();
      if (!doc) return;
      const result = nhejEdit(doc, guide, { seed: 1 });
      openDoc(result.doc);
      toast(`NHEJ: deleted ${result.change.deleted.length} bp`);
    },
    simulateHdr(guide, donor) {
      const doc = current();
      if (!doc) return;
      const result = hdrEdit(doc, guide, { donor });
      openDoc(result.doc);
      toast(`HDR: inserted ${donor.length} bp`);
    },
    crisprPrimers(guide) {
      const doc = current();
      if (!doc) return;
      try {
        const primers = editPrimers(doc, guide.cut, { arm: 40 });
        edit('Add edit primers', (d) => {
          d.primers.push({ id: newId(), name: 'edit forward', seq: primers.forward.sequence, description: `Flanks CRISPR cut at ${guide.cut + 1}`, extra: {} });
          d.primers.push({ id: newId(), name: 'edit reverse', seq: primers.reverse.sequence, description: `Flanks CRISPR cut at ${guide.cut + 1}`, extra: {} });
        });
        state.panel = 'Primers';
        render();
      } catch (e) { toast(e.message, 'error'); }
    },
    async runFold({ start, end, maxLength }) {
      const doc = current();
      if (!doc) return;
      const seq = start != null && end != null ? subseq(doc.seq, start, end, doc.circular) : doc.seq;
      state.fold = { docId: doc.id, busy: true, maxLength };
      render();
      try {
        if (!seq.length) throw new Error('Nothing to fold');
        const { result, worker } = await offload(
          ({ sequence, limit }) => fold(sequence, { maxLength: limit }),
          { sequence: seq, limit: maxLength },
          { workerUrl: new URL('../../workers/compute-worker.js', import.meta.url) },
        );
        if (result?.__error) throw new Error(result.__error);
        state.fold = { docId: doc.id, result, worker, maxLength, busy: false };
      } catch (e) {
        state.fold = { docId: doc.id, error: e.message, busy: false, maxLength };
      }
      render();
    },
    exportStructure() {
      const result = state.fold?.result;
      if (!result) return;
      const svg = structureToSvg(result.sequence, result.structure, { width: 900, perRow: 60 });
      const blob = new Blob([svg], { type: 'image/svg+xml' });
      const url = URL.createObjectURL(blob);
      const a = h('a', { href: url, download: 'structure.svg' });
      document.body.append(a); a.click(); a.remove();
      setTimeout(() => URL.revokeObjectURL(url), 2000);
    },
  };
}
