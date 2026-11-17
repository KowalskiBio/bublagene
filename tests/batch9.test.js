import test from 'node:test';
import assert from 'node:assert/strict';
import { History } from '../js/undo.js';
import { searchSequence, searchFeatures, filterFeatures, searchAll, featureTypes, qualifierKeys } from '../js/search.js';
import { saveDocument, openDocument, writeAutosave, readAutosave, supportsFileSystemAccess, supportsOPFS } from '../js/storage.js';
import { Autosave } from '../js/autosave.js';
import { offload, runLocal, runInWorker, canUseWorker } from '../js/compute.js';
import { newDocument, newId } from '../js/model.js';
import { revcomp } from '../js/seq.js';

test('undo/redo coalesces bursts and keeps a bounded history', () => {
  let now = 0;
  const h = new History({ coalesceMs: 100, limit: 3, clock: () => now });
  h.record({ v: 0 }, 'type', { coalesceKey: 'typing', now: 0 });
  now = 50;
  h.record({ v: 1 }, 'type', { coalesceKey: 'typing', now: 50 });
  assert.equal(h.depth, 1);
  now = 200;
  h.record({ v: 2 }, 'type', { coalesceKey: 'typing', now: 200 });
  assert.equal(h.depth, 2);
  const undone = h.undo({ v: 3 });
  assert.deepEqual(undone, { v: 2 });
  assert.equal(h.canRedo(), true);
  const redone = h.redo({ v: 2 });
  assert.deepEqual(redone, { v: 3 });
  assert.deepEqual(h.labels(), ['type', 'type']);
  h.record({ v: 3 }, 'a', { now: 300 });
  h.record({ v: 4 }, 'b', { now: 400 });
  h.record({ v: 5 }, 'c', { now: 500 });
  assert.equal(h.depth, 3);
  h.clear();
  assert.equal(h.canUndo(), false);
});

test('search covers sequence, features and qualifiers', () => {
  const doc = newDocument({ seq: 'AAAGAATTCAAAGAATTCAAA' });
  doc.features = [
    { id: newId(), name: 'promoter A', type: 'promoter', strand: 1, color: '#888', segments: [{ start: 0, end: 5 }],
      qualifiers: [{ name: 'note', values: [{ kind: 'text', value: 'strong expression' }] }], extra: {} },
    { id: newId(), name: 'CDS x', type: 'CDS', strand: -1, color: '#888', segments: [{ start: 6, end: 12 }],
      qualifiers: [{ name: 'gene', values: [{ kind: 'text', value: 'xylA' }] }], extra: {} },
  ];
  const seqHits = searchSequence(doc, 'GAATTC');
  assert.equal(seqHits.length, 2);
  assert.equal(searchSequence(doc, 'GAATTC', { bothStrands: false }).length, 2);
  assert.ok(searchSequence(doc, 'GAATTC').every((h) => h.kind === 'sequence'));
  const rev = searchSequence(newDocument({ seq: revcomp('GATCAG') }), 'GATCAG');
  assert.ok(rev.some((h) => h.strand === -1));
  const features = searchFeatures(doc, 'xylA');
  assert.equal(features.length, 1);
  assert.equal(features[0].kind, 'qualifier');
  assert.ok(searchFeatures(doc, 'CDS').some((h) => h.field === 'type'));
  assert.ok(searchAll(doc, 'promoter').some((h) => h.kind === 'feature'));
  assert.equal(searchSequence(doc, 'ZZZ').length, 0);
  const filtered = filterFeatures(doc, { type: 'CDS' });
  assert.equal(filtered.length, 1);
  assert.equal(filterFeatures(doc, { text: 'strong' }).length, 1);
  assert.equal(filterFeatures(doc, { hasQualifier: 'gene' }).length, 1);
  assert.deepEqual(featureTypes(doc), ['CDS', 'promoter']);
  assert.deepEqual(qualifierKeys(doc), ['gene', 'note']);
});

test('storage prefers the File System Access API, then OPFS, then download', async () => {
  let saved = null;
  const handle = { name: 'out.gb', createWritable: async () => ({ write: async (d) => { saved = d; }, close: async () => {} }) };
  const env = {
    showSaveFilePicker: async () => handle,
    showOpenFilePicker: async () => ({ name: 'in.gb', getFile: async () => ({ name: 'in.gb', arrayBuffer: async () => new TextEncoder().encode('ATG').buffer }) }),
  };
  assert.equal(supportsFileSystemAccess(env), true);
  const out = await saveDocument({ data: 'abc', name: 'out.gb', env });
  assert.equal(out.method, 'file-system-access');
  assert.equal(saved, 'abc');
  const opened = await openDocument({ env });
  assert.equal(opened.name, 'in.gb');
  assert.equal(new TextDecoder().decode(opened.data), 'ATG');

  let opfsSaved = null;
  const opfs = {
    navigator: { storage: { getDirectory: async () => ({ getFileHandle: async () => ({ createWritable: async () => ({ write: async (d) => { opfsSaved = d; }, close: async () => {} }) }) }) } },
  };
  assert.equal(supportsOPFS(opfs), true);
  const opfsOut = await saveDocument({ data: 'xyz', name: 'o.gb', env: opfs });
  assert.equal(opfsOut.method, 'opfs');
  assert.equal(opfsSaved, 'xyz');

  const dl = { URL: { createObjectURL: () => 'blob:1', revokeObjectURL: () => {} }, document: { createElement: () => ({ click: () => {} }) }, setTimeout: () => {} };
  const dlOut = await saveDocument({ data: 'q', name: 'd.gb', env: dl });
  assert.equal(dlOut.method, 'download');
});

test('autosave debounces writes and restores the last payload', async () => {
  let saved = null;
  let timerFn = null;
  const timer = { setTimeout: (fn) => { timerFn = fn; return 1; }, clearTimeout: () => {} };
  const auto = new Autosave({ write: async (p) => { saved = p; }, read: async () => saved, clear: async () => { saved = null; }, timer, delay: 10 });
  auto.schedule('one');
  auto.schedule('two');
  assert.equal(saved, null);
  await timerFn();
  assert.equal(saved, 'two');
  assert.equal((await auto.restore()), 'two');
  await auto.clearSaved();
  assert.equal(saved, null);
  const stored = await writeAutosave('k', { a: 1 }, { localStorage: { setItem: (k, v) => { saved = v; }, getItem: () => saved } });
  assert.equal(stored.method, 'localStorage');
  assert.equal(await readAutosave('k', { localStorage: { getItem: () => saved } }), '{"a":1}');
});

test('heavy work can be offloaded or run locally with the same result', async () => {
  const local = await runLocal((x) => x * 2, 3);
  assert.deepEqual(local, { result: 6, worker: false });
  const fallback = await offload((x) => x + 1, 4, {});
  assert.deepEqual(fallback, { result: 5, worker: false });
  class FakeWorker {
    constructor() { this.onmessage = null; this.onerror = null; }
    postMessage(p) { queueMicrotask(() => this.onmessage({ data: p * 3 })); }
    terminate() {}
  }
  const remote = await offload((x) => x, 4, { createWorker: () => new FakeWorker() });
  assert.deepEqual(remote, { result: 12, worker: true });
  assert.equal(canUseWorker({ Worker: FakeWorker }), true);
  const direct = await runInWorker(new FakeWorker(), 3);
  assert.deepEqual(direct, { result: 9, worker: true });
});
