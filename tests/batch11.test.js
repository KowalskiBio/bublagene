import test from 'node:test';
import assert from 'node:assert/strict';
import { seqDigest, SequenceCache, markAlignmentStale, alignmentIsCurrent, stampAlignment, getIndex } from '../js/cache.js';
import { createFastaParser, streamFasta, readFasta } from '../js/formats/genbank.js';
import { newDocument, replaceRange, rotateOrigin } from '../js/model.js';
import { gcContent, translate, findPattern } from '../js/seq.js';
import { performance } from 'node:perf_hooks';

test('alignment records go stale when the sequence changes', () => {
  const doc = newDocument({ seq: 'ATGC'.repeat(25) });
  stampAlignment(doc, { score: 10, identity: 100 });
  assert.equal(alignmentIsCurrent(doc), true);
  replaceRange(doc, 10, 14, 'GG');
  assert.equal(alignmentIsCurrent(doc), false);
  assert.equal(doc.alignment.stale, true);
  const before = doc.alignment.currentDigest;
  doc.alignment.digest = before;
  assert.equal(alignmentIsCurrent(doc), true);
  rotateOrigin(newDocument({ seq: 'ATGC', circular: true }), 1);
});

test('sequence caches key by digest and recompute after an edit', () => {
  const doc = newDocument({ seq: 'ATGCATGC' });
  const cache = new SequenceCache();
  assert.equal(cache.has('gc', doc.seq), false);
  const first = cache.get('gc', doc.seq, (s) => gcContent(s));
  assert.equal(cache.has('gc', doc.seq), true);
  assert.equal(cache.get('gc', doc.seq, () => -1), first);
  replaceRange(doc, 0, 0, 'GG');
  assert.equal(cache.has('gc', doc.seq), false);
  assert.notEqual(cache.get('gc', doc.seq, (s) => gcContent(s)), first);
  const idx = getIndex(doc, 'length', (s) => s.length);
  const before = doc.seq.length;
  doc.seq = doc.seq + 'A';
  assert.equal(getIndex(doc, 'length', (s) => s.length), before + 1);
  cache.invalidate('gc');
  assert.equal(cache.has('gc', doc.seq), false);
  cache.clear();
});

test('streaming FASTA handles chunk boundaries and a memory cap', async () => {
  const parser = createFastaParser();
  const out = [];
  out.push(...parser.push('>a desc\nATG'));
  out.push(...parser.push('C\n>b\nTT\n'));
  out.push(...parser.end());
  assert.equal(out.length, 2);
  assert.equal(out[0].name, 'a');
  assert.equal(out[0].seq, 'ATGC');
  assert.equal(out[1].seq, 'TT');
  assert.throws(() => createFastaParser({ maxRecordLength: 3 }).push('>x\nATGC\n'), /exceeds 3/);

  async function* chunks() { yield '>one\nAT'; yield 'GC\n>two\nTT\n'; }
  const streamed = [];
  for await (const doc of streamFasta(chunks())) streamed.push(doc);
  assert.deepEqual(streamed.map((d) => [d.name, d.seq]), [['one', 'ATGC'], ['two', 'TT']]);
  assert.equal(readFasta('>x\nATGC\n')[0].seq, 'ATGC');
});

test('a 5 Mb genome stays within the performance budget', () => {
  let x = 12345;
  let seq = '';
  const bases = 'ACGT';
  while (seq.length < 5_000_000) {
    x = (1103515245 * x + 12345) >>> 0;
    seq += bases[(x >>> 16) & 3];
  }
  const start = performance.now();
  const gc = gcContent(seq);
  const protein = translate(seq.slice(0, 300000));
  const hits = findPattern(seq, 'GAATTC', false);
  const ms = performance.now() - start;
  assert.ok(gc >= 0 && gc <= 100);
  assert.ok(protein.length > 0);
  assert.ok(hits.length >= 0);
  assert.ok(ms < 15000, `5 Mb operations took ${ms.toFixed(0)} ms`);
});
