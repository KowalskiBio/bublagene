// Round-trip tests against real .dna files. Point DNA_SAMPLES at a folder of
// .dna files.
import test from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync, readdirSync, statSync, existsSync } from 'node:fs';
import { join } from 'node:path';
import { readDna, writeDna, readSegments } from '../js/formats/dna.js';
import { replaceRange, codingSequence } from '../js/model.js';
import { findBindingSites } from '../js/primers.js';
import { translate } from '../js/seq.js';

const ROOT = process.env.DNA_SAMPLES ?? '';

function walk(dir, out = []) {
  for (const name of readdirSync(dir)) {
    const p = join(dir, name);
    const st = statSync(p);
    if (st.isDirectory()) walk(p, out);
    else if (name.toLowerCase().endsWith('.dna')) out.push(p);
  }
  return out;
}

const files = existsSync(ROOT) ? walk(ROOT) : [];

test('sample files are available', { skip: files.length ? false : `no .dna files under ${ROOT}` }, () => {
  assert.ok(files.length > 0);
});

test('every sample .dna file reads, writes and reads back identically', { skip: !files.length }, () => {
  let checked = 0;
  for (const f of files) {
    const doc = readDna(readFileSync(f), f.split('/').pop());
    assert.ok(doc.seq.length > 0, f);
    const bytes = writeDna(doc, { bindingSitesFor: (p) => findBindingSites(doc.seq, doc.circular, p.seq) });
    const back = readDna(bytes, f);
    assert.equal(back.seq, doc.seq, `sequence ${f}`);
    assert.equal(back.circular, doc.circular, `topology ${f}`);
    assert.equal(back.flags, doc.flags, `flags ${f}`);
    assert.equal(back.features.length, doc.features.length, `feature count ${f}`);
    for (let i = 0; i < doc.features.length; i++) {
      const a = doc.features[i];
      const b = back.features[i];
      assert.equal(b.name, a.name, f);
      assert.equal(b.type, a.type, f);
      assert.equal(b.strand, a.strand, f);
      assert.deepEqual(b.segments.map((s) => [s.start, s.end]), a.segments.map((s) => [s.start, s.end]), `${f} ${a.name}`);
      assert.deepEqual(b.qualifiers, a.qualifiers, `${f} ${a.name} qualifiers`);
    }
    assert.deepEqual(back.primers.map((p) => [p.name, p.seq]), doc.primers.map((p) => [p.name, p.seq]), f);
    assert.equal(back.notes.description ?? '', doc.notes.description ?? '', f);
    // segments we don't model (history, display options, ...) pass through byte-for-byte
    const keep = (segs) => segs.filter((s) => ![0x00, 0x05, 0x06, 0x08, 0x09, 0x0a].includes(s.type)).map((s) => [s.type, Buffer.from(s.data).toString('hex')]);
    assert.deepEqual(keep(readSegments(bytes)), keep(readSegments(new Uint8Array(readFileSync(f)))), `pass-through ${f}`);
    checked++;
  }
  assert.equal(checked, files.length);
});

// Checks splicing (gap segments) and the reading-frame offset against the
// own stored protein sequences. 874 of 877 agree. The three exceptions are the
// AmpR feature in the three pAd Gateway vectors, where the span drawn in the
// file is one codon longer than the protein its GenBank annotation carries;
// that is a disagreement inside those files, not a decoding difference.
const KNOWN_DIFFERENT = new Set([
  'pAd BLOCK-iT-DEST.dna/AmpR', 'pAd CMV V5-DEST.dna/AmpR', 'pAd PL-DEST.dna/AmpR',
]);

test('CDS translations match the protein stored in the file', { skip: !files.length }, () => {
  let matched = 0;
  const differ = [];
  for (const f of files) {
    const name = f.split('/').pop();
    const doc = readDna(readFileSync(f), name);
    for (const ft of doc.features.filter((x) => x.type === 'CDS')) {
      const q = ft.qualifiers.find((x) => x.name === 'translation');
      if (!q) continue;
      const want = q.values.map((v) => v.value).join('').replace(/[^A-Z*]/g, '');
      const got = translate(codingSequence(doc, ft));
      // some records store the protein without its stop codon
      if (got === want || got.replace(/\*$/, '') === want) matched++;
      else if (!KNOWN_DIFFERENT.has(`${name}/${ft.name}`)) {
        differ.push(`${name} ${ft.name}: expected ${want.slice(0, 20)}… got ${got.slice(0, 20)}…`);
      }
    }
  }
  assert.deepEqual(differ, []);
  assert.ok(matched >= 874, `only ${matched} of the expected 874 translations matched`);
});

test('editing a sample drops stale sequence caches and still round-trips', { skip: !files.length }, () => {
  const f = files.find((x) => x.endsWith('pIB2-SEC13-mEGFP.dna')) ?? files[0];
  const doc = readDna(readFileSync(f), 'x.dna');
  const featuresBefore = doc.features.length;
  replaceRange(doc, 0, 0, 'GATTACA');
  const bytes = writeDna(doc);
  const types = readSegments(bytes).map((s) => s.type);
  assert.ok(!types.includes(0x03), 'restriction-site cache dropped after edit');
  const back = readDna(bytes, 'x.dna');
  assert.equal(back.seq.slice(0, 7), 'GATTACA');
  assert.equal(back.features.length, featuresBefore);
});
