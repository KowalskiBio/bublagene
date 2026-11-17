import test from 'node:test';
import assert from 'node:assert/strict';
import { readAB1 } from '../js/formats/ab1.js';
import { readSCF } from '../js/formats/scf.js';
import {
  readTrace, traceFormat, trimTrace, trimFeatures, traceToDocument, sliceTrace,
  baseAtSample, qualitiesOf, qualityStats,
} from '../js/traces.js';

function writePString(str, buf, offset) {
  buf[offset] = str.length;
  for (let i = 0; i < str.length; i++) buf[offset + 1 + i] = str.charCodeAt(i);
}

function makeAB1({ bases = 'ACGTACGTAC', qualities = [5, 5, 5, 30, 30, 30, 30, 30, 5, 5], name = 'sample1' } = {}) {
  const n = bases.length;
  const tags = [];
  const dirOffset = 128;
  const entrySize = 28;
  let dataOffset = dirOffset + entrySize * 5;
  const add = (name, number, type, elemSize, count, dataSize) => {
    const tag = { name, number, type, elemSize, count, dataSize, dataOffset };
    tags.push(tag);
    dataOffset += dataSize;
  };
  const dataSize = n * 4 * 2;
  add('DATA', 9, 3, 2, n * 4, dataSize);
  add('PBAS', 2, 2, 1, n, n);
  add('PCON', 2, 1, 1, n, n);
  add('PLOC', 2, 3, 2, n, n * 2);
  add('SMPL', 1, 2, 1, name.length + 1, name.length + 1);
  const total = dataOffset;
  const buf = new Uint8Array(total);
  const dv = new DataView(buf.buffer);
  buf.set([0x41, 0x42, 0x49, 0x46], 0); // ABIF
  dv.setUint16(4, 101, false);
  dv.setUint32(10, dirOffset, false);
  dv.setUint32(14, tags.length, false);
  for (const [k, t] of tags.entries()) {
    const p = dirOffset + k * entrySize;
    for (let c = 0; c < 4; c++) buf[p + c] = t.name.charCodeAt(c);
    dv.setUint32(p + 4, t.number, false);
    dv.setUint16(p + 8, t.type, false);
    dv.setUint16(p + 10, t.elemSize, false);
    dv.setUint32(p + 12, t.count, false);
    dv.setUint32(p + 16, t.dataSize, false);
    dv.setUint32(p + 20, t.dataOffset, false);
    dv.setUint32(p + 24, 0, false);
  }
  const byKey = Object.fromEntries(tags.map((t) => [`${t.name}${t.number}`, t]));
  const data = byKey.DATA9;
  for (let i = 0; i < n; i++) {
    for (let c = 0; c < 4; c++) dv.setUint16(data.dataOffset + (c * n + i) * 2, 100 + c * 10 + i, false);
  }
  for (let i = 0; i < n; i++) {
    buf[byKey.PBAS2.dataOffset + i] = bases.charCodeAt(i);
    buf[byKey.PCON2.dataOffset + i] = qualities[i];
    dv.setUint16(byKey.PLOC2.dataOffset + i * 2, 10 + i * 5 + 1, false); // 1-based
  }
  writePString(name, buf, byKey.SMPL1.dataOffset);
  return buf;
}

function makeSCF({ bases = 'ACGTAC', samples = 30, name = 'scf1' } = {}) {
  const n = bases.length;
  const sampleOffset = 56;
  const basesOffset = sampleOffset + samples * 4;
  const total = basesOffset + n * 4;
  const buf = new Uint8Array(total);
  const dv = new DataView(buf.buffer);
  buf.set([0x2e, 0x73, 0x63, 0x66], 0); // .scf
  dv.setUint32(4, 3, false);
  dv.setUint32(8, samples, false);
  dv.setUint32(12, sampleOffset, false);
  dv.setUint32(16, n, false);
  dv.setUint32(20, basesOffset, false);
  dv.setUint32(24, 0, false);
  dv.setUint32(28, 0, false);
  dv.setUint32(32, 1, false);
  for (let i = 0; i < samples; i++) {
    for (let c = 0; c < 4; c++) buf[sampleOffset + (i * 4 + c)] = 50 + c * 10 + (i % 7);
  }
  for (let i = 0; i < n; i++) {
    const p = basesOffset + i * 4;
    buf[p] = Math.min(samples - 1, i * 5);
    buf[p + 1] = 40;
    buf[p + 2] = bases.charCodeAt(i);
    buf[p + 3] = 0;
  }
  return buf;
}

test('ABI traces parse channels, bases, positions, qualities and names', () => {
  const buffer = makeAB1();
  assert.equal(traceFormat(buffer), 'AB1');
  const trace = readAB1(buffer, 'x.ab1');
  assert.equal(trace.format, 'AB1');
  assert.equal(trace.name, 'sample1');
  assert.equal(trace.bases, 'ACGTACGTAC');
  assert.equal(trace.sampleCount, 10);
  assert.equal(trace.samples.A.length, 10);
  assert.equal(trace.samples.A[0], 100);
  assert.equal(trace.samples.T[0], 130);
  assert.equal(trace.positions[0], 10);
  assert.equal(trace.positions[1], 15);
  assert.equal(trace.confidence[3], 30);
  assert.ok(trace.tags.includes('DATA9'));
});

test('SCF traces parse channels, bases and Phred from probabilities', () => {
  const buffer = makeSCF();
  assert.equal(traceFormat(buffer), 'SCF');
  const trace = readSCF(buffer, 'x.scf');
  assert.equal(trace.format, 'SCF');
  assert.equal(trace.bases, 'ACGTAC');
  assert.equal(trace.sampleCount, 30);
  assert.equal(trace.samples.C[0], 60);
  assert.equal(trace.positions[1], 5);
  assert.ok(trace.confidence[0] > 0);
  assert.equal(trace.version, 3);
});

test('readTrace dispatches by magic', async () => {
  const ab1 = await readTrace(makeAB1(), 'a.ab1');
  assert.equal(ab1.format, 'AB1');
  const scf = await readTrace(makeSCF(), 's.scf');
  assert.equal(scf.format, 'SCF');
  assert.equal(traceFormat(new Uint8Array([1, 2, 3, 4])), null);
  await assert.rejects(() => readTrace(new Uint8Array([1, 2, 3, 4])), /Unrecognised/);
});

test('Mott and Phred trimming find the high-quality core', () => {
  const trace = readAB1(makeAB1());
  const phred = trimTrace(trace, { mode: 'phred', cutoff: 20 });
  assert.equal(phred.start, 3);
  assert.equal(phred.end, 8);
  assert.equal(phred.length, 5);
  assert.deepEqual(phred.removed, { start: 3, end: 2 });
  const mott = trimTrace(trace, { mode: 'mott', window: 2, cutoff: 20 });
  assert.ok(mott.start >= 3);
  assert.ok(mott.end <= 9);
  assert.ok(mott.length > 0);
  assert.ok(mott.mean >= 20);
  const stats = qualityStats(trace);
  assert.equal(stats.count, 10);
  assert.equal(stats.min, 5);
  assert.equal(stats.max, 30);
  assert.equal(qualitiesOf(trace).length, 10);
});

test('trim features, documents, slices and base navigation', () => {
  const trace = readAB1(makeAB1());
  const doc = traceToDocument(trace, { name: 'read 1', trimmed: true });
  assert.equal(doc.seq, trace.bases);
  assert.equal(doc.features.length, 2);
  assert.deepEqual(doc.features.map((f) => [f.segments[0].start, f.segments[0].end]), [[0, 3], [8, 10]]);
  const features = trimFeatures(doc, { start: 2, end: 9, mode: 'mott', cutoff: 20 });
  assert.equal(features.length, 2);
  const slice = sliceTrace(trace, 2, 6);
  assert.equal(slice.bases, 'GTAC');
  assert.equal(slice.positions.length, 4);
  assert.ok(slice.sampleCount > 0);
  const hit = baseAtSample(trace, 21);
  assert.equal(hit.base, 'G');
  assert.ok(hit.distance <= 3);
});
