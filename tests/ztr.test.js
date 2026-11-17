// There is no second implementation to compare a trace against, so the decode
// is judged on whether the result behaves like a chromatogram: smooth curves,
// and the base that was called standing tallest where it was called.
import test from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync, existsSync, readdirSync, statSync } from 'node:fs';
import { join } from 'node:path';
import { readTraces } from '../js/formats/ztr.js';
import { readDna } from '../js/formats/dna.js';

const ROOT = process.env.DNA_SAMPLES ?? '';

function walk(dir, out = []) {
  if (!existsSync(dir)) return out;
  for (const name of readdirSync(dir)) {
    const p = join(dir, name);
    statSync(p).isDirectory() ? walk(p, out) : name.toLowerCase().endsWith('.dna') && out.push(p);
  }
  return out;
}

const withTraces = [];
for (const f of walk(ROOT)) {
  let doc;
  try { doc = readDna(readFileSync(f), f.split('/').pop()); } catch { continue; }
  if ((doc.source?.segments ?? []).some((s) => s.type === 0x12)) withTraces.push({ file: f, doc });
}

const roughness = (a) => {
  let sum = 0;
  for (let i = 1; i < a.length; i++) sum += Math.abs(a[i] - a[i - 1]);
  return sum / Math.max(1, a.length - 1);
};

test('alignable-sequence segments are not mistaken for traces', async () => {
  const doc = { source: { segments: [{ type: 0x10, data: new Uint8Array([0, 0, 0, 0, 1, 2, 3]) }] } };
  assert.deepEqual(await readTraces(doc), []);
});

test('sample files contain a trace to read', { skip: withTraces.length ? false : `no .dna with traces under ${ROOT}` }, () => {
  assert.ok(withTraces.length > 0);
});

test('every trace decodes into a chromatogram', { skip: !withTraces.length }, async () => {
  let checkedTraces = 0;
  for (const { file, doc } of withTraces) {
    const traces = await readTraces(doc);
    assert.ok(traces.length, `${file}: no traces`);
    for (const t of traces) {
      assert.equal(t.error, undefined, `${file} trace ${t.id}: ${t.error}`);
      assert.ok(t.sampleCount > 100, `${file}: too few samples`);
      assert.ok(t.bases.length > 50, `${file}: too few bases`);
      assert.equal(t.positions.length, t.bases.length, `${file}: a position per called base`);
      assert.match(t.bases, /^[ACGTN]+$/, `${file}: bases should be nucleotides`);

      // real traces rise and fall smoothly; a bad decode looks like noise
      const rough = 'ACGT'.split('').reduce((s, c) => s + roughness(t.samples[c]), 0) / 4;
      assert.ok(rough < 200, `${file} trace ${t.id}: traces are not smooth (${Math.round(rough)})`);

      // and the base that was called should be the tallest channel there
      let agree = 0;
      let checked = 0;
      for (let k = 0; k < t.positions.length; k++) {
        const p = t.positions[k];
        if (p >= t.sampleCount) continue;
        let best = 'A';
        for (const c of 'CGT') if (t.samples[c][p] > t.samples[best][p]) best = c;
        if (best === t.bases[k]) agree++;
        checked++;
      }
      const pct = (100 * agree) / checked;
      assert.ok(pct > 95, `${file} trace ${t.id}: called base is the tallest channel only ${pct.toFixed(1)}% of the time`);
      checkedTraces++;
    }
  }
  assert.ok(checkedTraces > 0);
});

test('positions climb and stay inside the trace', { skip: !withTraces.length }, async () => {
  for (const { file, doc } of withTraces) {
    for (const t of await readTraces(doc)) {
      let last = -1;
      for (const p of t.positions) {
        assert.ok(p > last, `${file} trace ${t.id}: positions should increase`);
        assert.ok(p < t.sampleCount, `${file} trace ${t.id}: position past the end of the trace`);
        last = p;
      }
    }
  }
});
