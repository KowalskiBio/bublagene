// Bounded fuzzing for the binary and text parsers. Inputs are small and every
// failure is allowed to be an Error; hanging or a non-Error throw is not.

import test from 'node:test';
import assert from 'node:assert/strict';
import fc from 'fast-check';
import { readGenBank, readFasta } from '../js/formats/genbank.js';
import { readGFF3, readBED, readFeatureCSV } from '../js/formats/features.js';
import { readAB1 } from '../js/formats/ab1.js';
import { readSCF } from '../js/formats/scf.js';
import { traceFormat } from '../js/traces.js';

function attempt(fn) {
  try { fn(); return 'ok'; } catch (e) { assert.ok(e instanceof Error, `non-Error throw: ${e}`); return 'error'; }
}

test('random bytes never hang or throw non-Errors in the trace parsers', () => {
  fc.assert(fc.property(fc.uint8Array({ maxLength: 1500 }), (bytes) => {
    const withAbif = new Uint8Array(bytes.length + 4);
    withAbif.set([0x41, 0x42, 0x49, 0x46], 0);
    withAbif.set(bytes, 4);
    const withScf = new Uint8Array(bytes.length + 4);
    withScf.set([0x2e, 0x73, 0x63, 0x66], 0);
    withScf.set(bytes, 4);
    return ['ok', 'error'].includes(attempt(() => readAB1(withAbif))) &&
      ['ok', 'error'].includes(attempt(() => readSCF(withScf))) &&
      ['ok', 'error'].includes(attempt(() => traceFormat(bytes)));
  }), { numRuns: 120 });
});

test('truncated valid traces fail cleanly', () => {
  const base = new Uint8Array(300 + 28 * 4);
  base.set([0x41, 0x42, 0x49, 0x46], 0);
  fc.assert(fc.property(fc.integer({ min: 0, max: base.length - 1 }), (cut) => {
    try { readAB1(base.slice(0, cut)); } catch (e) { assert.ok(e instanceof Error); }
    return true;
  }), { numRuns: 60 });
});

test('random text never hangs the GenBank/FASTA/GFF readers', () => {
  fc.assert(fc.property(fc.string({ maxLength: 800 }), (text) => {
    return ['ok', 'error'].includes(attempt(() => readGenBank(text))) &&
      ['ok', 'error'].includes(attempt(() => readFasta(text))) &&
      ['ok', 'error'].includes(attempt(() => readGFF3(text))) &&
      ['ok', 'error'].includes(attempt(() => readBED(text))) &&
      ['ok', 'error'].includes(attempt(() => readFeatureCSV(text)));
  }), { numRuns: 120 });
});
