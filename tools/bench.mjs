#!/usr/bin/env node
// Performance budget for a 5 Mb genome: GC, translation, pattern search, sites.
import { performance } from 'node:perf_hooks';
import { gcContent, translate, findPattern, revcomp } from '../js/seq.js';
import { ENZYME_BY_NAME, findSites, digest } from '../js/enzymes.js';

const BUDGET_MS = Number(process.env.BENCH_BUDGET_MS ?? 8000);
let x = 12345;
let seq = '';
const bases = 'ACGT';
while (seq.length < 5_000_000) {
  x = (1103515245 * x + 12345) >>> 0;
  seq += bases[(x >>> 16) & 3];
}
const timings = {};
const time = (name, fn) => {
  const t = performance.now();
  const value = fn();
  timings[name] = performance.now() - t;
  return value;
};
const gc = time('gc', () => gcContent(seq));
const protein = time('translate-300kb', () => translate(seq.slice(0, 300000)));
const hits = time('pattern', () => findPattern(seq, 'GAATTC', false));
const cuts = time('sites', () => ['EcoRI', 'BamHI', 'HindIII', 'SmaI', 'BsaI'].flatMap((n) => findSites(seq, true, ENZYME_BY_NAME[n])));
const frags = time('digest', () => digest(seq, true, cuts));
time('revcomp', () => revcomp(seq.slice(0, 100000)));
const total = Object.values(timings).reduce((a, b) => a + b, 0);
console.log(`5 Mb genome benchmark (${BUDGET_MS} ms budget)`);
for (const [name, ms] of Object.entries(timings)) console.log(`  ${name.padEnd(16)} ${ms.toFixed(0)} ms`);
console.log(`  total           ${total.toFixed(0)} ms  (GC ${gc.toFixed(1)}%, ${protein.length} aa, ${hits.length} hits, ${frags.length} fragments)`);
if (total > BUDGET_MS) {
  console.error(`over budget by ${(total - BUDGET_MS).toFixed(0)} ms`);
  process.exit(1);
}
