// CLI parity: every engine is reachable from `bg`, emits JSON when asked, and
// the same seed reproduces the same report.

import test from 'node:test';
import assert from 'node:assert/strict';
import { spawnSync } from 'node:child_process';
import { mkdtempSync, writeFileSync, readFileSync, existsSync, rmSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';

const ROOT = new URL('..', import.meta.url).pathname;
const DIR = mkdtempSync(join(tmpdir(), 'bublagene-cli-'));

function randomSeq(n, seed) {
  let x = seed >>> 0;
  let out = '';
  for (let i = 0; i < n; i++) {
    x = (1103515245 * x + 12345) >>> 0;
    out += 'ACGT'[(x >>> 16) & 3];
  }
  return out;
}

function writeFasta(name, seq) {
  const path = join(DIR, `${name}.fa`);
  writeFileSync(path, `>${name}\n${seq}\n`);
  return path;
}

const run = (args) => spawnSync(process.execPath, ['tools/bg.mjs', ...args], { cwd: ROOT, encoding: 'utf8' });
const json = (args) => {
  const r = run(args);
  assert.equal(r.status, 0, `${args.join(' ')} failed: ${r.stderr}`);
  return JSON.parse(r.stdout);
};

const guide = 'GACCTGACCTGACCTGACCT';
const template = 'TTT' + guide + 'TGG' + 'A'.repeat(80);
const templatePath = writeFasta('template', template);
const randomTemplate = randomSeq(320, 5);
const designPath = writeFasta('design', randomTemplate);
const aBody = randomSeq(150, 1);
const bBody = randomSeq(120, 2);
const a = aBody;
const b = a.slice(-20) + bBody;
const aPath = writeFasta('a', a);
const bPath = writeFasta('b', b);
const g1seq = 'GGTCTC' + 'A' + 'CGTT' + 'ATGCATGCAT' + 'GGGGG' + 'GAGACC' + 'TTTT';
const g2seq = 'GGTCTC' + 'A' + 'GGGG' + 'CCCCAAAA' + 'CGTT' + 'A' + 'GAGACC' + 'TTTT';
const g1Path = writeFasta('g1', g1seq);
const g2Path = writeFasta('g2', g2seq);

test('info, design, pcr, digest and composition emit JSON', () => {
  const info = json(['info', templatePath, '--json']);
  assert.equal(info.schema, 'bublagene/report/1');
  assert.equal(info.length, template.length);
  const design = json(['design', designPath, '1..' + randomTemplate.length, '--json']);
  assert.ok(design.forward.sequence.length >= 18);
  const pcr = json(['pcr', designPath, design.forward.sequence, design.reverse.sequence, '--json']);
  assert.ok(pcr.product.length > 0);
  const digest = json(['digest', templatePath, 'EcoRI,BsaI', '--json']);
  assert.equal(digest.enzymes.length, 2);
  const composition = json(['composition', templatePath, '--json']);
  assert.equal(composition.length, template.length);
});

test('align, variants, contig, crispr and fold emit JSON', () => {
  const align = json(['align', aPath, bPath, '--json', '--band', '30']);
  assert.ok(align.score >= 0);
  assert.ok(align.identity >= 0);
  const variants = json(['variants', aPath, bPath, '--json']);
  assert.ok(Array.isArray(variants.variants));
  const contig = json(['contig', aPath, bPath, '--json', '--min-overlap', '15']);
  assert.equal(contig.length, 1);
  const crispr = json(['crispr', templatePath, '--json', '--offtargets', '1']);
  assert.ok(crispr.length >= 1);
  assert.equal(crispr[0].sequence, guide);
  const folded = json(['fold', writeFasta('foldseq', 'GGGGAAAACCCC'), '--json']);
  assert.equal(folded.structure, '((((....))))');
});

test('assembly, golden gate and exports emit JSON', () => {
  const assembled = json(['assemble', aPath, bPath, '--linear', '--json']);
  assert.equal(assembled.schema, 'bublagene/assembly/1');
  assert.equal(assembled.product.seqLen, a.length + b.length - 20);
  const gg = json(['golden-gate', g1Path, g2Path, '--enzyme', 'BsaI', '--json']);
  assert.equal(gg.method, 'golden-gate');
  assert.equal(gg.product.seqLen, 26);

  const svg = join(DIR, 'out.svg');
  const png = join(DIR, 'out.png');
  const gff = join(DIR, 'out.gff3');
  const pdf = join(DIR, 'out.pdf');
  for (const out of [svg, png, gff, pdf]) {
    const r = run(['export', templatePath, out, '--json']);
    assert.equal(r.status, 0, r.stderr);
    const meta = JSON.parse(r.stdout);
    assert.ok(meta.bytes > 20);
    assert.ok(existsSync(out));
  }
  assert.ok(readFileSync(svg, 'utf8').startsWith('<svg'));
  assert.deepEqual([...readFileSync(png).subarray(0, 8)], [137, 80, 78, 71, 13, 10, 26, 10]);
});

test('report, pipeline and seeded edits are reproducible', () => {
  const report = json(['report', templatePath, '--orfs', '10', '--enzymes', 'EcoRI,BsaI']);
  assert.equal(report.schema, 'bublagene/report/1');
  assert.ok(Array.isArray(report.orfs));
  assert.ok(Array.isArray(report.enzymes));
  const steps = join(DIR, 'steps.json');
  writeFileSync(steps, JSON.stringify({ seed: 7, steps: [{ op: 'info' }, { op: 'fold' }, { op: 'crispr', limit: 2 }, { op: 'export', out: join(DIR, 'pipe.svg') }] }));
  const pipe = json(['pipeline', templatePath, steps]);
  assert.equal(pipe.schema, 'bublagene/pipeline/1');
  assert.equal(pipe.seed, 7);
  assert.equal(pipe.results.length, 4);
  assert.ok(pipe.results.every((r) => !r.error));

  const one = json(['crispr-edit', templatePath, guide, '--seed', '9', '--json']);
  const two = json(['crispr-edit', templatePath, guide, '--seed', '9', '--json']);
  assert.deepEqual(one, two);
  const other = json(['crispr-edit', templatePath, guide, '--seed', '10', '--json']);
  assert.notDeepEqual(one.change, other.change);
});

test.after(() => { try { rmSync(DIR, { recursive: true, force: true }); } catch { /* ignore */ } });
