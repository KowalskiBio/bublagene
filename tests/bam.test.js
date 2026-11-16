// The BAM reader is checked against samtools' own output for the same file,
// taken from the alignment inside a real .dna document.
import test from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync, existsSync, readdirSync, statSync } from 'node:fs';
import { join } from 'node:path';
import { execFileSync } from 'node:child_process';
import { readBam, alignToReference, referenceSpan } from '../js/formats/bam.js';
import { readDna } from '../js/formats/dna.js';

const ROOT = process.env.DNA_SAMPLES ?? '';
const SAMTOOLS = process.env.SAMTOOLS;

function walk(dir, out = []) {
  if (!existsSync(dir)) return out;
  for (const name of readdirSync(dir)) {
    const p = join(dir, name);
    statSync(p).isDirectory() ? walk(p, out) : name.toLowerCase().endsWith('.dna') && out.push(p);
  }
  return out;
}

// every alignment segment found in the sample files
const alignments = [];
for (const f of walk(ROOT)) {
  let doc;
  try { doc = readDna(readFileSync(f), f.split('/').pop()); } catch { continue; }
  for (const seg of doc.source?.segments ?? []) {
    if (seg.type === 0x1b) alignments.push({ file: f, data: seg.data });
  }
}

test('sample files contain an alignment to read', { skip: alignments.length ? false : `no .dna with aligned reads under ${ROOT}` }, () => {
  assert.ok(alignments.length > 0);
});

test('reads every alignment in the sample files', { skip: !alignments.length }, async () => {
  for (const { file, data } of alignments) {
    const bam = await readBam(data);
    assert.ok(bam.refs.length > 0, `${file}: no reference listed`);
    assert.ok(bam.reads.length > 0, `${file}: no reads`);
    for (const r of bam.reads) {
      assert.equal(r.refLength, referenceSpan(r.cigar));
      const { chars } = alignToReference(r);
      assert.equal(chars.length, r.refLength, `${file} ${r.name}: aligned length should match the CIGAR`);
      assert.ok(r.pos >= 0);
    }
  }
});

// when samtools is available, compare field by field
const haveSamtools = existsSync(SAMTOOLS);

test('matches samtools field for field', { skip: !alignments.length || !haveSamtools ? 'samtools not available' : false }, async () => {
  const { writeFileSync, mkdtempSync } = await import('node:fs');
  const { tmpdir } = await import('node:os');
  const dir = mkdtempSync(join(tmpdir(), 'bublagene-bam-'));
  let compared = 0;
  for (const { data } of alignments) {
    const path = join(dir, 'a.bam');
    writeFileSync(path, data);
    const sam = execFileSync(SAMTOOLS, ['view', path], { encoding: 'utf8' }).trim();
    const mine = await readBam(data);
    const lines = sam ? sam.split('\n') : [];
    assert.equal(mine.reads.length, lines.length, 'read count');
    for (let i = 0; i < lines.length; i++) {
      const [name, flag, ref, pos, mapq, cigar, , , , seq] = lines[i].split('\t');
      const r = mine.reads[i];
      assert.equal(r.name, name, 'read name');
      assert.equal(r.flag, Number(flag), `${name}: flag`);
      assert.equal(r.reference, ref === '*' ? null : ref, `${name}: reference`);
      assert.equal(r.pos + 1, Number(pos), `${name}: position`);   // BAM counts from 0, SAM from 1
      assert.equal(r.mapq, Number(mapq), `${name}: mapping quality`);
      assert.equal(r.cigar.map(([n, op]) => `${n}${op}`).join(''), cigar, `${name}: CIGAR`);
      assert.equal(r.seq, seq, `${name}: sequence`);
      compared++;
    }
  }
  assert.ok(compared > 0, 'nothing was compared');
});
