import test from 'node:test';
import assert from 'node:assert/strict';
import { readGenBank, writeGenBank, parseLocation } from '../js/formats/genbank.js';

// A small record in the layout other programs write: no division, a molecule type
// of ds-DNA, colours as /ApEinfo qualifiers, and the non-standard keys "feature" and "primer".
const SEQ = 'acgt'.repeat(150);
const seqBlock = (s) => Array.from({ length: Math.ceil(s.length / 60) }, (_, i) => `${String(i * 60 + 1).padStart(9)} ${s.slice(i * 60, i * 60 + 60).match(/.{1,10}/g).join(' ')}`).join('\n');
const record = (extra = '') => `LOCUS       pTest                    600 bp ds-DNA     circular     18-NOV-2026
DEFINITION  .
ACCESSION   U00001
VERSION     U00001.1
KEYWORDS    .
SOURCE      Escherichia coli
  ORGANISM  Escherichia coli
            Bacteria; Pseudomonadota.
REFERENCE   1  (bases 1 to 600)
  TITLE     Direct Submission
FEATURES             Location/Qualifiers
     source          1..600
                     /organism="Escherichia coli"
                     /db_xref="taxon:562"
     CDS             101..400
                     /label="a name that is rather long for a feature, well past forty characters"
                     /ApEinfo_fwdcolor="#84b0dc"
                     /ApEinfo_revcolor="#84b0dc"
     feature         10..20
                     /label=ROP
     primer          complement(30..52)
                     /label=pGEX_3_primer
${extra}ORIGIN
${seqBlock(SEQ)}
//
`;

test('names are kept whole and the colours a file chose are used', () => {
  const [d] = readGenBank(record());
  const cds = d.features.find((f) => f.type === 'CDS');
  assert.equal(cds.name, 'a name that is rather long for a feature, well past forty characters');
  assert.equal(cds.color, '#84b0dc');
  assert.equal(cds.segments[0].color, '#84b0dc');
});

test('non-standard feature keys become the standard ones', () => {
  const [d] = readGenBank(record());
  assert.deepEqual(d.features.map((f) => f.type), ['CDS', 'misc_feature', 'primer_bind']);
});

test('the header is kept: molecule type, accession, version, references and the source qualifiers', () => {
  const [d] = readGenBank(record());
  const out = writeGenBank(d);
  assert.match(out, /^LOCUS\s+pTest\s+600 bp\s+ds-DNA\s+circular/);
  assert.match(out, /^ACCESSION   U00001$/m);
  assert.match(out, /^VERSION     U00001\.1$/m);
  assert.match(out, /Bacteria; Pseudomonadota\./);
  assert.match(out, /\/db_xref="taxon:562"/);
  assert.doesNotMatch(out, /BublaGene/);
});

test('a reference to "bases 1 to N" follows the sequence when it changes length', () => {
  const [d] = readGenBank(record());
  d.seq += 'ACGT';
  assert.match(writeGenBank(d), /\(bases 1 to 604\)/);
});

test('renaming a feature and changing its colour reach the saved file', () => {
  const [d] = readGenBank(record());
  const cds = d.features.find((f) => f.type === 'CDS');
  cds.name = 'My gene';
  cds.color = '#112233';
  const [back] = readGenBank(writeGenBank(d));
  const again = back.features.find((f) => f.type === 'CDS');
  assert.equal(again.name, 'My gene');
  assert.equal(again.color, '#112233');
});

test('a feature that was only ever given a default colour does not gain colour qualifiers', () => {
  const [d] = readGenBank(record());
  const out = writeGenBank(d);
  const rop = out.split('\n').findIndex((l) => /misc_feature\s+10\.\.20/.test(l));
  assert.doesNotMatch(out.split('\n').slice(rop, rop + 3).join('\n'), /ApEinfo/);
});

test('a byte-order mark or Windows line endings do not hide the LOCUS line', () => {
  const [d] = readGenBank(`﻿${record().replace(/\n/g, '\r\n')}`, 'x.gb');
  assert.equal(d.name, 'pTest');
  assert.equal(d.circular, true);
  assert.equal(d.features.length, 3);
});

test('a qualifier line that begins with a slash inside quotes continues the value', () => {
  const [d] = readGenBank(record(`     misc_feature    50..60
                     /note="alpha and
                     /or beta"
                     /codon_start=1
`));
  const f = d.features.find((x) => x.qualifiers.some((q) => q.name === 'note'));
  assert.equal(f.qualifiers.find((q) => q.name === 'note').values[0].value, 'alpha and /or beta');
  assert.ok(f.qualifiers.some((q) => q.name === 'codon_start'), 'the next qualifier is still its own');
  const [back] = readGenBank(writeGenBank(d));
  assert.equal(back.features.find((x) => x.qualifiers.some((q) => q.name === 'note')).qualifiers.find((q) => q.name === 'note').values[0].value, 'alpha and /or beta');
});

test('a flag stays a flag, and a long unbroken value survives a round trip', () => {
  const word = 'aVeryLongUnbrokenValue'.repeat(5);
  const [d] = readGenBank(record(`     misc_feature    70..80
                     /pseudo
                     /note="${word}"
`));
  const out = writeGenBank(d);
  assert.match(out, /^ {21}\/pseudo$/m);
  const [back] = readGenBank(out);
  const f = back.features.find((x) => x.qualifiers.some((q) => q.name === 'pseudo'));
  assert.equal(f.qualifiers.find((q) => q.name === 'note').values[0].value, word);
});

test('a quoted number stays text, an unquoted one stays a number', () => {
  const [d] = readGenBank(record(`     misc_feature    70..80
                     /codon_start=2
                     /note="12"
`));
  const f = d.features[d.features.length - 1];
  assert.equal(f.qualifiers.find((q) => q.name === 'codon_start').values[0].kind, 'int');
  assert.equal(f.qualifiers.find((q) => q.name === 'note').values[0].kind, 'text');
  assert.match(writeGenBank(d), /\/note="12"/);
});

test('uncertain positions such as (3.5)..8 are read, not dropped', () => {
  assert.deepEqual(parseLocation('(3.5)..8', 100).segments, [{ start: 2, end: 8 }]);
  assert.deepEqual(parseLocation('2..(8.10)', 100).segments, [{ start: 1, end: 10 }]);
  const [d] = readGenBank(record(`     misc_feature    (3.5)..8
                     /label=odd
`));
  assert.ok(d.features.some((f) => f.name === 'odd'));
});

test('a feature whose location cannot be read is reported, not silently lost', () => {
  const [d] = readGenBank(record(`     misc_feature    banana
                     /label=bad
`));
  assert.equal(d.source.warnings.length, 1);
  assert.match(d.source.warnings[0], /banana/);
});

test('a long locus name is not cut, and a long definition wraps inside 80 columns', () => {
  const [d] = readGenBank(record());
  d.name = 'a_rather_long_construct_name_for_a_plasmid';
  d.notes.description = 'word '.repeat(40).trim();
  const out = writeGenBank(d);
  assert.match(out, /^LOCUS\s+a_rather_long_construct_name_for_a_plasmid\s/);
  const lines = out.split('\n');
  const definition = lines.slice(lines.findIndex((l) => l.startsWith('DEFINITION')), lines.findIndex((l) => l.startsWith('ACCESSION')));
  assert.ok(definition.length > 1 && definition.every((l) => l.length <= 80), 'wrapped onto several lines, none over 80');
  assert.equal(readGenBank(out)[0].name, 'a_rather_long_construct_name_for_a_plasmid');
});

test('a long join wraps across lines and reads back the same', () => {
  const parts = Array.from({ length: 12 }, (_, i) => `${i * 40 + 1}..${i * 40 + 20}`).join(',');
  const [d] = readGenBank(record(`     misc_feature    join(${parts})
                     /label=many
`));
  const out = writeGenBank(d);
  assert.ok(out.split('\n').every((l) => l.length <= 80 || /^ {21}\//.test(l)));
  const [back] = readGenBank(out);
  assert.deepEqual(back.features.find((f) => f.name === 'many').segments.map((s) => [s.start, s.end]), d.features.find((f) => f.name === 'many').segments.map((s) => [s.start, s.end]));
});

test('RNA and linear records keep their molecule type and topology', () => {
  const text = `LOCUS       r1                        60 bp    RNA     linear   SYN 01-JAN-2020\nDEFINITION  an rna.\nORIGIN\n${seqBlock('acgu'.repeat(15))}\n//\n`;
  const [d] = readGenBank(text);
  const out = writeGenBank(d);
  assert.match(out, /^LOCUS\s+r1\s+60 bp\s+RNA\s+linear/);
  assert.match(out, /\/mol_type="other RNA"/);
});
