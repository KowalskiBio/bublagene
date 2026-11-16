import test from 'node:test';
import assert from 'node:assert/strict';
import { CODON_TABLES, TABLE_IDS, validateTable, validateAllTables, translateWithTable, getCodonTable } from '../js/codon-tables.js';
import { translate, sixFrameTranslation, nucleotideCount, dnaMolecularWeight, gcContent } from '../js/seq.js';
import {
  newDocument, newId, FEATURE_TYPES, FEATURE_CATALOG, featureIcon, colorFor, defaultQualifiers,
  codingSequence, featureSequence, frameOffset, phaseOffset, translateFeature, storedTranslation,
  recomputeFeature, recomputeFeatures, proteinMolecularWeight, featureProteinMass, featureGcContent,
  qualifierText, setQualifier, normalizeSegments, qualifierInt, featureTableId,
} from '../js/model.js';
import { parseLocation, parseLocationDetailed, formatLocation, readGenBank, writeGenBank, parseQualifiers } from '../js/formats/genbank.js';

const feat = (type, start, end, strand = 1, qualifiers = []) => ({
  id: newId(), name: type, type, strand, color: '#888',
  segments: [{ start, end }], qualifiers, extra: {},
});

test('every genetic code table is structurally valid and its anchors hold', () => {
  assert.deepEqual(TABLE_IDS, [1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 16, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 33]);
  for (const id of TABLE_IDS) {
    const v = validateTable(id);
    assert.ok(v.ok, `table ${id}: ${v.issues.join(', ')}`);
    assert.equal(CODON_TABLES[id].aas.length, 64);
  }
  assert.equal(validateAllTables().filter((v) => !v.ok).length, 0);
});

test('alternative codes translate their diagnostic codons', () => {
  assert.equal(translateWithTable('TGA', 1), '*');
  assert.equal(translateWithTable('TGA', 2), 'W');
  assert.equal(translateWithTable('AGA', 2), '*');
  assert.equal(translateWithTable('TAA', 6), 'Q');
  assert.equal(translateWithTable('TGA', 10), 'C');
  assert.equal(translateWithTable('CTG', 12), 'S');
  assert.equal(translateWithTable('TGA', 25), 'G');
  assert.equal(translateWithTable('AGA', 24), 'S');
  assert.equal(translateWithTable('AGG', 24), 'K');
  assert.equal(translateWithTable('TTA', 23), '*');
  assert.equal(translateWithTable('TCA', 22), '*');
  assert.equal(translateWithTable('TTG', 11, { initiator: true }), 'M');
  assert.equal(translateWithTable('TTG', 11), 'L');
  assert.equal(translate('ATGTAA', 2), 'M*');
});

test('codon_start / phase and partial codons', () => {
  const doc = newDocument({ seq: 'AATGGCCTAAG' });
  const f = feat('CDS', 0, 11, 1, [{ name: 'codon_start', values: [{ kind: 'int', value: '2' }] }]);
  assert.equal(frameOffset(f), 1);
  assert.equal(phaseOffset(f), 1);
  assert.equal(codingSequence(doc, f), 'AATGGCCTAAG'.slice(1));
  assert.equal(translateFeature(doc, f), 'MA*');
  assert.equal(storedTranslation(doc, f), 'MA');
  const partial = feat('CDS', 1, 6);
  assert.equal(translateFeature(doc, partial), 'M'); // dangling two bases dropped
  const phase = feat('CDS', 0, 6, 1, [{ name: 'codon_start', values: [{ kind: 'int', value: '3' }] }]);
  assert.equal(frameOffset(phase), 2);
  assert.equal(phaseOffset(phase), 2);
});

test('per-feature transl_table overrides the document default', () => {
  const doc = newDocument({ seq: 'ATGTAATGA' });
  doc.transl_table = 1;
  const f = feat('CDS', 0, 9, 1, [{ name: 'transl_table', values: [{ kind: 'int', value: '6' }] }]);
  assert.equal(featureTableId(f, doc), 6);
  assert.equal(translateFeature(doc, f), 'MQ*');
  const plain = feat('CDS', 0, 9);
  assert.equal(featureTableId(plain, doc), 1);
});

test('fuzzy ends, order(), nested locations and wrapping joins', () => {
  assert.deepEqual(parseLocation('<1..>50', 100), { strand: 1, segments: [{ start: 0, end: 50, partialStart: true, partialEnd: true }] });
  assert.deepEqual(parseLocation('order(5..8,1..2)', 100), { strand: 1, segments: [{ start: 4, end: 8 }, { start: 0, end: 2 }] });
  assert.equal(parseLocationDetailed('order(5..8,1..2)', 100).locationType, 'order');
  assert.equal(parseLocationDetailed('complement(join(90..100,1..5))', 100).locationType, 'join');
  assert.deepEqual(parseLocation('join(90..100,1..5)', 100), { strand: 1, segments: [{ start: 89, end: 5 }] });
  assert.deepEqual(parseLocation('100^101', 200), { strand: 1, segments: [{ start: 99, end: 100, between: true }] });
  const nested = parseLocation('complement(join(10..20,30..40))', 100);
  assert.deepEqual(nested, { strand: -1, segments: [{ start: 9, end: 20 }, { start: 29, end: 40 }] });
  const f = { strand: 1, segments: [{ start: 0, end: 50, partialStart: true, partialEnd: true }] };
  assert.equal(formatLocation(f, 100), '<1..>50');
  const ordered = { locationType: 'order', strand: 1, segments: [{ start: 4, end: 8 }, { start: 0, end: 2 }] };
  assert.equal(formatLocation(ordered, 100), 'order(5..8,1..2)');
  const wrap = { strand: -1, segments: [{ start: 89, end: 5 }] };
  assert.equal(formatLocation(wrap, 100), 'complement(join(90..100,1..5))');
});

test('multi-segment features: gaps, copy, edit and junction preservation', () => {
  const doc = newDocument({ seq: 'AAAA' + 'ATGGCC' + 'TTTT' + 'TAA' + 'CCCC' });
  doc.features = [{ id: newId(), name: 'spliced', type: 'CDS', strand: 1, color: '#888', qualifiers: [], extra: {},
    segments: [{ start: 4, end: 10 }, { start: 14, end: 17 }] }];
  assert.equal(featureSequence(doc, doc.features[0]), 'ATGGCC' + 'TAA');
  assert.equal(codingSequence(doc, doc.features[0]), 'ATGGCC' + 'TAA');
  const rc = { ...doc.features[0], strand: -1 };
  assert.equal(featureSequence({ ...doc, features: [rc] }, rc), 'TTA' + 'GGCCAT');
  assert.deepEqual(normalizeSegments([{ start: 1, end: 3 }, { start: 3, end: 5 }]).map((s) => [s.start, s.end]), [[1, 5]]);
  assert.equal(normalizeSegments([{ start: 1, end: 3, partialEnd: true }, { start: 3, end: 5 }]).length, 2);
});

test('qualifiers: order, multi-value keys and unknown keys are preserved', () => {
  const q = parseQualifiers(['gene="abc"', 'note="one"', 'note="two"', 'custom_thing="x y"', 'codon_start=2']);
  assert.deepEqual(q.map((x) => x.name), ['gene', 'note', 'custom_thing', 'codon_start']);
  assert.equal(q[1].values.length, 2);
  assert.deepEqual(q[1].values.map((v) => v.value), ['one', 'two']);
  assert.equal(q[3].values[0].kind, 'int');
  const f = feat('CDS', 0, 9, 1, q);
  assert.equal(qualifierInt(f, 'codon_start'), 2);
  assert.equal(qualifierText(f, 'note'), 'one; two');
  setQualifier(f, 'product', 'protein X');
  assert.equal(f.qualifiers.at(-1).name, 'product');
  assert.equal(f.qualifiers.filter((x) => x.name === 'note').length, 1);
});

test('GenBank parse -> write -> parse is semantically stable', () => {
  const text = `LOCUS       pRound              60 bp    DNA     circular SYN 01-JAN-2024
DEFINITION  round trip.
FEATURES             Location/Qualifiers
     source          1..60
     CDS             complement(join(<10..20,30..>40))
                     /gene="g1"
                     /note="first"
                     /note="second"
                     /transl_table=11
                     /codon_start=2
                     /custom_key="kept"
     misc_feature    order(1..2,5..6)
ORIGIN
        1 atgcatgcat gcatgcatgc atgcatgcat gcatgcatgc atgcatgcat gcatgcatgc
//
`;
  const [a] = readGenBank(text);
  const [b] = readGenBank(writeGenBank(a));
  const strip = (d) => d.features.map((f) => ({
    type: f.type, strand: f.strand, locationType: f.locationType ?? 'join',
    segments: f.segments.map((s) => ({ start: s.start, end: s.end, partialStart: !!s.partialStart, partialEnd: !!s.partialEnd })),
    qualifiers: f.qualifiers,
  }));
  assert.deepEqual(strip(b), strip(a));
  assert.equal(a.features[0].qualifiers.find((q) => q.name === 'note').values.length, 2);
  assert.equal(a.features[0].qualifiers.find((q) => q.name === 'custom_key').values[0].value, 'kept');
  assert.equal(b.features[0].segments[0].partialStart, true);
  assert.equal(b.features[0].segments[1].partialEnd, true);
});

test('translation, molecular weight and composition recompute deterministically', () => {
  const doc = newDocument({ seq: 'ATGGCCTAA' });
  const f = feat('CDS', 0, 9);
  doc.features = [f];
  recomputeFeature(doc, f);
  assert.equal(qualifierText(f, 'translation'), 'MA');
  doc.seq = 'ATGGCCTTTTAA';
  f.segments[0].end = 12;
  recomputeFeature(doc, f);
  assert.equal(qualifierText(f, 'translation'), 'MAF');
  assert.ok(proteinMolecularWeight('MA') > 200);
  assert.ok(featureProteinMass(doc, f) > 300);
  assert.equal(featureGcContent(doc, f), gcContent('ATGGCCTTTTAA'));
  const counts = nucleotideCount('AACCGGTTN');
  assert.equal(counts.A, 2);
  assert.equal(counts.N, 1);
  assert.equal(counts.length, 9);
  assert.ok(dnaMolecularWeight('ATGC') > 1200);
  recomputeFeatures(doc);
  assert.equal(qualifierText(f, 'translation'), 'MAF');
});

test('feature catalog covers the annotation keys with defaults and icons', () => {
  for (const key of ['misc_feature', 'regulatory', '-35_signal', '-10_signal', 'RBS', 'terminator', 'ncRNA', 'tRNA', 'rRNA', 'primer_bind', 'protein_bind', 'repeat_region', 'mobile_element']) {
    assert.ok(FEATURE_TYPES.includes(key), `missing ${key}`);
    assert.ok(FEATURE_CATALOG[key], `missing catalog entry ${key}`);
    assert.ok(featureIcon(key).length > 5);
    assert.match(colorFor(key), /^#[0-9a-f]{6}$/i);
  }
  assert.equal(defaultQualifiers('CDS')[0].name, 'transl_table');
  assert.deepEqual(defaultQualifiers('misc_feature'), []);
  assert.equal(featureIcon('unknown_key'), FEATURE_CATALOG.misc_feature.icon);
});

test('six-frame translation reports all frames', () => {
  const frames = sixFrameTranslation('ATGAAATAA');
  assert.equal(frames.length, 6);
  assert.equal(frames[0].protein, 'MK*');
  assert.deepEqual(frames.filter((f) => f.strand === -1).map((f) => f.frame), [1, 2, 3]);
  assert.equal(frames[3].protein, translate('TTATTTCAT'));
});
