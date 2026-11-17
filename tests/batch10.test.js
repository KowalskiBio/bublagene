import test from 'node:test';
import assert from 'node:assert/strict';
import { writeGFF3, readGFF3, writeBED, readBED, writeFeatureCSV, readFeatureCSV, writePlainMap } from '../js/formats/features.js';
import { exportSVG, exportPNG, exportPDF, exportDiagram } from '../js/export.js';
import {
  projectRecord, projectToJSON, projectFromJSON, projectDigest, migrateRecord,
  CURRENT_SCHEMAS, PROJECT_SCHEMA, ALIGNMENT_SCHEMA, alignmentRecord,
} from '../js/project.js';
import { newDocument, newId } from '../js/model.js';

function sampleDoc() {
  const doc = newDocument({ name: 'pTest', seq: 'ATGC'.repeat(25), circular: true });
  doc.features = [
    { id: 'f1', name: 'cds one', type: 'CDS', strand: 1, color: '#e0884c', segments: [{ start: 4, end: 30 }],
      qualifiers: [{ name: 'gene', values: [{ kind: 'text', value: 'xylA' }] }, { name: 'note', values: [{ kind: 'text', value: 'two words' }] }], extra: {} },
    { id: 'f2', name: 'spliced', type: 'misc_feature', strand: -1, color: '#888', segments: [{ start: 60, end: 70 }, { start: 20, end: 25 }], qualifiers: [], extra: {} },
  ];
  doc.primers = [{ id: 'p1', name: 'fw', seq: 'ATGCATGC', description: '', extra: {} }];
  doc.notes.description = 'a project test';
  return doc;
}

test('GFF3 and BED exports are valid and semantically round-trip', () => {
  const doc = sampleDoc();
  const gff = writeGFF3(doc);
  assert.match(gff, /^##gff-version 3/);
  assert.equal(gff.split('\n').filter((l) => l && !l.startsWith('#')).length, 3);
  const back = readGFF3(gff);
  assert.equal(back.features.length, 2);
  const cds = back.features.find((f) => f.type === 'CDS');
  assert.equal(cds.name, 'cds one');
  assert.deepEqual(cds.segments.map((s) => [s.start, s.end]), [[4, 30]]);
  assert.equal(cds.qualifiers.find((q) => q.name === 'gene').values[0].value, 'xylA');
  const spliced = back.features.find((f) => f.name === 'spliced');
  assert.equal(spliced.strand, -1);
  assert.equal(spliced.segments.length, 2);

  const bed = writeBED(doc);
  assert.match(bed, /^track/);
  const bedBack = readBED(bed);
  assert.equal(bedBack.features.length, 2);
  assert.deepEqual(bedBack.features.find((f) => f.name === 'spliced').segments.map((s) => [s.start, s.end]).sort(), [[20, 25], [60, 70]]);
});

test('CSV feature tables and plain-text maps', () => {
  const doc = sampleDoc();
  const csv = writeFeatureCSV(doc);
  assert.match(csv, /^name,type,strand,start,end,segments,qualifiers/);
  const back = readFeatureCSV(csv);
  assert.equal(back.features.length, 2);
  assert.equal(back.features[0].name, 'cds one');
  assert.equal(back.features[0].qualifiers.find((q) => q.name === 'gene').values[0].value, 'xylA');
  assert.ok(back.features[0].qualifiers.find((q) => q.name === 'note').values[0].value.includes('two words'));
  const map = writePlainMap(doc);
  assert.ok(map.includes('pTest'));
  assert.ok(map.includes('cds one'));
  assert.ok(map.includes('spliced'));
});

test('SVG, PNG and PDF diagrams have valid headers and content', () => {
  const doc = sampleDoc();
  const svg = exportSVG(doc);
  assert.ok(svg.startsWith('<svg'));
  assert.ok(svg.includes('pTest'));
  assert.ok(svg.includes('<rect'));
  const circ = exportSVG(doc, { circular: true });
  assert.ok(circ.includes('<circle'));
  const png = exportPNG(doc, { width: 300 });
  assert.deepEqual([...png.slice(0, 8)], [137, 80, 78, 71, 13, 10, 26, 10]);
  const dv = new DataView(png.buffer, png.byteOffset, png.byteLength);
  assert.equal(dv.getUint32(16, false), 300);
  assert.ok(png.length > 1000);
  const pdf = exportPDF(doc);
  assert.ok(pdf.startsWith('%PDF-1.4'));
  assert.ok(pdf.includes('/Type /Page'));
  assert.ok(pdf.includes('xref'));
  assert.ok(pdf.trimEnd().endsWith('%%EOF'));
  assert.equal(typeof exportDiagram(doc, 'svg'), 'string');
  assert.throws(() => exportDiagram(doc, 'bmp'), /Unknown diagram format/);
});

test('project archives round-trip losslessly with stable IDs', () => {
  const doc = sampleDoc();
  doc.cloning = { schema: CURRENT_SCHEMAS.cloning, type: 'restriction', operation: 'restriction-ligation', product: { name: 'p', seqLen: 100 }, inputs: [], junctions: [] };
  doc.assembly = { schema: CURRENT_SCHEMAS.assembly, method: 'homology', product: { name: 'a', seqLen: 100, digest: 'x' }, fragments: [], junctions: [] };
  doc.alignment = alignmentRecord({ schema: 'bublagene/alignment/0', target: 'x', score: 4, identity: 90 });
  const json = projectToJSON([doc], { settings: { theme: 'dark', view: 'map' }, alignments: new Map([[doc.id, [{ name: 'r1', seq: 'ATGC' }]]]), active: 0 });
  const loaded = projectFromJSON(json);
  assert.equal(loaded.docs.length, 1);
  const back = loaded.docs[0];
  assert.equal(back.id, doc.id);
  assert.equal(back.name, doc.name);
  assert.equal(back.seq, doc.seq);
  assert.equal(back.circular, true);
  assert.deepEqual(back.features, doc.features);
  assert.deepEqual(back.primers, doc.primers);
  assert.deepEqual(back.notes, doc.notes);
  assert.equal(back.cloning.schema, CURRENT_SCHEMAS.cloning);
  assert.equal(back.assembly.schema, CURRENT_SCHEMAS.assembly);
  assert.equal(back.alignment.schema, ALIGNMENT_SCHEMA);
  assert.equal(loaded.settings.theme, 'dark');
  assert.equal(loaded.alignments[0].reads[0].name, 'r1');
  assert.ok(projectDigest(loaded).length > 0);
});

test('older project and native records migrate forward', () => {
  const legacyProject = {
    schema: 'bublagene/project/0',
    documents: [{ id: 'd1', name: 'old', seq: 'ATGC', features: [{ name: 'f', type: 'misc_feature', segments: [{ start: 0, end: 2 }] }] }],
    settings: { theme: 'light' },
  };
  const loaded = projectFromJSON(legacyProject);
  assert.equal(loaded.schema, PROJECT_SCHEMA);
  assert.equal(loaded.docs[0].name, 'old');
  assert.ok(loaded.docs[0].features[0].id);
  assert.throws(() => projectFromJSON({ schema: 'bublagene/project/99' }), /Unsupported project schema/);

  const cloning = migrateRecord({ schema: 'bublagene/cloning/0', fragments: [{ name: 'a', seqLen: 10 }] }, 'cloning');
  assert.equal(cloning.schema, CURRENT_SCHEMAS.cloning);
  assert.equal(cloning.inputs.length, 1);
  assert.equal(cloning.operation, 'restriction');
  const assembly = migrateRecord({ schema: 'bublagene/assembly/0', name: 'x', seqLen: 5 }, 'assembly');
  assert.equal(assembly.schema, CURRENT_SCHEMAS.assembly);
  assert.equal(assembly.method, 'homology');
  const alignment = migrateRecord({ schema: 'bublagene/alignment/0', reference: 'acgt', read: 'acgt' }, 'alignment');
  assert.equal(alignment.schema, ALIGNMENT_SCHEMA);
  assert.equal(alignment.target, 'acgt');
  assert.equal(alignment.query, 'acgt');
});

test('duplicate document IDs are replaced rather than merged', () => {
  const a = newDocument({ name: 'a', seq: 'ATGC' });
  const b = newDocument({ name: 'b', seq: 'ATGC' });
  b.id = a.id;
  const loaded = projectFromJSON(projectToJSON([a, b]));
  assert.equal(loaded.docs.length, 2);
  assert.notEqual(loaded.docs[0].id, loaded.docs[1].id);
  assert.ok(loaded.warnings.some((w) => /duplicate document ID/i.test(w)));
});
