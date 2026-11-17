import test from 'node:test';
import assert from 'node:assert/strict';
import {
  ENZYMES, ENZYME_BY_NAME, enzymeInfo, isoschizomersOf, enzymeFromJSON, enzymesFromJSON,
  enzymesToJSON, mergeEnzymeSets, resolveEnzymes, findSites, digest, cutEnd, uniqueCutters,
  cuttersByCount, fragmentTable, typeIISOverhangs, overhangMatrix, enzymesProducingOverhang,
  ENZYME_SCHEMA,
} from '../js/enzymes.js';
import { revcomp } from '../js/seq.js';

test('enzyme metadata: types, overhangs, methylation and isoschizomers', () => {
  const bsa = enzymeInfo(ENZYME_BY_NAME.BsaI);
  assert.equal(bsa.type, 'IIS');
  assert.equal(bsa.overhang, '5');
  assert.equal(bsa.overhangLength, 4);
  assert.equal(bsa.outsideSite, true);
  const eco = enzymeInfo(ENZYME_BY_NAME.EcoRI);
  assert.equal(eco.type, 'IIP');
  assert.equal(eco.overhangLength, 4);
  assert.equal(eco.palindromic, true);
  assert.equal(enzymeInfo(ENZYME_BY_NAME.PstI).overhang, '3');
  assert.equal(enzymeInfo(ENZYME_BY_NAME.SmaI).overhang, 'blunt');
  assert.ok(ENZYME_BY_NAME.ClaI.methylation.includes('dam'));
  assert.ok(!ENZYME_BY_NAME.BamHI.methylation.includes('dam'), 'BamHI is not Dam-sensitive');
  assert.deepEqual(isoschizomersOf('BsmBI'), ['Esp3I']);
  assert.ok(isoschizomersOf('BspQI').includes('SapI'));
  assert.equal(ENZYMES.length, 94);
});

test('custom enzyme sets: JSON in, JSON out, validation and merging', () => {
  const custom = enzymeFromJSON({ name: 'XhoI2', site: 'CTGGAG', top: 1, bottom: 5 });
  assert.equal(custom.name, 'XhoI2');
  assert.equal(custom.overhang, '5');
  const set = enzymesFromJSON({
    schema: ENZYME_SCHEMA,
    enzymes: [
      { name: 'XhoI2', site: 'CTGGAG', top: 1, bottom: 5 },
      { name: 'Plain', site: 'GGCC', top: 2, bottom: 2 },
    ],
  });
  assert.equal(set.length, 2);
  assert.equal(set[1].overhang, 'blunt');
  const json = enzymesToJSON(set);
  assert.equal(json.schema, ENZYME_SCHEMA);
  const again = enzymesFromJSON(json);
  assert.deepEqual(again.map((e) => [e.name, e.site, e.top, e.bottom]), set.map((e) => [e.name, e.site, e.top, e.bottom]));
  assert.throws(() => enzymeFromJSON({ name: 'bad', site: 'ZZZ', top: 1, bottom: 2 }), /invalid site/);
  assert.throws(() => enzymesFromJSON({ schema: 'other/1', enzymes: [] }), /schema/);
  const merged = mergeEnzymeSets(set, [enzymeFromJSON({ name: 'XhoI2', site: 'CTCGAG', top: 1, bottom: 5 })]);
  assert.equal(merged.find((e) => e.name === 'XhoI2').site, 'CTCGAG');
  assert.equal(resolveEnzymes(['EcoRI', custom]).length, 2);
});

test('unique cutters, count filters and per-enzyme fragment tables', () => {
  const seq = 'A'.repeat(100) + 'GAATTC' + 'A'.repeat(200) + 'GGATCC' + 'A'.repeat(100) + 'GGATCC' + 'A'.repeat(200);
  const unique = uniqueCutters(seq, false, ['EcoRI', 'BamHI', 'HindIII']);
  assert.deepEqual(unique.map((r) => r.enzyme.name), ['EcoRI']);
  const twice = cuttersByCount(seq, false, ['EcoRI', 'BamHI'], 2);
  assert.deepEqual(twice.map((r) => r.enzyme.name), ['BamHI']);
  const table = fragmentTable(seq, false, ['EcoRI', 'BamHI']);
  const eco = table.find((r) => r.name === 'EcoRI');
  assert.equal(eco.count, 1);
  assert.equal(eco.fragments.length, 2);
  assert.equal(eco.fragments.reduce((n, f) => n + f.len, 0), seq.length);
  const bam = table.find((r) => r.name === 'BamHI');
  assert.equal(bam.fragments.length, 3);
  assert.equal(bam.fragments.reduce((n, f) => n + f.len, 0), seq.length);
});

test('digest of a known linear plasmid fixture', () => {
  let seq = 'A'.repeat(1000);
  seq = seq.slice(0, 100) + 'GAATTC' + seq.slice(106);
  seq = seq.slice(0, 400) + 'GGATCC' + seq.slice(406);
  seq = seq.slice(0, 700) + 'AAGCTT' + seq.slice(706);
  const cuts = ['EcoRI', 'BamHI', 'HindIII'].flatMap((n) => findSites(seq, false, ENZYME_BY_NAME[n]));
  const frags = digest(seq, false, cuts);
  assert.deepEqual(frags.map((f) => f.len), [101, 300, 300, 299]);
  assert.equal(frags.reduce((n, f) => n + f.len, 0), 1000);
  const circular = digest(seq, true, cuts.map((c) => ({ ...c })));
  assert.equal(circular.length, 3);
  assert.equal(circular.reduce((n, f) => n + f.len, 0), 1000);
  // custom enzyme works through the same digest path
  const custom = enzymeFromJSON({ name: 'XhoI2', site: 'CTGGAG', top: 1, bottom: 5 });
  const customSeq = 'AAAA' + 'CTGGAG' + 'AAAA';
  const site = findSites(customSeq, false, custom);
  assert.equal(site.length, 1);
  assert.deepEqual(cutEnd(customSeq, false, site[0]), { type: '5', oh: 'TGGA' });
});

test('Type IIS overhangs and compatibility matrix', () => {
  const seq = 'GGTCTC' + 'ACGT' + 'TTTT' + 'GGTCTC' + 'ACGT' + 'TTTT';
  const rows = typeIISOverhangs(seq, false, [ENZYME_BY_NAME.BsaI]);
  assert.equal(rows.length, 2);
  assert.equal(rows[0].overhang, 'CGTT');
  assert.equal(rows[0].type, '5');
  assert.equal(rows[0].length, 4);
  assert.equal(rows[1].overhang, 'CGTT');
  assert.equal(rows[1].reverseComplement, 'AACG');
  const matrix = overhangMatrix(seq, false, [ENZYME_BY_NAME.BsaI]);
  assert.deepEqual(matrix.overhangs, ['CGTT']);
  assert.equal(matrix.matrix[0][0], 'ligate');
  assert.equal(matrix.groups.CGTT.length, 2);
  assert.equal(enzymesProducingOverhang(seq, false, 'CGTT', [ENZYME_BY_NAME.BsaI]).length, 2);
  assert.equal(enzymesProducingOverhang(seq, false, 'AAAA', [ENZYME_BY_NAME.BsaI]).length, 0);

  // The reverse-strand site reports the reverse-complement overhang.
  const rc = revcomp(seq);
  const rev = typeIISOverhangs(rc, false, [ENZYME_BY_NAME.BsaI]).filter((r) => r.strand === -1);
  assert.ok(rev.length >= 1);
  assert.equal(rev[0].overhang, revcomp(rows[0].overhang));
});
