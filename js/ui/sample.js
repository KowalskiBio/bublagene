// A small example fragment so the app opens with something to look at.
// The sequence is made up — it has no biological meaning and is not taken
// from any real construct. It exists to show the views working.

import { newDocument, newId, colorFor } from '../model.js';

const SEQ =
  'AATATAAGGCCTAACCGTTTACTGTTTGCAGCATCAGAAATTTCTATGTCGAGGTCCATA' +
  'TAAGCGTTCTCAGGGCATACCAAGTATGAATAAACGCACTTAGTAGGAGGTAAAACAAGG' +
  'TATGGGACTAGTTGCTTCGCCAACTCGTACATGGCCAGGTCTGACAACCGAACGTTTAAA' +
  'CTTATTAGCCGCCCTTCTACTTATGATCACCCCATTTGCGATGGGGCACGGTTCTATTAA' +
  'TGACGTGGCTATGGTGTTAAGGAACATGAGTCGTATGGTGATTAATATACCAAAGTATAA' +
  'ATCTACAGGGAGTGTTCATCCTGACCGGGGGGGTGGACGCCATGCATACTGGCCTAATAC' +
  'ATGCGACGAGTGCTCATCATCGAATCGACGAAAGCTAGGCGCCTCCGGCGGAAACTCAGA' +
  'GCAGAGCAGGAACAGTGGTACAAGCGGCACTTAAGGGTACGAATTCGGATCCAAGCTTCT' +
  'CGAGCACATACGGGCAGTACGCCAGATCTG';

const PARTS = [
  ['Pbub promoter', 'promoter', 40, 100, 1, 'Example promoter region'],
  ['bubRBS', 'RBS', 104, 116, 1, 'Example ribosome binding site'],
  ['fizA', 'CDS', 121, 454, 1, 'Example 110-codon reading frame'],
  ['MCS', 'misc_feature', 460, 484, 0, 'EcoRI · BamHI · HindIII · XhoI, each cutting once'],
];

export function exampleDocument() {
  const doc = newDocument({ name: 'Example fragment', seq: SEQ, circular: false });
  doc.notes.description = 'A made-up example that ships with BublaGene. The sequence is invented and means nothing biologically.';
  doc.notes.organism = 'synthetic DNA construct';
  doc.features = PARTS.map(([name, type, start, end, strand, note]) => ({
    id: newId(), name, type, strand, color: colorFor(type),
    segments: [{ start, end, color: colorFor(type) }],
    qualifiers: [{ name: 'note', values: [{ kind: 'text', value: note }] }],
    extra: {},
  }));
  doc.primers = [
    { id: newId(), name: 'fizA-F', seq: SEQ.slice(121, 143), description: 'Anneals at the start of fizA', extra: {} },
    { id: newId(), name: 'MCS-R', seq: 'CAGATCTGGCGTACTGCCCGT', description: 'Anneals downstream of the cloning site', extra: {} },
  ];
  doc.history = [{ when: new Date().toISOString(), action: 'example', detail: 'Example fragment supplied with BublaGene' }];
  doc.source = { format: 'example' };
  doc.seqEdited = true;
  return doc;
}
