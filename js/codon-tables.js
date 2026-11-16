// NCBI genetic-code tables.
//
// Each table stores its 64 amino acids in the canonical NCBI order: first base T,
// C, A, G; then second base T, C, A, G; then third base T, C, A, G. `starts`
// lists the codons recognised as initiators; when translating with `initiator`,
// the first codon is rendered as M (as GenBank /translation records do).
//
// Tables included: 1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 16, 21, 22, 23, 24,
// 25, 26, 27, 28, 29, 30, 31 and 33. The strings are the NCBI `ncbieaa` lines;
// `validateTable()` re-derives anchors so a typo cannot pass silently.

const ORDER = 'TCAG';

function codonList() {
  const out = [];
  for (const a of ORDER) for (const b of ORDER) for (const c of ORDER) out.push(a + b + c);
  return out;
}

const CODONS = codonList();

const RAW = {
  1: ['Standard', 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG', ['ATG']],
  2: ['Vertebrate Mitochondrial', 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG', ['ATT', 'ATC', 'ATA', 'ATG', 'GTG']],
  3: ['Yeast Mitochondrial', 'FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG', ['ATA', 'ATG', 'GTG']],
  4: ['Mold/Protozoan/Coelenterate Mitochondrial; Mycoplasma/Spiroplasma', 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG', ['TTA', 'TTG', 'CTG', 'ATT', 'ATC', 'ATA', 'ATG', 'GTG']],
  5: ['Invertebrate Mitochondrial', 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG', ['TTG', 'ATT', 'ATC', 'ATA', 'ATG', 'GTG']],
  6: ['Ciliate/Dasycladacean/Hexamita Nuclear', 'FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG', ['ATG']],
  9: ['Echinoderm/Flatworm Mitochondrial', 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG', ['ATG', 'GTG']],
  10: ['Euplotid Nuclear', 'FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG', ['ATG']],
  11: ['Bacterial/Archaeal/Plant Plastid', 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG', ['TTG', 'CTG', 'ATT', 'ATC', 'ATA', 'ATG', 'GTG']],
  12: ['Alternative Yeast Nuclear', 'FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG', ['CTG', 'ATG']],
  13: ['Ascidian Mitochondrial', 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG', ['TTG', 'ATA', 'ATG', 'GTG']],
  14: ['Alternative Flatworm Mitochondrial', 'FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG', ['ATG']],
  16: ['Chlorophycean Mitochondrial', 'FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG', ['ATG']],
  21: ['Trematode Mitochondrial', 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG', ['ATG', 'GTG']],
  22: ['Scenedesmus obliquus Mitochondrial', 'FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG', ['ATG']],
  23: ['Thraustochytrium Mitochondrial', 'FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG', ['ATT', 'ATG', 'GTG']],
  24: ['Pterobranchia Mitochondrial', 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG', ['TTG', 'CTG', 'ATG', 'GTG']],
  25: ['Candidate Division SR1 and Gracilibacteria', 'FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG', ['TTG', 'ATG', 'GTG']],
  26: ['Pachysolen tannophilus Nuclear', 'FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG', ['CTG', 'ATG']],
  27: ['Karyorelict Nuclear', 'FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG', ['ATG']],
  28: ['Condylostoma Nuclear', 'FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG', ['ATG']],
  29: ['Mesodinium Nuclear', 'FFLLSSSSYYYYCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG', ['ATG']],
  30: ['Peritrich Nuclear', 'FFLLSSSSYYEECC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG', ['ATG']],
  31: ['Blastocrithidia Nuclear', 'FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG', ['ATG']],
  33: ['Cephalodiscidae Mitochondrial', 'FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG', ['TTG', 'CTG', 'ATG', 'GTG']],
};

const ALLOWED_AA = new Set('ACDEFGHIKLMNPQRSTVWY*X');

function build(id, [name, aas, starts]) {
  if (aas.length !== 64) throw new Error(`Codon table ${id} has ${aas.length} entries, expected 64`);
  const codons = {};
  const stops = [];
  for (let i = 0; i < 64; i++) {
    const aa = aas[i];
    codons[CODONS[i]] = aa;
    if (aa === '*') stops.push(CODONS[i]);
  }
  const startSet = new Set(starts);
  return { id, name, aas, codons, starts: [...startSet], stops };
}

export const CODON_TABLES = Object.freeze(
  Object.fromEntries(Object.entries(RAW).map(([id, raw]) => [Number(id), Object.freeze(build(Number(id), raw))])),
);

export const TABLE_IDS = Object.freeze(Object.keys(CODON_TABLES).map(Number).sort((a, b) => a - b));

export const DEFAULT_TABLE = 1;

export function getCodonTable(id = DEFAULT_TABLE) {
  const table = CODON_TABLES[Number(id)];
  if (!table) throw new Error(`Unknown genetic code table ${id}; known tables: ${TABLE_IDS.join(', ')}`);
  return table;
}

export function isKnownTable(id) {
  return Boolean(CODON_TABLES[Number(id)]);
}

// Structural and semantic validation. The semantic anchors below are hand-checked
// facts (not copied output): if a table string is edited incorrectly one of them
// changes and validation fails.
const ANCHORS = {
  1: { TGA: '*', TGG: 'W', ATG: 'M', TAA: '*', AGA: 'R' },
  2: { TGA: 'W', ATA: 'M', AGA: '*', AGG: '*', ATG: 'M' },
  3: { TGA: 'W', ATA: 'M', CTA: 'T', AGA: 'R' },
  4: { TGA: 'W', AGA: 'R', TTG: 'L' },
  5: { TGA: 'W', ATA: 'M', AGA: 'S', AGG: 'S' },
  6: { TAA: 'Q', TAG: 'Q', TGA: '*' },
  9: { TGA: 'W', AAA: 'N', AGA: 'S' },
  10: { TGA: 'C' },
  11: { TGA: '*', TTG: 'L' },
  12: { CTG: 'S', TGA: '*' },
  13: { TGA: 'W', ATA: 'M', AGA: 'G', AGG: 'G' },
  14: { TAA: 'Y', TAG: '*', TGA: 'W', AAA: 'N', AGA: 'S' },
  16: { TAG: 'L', TGA: '*' },
  21: { TGA: 'W', ATA: 'M', AAA: 'N', AGA: 'S' },
  22: { TCA: '*', TAG: 'L' },
  23: { TTA: '*', TGA: '*' },
  24: { TGA: 'W', AGA: 'S', AGG: 'K' },
  25: { TGA: 'G' },
  26: { CTG: 'A' },
  27: { TAA: 'Q', TAG: 'Q', TGA: 'W' },
  28: { TAA: 'Q', TAG: 'Q', TGA: 'W' },
  29: { TAA: 'Y', TAG: 'Y', TGA: '*' },
  30: { TAA: 'E', TAG: 'E', TGA: '*' },
  31: { TAA: 'E', TAG: 'E', TGA: 'W' },
  33: { TAA: 'Y', TAG: '*', TGA: 'W', AGA: 'S', AGG: 'K' },
};

export function validateTable(id) {
  const table = getCodonTable(id);
  const issues = [];
  if (table.aas.length !== 64) issues.push(`expected 64 codons, found ${table.aas.length}`);
  for (const aa of table.aas) if (!ALLOWED_AA.has(aa)) issues.push(`invalid amino acid ${JSON.stringify(aa)}`);
  for (const [codon, aa] of Object.entries(ANCHORS[id] ?? {})) {
    if (table.codons[codon] !== aa) issues.push(`${codon} is ${table.codons[codon]}, expected ${aa}`);
  }
  for (const codon of table.starts) if (!/^[ACGT]{3}$/.test(codon)) issues.push(`invalid start codon ${codon}`);
  return { id, name: table.name, ok: issues.length === 0, issues };
}

export function validateAllTables() {
  return TABLE_IDS.map(validateTable);
}

// Translate DNA/RNA with a table. `stopSymbol` defaults to `*`; `initiator`
// renders the first codon as M when it is a recognised start codon.
export function translateWithTable(dna, tableId = DEFAULT_TABLE, { initiator = false, stopSymbol = '*' } = {}) {
  const table = getCodonTable(tableId);
  let s = String(dna).toUpperCase().replace(/U/g, 'T');
  let out = '';
  for (let i = 0; i + 3 <= s.length; i += 3) {
    let aa = table.codons[s.slice(i, i + 3)] ?? 'X';
    if (i === 0 && initiator && table.starts.includes(s.slice(0, 3))) aa = 'M';
    out += aa === '*' ? stopSymbol : aa;
  }
  return out;
}
