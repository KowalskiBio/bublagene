// A small reader for PDB files: the Cα trace of each chain, with the secondary
// structure the file declares, which is all a backbone view needs.

const SS = { H: 'helix', G: 'helix', I: 'helix', E: 'sheet', B: 'sheet' };

export function parsePdb(text, { name = 'structure' } = {}) {
  const chains = new Map();
  const helices = [];
  const sheets = [];
  let title = '';
  let model = 0;
  for (const line of String(text).split(/\r?\n/)) {
    const rec = line.slice(0, 6).trim();
    if (rec === 'TITLE') title += `${line.slice(10).trim()} `;
    else if (rec === 'MODEL') model++;
    else if (rec === 'ENDMDL') break; // the first model is enough
    else if (rec === 'HELIX') helices.push({ chain: line[19], start: parseInt(line.slice(21, 25), 10), end: parseInt(line.slice(33, 37), 10) });
    else if (rec === 'SHEET') sheets.push({ chain: line[21], start: parseInt(line.slice(22, 26), 10), end: parseInt(line.slice(33, 37), 10) });
    else if (rec === 'ATOM' && line.slice(12, 16).trim() === 'CA' && (line[16] === ' ' || line[16] === 'A')) {
      const chain = line[21] || 'A';
      const x = parseFloat(line.slice(30, 38));
      const y = parseFloat(line.slice(38, 46));
      const z = parseFloat(line.slice(46, 54));
      if (![x, y, z].every(Number.isFinite)) continue;
      if (!chains.has(chain)) chains.set(chain, []);
      chains.get(chain).push({ number: parseInt(line.slice(22, 26), 10), residue: line.slice(17, 20).trim(), x, y, z, ss: 'coil' });
    }
  }
  const mark = (list, kind) => {
    for (const r of list) for (const a of chains.get(r.chain) ?? []) if (a.number >= r.start && a.number <= r.end) a.ss = kind;
  };
  mark(helices, 'helix');
  mark(sheets, 'sheet');
  const out = [...chains].map(([id, atoms]) => ({ id, atoms }));
  if (!out.length) throw new Error('No Cα atoms found: is this a PDB file?');
  return { name, title: title.trim(), chains: out, residues: out.reduce((n, c) => n + c.atoms.length, 0), models: Math.max(1, model) };
}

const THREE_TO_ONE = {
  ALA: 'A', ARG: 'R', ASN: 'N', ASP: 'D', CYS: 'C', GLN: 'Q', GLU: 'E', GLY: 'G', HIS: 'H', ILE: 'I',
  LEU: 'L', LYS: 'K', MET: 'M', PHE: 'F', PRO: 'P', SER: 'S', THR: 'T', TRP: 'W', TYR: 'Y', VAL: 'V',
};
export const chainSequence = (chain) => chain.atoms.map((a) => THREE_TO_ONE[a.residue] ?? 'X').join('');

// Centre the coordinates and return the radius that holds them all.
export function normalise(structure) {
  const all = structure.chains.flatMap((c) => c.atoms);
  const n = all.length || 1;
  const cx = all.reduce((s, a) => s + a.x, 0) / n;
  const cy = all.reduce((s, a) => s + a.y, 0) / n;
  const cz = all.reduce((s, a) => s + a.z, 0) / n;
  let radius = 1;
  for (const a of all) radius = Math.max(radius, Math.hypot(a.x - cx, a.y - cy, a.z - cz));
  return { centre: [cx, cy, cz], radius };
}
