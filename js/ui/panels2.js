// More inspector panels: history, plots, protein, contigs, search and libraries.

import { h, s, clear, fmt, toast } from './dom.js';

// ---- History ----------------------------------------------------------------

export function historyPanel(el, { doc, state, api }) {
  const tl = api.timeline();
  el.append(h('h3', {}, 'Undo history'));
  if (!tl.labels.length) {
    el.append(h('p', { class: 'empty' }, 'Nothing has been changed in this file yet.'));
    return;
  }
  el.append(h('p', { class: 'sub' }, 'Click a step to return to how the file looked right after it. Steps after the one you pick stay available to redo until you make a new edit.'));
  const row = (label, position, note) => h('div', {
    class: 'row', 'aria-current': String(tl.position === position), role: 'button', tabindex: '0',
    onclick: () => api.jumpHistory(position),
    onkeydown: (e) => { if (e.key === 'Enter') api.jumpHistory(position); },
  },
  h('span', { class: 'label' }, label),
  h('span', { class: 'meta' }, note));
  el.append(row('Original', 0, tl.position === 0 ? 'now' : ''));
  tl.labels.forEach((label, i) => {
    const position = i + 1;
    el.append(row(label, position, tl.position === position ? 'now' : position > tl.position ? 'redo' : ''));
  });
  el.append(h('div', { class: 'group', style: { marginTop: '10px' } },
    h('button', { class: 'btn', disabled: tl.position === 0, onclick: () => api.jumpHistory(0) }, 'Back to original'),
    h('button', { class: 'btn', disabled: tl.position === tl.labels.length, onclick: () => api.jumpHistory(tl.labels.length) }, 'Forward to latest')));
}

// ---- Plot -------------------------------------------------------------------

import { gcProfile, gcSkewProfile, cumulativeSkew, summarize } from '../tracks.js';
import { proteinProperties, hydropathyProfile, proteinMotifs, KYTE_DOOLITTLE } from '../protein.js';
import { chainSequence } from '../pdb.js';
import { proteinView } from './proteinview.js';
import { alignmentVariants } from '../msa.js';
import { formatVariant } from '../variants.js';
import { codingSequence } from '../model.js';
import { translate } from '../seq.js';

// A line plot of {position, value} points over [0, span), with the current
// selection shaded. Clicking selects a stretch around that position.
function linePlot(points, { span, width = 300, height = 130, yMin = null, yMax = null, zero = false, selection = null, onPick = null, label = '' }) {
  const pad = { l: 34, r: 8, t: 8, b: 20 };
  const root = s('svg', { viewBox: `0 0 ${width} ${height}`, width: '100%', role: 'img', 'aria-label': label, style: 'display:block;max-width:100%;touch-action:manipulation' });
  if (!points.length) return root;
  const vals = points.map((p) => p.value);
  const lo = yMin ?? Math.min(...vals);
  const hi = yMax ?? Math.max(...vals);
  const range = hi - lo || 1;
  const x = (p) => pad.l + (p / Math.max(1, span)) * (width - pad.l - pad.r);
  const y = (v) => height - pad.b - ((v - lo) / range) * (height - pad.t - pad.b);
  root.append(s('rect', { x: pad.l, y: pad.t, width: width - pad.l - pad.r, height: height - pad.t - pad.b, fill: 'var(--paper)', stroke: 'var(--line)' }));
  if (selection) root.append(s('rect', { x: x(selection.start), y: pad.t, width: Math.max(1.5, x(selection.end) - x(selection.start)), height: height - pad.t - pad.b, fill: 'var(--select)', opacity: 0.6 }));
  if (zero && lo < 0 && hi > 0) root.append(s('line', { x1: pad.l, x2: width - pad.r, y1: y(0), y2: y(0), stroke: 'var(--muted)', 'stroke-dasharray': '3 3' }));
  for (const v of [lo, (lo + hi) / 2, hi]) {
    root.append(s('text', { x: pad.l - 4, y: y(v) + 3, 'font-size': 9, fill: 'var(--muted)', 'text-anchor': 'end' }, Math.abs(v) >= 100 ? Math.round(v) : v.toFixed(Math.abs(hi - lo) < 3 ? 2 : 0)));
  }
  for (const f of [0, 0.5, 1]) root.append(s('text', { x: x(f * span), y: height - 6, 'font-size': 9, fill: 'var(--muted)', 'text-anchor': f === 0 ? 'start' : f === 1 ? 'end' : 'middle' }, fmt(Math.round(f * span))));
  const d = points.map((p, i) => `${i ? 'L' : 'M'}${x(p.position).toFixed(1)} ${y(p.value).toFixed(1)}`).join('');
  root.append(s('path', { d, fill: 'none', stroke: 'var(--accent)', 'stroke-width': 1.4 }));
  if (onPick) {
    root.style.cursor = 'crosshair';
    root.addEventListener('click', (e) => {
      const r = root.getBoundingClientRect();
      const frac = ((e.clientX - r.left) / r.width * width - pad.l) / (width - pad.l - pad.r);
      onPick(Math.max(0, Math.min(span - 1, Math.round(frac * span))));
    });
  }
  return root;
}

export function plotPanel(el, { doc, state, api }) {
  const p = state.plot;
  const L = doc.seq.length;
  el.append(h('h3', {}, 'Plots along the sequence'));
  if (!L) { el.append(h('p', { class: 'empty' }, 'There is no sequence to plot yet.')); return; }
  el.append(h('label', { class: 'field' }, h('span', {}, 'Track'),
    h('select', { id: 'plotTrack', onchange: (e) => api.setPlot({ track: e.target.value }) },
      [['gc', 'GC content (%)'], ['skew', 'GC skew (G−C)/(G+C)'], ['cumulative', 'Cumulative GC skew']].map(([v, t]) => h('option', { value: v, selected: p.track === v }, t)))));
  if (p.track !== 'cumulative') {
    el.append(h('label', { class: 'field' }, h('span', {}, `Window: ${fmt(p.window)} bp`),
      h('input', { id: 'plotWindow', type: 'range', min: 10, max: Math.max(50, Math.min(5000, L)), step: 10, value: Math.min(p.window, Math.max(50, Math.min(5000, L))), oninput: (e) => api.setPlot({ window: Number(e.target.value) }) })));
  }
  const step = Math.max(1, Math.floor(L / 600));
  const win = Math.min(p.window, L);
  const points = p.track === 'gc' ? gcProfile(doc.seq, win, step, doc.circular)
    : p.track === 'skew' ? gcSkewProfile(doc.seq, win, step, doc.circular)
      : cumulativeSkew(doc.seq, step);
  const pick = (pos) => {
    const half = Math.max(5, Math.floor(Math.min(win, L) / 2));
    api.setSelection({ start: Math.max(0, pos - half), end: Math.min(L, pos + half) });
  };
  el.append(linePlot(points, {
    span: L, selection: state.selection, onPick: pick, zero: p.track !== 'gc',
    yMin: p.track === 'gc' ? 0 : null, yMax: p.track === 'gc' ? 100 : null, label: `${p.track} along ${doc.name}`,
  }));
  const sum = summarize(points);
  const fixed = p.track === 'gc' ? 1 : p.track === 'skew' ? 3 : 0;
  el.append(h('p', { class: 'sub' }, `min ${sum.min.toFixed(fixed)} · mean ${sum.mean.toFixed(fixed)} · max ${sum.max.toFixed(fixed)}. Click the plot to select that stretch.`));
  if (p.track === 'cumulative') el.append(h('p', { class: 'sub' }, 'In bacterial genomes the minimum and maximum of this curve mark the origin and terminus of replication.'));
}

// ---- Protein ----------------------------------------------------------------

function hydropathyPlot(protein, selection) {
  const points = hydropathyProfile(protein, 9).map((q) => ({ position: q.position, value: q.value }));
  return linePlot(points, { span: protein.length, yMin: -4.5, yMax: 4.5, zero: true, selection, label: 'Hydropathy' });
}

export function proteinPanel(el, { doc, state, api }) {
  const cds = doc.features.filter((f) => f.type === 'CDS');
  const m = state.protein;
  el.append(h('h3', {}, 'Protein'));
  let protein = '';
  if (!cds.length) el.append(h('p', { class: 'empty' }, 'This file has no CDS feature. Add one from the Features panel, or use “+ CDS” beside an open reading frame.'));
  else {
    const feature = cds.find((f) => f.id === m.feature) ?? cds[0];
    el.append(h('label', { class: 'field' }, h('span', {}, 'Coding sequence'),
      h('select', { id: 'proteinCds', onchange: (e) => api.setProtein({ feature: e.target.value }) },
        cds.map((f) => h('option', { value: f.id, selected: f.id === feature.id }, f.name)))));
    protein = translate(codingSequence(doc, feature)).replace(/\*.*$/, '');
    const pr = proteinProperties(protein);
    const t = h('table', { class: 'data' });
    for (const [k, v] of [
      ['Length', `${fmt(pr.length)} aa`], ['Mass', `${(pr.mass / 1000).toFixed(2)} kDa`], ['Isoelectric point', pr.pI.toFixed(2)],
      ['Net charge at pH 7.4', pr.chargeAt7.toFixed(1)], ['Acidic / basic residues', `${pr.acidic} / ${pr.basic}`],
      ['Hydropathy (GRAVY)', pr.gravy.toFixed(2)], ['ε280, reduced', `${fmt(pr.extinction)} M⁻¹cm⁻¹`], ['ε280, cystines formed', `${fmt(pr.extinctionOxidised)} M⁻¹cm⁻¹`],
    ]) t.append(h('tr', {}, h('td', {}, k), h('td', {}, v)));
    el.append(t);
    el.append(h('h3', {}, 'Hydropathy'));
    el.append(hydropathyPlot(protein, null));
    el.append(h('p', { class: 'sub' }, 'Kyte–Doolittle, 9-residue window. Above zero is hydrophobic.'));
    const motifs = proteinMotifs(protein);
    el.append(h('h3', {}, `Motifs (${motifs.length})`));
    if (!motifs.length) el.append(h('p', { class: 'empty' }, 'No tags, cleavage sites or sequons found.'));
    for (const x of motifs) {
      el.append(h('div', { class: 'row', title: x.kind },
        h('span', { class: 'label' }, x.name),
        h('span', { class: 'meta' }, `${x.start + 1}–${x.end}`)));
    }
    el.append(h('p', { class: 'sub' }, 'A simple motif scan, not a domain prediction. Possible transmembrane helices use the classic hydropathy rule.'));
  }

  el.append(h('h3', {}, '3D structure'));
  el.append(h('p', { class: 'sub' }, 'Open a PDB file to see its backbone. Structures are not predicted here.'));
  el.append(h('button', { class: 'btn', onclick: () => api.pickFile() }, 'Open a PDB file…'));
  const structures = state.structures;
  if (!structures.length) return;
  const cur = structures.find((x) => x.id === m.structure) ?? structures[structures.length - 1];
  el.append(h('label', { class: 'field', style: { marginTop: '8px' } }, h('span', {}, 'Structure'),
    h('select', { onchange: (e) => api.setProtein({ structure: e.target.value }) },
      structures.map((x) => h('option', { value: x.id, selected: x.id === cur.id }, `${x.structure.name} · ${x.structure.residues} residues`)))));
  if (cur.structure.title) el.append(h('p', { class: 'sub' }, cur.structure.title));
  el.append(h('div', { class: 'group' },
    h('label', { class: 'field', style: { flex: '1 1 120px' } }, h('span', {}, 'Colour by'),
      h('select', { onchange: (e) => api.setProtein({ colour: e.target.value }) },
        [['structure', 'secondary structure'], ['rainbow', 'N → C terminus'], ['chain', 'chain']].map(([v, t]) => h('option', { value: v, selected: m.colour === v }, t)))),
    h('label', { class: 'field', style: { flex: '1 1 110px' } }, h('span', {}, 'Highlight residues'),
      h('input', { id: 'pdbHighlight', value: m.highlight, placeholder: 'e.g. 10-25', oninput: (e) => api.setProtein({ highlight: e.target.value }) }))));
  const range = /^\s*(\d+)\s*[-–]\s*(\d+)\s*$/.exec(m.highlight);
  el.append(proteinView(cur.structure, { colour: m.colour, highlight: range ? [Number(range[1]), Number(range[2])] : null }));
  el.append(h('p', { class: 'sub' }, 'Drag to rotate, scroll to zoom, double-click to reset.'));
  for (const chain of cur.structure.chains) {
    const seq = chainSequence(chain);
    const at = protein && seq.length > 5 ? protein.indexOf(seq.slice(0, Math.min(seq.length, 20))) : -1;
    el.append(h('p', { class: 'sub' }, `Chain ${chain.id}: ${chain.atoms.length} residues${at >= 0 ? ` · matches residue ${at + 1} onward of the selected CDS` : ''}`));
  }
  el.append(h('button', { class: 'btn', onclick: () => api.removeStructure(cur.id) }, 'Close this structure'));
}

// ---- Contigs ----------------------------------------------------------------

function layoutDiagram(contig) {
  const W = 300;
  const rowH = 7;
  const reads = contig.layout;
  const H = Math.min(260, reads.length * (rowH + 2) + 14);
  const root = s('svg', { viewBox: `0 0 ${W} ${H}`, width: '100%', role: 'img', 'aria-label': 'Read layout', style: 'display:block;max-width:100%' });
  const x = (p) => 4 + (p / contig.length) * (W - 8);
  reads.slice(0, 30).forEach((r, i) => {
    root.append(s('rect', { x: x(r.offset), y: 4 + i * (rowH + 2), width: Math.max(2, x(r.offset + r.length) - x(r.offset)), height: rowH, rx: 2, fill: r.reverse ? '#e0884c' : 'var(--accent)', opacity: 0.8 }, s('title', {}, `${r.name}${r.reverse ? ' (reverse)' : ''}`)));
  });
  return root;
}

export function contigsPanel(el, { doc, state, api, docs }) {
  const c = state.contigs;
  el.append(h('h3', {}, 'Assemble reads into contigs'));
  el.append(h('p', { class: 'sub' }, 'Pick the files that are Sanger reads (or fragments of one sequence). Reads in either orientation are joined where their ends overlap.'));
  for (const d of docs) {
    const on = c.ids.includes(d.id);
    el.append(h('div', {
      class: 'row', 'aria-current': String(on),
      onclick: () => api.setContigs({ ids: on ? c.ids.filter((x) => x !== d.id) : [...c.ids, d.id] }),
    }, h('input', { type: 'checkbox', checked: on, onclick: (e) => e.stopPropagation(), onchange: () => api.setContigs({ ids: on ? c.ids.filter((x) => x !== d.id) : [...c.ids, d.id] }) }),
    h('span', { class: 'label' }, d.name), h('span', { class: 'meta' }, `${fmt(d.seq.length)} bp`)));
  }
  el.append(h('label', { class: 'field', style: { marginTop: '8px' } }, h('span', {}, `Smallest overlap: ${c.minOverlap} bp`),
    h('input', { id: 'contigOverlap', type: 'range', min: 10, max: 150, value: c.minOverlap, oninput: (e) => api.setContigs({ minOverlap: Number(e.target.value) }) })));
  el.append(h('label', { class: 'field' }, h('span', {}, `Mismatches allowed in the overlap: ${c.maxMismatch}%`),
    h('input', { id: 'contigMismatch', type: 'range', min: 0, max: 20, value: c.maxMismatch, oninput: (e) => api.setContigs({ maxMismatch: Number(e.target.value) }) })));
  const traced = docs.filter((d) => api.traceCount(d) > 0);
  if (traced.length) {
    el.append(h('label', { class: 'sub', style: { display: 'flex', gap: '6px', alignItems: 'center', marginBottom: '8px' } },
      h('input', { type: 'checkbox', checked: c.useTraces, onchange: (e) => api.setContigs({ useTraces: e.target.checked }) }),
      'Use the chromatogram base calls (with their quality) where a file has them'));
  }
  el.append(h('button', { class: 'btn primary', disabled: c.ids.length < 1, onclick: () => api.runContigs() }, c.ids.length ? `Assemble ${c.ids.length} file${c.ids.length === 1 ? '' : 's'}` : 'Pick the reads first'));
  if (c.error) el.append(h('p', { class: 'empty' }, c.error));
  if (!c.result) return;
  el.append(h('h3', {}, `${c.result.contigs.length} contig${c.result.contigs.length === 1 ? '' : 's'} from ${c.result.readCount} reads`));
  c.result.contigs.forEach((ct, i) => {
    el.append(h('div', { class: 'row', style: { cursor: 'default' } },
      h('span', { class: 'label' }, `Contig ${i + 1}`),
      h('span', { class: 'meta' }, `${fmt(ct.length)} bp · ${ct.reads.length} read${ct.reads.length === 1 ? '' : 's'} · ${ct.depth}×`)));
    el.append(layoutDiagram(ct));
    el.append(h('div', { class: 'group', style: { margin: '4px 0 12px' } },
      h('button', { class: 'btn', onclick: () => api.openContig(i) }, 'Open as a file'),
      h('button', { class: 'btn', onclick: () => api.copy(ct.seq, 'bases') }, 'Copy sequence')));
  });
  if (c.result.unplaced?.length) el.append(h('p', { class: 'sub' }, `Not joined to anything: ${c.result.unplaced.join(', ')}`));
}

// ---- Several sequences (used by the Alignment panel) ---------------------------

const BASE_COLOR = { A: '#5aa58f', C: '#4f8fd6', G: '#d6a14a', T: '#c46f6f', U: '#c46f6f', '-': 'var(--muted)' };

function msaBlocks(result, width = 50) {
  const box = h('div', { class: 'mono', style: { fontSize: '11px', lineHeight: '1.35', maxHeight: '300px', overflow: 'auto', border: '1px solid var(--line)', borderRadius: '7px', padding: '6px', background: 'var(--paper)' } });
  const nameW = Math.min(14, Math.max(...result.names.map((n) => n.length)));
  const cols = result.consensus.length;
  for (let start = 0; start < cols; start += width) {
    const end = Math.min(cols, start + width);
    const block = h('div', { style: { whiteSpace: 'pre', marginBottom: '8px' } });
    result.rows.forEach((row, i) => {
      const line = h('div', {}, h('span', { style: { color: 'var(--muted)' } }, `${result.names[i].slice(0, nameW).padEnd(nameW)} `));
      for (let c = start; c < end; c++) {
        const ch = row[c];
        const agree = result.conservation[c] === 1;
        line.append(h('span', { style: { color: BASE_COLOR[ch] ?? 'var(--ink)', background: agree ? 'transparent' : 'var(--accent-soft)' } }, ch));
      }
      block.append(line);
    });
    block.append(h('div', { style: { color: 'var(--muted)' } }, `${'consensus'.padEnd(nameW)} ${result.consensus.slice(start, end)}`));
    block.append(h('div', { style: { color: 'var(--accent)' } }, `${''.padEnd(nameW)} ${[...result.conservation.slice(start, end)].map((v) => (v === 1 ? '*' : ' ')).join('')}`));
    box.append(block);
  }
  return box;
}

export function msaSection(el, { state, api, docs }) {
  const m = state.msa;
  el.append(h('h3', {}, 'Align several sequences'));
  el.append(h('p', { class: 'sub' }, 'Aligns every picked file against the longest one and merges the results, so insertions are kept. Good for variants and clones of one sequence.'));
  for (const d of docs) {
    const on = m.ids.includes(d.id);
    el.append(h('label', { class: 'row', style: { cursor: 'pointer' } },
      h('input', { type: 'checkbox', checked: on, onchange: () => api.setMsa({ ids: on ? m.ids.filter((x) => x !== d.id) : [...m.ids, d.id] }) }),
      h('span', { class: 'label' }, d.name), h('span', { class: 'meta' }, `${fmt(d.seq.length)} bp`)));
  }
  el.append(h('button', { class: 'btn primary', style: { marginTop: '6px' }, disabled: m.ids.length < 2, onclick: () => api.runMsa() }, m.ids.length < 2 ? 'Pick at least two files' : `Align ${m.ids.length} sequences`));
  if (m.error) el.append(h('p', { class: 'empty' }, m.error));
  if (!m.result) return;
  const r = m.result;
  el.append(h('p', { class: 'sub' }, `${r.consensus.length} columns · ${r.identity.toFixed(1)}% identical in every sequence`));
  el.append(msaBlocks(r));
  const vars = alignmentVariants(r.rows);
  el.append(h('h3', {}, `Differences from ${r.names[0]} (${vars.length})`));
  if (!vars.length) el.append(h('p', { class: 'empty' }, 'The sequences are identical.'));
  const t = h('table', { class: 'data' }, h('tr', {}, h('th', {}, 'Position'), h('th', {}, 'Type'), ...r.names.map((n) => h('th', { title: n }, n.slice(0, 6)))));
  for (const v of vars.slice(0, 80)) {
    t.append(h('tr', {}, h('td', {}, fmt(v.position + 1)), h('td', { class: 'sub' }, v.type), ...v.bases.map((b) => h('td', { class: 'mono' }, b))));
  }
  if (vars.length) el.append(t);
  if (vars.length > 80) el.append(h('p', { class: 'sub' }, `Showing 80 of ${vars.length}.`));
  el.append(h('button', { class: 'btn', style: { marginTop: '8px' }, onclick: () => api.copy(r.names.map((n, i) => `>${n}\n${r.rows[i]}`).join('\n'), 'alignment characters') }, 'Copy as FASTA'));
}

// ---- Search -----------------------------------------------------------------

export function searchPanel(el, { doc, state, api }) {
  const q = state.search;
  el.append(h('h3', {}, 'Search'));
  const run = () => api.runSearch(el.querySelector('#searchQuery').value);
  el.append(h('label', { class: 'field' }, h('span', {}, 'A sequence (either strand, IUPAC codes allowed), or a feature, primer or note'),
    h('input', {
      id: 'searchQuery', value: q.query, placeholder: 'GAATTC, lacZ, ampicillin …',
      onkeydown: (e) => { if (e.key === 'Enter') run(); },
    })));
  el.append(h('label', { class: 'field' }, h('span', {}, 'Where'),
    h('select', { id: 'searchScope', onchange: (e) => api.setSearch({ scope: e.target.value }) },
      [['file', 'This file'], ['open', 'All open files'], ['saved', 'Saved collections']].map(([v, t]) => h('option', { value: v, selected: q.scope === v }, t)))));
  el.append(h('button', { class: 'btn primary', onclick: run }, 'Search'));
  if (q.error) el.append(h('p', { class: 'empty' }, q.error));
  if (q.results) {
    const total = q.results.reduce((n, r) => n + r.items.length, 0);
    el.append(h('h3', {}, `${total} match${total === 1 ? '' : 'es'}${q.results.length > 1 ? ` in ${q.results.length} files` : ''}`));
    if (!total) el.append(h('p', { class: 'empty' }, 'Nothing found.'));
    for (const r of q.results) {
      if (q.results.length > 1 || q.scope !== 'file') el.append(h('p', { class: 'sub', style: { margin: '8px 0 2px', fontWeight: 600 } }, r.name));
      for (const it of r.items.slice(0, 60)) {
        el.append(h('div', { class: 'row', onclick: () => api.goTo(r.docId, it.start, it.end) },
          h('span', { class: 'label' }, it.label),
          h('span', { class: 'meta' }, `${fmt(it.start + 1)}–${fmt(it.end || 1)}${it.strand === -1 ? ' ←' : it.strand === 1 ? ' →' : ''}`)));
      }
      if (r.items.length > 60) el.append(h('p', { class: 'sub' }, `${r.items.length - 60} more in this file.`));
    }
  }
  if (q.saved) {
    el.append(h('h3', {}, `${q.saved.length} saved file${q.saved.length === 1 ? '' : 's'}`));
    if (!q.saved.length) el.append(h('p', { class: 'empty' }, 'No saved file matches.'));
    for (const r of q.saved) {
      el.append(h('div', { class: 'row', onclick: () => api.openSaved(r.id) },
        h('span', { class: 'label' }, h('span', {}, r.name), h('span', { class: 'sub', style: { display: 'block' } }, `${r.collection} · ${r.reasons.join('; ')}`)),
        h('span', { class: 'meta' }, 'Open')));
    }
  }

  el.append(h('h3', {}, 'Look it up elsewhere'));
  el.append(h('p', { class: 'sub' }, 'These open the other site in a new tab and send it what you typed or selected, only when you click.'));
  const selected = api.selectedSequence();
  const dna = /^[ACGTN]{12,}$/i.test(q.query.trim()) ? q.query.trim() : selected.length >= 12 ? selected : '';
  const link = (label, href, enabled) => h('a', {
    class: 'btn', href: enabled ? href : undefined, target: '_blank', rel: 'noopener noreferrer',
    style: { textDecoration: 'none', opacity: enabled ? 1 : 0.45, pointerEvents: enabled ? 'auto' : 'none' },
  }, label);
  el.append(h('div', { class: 'group' },
    link('NCBI BLAST', `https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastSearch&PROGRAM=blastn&QUERY=${encodeURIComponent(dna)}`, !!dna),
    link('Addgene', `https://www.addgene.org/search/catalog/plasmids/?q=${encodeURIComponent(q.query.trim())}`, !!q.query.trim())));
  el.append(h('p', { class: 'sub' }, dna ? `BLAST will receive ${fmt(dna.length)} bases.` : 'BLAST needs 12 or more bases: type them above or select them in the sequence.'));
}

// ---- Library ----------------------------------------------------------------

export function libraryPanel(el, { doc, state, api, docs }) {
  const lib = state.library;
  if (lib.records === null && !lib.loading) api.refreshLibrary();
  el.append(h('h3', {}, 'Saved collections'));
  el.append(h('p', { class: 'sub' }, 'Files kept in this browser so they are one click away, and searchable from the Search panel. They are never uploaded.'));
  el.append(h('label', { class: 'field' }, h('span', {}, 'Collection'),
    h('input', { id: 'libraryName', value: lib.name, list: 'libraryNames', oninput: (e) => api.setLibraryName(e.target.value) })));
  const names = [...new Set((lib.records ?? []).map((r) => r.collection))];
  el.append(h('datalist', { id: 'libraryNames' }, names.map((n) => h('option', { value: n }))));
  el.append(h('div', { class: 'group' },
    h('button', { class: 'btn primary', onclick: () => api.saveCurrentToCollection() }, `Save ${doc.name}`),
    docs.length > 1 ? h('button', { class: 'btn', onclick: () => api.saveAllToCollection() }, `Save all ${docs.length} open files`) : null));
  if (lib.error) el.append(h('p', { class: 'empty' }, lib.error));
  if (lib.loading) el.append(h('p', { class: 'empty' }, 'Loading…'));
  const records = lib.records ?? [];
  if (!lib.loading && !records.length && !lib.error) el.append(h('p', { class: 'empty' }, 'Nothing saved yet.'));
  let last = null;
  for (const r of records) {
    if (r.collection !== last) { last = r.collection; el.append(h('p', { class: 'sub', style: { margin: '12px 0 2px', fontWeight: 600 } }, r.collection)); }
    el.append(h('div', { class: 'row', onclick: () => api.openSaved(r.id) },
      h('span', { class: 'label' }, h('span', {}, r.name), h('span', { class: 'sub', style: { display: 'block' } }, `${fmt(r.bp)} bp · ${r.circular ? 'circular' : 'linear'} · ${r.featureNames.length} features`)),
      h('button', { class: 'btn', style: { padding: '1px 7px' }, onclick: (e) => { e.stopPropagation(); api.deleteSaved(r.id); } }, 'Remove')));
  }
  if (records.length) {
    el.append(h('div', { class: 'group', style: { marginTop: '12px' } },
      h('button', { class: 'btn', onclick: () => api.exportCollection() }, 'Export these as a project file…')));
    el.append(h('p', { class: 'sub' }, 'To share them, send that project file: anyone can open it with Open project in the File panel.'));
  }
  el.append(h('h3', {}, 'Also saved here'));
  const sets = Object.keys(api.enzymeSets()).length;
  const elements = api.userLibrary().length;
  el.append(h('p', { class: 'sub' }, `${sets} enzyme set${sets === 1 ? '' : 's'} (Enzymes panel) and ${elements} detection element${elements === 1 ? '' : 's'} (Detect panel).`));
}
