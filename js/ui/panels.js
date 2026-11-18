// @ts-nocheck — legacy DOM wiring; the engines it calls are type-checked.
// The inspector panels on the right.

import { h, s, clear, fmt, toast } from './dom.js';
import { ENZYMES, ENZYME_BY_NAME, COMMON_SET, findSites, digest, overhangLabel, migration, LADDERS } from '../enzymes.js';
import { findBindingSites, designPrimers, simulatePCR } from '../primers.js';
import { featureSpan, featureLength, codingSequence, colorFor, FEATURE_TYPES, featureLabel, newId, qualifierText, addHistory, recomputeFeature, setQualifier } from '../model.js';
import { gcContent, meltingTemp, translate, sixFrameTranslation, revcomp, subseq, rangeLength, cleanSequence } from '../seq.js';
import { filterFeatures, searchFeatures } from '../search.js';
import { alignPair, dotPlot } from '../align.js';
import { findGuides, guidesToFeatures, nhejEdit, hdrEdit, editPrimers } from '../crispr.js';
import { callVariants, variantsToFeatures, formatVariant } from '../variants.js';
import { structureToSvg } from '../fold.js';
import { describeSegments } from '../formats/dna.js';
import { cloningPanel } from './cloningpanel.js';
import { historyPanel, plotPanel, proteinPanel, contigsPanel, searchPanel, libraryPanel, msaSection } from './panels2.js';
import { aminoAcidChange } from '../mutagenesis.js';
import { orfsFor } from './overlays.js';
import { methylationSummary, METHYLATION_LABEL, HOST_NOTE } from '../methylation.js';
import { primerProperties, pairDimer } from '../thermo.js';

// Panels, grouped for the tab bar. Order within a group is the order shown.
export const PANEL_GROUPS = [
  ['Sequence', ['Features', 'Detect', 'Plot', 'Protein']],
  ['Primers', ['Primers', 'Mutate', 'CRISPR']],
  ['Enzymes', ['Enzymes', 'Digest']],
  ['Cloning', ['Cloning']],
  ['Analyse', ['Alignment', 'Contigs', 'Search', 'Structure']],
  ['File', ['File', 'Library', 'History']],
];
export const PANELS = PANEL_GROUPS.flatMap(([, names]) => names);

export function renderPanel(name, ctx) {
  const el = h('div');
  ({
    Features: features, Detect: detectPanel, Primers: primers, Mutate: mutatePanel, Enzymes: enzymes, Digest: digestPanel, Cloning: cloningPanel,
    Alignment: alignmentPanel, CRISPR: crisprPanel, Structure: structurePanel, File: file,
    History: historyPanel, Plot: plotPanel, Protein: proteinPanel, Contigs: contigsPanel, Search: searchPanel, Library: libraryPanel,
  })[name](el, ctx);
  return el;
}

// ---- Features ---------------------------------------------------------------

function features(el, { doc, state, api }) {
  const L = doc.seq.length;
  const sel = state.selection;
  el.append(h('h3', {}, `Features (${doc.features.length})`));
  if (!doc.features.length) el.append(h('p', { class: 'empty' }, 'No features yet. Select part of the sequence and add one.'));

  el.append(h('label', { class: 'field' }, h('span', {}, 'Filter'), h('input', {
    id: 'featureFilter', value: state.featureFilter ?? '', placeholder: 'name, type or qualifier',
    oninput: (e) => api.setFeatureFilter(e.target.value),
  })));
  const filtered = filterFeatures(doc, { text: state.featureFilter || null });
  if (state.featureFilter) el.append(h('p', { class: 'sub' }, `${filtered.length} of ${doc.features.length} feature(s)`));
  const sorted = [...filtered].sort((a, b) => featureSpan(a).start - featureSpan(b).start);
  for (const f of sorted) {
    const span = featureSpan(f);
    const picked = state.selectedFeatures?.has(f.id);
    const row = h('div', {
      class: 'row', 'aria-current': String(state.highlight === f.id), draggable: 'true',
      ondragstart: (e) => { e.dataTransfer.setData('text/plain', f.id); e.dataTransfer.effectAllowed = 'move'; },
      ondragover: (e) => e.preventDefault(),
      ondrop: (e) => { e.preventDefault(); const id = e.dataTransfer.getData('text/plain'); if (id && id !== f.id) api.reorderFeature(id, f.id); },
      onclick: (e) => { if (e.metaKey || e.ctrlKey) api.toggleFeatureSelected(f.id); else api.focusFeature(f); },
    },
      h('input', {
        type: 'checkbox', checked: picked, title: 'Select feature',
        onclick: (e) => { e.stopPropagation(); api.toggleFeatureSelected(f.id); },
      }),
      h('span', { class: 'swatch', style: { background: f.color } }),
      h('span', { class: 'label' },
        h('span', {}, f.name),
        h('span', { class: 'sub', style: { marginLeft: '6px' } }, f.type)),
      h('span', { class: 'meta' }, `${fmt(span.start + 1)}–${fmt(span.end)} ${f.strand === -1 ? '←' : f.strand === 1 ? '→' : ''}`));
    el.append(row);
    if (state.highlight === f.id) el.append(featureDetail(f, { doc, state, api }));
  }

  const selected = doc.features.filter((f) => state.selectedFeatures?.has(f.id));
  if (selected.length) {
    el.append(h('h3', {}, `Selected (${selected.length})`));
    el.append(h('div', { class: 'group' },
      h('button', { class: 'btn', onclick: () => api.copySelectedFeatures() }, 'Copy sequences'),
      h('button', { class: 'btn', title: 'Copy the annotations, to paste elsewhere', onclick: () => api.copyFeatureAnnotations() }, 'Copy features'),
      h('button', { class: 'btn', title: 'Copy the annotations and remove them here', onclick: () => api.cutFeatureAnnotations() }, 'Cut features'),
      h('button', { class: 'btn', onclick: () => api.setFeatureSelection(doc.features.map((f) => f.id)) }, 'Select all'),
      h('button', { class: 'btn', danger: true, onclick: () => api.deleteSelectedFeatures() }, 'Delete selected'),
      h('button', { class: 'btn', onclick: () => api.setFeatureSelection([]) }, 'Clear')));
    el.append(h('p', { class: 'sub' }, '⌘-click rows to add or remove them; drag a row onto another to reorder the list.'));
  }

  if (state.featureClip) {
    const n = state.featureClip.features.length;
    const where = sel ? fmt(sel.start + 1) : state.caret !== null && state.caret !== undefined ? fmt(state.caret + 1) : 'where they were copied from';
    el.append(h('div', { class: 'group', style: { marginTop: '8px' } },
      h('button', { class: 'btn', onclick: () => api.pasteFeatureAnnotations() }, `Paste ${n} feature${n === 1 ? '' : 's'} at ${where}`)));
  }

  el.append(h('h3', {}, 'Add'));
  el.append(h('button', {
    class: 'btn', disabled: !sel,
    onclick: () => {
      const type = el.querySelector('#newType').value;
      const name = el.querySelector('#newName').value.trim() || type;
      api.edit(`Add feature ${name}`, (d) => {
        d.features.push({
          id: newId(), name, type, strand: Number(el.querySelector('#newStrand').value),
          color: colorFor(type), segments: [{ start: sel.start, end: sel.end }], qualifiers: [], extra: {},
        });
      });
    },
  }, sel ? `Add feature over ${fmt(sel.start + 1)}–${fmt(sel.end)}` : 'Select a region first'));
  el.append(h('label', { class: 'field', style: { marginTop: '8px' } }, h('span', {}, 'Name'), h('input', { id: 'newName', placeholder: 'optional' })));
  el.append(h('label', { class: 'field' }, h('span', {}, 'Type'),
    h('select', { id: 'newType' }, FEATURE_TYPES.map((t) => h('option', { value: t }, t)))));
  el.append(h('label', { class: 'field' }, h('span', {}, 'Strand'),
    h('select', { id: 'newStrand' }, [['1', 'forward →'], ['-1', 'reverse ←'], ['0', 'none']].map(([v, t]) => h('option', { value: v }, t)))));

  el.append(h('h3', {}, 'Six-frame translation'));
  {
    const region = sel ? api.selectedSequence() : doc.seq.slice(0, 600);
    const frames = sixFrameTranslation(region);
    for (const f of frames) {
      if (!f.protein.length) continue;
      el.append(h('div', { class: 'row', title: `${f.strand === 1 ? '+' : '−'}${f.frame}` },
        h('span', { class: 'label mono', style: { fontSize: '11px' } }, `${f.strand === 1 ? '+' : '−'}${f.frame} ${f.protein.slice(0, 120)}`)));
    }
    if (region.length < doc.seq.length) el.append(h('p', { class: 'sub' }, 'Select a region to translate just that part.'));
  }

  const orfs = orfsFor(doc, 100);
  el.append(h('h3', {}, `Open reading frames (${orfs.length})`));
  if (!orfs.length) el.append(h('p', { class: 'empty' }, 'None of at least 100 codons.'));
  for (const o of orfs.slice(0, 12)) {
    el.append(h('div', {
      class: 'row', onclick: () => api.setSelection({ start: o.start, end: o.end }),
    },
      h('span', { class: 'label mono' }, `${o.strand === 1 ? '→' : '←'} ${fmt(o.start + 1)}–${fmt(o.end || L)}`),
      h('span', { class: 'meta' }, `${fmt(o.codons)} codons`),
      h('button', {
        class: 'btn', style: { padding: '1px 7px', fontSize: '12px' }, title: 'Add this reading frame as a CDS feature',
        onclick: (e) => { e.stopPropagation(); api.addOrfAsCds(o); },
      }, '+ CDS')));
  }
}

function qualifierLines(f) {
  return (f.qualifiers ?? []).flatMap((q) => q.values.map((v) => `${q.name}=${v.value}`)).join('\n');
}

function parseQualifierLines(text) {
  const out = [];
  for (const line of String(text).split(/\r?\n/)) {
    if (!line.trim()) continue;
    const eq = line.indexOf('=');
    const name = (eq < 0 ? line : line.slice(0, eq)).trim();
    const value = eq < 0 ? '' : line.slice(eq + 1).trim();
    if (!name) continue;
    const prev = out[out.length - 1];
    if (prev && prev.name === name) prev.values.push({ kind: /^-?\d+$/.test(value) ? 'int' : 'text', value });
    else out.push({ name, values: [{ kind: /^-?\d+$/.test(value) ? 'int' : 'text', value }] });
  }
  return out;
}

function featureDetail(f, { doc, api }) {
  const L = doc.seq.length;
  const box = h('div', { style: { padding: '4px 8px 12px', borderLeft: '2px solid var(--accent)', margin: '2px 0 8px 6px' } });
  const find = (d) => d.features.find((x) => x.id === f.id);
  box.append(h('label', { class: 'field' }, h('span', {}, 'Name'),
    h('input', { value: f.name, onchange: (e) => api.edit('Rename feature', (d) => { find(d).name = e.target.value; }) })));
  box.append(h('label', { class: 'field' }, h('span', {}, 'Type'),
    h('select', {
      onchange: (e) => api.edit('Change feature type', (d) => {
        const t = find(d);
        t.type = e.target.value;
        t.color = colorFor(e.target.value);
        for (const seg of t.segments) seg.color = t.color;
      }),
    }, FEATURE_TYPES.map((t) => h('option', { value: t, selected: t === f.type }, t)))));
  box.append(h('div', { class: 'group', style: { alignItems: 'flex-end' } },
    h('label', { class: 'field', style: { flex: '1 1 110px', marginBottom: '0' } },
      h('span', {}, 'Strand'),
      h('select', {
        onchange: (e) => api.edit('Change strand', (d) => { find(d).strand = Number(e.target.value); }),
      }, [['1', 'forward →'], ['-1', 'reverse ←'], ['0', 'none']].map(([v, t]) => h('option', { value: v, selected: Number(v) === f.strand }, t)))),
    h('label', { class: 'field', style: { flex: '0 0 70px', marginBottom: '0' } },
      h('span', {}, 'Colour'),
      h('input', {
        type: 'color', value: /^#[0-9a-f]{6}$/i.test(f.color) ? f.color : '#8c98a8',
        onchange: (e) => api.edit('Change feature colour', (d) => {
          const t = find(d);
          t.color = e.target.value;
          for (const seg of t.segments) seg.color = t.color;
        }),
      }))));

  const note = qualifierText(f, 'note') || qualifierText(f, 'product');
  if (note) box.append(h('p', { class: 'sub' }, note));
  box.append(h('p', { class: 'sub' }, `${fmt(featureLength(f, L))} bp${f.segments.length > 1 ? `, ${f.segments.length} segments` : ''} · ${featureLabel(f.type)}`));
  if (f.type === 'CDS') {
    const prot = translate(codingSequence(doc, f));
    box.append(h('p', { class: 'sub mono', style: { wordBreak: 'break-all', maxHeight: '90px', overflow: 'auto' } }, prot));
  }
  box.append(h('label', { class: 'field', style: { marginTop: '8px' } },
    h('span', {}, 'Qualifiers (one name=value per line)'),
    h('textarea', {
      style: { minHeight: '74px' },
      onchange: (e) => api.edit('Edit qualifiers', (d) => {
        const t = find(d);
        t.qualifiers = parseQualifierLines(e.target.value);
        recomputeFeature(d, t);
      }),
    }, qualifierLines(f))));
  box.append(h('div', { class: 'group' },
    h('button', { class: 'btn', onclick: () => api.setSelection({ start: featureSpan(f).start, end: featureSpan(f).end }) }, 'Select'),
    h('button', { class: 'btn', onclick: () => api.copy(api.featureSequence(f), 'bases') }, 'Copy sequence'),
    h('button', { class: 'btn', onclick: () => api.edit(`Delete ${f.name}`, (d) => { d.features = d.features.filter((x) => x.id !== f.id); }) }, 'Delete')));
  return box;
}

// ---- Detect -----------------------------------------------------------------

function detectPanel(el, { doc, state, api, others = [] }) {
  const run = state.detect?.docId === doc.id ? state.detect : null;
  const mine = api.userLibrary();
  el.append(h('h3', {}, 'Find known elements'));
  el.append(h('p', { class: 'sub' }, 'Scans both strands for promoters, tags, cleavage and recombination sites. Tags are found in any frame, whatever the codons.'));
  el.append(h('label', { class: 'field' }, h('span', {}, `Allowed mismatches: ${state.detectRate}%`),
    h('input', { id: 'detectRate', type: 'range', min: 0, max: 15, step: 1, value: state.detectRate, oninput: (e) => api.setDetectOptions({ detectRate: Number(e.target.value) }) })));
  if (others.length) {
    el.append(h('label', { class: 'sub', style: { display: 'flex', gap: '6px', alignItems: 'center', marginBottom: '8px' } },
      h('input', { type: 'checkbox', checked: state.detectUseDocs, onchange: (e) => api.setDetectOptions({ detectUseDocs: e.target.checked }) }),
      `Also use features annotated in the ${others.length} other open file${others.length === 1 ? '' : 's'}`));
  }
  el.append(h('button', { class: 'btn primary', disabled: !doc.seq.length || state.detectBusy, onclick: () => api.runDetect() }, state.detectBusy ? 'Detecting…' : 'Detect features'));

  if (run) {
    el.append(h('h3', {}, `${run.hits.length} new element${run.hits.length === 1 ? '' : 's'} found`));
    if (!run.hits.length) el.append(h('p', { class: 'empty' }, `Nothing new among ${run.libraryCount} known elements. Try allowing more mismatches.`));
    else {
      el.append(h('div', { class: 'group' },
        h('button', { class: 'btn primary', disabled: !run.hits.some((x) => x.checked), onclick: () => api.addDetected() }, `Add ${run.hits.filter((x) => x.checked).length} checked`),
        h('button', { class: 'btn', onclick: () => api.setDetectedAll(true) }, 'All'),
        h('button', { class: 'btn', onclick: () => api.setDetectedAll(false) }, 'None')));
      for (const x of run.hits) {
        el.append(h('div', { class: 'row', onclick: () => api.setSelection({ start: x.start, end: x.end }) },
          h('input', { type: 'checkbox', checked: x.checked, onclick: (e) => { e.stopPropagation(); api.toggleDetected(x.key); } }),
          h('span', { class: 'swatch', style: { background: colorFor(x.type) } }),
          h('span', { class: 'label' }, h('span', {}, x.name), h('span', { class: 'sub', style: { marginLeft: '6px' } }, x.source === 'document' ? 'from your files' : x.type)),
          h('span', { class: 'meta', title: x.note || '' }, `${fmt(x.start + 1)}–${fmt(x.end || doc.seq.length)} ${x.strand === -1 ? '←' : '→'}${x.identity < 100 ? ` ${x.identity.toFixed(0)}%` : ''}`)));
      }
    }
  }

  el.append(h('h3', {}, 'Your library'));
  el.append(h('p', { class: 'sub' }, 'Select an element in the sequence and save it here to have it found in future files. It stays in this browser.'));
  el.append(h('button', { class: 'btn', disabled: !state.selection, onclick: () => api.saveSelectionToLibrary() }, state.selection ? 'Save the selection as an element…' : 'Select bases to save an element'));
  for (const e of mine) {
    el.append(h('div', { class: 'row' },
      h('span', { class: 'label' }, e.name, h('span', { class: 'sub', style: { marginLeft: '6px' } }, `${e.seq.length} bp`)),
      h('button', { class: 'btn', style: { padding: '1px 7px' }, onclick: () => api.removeFromLibrary(e.name) }, 'Remove')));
  }
  if (!mine.length) el.append(h('p', { class: 'empty' }, 'Nothing saved yet.'));
}

// ---- Primers ----------------------------------------------------------------

function primers(el, { doc, state, api }) {
  const sel = state.selection;
  el.append(h('h3', {}, `Primers (${doc.primers.length})`));
  const tmSet = state.tm;
  const tmOpts = { na: tmSet.na / 1000, mg: tmSet.mg / 1000, oligo: tmSet.oligo * 1e-9 };
  const conditions = h('details', { style: { marginBottom: '8px' } },
    h('summary', { class: 'sub', style: { cursor: 'pointer' } }, `Melting-temperature conditions: ${tmSet.na} mM Na⁺, ${tmSet.mg} mM Mg²⁺, ${tmSet.oligo} nM primer`),
    h('div', { class: 'group', style: { marginTop: '6px' } },
      ...[['tmNa', 'Na⁺ (mM)', 'na', 0, 1000], ['tmMg', 'Mg²⁺ (mM)', 'mg', 0, 50], ['tmOligo', 'Primer (nM)', 'oligo', 1, 10000]].map(([id, label, key, min, max]) =>
        h('label', { class: 'field', style: { flex: '1 1 80px' } }, h('span', {}, label),
          h('input', { id, type: 'number', min, max, step: 'any', value: tmSet[key], oninput: (e) => { const v = Number(e.target.value); if (Number.isFinite(v) && v >= min && v <= max) api.setTm({ [key]: v }); } })))),
    h('p', { class: 'sub' }, 'Nearest-neighbour Tm with a salt correction. Use the conditions of your PCR buffer.'));
  el.append(conditions);
  if (!doc.primers.length) el.append(h('p', { class: 'empty' }, 'No primers yet.'));
  for (const p of doc.primers) {
    const sites = findBindingSites(doc.seq, doc.circular, p.seq);
    const props = primerProperties(p.seq, tmOpts);
    const hairpinBad = props.hairpin.dg <= -2;
    const selfBad = props.selfDimer.dg <= -6;
    el.append(h('div', {
      class: 'row',
      onclick: () => sites[0] && api.setSelection({ start: sites[0].start, end: sites[0].end }),
    },
      h('span', { class: 'label' }, h('span', {}, p.name),
        h('span', { class: 'sub mono', style: { display: 'block', fontSize: '11px' } }, p.seq.length > 34 ? `${p.seq.slice(0, 34)}…` : p.seq)),
      h('span', { class: 'meta' }, `${p.seq.length} nt`, h('br'), `${props.tm.toFixed(0)}°C`)));
    el.append(h('p', { class: 'sub', style: { marginLeft: '26px', color: hairpinBad || selfBad ? 'var(--mismatch)' : undefined } },
      `GC ${props.gc.toFixed(0)}% · hairpin ${props.hairpin.stem ? `ΔG ${props.hairpin.dg.toFixed(1)} (${props.hairpin.stem} bp stem)` : 'none'} · self-dimer ${props.selfDimer.length ? `ΔG ${props.selfDimer.dg.toFixed(1)} (${props.selfDimer.length} bp)` : 'none'}${hairpinBad || selfBad ? ' ⚠' : ''}`));
    if (!sites.length) el.append(h('p', { class: 'sub', style: { marginLeft: '26px' } }, 'No binding site on this sequence'));
    for (const b of sites) {
      el.append(h('p', { class: 'sub', style: { marginLeft: '26px' } },
        `${b.strand === 1 ? '→' : '←'} ${fmt(b.start + 1)}–${fmt(b.end)}, ${b.annealed.length} nt annealed, Tm ${b.tm.toFixed(0)}°C${b.tail ? `, ${b.tail.length} nt tail` : ''}`));
    }
  }

  if (doc.primers.length >= 2) {
    const a = doc.primers.find((x) => x.id === state.primerPair.a) ?? doc.primers[0];
    const b = doc.primers.find((x) => x.id === state.primerPair.b) ?? doc.primers[1];
    el.append(h('h3', {}, 'Check a pair'));
    const pick = (key, chosen) => h('select', { onchange: (e) => api.setPrimerPair({ [key]: e.target.value }) },
      doc.primers.map((x) => h('option', { value: x.id, selected: x.id === chosen.id }, x.name)));
    el.append(h('div', { class: 'group' }, pick('a', a), h('span', { class: 'sub' }, '×'), pick('b', b)));
    const cross = pairDimer(a.seq, b.seq);
    const tmA = primerProperties(a.seq, tmOpts).tm;
    const tmB = primerProperties(b.seq, tmOpts).tm;
    const bad = cross.dg <= -6;
    const warn = Math.abs(tmA - tmB) > 4;
    el.append(h('p', { class: 'sub', style: { color: bad || warn ? 'var(--mismatch)' : undefined } },
      `Tm ${tmA.toFixed(1)} and ${tmB.toFixed(1)} °C (difference ${Math.abs(tmA - tmB).toFixed(1)}${warn ? ', more than 4 ⚠' : ''}). Cross-dimer: ${cross.length ? `${cross.length} bp, ΔG ${cross.dg.toFixed(1)} kcal/mol${bad ? ' ⚠' : ''}` : 'none found'}.`));
    if (cross.structure) el.append(h('p', { class: 'sub mono', style: { fontSize: '11.5px' } }, cross.structure));
  }

  el.append(h('h3', {}, 'From the selection'));
  if (!sel) {
    el.append(h('p', { class: 'empty' }, 'Drag across the sequence, or right-click it, to use a region here.'));
  } else {
    const picked = api.selectedSequence();
    el.append(h('p', { class: 'sub' }, `${fmt(sel.start + 1)}–${fmt(sel.end)}, ${fmt(picked.length)} nt, Tm ${meltingTemp(picked).toFixed(1)}°C`));
    el.append(h('div', { class: 'group' },
      h('button', { class: 'btn', onclick: () => api.primerFromSelection(1) }, 'Forward primer'),
      h('button', { class: 'btn', onclick: () => api.primerFromSelection(-1) }, 'Reverse primer')));
  }

  el.append(h('h3', {}, 'Design a pair'));
  if (!sel) el.append(h('p', { class: 'empty' }, 'Select the region you want to amplify.'));
  else {
    const out = h('div');
    const run = () => {
      clear(out);
      try {
        const d = designPrimers(doc.seq, doc.circular, sel.start, sel.end, { targetTm: Number(el.querySelector('#tm').value) || 60 });
        for (const p of [d.forward, d.reverse]) {
          out.append(h('div', { style: { marginBottom: '8px' } },
            h('div', { class: 'sub' }, `${p.dir === 1 ? 'forward' : 'reverse'} · ${p.seq.length} nt · Tm ${p.tm.toFixed(1)}°C · GC ${p.gc.toFixed(0)}%`),
            h('div', { class: 'mono', style: { wordBreak: 'break-all', fontSize: '12px' } }, p.seq)));
        }
        out.append(h('div', { class: 'sub' }, `Product ${fmt(d.productLength)} bp`));
        out.append(h('button', {
          class: 'btn', style: { marginTop: '6px' },
          onclick: () => api.edit('Add designed primers', (doc2) => {
            doc2.primers.push({ id: newId(), name: 'forward', seq: d.forward.seq, description: 'Designed in BublaGene', extra: {} });
            doc2.primers.push({ id: newId(), name: 'reverse', seq: d.reverse.seq, description: 'Designed in BublaGene', extra: {} });
          }),
        }, 'Add both to this file'));
      } catch (e) { out.append(h('p', { class: 'empty' }, e.message)); }
    };
    el.append(h('label', { class: 'field' }, h('span', {}, 'Target melting temperature (°C)'),
      h('input', { id: 'tm', type: 'number', value: 60, min: 45, max: 75, oninput: run })));
    el.append(out);
    run();
  }

  el.append(h('h3', {}, 'PCR'));
  const fwd = h('input', { placeholder: 'forward primer 5′→3′', class: 'mono', style: { fontSize: '12px' }, value: doc.primers[0]?.seq ?? '' });
  const rev = h('input', { placeholder: 'reverse primer 5′→3′', class: 'mono', style: { fontSize: '12px' }, value: doc.primers[1]?.seq ?? '' });
  el.append(h('label', { class: 'field' }, h('span', {}, 'Forward'), fwd));
  el.append(h('label', { class: 'field' }, h('span', {}, 'Reverse'), rev));
  el.append(h('button', {
    class: 'btn primary',
    onclick: () => {
      try {
        const { doc: product, warnings } = simulatePCR(doc, fwd.value.trim().toUpperCase(), rev.value.trim().toUpperCase());
        api.openDoc(product);
        for (const w of warnings) toast(w);
      } catch (e) { toast(e.message, 'error'); }
    },
  }, 'Simulate PCR'));
  if (doc.primers.length >= 2) {
    el.append(h('p', { class: 'sub', style: { marginTop: '6px' } }, 'Tip: paste two of the primers above.'));
  }
}

// ---- Mutagenesis ------------------------------------------------------------

const AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY'.split('');

function mutatePanel(el, { doc, state, api }) {
  const m = state.mutate;
  const cds = doc.features.filter((f) => f.type === 'CDS');
  el.append(h('h3', {}, 'Site-directed mutagenesis'));
  el.append(h('p', { class: 'sub' }, 'Designs the primers for a substitution, insertion or deletion in a circular plasmid, and builds the expected mutant.'));
  if (!doc.circular) { el.append(h('p', { class: 'empty' }, 'This file is linear. Mutagenesis primers amplify the whole plasmid, so set the topology to circular in the File panel first.')); return; }

  el.append(h('label', { class: 'field' }, h('span', {}, 'What to change'),
    h('select', { id: 'mutMode', onchange: (e) => api.setMutate({ mode: e.target.value }) },
      [['aa', 'An amino acid of a CDS'], ['bases', 'The selected bases']].map(([v, t]) => h('option', { value: v, selected: m.mode === v }, t)))));

  let change = null;
  let problem = null;
  if (m.mode === 'aa') {
    if (!cds.length) problem = 'This file has no CDS feature. Add one, or edit the bases directly.';
    else {
      const feature = cds.find((f) => f.id === m.feature) ?? cds[0];
      el.append(h('label', { class: 'field' }, h('span', {}, 'Coding sequence'),
        h('select', { id: 'mutCds', onchange: (e) => api.setMutate({ feature: e.target.value }) },
          cds.map((f) => h('option', { value: f.id, selected: f.id === feature.id }, f.name)))));
      el.append(h('div', { class: 'group' },
        h('label', { class: 'field', style: { flex: '1 1 90px' } }, h('span', {}, 'Residue number'),
          h('input', { id: 'mutResidue', type: 'number', min: 1, value: m.residue, oninput: (e) => api.setMutate({ residue: Number(e.target.value) || 1 }) })),
        h('label', { class: 'field', style: { flex: '1 1 90px' } }, h('span', {}, 'Change it to'),
          h('select', { id: 'mutAa', onchange: (e) => api.setMutate({ newAa: e.target.value }) },
            [...AMINO_ACIDS, 'Stop'].map((a) => h('option', { value: a, selected: m.newAa === a }, a))))));
      try { change = aminoAcidChange(doc, feature, m.residue, m.newAa); } catch (e) { problem = e.message; }
    }
  } else {
    const sel = state.selection;
    const at = sel ?? (state.caret != null ? { start: state.caret, end: state.caret } : null);
    el.append(h('p', { class: 'sub' }, at ? (at.end > at.start ? `Replaces ${fmt(at.start + 1)}–${fmt(at.end)} (${fmt(at.end - at.start)} bp)` : `Inserts at ${fmt(at.start + 1)}`) : 'Select bases in the sequence, or click to place the cursor.'));
    el.append(h('label', { class: 'field' }, h('span', {}, 'New bases (leave empty to delete the selection)'),
      h('input', { id: 'mutBases', class: 'mono', value: m.bases, placeholder: 'e.g. GAA', oninput: (e) => api.setMutate({ bases: e.target.value }) })));
    if (at && at.end >= at.start) {
      const bases = cleanSequence(m.bases);
      if (!bases && at.end === at.start) problem = 'Nothing to change yet.';
      else change = { start: at.start, end: at.end, newBases: bases, oldBases: doc.seq.slice(at.start, at.end), label: `${at.end > at.start ? (bases ? 'sub' : 'del') : 'ins'}${at.start + 1}` };
    }
  }
  if (change) el.append(h('p', { class: 'sub mono' }, change.label && change.oldCodon ? `${change.label}: ${change.oldCodon} → ${change.newCodon} at ${fmt(change.start + 1)}–${fmt(change.end || doc.seq.length)}` : `${change.oldBases || '–'} → ${change.newBases || '–'}`));
  else if (problem) el.append(h('p', { class: 'empty' }, problem));

  el.append(h('div', { class: 'group' },
    h('label', { class: 'field', style: { flex: '2 1 140px' } }, h('span', {}, 'Primer design'),
      h('select', { id: 'mutMethod', onchange: (e) => api.setMutate({ method: e.target.value, targetTm: '' }) },
        [['overlap', 'Overlapping pair (QuikChange-style)'], ['backtoback', 'Back-to-back (Q5 SDM-style)']].map(([v, t]) => h('option', { value: v, selected: m.method === v }, t)))),
    h('label', { class: 'field', style: { flex: '1 1 70px' } }, h('span', {}, 'Target Tm °C'),
      h('input', { id: 'mutTm', type: 'number', min: 50, max: 85, value: m.targetTm, placeholder: m.method === 'overlap' ? '78' : '62', oninput: (e) => api.setMutate({ targetTm: e.target.value }) }))));
  el.append(h('button', { class: 'btn primary', disabled: !change, onclick: () => api.runMutagenesis(change) }, 'Design primers'));

  if (m.error) el.append(h('p', { class: 'empty' }, m.error));
  const r = m.result;
  if (!r) return;
  el.append(h('h3', {}, `${r.label} · ${r.kind}`));
  for (const p of [r.forward, r.reverse]) {
    el.append(h('div', { style: { marginBottom: '8px' } },
      h('div', { class: 'sub' }, `${p.name} · ${p.length} nt · Tm ${p.tm.toFixed(1)}°C · GC ${p.gc.toFixed(0)}%`),
      h('div', { class: 'mono', style: { wordBreak: 'break-all', fontSize: '12px' } }, p.seq)));
  }
  for (const n of r.notes) el.append(h('p', { class: 'sub' }, n));
  el.append(h('div', { class: 'group' },
    h('button', { class: 'btn primary', onclick: () => api.openMutant() }, `Create the mutant (${fmt(r.mutant.seq.length)} bp)`),
    h('button', { class: 'btn', onclick: () => api.addMutagenicPrimers() }, 'Add primers to this file')));
}

// ---- Enzymes ----------------------------------------------------------------

function enzymeRows(doc, names) {
  return names.map((n) => ({ name: n, sites: findSites(doc.seq, doc.circular, ENZYME_BY_NAME[n]) }));
}

function enzymes(el, { doc, state, api }) {
  el.append(h('h3', {}, 'Shown on the map'));
  el.append(h('p', { class: 'sub' }, 'Enzymes that cut this sequence are outlined; tap to show or hide.'));
  const rows = enzymeRows(doc, state.enzymeSet).sort((a, b) => a.sites.length - b.sites.length || a.name.localeCompare(b.name));
  const chips = h('div', { class: 'chips' });
  for (const r of rows) {
    chips.append(h('button', {
      class: `chip${r.sites.length ? ' cuts' : ' none'}`,
      'aria-pressed': String(state.shownEnzymes.has(r.name)),
      title: `${ENZYME_BY_NAME[r.name].site} · ${overhangLabel(ENZYME_BY_NAME[r.name])}`,
      onclick: () => api.toggleEnzyme(r.name),
    }, `${r.name}${r.sites.length ? ` ${r.sites.length}` : ''}`));
  }
  el.append(chips);
  el.append(h('div', { class: 'group', style: { marginTop: '10px' } },
    h('button', { class: 'btn', onclick: () => api.showEnzymes(rows.filter((r) => r.sites.length === 1).map((r) => r.name)) }, 'Unique cutters'),
    h('button', { class: 'btn', onclick: () => api.showEnzymes([]) }, 'None')));

  const unique = rows.filter((r) => r.sites.length === 1);
  el.append(h('h3', {}, `Unique cutters (${unique.length})`));
  const t = h('table', { class: 'data' }, h('tr', {}, h('th', {}, 'Enzyme'), h('th', {}, 'Site'), h('th', {}, 'Cuts after'), h('th', {}, 'Ends'), h('th', { title: 'Methylation that blocks or impairs this cut' }, 'Methylation')));
  for (const r of unique) {
    const meth = methylationSummary(doc, r.name, r.sites);
    t.append(h('tr', { style: { cursor: 'pointer' }, onclick: () => api.setSelection({ start: r.sites[0].site.start, end: r.sites[0].site.end }) },
      h('td', {}, r.name), h('td', { class: 'mono' }, ENZYME_BY_NAME[r.name].site),
      h('td', {}, fmt(r.sites[0].top)), h('td', { class: 'sub' }, overhangLabel(ENZYME_BY_NAME[r.name])),
      h('td', { class: 'sub', style: { color: meth.length ? 'var(--mismatch)' : undefined } },
        meth.length
          ? meth.map(({ kind, effect }) => h('span', { title: `${effect === 'blocked' ? 'Blocked' : 'Impaired by overlapping'} ${METHYLATION_LABEL[kind]} methylation: ${HOST_NOTE[kind]}`, style: { marginRight: '6px' } }, `${METHYLATION_LABEL[kind]}${effect === 'blocked' ? ' ✕' : ' ~'}`))
          : '–')));
  }
  el.append(t);
  el.append(h('p', { class: 'sub' }, '✕ blocked, ~ impaired when the methylated motif overlaps the site. Only well-established cases are listed; an enzyme not marked is not guaranteed insensitive.'));

  enzymeSetControls(el, { doc, state, api });
}

// Which enzymes are on offer: a saved set, or a hand-picked one.
const sameSet = (a, b) => a.length === b.length && [...a].sort().join() === [...b].sort().join();

function enzymeSetControls(el, { doc, state, api }) {
  el.append(h('h3', {}, 'Enzyme set'));
  const custom = api.enzymeSets();
  const all = ENZYMES.map((e) => e.name);
  const choices = [['Common', COMMON_SET], ['All', all], ...Object.entries(custom)];
  const here = choices.find(([, names]) => sameSet(names, state.enzymeSet));
  const select = h('select', {
    id: 'enzymeSetPick',
    onchange: (e) => {
      const picked = choices.find(([name]) => name === e.target.value);
      if (picked) api.setEnzymeSet([...picked[1]]);
    },
  },
  ...choices.map(([name, names]) => h('option', { value: name, selected: here?.[0] === name }, `${name} (${names.length})`)),
  here ? null : h('option', { value: '', selected: true }, `Hand-picked (${state.enzymeSet.length})`));
  el.append(h('label', { class: 'field' }, h('span', {}, 'Offered on the map and in the digest'), select));
  el.append(h('div', { class: 'group' },
    h('button', { class: 'btn', onclick: () => api.setEnzymeEdit(!state.enzymeEdit) }, state.enzymeEdit ? 'Done editing' : 'Edit this set…'),
    here && custom[here[0]] ? h('button', { class: 'btn', onclick: () => api.deleteEnzymeSet(here[0]) }, `Delete “${here[0]}”`) : null));
  if (!state.enzymeEdit) return;

  const q = (state.enzymeFilter ?? '').trim().toLowerCase();
  el.append(h('label', { class: 'field', style: { marginTop: '10px' } }, h('span', {}, 'Find an enzyme by name or site'),
    h('input', { id: 'enzymeFilter', value: state.enzymeFilter ?? '', placeholder: 'e.g. EcoRI or GAATTC', oninput: (e) => api.setEnzymeFilter(e.target.value) })));
  const inSet = new Set(state.enzymeSet);
  const matches = ENZYMES.filter((e) => !q || e.name.toLowerCase().includes(q) || e.site.toLowerCase().includes(q))
    // members first, so the current set is easy to see and prune
    .sort((a, b) => (inSet.has(b.name) - inSet.has(a.name)) || a.name.localeCompare(b.name));
  const LIMIT = 60;
  const list = h('div', { style: { maxHeight: '260px', overflow: 'auto', border: '1px solid var(--line)', borderRadius: '7px' } });
  for (const e of matches.slice(0, LIMIT)) {
    const cuts = findSites(doc.seq, doc.circular, e).length;
    list.append(h('label', { class: 'row', style: { cursor: 'pointer' } },
      h('input', { type: 'checkbox', checked: inSet.has(e.name), onchange: () => api.toggleEnzymeInSet(e.name) }),
      h('span', { class: 'label' }, e.name),
      h('span', { class: 'meta mono' }, e.site),
      h('span', { class: 'meta' }, cuts ? `${cuts}×` : '–')));
  }
  if (!matches.length) list.append(h('p', { class: 'empty' }, 'No enzyme matches.'));
  el.append(list);
  if (matches.length > LIMIT) el.append(h('p', { class: 'sub' }, `Showing ${LIMIT} of ${matches.length}; type to narrow it down.`));
  el.append(h('div', { class: 'group', style: { marginTop: '8px' } },
    h('button', {
      class: 'btn',
      onclick: () => { const name = window.prompt('Name this set:'); if (name) api.saveEnzymeSet(name); },
    }, 'Save as…'),
    h('button', { class: 'btn', onclick: () => api.setEnzymeSet([]) }, 'Clear'),
    h('button', { class: 'btn', onclick: () => api.setEnzymeSet(rowsThatCut(doc)) }, 'Only cutters of this sequence')));
}

const rowsThatCut = (doc) => ENZYMES.filter((e) => findSites(doc.seq, doc.circular, e).length).map((e) => e.name);

// ---- Digest -----------------------------------------------------------------

// Fragment sizes of one file cut with some enzymes; with none, it is left whole.
function laneSizes(d, names) {
  const cuts = names.flatMap((n) => findSites(d.seq, d.circular, ENZYME_BY_NAME[n]));
  return cuts.length ? digest(d.seq, d.circular, cuts).map((f) => f.len) : [d.seq.length];
}

// What methylation does to the enzymes picked, said in words, once per kind.
function methylationNotes(doc, names) {
  const notes = [];
  for (const n of names) {
    for (const { kind, effect } of methylationSummary(doc, n, findSites(doc.seq, doc.circular, ENZYME_BY_NAME[n]))) {
      notes.push(`${n} is ${effect === 'blocked' ? 'blocked by' : 'impaired by overlapping'} ${METHYLATION_LABEL[kind]} methylation: ${HOST_NOTE[kind]}.`);
    }
  }
  return notes;
}

function digestPanel(el, { doc, state, api, docs = [] }) {
  const chosen = state.digestWith;
  el.append(h('h3', {}, 'Digest with'));
  const chips = h('div', { class: 'chips' });
  for (const r of enzymeRows(doc, state.enzymeSet).filter((r) => r.sites.length)) {
    chips.append(h('button', {
      class: 'chip cuts', 'aria-pressed': String(chosen.includes(r.name)),
      onclick: () => api.toggleDigest(r.name),
    }, `${r.name} ${r.sites.length}`));
  }
  el.append(chips);

  const saved = state.gelLanes.map((l) => ({ l, d: docs.find((x) => x.id === l.docId) })).filter((x) => x.d);
  if (!chosen.length && !saved.length) { el.append(h('p', { class: 'empty' }, 'Pick one or more enzymes to cut with.')); return; }
  if (!chosen.length) el.append(h('p', { class: 'empty' }, 'Pick enzymes to add this file to the gel, or compare the lanes you saved.'));
  for (const note of methylationNotes(doc, chosen)) el.append(h('p', { class: 'sub' }, `⚠ ${note}`));

  let sizes = [];
  if (chosen.length) {
    const cuts = chosen.flatMap((n) => findSites(doc.seq, doc.circular, ENZYME_BY_NAME[n]));
    const frags = digest(doc.seq, doc.circular, cuts).map((f, i) => ({ ...f, i }));
    const bySize = [...frags].sort((a, b) => b.len - a.len);
    sizes = bySize.map((f) => f.len);

    el.append(h('h3', {}, `Fragments (${frags.length})`));
    const t = h('table', { class: 'data' }, h('tr', {}, h('th', {}, 'Size'), h('th', {}, 'From'), h('th', {}, 'Ends')));
    for (const f of bySize) {
      const end = (e) => (e.type === 'blunt' ? 'blunt' : e.type === 'natural' ? 'end' : `${e.type}′ ${e.oh}`);
      t.append(h('tr', {
        style: { cursor: 'pointer' },
        onclick: () => api.setSelection({ start: f.start, end: f.end === 0 && !doc.circular ? doc.seq.length : f.end }),
      },
        h('td', {}, `${fmt(f.len)} bp`),
        h('td', { class: 'sub' }, `${fmt(f.start + 1)}–${fmt(f.end || doc.seq.length)}`),
        h('td', { class: 'sub mono', style: { fontSize: '11px' } }, `${end(f.left.end)} / ${end(f.right.end)}`)));
    }
    el.append(t);
  }

  // every lane on the gel: this digest first, then the ones saved from any file
  const lanes = [
    ...(chosen.length ? [{ caption: 'digest', sizes }] : []),
    ...saved.map(({ l, d }) => ({ caption: l.label, sizes: laneSizes(d, l.enzymes) })),
  ];

  el.append(h('h3', {}, 'Gel'));
  // Redraw only the gel while the slider moves: rebuilding the whole panel would
  // replace the slider under the pointer and end the drag after one step.
  const agaroseLabel = h('span', {}, `Agarose ${state.agarose}%`);
  const gelBox = h('div', {}, gel(lanes, state.ladder, state.agarose));
  el.append(h('label', { class: 'field' }, agaroseLabel,
    h('input', {
      id: 'agarose', type: 'range', min: 0.5, max: 2, step: 0.1, value: state.agarose,
      oninput: (e) => {
        const v = Number(e.target.value);
        api.setAgarose(v, true);
        agaroseLabel.textContent = `Agarose ${v}%`;
        clear(gelBox).append(gel(lanes, state.ladder, v));
      },
    })));
  el.append(gelBox);
  el.append(h('label', { class: 'field', style: { marginTop: '8px' } }, h('span', {}, 'Ladder'),
    h('select', { onchange: (e) => api.setLadder(e.target.value) },
      Object.keys(LADDERS).map((k) => h('option', { value: k, selected: k === state.ladder }, k)))));

  el.append(h('div', { class: 'group' },
    h('button', {
      class: 'btn', disabled: !chosen.length, title: 'Keep this digest on the gel while you try another',
      onclick: () => api.addGelLane(`${doc.name}: ${chosen.join('+')}`, chosen),
    }, 'Keep this digest as a lane'),
    h('button', { class: 'btn', onclick: () => api.addGelLane(`${doc.name} uncut`, []) }, 'Add an uncut lane')));
  for (const { l, d } of saved) {
    el.append(h('div', { class: 'row', style: { cursor: 'default' } },
      h('span', { class: 'label', title: l.label }, l.label),
      h('span', { class: 'meta' }, `${laneSizes(d, l.enzymes).length} band${laneSizes(d, l.enzymes).length === 1 ? '' : 's'}`),
      h('button', { class: 'btn', style: { padding: '1px 7px' }, onclick: () => api.removeGelLane(l.id) }, 'Remove')));
  }
  if (saved.length) el.append(h('p', { class: 'sub' }, 'Lanes keep their file and enzymes, so you can compare digests of different files side by side.'));
}

// lanes: [{ caption, sizes }]
function gel(lanes, ladderName, agarose) {
  const box = h('div', { class: 'gel' });
  const H = 300;
  const lane = (caption, bands, isLadder) => {
    const l = h('div', { class: `lane${isLadder ? ' ladder' : ''}` }, h('div', { class: 'caption', title: caption }, caption));
    const counts = new Map();
    for (const bp of bands) counts.set(bp, (counts.get(bp) ?? 0) + 1);
    for (const [bp, n] of counts) {
      const y = 16 + migration(bp, agarose) * (H - 30);
      l.append(h('div', {
        class: 'band',
        style: { top: `${y}px`, opacity: String(Math.min(1, 0.45 + Math.log10(bp) / 5)), height: `${Math.min(6, 2 + n)}px` },
      }, h('span', {}, `${fmt(bp)}${n > 1 ? ` ×${n}` : ''}`)));
    }
    return l;
  };
  box.append(lane(ladderName, LADDERS[ladderName] ?? [], true));
  for (const l of lanes) box.append(lane(l.caption, l.sizes, false));
  return box;
}

// ---- Alignment -------------------------------------------------------------

function alignmentPanel(el, { doc, state, api, others = [], docs = [] }) {
  if (docs.length >= 2) msaSection(el, { state, api, docs });
  el.append(h('h3', {}, 'Align to another document'));
  if (!others.length) {
    el.append(h('p', { class: 'empty' }, 'Open a second document to align against this one.'));
    return;
  }
  const target = others.find((d) => d.id === state.alignTarget) ?? others[0];
  el.append(h('label', { class: 'field' }, h('span', {}, 'Reference'),
    h('select', { onchange: (e) => api.setAlignTarget(e.target.value) },
      others.map((d) => h('option', { value: d.id, selected: d.id === target.id }, `${d.name} (${fmt(d.seq.length)} bp)`)))));
  el.append(h('div', { class: 'group' },
    h('label', { class: 'sub', style: { display: 'flex', gap: '5px', alignItems: 'center' } },
      h('input', { type: 'checkbox', checked: state.alignLocal, onchange: (e) => api.setAlignLocal(e.target.checked) }), 'Local'),
    h('label', { class: 'sub', style: { display: 'flex', gap: '5px', alignItems: 'center' } },
      'Band', h('input', { type: 'number', min: 0, value: state.alignBand ?? 0, style: { width: '62px' }, onchange: (e) => api.setAlignBand(Number(e.target.value)) }))));
  el.append(h('div', { class: 'group', style: { marginTop: '8px' } },
    h('button', { class: 'btn primary', disabled: state.alignBusy, onclick: () => api.runAlignment() }, state.alignBusy ? 'Aligning…' : 'Align'),
    h('button', { class: 'btn', disabled: !state.alignResult || state.alignResult.docId !== doc.id, onclick: () => api.alignmentVariants() }, 'Make variant features')));

  const result = state.alignResult?.docId === doc.id ? state.alignResult : null;
  if (!result) { el.append(h('p', { class: 'empty' }, 'No alignment yet.')); return; }
  if (result.error) { el.append(h('p', { class: 'empty' }, result.error)); return; }
  const r = result.result;
  el.append(h('p', { class: 'sub' }, `score ${fmt(r.score)} · identity ${r.identity.toFixed(1)}% · ${fmt(r.matches)} matches · ${fmt(r.mismatches)} mismatches · ${fmt(r.gaps)} gap columns`));
  el.append(dotPlotSvg(doc.seq, target.seq, 260));
  el.append(h('h3', {}, 'Aligned'));
  el.append(alignmentText(r));
  if (r.alignedA.length > 2000) el.append(h('p', { class: 'sub' }, 'Showing the first 2,000 alignment columns.'));
}

function dotPlotSvg(a, b, size) {
  const step = Math.max(1, Math.floor(Math.max(a.length, b.length) / 180));
  const plot = dotPlot(a, b, { window: 10, step, threshold: 0.7, maxPoints: 6000 });
  const root = s('svg', { viewBox: `0 0 ${size} ${size}`, width: size, height: size, role: 'img', 'aria-label': 'Dot plot' });
  root.append(s('rect', { x: 0, y: 0, width: size, height: size, fill: 'var(--paper)', stroke: 'var(--line)' }));
  const sx = (p) => (p / Math.max(1, plot.width)) * size;
  const sy = (p) => size - (p / Math.max(1, plot.height)) * size;
  for (const [i, j, score] of plot.points) {
    root.append(s('rect', { x: sx(i), y: sy(j), width: 1.4, height: 1.4, fill: score > 0.95 ? 'var(--accent)' : 'var(--muted)', opacity: 0.7 }));
  }
  if (plot.truncated) root.append(s('text', { x: 6, y: size - 6, 'font-size': 9, fill: 'var(--muted)' }, 'points capped'));
  return root;
}

function alignmentText(r) {
  const box = h('div', { class: 'mono', style: { fontSize: '11px', lineHeight: '1.35', whiteSpace: 'pre-wrap', wordBreak: 'break-all', maxHeight: '260px', overflow: 'auto' } });
  const limit = Math.min(r.alignedA.length, 2000);
  const addRow = (text, other, side) => {
    let run = { cls: 'match', text: '' };
    const flush = () => { if (run.text) box.append(h('span', { style: { color: run.cls === 'match' ? 'var(--muted)' : run.cls === 'gap' ? 'var(--danger)' : 'var(--mismatch)' } }, run.text)); run = { cls: run.cls, text: '' }; };
    for (let i = 0; i < limit; i++) {
      const c = text[i];
      const cls = c === '-' ? 'gap' : c === other[i] ? 'match' : 'mismatch';
      if (cls !== run.cls) { flush(); run.cls = cls; }
      run.text += c;
    }
    flush();
    box.append('\n');
  };
  addRow(r.alignedA, r.alignedB, 'target');
  addRow(r.alignedB, r.alignedA, 'query');
  return box;
}

// ---- CRISPR ----------------------------------------------------------------

function crisprPanel(el, { doc, state, api }) {
  const current = state.crispr?.docId === doc.id ? state.crispr : null;
  el.append(h('h3', {}, 'Find guides'));
  el.append(h('label', { class: 'field' }, h('span', {}, 'Enzyme / PAM'),
    h('select', { id: 'crisprEnzyme' }, ['SpCas9', 'SaCas9', 'Cas12a', 'custom'].map((x) => h('option', { value: x, selected: (current?.enzyme ?? 'SpCas9') === x }, x)))));
  el.append(h('div', { class: 'group' },
    h('label', { class: 'field', style: { flex: '1 1 70px' } }, h('span', {}, 'Guide nt'),
      h('input', { id: 'crisprLen', type: 'number', min: 10, max: 30, value: current?.guideLength ?? 20 })),
    h('label', { class: 'field', style: { flex: '1 1 70px' } }, h('span', {}, 'Min GC'),
      h('input', { id: 'crisprMinGc', type: 'number', min: 0, max: 100, value: current?.minGC ?? 20 })),
    h('label', { class: 'field', style: { flex: '1 1 70px' } }, h('span', {}, 'Max GC'),
      h('input', { id: 'crisprMaxGc', type: 'number', min: 0, max: 100, value: current?.maxGC ?? 80 }))));
  el.append(h('div', { class: 'group' },
    h('label', { class: 'sub', style: { display: 'flex', gap: '5px', alignItems: 'center' } },
      h('input', { id: 'crisprOff', type: 'checkbox', checked: current?.offTargets ?? false }), 'Off-target scan (heuristic)'),
    h('button', { class: 'btn primary', onclick: () => api.runCrispr({
      enzyme: el.querySelector('#crisprEnzyme').value,
      guideLength: Number(el.querySelector('#crisprLen').value) || 20,
      minGC: Number(el.querySelector('#crisprMinGc').value),
      maxGC: Number(el.querySelector('#crisprMaxGc').value),
      offTargets: el.querySelector('#crisprOff').checked,
    }) }, 'Find guides')));

  if (!current) { el.append(h('p', { class: 'empty' }, 'No guide search yet.')); return; }
  if (current.error) { el.append(h('p', { class: 'empty' }, current.error)); return; }
  el.append(h('h3', {}, `${current.guides.length} guide(s)`));
  if (current.guides.length) {
    el.append(h('div', { class: 'group' },
      h('button', { class: 'btn', onclick: () => api.addCrisprFeatures() }, 'Add guide features'),
      h('button', { class: 'btn', onclick: () => api.clearCrispr() }, 'Clear results')));
  }
  for (const g of current.guides.slice(0, 80)) {
    const row = h('div', {
      class: 'row', title: `${g.strand === 1 ? '+' : '−'} ${g.start + 1}–${g.end}, PAM ${g.pam}`,
      onclick: () => api.setSelection({ start: g.start, end: g.end }),
    },
      h('span', { class: 'label mono', style: { fontSize: '12px' } }, `${g.strand === 1 ? '→' : '←'} ${g.sequence}`),
      h('span', { class: 'meta' }, `${fmt(g.start + 1)} · ${g.gc.toFixed(0)}%`));
    el.append(row);
    el.append(h('div', { class: 'group', style: { margin: '0 0 6px 20px' } },
      h('button', { class: 'btn', onclick: () => api.simulateNhej(g) }, 'NHEJ'),
      h('button', { class: 'btn', onclick: () => {
        const donor = window.prompt('HDR donor sequence (5′→3′):');
        if (donor) api.simulateHdr(g, cleanSequence(donor));
      } }, 'HDR…'),
      h('button', { class: 'btn', onclick: () => api.crisprPrimers(g) }, 'Primers'),
      g.offTargets ? h('span', { class: 'sub' }, `${g.offTargets.length} near-match(es)`) : null));
  }
  el.append(h('p', { class: 'sub' }, 'Off-target counts are a seed/mismatch heuristic, not a cleavage prediction.'));
}

// ---- Structure -------------------------------------------------------------

function structurePanel(el, { doc, state, api }) {
  const result = state.fold?.docId === doc.id ? state.fold : null;
  const sel = state.selection;
  el.append(h('h3', {}, 'RNA secondary structure'));
  el.append(h('p', { class: 'sub' }, 'Bounded Turner subset: WC + G·U, stacking, hairpin and simple interior loops.'));
  el.append(h('div', { class: 'group' },
    h('label', { class: 'field', style: { flex: '1 1 90px' } }, h('span', {}, 'Max length'),
      h('input', { id: 'foldMax', type: 'number', min: 20, max: 400, value: result?.maxLength ?? 250 })),
    h('button', { class: 'btn primary', disabled: state.fold?.busy, onclick: () => api.runFold({
      start: sel?.start ?? null, end: sel?.end ?? null,
      maxLength: Number(el.querySelector('#foldMax').value) || 250,
    }) }, state.fold?.busy ? 'Folding…' : (sel ? 'Fold selection' : 'Fold sequence'))));
  if (!result) { el.append(h('p', { class: 'empty' }, 'No fold yet.')); return; }
  if (result.busy) { el.append(h('p', { class: 'empty' }, 'Folding…')); return; }
  if (result.error) { el.append(h('p', { class: 'empty' }, result.error)); return; }
  const r = result.result;
  el.append(h('p', { class: 'sub' }, `ΔG ${r.energy} kcal/mol · ${r.length} nt${result.worker ? ' · computed in a worker' : ''}`));
  el.append(h('p', { class: 'mono', style: { fontSize: '11px', whiteSpace: 'pre-wrap', wordBreak: 'break-all' } }, r.structure));
  const svgBox = h('div', { style: { overflow: 'auto', border: '1px solid var(--line)', borderRadius: '7px', background: 'var(--paper)' } });
  svgBox.innerHTML = structureToSvg(r.sequence, r.structure, { width: 640, perRow: 50, rowHeight: 110 });
  el.append(svgBox);
  el.append(h('div', { class: 'group', style: { marginTop: '8px' } },
    h('button', { class: 'btn', onclick: () => api.exportStructure() }, 'Export SVG'),
    h('button', { class: 'btn', onclick: () => api.copy(r.structure, 'structure characters') }, 'Copy dot-bracket')));
}

// ---- File -------------------------------------------------------------------

// What else a .dna file carries, in the reader's terms rather than the format's.
// Anything not listed here is either shown in BublaGene already or a search
// cache the format rebuilds, and neither is worth mentioning.
const CARRIED = {
  0x07: ['the cloning history', 'cloning histories'],
  0x0b: ['an earlier version of the sequence', 'earlier versions of the sequence'],
  0x0c: ['extra annotations', 'sets of extra annotations'],
  0x0d: ['display settings', 'sets of display settings'],
  0x0e: ['custom enzyme sets', 'custom enzyme sets'],
  0x0f: ['unsynced case changes', 'unsynced case changes'],
  // 0x10-0x1b are the sequencing data, drawn in the chromatogram and sequence
  // views rather than listed here as unviewable
  0x13: ['marked uracil positions', 'sets of marked uracil positions'],
  0x14: ['custom strand colours', 'sets of custom strand colours'],
  0x16: ['custom numbering', 'sets of custom numbering'],
  0x17: ['an embedded file', 'embedded files'],
  0x19: ['an alignment', 'alignments'],
  0x1a: ['an edited alignment', 'edited alignments'],
  // 0x1b, the aligned reads, is shown in the sequence view, so it is not listed
  // here among the things BublaGene cannot display
  0x1c: ['which enzymes are shown', 'records of which enzymes are shown'],
  0x1d: ['an edit to an earlier version', 'edits to earlier versions'],
  0x22: ['a secondary-structure prediction', 'secondary-structure predictions'],
  0x23: ['saved CRISPR searches', 'sets of saved CRISPR searches'],
};

function fileContents(el, doc) {
  const segs = describeSegments(doc);
  if (!segs.length) return;

  const counts = new Map();
  for (const s of segs) {
    if (CARRIED[s.type]) counts.set(s.type, (counts.get(s.type) ?? 0) + 1);
  }

  el.append(h('h3', {}, 'Also in this file'));
  if (!counts.size) {
    el.append(h('p', { class: 'sub' }, 'Nothing beyond what you can see here.'));
  } else {
    el.append(h('p', { class: 'sub' }, 'BublaGene cannot show these, but keeps them exactly as they are when you save:'));
    const list = h('ul', { class: 'sub', style: { margin: '6px 0 0', paddingLeft: '18px' } });
    for (const [type, n] of [...counts].sort((a, b) => b[1] - a[1])) {
      const [one, many] = CARRIED[type];
      list.append(h('li', {}, n > 1 ? `${fmt(n)} ${many}` : one));
    }
    el.append(list);
  }

  // the format-level view, for anyone who wants it
  const details = h('details', { style: { marginTop: '10px' } });
  details.append(h('summary', { class: 'sub', style: { cursor: 'pointer' } }, `All ${segs.length} parts of the file`));
  const st = h('table', { class: 'data', style: { marginTop: '6px' } });
  for (const sgm of segs) {
    st.append(h('tr', {},
      h('td', { class: 'mono sub' }, `0x${sgm.type.toString(16).padStart(2, '0')}`),
      h('td', {}, sgm.name),
      h('td', {}, `${fmt(sgm.bytes)} B`)));
  }
  details.append(st);
  el.append(details);
}

function file(el, { doc, state, api }) {
  const L = doc.seq.length;
  el.append(h('h3', {}, 'Save and projects'));
  el.append(h('p', { class: 'sub' }, 'On a browser with the File System Access API this asks where to save; otherwise it downloads the file.'));
  el.append(h('div', { class: 'group' },
    h('button', { class: 'btn', onclick: () => api.save('dna') }, 'Save .dna'),
    h('button', { class: 'btn', onclick: () => api.save('gb') }, 'Save GenBank'),
    h('button', { class: 'btn', onclick: () => api.save('fasta') }, 'Save FASTA')));
  el.append(h('div', { class: 'group', style: { marginTop: '6px' } },
    h('button', { class: 'btn', onclick: () => api.saveProject() }, 'Save project…'),
    h('button', { class: 'btn', onclick: () => api.openProject() }, 'Open project…')));
  el.append(h('p', { class: 'sub' }, 'Projects keep every open document, its features, primers, cloning and assembly records, alignments and settings.'));
  el.append(h('label', { class: 'field' }, h('span', {}, 'Name'),
    h('input', { value: doc.name, onchange: (e) => api.edit('Rename', (d) => { d.name = e.target.value; }) })));
  el.append(h('label', { class: 'field' }, h('span', {}, 'Topology'),
    h('select', {
      onchange: (e) => api.edit('Change topology', (d) => { d.circular = e.target.value === 'circular'; d.flags = d.circular ? d.flags | 1 : d.flags & ~1; }),
    }, [['circular', 'circular'], ['linear', 'linear']].map(([v, t]) => h('option', { value: v, selected: doc.circular === (v === 'circular') }, t)))));
  el.append(h('label', { class: 'field' }, h('span', {}, 'Description'),
    h('textarea', { onchange: (e) => api.edit('Edit description', (d) => { d.notes.description = e.target.value; }) }, doc.notes.description ?? '')));

  el.append(h('h3', {}, 'Composition'));
  const counts = { A: 0, C: 0, G: 0, T: 0, other: 0 };
  for (const ch of doc.seq) counts[ch] === undefined ? counts.other++ : counts[ch]++;
  const t = h('table', { class: 'data' });
  t.append(h('tr', {}, h('td', {}, 'Length'), h('td', {}, `${fmt(L)} bp`)));
  t.append(h('tr', {}, h('td', {}, 'GC'), h('td', {}, `${gcContent(doc.seq).toFixed(1)}%`)));
  for (const b of ['A', 'C', 'G', 'T']) {
    const pct = L ? ` (${((100 * counts[b]) / L).toFixed(1)}%)` : '';
    t.append(h('tr', {}, h('td', {}, b), h('td', {}, `${fmt(counts[b])}${pct}`)));
  }
  if (counts.other) t.append(h('tr', {}, h('td', {}, 'other'), h('td', {}, fmt(counts.other))));
  el.append(t);

  if (doc.notes.organism) el.append(h('p', { class: 'sub' }, `Organism: ${doc.notes.organism}`));

  fileContents(el, doc);
  for (const w of doc.source?.warnings ?? []) el.append(h('p', { class: 'sub' }, w));

  if (doc.history?.length) {
    el.append(h('h3', {}, 'What you have done here'));
    for (const ev of [...doc.history].reverse().slice(0, 25)) {
      el.append(h('p', { class: 'sub' }, `${ev.detail || ev.action}`));
    }
  }
}
