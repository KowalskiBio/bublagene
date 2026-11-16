// The inspector panels on the right.

import { h, clear, fmt, toast } from './dom.js';
import { ENZYMES, ENZYME_BY_NAME, COMMON_SET, findSites, digest, overhangLabel, migration, LADDERS } from '../enzymes.js';
import { findBindingSites, designPrimers, simulatePCR } from '../primers.js';
import { featureSpan, featureLength, codingSequence, colorFor, FEATURE_TYPES, newId, qualifierText, addHistory } from '../model.js';
import { gcContent, meltingTemp, translate, revcomp, subseq, rangeLength, findOrfs } from '../seq.js';
import { describeSegments } from '../formats/dna.js';

export const PANELS = ['Features', 'Primers', 'Enzymes', 'Digest', 'File'];

export function renderPanel(name, ctx) {
  const el = h('div');
  ({ Features: features, Primers: primers, Enzymes: enzymes, Digest: digestPanel, File: file })[name](el, ctx);
  return el;
}

// ---- Features ---------------------------------------------------------------

function features(el, { doc, state, api }) {
  const L = doc.seq.length;
  const sel = state.selection;
  el.append(h('h3', {}, `Features (${doc.features.length})`));
  if (!doc.features.length) el.append(h('p', { class: 'empty' }, 'No features yet. Select part of the sequence and add one.'));

  const sorted = [...doc.features].sort((a, b) => featureSpan(a).start - featureSpan(b).start);
  for (const f of sorted) {
    const span = featureSpan(f);
    const row = h('div', {
      class: 'row', 'aria-current': String(state.highlight === f.id),
      onclick: () => api.focusFeature(f),
    },
      h('span', { class: 'swatch', style: { background: f.color } }),
      h('span', { class: 'label' },
        h('span', {}, f.name),
        h('span', { class: 'sub', style: { marginLeft: '6px' } }, f.type)),
      h('span', { class: 'meta' }, `${fmt(span.start + 1)}–${fmt(span.end)} ${f.strand === -1 ? '←' : f.strand === 1 ? '→' : ''}`));
    el.append(row);
    if (state.highlight === f.id) el.append(featureDetail(f, { doc, state, api }));
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

  const orfs = findOrfs(doc.seq, doc.circular, 100);
  el.append(h('h3', {}, `Open reading frames (${orfs.length})`));
  if (!orfs.length) el.append(h('p', { class: 'empty' }, 'None of at least 100 codons.'));
  for (const o of orfs.slice(0, 12)) {
    el.append(h('div', {
      class: 'row', onclick: () => api.setSelection({ start: o.start, end: o.end }),
    },
      h('span', { class: 'label mono' }, `${o.strand === 1 ? '→' : '←'} ${fmt(o.start + 1)}–${fmt(o.end)}`),
      h('span', { class: 'meta' }, `${fmt(o.codons)} codons`)));
  }
}

function featureDetail(f, { doc, api }) {
  const L = doc.seq.length;
  const box = h('div', { style: { padding: '4px 8px 12px', borderLeft: '2px solid var(--accent)', margin: '2px 0 8px 6px' } });
  box.append(h('label', { class: 'field' }, h('span', {}, 'Name'),
    h('input', { value: f.name, onchange: (e) => api.edit('Rename feature', (d) => { d.features.find((x) => x.id === f.id).name = e.target.value; }) })));
  box.append(h('label', { class: 'field' }, h('span', {}, 'Type'),
    h('select', {
      onchange: (e) => api.edit('Change feature type', (d) => {
        const t = d.features.find((x) => x.id === f.id);
        t.type = e.target.value;
        t.color = colorFor(e.target.value);
        for (const s of t.segments) s.color = t.color;
      }),
    }, FEATURE_TYPES.map((t) => h('option', { value: t, selected: t === f.type }, t)))));

  const note = qualifierText(f, 'note') || qualifierText(f, 'product');
  if (note) box.append(h('p', { class: 'sub' }, note));
  box.append(h('p', { class: 'sub' }, `${fmt(featureLength(f, L))} bp${f.segments.length > 1 ? `, ${f.segments.length} segments` : ''}`));
  if (f.type === 'CDS') {
    const prot = translate(codingSequence(doc, f));
    box.append(h('p', { class: 'sub mono', style: { wordBreak: 'break-all', maxHeight: '90px', overflow: 'auto' } }, prot));
  }
  box.append(h('div', { class: 'group' },
    h('button', { class: 'btn', onclick: () => api.setSelection({ start: featureSpan(f).start, end: featureSpan(f).end }) }, 'Select'),
    h('button', { class: 'btn', onclick: () => api.edit(`Delete ${f.name}`, (d) => { d.features = d.features.filter((x) => x.id !== f.id); }) }, 'Delete')));
  return box;
}

// ---- Primers ----------------------------------------------------------------

function primers(el, { doc, state, api }) {
  const sel = state.selection;
  el.append(h('h3', {}, `Primers (${doc.primers.length})`));
  if (!doc.primers.length) el.append(h('p', { class: 'empty' }, 'No primers yet.'));
  for (const p of doc.primers) {
    const sites = findBindingSites(doc.seq, doc.circular, p.seq);
    el.append(h('div', {
      class: 'row',
      onclick: () => sites[0] && api.setSelection({ start: sites[0].start, end: sites[0].end }),
    },
      h('span', { class: 'label' }, h('span', {}, p.name),
        h('span', { class: 'sub mono', style: { display: 'block', fontSize: '11px' } }, p.seq.length > 34 ? `${p.seq.slice(0, 34)}…` : p.seq)),
      h('span', { class: 'meta' }, `${p.seq.length} nt`, h('br'), `${meltingTemp(p.seq).toFixed(0)}°C`)));
    if (!sites.length) el.append(h('p', { class: 'sub', style: { marginLeft: '26px' } }, 'No binding site on this sequence'));
    for (const b of sites) {
      el.append(h('p', { class: 'sub', style: { marginLeft: '26px' } },
        `${b.strand === 1 ? '→' : '←'} ${fmt(b.start + 1)}–${fmt(b.end)}, ${b.annealed.length} nt annealed, Tm ${b.tm.toFixed(0)}°C${b.tail ? `, ${b.tail.length} nt tail` : ''}`));
    }
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
  const t = h('table', { class: 'data' }, h('tr', {}, h('th', {}, 'Enzyme'), h('th', {}, 'Site'), h('th', {}, 'Cuts after'), h('th', {}, 'Ends')));
  for (const r of unique) {
    t.append(h('tr', { style: { cursor: 'pointer' }, onclick: () => api.setSelection({ start: r.sites[0].site.start, end: r.sites[0].site.end }) },
      h('td', {}, r.name), h('td', { class: 'mono' }, ENZYME_BY_NAME[r.name].site),
      h('td', {}, fmt(r.sites[0].top)), h('td', { class: 'sub' }, overhangLabel(ENZYME_BY_NAME[r.name]))));
  }
  el.append(t);

  el.append(h('h3', {}, 'Enzyme set'));
  el.append(h('div', { class: 'group' },
    h('button', { class: 'btn', onclick: () => api.setEnzymeSet(COMMON_SET) }, `Common (${COMMON_SET.length})`),
    h('button', { class: 'btn', onclick: () => api.setEnzymeSet(ENZYMES.map((e) => e.name)) }, `All (${ENZYMES.length})`)));
}

// ---- Digest -----------------------------------------------------------------

function digestPanel(el, { doc, state, api }) {
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
  if (!chosen.length) { el.append(h('p', { class: 'empty' }, 'Pick one or more enzymes to cut with.')); return; }

  const cuts = chosen.flatMap((n) => findSites(doc.seq, doc.circular, ENZYME_BY_NAME[n]));
  const frags = digest(doc.seq, doc.circular, cuts).map((f, i) => ({ ...f, i }));
  const bySize = [...frags].sort((a, b) => b.len - a.len);

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

  el.append(h('h3', {}, 'Gel'));
  const pct = state.agarose;
  el.append(h('label', { class: 'field' }, h('span', {}, `Agarose ${pct}%`),
    h('input', { type: 'range', min: 0.5, max: 2, step: 0.1, value: pct, oninput: (e) => api.setAgarose(Number(e.target.value)) })));
  el.append(gel(bySize.map((f) => f.len), state.ladder, pct));
  el.append(h('label', { class: 'field', style: { marginTop: '8px' } }, h('span', {}, 'Ladder'),
    h('select', { onchange: (e) => api.setLadder(e.target.value) },
      Object.keys(LADDERS).map((k) => h('option', { value: k, selected: k === state.ladder }, k)))));
}

function gel(sizes, ladderName, agarose) {
  const box = h('div', { class: 'gel' });
  const H = 300;
  const lane = (caption, bands, isLadder) => {
    const l = h('div', { class: `lane${isLadder ? ' ladder' : ''}` }, h('div', { class: 'caption' }, caption));
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
  box.append(lane('digest', sizes, false));
  return box;
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
  for (const b of ['A', 'C', 'G', 'T']) t.append(h('tr', {}, h('td', {}, b), h('td', {}, `${fmt(counts[b])} (${((100 * counts[b]) / L).toFixed(1)}%)`)));
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
