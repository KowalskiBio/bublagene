// The Cloning panel. Everything here is presentation: choosing which fragments
// to join and showing what comes back. The biology is in js/cloning.js.

import { h, clear, fmt, toast } from './dom.js';
import { ENZYME_BY_NAME, COMMON_SET, findSites, digest } from '../enzymes.js';
import { restrictionClone, gibsonAssemble, infusionAssemble, checkLigation, pieceFromFragment, describeEnd } from '../cloning.js';
import { goldenGateAssemble, goldenGatePieces } from '../assembly.js';
import { findRecombinationSites, recombine, topoClone } from '../recombination.js';
import { ligationPreview } from './ligation-preview.js';
import { historyTree } from './historytree.js';

const MODES = [
  ['restriction', 'Restriction + ligation'], ['goldengate', 'Golden Gate'], ['gibson', 'Gibson'],
  ['infusion', 'In-Fusion'], ['recombination', 'Recombination (att, lox, FRT)'], ['topo', 'TOPO / TA / blunt'],
];
const TYPE_IIS = ['BsaI', 'BsmBI', 'BbsI', 'SapI', 'Esp3I'].filter((n) => ENZYME_BY_NAME[n]);

// Enzymes worth offering for a cut: those that cut this sequence a handful of
// times, rarest first, since a cloning site is usually a unique cutter.
function cuttersFor(doc, limit = 18) {
  const rows = [];
  for (const name of COMMON_SET) {
    const n = findSites(doc.seq, doc.circular, ENZYME_BY_NAME[name]).length;
    if (n >= 1 && n <= 4) rows.push({ name, n });
  }
  return rows.sort((a, b) => a.n - b.n || a.name.localeCompare(b.name)).slice(0, limit);
}

function fragmentsOf(doc, names) {
  if (!names.length) return [];
  const cuts = names.flatMap((n) => findSites(doc.seq, doc.circular, ENZYME_BY_NAME[n]));
  if (!cuts.length) return [];
  return digest(doc.seq, doc.circular, cuts).map((f, i) => ({ ...f, index: i }));
}

const fragLabel = (f, doc) => {
  const ends = [f.left.cut?.enzyme?.name, f.right.cut?.enzyme?.name].filter(Boolean);
  return `${fmt(f.len)} bp · ${fmt(f.start + 1)}–${fmt(f.end || doc.seq.length)}${ends.length ? ` · ${[...new Set(ends)].join('/')}` : ''}`;
};

// One side of a restriction cloning: which document, cut with what, which piece.
function sidePicker(el, { label, side, state, api, docs }) {
  const sel = state.cloning[side];
  // until a choice is made, start the two sides on different files
  const fallback = side === 'insert'
    ? docs.find((d) => d.id !== state.cloning.vector.docId) ?? docs[docs.length - 1]
    : docs[0];
  const doc = docs.find((d) => d.id === sel.docId) ?? fallback;

  el.append(h('h3', {}, label));
  el.append(h('label', { class: 'field' },
    h('select', {
      onchange: (e) => api.setCloning({ [side]: { ...sel, docId: e.target.value, enzymes: [], fragment: null } }),
    }, docs.map((d) => h('option', { value: d.id, selected: d.id === doc?.id }, `${d.name} (${fmt(d.seq.length)} bp)`)))));

  if (!doc) return null;

  const chips = h('div', { class: 'chips' });
  for (const { name, n } of cuttersFor(doc)) {
    chips.append(h('button', {
      class: 'chip cuts', 'aria-pressed': String(sel.enzymes.includes(name)),
      title: `${ENZYME_BY_NAME[name].site}, cuts ${n}×`,
      onclick: () => {
        const next = sel.enzymes.includes(name) ? sel.enzymes.filter((x) => x !== name) : [...sel.enzymes, name];
        api.setCloning({ [side]: { ...sel, enzymes: next, fragment: null } });
      },
    }, `${name}${n > 1 ? ` ${n}` : ''}`));
  }
  el.append(chips);

  const frags = fragmentsOf(doc, sel.enzymes);
  if (!sel.enzymes.length) {
    el.append(h('p', { class: 'empty' }, 'Choose the enzymes to cut with.'));
    return null;
  }
  if (!frags.length) {
    el.append(h('p', { class: 'empty' }, 'Those enzymes do not cut this sequence.'));
    return null;
  }

  // Default to the largest piece cut at both ends: on a vector that is the
  // backbone, on an insert the insert itself. A piece with an uncut end is the
  // leftover from a linear sequence and cannot be ligated.
  const usable = frags.filter((f) => f.left.end.type !== 'natural' && f.right.end.type !== 'natural');
  const best = [...(usable.length ? usable : frags)].sort((a, b) => b.len - a.len)[0];
  const chosen = frags.find((f) => f.index === sel.fragment) ?? best;
  el.append(h('label', { class: 'field' },
    h('select', {
      onchange: (e) => api.setCloning({ [side]: { ...sel, fragment: Number(e.target.value) } }),
    }, frags.map((f) => {
      const uncut = f.left.end.type === 'natural' || f.right.end.type === 'natural';
      return h('option', { value: f.index, selected: f.index === chosen.index }, `${fragLabel(f, doc)}${uncut ? ' · has an uncut end' : ''}`);
    }))));
  el.append(h('p', { class: 'sub' }, `Ends: ${describeEnd({ type: chosen.left.end.type, oh: chosen.left.end.oh })} / ${describeEnd({ type: chosen.right.end.type, oh: chosen.right.end.oh })}`));
  return { doc, fragment: chosen };
}

export function cloningPanel(el, { doc, state, api, docs }) {
  const c = state.cloning;
  if (!docs.length) {
    el.append(h('p', { class: 'empty' }, 'Open a vector and an insert first.'));
    return;
  }
  if (!c.vector.docId && docs.length) c.vector = { ...c.vector, docId: docs[0].id };
  if (docs.length < 2 && c.mode === 'restriction') {
    el.append(h('p', { class: 'sub' }, 'Open the vector and the insert as two files to clone between them. You can open the same file twice.'));
  }

  el.append(h('label', { class: 'field' }, h('span', {}, 'Method'),
    h('select', { id: 'cloneMode', onchange: (e) => api.setCloning({ mode: e.target.value }) },
      MODES.map(([v, t]) => h('option', { value: v, selected: c.mode === v }, t)))));

  const done = () => { if (doc?.cloning) el.append(historyTree(doc, docs, (id) => api.focusDoc(id))); };
  if (c.mode === 'gibson') { gibsonSection(el, { state, api, docs }, { title: 'Gibson', assemble: gibsonAssemble, overlapKey: 'minOverlap' }); return done(); }
  if (c.mode === 'infusion') { gibsonSection(el, { state, api, docs }, { title: 'In-Fusion', assemble: infusionAssemble, overlapKey: 'infusionOverlap' }); return done(); }
  if (c.mode === 'goldengate') { goldenGateSection(el, { state, api, docs }); return done(); }
  if (c.mode === 'recombination') { recombinationSection(el, { state, api, docs }); return done(); }
  if (c.mode === 'topo') { topoSection(el, { doc, state, api, docs }); return done(); }

  const vector = sidePicker(el, { label: 'Vector', side: 'vector', state, api, docs });
  const insert = sidePicker(el, { label: 'Insert', side: 'insert', state, api, docs });

  el.append(h('h3', {}, 'Product'));
  el.append(h('label', { class: 'field' }, h('span', {}, 'Shape'),
    h('select', { onchange: (e) => api.setCloning({ circular: e.target.value === 'circular' }) },
      [['circular', 'circular plasmid'], ['linear', 'linear fragment']].map(([v, t]) => h('option', { value: v, selected: c.circular === (v === 'circular') }, t)))));
  el.append(h('label', { class: 'field' }, h('span', {}, 'Insert orientation'),
    h('select', { onchange: (e) => api.setCloning({ orientation: e.target.value }) },
      [['auto', 'whichever fits'], ['forward', 'forward'], ['reverse', 'reverse']].map(([v, t]) => h('option', { value: v, selected: c.orientation === v }, t)))));

  if (!vector || !insert) return done();

  // tell the user whether it will ligate before they commit to it
  let check = null;
  try {
    check = checkLigation({
      backbone: { doc: vector.doc, fragment: vector.fragment },
      insert: { doc: insert.doc, fragment: insert.fragment },
      circular: c.circular,
    });
  } catch (e) {
    el.append(h('p', { class: 'empty' }, e.message));
    return done();
  }

  const box = h('div', { class: 'checkbox' });
  for (const o of ['forward', 'reverse']) {
    const v = check[o];
    box.append(h('div', { class: `checkrow ${v.compatible ? 'ok' : 'no'}` },
      h('span', {}, v.compatible ? '✓' : '✕'),
      h('span', { class: 'label' }, `${o}: ${v.compatible ? `${fmt(v.productSize)} bp product` : v.problems[0] ?? 'ends do not meet'}`)));
  }
  el.append(box);
  for (const w of check.warnings) el.append(h('p', { class: 'sub' }, w));

  const can = c.orientation === 'auto' ? check.orientations.length > 0 : check[c.orientation]?.compatible;
  if (can) {
    const preview = ligationPreview({ vector, insert, circular: c.circular, orientation: c.orientation });
    if (preview) el.append(preview);
  }
  el.append(h('button', {
    class: 'btn primary', disabled: !can, style: { marginTop: '10px' },
    onclick: () => {
      try {
        const result = restrictionClone({
          backbone: { doc: vector.doc, fragment: vector.fragment },
          insert: { doc: insert.doc, fragment: insert.fragment },
          circular: c.circular,
          orientation: c.orientation,
        });
        api.openDoc(result.doc);
        for (const w of result.warnings) toast(w);
        toast(`${result.doc.name}: ${fmt(result.doc.seq.length)} bp, insert ${result.orientation}`);
      } catch (e) {
        toast(e.message, 'error');
      }
    },
  }, can ? 'Ligate' : 'These ends do not ligate'));

  done();
}

function gibsonSection(el, { state, api, docs }, { title, assemble, overlapKey }) {
  const c = state.cloning;
  const overlap = c[overlapKey] ?? 15;
  el.append(h('h3', {}, 'Fragments to assemble'));
  el.append(h('p', { class: 'sub' }, `${title} joins linear fragments by their shared end sequences. Pick them in the order they should join; each needs an end overlapping the next.`));
  for (const d of docs) {
    const on = c.gibson.includes(d.id);
    el.append(h('div', {
      class: 'row', 'aria-current': String(on),
      onclick: () => api.setCloning({ gibson: on ? c.gibson.filter((x) => x !== d.id) : [...c.gibson, d.id] }),
    },
      h('span', { class: 'label' }, d.name),
      h('span', { class: 'meta' }, on ? `#${c.gibson.indexOf(d.id) + 1} · ${fmt(d.seq.length)} bp` : `${fmt(d.seq.length)} bp`)));
  }
  el.append(h('label', { class: 'field', style: { marginTop: '10px' } }, h('span', {}, `Smallest overlap: ${overlap} bp`),
    h('input', { id: `overlap-${overlapKey}`, type: 'range', min: 10, max: 60, value: overlap, oninput: (e) => api.setCloning({ [overlapKey]: Number(e.target.value) }) })));
  el.append(h('label', { class: 'field' }, h('span', {}, 'Shape'),
    h('select', { onchange: (e) => api.setCloning({ circular: e.target.value === 'circular' }) },
      [['circular', 'circular plasmid'], ['linear', 'linear fragment']].map(([v, t]) => h('option', { value: v, selected: c.circular === (v === 'circular') }, t)))));

  const chosen = c.gibson.map((id) => docs.find((d) => d.id === id)).filter(Boolean);
  el.append(h('button', {
    class: 'btn primary', disabled: chosen.length < 2,
    onclick: () => {
      try {
        const result = assemble(chosen.map((d) => ({ doc: d })), {
          circular: c.circular, minOverlap: overlap,
        });
        api.openDoc(result.doc);
        for (const w of result.warnings) toast(w);
        toast(`${result.doc.name}: ${fmt(result.doc.seq.length)} bp from ${chosen.length} fragments`);
      } catch (e) {
        toast(e.message, 'error');
      }
    },
  }, chosen.length < 2 ? 'Pick at least two fragments' : `Assemble ${chosen.length} fragments`));
}

// ---- Golden Gate -------------------------------------------------------------

function goldenGateSection(el, { state, api, docs }) {
  const c = state.cloning;
  const gg = c.goldenGate;
  el.append(h('h3', {}, 'Golden Gate'));
  el.append(h('p', { class: 'sub' }, 'Each file carries a Type IIS site at both ends, pointing inward. The overhangs they leave decide the order; the first file is the backbone.'));
  el.append(h('label', { class: 'field' }, h('span', {}, 'Enzyme'),
    h('select', { onchange: (e) => api.setCloning({ goldenGate: { ...gg, enzyme: e.target.value } }) },
      TYPE_IIS.map((n) => h('option', { value: n, selected: gg.enzyme === n }, `${n} (${ENZYME_BY_NAME[n].site})`)))));

  for (const d of docs) {
    const at = gg.docs.indexOf(d.id);
    el.append(h('div', {
      class: 'row', 'aria-current': String(at >= 0),
      onclick: () => api.setCloning({ goldenGate: { ...gg, docs: at >= 0 ? gg.docs.filter((x) => x !== d.id) : [...gg.docs, d.id] } }),
    },
      h('span', { class: 'label' }, d.name),
      h('span', { class: 'meta' }, at >= 0 ? `#${at + 1}${at === 0 ? ' backbone' : ''} · ${fmt(d.seq.length)} bp` : `${fmt(d.seq.length)} bp`)));
  }

  const chosen = gg.docs.map((id) => docs.find((d) => d.id === id)).filter(Boolean);
  const pieces = chosen.map((d) => {
    try { return { doc: d, ...goldenGatePieces(d, { enzyme: gg.enzyme }) }; } catch (e) { return { doc: d, error: e.message }; }
  });
  if (pieces.length) {
    el.append(h('h3', {}, 'Overhangs'));
    const table = h('table', { class: 'data' }, h('tr', {}, h('th', {}, 'File'), h('th', {}, 'Left'), h('th', {}, 'Right'), h('th', {}, 'Body')));
    for (const p of pieces) {
      table.append(p.error
        ? h('tr', {}, h('td', {}, p.doc.name), h('td', { colSpan: 3, style: { color: 'var(--danger)' } }, p.error))
        : h('tr', {}, h('td', {}, p.doc.name), h('td', { class: 'mono' }, p.left), h('td', { class: 'mono' }, p.right), h('td', {}, `${fmt(p.seq.length)} bp`)));
    }
    el.append(table);
    const ok = pieces.filter((p) => !p.error);
    if (ok.length === pieces.length && pieces.length > 1) {
      const rows = [];
      pieces.forEach((p, i) => {
        const next = pieces[i + 1] ?? (c.circular ? pieces[0] : null);
        if (!next) return;
        const fits = p.right === next.left;
        rows.push(h('div', { class: `checkrow ${fits ? 'ok' : 'no'}` }, h('span', {}, fits ? '✓' : '✕'),
          h('span', { class: 'label' }, `${p.doc.name} → ${next.doc.name}: ${p.right} ${fits ? '=' : '≠'} ${next.left}`)));
      });
      el.append(h('div', { class: 'checkbox' }, rows));
    }
  }
  el.append(h('label', { class: 'field', style: { marginTop: '8px' } }, h('span', {}, 'Shape'),
    h('select', { onchange: (e) => api.setCloning({ circular: e.target.value === 'circular' }) },
      [['circular', 'circular plasmid'], ['linear', 'linear fragment']].map(([v, t]) => h('option', { value: v, selected: c.circular === (v === 'circular') }, t)))));
  el.append(h('button', {
    class: 'btn primary', disabled: chosen.length < 2,
    onclick: () => {
      try {
        const result = goldenGateAssemble(chosen.map((d) => ({ doc: d })), { circular: c.circular, enzyme: gg.enzyme });
        api.openDoc(result.doc);
        for (const e of result.errors) toast(e.message, 'error');
        for (const w of result.warnings) toast(w);
        toast(`${result.doc.name}: ${fmt(result.doc.seq.length)} bp from ${chosen.length} fragments`);
      } catch (e) { toast(e.message, 'error'); }
    },
  }, chosen.length < 2 ? 'Pick at least two files' : `Assemble ${chosen.length} files`));
}

// ---- Recombination -----------------------------------------------------------

function sitesLine(doc) {
  const sites = findRecombinationSites(doc);
  return sites.length ? sites.map((x) => `${x.name}${x.strand === -1 ? ' ←' : ''}`).join(', ') : 'no att, lox or FRT sites';
}

function recombinationSection(el, { state, api, docs }) {
  const c = state.cloning.recomb;
  el.append(h('h3', {}, 'Site-specific recombination'));
  el.append(h('p', { class: 'sub' }, 'Swaps the stretch between two matching sites of one file into the other, inside the sites\' shared core. Works for Gateway att sites (LR and BP), loxP and FRT. Sites come from the files\' own annotations; attB1 and attB2 are also found by sequence.'));
  const pick = (key, label) => {
    const chosen = docs.find((d) => d.id === c[key]) ?? (key === 'insert' ? docs[0] : docs[docs.length - 1]);
    el.append(h('label', { class: 'field' }, h('span', {}, label),
      h('select', { onchange: (e) => api.setCloning({ recomb: { ...c, [key]: e.target.value } }) },
        docs.map((d) => h('option', { value: d.id, selected: d.id === chosen?.id }, `${d.name} (${fmt(d.seq.length)} bp)`)))));
    if (chosen) el.append(h('p', { class: 'sub' }, `Sites: ${sitesLine(chosen)}`));
    return chosen;
  };
  const insert = pick('insert', 'Insert from (entry clone, or attB-flanked fragment)');
  const acceptor = pick('acceptor', 'Into (destination, or donor vector)');
  el.append(h('button', {
    class: 'btn primary', disabled: !insert || !acceptor,
    onclick: () => {
      try {
        const result = recombine(insert, acceptor);
        state.cloning.recomb.byproduct = result.byproduct ? { doc: result.byproduct } : null;
        api.openDoc(result.doc);
        for (const w of result.warnings) toast(w);
        toast(`${result.doc.name}: ${fmt(result.doc.seq.length)} bp`);
      } catch (e) { toast(e.message, 'error'); }
    },
  }, 'Recombine'));
  if (c.byproduct) {
    el.append(h('button', { class: 'btn', style: { marginLeft: '6px' }, onclick: () => { api.openDoc(c.byproduct.doc); } }, 'Open the by-product'));
  }
}

// ---- TOPO / TA / blunt -------------------------------------------------------

function topoSection(el, { doc, state, api, docs }) {
  const c = state.cloning.topo;
  const insertChoices = docs.filter((d) => d.id !== doc.id);
  const sel = state.selection;
  const site = sel ? { start: sel.start, end: sel.end } : state.caret != null ? { start: state.caret, end: state.caret } : null;
  el.append(h('h3', {}, 'TOPO, TA and blunt cloning'));
  el.append(h('p', { class: 'sub' }, `The vector is the file you are looking at (${doc.name}). Select the cloning site in its sequence, or click to put the cursor where the insert goes; a selection is replaced.`));
  el.append(h('p', { class: 'sub' }, site ? (site.end > site.start ? `Site: replaces ${fmt(site.start + 1)}–${fmt(site.end)}` : `Site: insert at ${fmt(site.start + 1)}`) : 'Site: not chosen yet'));
  if (!insertChoices.length) { el.append(h('p', { class: 'empty' }, 'Open the insert as a second file.')); return; }
  const chosen = insertChoices.find((d) => d.id === c.insert) ?? insertChoices[0];
  el.append(h('label', { class: 'field' }, h('span', {}, 'Insert'),
    h('select', { onchange: (e) => api.setCloning({ topo: { ...c, insert: e.target.value } }) },
      insertChoices.map((d) => h('option', { value: d.id, selected: d.id === chosen.id }, `${d.name} (${fmt(d.seq.length)} bp)`)))));
  el.append(h('label', { class: 'field' }, h('span', {}, 'Kind'),
    h('select', { onchange: (e) => api.setCloning({ topo: { ...c, kind: e.target.value } }) },
      [['TA', 'TA (3′ A overhang insert)'], ['blunt', 'Blunt'], ['directional', 'Directional (insert starts CACC)']].map(([v, t]) => h('option', { value: v, selected: c.kind === v }, t)))));
  el.append(h('label', { class: 'field' }, h('span', {}, 'Insert orientation'),
    h('select', { onchange: (e) => api.setCloning({ topo: { ...c, orientation: e.target.value } }) },
      [['forward', 'forward'], ['reverse', 'reverse']].map(([v, t]) => h('option', { value: v, selected: c.orientation === v }, t)))));
  el.append(h('button', {
    class: 'btn primary', disabled: !site,
    onclick: () => {
      try {
        const result = topoClone(doc, chosen, { kind: c.kind, range: site, orientation: c.orientation });
        api.openDoc(result.doc);
        for (const w of result.warnings) toast(w);
        toast(`${result.doc.name}: ${fmt(result.doc.seq.length)} bp`);
      } catch (e) { toast(e.message, 'error'); }
    },
  }, site ? 'Clone' : 'Choose a site first'));
}
