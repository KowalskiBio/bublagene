// The Cloning panel. Everything here is presentation: choosing which fragments
// to join and showing what comes back. The biology is in js/cloning.js.

import { h, clear, fmt, toast } from './dom.js';
import { ENZYME_BY_NAME, COMMON_SET, findSites, digest } from '../enzymes.js';
import { restrictionClone, gibsonAssemble, checkLigation, pieceFromFragment, describeEnd, cloningSummary } from '../cloning.js';

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

  el.append(h('div', { class: 'seg', style: { marginBottom: '12px' } },
    h('button', { 'aria-pressed': String(c.mode === 'restriction'), onclick: () => api.setCloning({ mode: 'restriction' }) }, 'Restriction'),
    h('button', { 'aria-pressed': String(c.mode === 'gibson'), onclick: () => api.setCloning({ mode: 'gibson' }) }, 'Gibson')));

  if (c.mode === 'gibson') return gibsonSection(el, { state, api, docs });

  const vector = sidePicker(el, { label: 'Vector', side: 'vector', state, api, docs });
  const insert = sidePicker(el, { label: 'Insert', side: 'insert', state, api, docs });

  el.append(h('h3', {}, 'Product'));
  el.append(h('label', { class: 'field' }, h('span', {}, 'Shape'),
    h('select', { onchange: (e) => api.setCloning({ circular: e.target.value === 'circular' }) },
      [['circular', 'circular plasmid'], ['linear', 'linear fragment']].map(([v, t]) => h('option', { value: v, selected: c.circular === (v === 'circular') }, t)))));
  el.append(h('label', { class: 'field' }, h('span', {}, 'Insert orientation'),
    h('select', { onchange: (e) => api.setCloning({ orientation: e.target.value }) },
      [['auto', 'whichever fits'], ['forward', 'forward'], ['reverse', 'reverse']].map(([v, t]) => h('option', { value: v, selected: c.orientation === v }, t)))));

  if (!vector || !insert) return;

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
    return;
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

  if (doc?.cloning) el.append(provenance(doc));
}

function gibsonSection(el, { state, api, docs }) {
  const c = state.cloning;
  el.append(h('h3', {}, 'Fragments to assemble'));
  el.append(h('p', { class: 'sub' }, 'Pick them in the order they should join. Each needs an end overlapping the next.'));
  for (const d of docs) {
    const on = c.gibson.includes(d.id);
    el.append(h('div', {
      class: 'row', 'aria-current': String(on),
      onclick: () => api.setCloning({ gibson: on ? c.gibson.filter((x) => x !== d.id) : [...c.gibson, d.id] }),
    },
      h('span', { class: 'label' }, d.name),
      h('span', { class: 'meta' }, on ? `#${c.gibson.indexOf(d.id) + 1} · ${fmt(d.seq.length)} bp` : `${fmt(d.seq.length)} bp`)));
  }
  el.append(h('label', { class: 'field', style: { marginTop: '10px' } }, h('span', {}, `Smallest overlap: ${c.minOverlap} bp`),
    h('input', { type: 'range', min: 10, max: 60, value: c.minOverlap, oninput: (e) => api.setCloning({ minOverlap: Number(e.target.value) }) })));
  el.append(h('label', { class: 'field' }, h('span', {}, 'Shape'),
    h('select', { onchange: (e) => api.setCloning({ circular: e.target.value === 'circular' }) },
      [['circular', 'circular plasmid'], ['linear', 'linear fragment']].map(([v, t]) => h('option', { value: v, selected: c.circular === (v === 'circular') }, t)))));

  const chosen = c.gibson.map((id) => docs.find((d) => d.id === id)).filter(Boolean);
  el.append(h('button', {
    class: 'btn primary', disabled: chosen.length < 2,
    onclick: () => {
      try {
        const result = gibsonAssemble(chosen.map((d) => ({ doc: d })), {
          circular: c.circular, minOverlap: c.minOverlap,
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

// What this document was made from, once it is a cloning product.
function provenance(doc) {
  const s = cloningSummary(doc);
  const box = h('div');
  if (!s) return box;
  box.append(h('h3', {}, 'How this was made'));
  box.append(h('p', { class: 'sub' }, `${s.operation} · ${fmt(s.size)} bp · ${s.circular ? 'circular' : 'linear'}`));
  for (const i of s.inputs) {
    box.append(h('div', { class: 'row' },
      h('span', { class: 'label' }, i.name),
      h('span', { class: 'meta' }, `${i.role}${i.orientation && i.orientation !== 'forward' ? ` · ${i.orientation}` : ''} · ${fmt(i.size)} bp`)));
  }
  for (const j of s.junctions) {
    box.append(h('p', { class: 'sub mono', style: { fontSize: '11.5px' } },
      `${fmt(j.position)}: ${j.kind}${j.enzymes.length ? ` (${j.enzymes.join('/')})` : ''}${j.regenerated.length ? ` · remakes ${j.regenerated.join(', ')}` : ''}`));
  }
  return box;
}
