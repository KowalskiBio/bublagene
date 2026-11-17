// @ts-nocheck — legacy DOM wiring; the engines it calls are type-checked.
// The page around the view: document tabs, view options, status bar and inspector tabs.

import { reverseComplementDoc, rotateOrigin } from '../model.js';
import { gcContent, meltingTemp, rangeLength, revcomp, subseq } from '../seq.js';
import { api } from './api/index.js';
import { copy } from './clipboard.js';
import { alignments, el, isDirty, render, state } from './core.js';
import { closeDoc, deleteSelection, edit, replaceSelection } from './documents.js';
import { clear, fmt, h } from './dom.js';
import { PANEL_GROUPS, renderPanel } from './panels.js';

export function renderTabs() {
  const box = clear(el('tabs'));
  state.docs.forEach((doc, i) => {
    const dirty = isDirty(doc);
    box.append(h('div', { class: 'tab', role: 'tab', 'aria-selected': String(i === state.active), onclick: () => { state.active = i; state.selection = null; state.caret = null; state.highlight = null; state.selectedFeatures = new Set(); state.traceBase = null; render(); el('panel').scrollTop = 0; } },
      h('span', { class: 'name' }, doc.name),
      dirty ? h('span', { class: 'dirty', title: 'edited' }, '•') : null,
      h('button', { class: 'close', title: 'Close', onclick: (e) => { e.stopPropagation(); closeDoc(i); } }, '×')));
  });
}

export function renderViewBar(doc) {
  const opts = clear(el('viewOpts'));
  const ops = clear(el('editOps'));
  if (!doc) return;
  const toggles = [['showComplement', 'Both strands'], ['showTranslation', 'Translation'], ['showSixFrame', 'Six frames'], ['showOrfs', 'ORFs']];
  if (doc.primers.length) toggles.push(['showPrimers', 'Primers']);
  if (alignments.get(doc.id)?.length) toggles.push(['showReads', `Reads (${alignments.get(doc.id).length})`]);
  for (const [key, label] of toggles) {
    opts.append(h('label', { class: 'sub', style: { display: 'flex', gap: '4px', alignItems: 'center' } },
      h('input', { type: 'checkbox', checked: state[key], onchange: (e) => { state[key] = e.target.checked; render(); } }), label));
  }
  if (state.showOrfs) {
    opts.append(h('label', { class: 'sub', style: { display: 'flex', gap: '4px', alignItems: 'center' }, title: 'Shortest reading frame to draw, in codons' },
      'min',
      h('input', {
        id: 'orfMin', type: 'number', min: 20, max: 2000, step: 10, value: state.orfMin, style: { width: '64px', padding: '2px 5px' },
        onchange: (e) => { state.orfMin = Math.max(20, Number(e.target.value) || 100); render(); },
      }), 'aa'));
  }
  if (state.view === 'map' || state.view === 'split') {
    opts.append(h('label', { class: 'sub', style: { display: 'flex', gap: '4px', alignItems: 'center' } }, 'Map style',
      h('select', { id: 'mapStyle', onchange: (e) => api.setMapStyle(e.target.value), style: { padding: '2px 4px' } },
        [['auto', 'as the file'], ['circular', 'circular'], ['linear', 'linear']].map(([v, t]) => h('option', { value: v, selected: state.mapStyle === v }, t)))));
  }
  const sel = state.selection;
  ops.append(h('button', { class: 'btn', disabled: !sel, onclick: deleteSelection }, 'Delete'));
  ops.append(h('button', { class: 'btn', disabled: !sel, onclick: replaceSelection }, 'Replace…'));
  ops.append(h('button', {
    class: 'btn', disabled: !sel || !doc.circular,
    title: 'Make the selection start base 1',
    onclick: () => edit('Set origin', (d) => { rotateOrigin(d, sel.start); }),
  }, 'Set origin'));
  ops.append(h('button', { class: 'btn', onclick: () => edit('Reverse complement', (d) => reverseComplementDoc(d)) }, 'Rev. comp.'));
  ops.append(h('select', {
    id: 'exportMenu', class: 'btn', title: 'Export or print',
    onchange: (e) => { const v = e.target.value; e.target.value = ''; if (v === 'print') api.print(); else if (v) api.exportMap(v); },
  },
  h('option', { value: '' }, 'Export…'),
  h('option', { value: 'svg' }, 'Map as SVG'), h('option', { value: 'png' }, 'Map as PNG'), h('option', { value: 'pdf' }, 'Map as PDF'),
  h('option', { value: 'print' }, 'Print…')));
}

export function renderStatus(doc) {
  const bar = clear(el('status'));
  if (!doc) return;
  bar.append(h('span', {}, h('b', {}, fmt(doc.seq.length)), ' bp · ', doc.circular ? 'circular' : 'linear'));
  bar.append(h('span', {}, 'GC ', h('b', {}, `${gcContent(doc.seq).toFixed(1)}%`)));
  const sel = state.selection;
  if (sel) {
    const len = rangeLength(sel, doc.seq.length);
    const s = subseq(doc.seq, sel.start, sel.end, doc.circular);
    bar.append(h('span', {}, 'Selection ', h('b', {}, `${fmt(sel.start + 1)}–${fmt(sel.end)}`), ` · ${fmt(len)} bp · GC ${gcContent(s).toFixed(1)}%`));
    if (len >= 8 && len <= 60) bar.append(h('span', {}, 'Tm ', h('b', {}, `${meltingTemp(s).toFixed(1)}°C`)));
    bar.append(h('button', { class: 'btn', style: { padding: '2px 8px' }, onclick: () => copy(s) }, 'Copy'));
    bar.append(h('button', { class: 'btn', style: { padding: '2px 8px' }, onclick: () => copy(revcomp(s)) }, 'Copy rev. comp.'));
  } else {
    bar.append(h('span', {}, 'Drag to select · right-click for everything you can do with it'));
  }
}

export function renderPanels(doc) {
  const tabs = clear(el('panelTabs'));
  for (const [group, names] of PANEL_GROUPS) {
    const box = h('div', { class: 'tabgroup', role: 'presentation' }, names.length > 1 ? h('span', { class: 'grouplabel' }, group) : null);
    for (const name of names) {
      box.append(h('button', {
        role: 'tab', 'aria-selected': String(state.panel === name),
        onclick: () => { state.panel = name; render(); el('panel').scrollTop = 0; },
      }, name));
    }
    tabs.append(box);
  }
  const panel = clear(el('panel'));
  if (!doc) { panel.append(h('p', { class: 'empty' }, 'Open a file to see its details.')); return; }
  panel.append(renderPanel(state.panel, {
    doc, state, api,
    others: state.docs.filter((d) => d.id !== doc.id),
    docs: state.docs,
  }));
}
