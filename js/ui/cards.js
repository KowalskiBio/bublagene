// @ts-nocheck — legacy DOM wiring; the engines it calls are type-checked.
// The cards shown in place of a view: welcome, restore and empty document.

import { newDocument, replaceRange } from '../model.js';
import { cleanSequence } from '../seq.js';
import { el, state } from './core.js';
import { edit, openDoc } from './documents.js';
import { fmt, h, toast } from './dom.js';
import { restoreSession, startFresh } from './persist.js';
import { exampleDocument } from './sample.js';
import { timeAgo } from './session.js';

// A new document has no bases, and nothing to select, so give it a way in.
export function emptyCard(doc) {
  const box = h('textarea', { id: 'pasteSeq', rows: 6, placeholder: 'ACGT… or a FASTA record', class: 'mono', style: { width: '100%' } });
  const add = () => {
    const text = box.value.replace(/^>.*$/gm, '');
    const seq = cleanSequence(text);
    if (!seq) { toast('No bases found in what was pasted', 'error'); return; }
    edit('Add sequence', (d) => { replaceRange(d, 0, 0, seq); });
  };
  return h('div', { class: 'welcome' },
    h('h1', {}, `${doc.name} is empty`),
    h('p', {}, 'Paste a sequence to start, or open a file. Spaces, digits and line breaks are ignored.'),
    h('label', { class: 'field' }, h('span', {}, 'Sequence'), box),
    h('div', { class: 'group' }, h('button', { class: 'btn primary', onclick: add }, 'Add sequence')));
}

export function welcome() {
  if (state.restorable) return restoreCard(state.restorable);
  const drop = h('div', { class: 'drop' },
    h('p', { style: { margin: '0 0 12px' } }, 'Drop a sequence file here'),
    h('button', { class: 'btn primary', onclick: () => el('fileInput').click() }, 'Choose a file…'));
  const box = h('div', { class: 'welcome' },
    h('h1', {}, 'A plasmid map in the browser'),
    h('p', {}, 'Open a .dna, GenBank or FASTA file to see its map, sequence, features and cut sites. Files stay on this computer.'),
    drop,
    h('div', { class: 'group' },
      h('button', { class: 'btn', onclick: () => openDoc(exampleDocument()) }, 'Open the example fragment'),
      h('button', { class: 'btn', onclick: () => openDoc(newDocument({ name: 'Untitled', seq: '', circular: true })) }, 'Start an empty sequence')),
    h('ul', { style: { marginTop: '22px' } },
      h('li', {}, 'Click a feature or an enzyme label to select it'),
      h('li', {}, 'Drag across the sequence to select a region, then ⌘C to copy it'),
      h('li', {}, 'Right-click anywhere for copying, primers and editing'),
      h('li', {}, 'Saving a .dna file keeps the parts BublaGene does not edit, such as the cloning history')),
    h('p', { class: 'sub', style: { marginTop: '18px' } },
      'Free software under the AGPL-3.0. ',
      h('a', { href: 'https://github.com/KowalskiBio/bublagene', target: '_blank', rel: 'noopener', style: { color: 'var(--accent)' } }, 'Source code')));
  return box;
}

export function restoreCard(saved) {
  const docs = saved.docs ?? [];
  const list = h('ul', { style: { listStyle: 'none', padding: '0', margin: '0 0 20px' } });
  for (const d of docs) {
    list.append(h('li', {
      style: { display: 'flex', gap: '10px', alignItems: 'baseline', padding: '7px 0', borderBottom: '1px solid var(--line)' },
    },
      h('span', { style: { fontWeight: '500' } }, d.name),
      h('span', { class: 'sub' }, `${fmt((d.seq ?? '').length)} bp · ${d.circular ? 'circular' : 'linear'}`),
      d.restoredDirty ? h('span', { class: 'sub', style: { color: 'var(--accent)' } }, '· unsaved changes') : null));
  }
  return h('div', { class: 'welcome' },
    h('h1', {}, 'Pick up where you left off?'),
    h('p', {}, `${docs.length === 1 ? 'One file was' : `${docs.length} files were`} open ${timeAgo(saved.savedAt)}.`),
    list,
    h('div', { class: 'group' },
      h('button', { class: 'btn primary', onclick: () => restoreSession(saved) }, docs.length === 1 ? 'Reopen it' : 'Reopen them'),
      h('button', { class: 'btn', onclick: startFresh }, 'Start fresh')),
    h('p', { class: 'sub', style: { marginTop: '18px' } },
      'Kept in this browser only, never uploaded. “Start fresh” forgets it.'));
}
