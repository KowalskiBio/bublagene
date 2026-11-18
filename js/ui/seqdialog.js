// @ts-nocheck — DOM wiring.
// The dialog for pasting a sequence: swapping a feature's bases, replacing a
// selection, or entering a donor. The current sequence (when there is one) sits
// beside a paste box, and a formatted preview shows what will go in, both
// scrollable however long the sequence is.

import { readPastedSequence } from '../seq.js';
import { fmt, h } from './dom.js';

const PER_LINE = 60;
const BLOCK = 10;

// A sequence as numbered lines of blocks, the way sequence records print it.
function formatted(seq, offset = 0) {
  const lines = [];
  const width = String(offset + seq.length).length;
  for (let i = 0; i < seq.length; i += PER_LINE) {
    const line = seq.slice(i, i + PER_LINE).match(new RegExp(`.{1,${BLOCK}}`, 'g')).join(' ');
    lines.push(h('div', { class: 'seqline' }, h('span', { class: 'seqnum' }, String(offset + i + 1).padStart(width, ' ')), line));
  }
  return lines;
}

/**
 * @param {{
 *   title: string, intro?: string, confirmLabel?: string, placeholder?: string,
 *   current?: { label: string, seq: string, offset?: number } | null,
 *   notes?: (parsed: { seq: string, ignored: number }) => string[],
 *   onConfirm: (seq: string) => void,
 * }} opts
 */
export function showSequenceDialog({ title, intro = '', confirmLabel = 'Apply', placeholder = 'Paste the sequence here (plain text or FASTA)', current = null, notes = () => [], onConfirm }) {
  const text = h('textarea', {
    class: 'seqtext', spellcheck: false, autocomplete: 'off', autocapitalize: 'off', wrap: 'soft',
    placeholder, 'aria-label': 'New sequence',
  });
  const preview = h('div', { class: 'seqpreview', hidden: true, tabindex: '0', 'aria-label': 'Preview of the new sequence' });
  const count = h('span', { class: 'sub' }, '');
  const pasteTab = h('button', { class: 'btn', type: 'button', 'aria-pressed': 'true' }, 'Paste');
  const previewTab = h('button', { class: 'btn', type: 'button', 'aria-pressed': 'false' }, 'Preview');
  const info = h('div', { class: 'seqinfo', role: 'status', 'aria-live': 'polite' });
  const apply = h('button', { class: 'btn primary', type: 'button', disabled: true }, confirmLabel);
  const cancel = h('button', { class: 'btn', type: 'button' }, 'Cancel');
  const closeX = h('button', { class: 'btn icon', type: 'button', title: 'Close', 'aria-label': 'Close' }, '×');

  const newPane = h('section', { class: 'seqpane' },
    h('div', { class: 'seqpanehead' }, h('strong', {}, 'New sequence'), count, h('span', { class: 'spacer' }), pasteTab, previewTab),
    h('div', { class: 'seqpanebody' }, text, preview));
  const panes = [];
  if (current) {
    panes.push(h('section', { class: 'seqpane' },
      h('div', { class: 'seqpanehead' }, h('strong', {}, current.label), h('span', { class: 'sub' }, `${fmt(current.seq.length)} bp`)),
      h('div', { class: 'seqpanebody seqpreview', tabindex: '0', 'aria-label': current.label }, ...formatted(current.seq, current.offset ?? 0))));
  }
  panes.push(newPane);

  const box = h('div', { class: 'modal seqmodal', role: 'dialog', 'aria-modal': 'true', 'aria-label': title },
    h('header', { class: 'modalhead' }, h('div', {}, h('h2', {}, title), intro ? h('p', { class: 'sub' }, intro) : null), closeX),
    h('div', { class: `seqpanes${current ? '' : ' single'}` }, ...panes),
    info,
    h('footer', { class: 'modalbtns' }, h('span', { class: 'sub hint' }, '⌘↵ to apply · Esc to cancel'), h('span', { class: 'spacer' }), cancel, apply));
  const backdrop = h('div', { class: 'modalback' }, box);

  let parsed = { seq: '', ignored: 0 };
  let showing = 'paste';
  const show = (which) => {
    showing = which;
    text.hidden = which !== 'paste';
    preview.hidden = which !== 'preview';
    pasteTab.setAttribute('aria-pressed', String(which === 'paste'));
    previewTab.setAttribute('aria-pressed', String(which === 'preview'));
    if (which === 'preview') {
      preview.replaceChildren(...(parsed.seq ? formatted(parsed.seq) : [h('div', { class: 'sub' }, 'Nothing pasted yet.')]));
    }
  };
  const update = () => {
    parsed = readPastedSequence(text.value);
    const n = parsed.seq.length;
    count.textContent = n ? `${fmt(n)} bp${current ? ` (${n - current.seq.length >= 0 ? '+' : '−'}${fmt(Math.abs(n - current.seq.length))})` : ''}` : '';
    const lines = [...(parsed.ignored ? [`${fmt(parsed.ignored)} letter${parsed.ignored === 1 ? '' : 's'} that are not bases were left out.`] : []), ...(n ? notes(parsed) : [])];
    info.replaceChildren(...lines.map((t) => h('div', {}, t)));
    apply.disabled = !n;
    if (showing === 'preview') show('preview');
  };
  text.addEventListener('input', update);
  pasteTab.addEventListener('click', () => { show('paste'); text.focus(); });
  previewTab.addEventListener('click', () => show('preview'));
  update();

  const focusable = () => [...box.querySelectorAll('button:not([disabled]), textarea:not([hidden]), [tabindex="0"]')];
  const close = () => { backdrop.remove(); removeEventListener('keydown', onKey, true); previouslyFocused?.focus?.(); };
  const confirm = () => { if (!parsed.seq) return; close(); onConfirm(parsed.seq); };
  const previouslyFocused = document.activeElement;
  const onKey = (e) => {
    if (e.key === 'Escape') { e.preventDefault(); e.stopPropagation(); close(); return; }
    if (e.key === 'Enter' && (e.metaKey || e.ctrlKey)) { e.preventDefault(); confirm(); return; }
    if (e.key === 'Tab') { // keep the focus inside the dialog
      const items = focusable();
      if (!items.length) return;
      const i = items.indexOf(document.activeElement);
      const next = e.shiftKey ? (i <= 0 ? items.length - 1 : i - 1) : (i === items.length - 1 ? 0 : i + 1);
      e.preventDefault();
      items[next].focus();
    }
  };
  addEventListener('keydown', onKey, true);
  cancel.addEventListener('click', close);
  closeX.addEventListener('click', close);
  apply.addEventListener('click', confirm);
  backdrop.addEventListener('pointerdown', (e) => { if (e.target === backdrop && !text.value.trim()) close(); });
  document.body.append(backdrop);
  text.focus();
}
