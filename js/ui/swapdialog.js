// @ts-nocheck — DOM wiring.
// The dialog for swapping a feature's sequence: paste the new bases, see what
// the swap will do, confirm.

import { featureLength, featureSpan } from '../model.js';
import { readPastedSequence } from '../seq.js';
import { fmt, h } from './dom.js';

export function showSwapDialog(doc, feature, onConfirm) {
  const span = featureSpan(feature);
  const length = featureLength(feature, doc.seq.length);
  // annotations lying wholly inside the feature go with the bases they described
  const inside = doc.features.filter((f) => f.id !== feature.id
    && f.segments.every((s) => s.start >= span.start && s.end <= span.end));

  const text = h('textarea', {
    class: 'swaptext', rows: 9, spellcheck: false, autocomplete: 'off', autocapitalize: 'off',
    placeholder: 'Paste the new sequence (plain text or FASTA)',
    'aria-label': 'New sequence',
  });
  const info = h('div', { class: 'swapinfo', role: 'status' });
  const swap = h('button', { class: 'btn primary', type: 'button', disabled: true }, 'Swap');
  const cancel = h('button', { class: 'btn', type: 'button' }, 'Cancel');

  const direction = feature.strand === -1
    ? 'It is on the reverse strand, so paste it as the feature reads (5′→3′ along the feature), not as the top strand shows it.'
    : 'Paste it in the same direction as the feature (5′→3′).';
  const box = h('div', { class: 'modal', role: 'dialog', 'aria-modal': 'true', 'aria-label': `Swap ${feature.name}` },
    h('h2', {}, `Swap “${feature.name}”`),
    h('p', { class: 'sub' }, `Replaces its ${fmt(length)} bp (${fmt(span.start + 1)}–${fmt(span.end)}) with the sequence you paste. ${direction}`),
    text, info,
    h('div', { class: 'modalbtns' }, cancel, swap));
  const backdrop = h('div', { class: 'modalback' }, box);

  let parsed = { seq: '', ignored: 0 };
  const update = () => {
    parsed = readPastedSequence(text.value);
    const n = parsed.seq.length;
    const notes = [];
    if (n) {
      const d = n - length;
      notes.push(`${fmt(n)} bp (${d === 0 ? 'same length' : `${d > 0 ? '+' : '−'}${fmt(Math.abs(d))} bp`})`);
      if (feature.type === 'CDS' && n % 3) notes.push('Not a multiple of 3: the reading frame will shift.');
      if (feature.type === 'CDS' && n % 3 === 0 && length % 3 === 0 && n !== length) notes.push('Length changes in whole codons, so the frame after it is kept.');
    }
    if (parsed.ignored) notes.push(`${fmt(parsed.ignored)} letter${parsed.ignored === 1 ? '' : 's'} that are not bases were left out.`);
    if (inside.length) {
      const names = inside.slice(0, 4).map((f) => f.name).join(', ');
      notes.push(`${inside.length} annotation${inside.length === 1 ? '' : 's'} inside it will be removed: ${names}${inside.length > 4 ? '…' : ''}.`);
    }
    info.replaceChildren(...notes.map((t) => h('div', {}, t)));
    swap.disabled = !n;
  };
  text.addEventListener('input', update);
  update();

  const close = () => { backdrop.remove(); removeEventListener('keydown', onKey, true); };
  const confirm = () => { if (!parsed.seq) return; close(); onConfirm(parsed.seq); };
  const onKey = (e) => {
    if (e.key === 'Escape') { e.preventDefault(); e.stopPropagation(); close(); }
    else if (e.key === 'Enter' && (e.metaKey || e.ctrlKey)) { e.preventDefault(); confirm(); }
  };
  addEventListener('keydown', onKey, true);
  cancel.addEventListener('click', close);
  swap.addEventListener('click', confirm);
  backdrop.addEventListener('pointerdown', (e) => { if (e.target === backdrop) close(); });
  document.body.append(backdrop);
  text.focus();
}
