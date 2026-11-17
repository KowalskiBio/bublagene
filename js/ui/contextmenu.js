// @ts-nocheck — legacy DOM wiring; the engines it calls are type-checked.
// The right-click menu.

import { codingSequence, featureSequence, rotateOrigin } from '../model.js';
import { revcomp, translate } from '../seq.js';
import { api } from './api/index.js';
import { addFeatureFromSelection, addPrimerFromSelection, copy, selectedSeq } from './clipboard.js';
import { current, render, state } from './core.js';
import { deleteSelection, edit, replaceSelection } from './documents.js';
import { fmt } from './dom.js';
import { showMenu } from './menu.js';
import { isOrf, isPrimer } from './overlays.js';
import { overlaysFor } from './viewdata.js';

function featureFromEvent(ev) {
  const id = ev.target?.closest?.('[data-feature]')?.getAttribute('data-feature');
  if (!id) return null;
  const doc = current();
  return doc?.features.find((f) => f.id === id) ?? overlaysFor(doc).find((f) => f.id === id) ?? null;
}

export function openContextMenu(ev) {
  const doc = current();
  if (!doc) return;
  ev.preventDefault();
  const feature = featureFromEvent(ev);
  const sel = state.selection;
  const picked = selectedSeq();
  const items = [];

  if (feature && isOrf(feature)) {
    const seg = feature.segments[0];
    items.push({ label: `Select “${feature.name}”`, onClick: () => api.setSelection({ start: seg.start, end: seg.end }) });
    items.push({ label: 'Add as a CDS feature', onClick: () => api.addOrfAsCds(feature.orfInfo) });
    items.push('-');
  } else if (feature && isPrimer(feature)) {
    const seg = feature.segments[0];
    items.push({ label: `Select “${feature.name}” binding site`, onClick: () => api.setSelection({ start: seg.start, end: seg.end }) });
    items.push('-');
  } else if (feature) {
    items.push({ label: `Select “${feature.name}”`, onClick: () => api.focusFeature(feature) });
    items.push({ label: 'Copy feature sequence', onClick: () => copy(featureSequence(doc, feature)) });
    items.push({ label: 'Copy feature (annotation only)', onClick: () => api.copyFeatureAnnotations([feature.id]) });
    if (feature.type === 'CDS') {
      items.push({ label: 'Copy protein', onClick: () => copy(translate(codingSequence(doc, feature)), 'residues') });
    }
    items.push({ label: 'Delete feature', danger: true, onClick: () => edit(`Delete ${feature.name}`, (d) => { d.features = d.features.filter((x) => x.id !== feature.id); }) });
    items.push('-');
  }

  items.push({ label: 'Copy', hint: '⌘C', disabled: !picked, onClick: () => copy(picked) });
  items.push({ label: 'Copy reverse complement', disabled: !picked, onClick: () => copy(revcomp(picked)) });
  items.push({ label: 'Copy translation', disabled: !picked, onClick: () => copy(translate(picked), 'residues') });
  items.push({
    label: 'Copy as FASTA', disabled: !picked,
    onClick: () => copy(`>${doc.name} ${sel.start + 1}-${sel.end}\n${picked.match(/.{1,70}/g).join('\n')}\n`),
  });
  items.push({ label: 'Copy fragment with features', disabled: !sel, onClick: () => api.copyFragment() });
  items.push({ label: 'Paste fragment here', disabled: !state.clipboard, onClick: () => api.pasteFragment() });
  items.push({ label: 'Paste copied features here', disabled: !state.featureClip, onClick: () => api.pasteFeatureAnnotations() });
  items.push('-');
  items.push({ label: 'Make a forward primer', disabled: !picked, onClick: () => addPrimerFromSelection(1) });
  items.push({ label: 'Make a reverse primer', disabled: !picked, onClick: () => addPrimerFromSelection(-1) });
  items.push({ label: 'Design a primer pair…', disabled: !picked, onClick: () => { state.panel = 'Primers'; render(); } });
  items.push('-');
  items.push({ label: 'Add a feature here', disabled: !sel, onClick: addFeatureFromSelection });
  items.push({ label: 'Delete selected bases', danger: true, disabled: !sel, onClick: deleteSelection });
  items.push({ label: 'Replace selected bases…', disabled: !sel, onClick: replaceSelection });
  items.push('-');
  if (doc.circular && sel) items.push({ label: `Start the sequence at ${fmt(sel.start + 1)}`, onClick: () => edit('Set origin', (d) => rotateOrigin(d, sel.start)) });
  items.push({ label: 'Select everything', hint: '⌘A', onClick: () => api.setSelection({ start: 0, end: doc.seq.length }) });
  items.push({ label: 'Clear the selection', disabled: !sel, onClick: () => api.setSelection(null) });
  showMenu(ev.clientX, ev.clientY, items);
}
