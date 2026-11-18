// @ts-nocheck — DOM wiring.
// Swapping a feature's sequence: what the dialog says about the feature and the swap.

import { featureLength, featureSequence, featureSpan } from '../model.js';
import { fmt } from './dom.js';
import { showSequenceDialog } from './seqdialog.js';

export function swapDialog(doc, feature, onConfirm) {
  const span = featureSpan(feature);
  const length = featureLength(feature, doc.seq.length);
  // annotations lying wholly inside the feature go with the bases they described
  const inside = doc.features.filter((f) => f.id !== feature.id
    && f.segments.every((s) => s.start >= span.start && s.end <= span.end));
  const reverse = feature.strand === -1;
  showSequenceDialog({
    title: `Swap “${feature.name}”`,
    intro: `${feature.type} · ${fmt(span.start + 1)}–${fmt(span.end)} · ${reverse ? 'reverse strand' : 'forward strand'}. `
      + (reverse
        ? "Paste the new sequence as the feature reads (5′→3′ along the feature), not as the top strand shows it."
        : 'Paste the new sequence in the same direction as the feature (5′→3′).'),
    current: { label: reverse ? 'Current (as the feature reads)' : 'Current', seq: featureSequence(doc, feature), offset: 0 },
    confirmLabel: 'Swap',
    notes: ({ seq }) => {
      const n = seq.length;
      const out = [];
      if (feature.type === 'CDS' && n % 3) out.push('Not a multiple of 3: the reading frame will shift.');
      else if (feature.type === 'CDS' && n !== length && length % 3 === 0) out.push('Length changes in whole codons, so the frame after it is kept.');
      if (inside.length) {
        const names = inside.slice(0, 4).map((f) => f.name).join(', ');
        out.push(`${inside.length} annotation${inside.length === 1 ? '' : 's'} inside it will be removed: ${names}${inside.length > 4 ? '…' : ''}.`);
      }
      return out;
    },
    onConfirm,
  });
}
