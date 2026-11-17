// @ts-nocheck — legacy DOM wiring; the engines it calls are type-checked.
// Whole-sequence scans (ORFs, enzyme sites) are cheap on a plasmid and slow on a
// genome. While a long sequence is being edited, keep showing the last answer and
// recompute once typing pauses, so each keystroke costs one redraw and nothing more.

import { render } from './core.js';

const answers = new Map();
const timers = new Map();

export const LONG_SEQUENCE = 20000;

export function deferredScan(key, seq, compute, { threshold = LONG_SEQUENCE, delay = 250 } = {}) {
  const hit = answers.get(key);
  if (hit && hit.seq === seq) return hit.value;
  if (!hit || seq.length < threshold) {
    const value = compute();
    answers.delete(key);
    answers.set(key, { seq, value });
    if (answers.size > 200) answers.delete(answers.keys().next().value); // keys vary with the shown enzymes
    return value;
  }
  clearTimeout(timers.get(key));
  timers.set(key, setTimeout(() => {
    timers.delete(key);
    answers.set(key, { seq, value: compute() });
    render();
  }, delay));
  return hit.value;
}
