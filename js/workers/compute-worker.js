// Worker for bounded heavy work. The main thread falls back to the same
// algorithm when workers are unavailable.
import { fold } from '../fold.js';

self.onmessage = (event) => {
  try {
    const { sequence, limit } = event.data ?? {};
    self.postMessage(fold(sequence, { maxLength: limit }));
  } catch (e) {
    self.postMessage({ __error: e.message });
  }
};
