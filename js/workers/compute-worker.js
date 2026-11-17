// Worker for bounded heavy work. The main thread falls back to the same
// algorithms when workers are unavailable.
import { fold } from '../fold.js';
import { alignPair } from '../align.js';

self.onmessage = (event) => {
  try {
    const payload = event.data ?? {};
    if (payload.kind === 'align') self.postMessage(alignPair(payload.a, payload.b, payload.options ?? {}));
    else self.postMessage(fold(payload.sequence, { maxLength: payload.limit }));
  } catch (e) {
    self.postMessage({ __error: e.message });
  }
};
