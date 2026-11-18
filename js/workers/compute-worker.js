// Worker for bounded heavy work. The main thread falls back to the same
// algorithms when workers are unavailable.
import { fold } from '../fold.js';
import { alignPair } from '../align.js';
import { detectFeatures } from '../detect.js';

self.onmessage = (event) => {
  try {
    const payload = event.data ?? {};
    if (payload.kind === 'align') self.postMessage(alignPair(payload.a, payload.b, payload.options ?? {}));
    else if (payload.kind === 'detect') self.postMessage(detectFeatures(payload.doc, payload.library, payload.options ?? {}));
    else self.postMessage(fold(payload.sequence, { maxLength: payload.limit }));
  } catch (e) {
    self.postMessage({ __error: e.message });
  }
};
