// Snapshot undo/redo with coalescing.
//
// `record(before, label, { coalesceKey, now })` stores the state *before* an
// edit. Consecutive records with the same coalesce key inside `coalesceMs`
// keep the earliest snapshot, so a burst of typing undoes as one action.

export class History {
  constructor({ limit = 100, coalesceMs = 750, clock = () => Date.now() } = {}) {
    this.past = [];
    this.future = [];
    this.limit = limit;
    this.coalesceMs = coalesceMs;
    this.clock = clock;
    this.last = null;
  }

  get depth() { return this.past.length; }
  canUndo() { return this.past.length > 0; }
  canRedo() { return this.future.length > 0; }

  record(before, label = 'edit', { coalesceKey = null, now = this.clock() } = {}) {
    const coalesce = coalesceKey != null && this.last && this.last.key === coalesceKey && now - this.last.time <= this.coalesceMs;
    if (!coalesce) {
      this.past.push({ state: before, label });
      if (this.past.length > this.limit) this.past.shift();
      this.last = coalesceKey != null ? { key: coalesceKey, time: now } : null;
    } else {
      this.last.time = now;
    }
    this.future.length = 0;
    return this;
  }

  undo(current) {
    const entry = this.past.pop();
    if (!entry) return null;
    this.future.push({ state: current, label: entry.label });
    this.last = null;
    return entry.state;
  }

  redo(current) {
    const entry = this.future.pop();
    if (!entry) return null;
    this.past.push({ state: current, label: entry.label });
    this.last = null;
    return entry.state;
  }

  clear() {
    this.past.length = 0;
    this.future.length = 0;
    this.last = null;
    return this;
  }

  labels() {
    return this.past.map((e) => e.label);
  }
}
