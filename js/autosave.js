// Debounced autosave with an explicit crash-recovery read.
//
// The app still keeps its IndexedDB session; this helper is for periodic
// snapshots and for hosts that prefer a small JSON payload.

export class Autosave {
  /** @param {{write?:Function,read?:Function,clear?:Function,delay?:number,clock?:Function,timer?:any}} [options] */
  constructor({ write, read, clear, delay = 1500, clock = () => Date.now(), timer = null } = {}) {
    this.write = write;
    this.read = read;
    this.clear = clear;
    this.delay = delay;
    this.clock = clock;
    this.timer = timer;
    this.pending = null;
    this.lastSaved = null;
    this._timeout = null;
    this._setTimeout = timer?.setTimeout ?? globalThis.setTimeout?.bind(globalThis);
    this._clearTimeout = timer?.clearTimeout ?? globalThis.clearTimeout?.bind(globalThis);
  }

  schedule(payload) {
    this.pending = payload;
    if (!this._setTimeout) return this;
    if (this._timeout != null) this._clearTimeout?.(this._timeout);
    this._timeout = this._setTimeout(() => { this._timeout = null; this.flush(); }, this.delay);
    return this;
  }

  async flush() {
    if (this.pending == null) return this.lastSaved;
    const payload = this.pending;
    this.pending = null;
    try {
      await this.write(payload);
      this.lastSaved = { payload, at: this.clock() };
      return this.lastSaved;
    } catch (e) {
      this.pending = payload;
      throw e;
    }
  }

  async restore() {
    return (await this.read?.()) ?? null;
  }

  async clearSaved() {
    this.pending = null;
    if (this._timeout != null) this._clearTimeout?.(this._timeout);
    this._timeout = null;
    this.lastSaved = null;
    return (await this.clear?.()) ?? true;
  }
}
