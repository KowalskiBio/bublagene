// Offload heavy work to a Worker when one is available, otherwise run it here.
//
// The caller supplies a Worker URL (or a Worker factory) and a fallback
// function. Results keep the same shape either way, with `worker: false` when
// the work stayed on the main thread.

export function canUseWorker(env = globalThis) {
  return typeof env.Worker === 'function';
}

// Fallback: run a pure function with structured-clone-compatible arguments.
export async function runLocal(task, payload) {
  return { result: await task(payload), worker: false };
}

// Worker protocol: one message in, one message out, then terminate.
export function runInWorker(worker, payload, { timeoutMs = 30000 } = {}) {
  return new Promise((resolve, reject) => {
    let done = false;
    const timer = setTimeout(() => {
      if (done) return;
      done = true;
      worker.terminate?.();
      reject(new Error(`Worker timed out after ${timeoutMs} ms`));
    }, timeoutMs);
    worker.onmessage = (event) => {
      if (done) return;
      done = true;
      clearTimeout(timer);
      worker.terminate?.();
      resolve({ result: event.data, worker: true });
    };
    worker.onerror = (event) => {
      if (done) return;
      done = true;
      clearTimeout(timer);
      worker.terminate?.();
      reject(new Error(event.message ?? 'Worker failed'));
    };
    worker.postMessage(payload);
  });
}

export async function offload(task, payload, {
  workerUrl = null, createWorker = null, env = globalThis, timeoutMs = 30000,
} = {}) {
  const factory = createWorker ?? (workerUrl && canUseWorker(env) ? () => new env.Worker(workerUrl, { type: 'module' }) : null);
  if (!factory) return runLocal(task, payload);
  try {
    const worker = factory();
    return await runInWorker(worker, payload, { timeoutMs });
  } catch {
    return runLocal(task, payload);
  }
}
