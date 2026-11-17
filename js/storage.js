// Saving and opening files in the browser, with graceful fallbacks.
//
// Preference order: File System Access API → OPFS → a normal download. The
// functions take an `env` so they can be exercised without a browser.

export function supportsFileSystemAccess(env = globalThis) {
  return typeof env.showSaveFilePicker === 'function' && typeof env.showOpenFilePicker === 'function';
}

export function supportsOPFS(env = globalThis) {
  return Boolean(env.navigator?.storage?.getDirectory);
}

function download(data, name, env = globalThis) {
  const blob = new Blob([data]);
  const url = env.URL.createObjectURL(blob);
  const a = env.document?.createElement?.('a');
  if (a) {
    a.href = url;
    a.download = name;
    a.click();
  }
  env.setTimeout?.(() => env.URL.revokeObjectURL(url), 1000);
  return { method: 'download', name };
}

export async function saveDocument({ data, name = 'sequence.gb', types = null, env = globalThis } = {}) {
  if (supportsFileSystemAccess(env)) {
    const handle = await env.showSaveFilePicker({ suggestedName: name, types: types ?? undefined });
    const writable = await handle.createWritable();
    await writable.write(data);
    await writable.close();
    return { method: 'file-system-access', name: handle.name ?? name, handle };
  }
  if (supportsOPFS(env)) {
    const root = await env.navigator.storage.getDirectory();
    const handle = await root.getFileHandle(name, { create: true });
    const writable = await handle.createWritable();
    await writable.write(data);
    await writable.close();
    return { method: 'opfs', name, handle };
  }
  return download(data, name, env);
}

export async function openDocument({ types = null, multiple = false, env = globalThis } = {}) {
  if (supportsFileSystemAccess(env)) {
    const handles = await env.showOpenFilePicker({ multiple, types: types ?? undefined });
    const list = Array.isArray(handles) ? handles : [handles];
    const files = [];
    for (const handle of list) {
      const file = await handle.getFile();
      files.push({ name: handle.name ?? file.name, data: new Uint8Array(await file.arrayBuffer()), handle });
    }
    return multiple ? files : files[0] ?? null;
  }
  if (env.document?.createElement) {
    return new Promise((resolve) => {
      const input = env.document.createElement('input');
      input.type = 'file';
      input.multiple = multiple;
      input.onchange = async () => {
        const files = [];
        for (const file of input.files ?? []) files.push({ name: file.name, data: new Uint8Array(await file.arrayBuffer()) });
        resolve(multiple ? files : files[0] ?? null);
      };
      input.click();
    });
  }
  return null;
}

// Small key/value helpers for autosave payloads (OPFS first, localStorage next).
export async function writeAutosave(key, payload, env = globalThis) {
  try {
    if (supportsOPFS(env)) {
      const root = await env.navigator.storage.getDirectory();
      const handle = await root.getFileHandle(`${key}.json`, { create: true });
      const writable = await handle.createWritable();
      await writable.write(typeof payload === 'string' ? payload : JSON.stringify(payload));
      await writable.close();
      return { method: 'opfs' };
    }
    env.localStorage?.setItem(key, typeof payload === 'string' ? payload : JSON.stringify(payload));
    return { method: 'localStorage' };
  } catch (e) {
    return { method: 'failed', error: e.message };
  }
}

export async function readAutosave(key, env = globalThis) {
  try {
    if (supportsOPFS(env)) {
      const root = await env.navigator.storage.getDirectory();
      const handle = await root.getFileHandle(`${key}.json`);
      const file = await handle.getFile();
      return await file.text();
    }
    return env.localStorage?.getItem(key) ?? null;
  } catch {
    return null;
  }
}

export async function clearAutosave(key, env = globalThis) {
  try {
    env.localStorage?.removeItem?.(key);
    if (supportsOPFS(env)) {
      const root = await env.navigator.storage.getDirectory();
      await root.removeEntry?.(`${key}.json`);
    }
    return true;
  } catch {
    return false;
  }
}
