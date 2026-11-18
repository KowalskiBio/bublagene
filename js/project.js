// Project archives and native-record migrations.
//
// A project is one JSON document holding every open sequence (with stable
// document and feature IDs), notes, primers, cloning/assembly records,
// alignments and UI settings. Readers migrate older schemas forward.

import { newDocument, newId } from './model.js';
import { contentHash } from './assembly.js';

export const PROJECT_SCHEMA = 'bublagene/project/1';
export const ALIGNMENT_SCHEMA = 'bublagene/alignment/1';

export const CURRENT_SCHEMAS = {
  project: PROJECT_SCHEMA,
  cloning: 'bublagene/cloning/1',
  assembly: 'bublagene/assembly/1',
  alignment: ALIGNMENT_SCHEMA,
};

const clone = (x) => (x == null ? x : structuredClone(x));

function serialiseDoc(doc) {
  return {
    id: doc.id,
    name: doc.name,
    seq: doc.seq,
    circular: !!doc.circular,
    flags: doc.flags,
    features: clone(doc.features ?? []),
    primers: clone(doc.primers ?? []),
    notes: clone(doc.notes ?? {}),
    history: clone(doc.history ?? []),
    cloning: clone(doc.cloning ?? null),
    assembly: clone(doc.assembly ?? null),
    alignment: clone(doc.alignment ?? null),
    seqEdited: !!doc.seqEdited,
  };
}

// A trace's dye curves are four long runs of 16-bit numbers; as JSON arrays they
// would dwarf the rest of the project, so they travel as base64.
const toBase64 = (u16) => {
  const bytes = new Uint8Array(u16.buffer, u16.byteOffset, u16.byteLength);
  let bin = '';
  for (let i = 0; i < bytes.length; i += 0x8000) bin += String.fromCharCode(...bytes.subarray(i, i + 0x8000));
  return btoa(bin);
};
const fromBase64 = (text) => {
  const bin = atob(text);
  const bytes = new Uint8Array(bin.length);
  for (let i = 0; i < bin.length; i++) bytes[i] = bin.charCodeAt(i);
  return new Uint16Array(bytes.buffer, 0, bytes.length >> 1);
};

export function encodeTrace(trace) {
  const samples = {};
  for (const c of 'ACGT') samples[c] = toBase64(Uint16Array.from(trace.samples?.[c] ?? []));
  return { ...clone({ ...trace, samples: undefined }), samples, encoded: 'u16-base64' };
}

export function decodeTrace(record) {
  if (record.encoded !== 'u16-base64') return record;
  const samples = {};
  for (const c of 'ACGT') samples[c] = fromBase64(record.samples?.[c] ?? '');
  const { encoded, ...rest } = record;
  return { ...rest, samples };
}

export function projectRecord(docs, {
  settings = {}, alignments = {}, traces = {}, active = 0, created = new Date().toISOString(), name = null,
} = {}) {
  const list = Array.isArray(docs) ? docs : [docs];
  const align = alignments instanceof Map
    ? [...alignments.entries()].map(([docId, reads]) => ({ docId, reads: clone(reads ?? []) }))
    : clone(alignments ?? []);
  const tr = traces instanceof Map
    ? [...traces.entries()].map(([docId, list]) => ({ docId, traces: (list ?? []).map(encodeTrace) }))
    : (Array.isArray(traces) ? traces : []).map((e) => ({ ...e, traces: (e.traces ?? []).map((t) => (t.encoded ? t : encodeTrace(t))) }));
  return {
    schema: PROJECT_SCHEMA,
    created,
    name: name ?? list[0]?.name ?? 'project',
    active,
    settings: clone(settings),
    docs: list.map(serialiseDoc),
    alignments: Array.isArray(align) ? align : [],
    traces: Array.isArray(tr) ? tr : [],
  };
}

export function projectToJSON(docs, options = {}) {
  return JSON.stringify(projectRecord(docs, options), null, 2);
}

// ---- migrations -------------------------------------------------------------

const MIGRATIONS = {
  cloning: {
    'bublagene/cloning/0': (r) => ({
      ...r,
      schema: CURRENT_SCHEMAS.cloning,
      type: r.type ?? r.operation ?? 'restriction',
      operation: r.operation ?? r.type ?? 'restriction',
      product: r.product ?? { name: r.name ?? 'construct', seqLen: r.seqLen ?? 0 },
      inputs: r.inputs ?? (r.fragments ?? []).map((f) => ({ name: f.name, seqLen: f.seqLen ?? 0, role: f.role ?? 'insert', orientation: f.orientation ?? 'forward' })),
      junctions: r.junctions ?? [],
    }),
  },
  assembly: {
    'bublagene/assembly/0': (r) => ({
      ...r,
      schema: CURRENT_SCHEMAS.assembly,
      method: r.method ?? 'homology',
      fragments: r.fragments ?? [],
      junctions: r.junctions ?? [],
      product: r.product ?? { name: r.name ?? 'assembly', seqLen: r.seqLen ?? 0, digest: r.digest ?? '' },
    }),
  },
  alignment: {
    'bublagene/alignment/0': (r) => ({
      ...r,
      schema: ALIGNMENT_SCHEMA,
      target: r.target ?? r.reference ?? null,
      query: r.query ?? r.read ?? null,
      score: r.score ?? 0,
      identity: r.identity ?? 0,
      alignedA: r.alignedA ?? r.targetAligned ?? '',
      alignedB: r.alignedB ?? r.queryAligned ?? '',
    }),
  },
  project: {
    'bublagene/project/0': (r) => ({
      ...r,
      schema: PROJECT_SCHEMA,
      docs: r.docs ?? r.documents ?? [],
      alignments: r.alignments ?? [],
      traces: r.traces ?? [],
      settings: r.settings ?? {},
      active: r.active ?? 0,
    }),
  },
};

export function migrateRecord(record, kind) {
  if (!record || typeof record !== 'object') return record;
  const table = MIGRATIONS[kind];
  if (!table) return record;
  let current = record;
  let guard = 0;
  while (current.schema && current.schema !== CURRENT_SCHEMAS[kind] && table[current.schema] && guard++ < 10) {
    current = table[current.schema](current);
  }
  return current;
}

export function migrateNativeRecords(doc) {
  if (doc.cloning) doc.cloning = migrateRecord(doc.cloning, 'cloning');
  if (doc.assembly) doc.assembly = migrateRecord(doc.assembly, 'assembly');
  if (doc.alignment) doc.alignment = migrateRecord(doc.alignment, 'alignment');
  return doc;
}

function stableFeature(feature) {
  const f = clone(feature);
  if (!f.id) f.id = newId();
  f.segments = (f.segments ?? []).map((s) => ({ ...s }));
  f.qualifiers = f.qualifiers ?? [];
  f.extra = f.extra ?? {};
  return f;
}

export function projectFromJSON(value) {
  const parsed = typeof value === 'string' ? JSON.parse(value) : clone(value);
  const record = migrateRecord(parsed, 'project');
  if (!record || record.schema !== PROJECT_SCHEMA) throw new Error(`Unsupported project schema ${record?.schema ?? '(none)'}`);
  const warnings = [];
  const docs = (record.docs ?? []).map((d) => {
    const doc = newDocument({ name: d.name ?? 'Untitled', seq: d.seq ?? '', circular: !!d.circular });
    doc.id = d.id ?? doc.id;
    doc.flags = d.flags ?? doc.flags;
    doc.features = (d.features ?? []).map(stableFeature);
    doc.primers = clone(d.primers ?? []);
    doc.notes = { ...doc.notes, ...(d.notes ?? {}) };
    doc.history = clone(d.history ?? []);
    doc.cloning = d.cloning ? migrateRecord(d.cloning, 'cloning') : null;
    doc.assembly = d.assembly ? migrateRecord(d.assembly, 'assembly') : null;
    doc.alignment = d.alignment ? migrateRecord(d.alignment, 'alignment') : null;
    doc.seqEdited = !!d.seqEdited;
    doc.source = { format: 'project' };
    if (d.cloning && d.cloning.schema !== CURRENT_SCHEMAS.cloning) warnings.push(`${doc.name}: cloning record was migrated`);
    if (d.assembly && d.assembly.schema !== CURRENT_SCHEMAS.assembly) warnings.push(`${doc.name}: assembly record was migrated`);
    return doc;
  });
  const ids = new Set();
  for (const doc of docs) {
    if (ids.has(doc.id)) { doc.id = newId(); warnings.push('A duplicate document ID was replaced'); }
    ids.add(doc.id);
  }
  return {
    docs,
    settings: clone(record.settings ?? {}),
    alignments: clone(record.alignments ?? []),
    traces: (record.traces ?? []).map((e) => ({ docId: e.docId, traces: (e.traces ?? []).map(decodeTrace) })),
    active: Math.max(0, Math.min(Number(record.active) || 0, Math.max(0, docs.length - 1))),
    name: record.name ?? docs[0]?.name ?? 'project',
    created: record.created ?? null,
    warnings,
    schema: PROJECT_SCHEMA,
  };
}

export function projectDigest(record) {
  return contentHash((record.docs ?? []).map((d) => `${d.id}:${d.seq}`).join('|'));
}

// Validate alignment records and give them a schema.
export function alignmentRecord(value) {
  const record = migrateRecord(value, 'alignment');
  if (!record) return null;
  return { ...record, schema: ALIGNMENT_SCHEMA };
}
