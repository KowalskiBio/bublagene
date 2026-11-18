// @ts-nocheck — legacy DOM wiring; the engines it calls are type-checked.
// Opening and saving files and projects.

import { readDna, writeDna } from '../formats/dna.js';
import { readFasta, readGenBank, writeFasta, writeGenBank } from '../formats/genbank.js';
import { newId } from '../model.js';
import { parsePdb } from '../pdb.js';
import { findBindingSites } from '../primers.js';
import { projectFromJSON, projectToJSON } from '../project.js';
import { openDocument, saveDocument } from '../storage.js';
import { collectAttached, restoreAttached } from './attached.js';
import { current, render, state } from './core.js';
import { openDoc, openTraceFile } from './documents.js';
import { traceFormat } from '../traces.js';
import { toast } from './dom.js';

export async function readFileData(file) {
  if (file.data) return file.data instanceof Uint8Array ? file.data : new Uint8Array(file.data);
  return new Uint8Array(await file.arrayBuffer());
}

export async function readFileText(file) {
  if (file.text) return file.text();
  return new TextDecoder().decode(await readFileData(file));
}

const TRACE_EXT = /\.(ab1|abi|scf|ztr)$/i;

export async function openFiles(files) {
  // references first, so a read opened in the same batch has something to align to
  const ordered = [...files].sort((a, b) => TRACE_EXT.test(a.name) - TRACE_EXT.test(b.name));
  for (const file of ordered) {
    try {
      const lower = file.name.toLowerCase();
      if (lower.endsWith('.pdb') || lower.endsWith('.ent')) {
        const structure = parsePdb(await readFileText(file), { name: file.name.replace(/\.(pdb|ent)$/i, '') });
        state.structures.push({ id: newId(), structure });
        state.protein = { ...state.protein, structure: state.structures[state.structures.length - 1].id };
        state.panel = 'Protein';
        render();
        toast(`${structure.name}: ${structure.residues} residues in ${structure.chains.length} chain${structure.chains.length === 1 ? '' : 's'}`);
      } else if (lower.endsWith('.dna')) {
        const doc = readDna(await readFileData(file), file.name);
        openDoc(doc);
        for (const w of doc.source.warnings ?? []) toast(w);
      } else if (TRACE_EXT.test(lower)) {
        await openTraceFile(await readFileData(file), file.name);
      } else if (lower.endsWith('.json')) {
        const project = projectFromJSON(await readFileText(file));
        for (const d of project.docs) openDoc(d);
        if (project.warnings.length) toast(project.warnings[0]);
      } else {
        // a chromatogram with an unexpected extension is still recognised by its first bytes
        if (file.data || typeof file.arrayBuffer === 'function') {
          const bytes = await readFileData(file);
          if (traceFormat(bytes)) { await openTraceFile(bytes, file.name); continue; }
        }
        const text = (await readFileText(file)).replace(/^\uFEFF/, '');
        const docs = /^LOCUS/m.test(text) ? readGenBank(text, file.name)
          : text.trimStart().startsWith('>') ? readFasta(text, file.name)
            : readFasta(`>${file.name}\n${text}`, file.name);
        for (const d of docs) openDoc(d);
        const skipped = docs.flatMap((d) => d.source?.warnings ?? []);
        if (skipped.length) toast(`${file.name}: ${skipped.length} feature${skipped.length === 1 ? '' : 's'} skipped, location not understood (${skipped[0].slice(0, 50)})`, 'error');
      }
    } catch (e) {
      toast(`${file.name}: ${e.message}`, 'error');
    }
  }
}

export async function save(kind) {
  const doc = current();
  if (!doc) return;
  let data;
  let name = doc.name.replace(/[^\w.-]+/g, '_');
  let types;
  if (kind === 'dna') {
    data = writeDna(doc, { bindingSitesFor: (p) => findBindingSites(doc.seq, doc.circular, p.seq) });
    name += '.dna';
    types = [{ description: 'BublaGene sequence', accept: { 'application/octet-stream': ['.dna'] } }];
  } else if (kind === 'gb') {
    data = writeGenBank(doc);
    name += '.gb';
    types = [{ description: 'GenBank', accept: { 'chemical/seq-na-genbank': ['.gb', '.gbk'] } }];
  } else {
    data = writeFasta(doc);
    name += '.fasta';
    types = [{ description: 'FASTA', accept: { 'text/plain': ['.fasta', '.fa'] } }];
  }
  try {
    const result = await saveDocument({ data, name, types });
    toast(result.method === 'download' ? `Saved ${result.name}` : `${result.name} saved`);
  } catch (e) {
    toast(`Could not save: ${e.message}`, 'error');
  }
}

export async function saveProject() {
  if (!state.docs.length) return;
  // every read and chromatogram opened against a document, whether the reader attached it or a .dna file carried it
  const kept = collectAttached({ plain: true, onlyAttached: false });
  const json = projectToJSON(state.docs, {
    active: state.active,
    settings: {
      view: state.view, panel: state.panel, enzymeSet: state.enzymeSet,
      shownEnzymes: [...state.shownEnzymes], digestWith: state.digestWith,
      agarose: state.agarose, ladder: state.ladder,
      showTranslation: state.showTranslation, showSixFrame: state.showSixFrame, showComplement: state.showComplement,
      featureFilter: state.featureFilter,
    },
    alignments: kept.alignments,
    traces: kept.traces,
  });
  try {
    const result = await saveDocument({ data: json, name: `${state.docs[state.active]?.name ?? 'project'}${state.docs.length > 1 ? `+${state.docs.length - 1}` : ''}.bublagene.json`, types: [{ description: 'BublaGene project', accept: { 'application/json': ['.json'] } }] });
    toast(`${result.method === 'download' ? 'Downloaded' : 'Saved'} ${state.docs.length} tab${state.docs.length === 1 ? '' : 's'}${kept.traces.length ? ' with their chromatograms' : ''}`);
  } catch (e) { toast(`Could not save the project: ${e.message}`, 'error'); }
}

export async function openProject() {
  try {
    const file = await openDocument({ types: [{ description: 'BublaGene project', accept: { 'application/json': ['.json'] } }] });
    if (!file) return;
    const project = projectFromJSON(new TextDecoder().decode(file.data));
    const first = state.docs.length;
    // opening the same project twice must not give two tabs the same identity
    const open = new Set(state.docs.map((d) => d.id));
    const renamed = new Map();
    for (const d of project.docs) if (open.has(d.id)) { const was = d.id; d.id = newId(); renamed.set(was, d.id); }
    const remap = (e) => ({ ...e, docId: renamed.get(e.docId) ?? e.docId });
    restoreAttached({ alignments: project.alignments.map(remap), traces: project.traces.map(remap) });
    for (const d of project.docs) openDoc(d);
    // reopen on the tab that was in front when it was saved
    state.active = Math.min(first + project.active, state.docs.length - 1);
    const s = project.settings ?? {};
    state.view = s.view ?? state.view;
    state.panel = s.panel ?? state.panel;
    state.enzymeSet = s.enzymeSet ?? state.enzymeSet;
    state.shownEnzymes = new Set(s.shownEnzymes ?? state.shownEnzymes);
    state.digestWith = s.digestWith ?? state.digestWith;
    state.showTranslation = s.showTranslation ?? state.showTranslation;
    state.showSixFrame = s.showSixFrame ?? state.showSixFrame;
    state.showComplement = s.showComplement ?? state.showComplement;
    for (const w of project.warnings) toast(w);
    toast(`Opened ${project.docs.length} tab${project.docs.length === 1 ? '' : 's'}`);
    render();
  } catch (e) { toast(`Could not open the project: ${e.message}`, 'error'); }
}
