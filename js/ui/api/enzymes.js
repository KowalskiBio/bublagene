// @ts-nocheck — legacy DOM wiring; the engines it calls are type-checked.
// The panels' view of the app: Enzymes, digests and gels, primers, and mutagenesis.

import { newId } from '../../model.js';
import { designMutagenesis } from '../../mutagenesis.js';
import { current, render, state } from '../core.js';
import { edit } from '../documents.js';
import { fmt, toast } from '../dom.js';
import { loadEnzymeSets, storeEnzymeSets } from '../local-libraries.js';
import { rememberSoon } from '../persist.js';

export function enzymesApi(api) {
  return {
    toggleEnzyme(name) {
      state.shownEnzymes.has(name) ? state.shownEnzymes.delete(name) : state.shownEnzymes.add(name);
      render();
    },
    showEnzymes(names) { state.shownEnzymes = new Set(names); render(); },
    setEnzymeSet(names) {
      state.enzymeSet = names;
      // an enzyme that left the set cannot be hidden from the Enzymes panel any more, so stop drawing it
      state.shownEnzymes = new Set([...state.shownEnzymes].filter((n) => names.includes(n)));
      state.digestWith = state.digestWith.filter((n) => names.includes(n));
      render();
    },
    enzymeSets() { return loadEnzymeSets(); },
    toggleEnzymeInSet(name) {
      const has = state.enzymeSet.includes(name);
      api.setEnzymeSet(has ? state.enzymeSet.filter((n) => n !== name) : [...state.enzymeSet, name]);
    },
    setEnzymeEdit(on) { state.enzymeEdit = on; render(); },
    setEnzymeFilter(text) { state.enzymeFilter = text; render(); },
    saveEnzymeSet(name) {
      const clean = String(name ?? '').trim();
      if (!clean || /^(common|all)\b/i.test(clean)) { toast('Give the set a name other than Common or All', 'error'); return; }
      if (!state.enzymeSet.length) { toast('Pick at least one enzyme first', 'error'); return; }
      const sets = loadEnzymeSets();
      sets[clean] = [...state.enzymeSet];
      if (!storeEnzymeSets(sets)) { toast('This browser would not keep the set', 'error'); return; }
      toast(`Saved “${clean}” (${state.enzymeSet.length} enzymes)`);
      render();
    },
    deleteEnzymeSet(name) {
      const sets = loadEnzymeSets();
      if (!(name in sets)) return;
      delete sets[name];
      storeEnzymeSets(sets);
      render();
    },
    toggleDigest(name) {
      state.digestWith = state.digestWith.includes(name) ? state.digestWith.filter((x) => x !== name) : [...state.digestWith, name];
      render();
    },
    // quiet: the caller redraws the gel itself, so a slider being dragged is not rebuilt under the pointer
    setAgarose(v, quiet = false) { state.agarose = v; if (quiet) rememberSoon(); else render(); },
    setLadder(v) { state.ladder = v; render(); },
    // ---- gel, Tm, display ----
    addGelLane(label, enzymes) {
      const doc = current();
      if (!doc) return;
      if (state.gelLanes.length >= 8) { toast('The gel has room for eight saved lanes', 'error'); return; }
      state.gelLanes = [...state.gelLanes, { id: newId(), docId: doc.id, label, enzymes: [...enzymes] }];
      render();
    },
    removeGelLane(id) { state.gelLanes = state.gelLanes.filter((l) => l.id !== id); render(); },
    setTm(patch) { state.tm = { ...state.tm, ...patch }; render(); },
    setPrimerPair(patch) { state.primerPair = { ...state.primerPair, ...patch }; render(); },
    setMutate(patch) { state.mutate = { ...state.mutate, result: null, error: null, ...patch }; render(); },
    runMutagenesis(change) {
      const doc = current();
      if (!doc || !change) return;
      try {
        const result = designMutagenesis(doc, change, { method: state.mutate.method, targetTm: Number(state.mutate.targetTm) || null });
        state.mutate = { ...state.mutate, result, error: null };
      } catch (e) {
        state.mutate = { ...state.mutate, result: null, error: e.message };
      }
      render();
    },
    openMutant() {
      const r = state.mutate.result;
      if (!r) return;
      api.openDoc(r.mutant);
      toast(`${r.mutant.name}: ${fmt(r.mutant.seq.length)} bp`);
    },
    addMutagenicPrimers() {
      const r = state.mutate.result;
      if (!r) return;
      const keep = state.mutate;
      edit(`Add primers for ${r.label}`, (d) => {
        for (const p of [r.forward, r.reverse]) d.primers.push({ id: newId(), name: p.name, seq: p.seq, description: `Mutagenic primer for ${r.label}`, extra: {} });
      });
      state.mutate = keep;
      state.panel = 'Primers';
      render();
    },
  };
}
