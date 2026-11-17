// @ts-nocheck — legacy DOM wiring; the engines it calls are type-checked.
// BublaGene: the object handed to every inspector panel. Panels call these
// methods instead of reaching into the page state, so what a panel can do is
// listed in the files beside this one.

import { editApi } from './edit.js';
import { enzymesApi } from './enzymes.js';
import { designApi } from './design.js';
import { analysisApi } from './analysis.js';

// Methods call each other (api.setEnzymeSet, api.refreshLibrary, ...), so each
// group is built with the finished object in hand.
const api = {};
Object.assign(api, editApi(api), enzymesApi(api), designApi(api), analysisApi(api));

export { api };
