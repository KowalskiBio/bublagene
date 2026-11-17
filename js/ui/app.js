// @ts-nocheck — legacy DOM wiring; the engines it calls are type-checked.
// BublaGene: the entry point. The work lives in the modules beside this one:
//
//   core.js        page state, undo history, and the registered renderer
//   documents.js   open, close, edit and undo
//   typing.js      typing into the sequence view
//   clipboard.js   copying, and selection to primer or feature
//   files.js       opening and saving
//   persist.js     the session kept between visits
//   api/           what the inspector panels may do
//   render.js      redraws the page from the state
//   chrome.js      tabs, view options, status bar; viewpane.js the main view
//   cards.js       welcome, restore and empty-document cards
//   contextmenu.js the right-click menu
//   events.js      every listener

import './render.js'; // registers the renderer
import { render } from './core.js';
import { bindEvents } from './events.js';
import { offerPreviousSession } from './persist.js';

try {
  const t = localStorage.getItem('bublagene.theme');
  if (t) document.documentElement.dataset.theme = t;
} catch { /* ignore */ }

bindEvents();
render();
offerPreviousSession();
