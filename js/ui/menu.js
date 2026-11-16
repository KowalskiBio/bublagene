// A small right-click menu.
// items: [{ label, hint, disabled, danger, onClick } | '-']

import { h } from './dom.js';

let open = null;

export function closeMenu() {
  open?.remove();
  open = null;
}

export function menuIsOpen() {
  return !!open;
}

export function showMenu(x, y, items) {
  closeMenu();
  const menu = h('div', { class: 'ctxmenu', role: 'menu' });
  let lastWasSep = true;
  for (const it of items) {
    if (it === '-') {
      if (!lastWasSep) { menu.append(h('div', { class: 'sep' })); lastWasSep = true; }
      continue;
    }
    if (!it) continue;
    lastWasSep = false;
    menu.append(h('button', {
      class: `item${it.danger ? ' danger' : ''}`, role: 'menuitem', type: 'button', disabled: !!it.disabled,
      onclick: () => { closeMenu(); it.onClick?.(); },
    }, h('span', {}, it.label), it.hint ? h('kbd', {}, it.hint) : null));
  }
  if (!menu.childElementCount) return;
  document.body.append(menu);
  const r = menu.getBoundingClientRect();
  menu.style.left = `${Math.max(6, Math.min(x, window.innerWidth - r.width - 6))}px`;
  menu.style.top = `${Math.max(6, Math.min(y, window.innerHeight - r.height - 6))}px`;
  open = menu;
  menu.querySelector('.item:not([disabled])')?.focus({ preventScroll: true });
}

// one set of dismissers for the life of the page
addEventListener('pointerdown', (e) => { if (open && !open.contains(e.target)) closeMenu(); }, true);
addEventListener('keydown', (e) => {
  if (!open) return;
  if (e.key === 'Escape') { e.preventDefault(); closeMenu(); return; }
  const items = [...open.querySelectorAll('.item:not([disabled])')];
  const i = items.indexOf(document.activeElement);
  if (e.key === 'ArrowDown') { e.preventDefault(); items[(i + 1) % items.length]?.focus(); }
  if (e.key === 'ArrowUp') { e.preventDefault(); items[(i - 1 + items.length) % items.length]?.focus(); }
}, true);
addEventListener('resize', closeMenu);
addEventListener('wheel', () => closeMenu(), { passive: true });
