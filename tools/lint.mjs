#!/usr/bin/env node
// Syntax lint: parse every source, test and tool file with `node --check`.
import { readdirSync, statSync } from 'node:fs';
import { join } from 'node:path';
import { spawnSync } from 'node:child_process';

const roots = ['js', 'tests', 'tools'];
const files = [];
const walk = (dir) => {
  for (const name of readdirSync(dir)) {
    const p = join(dir, name);
    if (statSync(p).isDirectory()) { if (name !== 'node_modules') walk(p); }
    else if (/\.(m?js)$/.test(name)) files.push(p);
  }
};
for (const r of roots) walk(r);
let failed = 0;
for (const f of files) {
  const r = spawnSync(process.execPath, ['--check', f], { encoding: 'utf8' });
  if (r.status !== 0) { failed++; console.error(`FAIL ${f}\n${r.stderr}`); }
}
console.log(`${files.length - failed}/${files.length} file(s) parsed`);
process.exit(failed ? 1 : 0);
