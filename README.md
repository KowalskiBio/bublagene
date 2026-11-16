# BublaGene

**[Try it here](https://kowalskibio.github.io/bublagene/)** — it runs entirely in your browser.

A browser-based sequence viewer and editor for `.dna`, GenBank and FASTA
files. Everything runs locally; no file leaves the machine.

## Run it

```bash
python3 -m http.server 8642 --directory bublagene
```

Then open <http://localhost:8642>. Drop a sequence file on the page, or start
from the built-in example.

- **Map** — circular or linear, features as arrows stacked into lanes, enzyme
  sites labelled with their cut positions, click anything to select it.
- **Sequence** — both strands, feature bars, amino acids under their codons,
  cut marks. Drag across the bases to select a region and press ⌘C to copy it;
  only the lines in view are built, so a 35 kb plasmid stays responsive.
- **Right-click** anywhere in the map or sequence for everything you can do
  with a region: copy it as bases, reverse complement, protein or FASTA; turn
  it into a forward or reverse primer; design a pair; add a feature; delete or
  replace the bases; or make it the start of a circular sequence. Right-clicking
  a feature adds its own actions on top.
- **Sequencing data** — reads stored in a `.dna` file are decoded and drawn
  under the reference, matches dimmed to dots so mismatches and insertions stand
  out. The BAM reader is checked field by field against samtools. Trace
  chromatograms are decoded too and drawn as the four dye curves with the called
  bases above them.
- **Panels** — features (add, edit, delete, ORFs), primers (binding sites, pair
  design, PCR), enzymes, digest with a simulated gel, and file details.
- **Cloning** — cut a vector and an insert with the same enzymes, check that the
  ends actually ligate (and whether the insert also fits backwards), and join
  them into a new construct carrying features from both pieces. Gibson assembly
  joins fragments that share terminal overlaps. Products keep a structured
  cloning record in portable JSON.
- Undo/redo, several files open at once, light and dark themes.
- **Your files come back.** What you had open is kept in the browser's own
  storage, so closing the tab and returning offers to reopen it, listing the
  files and whether any had unsaved changes. "Start fresh" forgets it. Nothing
  is uploaded, and the stored copy keeps a `.dna` file's untouched bytes, so a
  reopened file still saves back correctly.

There is also a command line, described below.

## Requirements

Node 20 or newer. On this machine the default `node` (25.9) is broken — it fails
with a missing `libsimdjson` library after a Homebrew upgrade — so use the
Node 24 installed alongside it:

```bash
export PATH=/opt/homebrew/opt/node@24/bin:$PATH
```

`brew reinstall node` should repair the default one.

## Run the tests

```bash
cd bublagene && npm test
```

37 tests. Most run on their own; the `.dna` tests additionally read every
`.dna` file in a folder you point them at with `DNA_SAMPLES`:

```bash
DNA_SAMPLES=~/my-plasmids npm test
```

Without samples those tests report themselves as skipped and the rest still run.

What they check: reverse complement, translation, melting temperature, IUPAC
searching across the origin, enzyme sites (including enzymes that cut outside
their recognition site), digests, edits that keep features on the same bases,
origin rotation, reverse complementing, ORF finding, primer binding, PCR,
restriction cloning and Gibson assembly, XML, GenBank and FASTA round trips,
and `.dna` round trips over all 206 sample files. One test translates all 877
coding sequences in those files and compares each against the protein stored in
the file itself; 874 agree exactly and the three that don't are listed in the
test with the reason.

## Use it

```bash
node tools/bg.mjs info plasmid.dna
```

| Command | What it does |
|---|---|
| `bg info <file>` | summary, features with protein lengths, primers with binding sites, and the raw segment list for `.dna` files |
| `bg seq <file> [start..end]` | the sequence, numbered, optionally for one region |
| `bg cut <file> [enzymes] [digest-with]` | which enzymes cut where, unique cutters first; add a second list to see fragment sizes and ends |
| `bg orfs <file> [minCodons]` | open reading frames on both strands |
| `bg primers <file> <start..end> [Tm]` | design a primer pair for a region |
| `bg pcr <file> <fwd> <rev>` | simulate PCR and describe the product |
| `bg clone <vector> <insert> <enzyme,...> [out]` | cut both, ligate the insert, describe or save the construct |
| `bg gibson <f1> <f2> [...] [-o out] [--linear]` | assemble fragments by shared terminal overlaps |
| `bg history <file>` | print the cloning record stored beside a file |
| `bg convert <file> <out.gb\|.fasta\|.dna>` | convert between formats |

Positions are 1-based and inclusive, as in GenBank. Files may be `.dna`,
`.gb`/`.gbk`/`.genbank`, or `.fasta`/`.fa`.

Examples:

```bash
node tools/bg.mjs cut plasmid.dna EcoRI,BamHI,HindIII,NotI
```

```bash
node tools/bg.mjs primers plasmid.dna 1653..2372
```

```bash
node tools/bg.mjs clone vector.dna insert.dna SalI,KpnI construct.dna
```

## Cloning

`js/cloning.js` turns digested fragments into constructs. A fragment from
`digest()` becomes a piece: its ends are checked, the insert is tried in both
orientations, and the product carries every feature that survives the cuts —
split at the junction.

```js
import { checkLigation, restrictionClone } from './js/cloning.js';

const backbone = { doc: vector, fragment: vectorFragments[0] };
const insert = { doc: insertDoc, fragment: insertFragments[1] };

const check = checkLigation({ backbone, insert });   // orientations, junctions, size
const { doc } = restrictionClone({ backbone, insert });
```

Gibson assembly uses `gibsonAssemble([doc1, doc2, …])`, which infers the order
and orientation from shared terminal overlaps (15 bp by default).

Products keep `doc.cloning` (operation, parents, fragments, junctions and
regenerated sites). `cloningRecord(doc)` returns it as plain JSON under the
schema `bublagene/cloning/1`, and the CLI writes it beside the construct as
`<name>.cloning.json`; `bg history <file>` prints it. The product's Notes
description also records the parents and enzymes, so a saved `.dna` keeps a
readable provenance line. Compressed history blocks in existing files are
passed through untouched and are not decoded.

## Working with .dna files

Reading keeps every segment of the original file. Saving regenerates only the
parts BublaGene understands (sequence, features, primers, notes, sequence
properties) and passes everything else through byte for byte, so the history
blocks, display options and alignments in a file survive a round trip. Caches
the format rebuilds are dropped once the sequence has been edited, so they can
never go stale. A segment that is xz-compressed is passed through untouched
rather than guessed at, and `info` reports it.

## Layout

```
js/seq.js              sequence maths: complement, translation, Tm, ORFs, IUPAC
js/model.js            the document model and edits that keep features in place
js/xml.js              a small XML parser and writer
js/enzymes.js          ~90 restriction enzymes, site finding, digests, gel model
js/primers.js          primer binding, design, PCR
js/cloning.js          restriction cloning, Gibson assembly, cloning history
js/formats/dna.js      binary .dna
js/formats/genbank.js  GenBank and FASTA
js/ui/dom.js           DOM helpers (for the interface that isn't built yet)
tools/bg.mjs           the command line
```

## Not built

- The digest/gel and primer panels work but have had less testing than the map and sequence views.

## Licence

AGPL-3.0-or-later. See [LICENSE](LICENSE). You may use, study, change and share
this; anything you distribute or run as a network service that is built on it
must be offered under the same licence, with its source available.
