# BublaGene

![CI](https://github.com/KowalskiBio/bublagene/actions/workflows/ci.yml/badge.svg)
![Licence: AGPL-3.0-or-later](https://img.shields.io/badge/licence-AGPL--3.0--or--later-blue)

**A local-first DNA workbench for the browser and the command line.**

BublaGene reads `.dna`, GenBank and FASTA files and gives you maps, a typeable
sequence view, editable annotations, automatic feature detection, primers and
mutagenesis, PCR, six cloning methods, assembly, alignment, variants,
chromatograms, CRISPR guides, RNA folding and protein properties. Everything
runs on your own machine — no sequence is uploaded.

**[Open BublaGene](https://kowalskibio.github.io/bublagene/)** · **[Source](https://github.com/KowalskiBio/bublagene)**

---

## Why it exists

Most sequence tools are either a viewer with no analysis or a heavy desktop
package. BublaGene keeps the useful engines small, testable and scriptable:
the browser and `bg` call the same JavaScript modules, and every record it
writes has a stable JSON schema.

- **Local-first.** Files stay in the browser or on your disk.
- **Round-trip safe.** Unknown GenBank qualifiers, fuzzy locations, segments,
  `.dna` pass-through bytes and unsupported records survive a save.
- **Scriptable.** Every engine has a CLI command and a JSON report.
- **Honest about models.** CRISPR off-target scores and RNA folding are
  explicitly bounded; they are not sold as predictions they cannot make.

## Quick start

### In the browser

```bash
git clone https://github.com/KowalskiBio/bublagene.git
cd bublagene
npm run serve
```

Open <http://localhost:8642>. Drop a file on the page or start from the built-in
example.

No build step, no bundler, no server-side component. (If you edit the code while
serving it, hard-reload: browsers cache ES modules aggressively.) Node 20 or newer is needed
for the tests and the CLI; the browser app runs from static files.

### Command line

```bash
node tools/bg.mjs info plasmid.dna
node tools/bg.mjs cut plasmid.dna EcoRI,BamHI,HindIII
node tools/bg.mjs primers plasmid.dna 1653..2372
node tools/bg.mjs pcr plasmid.dna <forward> <reverse> --json
```

`npm link` once and the command is just `bg`.

## What you can do

The browser app is the whole toolbox; the [command line](#the-command-line)
covers the engines that make sense without a screen. Panels are grouped in the
inspector: **Sequence** (Features, Detect, Plot, Protein), **Primers** (Primers,
Mutate, CRISPR), **Enzymes** (Enzymes, Digest), **Cloning**, **Analyse**
(Alignment, Contigs, Search, Structure) and **File** (File, Library, History).

### Maps, sequence and editing

- Circular and linear maps with stacked feature lanes, enzyme sites and cut
  positions. Force either style for any file, or show the map above the
  sequence (**Map + sequence**).
- Sequence view with both strands, feature bars, codon-aligned amino acids and
  a six-frame translation track. Optional dashed bars draw open reading frames
  and every place a primer anneals, on the map and in the sequence.
- Click to put a caret and type bases; Backspace, Delete, the arrow keys (with
  Shift to select), Home/End and paste all edit at the caret. Drag to select;
  right-click to copy as bases, reverse complement, protein or FASTA, make
  primers, add a feature, delete or replace bases, or set a new origin.
- Feature editing covers name, type, strand, colour and arbitrary qualifiers.
  Features can be multi-selected, reordered by drag, and copied, cut and pasted
  as annotations (positions kept relative to each other).
- Fragments copy and paste across documents. Undo and redo coalesce bursts of
  typing, and the **History** panel lists every step so you can jump to any of
  them. Autosave keeps an IndexedDB session plus a small recovery mirror.
- **Export** the map as SVG, PNG or a one-page PDF (always dark-on-white), or
  **Print** a sheet built for paper: the map, a feature table, primers and the
  whole numbered sequence, however little of it is on screen.

### Annotation and analysis

- **Detect** scans both strands for known elements: promoters, epitope and
  purification tags, protease sites, recombination sites and common primer
  sites. Tags are found in any reading frame whatever the codons; DNA elements
  tolerate a chosen percentage of mismatches and wrap the origin of a circle.
  It can also learn from the features annotated in your other open files, and
  you can save any selection as your own element.
- **Plot** shows GC content, GC skew and cumulative skew along the sequence;
  click the plot to select that stretch.
- **Protein** gives mass, pI, charge, hydropathy and a motif scan (tags,
  cleavage sites, N-glycosylation sequons, possible membrane helices) for any
  CDS, and a rotatable 3D backbone for a PDB file you open. Structures are not
  predicted.
- **Search** looks through this file, all open files or your saved collections,
  by sequence (either strand, IUPAC codes), feature, primer or note. Links to
  NCBI BLAST and Addgene open the other site in a new tab, and send it your
  query only when you click.
- **Library** keeps files in this browser, grouped into named collections and
  searchable without opening them. Share a collection as a project file.

### Primers, PCR and mutagenesis

- SantaLucia nearest-neighbour Tm with salt/Mg correction, GC, 3′ clamp and
  degenerate-primer expansion. Set the Na⁺, Mg²⁺ and primer concentrations of
  your own buffer.
- Hairpin, self-dimer and cross-dimer ΔG for every primer and any pair, with
  warnings.
- Ranked primer pairs with product-size, Tm-window, specificity, tail and
  restriction-site options.
- Mismatch-aware PCR with a strict 3′ rule, circular templates, primer-binding
  features on the product, and a simulated gel.
- **Mutate**: pick an amino acid change in a CDS (the fewest-base codon is
  chosen) or edit the selected bases, then get overlapping (QuikChange-style) or
  back-to-back (Q5 SDM-style) primers and the expected mutant, with insertions
  and deletions supported.

### Enzymes and digests

- 94 curated enzymes with cut offsets, overhang type and length, and
  isoschizomers. Custom enzyme sets can be edited and saved in the browser, or
  exchanged as plain JSON (`bublagene/enzymes/1`).
- Unique cutters, cut-count filters, per-enzyme fragment tables, and a Type IIS
  overhang compatibility matrix.
- Methylation: Dam, Dcm and CpG sensitivity is judged from the sequence around
  each site, so an enzyme is flagged only where the motif actually overlaps.
  Only well-established cases are listed.
- A simulated gel with a ladder, adjustable agarose, and several lanes at once:
  keep a digest as a lane and try another, or compare digests of different
  files.

### Cloning

The Cloning panel has six methods, all working from the open documents and
opening the product as a new tab:

- **Restriction + ligation**, with a picture of the digested pieces and their
  real sticky ends, the ligation verdict for both orientations, and a map of the
  product.
- **Golden Gate** with a Type IIS enzyme, the overhang each file leaves, and a
  junction-by-junction match check.
- **Gibson** and **In-Fusion** (the same homology assembly, with a 15 bp default
  overlap for In-Fusion), circular or linear: up to eight fragments in the
  panel, twelve with `bg assemble`.
- **Site-specific recombination**: Gateway LR and BP, Cre/lox and Flp/FRT. The
  stretch between two matching sites is swapped through their shared core, using
  the sites annotated on the files (attB1 and attB2 are also found by sequence);
  the by-product is built too.
- **TOPO, TA and blunt** insertion at your selection or caret, in either
  orientation; directional cloning checks for the CACC end.

Every construct keeps a cloning record, shown as a tree: the product on top,
the pieces beneath, and, when those are open and were cloned themselves, what
they were made from in turn. Click a node to jump to that file. Gibson and
Golden Gate assembly also write a replayable plan (`bublagene/assembly/1`), and
a fragment library can be batch-assembled from the command line.

### Alignment, variants and contigs

- Affine-gap global, local and banded pairwise alignment with identity and
  gap statistics, dot plots, and mismatch/indel colouring.
- Multiple alignment of several files (centre-star) with a difference table.
- **Contigs**: quality-weighted assembly from Sanger reads in either
  orientation, using chromatogram base calls and qualities where a file has
  them; open a contig as a file with its reads annotated.
- Normalised SNPs, indels and MNVs with reference-anchored coordinates that
  round-trip through `variation` features.

### Traces

- ABI (`.ab1`, `.abi`), SCF and ZTR chromatograms with four dye channels, called
  bases, peak positions and Phred qualities.
- Open a read with a reference (`.gb`, `.dna`, FASTA) already open and it is
  aligned to that reference in either orientation: it appears as a row under the
  sequence (dots for matches, letters for mismatches, a marker for insertions)
  and its chromatogram is behind the Chromatogram view. A read that does not
  match the reference opens as its own document. Files are recognised by their
  first bytes, so the extension does not matter.
- The **Alignment** tab shows the reference as a feature map, then each aligned
  sequence column by column under the matching reference bases (mismatches
  filled in, gaps and insertions marked), then the chromatogram of the selected
  one stretched onto the same columns, so each peak sits under the base it
  called. Align a document to another from the Alignment panel and it opens
  there too; when that document is a trace, the chromatogram comes with it.
- Mott and Phred trimming, trim features and synchronised trace-to-reference
  navigation.
- Trace bases map back to the reference through an alignment; “Show in sequence”
  jumps to that base.

### CRISPR and RNA

- PAM-aware guides for SpCas9, SaCas9, Cas12a or a custom PAM, on both strands,
  with GC bounds.
- Labelled seed/mismatch off-target heuristic — **not** a cleavage prediction.
- Deterministic NHEJ and HDR simulation with frame-effect reports and edit
  primers.
- Bounded RNA minimum-free-energy folding over a documented Turner subset, with
  dot-bracket and SVG output. Sequences past the limit are refused, not silently
  truncated.

### Formats and persistence

- GenBank, FASTA, GFF3, BED, CSV feature tables, plain-text maps.
- SVG, PNG and PDF diagrams (from the browser and the command line), and PDB
  files for the 3D backbone view.
- Project archives (`bublagene/project/1`) keep open documents, features,
  primers, cloning/assembly records, alignments and settings with stable IDs.
- Native records carry schema versions; older `project`, `cloning`, `assembly`
  and `alignment` records migrate forward on read.

## The command line

Positions are 1-based and inclusive, as in GenBank. Files may be `.dna`,
`.gb`/`.gbk`/`.genbank`, or `.fasta`/`.fa`.

### Inspect and edit

| Command | What it does |
|---|---|
| `bg info <file>` | summary, features, primers, protein lengths and `.dna` segments |
| `bg seq <file> [start..end]` | numbered sequence, optionally one region |
| `bg composition <file>` | base counts, GC and length |
| `bg search <file> <query>` | sequence, feature and qualifier search |
| `bg frames <file> [start..end]` | six-frame translation |
| `bg translate <file> [start..end] [table]` | translate with an NCBI genetic code |
| `bg codons` | list and validate supported genetic codes |
| `bg annotate <file> [out]` | recompute `translation` / `gc_content` after edits |
| `bg orfs <file> [minCodons]` | open reading frames on both strands |
| `bg export <file> <out>` | GFF3, BED, CSV, text map, SVG, PNG or PDF |
| `bg convert <file> <out.gb\|.fasta\|.dna>` | convert between formats |

### Primers and PCR

| Command | What it does |
|---|---|
| `bg primers <file> <start..end> [Tm]` | design a primer pair |
| `bg design <file> <start..end> [flags]` | Tm/product-size/tail/restriction options |
| `bg pcr <file> <fwd> <rev> [--json]` | simulate PCR and describe the product |
| `bg gel <bp,bp,...> [flags]` | fragment sizes on a simulated gel |

### Enzymes and cloning

| Command | What it does |
|---|---|
| `bg cut <file> [enzymes] [digest-with]` | which enzymes cut where; optional fragment sizes |
| `bg digest <file> <enzyme,...>` | per-enzyme fragment tables |
| `bg overhangs <file> [enzymes]` | Type IIS overhangs and compatibility |
| `bg enzymes list\|export\|import` | curated metadata or custom JSON sets |
| `bg clone <vector> <insert> <enzyme,...> [out]` | restriction cloning |
| `bg gibson <f1> <f2> [...] [-o out] [--linear]` | overlap assembly |
| `bg assemble <f1> <f2> ...` | N-fragment homology assembly with a plan |
| `bg golden-gate <f1> <f2> ...` | Type IIS modular assembly |
| `bg assembly-primers <f1> <f2> ...` | PCR primers carrying assembly overlaps |
| `bg plan <file>` / `bg replay <plan.json>` | inspect or rebuild an assembly plan |
| `bg batch <library.json> <designs.json>` | many constructs from one library |
| `bg history <file>` | print the cloning record beside a file |

### Alignment, variants and traces

| Command | What it does |
|---|---|
| `bg align <a> <b> [--local] [--band N]` | affine-gap pairwise alignment |
| `bg dotplot <a> <b>` | dot-plot grid |
| `bg contig <r1> <r2> ...` | quality-aware contig assembly |
| `bg variants <ref> <query>` | normalised SNPs/indels between sequences |
| `bg trace <trace.ab1\|.scf>` | chromatogram summary and quality |
| `bg trim <trace> [flags]` | Mott/Phred trimming, optionally FASTA |

### CRISPR, structure and automation

| Command | What it does |
|---|---|
| `bg crispr <file>` | PAM-aware guide finding |
| `bg crispr-edit <file> <guide>` | deterministic NHEJ/HDR simulation |
| `bg crispr-primers <file> <guide>` | primers flanking the cut |
| `bg fold <file> [start..end]` | bounded MFE structure, dot-bracket or SVG |
| `bg report <file>` | comprehensive JSON report |
| `bg pipeline <file> <steps.json>` | reproducible JSON pipeline |
| `bg project save\|open ...` | lossless project archive with migrations |

### JSON pipelines

```bash
cat > steps.json <<'JSON'
{
  "seed": 7,
  "steps": [
    { "op": "info" },
    { "op": "digest", "enzymes": ["EcoRI", "BamHI"] },
    { "op": "fold" },
    { "op": "export", "out": "map.svg" }
  ]
}
JSON

node tools/bg.mjs pipeline plasmid.dna steps.json
```

Pipeline steps run in order and produce one JSON report. Supported operations:
`info`, `cut`, `digest`, `design`, `pcr`, `variants`, `align`, `crispr`, `fold`,
`translate`, `export`.

## Library examples

```js
import { translate, sixFrameTranslation, revcomp } from './js/seq.js';
import { findSites, digest, ENZYME_BY_NAME } from './js/enzymes.js';
import { designPrimerPairs, simulatePCR } from './js/primers.js';
import { gibsonAssemble } from './js/cloning.js';
import { callVariants, variantsToFeatures } from './js/variants.js';
import { detectFeatures, hitToFeature } from './js/detect.js';
import { aminoAcidChange, designMutagenesis } from './js/mutagenesis.js';
import { recombine } from './js/recombination.js';
import { multipleAlign } from './js/msa.js';

translate('ATGGCCTAA', 11);                 // NCBI table 11
sixFrameTranslation('ATGGCCTAA');
revcomp('GAATTC');

// These values come from readDna()/readGenBank()/readFasta() or the clipboard.
const cuts = findSites(seq, true, ENZYME_BY_NAME.EcoRI);
const fragments = digest(seq, true, cuts);

const pair = designPrimerPairs(seq, false, 100, 600, { targetTm: 60 });
const { doc: product } = simulatePCR(templateDoc, pair.forward.seq, pair.reverse.seq);

const { doc: construct } = gibsonAssemble([vector, insert], { circular: true });

const calls = callVariants(reference, edited);
doc.features.push(...variantsToFeatures(doc, calls.variants));

// Find known elements, and add the ones you accept as features.
const hits = detectFeatures(doc, undefined, { maxMismatchRate: 0.05 });
doc.features.push(...hits.map(hitToFeature));

// Design primers for A5E in the first CDS, and get the expected mutant.
const change = aminoAcidChange(doc, doc.features.find((f) => f.type === 'CDS'), 5, 'E');
const { forward, reverse, mutant } = designMutagenesis(doc, change, { method: 'overlap' });

// Swap the stretch between two att sites of an entry clone into a destination.
const { doc: expression } = recombine(entryClone, destination);

const msa = multipleAlign([{ name: 'a', seq: a }, { name: 'b', seq: b }]);
```

## How it is built

```
js/
  seq.js              sequence maths: complement, translation, Tm, ORFs, IUPAC
  codon-tables.js     25 NCBI genetic codes, validation, per-table translation
  model.js            document model, feature catalog, edits that keep features in place
  thermo.js           nearest-neighbour thermodynamics and secondary-structure ΔG
  primers.js          binding, ranked design, mismatch PCR, gel model
  mutagenesis.js      amino-acid and base changes, overlapping and back-to-back primers
  enzymes.js          enzyme dataset, custom sets, digests, overhangs
  methylation.js      Dam, Dcm and CpG sensitivity from sequence context
  cloning.js          restriction cloning, Gibson and In-Fusion assembly
  assembly.js         homology/Golden Gate assembly, plans, batch runs
  recombination.js    att, lox and FRT recombination; TOPO-style insertion
  detect.js           known-element detection; feature-library.js is its built-in library
  align.js            affine-gap global/local/banded alignment and dot plots
  msa.js              centre-star multiple alignment and difference table
  contig.js           quality-weighted contig assembly
  variants.js         normalised variant calling and feature round-trip
  crispr.js           guide finding, off-target heuristic, edit simulation
  fold.js             bounded RNA MFE folding
  protein.js          mass, pI, hydropathy and motif scan; pdb.js reads Cα traces
  tracks.js           GC content and skew along a sequence
  traces.js           trace trimming, slicing and reference mapping
  readalign.js        place a Sanger read on a reference (k-mer seed, local alignment)
  ui/alignview.js     the Alignment tab: reference map, aligned rows, chromatogram
  fragment.js         fragment and feature copy/paste, feature reordering
  cache.js            digest-keyed caches and stale-alignment invalidation
  undo.js             coalescing undo/redo with a jumpable timeline
  search.js           sequence/feature/qualifier search and filters
  collection.js       saved collections in IndexedDB, and searching them
  storage.js          File System Access, OPFS and download fallbacks
  compute.js          worker offload with local fallback
  project.js          project archives and schema migrations
  export.js           SVG, PNG and PDF diagrams (command line)
  pdf.js              a minimal one-page PDF writer
  xml.js              small XML parser/writer
  formats/            dna, genbank/fasta, gff3/bed/csv, ab1, scf, ztr, bam
  workers/            bounded heavy work
  ui/                 the browser app, below
tools/
  bg.mjs              command line
  bench.mjs           5 Mb performance budget
  lint.mjs            syntax lint
tests/                unit, property, fuzz, CLI and jsdom UI tests
```

The browser and CLI share the same engine modules. There is no duplicated
analysis logic.

### The browser app

`js/ui/app.js` is a short entry point; the work is split along the lines of what
changes together:

```
js/ui/
  app.js            entry: binds events, first render, offers the last session
  core.js           page state, undo histories, and the registered renderer
  documents.js      open, close, edit and undo
  typing.js         typing into the sequence view
  clipboard.js      copying, and selection to primer or feature
  files.js          opening and saving files and projects
  persist.js        the session kept between visits
  local-libraries.js  saved enzyme sets and detection elements
  viewdata.js       what the map and sequence draw (cut sites, ORFs, primer sites)
  render.js         the one function that redraws the page
  chrome.js         document tabs, view options, status bar, inspector tabs
  viewpane.js       the main view: map, sequence, both, or chromatogram
  cards.js          welcome, restore and empty-document cards
  contextmenu.js    the right-click menu
  events.js         buttons, keys, drag-and-drop, clipboard and print
  api/              what the inspector panels may do (edit, enzymes, design, analysis)
  panels.js, panels2.js, cloningpanel.js   the inspector panels
  map.js, seqview.js, traceview.js, proteinview.js   the drawn views
  export.js         map export and the print sheet
```

Everything that changes the page calls `render()` from `core.js`. The renderer
registers itself there rather than being imported, so the code that draws the
page and the code that changes it never import each other. Panels never touch
the state directly: they call methods on the `api` object, so what a panel can
do is listed in one place.

## Testing

```bash
npm test              # 166 tests
npm run lint          # parse every source file with node --check
npm run typecheck     # tsc --checkJs over engines and CLI
npm run test:coverage # line/function thresholds
npm run bench         # 5 Mb genome budget
npm run check         # lint + typecheck + tests
```

Without external fixtures the suite reports **156 pass, 10 skip, 0 fail**. The
skips are environment-gated: tests that need `.dna` files and tests that compare
the BAM reader with `samtools`.

Point the sample tests at a folder of `.dna` files:

```bash
DNA_SAMPLES=~/plasmids npm test
```

To run the samtools comparison too:

```bash
SAMTOOLS="$(which samtools)" DNA_SAMPLES=~/plasmids npm test
```

The suite includes property tests (fast-check), bounded fuzzing for the binary
and text parsers, a jsdom smoke test that boots the real browser app, a 5 Mb
genome benchmark, and three consecutive test runs in CI to catch flakes.

## Scope and limits

- Genetic-code, PAM and RNA-folding behaviour is limited to the models
  documented in their modules. CRISPR off-target scores are heuristics, not
  cleavage predictions.
- Folding refuses sequences past its stated length instead of truncating.
- The ABI/SCF readers interpret the common layouts, not every writer variant.
- Restriction cloning and assembly are deterministic calculators; they do not
  model ligation efficiency, star activity or cellular repair outcomes.
- Feature detection ships only short, well-known elements (promoters, tags,
  cleavage and recombination sites, primer sites). Long elements such as
  resistance genes and origins come from your own saved elements and from the
  features annotated in your other open files; none are guessed at.
- Recombination is a crossover inside the sites' shared core, using sites
  annotated on your files; only attB1 and attB2 are recognised by sequence.
- Multiple alignment is centre-star: good for variants and clones of one
  sequence, not for distant homologues.
- Methylation lists only well-established cases; an enzyme that is not flagged
  is not guaranteed insensitive.
- The protein view draws the Cα trace of a PDB file you open. It does not
  predict structure.
- BLAST and Addgene are links that open in a new tab when you click; nothing is
  searched online from inside the app, and there is no sync between devices
  (share a collection as a project file).
- `.dna` history blocks that use an unsupported compression are preserved
  byte-for-byte and left undecoded.

## Licence

AGPL-3.0-or-later. See [LICENSE](LICENSE). You may use, study, change and share
this; anything you distribute or run as a network service built on it must be
offered under the same licence, with its source available.
