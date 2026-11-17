# BublaGene

![CI](https://github.com/KowalskiBio/bublagene/actions/workflows/ci.yml/badge.svg)
![Licence: AGPL-3.0-or-later](https://img.shields.io/badge/licence-AGPL--3.0--or--later-blue)

**A local-first DNA workbench for the browser and the command line.**

BublaGene reads `.dna`, GenBank and FASTA files and gives you maps, editable
annotations, primers, PCR, restriction cloning, assembly, alignment, variants,
chromatograms, CRISPR guides, and RNA folding. Everything runs on your own
machine — no sequence is uploaded.

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

No build step, no bundler, no server-side component. Node 20 or newer is needed
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

### Maps, sequence and editing

- Circular and linear maps with stacked feature lanes, enzyme sites and cut
  positions.
- Both strands, feature bars, codon-aligned amino acids, a six-frame translation
  track and a live ORF list.
- Drag to select; right-click to copy as bases, reverse complement, protein or
  FASTA, make primers, add a feature, delete/replace bases, or set a new origin.
- Feature editing includes name, type, strand, colour and arbitrary
  qualifiers. Features can be multi-selected, reordered by drag, and copied
  with their annotations.
- Fragments copy and paste across documents; undo/redo coalesces bursts of
  typing; autosave keeps an IndexedDB session plus a small recovery mirror.

### Primers and PCR

- SantaLucia nearest-neighbour Tm with salt/Mg correction, GC, 3′ clamp and
  degenerate-primer expansion.
- Hairpin, self-dimer and pair-dimer ΔG for ranking candidates.
- Ranked primer pairs with product-size, Tm-window, specificity, tail and
  restriction-site options.
- Mismatch-aware PCR with a strict 3′ rule, circular templates, primer-binding
  features on the product, and a simulated gel.

### Enzymes and cloning

- 94 curated enzymes with cut offsets, overhang type/length, isoschizomers and
  methylation sensitivity.
- Custom enzyme sets as plain JSON (`bublagene/enzymes/1`).
- Unique cutters, cut-count filters, per-enzyme fragment tables, simulated gels,
  and a Type IIS overhang compatibility matrix.
- Restriction cloning checks both insert orientations and carries surviving
  features into the construct.
- Gibson/homology assembly joins up to 12 fragments (circular or linear);
  Golden Gate assembly matches Type IIS overhangs. Every assembly has a
  self-contained replayable plan (`bublagene/assembly/1`).
- Batch-assemble many constructs from one fragment library.

### Alignment, variants and contigs

- Affine-gap global, local and banded pairwise alignment with identity and
  gap statistics.
- Dot plots and an alignment panel with mismatch/indel colouring.
- Quality-weighted contig assembly from Sanger reads.
- Normalised SNPs, indels and MNVs with reference-anchored coordinates that
  round-trip through `variation` features.

### Traces

- ABI (`.ab1`), SCF and ZTR chromatograms with four dye channels, called bases,
  peak positions and Phred qualities.
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
- SVG, PNG and PDF diagrams.
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
```

## How it is built

```
js/
  seq.js              sequence maths: complement, translation, Tm, ORFs, IUPAC
  codon-tables.js     25 NCBI genetic codes, validation, per-table translation
  model.js            document model, feature catalog, edits that keep features in place
  thermo.js           nearest-neighbour thermodynamics and secondary-structure ΔG
  primers.js          binding, ranked design, mismatch PCR, gel model
  enzymes.js          enzyme dataset, custom sets, digests, overhangs
  cloning.js          restriction cloning and Gibson assembly
  assembly.js         homology/Golden Gate assembly, plans, batch runs
  align.js            affine-gap global/local/banded alignment and dot plots
  contig.js           quality-weighted contig assembly
  variants.js         normalised variant calling and feature round-trip
  crispr.js           guide finding, off-target heuristic, edit simulation
  fold.js             bounded RNA MFE folding
  traces.js           trace trimming, slicing and reference mapping
  fragment.js         fragment copy/paste and feature reordering
  cache.js            digest-keyed caches and stale-alignment invalidation
  undo.js             coalescing undo/redo
  search.js           sequence/feature/qualifier search and filters
  storage.js          File System Access, OPFS and download fallbacks
  compute.js          worker offload with local fallback
  project.js          project archives and schema migrations
  export.js           SVG, PNG and PDF diagrams
  xml.js              small XML parser/writer
  formats/            dna, genbank/fasta, gff3/bed/csv, ab1, scf, ztr, bam
  ui/                 browser views and panels
  workers/            bounded heavy work
tools/
  bg.mjs              command line
  bench.mjs           5 Mb performance budget
  lint.mjs            syntax lint
tests/                unit, property, fuzz, CLI and jsdom UI tests
```

The browser and CLI share the same engine modules. There is no duplicated
analysis logic.

## Testing

```bash
npm test              # 128 tests
npm run lint          # parse every source file with node --check
npm run typecheck     # tsc --checkJs over engines and CLI
npm run test:coverage # line/function thresholds
npm run bench         # 5 Mb genome budget
npm run check         # lint + typecheck + tests
```

Without external fixtures the suite reports **118 pass, 10 skip, 0 fail**. The
skips are environment-gated: tests that need `.dna` files and tests that compare
the BAM reader with `samtools`.

Point the sample tests at a folder of `.dna` files:

```bash
DNA_SAMPLES=~/plasmids npm test
```

With a folder containing alignments but no traces or `samtools`, the suite
reports **124 pass, 4 skip, 0 fail**. To run the samtools comparison too:

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
- `.dna` history blocks that use an unsupported compression are preserved
  byte-for-byte and left undecoded.

## Licence

AGPL-3.0-or-later. See [LICENSE](LICENSE). You may use, study, change and share
this; anything you distribute or run as a network service built on it must be
offered under the same licence, with its source available.
