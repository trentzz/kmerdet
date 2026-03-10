# kmerdet — Development Guide

## What This Project Is

kmerdet is a Rust reimplementation of the km + kmtools + kam ecosystem for k-mer-based variant detection in liquid biopsy ctDNA monitoring. It must function as a **drop-in replacement** for `km find_mutation` and `kmtools chunk`, while adding features from `kmtools filter`, `stats`, `plot`, and more.

**Primary users**: Bioinformatics researchers running targeted variant detection on UMI-deduplicated sequencing data from liquid biopsy panels.

## Guiding Principles

### 1. Correctness Over Cleverness
Every algorithm must produce results that match km's output on identical inputs. When in doubt, match km's behavior. Deviations from km must be intentional, documented, and tested.

### 2. Types Are Documentation
Use Rust's type system to encode domain knowledge. A `Kmer` is not a `u64` — it's a `Kmer` struct. A `VariantCall` is not a bag of strings — it's a typed struct with `VariantType`, `f64` rVAF, and `Option<String>` fields. Newtypes, enums, and domain-specific types prevent bugs and make the code self-documenting.

### 3. Test Everything
- Unit tests for every public function
- Property-based tests for encoding/decoding roundtrips
- Integration tests comparing output against km on reference data
- Tests for edge cases: empty inputs, single k-mer targets, zero-count k-mers, homopolymer sequences
- The `MockDb` in `tests/common/mod.rs` enables testing without jellyfish installed

### 4. Performance Matters
This tool runs in clinical pipelines where turnaround time affects patient care. The thesis achieved ~6 minutes end-to-end; kmerdet should be faster. Use rayon for parallel target processing. Avoid unnecessary allocations in hot paths (walking, graph construction). Profile before optimizing.

### 5. Researcher-Friendly CLI
Researchers using this tool are accustomed to `km find_mutation` and `kmtools chunk`. The CLI must feel familiar:
- `kmerdet detect` replaces `km find_mutation` + `kmtools chunk` (parallel execution built-in)
- `kmerdet filter` replaces `kmtools filter`
- `kmerdet merge` replaces `kmtools merge`
- `kmerdet stats` replaces `kmtools stats`
- `kmerdet run` orchestrates the full pipeline in one command

Default parameters match km's liquid biopsy defaults: `--count 2 --ratio 0.05`. TSV output matches km's column format for backward compatibility.

### 6. Fail Loudly at Boundaries, Recover Gracefully Inside
- Validate ALL inputs eagerly (file existence, format, parameter ranges) before processing
- Per-target errors are non-fatal by default — log the error, skip the target, continue
- Use `anyhow` at CLI boundaries, `thiserror` for typed library errors
- Never silently produce wrong results — prefer an error over a wrong answer

## Documentation Structure

### `docs/research/initial/` — Phase 0 Research (Foundation)
10 documents covering the original tools (km, kmtools, kam), jellyfish internals, the k-mer variant detection algorithm, liquid biopsy context, and the thesis findings. **Read these first** to understand the domain.

### `docs/research/thesis-context/` — Thesis Analysis
- `thesis-summary.md` — What the thesis found: 77% SNV sensitivity, 38% INDEL sensitivity, ~0.1% VAF threshold, 6-minute pipeline
- `room-for-improvement.md` — Every limitation identified, with severity ratings and fix approaches

### `docs/research/sensitivity-and-confidence/` — Core Research Question
How to detect more true variants (sensitivity) and trust the ones we find (confidence):
- `sensitivity-landscape.md` — Where variants are lost, by type/VAF/pipeline stage
- `confidence-metrics.md` — Statistical frameworks: binomial, Poisson, Bayesian, Fisher's exact
- `false-positive-analysis.md` — FP sources and suppression strategies

### `docs/research/kmer-length-optimization/` — Multi-k Approaches
- `kmer-length-tradeoffs.md` — Theory and practice of k-mer length selection
- `multi-k-strategy.md` — Concrete multi-k detection approaches

### `docs/research/error-correction/` — Upstream Improvements
- `humid-analysis.md` — HUMID strengths, limitations, alternatives
- `error-correction-alternatives.md` — K-mer error correction, UMI-aware counting

### `docs/research/kmer-counting/` — Counting Infrastructure
- `jellyfish-deep-dive.md` — JF internals: counters, canonical mode, hash behavior
- `counting-alternatives.md` — KMC3, DSK, native Rust counting

### `docs/research/graph-and-walking/` — Algorithm Improvements
- `walking-improvements.md` — Adaptive thresholds, bidirectional walking, branch pruning
- `graph-improvements.md` — Edge weighting, graph pruning, alternative pathfinding

### `docs/research/filtering-strategies/` — Smarter Filtering
- `adaptive-filtering.md` — Coverage-based thresholds, panel-of-normals
- `indel-filtering.md` — INDEL normalization, representation equivalence

### `docs/features/` — Feature Specifications
One file per feature with: what, why, research backing, design, acceptance criteria. These are the **implementation specs** — read the relevant feature doc before implementing.

### `docs/planning/` — Roadmaps
- `research-roadmap.md` — 10 experiments to validate improvements
- `development-roadmap.md` — 6-phase implementation plan
- `open-questions.md` — Unresolved questions and investigation approaches

## Architecture

```
src/
├── main.rs              # CLI entry, logging init, subcommand dispatch
├── lib.rs               # Public module declarations, run() entry point
├── config.rs            # TOML config loading (DetectConfig, FilterConfig, etc.)
├── cli/                 # Clap derive subcommands (detect, filter, merge, stats, plot, coverage, run)
├── kmer/                # 2-bit packed Kmer struct, encoding, canonical form
├── jellyfish/           # KmerDatabase trait + conditional FFI to libjellyfish
├── walker/              # Iterative DFS k-mer walking with extension thresholding
├── graph/               # Directed weighted k-mer graph, Dijkstra, all-shortest-paths
├── sequence/            # Target FASTA loading (needletail), RefSeq decomposition, KmerPath
├── variant/             # VariantType, classifier (diff_path_without_overlap), NNLS quantifier
├── filter/              # Tumor-informed filtering (reference mode + alt-sequence mode)
└── output/              # Multi-format writers: TSV, CSV, VCF 4.3, JSON/JSONL, Excel
```

## Key Types

- `Kmer` — 2-bit packed k-mer in u64, with canonical form and extension methods
- `KmerDatabase` trait — Abstraction over jellyfish DB (FFI or mock)
- `WalkResult` — Nodes (k-mer→count map) + reference k-mer set from DFS walking
- `KmerGraph` — Directed weighted graph with BigBang/BigCrunch virtual nodes
- `KmerPath` — Ordered k-mer sequence reconstructible to a DNA string
- `Classification` — Variant type + name + alleles from path comparison
- `Quantification` — NNLS coefficients, rVAFs, min coverages
- `VariantCall` — Complete variant record ready for output
- `FilterResult` — Filtering outcome with match status

## The km Algorithm (Reference Implementation)

The core algorithm that `detect` implements:

1. **Load targets**: Parse FASTA files into `RefSeq` (target + overlapping k-mers)
2. **Walk**: For each reference k-mer, DFS-extend trying all 4 bases, keeping extensions above `max(count, ratio * sibling_total)`
3. **Build graph**: K-mers as nodes, (k-1)-overlap edges, ref edges weighted 0.01, alt edges 1.0
4. **Find paths**: Dijkstra from source/sink, remove ref edges, enumerate all shortest paths through non-ref edges
5. **Classify**: Compare ref path vs alt paths using `diff_path_without_overlap` → SNV, insertion, deletion, ITD, complex
6. **Quantify**: Build contribution matrix, solve via NNLS least squares, compute rVAF = coef/sum(coefs)

## Testing Approach

```bash
cargo test                    # All unit + integration tests
cargo test --test test_kmer   # Just k-mer tests
cargo test -- --ignored       # Run ignored tests (require implementation)
```

The `MockDb` in `tests/common/mod.rs` is the primary testing tool:
```rust
let mut db = MockDb::new(5);
db.set_sequence("ACGTACGT", 100);  // Set all 5-mers from this sequence to count 100
db.set("ACGTA", 5);                // Override specific k-mer count
```

## Build Notes

- `build.rs` conditionally compiles jellyfish FFI — the project builds without jellyfish installed
- `cfg(has_jellyfish)` gates FFI code — set by build.rs when pkg-config finds jellyfish
- All tests pass without jellyfish via `MockDb`
