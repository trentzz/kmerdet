# Feature: `detect` -- Core Variant Detection

## What It Does

The `detect` subcommand is the primary entry point for k-mer variant detection. It takes
a jellyfish k-mer database (.jf file) and one or more target FASTA sequences, then
executes the full detection pipeline: k-mer walking, graph construction, pathfinding,
variant classification, and NNLS quantification. The output is a set of variant calls
with relative VAF estimates, expression coefficients, and coverage metrics.

```
kmerdet detect -d sample.jf -t targets/ -o results.tsv
```

### Input Specification

| Input | Flag | Format | Required |
|-------|------|--------|----------|
| Jellyfish database | `-d` / `--database` | `.jf` binary (jellyfish 2.x) | Yes |
| Target sequences | `-t` / `--targets` | Directory of `.fa` files or a single multi-FASTA | Yes |
| Config file | `--config` | TOML (kmerdet.toml) | No |
| Output path | `-o` / `--output` | File path or `-` for stdout | No (defaults to stdout) |

### Output Specification

Default output is TSV with the following columns, preserving backward compatibility
with km while adding new fields:

| Column | Type | Description | Source |
|--------|------|-------------|--------|
| Database | string | Path to the .jf file queried | Input |
| Query | string | Target FASTA filename | Input |
| Type | string | Variant type (Substitution, Insertion, Deletion, ITD, Complex, Reference) | Classifier |
| Variant_name | string | Human-readable variant descriptor (e.g., chr17:7577120:C/T) | Classifier |
| rVAF | float | Relative variant allele frequency (0.0-1.0) | NNLS quantifier |
| Expression | float | Absolute expression coefficient for the variant path | NNLS quantifier |
| Min_coverage | u64 | Minimum k-mer count along the variant path | Walker |
| Start_offset | u64 | Offset from the target start to the variant | Classifier |
| Sequence | string | Variant path sequence | Walker |
| Reference_expression | float | Expression coefficient for the reference path | NNLS quantifier |
| Reference_sequence | string | Reference path sequence | Walker |
| Info | string | Semicolon-delimited key=value annotation pairs | Pipeline |
| Confidence | float | Phred-scaled quality score (kmerdet extension) | Statistical model |
| Strand_bias | float | Fisher strand bias metric, 0=balanced (kmerdet extension) | Strand-aware counting |

Additional output formats: CSV, VCF 4.3, JSON, JSONL, Excel (.xlsx). Selected via
`--format`.

---

## Why It Matters

The `detect` subcommand replaces three tools from the original pipeline: `km find_mutation`
(the Python k-mer walker), `kmtools chunk` (the parallel executor), and `refolder`
(the file organizer). Consolidating these into a single Rust binary eliminates subprocess
spawning overhead, inter-process serialization, and file-based partitioning. The thesis
pipeline spent approximately 2 minutes on `kmtools chunk` (running km as subprocesses
across 4 threads) for a typical 50-target panel. The performance target for kmerdet
`detect` is under 30 seconds for the same workload -- a 4x improvement from native
execution alone, before algorithmic improvements.

In the clinical workflow, `detect` is called once per patient sample per timepoint. For
MRD monitoring with serial blood draws, a patient may have 10-20 timepoints, each
producing a separate .jf database. Fast detection enables same-day reporting, which the
thesis identified as a critical advantage of k-mer methods over alignment-based pipelines
(6 minutes end-to-end vs. 1-2 hours).

---

## Algorithm Pipeline Stages

### Stage 1: Database and Target Loading

Open the jellyfish database via FFI (or pure Rust reader in Phase 10). Load all target
FASTA files using needletail for zero-copy parsing. Validate that anchor k-mers (first
and last k-mer of each target) are present in the database -- targets with missing anchors
are skipped with a warning.

### Stage 2: K-mer Walking (Iterative DFS)

For each target, perform iterative depth-first search starting from the first k-mer
(source/"BigBang"). At each node, try all four single-base extensions. Keep any extension
whose count in the database exceeds the threshold:

```
threshold = max(sum_of_sibling_counts * ratio, count)
```

Walking is bounded by three parameters:
- `max_stack` (default: 500) -- maximum DFS depth
- `max_break` (default: 10) -- maximum branching points explored
- `max_node` (default: 10000) -- maximum total nodes visited

When any limit is reached, the walk terminates for that target. The result is a set of
discovered k-mers and their extension relationships.

### Stage 3: Graph Construction

Build a directed weighted graph from the walking results:
- Nodes: discovered k-mers
- Edges: k-1 suffix-prefix overlaps
- Weights: reference edges = 0.01, non-reference edges = 1.0

The source node is the first k-mer of the target; the sink is the last k-mer.

### Stage 4: Pathfinding

Run Dijkstra's algorithm from source to find shortest paths to all reachable nodes.
Run Dijkstra from sink (reversed graph) to find shortest paths back. The reference
path is the shortest source-to-sink path (dominated by 0.01-weight reference edges).

Remove all reference edges. For each remaining non-reference edge, reconstruct the
shortest path through that edge using the forward and backward Dijkstra results.
Deduplicate the resulting alternative paths.

### Stage 5: Variant Classification

Compare each alternative path to the reference path using the `diff_path_without_overlap`
algorithm. Walk from both ends to find the divergence region, then classify:

| Condition | Classification |
|-----------|---------------|
| Paths identical | Reference |
| Same length, bases differ | Substitution (SNV/MNP) |
| Alt path longer, no deleted bases | Insertion |
| Ref path longer, no inserted bases | Deletion |
| Variant retraces reference segment | Internal Tandem Duplication |
| Mixed insertion and deletion | Complex Indel |

### Stage 6: NNLS Quantification

Build a contribution matrix where entry (i, j) = number of times k-mer i appears in
path j. Solve the non-negative least squares problem:

```
minimize ||contrib @ coef - counts||^2  subject to coef >= 0
```

Compute rVAF as `coef[j] / sum(coef)` for each path j. The Lawson-Hanson NNLS algorithm
is used directly, replacing km's approach of unconstrained least squares followed by
gradient descent correction for negative coefficients.

---

## Key Parameters

| Parameter | Flag | Default | Range | Effect |
|-----------|------|---------|-------|--------|
| k-mer length | `-k` / `--kmer-length` | 31 | 15-63 (odd) | Specificity vs. INDEL sensitivity |
| Extension ratio | `--ratio` | 0.05 | 0.0-1.0 | Lower = more sensitive, more noise |
| Minimum count | `--count` | 2 | 1-1000 | Absolute floor for k-mer extension |
| Max stack depth | `--max-stack` | 500 | 100-10000 | Walking depth limit |
| Max branch points | `--max-break` | 10 | 1-100 | Walking breadth limit |
| Max nodes | `--max-node` | 10000 | 1000-1000000 | Total walking budget |
| Canonical mode | `--canonical` | true | bool | Count both strand orientations together |

The thesis validation used `k=31, ratio=0.00001, count=2` for maximum sensitivity in
liquid biopsy. The recommended balanced configuration is `ratio=0.0001, count=3` per
thesis Chapter 5. These defaults can be overridden per-run via CLI flags or per-project
via the TOML config file, following the three-tier precedence: CLI > config file > defaults.

---

## Parallelism Model

Detection uses rayon's work-stealing thread pool for target-level parallelism. Each
target is an independent unit of work: it reads from the shared (read-only) jellyfish
database, performs walking and classification, and produces a result struct. No
synchronization is needed between targets during processing.

```
[Target 1] ──> [Walk] ──> [Graph] ──> [Classify] ──> [Quantify] ──> Result 1
[Target 2] ──> [Walk] ──> [Graph] ──> [Classify] ──> [Quantify] ──> Result 2
   ...                                                                  ...
[Target N] ──> [Walk] ──> [Graph] ──> [Classify] ──> [Quantify] ──> Result N
              └──────────────── rayon work-stealing ────────────────┘
```

Thread count defaults to the number of logical CPUs but can be set via `--threads`.
Results are collected into a Vec and sorted deterministically (by target filename) before
output, ensuring reproducible output regardless of thread scheduling.

---

## Error Recovery

Per-target errors are non-fatal by default. If walking exceeds `max_node`, if the graph
has no source-to-sink path, or if quantification produces degenerate results, the target
is logged as failed and processing continues with remaining targets.

| Error condition | Default behavior | `--strict` behavior |
|-----------------|-----------------|---------------------|
| Walking timeout (max_node) | Warn, skip target | Fatal error |
| No paths found | Warn, skip target | Fatal error |
| Quantification failure | Warn, report as failed | Fatal error |
| Anchor k-mer missing | Warn, skip target | Fatal error |
| Thread panic | Catch at rayon boundary, continue | Fatal error |
| Database read error | Fatal error | Fatal error |
| Output write failure | Fatal error | Fatal error |

The exit code is 0 if all targets succeed, 2 if some targets fail (partial success),
and 1 for fatal errors. A summary of failed targets is printed to stderr at completion.

---

## Research Backing

The core algorithm is from Audoux et al. (2017), as implemented in km and validated
in the thesis on a 10-patient cohort with UMI-based duplex sequencing. The thesis
established 77% SNV sensitivity and 38% INDEL sensitivity with k=31 and permissive
parameters.

Key improvements over the original km implementation incorporated into `detect`:

1. **Iterative DFS instead of recursive**: Eliminates stack overflow risk for deep walks
   and enables precise control of the stack depth limit. Research in walking-improvements.md
   details the conversion from Python recursion to Rust iterative DFS with an explicit
   stack.

2. **Lawson-Hanson NNLS**: Replaces the two-phase approach (unconstrained least squares
   then gradient descent correction) with a proper NNLS solver. This produces provably
   non-negative coefficients without iterative correction, improving quantification
   accuracy at low VAF. See room-for-improvement.md limitation #9.

3. **Confidence scoring**: Each variant call receives a Phred-scaled quality score based
   on a binomial test of the variant k-mer count against the local error rate. This
   addresses the thesis finding that k-mer calls lacked statistical quality annotation
   (confidence-metrics.md Section 1.2).

4. **Strand bias reporting**: When strand-aware k-mer counting is available, report
   the Fisher strand bias metric. Guo et al. (2012) showed that strand bias filtering
   removes 70-90% of false positive variant calls in Illumina data.

5. **Native parallelism**: Rayon work-stealing replaces the file-based partitioning of
   kmtools chunk + refolder, eliminating subprocess overhead and enabling dynamic load
   balancing.

---

## Performance Target

| Metric | Thesis pipeline (km + kmtools) | kmerdet detect target |
|--------|-------------------------------|----------------------|
| 50 targets, 4 threads | ~2 minutes | <30 seconds |
| Per-target latency | ~2.4 seconds | <0.6 seconds |
| Memory (peak RSS) | ~400 MB (Python + jellyfish) | <300 MB |
| Startup overhead | ~1s (Python interpreter) | <50ms |

The 4x speedup comes from: eliminating subprocess spawning (~0.5s per km invocation x 50
targets / 4 threads = ~6s), eliminating Python interpreter overhead, native in-process
parallelism with rayon work-stealing, and Rust's zero-cost abstractions for the walking
inner loop (2-bit packed k-mers, branch-free extension).

---

## Acceptance Criteria

### Unit Tests

- [ ] K-mer encoding: round-trip encode/decode for all 4^k k-mers at k=3,5 (exhaustive)
- [ ] K-mer extension: verify all 4 extensions produce correct k-mers
- [ ] Walking: known target with known variants produces expected paths
- [ ] Graph construction: correct nodes, edges, and weights from walking output
- [ ] Pathfinding: Dijkstra finds shortest paths; alternative paths correctly extracted
- [ ] Classification: SNV, insertion, deletion, ITD, complex each correctly classified
- [ ] Quantification: NNLS produces correct rVAF for synthetic contribution matrices
- [ ] Confidence score: binomial p-value matches hand-calculated value

### Integration Tests

- [ ] Full pipeline on reference test data matches km output (variant type, rVAF within
      tolerance of +/- 0.01)
- [ ] Multi-target parallel execution produces deterministic output
- [ ] `--strict` mode fails on first target error
- [ ] All output formats (TSV, CSV, VCF, JSON, JSONL) are well-formed
- [ ] Empty target directory produces empty output (not an error)
- [ ] Target with non-unique anchor k-mers is skipped with warning

### Performance Tests

- [ ] 50-target panel completes in <30 seconds on 4 cores
- [ ] Memory usage stays below 300 MB for typical panel
- [ ] Thread scaling: 2x threads yields >1.5x speedup (up to core count)
