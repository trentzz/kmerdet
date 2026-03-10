# Development Roadmap: Research-Informed Phased Implementation

This document replaces the original Phase 1-10 roadmap with a research-informed plan. Each phase has clear milestones, lists files to modify or create, and references the research findings that motivate the design decisions.

---

## Phase Overview

| Phase | Name | Scope | Milestone |
|-------|------|-------|-----------|
| 1 | Core Detection | K-mer primitives through variant classification | Match km output on test data |
| 2 | Filtering + Output | Filter, merge, output formats, stats | Full detect-filter-stats pipeline |
| 3 | Confidence + Adaptive | Statistical scoring, adaptive thresholds, walking improvements | Measurable sensitivity improvement |
| 4 | Multi-k + Advanced | Multi-k detection, PoN, plot, run subcommands | Multi-k improves INDEL sensitivity >20% |
| 5 | Pure Rust + Performance | Eliminate C++ FFI, native counting, SIMD | No external dependencies |
| 6 | UMI Integration + Stretch | UMI-aware counting, ML filtering, SV extensions | Detection without HUMID preprocessing |

---

## Phase 1: Core Detection (Reproduce Thesis Results)

### Objective
Build the complete detection pipeline from k-mer primitives through variant classification and quantification. The end-to-end `detect` subcommand must produce output matching km on the validation data.

### Motivation
The thesis achieved 77% SNV sensitivity and 38% INDEL sensitivity with this pipeline (`thesis-summary.md`). Reproducing these numbers validates the Rust reimplementation before any improvements are attempted. This corresponds to Research Experiment 1 (Baseline Reproduction).

### Components

#### 1.1 K-mer Primitives
- **Files**: `src/kmer/encoding.rs`, `src/kmer/mod.rs`
- **Work**: Complete 2-bit packed encoding (A=00, C=01, G=10, T=11) in u64. Implement canonical form (lexicographic minimum of k-mer and reverse complement). Support k=21 through k=43 (odd values only, per `kmer-length-tradeoffs.md` convention for avoiding palindromes).
- **Reference**: `counting-alternatives.md` Section "2-Bit Encoding and Canonical Form" provides the encoding and canonical k-mer iterator implementation.
- **Validation**: Unit tests comparing canonical forms against jellyfish query output.

#### 1.2 Jellyfish FFI Reader
- **Files**: `src/jellyfish/ffi.rs`, `src/jellyfish/database.rs`, `jf_wrapper/`
- **Work**: Complete the C++ FFI wrapper for reading .jf files. Implement the `KmerDatabase` trait with `query(kmer) -> count` and `check(kmer) -> bool` methods. Handle canonical mode matching jellyfish's `-C` flag. Support memory-mapped access for .jf files per `jellyfish-deep-dive.md` recommendation to use `mmap` with `populate()` for targeted panel data.
- **Reference**: `jellyfish-deep-dive.md` for hash table internals, reversible hash function, variable-length counters, and memory-mapping strategy.
- **Validation**: Query known k-mers from test .jf databases and verify counts match `jellyfish query` output.

#### 1.3 Target Sequence Loading
- **Files**: `src/sequence/target.rs`, `src/sequence/refseq.rs`
- **Work**: Parse target FASTA files using needletail. Extract reference k-mer sequence for each target. Implement `RefSeq` decomposition that identifies the reference k-mer path through each target. Handle the `--flank 35` convention from `thesis-summary.md`.
- **Validation**: Decomposed reference k-mers match jellyfish counts at expected coverage levels.

#### 1.4 K-mer Walking (Iterative DFS)
- **Files**: `src/walker/mod.rs`, `src/walker/extension.rs`, `src/walker/state.rs`
- **Work**: Implement iterative DFS walking with explicit stack. The extension threshold mirrors km's logic: `threshold = max(sum_of_sibling_counts * ratio, n_cutoff)` per `walking-improvements.md` Section 1. Support configurable `max_stack` (default 500), `max_break` (default 10), and `max_node` (default 10000). Forward extension queries all four single-base extensions, retains those meeting threshold.
- **Design decisions informed by research**:
  - Use iterative (not recursive) DFS for stack safety and pruning control (`walking-improvements.md` Section 7)
  - Implement `extend_backward` stub for Phase 3 bidirectional walking (`walking-improvements.md` Section 3)
  - Store walk results as `HashMap<String, u64>` (k-mer string -> count) initially; optimize to u64 keys in Phase 5 (`walking-improvements.md` Section 7 on memory efficiency)
- **Reference**: `walking-improvements.md` for complete threshold analysis, `sensitivity-landscape.md` Section 3.3 for walking as the primary false negative source.
- **Validation**: Walk results on test targets match km `find_mutation` output.

#### 1.5 Graph Construction + Pathfinding
- **Files**: `src/graph/builder.rs`, `src/graph/dijkstra.rs`, `src/graph/paths.rs`, `src/graph/types.rs`
- **Work**: Build directed weighted graph from walk results. Add "BigBang" source and "BigCrunch" sink nodes. Use binary edge weights (ref=0.01, alt=1.0) matching km per `graph-improvements.md` Section 1. Implement Dijkstra's algorithm from both source and sink. Enumerate all shortest paths through non-reference edges.
- **Design decisions informed by research**:
  - Retain binary weights for Phase 1 compatibility; plan coverage-weighted alternative for Phase 3 (`graph-improvements.md` Section 2)
  - Implement dead-end removal as a baseline pruning step (`graph-improvements.md` Section 4)
  - Add `--debug-graph` DOT format export for debugging from the start (`graph-improvements.md` Section 6)
- **Reference**: `graph-improvements.md` for current algorithm analysis and planned improvements.
- **Validation**: Paths discovered match km's path output for all test targets.

#### 1.6 Variant Classification
- **Files**: `src/variant/classifier.rs`, `src/variant/types.rs`
- **Work**: Implement `diff_path_without_overlap` logic from km. Compare reference and alternative paths to determine variant type (Reference, SNV, Insertion, Deletion, ITD, Complex). Walk from both ends to find divergence region. Handle edge cases documented in `sensitivity-landscape.md` Section 3.5 (variants at target boundaries, multiple variants on same path, compound indel+SNV).
- **Reference**: `sensitivity-landscape.md` Section 3.5 for classification failure modes, `room-for-improvement.md` Section 2 for INDEL classification edge cases in `diff_path_without_overlap`.
- **Validation**: Classification matches km output on all test targets.

#### 1.7 NNLS Quantification
- **Files**: `src/variant/quantifier.rs`, `src/variant/cluster.rs`
- **Work**: Implement NNLS (non-negative least squares) decomposition. Build contribution matrix mapping k-mers to paths. Solve `||contrib @ coef - counts||^2` subject to `coef >= 0`. Compute rVAF as `coef[j] / sum(coef)`. Handle cluster mode for overlapping variants.
- **Design decisions informed by research**:
  - Use a proper NNLS solver (Lawson-Hanson) rather than least squares + gradient descent correction per `room-for-improvement.md` Section 9
  - Plan for bootstrap CI capability in Phase 3 per `confidence-metrics.md` Section 3.4
- **Reference**: `thesis-summary.md` Section "Quantification Approach", `room-for-improvement.md` Section 9 for NNLS limitations and alternatives.
- **Validation**: rVAF values match km within 1% absolute tolerance.

#### 1.8 Detect Subcommand + Coverage Subcommand
- **Files**: `src/cli/detect.rs`, `src/cli/coverage.rs`, `src/pipeline/detect.rs`
- **Work**: Wire together all components into the `detect` subcommand. Process targets in parallel via rayon. Output per-target results with variant name, type, rVAF, expression, min_coverage. Implement `coverage` subcommand to query reference k-mer coverage statistics per target.
- **Validation**: End-to-end test on validation data matches km output.

### Dependencies
- Phase 0 scaffolding (complete)
- Jellyfish installed for .jf file generation (build.rs conditional compilation)
- Test data (10-patient .jf files, target FASTA, reference mutations)

### Estimated Scope
- ~15-20 source files with substantive implementation
- ~3000-5000 lines of Rust code
- 4-6 weeks for a single developer

### Milestone
**Match km output on the 10-patient test data within 1% rVAF tolerance. SNV sensitivity matches 77% (+/- 2pp). INDEL sensitivity matches 38% (+/- 3pp).**

---

## Phase 2: Filtering + Output

### Objective
Complete the filtering pipeline and all output formats. The full detect-filter-stats workflow should be functional, producing clinically interpretable output.

### Motivation
The thesis pipeline requires filtering to achieve zero false positives (`false-positive-analysis.md` Section 5). INDEL normalization during filtering is critical for matching INDELs correctly (`indel-filtering.md`). Multiple output formats serve different downstream needs (VCF for clinical, TSV for km compatibility, JSON for programmatic access).

### Components

#### 2.1 Filter Subcommand -- Reference Mode
- **Files**: `src/filter/reference_mode.rs`, `src/filter/mod.rs`, `src/cli/filter.rs`
- **Work**: Implement coordinate-based matching of detected variants against a reference mutation file. Apply min_coverage, min_vaf, and min_expression thresholds. Support variant type filtering.
- **Critical addition**: Implement INDEL left-alignment normalization as the default matching behavior per `indel-filtering.md` Section 2. The Tan et al. (2015) algorithm: trim suffix, trim prefix, left-align using reference sequence, re-trim. This addresses the representation mismatch problem identified in `indel-filtering.md` Section 1 where different tools produce different INDEL coordinates for the same biological event.
- **Fallback**: ALT sequence comparison when normalized coordinates do not match per `indel-filtering.md` Section 3. The recommended strategy is normalize-then-alt-fallback with the match source reported in output.
- **Reference**: `indel-filtering.md` for the complete normalization algorithm with Rust implementation, test cases, and validation against bcftools norm.
- **Validation**: Filter results match kmtools filter output. INDEL normalization matches bcftools norm on test cases.

#### 2.2 Filter Subcommand -- ALT-Sequence Mode
- **Files**: `src/filter/alt_mode.rs`
- **Work**: Implement alternative sequence matching using full path sequences rather than coordinates. Substring and overlap matching per `indel-filtering.md` Section 3. This mode handles complex INDELs and delins where coordinate matching fails even after normalization.
- **Reference**: `indel-filtering.md` Section 3 for ALT sequence matching algorithm.

#### 2.3 INDEL Deduplication
- **Files**: `src/variant/normalize.rs`
- **Work**: Implement INDEL deduplication by normalizing all detected INDELs before output per `indel-filtering.md` Section 5. Multiple k-mer graph paths representing the same deletion in a repeat (e.g., deletion of one T from TTTTT reported at three different positions) are collapsed to a single canonical representation.
- **Reference**: `indel-filtering.md` Section 5.

#### 2.4 Merge Subcommand
- **Files**: `src/cli/merge.rs`, `src/pipeline/merge.rs`
- **Work**: Combine results from multiple detection runs (e.g., different k values, or chunked execution). Deduplicate variants, merge rVAF estimates, preserve per-run metadata. This supports the multi-k workflow planned for Phase 4.

#### 2.5 Output Formats
- **Files**: `src/output/tsv.rs`, `src/output/vcf.rs`, `src/output/json.rs`, `src/output/csv.rs`, `src/output/excel.rs`
- **Work**:
  - **TSV** (km-compatible): Tab-separated format matching km output columns for backward compatibility.
  - **VCF 4.3**: Standard variant call format with INFO fields for rVAF, expression, min_coverage. FILTER field with PASS/LowQual/etc. per `confidence-metrics.md` Section 3.5. Normalized INDEL coordinates.
  - **JSON/JSONL**: Machine-readable format with full variant details, graph metadata, and filter annotations.
  - **CSV**: Spreadsheet-compatible format.
  - **Excel**: Optional xlsx output via the `rust_xlsxwriter` crate.
- **Reference**: `confidence-metrics.md` Section 4.1 for VCF annotation design, including k-mer equivalents for standard VCF fields (DP, AD, AF, FS).

#### 2.6 Stats Subcommand
- **Files**: `src/cli/stats.rs`, `src/pipeline/stats.rs`
- **Work**: Compute summary statistics from detection and filter results. Per-sample: coverage statistics (median, CV, per-target), variant counts by type, sensitivity estimates (if ground truth provided). Per-target: coverage, number of variants, filter pass rate. This supports Research Experiment 2 (Bottleneck Analysis).

### Dependencies
- Phase 1 complete (working detect subcommand)

### Estimated Scope
- ~10-12 source files
- ~2000-3000 lines of Rust code
- 3-4 weeks

### Milestone
**Full detect -> filter -> stats pipeline works end-to-end. VCF output passes bcftools validation. INDEL normalization matches bcftools norm. Filter reference mode with normalization achieves higher INDEL match rate than without.**

---

## Phase 3: Confidence + Adaptive (Research-Driven Improvements)

### Objective
Implement the research-driven improvements that are expected to provide the largest sensitivity gains: confidence scoring, adaptive thresholds, walking improvements, and graph improvements. This phase corresponds to Research Experiments 2-9.

### Motivation
The thesis identified fixed thresholds (`room-for-improvement.md` Section 6), lack of quality metrics (`confidence-metrics.md`), and walking limitations (`walking-improvements.md`) as the primary barriers to higher sensitivity. Phase 3 addresses all three simultaneously.

### Components

#### 3.1 Statistical Confidence Scoring
- **Files**: `src/confidence/pvalue.rs`, `src/confidence/qual.rs`, `src/confidence/mod.rs`
- **Work**:
  - **Binomial p-value**: Per variant k-mer, compute one-sided binomial p-value given local coverage and estimated error rate per `confidence-metrics.md` Section 2.1. Formula: `p = 1 - BinomCDF(c-1, D, e_k)`.
  - **Combined p-value**: Fisher's method across k-mers along the variant path with autocorrelation correction (effective m_eff = m * k / read_length) per `confidence-metrics.md` Section 6.1.
  - **Phred-scaled QUAL**: `QUAL = -10 * log10(p_value)` per `confidence-metrics.md` Section 3.5. Add to VCF output.
  - **Negative binomial model**: Estimate overdispersion from reference k-mer count variance per `confidence-metrics.md` Section 2.4. Provide NB p-values as an alternative for high-coverage regions where k-mer counts show overdispersion (typical alpha=0.01-0.1).
- **Error rate estimation**: Implement the empirical approach from `walking-improvements.md` Section 2: query all four extensions of reference k-mers, use non-reference counts to estimate per-sample error rate.
- **Reference**: `confidence-metrics.md` Sections 2-3 and 6 for mathematical details.

#### 3.2 Positional Uniformity + Strand Bias Metrics
- **Files**: `src/confidence/strand_bias.rs`, `src/confidence/positional.rs`
- **Work**:
  - **Positional uniformity**: CV of k-mer counts along variant path per `confidence-metrics.md` Section 3.2. Expected CV ~0.3-0.5 for true variants; >1.0 suggests positional bias.
  - **Strand bias** (partial, full in Phase 5): If strand information is available (from non-canonical counting or separate strand tracking), compute Fisher strand bias per `confidence-metrics.md` Section 3.1. Otherwise, defer to Phase 5 when native Rust counting provides strand-aware data.
- **Reference**: `confidence-metrics.md` Section 3.1-3.2, `jellyfish-deep-dive.md` Section "Loss of Strand-Specific Information" on workarounds for canonical counting.

#### 3.3 Bootstrap Confidence Intervals on rVAF
- **Files**: `src/variant/bootstrap.rs`
- **Work**: Implement 1000-replicate bootstrap for NNLS decomposition per `confidence-metrics.md` Section 3.4. For each replicate: resample k-mer counts from a Poisson distribution with observed counts as parameters, solve NNLS, compute rVAF. Report 2.5th and 97.5th percentiles as 95% CI. NNLS for 31-variable system is microseconds per solve, so 1000 replicates per variant is negligible overhead.
- **Reference**: `confidence-metrics.md` Section 3.4 for bootstrap methodology, `room-for-improvement.md` Section 9 for quantification model limitations motivating CIs.

#### 3.4 Adaptive Walking Thresholds
- **Files**: `src/walker/adaptive.rs`, modify `src/walker/extension.rs`
- **Work**:
  - **Per-sample adaptive threshold**: Compute median coverage and error rate from reference k-mers before walking. Set threshold at Poisson 99th percentile of error distribution per `walking-improvements.md` Section 2. Config: `AdaptiveThresholdConfig { false_extension_rate: 0.01, error_rate: 0.001, coverage_window: 15, min_threshold: 2 }`.
  - **Per-target adaptive threshold**: Adjust based on local GC content, linguistic complexity, and homopolymer length per `adaptive-filtering.md` Section 4. Precompute target annotations (GC content, complexity, max homopolymer) and store in `TargetAnnotation` struct.
  - **Depth-tier auto-detection**: Automatically classify sample into coverage tier (ultra-low <500x, standard 500-2000x, high-depth 2000-10000x, ultra-deep >10000x) and apply tier-appropriate parameters per `adaptive-filtering.md` Section 3.
  - **Backward compatibility**: Retain fixed ratio/count parameters as `--legacy-thresholds` option.
- **Reference**: `adaptive-filtering.md` Sections 2-4, `walking-improvements.md` Section 2, `room-for-improvement.md` Section 6.

#### 3.5 Walking Improvements
- **Files**: `src/walker/bidirectional.rs`, `src/walker/pruning.rs`, modify `src/walker/mod.rs`
- **Work**:
  - **Bidirectional walking**: Walk forward from each reference k-mer AND backward. Merge results (set union of discovered k-mers). This addresses deletions at target boundaries and variants near anchors per `walking-improvements.md` Section 3.
  - **Branch pruning**: Implement tip removal (dead-end branches < k/2 extensions) and bubble smoothing (parallel paths < k with >100:1 coverage ratio) per `walking-improvements.md` Section 5 and `error-correction-alternatives.md` Section "Lightweight Error Correction During Walking". Safety: never remove branches with >= 2 UMI family support or reads from both strands.
  - **Relaxed thresholds during non-reference traversal**: When >= 5 consecutive non-reference extensions (suggesting INDEL traversal), reduce threshold by 50% per `walking-improvements.md` Section 4. This addresses the long INDEL problem where each junction k-mer must independently pass threshold.
  - **Dynamic max_stack**: For targets with expected large INDELs, increase max_stack proportionally per `walking-improvements.md` Section 4.

#### 3.6 Coverage-Weighted Graph Edges (Hybrid)
- **Files**: modify `src/graph/builder.rs`, `src/graph/paths.rs`
- **Work**: Implement the recommended hybrid approach from `graph-improvements.md` Section 2: retain binary weights (ref=0.01, alt=1.0) for structural pathfinding (steps 1-3), then rank discovered alternative paths by coverage-based `path_score = min(variant_kmer_counts)`. This provides meaningful path ordering without changing the pathfinding algorithm.
- **Reference**: `graph-improvements.md` Section 2 for hybrid approach rationale.

#### 3.7 Graph Pruning
- **Files**: `src/graph/pruning.rs`, modify `src/graph/builder.rs`
- **Work**: Implement pruning pipeline from `graph-improvements.md` Section 4, applied between graph construction and pathfinding:
  1. Dead-end removal (BFS forward from source + backward from sink, remove unreachable nodes)
  2. Low-coverage path removal (nodes with count below `max(2, median_coverage * error_rate)`)
  3. Dead-end removal again (clean up after step 2)
  4. Tip clipping (branches < k/2 extensions)
  5. Bubble collapsing (parallel paths < k with < 10% coverage of main path)
- **Reference**: `graph-improvements.md` Section 4.

#### 3.8 Composite Filtering Pipeline
- **Files**: `src/filter/adaptive.rs`, `src/filter/context.rs`, modify `src/filter/mod.rs`
- **Work**: Implement the multi-stage filter pipeline from `adaptive-filtering.md` Section 7:
  1. Adaptive hard thresholds (coverage-proportional count and VAF)
  2. Type filtering
  3. Context filtering (homopolymer adjustment, GC adjustment)
  4. Sequence context annotation (flag GGC, homopolymer, low-complexity)
  Stage 4 (PoN) and Stage 5 (ML) deferred to Phase 4 and 6 respectively.
- Add filter annotations to output per `adaptive-filtering.md` Section 7: FILTER_STAGES, ADAPTIVE_THRESHOLD, CONTEXT_FLAGS columns.

### Dependencies
- Phase 2 complete (filtering + output infrastructure)
- Research Experiments 1-2 complete (baseline established, bottleneck identified)

### Estimated Scope
- ~12-15 source files
- ~4000-6000 lines of Rust code
- 6-8 weeks

### Milestone
**Measurable sensitivity improvement over km baseline. Target: SNV sensitivity >= 82% (from 77%), INDEL sensitivity >= 45% (from 38%). QUAL scores with AUROC >= 0.90 for TP/FP discrimination. Bootstrap 95% CI on rVAF with 90-98% coverage.**

---

## Phase 4: Multi-k + Advanced Features

### Objective
Implement multi-k detection, panel-of-normals, and remaining subcommands. Multi-k is the single most impactful improvement for INDEL sensitivity per `multi-k-strategy.md`.

### Motivation
The thesis identified INDEL sensitivity as the single largest clinical limitation (`room-for-improvement.md` Section 2). Shorter k values (k=21, k=25) dramatically improve junction coverage for insertions >15 bp: at k=21, a 15 bp insertion has 26 spanning k-mers (84% of full signal) vs. 16 at k=31 (52%) per `kmer-length-tradeoffs.md`. Multi-k detection exploits this while maintaining k=31 specificity for SNVs.

### Components

#### 4.1 Per-Target K Selection
- **Files**: `src/config.rs` (extend target catalog), `src/pipeline/detect.rs` (routing)
- **Work**: Implement Approach 3 from `multi-k-strategy.md`: assign each target an optimal k based on variant type, repeat content, and anchor uniqueness. Add `recommended_k` field to target metadata. Routing logic directs each target to the appropriate .jf database during walking.
- **Selection rules** per `multi-k-strategy.md` Section "Per-Target k Assignment Rules":
  - SNV in non-repetitive: k=31
  - Short INDEL (1-10 bp): k=31
  - Medium INDEL (11-20 bp): k=25
  - Long INDEL (>20 bp): k=21
  - ITD: k=21 or k=25
  - Target in segmental duplication: k=43
- **Reference**: `multi-k-strategy.md` Approach 3.

#### 4.2 Multi-k Detection with Consensus
- **Files**: `src/pipeline/multik.rs`, `src/variant/consensus.rs`
- **Work**: Implement Approach 1+2 from `multi-k-strategy.md`: run detection independently at k=21, 31, 41 (or configured set), then merge results with consensus voting.
  - **Variant matching across k values**: Normalize coordinates (especially for INDELs which may have different breakpoints at different k) per `indel-filtering.md`.
  - **Weighted consensus voting**: Weight by variant type per `multi-k-strategy.md` Section "Weighted Voting":
    - SNVs: w(k=21)=0.7, w(k=31)=1.0, w(k=41)=0.9
    - Short INDELs: w(k=21)=1.0, w(k=31)=0.9, w(k=41)=0.7
    - Long INDELs: w(k=21)=1.0, w(k=31)=0.5, w(k=41)=0.2
  - **Confidence tiers**: Tier 1 (all k), Tier 2 (2/3 k), Tier 3 (k=31 only), Tier 4 (single non-standard k).
  - **rVAF reconciliation**: Coverage-weighted average across k values per `multi-k-strategy.md` Section "Handling Conflicting rVAF Estimates".
- **Reference**: `multi-k-strategy.md` Approaches 1-2.

#### 4.3 Panel-of-Normals
- **Files**: `src/filter/pon.rs`, `src/cli/pon.rs`
- **Work**: Implement PoN construction and filtering per `adaptive-filtering.md` Section 5 and `false-positive-analysis.md` Section 5.3.
  - **Build command**: `kmerdet pon build --normals dir/ --targets panel/ --output pon.db`
  - **Storage**: Hash map `(target, normalized_variant) -> PonEntry { frequency, max_vaf, mean_vaf }`. Serialize with bincode. Storage is trivial (~50 bytes per entry per `false-positive-analysis.md`).
  - **Filter integration**: Stage 4 in the composite filter pipeline. Default frequency threshold 5%.
  - **What it catches**: Systematic sequencing errors (GGC motifs), reference genome errors, common CHIP variants, PCR stutter at microsatellites per `false-positive-analysis.md` Section 5.3.
- **Reference**: `adaptive-filtering.md` Section 5, `false-positive-analysis.md` Sections 5.3 and 6.

#### 4.4 Plot Subcommand
- **Files**: `src/cli/plot.rs`, `src/pipeline/plot.rs`
- **Work**: Generate VAF plots for longitudinal MRD monitoring. Output SVG or PNG using a Rust plotting library. Show rVAF trend over timepoints for each variant, with confidence intervals from Phase 3 bootstrap implementation.

#### 4.5 Run Subcommand (Pipeline Orchestration)
- **Files**: `src/cli/run.rs`, `src/pipeline/orchestrator.rs`
- **Work**: Full pipeline orchestration from HUMID-deduplicated FASTQ through final filtered output. Sequence: jellyfish count (at required k values) -> detect (per-target k routing) -> filter (all stages) -> stats -> optional plot. Support `--threads` for Nextflow compatibility per `10-future-improvements.md`.

### Dependencies
- Phase 3 complete (adaptive thresholds, confidence scoring)
- Research Experiment 3 (K-mer Length Sweep) results to validate multi-k parameters

### Estimated Scope
- ~10-12 source files
- ~3000-4000 lines of Rust code
- 5-7 weeks

### Milestone
**Multi-k detection improves INDEL sensitivity by >20% (from ~45% after Phase 3 to >55%). Panel-of-normals reduces false positives when using more permissive detection thresholds. `run` subcommand provides single-command pipeline execution.**

---

## Phase 5: Pure Rust + Performance

### Objective
Eliminate the C++ FFI dependency on jellyfish. Implement native Rust k-mer counting for targeted panels. Add SIMD acceleration and comprehensive benchmarks.

### Motivation
The C++ FFI creates build complexity (requires jellyfish dev headers and pkg-config), limits portability, and prevents UMI-aware extensions. A pure Rust implementation enables single-binary distribution with no external dependencies per `10-future-improvements.md` Section "Pure Rust Jellyfish Reader" and `counting-alternatives.md` Section "Counting K-mers Directly in Rust from FASTQ".

### Components

#### 5.1 Pure Rust .jf File Reader
- **Files**: `src/jellyfish/reader.rs`, `src/jellyfish/hash.rs`
- **Work**: Parse the jellyfish binary format in pure Rust:
  1. Parse JSON header (k, table_size, counter_len, canonical flag, max_reprobe) per `jellyfish-deep-dive.md` Section "Jellyfish Hash Function"
  2. Replicate the reversible hash function from C++ source (`large_hash_array.hpp`)
  3. Implement probing sequence (double hashing with reprobe offset)
  4. Bit-level extraction of key remainder and variable-length counters
  5. Reconstruct full k-mer from table position + key remainder via inverse hash
- **Challenge**: The hash function is not formally documented; it must be reverse-engineered from the C++ source per `jellyfish-deep-dive.md`. This is the primary implementation challenge.
- **Memory mapping**: Use `memmap2` with `MmapOptions::new().populate()` for targeted panel data per `jellyfish-deep-dive.md` recommendation.
- **Reference**: `jellyfish-deep-dive.md` Sections on hash function, counter mechanism, and memory-mapped access. `10-future-improvements.md` Section "Pure Rust Jellyfish Reader".
- **Validation**: Query results must match jellyfish CLI output bit-for-bit on all test .jf files.

#### 5.2 Native Rust K-mer Counting
- **Files**: `src/counting/mod.rs`, `src/counting/hashtable.rs`, `src/counting/stranded.rs`
- **Work**: Implement in-memory k-mer counting from FASTQ per `counting-alternatives.md` Section "Counting K-mers Directly in Rust from FASTQ":
  - **Phase A**: Use DashMap for simplicity (~3.8 GB for 50M k-mers). Parse FASTQ with needletail, extract canonical k-mers, count concurrently with rayon.
  - **Phase B**: Optimize with custom lock-free open-addressing hash table (~750 MB for 50M k-mers) using AtomicU64 keys and AtomicU32 counts per `counting-alternatives.md` Section "Custom hash table".
  - **Strand-aware counting**: Maintain separate forward/reverse counts per canonical k-mer per `jellyfish-deep-dive.md` Section "Workarounds for Strand Information" Approach C. This enables strand bias computation without running jellyfish twice.
  ```rust
  struct StrandedCount {
      forward: u32,
      reverse: u32,
  }
  ```
- **Reference**: `counting-alternatives.md` for architecture and performance comparison, `jellyfish-deep-dive.md` for strand-aware counting design.

#### 5.3 SIMD-Accelerated K-mer Operations
- **Files**: `src/kmer/simd.rs`
- **Work**: SIMD-accelerate hot paths:
  - K-mer encoding (batch 2-bit encoding of nucleotide sequences)
  - Canonical form computation (parallel comparison of forward and reverse complement)
  - Hash function (if applicable to the chosen hash)
- Use `std::simd` (nightly) or portable SIMD via `packed_simd2` crate. Provide non-SIMD fallback for portability.

#### 5.4 Benchmark Suite
- **Files**: `benches/counting.rs`, `benches/walking.rs`, `benches/pathfinding.rs`
- **Work**: Criterion-based benchmarks for all hot paths. Compare native Rust counting against jellyfish on equivalent inputs. Measure walking throughput (targets per second) and pathfinding time.

### Dependencies
- Phase 4 complete (multi-k infrastructure requires counting at multiple k values)
- Access to jellyfish C++ source for hash function reverse engineering

### Estimated Scope
- ~8-10 source files
- ~3000-5000 lines of Rust code
- 6-8 weeks

### Milestone
**No external dependencies required. `kmerdet count` produces counts matching jellyfish within 0.01% tolerance. Strand-aware counting enables Fisher strand bias filter. Benchmark shows native counting within 1.5x of jellyfish speed on targeted panel data.**

---

## Phase 6: UMI Integration + Stretch

### Objective
Integrate UMI-aware k-mer counting, eliminating the need for external HUMID preprocessing. Implement ML-based filtering and structural variant extensions as stretch goals.

### Motivation
UMI-aware counting is the single most impactful improvement for lowering the VAF detection threshold per `room-for-improvement.md` Section 5 and `error-correction-alternatives.md`. Counting molecules instead of reads eliminates PCR bias at its source. A k-mer supported by 3 unique UMI families is far more credible than one with 20 reads from a single PCR family (both have the same read count but very different molecular evidence). Implementing this in kmerdet simultaneously solves three limitations: no UMI integration (#5), HUMID bottleneck (#7), and VAF threshold (#3, partially) per `room-for-improvement.md` cross-limitation interaction analysis.

### Components

#### 6.1 UMI-Aware K-mer Counting
- **Files**: `src/counting/umi.rs`, `src/counting/umi_counter.rs`
- **Work**: Implement UMI-aware k-mer counting directly from raw FASTQ per `error-correction-alternatives.md` Section "UMI-Aware K-mer Counting":
  - **Hybrid counter**: Exact counting for k-mers with <= 8 UMI families (SmallVec inline), HyperLogLog for high-count k-mers per `error-correction-alternatives.md` Option C. This keeps exact counts where precision matters most (variant k-mers at low VAF) while saving memory for reference k-mers.
  - **Dual count output**: Each k-mer has both `read_count` and `molecule_count`. Walking uses molecule_count for thresholding; quantification uses read_count.
  - **UMI clustering**: Port directional adjacency algorithm from rumi or RUMINA per `humid-analysis.md` Section "Alternatives to HUMID". Handle UMI errors (Hamming distance-based clustering with configurable -m parameter).
  - **Pipeline simplification**:
    ```
    Current:  FASTQ -> HUMID dedup -> dedup FASTQ -> jellyfish -> .jf -> kmerdet
    New:      FASTQ -> kmerdet count-and-detect -> results
    ```
  - **Benefits**: Eliminates ~2 min HUMID step, eliminates 10-15 GB intermediate files, enables molecule-level thresholding per `humid-analysis.md` Section "Could kmerdet Incorporate Its Own Dedup?"
- **Reference**: `error-correction-alternatives.md` Section "UMI-Aware K-mer Counting", `humid-analysis.md` for HUMID limitations and RUMINA as Rust-native alternative, `room-for-improvement.md` Section 5.

#### 6.2 K-mer-Level Duplex Evidence
- **Files**: `src/counting/duplex.rs`
- **Work**: Implement k-mer-level duplex evidence tracking per `error-correction-alternatives.md` Section "K-mer-Level Duplex Evidence":
  - For each k-mer, track forward_umi_count and reverse_umi_count separately
  - `duplex_support = min(forward_umi_count, reverse_umi_count)`
  - A k-mer with duplex_support >= 1 has cross-strand validation, dramatically reducing false positive probability
  - Integrate duplex_support into confidence scoring (add as feature to composite QUAL score)
- **Reference**: `error-correction-alternatives.md` Section "Duplex Consensus in an Alignment-Free Context".

#### 6.3 ML-Based Filtering (Optional)
- **Files**: `src/filter/ml.rs`, `src/filter/features.rs`
- **Work**: Implement optional ML-based post-filter per `adaptive-filtering.md` Section 6:
  - **Feature extraction**: rVAF, min_coverage, expression, variant_type, gc_content, homopolymer_length, linguistic_complexity, path_count, coverage_uniformity, pon_frequency, QUAL score, strand_bias, positional_cv, molecule_count.
  - **Model**: Gradient-boosted trees trained in Python (scikit-learn/XGBoost), exported as serialized decision tree for Rust evaluation. Pre-trained model shipped with kmerdet.
  - **Integration**: Optional stage 5 in composite filter pipeline. `--ml-filter` flag with configurable `--ml-threshold` (default 0.5, lower for high sensitivity). Always optional; rule-based filter remains default.
  - **Training data**: Simulated data + transfer learning from alignment-based callers per `adaptive-filtering.md` Section 6.
- **Reference**: `adaptive-filtering.md` Section 6, `10-future-improvements.md` Section "Machine Learning for Adaptive Thresholding".

#### 6.4 Structural Variant Extensions
- **Files**: `src/variant/structural.rs`
- **Work**: Extend variant detection to additional SV types per `10-future-improvements.md` Section "Structural Variant Detection Extensions":
  - **Gene fusions**: Targets designed across fusion partners; walking discovers fusion junction k-mers. Highest clinical priority.
  - **Large deletions (>50 bp)**: Extended walking with higher max_stack and bidirectional approach from Phase 3.
  - **Copy number changes**: Coverage-based detection from k-mer count profiles across targets.

#### 6.5 Longitudinal Multi-Sample Analysis
- **Files**: `src/pipeline/longitudinal.rs`
- **Work**: Support serial blood draw analysis per `10-future-improvements.md` Section "Multi-Sample Joint Calling":
  - Track rVAF trends over timepoints for each variant
  - Statistical test for rising vs. stable vs. declining VAF (relapse detection)
  - Borrow strength across timepoints: a variant consistently detected at very low VAF across 5 timepoints is more credible than one detected once

### Dependencies
- Phase 5 complete (native Rust counting infrastructure)
- ML model training data from Research Experiments

### Estimated Scope
- ~10-14 source files
- ~5000-7000 lines of Rust code
- 8-12 weeks

### Milestone
**Detection without HUMID preprocessing. Molecule-level k-mer counting with dual (read/molecule) count output. VAF detection threshold pushed from ~0.1% to ~0.05% through molecule counting + duplex evidence. ML filter optional but demonstrates improved FP/FN tradeoff on validation data.**

---

## Phase Dependency Summary

```
Phase 0 (Complete)
    |
    v
Phase 1 (Core Detection)
    |
    v
Phase 2 (Filtering + Output)
    |
    v
Phase 3 (Confidence + Adaptive)  <-- Research Experiments 1-9 run here
    |
    v
Phase 4 (Multi-k + Advanced)     <-- Research Experiment 3, 10 inform this
    |
    v
Phase 5 (Pure Rust + Performance)
    |
    v
Phase 6 (UMI Integration + Stretch)
```

Each phase depends on the previous phase. Research experiments run in parallel with development phases 3-4, informing parameter choices and validating improvements. The estimated total timeline is 8-12 months for a single developer working full-time.

---

## References

All research documents in `docs/research/`:
- `thesis-context/thesis-summary.md`, `thesis-context/room-for-improvement.md`
- `sensitivity-and-confidence/sensitivity-landscape.md`, `confidence-metrics.md`, `false-positive-analysis.md`
- `kmer-length-optimization/kmer-length-tradeoffs.md`, `multi-k-strategy.md`
- `error-correction/humid-analysis.md`, `error-correction-alternatives.md`
- `kmer-counting/jellyfish-deep-dive.md`, `counting-alternatives.md`
- `graph-and-walking/walking-improvements.md`, `graph-improvements.md`
- `filtering-strategies/adaptive-filtering.md`, `indel-filtering.md`
- `initial/10-future-improvements.md`
