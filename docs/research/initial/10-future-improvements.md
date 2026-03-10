# Future Improvements

## Overview

This document covers future improvements and stretch goals for kmerdet, drawn from thesis Chapter 6.5 recommendations, pipeline integration needs, and opportunities unique to the Rust reimplementation.

## From Thesis Ch6.5: Immediate Recommendations

### Parallel HUMID Integration

Currently, HUMID deduplication runs as a separate preprocessing step before jellyfish counting. The thesis recommends tighter integration:

- **Current**: `FASTQ -> HUMID dedup -> deduplicated FASTQ -> jellyfish count -> .jf`
- **Proposed**: `FASTQ -> kmerdet preprocess (HUMID + jellyfish in one pass) -> .jf`

Implementation approaches:
1. **Streaming integration**: Pipe HUMID output directly into jellyfish count via Unix named pipes (FIFOs)
2. **In-process deduplication**: Implement UMI-aware dedup in Rust, feed deduplicated reads directly to the k-mer counter
3. Eliminates disk I/O for intermediate deduplicated FASTQ files (can be >10 GB)
4. Estimated ~30% wall-clock improvement from overlapping dedup and counting

### Adaptive Filtering Parameters

Instead of fixed `ratio` and `count` values, automatically adjust based on sample characteristics:

- **Median k-mer coverage**: Higher coverage allows higher thresholds; lower coverage needs lower thresholds
- **Coverage coefficient of variation**: High variability warrants more permissive ratio
- **Deduplication rate**: Proxy for library complexity
- **GC bias profile**: Sequence-context-dependent adjustments

Example: At >5000x median coverage, set `count >= 5`; at <500x, set `count <= 2`. High coverage CV (>0.5) halves the effective ratio.

## K-mer Length Optimization

### Automatic Determination of Optimal k

| Factor | Shorter k (21-25) | Longer k (31-35) |
|--------|-------------------|-------------------|
| Specificity | Lower (more false matches) | Higher (more unique k-mers) |
| Large indel sensitivity | Better (more spanning k-mers) | Worse |
| Fragment length compatibility | Better for short cfDNA | May lose terminal k-mers |
| Homopolymer tolerance | Worse | Better |

**Algorithm**: For each target, find the minimum k (odd values only: 21, 23, ..., 35) such that anchor k-mers (first and last) are unique in the reference genome. Use shorter k for targets with known large indels, longer k for repetitive regions.

**Practical first step**: Allow per-target k specification in the target catalog, rather than a global k value.

## Nextflow Integration

Nextflow enables scalable, reproducible pipeline execution across local machines, HPC clusters, and cloud environments. kmerdet's role: a single binary handling detect + filter in one invocation (`kmerdet run`), minimizing process transitions.

**Requirements for Nextflow compatibility**:
- Clean exit codes (0 = success, non-zero = failure)
- All output written to files (not just stdout) for Nextflow's output channel collection
- Deterministic output (same inputs produce same outputs) for Nextflow's caching
- `--threads` support to work within Nextflow's resource allocation

## Pure Rust Jellyfish Reader

### Eliminating the C++ FFI Dependency

Replace FFI bindings to libjellyfish with a pure Rust reader using `memmap2` for memory-mapped I/O:

1. Parse the jellyfish binary header (JSON metadata: k-mer length, hash size, canonical flag, counter length)
2. Memory-map the hash table portion of the `.jf` file
3. Implement jellyfish's reversible hash function for key lookup
4. Handle variable-length counters (bit-level parsing for compact entries)
5. Linear probing with reprobe limit for collision resolution

**Challenges**:
- Jellyfish's hash function must be replicated exactly (not formally documented; reverse-engineer from C++ source)
- Variable-length counters pack small values tightly at the bit level
- Canonical k-mer representation must match jellyfish bit-for-bit

**Benefits**:
- No C++ toolchain required for building kmerdet
- No runtime dependency on libjellyfish shared library
- Memory-mapped I/O is optimal for random-access k-mer queries
- Potential for SIMD-accelerated hash lookup

## UMI-Aware K-mer Counting

Standard jellyfish counting treats each read independently. With UMI-based sequencing, reads sharing a UMI originate from the same original molecule.

**Proposed**: Count k-mers at the molecule level:
- `molecule_count[K]` = number of distinct UMI families containing k-mer K
- `read_count[K]` = total reads containing K (traditional count)

A k-mer at a true variant site supported by 3 unique molecules is more credible than one supported by 20 reads from a single PCR duplicate family. This directly addresses PCR amplification bias.

**Implementation options**:
- Custom k-mer counting data structure storing UMI-to-k-mer associations
- Bloom filter or HyperLogLog for approximate unique molecule counting per k-mer
- Requires more memory than standard counting, but eliminates separate dedup step

## Multi-Sample Joint Calling

Process multiple samples together to improve sensitivity and specificity:

| Mode | Use Case | Approach |
|------|----------|----------|
| **Longitudinal** | Serial blood draws, same patient | Borrow strength across timepoints |
| **Cohort** | Multiple patients, shared panel | Identify recurrent variants vs. artifacts |
| **Tumor-Normal** | Paired samples | Subtract normal k-mer counts to filter germline |

For tumor-normal mode, subtract normal k-mer counts from tumor counts before thresholding — effectively filtering germline variants and common artifacts without a separate reference file.

## Panel-of-Normals for Filtering

Build a database of k-mer paths from healthy control samples (N >= 20 recommended) to filter recurrent false positives:

1. Process healthy controls through the same pipeline
2. For each target, record non-reference k-mer paths found in normals
3. Store as a compact database: `(target, variant_path_hash) -> frequency_in_normals`
4. At detection time, flag any variant found in >5% of normals as a recurrent artifact (configurable via `--pon-max-freq`)

This catches systematic artifacts that are not sample-specific: sequencer-specific error patterns, reference genome errors, and common CHIP variants.

## Structural Variant Detection Extensions

| SV Type | K-mer Signature | Detection Approach |
|---------|----------------|--------------------|
| Large deletions (>50 bp) | Junction k-mers spanning breakpoint | Extended walking with higher max_stack |
| Tandem duplications | Repeated k-mer segments with elevated counts | ITD detection generalized |
| Inversions | Reverse-complement junction k-mers | Bidirectional walking at breakpoint |
| Gene fusions | k-mers spanning fusion junction | Targets designed across fusion partners |
| Copy number changes | Uniform shift in k-mer counts | Coverage-based detection from count profiles |

**Priority**: Gene fusions first (high clinical relevance, straightforward target design), then large deletions and duplications.

## Cloud/HPC Deployment Considerations

| Environment | Approach | Notes |
|-------------|----------|-------|
| Local workstation | Single binary | Default — no special configuration |
| HPC cluster (SLURM) | Binary + Nextflow | Nextflow handles job submission |
| AWS/GCP/Azure | Docker container | Minimal image: kmerdet binary + jellyfish |
| Kubernetes | Container + Argo/Tower | Workflow-engine managed |

**Cloud-specific features**:
- Support `s3://` and `gs://` paths for input jellyfish databases
- Checkpointing for spot/preemptible instances (save partial results to resume)
- Typical liquid biopsy analysis fits on a single 8-vCPU node in <6 minutes

## Real-Time Monitoring Dashboard

Optional web-based dashboard for pipeline monitoring (stretch goal):

- **Status endpoint**: kmerdet exposes lightweight HTTP endpoint (`--dashboard` flag) showing progress
- **Results viewer**: Browser interface for detection results, VAF plots, per-target status
- **Longitudinal tracking**: VAF trend plots over time for MRD monitoring
- **Technology**: `axum` HTTP server, minimal bundled HTML/JS (optional compile-time feature)

## Machine Learning for Adaptive Thresholding

Train a classifier to distinguish true variants from artifacts using features beyond raw k-mer count:

**Feature set**: rVAF, min_coverage, expression, ref_expression, variant_type, gc_content, homopolymer_length, kmer_complexity, path_count, graph_size, coverage_uniformity

**Approach**: Gradient-boosted trees via `smartcore` or `linfa` crate — lightweight and interpretable. Pre-trained model shipped with kmerdet, applied as optional post-filter (`--ml-filter`). Each call annotated with confidence score (0.0-1.0).

**Practical constraints**:
- 10-patient training set is small — augment with simulated data
- Model requires recalibration for different sequencing platforms/panels
- ML filter is always optional; deterministic rule-based filter remains default

## Implementation Priority

| Priority | Feature | Effort | Impact |
|----------|---------|--------|--------|
| P1 | Adaptive filtering parameters | Low | High |
| P1 | Panel-of-normals | Medium | High |
| P2 | Pure Rust jellyfish reader | High | High |
| P2 | K-mer length optimization | Medium | Medium |
| P2 | Nextflow integration | Low | Medium |
| P3 | UMI-aware k-mer counting | High | High |
| P3 | Multi-sample joint calling | Medium | Medium |
| P3 | Structural variant extensions | High | Medium |
| P4 | Real-time dashboard | Medium | Low |
| P4 | ML adaptive thresholding | High | Medium |

## References

- Thesis Chapter 6.5: Recommendations for Future Work
- memmap2 crate: https://docs.rs/memmap2/
- Nextflow documentation: https://www.nextflow.io/docs/latest/
- smartcore (Rust ML): https://docs.rs/smartcore/
- linfa (Rust ML): https://docs.rs/linfa/
- HyperLogLog: Flajolet et al. (2007)
