# HUMID Analysis: Strengths, Limitations, and Alternatives

## Overview of HUMID

HUMID (High-performance UMI Deduplication) is a UMI-based deduplication tool developed by Pockrandt et al. that clusters reads by their UMI (Unique Molecular Identifier) sequences and outputs one representative read per UMI family. It was designed for speed and low memory usage in high-throughput sequencing workflows.

In the kam pipeline for ctDNA monitoring, HUMID is the first processing step after sequencing, responsible for collapsing PCR duplicates before k-mer counting with jellyfish.

## HUMID's Approach

### UMI Clustering via Hamming Distance

HUMID clusters reads based on their UMI sequences using Hamming distance as the similarity metric:

1. **Extract UMIs**: Parse UMI sequences from read headers or inline positions
2. **Cluster by Hamming distance**: Group UMIs within a configurable Hamming distance threshold (the `-m` flag)
3. **Select representative**: For each cluster, output one representative read (typically the first encountered or highest quality)
4. **Write output**: Produce a deduplicated FASTQ file

The Hamming distance approach is efficient (O(1) per comparison for fixed-length UMIs) but only accounts for substitution errors. It does not handle insertions or deletions within UMIs, which can occur during sequencing or library preparation.

### Error Tolerance Parameter (-m flag)

The `-m` flag controls the maximum Hamming distance between UMIs in the same cluster:

| -m value | Behavior | Trade-off |
|----------|----------|-----------|
| 0 | Exact UMI matching only | Misses UMIs with sequencing errors (under-dedup) |
| 1 | Allow 1 mismatch | Good balance for short UMIs (6-8 bp) |
| 2 | Allow 2 mismatches | May over-cluster for short UMIs |
| 3+ | Very permissive | Risk of merging distinct molecules (over-dedup) |

For a 8-bp UMI (4^8 = 65,536 possible UMIs), setting -m=1 means each UMI has 8*3 = 24 neighbors (each position can mutate to 3 other bases). At -m=2, each UMI has ~24 + 24*23 = 576 neighbors, significantly increasing collision risk.

### Output: One Read Per Family

HUMID outputs a single representative read per UMI family. This is deduplication, not consensus calling -- the representative read carries whatever errors it happened to contain. No information from other reads in the family is used to correct errors.

## What HUMID Does NOT Do

### No Duplex Consensus

True duplex sequencing uses complementary UMI pairs to tag both strands of a double-stranded DNA molecule. The protocol typically uses a dual-index scheme where:

- **Alpha UMI**: Tags the top strand
- **Beta UMI**: Tags the bottom strand
- **Duplex pair**: Reads with complementary alpha/beta UMIs originate from opposite strands of the same molecule

Duplex consensus requires:
1. Grouping reads by alpha UMI (forward strand family)
2. Grouping reads by beta UMI (reverse strand family)
3. Matching alpha and beta families (same original molecule)
4. Calling consensus separately within each single-strand family
5. Requiring agreement between forward and reverse strand consensus

HUMID performs only step 1 (and a simplified version at that). It does not:
- Pair alpha and beta UMIs to identify duplex families
- Call consensus within families
- Require cross-strand agreement

**Impact**: Without duplex consensus, single-strand errors (from PCR, oxidative damage during library prep, or deamination) survive deduplication. The fgbio best-practice pipeline (Fulcrum Genomics) demonstrates that duplex consensus suppresses errors by approximately 2 additional orders of magnitude compared to single-strand dedup alone. For ctDNA detection at 0.1% VAF, this difference is critical -- single-strand errors at ~0.1% frequency are indistinguishable from true variants without duplex filtering.

### No Consensus Calling Within Families

Even without duplex pairing, consensus calling within a single UMI family can correct random sequencing errors:

```
Family of 5 reads at same UMI:
  Read 1: ACGTACGT
  Read 2: ACGTACGT
  Read 3: ACGTACGT
  Read 4: ACCTACGT  <- sequencing error at position 3 (G->C)
  Read 5: ACGTACGT

Majority-rule consensus: ACGTACGT (error corrected)
HUMID output: Read 1 (error not corrected, but error in read 4 is discarded)
```

In this example, HUMID discards reads 2-5 entirely. If read 1 happened to contain the error, that error would persist. Consensus calling would correct it by majority vote.

The fgbio pipeline's `CallMolecularConsensusReads` tool performs this consensus calling, producing a per-base quality score that reflects the agreement across reads in the family. The `CallDuplexConsensusReads` tool extends this to require cross-strand agreement.

### No Positional Information

HUMID operates on FASTQ files (unaligned reads) and clusters solely by UMI sequence. It does not consider:
- Mapping position (reads with the same UMI but different genomic origins)
- Insert size
- Read pair concordance

Tools like UMI-tools (`dedup` command) use both UMI and mapping position to cluster reads, which reduces over-clustering when different molecules happen to share the same UMI (UMI collision).

## UMI Collision Problem in Targeted Panels

### The Mathematics of UMI Diversity

For a UMI of length N bases, there are 4^N possible unique UMI sequences:

| UMI length | Distinct UMIs | Sufficient for |
|------------|---------------|----------------|
| 6 bp | 4,096 | Low-depth WGS |
| 8 bp | 65,536 | Moderate-depth targeted |
| 10 bp | 1,048,576 | High-depth targeted |
| 12 bp | 16,777,216 | Ultra-deep sequencing |

### Collision Probability (Birthday Paradox)

The probability of at least one UMI collision among M molecules with U possible UMIs follows the birthday paradox approximation:

```
P(collision) ~= 1 - exp(-M^2 / (2 * U))
```

For a targeted panel at 5000x unique molecular depth per target region, with a typical target panel covering ~50 targets each ~200 bp:

- Total molecules: ~5000 * 50 = 250,000 unique molecules (rough upper bound per locus = 5000)
- With 8-bp UMI (U = 65,536): P(collision at a single locus) ~= 1 - exp(-5000^2 / (2 * 65536)) = ~1.0 (essentially certain)
- With 12-bp UMI (U = 16.7M): P(collision) ~= 1 - exp(-5000^2 / (2 * 16.7M)) = ~0.52 (50% chance)

At 5000x depth with 8-bp UMIs, collisions are almost guaranteed. Multiple distinct molecules will share the same UMI, and HUMID will incorrectly merge them into a single family. This has two consequences:

1. **Under-counting**: True molecular diversity is underestimated (inflates deduplication rate)
2. **Chimeric representatives**: The representative read may not accurately represent any single molecule

### Impact on Variant Detection

UMI collisions directly affect k-mer counting because HUMID outputs one read per (collided) family. If a UMI family contains both reference and variant molecules (because two different molecules share a UMI), the representative read will carry either the reference or variant allele, but not both. This effectively drops the variant molecule's k-mers from the count, reducing sensitivity.

Expected count loss at 0.1% VAF with 5000x depth:
```
True variant molecules: 5000 * 0.001 = 5
With 8-bp UMI collisions (~17% collision rate at per-molecule level):
  Expected variant molecules lost to collision: ~0.85
  Post-collision variant molecules: ~4.15
  Sensitivity reduction: ~17%
```

This 17% sensitivity loss is significant when operating near the detection threshold.

## Single-Threaded Performance Bottleneck

### Runtime Profile in the kam Pipeline

The thesis timing breakdown shows HUMID consuming ~2 minutes of the ~6 minute total pipeline. However, this understates the bottleneck because:

1. **HUMID is single-threaded**: It processes reads sequentially, using only one CPU core
2. **Jellyfish is multi-threaded**: It uses all available cores (typically 4-8)
3. **Walking is parallelized**: Multiple targets processed simultaneously

On a modern 8-core machine, HUMID uses 1/8 of available compute while consuming 1/3 of total wall-clock time. If HUMID were parallelized to use all 8 cores, its share would drop to ~15 seconds (assuming linear scaling), reducing total pipeline time to ~4.25 minutes.

### Disk I/O Overhead

HUMID reads the full FASTQ, processes it, and writes a deduplicated FASTQ. For a typical liquid biopsy sample:
- Input FASTQ: ~10 GB (compressed) / ~30 GB (uncompressed)
- Output deduplicated FASTQ: ~3-5 GB (compressed) / ~10-15 GB (uncompressed)

This I/O overhead could be eliminated by streaming deduplication directly into k-mer counting, avoiding the intermediate file entirely.

## Alternatives to HUMID

### fgbio (Fulcrum Genomics)

**Approach**: Alignment-based UMI grouping with molecular and duplex consensus calling.

| Feature | HUMID | fgbio |
|---------|-------|-------|
| Alignment required | No | Yes (BWA) |
| Duplex consensus | No | Yes (`CallDuplexConsensusReads`) |
| Single-strand consensus | No | Yes (`CallMolecularConsensusReads`) |
| UMI error correction | Hamming distance | Adjacency-based directional |
| Language | C++ | Scala (JVM) |
| Speed | Fast | Slow (alignment + JVM overhead) |
| Memory | Low | High |
| Per-base quality | No | Yes (consensus quality) |

fgbio's best-practice consensus pipeline (documented at `github.com/fulcrumgenomics/fgbio/docs/best-practice-consensus-pipeline.md`) is the gold standard for duplex sequencing analysis. The workflow:
1. Extract UMIs (`FastqToBam`)
2. Align raw reads (`bwa mem`)
3. Group by position + UMI (`GroupReadsByUmi`)
4. Call duplex consensus (`CallDuplexConsensusReads`)
5. Re-align consensus reads (`bwa mem`)
6. Filter consensus reads (`FilterConsensusReads`)

**Limitation for kmerdet**: fgbio requires alignment, which defeats the purpose of alignment-free k-mer analysis. However, its consensus-calling algorithms could be adapted for alignment-free use.

### UMI-tools

**Approach**: Alignment-based UMI deduplication with network-based error correction (Smith et al., 2017).

UMI-tools uses a directional adjacency method that accounts for UMI amplification patterns: a UMI with count N could have spawned UMIs with count <= N/2 + 1 that differ by 1 edit. This models PCR error more accurately than simple Hamming distance.

**Limitation**: Single-threaded and alignment-dependent. Recent benchmarks excluded UMI-tools due to its scalability issues at high depth.

### UMICollapse

**Approach**: Optimized UMI clustering using BK-trees and n-gram indexing (Liu, 2019).

UMICollapse achieves orders-of-magnitude speedup over UMI-tools by using efficient data structures for Hamming distance queries. It can deduplicate over 1 million unique UMIs at a single position in ~26 seconds.

**Relevance**: Could replace HUMID for the clustering step, but still does not provide consensus calling.

### gencore

**Approach**: Alignment-based UMI deduplication with consensus sequence generation and quality correction.

gencore generates consensus reads from UMI families, incorporating per-base quality scores. It is implemented in C++ and is faster than fgbio. However, it does not support duplex consensus.

### UMI-nea

**Approach**: Reference-free UMI deduplication that handles both substitution and indel errors in UMIs (published 2025).

UMI-nea is specifically designed to address the limitation of Hamming-distance-only approaches. It uses an edit-distance-based clustering that catches UMI indel errors that HUMID misses. This is relevant for sequencing platforms where indel errors in UMIs are non-negligible.

### RUMINA

**Approach**: Rust-based UMI deduplication with enhanced error correction (published 2026 in Bioinformatics).

RUMINA is particularly relevant because:
1. **Implemented in Rust**: Natural integration with kmerdet
2. **Supports multiple clustering strategies**: Directional, adjacency, and simple
3. **Majority-rule read selection**: Independent of mapping quality
4. **Performance**: Up to 10x faster than existing tools
5. **Improved low-frequency detection**: Benchmarked specifically on ultra-low frequency SNV detection (0.01-1%), the exact regime relevant to ctDNA monitoring

RUMINA could serve as the basis for kmerdet's built-in deduplication, eliminating the HUMID dependency entirely.

### rumi

**Approach**: Minimal Rust implementation of UMI-tools' directional adjacency algorithm (sstadick/rumi on GitHub).

rumi is a lightweight Rust port focused on the specific directional adjacency dedup algorithm from UMI-tools. It demonstrates that the core UMI clustering algorithms can be efficiently implemented in Rust with minimal dependencies.

## Could kmerdet Incorporate Its Own Dedup?

### Feasibility Assessment

Building UMI deduplication into kmerdet is feasible and potentially advantageous:

**Advantages**:
1. Eliminate external dependency (HUMID binary)
2. Stream dedup directly into k-mer counting (no intermediate FASTQ)
3. Enable UMI-aware k-mer counting (count molecules, not reads)
4. Full control over clustering parameters and error models
5. Parallelized from the start (Rust + rayon)

**Implementation sketch**:
```rust
// Phase 1: Read FASTQ, extract UMIs, cluster
let umi_families: HashMap<UmiCluster, Vec<ReadPair>> = cluster_by_umi(fastq_reader, max_dist);

// Phase 2: For each family, either select representative or call consensus
let consensus_reads: Vec<ConsensusRead> = umi_families.par_iter()
    .map(|(cluster, reads)| call_consensus(reads))
    .collect();

// Phase 3: Count k-mers directly from consensus reads (skip jellyfish)
let kmer_counts: DashMap<u64, u32> = count_kmers_parallel(&consensus_reads, k);
```

**Key library support in Rust**:
- `needletail` for FASTQ parsing (already used in kmerdet)
- `dashmap` for concurrent k-mer counting
- `rayon` for parallel processing
- Custom UMI clustering (port from RUMINA or rumi)

### Timeline and Priority

This is a Phase 3+ feature (see `10-future-improvements.md`). The immediate path uses HUMID as-is, with kmerdet reading the deduplicated FASTQ. A built-in dedup module would be developed alongside the UMI-aware k-mer counting feature, since both benefit from tight integration.

## References

- Pockrandt, C. et al. HUMID: High-performance UMI Deduplication. GitHub repository.
- Smith, T. et al. (2017). UMI-tools: modeling sequencing errors in Unique Molecular Identifiers to improve quantification accuracy. *Genome Research*, 27(3), 491-499.
- Fulcrum Genomics. fgbio Best Practice Consensus Pipeline. https://github.com/fulcrumgenomics/fgbio/blob/main/docs/best-practice-consensus-pipeline.md
- Liu, D. (2019). UMICollapse: Accelerating the deduplication and collapsing process for reads with Unique Molecular Identifiers. GitHub: Daniel-Liu-c0deb0t/UMICollapse.
- RUMINA (2026). High-throughput UMI deduplication for amplicon and whole-genome sequencing with enhanced error correction. *Bioinformatics*.
- sstadick/rumi. Rust UMI Directional Adjacency Deduplicator. GitHub repository.
- Salk, J.J. et al. (2018). Enhancing the accuracy of next-generation sequencing for detecting rare and subclonal mutations. *Nature Reviews Genetics*, 19, 269-285.
