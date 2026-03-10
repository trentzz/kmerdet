# Error Correction Alternatives for K-mer-Based Variant Detection

## Context

The kmerdet pipeline operates on sequencing data that has been deduplicated by HUMID but not error-corrected. Residual errors -- from PCR amplification, oxidative damage during library preparation, and sequencing itself -- create false positive k-mers that can produce spurious variant calls. This document evaluates error correction strategies applicable to the alignment-free k-mer detection workflow.

## Duplex Consensus in an Alignment-Free Context

### The Challenge

Traditional duplex consensus requires alignment to identify read pairs originating from the same genomic position on opposite strands. Without alignment, we need alternative methods to:

1. **Pair alpha/beta UMIs**: Identify which forward-strand UMI family corresponds to which reverse-strand family
2. **Determine strandedness**: Know which reads come from the Watson strand vs the Crick strand
3. **Call consensus**: Agree on a sequence for each strand, then require cross-strand agreement

### Alignment-Free UMI Pairing Approaches

**Approach A: Sequence-content clustering**

Instead of using genomic coordinates to group reads, cluster by sequence similarity:

```
1. Extract UMIs from all reads
2. For each UMI family, compute a "signature" from the read sequence
   (e.g., a MinHash sketch of the k-mers in the read)
3. Match UMI families whose signatures indicate overlapping genomic origin
4. Pair matched families as duplex partners
```

This approach works for targeted panels where reads from different targets are easily distinguishable by sequence content. Reads from the same target region will share most k-mers, enabling clustering without alignment.

**Approach B: K-mer overlap pairing**

```
1. For read pair R1/R2 (forward/reverse), the insert sequence is R1...R2_revcomp
2. Two duplex partners from the same molecule have complementary inserts
3. Compute a canonical k-mer set for each UMI family
4. Families with high Jaccard similarity in canonical k-mer space are duplex partners
```

Since jellyfish already uses canonical k-mers (storing min(kmer, reverse_complement)), reads from opposite strands produce the same canonical k-mers. This means duplex partners naturally share k-mer content and can be identified by k-mer set similarity.

**Practical feasibility**: For targeted panels with 50 targets, sequence-content clustering is straightforward because each target's reads are highly distinct. The challenge is computational cost: comparing all UMI family pairs is O(N^2). However, using MinHash or locality-sensitive hashing reduces this to near-linear time.

### K-mer-Level Duplex Evidence

Rather than forming duplex consensus at the read level, we can extract duplex-like evidence at the k-mer level:

```
For each k-mer K:
  forward_umi_count = number of UMI families where K appears in forward reads
  reverse_umi_count = number of UMI families where K appears in reverse reads

  duplex_support = min(forward_umi_count, reverse_umi_count)
  simplex_only = max(forward_umi_count, reverse_umi_count) - duplex_support
```

A k-mer supported by both forward and reverse UMI families has duplex-level evidence -- it was observed on both strands of at least one molecule. A k-mer supported on only one strand may be an artifact (e.g., oxidative damage creates G>T errors preferentially on one strand).

This approach provides a per-k-mer strand bias metric without requiring explicit duplex consensus calling, and integrates naturally with the k-mer counting paradigm.

## K-mer-Level Error Correction Tools

### Lighter (Song et al., 2014)

**Algorithm**: Lighter uses a pair of Bloom filters to identify and correct sequencing errors without explicit k-mer counting.

**Three-pass approach**:
1. **Sampling pass**: Subsample input k-mers at rate alpha = c/D (c = constant, D = estimated depth). Store sampled k-mers in Bloom filter A.
2. **Trusted k-mer identification**: A k-mer seen in Bloom filter A at the expected sampling rate is "trusted." Store trusted k-mers in Bloom filter B.
3. **Correction pass**: For each read, check all k-mers against Bloom filter B. If a k-mer is untrusted, try all single-base substitutions. If exactly one substitution produces a trusted k-mer, apply the correction.

**Performance**: Lighter is parallelized, uses no secondary storage, and requires only ~1-2 GB memory regardless of dataset size (because Bloom filter size is held constant by adjusting sampling rate).

**Applicability to targeted panels**: Lighter was designed for WGS data where the k-mer spectrum has a clear bimodal distribution (error k-mers at low counts, true k-mers at high counts). For targeted panel data:
- **Advantage**: Targeted data has very high depth per target, making the true/error distinction clearer
- **Caveat**: The sampling rate must be adjusted for the much higher depth (5000x vs 30x for WGS)
- **Caveat**: At very low VAF (0.1%), true variant k-mers have counts of ~5, which may be below the "trusted" threshold if the threshold is calibrated for the dominant reference coverage

**Risk**: Lighter could "correct" true variant k-mers to reference if the variant k-mer count is below the trusted threshold. For ctDNA detection, this would be catastrophic -- the error corrector would remove the signal we are trying to detect.

### BFC (Li, 2015)

**Algorithm**: BFC (Bloom Filter Corrector) uses a BWT (Burrows-Wheeler Transform) approach combined with Bloom filters for error correction. It constructs a k-mer count table using a blocked Bloom filter, then applies a correction strategy based on the k-mer spectrum.

**Key features**:
- Supports both substitution and indel error correction
- Uses a 2-pass approach: counting then correction
- High accuracy on benchmark datasets

**Applicability**: Like Lighter, BFC risks correcting true low-VAF variant k-mers. BFC has a minimum count threshold for "solid" k-mers; variants below this threshold would be eliminated.

### Musket (Liu et al., 2013)

**Algorithm**: Multi-stage k-mer spectrum-based error correction. Musket operates in three stages:
1. **Stage 1**: Two-sided k-mer extension correction (considers context on both sides of a suspect k-mer)
2. **Stage 2**: One-sided correction for errors near read ends
3. **Stage 3**: Voting-based correction for remaining errors

**Performance**: In a benchmark by AlEisa et al. (2022), Musket exhibited the best total effectiveness across six databases with varied error profiles. However, it is designed for moderate-coverage WGS, not ultra-high-depth targeted sequencing.

**Applicability**: The multi-stage approach is relevant to kmerdet because it uses contextual information (surrounding k-mers) to distinguish errors from true variation. This principle could be adapted for the walking algorithm: a k-mer that appears erroneous in the context of its neighbors in the de Bruijn graph can be flagged.

### CARE 2.0 (2022)

**Algorithm**: Uses machine learning (random forests) to reduce false-positive corrections. CARE 2.0 trains a model to predict whether a proposed correction is genuine, reducing the over-correction problem that afflicts simpler methods.

**Relevance**: The ML-guided approach is appealing for ctDNA analysis because it could be trained to preserve low-VAF variant k-mers while correcting true errors. However, it requires training data with known ground truth.

### Applicability Summary for Targeted Panel Data

| Tool | Designed for | Risk of correcting variants | Applicable to targeted panels? |
|------|-------------|---------------------------|-------------------------------|
| Lighter | WGS (30x) | High at low VAF | With caution; adjust thresholds |
| BFC | WGS (30x) | High at low VAF | With caution |
| Musket | WGS (moderate coverage) | Moderate | Partially; context-based approach is useful |
| CARE 2.0 | WGS (various) | Low (ML-guided) | Promising; needs retraining |

**Critical insight**: Standard error correction tools are designed for WGS at 30-50x coverage, where errors appear at count ~1 and true k-mers at count ~30. In targeted panels at 5000x, true variant k-mers at 0.1% VAF have counts of ~5 -- indistinguishable from moderate-frequency errors. Running these tools naively on targeted panel data would almost certainly eliminate true variant signal.

## UMI-Aware K-mer Counting

### Concept

The most promising error correction strategy for kmerdet is not to correct errors in reads, but to count k-mers by the number of distinct UMI families rather than by raw read count. This directly addresses PCR amplification bias and indirectly reduces error noise.

```
Traditional counting (jellyfish):
  kmer_count[K] = total number of reads containing k-mer K

UMI-aware counting:
  kmer_count[K] = number of distinct UMI families containing k-mer K
```

### Why This Works

Consider a scenario with 5000 raw reads covering a locus, collapsing to 500 UMI families after dedup:

- **True variant at 0.1% VAF**: Present in ~0.5 UMI families (expectation), but actually in 0 or 1 family (Poisson)
- **PCR error in one molecule amplified 20x**: Appears in 20 raw reads but only 1 UMI family
- **Sequencing error**: Appears in ~1 raw read, 1 UMI family

With raw counting: PCR error has count=20, far above typical thresholds. It would be called as a variant.
With UMI-aware counting: PCR error has count=1, same as sequencing errors. The threshold can be set based on molecular evidence rather than read evidence.

### Implementation Options

**Option A: Exact counting with HashMap**

```rust
use std::collections::{HashMap, HashSet};

// Key: k-mer (u64), Value: set of UMI hashes
let kmer_umis: HashMap<u64, HashSet<u64>> = HashMap::new();

for read in fastq_reader {
    let umi = extract_umi(&read);
    let umi_hash = hash_umi(umi);
    for kmer in extract_kmers(&read.sequence, k) {
        kmer_umis.entry(kmer)
            .or_insert_with(HashSet::new)
            .insert(umi_hash);
    }
}

// Final count: number of distinct UMIs per k-mer
let kmer_counts: HashMap<u64, u32> = kmer_umis.iter()
    .map(|(kmer, umis)| (*kmer, umis.len() as u32))
    .collect();
```

**Memory cost**: Each distinct k-mer stores a HashSet of UMI hashes. For a targeted panel with ~45M distinct k-mers and average ~10 UMI families per k-mer:
- K-mer key: 8 bytes
- HashSet overhead: ~64 bytes (for small sets)
- UMI hashes: 10 * 8 = 80 bytes
- Total per k-mer: ~152 bytes
- Total: 45M * 152 = ~6.8 GB

This is feasible but memory-heavy. Optimization: use a more compact set representation for small UMI counts (inline array for count <= 8, HashSet for larger).

**Option B: Approximate counting with HyperLogLog**

HyperLogLog (Flajolet et al., 2007) estimates cardinality (number of distinct elements) using only O(log log N) space per counter, with typical error of ~2% using 64 registers.

```rust
use hyperloglog::HyperLogLog;

// Key: k-mer (u64), Value: HyperLogLog counter for UMIs
let kmer_hll: HashMap<u64, HyperLogLog> = HashMap::new();

for read in fastq_reader {
    let umi = extract_umi(&read);
    for kmer in extract_kmers(&read.sequence, k) {
        kmer_hll.entry(kmer)
            .or_insert_with(|| HyperLogLog::new(6))  // 2^6 = 64 registers
            .insert(&umi);
    }
}

let kmer_counts: HashMap<u64, u32> = kmer_hll.iter()
    .map(|(kmer, hll)| (*kmer, hll.len() as u32))
    .collect();
```

**Memory cost**: Each HyperLogLog with 64 registers uses 64 bytes. For 45M distinct k-mers:
- K-mer key: 8 bytes
- HLL counter: 64 bytes
- Total per k-mer: 72 bytes
- Total: 45M * 72 = ~3.2 GB

This is approximately half the memory of exact counting, with ~2-5% error in the count. For variant detection where counts are small (1-50 molecules), the absolute error is less than 1 molecule, which is acceptable.

**However**, for very low counts (1-5 molecules at low VAF), HyperLogLog's relative error becomes significant. A count of 3 +/- 2% = 3.06, which rounds to 3 -- fine. But the quantization from HyperLogLog at these small cardinalities may introduce systematic bias. For the critical low-count regime, exact counting may be necessary.

**Option C: Hybrid approach**

Use exact counting for k-mers with low UMI counts (< threshold) and HyperLogLog for high-count k-mers:

```rust
enum UmiCounter {
    Exact(SmallVec<[u64; 8]>),  // Up to 8 UMIs stored exactly
    Approximate(HyperLogLog),    // Switch to HLL above threshold
}
```

This keeps exact counts where precision matters most (variant k-mers at low VAF) while saving memory for high-count reference k-mers.

### Replacing HUMID Entirely

If kmerdet performs UMI-aware k-mer counting directly from raw FASTQ, HUMID becomes unnecessary:

```
Current pipeline:
  FASTQ -> HUMID dedup -> deduped FASTQ -> jellyfish count -> .jf -> kmerdet detect

Proposed pipeline:
  FASTQ -> kmerdet count-and-detect -> results
```

Benefits:
1. **No intermediate files**: Eliminates ~10-15 GB intermediate FASTQ
2. **Single pass possible**: Read FASTQ once, extract UMIs and k-mers simultaneously
3. **Integrated UMI error correction**: Apply UMI clustering during counting
4. **Parallel from the start**: Use rayon for multithreaded counting (HUMID is single-threaded)
5. **True molecule counts**: More accurate than dedup-then-count

## Hybrid Approaches: Dedup + K-mer Error Correction

### Strategy: HUMID + Lighter/BFC + Jellyfish

Apply error correction between deduplication and k-mer counting:

```
FASTQ -> HUMID dedup -> Lighter correction -> jellyfish count -> kmerdet detect
```

**Expected improvement**: Lighter would correct random sequencing errors in the deduplicated reads before counting, reducing the noise floor of k-mer counts. The key question is whether Lighter's correction threshold can be tuned to avoid correcting true variant k-mers.

**Threshold tuning**: For targeted panel data at ~500x post-dedup coverage (from ~5000x raw depth with ~90% dedup rate):
- True reference k-mers: count ~500
- True variant at 0.1% VAF: count ~0.5 (expected) -> 0 or 1 (actual, Poisson)
- Sequencing error k-mers: count ~1-2

The overlap between variant and error k-mer counts is severe. Lighter cannot distinguish them without additional information (like UMI support or strand bias).

**Conclusion**: Lighter-style error correction is too aggressive for low-VAF detection. It may be useful as a first pass to reduce the most obvious errors (count=1 error k-mers that have a count=500 1-base-away neighbor), but must be parameterized conservatively to avoid removing variant signal.

### Strategy: HUMID + Context-Aware Correction

Rather than using coverage-based correction, use sequence context to identify likely errors:

1. For each low-count k-mer K, check its 1-Hamming-distance neighbors
2. If exactly one neighbor N has count >> count(K), and the base change matches a known error profile (e.g., G>T for oxidative damage), flag K as a probable error
3. Only correct/remove K if the evidence is overwhelming (count ratio > 100:1)

This is more conservative than Lighter and specifically targets error modes known to affect the sequencing chemistry in use.

## Lightweight Error Correction During Walking

### Tip Removal in de Bruijn Graphs

The de Bruijn graph assembly literature provides well-studied algorithms for removing sequencing errors from k-mer graphs. The key insight: errors create characteristic structures in the graph:

1. **Tips**: Short dead-end branches caused by errors near read ends
2. **Bulges**: Short alternative paths caused by internal read errors
3. **Chimeric edges**: Spurious connections caused by chimeric reads or contamination

The Velvet assembler (Zerbino & Birney, 2008) introduced systematic tip removal and bulge smoothing. SPAdes uses iteratively increasing coverage thresholds for these operations.

### Adapting Tip Removal for kmerdet

During k-mer walking, kmerdet already constructs a local de Bruijn graph around each target. Error-induced structures in this graph can be identified and removed:

```rust
/// Remove tips (dead-end branches) shorter than max_tip_length
/// with coverage below tip_coverage_threshold
fn remove_tips(graph: &mut KmerGraph, max_tip_length: usize, tip_coverage_ratio: f64) {
    let tips: Vec<NodeId> = graph.nodes()
        .filter(|node| {
            // Tip: has in-degree > 0 but out-degree = 0 (dead end)
            graph.out_degree(*node) == 0
            && graph.in_degree(*node) > 0
        })
        .collect();

    for tip_end in tips {
        // Walk backward from tip end to find the branch point
        let (branch_point, tip_path) = walk_back_to_branch(graph, tip_end);

        // Check tip length and coverage
        if tip_path.len() <= max_tip_length {
            let tip_coverage = min_coverage(&tip_path, graph);
            let main_path_coverage = max_sibling_coverage(branch_point, graph);

            if tip_coverage < main_path_coverage * tip_coverage_ratio {
                // Remove the tip
                for node in &tip_path {
                    graph.remove_node(*node);
                }
            }
        }
    }
}
```

### Bulge Smoothing

Bulges are short parallel paths between the same two nodes, where one path has much higher coverage:

```
                 (error path, low coverage)
    NODE_A  ──────────────────────────>  NODE_B
       \                                   /
        ──────────────────────────────────
              (correct path, high coverage)
```

For variant detection, care must be taken: a true low-VAF variant also creates a parallel path with lower coverage. The key distinguishing features:

| Feature | Error bulge | True variant |
|---------|------------|--------------|
| Path length difference | Often 0 (SNV-like error) | Varies (0 for SNV, >0 for INDEL) |
| Coverage ratio | Very high (>100:1) | Moderate (depends on VAF) |
| Base change profile | Matches known error modes | Random |
| UMI family support | Usually 1 family | Multiple families |
| Strand bias | Often single-strand | Both strands (for true variants) |

By checking UMI family support and strand representation, error bulges can be distinguished from true variants with high confidence. This integration of UMI-level and k-mer-level information is a unique advantage of the kmerdet approach.

### Error Branch Collapsing

A more aggressive approach: during graph construction, collapse low-coverage branches into their high-coverage parent if the branch differs by exactly one base (consistent with a substitution error):

```rust
fn collapse_error_branches(graph: &mut KmerGraph, min_ratio: f64) {
    for node in graph.nodes() {
        let children = graph.successors(node);
        if children.len() <= 1 { continue; }

        // Find the highest-coverage child
        let max_child = children.iter()
            .max_by_key(|c| graph.edge_weight(node, **c))
            .unwrap();
        let max_weight = graph.edge_weight(node, *max_child);

        // Collapse children with much lower coverage
        for child in &children {
            if child == max_child { continue; }
            let child_weight = graph.edge_weight(node, *child);

            if (child_weight as f64) < (max_weight as f64) * min_ratio {
                // This child is likely an error
                // Check: does it differ from max_child by exactly 1 base?
                let hamming = hamming_distance(
                    graph.kmer(*child), graph.kmer(*max_child)
                );
                if hamming == 1 {
                    // Collapse: remove error edge, add weight to correct edge
                    graph.add_edge_weight(node, *max_child, child_weight);
                    graph.remove_edge(node, *child);
                }
            }
        }
    }
}
```

### Expected Impact

Tip removal and bulge smoothing during walking would:
- Reduce the number of spurious variant paths found during graph exploration
- Lower the false positive rate without affecting the false negative rate (true variants have higher coverage than tips)
- Speed up pathfinding (smaller graph = faster Dijkstra)

The risk is removing true ultra-low-VAF variants that happen to look like tips. Safeguards:
- Never remove a branch with UMI family support >= 2
- Never remove a branch with reads from both strands
- Log all removed branches for review

## Recommended Error Correction Strategy for kmerdet

### Phase 1 (Current): No Correction

Continue using HUMID dedup + raw jellyfish counting. Focus on getting the core pipeline working correctly before adding error correction complexity. Baseline performance establishes what improvement error correction provides.

### Phase 2: Walking-Phase Error Correction

Implement tip removal and bulge smoothing during graph construction. This is self-contained within kmerdet (no external tools), low risk (conservative thresholds), and directly reduces false positives.

### Phase 3: UMI-Aware K-mer Counting

Replace HUMID + jellyfish with integrated UMI-aware k-mer counting in Rust. This provides the most fundamental improvement: counting molecules instead of reads eliminates PCR bias at its source. Combine with strand-aware counting for duplex-like evidence.

### Phase 4: ML-Guided Filtering

Train a classifier on features including: k-mer count, UMI family count, strand ratio, sequence context, neighboring k-mer coverage profile, and walking quality metrics. This replaces heuristic thresholds with learned decision boundaries, adapting to the specific error profile of each sequencing platform and library prep.

## References

- Song, L., Florea, L. & Langmead, B. (2014). Lighter: fast and memory-efficient sequencing error correction without counting. *Genome Biology*, 15, 509.
- Li, H. (2015). BFC: correcting Illumina sequencing errors. *Bioinformatics*, 31(17), 2885-2887.
- Liu, Y. et al. (2013). Musket: a multistage k-mer spectrum-based error corrector for Illumina sequence data. *Bioinformatics*, 29(3), 308-315.
- CARE 2.0 (2022). Reducing false-positive sequencing error corrections using machine learning. *Bioinformatics*, 38(12).
- Zerbino, D.R. & Birney, E. (2008). Velvet: Algorithms for de novo short read assembly using de Bruijn graphs. *Genome Research*, 18(5), 821-829.
- Flajolet, P. et al. (2007). HyperLogLog: the analysis of a near-optimal cardinality estimation algorithm. *AofA*, 127-146.
- Smith, T. et al. (2017). UMI-tools: modeling sequencing errors in Unique Molecular Identifiers. *Genome Research*, 27(3), 491-499.
- AlEisa, H.N. et al. (2022). K-mer Spectrum-Based Error Correction Algorithm for Next-Generation Sequencing Data. *Computational Intelligence and Neuroscience*, 2022.
- Salk, J.J. et al. (2018). Enhancing the accuracy of next-generation sequencing. *Nature Reviews Genetics*, 19, 269-285.
