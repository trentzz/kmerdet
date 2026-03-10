# Walking Algorithm Improvements

## 1. Current Thresholding Analysis

### How the Threshold Works Today

The k-mer extension threshold determines which child k-mers are followed during DFS walking. The current implementation in `src/walker/extension.rs` computes:

```
threshold = max(sum_of_sibling_counts * ratio, n_cutoff)
```

Where all four possible single-base extensions are queried, their counts summed, and any child whose count meets or exceeds the threshold is retained. This mirrors the original km Python logic in `Jellyfish.get_child`.

Two parameter regimes exist:

| Context | ratio | n_cutoff | Purpose |
|---------|-------|----------|---------|
| km standard (RNA-seq, high-VAF) | 0.30 | 500 | Aggressively prune noise at high coverage |
| Liquid biopsy (ctDNA, low-VAF) | 0.00001 | 2 | Detect variants at VAF 0.1% or below |

### Problem: Ratio Threshold at Low VAF

The ratio threshold is computed from the sum of sibling counts, which is dominated by the reference extension. At a locus with 5,000x coverage (typical for targeted ctDNA panels):

```
Reference child count: ~5000
Variant child count:   ~5  (at 0.1% VAF)
Error children:        ~1-2 each

sum_of_siblings = 5000 + 5 + 1 + 1 = 5007
ratio_threshold = 5007 * 0.00001 = 0.05 (rounds to 0)
threshold = max(0, 2) = 2
```

The variant k-mer at count=5 passes, but this is precarious. If coverage drops to 2,000x and VAF is 0.05%:

```
Variant count: 2000 * 0.0005 = 1
threshold = max(2000 * 0.00001, 2) = 2
```

The variant k-mer at count=1 is filtered. The fixed absolute minimum of 2 becomes the binding constraint, and genuine low-VAF variants are lost.

### Problem: Ratio Threshold in High-Coverage Hotspots

Some genomic regions have locally elevated coverage due to GC content, amplification bias, or mappability. If a reference k-mer has count 50,000 (10x the median), even a generous ratio of 0.00001 yields:

```
ratio_threshold = 50000 * 0.00001 = 0.5 (rounds to 0)
threshold = max(0, 2) = 2
```

This is fine numerically, but the problem is that high-coverage regions also have proportionally more sequencing errors. Error k-mers at count 3-4 now pass the threshold and pollute the graph with false branches.

The core issue: the threshold does not account for the expected error rate, which scales with coverage.

## 2. Adaptive Thresholds Based on Local Coverage

### Statistical Model

Instead of a fixed ratio, compute the threshold from a local coverage model that accounts for sequencing error rates. The key insight is that error k-mers follow a Poisson-like distribution, while genuine variant k-mers have counts proportional to coverage and VAF.

**Error k-mer count model:**

For a k-mer at a position with local coverage `C` and per-base error rate `e` (typically 0.1-1% for Illumina after UMI dedup):

```
Expected error k-mer count = C * e / 3
```

The division by 3 accounts for three possible erroneous bases. The actual count follows approximately a Poisson distribution:

```
P(error_count = n) = Poisson(n; lambda = C * e / 3)
```

**Threshold from the Poisson model:**

Set the threshold at the value where the probability of an error k-mer reaching that count is below a desired false extension rate (FER):

```
threshold(C, e, FER) = min { t : P(Poisson(C*e/3) >= t) < FER }
```

For example, with C=5000, e=0.001, FER=0.01:
```
lambda = 5000 * 0.001 / 3 = 1.67
P(X >= 5) = 1 - CDF_Poisson(4; 1.67) = 0.0038 < 0.01
threshold = 5
```

This approach draws from the work of Kelley et al. on "statistically solid k-mers" for error correction, where a Gamma distribution models erroneous k-mer frequencies and the 99.9th percentile defines the error-vs-solid boundary. The principle is the same: use the expected noise distribution to set a statistically principled threshold rather than an arbitrary ratio.

### Computing Local Coverage

Local coverage for a given k-mer position can be estimated from the reference k-mers in its neighborhood:

```rust
fn local_coverage(ref_counts: &[u64], position: usize, window: usize) -> f64 {
    let start = position.saturating_sub(window);
    let end = (position + window + 1).min(ref_counts.len());
    let slice = &ref_counts[start..end];
    // Use median rather than mean to resist outliers
    median(slice)
}
```

A window of 10-20 k-mer positions (covering 40-50 bp at k=31) provides a smooth local estimate while remaining responsive to genuine coverage variation.

### Estimating Error Rate

The per-base error rate can be estimated per-sample from the k-mer spectrum:

1. For each reference k-mer, query all four extensions
2. The three non-reference extensions are (mostly) errors
3. Compute: `e_estimate = sum(error_child_counts) / (3 * sum(reference_child_counts))`
4. Average across all reference k-mers to get a global `e` estimate

Alternatively, for UMI-deduplicated data where most PCR errors are removed, a fixed estimate of `e = 0.001` (0.1%) is reasonable for Illumina platforms. Post-deduplication error rates can be as low as 0.0001 for duplex-consensus reads.

### Implementation in kmerdet

The adaptive threshold replaces the fixed ratio/count parameters:

```rust
pub struct AdaptiveThresholdConfig {
    /// False extension rate (probability of following an error k-mer)
    pub false_extension_rate: f64,  // default: 0.01
    /// Per-base error rate (estimated or fixed)
    pub error_rate: f64,            // default: 0.001
    /// Window size for local coverage estimation (in k-mer positions)
    pub coverage_window: usize,     // default: 15
    /// Absolute minimum threshold (safety floor)
    pub min_threshold: u32,         // default: 2
}
```

The fixed `ratio` and `count` parameters remain available as a fallback for users who prefer deterministic behavior or backward compatibility with km parameters.

## 3. Bidirectional Walking

### Current: Forward-Only Walking

The current algorithm walks forward from each reference k-mer: given k-mer `S[1..k]`, try all extensions `S[2..k]+{A,C,G,T}`. This means variants are discovered by branching off the reference path to the right.

The `extend_backward` function already exists in `src/walker/extension.rs` but is not yet used in the walking loop. Backward extension prepends a base: try `{A,C,G,T}+S[1..k-1]`.

### Why Bidirectional Walking Helps

**Deletions at target boundaries:** If a deletion removes bases near the right end of the target, forward walking from the last reference k-mer before the deletion may not have enough remaining reference to rejoin. Walking backward from the first reference k-mer after the deletion can bridge the gap from the other direction.

```
Target:   =====[ref_kmers]=====[DELETED REGION]=====[ref_kmers]=====
Forward:  →→→→→→→→→→→→→→→→→→→→→ (can't bridge)
Backward:                         ←←←←←←←←←←←←←←←← (bridges the gap)
```

**Variants near anchors:** Variants within k positions of the target start are missed by forward walking because there are insufficient upstream reference k-mers to serve as starting points. Backward walking from the target end provides coverage.

**Asymmetric error profiles:** Some sequence contexts (e.g., GGC motifs on Illumina) have directional error patterns. Walking from both directions provides redundancy.

### Implementation Strategy

```rust
pub fn walk_bidirectional(
    db: &dyn KmerDatabase,
    ref_kmers: &[String],
    config: &WalkerConfig,
) -> WalkResult {
    // Forward pass: walk right from each reference k-mer
    let forward_nodes = walk_forward(db, ref_kmers, config);

    // Backward pass: walk left from each reference k-mer
    let backward_nodes = walk_backward(db, ref_kmers, config);

    // Merge: union of discovered nodes, keep max count for duplicates
    merge_walk_results(forward_nodes, backward_nodes)
}
```

The merged result feeds into graph construction identically to the current forward-only result. No changes to the graph or pathfinding code are needed.

### Deduplication After Bidirectional Walking

When both directions discover the same variant k-mers, the merge step must deduplicate. Since k-mer counts come from the jellyfish database (not from walking), the count for a given k-mer is the same regardless of which direction discovered it. The merge is simply a set union of k-mer strings.

## 4. Walking Strategies for Long INDELs

### The Long INDEL Problem

For insertions longer than ~15 bp (approaching k/2 for k=31), the variant path consists of many consecutive non-reference k-mers. Each must pass the threshold independently, and the walking algorithm must traverse all of them without returning to the reference.

With k=31, an insertion of length `L` introduces `L + k - 1` novel k-mers at the insertion junction. For L=50, that is 80 consecutive non-reference extensions that must all succeed. At low VAF, any one of these k-mers could drop below threshold and terminate the walk.

### Current Limitations

- `max_stack=500` limits DFS depth. For very long insertions (>100 bp) with branching, the stack may be exhausted before the path rejoins the reference.
- Each extension is evaluated independently against the threshold. There is no "momentum" - a single low-count k-mer breaks the chain.
- The threshold is computed from sibling counts that include the reference extension. During non-reference traversal, there may be no reference sibling, altering the threshold calculation.

### Proposed: Relaxed Thresholds During Non-Reference Traversal

When the walker is in a region of consecutive non-reference k-mers (suggesting an INDEL or structural variant), relax the threshold:

```rust
fn effective_threshold(
    config: &WalkerConfig,
    consecutive_non_ref: usize,
) -> (f64, u32) {
    if consecutive_non_ref > 5 {
        // Likely inside an INDEL - relax thresholds
        let relaxation = 0.5; // use half the normal threshold
        (config.ratio * relaxation, config.count.max(1))
    } else {
        (config.ratio, config.count)
    }
}
```

The relaxation engages only after 5+ consecutive non-reference extensions, reducing the risk of following errors at reference positions while allowing longer variant traversals.

### Proposed: Split Walking for Large INDELs

For targets with known large INDELs (from the variant catalog), walk from both ends of the target toward the middle:

```
Target:   [anchor_left]====[insertion site]====[anchor_right]
Walk 1:   →→→→→→→→→→→→→→→→→ (from left anchor, rightward)
Walk 2:   ←←←←←←←←←←←←←←←←← (from right anchor, leftward)
```

If both walks discover non-reference k-mers that overlap, the insertion path can be reconstructed even if neither walk alone would traverse the full insertion.

### Proposed: Guided Walking Using Partial Sequence Matching

When a target catalog specifies the expected variant sequence, use partial matches to guide the walk:

1. Decompose the expected variant sequence into k-mers
2. During DFS, if a candidate extension matches an expected variant k-mer, prioritize it (e.g., explore it first in the DFS order)
3. Do not use expected k-mers to bypass the threshold - this would bias toward expected variants and miss novel ones

This is optional guidance, not a filter. The walker still discovers all paths above threshold, but explores expected paths first for efficiency.

### Increased Stack Depth for INDEL-Suspected Targets

If target metadata indicates a large INDEL (e.g., FLT3-ITD can be 3-400 bp), increase `max_stack` dynamically:

```rust
let effective_max_stack = if target.expected_indel_size > 50 {
    config.max_stack + target.expected_indel_size * 2
} else {
    config.max_stack
};
```

## 5. Branch Pruning During Walking

### The Problem: Combinatorial Explosion in Repetitive Regions

Repetitive sequences generate many branching points. At homopolymer runs (e.g., AAAAAAA), each A->A extension also has plausible A->C, A->G, A->T children from sequencing errors, and the k-mers after returning to the reference are identical regardless of where the "error" occurred. This creates an exponential number of paths, most of which represent the same sequencing error at different positions.

The current `max_break=10` and `max_node=10000` limits provide hard caps, but these are blunt instruments that may terminate walking prematurely in complex but genuine variant regions.

### Tip Removal

Tips are short dead-end branches in the k-mer graph - a chain of nodes that extends a few steps then terminates without rejoining the reference or reaching another path.

From de Bruijn graph assembly (Velvet, SPAdes), tips shorter than `2k` are considered errors. In the kmerdet context, a tip of 2-3 extensions is almost certainly a sequencing error:

```rust
fn is_tip(branch: &[String], ref_kmers: &HashSet<String>, k: usize) -> bool {
    // A tip is a branch that terminates without rejoining reference
    // and is shorter than k/2 extensions
    let rejoins = branch.iter().any(|kmer| ref_kmers.contains(kmer));
    !rejoins && branch.len() < k / 2
}
```

Tips can be identified during walking (prune immediately) or after walking (remove from the node set before graph construction). Post-walking removal is simpler and allows the pruning criteria to use global information.

### Bubble Collapsing

Bubbles are pairs of paths that diverge from and rejoin the same nodes, differing by one or a few bases. In genome assembly, bubbles represent either heterozygous variants or sequencing errors. In our context, a bubble that differs by exactly one base and where one path has much lower coverage than the other is likely a sequencing error.

From BubbleGun and Velvet's Tour Bus algorithm, bubble detection works by:
1. Finding pairs of nodes with the same predecessors and successors
2. Comparing the sequences of the two paths
3. If paths differ by 1 base and one has <10% the coverage of the other, collapse to the higher-coverage path

```
      ┌─ ACGT (count=5000) ──┐
  ○───┤                      ├───○
      └─ ACGC (count=3)    ──┘
                ↑ error bubble → remove lower path
```

However, in variant detection we must be careful not to collapse genuine variant bubbles. The key distinguishing feature: error bubbles rejoin within 1-2 k-mer positions, while variant bubbles span `k` or more positions (at least `k` k-mers differ for an SNV).

**Rule:** Only collapse bubbles shorter than `k` k-mers. Longer bubbles are candidate variants and must be preserved.

### Quick-Rejoin Pruning

Branches that leave the reference and rejoin within `k` positions likely represent single-base errors that create a "shortcut" through k-mer space:

```
Reference: K1 → K2 → K3 → K4 → K5
Error:     K1 → Kx → K3 (rejoins after 1 step)
```

This pattern cannot represent a real variant (an SNV changes `k` consecutive k-mers). Prune any branch that rejoins within fewer than `k` positions.

## 6. Walking with Quality-Weighted Extensions

### The Information Deficit

Jellyfish databases store only k-mer counts, not quality information. A k-mer with count=10 could represent 10 high-quality observations or 10 low-quality observations. In the liquid biopsy setting, this distinction matters: UMI-deduplicated consensus reads have very high effective quality, but some consensus families are better supported than others.

### Count Ratio as Quality Proxy

Without per-base quality, the ratio of a child's count to the total sibling count serves as a proxy for confidence:

```
confidence(child) = count(child) / sum_of_siblings
```

A child with confidence 0.001 (1 in 1000 reads support it) is much less likely to be a true variant than one with confidence 0.01 (1 in 100). This is already implicitly used by the ratio threshold, but could be used more granularly:

- During DFS, explore higher-confidence branches first (best-first search element)
- When the graph becomes too large, prune lowest-confidence branches first (instead of hitting max_node and stopping entirely)

### Future: UMI-Aware Walking

If kmerdet implements UMI-aware k-mer counting (see `docs/research/initial/10-future-improvements.md`), each k-mer would have two counts: read count and molecule count. Walking decisions could use molecule count for thresholding while reporting read count for quantification:

```rust
struct KmerCounts {
    read_count: u64,
    molecule_count: u64,  // number of unique UMI families
}
```

A k-mer supported by 3 unique molecules (molecule_count=3) is more credible than one supported by 20 reads from a single UMI family (molecule_count=1, read_count=20). The molecule count is a direct measure of independent observations.

## 7. Performance Considerations

### Iterative vs. Recursive DFS

The km Python implementation uses recursive DFS (`__extend` calls itself), which is limited by Python's recursion limit and creates stack frames for each extension. The kmerdet implementation uses iterative DFS with an explicit stack (already planned in `src/walker/mod.rs`), which:

- Avoids stack overflow for deep traversals
- Allows easy implementation of max_stack as a stack size check
- Enables mid-traversal pruning (pop and skip)
- Is more cache-friendly than recursive calls

### Parallelizing Walking Across Targets

Each target is independently walkable (different region of the reference, same jellyfish database). Rayon's `par_iter` over targets provides natural parallelism:

```rust
targets.par_iter().map(|target| {
    let ref_kmers = target.decompose(k);
    walk_bidirectional(&db, &ref_kmers, &config)
}).collect()
```

Walking within a single target is inherently sequential (DFS), so parallelism is at the target level, not within a single walk.

### Memory Efficiency

The current `WalkResult` stores k-mer strings in a `HashMap<String, u64>`. For 31-mers, each entry is:
- String: 24 bytes (heap) + 31 bytes (content) + padding = ~64 bytes
- u64: 8 bytes
- HashMap overhead: ~32 bytes per entry

Total: ~104 bytes per discovered k-mer. For `max_node=10000`, this is ~1 MB per target, acceptable even with hundreds of targets.

For future optimization, 2-bit encoded k-mers (8 bytes per 31-mer in a u64) could replace String keys, reducing per-node overhead to ~48 bytes. This optimization is tracked in `src/kmer/` and is a Phase 1 task.

## References

- Audoux et al. (2017) -- km: RNA-seq investigation using k-mer decomposition
- Kelley et al. (2010) -- Quake: quality-aware detection and correction of sequencing errors (Poisson threshold model)
- Zerbino & Birney (2008) -- Velvet: de novo short read assembly using de Bruijn graphs (tip clipping, Tour Bus bubble merging)
- Onodera et al. (2016) -- Bidirectional Variable-Order de Bruijn Graphs
- Guo et al. (2024) -- Mining statistically-solid k-mers for accurate NGS error correction (Gamma distribution for error k-mer modeling)
- Bankevich et al. (2012) -- SPAdes: genome assembler using de Bruijn graphs (graph cleaning heuristics)
- Pockrandt et al. -- HUMID: UMI-aware deduplication for error-corrected sequencing
