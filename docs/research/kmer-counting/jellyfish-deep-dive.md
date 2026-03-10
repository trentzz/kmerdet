# Jellyfish Deep Dive: Internals Relevant to Sensitivity

## Counter Saturation and Overflow

### Default Counter Size

Jellyfish uses variable-length counters controlled by the `-c` (or `--counter-len`) flag. The default counter length is 7 bits, supporting values 0 to 127 in the compact representation. However, this does not impose a maximum count -- jellyfish uses a variable-length encoding that spills into additional hash table entries when the counter overflows its initial allocation.

From the jellyfish manual:
> There is no maximum value in the hash. Even if the counting field uses 5 bits, a k-mer occurring 2 million times will have a value reported of 2 million.

### Variable-Length Counter Mechanism

The counter for each k-mer starts as a fixed-width field of `-c` bits stored alongside the k-mer key in a single hash table entry. When the count exceeds 2^c - 1, the counter spills into additional "overflow" entries in the hash table. Each overflow entry provides additional bits of counter value.

The process works as follows:

1. **Initial state**: K-mer and counter stored in one entry. Counter has c bits. Max value = 2^c - 1.
2. **Overflow**: When count reaches 2^c, an additional entry is allocated in the hash table to store the higher-order bits of the counter.
3. **Chaining**: Multiple overflow entries can be chained for very large counts.

The trade-off from the `-c` parameter:
- **Low c (e.g., 3-5 bits)**: Saves space per entry for low-count k-mers, but high-count k-mers consume more entries (overflow entries).
- **High c (e.g., 12-16 bits)**: Wastes space on low-count k-mers (which are the majority), but high-count k-mers fit in a single entry.

For targeted panel sequencing at 5000x depth, a reference k-mer might have count ~5000. With the default c=7 (max 127), this requires overflow entries. With c=13 (max 8191), most k-mers fit in a single entry.

### Memory Formula

The memory usage per hash table entry (in bits) is:

```
bits_per_entry = 2k - l + r + 1
```

Where:
- k = k-mer length (in bases; each base = 2 bits, so 2k bits for the key)
- l = log2(table_size) -- bits saved by the table size itself (position in table encodes l bits of the key)
- r = log2(max_reprobe) -- bits used to store the reprobe distance
- 1 = extra bit for "busy" flag (used for lock-free synchronization)

For k=31, table size = 2^26 (~67M entries), max reprobe = 126 (7 bits):
```
bits_per_entry = 62 - 26 + 7 + 1 = 44 bits = 5.5 bytes
```

Plus the counter field: 7 bits default = 51 bits total per entry, or ~6.4 bytes.

### Impact on High-Coverage Regions

In targeted panel sequencing, certain regions may have extremely high coverage (>10,000x) due to capture efficiency variation. At these coverage levels:

1. **Counter overflow is guaranteed**: Counts of 5000-10,000 overflow a 7-bit counter (max 127), requiring 1-2 additional entries per k-mer.
2. **Hash table fill increases**: Overflow entries consume hash table slots, increasing the effective fill factor.
3. **No precision loss**: Counts remain exact regardless of overflow. There is no saturation or capping.
4. **Recommendation**: For targeted panel data, use `-c 13` or `-c 16` to avoid overflow overhead:
   ```bash
   jellyfish count -m 31 -s 100M -t 4 -C -c 16 -L 2 reads.fastq
   ```
   With c=16 (max 65,535), counts up to 65K fit in a single entry. This wastes ~9 bits per low-count k-mer but avoids the overhead of overflow entries for the many high-count reference k-mers in targeted data.

## Canonical Counting Implications

### How Canonical Mode Works

With the `-C` flag (the standard and recommended mode), jellyfish stores each k-mer in its canonical form: the lexicographically smaller of the k-mer and its reverse complement.

```
K-mer:              ACGTACGTACGTACGTACGTACGTACGTACG  (31-mer)
Reverse complement: CGTACGTACGTACGTACGTACGTACGTACGT

Canonical = min(K-mer, RC) = ACGTACGTACGTACGTACGTACGTACGTACG
(stored as this sequence, regardless of which strand was sequenced)
```

Both forward and reverse strand reads contribute to the same counter. This means:

- The count for a canonical k-mer reflects **total coverage from both strands**
- A locus with 2500 forward reads and 2500 reverse reads produces a count of 5000 for the canonical k-mer

### Advantages for Variant Detection

1. **Higher effective count**: Both strands contribute, doubling the count compared to strand-specific counting. Higher counts improve the signal-to-noise ratio.
2. **Simplified analysis**: One counter per genomic position (modulo repeat ambiguity) rather than two.
3. **Memory savings**: Hash table is approximately half the size.

### Loss of Strand-Specific Information

**The critical limitation**: Canonical counting discards strand identity. This means:

1. **Cannot compute strand bias from JF alone**: Strand bias -- the ratio of forward to reverse strand support -- is a powerful quality metric for variant calling. True variants should have roughly equal support from both strands (strand ratio ~0.5). Single-strand artifacts (e.g., oxidative damage creating G>T on one strand) show extreme strand bias (ratio near 0 or 1).

2. **No duplex-level evidence**: As discussed in the error correction documents, duplex sequencing requires separate strand information. Canonical counting merges the strands.

3. **Cannot detect strand-specific errors**: Some library preparation protocols introduce strand-specific artifacts (e.g., FFPE-derived DNA has cytosine deamination artifacts preferentially on one strand). These cannot be identified from canonical counts.

### Workarounds for Strand Information

To recover strand information while using jellyfish, two approaches are possible:

**Approach A: Run jellyfish twice (canonical and non-canonical)**:
```bash
# Canonical (both strands merged)
jellyfish count -m 31 -s 100M -t 4 -C reads.fastq -o counts_canonical.jf

# Non-canonical (strand-specific)
jellyfish count -m 31 -s 100M -t 4 reads.fastq -o counts_noncanonical.jf
```

Then: forward_count = query(kmer, noncanonical); total_count = query(kmer, canonical); reverse_count = total_count - forward_count.

**Approach B: Count strands in kmerdet directly**:
When kmerdet reads the FASTQ for k-mer walking, it can track which strand each k-mer came from. This requires parsing the FASTQ rather than relying solely on the .jf database, but provides per-k-mer strand counts for quality scoring without running jellyfish twice.

**Approach C: Strand-aware k-mer counting in Rust**:
If kmerdet implements its own k-mer counting (see `counting-alternatives.md`), it can maintain separate forward and reverse counters per canonical k-mer:
```rust
struct StrandedCount {
    forward: u32,
    reverse: u32,
}
// Total count = forward + reverse (equivalent to canonical counting)
// Strand bias = forward / (forward + reverse)
```

### Implication for rVAF Quantification

The NNLS (non-negative least squares) quantification uses k-mer counts to estimate the contribution of each path. With canonical counting, the counts are higher (both strands) but the contribution matrix must also account for both strands. In practice, canonical counting works correctly for rVAF estimation because the contribution matrix is constructed from canonical k-mers. However, strand bias information, if available, could improve the quantification by down-weighting k-mers with extreme strand bias.

## Lower-Count Cutoff (-L Flag)

### How -L Works

The `-L` flag sets a lower count threshold during jellyfish counting. K-mers with count below this threshold are discarded from the database:

```bash
jellyfish count -m 31 -s 100M -t 4 -C -L 2 reads.fastq
```

With `-L 2`, all k-mers appearing only once (singletons) are discarded. This has profound implications for the database:

- **WGS at 30x**: ~70-80% of distinct k-mers are singletons (sequencing errors). `-L 2` removes most of these, reducing database size by ~70%.
- **Targeted panel at 5000x**: Fewer singleton k-mers in the target regions (most are counted many times), but off-target k-mers may be singletons.

### Impact on Low-VAF Detection

The critical question: what is the expected count of a true variant k-mer at very low VAF?

```
expected_count(VAF, depth, k, read_length) = depth * VAF * (1 - (k-1) / read_length)
```

For 150bp reads, k=31:

| VAF   | Depth (post-dedup) | Expected k-mer count | Survives -L 2? |
|-------|--------------------|---------------------|-----------------|
| 1%    | 500                | 500 * 0.01 * 0.80 = 4.0  | Yes |
| 0.5%  | 500                | 500 * 0.005 * 0.80 = 2.0 | Borderline |
| 0.1%  | 500                | 500 * 0.001 * 0.80 = 0.4 | **No** (expected < 1) |
| 0.1%  | 5000               | 5000 * 0.001 * 0.80 = 4.0 | Yes |
| 0.05% | 5000               | 5000 * 0.0005 * 0.80 = 2.0 | Borderline |
| 0.01% | 5000               | 5000 * 0.0001 * 0.80 = 0.4 | **No** |

At 500x post-dedup depth (typical after HUMID), variants at 0.1% VAF have expected k-mer counts below 1. With `-L 2`, these variants are discarded before detection even begins.

At 5000x raw depth (before HUMID), the situation is better: 0.1% VAF yields expected counts of ~4, well above `-L 2`. But after HUMID deduplication (which typically collapses ~90% of reads), the effective depth drops to ~500x, pushing low-VAF signals dangerously close to the threshold.

### Recommendations

1. **Use `-L 1` (or omit -L entirely) for ctDNA detection**: The noise floor from singletons can be handled during the walking phase (thresholding during extension) rather than at the counting phase. Removing singletons at counting time risks discarding true signal.

2. **Account for dedup rate when choosing -L**: If counting is performed on deduplicated reads (post-HUMID), the effective depth is lower, and `-L` should be more conservative (lower or disabled).

3. **Use -L for WGS but not for targeted panels**: For WGS at 30x, `-L 2` removes the vast majority of error k-mers with minimal signal loss. For targeted panels at 5000x, the benefit of `-L 2` is smaller (fewer error singletons as a fraction) and the risk is higher (variant k-mers near count 1-2).

### Alternative: Post-Counting Filtering

Rather than discarding low-count k-mers during counting, keep all k-mers and filter during the walking phase:

```rust
fn get_children(&self, kmer: Kmer, db: &KmerDatabase) -> Vec<(Kmer, u32)> {
    let children: Vec<(Kmer, u32)> = kmer.extend_all()
        .iter()
        .map(|child| (*child, db.query(*child)))
        .filter(|(_, count)| *count > 0)  // Keep all non-zero counts
        .collect();

    // Apply relative threshold during walking, not absolute threshold
    let sum: u32 = children.iter().map(|(_, c)| *c).sum();
    let threshold = max(
        (sum as f64 * self.ratio) as u32,
        self.min_count  // This is the effective -L, applied contextually
    );

    children.into_iter()
        .filter(|(_, count)| *count >= threshold)
        .collect()
}
```

This approach defers the count threshold to the walking phase, where it can be applied relative to the local coverage (via the `ratio` parameter) rather than as a global absolute cutoff.

## Hash Collision Behavior at High Load

### Open Addressing with Linear Probing

Jellyfish uses open addressing with a probing strategy to handle hash collisions. When a k-mer hashes to an occupied slot, jellyfish probes subsequent slots until it finds the k-mer (existing entry) or an empty slot (new entry).

The probe sequence is not simple linear probing -- jellyfish uses a "reprobing" scheme where the step size is derived from the k-mer itself, providing a form of double hashing.

### Reprobe Limit

The `--reprobe` parameter (or compiled default) sets the maximum number of probes before declaring a hash table miss. If a k-mer requires more than `reprobe_limit` probes, the insertion fails.

The default reprobe limit in jellyfish 2 is 126. This means:
- If the hash table is sparsely loaded, most insertions succeed in 1-3 probes
- As the table fills, average probe count increases
- When the table is >80% full, probe counts can reach the limit, causing failures

### What Happens When the Table is Undersized (-s Flag)

The `-s` flag sets the initial hash table size. If the actual number of distinct k-mers exceeds the table capacity (accounting for the ~80% load factor), jellyfish's behavior depends on the version and configuration:

**Jellyfish 1.x**: Creates multiple intermediate files and requires a separate `jellyfish merge` step to combine them. Each file contains a full hash table; when one fills, a new one is started.

**Jellyfish 2.x**: Uses a more sophisticated approach with dynamic resizing. However, undersizing still causes:
1. **Higher probe counts**: More collisions degrade insertion and query performance
2. **Overflow files**: Excess k-mers may be written to secondary files
3. **Potential count inaccuracy**: If the same k-mer ends up in multiple files before merging, its count may be split

**Impact on count accuracy**: In practice, jellyfish 2 handles undersized tables gracefully by creating multiple output files that can be merged. The final counts are accurate, but:
- Wall-clock time increases due to merge overhead
- Memory usage spikes during merge
- For the walking algorithm, querying a merged database vs. multiple unmerged databases may produce different results if the merge is incomplete

**Recommendation for targeted panel data**:
```bash
# Estimate distinct k-mers conservatively
# For 50 targets, each ~200bp, at 5000x depth:
# ~50 * 200 * 5000 = 50M reads * 120 k-mers each = 6B total k-mers
# But distinct k-mers are much fewer: ~50M for targeted data
# Use 2x margin: -s 100M

jellyfish count -m 31 -s 100M -t 4 -C reads.fastq
```

Oversizing the table by 2-3x wastes ~30% memory but ensures optimal performance with minimal probing.

### Cache Behavior for K-mer Walking

K-mer walking involves a specific access pattern: starting from an anchor k-mer, extending one base at a time (4 queries per extension), following a path through the hash table.

The access pattern has two relevant properties:

1. **Spatial locality**: Consecutive k-mers in a walk differ by one base (k-1 overlap). Their hash values are completely different (hash functions are designed to spread similar inputs). This means **no spatial locality** in the hash table -- each step in the walk jumps to a random location.

2. **Temporal locality**: The same k-mers may be queried multiple times (once during forward walking, once during backward walking, once during pathfinding). This provides some temporal locality, but the time between re-queries may be long (thousands of other queries in between).

For a hash table of 100M entries at 6 bytes each = 600 MB, the entire table does not fit in L3 cache (typically 8-32 MB). Each query in the walking phase is essentially a random access to a 600 MB structure, resulting in cache misses.

**Optimization opportunities**:
- **Memory-mapped I/O**: The OS page cache handles caching transparently. Frequently accessed pages stay in memory.
- **Prefetching**: If the walking algorithm can predict which k-mers it will query next (which it can -- the 4 possible extensions), it can issue prefetch instructions for those hash table locations while processing the current k-mer.
- **Batch queries**: Group k-mer queries from multiple concurrent walks to amortize cache miss latency.

## Memory-Mapped vs In-Memory Access

### JF File Sizes

For targeted panel sequencing after HUMID deduplication, typical .jf file sizes:

| Data type | Distinct k-mers | JF file size |
|-----------|-----------------|--------------|
| Targeted panel (50 targets, 5000x) | ~20-50M | 200-500 MB |
| Whole exome (5000 genes) | ~200-500M | 2-5 GB |
| WGS (30x) | ~2-3 billion | 15-25 GB |

For targeted panels, the entire .jf file fits comfortably in RAM (most workstations have 16+ GB). For WGS, memory mapping is essential.

### Memory Mapping for Random Access

Using `memmap2` in Rust:

```rust
use memmap2::Mmap;
use std::fs::File;

let file = File::open("counts.jf")?;
let mmap = unsafe { Mmap::map(&file)? };

// Random access via byte offset
let value = &mmap[offset..offset + entry_size];
```

Memory mapping provides:
- **Lazy loading**: Only pages accessed are loaded into memory
- **OS-managed caching**: The kernel manages page cache efficiently
- **Shared access**: Multiple processes can share the same mapping
- **No explicit I/O**: Reads look like memory accesses

### Access Pattern Analysis

For the k-mer walking use case:

```
Walking step 1: Query kmer_A -> hash_A -> page_X
Walking step 2: Query kmer_B -> hash_B -> page_Y (random, different page)
Walking step 3: Query kmer_C -> hash_C -> page_Z (random, different page)
...
```

Each walking step accesses a different, effectively random page. For a 500 MB .jf file with 4 KB pages, there are ~125,000 pages. A typical walk visits ~200 k-mers (for a 200 bp target), accessing ~200 pages (assuming all different). This represents 200/125,000 = 0.16% of pages, which fits easily in the page cache.

For 50 targets walking concurrently, ~10,000 page accesses are needed, still a small fraction of total pages. The working set is manageable.

### Recommendation for kmerdet

For targeted panel data, load the entire .jf file into memory (via mmap with `madvise(MADV_WILLNEED)` to pre-fault all pages):

```rust
use memmap2::MmapOptions;

let file = File::open("counts.jf")?;
let mmap = unsafe {
    MmapOptions::new()
        .populate()  // Pre-fault all pages
        .map(&file)?
};
```

This ensures all queries are served from RAM with no page faults during the walking phase. For a 500 MB file, this is a trivial memory cost on modern systems.

For WGS or very large databases, use standard mmap without pre-population and rely on the OS page cache.

## Jellyfish Hash Function

### Reversible Hash

Jellyfish uses a reversible hash function for its hash table. The key insight: because the hash function is reversible (bijective), the original k-mer can be recovered from the hash value plus the table position. This means the hash table does not need to store the full k-mer key -- it stores only the bits of the k-mer that are not encoded by the table position.

For a table of size 2^l, the table position encodes l bits of the k-mer. The remaining (2k - l) bits must be stored explicitly. The reprobe offset uses r additional bits. This gives the memory formula:

```
bits_per_entry = (2k - l) + r + 1 + c
               = key_remainder + reprobe_offset + busy_bit + counter
```

Where c is the counter bits (-c flag).

### Lock-Free CAS Protocol

The insertion protocol uses Compare-and-Swap (CAS) atomics:

1. Hash the k-mer to get a table position
2. Read the entry at that position atomically
3. If empty: CAS the entry from empty to (key_remainder | counter=1 | busy=1)
4. If occupied with same key: CAS to increment counter
5. If occupied with different key: reprobe (try next position in probe sequence)
6. If CAS fails (another thread modified the entry): retry from step 2

This lock-free protocol enables near-linear scaling with CPU threads, which is jellyfish's primary performance advantage over sequential counters.

### Implications for Rust Reimplementation

To read .jf files in pure Rust (without the C++ FFI), kmerdet must:

1. Parse the JSON header to extract k, table_size, counter_len, canonical flag, and max_reprobe
2. Implement the same reversible hash function (must match the C++ implementation bit-for-bit)
3. Implement the probing sequence
4. Extract the key remainder and counter from each entry (bit-level parsing)
5. Reconstruct the full k-mer from table position + key remainder via the inverse hash

The hash function is not formally documented in the jellyfish paper. It must be reverse-engineered from the C++ source code (`include/jellyfish/large_hash_array.hpp` and related headers). This is the primary challenge for a pure Rust .jf reader.

## Summary of Sensitivity-Relevant Parameters

| Parameter | Default | Sensitivity impact | Recommendation for ctDNA |
|-----------|---------|--------------------|--------------------------|
| `-m` (k-mer length) | 31 | Larger k = more specific, less INDEL sensitive | 31 for SNV, 21 for large INDEL (see multi-k docs) |
| `-C` (canonical) | Off | Canonical doubles count, loses strand info | Use -C for counting; track strands separately |
| `-L` (lower count) | 0 | Higher -L discards low-count variant k-mers | Use -L 1 or omit for low-VAF detection |
| `-s` (hash size) | Required | Undersized table may lose counts | 2-3x estimated distinct k-mers |
| `-c` (counter bits) | 7 | Low -c causes overflow at high coverage | Use -c 16 for targeted panel (5000x) |
| `-t` (threads) | 1 | Performance only, no accuracy impact | Use all available cores |

## References

- Marcais, G. & Kingsford, C. (2011). A fast, lock-free approach for efficient parallel counting of occurrences of k-mers. *Bioinformatics*, 27(6), 764-770.
- Jellyfish documentation. https://github.com/gmarcais/Jellyfish/blob/master/doc/Readme.md
- Jellyfish manual v1.1. https://www.cbcb.umd.edu/software/jellyfish/jellyfish-manual-1.1.pdf
- Jellyfish source code. https://github.com/gmarcais/Jellyfish (include/jellyfish/large_hash_array.hpp)
