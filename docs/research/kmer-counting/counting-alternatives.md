# K-mer Counting Alternatives and Their Tradeoffs

## Overview

Jellyfish is the current k-mer counter used in the kam/kmerdet pipeline. This document evaluates alternative counting approaches, from established tools (KMC3, DSK, Squeakr) to novel strategies (streaming counting, approximate sketches, and native Rust counting). The evaluation focuses on the specific requirements of targeted panel sequencing for ctDNA detection: exact counts at low frequencies, moderate data volumes, and integration with an alignment-free variant detection pipeline.

## KMC3

### How It Works

KMC3 (Kokot et al., 2017) uses a disk-based approach that partitions the k-mer space to enable counting with bounded memory:

1. **Signature-based partitioning**: Extract a short "signature" (prefix) from each k-mer. Partition k-mers into bins based on their signature.
2. **Disk buffering**: Write k-mers to disk-based bins (temporary files) organized by signature partition.
3. **Per-partition counting**: Load one partition at a time into a compact in-memory hash table, count all k-mers in that partition, and output the counts.
4. **Merge**: Concatenate partition outputs into the final k-mer count database.

The key innovation is the two-stage pipeline: a "splitter" stage reads the input and distributes k-mers to disk partitions, while a "counter" stage processes partitions independently. Both stages run concurrently, with disk I/O buffered between them.

### Performance Characteristics

| Metric | KMC3 | Jellyfish 2 |
|--------|------|-------------|
| Counting speed (human WGS, 30x) | ~30-45 min | ~45-60 min |
| Peak memory | 4-12 GB (configurable) | Proportional to distinct k-mers |
| Disk usage | Significant (temporary files) | Minimal |
| Thread scaling | Near-linear | Near-linear |
| Max k supported | Arbitrary (>200) | Limited by hash entry size |
| Accuracy | Exact | Exact |

In comprehensive benchmarks (Manekar & Sathe, 2018), KMC3 was the fastest tool for most k values and dataset sizes, with DSK as the only other tool that consistently produced correct results across all tested configurations.

### Advantages Over Jellyfish for kmerdet

1. **Lower memory ceiling**: KMC3 can count k-mers using a fixed, user-specified memory budget. For resource-constrained environments, this is advantageous.
2. **Arbitrary k support**: KMC3 handles k>32 natively, while jellyfish requires multiple machine words for k>31 (62 bits > 64-bit word). If kmerdet needs k=41 or k=43 for multi-k strategies, KMC3 handles this more naturally.
3. **Set operations**: KMC3 includes tools for set operations (union, intersection, subtraction) on k-mer databases. This could support panel-of-normals filtering.

### Disadvantages

1. **Disk I/O dependency**: KMC3 writes significant temporary data to disk. On systems with slow storage, this can become a bottleneck. For targeted panel data (small input), this overhead is minimal.
2. **Output format**: KMC3 produces its own binary format, not .jf files. kmerdet would need a KMC3 reader in addition to (or instead of) the jellyfish reader.
3. **Query interface**: KMC3's query API is different from jellyfish's. Random-access queries require loading the KMC database, which has different performance characteristics.

### Verdict for kmerdet

KMC3 is a viable alternative to jellyfish, especially for multi-k workflows (k>31) and memory-constrained environments. However, for the targeted panel use case (small data, high depth), jellyfish's in-memory hash table is simpler and fast enough. KMC3 becomes more attractive if kmerdet expands to WES or WGS applications.

## DSK

### How It Works

DSK (Rizk et al., 2013) uses a streaming disk-based approach to count k-mers with very low memory:

1. **Partition the k-mer space**: Divide all possible k-mers into P partitions based on a hash of the k-mer.
2. **Multi-pass reading**: For each partition, make a pass over the entire input, extracting only k-mers belonging to that partition.
3. **In-memory counting**: For each partition, load the extracted k-mers into a hash table, count, and output.
4. **Iterate**: Repeat for all P partitions.

The number of partitions P determines the memory-time tradeoff:
- More partitions = less memory per partition, but more passes over the input
- Fewer partitions = more memory needed, but fewer passes

### Key Properties

- **Constant memory**: DSK requires only `total_distinct_kmers / P` entries in memory at any time. For a human genome with ~3 billion distinct k-mers and P=100 partitions, each partition has ~30M k-mers, requiring ~500 MB.
- **Disk usage**: Temporary disk space for buffering k-mers during partitioning. For a human genome at 30x, DSK used ~160 GB of temporary disk.
- **Wall-clock time**: 17.9 hours for human genome 27-mers with 4 GB RAM (DSK paper). Much slower than jellyfish or KMC3.

### Handling of Low-Count K-mers

DSK counts all k-mers exactly, including singletons. Unlike jellyfish with `-L 2`, DSK does not discard low-count k-mers during counting (though it can filter during output). This is advantageous for low-VAF variant detection where true variant k-mers may have count=1.

However, retaining all k-mers significantly increases the output database size. For targeted panel data, where the vast majority of k-mers are high-count (target regions at 5000x), the singleton fraction is small, and this is not a concern.

### Verdict for kmerdet

DSK's extreme memory efficiency is unnecessary for targeted panel data (which has modest k-mer counts). Its multi-pass design makes it slower than jellyfish for small inputs. DSK would only be relevant if kmerdet were applied to WGS data on memory-constrained systems. For the current use case, jellyfish or KMC3 is preferred.

## Squeakr

### How It Works

Squeakr (Pandey et al., 2018) uses a counting quotient filter (CQF), a space-efficient probabilistic data structure that supports exact counting:

1. **Quotient filter basis**: Hash each k-mer to a fingerprint. Store the fingerprint in a quotient filter, which splits the fingerprint into a quotient (used as the table index) and a remainder (stored at that index).
2. **Counting extension**: Extend the quotient filter with per-element counters. Unlike Bloom filter-based approaches, the CQF supports exact counting of individual elements.
3. **Compact storage**: The CQF uses ~2-3 bits per element for the filter structure, plus counter bits.

### Performance

From Pandey et al. (2018):
- Squeakr takes less time for counting and random-point queries than KMC2
- Uses considerably less memory than KMC2 and Jellyfish2
- Supports exact counting (unlike count-min sketch)

### Advantages for Low-Count Regime

Squeakr's CQF structure is particularly efficient for skewed count distributions (many low-count, few high-count), which is exactly the distribution in sequencing data. The compact representation means:
- More distinct k-mers can be stored in the same memory
- Random queries are fast (single hash computation + short probe)
- Exact counts are maintained (no approximation error)

### Verdict for kmerdet

Squeakr offers an interesting middle ground between exact counting (jellyfish) and approximate counting (count-min sketch). Its memory efficiency could be valuable if kmerdet needs to hold multiple k-mer databases simultaneously (multi-k strategy). However, it lacks the ecosystem maturity and documentation of jellyfish and KMC3. Worth monitoring but not a priority for adoption.

## Streaming vs Batch Counting

### Current Batch Approach

The current pipeline uses batch counting:
```
FASTQ -> jellyfish count -> .jf file -> kmerdet detect (queries .jf)
```

All k-mers are counted before any detection begins. The entire count database exists as a file on disk, queryable at any time.

### Streaming Alternative

A streaming approach would count k-mers and detect variants simultaneously:
```
FASTQ -> kmerdet stream-detect -> results
```

In this model:
1. Read FASTQ incrementally
2. Update an in-memory k-mer count table as reads arrive
3. Periodically run detection on the current counts
4. Finalize detection when all reads have been processed

### Advantages of Streaming

1. **No intermediate file**: Eliminates the .jf file (200-500 MB for targeted panels)
2. **Lower latency**: Results begin emerging as soon as sufficient depth is reached
3. **Single binary**: No dependency on jellyfish installation
4. **Integrated UMI handling**: Can incorporate UMI-aware counting natively (see `error-correction-alternatives.md`)

### Disadvantages of Streaming

1. **No re-querying**: The batch .jf file can be queried multiple times (e.g., with different parameters, different targets). Streaming requires re-processing from FASTQ.
2. **Memory footprint**: The entire count table must be held in memory for the duration of processing. For targeted panels (~50M distinct k-mers at ~16 bytes each = ~800 MB), this is feasible. For WGS (~3B distinct k-mers = ~48 GB), this is prohibitive.
3. **Complexity**: Walking requires stable, final counts. If counts change as more reads arrive, the walking algorithm must handle incremental updates or wait for completion.
4. **Reproducibility**: Batch counting is deterministic (same FASTQ always produces same .jf). Streaming may introduce ordering dependencies if parallelized.

### Practical Streaming Design for kmerdet

A pragmatic streaming approach for targeted panels:

```rust
// Phase 1: Stream all reads, build count table in memory
let kmer_counts: DashMap<u64, AtomicU32> = DashMap::with_capacity(50_000_000);

fastq_reader.par_records().for_each(|record| {
    for kmer in extract_canonical_kmers(&record.seq, k) {
        kmer_counts.entry(kmer)
            .or_insert_with(|| AtomicU32::new(0))
            .fetch_add(1, Ordering::Relaxed);
    }
});

// Phase 2: Detect variants using the in-memory count table
let results = detect_variants(&targets, &kmer_counts, &config);
```

Phase 1 is a single pass over the FASTQ, building the full count table. Phase 2 uses the table for detection. This is technically "batch in memory" rather than true streaming, but it eliminates the .jf file and the jellyfish dependency.

### Verdict

For targeted panel data, the "batch in memory" approach is optimal: read all k-mers into memory once, then detect. This combines the simplicity of batch counting with the integration benefits of eliminating jellyfish. True streaming (detecting before all reads are counted) adds complexity with minimal benefit for a pipeline that completes in 6 minutes.

## Count-Min Sketch for Approximate Counting

### How Count-Min Sketch Works

The count-min sketch (Cormode & Muthukrishnan, 2005) is a probabilistic data structure for frequency estimation:

1. Allocate a 2D array of d rows and w columns, initialized to zero
2. Use d independent hash functions, one per row
3. To increment k-mer K: for each row i, increment array[i][hash_i(K) % w]
4. To query k-mer K: return min over all rows of array[i][hash_i(K) % w]

The structure uses O(d * w) space regardless of the number of distinct k-mers. The error is one-sided: counts are never underestimated, only overestimated (due to hash collisions).

### Error Bounds

For a sketch with w columns and d rows:
- Space: d * w counters
- Error probability: delta = (1/e)^d (probability that error exceeds epsilon * N)
- Error magnitude: epsilon = e / w (maximum overcount as fraction of total count N)

For typical parameters (d=4, w=2^20 ~= 1M):
- Space: 4 * 1M * 4 bytes = 16 MB
- Error probability: (1/e)^4 ~= 1.8%
- Error magnitude per query: e / 2^20 * N_total

For a targeted panel with N_total = 5 billion total k-mers (5000x depth, 50 targets, 200 bp each, 120 k-mers per read):
- Expected overcount: e / 2^20 * 5 billion = 2.72 / 1048576 * 5e9 = ~12,974

An expected overcount of ~13,000 is catastrophic for variant detection where true variant k-mers have counts of 5-50. The count-min sketch would make every k-mer appear to have count >= 13,000, obliterating the signal.

### Can We Fix This?

**Conservative update**: Instead of incrementing all d hash positions, only increment positions where the current value equals the minimum. This significantly reduces overcount but does not eliminate it.

**Count-min-log**: Use logarithmic counters to reduce the impact of overcounting. Still fundamentally unsuitable for the low-count regime.

**Larger sketch**: Increase w to reduce epsilon. To get epsilon * N_total < 1 (overcount less than 1):
- w = e * N_total = 2.72 * 5e9 = 13.6 billion columns
- Space: 4 * 13.6B * 4 bytes = 217 GB
- This defeats the purpose of approximate counting.

### Verdict for kmerdet

**Count-min sketch is not suitable for low-VAF variant detection.** The fundamental problem is that overcounting is proportional to total k-mer count (N_total), which is enormous for high-depth targeted sequencing. True variant signals at 0.1% VAF produce counts of ~5, while expected overcount is ~13,000. No amount of tuning makes this work.

Count-min sketch could potentially serve as a **pre-filter**: quickly identify k-mers with count > threshold (using the sketch's property that counts are never underestimated). All k-mers that pass the sketch filter would then be counted exactly. However, for targeted panel data, exact counting is already fast and memory-efficient enough that a pre-filter adds complexity without meaningful benefit.

### khmer Software

The khmer package (Brown et al., 2015) used count-min sketches for k-mer counting in sequencing analysis. However, khmer was designed for WGS at moderate depth where the count distribution is broad and approximate counts are acceptable for applications like coverage normalization. It explicitly acknowledges that the data structure "can also blow up the memory usage for skewed data distributions, as often occur with k-mers in sequencing datasets."

## Counting K-mers Directly in Rust from FASTQ

### Architecture

For targeted panel data, the most compelling alternative to jellyfish is counting k-mers directly in Rust within the kmerdet binary:

```rust
use dashmap::DashMap;
use rayon::prelude::*;
use needletail::parse_fastx_file;

pub struct RustKmerCounter {
    counts: DashMap<u64, u32>,
    k: usize,
}

impl RustKmerCounter {
    pub fn new(k: usize, estimated_kmers: usize) -> Self {
        Self {
            counts: DashMap::with_capacity(estimated_kmers),
            k,
        }
    }

    pub fn count_fastq(&self, path: &Path) -> Result<()> {
        let reader = parse_fastx_file(path)?;

        // Process records in parallel using rayon
        reader.par_bridge().for_each(|record| {
            let record = record.expect("invalid FASTQ record");
            let seq = record.normalize(false);

            // Extract canonical k-mers using 2-bit encoding
            for kmer_u64 in canonical_kmers(&seq, self.k) {
                self.counts
                    .entry(kmer_u64)
                    .and_modify(|c| *c += 1)
                    .or_insert(1);
            }
        });

        Ok(())
    }

    pub fn query(&self, kmer: u64) -> u32 {
        self.counts.get(&kmer).map(|v| *v).unwrap_or(0)
    }
}
```

### 2-Bit Encoding and Canonical Form

Each nucleotide is encoded in 2 bits (A=00, C=01, G=10, T=11), allowing a 31-mer to fit in a single u64 (62 bits):

```rust
fn encode_base(b: u8) -> u64 {
    match b {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => panic!("invalid base"),
    }
}

fn canonical_kmers(seq: &[u8], k: usize) -> impl Iterator<Item = u64> + '_ {
    let mask: u64 = (1u64 << (2 * k)) - 1;
    let mut forward: u64 = 0;
    let mut reverse: u64 = 0;

    // Initialize first k-1 bases
    for i in 0..k-1 {
        let base = encode_base(seq[i]);
        forward = (forward << 2) | base;
        reverse = reverse | ((3 - base) << (2 * i));
    }

    (k-1..seq.len()).map(move |i| {
        let base = encode_base(seq[i]);
        forward = ((forward << 2) | base) & mask;
        reverse = (reverse >> 2) | ((3 - base) << (2 * (k - 1)));
        std::cmp::min(forward, reverse)  // canonical form
    })
}
```

### Memory Requirements

For a targeted panel with ~50M distinct k-mers:

| Component | Per entry | Total |
|-----------|-----------|-------|
| DashMap overhead | ~64 bytes (bucket + metadata) | ~3.2 GB |
| K-mer key (u64) | 8 bytes | 400 MB |
| Count value (u32) | 4 bytes | 200 MB |
| **Total** | | **~3.8 GB** |

This is higher than jellyfish (~500 MB for the same data) because DashMap has significant per-entry overhead for its concurrent access support (sharded locks, hash metadata). Alternatives:

**Custom hash table**: A purpose-built open-addressing hash table with 2-bit k-mer keys and 32-bit counters:
```rust
struct CompactKmerTable {
    keys: Vec<u64>,      // 2-bit encoded k-mers (0 = empty)
    counts: Vec<u32>,    // Parallel array of counts
    mask: usize,         // Table size - 1 (power of 2)
}
```

Memory: ~12 bytes per slot. At 80% load with 50M k-mers: 62.5M slots * 12 bytes = ~750 MB. Much closer to jellyfish.

**CHTKC approach**: The CHTKC tool (Wang et al., 2021) demonstrated a lock-free chaining hash table for k-mer counting that achieves competitive performance with jellyfish. Its Rust equivalent would use atomics for lock-free concurrent access:

```rust
use std::sync::atomic::{AtomicU64, AtomicU32, Ordering};

struct LockFreeKmerTable {
    keys: Vec<AtomicU64>,
    counts: Vec<AtomicU32>,
    size: usize,
}

impl LockFreeKmerTable {
    fn insert_or_increment(&self, kmer: u64) {
        let mut pos = hash(kmer) & (self.size - 1);
        loop {
            let current = self.keys[pos].load(Ordering::Relaxed);
            if current == 0 {
                // Empty slot: try to claim it
                match self.keys[pos].compare_exchange(
                    0, kmer, Ordering::AcqRel, Ordering::Relaxed
                ) {
                    Ok(_) => {
                        self.counts[pos].store(1, Ordering::Release);
                        return;
                    }
                    Err(_) => continue, // Another thread claimed it; retry
                }
            } else if current == kmer {
                // Our k-mer: increment count
                self.counts[pos].fetch_add(1, Ordering::Relaxed);
                return;
            } else {
                // Collision: linear probe
                pos = (pos + 1) & (self.size - 1);
            }
        }
    }
}
```

### Advantages of Native Rust Counting

1. **No external dependency**: Eliminates the jellyfish binary and .jf file format dependency
2. **Integrated pipeline**: Count and detect in a single binary invocation
3. **UMI-aware counting**: Extend the counter to track UMI families per k-mer (see `error-correction-alternatives.md`)
4. **Strand-aware counting**: Maintain forward/reverse counts per canonical k-mer
5. **Streaming**: No intermediate .jf file; count directly into memory
6. **Custom optimizations**: Tune the hash table for the specific access patterns of k-mer walking (e.g., prefetching for extension queries)

### Feasibility for Targeted Panels

For targeted panel data (the primary use case), native Rust counting is highly feasible:

| Metric | Jellyfish | Native Rust (custom hash) |
|--------|-----------|---------------------------|
| Distinct k-mers | ~50M | ~50M |
| Memory | ~500 MB | ~750 MB |
| Counting time (estimate) | ~1.5 min | ~1-2 min |
| Query time per k-mer | O(1) amortized | O(1) amortized |
| Dependencies | jellyfish binary + FFI | None (pure Rust) |
| UMI-aware | No | Yes (extensible) |
| Strand-aware | Only with -C workaround | Native |

The slight memory increase (~50%) is negligible on modern workstations. The counting speed should be comparable, as the bottleneck is FASTQ I/O and hashing, both of which Rust handles efficiently.

### Disadvantages

1. **Implementation effort**: Building a correct, performant concurrent hash table is non-trivial. Edge cases in lock-free programming are subtle.
2. **Testing burden**: Must verify that counts match jellyfish exactly for validation.
3. **Not suitable for WGS**: 3 billion distinct k-mers at 12 bytes each = 36 GB. This requires disk-based partitioning (KMC3-style), which is a significantly larger implementation effort.

### Recommended Implementation Path

1. **Phase 1 (current)**: Use jellyfish via FFI for .jf reading. This works and is validated.
2. **Phase 2**: Implement a pure Rust .jf file reader (parse jellyfish format without the C++ library). Eliminates compile-time dependency on libjellyfish.
3. **Phase 3**: Implement native Rust k-mer counting for targeted panels. Use DashMap initially (simpler, correct), then optimize with a custom lock-free table if profiling shows DashMap overhead is significant.
4. **Phase 4**: Extend native counting with UMI-aware and strand-aware features.

## Comparison Summary

| Tool/Approach | Memory | Speed | Accuracy | k range | UMI support | kmerdet fit |
|---------------|--------|-------|----------|---------|-------------|-------------|
| Jellyfish 2 | Proportional to distinct k-mers | Fast (lock-free) | Exact | k<=31 natively | No | Current choice |
| KMC3 | Configurable (disk-backed) | Fastest for large data | Exact | Arbitrary | No | Good for WGS/multi-k |
| DSK | Very low (multi-pass) | Slow | Exact | Arbitrary | No | Low-memory environments only |
| Squeakr | Low (quotient filter) | Fast | Exact | Moderate | No | Interesting but immature |
| Count-min sketch | Very low | Very fast | Approximate (overcounts) | Any | No | **Not suitable** for low-VAF |
| Native Rust (DashMap) | ~3.8 GB for 50M k-mers | Fast | Exact | k<=31 (u64) | Extensible | Good for targeted panels |
| Native Rust (custom) | ~750 MB for 50M k-mers | Fast | Exact | k<=31 (u64) | Extensible | **Best long-term option** |

## References

- Kokot, M., Dlugosz, M. & Deorowicz, S. (2017). KMC 3: counting and manipulating k-mer statistics. *Bioinformatics*, 33(17), 2759-2761.
- Rizk, G., Lavenier, D. & Chikhi, R. (2013). DSK: k-mer counting with very low memory usage. *Bioinformatics*, 29(5), 652-653.
- Pandey, P. et al. (2018). Squeakr: an exact and approximate k-mer counting system. *Bioinformatics*, 34(4), 568-575.
- Cormode, G. & Muthukrishnan, S. (2005). An improved data stream summary: the count-min sketch and its applications. *Journal of Algorithms*, 55(1), 58-75.
- Manekar, S.C. & Sathe, S.R. (2018). A benchmark study of k-mer counting methods for high-throughput sequencing. *GigaScience*, 7(12), giy125.
- Wang, Y. et al. (2021). CHTKC: a robust and efficient k-mer counting algorithm based on a lock-free chaining hash table. *Briefings in Bioinformatics*, 22(3), bbaa063.
- Brown, C.T. et al. (2015). khmer: Working with Big Data in Bioinformatics. https://github.com/dib-lab/khmer
- Marcais, G. & Kingsford, C. (2011). A fast, lock-free approach for efficient parallel counting of occurrences of k-mers. *Bioinformatics*, 27(6), 764-770.
- DashMap: Blazing fast concurrent HashMap for Rust. https://github.com/xacrimon/dashmap
- Bauer, S. et al. (2024). Hyper-k-mers: efficient streaming k-mers representation. *bioRxiv*, 2024.11.06.620789.
