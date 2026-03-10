# Feature: UMI-Aware K-mer Counting (Stretch Goal)

## What It Does

UMI-aware k-mer counting replaces standard read-level k-mer counting with
molecule-level counting: for each k-mer, the count represents the number of
distinct UMI families (original DNA molecules) containing that k-mer, rather
than the total number of sequencing reads. This eliminates PCR amplification
bias at its source, produces more accurate rVAF estimates, and may eliminate
the need for external UMI deduplication (HUMID) entirely.

```
Standard counting (jellyfish):
    kmer_count[K] = total reads containing k-mer K

UMI-aware counting:
    kmer_count[K] = distinct UMI families containing k-mer K
```

## Why It Matters

### PCR Amplification Bias

PCR amplification is inherently stochastic: some molecules are amplified 50x
while others are amplified 5x. This creates a systematic distortion in k-mer
counts that directly biases rVAF estimates.

Consider a locus with 500 unique molecules (after deduplication), 499 carrying
the reference allele and 1 carrying a variant (true VAF = 0.2%). If the variant
molecule is amplified 30x and reference molecules are amplified an average of
10x:

| Counting method | Variant count | Ref count | Observed VAF |
|-----------------|---------------|-----------|--------------|
| Read-level (raw) | 30 | 4990 | 0.60% |
| Read-level (post-HUMID dedup) | 1 | 499 | 0.20% |
| Molecule-level (UMI-aware) | 1 | 499 | 0.20% |

HUMID dedup and UMI-aware counting give the same result in this case. But
consider a PCR error that occurs in the first amplification cycle of a reference
molecule:

| Counting method | Error k-mer count | True signal? |
|-----------------|-------------------|--------------|
| Read-level (raw) | 20 (amplified from early-cycle error) | No |
| Read-level (post-HUMID dedup) | 1 (dedup to single representative) | Ambiguous |
| Molecule-level (UMI-aware) | 1 (one UMI family) | Clearly single-molecule |

With read-level counting, the early-cycle PCR error produces a count of 20,
well above typical thresholds. After HUMID dedup, the count drops to 1, but
HUMID outputs a single representative read that may or may not carry the error.
With UMI-aware counting, the error is unambiguously identified as a
single-molecule event, and a dual-threshold model (require both read_count >= 2
AND molecule_count >= 2) would correctly reject it.

### Eliminating the HUMID Bottleneck

HUMID consumes approximately 2 minutes of the 6-minute pipeline (33% of total
runtime). It is single-threaded, writes ~10-15 GB of intermediate FASTQ to disk,
and loses UMI group membership information at the output boundary.

If kmerdet counts k-mers directly from raw FASTQ with UMI awareness, HUMID
becomes unnecessary:

```
Current pipeline:
    FASTQ -> HUMID dedup -> deduped FASTQ -> jellyfish count -> .jf -> kmerdet

Proposed pipeline:
    FASTQ -> kmerdet count-and-detect -> results
```

This eliminates: the HUMID dependency, the intermediate FASTQ files (10-15 GB
disk), the jellyfish dependency, and the .jf file (200-500 MB disk). The
pipeline becomes a single Rust binary reading raw FASTQ.

### Lowering the VAF Detection Threshold

The current practical detection threshold is ~0.1% VAF, limited partly by the
inability to distinguish true low-VAF variants from PCR artifacts at the k-mer
level. Molecule-level counting provides a natural noise model:

- A k-mer supported by 3 distinct UMI families is far more credible than one
  supported by 30 reads from a single family
- A dual-threshold model (molecule_count >= M AND read_count >= R) can be more
  permissive on read counts when molecule evidence is strong
- At 0.05% VAF with 5000 unique molecules, expected variant molecules = 2.5,
  which is detectable with molecule_count >= 2

This could push the detection threshold toward 0.05%, closer to competing
technologies like duplex sequencing with error correction (0.01-0.02% VAF
claimed).

## Research Backing

### Thesis Limitations (room-for-improvement.md)

The thesis identified "No UMI Integration into K-mer Counting" as Limitation #5
(Severity: High). The thesis noted that UMI information is discarded after
HUMID deduplication and is not available during k-mer counting or walking. This
was identified as a "significant missed opportunity for improved artifact
detection."

The cross-limitation analysis shows that solving UMI-aware counting (Limitation
#5) also solves the HUMID bottleneck (Limitation #7) and partially addresses the
VAF detection threshold (Limitation #3).

### Error Correction Alternatives (error-correction-alternatives.md)

The document evaluates standard k-mer error correction tools (Lighter, BFC,
Musket) and finds them unsuitable for low-VAF variant detection because they
risk correcting true variant k-mers to reference. The recommended strategy is
molecule-level counting rather than error correction.

The document also describes k-mer-level duplex evidence as an alternative to
read-level duplex consensus: for each k-mer, track forward-strand and
reverse-strand UMI family counts. A k-mer with support from both strands has
duplex-level evidence without requiring explicit duplex consensus calling.

### HUMID Analysis (humid-analysis.md)

The analysis documents HUMID's limitations:

- **No duplex consensus**: HUMID does not pair alpha/beta UMIs for cross-strand
  validation. Single-strand errors survive deduplication.
- **No consensus calling**: HUMID selects one representative read per family.
  If that read happens to contain an error, the error persists.
- **UMI collision**: At 5000x depth with 8-bp UMIs (65,536 possible), collisions
  are near-certain. Multiple distinct molecules sharing a UMI are incorrectly
  merged, losing ~17% of variant molecules at the margins.
- **Single-threaded**: Uses one CPU core while the rest of the pipeline is
  multi-threaded.

The analysis also surveys alternatives: fgbio (alignment-required duplex
consensus), UMI-tools (single-threaded, alignment-dependent), UMICollapse
(fast clustering but no consensus), RUMINA (Rust-based, up to 10x faster,
could serve as a basis for kmerdet's built-in dedup).

### Counting Alternatives (counting-alternatives.md)

The document evaluates native Rust k-mer counting as the long-term alternative
to jellyfish. A custom lock-free hash table in Rust achieves ~750 MB for 50M
distinct k-mers (comparable to jellyfish) while being extensible to UMI-aware
and strand-aware counting. The DashMap-based approach is simpler but uses more
memory (~3.8 GB for 50M k-mers).

## Design Considerations

### UMI Extraction

UMIs are extracted from FASTQ read headers. Standard UMI formats:

```
@READ_NAME:UMI_SEQUENCE 1:N:0:AGATCTCG
@READ_NAME:ACGTACGT+TGCATGCA 1:N:0:AGATCTCG   (dual-index UMI)
```

kmerdet supports configurable UMI parsing via regex or fixed-position extraction:

```toml
[umi]
source = "header"          # "header" or "inline"
pattern = ":([ACGTN]+)$"  # regex to extract UMI from header
# OR
position = "0:8"           # first 8 bases of read are UMI (inline)
```

### Data Structure: Hybrid Exact/Approximate Counting

Three implementation options were evaluated:

**Option A: Exact counting with HashMap<Kmer, HashSet<UMI>>**

Each k-mer stores the full set of UMI hashes. Memory: ~152 bytes per k-mer
(8 bytes key + 64 bytes HashSet overhead + 10 * 8 bytes for average 10 UMIs).
For 45M distinct k-mers: ~6.8 GB. Exact but memory-heavy.

**Option B: Approximate counting with HashMap<Kmer, HyperLogLog>**

HyperLogLog (Flajolet et al., 2007) estimates distinct element count using
~64 bytes per counter with ~2% error. For 45M k-mers: ~3.2 GB. Memory-efficient
but introduces estimation error at low counts (1-5 molecules), which is exactly
the regime where precision matters most.

**Option C: Hybrid approach (recommended)**

Use exact counting for k-mers with few UMI observations and HyperLogLog for
high-count k-mers:

```rust
enum UmiCounter {
    /// Up to 8 UMIs stored exactly in a small inline array
    Exact(SmallVec<[u64; 8]>),
    /// Above threshold, switch to approximate counting
    Approximate(HyperLogLog),
}

impl UmiCounter {
    fn add(&mut self, umi_hash: u64) {
        match self {
            UmiCounter::Exact(vec) => {
                if !vec.contains(&umi_hash) {
                    if vec.len() < 8 {
                        vec.push(umi_hash);
                    } else {
                        // Transition to HLL
                        let mut hll = HyperLogLog::new(6);
                        for existing in vec.iter() {
                            hll.insert(existing);
                        }
                        hll.insert(&umi_hash);
                        *self = UmiCounter::Approximate(hll);
                    }
                }
            }
            UmiCounter::Approximate(hll) => {
                hll.insert(&umi_hash);
            }
        }
    }

    fn count(&self) -> u32 {
        match self {
            UmiCounter::Exact(vec) => vec.len() as u32,
            UmiCounter::Approximate(hll) => hll.len() as u32,
        }
    }
}
```

Memory analysis:
- Variant k-mers (low count, ~1% of k-mers): SmallVec with 1-5 entries =
  8 + 5*8 = 48 bytes. Exact counts where they matter most.
- Reference k-mers (high count, ~99%): HyperLogLog = 64 bytes.
  Approximate counts where precision is less critical.
- Total for 45M k-mers: ~2.9 GB (weighted average of exact and approximate).

### Strand-Aware UMI Counting

Extend the molecule counter to track forward and reverse strand UMI families
separately:

```rust
struct StrandAwareUmiCounter {
    forward: UmiCounter,
    reverse: UmiCounter,
}

impl StrandAwareUmiCounter {
    fn molecule_count(&self) -> u32 {
        self.forward.count() + self.reverse.count()
    }

    fn duplex_count(&self) -> u32 {
        // Minimum of forward and reverse = duplex-supported molecules
        std::cmp::min(self.forward.count(), self.reverse.count())
    }

    fn strand_bias(&self) -> f64 {
        let f = self.forward.count() as f64;
        let r = self.reverse.count() as f64;
        (f - r).abs() / (f + r + 1e-10)
    }
}
```

This provides duplex-like evidence at the k-mer level without requiring
explicit duplex consensus calling. A k-mer with `duplex_count >= 2` has been
observed on both strands of at least 2 molecules, providing strong evidence
against artifact.

### Dual-Threshold Detection Model

With both read counts and molecule counts available, the extension threshold
becomes two-dimensional:

```
pass_extension = (read_count >= R) AND (molecule_count >= M)
```

Recommended defaults:
- R = adaptive read count threshold (from coverage estimation)
- M = 2 (require at least 2 distinct molecules)

The molecule count threshold M = 2 is the key innovation: a k-mer observed in
2 or more independent UMI families is extremely unlikely to be a PCR artifact
(which would require the same error occurring independently in 2 molecules).
Sequencing errors can produce count-2 k-mers from a single molecule (if 2 reads
from the same UMI family both contain the same error), but molecule_count = 2
requires true independent observations.

### Integration with Existing Pipeline

UMI-aware counting replaces both HUMID and jellyfish. The walking algorithm,
graph construction, pathfinding, variant classification, and NNLS quantification
are unchanged -- they operate on k-mer counts regardless of how those counts
were obtained.

The KmerDatabase trait in `src/jellyfish/mod.rs` abstracts the count source:

```rust
pub trait KmerDatabase: Send + Sync {
    fn query(&self, kmer: &Kmer) -> u64;
    fn k(&self) -> usize;
}
```

A UMI-aware counter implements this trait, returning molecule counts instead of
read counts. No changes to the walking or graph code are needed.

### Processing Architecture

```rust
// Single-pass processing from raw FASTQ
fn count_with_umis(fastq_paths: &[PathBuf], k: usize, umi_config: &UmiConfig)
    -> UmiAwareKmerDatabase
{
    let db = UmiAwareKmerDatabase::new(k);

    for path in fastq_paths {
        let reader = parse_fastx_file(path).unwrap();
        reader.par_bridge().for_each(|record| {
            let record = record.unwrap();
            let umi = extract_umi(&record, umi_config);
            let umi_hash = hash_umi(&umi);
            let strand = determine_strand(&record);

            for kmer in canonical_kmers(&record.seq(), k) {
                db.insert(kmer, umi_hash, strand);
            }
        });
    }

    db
}
```

The `par_bridge()` call enables multi-threaded processing via rayon, utilizing
all available CPU cores. This alone addresses the single-threaded HUMID
bottleneck.

### Memory Optimization

For targeted panel data with ~45M distinct k-mers, the hybrid counter uses
~2.9 GB. This is higher than jellyfish (~500 MB) but acceptable on modern
workstations. Optimization strategies if memory is constrained:

1. **Two-pass approach**: First pass counts total occurrences (8 bytes/k-mer,
   ~360 MB). Second pass adds UMI tracking only for k-mers with count > 1
   (most error k-mers are singletons and do not need UMI tracking).

2. **Region-restricted counting**: Only count k-mers that appear in or near
   target sequences. For a 50-target panel, this reduces the k-mer space from
   ~45M to ~500K relevant k-mers, with trivial memory requirements.

3. **Streaming with disk spill**: For very large datasets (WGS), partition
   k-mers by hash prefix and process one partition at a time (KMC3-style).

### Benefits Summary

| Benefit | Impact | Replaces |
|---------|--------|----------|
| Eliminate PCR bias | More accurate rVAF | HUMID dedup (partial) |
| Lower VAF threshold | Detect 0.05% VAF | Conservative thresholds |
| Strand-aware counting | Artifact detection | N/A (not available today) |
| Single binary | No external deps | HUMID + jellyfish |
| Multi-threaded | ~2x faster preprocess | Single-threaded HUMID |
| No intermediate files | Save 10-15 GB disk | deduped FASTQ + .jf file |

## Acceptance Criteria

### Counting Correctness

1. For a k-mer present in N distinct UMI families, `molecule_count(K) == N`
   (exact mode) or `|molecule_count(K) - N| <= 0.05 * N` (approximate mode).

2. Molecule counts differ from read counts by >20% in regions with high PCR
   amplification (expected: molecules << reads in amplified regions).

3. The UMI extraction correctly parses standard UMI header formats
   (Illumina, 10x, custom regex).

### Functional Equivalence

4. When all UMI families have exactly 1 read (no PCR duplicates), molecule
   counts equal read counts. UMI-aware counting degenerates to standard
   counting.

5. The walking algorithm produces identical paths when given molecule counts
   as when given deduplicated read counts (same k-mer count, different source).

### Strand Awareness

6. Forward and reverse strand molecule counts are tracked separately.

7. The strand bias metric correctly identifies known single-strand artifacts
   (e.g., oxidative damage G>T enriched on one strand).

8. Duplex molecule count (min of forward and reverse) provides stronger
   evidence than simplex count for variant calls.

### Performance

9. UMI-aware counting from raw FASTQ completes in less time than
   HUMID + jellyfish combined (target: <3 minutes for a typical panel sample).

10. Memory usage for a 50-target panel does not exceed 4 GB.

### Integration

11. The UMI-aware database implements the KmerDatabase trait and is
    transparent to the walking/graph/classification code.

12. The `run` subcommand supports a `--umi` flag that enables UMI-aware
    counting from raw FASTQ, bypassing the need for pre-existing .jf databases.

### Validation

13. On a sample with known variants at 0.1% VAF, UMI-aware counting achieves
    equal or better sensitivity compared to HUMID + jellyfish counting.

14. On a sample with known PCR artifacts, UMI-aware counting reduces false
    positive rate compared to HUMID + jellyfish counting.
