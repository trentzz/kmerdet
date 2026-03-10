# Rust Implementation Considerations

## Why Rust?

1. **Performance**: The km algorithm involves millions of k-mer lookups, graph construction, and path finding — all CPU-bound operations that benefit from Rust's zero-cost abstractions
2. **Parallelism**: Rust's ownership model makes concurrent k-mer walking safe and efficient (rayon, crossbeam)
3. **Memory efficiency**: 2-bit k-mer encoding, compact graph representation
4. **Single binary**: No Python runtime, no jellyfish Python bindings — simpler deployment
5. **Bioinformatics ecosystem**: rust-bio, noodles, needletail for sequence parsing

## Core Components to Implement

### 1. Jellyfish Database Reader

**Options**:
- **FFI bindings to libjellyfish**: Use `bindgen` to create Rust bindings to jellyfish's C++ API
  - Pro: Exact compatibility with .jf format
  - Con: C++ dependency, complex build
- **Pure Rust reader**: Parse the .jf binary format directly
  - Pro: No external dependencies, full control
  - Con: Must reverse-engineer the binary format
- **Shell out to jellyfish query**: Use `std::process::Command`
  - Pro: Simple
  - Con: Very slow for millions of queries (process spawn overhead)

**Recommendation**: FFI bindings to libjellyfish for compatibility, with a pure-Rust reader as a stretch goal.

### 2. K-mer Representation

```rust
/// 2-bit encoded k-mer (up to 32 bases in a u64)
struct Kmer {
    data: u64,
    k: u8,
}

impl Kmer {
    fn canonicalize(&self) -> Self { /* min of self and reverse complement */ }
    fn extend_right(&self, base: u8) -> Self { /* shift left, add base */ }
    fn extend_left(&self, base: u8) -> Self { /* shift right, add base */ }
    fn to_string(&self) -> String { /* decode to ACGT string */ }
}
```

### 3. K-mer Graph

```rust
struct KmerGraph {
    nodes: Vec<KmerNode>,           // All discovered k-mers
    edges: HashMap<usize, Vec<(usize, f32)>>,  // Adjacency list with weights
    node_index: HashMap<u64, usize>,  // Kmer hash → node index
}

struct KmerNode {
    kmer: Kmer,
    count: u32,        // From jellyfish DB
    is_reference: bool,
}
```

### 4. Path Finding

Implement Dijkstra's algorithm with:
- Forward pass from source
- Backward pass from sink
- Reference edge removal
- All-shortest-paths enumeration through non-reference edges

### 5. Path Quantification

Linear algebra for VAF estimation:
- Use `ndarray` + `ndarray-linalg` for least squares
- Gradient descent refinement for non-negative coefficients

### 6. Variant Classification

Port the `diff_path_without_overlap` logic:
- Left-to-right scan for divergence point
- Right-to-left scan for convergence point
- Classification based on path length differences

## Potential Rust Crates

| Crate | Purpose |
|-------|---------|
| `clap` | CLI argument parsing |
| `rayon` | Data parallelism (parallel target processing) |
| `ndarray` | Numerical arrays for path quantification |
| `ndarray-linalg` | Linear algebra (least squares) |
| `petgraph` | Graph data structures and algorithms |
| `needletail` | Fast FASTA/FASTQ parsing |
| `noodles` | Bioinformatics file format support |
| `rust-bio` | Biological sequence utilities |
| `serde` / `serde_json` | JSON parsing (jellyfish header) |
| `dashmap` | Concurrent hash map for parallel k-mer collection |
| `crossbeam` | Concurrent utilities |
| `anyhow` / `thiserror` | Error handling |
| `tracing` | Logging and diagnostics |

## Architecture Sketch

```
kmerdet/
├── src/
│   ├── main.rs              # CLI entry point
│   ├── lib.rs               # Library root
│   ├── cli/
│   │   ├── mod.rs
│   │   ├── find_mutation.rs  # find_mutation subcommand
│   │   ├── find_report.rs   # find_report subcommand
│   │   └── min_cov.rs       # min_cov subcommand
│   ├── jellyfish/
│   │   ├── mod.rs
│   │   ├── reader.rs         # .jf file reader
│   │   └── query.rs          # K-mer query interface
│   ├── kmer/
│   │   ├── mod.rs
│   │   ├── encoding.rs       # 2-bit k-mer encoding
│   │   ├── walker.rs         # K-mer extension/walking
│   │   └── canonical.rs      # Canonical form handling
│   ├── graph/
│   │   ├── mod.rs
│   │   ├── builder.rs        # Graph construction from k-mers
│   │   └── pathfind.rs       # Dijkstra + shortest paths
│   ├── variant/
│   │   ├── mod.rs
│   │   ├── classifier.rs     # Variant type classification
│   │   └── quantifier.rs     # Path quantification (linear regression)
│   ├── sequence/
│   │   ├── mod.rs
│   │   ├── target.rs         # Target sequence handling
│   │   └── reference.rs      # Reference path representation
│   └── output/
│       ├── mod.rs
│       ├── tsv.rs            # TSV output format
│       └── vcf.rs            # VCF output format (stretch goal)
├── tests/
│   ├── integration/
│   └── data/                 # Test jellyfish DBs and targets
├── benches/                  # Benchmarks
├── Cargo.toml
└── README.md
```

## Performance Opportunities

1. **Parallel target processing**: Process each target sequence independently via rayon
2. **Batch k-mer queries**: Read jellyfish DB into memory-mapped structure for O(1) lookups
3. **Compact graph representation**: Use CSR (compressed sparse row) format instead of HashMap
4. **SIMD k-mer operations**: Use SIMD for k-mer encoding/comparison
5. **Streaming output**: Write results as they're computed, not after all targets complete

## Compatibility Considerations

1. **Output format**: Must match km's TSV output format for downstream tool compatibility
2. **Jellyfish version**: Target jellyfish 2.2.x .jf format
3. **Canonical k-mer handling**: Must exactly match jellyfish's canonicalization
4. **Default parameters**: Match km's defaults for equivalent results
5. **Target FASTA format**: Same input format as km

## Testing Strategy

1. **Unit tests**: K-mer encoding, canonical form, graph operations
2. **Integration tests**: Compare output against km's output on same inputs
3. **Benchmark**: Performance comparison against km on real datasets
4. **Test data**: Use km's example .jf files and target catalogs
