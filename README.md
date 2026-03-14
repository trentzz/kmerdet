# kmerdet

K-mer based variant detection for targeted liquid biopsy ctDNA monitoring.

kmerdet is a Rust reimplementation of the [km](https://github.com/iric-soft/km) + kmtools + kam ecosystem, designed for detecting somatic variants in UMI-deduplicated sequencing data from liquid biopsy panels. It functions as a drop-in replacement for `km find_mutation` and `kmtools chunk`, while adding configurable sensitivity, statistical confidence scoring, and multi-k detection.

## Key Features

- **K-mer graph walking** — DFS-based exploration of variant branches with adaptive thresholds
- **NNLS quantification** — Non-negative least squares estimation of variant allele frequencies (rVAF)
- **Statistical confidence** — Binomial p-values, Fisher's method, Phred-scaled QUAL scores, bootstrap CIs
- **Configurable sensitivity** — Preset profiles (ultra/high/standard/strict) or fine-grained numeric control via config file
- **Multi-k detection** — Run at multiple k-mer lengths and merge with consensus voting for improved sensitivity across variant types
- **Graph pruning** — Dead-end removal, tip clipping, bubble collapsing to reduce false positives
- **Bidirectional walking** — Walk both forward and backward from reference k-mers
- **Multiple output formats** — TSV, CSV, VCF 4.3, JSON, JSONL, Excel
- **Tumor-informed filtering** — Reference-mode and alt-sequence matching against expected variants
- **Panel-of-normals** — Build and filter against a PoN to suppress recurrent artifacts
- **Ground truth benchmarking** — Built-in `benchmark` subcommand with sensitivity/specificity/F1 metrics

## Installation

```bash
git clone https://github.com/<your-org>/kmerdet.git
cd kmerdet
cargo build --release
```

The binary will be at `target/release/kmerdet`.

### Dependencies

- **Rust** >= 1.70
- **jellyfish** >= 2.3 (for building k-mer databases from FASTQ/BAM)

kmerdet reads jellyfish `.jf` databases using a pure Rust reader — no C/C++ jellyfish library is needed at runtime.

## Quick Start

```bash
# 1. Build a jellyfish database from your sequencing data
jellyfish count -m 31 -s 100M -t 8 -C reads.fastq -o sample.jf

# 2. Run variant detection
kmerdet detect -d sample.jf -T targets/ --format tsv -o results.tsv

# 3. Filter against expected variants
kmerdet filter -i results.tsv --targets expected_variants.tsv -o filtered.tsv

# 4. Or run the full pipeline in one command
kmerdet run -d sample.jf -T targets/ --expected expected_variants.tsv -o final.tsv
```

## Subcommands

| Command | Description |
|---------|-------------|
| `detect` | Run variant detection on targets against a jellyfish database |
| `filter` | Apply tumor-informed filtering to detection results |
| `merge` | Merge results from multiple detection runs |
| `stats` | Generate summary statistics from results |
| `plot` | Generate visualization charts |
| `coverage` | Report k-mer coverage statistics for a database |
| `run` | Full pipeline: detect -> merge -> filter |
| `benchmark` | Compare detection results against ground truth |
| `pon` | Build or filter against a panel-of-normals |

## Sensitivity Configuration

kmerdet provides a configurable sensitivity dial for tuning the sensitivity/specificity tradeoff. This can be set via CLI flags, named presets, or a TOML config file.

### Presets

```bash
kmerdet detect -d sample.jf -T targets/ --sensitivity ultra    # Maximum sensitivity
kmerdet detect -d sample.jf -T targets/ --sensitivity high     # High sensitivity
kmerdet detect -d sample.jf -T targets/ --sensitivity standard # Balanced (default)
kmerdet detect -d sample.jf -T targets/ --sensitivity strict   # High specificity
```

| Preset | count | ratio | min_qual | min_rvaf | Use case |
|--------|-------|-------|----------|----------|----------|
| ultra | 1 | 1e-6 | 0 | 0 | MRD monitoring, find everything |
| high | 1 | 1e-4 | 5 | 0.00001 | Low-VAF detection |
| standard | 2 | 0.05 | 10 | 0.0001 | General use (default) |
| strict | 3 | 0.05 | 20 | 0.001 | Treatment selection, high confidence |

### Fine-Grained Control

Individual parameters can be set on the CLI or in a config file. CLI flags always take precedence.

```bash
# CLI flags
kmerdet detect -d sample.jf -T targets/ \
    --count 1 --ratio 0.0001 --min-qual 8 --min-rvaf 0.00005

# Config file
kmerdet --config my_settings.toml detect -d sample.jf -T targets/
```

Example config file (`my_settings.toml`):

```toml
[sensitivity]
preset = "high"          # start from "high" defaults
min_qual = 8.0           # override just the QUAL threshold
min_rvaf = 0.00005       # and the rVAF threshold

[detect]
max_node = 20000

[output]
format = "tsv"
```

## Benchmarking

kmerdet includes a built-in benchmarking framework for evaluating detection accuracy against ground truth.

```bash
# Compare detection results against known truth
kmerdet benchmark --results detect.tsv --truth ground_truth.tsv

# With VAF bins and threshold sweep
kmerdet benchmark --results detect.tsv --truth ground_truth.tsv \
    --vaf-bins "0,0.001,0.01,0.1,1.0" \
    --sweep-vaf 0.0001,0.001,0.005,0.01,0.05
```

See `docs/benchmarking/` for the full benchmarking framework, including simulated data generation, analysis scripts, and plotting tools.

### Ground Truth Philosophy

kmerdet is designed to **exceed** km's detection capabilities. km output must never be used as ground truth. Ground truth must come from independent sources:

1. **Orthogonal validation** — ddPCR, Sanger sequencing
2. **Gold-standard alignment workflow** — BWA-MEM2 + GATK/Mutect2
3. **Simulated data** — Synthetic variants with known VAF

## Diagnostic Reports

Use `--report-dir` to generate detailed diagnostic output for debugging and understanding detection decisions:

```bash
kmerdet detect -d sample.jf -T targets/ \
    --report-dir /tmp/report --report-level detailed
```

Report levels: `summary`, `standard`, `detailed`, `full`.

## Project Structure

```
src/
├── cli/         # Subcommand definitions (detect, filter, merge, stats, etc.)
├── kmer/        # 2-bit packed Kmer struct, encoding, canonical form
├── jellyfish/   # KmerDatabase trait + pure Rust .jf reader
├── walker/      # DFS k-mer walking, adaptive thresholds, bidirectional
├── graph/       # Directed weighted graph, Dijkstra, pruning pipeline
├── sequence/    # Target FASTA loading, RefSeq decomposition
├── variant/     # Classifier, NNLS quantifier, bootstrap CIs, consensus
├── confidence/  # Binomial p-values, Fisher's method, Phred QUAL
├── filter/      # Tumor-informed filtering (reference + alt-sequence mode)
├── output/      # Multi-format writers (TSV, CSV, VCF, JSON, Excel)
└── report/      # Diagnostic report system
```

## Testing

```bash
cargo test                     # All tests
cargo test -- --ignored        # Run ignored tests
```

Tests use a `MockDb` (in `tests/common/mod.rs`) that simulates a jellyfish database without requiring external dependencies.

## Documentation

- `docs/features/` — Feature specifications with acceptance criteria
- `docs/research/` — Research documents on algorithms, sensitivity, and error correction
- `docs/benchmarking/` — Benchmarking framework and dataset generation
- `docs/planning/` — Development roadmap and open questions

## License

MIT
