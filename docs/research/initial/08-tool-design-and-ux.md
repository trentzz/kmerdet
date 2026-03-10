# Tool Design and UX

## CLI Design Philosophy

kmerdet follows a **verb-based subcommand** pattern, where each subcommand represents an action the user performs on their data. The CLI is designed for both interactive use (human at a terminal) and pipeline integration (scripted workflows, Nextflow).

### Subcommand Structure

```
kmerdet <SUBCOMMAND> [OPTIONS] [ARGS]
```

| Subcommand | Purpose | Primary Input | Primary Output |
|------------|---------|---------------|----------------|
| `detect` | Run k-mer walking and variant detection | .jf + target FASTA | Variant calls (TSV/VCF) |
| `filter` | Tumor-informed filtering of results | Detection output + reference | Filtered variant calls |
| `merge` | Combine results from multiple runs | Multiple detection outputs | Single merged output |
| `stats` | Compute detection statistics | Detection or filtered output | Summary statistics |
| `plot` | Generate visualization plots | Detection or filtered output | SVG/PNG images |
| `coverage` | Compute k-mer coverage statistics | .jf + target FASTA | Coverage report |
| `run` | Execute full pipeline (detect + filter) | .jf + targets + reference | Filtered variant calls |

### Naming Rationale

- **Verbs, not nouns**: `detect` not `detection`, `filter` not `filtering` — commands describe actions
- **Short and memorable**: Single words where possible
- **Composable**: Each subcommand does one thing; `run` composes them for convenience

## Configuration File Format

### TOML Configuration

kmerdet uses **TOML** as its configuration format, chosen for readability and alignment with the Rust ecosystem (Cargo.toml precedent).

```toml
# kmerdet.toml — example configuration

[detect]
kmer_length = 31
ratio = 0.00001
count = 2
max_stack = 500
max_break = 10
max_node = 10000
canonical = true

[filter]
mode = "reference"           # "reference" or "use-alt"
reference_file = "reference-info.tsv"
min_vaf = 0.001
exclude_types = ["Reference"]

[output]
format = "tsv"               # tsv, csv, vcf, json, jsonl, xlsx
output_dir = "./results"
prefix = "sample01"
include_reference_calls = false
decimal_precision = 6

[runtime]
threads = 4
verbose = 1                  # 0=warn, 1=info, 2=debug, 3=trace
progress = true
temp_dir = "/tmp/kmerdet"
```

### Configuration Precedence

Three-tier precedence model:

```
CLI arguments  >  Config file  >  Built-in defaults
   (highest)      (middle)        (lowest)
```

- **Built-in defaults**: Hardcoded in the application (e.g., k=31, ratio=0.05, format=tsv)
- **Config file**: Loaded from `kmerdet.toml` (explicit `--config` flag or auto-discovery in cwd)
- **CLI arguments**: Always win — any flag on the command line overrides everything else

## Output Format Strategy

| Format | Extension | Use Case | Crate |
|--------|-----------|----------|-------|
| **TSV** | `.tsv` | Default — backwards compatible with km output | `csv` (tab delimiter) |
| **CSV** | `.csv` | Spreadsheet import, general data exchange | `csv` |
| **VCF 4.3** | `.vcf` | Bioinformatics tools (bcftools, IGV, GATK) | `noodles-vcf` or manual |
| **JSON** | `.json` | Programmatic consumption, API responses | `serde_json` |
| **JSONL** | `.jsonl` | Streaming consumption, log aggregation | `serde_json` (line-delimited) |
| **Excel** | `.xlsx` | Clinical reporting, non-technical stakeholders | `rust_xlsxwriter` |

### TSV Format (Default)

Preserves km's original output columns for backwards compatibility:

```
Database	Query	Type	Variant_name	rVAF	Expression	Min_coverage	Start_offset	Sequence	Reference_expression	Reference_sequence	Info
```

### VCF 4.3 Output

Enables integration with standard bioinformatics tooling. Maps kmerdet fields to VCF INFO tags: `RVAF` (relative VAF), `EXPR` (expression coefficient), `MINCOV` (minimum k-mer coverage), `KMER` (k-mer length).

### JSON/JSONL Output

JSON provides structured output for programmatic consumers. JSONL writes one JSON object per line — suitable for streaming with `jq`, log ingestion, and piping without buffering entire arrays.

## Error Handling Patterns

### Library vs CLI Error Strategy

1. **Library layer** (`thiserror`): Typed, structured errors for programmatic handling
2. **CLI layer** (`anyhow`): Context-rich errors for human-readable messages

Library errors are defined as a `KmerdetError` enum with variants for each failure mode: `JellyfishDbNotFound`, `InvalidKmerLength`, `TargetParseError`, `MaxNodesExceeded`, `NoPathsFound`, `QuantificationError`.

CLI code wraps these with `anyhow::Context` to add file paths and operation context.

### Error Recovery Strategy

| Error Type | Behavior |
|------------|----------|
| Missing input file | Fail immediately with clear path in message |
| Invalid k-mer length | Fail immediately with valid range hint |
| Walking timeout (max_node) | Warn and skip target, continue with remaining targets |
| Quantification failure | Warn and report target as failed, continue |
| Output write failure | Fail immediately (partial output is worse than none) |
| Thread panic | Catch at rayon boundary, report, continue other targets |

Per-target errors are **non-fatal** by default — the tool processes all targets and reports failures in a summary. The `--strict` flag makes per-target errors fatal.

## Progress Reporting

Long-running operations use `indicatif` progress bars showing elapsed time, completion percentage, targets processed / total, and estimated time remaining.

Progress bars are automatically suppressed when stderr is not a TTY (piped or redirected output).

## Verbose/Quiet Modes

kmerdet uses the `tracing` crate for structured logging with escalating verbosity:

| Flag | Level | What is Logged |
|------|-------|----------------|
| (default) | `warn` | Warnings and errors only |
| `-v` | `info` | Summary statistics, file paths, parameter values |
| `-vv` | `debug` | Per-target progress, k-mer counts, graph sizes |
| `-vvv` | `trace` | Individual k-mer queries, extension decisions, path details |
| `-q` | `error` | Errors only (suppress warnings) |

Verbose counting uses `clap::ArgAction::Count` on the `-v` flag. All log output goes to stderr via `tracing_subscriber::fmt().with_writer(std::io::stderr)`.

## Exit Codes and I/O Conventions

### Exit Codes

| Code | Meaning |
|------|---------|
| `0` | Success — all targets processed, output written |
| `1` | General error — invalid arguments, missing files, fatal errors |
| `2` | Partial failure — some targets failed (logged to stderr), results written for successful targets |

### Stdout / Stderr Convention

- **stdout**: Data output only (variant calls, statistics, results)
- **stderr**: Human-facing messages (progress bars, logs, warnings, errors)

When `--output` / `-o` is specified, results are written to the named file and stdout remains available.

## Shell Completion Support

kmerdet generates shell completions via `clap_complete`, available through a `completions` subcommand:

```bash
kmerdet completions bash > /etc/bash_completion.d/kmerdet
kmerdet completions zsh > ~/.zfunc/_kmerdet
kmerdet completions fish > ~/.config/fish/completions/kmerdet.fish
```

Completions cover subcommand names, flag names/aliases, file path arguments (with extensions like `.jf`, `.fa`, `.tsv`), and enum values (output formats, filter modes).

## Example CLI Invocations

```bash
# Basic detection
kmerdet detect -d sample.jf -t targets/ -o results.tsv

# Detection with custom parameters and VCF output
kmerdet detect -d sample.jf -t targets/ --count 2 --ratio 0.00001 --format vcf -o results.vcf

# Full pipeline with config file
kmerdet run --config kmerdet.toml -d sample.jf -t targets/ -r reference-info.tsv

# Filter existing results
kmerdet filter -i km_results.tsv -r reference-info.tsv --mode reference -o filtered.tsv

# Generate statistics and plots
kmerdet stats -i filtered.tsv
kmerdet plot -i filtered.tsv -o plots/ --format svg

# Merge multiple result files
kmerdet merge -i run1.tsv run2.tsv run3.tsv -o merged.tsv
```

## References

- clap documentation: https://docs.rs/clap/
- indicatif documentation: https://docs.rs/indicatif/
- tracing documentation: https://docs.rs/tracing/
- thiserror documentation: https://docs.rs/thiserror/
- VCF 4.3 specification: https://samtools.github.io/hts-specs/VCFv4.3.pdf
