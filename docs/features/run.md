# Feature: `run` -- Full Pipeline Orchestration

## What It Does

The `run` subcommand executes the complete kmerdet analysis pipeline in a single
invocation: detect variants, filter results against a tumor-informed reference,
and optionally compute statistics and generate plots. It is the convenience
entry point for end-to-end clinical analysis, replacing the need to manually
chain `detect | filter | stats | plot` as separate commands.

```
kmerdet run --config kmerdet.toml \
    -d sample.jf \
    -t targets/ \
    -r reference-info.tsv \
    -o results/
```

This produces filtered variant calls, a summary statistics report, and
(optionally) visualization plots -- all from a single command.

## Why It Matters

### Clinical Workflow Simplicity

The thesis validation study demonstrated a complete pipeline execution time of
approximately 6 minutes per sample. In a clinical laboratory processing 20-50
samples per sequencing run, operators need a single reliable command rather than
a multi-step script. Every manual step introduces opportunities for error:
mis-ordered arguments, forgotten flags, accidentally piping unfiltered results
to the wrong output file.

A single `run` command with a validated TOML configuration file eliminates this
class of operator error. The configuration file can be version-controlled per
assay panel, ensuring reproducibility across runs and across sites.

### Automation and Pipeline Integration

Clinical sequencing pipelines (Nextflow, Snakemake, CWL) benefit from tools that
expose a single entry point with well-defined inputs and outputs. The `run`
subcommand serves as a self-contained processing node:

- **Input**: jellyfish database (.jf), target FASTA directory, reference info file
- **Output**: filtered results (TSV/VCF), statistics summary, optional plots
- **Exit code**: 0 (success), 1 (fatal error), 2 (partial failure with results)

Pipeline frameworks can wrap this single command without needing to understand
kmerdet's internal stage decomposition.

### Reproducibility

When all parameters live in a single TOML config file, the exact analysis
configuration is captured alongside results. Re-running the same config on the
same input produces identical output. This is essential for clinical validation
(CAP/CLIA requirements) and for debugging discrepancies between runs.

## Research Backing

The thesis (Ch5, Ch6) documented the 6-minute end-to-end workflow for the kam
pipeline: HUMID dedup (~2 min) + jellyfish count (~1.5 min) + km detect
(~2 min) + kamiltool filter (~0.5 min). The `run` subcommand covers the
detect + filter portion (~2.5 min), which is the analytically complex part.
The thesis recommended single-command execution for clinical deployment
(Ch6.5), noting that multi-step pipelines with manual parameter passing were
a source of errors during the validation study.

The tool design research (`08-tool-design-and-ux.md`) established the verb-based
subcommand pattern and the three-tier configuration precedence (CLI > config
file > defaults). The `run` subcommand is the composition point where all
per-stage configurations are unified in a single TOML file.

## Design Considerations

### Pipeline Stages and Ordering

The `run` subcommand executes stages in a fixed order:

```
Stage 1: VALIDATE
  - Verify all input files exist and are readable
  - Parse and validate the TOML config
  - Check jellyfish database integrity (header read)
  - Validate target FASTA files (parseable, non-empty)
  - Validate reference info file format

Stage 2: COVERAGE (optional, if adaptive thresholds enabled)
  - Query reference k-mer counts from the jellyfish database
  - Compute per-sample and per-target coverage statistics
  - Derive adaptive count/ratio thresholds

Stage 3: DETECT
  - For each target: k-mer walking, graph construction, pathfinding
  - Variant classification and NNLS quantification
  - Produces raw (unfiltered) variant calls

Stage 4: FILTER
  - Apply tumor-informed filtering (reference mode or alt-sequence mode)
  - Apply count/VAF/expression thresholds
  - Apply variant type filtering

Stage 5: STATS (optional, --stats flag)
  - Compute summary statistics: detection rate, variant type distribution,
    rVAF distribution, per-target success rate

Stage 6: PLOT (optional, --plot flag)
  - Generate rVAF scatter plots, coverage heatmaps, detection summary bar charts
```

Stages 1-4 always execute. Stages 5-6 are opt-in via flags or config.

### Configuration via TOML

All subcommand parameters are consolidated into a single TOML file. The `run`
subcommand reads sections `[detect]`, `[filter]`, `[output]`, and `[runtime]`:

```toml
[detect]
count = 2
ratio = 0.00001
max_stack = 500
max_break = 10
max_node = 10000

[filter]
mode = "reference"
min_vaf = 0.001
min_coverage = 2
exclude_types = ["Reference"]

[output]
format = "tsv"
output_dir = "./results"
prefix = "sample01"

[runtime]
threads = 4
verbose = 1
progress = true
```

CLI arguments override any config file value. For example,
`kmerdet run --config kmerdet.toml --count 5` uses count=5 regardless of the
config file's `[detect].count` setting.

### Input Requirements

| Input | Flag | Required | Description |
|-------|------|----------|-------------|
| Jellyfish DB | `-d` / `--db` | Yes | Path to .jf file (k-mer count database) |
| Targets | `-t` / `--targets` | Yes | Directory of target FASTA files or single multi-FASTA |
| Reference info | `-r` / `--reference` | Yes | TSV mapping targets to expected variants (tumor biopsy) |
| Config file | `--config` | No | TOML config (auto-discovered as `kmerdet.toml` in cwd) |

### Output Structure

When `--output-dir` is specified, `run` creates a structured output directory:

```
results/
  sample01_detect.tsv       # Raw detection results (all targets)
  sample01_filtered.tsv     # Filtered results (passed filter)
  sample01_stats.json       # Summary statistics (if --stats)
  sample01_timing.json      # Per-stage timing breakdown
  plots/                    # Visualization directory (if --plot)
    rvaf_scatter.svg
    coverage_heatmap.svg
    detection_summary.svg
```

### Intermediate Results

By default, the `run` subcommand streams data between stages in memory without
writing intermediate files. The `--save-intermediates` flag writes per-stage
output for debugging:

- `*_detect.tsv` -- raw unfiltered detection results
- `*_coverage.json` -- coverage statistics (if adaptive thresholds used)

When intermediates are not saved, data flows through an in-memory pipeline using
Rust's ownership model: the detect stage produces a `Vec<VariantCall>`, which is
consumed by the filter stage, which produces a `Vec<FilteredCall>`.

### Error Handling Modes

**Default (continue-on-error)**: Per-target failures (walking timeout, graph
construction failure, quantification error) are logged as warnings. The pipeline
continues with remaining targets. The exit code is 2 if any targets failed, 0 if
all succeeded.

**Strict mode (`--strict`)**: Any per-target failure is fatal. The pipeline
terminates immediately and returns exit code 1. Use this when every target must
succeed (e.g., a minimal MRD panel where missing any target is clinically
significant).

### Progress Reporting

The `run` subcommand displays progress across all stages using `indicatif`
multi-progress bars:

```
[Stage 2/4] Detecting variants
  [=========>           ] 23/50 targets  [00:42 / ~01:30]
```

Progress is suppressed when stderr is not a TTY (piped output). The `-q` flag
suppresses all non-error output.

### Timing Breakdown

Every `run` invocation produces a timing breakdown, either to stderr (at
`-v` verbosity) or to a JSON file (with `--save-intermediates`):

```json
{
  "total_seconds": 152.3,
  "stages": {
    "validate": 0.2,
    "coverage": 3.1,
    "detect": 128.7,
    "filter": 0.4,
    "stats": 1.2,
    "plot": 18.7
  },
  "targets_processed": 50,
  "targets_failed": 2,
  "targets_skipped": 0
}
```

This enables performance monitoring and identification of bottlenecks. The
detect stage dominates runtime; if its share exceeds 95%, the pipeline is
compute-bound. If validate or coverage stages are slow, there may be I/O issues
with the jellyfish database.

## Acceptance Criteria

### Functional Equivalence

1. **Output parity**: `kmerdet run --config c.toml -d s.jf -t targets/ -r ref.tsv`
   produces identical filtered output to running the stages separately:
   ```
   kmerdet detect --config c.toml -d s.jf -t targets/ -o raw.tsv
   kmerdet filter --config c.toml -i raw.tsv -r ref.tsv -o filtered.tsv
   ```
   Byte-for-byte identical output (given deterministic ordering).

2. **Config file loading**: All `[detect]`, `[filter]`, `[output]`, and
   `[runtime]` sections are correctly parsed and applied.

3. **CLI override**: Any CLI flag overrides the corresponding config file value.

### Error Handling

4. **Graceful degradation**: In default mode, failure on target N does not
   prevent processing of target N+1. Failed targets are listed in the summary.

5. **Strict mode**: `--strict` causes immediate termination on the first
   per-target failure, with exit code 1.

6. **Input validation**: Missing or unreadable input files produce a clear error
   message with the file path before any processing begins.

### Progress and Timing

7. **Timing output**: A timing breakdown is emitted at `-v` verbosity or higher,
   showing wall-clock seconds per stage.

8. **Progress bars**: Displayed when stderr is a TTY and `--progress` is not
   disabled. Show current stage, targets completed, and ETA.

### Output Completeness

9. **All stages produce output**: When `--stats` and `--plot` are enabled,
   the corresponding output files are created in the output directory.

10. **Intermediate files**: When `--save-intermediates` is set, raw detection
    results are written alongside filtered results.

### Performance

11. **No overhead vs separate invocation**: The `run` subcommand should not be
    measurably slower than running stages separately (within 5% wall-clock),
    since it avoids file I/O between stages.
