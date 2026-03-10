# kmtools Features to Incorporate

## Overview

kmtools is the Python companion toolset within the kam pipeline that extends km's functionality with parallelism, filtering, merging, statistics, and plotting. These features need to be reimplemented as native subcommands in kmerdet.

- **Language**: Python
- **Location**: Part of the kam repository (https://github.com/trentzz/kam)
- **Entry point**: `kmtools` CLI with subcommands (`chunk`, `filter`, `merge`, `plot`, `stats`)
- **Dependencies**: pandas, matplotlib, openpyxl, click (CLI framework)

## kmtools chunk: Parallel km Execution

### Original Behavior

`kmtools chunk` splits km execution across multiple threads by distributing target FASTA files across workers:

1. Discovers all subdirectories within `--km-target-directory` (created by `refolder`)
2. Assigns subdirectories to threads in a round-robin fashion
3. Each thread calls `km find_mutation` as a subprocess for all `.fa` files in its subdirectory
4. Collects stdout from each km invocation
5. Optionally merges results into a single output file (`--merge`)

### Rust Replacement Strategy

In kmerdet, `kmtools chunk` is **not needed as a separate tool**. The `detect` subcommand natively parallelizes using rayon:

- **No file-based partitioning**: rayon's work-stealing scheduler distributes targets dynamically
- **No subprocess spawning**: Detection runs in-process, avoiding process creation overhead
- **Shared database access**: The jellyfish database handle is shared (read-only) across threads
- **Dynamic load balancing**: rayon automatically balances work; no need for `refolder`

The `refolder` tool is also unnecessary — kmerdet reads all `.fa` files from a directory (optionally recursive) without requiring pre-organized subdirectories.

## kmtools filter: Tumor-Informed Filtering

### Reference Mode

Matches variants by genomic coordinates and allele identity:

**Matching columns** (all must match for a positive hit):
- `CHROM` — chromosome
- `POS` — genomic position
- `REF` — reference allele
- `ALT` — alternate allele
- `TYPE` — variant type (Substitution, Insertion, Deletion, etc.)

**Reference information file format**:

```tsv
CHROM	POS	REF	ALT	TYPE	GENE	ANNOTATION
chr17	7577120	C	T	Substitution	TP53	R248W
chr5	170837543	TCTG	T	Deletion	NPM1	W288fs
```

Filtering logic: For each km result row, check if (CHROM, POS, REF, ALT, TYPE) exists in the reference. Mark matching rows as `PASS`, non-matching as `FILTERED`.

### Use-Alt Mode

Matches variants by the full alternative sequence (`ALT_SEQUENCE`) rather than genomic coordinates. Used when coordinate-based representation is ambiguous, particularly for complex indels where km and the reference may represent the same variant differently.

### Both Modes Needed in kmerdet

```rust
enum FilterMode {
    Reference,  // Match by CHROM/POS/REF/ALT/TYPE
    UseAlt,     // Match by full ALT_SEQUENCE
}
```

### Additional Filter Features to Implement

- **VAF threshold filtering**: Remove calls below a minimum rVAF (`--min-vaf 0.001`)
- **Type exclusion**: Remove specific variant types (`--exclude-type Reference`)
- **Expression threshold**: Remove calls with expression below a minimum (`--min-expression 5.0`)
- **Multi-reference support**: Accept multiple reference files for joint filtering

## kmtools merge: Combining Results

### Original Behavior

`kmtools merge` combines results from multiple parallel km runs:

1. **File discovery**: Find all result files in the input directory (recursive)
2. **Header validation**: Verify all files have the same column headers
3. **Header handling**: Write header once from the first file, skip headers from subsequent files
4. **Concatenation**: Append rows from all files in deterministic order (sorted by filename)
5. **Deduplication**: Remove exact duplicate rows (same Database + Query + Variant_name)

### Rust Implementation Notes

- Deduplication uses a HashSet of `(database, query, variant_name)` tuples; first occurrence wins
- Header mismatch produces a typed `KmerdetError::HeaderMismatch` with expected vs found columns
- In kmerdet, merge is less critical since `detect` already produces unified output from parallel target processing, but the subcommand is still needed for merging results across different samples or runs

## kmtools plot: Visualization

### Original Behavior

`kmtools plot` generates matplotlib-based visualizations from km output:

| Plot | Description | Chart Type |
|------|-------------|------------|
| VAF histogram | Distribution of rVAF values across all detected variants | Histogram |
| Detection bar chart | Count of detected vs. not-detected targets | Grouped bar |
| Type distribution | Breakdown by variant type (Substitution, Insertion, etc.) | Bar chart |
| Summary pie | Overall detection rate as pass/fail/filtered | Pie chart |

### Rust Replacement with plotters

Replace matplotlib with the `plotters` crate for SVG and PNG output:

- **SVG**: Default — scalable, editable, web-friendly
- **PNG**: For reports and presentations (via `plotters` bitmap backend)

**Plots to implement** (matching kmtools):
1. `vaf_histogram.svg` — rVAF distribution with configurable bin width
2. `detection_summary.svg` — detected vs. not-detected count per target
3. `type_distribution.svg` — variant type breakdown
4. `summary_pie.svg` — overall detection rate

## kmtools stats: Detection Statistics

### Statistics Computed

| Statistic | Description |
|-----------|-------------|
| Total targets | Number of target sequences queried |
| Detected variants | Number of non-reference calls |
| Detection rate | Detected / total (percentage) |
| VAF mean/median/min/max | Summary statistics for rVAF values |
| VAF percentiles | p5, p25, p50, p75, p95 |
| Per-type counts | Breakdown by variant type |
| Per-target summary | Detection status for each target |
| Expression summary | Mean/median expression across variants |

Output as formatted text to stdout (default) or JSON (`--format json`) for programmatic consumption.

## Custom Exceptions and Error Handling

### kmtools Error Patterns

| Exception | When Raised | Rust Equivalent |
|-----------|-------------|-----------------|
| `FileNotFoundError` | Input file does not exist | `KmerdetError::FileNotFound` |
| `InvalidHeaderError` | Column headers do not match expected | `KmerdetError::HeaderMismatch` |
| `EmptyFileError` | Input file has no data rows | `KmerdetError::EmptyFile` |
| `FilterModeError` | Invalid filter mode specified | `KmerdetError::InvalidFilterMode` |
| `MergeError` | Header mismatch during merge | `KmerdetError::HeaderMismatch` |

## File Validation

### Validation Checks from kmtools

kmtools performs several validation steps before processing:

1. **File existence**: Check all input files exist before starting
2. **Column header validation**: Verify expected columns are present and in correct order
3. **Format validation**: Check that numeric columns (rVAF, Expression, Min_coverage) contain valid numbers
4. **Empty file detection**: Warn or error on files with headers but no data rows
5. **Jellyfish DB validation**: Verify the `.jf` file is a valid jellyfish database (check header)

Validation runs **eagerly** — all inputs are validated before any processing begins. This avoids partial output from runs that would fail mid-way due to a bad input file.

## Development Branch vs Main Branch

### Feature Divergence

| Feature | Main Branch | Development Branch |
|---------|-------------|-------------------|
| `chunk` subcommand | Yes | Yes |
| `filter` subcommand | Yes | Yes (enhanced) |
| `merge` subcommand | Yes | Yes |
| `plot` subcommand | No | Yes |
| `stats` subcommand | No | Yes |
| Use-alt filter mode | No | Yes |
| Multi-reference filter | No | Yes |
| Excel output | No | Yes (via vcf2xlsx) |

For kmerdet, implement the **development branch feature set** — it represents the most complete version of the intended functionality.

### Feature Priority for kmerdet

| Priority | Feature | Rationale |
|----------|---------|-----------|
| P0 (core) | detect (replaces km + chunk) | Core functionality |
| P0 (core) | filter (both modes) | Required for clinical use |
| P1 (important) | merge | Needed for multi-run workflows |
| P1 (important) | stats | Essential for result interpretation |
| P2 (nice to have) | plot | Useful but not blocking |
| P2 (nice to have) | coverage | Diagnostic utility |

## References

- kam repository: https://github.com/trentzz/kam
- plotters crate: https://docs.rs/plotters/
- click (Python CLI): https://click.palletsprojects.com/
- pandas DataFrame API: https://pandas.pydata.org/docs/
