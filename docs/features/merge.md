# Feature: `merge` -- Result Merging and Deduplication

## What It Does

The `merge` subcommand combines variant detection results from multiple files into a
single unified output. It validates headers across input files, deduplicates by a
composite key, and produces deterministic sorted output. The primary use cases are
combining results from different samples, different detection runs, or different
k-mer lengths.

```
kmerdet merge -i run1.tsv run2.tsv run3.tsv -o merged.tsv
kmerdet merge -i results_dir/ -o merged.tsv --recursive
```

### Input Specification

| Input | Flag | Format | Required |
|-------|------|--------|----------|
| Result files | `-i` / `--input` | One or more TSV/CSV files, or a directory | Yes |
| Recursive scan | `--recursive` | Scan subdirectories for result files | No |
| Output path | `-o` / `--output` | File path or `-` for stdout | No (stdout) |

### Output Specification

Same format as input. One header row followed by deduplicated, sorted data rows.
The output includes all columns from the input files. When input files have different
column sets (e.g., one has Confidence and another does not), the output includes the
union of all columns with empty values for missing fields.

---

## Why It Matters

### Multi-Sample Workflows

In longitudinal MRD monitoring, each blood draw produces a separate .jf database and a
separate set of detection results. Clinicians need a consolidated view across all
timepoints for a patient to track VAF trends. The merge subcommand combines results
from all timepoints into a single file that can be analyzed by `stats` or `plot`.

### Multi-Run Workflows

When the same sample is processed multiple times (e.g., with different parameters,
or during pipeline development and validation), results from each run need to be
compared or combined. Merging with deduplication ensures that repeated runs do not
inflate variant counts.

### Multi-k Workflows

The multi-k-strategy.md research document describes running detection at multiple k-mer
lengths (k=21, k=31, k=41) and merging the results. Each k value may detect different
variants or the same variant with different rVAF estimates. The merge subcommand handles
cross-k merging with specific conflict resolution rules.

---

## Header Validation

Before processing any data, `merge` validates that all input files have compatible
column headers. This catches accidental mixing of files from different pipeline versions
or different output formats.

### Validation Rules

1. **Required columns**: All input files must contain at minimum: Database, Query, Type,
   Variant_name, rVAF, Expression, Min_coverage.
2. **Column order**: Columns need not be in the same order across files. Merging
   reorders to the canonical column order.
3. **Extra columns**: Files may have additional columns (e.g., Confidence from newer
   kmerdet versions). These are preserved and appear in the output.
4. **Missing columns**: If a required column is missing from any file, the merge fails
   with `KmerdetError::HeaderMismatch` listing expected vs. found columns.

### Strict vs. Permissive Mode

| Mode | Flag | Behavior |
|------|------|----------|
| Strict (default) | none | All files must have identical column sets |
| Permissive | `--permissive` | Union of columns; missing values filled with empty |

---

## Deduplication

### Composite Key

Deduplication is based on the composite key `(Database, Query, Variant_name)`. Two
rows with the same key are considered duplicates. The first occurrence (in file order,
then row order) wins.

```
Key = (database_path, query_filename, variant_name_string)
```

### Why This Key

- **Database**: Identifies the sample (each .jf comes from one sample). Including
  the database path ensures that the same variant detected in different samples is
  not deduplicated.
- **Query**: Identifies the target. The same variant at different targets would have
  different query filenames.
- **Variant_name**: Identifies the specific variant. Encodes chromosome, position,
  reference allele, and alternate allele.

### First-Occurrence Wins

When duplicates are encountered, the first occurrence is retained. This means file
order matters: files listed first on the command line have priority. For multi-k merges,
this allows the user to prioritize the primary k value:

```
# k=31 results take priority over k=21 for duplicate variants
kmerdet merge -i results_k31.tsv results_k21.tsv -o merged.tsv
```

---

## Sorting

Output rows are sorted deterministically to ensure reproducibility:

### Default Sort Order

1. Database (lexicographic)
2. Query (lexicographic)
3. Variant_name (lexicographic, which produces chromosome-then-position ordering
   for standard variant name formats like `chr17:7577120:C/T`)

### Custom Sort

```
kmerdet merge -i files/ --sort-by rVAF --sort-order desc -o merged.tsv
```

Supported sort keys: Database, Query, Type, Variant_name, rVAF, Expression,
Min_coverage.

---

## Multi-Sample Merge

When combining results from different .jf databases (different samples), each row's
Database field identifies its source sample. The merged output contains all variants
from all samples, with no cross-sample deduplication (since the same variant in
different samples represents independent observations).

### Sample Identification

By default, the Database column contains the full path to the .jf file. For cleaner
output, `--strip-database-path` retains only the filename, and `--sample-name` maps
database paths to sample identifiers:

```
kmerdet merge -i results/ --sample-name sample01=path/to/s01.jf \
                          --sample-name sample02=path/to/s02.jf -o merged.tsv
```

---

## Multi-k Merge

When combining results from different k-mer lengths, additional handling is needed
because the same variant may be detected at multiple k values with different rVAF
estimates and coverage metrics.

### Variant Matching Across k Values

Two calls from different k values match if they describe the same genomic variant after
INDEL normalization. The matching uses the same normalization logic as the filter
subcommand:

1. Normalize both variants (left-align INDELs)
2. Compare (CHROM, POS, REF, ALT, TYPE)
3. If all five match, the variants are the same

### Conflict Resolution for Duplicate rVAF

When the same variant is detected at multiple k values with different rVAF estimates,
the merge must choose or combine:

| Strategy | Flag | Behavior |
|----------|------|----------|
| First wins (default) | none | Keep rVAF from the first file (primary k) |
| Best coverage | `--conflict best-coverage` | Keep the call with highest Min_coverage |
| Weighted average | `--conflict weighted-avg` | Weighted average of rVAF by coverage |
| Keep all | `--conflict keep-all` | Retain all calls, annotated with k value |

The `keep-all` strategy is useful for analysis: it shows how rVAF varies across k values,
which the multi-k-strategy.md document identifies as a quality signal (variants detected
consistently across k values are higher confidence).

### Multi-k Annotations

When merging multi-k results, the output includes additional annotation columns:

| Column | Description |
|--------|-------------|
| k_values_detected | Comma-separated list of k values that detected this variant |
| n_k_detected | Number of k values detecting this variant |
| consensus_tier | Confidence tier (1-4) based on detection pattern |

The consensus tier follows the voting framework from multi-k-strategy.md:
- Tier 1: Detected at all k values (highest confidence)
- Tier 2: Detected at 2 of 3 k values
- Tier 3: Detected at k=31 only (standard)
- Tier 4: Detected at only k=21 or k=41 (review recommended)

---

## Research Backing

### Multi-k Strategy

The multi-k-strategy.md research document provides the theoretical and practical
framework for combining detection results across k-mer lengths. The key insight is that
no single k value is optimal for all variant types: k=31 works well for SNVs, but shorter
k values (21-25) provide better sensitivity for large INDELs because more k-mers span
the insertion junction.

The consensus voting approach (Approach 2 in the research document) shows that a variant
detected at all three k values receives a 40% confidence bonus, while single-k detections
are flagged for review. Expected sensitivity improvements: SNVs +5-10%, short INDELs
+15-25%, long INDELs +20-40%.

### BAYSIC Precedent

Cantarel et al. (2014) demonstrated that Bayesian combination of variant calls from
multiple callers improves both sensitivity and specificity. The multi-k merge is
analogous: each k value acts as a different "caller" with different strengths, and
combining results extracts maximum information.

### INDEL Normalization for Cross-k Matching

Different k values may produce slightly different breakpoint coordinates for the same
INDEL, because the k-mer that first diverges from the reference depends on k. The Tan
et al. (2015) normalization algorithm, applied before cross-k matching, ensures that
equivalent representations converge to the same canonical form.

---

## Design Considerations

### Memory Efficiency

For small merges (a few files, <100K total rows), all rows are loaded into memory for
sorting and deduplication. For large merges (millions of rows), an external sort
strategy can be used: sort each file individually, then perform a k-way merge using
a min-heap.

In practice, kmerdet detection results are small (50-500 rows per sample for a typical
panel), so in-memory merging is always sufficient for clinical use cases.

### Determinism

The merge output must be byte-identical across runs with the same inputs in the same
order. This requires:
- Deterministic sort (stable sort with total ordering)
- Deterministic deduplication (first-occurrence wins, not arbitrary)
- Deterministic floating-point formatting (fixed decimal precision via `--precision`)

### Composability

Merge is designed to be chained:

```
# Merge per-timepoint results, then filter
kmerdet merge -i timepoint_*/ | kmerdet filter -r reference.tsv -o final.tsv
```

---

## Acceptance Criteria

### Unit Tests

- [ ] Two files with identical headers merge correctly
- [ ] Deduplication removes exact duplicates by (Database, Query, Variant_name)
- [ ] First-occurrence wins: first file's row is kept for duplicates
- [ ] Header mismatch in strict mode produces KmerdetError::HeaderMismatch
- [ ] Permissive mode fills missing columns with empty values
- [ ] Sorting: output is deterministically sorted by default sort order
- [ ] Custom sort: --sort-by rVAF produces correct ordering
- [ ] Empty files (header only, no data) handled gracefully
- [ ] Single input file passes through unchanged

### Integration Tests

- [ ] Multi-sample merge: variants from different databases are all preserved
- [ ] Multi-k merge: same variant at k=21 and k=31 is deduplicated
- [ ] Multi-k annotations: k_values_detected and consensus_tier are correct
- [ ] Conflict resolution: all strategies produce expected rVAF
- [ ] Directory input with --recursive finds all .tsv files
- [ ] Pipe compatibility: detect | merge works for streaming

### Edge Cases

- [ ] Files with different column orders but same column names
- [ ] Files with zero data rows (header only)
- [ ] Merging a single file (no deduplication needed)
- [ ] Duplicate rows within a single file (should deduplicate)
- [ ] Very long Variant_name strings (complex indels)
- [ ] Unicode in file paths and sample names
