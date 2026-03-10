# Feature: `coverage` -- K-mer Coverage Reporting

## What It Does

The `coverage` subcommand reports k-mer coverage statistics for target sequences without
performing variant detection. It queries the jellyfish database for reference k-mer
counts across each target, computes summary statistics, assesses GC bias, and optionally
produces a k-mer spectrum histogram. This serves as a quality control step before
detection and as a diagnostic tool for understanding sample characteristics.

```
kmerdet coverage -d sample.jf -t targets/ -o coverage_report.tsv
kmerdet coverage -d sample.jf -t targets/ --spectrum -o spectrum.tsv
```

### Input Specification

| Input | Flag | Format | Required |
|-------|------|--------|----------|
| Jellyfish database | `-d` / `--database` | `.jf` binary | Yes |
| Target sequences | `-t` / `--targets` | Directory of .fa or multi-FASTA | Yes |
| Output path | `-o` / `--output` | File path or stdout | No (stdout) |
| Format | `--format` | tsv / csv / json | No (tsv) |

### Output Specification

#### Per-Target Coverage Report

| Column | Type | Description |
|--------|------|-------------|
| Target | string | Target FASTA filename |
| total_kmers | u64 | Number of reference k-mers in the target |
| mean_count | f64 | Mean k-mer count across the target |
| median_count | u64 | Median k-mer count |
| min_count | u64 | Minimum k-mer count (coverage bottleneck) |
| max_count | u64 | Maximum k-mer count |
| std_dev | f64 | Standard deviation of counts |
| cv | f64 | Coefficient of variation (std_dev / mean) |
| gc_content | f64 | GC fraction of the target sequence |
| zero_count_kmers | u64 | Number of reference k-mers with count 0 (missing) |
| anchor_first_count | u64 | Count of the first k-mer (source anchor) |
| anchor_last_count | u64 | Count of the last k-mer (sink anchor) |
| status | string | OK, LOW_COVERAGE, MISSING_ANCHOR, HIGH_CV |

#### Sample-Level Summary

| Statistic | Description |
|-----------|-------------|
| sample_median_coverage | Median of per-target median counts |
| sample_mean_coverage | Mean of per-target mean counts |
| sample_cv | CV across all reference k-mers |
| targets_above_threshold | Count of targets with median > threshold |
| targets_below_threshold | Count of targets with median < threshold |
| estimated_depth | Estimated sequencing depth from k-mer counts |

---

## Why It Matters

### Pre-Detection Quality Control

Running `coverage` before `detect` reveals problems that would cause detection failures
or unreliable results:

1. **Missing anchors**: If a target's first or last k-mer has count 0, the walking
   algorithm cannot start or cannot reach the sink. The target will fail during
   detection. Identifying these targets in advance allows the user to adjust target
   design or flag expected failures.

2. **Low coverage targets**: Targets with median k-mer counts below the walking
   threshold (count=2) will not produce variant calls even if variants are present.
   These targets are effectively invisible to the detection algorithm.

3. **High variability**: Targets with high coefficient of variation (CV > 1.0) have
   pockets of very low coverage where variants may be missed, even though the median
   coverage appears adequate.

4. **GC bias**: Systematic coverage differences correlated with GC content indicate
   library preparation or capture bias that should be accounted for in the detection
   thresholds.

### Adaptive Threshold Computation

The adaptive-filtering.md research document proposes setting detection thresholds based
on coverage statistics rather than using fixed values. The `coverage` subcommand provides
the measurements needed for this adaptation:

| Coverage metric | Used for |
|----------------|----------|
| Sample median coverage | Setting depth-proportional count threshold |
| Per-target median | Setting per-target adaptive thresholds |
| Coverage CV | Adjusting ratio threshold (high CV -> lower ratio) |
| GC bias slope | Adjusting per-target expected coverage |
| K-mer spectrum shape | Estimating noise floor for adaptive thresholding |

The recommended depth tiers from the research document:

| Tier | Median Coverage | Recommended count | Recommended ratio |
|------|----------------|-------------------|-------------------|
| Ultra-low | <500x | 2 | 0.00001 |
| Standard | 500-2000x | 2 | 0.0001 |
| High-depth | 2000-10000x | 5 | 0.0003 |
| Ultra-deep | >10000x | 10 | 0.001 |

The `coverage` subcommand can output the recommended tier and parameters:

```
kmerdet coverage -d sample.jf -t targets/ --recommend-params
```

```
Sample coverage tier: High-depth (median: 3,450x)
Recommended parameters:
  --count 5
  --ratio 0.0003
  --min-vaf 0.001
```

---

## Coverage Depth Estimation

### From K-mer Counts to Sequencing Depth

The relationship between k-mer count and sequencing depth depends on the read length
and k-mer length:

```
expected_kmer_count = depth * (read_length - k + 1) / read_length
```

For 150 bp paired-end reads with k=31:
```
expected_kmer_count = depth * 120/150 = depth * 0.8
```

So a median k-mer count of 4,000 corresponds to an estimated sequencing depth of
~5,000x. This relationship is approximate because it assumes uniform coverage and
does not account for paired-end overlap, but it provides a useful ballpark.

### Factors Affecting the Estimate

| Factor | Effect | Correction |
|--------|--------|-----------|
| Paired-end overlap | Inflates count for short fragments | Subtract overlap contribution |
| PCR duplication | Inflates count | Use post-dedup database |
| GC bias | Variable count by target | Normalize by GC model |
| Capture efficiency | Variable count by target | Compare within-panel |
| Off-target reads | Dilutes count for targeted regions | Not applicable (jellyfish counts all) |

For clinical reporting, the estimated depth should be presented alongside the k-mer
count to give context: "median k-mer count: 4,000 (estimated depth: ~5,000x)" is more
interpretable than the raw k-mer count alone.

---

## GC Bias Assessment

### Per-Target GC Analysis

For each target, compute the GC content and plot it against the median k-mer count.
Ideally, all targets should have similar coverage regardless of GC content. In practice,
library preparation and capture introduce GC bias.

```
kmerdet coverage -d sample.jf -t targets/ --gc-analysis -o gc_bias.tsv
```

Output columns:

| Column | Description |
|--------|-------------|
| target | Target name |
| gc_content | GC fraction (0.0-1.0) |
| median_count | Median reference k-mer count |
| expected_count | Expected count from GC bias model |
| gc_ratio | median_count / expected_count (1.0 = no bias) |

### GC Bias Model

Fit a LOESS (locally estimated scatterplot smoothing) curve to the (GC content, median
count) data points across all targets. The residuals from this curve represent
target-specific coverage deviations not explained by GC content.

The GC bias model can be exported for use by the `detect` subcommand when adaptive
thresholds are enabled:

```
kmerdet coverage -d sample.jf -t targets/ --export-gc-model gc_model.json
kmerdet detect -d sample.jf -t targets/ --gc-model gc_model.json --adaptive
```

This connects coverage QC to detection parameter tuning, implementing the per-target
threshold adjustment described in adaptive-filtering.md Section 4.

---

## K-mer Spectrum Analysis

### Count Distribution Histogram

The k-mer spectrum is the distribution of k-mer counts across all queried k-mers. For
targeted sequencing, the spectrum of reference k-mers shows the coverage distribution:

```
kmerdet coverage -d sample.jf -t targets/ --spectrum -o spectrum.tsv
```

Output: a two-column TSV with (count, frequency) pairs:

```tsv
count	frequency
1	450
2	1230
3	890
...
3000	12
3001	15
...
```

### Interpreting the Spectrum

The k-mer spectrum reveals the noise structure of the sample. A typical targeted
sequencing spectrum has two peaks:

1. **Low-count peak (count 1-5)**: Error k-mers from sequencing errors. Each error
   creates k novel k-mers, most appearing once or twice. After UMI deduplication,
   this peak is reduced but not eliminated.

2. **Main peak (at coverage depth)**: Genuine genomic k-mers centered at the
   sequencing depth. The width of this peak reflects coverage variability.

The **valley** between these peaks represents the natural threshold separating noise
from signal. Tools from genome assembly (KmerGenie, GenomeScope) use this valley to
automatically set count thresholds. For variant detection, the valley position indicates
the noise floor: variant k-mers with counts below the valley are indistinguishable
from noise.

### Noise Floor Estimation

From the error k-mer distribution, estimate the noise floor:

```
noise_floor = percentile_99(error_kmer_counts)
```

Where error k-mers are identified as non-reference extensions at reference positions
(the three incorrect nucleotide extensions at each reference k-mer). This noise floor
estimate feeds into the adaptive thresholding described in adaptive-filtering.md
Section 2.

---

## Target Quality Annotations

Beyond coverage statistics, the `coverage` subcommand computes additional quality
annotations per target:

| Annotation | Description | Source |
|------------|-------------|--------|
| linguistic_complexity | Ratio of distinct k-mers to total k-mers | Target sequence |
| max_homopolymer | Length of longest homopolymer run | Target sequence |
| anchor_unique | Whether first/last k-mers appear unique in database | Database query |
| estimated_gc_bias | Coverage deviation attributable to GC content | GC model |

These annotations can be exported as a target catalog sidecar file:

```
kmerdet coverage -d sample.jf -t targets/ --annotate-targets -o target_annotations.json
```

This sidecar file can be reused across samples (the target-intrinsic annotations like
GC content and homopolymer length do not change between samples) and consumed by
`detect` for per-target parameter adjustment.

---

## Research Backing

### Adaptive Filtering Foundation

The adaptive-filtering.md research document describes a multi-stage adaptive filtering
pipeline where the first stage is "noise estimation from coverage distribution." The
`coverage` subcommand implements this first stage: it measures the coverage distribution,
estimates the noise floor, and recommends depth-proportional thresholds.

The document's Section 2 (Estimating Noise from Reference K-mer Coverage) outlines the
approach of querying non-reference extensions to characterize the error k-mer
distribution. The coverage subcommand implements this by sampling non-reference
extensions at reference k-mer positions and fitting a gamma distribution to their count
distribution.

### Room for Improvement

The room-for-improvement.md document identifies "Fixed Parameters Not Adapting to Sample
Characteristics" as limitation #6, severity High. Coverage-based parameter adaptation
is the recommended fix (P1 priority). The `coverage` subcommand provides the measurement
infrastructure for this adaptation.

### False Positive Analysis

The false-positive-analysis.md document shows that at high coverage (>5000x), error
k-mers at systematic error positions (GGC motifs) can have counts of 25-50, well above
typical detection thresholds. The coverage subcommand's noise floor estimation
identifies these high-count error positions, enabling position-specific thresholding
that maintains specificity at high depth.

### Sensitivity Landscape

The sensitivity-landscape.md document models the detection probability as a function of
VAF, sequencing depth, and k-mer length. The coverage subcommand provides the depth
measurement needed to evaluate where a sample falls on this landscape and what minimum
VAF is detectable.

---

## Design Considerations

### Performance

The coverage subcommand queries the jellyfish database for every reference k-mer in every
target. For a 50-target panel with ~200 bp per target and k=31, this is approximately
50 * 170 = 8,500 queries. At ~100 ns per jellyfish query (hash table lookup), the total
query time is <1 ms. The dominant cost is database loading (~100 ms) and output
formatting.

For a very large panel (10,000 targets), the query count is ~1.7 million, still
completing in <200 ms. Coverage computation is always fast relative to detection.

### Parallelism

Target-level parallelism via rayon is available but rarely needed given the sub-second
runtime. It is implemented for consistency with the `detect` subcommand and becomes
relevant only for extremely large panels.

### Integration with `detect`

The `coverage` subcommand can be invoked implicitly by `detect` when `--adaptive` is
enabled. In this mode, `detect` runs coverage analysis first, computes adaptive
thresholds, then proceeds with walking using the computed thresholds. This eliminates
the need for a separate coverage step in automated pipelines:

```
kmerdet detect -d sample.jf -t targets/ --adaptive -o results.tsv
```

### Comparison with Manual Jellyfish Queries

Coverage statistics can also be computed manually:

```bash
jellyfish query sample.jf ACGTACGTACGTACGTACGTACGTACGTACG
```

The `coverage` subcommand automates this for all reference k-mers across all targets,
computes summary statistics, and produces structured output. Acceptance testing validates
that kmerdet coverage statistics match manual jellyfish queries.

---

## Acceptance Criteria

### Unit Tests

- [ ] Mean, median, min, max, std_dev, CV computed correctly for known count vectors
- [ ] GC content computation matches hand calculation
- [ ] Zero-count k-mer detection: targets with missing k-mers are flagged
- [ ] Anchor validation: missing anchor k-mers flagged as MISSING_ANCHOR
- [ ] Depth estimation formula produces correct values
- [ ] Noise floor estimation from error k-mer distribution
- [ ] Status classification: LOW_COVERAGE, HIGH_CV, MISSING_ANCHOR thresholds

### Integration Tests

- [ ] Coverage statistics match manual jellyfish query results for test targets
- [ ] GC bias assessment produces reasonable bias curves on real data
- [ ] K-mer spectrum output matches jellyfish histo for the same database
- [ ] Recommended parameters are consistent with the depth tier table
- [ ] Target annotations are correct for targets with known properties
- [ ] Coverage report is consistent across multiple runs (deterministic)

### Performance Tests

- [ ] 50-target panel completes in <1 second
- [ ] 10,000-target panel completes in <5 seconds
- [ ] Memory usage proportional to number of targets, not database size

### Edge Cases

- [ ] Target with all k-mers having count 0 (completely missing)
- [ ] Target with uniform counts (CV = 0)
- [ ] Target consisting entirely of a homopolymer (all k-mers identical)
- [ ] Single-target input
- [ ] Very short target (fewer than k+1 bases, no valid k-mers)
- [ ] Database built with different k than specified (error with helpful message)
