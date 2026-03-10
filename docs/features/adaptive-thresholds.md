# Feature: Adaptive Thresholds -- Sample-Aware Parameter Tuning

## What It Does

Adaptive thresholds automatically adjust variant detection and filtering
parameters based on the characteristics of each sample and each target region.
Instead of applying fixed `count=2, ratio=0.00001` across all samples and
targets, kmerdet estimates the coverage depth, noise floor, and genomic context
for each analysis unit, then sets thresholds that optimize the
sensitivity-specificity tradeoff for that specific context.

## Why It Matters

### The Fixed Threshold Problem

The thesis validation study used `count=2, ratio=0.00001` for all samples. These
values were chosen as a compromise: low enough to detect variants at the lowest
expected VAF (0.1% at ~3000x coverage), but high enough to avoid overwhelming
false positives. This compromise fails at both extremes:

**High-coverage samples (>5000x)**: At 10,000x coverage, a sequencing error
k-mer (per-base error rate ~0.1%) has an expected count of ~10. The fixed
count=2 threshold passes these errors freely. A sample at 10,000x has an error
noise floor roughly 5x higher than at 2000x, but the threshold does not adapt.
The result: dramatically more false positive variant calls at high depth.

**Low-coverage samples (<1000x)**: At 500x coverage, a true variant at 0.1% VAF
produces a k-mer with expected count ~0.5. The count=2 threshold eliminates
this signal. A more permissive count=1 with additional contextual evidence
could recover it, but the fixed threshold cannot express "this count-1 k-mer is
credible because its neighbors all have count 2-3."

**Between targets**: Within a single sample, capture efficiency and GC content
cause 3-5x coverage variation across targets. A target at 1500x effective
coverage operates in a fundamentally different regime than one at 8000x. The
thesis found that coverage CV within a sample significantly affected per-target
detection sensitivity.

### Impact on Clinical Performance

In a clinical lab processing samples with variable cfDNA input (5-50 ng),
extraction efficiency (30-80%), and library complexity, coverage varies 5-10x
across samples. Fixed thresholds tuned for the median sample produce:

- Excess false positives on high-quality, high-coverage samples (wasted follow-up)
- Missed true variants on low-quality, low-coverage samples (delayed detection)
- Inconsistent per-target sensitivity within a sample

Adaptive thresholds address all three failure modes by matching the threshold to
the data.

## Research Backing

### Thesis Recommendations (room-for-improvement.md)

The thesis explicitly recommended adaptive parameter adjustment (Limitation #6,
Priority P1). It provided a recommended parameter table:

| Median Coverage | Recommended count | Recommended ratio | Rationale |
|----------------|-------------------|-------------------|-----------|
| <500x | 2 | 0.00001 | Maximize sensitivity at low depth |
| 500-2000x | 2-3 | 0.00005 | Balanced |
| 2000-5000x | 3-5 | 0.0001 | Reduce noise at moderate depth |
| >5000x | 5-10 | 0.001 | Can afford strict thresholds |

Additional guidance: if coverage CV > 0.5, halve the effective ratio to avoid
losing low-coverage targets. If dedup rate > 80%, use more conservative
thresholds (low library complexity means effective molecule count is low).

### K-mer Spectrum Analysis (adaptive-filtering.md)

The k-mer frequency spectrum reveals the noise structure of each sample. The
low-count peak (error k-mers) and high-count peak (true k-mers) are separated
by a valley whose position varies by sample. Fitting a Gamma distribution to
the error k-mer counts and setting the threshold at the 99.9th percentile
automatically adapts to the sample's error rate.

This approach follows the precedent of Kelley et al. (Quake, 2010) for error
detection and Guo et al. (2024) for statistically-solid k-mer identification.

### Genomic Context Effects (adaptive-filtering.md)

Per-target context features that affect threshold selection:

- **GC content**: Extreme GC (>70% or <30%) reduces coverage 2-3x. Thresholds
  must be lowered proportionally. GCfix (Zheng et al., 2024) demonstrated that
  fragment-length-specific GC correction is needed for cfDNA data.

- **Homopolymer runs**: INDEL error rate increases exponentially with
  homopolymer length: ~0.01% at 1-3 bp, ~1% at 6-7 bp, ~5-10% at 8+ bp.
  Targets with long homopolymers need elevated INDEL thresholds.

- **Linguistic complexity**: Repetitive regions produce more background noise.
  Targets with complexity < 0.8 (ratio of distinct k-mers to total) should
  have elevated thresholds.

## Design Considerations

### Coverage Estimation

Coverage is estimated from reference k-mer counts -- the counts of k-mers
derived from the known target reference sequences. These k-mers are expected
to be present at the sequencing depth.

```rust
pub struct CoverageStats {
    pub median: u64,
    pub mean: f64,
    pub p5: u64,          // 5th percentile
    pub p95: u64,         // 95th percentile
    pub cv: f64,          // coefficient of variation
    pub dedup_rate: f64,  // estimated from coverage vs raw read count
}
```

The `coverage` subcommand already computes these statistics. For the `run`
subcommand, coverage estimation is integrated as an automatic first step before
detection.

### Depth-Proportional Count Threshold

The core formula:

```
adaptive_count = max(min_threshold, floor(median_coverage * min_detectable_vaf))
```

This ensures that at any coverage level, the threshold is set just high enough
to detect variants at the target VAF while filtering the majority of error
k-mers. Concrete examples:

| Median Coverage | min_vaf=0.001 | min_vaf=0.0005 |
|----------------|---------------|----------------|
| 500x | max(2, 0.5) = 2 | max(2, 0.25) = 2 |
| 2000x | max(2, 2) = 2 | max(2, 1) = 2 |
| 5000x | max(2, 5) = 5 | max(2, 2.5) = 3 |
| 10000x | max(2, 10) = 10 | max(2, 5) = 5 |
| 50000x | max(2, 50) = 50 | max(2, 25) = 25 |

### Noise Floor Estimation

Rather than relying on a fixed error rate, estimate the actual noise floor from
the data by querying non-reference extensions at reference k-mer positions:

```
For each reference k-mer R:
    Query all 4 possible forward extensions (RA, RC, RG, RT)
    The 3 non-reference extensions are predominantly errors
    Collect their counts as error_counts

Fit Gamma(shape, rate) to error_counts
Threshold = Gamma_quantile(0.999, shape, rate)
```

This approach adapts to platform-specific error rates, library prep artifacts,
and per-sample noise characteristics simultaneously. It is computed once per
sample and adds negligible runtime (a few thousand jellyfish queries).

### Per-Target Context Adjustment

After computing the sample-level threshold, adjust per target:

```rust
pub struct TargetAnnotation {
    pub gc_content: f64,
    pub linguistic_complexity: f64,
    pub max_homopolymer: usize,
    pub total_homopolymer_bases: usize,
    pub local_median_coverage: u64,
    pub adjusted_count_threshold: u32,
    pub adjusted_ratio_threshold: f64,
}
```

Adjustments:

1. **GC bias correction**: Fit a GC bias model from reference k-mer counts
   grouped by GC content. Low-coverage targets (due to GC bias) get lower
   thresholds; high-coverage targets get higher thresholds.

2. **Homopolymer penalty**: For targets with homopolymer runs >= 6 bp,
   multiply the INDEL-specific count threshold by 5. Flag INDEL calls near
   homopolymers.

3. **Complexity penalty**: For targets with linguistic complexity < 0.8,
   multiply the count threshold by 2.

4. **Local coverage override**: If a target's local reference k-mer coverage
   differs by more than 2x from the sample median, use the local coverage for
   threshold computation.

These annotations are precomputed once per panel design and stored as a sidecar
file (e.g., `panel_annotations.json`), so they do not need to be recomputed for
every sample.

### Walking Threshold Adaptation

The adaptive threshold applies not just to final filtering but also to the
k-mer walking extension step. The walking extension rule is:

```
threshold = max(sum_of_sibling_counts * ratio, count)
```

With adaptive parameters:

```
adaptive_ratio = max(estimated_error_rate / 3, base_ratio)
adaptive_count = max(min_threshold, floor(local_coverage * min_vaf))

threshold = max(sum_of_sibling_counts * adaptive_ratio, adaptive_count)
```

This means high-coverage regions use stricter walking thresholds (reducing
spurious branches), while low-coverage regions use more permissive thresholds
(preserving true variant paths).

### Coverage Tier Auto-Detection

For clinical deployment, define validated coverage tiers with pre-set
parameters. The auto-detected tier is logged and reported:

```rust
pub enum CoverageTier {
    UltraLow,   // <500x:    count=2, ratio=0.00001, min_vaf=0.005
    Standard,    // 500-2000x: count=2, ratio=0.0001, min_vaf=0.002
    HighDepth,   // 2000-10000x: count=5, ratio=0.0003, min_vaf=0.001
    UltraDeep,   // >10000x:  count=10, ratio=0.001, min_vaf=0.0005
}
```

Users can override the auto-detected tier with `--coverage-tier` or individual
parameter flags.

### Poisson Model for Extension Thresholds

At each walking extension, the expected number of error k-mers follows a
Poisson distribution with rate:

```
lambda_error = coverage * per_base_error_rate * (number of error possibilities)
```

For a single-base extension with 3 possible error outcomes:

```
lambda_error = coverage * 0.001 * 3 = coverage * 0.003
```

At 5000x coverage, lambda_error = 15. The 99.9th percentile of Poisson(15) is
~27. Setting the extension threshold at this level ensures that fewer than
0.1% of error extensions pass, while any true variant extension at VAF >= 0.1%
(expected count = 5) would need additional support from neighboring k-mers in
the walk.

This Poisson model replaces the heuristic ratio parameter with a statistically
grounded threshold.

## Acceptance Criteria

### Coverage Estimation

1. The `coverage` subcommand (and the coverage estimation step within `run`)
   correctly computes median, mean, P5, P95, and CV of reference k-mer counts.

2. Coverage statistics are computed per-sample and per-target.

3. The coverage tier is auto-detected and logged.

### Threshold Computation

4. The adaptive count threshold follows the formula:
   `max(min_threshold, floor(median_coverage * min_vaf))`.

5. The adaptive ratio threshold incorporates the estimated error rate from the
   k-mer spectrum.

6. Per-target adjustments for GC content, homopolymer runs, and linguistic
   complexity are applied when target annotations are available.

### Sensitivity Improvement

7. At low coverage (<1000x), adaptive thresholds achieve equal or better
   sensitivity compared to fixed thresholds, without increasing the false
   positive rate.

8. At high coverage (>5000x), adaptive thresholds reduce the false positive
   rate by at least 50% compared to fixed count=2, without reducing sensitivity
   for variants above the minimum detectable VAF.

### Clinical Safety

9. Adaptive thresholds never produce a threshold below the configured
   `min_threshold` floor (default: 2). This prevents ultra-permissive
   thresholds from noise model fitting failures.

10. The computed threshold and the rationale (coverage tier, local coverage,
    context adjustments) are reported in the output for clinical traceability.

### Integration

11. Adaptive thresholds work correctly with both the `detect` subcommand
    (when `--adaptive` is enabled) and the `run` subcommand (where coverage
    estimation is automatic).

12. CLI flags `--count` and `--ratio` override adaptive computation when
    explicitly provided, allowing users to force specific thresholds.
