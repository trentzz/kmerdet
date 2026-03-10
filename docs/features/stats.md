# Feature: `stats` -- Summary Statistics

## What It Does

The `stats` subcommand computes detection statistics from variant call result files.
It produces summary metrics covering detection rates, VAF distributions, per-type
breakdowns, and expression statistics. Output can be formatted text (for human
consumption), JSON (for programmatic consumption), or TSV (for downstream analysis).

```
kmerdet stats -i filtered.tsv
kmerdet stats -i filtered.tsv --format json -o stats.json
kmerdet stats -i before.tsv --compare after.tsv
```

### Input Specification

| Input | Flag | Format | Required |
|-------|------|--------|----------|
| Result file | `-i` / `--input` | TSV/CSV from detect or filter | Yes |
| Comparison file | `--compare` | Second result file for comparison | No |
| Output path | `-o` / `--output` | File path or stdout | No (stdout) |
| Format | `--format` | text / json / tsv | No (text) |

---

## Why It Matters

### Clinical Reporting

For each patient timepoint, clinicians need a concise summary: how many targets were
queried, how many variants were detected, what is the overall detection rate, and what
are the VAF values. These numbers go directly into clinical reports. The thesis pipeline
computed these manually or via ad-hoc scripts; `stats` standardizes and automates
this reporting.

### Longitudinal Monitoring

MRD monitoring tracks patients over serial blood draws. The critical question is whether
the detection rate or mean VAF is changing over time. The `stats` subcommand, applied
to each timepoint's results, produces the numbers needed for longitudinal trend analysis.
A rising mean VAF or increasing detection rate may indicate relapse; a declining trend
confirms treatment response.

### Quality Control

Detection statistics serve as QC metrics for the sequencing run and the analysis pipeline.
Anomalous statistics (e.g., 0% detection rate, extremely high mean VAF, or unexpected
variant type distribution) flag problems with the sample, the library prep, or the
pipeline parameters.

---

## Statistics Computed

### Summary Statistics

| Statistic | Description | Computation |
|-----------|-------------|-------------|
| total_targets | Number of unique target sequences queried | Count of distinct Query values |
| total_variants | Number of non-reference variant calls | Count of rows where Type != Reference |
| detection_rate | Fraction of targets with at least one variant | targets_with_variants / total_targets |
| reference_calls | Number of reference-only calls | Count of rows where Type == Reference |
| failed_targets | Number of targets with no result (not even Reference) | total_targets - targets_with_any_result |

### VAF Distribution

| Statistic | Description |
|-----------|-------------|
| vaf_mean | Arithmetic mean of rVAF values (excluding Reference calls) |
| vaf_median | Median rVAF |
| vaf_min | Minimum rVAF |
| vaf_max | Maximum rVAF |
| vaf_std | Standard deviation of rVAF |
| vaf_p5 | 5th percentile |
| vaf_p25 | 25th percentile (Q1) |
| vaf_p50 | 50th percentile (median, same as vaf_median) |
| vaf_p75 | 75th percentile (Q3) |
| vaf_p95 | 95th percentile |
| vaf_iqr | Interquartile range (p75 - p25) |

Percentiles are computed using linear interpolation between adjacent ranks, matching
the numpy `percentile` function with `method='linear'` for reproducibility.

### Per-Type Breakdown

| Statistic | Description |
|-----------|-------------|
| substitution_count | Number of Substitution (SNV/MNP) calls |
| insertion_count | Number of Insertion calls |
| deletion_count | Number of Deletion calls |
| itd_count | Number of Internal Tandem Duplication calls |
| complex_count | Number of Complex Indel calls |
| reference_count | Number of Reference calls |
| substitution_mean_vaf | Mean rVAF for Substitutions |
| insertion_mean_vaf | Mean rVAF for Insertions |
| deletion_mean_vaf | Mean rVAF for Deletions |

### Per-Target Summary

For each target, report:
- Whether a variant was detected (yes/no)
- The variant type if detected
- The rVAF if detected
- The expression coefficient

This per-target view is essential for clinical reporting: the clinician needs to know
the status of each tracked mutation individually, not just aggregate statistics.

### Expression Summary

| Statistic | Description |
|-----------|-------------|
| expression_mean | Mean expression coefficient across variant calls |
| expression_median | Median expression |
| ref_expression_mean | Mean reference expression (proxy for coverage depth) |
| ref_expression_median | Median reference expression |
| expression_ratio_mean | Mean ratio of variant / reference expression |

---

## Output Formats

### Formatted Text (Default)

Human-readable summary printed to stdout:

```
Detection Summary
=================
Total targets:          50
Detected variants:      38
Detection rate:         76.0%
Reference-only:         12
Failed targets:         0

VAF Distribution
================
Mean:       0.0234
Median:     0.0189
Min:        0.0011
Max:        0.1520
Std dev:    0.0285
IQR:        0.0087 - 0.0342

Per-Type Counts
===============
Substitution:   32  (mean VAF: 0.0256)
Insertion:       3  (mean VAF: 0.0145)
Deletion:        2  (mean VAF: 0.0098)
ITD:             1  (mean VAF: 0.0340)

Per-Target Details
==================
TP53_R175H:     Substitution    rVAF=0.0234     Expression=23.4
NPM1_W288fs:    Insertion       rVAF=0.0145     Expression=14.5
FLT3_ITD:       ITD             rVAF=0.0340     Expression=34.0
BRCA1_ex11:     Reference       -               Expression=0.0
...
```

### JSON

Structured output for programmatic consumption:

```json
{
  "summary": {
    "total_targets": 50,
    "total_variants": 38,
    "detection_rate": 0.76,
    "reference_calls": 12,
    "failed_targets": 0
  },
  "vaf": {
    "mean": 0.0234,
    "median": 0.0189,
    "min": 0.0011,
    "max": 0.1520,
    "percentiles": { "p5": 0.0015, "p25": 0.0087, "p75": 0.0342, "p95": 0.0890 }
  },
  "per_type": { ... },
  "per_target": [ ... ]
}
```

### TSV

One row per statistic, suitable for appending to a longitudinal tracking spreadsheet:

```tsv
sample	timepoint	total_targets	detected	rate	vaf_mean	vaf_median	...
sample01	T1	50	38	0.76	0.0234	0.0189	...
```

---

## Comparison Mode

Compare two result files side by side to assess the effect of filtering, parameter
changes, or longitudinal changes:

```
kmerdet stats -i before_filter.tsv --compare after_filter.tsv
```

### Comparison Output

```
Comparison: before_filter.tsv vs after_filter.tsv
=================================================
                        Before      After       Delta
Total targets:          50          50          0
Detected variants:      42          38          -4
Detection rate:         84.0%       76.0%       -8.0%
Mean rVAF:              0.0256      0.0234      -0.0022

Variants removed by filter:     4
  TP53_art1:    Substitution    rVAF=0.0012    (below min_vaf threshold)
  KRAS_noise:   Substitution    rVAF=0.0008    (below min_vaf threshold)
  ...

Variants retained:              38
Variants gained:                0
```

This comparison is particularly useful for evaluating filter parameter tuning: how many
true positives are lost for each false positive removed.

---

## Per-Patient Longitudinal Tracking

When processing results from serial timepoints for the same patient, `stats` can
produce a longitudinal summary:

```
kmerdet stats --longitudinal -i T1.tsv T2.tsv T3.tsv T4.tsv \
              --timepoint-labels "Baseline,Week4,Week8,Week12"
```

### Longitudinal Output

```
Longitudinal Summary: Patient P001
===================================
Timepoint   Detection Rate   Mean VAF    Detected Mutations
Baseline    82.0%           0.0340       41/50
Week 4      64.0%           0.0120       32/50
Week 8      44.0%           0.0045       22/50
Week 12     26.0%           0.0018       13/50

VAF Trend: Declining (linear regression slope = -0.0027/timepoint, p < 0.001)
Interpretation: Treatment response -- consistent decline in ctDNA burden
```

This directly supports the clinical workflow described in the thesis: tracking
variant detection rates and VAF over time to monitor treatment response and detect
relapse.

---

## Research Backing

### Clinical Validation Context

The thesis validation study used summary statistics (detection rate, mean VAF) as the
primary metrics for evaluating pipeline performance. The 77% SNV sensitivity figure is
itself a detection rate statistic. Standardizing this computation in a dedicated
subcommand ensures consistency across analyses and eliminates manual calculation errors.

### Confidence Metrics

The confidence-metrics.md research document emphasizes that rVAF point estimates without
uncertainty are insufficient for clinical use. When bootstrap confidence intervals are
available (from the `detect` subcommand), `stats` includes CI summary statistics:

| Statistic | Description |
|-----------|-------------|
| vaf_ci_width_mean | Mean width of 95% CI across variants |
| vaf_ci_width_median | Median CI width |
| high_confidence_count | Variants with CI width < 0.01 |
| low_confidence_count | Variants with CI width > 0.05 |

### MRD Monitoring

The thesis proposed using k-mer detection as a rapid screening tool for MRD monitoring.
The longitudinal tracking feature in `stats` directly implements the clinical monitoring
workflow: track detection rate and mean VAF across serial blood draws to assess treatment
response and detect early relapse. A rising detection rate or increasing mean VAF after
initial decline (molecular relapse) triggers escalation to full alignment-based analysis.

---

## Design Considerations

### Handling Empty Input

Empty input (header only, no data rows) produces a valid statistics output with zero
counts and NaN/null for distribution statistics. This is not an error: a sample with
no detected variants is a meaningful clinical result (no evidence of residual disease).

### Handling Mixed Formats

The stats subcommand auto-detects the input format (TSV vs CSV) from the file extension
and the first line. It can also accept input from stdin (pipe from `detect` or `filter`),
defaulting to TSV.

### Precision

Floating-point statistics are reported with 4 decimal places by default (sufficient for
VAF values in the 0.001-0.100 range). The `--precision` flag overrides this for
applications requiring more or fewer digits.

---

## Acceptance Criteria

### Unit Tests

- [ ] Correct percentile calculation (linear interpolation) for small datasets
- [ ] Correct percentile calculation for large datasets
- [ ] Detection rate: correctly counts targets with at least one variant
- [ ] Per-type breakdown: correct counts for each variant type
- [ ] Empty input: produces valid output with zero counts
- [ ] Single-row input: all statistics computed correctly (no division by zero)
- [ ] All NaN values handled (e.g., std dev of single value)
- [ ] JSON output is valid JSON (parseable by serde_json)
- [ ] TSV output has correct column count and types

### Integration Tests

- [ ] Stats on detect output matches expected values for test dataset
- [ ] Stats on filtered output matches expected values
- [ ] Comparison mode correctly identifies added/removed/retained variants
- [ ] Longitudinal mode produces correct trend across timepoints
- [ ] Pipe compatibility: detect | filter | stats works end-to-end
- [ ] All three output formats (text, json, tsv) produce consistent numbers

### Edge Cases

- [ ] All variants are Reference type (detection rate = 0%)
- [ ] All variants are the same type (per-type has one entry)
- [ ] rVAF = 0.0 for some variants (included in distribution)
- [ ] Very large rVAF (>1.0, theoretically impossible but should not crash)
- [ ] Duplicate variant names in input (counted separately)
- [ ] Missing Expression column (older format) handled gracefully
