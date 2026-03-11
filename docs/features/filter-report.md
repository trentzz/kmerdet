# Filter Report — Detailed Filtering Analysis

## What It Does

Generates a comprehensive report analyzing filtering outcomes: what was filtered out, why, threshold sensitivity, and per-criterion breakdown. Helps researchers understand why specific variants were or weren't found, and tune parameters.

## Why It Matters

Filtering is the most tunable part of the pipeline. Researchers need to understand:
- Which variants narrowly missed thresholds (could be rescued with parameter adjustment)
- Which filter criterion is most restrictive (bottleneck analysis)
- Whether thresholds are appropriate for this sample's depth/noise profile
- How detection rates change across threshold sweeps

## CLI Interface

```bash
# As part of filter command
kmerdet filter --results detect.tsv --expected expected.tsv --report report.json

# As part of run pipeline
kmerdet run --db sample.jf --targets-dir targets/ --expected expected.tsv --filter-report report.json

# Report formats
kmerdet filter --results detect.tsv --expected expected.tsv --report report.json   # JSON (default)
kmerdet filter --results detect.tsv --expected expected.tsv --report report.txt    # Human-readable text
```

## Report Contents

### Summary Section
- Total expected variants
- Found count / not-found count / detection rate
- Filter configuration used (all threshold values)

### Per-Criterion Breakdown
For each filter criterion (coverage, VAF, expression, type):
- How many variants were rejected by this criterion alone
- How many variants were rejected by this criterion in combination
- Distribution of values for found vs not-found variants
- Suggested threshold adjustments (if applicable)

### Near-Miss Analysis
Variants that were detected but failed filtering by a narrow margin:
- Coverage within 2x of threshold
- VAF within 2x of threshold
- These are candidates for threshold relaxation

### Per-Variant Detail
For each expected variant:
- Match status (found / not found / detected but filtered)
- If detected: all metric values (rVAF, coverage, expression)
- If filtered: which criteria failed, by how much
- If not detected: was the target processed? Did walking find any paths?
- Closest matching detected variant (if any partial match)

### Threshold Sensitivity
- Detection rate at threshold multiples: 0.5x, 1x (current), 2x, 5x
- For each criterion independently
- Helps identify the optimal operating point

## Data Structures

```rust
pub struct FilterReport {
    pub summary: FilterSummary,
    pub per_criterion: Vec<CriterionAnalysis>,
    pub near_misses: Vec<NearMiss>,
    pub per_variant: Vec<VariantDetail>,
    pub threshold_sensitivity: ThresholdSensitivity,
}

pub struct FilterSummary {
    pub total_expected: usize,
    pub found: usize,
    pub not_found: usize,
    pub detected_but_filtered: usize,
    pub not_detected: usize,
    pub detection_rate: f64,
    pub config: FilterConfig,
}

pub struct CriterionAnalysis {
    pub name: String,           // "coverage", "vaf", "expression", "type"
    pub threshold: String,      // e.g., "3" or "0.001"
    pub rejected_alone: usize,  // failed ONLY this criterion
    pub rejected_any: usize,    // failed this criterion (possibly with others)
    pub found_values: Vec<f64>, // metric values for found variants
    pub filtered_values: Vec<f64>, // metric values for filtered variants
}

pub struct NearMiss {
    pub variant: ExpectedVariant,
    pub criterion: String,
    pub actual_value: f64,
    pub threshold: f64,
    pub margin: f64,  // how close (ratio: actual/threshold)
}
```

## Acceptance Criteria

- [ ] `--report <path>` flag on filter and run subcommands
- [ ] JSON report with all sections (summary, per-criterion, near-misses, per-variant)
- [ ] Human-readable text report when path ends in .txt
- [ ] Near-miss detection with configurable margin (default 2x threshold)
- [ ] Threshold sensitivity sweep at 0.5x, 1x, 2x, 5x multipliers
- [ ] Per-criterion rejection counts (alone and combined)
- [ ] Correctly distinguishes "not detected" vs "detected but filtered"
- [ ] Report includes the filter configuration used for reproducibility
