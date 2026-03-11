# Ground Truth Benchmarking

## What It Does

Compares detection results against a known ground truth variant set, computing standard classification metrics (sensitivity, specificity, PPV, NPV, F1) with breakdowns by variant type and VAF range. This is separate from the expected variants used in filtering — ground truth represents what is truly present in the sample.

## Why It Matters

Researchers need to evaluate detection accuracy against validated truth sets:
- Synthetic spike-in experiments with known variants at known VAFs
- Cell-line dilution series with characterized mutations
- Clinical samples with orthogonal validation (e.g., ddPCR confirmed)
- Benchmarking against km/kmtools for regression testing

The ground truth is fundamentally different from the expected variants file:
- **Expected variants** = "look for these specific variants" (drives filtering)
- **Ground truth** = "these variants are truly present at these VAFs" (drives evaluation)

## CLI Interface

```bash
# New subcommand
kmerdet benchmark --results detect.tsv --truth ground_truth.tsv

# With VAF bins
kmerdet benchmark --results detect.tsv --truth ground_truth.tsv --vaf-bins "0,0.001,0.01,0.1,1.0"

# JSON output
kmerdet benchmark --results detect.tsv --truth ground_truth.tsv -o benchmark.json -f json

# Threshold sweep (ROC-like analysis)
kmerdet benchmark --results detect.tsv --truth ground_truth.tsv --sweep-vaf 0.0001,0.001,0.005,0.01,0.05
```

## Ground Truth File Format

TSV with header:
```
chrom	pos	ref	alt	type	true_vaf	category
chr17	7577120	C	T	Substitution	0.05	somatic
chr12	25398284	C	A	Substitution	0.01	somatic
chr13	28610183	A	AGTG	Insertion	0.005	somatic
chr17	7577120	C	T	Substitution	0.0	absent
```

- `true_vaf`: the known VAF (0.0 means the variant is confirmed absent — true negative)
- `category`: optional grouping (somatic, germline, artifact, absent)

## Report Contents

### Confusion Matrix
```
                 Detected    Not Detected
Present (truth)    TP            FN
Absent (truth)     FP            TN
```

### Summary Metrics
- Sensitivity (recall): TP / (TP + FN)
- Specificity: TN / (TN + FP)
- Precision (PPV): TP / (TP + FP)
- NPV: TN / (TN + FN)
- F1 score: 2 * (precision * recall) / (precision + recall)
- Accuracy: (TP + TN) / total
- Detection rate: TP / total_present

### Per-Type Breakdown
| Type | Present | TP | FN | Sensitivity |
|------|---------|----|----|-------------|
| Substitution | 15 | 12 | 3 | 80.0% |
| Insertion | 5 | 2 | 3 | 40.0% |
| Deletion | 3 | 1 | 2 | 33.3% |

### Per-VAF-Bin Breakdown
| VAF Range | Present | TP | FN | Sensitivity |
|-----------|---------|----|----|-------------|
| 0.1-1.0 | 5 | 5 | 0 | 100% |
| 0.01-0.1 | 8 | 6 | 2 | 75% |
| 0.001-0.01 | 7 | 2 | 5 | 28.6% |

### Per-Variant Detail
For each ground truth variant:
- Match status (TP/FP/FN/TN)
- If detected: measured rVAF vs true VAF, VAF error ratio
- If missed: closest detected variant (if any), partial match info

### Threshold Sweep (optional --sweep-vaf)
For each VAF threshold:
- Sensitivity and specificity at that threshold
- Enables ROC-like analysis without formal ROC curves

## Data Structures

```rust
pub struct BenchmarkReport {
    pub confusion: ConfusionMatrix,
    pub summary: BenchmarkSummary,
    pub per_type: Vec<TypeBreakdown>,
    pub per_vaf_bin: Vec<VafBinBreakdown>,
    pub per_variant: Vec<VariantBenchmark>,
    pub threshold_sweep: Option<Vec<ThresholdPoint>>,
}

pub struct ConfusionMatrix {
    pub tp: usize,
    pub fp: usize,
    pub fn_: usize,
    pub tn: usize,
}

pub struct BenchmarkSummary {
    pub sensitivity: f64,
    pub specificity: f64,
    pub precision: f64,
    pub npv: f64,
    pub f1: f64,
    pub accuracy: f64,
}

pub struct GroundTruthVariant {
    pub chrom: String,
    pub pos: u64,
    pub ref_allele: String,
    pub alt_allele: String,
    pub variant_type: String,
    pub true_vaf: f64,
    pub category: Option<String>,
}
```

## Acceptance Criteria

- [ ] `kmerdet benchmark` subcommand with --results and --truth
- [ ] Ground truth TSV parsing with true_vaf and optional category
- [ ] Confusion matrix computation (TP/FP/FN/TN)
- [ ] Standard metrics: sensitivity, specificity, PPV, NPV, F1, accuracy
- [ ] Per-variant-type breakdown
- [ ] Per-VAF-bin breakdown with configurable bins (--vaf-bins)
- [ ] Per-variant detail: match status, measured vs true VAF
- [ ] Threshold sweep with --sweep-vaf for ROC-like analysis
- [ ] Text output (human-readable) and JSON output
- [ ] Handles edge cases: no true positives, no true negatives, all same type
- [ ] Variant matching uses same logic as filter (reference mode + chr normalization)
