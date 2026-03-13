# kmerdet vs km: Performance Comparison

**Analysis date**: YYYY-MM-DD
**kmerdet version**: vX.Y.Z (git hash: XXXXXXX)
**km version**: X.X.X
**Dataset**: [dataset name and description]
**Analyst**: [name]

---

## Executive Summary

[2-4 sentences summarizing the key findings. Example:
"kmerdet v0.2.0 achieves 82% SNV sensitivity on the simulated_snv_indel
benchmark, a 5 percentage-point improvement over the km baseline (77%).
INDEL sensitivity improved from 38% to 51%, driven primarily by the
per-target k-mer length optimization. Runtime decreased from 120 s to 45 s
(2.7x speedup) on the 50-target dataset."]

---

## 1. Dataset Description

| Parameter | Value |
|-----------|-------|
| Dataset name | [e.g., simulated_snv_indel / patient_cohort_10] |
| N targets | |
| N variants (present) | |
| N variants (absent, specificity controls) | |
| Variant type mix | [e.g., 30 SNVs, 15 insertions, 15 deletions, 5 ITDs] |
| VAF range | [e.g., 0.1% – 50%] |
| Sequencing depth | [e.g., 3000x simulated / real clinical samples] |
| Sample type | [simulated / clinical cfDNA / cell line dilution] |

---

## 2. Detection Performance: Overall

| Metric | kmerdet vX.Y.Z | km (thesis baseline) | Delta | p-value (McNemar) |
|--------|---------------|---------------------|-------|-------------------|
| Sensitivity | | 0.77 (SNV), 0.38 (INDEL) | | |
| Specificity | | N/A | | |
| Precision (PPV) | | N/A | | |
| NPV | | N/A | | |
| F1-score | | N/A | | |
| MCC | | N/A | | |
| TP | | | | |
| FP | | | | |
| FN | | | | |
| TN | | | | |

**McNemar's test interpretation**: [e.g., "Significant improvement (p=0.023); kmerdet
detected 8 additional variants that km missed, while km detected 1 that kmerdet missed."]

---

## 3. Sensitivity by Variant Type

| Variant Type | kmerdet | km Baseline | Delta | vs Alignment-Based |
|-------------|---------|------------|-------|--------------------|
| Substitution (SNV) | | 77% | | 90-95% |
| Insertion | | 38% | | 80-95% |
| Deletion | | 38% | | 80-95% |
| ITD | | ~3-20% | | varies |
| Complex | | unknown | | varies |

**Figure 1**: `variant_type_breakdown.png` — Stacked bar chart of TP vs FN by type.

**Observations**:
- [Key finding 1: e.g., "SNV improvement driven by adaptive threshold on high-coverage targets"]
- [Key finding 2: e.g., "Large INDEL sensitivity remains low at X% due to k-mer length constraint"]
- [Key finding 3]

---

## 4. Sensitivity by VAF Bin

| VAF Range | kmerdet | km Baseline | Clinical Context |
|-----------|---------|------------|-----------------|
| 0–0.1% | | ~15% estimated | Below practical threshold |
| 0.1–1% | | ~60% estimated | Critical MRD monitoring range |
| 1–5% | | ~88% estimated | Elevated ctDNA |
| 5–100% | | ~98% estimated | High tumor burden |

**Figure 2**: `sensitivity_by_vaf.png` — Bar chart with thesis baseline overlaid.

**VAF detection threshold**: [estimated LOD, e.g., "kmerdet detects reliably above 0.08% VAF
on this dataset with the adaptive threshold feature enabled."]

**Observations**:
- [Key finding: e.g., "The 0.1–1% bin improved from 60% to 72%, most important for MRD"]
- [Key finding: e.g., "Ultra-low-VAF (<0.1%) detection remains at ~15% without duplex UMI support"]

---

## 5. rVAF Quantification Accuracy

| Metric | Value |
|--------|-------|
| Pearson r (rVAF vs true VAF) | |
| Mean absolute error (all variants) | |
| Mean absolute error (VAF > 0.1%) | |
| Median rVAF bias | |
| 95th percentile absolute error | |

**Figure 3**: Scatter plot of measured rVAF vs true VAF (not included in standard benchmark suite;
generate manually using per_variant.tsv).

**Observations**:
- [e.g., "rVAF estimates are well-calibrated above 0.1% VAF (r=0.97), but noisy below (r=0.61)"]

---

## 6. Runtime Performance

| Step | kmerdet vX.Y.Z | km (thesis) | Speedup |
|------|---------------|------------|---------|
| detect (50 targets, 4 threads) | | ~120 s | |
| filter | | ~15 s | |
| merge | | ~5 s | |
| total pipeline | | ~360 s (incl. HUMID+JF) | |

**Figure 4**: `runtime_scaling.png` — Runtime vs number of targets.

**Peak memory**: [e.g., "312 MB peak RSS vs not measured for km"]

**Thread scaling efficiency**:

| Threads | Runtime (s) | Speedup vs 1 thread | Efficiency |
|---------|-------------|---------------------|------------|
| 1 | | 1.0x | 100% |
| 2 | | | |
| 4 | | | |
| 8 | | | |

**Observations**:
- [e.g., "kmerdet achieves 3.8x speedup with 4 threads (95% efficiency)"]
- [e.g., "Thread scaling degrades above 8 threads due to graph construction bottleneck"]

---

## 7. Multi-k vs Single-k Comparison (if applicable)

*Only applicable when multi-k feature is enabled.*

| Variant Type | Single-k (k=31) | Multi-k (k=21+31) | Improvement |
|-------------|----------------|-------------------|-------------|
| Substitution | | | |
| Insertion (≤7 bp) | | | |
| Insertion (>7 bp) | | | |
| Deletion (≤7 bp) | | | |
| Deletion (>7 bp) | | | |
| ITD | | | |

**Figure 5**: `multi_k_comparison.png` — Grouped bar chart.

---

## 8. False Positive Analysis

| FP Category | Count | % of all FPs | Root Cause |
|-------------|-------|--------------|------------|
| Unmatched detected calls | | | No truth entry |
| Called in absent samples | | | Artifact or CHIP |
| rVAF near threshold | | | Borderline call |

**False positive rate**: [TP / (TP + FP)]

**Top false positive contexts** (if available):
- [e.g., "3/7 FPs occurred in homopolymer regions (≥5 A/T run)"]
- [e.g., "2/7 FPs had very low rVAF (<0.001) and borderline k-mer counts"]

---

## 9. Parameter Sensitivity Analysis

*From the accuracy sweep in run_benchmarks.sh.*

| count | ratio | SNV sens | INDEL sens | Precision | F1 |
|-------|-------|----------|------------|-----------|-----|
| 2 | 0.00001 | | | | |
| 2 | 0.0001 | | | | |
| 3 | 0.00001 | | | | |
| 5 | 0.00001 | | | | |
| 5 | 0.0001 | | | | |

**Recommended parameters**: [e.g., "count=2, ratio=0.0001 gives the best F1 on this dataset
while maintaining precision above 0.90."]

---

## 10. Statistical Testing

### McNemar's Test (kmerdet vs km)

This test evaluates whether the two methods differ significantly in which
specific variants they detect, accounting for the paired structure of the data
(same truth variants, different detectors).

```
                km detects    km misses
kmerdet detects      a             b
kmerdet misses       c             d

b = kmerdet only: X
c = km only:      Y
McNemar chi^2 = (|b - c| - 1)^2 / (b + c) = Z
p-value = W
```

**Interpretation**: [Copy from analyze_results.py output]

### Effect size

If McNemar's test is significant, report:
- Absolute improvement: [b - c] additional variants detected by kmerdet
- Relative improvement: [(b - c) / (a + c)] × 100%

---

## 11. Key Findings and Recommendations

### Confirmed improvements

- [List specific improvements verified by this analysis]

### Confirmed regressions

- [Any metrics that decreased — must be investigated before release]

### Areas still needing improvement

- [ ] INDEL sensitivity for >7 bp INDELs (currently X%, thesis ~3%)
- [ ] Ultra-low-VAF detection (<0.1% VAF)
- [ ] [Other specific gaps]

### Recommended parameter settings

| Use case | count | ratio | Notes |
|----------|-------|-------|-------|
| MRD screening (maximize sensitivity) | 2 | 0.00001 | May have FPs at 0.1-0.3% |
| Clinical reporting (balanced) | 2 | 0.0001 | Recommended |
| Confirmatory (maximize specificity) | 5 | 0.001 | Miss borderline calls |

---

## 12. Appendix: Run Metadata

```json
{
  "run_date": "YYYY-MM-DD",
  "kmerdet_version": "vX.Y.Z",
  "git_hash": "XXXXXXX",
  "dataset": "...",
  "n_targets": ...,
  "threads": 4,
  "vaf_bins": [0.0, 0.001, 0.01, 0.05, 0.1, 0.5, 1.0],
  "sweep_thresholds": [0.0, 0.0001, 0.001, 0.01, 0.05, 0.1],
  "files": {
    "detection_results": "...",
    "ground_truth": "...",
    "benchmark_json": "...",
    "metrics_json": "..."
  }
}
```

---

## 13. Figures

| Figure | File | Description |
|--------|------|-------------|
| 1 | `variant_type_breakdown.png` | TP vs FN stacked bar by variant type |
| 2 | `sensitivity_by_vaf.png` | Sensitivity by VAF bin with thesis baseline |
| 3 | `sensitivity_heatmap.png` | Heatmap: variant type × VAF bin |
| 4 | `threshold_sweep.png` | Sensitivity/precision vs rVAF threshold |
| 5 | `runtime_scaling.png` | Runtime scaling (if performance benchmark run) |
| 6 | `multi_k_comparison.png` | Multi-k vs single-k (if multi-k benchmark run) |

*Figures generated by: `python3 docs/benchmarking/framework/plot_results.py --results-dir <dir>`*
