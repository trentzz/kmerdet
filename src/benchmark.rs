//! Ground truth benchmarking for variant detection results.
//!
//! Compares detected [`VariantCall`]s against known ground truth variants
//! to compute sensitivity, specificity, precision, F1-score, and other
//! performance metrics. Supports breakdown by variant type, VAF bin, and
//! optional VAF threshold sweep for ROC-like analysis.

use std::fmt;
use std::io::Write;

use anyhow::Result;
use serde::Serialize;

use crate::variant::VariantCall;

// ---------------------------------------------------------------------------
// Ground truth input
// ---------------------------------------------------------------------------

/// A known ground truth variant with expected VAF.
#[derive(Debug, Clone)]
pub struct GroundTruthVariant {
    /// Chromosome name (e.g., "chr1" or "1").
    pub chrom: String,
    /// 1-based genomic position.
    pub pos: u64,
    /// Reference allele.
    pub ref_allele: String,
    /// Alternative allele.
    pub alt_allele: String,
    /// Variant type label (e.g., "Substitution", "Insertion").
    pub variant_type: String,
    /// True variant allele frequency (0.0 means expected absent).
    pub true_vaf: f64,
    /// Optional category label for grouping (e.g., "SNV", "tier1").
    pub category: Option<String>,
}

// ---------------------------------------------------------------------------
// Classification status
// ---------------------------------------------------------------------------

/// Classification outcome for a single variant.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
pub enum ClassificationStatus {
    /// True positive: present in truth and detected.
    TP,
    /// False positive: not in truth (or truth VAF == 0) but detected.
    FP,
    /// False negative: present in truth but not detected.
    FN,
    /// True negative: absent in truth and not detected.
    TN,
}

impl fmt::Display for ClassificationStatus {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::TP => write!(f, "TP"),
            Self::FP => write!(f, "FP"),
            Self::FN => write!(f, "FN"),
            Self::TN => write!(f, "TN"),
        }
    }
}

// ---------------------------------------------------------------------------
// Confusion matrix
// ---------------------------------------------------------------------------

/// Standard 2x2 confusion matrix counts.
#[derive(Debug, Clone, Copy, Default, Serialize)]
pub struct ConfusionMatrix {
    /// True positives.
    pub tp: usize,
    /// False positives.
    pub fp: usize,
    /// False negatives.
    pub fn_count: usize,
    /// True negatives.
    pub tn: usize,
}

impl ConfusionMatrix {
    /// Sensitivity (recall, true positive rate). TP / (TP + FN).
    pub fn sensitivity(&self) -> f64 {
        let denom = self.tp + self.fn_count;
        if denom == 0 {
            f64::NAN
        } else {
            self.tp as f64 / denom as f64
        }
    }

    /// Specificity (true negative rate). TN / (TN + FP).
    pub fn specificity(&self) -> f64 {
        let denom = self.tn + self.fp;
        if denom == 0 {
            f64::NAN
        } else {
            self.tn as f64 / denom as f64
        }
    }

    /// Precision (positive predictive value). TP / (TP + FP).
    pub fn precision(&self) -> f64 {
        let denom = self.tp + self.fp;
        if denom == 0 {
            f64::NAN
        } else {
            self.tp as f64 / denom as f64
        }
    }

    /// Negative predictive value. TN / (TN + FN).
    pub fn npv(&self) -> f64 {
        let denom = self.tn + self.fn_count;
        if denom == 0 {
            f64::NAN
        } else {
            self.tn as f64 / denom as f64
        }
    }

    /// F1-score: harmonic mean of precision and sensitivity.
    pub fn f1(&self) -> f64 {
        let p = self.precision();
        let s = self.sensitivity();
        if p.is_nan() || s.is_nan() || (p + s) == 0.0 {
            f64::NAN
        } else {
            2.0 * p * s / (p + s)
        }
    }

    /// Accuracy: (TP + TN) / (TP + TN + FP + FN).
    pub fn accuracy(&self) -> f64 {
        let total = self.tp + self.tn + self.fp + self.fn_count;
        if total == 0 {
            f64::NAN
        } else {
            (self.tp + self.tn) as f64 / total as f64
        }
    }
}

// ---------------------------------------------------------------------------
// Benchmark report structures
// ---------------------------------------------------------------------------

/// Top-level numeric summary of benchmark metrics.
#[derive(Debug, Clone, Serialize)]
pub struct BenchmarkSummary {
    pub sensitivity: f64,
    pub specificity: f64,
    pub precision: f64,
    pub npv: f64,
    pub f1: f64,
    pub accuracy: f64,
    pub total_truth: usize,
    pub total_present: usize,
    pub total_absent: usize,
    pub total_detected: usize,
}

/// Per-variant-type breakdown.
#[derive(Debug, Clone, Serialize)]
pub struct TypeBreakdown {
    pub variant_type: String,
    pub present: usize,
    pub tp: usize,
    pub fn_count: usize,
    pub sensitivity: f64,
}

/// Per-VAF-bin breakdown.
#[derive(Debug, Clone, Serialize)]
pub struct VafBinBreakdown {
    pub bin_label: String,
    pub bin_low: f64,
    pub bin_high: f64,
    pub present: usize,
    pub tp: usize,
    pub fn_count: usize,
    pub sensitivity: f64,
}

/// Per-variant detailed result.
#[derive(Debug, Clone, Serialize)]
pub struct VariantBenchmark {
    pub chrom: String,
    pub pos: u64,
    pub ref_allele: String,
    pub alt_allele: String,
    pub variant_type: String,
    pub true_vaf: f64,
    pub category: Option<String>,
    pub status: ClassificationStatus,
    pub measured_rvaf: Option<f64>,
    pub vaf_error: Option<f64>,
}

/// One point in a VAF-threshold sweep (for ROC-like analysis).
#[derive(Debug, Clone, Serialize)]
pub struct ThresholdPoint {
    pub vaf_threshold: f64,
    pub sensitivity: f64,
    pub specificity: f64,
    pub precision: f64,
    pub f1: f64,
}

/// Complete benchmark report.
#[derive(Debug, Clone, Serialize)]
pub struct BenchmarkReport {
    pub confusion: ConfusionMatrix,
    pub summary: BenchmarkSummary,
    pub per_type: Vec<TypeBreakdown>,
    pub per_vaf_bin: Vec<VafBinBreakdown>,
    pub per_variant: Vec<VariantBenchmark>,
    pub threshold_sweep: Option<Vec<ThresholdPoint>>,
}

// ---------------------------------------------------------------------------
// Chromosome normalization (same logic as filter/reference_mode.rs)
// ---------------------------------------------------------------------------

/// Strip the "chr" prefix for case-insensitive chromosome comparison.
fn normalize_chrom(chrom: &str) -> &str {
    chrom.strip_prefix("chr").unwrap_or(chrom)
}

// ---------------------------------------------------------------------------
// Matching logic
// ---------------------------------------------------------------------------

/// Check whether a detected call matches a ground truth variant by
/// chromosome, position, reference allele, and alternative allele.
///
/// Chromosome comparison uses `normalize_chrom` so "chr1" matches "1".
fn call_matches_truth(call: &VariantCall, truth: &GroundTruthVariant) -> bool {
    let chrom_match = call
        .chrom
        .as_deref()
        .map(|c| normalize_chrom(c) == normalize_chrom(&truth.chrom))
        .unwrap_or(false);

    let pos_match = call.pos == Some(truth.pos);
    let ref_match = call.ref_allele.as_deref() == Some(truth.ref_allele.as_str());
    let alt_match = call.alt_allele.as_deref() == Some(truth.alt_allele.as_str());

    chrom_match && pos_match && ref_match && alt_match
}

// ---------------------------------------------------------------------------
// Core benchmarking
// ---------------------------------------------------------------------------

/// Run a benchmark comparing detected variant calls against ground truth.
///
/// # Arguments
///
/// * `calls` - Detected variant calls from kmerdet.
/// * `truth` - Ground truth variants with known VAFs.
/// * `vaf_bins` - Ordered VAF bin boundaries (e.g., `[0.0, 0.001, 0.01, 0.1, 1.0]`).
///   Adjacent pairs define bins: `[0, 0.001)`, `[0.001, 0.01)`, etc.
/// * `sweep_thresholds` - Optional list of VAF thresholds for threshold sweep analysis.
///
/// # Returns
///
/// A [`BenchmarkReport`] with confusion matrix, summary metrics, and breakdowns.
pub fn run_benchmark(
    calls: &[VariantCall],
    truth: &[GroundTruthVariant],
    vaf_bins: &[f64],
    sweep_thresholds: Option<&[f64]>,
) -> BenchmarkReport {
    // Track which calls have been matched to a truth variant.
    let mut call_matched = vec![false; calls.len()];

    let mut per_variant = Vec::with_capacity(truth.len());
    let mut tp: usize = 0;
    let mut fn_count: usize = 0;
    let mut tn: usize = 0;

    // --- Classify each truth variant ---
    for tv in truth {
        let is_present = tv.true_vaf > 0.0;

        // Find the best matching call (prefer highest rVAF if multiple match).
        let best_match = calls
            .iter()
            .enumerate()
            .filter(|(_, c)| call_matches_truth(c, tv))
            .max_by(|(_, a), (_, b)| {
                a.rvaf
                    .partial_cmp(&b.rvaf)
                    .unwrap_or(std::cmp::Ordering::Equal)
            });

        let (status, measured_rvaf, vaf_error) = match best_match {
            Some((idx, call)) => {
                call_matched[idx] = true;
                let measured = call.rvaf;
                if is_present {
                    tp += 1;
                    let err = measured - tv.true_vaf;
                    (ClassificationStatus::TP, Some(measured), Some(err))
                } else {
                    // Truth says absent but we detected it: FP is counted
                    // below in the unmatched-call sweep. Here we record the
                    // truth entry as FP since the absent variant was called.
                    (ClassificationStatus::FP, Some(measured), None)
                }
            }
            None => {
                if is_present {
                    fn_count += 1;
                    (ClassificationStatus::FN, None, None)
                } else {
                    tn += 1;
                    (ClassificationStatus::TN, None, None)
                }
            }
        };

        per_variant.push(VariantBenchmark {
            chrom: tv.chrom.clone(),
            pos: tv.pos,
            ref_allele: tv.ref_allele.clone(),
            alt_allele: tv.alt_allele.clone(),
            variant_type: tv.variant_type.clone(),
            true_vaf: tv.true_vaf,
            category: tv.category.clone(),
            status,
            measured_rvaf,
            vaf_error,
        });
    }

    // --- Count FP from detected calls not matching any truth variant ---
    let mut fp: usize = 0;
    // Also count FPs from truth entries that were absent but detected
    // (already recorded in per_variant with status FP).
    let truth_fps = per_variant
        .iter()
        .filter(|v| v.status == ClassificationStatus::FP)
        .count();
    fp += truth_fps;

    for (idx, call) in calls.iter().enumerate() {
        if call_matched[idx] {
            continue;
        }
        // Skip Reference-type calls: they represent "no variant found"
        // and should not count as false positives.
        if call.variant_type == crate::variant::VariantType::Reference {
            continue;
        }
        fp += 1;
        per_variant.push(VariantBenchmark {
            chrom: call.chrom.clone().unwrap_or_default(),
            pos: call.pos.unwrap_or(0),
            ref_allele: call.ref_allele.clone().unwrap_or_default(),
            alt_allele: call.alt_allele.clone().unwrap_or_default(),
            variant_type: call.variant_type.to_string(),
            true_vaf: 0.0,
            category: None,
            status: ClassificationStatus::FP,
            measured_rvaf: Some(call.rvaf),
            vaf_error: None,
        });
    }

    let confusion = ConfusionMatrix {
        tp,
        fp,
        fn_count,
        tn,
    };

    let total_present = truth.iter().filter(|t| t.true_vaf > 0.0).count();
    let total_absent = truth.len() - total_present;
    let total_detected = calls
        .iter()
        .filter(|c| c.variant_type != crate::variant::VariantType::Reference)
        .count();

    let summary = BenchmarkSummary {
        sensitivity: confusion.sensitivity(),
        specificity: confusion.specificity(),
        precision: confusion.precision(),
        npv: confusion.npv(),
        f1: confusion.f1(),
        accuracy: confusion.accuracy(),
        total_truth: truth.len(),
        total_present,
        total_absent,
        total_detected,
    };

    // --- Per-type breakdown (truth variants with true_vaf > 0 only) ---
    let per_type = compute_type_breakdown(&per_variant);

    // --- Per-VAF-bin breakdown ---
    let per_vaf_bin = compute_vaf_bin_breakdown(&per_variant, vaf_bins);

    // --- Threshold sweep ---
    let threshold_sweep = sweep_thresholds.map(|thresholds| {
        compute_threshold_sweep(calls, truth, thresholds)
    });

    BenchmarkReport {
        confusion,
        summary,
        per_type,
        per_vaf_bin,
        per_variant,
        threshold_sweep,
    }
}

/// Compute per-variant-type breakdown from classified variants.
fn compute_type_breakdown(per_variant: &[VariantBenchmark]) -> Vec<TypeBreakdown> {
    use std::collections::BTreeMap;

    let mut by_type: BTreeMap<String, (usize, usize, usize)> = BTreeMap::new();

    for v in per_variant {
        // Only consider truth entries that are expected present (true_vaf > 0)
        // for the per-type sensitivity breakdown.
        if v.true_vaf <= 0.0 && v.status != ClassificationStatus::FP {
            continue;
        }
        // FP from unmatched calls don't belong to a truth type breakdown
        if v.true_vaf == 0.0 && v.status == ClassificationStatus::FP {
            continue;
        }
        let entry = by_type.entry(v.variant_type.clone()).or_insert((0, 0, 0));
        match v.status {
            ClassificationStatus::TP => {
                entry.0 += 1; // present
                entry.1 += 1; // tp
            }
            ClassificationStatus::FN => {
                entry.0 += 1; // present
                entry.2 += 1; // fn
            }
            _ => {}
        }
    }

    by_type
        .into_iter()
        .map(|(variant_type, (present, tp, fn_count))| {
            let sensitivity = if present == 0 {
                f64::NAN
            } else {
                tp as f64 / present as f64
            };
            TypeBreakdown {
                variant_type,
                present,
                tp,
                fn_count,
                sensitivity,
            }
        })
        .collect()
}

/// Compute per-VAF-bin breakdown from classified variants.
fn compute_vaf_bin_breakdown(
    per_variant: &[VariantBenchmark],
    vaf_bins: &[f64],
) -> Vec<VafBinBreakdown> {
    if vaf_bins.len() < 2 {
        return Vec::new();
    }

    let mut bins = Vec::with_capacity(vaf_bins.len() - 1);

    for window in vaf_bins.windows(2) {
        let low = window[0];
        let high = window[1];
        let bin_label = format!("[{}, {})", format_vaf(low), format_vaf(high));

        let mut present: usize = 0;
        let mut tp: usize = 0;
        let mut fn_count: usize = 0;

        for v in per_variant {
            // Only truth entries with positive true_vaf
            if v.true_vaf <= 0.0 {
                continue;
            }
            // Check if this variant falls in the current bin.
            // Last bin is inclusive on the right: [low, high].
            let in_bin = if high == vaf_bins[vaf_bins.len() - 1] {
                v.true_vaf >= low && v.true_vaf <= high
            } else {
                v.true_vaf >= low && v.true_vaf < high
            };
            if !in_bin {
                continue;
            }
            present += 1;
            match v.status {
                ClassificationStatus::TP => tp += 1,
                ClassificationStatus::FN => fn_count += 1,
                _ => {}
            }
        }

        let sensitivity = if present == 0 {
            f64::NAN
        } else {
            tp as f64 / present as f64
        };

        bins.push(VafBinBreakdown {
            bin_label,
            bin_low: low,
            bin_high: high,
            present,
            tp,
            fn_count,
            sensitivity,
        });
    }

    bins
}

/// Format a VAF value for display in bin labels.
fn format_vaf(v: f64) -> String {
    if v == 0.0 {
        "0".to_string()
    } else if v >= 0.01 {
        format!("{:.2}", v)
    } else {
        format!("{:.4}", v)
    }
}

/// Compute metrics at each VAF threshold for a threshold sweep.
///
/// At each threshold, calls with `rvaf < threshold` are excluded from the
/// detected set, and the confusion matrix is recomputed.
fn compute_threshold_sweep(
    calls: &[VariantCall],
    truth: &[GroundTruthVariant],
    thresholds: &[f64],
) -> Vec<ThresholdPoint> {
    thresholds
        .iter()
        .map(|&threshold| {
            let filtered_calls: Vec<&VariantCall> = calls
                .iter()
                .filter(|c| {
                    c.variant_type != crate::variant::VariantType::Reference
                        && c.rvaf >= threshold
                })
                .collect();

            let mut tp: usize = 0;
            let mut fn_count: usize = 0;
            let mut tn: usize = 0;
            let mut fp_from_truth: usize = 0;

            let mut call_used = vec![false; filtered_calls.len()];

            for tv in truth {
                let is_present = tv.true_vaf > 0.0;

                let best = filtered_calls
                    .iter()
                    .enumerate()
                    .filter(|(_, c)| call_matches_truth(c, tv))
                    .max_by(|(_, a), (_, b)| {
                        a.rvaf
                            .partial_cmp(&b.rvaf)
                            .unwrap_or(std::cmp::Ordering::Equal)
                    });

                match best {
                    Some((idx, _)) => {
                        call_used[idx] = true;
                        if is_present {
                            tp += 1;
                        } else {
                            fp_from_truth += 1;
                        }
                    }
                    None => {
                        if is_present {
                            fn_count += 1;
                        } else {
                            tn += 1;
                        }
                    }
                }
            }

            let fp_unmatched = call_used.iter().filter(|&&used| !used).count();
            let fp = fp_from_truth + fp_unmatched;

            let cm = ConfusionMatrix {
                tp,
                fp,
                fn_count,
                tn,
            };

            ThresholdPoint {
                vaf_threshold: threshold,
                sensitivity: cm.sensitivity(),
                specificity: cm.specificity(),
                precision: cm.precision(),
                f1: cm.f1(),
            }
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Output: human-readable text
// ---------------------------------------------------------------------------

/// Write benchmark results as human-readable text.
pub fn write_benchmark_text(report: &BenchmarkReport, writer: &mut dyn Write) -> Result<()> {
    writeln!(writer, "=== kmerdet Benchmark Report ===")?;
    writeln!(writer)?;

    // Confusion matrix
    writeln!(writer, "Confusion Matrix:")?;
    writeln!(writer, "  TP: {:>6}    FP: {:>6}", report.confusion.tp, report.confusion.fp)?;
    writeln!(
        writer,
        "  FN: {:>6}    TN: {:>6}",
        report.confusion.fn_count, report.confusion.tn
    )?;
    writeln!(writer)?;

    // Summary metrics
    writeln!(writer, "Summary Metrics:")?;
    writeln!(
        writer,
        "  Sensitivity (recall):  {}",
        format_metric(report.summary.sensitivity)
    )?;
    writeln!(
        writer,
        "  Specificity:           {}",
        format_metric(report.summary.specificity)
    )?;
    writeln!(
        writer,
        "  Precision (PPV):       {}",
        format_metric(report.summary.precision)
    )?;
    writeln!(
        writer,
        "  NPV:                   {}",
        format_metric(report.summary.npv)
    )?;
    writeln!(
        writer,
        "  F1-score:              {}",
        format_metric(report.summary.f1)
    )?;
    writeln!(
        writer,
        "  Accuracy:              {}",
        format_metric(report.summary.accuracy)
    )?;
    writeln!(writer)?;
    writeln!(
        writer,
        "  Truth variants:        {}",
        report.summary.total_truth
    )?;
    writeln!(
        writer,
        "  Expected present:      {}",
        report.summary.total_present
    )?;
    writeln!(
        writer,
        "  Expected absent:       {}",
        report.summary.total_absent
    )?;
    writeln!(
        writer,
        "  Total detected:        {}",
        report.summary.total_detected
    )?;
    writeln!(writer)?;

    // Per-type breakdown
    if !report.per_type.is_empty() {
        writeln!(writer, "Per-Type Breakdown:")?;
        writeln!(
            writer,
            "  {:15} {:>8} {:>6} {:>6} {:>12}",
            "Type", "Present", "TP", "FN", "Sensitivity"
        )?;
        writeln!(writer, "  {}", "-".repeat(53))?;
        for t in &report.per_type {
            writeln!(
                writer,
                "  {:15} {:>8} {:>6} {:>6} {:>12}",
                t.variant_type,
                t.present,
                t.tp,
                t.fn_count,
                format_metric(t.sensitivity)
            )?;
        }
        writeln!(writer)?;
    }

    // Per-VAF-bin breakdown
    if !report.per_vaf_bin.is_empty() {
        writeln!(writer, "Per-VAF-Bin Breakdown:")?;
        writeln!(
            writer,
            "  {:20} {:>8} {:>6} {:>6} {:>12}",
            "Bin", "Present", "TP", "FN", "Sensitivity"
        )?;
        writeln!(writer, "  {}", "-".repeat(58))?;
        for b in &report.per_vaf_bin {
            writeln!(
                writer,
                "  {:20} {:>8} {:>6} {:>6} {:>12}",
                b.bin_label,
                b.present,
                b.tp,
                b.fn_count,
                format_metric(b.sensitivity)
            )?;
        }
        writeln!(writer)?;
    }

    // Per-variant detail
    writeln!(writer, "Per-Variant Detail:")?;
    writeln!(
        writer,
        "  {:8} {:>10} {:6} {:6} {:15} {:>10} {:>4} {:>12} {:>12}",
        "Chrom", "Pos", "Ref", "Alt", "Type", "TrueVAF", "Stat", "MeasuredVAF", "VAF Error"
    )?;
    writeln!(writer, "  {}", "-".repeat(95))?;
    for v in &report.per_variant {
        writeln!(
            writer,
            "  {:8} {:>10} {:6} {:6} {:15} {:>10.6} {:>4} {:>12} {:>12}",
            v.chrom,
            v.pos,
            v.ref_allele,
            v.alt_allele,
            v.variant_type,
            v.true_vaf,
            v.status,
            v.measured_rvaf
                .map(|r| format!("{:.6}", r))
                .unwrap_or_else(|| "-".to_string()),
            v.vaf_error
                .map(|e| format!("{:+.6}", e))
                .unwrap_or_else(|| "-".to_string()),
        )?;
    }
    writeln!(writer)?;

    // Threshold sweep
    if let Some(sweep) = &report.threshold_sweep {
        writeln!(writer, "Threshold Sweep:")?;
        writeln!(
            writer,
            "  {:>12} {:>12} {:>12} {:>12} {:>12}",
            "Threshold", "Sensitivity", "Specificity", "Precision", "F1"
        )?;
        writeln!(writer, "  {}", "-".repeat(64))?;
        for pt in sweep {
            writeln!(
                writer,
                "  {:>12.6} {:>12} {:>12} {:>12} {:>12}",
                pt.vaf_threshold,
                format_metric(pt.sensitivity),
                format_metric(pt.specificity),
                format_metric(pt.precision),
                format_metric(pt.f1),
            )?;
        }
        writeln!(writer)?;
    }

    Ok(())
}

/// Format a metric value, showing "N/A" for NaN.
fn format_metric(v: f64) -> String {
    if v.is_nan() {
        "N/A".to_string()
    } else {
        format!("{:.4}", v)
    }
}

// ---------------------------------------------------------------------------
// Output: JSON
// ---------------------------------------------------------------------------

/// Write benchmark results as JSON.
pub fn write_benchmark_json(report: &BenchmarkReport, writer: &mut dyn Write) -> Result<()> {
    serde_json::to_writer_pretty(&mut *writer, report)?;
    writeln!(writer)?;
    Ok(())
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::variant::VariantType;

    /// Helper to build a VariantCall for testing.
    fn make_call(
        chrom: &str,
        pos: u64,
        ref_allele: &str,
        alt_allele: &str,
        vtype: VariantType,
        rvaf: f64,
    ) -> VariantCall {
        VariantCall {
            sample: "sample1".to_string(),
            target: "target1".to_string(),
            variant_type: vtype,
            variant_name: "test".to_string(),
            rvaf,
            expression: 1.0,
            min_coverage: 10,
            start_kmer_count: 20,
            ref_sequence: "ACGT".to_string(),
            alt_sequence: "TCGT".to_string(),
            info: "vs_ref".to_string(),
            chrom: Some(chrom.to_string()),
            pos: Some(pos),
            ref_allele: Some(ref_allele.to_string()),
            alt_allele: Some(alt_allele.to_string()),
            pvalue: None,
            qual: None,
            ci_lower: None,
            ci_upper: None,
        }
    }

    /// Helper to build a GroundTruthVariant.
    fn make_truth(
        chrom: &str,
        pos: u64,
        ref_allele: &str,
        alt_allele: &str,
        variant_type: &str,
        true_vaf: f64,
    ) -> GroundTruthVariant {
        GroundTruthVariant {
            chrom: chrom.to_string(),
            pos,
            ref_allele: ref_allele.to_string(),
            alt_allele: alt_allele.to_string(),
            variant_type: variant_type.to_string(),
            true_vaf,
            category: None,
        }
    }

    // -----------------------------------------------------------------------
    // ConfusionMatrix unit tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_confusion_matrix_perfect() {
        let cm = ConfusionMatrix {
            tp: 10,
            fp: 0,
            fn_count: 0,
            tn: 5,
        };
        assert_eq!(cm.sensitivity(), 1.0);
        assert_eq!(cm.specificity(), 1.0);
        assert_eq!(cm.precision(), 1.0);
        assert_eq!(cm.npv(), 1.0);
        assert_eq!(cm.f1(), 1.0);
        assert_eq!(cm.accuracy(), 1.0);
    }

    #[test]
    fn test_confusion_matrix_no_positives() {
        let cm = ConfusionMatrix {
            tp: 0,
            fp: 0,
            fn_count: 0,
            tn: 10,
        };
        assert!(cm.sensitivity().is_nan());
        assert_eq!(cm.specificity(), 1.0);
        assert!(cm.precision().is_nan());
        assert_eq!(cm.npv(), 1.0);
    }

    #[test]
    fn test_confusion_matrix_no_negatives() {
        let cm = ConfusionMatrix {
            tp: 5,
            fp: 3,
            fn_count: 2,
            tn: 0,
        };
        assert_eq!(cm.sensitivity(), 5.0 / 7.0);
        // TN=0, FP=3: specificity = 0 / (0 + 3) = 0.0
        assert_eq!(cm.specificity(), 0.0);
        assert_eq!(cm.precision(), 5.0 / 8.0);
    }

    #[test]
    fn test_confusion_matrix_no_negatives_at_all() {
        // When both TN=0 and FP=0, specificity is undefined
        let cm = ConfusionMatrix {
            tp: 5,
            fp: 0,
            fn_count: 2,
            tn: 0,
        };
        assert!(cm.specificity().is_nan());
    }

    #[test]
    fn test_confusion_matrix_empty() {
        let cm = ConfusionMatrix::default();
        assert!(cm.sensitivity().is_nan());
        assert!(cm.specificity().is_nan());
        assert!(cm.precision().is_nan());
        assert!(cm.npv().is_nan());
        assert!(cm.f1().is_nan());
        assert!(cm.accuracy().is_nan());
    }

    #[test]
    fn test_confusion_matrix_known_values() {
        // Classic example: 50 TP, 10 FP, 5 FN, 100 TN
        let cm = ConfusionMatrix {
            tp: 50,
            fp: 10,
            fn_count: 5,
            tn: 100,
        };
        let sens = 50.0 / 55.0;
        let spec = 100.0 / 110.0;
        let prec = 50.0 / 60.0;
        let npv = 100.0 / 105.0;
        let f1 = 2.0 * prec * sens / (prec + sens);
        let acc = 150.0 / 165.0;

        assert!((cm.sensitivity() - sens).abs() < 1e-10);
        assert!((cm.specificity() - spec).abs() < 1e-10);
        assert!((cm.precision() - prec).abs() < 1e-10);
        assert!((cm.npv() - npv).abs() < 1e-10);
        assert!((cm.f1() - f1).abs() < 1e-10);
        assert!((cm.accuracy() - acc).abs() < 1e-10);
    }

    // -----------------------------------------------------------------------
    // Matching tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_call_matches_truth_exact() {
        let call = make_call("chr1", 100, "A", "T", VariantType::Substitution, 0.05);
        let truth = make_truth("chr1", 100, "A", "T", "Substitution", 0.05);
        assert!(call_matches_truth(&call, &truth));
    }

    #[test]
    fn test_call_matches_truth_chr_normalization() {
        let call = make_call("1", 100, "A", "T", VariantType::Substitution, 0.05);
        let truth = make_truth("chr1", 100, "A", "T", "Substitution", 0.05);
        assert!(call_matches_truth(&call, &truth));

        // Reverse direction
        let call2 = make_call("chr1", 100, "A", "T", VariantType::Substitution, 0.05);
        let truth2 = make_truth("1", 100, "A", "T", "Substitution", 0.05);
        assert!(call_matches_truth(&call2, &truth2));
    }

    #[test]
    fn test_call_matches_truth_mismatch_pos() {
        let call = make_call("chr1", 101, "A", "T", VariantType::Substitution, 0.05);
        let truth = make_truth("chr1", 100, "A", "T", "Substitution", 0.05);
        assert!(!call_matches_truth(&call, &truth));
    }

    #[test]
    fn test_call_matches_truth_mismatch_allele() {
        let call = make_call("chr1", 100, "A", "C", VariantType::Substitution, 0.05);
        let truth = make_truth("chr1", 100, "A", "T", "Substitution", 0.05);
        assert!(!call_matches_truth(&call, &truth));
    }

    #[test]
    fn test_call_matches_truth_missing_chrom() {
        let mut call = make_call("chr1", 100, "A", "T", VariantType::Substitution, 0.05);
        call.chrom = None;
        let truth = make_truth("chr1", 100, "A", "T", "Substitution", 0.05);
        assert!(!call_matches_truth(&call, &truth));
    }

    // -----------------------------------------------------------------------
    // run_benchmark integration tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_benchmark_all_tp() {
        let calls = vec![
            make_call("chr1", 100, "A", "T", VariantType::Substitution, 0.05),
            make_call("chr2", 200, "G", "GACC", VariantType::Insertion, 0.10),
        ];
        let truth = vec![
            make_truth("chr1", 100, "A", "T", "Substitution", 0.05),
            make_truth("chr2", 200, "G", "GACC", "Insertion", 0.10),
        ];

        let report = run_benchmark(&calls, &truth, &[0.0, 0.1, 1.0], None);

        assert_eq!(report.confusion.tp, 2);
        assert_eq!(report.confusion.fp, 0);
        assert_eq!(report.confusion.fn_count, 0);
        assert_eq!(report.confusion.tn, 0);
        assert_eq!(report.summary.sensitivity, 1.0);
        assert_eq!(report.summary.precision, 1.0);
        assert_eq!(report.summary.f1, 1.0);
    }

    #[test]
    fn test_benchmark_all_fn() {
        let calls: Vec<VariantCall> = Vec::new();
        let truth = vec![
            make_truth("chr1", 100, "A", "T", "Substitution", 0.05),
            make_truth("chr2", 200, "G", "GACC", "Insertion", 0.10),
        ];

        let report = run_benchmark(&calls, &truth, &[0.0, 1.0], None);

        assert_eq!(report.confusion.tp, 0);
        assert_eq!(report.confusion.fn_count, 2);
        assert_eq!(report.summary.sensitivity, 0.0);
    }

    #[test]
    fn test_benchmark_with_fp_unmatched() {
        // One call that matches truth, one call that does not
        let calls = vec![
            make_call("chr1", 100, "A", "T", VariantType::Substitution, 0.05),
            make_call("chr3", 300, "C", "G", VariantType::Substitution, 0.02),
        ];
        let truth = vec![make_truth("chr1", 100, "A", "T", "Substitution", 0.05)];

        let report = run_benchmark(&calls, &truth, &[0.0, 1.0], None);

        assert_eq!(report.confusion.tp, 1);
        assert_eq!(report.confusion.fp, 1);
        assert_eq!(report.confusion.fn_count, 0);
        assert_eq!(report.summary.precision, 0.5);
    }

    #[test]
    fn test_benchmark_with_tn() {
        // Truth says variant is absent (true_vaf=0), and we don't detect it
        let calls = vec![make_call(
            "chr1",
            100,
            "A",
            "T",
            VariantType::Substitution,
            0.05,
        )];
        let truth = vec![
            make_truth("chr1", 100, "A", "T", "Substitution", 0.05),
            make_truth("chr2", 200, "G", "C", "Substitution", 0.0), // absent
        ];

        let report = run_benchmark(&calls, &truth, &[0.0, 1.0], None);

        assert_eq!(report.confusion.tp, 1);
        assert_eq!(report.confusion.tn, 1);
        assert_eq!(report.confusion.fp, 0);
        assert_eq!(report.confusion.fn_count, 0);
    }

    #[test]
    fn test_benchmark_absent_but_detected() {
        // Truth says absent (true_vaf=0), but we detect it => FP
        let calls = vec![make_call(
            "chr2",
            200,
            "G",
            "C",
            VariantType::Substitution,
            0.03,
        )];
        let truth = vec![make_truth("chr2", 200, "G", "C", "Substitution", 0.0)];

        let report = run_benchmark(&calls, &truth, &[0.0, 1.0], None);

        assert_eq!(report.confusion.tp, 0);
        assert_eq!(report.confusion.fp, 1);
        assert_eq!(report.confusion.fn_count, 0);
        assert_eq!(report.confusion.tn, 0);
    }

    #[test]
    fn test_benchmark_reference_calls_not_fp() {
        // Reference-type calls should not count as FP
        let calls = vec![
            make_call("chr1", 100, "A", "T", VariantType::Reference, 0.0),
        ];
        let truth: Vec<GroundTruthVariant> = Vec::new();

        let report = run_benchmark(&calls, &truth, &[0.0, 1.0], None);

        assert_eq!(report.confusion.fp, 0);
    }

    #[test]
    fn test_benchmark_vaf_error_computed() {
        let calls = vec![make_call(
            "chr1",
            100,
            "A",
            "T",
            VariantType::Substitution,
            0.048,
        )];
        let truth = vec![make_truth("chr1", 100, "A", "T", "Substitution", 0.05)];

        let report = run_benchmark(&calls, &truth, &[0.0, 1.0], None);

        let v = &report.per_variant[0];
        assert_eq!(v.status, ClassificationStatus::TP);
        assert!((v.measured_rvaf.unwrap() - 0.048).abs() < 1e-10);
        assert!((v.vaf_error.unwrap() - (-0.002)).abs() < 1e-10);
    }

    // -----------------------------------------------------------------------
    // Per-type breakdown tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_per_type_breakdown() {
        let calls = vec![
            make_call("chr1", 100, "A", "T", VariantType::Substitution, 0.05),
            make_call("chr2", 200, "G", "GACC", VariantType::Insertion, 0.10),
        ];
        let truth = vec![
            make_truth("chr1", 100, "A", "T", "Substitution", 0.05),
            make_truth("chr2", 200, "G", "GACC", "Insertion", 0.10),
            make_truth("chr3", 300, "C", "T", "Substitution", 0.03), // FN
        ];

        let report = run_benchmark(&calls, &truth, &[0.0, 1.0], None);

        let sub_type = report
            .per_type
            .iter()
            .find(|t| t.variant_type == "Substitution")
            .unwrap();
        assert_eq!(sub_type.present, 2);
        assert_eq!(sub_type.tp, 1);
        assert_eq!(sub_type.fn_count, 1);
        assert_eq!(sub_type.sensitivity, 0.5);

        let ins_type = report
            .per_type
            .iter()
            .find(|t| t.variant_type == "Insertion")
            .unwrap();
        assert_eq!(ins_type.present, 1);
        assert_eq!(ins_type.tp, 1);
        assert_eq!(ins_type.fn_count, 0);
        assert_eq!(ins_type.sensitivity, 1.0);
    }

    // -----------------------------------------------------------------------
    // Per-VAF-bin breakdown tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_per_vaf_bin_breakdown() {
        let calls = vec![
            make_call("chr1", 100, "A", "T", VariantType::Substitution, 0.05),
        ];
        let truth = vec![
            make_truth("chr1", 100, "A", "T", "Substitution", 0.05),
            make_truth("chr2", 200, "G", "C", "Substitution", 0.005), // FN, bin [0.001, 0.01)
        ];

        let bins = vec![0.0, 0.001, 0.01, 0.1, 1.0];
        let report = run_benchmark(&calls, &truth, &bins, None);

        // The 0.05 variant falls in [0.01, 0.1)
        let bin_01 = report
            .per_vaf_bin
            .iter()
            .find(|b| b.bin_low == 0.01 && b.bin_high == 0.1)
            .unwrap();
        assert_eq!(bin_01.present, 1);
        assert_eq!(bin_01.tp, 1);

        // The 0.005 variant falls in [0.001, 0.01)
        let bin_001 = report
            .per_vaf_bin
            .iter()
            .find(|b| b.bin_low == 0.001 && b.bin_high == 0.01)
            .unwrap();
        assert_eq!(bin_001.present, 1);
        assert_eq!(bin_001.fn_count, 1);
        assert_eq!(bin_001.sensitivity, 0.0);
    }

    #[test]
    fn test_vaf_bin_last_bin_inclusive() {
        // A variant with true_vaf = 1.0 should fall in the last bin [0.1, 1.0]
        let calls = vec![make_call(
            "chr1",
            100,
            "A",
            "T",
            VariantType::Substitution,
            0.95,
        )];
        let truth = vec![make_truth("chr1", 100, "A", "T", "Substitution", 1.0)];

        let bins = vec![0.0, 0.1, 1.0];
        let report = run_benchmark(&calls, &truth, &bins, None);

        let last_bin = report
            .per_vaf_bin
            .iter()
            .find(|b| b.bin_high == 1.0)
            .unwrap();
        assert_eq!(last_bin.present, 1);
        assert_eq!(last_bin.tp, 1);
    }

    #[test]
    fn test_vaf_bin_too_few_boundaries() {
        let report = run_benchmark(&[], &[], &[0.5], None);
        assert!(report.per_vaf_bin.is_empty());
    }

    // -----------------------------------------------------------------------
    // Threshold sweep tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_threshold_sweep() {
        let calls = vec![
            make_call("chr1", 100, "A", "T", VariantType::Substitution, 0.05),
            make_call("chr2", 200, "G", "C", VariantType::Substitution, 0.001),
        ];
        let truth = vec![
            make_truth("chr1", 100, "A", "T", "Substitution", 0.05),
            make_truth("chr2", 200, "G", "C", "Substitution", 0.001),
        ];

        let thresholds = vec![0.0, 0.01, 0.1];
        let report = run_benchmark(&calls, &truth, &[0.0, 1.0], Some(&thresholds));

        let sweep = report.threshold_sweep.as_ref().unwrap();
        assert_eq!(sweep.len(), 3);

        // At threshold 0.0: both detected => sens=1.0
        assert_eq!(sweep[0].sensitivity, 1.0);

        // At threshold 0.01: only chr1 (0.05) detected => sens=0.5
        assert_eq!(sweep[1].sensitivity, 0.5);

        // At threshold 0.1: neither detected => sens=0.0
        assert_eq!(sweep[2].sensitivity, 0.0);
    }

    #[test]
    fn test_threshold_sweep_none() {
        let report = run_benchmark(&[], &[], &[0.0, 1.0], None);
        assert!(report.threshold_sweep.is_none());
    }

    // -----------------------------------------------------------------------
    // Output tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_write_benchmark_text_output() {
        let calls = vec![
            make_call("chr1", 100, "A", "T", VariantType::Substitution, 0.05),
        ];
        let truth = vec![
            make_truth("chr1", 100, "A", "T", "Substitution", 0.05),
            make_truth("chr2", 200, "G", "C", "Substitution", 0.0),
        ];

        let report = run_benchmark(&calls, &truth, &[0.0, 0.1, 1.0], None);

        let mut buf = Vec::new();
        write_benchmark_text(&report, &mut buf).unwrap();
        let output = String::from_utf8(buf).unwrap();

        assert!(output.contains("kmerdet Benchmark Report"));
        assert!(output.contains("Confusion Matrix:"));
        assert!(output.contains("TP:"));
        assert!(output.contains("Sensitivity (recall):"));
        assert!(output.contains("Per-Type Breakdown:"));
        assert!(output.contains("Per-VAF-Bin Breakdown:"));
        assert!(output.contains("Per-Variant Detail:"));
    }

    #[test]
    fn test_write_benchmark_json_output() {
        let calls = vec![make_call(
            "chr1",
            100,
            "A",
            "T",
            VariantType::Substitution,
            0.05,
        )];
        let truth = vec![make_truth("chr1", 100, "A", "T", "Substitution", 0.05)];

        let report = run_benchmark(&calls, &truth, &[0.0, 1.0], None);

        let mut buf = Vec::new();
        write_benchmark_json(&report, &mut buf).unwrap();
        let output = String::from_utf8(buf).unwrap();

        // Verify valid JSON
        let parsed: serde_json::Value = serde_json::from_str(&output).unwrap();
        assert!(parsed.get("confusion").is_some());
        assert!(parsed.get("summary").is_some());
        assert!(parsed.get("per_type").is_some());
        assert!(parsed.get("per_vaf_bin").is_some());
        assert!(parsed.get("per_variant").is_some());
        assert_eq!(parsed["confusion"]["tp"], 1);
    }

    #[test]
    fn test_write_benchmark_json_with_sweep() {
        let calls = vec![make_call(
            "chr1",
            100,
            "A",
            "T",
            VariantType::Substitution,
            0.05,
        )];
        let truth = vec![make_truth("chr1", 100, "A", "T", "Substitution", 0.05)];
        let thresholds = vec![0.0, 0.01];

        let report = run_benchmark(&calls, &truth, &[0.0, 1.0], Some(&thresholds));

        let mut buf = Vec::new();
        write_benchmark_json(&report, &mut buf).unwrap();
        let output = String::from_utf8(buf).unwrap();

        let parsed: serde_json::Value = serde_json::from_str(&output).unwrap();
        assert!(parsed["threshold_sweep"].is_array());
        assert_eq!(parsed["threshold_sweep"].as_array().unwrap().len(), 2);
    }

    // -----------------------------------------------------------------------
    // Edge case / mixed scenario tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_benchmark_empty_inputs() {
        let report = run_benchmark(&[], &[], &[0.0, 1.0], None);

        assert_eq!(report.confusion.tp, 0);
        assert_eq!(report.confusion.fp, 0);
        assert_eq!(report.confusion.fn_count, 0);
        assert_eq!(report.confusion.tn, 0);
        assert!(report.per_variant.is_empty());
        assert!(report.per_type.is_empty());
        assert_eq!(report.summary.total_truth, 0);
    }

    #[test]
    fn test_benchmark_multiple_calls_same_truth() {
        // Two calls match the same truth; the one with higher rVAF should be used
        let calls = vec![
            make_call("chr1", 100, "A", "T", VariantType::Substitution, 0.03),
            make_call("chr1", 100, "A", "T", VariantType::Substitution, 0.07),
        ];
        let truth = vec![make_truth("chr1", 100, "A", "T", "Substitution", 0.05)];

        let report = run_benchmark(&calls, &truth, &[0.0, 1.0], None);

        assert_eq!(report.confusion.tp, 1);
        // The other unmatched call is FP
        assert_eq!(report.confusion.fp, 1);

        let tp_variant = report
            .per_variant
            .iter()
            .find(|v| v.status == ClassificationStatus::TP)
            .unwrap();
        // Should have picked the higher-rVAF call (0.07)
        assert!((tp_variant.measured_rvaf.unwrap() - 0.07).abs() < 1e-10);
    }

    #[test]
    fn test_benchmark_mixed_scenario() {
        // Complex scenario with multiple types and outcomes
        let calls = vec![
            make_call("chr1", 100, "A", "T", VariantType::Substitution, 0.05),
            make_call("chr2", 200, "G", "GACC", VariantType::Insertion, 0.10),
            make_call("chr5", 500, "T", "A", VariantType::Substitution, 0.01), // FP (not in truth)
        ];
        let truth = vec![
            make_truth("chr1", 100, "A", "T", "Substitution", 0.05),  // TP
            make_truth("chr2", 200, "G", "GACC", "Insertion", 0.10),   // TP
            make_truth("chr3", 300, "C", "T", "Deletion", 0.02),       // FN
            make_truth("chr4", 400, "A", "G", "Substitution", 0.0),    // TN (absent)
        ];

        let report = run_benchmark(&calls, &truth, &[0.0, 0.01, 0.1, 1.0], None);

        assert_eq!(report.confusion.tp, 2);
        assert_eq!(report.confusion.fp, 1);
        assert_eq!(report.confusion.fn_count, 1);
        assert_eq!(report.confusion.tn, 1);
        assert_eq!(report.summary.total_truth, 4);
        assert_eq!(report.summary.total_present, 3);
        assert_eq!(report.summary.total_absent, 1);

        // Sensitivity: 2/3
        assert!((report.summary.sensitivity - 2.0 / 3.0).abs() < 1e-10);
        // Precision: 2/3
        assert!((report.summary.precision - 2.0 / 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_classification_status_display() {
        assert_eq!(ClassificationStatus::TP.to_string(), "TP");
        assert_eq!(ClassificationStatus::FP.to_string(), "FP");
        assert_eq!(ClassificationStatus::FN.to_string(), "FN");
        assert_eq!(ClassificationStatus::TN.to_string(), "TN");
    }

    #[test]
    fn test_format_vaf() {
        assert_eq!(format_vaf(0.0), "0");
        assert_eq!(format_vaf(0.1), "0.10");
        assert_eq!(format_vaf(0.001), "0.0010");
        assert_eq!(format_vaf(1.0), "1.00");
        assert_eq!(format_vaf(0.05), "0.05");
    }

    #[test]
    fn test_format_metric() {
        assert_eq!(format_metric(1.0), "1.0000");
        assert_eq!(format_metric(0.5), "0.5000");
        assert_eq!(format_metric(f64::NAN), "N/A");
    }

    #[test]
    fn test_normalize_chrom() {
        assert_eq!(normalize_chrom("chr1"), "1");
        assert_eq!(normalize_chrom("1"), "1");
        assert_eq!(normalize_chrom("chrX"), "X");
        assert_eq!(normalize_chrom("X"), "X");
        assert_eq!(normalize_chrom("chrMT"), "MT");
    }
}
