//! Filter report generation for tumor-informed variant filtering.
//!
//! Produces a structured [`FilterReport`] summarizing detection outcomes,
//! per-criterion rejection counts, near-miss analysis, per-variant detail,
//! and threshold sensitivity analysis.

use std::io::Write;

use anyhow::Result;
use serde::Serialize;

use super::{FilterConfig, FilterResult};

// ---------------------------------------------------------------------------
// Report data structures
// ---------------------------------------------------------------------------

/// Top-level filter report containing all analysis sections.
#[derive(Debug, Clone, Serialize)]
pub struct FilterReport {
    /// High-level summary statistics.
    pub summary: FilterSummary,
    /// Breakdown of rejections per filter criterion.
    pub per_criterion: Vec<CriterionAnalysis>,
    /// Variants that narrowly missed one or more thresholds.
    pub near_misses: Vec<NearMiss>,
    /// Detailed per-variant match status and metrics.
    pub per_variant: Vec<VariantDetail>,
    /// Detection rate at various threshold multipliers.
    pub threshold_sensitivity: Vec<ThresholdPoint>,
}

/// Aggregate detection statistics.
#[derive(Debug, Clone, Serialize)]
pub struct FilterSummary {
    /// Total expected variants evaluated.
    pub total_expected: usize,
    /// Variants that passed all filter criteria ("Found").
    pub found: usize,
    /// Variants that did not pass ("Not Found"), irrespective of reason.
    pub not_found: usize,
    /// Subset of not-found: a variant was detected but failed one or more criteria.
    pub detected_but_filtered: usize,
    /// Subset of not-found: no matching variant was detected at all.
    pub not_detected: usize,
    /// Fraction of expected variants that were found (found / total_expected).
    pub detection_rate: f64,
    /// Human-readable description of the filter configuration used.
    pub config_description: String,
}

/// Rejection statistics for a single filter criterion.
#[derive(Debug, Clone, Serialize)]
pub struct CriterionAnalysis {
    /// Criterion name (e.g. "coverage", "VAF", "expression", "type").
    pub name: String,
    /// Threshold value as a human-readable string.
    pub threshold: String,
    /// Number of variants rejected *only* by this criterion (no other failures).
    pub rejected_alone: usize,
    /// Number of variants where this criterion was among the failures.
    pub rejected_any: usize,
}

/// A variant that narrowly missed a filter threshold.
#[derive(Debug, Clone, Serialize)]
pub struct NearMiss {
    /// Genomic coordinates and alleles.
    pub chrom: String,
    pub pos: u64,
    pub ref_allele: String,
    pub alt_allele: String,
    /// Which criterion was narrowly missed.
    pub criterion: String,
    /// Observed metric value.
    pub observed_value: f64,
    /// Threshold that was not met.
    pub threshold: f64,
    /// Ratio of observed / threshold (values in (0.5, 1.0) are near-misses).
    pub ratio: f64,
}

/// Detailed per-variant information for the report.
#[derive(Debug, Clone, Serialize)]
pub struct VariantDetail {
    pub chrom: String,
    pub pos: u64,
    pub ref_allele: String,
    pub alt_allele: String,
    pub variant_type: String,
    /// "found", "detected_but_filtered", or "not_detected".
    pub status: String,
    /// Observed metric values (None when variant was not detected).
    pub kmer_vaf: Option<f64>,
    pub kmer_min_coverage: Option<u64>,
    pub kmer_expression: Option<f64>,
    /// Which criteria failed, if any.
    pub failed_criteria: Vec<String>,
}

/// Detection rate at a simulated threshold multiplier.
#[derive(Debug, Clone, Serialize)]
pub struct ThresholdPoint {
    /// Multiplier applied to all numeric thresholds (e.g. 0.5, 1.0, 2.0).
    pub multiplier: f64,
    /// Label for the multiplier (e.g. "0.5x", "1x").
    pub label: String,
    /// Number of variants that would be found at this threshold.
    pub found: usize,
    /// Detection rate at this threshold (found / total with detection data).
    pub detection_rate: f64,
}

// ---------------------------------------------------------------------------
// Filter-notes parsing helpers
// ---------------------------------------------------------------------------

/// Parsed criterion failure extracted from a filter_notes string.
#[derive(Debug, Clone, PartialEq)]
struct ParsedFailure {
    criterion: String,
    observed: f64,
    threshold: f64,
}

/// Parse the filter_notes field into individual criterion failures.
///
/// Expected formats:
///   - "PASS"
///   - "No matching variant detected"
///   - "coverage 5 < 10; VAF 0.001200 < 0.01; expression 0.50 < 1"
///   - "type snv not in allowed types"
fn parse_filter_notes(notes: &str) -> Vec<ParsedFailure> {
    if notes == "PASS" || notes == "No matching variant detected" {
        return Vec::new();
    }

    let mut failures = Vec::new();
    for part in notes.split("; ") {
        let part = part.trim();
        if part.starts_with("coverage ") {
            if let Some(f) = parse_numeric_failure("coverage", part) {
                failures.push(f);
            }
        } else if part.starts_with("VAF ") {
            if let Some(f) = parse_numeric_failure("VAF", part) {
                failures.push(f);
            }
        } else if part.starts_with("expression ") {
            if let Some(f) = parse_numeric_failure("expression", part) {
                failures.push(f);
            }
        } else if part.starts_with("type ") {
            // Type failures don't have numeric observed/threshold values,
            // but we still track them for criterion counting.
            failures.push(ParsedFailure {
                criterion: "type".to_string(),
                observed: f64::NAN,
                threshold: f64::NAN,
            });
        }
    }
    failures
}

/// Parse a "criterion OBSERVED < THRESHOLD" fragment.
fn parse_numeric_failure(criterion: &str, text: &str) -> Option<ParsedFailure> {
    // Strip the criterion prefix and parse "OBSERVED < THRESHOLD"
    let rest = text.strip_prefix(criterion)?.trim();
    let parts: Vec<&str> = rest.split('<').collect();
    if parts.len() != 2 {
        return None;
    }
    let observed: f64 = parts[0].trim().parse().ok()?;
    let threshold: f64 = parts[1].trim().parse().ok()?;
    Some(ParsedFailure {
        criterion: criterion.to_string(),
        observed,
        threshold,
    })
}

/// Classify a FilterResult into one of three categories.
fn classify_result(result: &FilterResult) -> &'static str {
    if result.found == "Found" {
        "found"
    } else if result.filter_notes == "No matching variant detected" {
        "not_detected"
    } else {
        "detected_but_filtered"
    }
}

// ---------------------------------------------------------------------------
// Report generation
// ---------------------------------------------------------------------------

/// Build a comprehensive [`FilterReport`] from a set of filter results and the
/// configuration that produced them.
pub fn generate_report(results: &[FilterResult], config: &FilterConfig) -> FilterReport {
    let total = results.len();
    let found_count = results.iter().filter(|r| r.found == "Found").count();
    let not_found_count = total - found_count;
    let detected_but_filtered = results
        .iter()
        .filter(|r| classify_result(r) == "detected_but_filtered")
        .count();
    let not_detected = results
        .iter()
        .filter(|r| classify_result(r) == "not_detected")
        .count();

    let detection_rate = if total > 0 {
        found_count as f64 / total as f64
    } else {
        0.0
    };

    let config_description = format_config_description(config);

    let summary = FilterSummary {
        total_expected: total,
        found: found_count,
        not_found: not_found_count,
        detected_but_filtered,
        not_detected,
        detection_rate,
        config_description,
    };

    let per_criterion = build_per_criterion(results, config);
    let near_misses = build_near_misses(results, config);
    let per_variant = build_per_variant(results);
    let threshold_sensitivity = build_threshold_sensitivity(results, config);

    FilterReport {
        summary,
        per_criterion,
        near_misses,
        per_variant,
        threshold_sensitivity,
    }
}

fn format_config_description(config: &FilterConfig) -> String {
    let mut parts = vec![
        format!("min_coverage={}", config.min_coverage),
        format!("min_vaf={}", config.min_vaf),
        format!("min_expression={}", config.min_expression),
    ];
    if config.use_alt {
        parts.push("mode=alt".to_string());
    } else {
        parts.push("mode=reference".to_string());
    }
    if !config.types.is_empty() {
        parts.push(format!("types=[{}]", config.types.join(",")));
    }
    parts.join(", ")
}

// ---------------------------------------------------------------------------
// Per-criterion analysis
// ---------------------------------------------------------------------------

fn build_per_criterion(results: &[FilterResult], config: &FilterConfig) -> Vec<CriterionAnalysis> {
    // Only look at results that were detected but filtered.
    let filtered: Vec<_> = results
        .iter()
        .filter(|r| classify_result(r) == "detected_but_filtered")
        .collect();

    // For each filtered result, parse its failure set.
    let parsed: Vec<Vec<ParsedFailure>> = filtered
        .iter()
        .map(|r| parse_filter_notes(&r.filter_notes))
        .collect();

    let criteria = [
        ("coverage", format!("{}", config.min_coverage)),
        ("VAF", format!("{}", config.min_vaf)),
        ("expression", format!("{}", config.min_expression)),
        (
            "type",
            if config.types.is_empty() {
                "all allowed".to_string()
            } else {
                config.types.join(",")
            },
        ),
    ];

    criteria
        .iter()
        .map(|(name, threshold)| {
            let rejected_any = parsed
                .iter()
                .filter(|failures| failures.iter().any(|f| f.criterion == *name))
                .count();

            let rejected_alone = parsed
                .iter()
                .filter(|failures| {
                    failures.len() == 1
                        && failures.first().map(|f| f.criterion.as_str()) == Some(*name)
                })
                .count();

            CriterionAnalysis {
                name: name.to_string(),
                threshold: threshold.clone(),
                rejected_alone,
                rejected_any,
            }
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Near-miss detection
// ---------------------------------------------------------------------------

/// A variant is a near-miss for a criterion when:
/// - The variant was *not* found (detected but filtered or not detected with metrics),
/// - The observed metric value is within a factor of 2 of the threshold
///   (i.e. ratio = observed / threshold is in [0.5, 1.0) for criteria where
///   higher is better).
fn build_near_misses(results: &[FilterResult], config: &FilterConfig) -> Vec<NearMiss> {
    let mut misses = Vec::new();

    for r in results {
        // Only consider results that were not found.
        if r.found == "Found" {
            continue;
        }

        // Coverage near-miss.
        if let Some(cov) = r.kmer_min_coverage {
            let threshold = config.min_coverage as f64;
            if threshold > 0.0 {
                let observed = cov as f64;
                let ratio = observed / threshold;
                if ratio >= 0.5 && ratio < 1.0 {
                    misses.push(NearMiss {
                        chrom: r.chrom.clone(),
                        pos: r.pos,
                        ref_allele: r.ref_allele.clone(),
                        alt_allele: r.alt_allele.clone(),
                        criterion: "coverage".to_string(),
                        observed_value: observed,
                        threshold,
                        ratio,
                    });
                }
            }
        }

        // VAF near-miss.
        if let Some(vaf) = r.kmer_vaf {
            let threshold = config.min_vaf;
            if threshold > 0.0 {
                let ratio = vaf / threshold;
                if ratio >= 0.5 && ratio < 1.0 {
                    misses.push(NearMiss {
                        chrom: r.chrom.clone(),
                        pos: r.pos,
                        ref_allele: r.ref_allele.clone(),
                        alt_allele: r.alt_allele.clone(),
                        criterion: "VAF".to_string(),
                        observed_value: vaf,
                        threshold,
                        ratio,
                    });
                }
            }
        }

        // Expression near-miss.
        if let Some(expr) = r.kmer_expression {
            let threshold = config.min_expression;
            if threshold > 0.0 {
                let ratio = expr / threshold;
                if ratio >= 0.5 && ratio < 1.0 {
                    misses.push(NearMiss {
                        chrom: r.chrom.clone(),
                        pos: r.pos,
                        ref_allele: r.ref_allele.clone(),
                        alt_allele: r.alt_allele.clone(),
                        criterion: "expression".to_string(),
                        observed_value: expr,
                        threshold,
                        ratio,
                    });
                }
            }
        }
    }

    misses
}

// ---------------------------------------------------------------------------
// Per-variant detail
// ---------------------------------------------------------------------------

fn build_per_variant(results: &[FilterResult]) -> Vec<VariantDetail> {
    results
        .iter()
        .map(|r| {
            let status = classify_result(r).to_string();
            let failed_criteria = if status == "detected_but_filtered" {
                parse_filter_notes(&r.filter_notes)
                    .into_iter()
                    .map(|f| f.criterion)
                    .collect()
            } else {
                Vec::new()
            };

            VariantDetail {
                chrom: r.chrom.clone(),
                pos: r.pos,
                ref_allele: r.ref_allele.clone(),
                alt_allele: r.alt_allele.clone(),
                variant_type: r.variant_type.clone(),
                status,
                kmer_vaf: r.kmer_vaf,
                kmer_min_coverage: r.kmer_min_coverage,
                kmer_expression: r.kmer_expression,
                failed_criteria,
            }
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Threshold sensitivity
// ---------------------------------------------------------------------------

/// Simulate detection rates at various threshold multipliers.
///
/// For each multiplier, re-evaluate every result that has detection data
/// (i.e. not "not_detected") against the scaled thresholds.
fn build_threshold_sensitivity(
    results: &[FilterResult],
    config: &FilterConfig,
) -> Vec<ThresholdPoint> {
    let multipliers = [0.5, 1.0, 2.0, 5.0];

    // Only results that were actually detected (have metric values).
    let detectable: Vec<_> = results
        .iter()
        .filter(|r| classify_result(r) != "not_detected")
        .collect();

    let total = results.len();

    // Results with no detection never pass regardless of threshold — they are
    // excluded from `detectable` above and implicitly reduce the detection rate.

    multipliers
        .iter()
        .map(|&mult| {
            let scaled_coverage = (config.min_coverage as f64 * mult) as u64;
            let scaled_vaf = config.min_vaf * mult;
            let scaled_expression = config.min_expression * mult;

            let found_at = detectable
                .iter()
                .filter(|r| {
                    let cov_ok = r
                        .kmer_min_coverage
                        .map_or(false, |c| c >= scaled_coverage);
                    let vaf_ok = r.kmer_vaf.map_or(false, |v| v >= scaled_vaf);
                    let expr_ok = r
                        .kmer_expression
                        .map_or(false, |e| e >= scaled_expression);
                    // Type filter does not scale with multiplier.
                    cov_ok && vaf_ok && expr_ok
                })
                .count();

            let rate = if total > 0 {
                found_at as f64 / total as f64
            } else {
                0.0
            };

            let label = if mult == 1.0 {
                "1x (current)".to_string()
            } else {
                format!("{}x", mult)
            };

            ThresholdPoint {
                multiplier: mult,
                label,
                found: found_at,
                detection_rate: rate,
            }
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Report writers
// ---------------------------------------------------------------------------

/// Write the filter report as pretty-printed JSON to the given writer.
pub fn write_report_json<W: Write>(report: &FilterReport, writer: &mut W) -> Result<()> {
    serde_json::to_writer_pretty(&mut *writer, report)?;
    writeln!(writer)?;
    Ok(())
}

/// Write a human-readable text report to the given writer.
pub fn write_report_text<W: Write>(report: &FilterReport, writer: &mut W) -> Result<()> {
    writeln!(writer, "==========================================")?;
    writeln!(writer, "  Filter Report")?;
    writeln!(writer, "==========================================")?;
    writeln!(writer)?;

    // -- Summary --
    let s = &report.summary;
    writeln!(writer, "--- Summary ---")?;
    writeln!(writer, "Config: {}", s.config_description)?;
    writeln!(writer, "Total expected variants:  {}", s.total_expected)?;
    writeln!(writer, "Found (pass all filters): {}", s.found)?;
    writeln!(writer, "Not found:                {}", s.not_found)?;
    writeln!(writer, "  Detected but filtered:  {}", s.detected_but_filtered)?;
    writeln!(writer, "  Not detected at all:    {}", s.not_detected)?;
    writeln!(writer, "Detection rate:           {:.1}%", s.detection_rate * 100.0)?;
    writeln!(writer)?;

    // -- Per-criterion --
    writeln!(writer, "--- Per-Criterion Rejection Analysis ---")?;
    if report.per_criterion.is_empty() || report.per_criterion.iter().all(|c| c.rejected_any == 0) {
        writeln!(writer, "(no criterion rejections)")?;
    } else {
        writeln!(
            writer,
            "{:<12} {:<16} {:>14} {:>14}",
            "Criterion", "Threshold", "Rejected Alone", "Rejected Any"
        )?;
        for c in &report.per_criterion {
            writeln!(
                writer,
                "{:<12} {:<16} {:>14} {:>14}",
                c.name, c.threshold, c.rejected_alone, c.rejected_any
            )?;
        }
    }
    writeln!(writer)?;

    // -- Near misses --
    writeln!(writer, "--- Near Misses (within 2x of threshold) ---")?;
    if report.near_misses.is_empty() {
        writeln!(writer, "(none)")?;
    } else {
        for nm in &report.near_misses {
            writeln!(
                writer,
                "  {}:{} {}/{} — {} {:.4} vs threshold {:.4} (ratio {:.3})",
                nm.chrom,
                nm.pos,
                nm.ref_allele,
                nm.alt_allele,
                nm.criterion,
                nm.observed_value,
                nm.threshold,
                nm.ratio,
            )?;
        }
    }
    writeln!(writer)?;

    // -- Per-variant detail --
    writeln!(writer, "--- Per-Variant Detail ---")?;
    for v in &report.per_variant {
        let metrics = match (v.kmer_vaf, v.kmer_min_coverage, v.kmer_expression) {
            (Some(vaf), Some(cov), Some(expr)) => {
                format!("VAF={:.6} cov={} expr={:.2}", vaf, cov, expr)
            }
            _ => "no metrics".to_string(),
        };
        let failures = if v.failed_criteria.is_empty() {
            String::new()
        } else {
            format!(" [failed: {}]", v.failed_criteria.join(", "))
        };
        writeln!(
            writer,
            "  {}:{} {}/{} ({}) — {}: {}{}",
            v.chrom,
            v.pos,
            v.ref_allele,
            v.alt_allele,
            v.variant_type,
            v.status,
            metrics,
            failures,
        )?;
    }
    writeln!(writer)?;

    // -- Threshold sensitivity --
    writeln!(writer, "--- Threshold Sensitivity ---")?;
    writeln!(
        writer,
        "{:<16} {:>6} {:>14}",
        "Threshold", "Found", "Detection Rate"
    )?;
    for tp in &report.threshold_sensitivity {
        writeln!(
            writer,
            "{:<16} {:>6} {:>13.1}%",
            tp.label,
            tp.found,
            tp.detection_rate * 100.0,
        )?;
    }
    writeln!(writer)?;

    Ok(())
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    /// Helper to build a default FilterConfig for testing.
    fn test_config() -> FilterConfig {
        FilterConfig {
            min_coverage: 10,
            min_vaf: 0.01,
            min_expression: 1.0,
            use_alt: false,
            types: vec![],
        }
    }

    /// Helper to build a FilterResult that passes all filters.
    fn make_found(chrom: &str, pos: u64, ref_a: &str, alt_a: &str) -> FilterResult {
        FilterResult {
            sample: "sample1".into(),
            chrom: chrom.into(),
            pos,
            ref_allele: ref_a.into(),
            alt_allele: alt_a.into(),
            variant_type: "SNV".into(),
            found: "Found".into(),
            filter_notes: "PASS".into(),
            kmer_vaf: Some(0.05),
            kmer_min_coverage: Some(100),
            kmer_expression: Some(5.0),
            ref_sequence: Some("ATCG".into()),
            variant_sequence: Some("ATGG".into()),
        }
    }

    /// Helper: detected but filtered on coverage.
    fn make_filtered_coverage(chrom: &str, pos: u64, cov: u64) -> FilterResult {
        FilterResult {
            sample: "sample1".into(),
            chrom: chrom.into(),
            pos,
            ref_allele: "A".into(),
            alt_allele: "T".into(),
            variant_type: "SNV".into(),
            found: "Not Found".into(),
            filter_notes: format!("coverage {} < 10", cov),
            kmer_vaf: Some(0.05),
            kmer_min_coverage: Some(cov),
            kmer_expression: Some(5.0),
            ref_sequence: None,
            variant_sequence: None,
        }
    }

    /// Helper: detected but filtered on multiple criteria.
    fn make_filtered_multi(
        chrom: &str,
        pos: u64,
        cov: u64,
        vaf: f64,
    ) -> FilterResult {
        FilterResult {
            sample: "sample1".into(),
            chrom: chrom.into(),
            pos,
            ref_allele: "C".into(),
            alt_allele: "G".into(),
            variant_type: "SNV".into(),
            found: "Not Found".into(),
            filter_notes: format!("coverage {} < 10; VAF {:.6} < 0.01", cov, vaf),
            kmer_vaf: Some(vaf),
            kmer_min_coverage: Some(cov),
            kmer_expression: Some(5.0),
            ref_sequence: None,
            variant_sequence: None,
        }
    }

    /// Helper: not detected at all.
    fn make_not_detected(chrom: &str, pos: u64) -> FilterResult {
        FilterResult {
            sample: String::new(),
            chrom: chrom.into(),
            pos,
            ref_allele: "G".into(),
            alt_allele: "A".into(),
            variant_type: "SNV".into(),
            found: "Not Found".into(),
            filter_notes: "No matching variant detected".into(),
            kmer_vaf: None,
            kmer_min_coverage: None,
            kmer_expression: None,
            ref_sequence: None,
            variant_sequence: None,
        }
    }

    // -- generate_report with mixed results --

    #[test]
    fn test_generate_report_summary_counts() {
        let results = vec![
            make_found("chr1", 100, "A", "T"),
            make_found("chr1", 200, "C", "G"),
            make_filtered_coverage("chr2", 300, 5),
            make_not_detected("chr3", 400),
        ];
        let config = test_config();
        let report = generate_report(&results, &config);

        assert_eq!(report.summary.total_expected, 4);
        assert_eq!(report.summary.found, 2);
        assert_eq!(report.summary.not_found, 2);
        assert_eq!(report.summary.detected_but_filtered, 1);
        assert_eq!(report.summary.not_detected, 1);
        assert!((report.summary.detection_rate - 0.5).abs() < 1e-9);
    }

    #[test]
    fn test_generate_report_empty_input() {
        let config = test_config();
        let report = generate_report(&[], &config);

        assert_eq!(report.summary.total_expected, 0);
        assert_eq!(report.summary.found, 0);
        assert_eq!(report.summary.not_found, 0);
        assert_eq!(report.summary.detection_rate, 0.0);
        assert!(report.per_variant.is_empty());
        assert!(report.near_misses.is_empty());
    }

    #[test]
    fn test_generate_report_all_found() {
        let results = vec![
            make_found("chr1", 100, "A", "T"),
            make_found("chr1", 200, "C", "G"),
        ];
        let config = test_config();
        let report = generate_report(&results, &config);

        assert_eq!(report.summary.found, 2);
        assert_eq!(report.summary.not_found, 0);
        assert!((report.summary.detection_rate - 1.0).abs() < 1e-9);
    }

    // -- Per-criterion counting --

    #[test]
    fn test_per_criterion_single_failure() {
        let results = vec![
            make_found("chr1", 100, "A", "T"),
            make_filtered_coverage("chr2", 200, 5),
        ];
        let config = test_config();
        let report = generate_report(&results, &config);

        let cov_criterion = report
            .per_criterion
            .iter()
            .find(|c| c.name == "coverage")
            .unwrap();
        assert_eq!(cov_criterion.rejected_alone, 1);
        assert_eq!(cov_criterion.rejected_any, 1);

        // VAF should have 0 rejections.
        let vaf_criterion = report
            .per_criterion
            .iter()
            .find(|c| c.name == "VAF")
            .unwrap();
        assert_eq!(vaf_criterion.rejected_alone, 0);
        assert_eq!(vaf_criterion.rejected_any, 0);
    }

    #[test]
    fn test_per_criterion_multi_failure() {
        // Variant fails both coverage and VAF.
        let results = vec![make_filtered_multi("chr1", 100, 3, 0.005)];
        let config = test_config();
        let report = generate_report(&results, &config);

        let cov = report
            .per_criterion
            .iter()
            .find(|c| c.name == "coverage")
            .unwrap();
        assert_eq!(cov.rejected_any, 1);
        assert_eq!(cov.rejected_alone, 0); // not alone, VAF also failed

        let vaf = report
            .per_criterion
            .iter()
            .find(|c| c.name == "VAF")
            .unwrap();
        assert_eq!(vaf.rejected_any, 1);
        assert_eq!(vaf.rejected_alone, 0);
    }

    #[test]
    fn test_per_criterion_mixed_alone_and_multi() {
        let results = vec![
            make_filtered_coverage("chr1", 100, 7),       // coverage alone
            make_filtered_multi("chr2", 200, 4, 0.008),   // coverage + VAF
        ];
        let config = test_config();
        let report = generate_report(&results, &config);

        let cov = report
            .per_criterion
            .iter()
            .find(|c| c.name == "coverage")
            .unwrap();
        assert_eq!(cov.rejected_any, 2);
        assert_eq!(cov.rejected_alone, 1);

        let vaf = report
            .per_criterion
            .iter()
            .find(|c| c.name == "VAF")
            .unwrap();
        assert_eq!(vaf.rejected_any, 1);
        assert_eq!(vaf.rejected_alone, 0);
    }

    // -- Near-miss detection --

    #[test]
    fn test_near_miss_coverage() {
        // Coverage 7 vs threshold 10 → ratio 0.7, within [0.5, 1.0).
        let results = vec![make_filtered_coverage("chr1", 100, 7)];
        let config = test_config();
        let report = generate_report(&results, &config);

        assert_eq!(report.near_misses.len(), 1);
        let nm = &report.near_misses[0];
        assert_eq!(nm.criterion, "coverage");
        assert!((nm.observed_value - 7.0).abs() < 1e-9);
        assert!((nm.threshold - 10.0).abs() < 1e-9);
        assert!((nm.ratio - 0.7).abs() < 1e-9);
    }

    #[test]
    fn test_near_miss_vaf() {
        // VAF 0.007 vs threshold 0.01 → ratio 0.7.
        let mut r = make_found("chr1", 100, "A", "T");
        r.found = "Not Found".into();
        r.filter_notes = "VAF 0.007000 < 0.01".into();
        r.kmer_vaf = Some(0.007);

        let config = test_config();
        let report = generate_report(&[r], &config);

        let vaf_miss: Vec<_> = report
            .near_misses
            .iter()
            .filter(|nm| nm.criterion == "VAF")
            .collect();
        assert_eq!(vaf_miss.len(), 1);
        assert!((vaf_miss[0].ratio - 0.7).abs() < 1e-9);
    }

    #[test]
    fn test_near_miss_not_triggered_below_half() {
        // Coverage 3 vs threshold 10 → ratio 0.3, below 0.5, not a near-miss.
        let results = vec![make_filtered_coverage("chr1", 100, 3)];
        let config = test_config();
        let report = generate_report(&results, &config);

        let cov_misses: Vec<_> = report
            .near_misses
            .iter()
            .filter(|nm| nm.criterion == "coverage")
            .collect();
        assert!(cov_misses.is_empty());
    }

    #[test]
    fn test_near_miss_not_triggered_when_found() {
        // Found variants should not appear as near-misses, even if they'd
        // otherwise qualify (they passed the threshold).
        let results = vec![make_found("chr1", 100, "A", "T")];
        let config = test_config();
        let report = generate_report(&results, &config);

        assert!(report.near_misses.is_empty());
    }

    #[test]
    fn test_near_miss_not_detected_variant() {
        // Not-detected variants have no metrics → no near-miss possible.
        let results = vec![make_not_detected("chr1", 100)];
        let config = test_config();
        let report = generate_report(&results, &config);

        assert!(report.near_misses.is_empty());
    }

    // -- Threshold sensitivity --

    #[test]
    fn test_threshold_sensitivity_multipliers() {
        let results = vec![
            make_found("chr1", 100, "A", "T"),        // cov=100, vaf=0.05, expr=5.0
            make_filtered_coverage("chr2", 200, 5),    // cov=5, vaf=0.05, expr=5.0
        ];
        let config = test_config();
        let report = generate_report(&results, &config);

        assert_eq!(report.threshold_sensitivity.len(), 4);

        // At 0.5x: thresholds are cov>=5, vaf>=0.005, expr>=0.5
        // Result 1: cov 100>=5 ✓, vaf 0.05>=0.005 ✓, expr 5.0>=0.5 ✓ → found
        // Result 2: cov 5>=5 ✓, vaf 0.05>=0.005 ✓, expr 5.0>=0.5 ✓ → found
        let half = &report.threshold_sensitivity[0];
        assert!((half.multiplier - 0.5).abs() < 1e-9);
        assert_eq!(half.found, 2);

        // At 1x: thresholds are cov>=10, vaf>=0.01, expr>=1.0
        // Result 1: found; Result 2: cov 5<10 → not found
        let one = &report.threshold_sensitivity[1];
        assert!((one.multiplier - 1.0).abs() < 1e-9);
        assert_eq!(one.found, 1);

        // At 5x: thresholds are cov>=50, vaf>=0.05, expr>=5.0
        // Result 1: cov 100>=50 ✓, vaf 0.05>=0.05 ✓, expr 5.0>=5.0 ✓ → found
        // Result 2: cov 5<50 → not found
        let five = &report.threshold_sensitivity[3];
        assert!((five.multiplier - 5.0).abs() < 1e-9);
        assert_eq!(five.found, 1);
    }

    #[test]
    fn test_threshold_sensitivity_not_detected_never_pass() {
        // A not-detected variant should never count as found at any multiplier.
        let results = vec![make_not_detected("chr1", 100)];
        let config = test_config();
        let report = generate_report(&results, &config);

        for tp in &report.threshold_sensitivity {
            assert_eq!(tp.found, 0);
        }
    }

    #[test]
    fn test_threshold_sensitivity_labels() {
        let config = test_config();
        let report = generate_report(&[], &config);

        let labels: Vec<&str> = report
            .threshold_sensitivity
            .iter()
            .map(|tp| tp.label.as_str())
            .collect();
        assert_eq!(labels, vec!["0.5x", "1x (current)", "2x", "5x"]);
    }

    // -- Per-variant detail --

    #[test]
    fn test_per_variant_detail_statuses() {
        let results = vec![
            make_found("chr1", 100, "A", "T"),
            make_filtered_coverage("chr2", 200, 5),
            make_not_detected("chr3", 300),
        ];
        let config = test_config();
        let report = generate_report(&results, &config);

        assert_eq!(report.per_variant.len(), 3);
        assert_eq!(report.per_variant[0].status, "found");
        assert_eq!(report.per_variant[1].status, "detected_but_filtered");
        assert_eq!(report.per_variant[2].status, "not_detected");
    }

    #[test]
    fn test_per_variant_failed_criteria() {
        let results = vec![make_filtered_multi("chr1", 100, 3, 0.005)];
        let config = test_config();
        let report = generate_report(&results, &config);

        let detail = &report.per_variant[0];
        assert_eq!(detail.failed_criteria.len(), 2);
        assert!(detail.failed_criteria.contains(&"coverage".to_string()));
        assert!(detail.failed_criteria.contains(&"VAF".to_string()));
    }

    #[test]
    fn test_per_variant_found_has_no_failures() {
        let results = vec![make_found("chr1", 100, "A", "T")];
        let config = test_config();
        let report = generate_report(&results, &config);

        assert!(report.per_variant[0].failed_criteria.is_empty());
    }

    // -- Report serialization --

    #[test]
    fn test_write_report_json_roundtrip() {
        let results = vec![
            make_found("chr1", 100, "A", "T"),
            make_filtered_coverage("chr2", 200, 7),
            make_not_detected("chr3", 300),
        ];
        let config = test_config();
        let report = generate_report(&results, &config);

        let mut buf = Vec::new();
        write_report_json(&report, &mut buf).unwrap();

        let json_str = String::from_utf8(buf).unwrap();
        let parsed: serde_json::Value = serde_json::from_str(&json_str).unwrap();

        // Verify top-level structure.
        assert!(parsed.get("summary").is_some());
        assert!(parsed.get("per_criterion").is_some());
        assert!(parsed.get("near_misses").is_some());
        assert!(parsed.get("per_variant").is_some());
        assert!(parsed.get("threshold_sensitivity").is_some());

        // Verify summary values.
        let summary = &parsed["summary"];
        assert_eq!(summary["total_expected"], 3);
        assert_eq!(summary["found"], 1);
        assert_eq!(summary["not_found"], 2);
    }

    #[test]
    fn test_write_report_text_contains_sections() {
        let results = vec![
            make_found("chr1", 100, "A", "T"),
            make_filtered_coverage("chr2", 200, 7),
            make_not_detected("chr3", 300),
        ];
        let config = test_config();
        let report = generate_report(&results, &config);

        let mut buf = Vec::new();
        write_report_text(&report, &mut buf).unwrap();
        let text = String::from_utf8(buf).unwrap();

        assert!(text.contains("Filter Report"));
        assert!(text.contains("--- Summary ---"));
        assert!(text.contains("--- Per-Criterion Rejection Analysis ---"));
        assert!(text.contains("--- Near Misses"));
        assert!(text.contains("--- Per-Variant Detail ---"));
        assert!(text.contains("--- Threshold Sensitivity ---"));

        // Verify specific content.
        assert!(text.contains("Total expected variants:  3"));
        assert!(text.contains("Found (pass all filters): 1"));
        assert!(text.contains("Detection rate:           33.3%"));
    }

    #[test]
    fn test_write_report_text_near_miss_shown() {
        // Coverage 7 vs threshold 10 → near-miss.
        let results = vec![make_filtered_coverage("chr2", 200, 7)];
        let config = test_config();
        let report = generate_report(&results, &config);

        let mut buf = Vec::new();
        write_report_text(&report, &mut buf).unwrap();
        let text = String::from_utf8(buf).unwrap();

        assert!(text.contains("chr2:200"));
        assert!(text.contains("coverage"));
        assert!(text.contains("0.700"));
    }

    // -- parse_filter_notes --

    #[test]
    fn test_parse_filter_notes_pass() {
        assert!(parse_filter_notes("PASS").is_empty());
    }

    #[test]
    fn test_parse_filter_notes_not_detected() {
        assert!(parse_filter_notes("No matching variant detected").is_empty());
    }

    #[test]
    fn test_parse_filter_notes_single_criterion() {
        let failures = parse_filter_notes("coverage 5 < 10");
        assert_eq!(failures.len(), 1);
        assert_eq!(failures[0].criterion, "coverage");
        assert!((failures[0].observed - 5.0).abs() < 1e-9);
        assert!((failures[0].threshold - 10.0).abs() < 1e-9);
    }

    #[test]
    fn test_parse_filter_notes_multiple_criteria() {
        let failures = parse_filter_notes("coverage 3 < 10; VAF 0.005000 < 0.01");
        assert_eq!(failures.len(), 2);
        assert_eq!(failures[0].criterion, "coverage");
        assert_eq!(failures[1].criterion, "VAF");
    }

    #[test]
    fn test_parse_filter_notes_type_failure() {
        let failures = parse_filter_notes("type snv not in allowed types");
        assert_eq!(failures.len(), 1);
        assert_eq!(failures[0].criterion, "type");
    }

    #[test]
    fn test_parse_filter_notes_expression() {
        let failures = parse_filter_notes("expression 0.50 < 1");
        assert_eq!(failures.len(), 1);
        assert_eq!(failures[0].criterion, "expression");
        assert!((failures[0].observed - 0.50).abs() < 1e-9);
        assert!((failures[0].threshold - 1.0).abs() < 1e-9);
    }

    // -- Config description --

    #[test]
    fn test_config_description_default() {
        let config = test_config();
        let report = generate_report(&[], &config);

        assert!(report.summary.config_description.contains("min_coverage=10"));
        assert!(report.summary.config_description.contains("min_vaf=0.01"));
        assert!(report.summary.config_description.contains("min_expression=1"));
        assert!(report.summary.config_description.contains("mode=reference"));
    }

    #[test]
    fn test_config_description_with_types_and_alt() {
        let config = FilterConfig {
            min_coverage: 5,
            min_vaf: 0.005,
            min_expression: 0.5,
            use_alt: true,
            types: vec!["SNV".into(), "DEL".into()],
        };
        let report = generate_report(&[], &config);

        assert!(report.summary.config_description.contains("mode=alt"));
        assert!(report.summary.config_description.contains("types=[SNV,DEL]"));
    }
}
