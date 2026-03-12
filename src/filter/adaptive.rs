//! Multi-stage adaptive filtering pipeline for variant calls.
//!
//! This module provides a composite filter that applies multiple stages of
//! quality-based and context-based filtering to reduce false positives.
//! It complements the existing coordinate/sequence matching filter by adding:
//!
//! - **Stage 1**: Hard thresholds on coverage, VAF, and expression
//! - **Stage 2**: Variant type filtering (optional allow-list)
//! - **Stage 3**: Quality filtering on QUAL score and p-value
//! - **Stage 4**: Sequence context annotation (homopolymers, GC content, motifs)
//!
//! Stage 4 annotates but does not reject -- it produces flags that downstream
//! consumers can use for soft filtering or manual review.

use crate::variant::VariantCall;

/// Configuration for adaptive filtering stages.
#[derive(Debug, Clone)]
pub struct AdaptiveFilterConfig {
    /// Stage 1: Min coverage (absolute)
    pub min_coverage: u64,
    /// Stage 1: Min VAF
    pub min_vaf: f64,
    /// Stage 1: Min expression
    pub min_expression: f64,
    /// Stage 2: Variant types to keep (empty = all)
    pub allowed_types: Vec<String>,
    /// Stage 3: Min QUAL score (Phred-scaled)
    pub min_qual: Option<f64>,
    /// Stage 3: Max p-value
    pub max_pvalue: Option<f64>,
    /// Stage 4: Enable context filtering
    pub context_filter: bool,
    /// Stage 4: Homopolymer length to flag
    pub homopolymer_threshold: usize,
    /// Stage 4: GC content range to flag
    pub gc_low: f64,
    pub gc_high: f64,
}

impl Default for AdaptiveFilterConfig {
    fn default() -> Self {
        Self {
            min_coverage: 3,
            min_vaf: 0.0,
            min_expression: 0.0,
            allowed_types: vec![],
            min_qual: None,
            max_pvalue: None,
            context_filter: false,
            homopolymer_threshold: 6,
            gc_low: 0.2,
            gc_high: 0.8,
        }
    }
}

/// Result of multi-stage filtering for a single call.
#[derive(Debug, Clone)]
pub struct StageFilterResult {
    /// Whether the call passed all stages.
    pub passed: bool,
    /// Which stages were applied and their results.
    pub stage_results: Vec<StageResult>,
    /// Context flags (if context filtering enabled).
    pub context_flags: Vec<String>,
}

#[derive(Debug, Clone)]
pub struct StageResult {
    pub stage: String,
    pub passed: bool,
    pub reason: Option<String>,
}

/// Apply multi-stage filtering to a variant call.
///
/// Stages are evaluated in order:
/// 1. Hard thresholds (always applied)
/// 2. Type filter (only if `allowed_types` is non-empty)
/// 3. Quality filter (only if `min_qual` or `max_pvalue` is set)
/// 4. Context annotation (only if `context_filter` is enabled; flags only, no rejection)
pub fn apply_stages(call: &VariantCall, config: &AdaptiveFilterConfig) -> StageFilterResult {
    let mut results = Vec::new();
    let mut all_passed = true;

    // Stage 1: Hard thresholds (coverage, VAF, expression)
    let s1 = stage_hard_thresholds(call, config);
    if !s1.passed {
        all_passed = false;
    }
    results.push(s1);

    // Stage 2: Type filter
    if !config.allowed_types.is_empty() {
        let s2 = stage_type_filter(call, &config.allowed_types);
        if !s2.passed {
            all_passed = false;
        }
        results.push(s2);
    }

    // Stage 3: Quality filter (QUAL, p-value)
    if config.min_qual.is_some() || config.max_pvalue.is_some() {
        let s3 = stage_quality_filter(call, config);
        if !s3.passed {
            all_passed = false;
        }
        results.push(s3);
    }

    // Stage 4: Context annotation (flags, doesn't hard-filter)
    let context_flags = if config.context_filter {
        annotate_context(call, config)
    } else {
        vec![]
    };

    StageFilterResult {
        passed: all_passed,
        stage_results: results,
        context_flags,
    }
}

fn stage_hard_thresholds(call: &VariantCall, config: &AdaptiveFilterConfig) -> StageResult {
    if call.min_coverage < config.min_coverage {
        return StageResult {
            stage: "hard_threshold".to_string(),
            passed: false,
            reason: Some(format!(
                "min_coverage {} < {}",
                call.min_coverage, config.min_coverage
            )),
        };
    }
    if call.rvaf < config.min_vaf {
        return StageResult {
            stage: "hard_threshold".to_string(),
            passed: false,
            reason: Some(format!("vaf {:.4} < {:.4}", call.rvaf, config.min_vaf)),
        };
    }
    if call.expression < config.min_expression {
        return StageResult {
            stage: "hard_threshold".to_string(),
            passed: false,
            reason: Some(format!(
                "expression {:.2} < {:.2}",
                call.expression, config.min_expression
            )),
        };
    }
    StageResult {
        stage: "hard_threshold".to_string(),
        passed: true,
        reason: None,
    }
}

fn stage_type_filter(call: &VariantCall, allowed: &[String]) -> StageResult {
    let type_str = format!("{}", call.variant_type);
    let passed = allowed
        .iter()
        .any(|t| t.eq_ignore_ascii_case(&type_str));
    StageResult {
        stage: "type_filter".to_string(),
        passed,
        reason: if !passed {
            Some(format!("type {} not in allowed list", type_str))
        } else {
            None
        },
    }
}

fn stage_quality_filter(call: &VariantCall, config: &AdaptiveFilterConfig) -> StageResult {
    if let Some(min_qual) = config.min_qual {
        if let Some(qual) = call.qual {
            if qual < min_qual {
                return StageResult {
                    stage: "quality".to_string(),
                    passed: false,
                    reason: Some(format!("QUAL {:.1} < {:.1}", qual, min_qual)),
                };
            }
        }
    }
    if let Some(max_p) = config.max_pvalue {
        if let Some(pval) = call.pvalue {
            if pval > max_p {
                return StageResult {
                    stage: "quality".to_string(),
                    passed: false,
                    reason: Some(format!("pvalue {:.2e} > {:.2e}", pval, max_p)),
                };
            }
        }
    }
    StageResult {
        stage: "quality".to_string(),
        passed: true,
        reason: None,
    }
}

/// Annotate sequence context (does not filter, just flags).
fn annotate_context(call: &VariantCall, config: &AdaptiveFilterConfig) -> Vec<String> {
    let mut flags = Vec::new();
    let seq = &call.ref_sequence;

    // Check for homopolymer
    if let Some(max_run) = max_homopolymer_run(seq) {
        if max_run >= config.homopolymer_threshold {
            flags.push(format!("HOMOPOLYMER_{}", max_run));
        }
    }

    // Check GC content
    let gc = gc_content(seq);
    if gc < config.gc_low {
        flags.push(format!("LOW_GC_{:.2}", gc));
    } else if gc > config.gc_high {
        flags.push(format!("HIGH_GC_{:.2}", gc));
    }

    // Check for GGC motif (known error-prone in Illumina sequencing)
    if seq.contains("GGC") || seq.contains("GCC") {
        flags.push("GGC_MOTIF".to_string());
    }

    // Check low complexity
    if is_low_complexity(seq) {
        flags.push("LOW_COMPLEXITY".to_string());
    }

    flags
}

/// Find the longest homopolymer run in a sequence.
pub fn max_homopolymer_run(seq: &str) -> Option<usize> {
    if seq.is_empty() {
        return None;
    }
    let bytes = seq.as_bytes();
    let mut max_run = 1;
    let mut current_run = 1;
    for i in 1..bytes.len() {
        if bytes[i].to_ascii_uppercase() == bytes[i - 1].to_ascii_uppercase() {
            current_run += 1;
            max_run = max_run.max(current_run);
        } else {
            current_run = 1;
        }
    }
    Some(max_run)
}

/// Compute GC content of a sequence.
pub fn gc_content(seq: &str) -> f64 {
    if seq.is_empty() {
        return 0.0;
    }
    let gc = seq
        .bytes()
        .filter(|&b| b == b'G' || b == b'C' || b == b'g' || b == b'c')
        .count();
    gc as f64 / seq.len() as f64
}

/// Check if a sequence has low linguistic complexity (many repeated k-mers).
///
/// Uses a simple heuristic: count distinct 3-mers vs possible 3-mers.
/// Low complexity if fewer than 50% of possible 3-mers are distinct.
pub fn is_low_complexity(seq: &str) -> bool {
    if seq.len() < 10 {
        return false;
    }
    // Simple check: count distinct 3-mers vs possible 3-mers
    let mut trimers = std::collections::HashSet::new();
    let bytes = seq.as_bytes();
    for i in 0..bytes.len().saturating_sub(2) {
        trimers.insert(&bytes[i..i + 3]);
    }
    let possible = seq.len() - 2;
    let distinct = trimers.len();
    // Low complexity if < 50% of possible 3-mers are distinct
    (distinct as f64) < (possible as f64 * 0.5)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::variant::{VariantCall, VariantType};

    /// Build a test VariantCall with the given parameters.
    fn make_call(
        vtype: VariantType,
        vaf: f64,
        cov: u64,
        expr: f64,
        ref_seq: &str,
    ) -> VariantCall {
        VariantCall {
            sample: "test".into(),
            target: "target".into(),
            variant_type: vtype,
            variant_name: "test_variant".into(),
            rvaf: vaf,
            expression: expr,
            min_coverage: cov,
            start_kmer_count: 100,
            ref_sequence: ref_seq.into(),
            alt_sequence: "ACGTACGT".into(),
            info: "vs_ref".into(),
            chrom: Some("chr1".into()),
            pos: Some(100),
            ref_allele: Some("A".into()),
            alt_allele: Some("T".into()),
            pvalue: None,
            qual: None,
            ci_lower: None,
            ci_upper: None,
            path_score: 0,
        }
    }

    fn make_call_with_qual(
        vaf: f64,
        cov: u64,
        qual: Option<f64>,
        pvalue: Option<f64>,
    ) -> VariantCall {
        VariantCall {
            sample: "test".into(),
            target: "target".into(),
            variant_type: VariantType::Substitution,
            variant_name: "test_variant".into(),
            rvaf: vaf,
            expression: 10.0,
            min_coverage: cov,
            start_kmer_count: 100,
            ref_sequence: "ACGTACGT".into(),
            alt_sequence: "TCGTACGT".into(),
            info: "vs_ref".into(),
            chrom: Some("chr1".into()),
            pos: Some(100),
            ref_allele: Some("A".into()),
            alt_allele: Some("T".into()),
            pvalue,
            qual,
            ci_lower: None,
            ci_upper: None,
            path_score: 0,
        }
    }

    // ── Stage 1: Hard thresholds ─────────────────────────────────────────

    #[test]
    fn test_hard_thresholds_pass() {
        let call = make_call(VariantType::Substitution, 0.1, 50, 10.0, "ACGTACGT");
        let config = AdaptiveFilterConfig {
            min_coverage: 10,
            min_vaf: 0.05,
            min_expression: 5.0,
            ..Default::default()
        };
        let result = stage_hard_thresholds(&call, &config);
        assert!(result.passed);
        assert!(result.reason.is_none());
    }

    #[test]
    fn test_hard_thresholds_fail_coverage() {
        let call = make_call(VariantType::Substitution, 0.1, 5, 10.0, "ACGTACGT");
        let config = AdaptiveFilterConfig {
            min_coverage: 10,
            ..Default::default()
        };
        let result = stage_hard_thresholds(&call, &config);
        assert!(!result.passed);
        assert!(result.reason.as_ref().unwrap().contains("min_coverage"));
        assert!(result.reason.as_ref().unwrap().contains("5"));
    }

    #[test]
    fn test_hard_thresholds_fail_vaf() {
        let call = make_call(VariantType::Substitution, 0.001, 50, 10.0, "ACGTACGT");
        let config = AdaptiveFilterConfig {
            min_vaf: 0.01,
            ..Default::default()
        };
        let result = stage_hard_thresholds(&call, &config);
        assert!(!result.passed);
        assert!(result.reason.as_ref().unwrap().contains("vaf"));
    }

    #[test]
    fn test_hard_thresholds_fail_expression() {
        let call = make_call(VariantType::Substitution, 0.1, 50, 0.5, "ACGTACGT");
        let config = AdaptiveFilterConfig {
            min_expression: 5.0,
            ..Default::default()
        };
        let result = stage_hard_thresholds(&call, &config);
        assert!(!result.passed);
        assert!(result.reason.as_ref().unwrap().contains("expression"));
    }

    #[test]
    fn test_hard_thresholds_boundary_exact_values() {
        // Values exactly at thresholds should pass
        let call = make_call(VariantType::Substitution, 0.05, 10, 5.0, "ACGTACGT");
        let config = AdaptiveFilterConfig {
            min_coverage: 10,
            min_vaf: 0.05,
            min_expression: 5.0,
            ..Default::default()
        };
        let result = stage_hard_thresholds(&call, &config);
        assert!(result.passed);
    }

    // ── Stage 2: Type filter ─────────────────────────────────────────────

    #[test]
    fn test_type_filter_allowed() {
        let call = make_call(VariantType::Substitution, 0.1, 50, 10.0, "ACGTACGT");
        let allowed = vec!["Substitution".to_string(), "Insertion".to_string()];
        let result = stage_type_filter(&call, &allowed);
        assert!(result.passed);
    }

    #[test]
    fn test_type_filter_not_allowed() {
        let call = make_call(VariantType::Deletion, 0.1, 50, 10.0, "ACGTACGT");
        let allowed = vec!["Substitution".to_string(), "Insertion".to_string()];
        let result = stage_type_filter(&call, &allowed);
        assert!(!result.passed);
        assert!(result.reason.as_ref().unwrap().contains("Deletion"));
        assert!(result.reason.as_ref().unwrap().contains("not in allowed list"));
    }

    #[test]
    fn test_type_filter_case_insensitive() {
        let call = make_call(VariantType::Substitution, 0.1, 50, 10.0, "ACGTACGT");
        let allowed = vec!["substitution".to_string()];
        let result = stage_type_filter(&call, &allowed);
        assert!(result.passed);
    }

    #[test]
    fn test_type_filter_itd() {
        let call = make_call(VariantType::Itd, 0.1, 50, 10.0, "ACGTACGT");
        let allowed = vec!["ITD".to_string()];
        let result = stage_type_filter(&call, &allowed);
        assert!(result.passed);
    }

    // ── Stage 3: Quality filter ──────────────────────────────────────────

    #[test]
    fn test_quality_filter_pass_qual() {
        let call = make_call_with_qual(0.1, 50, Some(30.0), None);
        let config = AdaptiveFilterConfig {
            min_qual: Some(20.0),
            ..Default::default()
        };
        let result = stage_quality_filter(&call, &config);
        assert!(result.passed);
    }

    #[test]
    fn test_quality_filter_fail_qual() {
        let call = make_call_with_qual(0.1, 50, Some(10.0), None);
        let config = AdaptiveFilterConfig {
            min_qual: Some(20.0),
            ..Default::default()
        };
        let result = stage_quality_filter(&call, &config);
        assert!(!result.passed);
        assert!(result.reason.as_ref().unwrap().contains("QUAL"));
        assert!(result.reason.as_ref().unwrap().contains("10.0"));
    }

    #[test]
    fn test_quality_filter_pass_pvalue() {
        let call = make_call_with_qual(0.1, 50, None, Some(0.001));
        let config = AdaptiveFilterConfig {
            max_pvalue: Some(0.01),
            ..Default::default()
        };
        let result = stage_quality_filter(&call, &config);
        assert!(result.passed);
    }

    #[test]
    fn test_quality_filter_fail_pvalue() {
        let call = make_call_with_qual(0.1, 50, None, Some(0.1));
        let config = AdaptiveFilterConfig {
            max_pvalue: Some(0.01),
            ..Default::default()
        };
        let result = stage_quality_filter(&call, &config);
        assert!(!result.passed);
        assert!(result.reason.as_ref().unwrap().contains("pvalue"));
    }

    #[test]
    fn test_quality_filter_no_qual_data_passes() {
        // If the call has no QUAL or p-value, the quality stage still passes
        // (there's nothing to fail on).
        let call = make_call_with_qual(0.1, 50, None, None);
        let config = AdaptiveFilterConfig {
            min_qual: Some(20.0),
            max_pvalue: Some(0.01),
            ..Default::default()
        };
        let result = stage_quality_filter(&call, &config);
        assert!(result.passed);
    }

    #[test]
    fn test_quality_filter_qual_takes_priority() {
        // QUAL fails, p-value would also fail, but QUAL check returns first.
        let call = make_call_with_qual(0.1, 50, Some(5.0), Some(0.5));
        let config = AdaptiveFilterConfig {
            min_qual: Some(20.0),
            max_pvalue: Some(0.01),
            ..Default::default()
        };
        let result = stage_quality_filter(&call, &config);
        assert!(!result.passed);
        assert!(result.reason.as_ref().unwrap().contains("QUAL"));
    }

    // ── Context annotation ───────────────────────────────────────────────

    #[test]
    fn test_annotate_context_homopolymer() {
        let call = make_call(
            VariantType::Substitution,
            0.1,
            50,
            10.0,
            "ACGTAAAAAAAACGT",
        );
        let config = AdaptiveFilterConfig {
            context_filter: true,
            homopolymer_threshold: 6,
            ..Default::default()
        };
        let flags = annotate_context(&call, &config);
        assert!(flags.iter().any(|f| f.starts_with("HOMOPOLYMER_")));
        assert!(flags.iter().any(|f| f.contains("8"))); // 8 A's
    }

    #[test]
    fn test_annotate_context_no_homopolymer_below_threshold() {
        let call = make_call(VariantType::Substitution, 0.1, 50, 10.0, "ACGTAAACGT");
        let config = AdaptiveFilterConfig {
            context_filter: true,
            homopolymer_threshold: 6,
            ..Default::default()
        };
        let flags = annotate_context(&call, &config);
        assert!(!flags.iter().any(|f| f.starts_with("HOMOPOLYMER_")));
    }

    #[test]
    fn test_annotate_context_low_gc() {
        // All A/T => GC = 0.0
        let call = make_call(VariantType::Substitution, 0.1, 50, 10.0, "AAAAAATTTTTT");
        let config = AdaptiveFilterConfig {
            context_filter: true,
            gc_low: 0.2,
            gc_high: 0.8,
            ..Default::default()
        };
        let flags = annotate_context(&call, &config);
        assert!(flags.iter().any(|f| f.starts_with("LOW_GC_")));
    }

    #[test]
    fn test_annotate_context_high_gc() {
        // All G/C => GC = 1.0
        let call = make_call(VariantType::Substitution, 0.1, 50, 10.0, "GGGGGCCCCCCC");
        let config = AdaptiveFilterConfig {
            context_filter: true,
            gc_low: 0.2,
            gc_high: 0.8,
            ..Default::default()
        };
        let flags = annotate_context(&call, &config);
        assert!(flags.iter().any(|f| f.starts_with("HIGH_GC_")));
    }

    #[test]
    fn test_annotate_context_normal_gc() {
        // 50% GC
        let call = make_call(VariantType::Substitution, 0.1, 50, 10.0, "ACGTACGT");
        let config = AdaptiveFilterConfig {
            context_filter: true,
            gc_low: 0.2,
            gc_high: 0.8,
            ..Default::default()
        };
        let flags = annotate_context(&call, &config);
        assert!(!flags.iter().any(|f| f.contains("GC_")));
    }

    #[test]
    fn test_annotate_context_ggc_motif() {
        let call = make_call(VariantType::Substitution, 0.1, 50, 10.0, "ATGGCATG");
        let config = AdaptiveFilterConfig {
            context_filter: true,
            ..Default::default()
        };
        let flags = annotate_context(&call, &config);
        assert!(flags.iter().any(|f| f == "GGC_MOTIF"));
    }

    #[test]
    fn test_annotate_context_gcc_motif() {
        let call = make_call(VariantType::Substitution, 0.1, 50, 10.0, "ATGCCATG");
        let config = AdaptiveFilterConfig {
            context_filter: true,
            ..Default::default()
        };
        let flags = annotate_context(&call, &config);
        assert!(flags.iter().any(|f| f == "GGC_MOTIF"));
    }

    #[test]
    fn test_annotate_context_low_complexity() {
        // "ATATATATAT" has only 2 distinct 3-mers (ATA, TAT) out of 8 possible
        let call = make_call(VariantType::Substitution, 0.1, 50, 10.0, "ATATATATAT");
        let config = AdaptiveFilterConfig {
            context_filter: true,
            ..Default::default()
        };
        let flags = annotate_context(&call, &config);
        assert!(flags.iter().any(|f| f == "LOW_COMPLEXITY"));
    }

    // ── max_homopolymer_run ──────────────────────────────────────────────

    #[test]
    fn test_max_homopolymer_run_empty() {
        assert_eq!(max_homopolymer_run(""), None);
    }

    #[test]
    fn test_max_homopolymer_run_single_base() {
        assert_eq!(max_homopolymer_run("A"), Some(1));
    }

    #[test]
    fn test_max_homopolymer_run_no_repeat() {
        assert_eq!(max_homopolymer_run("ACGT"), Some(1));
    }

    #[test]
    fn test_max_homopolymer_run_all_same() {
        assert_eq!(max_homopolymer_run("AAAAAAA"), Some(7));
    }

    #[test]
    fn test_max_homopolymer_run_mixed() {
        assert_eq!(max_homopolymer_run("AACCCGTTT"), Some(3));
    }

    #[test]
    fn test_max_homopolymer_run_at_end() {
        assert_eq!(max_homopolymer_run("ACGTTTTT"), Some(5));
    }

    #[test]
    fn test_max_homopolymer_run_at_start() {
        assert_eq!(max_homopolymer_run("GGGGGACGT"), Some(5));
    }

    // ── gc_content ───────────────────────────────────────────────────────

    #[test]
    fn test_gc_content_empty() {
        assert_eq!(gc_content(""), 0.0);
    }

    #[test]
    fn test_gc_content_all_gc() {
        assert!((gc_content("GCGCGC") - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_gc_content_no_gc() {
        assert!((gc_content("ATATAT") - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_gc_content_half() {
        assert!((gc_content("ACGT") - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_gc_content_lowercase() {
        assert!((gc_content("acgt") - 0.5).abs() < 1e-10);
    }

    // ── is_low_complexity ────────────────────────────────────────────────

    #[test]
    fn test_low_complexity_short_seq() {
        // Sequences shorter than 10 bases are never flagged
        assert!(!is_low_complexity("ATAT"));
        assert!(!is_low_complexity("ATATATAAA"));
    }

    #[test]
    fn test_low_complexity_repetitive() {
        // "ATATATATAT" => 3-mers: ATA, TAT repeated => 2 distinct out of 8
        assert!(is_low_complexity("ATATATATAT"));
    }

    #[test]
    fn test_low_complexity_diverse() {
        // Diverse sequence should not be flagged
        assert!(!is_low_complexity("ACGTACGATCAGT"));
    }

    #[test]
    fn test_low_complexity_homopolymer() {
        // "AAAAAAAAAA" => 3-mers: all "AAA" => 1 distinct out of 8
        assert!(is_low_complexity("AAAAAAAAAA"));
    }

    // ── apply_stages (full pipeline) ─────────────────────────────────────

    #[test]
    fn test_apply_stages_all_pass_default() {
        let call = make_call(VariantType::Substitution, 0.1, 50, 10.0, "ACGTACGT");
        let config = AdaptiveFilterConfig::default();
        let result = apply_stages(&call, &config);
        assert!(result.passed);
        // Only hard_threshold stage with default config (no type filter, no quality filter)
        assert_eq!(result.stage_results.len(), 1);
        assert_eq!(result.stage_results[0].stage, "hard_threshold");
        assert!(result.context_flags.is_empty());
    }

    #[test]
    fn test_apply_stages_hard_threshold_fail() {
        let call = make_call(VariantType::Substitution, 0.1, 1, 10.0, "ACGTACGT");
        let config = AdaptiveFilterConfig {
            min_coverage: 10,
            ..Default::default()
        };
        let result = apply_stages(&call, &config);
        assert!(!result.passed);
        assert!(!result.stage_results[0].passed);
    }

    #[test]
    fn test_apply_stages_with_type_filter() {
        let call = make_call(VariantType::Substitution, 0.1, 50, 10.0, "ACGTACGT");
        let config = AdaptiveFilterConfig {
            allowed_types: vec!["Insertion".to_string()],
            ..Default::default()
        };
        let result = apply_stages(&call, &config);
        assert!(!result.passed);
        // Should have 2 stages: hard_threshold + type_filter
        assert_eq!(result.stage_results.len(), 2);
        assert!(result.stage_results[0].passed); // hard thresholds pass
        assert!(!result.stage_results[1].passed); // type filter fails
    }

    #[test]
    fn test_apply_stages_with_quality_filter() {
        let mut call = make_call(VariantType::Substitution, 0.1, 50, 10.0, "ACGTACGT");
        call.qual = Some(5.0);
        call.pvalue = Some(0.3);
        let config = AdaptiveFilterConfig {
            min_qual: Some(20.0),
            max_pvalue: Some(0.01),
            ..Default::default()
        };
        let result = apply_stages(&call, &config);
        assert!(!result.passed);
        // Should have 2 stages: hard_threshold + quality
        assert_eq!(result.stage_results.len(), 2);
        assert!(result.stage_results[0].passed);
        assert!(!result.stage_results[1].passed);
    }

    #[test]
    fn test_apply_stages_with_context_flags() {
        let call = make_call(
            VariantType::Substitution,
            0.1,
            50,
            10.0,
            "AAAAAAAACGT",
        );
        let config = AdaptiveFilterConfig {
            context_filter: true,
            homopolymer_threshold: 6,
            ..Default::default()
        };
        let result = apply_stages(&call, &config);
        assert!(result.passed); // Context flags don't cause failure
        assert!(!result.context_flags.is_empty());
        assert!(result.context_flags.iter().any(|f| f.starts_with("HOMOPOLYMER_")));
    }

    #[test]
    fn test_apply_stages_all_stages_enabled() {
        let mut call = make_call(VariantType::Substitution, 0.1, 50, 10.0, "ACGTACGT");
        call.qual = Some(30.0);
        call.pvalue = Some(0.001);
        let config = AdaptiveFilterConfig {
            min_coverage: 10,
            min_vaf: 0.05,
            min_expression: 5.0,
            allowed_types: vec!["Substitution".to_string()],
            min_qual: Some(20.0),
            max_pvalue: Some(0.01),
            context_filter: true,
            homopolymer_threshold: 6,
            gc_low: 0.2,
            gc_high: 0.8,
        };
        let result = apply_stages(&call, &config);
        assert!(result.passed);
        // Should have 3 stages: hard_threshold + type_filter + quality
        assert_eq!(result.stage_results.len(), 3);
        for sr in &result.stage_results {
            assert!(sr.passed);
        }
    }

    #[test]
    fn test_apply_stages_multiple_failures() {
        let mut call = make_call(VariantType::Deletion, 0.001, 1, 0.1, "ACGTACGT");
        call.qual = Some(5.0);
        let config = AdaptiveFilterConfig {
            min_coverage: 10,
            min_vaf: 0.05,
            min_expression: 5.0,
            allowed_types: vec!["Substitution".to_string()],
            min_qual: Some(20.0),
            ..Default::default()
        };
        let result = apply_stages(&call, &config);
        assert!(!result.passed);
        // All three stages should be present and all should fail
        assert_eq!(result.stage_results.len(), 3);
        assert!(!result.stage_results[0].passed); // hard threshold
        assert!(!result.stage_results[1].passed); // type filter
        assert!(!result.stage_results[2].passed); // quality
    }
}
