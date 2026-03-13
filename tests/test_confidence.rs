// Confidence scoring integration tests.
//
// Tests the statistical modules used for variant confidence assessment:
// - binomial_sf: returns correct p-values for known inputs
// - fisher_combined_pvalue: works with multiple p-values
// - phred_qual: conversion is correct (p-value <-> QUAL)
// - bootstrap CI: produces reasonable intervals with MockDb
// - compute_variant_pvalue: end-to-end variant p-value via MockDb
// - estimate_error_rate: returns sensible rates from MockDb reference data

mod common;

use common::MockDb;
use kmerdet::confidence::pvalue::{
    binomial_pvalue, binomial_sf, compute_variant_pvalue, estimate_error_rate,
    fisher_combined_pvalue,
};
use kmerdet::confidence::qual::{phred_qual, pvalue_from_qual};
use kmerdet::sequence::path::KmerPath;
use kmerdet::variant::bootstrap::{bootstrap_confidence_intervals, BootstrapConfig};

// ---------------------------------------------------------------------------
// binomial_sf: known-value tests
// ---------------------------------------------------------------------------

#[test]
fn test_binomial_sf_fair_coin_known() {
    // P(X >= 50 | n=100, p=0.5) ~= 0.5398 (cumulative from 50 upward)
    let p = binomial_sf(50, 100, 0.5);
    assert!(
        (p - 0.5398).abs() < 0.01,
        "P(X>=50|100,0.5) should be ~0.54, got {}",
        p
    );
}

#[test]
fn test_binomial_sf_extreme_tail() {
    // P(X >= 90 | n=100, p=0.5) should be astronomically small
    let p = binomial_sf(90, 100, 0.5);
    assert!(p < 1e-15, "P(X>=90|100,0.5) should be tiny, got {}", p);
}

#[test]
fn test_binomial_sf_boundary_k_equals_one() {
    // P(X >= 1 | n=10, p=0.5) = 1 - P(X=0) = 1 - 0.5^10
    let p = binomial_sf(1, 10, 0.5);
    let expected = 1.0 - 0.5f64.powi(10);
    assert!(
        (p - expected).abs() < 1e-10,
        "P(X>=1|10,0.5) should be {}, got {}",
        expected,
        p
    );
}

#[test]
fn test_binomial_sf_boundary_k_equals_n() {
    // P(X >= n | n, p) = p^n
    let p = binomial_sf(5, 5, 0.3);
    let expected = 0.3f64.powi(5);
    assert!(
        (p - expected).abs() < 1e-10,
        "P(X>=5|5,0.3) should be {}, got {}",
        expected,
        p
    );
}

#[test]
fn test_binomial_sf_edge_case_k_zero() {
    // P(X >= 0) = 1.0 always
    assert_eq!(binomial_sf(0, 1000, 0.99), 1.0);
    assert_eq!(binomial_sf(0, 0, 0.5), 1.0);
}

#[test]
fn test_binomial_sf_edge_case_k_greater_than_n() {
    // P(X >= k) = 0 when k > n
    assert_eq!(binomial_sf(200, 100, 0.5), 0.0);
}

#[test]
fn test_binomial_sf_edge_case_p_zero() {
    // No successes possible
    assert_eq!(binomial_sf(1, 100, 0.0), 0.0);
}

#[test]
fn test_binomial_sf_edge_case_p_one() {
    // All successes guaranteed
    assert_eq!(binomial_sf(50, 100, 1.0), 1.0);
}

#[test]
fn test_binomial_sf_large_n_triggers_incomplete_beta() {
    // n > 1000 triggers the regularized incomplete beta code path
    // n > 1000 triggers the regularized incomplete beta code path
    // E[X] = n*p = 5000, so 5100 is slightly above mean
    let p = binomial_sf(5100, 10000, 0.5);
    // 5100 is slightly above mean=5000. P(X >= 5100) should be moderate but < 0.5
    assert!(p < 0.5, "P(X>=5100|10000,0.5) should be < 0.5, got {}", p);
    assert!(p > 0.01, "P(X>=5100|10000,0.5) should be > 0.01, got {}", p);
}

#[test]
fn test_binomial_sf_variant_detection_scenario() {
    // Realistic scenario: 10 variant k-mers out of 1000 total, error rate 0.001
    // E[errors] = 1, observing 10 is very unlikely under null
    let p = binomial_sf(10, 1000, 0.001);
    assert!(p < 1e-5, "10 out of 1000 at 0.001 error rate should be very significant, got {}", p);
}

// ---------------------------------------------------------------------------
// binomial_pvalue: wrapper tests
// ---------------------------------------------------------------------------

#[test]
fn test_binomial_pvalue_zero_count_returns_one() {
    assert_eq!(binomial_pvalue(0, 1000, 0.001), 1.0);
}

#[test]
fn test_binomial_pvalue_zero_coverage_returns_one() {
    assert_eq!(binomial_pvalue(5, 0, 0.001), 1.0);
}

#[test]
fn test_binomial_pvalue_strong_variant() {
    // 50 counts at 1000x with 0.001 error rate -> very significant
    let p = binomial_pvalue(50, 1000, 0.001);
    assert!(p < 1e-10, "50/1000 at 0.001 should be extremely significant, got {}", p);
}

#[test]
fn test_binomial_pvalue_borderline_variant() {
    // 2 counts at 1000x with 0.001 error rate -> E[X]=1, observing 2 is borderline
    let p = binomial_pvalue(2, 1000, 0.001);
    // Should not be very significant
    assert!(p > 0.05, "2/1000 at 0.001 should not be very significant, got {}", p);
}

// ---------------------------------------------------------------------------
// fisher_combined_pvalue: combining multiple p-values
// ---------------------------------------------------------------------------

#[test]
fn test_fisher_combined_empty_returns_one() {
    let p = fisher_combined_pvalue(&[], 31, 150);
    assert_eq!(p, 1.0);
}

#[test]
fn test_fisher_combined_single_pvalue_preserved() {
    // Single p-value with effective k clamped to 1: df=2
    // chi2 = -2*ln(0.05) = 5.99, P(chi2(2) > 5.99) ~= 0.05
    let p = fisher_combined_pvalue(&[0.05], 31, 150);
    assert!(
        (p - 0.05).abs() < 0.02,
        "Single p-value should be roughly preserved, got {}",
        p
    );
}

#[test]
fn test_fisher_combined_multiple_significant() {
    // Multiple very small p-values should combine to even smaller
    let pvals = vec![0.001, 0.001, 0.001, 0.001, 0.001];
    let combined = fisher_combined_pvalue(&pvals, 31, 150);
    assert!(combined < 0.0001, "Combined small p-values should be very significant, got {}", combined);
    // Should be more significant than any individual
    assert!(combined < 0.001, "Combined should be more significant than individual");
}

#[test]
fn test_fisher_combined_all_nonsignificant() {
    // Multiple non-significant p-values should remain non-significant
    let pvals = vec![0.5, 0.6, 0.7, 0.8, 0.9];
    let combined = fisher_combined_pvalue(&pvals, 31, 150);
    assert!(combined > 0.1, "Non-significant p-values should stay non-significant, got {}", combined);
}

#[test]
fn test_fisher_combined_mixed_pvalues() {
    // Mix of significant and non-significant
    let pvals = vec![0.001, 0.5, 0.001, 0.5];
    let combined = fisher_combined_pvalue(&pvals, 31, 150);
    // Should be somewhat significant due to the two 0.001 values
    assert!(combined < 0.05, "Mixed p-values with two significant ones should be significant, got {}", combined);
}

#[test]
fn test_fisher_combined_autocorrelation_correction() {
    // Same p-values but different correction factors should yield different results.
    let pvals = vec![0.01, 0.01, 0.01, 0.01, 0.01];

    // No correction (read_length == 0 => correction = 1.0, effective_k = 5, df = 10)
    let p_no_corr = fisher_combined_pvalue(&pvals, 31, 0);
    // With correction (31/150 ~ 0.207, effective_k ~ 1.03, df ~ 2.07)
    let p_corr = fisher_combined_pvalue(&pvals, 31, 150);

    // The two should be different (autocorrelation correction has an effect).
    assert!(
        (p_corr - p_no_corr).abs() > 1e-15,
        "Autocorrelation correction should change the result: corrected={}, uncorrected={}",
        p_corr,
        p_no_corr
    );

    // With correction, effective df is lower. For the same chi2 statistic,
    // a chi-squared distribution with fewer df has a heavier tail at high values,
    // but the key point is that both should be significant.
    assert!(p_corr < 0.01, "Corrected should still be significant, got {}", p_corr);
    assert!(p_no_corr < 0.01, "Uncorrected should be significant, got {}", p_no_corr);
}

#[test]
fn test_fisher_combined_very_small_pvalues_no_overflow() {
    // Extremely small p-values should not cause inf/NaN
    let pvals = vec![1e-100, 1e-100, 1e-100];
    let combined = fisher_combined_pvalue(&pvals, 31, 150);
    assert!(combined.is_finite(), "Should handle very small p-values without overflow");
    assert!(combined >= 0.0, "Should be non-negative");
    assert!(combined <= 1.0, "Should be <= 1.0");
}

// ---------------------------------------------------------------------------
// phred_qual: p-value <-> QUAL conversion
// ---------------------------------------------------------------------------

#[test]
fn test_phred_qual_standard_values() {
    // p = 0.1 -> Q10
    assert!((phred_qual(0.1) - 10.0).abs() < 1e-10);
    // p = 0.01 -> Q20
    assert!((phred_qual(0.01) - 20.0).abs() < 1e-10);
    // p = 0.001 -> Q30
    assert!((phred_qual(0.001) - 30.0).abs() < 1e-10);
    // p = 0.0001 -> Q40
    assert!((phred_qual(0.0001) - 40.0).abs() < 1e-10);
}

#[test]
fn test_phred_qual_zero_pvalue_capped() {
    assert_eq!(phred_qual(0.0), 999.0);
}

#[test]
fn test_phred_qual_negative_pvalue_capped() {
    assert_eq!(phred_qual(-1.0), 999.0);
}

#[test]
fn test_phred_qual_one_returns_zero() {
    assert_eq!(phred_qual(1.0), 0.0);
}

#[test]
fn test_phred_qual_greater_than_one_returns_zero() {
    assert_eq!(phred_qual(1.5), 0.0);
}

#[test]
fn test_phred_qual_very_small_pvalue_capped_at_999() {
    let q = phred_qual(1e-200);
    assert_eq!(q, 999.0);
}

#[test]
fn test_phred_qual_roundtrip() {
    // Convert p-value to qual and back; should be identical
    let original = 0.00123;
    let qual = phred_qual(original);
    let recovered = pvalue_from_qual(qual);
    assert!(
        (original - recovered).abs() < 1e-12,
        "Roundtrip failed: {} -> {} -> {}",
        original,
        qual,
        recovered
    );
}

#[test]
fn test_phred_qual_roundtrip_multiple_values() {
    for &pval in &[0.5, 0.1, 0.01, 0.001, 1e-5, 1e-10, 1e-50] {
        let qual = phred_qual(pval);
        if qual >= 999.0 {
            continue; // Capped values can't roundtrip
        }
        let recovered = pvalue_from_qual(qual);
        assert!(
            (pval - recovered).abs() / pval < 1e-10,
            "Roundtrip failed for p={}: qual={}, recovered={}",
            pval,
            qual,
            recovered
        );
    }
}

#[test]
fn test_pvalue_from_qual_standard_values() {
    assert!((pvalue_from_qual(10.0) - 0.1).abs() < 1e-12);
    assert!((pvalue_from_qual(20.0) - 0.01).abs() < 1e-12);
    assert!((pvalue_from_qual(30.0) - 0.001).abs() < 1e-12);
}

#[test]
fn test_pvalue_from_qual_zero_returns_one() {
    assert_eq!(pvalue_from_qual(0.0), 1.0);
}

// ---------------------------------------------------------------------------
// estimate_error_rate: error estimation from reference k-mers via MockDb
// ---------------------------------------------------------------------------

#[test]
fn test_estimate_error_rate_default_when_insufficient_data() {
    // Fewer than 2 k-mers -> default 0.001
    let db = MockDb::new(5);
    let rate = estimate_error_rate(&db, &["ACGTA".to_string()]);
    assert!(
        (rate - 0.001).abs() < 1e-10,
        "Should return default 0.001 with insufficient data, got {}",
        rate
    );
}

#[test]
fn test_estimate_error_rate_default_when_no_coverage() {
    // k-mers with zero coverage -> default
    let db = MockDb::new(5);
    let ref_kmers: Vec<String> = (0..20).map(|i| format!("ACGT{}", "A".repeat(i % 4 + 1))).collect();
    // None of these k-mers have counts in the empty db
    let rate = estimate_error_rate(&db, &ref_kmers);
    assert!(
        (rate - 0.001).abs() < 1e-10,
        "Should return default with zero coverage, got {}",
        rate
    );
}

#[test]
fn test_estimate_error_rate_clean_reference() {
    // Build a clean reference: high counts for correct extensions, low for errors
    let mut db = MockDb::new(5);

    // Create an overlapping k-mer chain from a known sequence
    let sequence = "ACGTACGATCAGTACGATCA";
    let k = 5;
    let ref_kmers: Vec<String> = sequence
        .as_bytes()
        .windows(k)
        .map(|w| std::str::from_utf8(w).unwrap().to_string())
        .collect();

    // For each consecutive pair of k-mers, set extensions
    for i in 0..ref_kmers.len() - 1 {
        let kmer = &ref_kmers[i];
        let suffix = &kmer[1..];
        let next_base = ref_kmers[i + 1].as_bytes()[k - 1];

        for &base in &[b'A', b'C', b'G', b'T'] {
            let mut child = String::from(suffix);
            child.push(base as char);
            if base == next_base {
                db.set(&child, 1000); // Correct extension
            } else {
                db.set(&child, 1); // Error extension
            }
        }
    }

    let rate = estimate_error_rate(&db, &ref_kmers);
    // With 1000 correct and 1*3 errors per position: error ~ 3/1003 ~ 0.003
    assert!(rate < 0.01, "Clean reference should have low error rate, got {}", rate);
    assert!(rate > 0.0, "Error rate should be positive, got {}", rate);
}

#[test]
fn test_estimate_error_rate_noisy_reference() {
    // Build a noisy reference: moderate error counts
    let mut db = MockDb::new(5);

    let sequence = "ACGTACGATCAGTACGATCA";
    let k = 5;
    let ref_kmers: Vec<String> = sequence
        .as_bytes()
        .windows(k)
        .map(|w| std::str::from_utf8(w).unwrap().to_string())
        .collect();

    for i in 0..ref_kmers.len() - 1 {
        let kmer = &ref_kmers[i];
        let suffix = &kmer[1..];
        let next_base = ref_kmers[i + 1].as_bytes()[k - 1];

        for &base in &[b'A', b'C', b'G', b'T'] {
            let mut child = String::from(suffix);
            child.push(base as char);
            if base == next_base {
                db.set(&child, 100); // Correct extension
            } else {
                db.set(&child, 20); // Significant error
            }
        }
    }

    let rate = estimate_error_rate(&db, &ref_kmers);
    // 3*20 errors out of 100 + 3*20 = 160 total per k-mer: ~0.375
    // But clamped to max 0.1
    assert!(rate > 0.01, "Noisy reference should have higher error rate, got {}", rate);
    assert!(rate <= 0.1, "Error rate should be clamped to 0.1, got {}", rate);
}

// ---------------------------------------------------------------------------
// compute_variant_pvalue: end-to-end variant p-value with MockDb
// ---------------------------------------------------------------------------

#[test]
fn test_compute_variant_pvalue_reference_call() {
    // When alt path == ref path, there are no variant-specific k-mers -> p=1.0
    let mut db = MockDb::new(5);
    db.set("ACGTA", 1000);
    db.set("CGTAC", 1000);

    let ref_path = KmerPath {
        kmers: vec!["ACGTA".to_string(), "CGTAC".to_string()],
        is_reference: true,
    };
    let alt_path = KmerPath {
        kmers: vec!["ACGTA".to_string(), "CGTAC".to_string()],
        is_reference: false,
    };

    let p = compute_variant_pvalue(&db, &alt_path, &ref_path, 0.001, 150);
    assert_eq!(p, 1.0, "Reference call should have p-value of 1.0");
}

#[test]
fn test_compute_variant_pvalue_strong_variant() {
    let mut db = MockDb::new(5);

    // Reference k-mers at high coverage
    db.set("ACGTA", 1000);
    db.set("CGTAC", 1000);
    db.set("GTACG", 1000);

    // Variant-specific k-mers with moderate counts (real variant)
    db.set("CGTTT", 50);
    db.set("GTTTC", 50);

    let ref_path = KmerPath {
        kmers: vec![
            "ACGTA".to_string(),
            "CGTAC".to_string(),
            "GTACG".to_string(),
        ],
        is_reference: true,
    };
    let alt_path = KmerPath {
        kmers: vec![
            "ACGTA".to_string(),
            "CGTTT".to_string(),
            "GTTTC".to_string(),
        ],
        is_reference: false,
    };

    let p = compute_variant_pvalue(&db, &alt_path, &ref_path, 0.001, 150);
    assert!(
        p < 1e-10,
        "Strong variant should have very small p-value, got {}",
        p
    );
}

#[test]
fn test_compute_variant_pvalue_weak_variant() {
    let mut db = MockDb::new(5);

    // Reference k-mers at high coverage
    db.set("ACGTA", 1000);
    db.set("CGTAC", 1000);
    db.set("GTACG", 1000);

    // Variant-specific k-mers with very low counts (could be error)
    db.set("CGTTT", 1);
    db.set("GTTTC", 1);

    let ref_path = KmerPath {
        kmers: vec![
            "ACGTA".to_string(),
            "CGTAC".to_string(),
            "GTACG".to_string(),
        ],
        is_reference: true,
    };
    let alt_path = KmerPath {
        kmers: vec![
            "ACGTA".to_string(),
            "CGTTT".to_string(),
            "GTTTC".to_string(),
        ],
        is_reference: false,
    };

    let p = compute_variant_pvalue(&db, &alt_path, &ref_path, 0.001, 150);
    // 1 count at 1000x coverage with 0.001 error rate -> E[X]=1, so observing 1 is expected
    assert!(
        p > 0.1,
        "Weak variant (1 count) should not be significant, got {}",
        p
    );
}

#[test]
fn test_compute_variant_pvalue_empty_alt_path() {
    let db = MockDb::new(5);
    let ref_path = KmerPath {
        kmers: vec!["ACGTA".to_string()],
        is_reference: true,
    };
    let alt_path = KmerPath {
        kmers: vec![],
        is_reference: false,
    };

    let p = compute_variant_pvalue(&db, &alt_path, &ref_path, 0.001, 150);
    assert_eq!(p, 1.0, "Empty alt path should return 1.0");
}

#[test]
fn test_compute_variant_pvalue_zero_count_variant_kmers() {
    // Variant k-mers with zero counts get skipped
    let mut db = MockDb::new(5);
    db.set("ACGTA", 1000);
    db.set("CGTAC", 1000);
    // Variant k-mers not set -> 0 counts

    let ref_path = KmerPath {
        kmers: vec!["ACGTA".to_string(), "CGTAC".to_string()],
        is_reference: true,
    };
    let alt_path = KmerPath {
        kmers: vec!["ACGTA".to_string(), "TTTTT".to_string()],
        is_reference: false,
    };

    let p = compute_variant_pvalue(&db, &alt_path, &ref_path, 0.001, 150);
    assert_eq!(p, 1.0, "Zero-count variant k-mers should result in p=1.0");
}

// ---------------------------------------------------------------------------
// bootstrap_confidence_intervals: integration tests with MockDb
// ---------------------------------------------------------------------------

#[test]
fn test_bootstrap_single_path_tight_ci() {
    // A single path always has rVAF = 1.0, so CI should be [1.0, 1.0]
    let mut db = MockDb::new(4);
    db.set("ACGT", 1000);
    db.set("CGTA", 1000);
    db.set("GTAC", 1000);

    let path = KmerPath {
        kmers: vec![
            "ACGT".to_string(),
            "CGTA".to_string(),
            "GTAC".to_string(),
        ],
        is_reference: true,
    };

    let config = BootstrapConfig {
        n_replicates: 500,
        confidence_level: 0.95,
        seed: Some(42),
    };

    let cis = bootstrap_confidence_intervals(&[path], &db, &config);
    assert_eq!(cis.len(), 1);
    assert!(
        (cis[0].ci_lower - 1.0).abs() < 1e-9,
        "Single path CI lower should be 1.0, got {}",
        cis[0].ci_lower
    );
    assert!(
        (cis[0].ci_upper - 1.0).abs() < 1e-9,
        "Single path CI upper should be 1.0, got {}",
        cis[0].ci_upper
    );
    assert!(
        cis[0].std_error < 1e-9,
        "Single path std_error should be ~0, got {}",
        cis[0].std_error
    );
}

#[test]
fn test_bootstrap_two_paths_equal_coverage() {
    // Two paths with equal coverage should have rVAF ~0.5 each
    let mut db = MockDb::new(4);
    db.set("AAAA", 2000); // Shared
    db.set("AAAC", 1000); // Ref-only
    db.set("AACG", 1000); // Ref-only
    db.set("AAAT", 1000); // Alt-only
    db.set("AATG", 1000); // Alt-only

    let ref_path = KmerPath {
        kmers: vec![
            "AAAA".to_string(),
            "AAAC".to_string(),
            "AACG".to_string(),
        ],
        is_reference: true,
    };
    let alt_path = KmerPath {
        kmers: vec![
            "AAAA".to_string(),
            "AAAT".to_string(),
            "AATG".to_string(),
        ],
        is_reference: false,
    };

    let config = BootstrapConfig {
        n_replicates: 1000,
        confidence_level: 0.95,
        seed: Some(42),
    };

    let cis = bootstrap_confidence_intervals(&[ref_path, alt_path], &db, &config);
    assert_eq!(cis.len(), 2);

    // Both CIs should contain 0.5
    for (i, ci) in cis.iter().enumerate() {
        assert!(
            ci.ci_lower <= 0.55 && ci.ci_upper >= 0.45,
            "CI for path {} should contain ~0.5: [{}, {}]",
            i,
            ci.ci_lower,
            ci.ci_upper
        );
    }

    // CI width should be narrow at high coverage
    let width0 = cis[0].ci_upper - cis[0].ci_lower;
    assert!(
        width0 < 0.10,
        "High coverage CI width should be <0.10, got {}",
        width0
    );
}

#[test]
fn test_bootstrap_low_vs_high_coverage_ci_width() {
    // Low coverage should produce wider CI than high coverage
    let config = BootstrapConfig {
        n_replicates: 500,
        confidence_level: 0.95,
        seed: Some(42),
    };

    // High coverage
    let mut db_high = MockDb::new(4);
    db_high.set("AAAA", 2000);
    db_high.set("AAAC", 1000);
    db_high.set("AACG", 1000);
    db_high.set("AAAT", 1000);
    db_high.set("AATG", 1000);

    let ref_path = KmerPath {
        kmers: vec!["AAAA".to_string(), "AAAC".to_string(), "AACG".to_string()],
        is_reference: true,
    };
    let alt_path = KmerPath {
        kmers: vec!["AAAA".to_string(), "AAAT".to_string(), "AATG".to_string()],
        is_reference: false,
    };

    let cis_high = bootstrap_confidence_intervals(
        &[ref_path.clone(), alt_path.clone()],
        &db_high,
        &config,
    );

    // Low coverage
    let mut db_low = MockDb::new(4);
    db_low.set("AAAA", 20);
    db_low.set("AAAC", 10);
    db_low.set("AACG", 10);
    db_low.set("AAAT", 10);
    db_low.set("AATG", 10);

    let cis_low = bootstrap_confidence_intervals(&[ref_path, alt_path], &db_low, &config);

    let width_high = cis_high[0].ci_upper - cis_high[0].ci_lower;
    let width_low = cis_low[0].ci_upper - cis_low[0].ci_lower;

    assert!(
        width_low > width_high,
        "Low coverage CI ({}) should be wider than high coverage CI ({})",
        width_low,
        width_high
    );
}

#[test]
fn test_bootstrap_reproducibility_same_seed() {
    let mut db = MockDb::new(4);
    db.set("AAAA", 200);
    db.set("AAAC", 100);
    db.set("AACG", 100);
    db.set("AAAT", 100);
    db.set("AATG", 100);

    let ref_path = KmerPath {
        kmers: vec!["AAAA".to_string(), "AAAC".to_string(), "AACG".to_string()],
        is_reference: true,
    };
    let alt_path = KmerPath {
        kmers: vec!["AAAA".to_string(), "AAAT".to_string(), "AATG".to_string()],
        is_reference: false,
    };

    let config = BootstrapConfig {
        n_replicates: 200,
        confidence_level: 0.95,
        seed: Some(12345),
    };

    let cis1 = bootstrap_confidence_intervals(
        &[ref_path.clone(), alt_path.clone()],
        &db,
        &config,
    );
    let cis2 = bootstrap_confidence_intervals(&[ref_path, alt_path], &db, &config);

    for (a, b) in cis1.iter().zip(cis2.iter()) {
        assert!(
            (a.ci_lower - b.ci_lower).abs() < 1e-12,
            "Same seed should produce identical ci_lower"
        );
        assert!(
            (a.ci_upper - b.ci_upper).abs() < 1e-12,
            "Same seed should produce identical ci_upper"
        );
        assert!(
            (a.std_error - b.std_error).abs() < 1e-12,
            "Same seed should produce identical std_error"
        );
    }
}

#[test]
fn test_bootstrap_empty_paths() {
    let db = MockDb::new(4);
    let config = BootstrapConfig::default();
    let cis = bootstrap_confidence_intervals(&[], &db, &config);
    assert!(cis.is_empty(), "Empty paths should return empty CIs");
}

#[test]
fn test_bootstrap_zero_replicates() {
    let mut db = MockDb::new(4);
    db.set("ACGT", 100);

    let path = KmerPath {
        kmers: vec!["ACGT".to_string()],
        is_reference: true,
    };

    let config = BootstrapConfig {
        n_replicates: 0,
        confidence_level: 0.95,
        seed: Some(42),
    };

    let cis = bootstrap_confidence_intervals(&[path], &db, &config);
    assert!(cis.is_empty(), "Zero replicates should return empty CIs");
}

#[test]
fn test_bootstrap_asymmetric_paths() {
    // Ref path dominates: 90% ref, 10% alt -> rVAF_alt ~0.1
    let mut db = MockDb::new(4);
    db.set("AAAA", 1000); // Shared
    db.set("AAAC", 900);  // Ref-only
    db.set("AACG", 900);  // Ref-only
    db.set("AAAT", 100);  // Alt-only
    db.set("AATG", 100);  // Alt-only

    let ref_path = KmerPath {
        kmers: vec!["AAAA".to_string(), "AAAC".to_string(), "AACG".to_string()],
        is_reference: true,
    };
    let alt_path = KmerPath {
        kmers: vec!["AAAA".to_string(), "AAAT".to_string(), "AATG".to_string()],
        is_reference: false,
    };

    let config = BootstrapConfig {
        n_replicates: 1000,
        confidence_level: 0.95,
        seed: Some(42),
    };

    let cis = bootstrap_confidence_intervals(&[ref_path, alt_path], &db, &config);
    assert_eq!(cis.len(), 2);

    // Alt path CI should be centered around a low rVAF (roughly 0.1)
    // The exact value depends on NNLS behavior with the shared k-mer
    let alt_ci = &cis[1];
    assert!(
        alt_ci.ci_upper < 0.5,
        "Alt path (low frequency) CI upper should be < 0.5, got {}",
        alt_ci.ci_upper
    );
    assert!(
        alt_ci.ci_lower >= 0.0,
        "CI lower should be non-negative, got {}",
        alt_ci.ci_lower
    );
}

// ---------------------------------------------------------------------------
// End-to-end: binomial_sf -> phred_qual pipeline
// ---------------------------------------------------------------------------

#[test]
fn test_pvalue_to_qual_pipeline() {
    // Simulate: variant k-mer count=50, coverage=1000, error_rate=0.001
    let p = binomial_pvalue(50, 1000, 0.001);
    let qual = phred_qual(p);

    // Should be extremely significant -> very high QUAL
    assert!(qual > 100.0, "Strong variant should have QUAL > 100, got {}", qual);
    assert!(qual <= 999.0, "QUAL should not exceed 999, got {}", qual);

    // Roundtrip check
    let recovered_p = pvalue_from_qual(qual);
    assert!(
        (p - recovered_p).abs() / p.max(1e-300) < 1e-6,
        "Roundtrip failed: p={}, qual={}, recovered={}",
        p,
        qual,
        recovered_p
    );
}

#[test]
fn test_variant_pvalue_to_qual_pipeline_with_mockdb() {
    // Full pipeline: create MockDb, compute variant p-value, convert to QUAL
    let mut db = MockDb::new(5);

    // Set up reference k-mers at 1000x
    db.set("ACGTA", 1000);
    db.set("CGTAC", 1000);
    db.set("GTACG", 1000);

    // Set up variant k-mers at 100x (true variant at ~10% VAF)
    db.set("CGTTT", 100);
    db.set("GTTTC", 100);

    let ref_path = KmerPath {
        kmers: vec![
            "ACGTA".to_string(),
            "CGTAC".to_string(),
            "GTACG".to_string(),
        ],
        is_reference: true,
    };
    let alt_path = KmerPath {
        kmers: vec![
            "ACGTA".to_string(),
            "CGTTT".to_string(),
            "GTTTC".to_string(),
        ],
        is_reference: false,
    };

    let pvalue = compute_variant_pvalue(&db, &alt_path, &ref_path, 0.001, 150);
    let qual = phred_qual(pvalue);

    // This should be a very confident variant call
    assert!(qual > 50.0, "True variant at 10% VAF should have high QUAL, got {}", qual);
}

// ---------------------------------------------------------------------------
// Confidence level affects CI width
// ---------------------------------------------------------------------------

#[test]
fn test_bootstrap_wider_ci_at_higher_confidence() {
    let mut db = MockDb::new(4);
    db.set("AAAA", 200);
    db.set("AAAC", 100);
    db.set("AACG", 100);
    db.set("AAAT", 100);
    db.set("AATG", 100);

    let ref_path = KmerPath {
        kmers: vec!["AAAA".to_string(), "AAAC".to_string(), "AACG".to_string()],
        is_reference: true,
    };
    let alt_path = KmerPath {
        kmers: vec!["AAAA".to_string(), "AAAT".to_string(), "AATG".to_string()],
        is_reference: false,
    };

    let config_90 = BootstrapConfig {
        n_replicates: 1000,
        confidence_level: 0.90,
        seed: Some(42),
    };
    let config_99 = BootstrapConfig {
        n_replicates: 1000,
        confidence_level: 0.99,
        seed: Some(42),
    };

    let cis_90 = bootstrap_confidence_intervals(
        &[ref_path.clone(), alt_path.clone()],
        &db,
        &config_90,
    );
    let cis_99 = bootstrap_confidence_intervals(&[ref_path, alt_path], &db, &config_99);

    let width_90 = cis_90[0].ci_upper - cis_90[0].ci_lower;
    let width_99 = cis_99[0].ci_upper - cis_99[0].ci_lower;

    assert!(
        width_99 > width_90,
        "99% CI ({}) should be wider than 90% CI ({})",
        width_99,
        width_90
    );
}
