// Filter pipeline integration tests.
//
// Tests the filter subcommand logic by:
// 1. Creating detection TSV output (temp files)
// 2. Parsing detection results back
// 3. Running filter logic against expected variants
// 4. Verifying found/not-found classification and filter notes

use kmerdet::filter::{self, ExpectedVariant, FilterConfig};
use kmerdet::output;
use kmerdet::variant::{VariantCall, VariantType};

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Build a VariantCall with all fields populated.
fn make_call(
    sample: &str,
    target: &str,
    chrom: &str,
    pos: u64,
    ref_a: &str,
    alt_a: &str,
    vtype: VariantType,
    rvaf: f64,
    cov: u64,
    expr: f64,
) -> VariantCall {
    VariantCall {
        sample: sample.to_string(),
        target: target.to_string(),
        variant_type: vtype,
        variant_name: format!("{}:{}/{}:{}", pos, ref_a, alt_a, pos + ref_a.len() as u64),
        rvaf,
        expression: expr,
        min_coverage: cov,
        path_score: cov,
        start_kmer_count: 200,
        ref_sequence: "ACGTACGTACGTACGT".to_string(),
        alt_sequence: format!("ACGT{}ACGT", alt_a),
        info: "vs_ref".to_string(),
        chrom: Some(chrom.to_string()),
        pos: Some(pos),
        ref_allele: Some(ref_a.to_string()),
        alt_allele: Some(alt_a.to_string()),
        pvalue: Some(1e-10),
        qual: Some(100.0),
        ci_lower: Some(rvaf * 0.8),
        ci_upper: Some(rvaf * 1.2),
    }
}

fn make_expected(
    chrom: &str,
    pos: u64,
    ref_a: &str,
    alt_a: &str,
    vtype: &str,
) -> ExpectedVariant {
    ExpectedVariant {
        chrom: chrom.to_string(),
        pos,
        ref_allele: ref_a.to_string(),
        alt_allele: alt_a.to_string(),
        variant_type: vtype.to_string(),
    }
}

// ---------------------------------------------------------------------------
// TSV round-trip tests: write calls to TSV, parse back, verify
// ---------------------------------------------------------------------------

#[test]
fn test_tsv_roundtrip_preserves_variant_calls() {
    let calls = vec![
        make_call("s1", "TP53", "chr17", 7577120, "C", "T", VariantType::Substitution, 0.15, 42, 150.0),
        make_call("s1", "NPM1", "chr5", 170837543, "T", "TCTG", VariantType::Insertion, 0.25, 100, 250.0),
        make_call("s1", "FLT3", "chr13", 28608250, "GATA", "", VariantType::Deletion, 0.08, 20, 80.0),
    ];

    // Write to temp file
    let dir = tempfile::tempdir().unwrap();
    let tsv_path = dir.path().join("detection.tsv");
    let mut file = std::fs::File::create(&tsv_path).unwrap();
    output::tsv::write(&calls, &mut file, false).unwrap();
    drop(file);

    // Parse back
    let parsed = output::parse_detection_tsv(&tsv_path).unwrap();

    assert_eq!(parsed.len(), 3, "Should parse 3 calls from TSV");

    // Verify first call
    assert_eq!(parsed[0].sample, "s1");
    assert_eq!(parsed[0].target, "TP53");
    assert_eq!(parsed[0].variant_type, VariantType::Substitution);
    assert!((parsed[0].rvaf - 0.15).abs() < 1e-4);
    assert_eq!(parsed[0].min_coverage, 42);
    assert_eq!(parsed[0].chrom.as_deref(), Some("chr17"));
    assert_eq!(parsed[0].pos, Some(7577120));
    assert_eq!(parsed[0].ref_allele.as_deref(), Some("C"));
    assert_eq!(parsed[0].alt_allele.as_deref(), Some("T"));

    // Verify second call (insertion)
    assert_eq!(parsed[1].variant_type, VariantType::Insertion);
    assert_eq!(parsed[1].alt_allele.as_deref(), Some("TCTG"));

    // Verify third call (deletion)
    assert_eq!(parsed[2].variant_type, VariantType::Deletion);
}

#[test]
fn test_tsv_roundtrip_empty_calls() {
    let calls: Vec<VariantCall> = vec![];

    let dir = tempfile::tempdir().unwrap();
    let tsv_path = dir.path().join("empty.tsv");
    let mut file = std::fs::File::create(&tsv_path).unwrap();
    output::tsv::write(&calls, &mut file, false).unwrap();
    drop(file);

    let parsed = output::parse_detection_tsv(&tsv_path).unwrap();
    assert!(parsed.is_empty(), "Empty TSV should parse to empty vec");
}

#[test]
fn test_tsv_roundtrip_preserves_optional_fields() {
    // Test that pvalue, qual, ci_lower, ci_upper survive the round-trip
    let calls = vec![make_call(
        "s1", "target1", "chr1", 100, "A", "T",
        VariantType::Substitution, 0.5, 100, 50.0,
    )];

    let dir = tempfile::tempdir().unwrap();
    let tsv_path = dir.path().join("with_pvalue.tsv");
    let mut file = std::fs::File::create(&tsv_path).unwrap();
    output::tsv::write(&calls, &mut file, false).unwrap();
    drop(file);

    let parsed = output::parse_detection_tsv(&tsv_path).unwrap();
    assert_eq!(parsed.len(), 1);

    // pvalue should be preserved (within formatting precision)
    assert!(parsed[0].pvalue.is_some(), "pvalue should survive round-trip");
    assert!(parsed[0].qual.is_some(), "qual should survive round-trip");
    assert!(parsed[0].ci_lower.is_some(), "ci_lower should survive round-trip");
    assert!(parsed[0].ci_upper.is_some(), "ci_upper should survive round-trip");
}

// ---------------------------------------------------------------------------
// Filter pipeline: write TSV, parse, filter against expected variants
// ---------------------------------------------------------------------------

#[test]
fn test_filter_pipeline_found_variants() {
    // Create detection results with known variants
    let calls = vec![
        make_call("s1", "TP53", "chr17", 100, "A", "T", VariantType::Substitution, 0.15, 50, 150.0),
        make_call("s1", "BRCA1", "chr17", 200, "GG", "AA", VariantType::Substitution, 0.10, 30, 100.0),
        make_call("s1", "NPM1", "chr5", 300, "T", "TCTG", VariantType::Insertion, 0.25, 80, 250.0),
    ];

    // Expected variants that should all be found
    let expected = vec![
        make_expected("chr17", 100, "A", "T", "Substitution"),
        make_expected("chr17", 200, "GG", "AA", "Substitution"),
        make_expected("chr5", 300, "T", "TCTG", "Insertion"),
    ];

    let config = FilterConfig {
        min_coverage: 0,
        min_vaf: 0.0,
        min_expression: 0.0,
        use_alt: false,
        types: vec![],
    };

    let results = filter::filter_results(&calls, &expected, &config).unwrap();

    assert_eq!(results.len(), 3);
    for (i, r) in results.iter().enumerate() {
        assert_eq!(
            r.found, "Found",
            "Expected variant {} should be found, got: {}",
            i,
            r.filter_notes
        );
        assert_eq!(r.filter_notes, "PASS");
    }
}

#[test]
fn test_filter_pipeline_mixed_found_and_missing() {
    let calls = vec![
        make_call("s1", "TP53", "chr17", 100, "A", "T", VariantType::Substitution, 0.15, 50, 150.0),
    ];

    let expected = vec![
        make_expected("chr17", 100, "A", "T", "Substitution"),  // found
        make_expected("chr2", 500, "C", "G", "Substitution"),   // not in calls
        make_expected("chrX", 1000, "AT", "", "Deletion"),       // not in calls
    ];

    let config = FilterConfig {
        min_coverage: 0,
        min_vaf: 0.0,
        min_expression: 0.0,
        use_alt: false,
        types: vec![],
    };

    let results = filter::filter_results(&calls, &expected, &config).unwrap();

    assert_eq!(results.len(), 3);
    assert_eq!(results[0].found, "Found");
    assert_eq!(results[0].filter_notes, "PASS");
    assert_eq!(results[0].kmer_vaf, Some(0.15));

    assert_eq!(results[1].found, "Not Found");
    assert_eq!(results[1].filter_notes, "No matching variant detected");
    assert!(results[1].kmer_vaf.is_none());

    assert_eq!(results[2].found, "Not Found");
    assert_eq!(results[2].filter_notes, "No matching variant detected");
}

#[test]
fn test_filter_pipeline_threshold_enforcement() {
    // Variant is detected but fails thresholds
    let calls = vec![
        make_call("s1", "TP53", "chr17", 100, "A", "T", VariantType::Substitution, 0.02, 3, 0.5),
    ];

    let expected = vec![
        make_expected("chr17", 100, "A", "T", "Substitution"),
    ];

    let config = FilterConfig {
        min_coverage: 10,
        min_vaf: 0.05,
        min_expression: 5.0,
        use_alt: false,
        types: vec![],
    };

    let results = filter::filter_results(&calls, &expected, &config).unwrap();

    assert_eq!(results.len(), 1);
    assert_eq!(results[0].found, "Not Found");
    // All three failures should be reported
    assert!(results[0].filter_notes.contains("coverage"), "Should report coverage failure");
    assert!(results[0].filter_notes.contains("VAF"), "Should report VAF failure");
    assert!(results[0].filter_notes.contains("expression"), "Should report expression failure");
}

#[test]
fn test_filter_pipeline_type_restriction() {
    let calls = vec![
        make_call("s1", "TP53", "chr17", 100, "A", "T", VariantType::Substitution, 0.5, 100, 50.0),
        make_call("s1", "NPM1", "chr5", 200, "T", "TCTG", VariantType::Insertion, 0.3, 80, 40.0),
    ];

    let expected = vec![
        make_expected("chr17", 100, "A", "T", "Substitution"),
        make_expected("chr5", 200, "T", "TCTG", "Insertion"),
    ];

    // Only allow Insertion type
    let config = FilterConfig {
        min_coverage: 0,
        min_vaf: 0.0,
        min_expression: 0.0,
        use_alt: false,
        types: vec!["Insertion".to_string()],
    };

    let results = filter::filter_results(&calls, &expected, &config).unwrap();

    assert_eq!(results.len(), 2);
    // Substitution should fail type check
    assert_eq!(results[0].found, "Not Found");
    assert!(results[0].filter_notes.contains("type"));
    // Insertion should pass
    assert_eq!(results[1].found, "Found");
    assert_eq!(results[1].filter_notes, "PASS");
}

#[test]
fn test_filter_pipeline_alt_mode() {
    // In alt mode, matching is done by checking if the expected alt_allele
    // appears in the call's alt_sequence field.
    let calls = vec![
        make_call("s1", "NPM1", "chr5", 170837543, "T", "TCTG", VariantType::Insertion, 0.25, 100, 250.0),
    ];

    // Expected variant with different coordinates but matching alt allele
    let expected = vec![
        make_expected("chr17", 99999, "X", "TCTG", "Insertion"),
    ];

    // Reference mode should fail (different coordinates)
    let ref_config = FilterConfig {
        use_alt: false,
        ..FilterConfig {
            min_coverage: 0,
            min_vaf: 0.0,
            min_expression: 0.0,
            use_alt: false,
            types: vec![],
        }
    };
    let ref_results = filter::filter_results(&calls, &expected, &ref_config).unwrap();
    assert_eq!(ref_results[0].found, "Not Found");

    // Alt mode should succeed (alt_sequence contains "TCTG")
    let alt_config = FilterConfig {
        use_alt: true,
        min_coverage: 0,
        min_vaf: 0.0,
        min_expression: 0.0,
        types: vec![],
    };
    let alt_results = filter::filter_results(&calls, &expected, &alt_config).unwrap();
    assert_eq!(alt_results[0].found, "Found");
    assert_eq!(alt_results[0].filter_notes, "PASS");
}

// ---------------------------------------------------------------------------
// Filter result field verification
// ---------------------------------------------------------------------------

#[test]
fn test_filter_result_fields_populated() {
    let calls = vec![
        make_call("sample_001", "FLT3_ITD", "chr13", 28608250, "A", "AGATAGAT", VariantType::Itd, 0.35, 120, 350.0),
    ];

    let expected = vec![
        make_expected("chr13", 28608250, "A", "AGATAGAT", "ITD"),
    ];

    let config = FilterConfig {
        min_coverage: 0,
        min_vaf: 0.0,
        min_expression: 0.0,
        use_alt: false,
        types: vec![],
    };

    let results = filter::filter_results(&calls, &expected, &config).unwrap();

    assert_eq!(results.len(), 1);
    let r = &results[0];

    assert_eq!(r.sample, "sample_001");
    assert_eq!(r.chrom, "chr13");
    assert_eq!(r.pos, 28608250);
    assert_eq!(r.ref_allele, "A");
    assert_eq!(r.alt_allele, "AGATAGAT");
    assert_eq!(r.variant_type, "ITD");
    assert_eq!(r.found, "Found");
    assert_eq!(r.filter_notes, "PASS");
    assert_eq!(r.kmer_vaf, Some(0.35));
    assert_eq!(r.kmer_min_coverage, Some(120));
    assert_eq!(r.kmer_expression, Some(350.0));
    assert!(r.ref_sequence.is_some());
    assert!(r.variant_sequence.is_some());
}

#[test]
fn test_filter_result_not_found_fields() {
    let calls: Vec<VariantCall> = vec![];
    let expected = vec![
        make_expected("chr1", 100, "A", "T", "Substitution"),
    ];

    let config = FilterConfig {
        min_coverage: 0,
        min_vaf: 0.0,
        min_expression: 0.0,
        use_alt: false,
        types: vec![],
    };

    let results = filter::filter_results(&calls, &expected, &config).unwrap();

    assert_eq!(results.len(), 1);
    let r = &results[0];

    assert_eq!(r.found, "Not Found");
    assert_eq!(r.filter_notes, "No matching variant detected");
    assert!(r.kmer_vaf.is_none());
    assert!(r.kmer_min_coverage.is_none());
    assert!(r.kmer_expression.is_none());
    assert!(r.ref_sequence.is_none());
    assert!(r.variant_sequence.is_none());
    // Fields from expected variant should still be populated
    assert_eq!(r.chrom, "chr1");
    assert_eq!(r.pos, 100);
    assert_eq!(r.ref_allele, "A");
    assert_eq!(r.alt_allele, "T");
    assert_eq!(r.variant_type, "Substitution");
}

// ---------------------------------------------------------------------------
// Chromosome normalization in filter pipeline
// ---------------------------------------------------------------------------

#[test]
fn test_filter_pipeline_chr_normalization() {
    // Call uses "chr1", expected uses "1" (no prefix)
    let calls = vec![
        make_call("s1", "TP53", "chr1", 100, "A", "T", VariantType::Substitution, 0.5, 100, 50.0),
    ];
    let expected = vec![
        make_expected("1", 100, "A", "T", "Substitution"),
    ];

    let config = FilterConfig {
        min_coverage: 0,
        min_vaf: 0.0,
        min_expression: 0.0,
        use_alt: false,
        types: vec![],
    };

    let results = filter::filter_results(&calls, &expected, &config).unwrap();
    assert_eq!(results[0].found, "Found", "chr normalization should match chr1 == 1");
}
