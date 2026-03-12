// Comprehensive tests for the tumor-informed filter engine.

use kmerdet::filter::{self, ExpectedVariant, FilterConfig};
use kmerdet::variant::{VariantCall, VariantType};

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn make_call(
    chrom: &str,
    pos: u64,
    ref_a: &str,
    alt_a: &str,
    vtype: VariantType,
    vaf: f64,
    cov: u64,
    expr: f64,
) -> VariantCall {
    VariantCall {
        sample: "test_sample".into(),
        target: "test_target".into(),
        variant_type: vtype,
        variant_name: format!("{}:{}/{}:{}", pos, ref_a, alt_a, pos + ref_a.len() as u64),
        rvaf: vaf,
        expression: expr,
        min_coverage: cov,
        start_kmer_count: 100,
        ref_sequence: "ACGTACGT".into(),
        alt_sequence: format!("ACGT{}ACGT", alt_a),
        info: "vs_ref".into(),
        chrom: Some(chrom.into()),
        pos: Some(pos),
        ref_allele: Some(ref_a.into()),
        alt_allele: Some(alt_a.into()),
        ci_lower: None,
        ci_upper: None,
    }
}

fn make_expected(chrom: &str, pos: u64, ref_a: &str, alt_a: &str, vtype: &str) -> ExpectedVariant {
    ExpectedVariant {
        chrom: chrom.into(),
        pos,
        ref_allele: ref_a.into(),
        alt_allele: alt_a.into(),
        variant_type: vtype.into(),
    }
}

fn default_config() -> FilterConfig {
    FilterConfig {
        min_coverage: 0,
        min_vaf: 0.0,
        min_expression: 0.0,
        use_alt: false,
        types: vec![],
    }
}

// ---------------------------------------------------------------------------
// Reference mode tests
// ---------------------------------------------------------------------------

#[test]
fn test_reference_mode_exact_match() {
    let calls = vec![make_call(
        "chr1",
        12345,
        "A",
        "T",
        VariantType::Substitution,
        0.5,
        100,
        50.0,
    )];
    let expected = make_expected("chr1", 12345, "A", "T", "Substitution");

    let result = filter::reference_mode::find_match(&expected, &calls);
    assert!(result.is_some());
    let matched = result.unwrap();
    assert_eq!(matched.pos, Some(12345));
    assert_eq!(matched.ref_allele.as_deref(), Some("A"));
    assert_eq!(matched.alt_allele.as_deref(), Some("T"));
}

#[test]
fn test_reference_mode_no_match() {
    let calls = vec![make_call(
        "chr1",
        12345,
        "A",
        "T",
        VariantType::Substitution,
        0.5,
        100,
        50.0,
    )];
    // Different position
    let expected = make_expected("chr1", 99999, "A", "T", "Substitution");
    assert!(filter::reference_mode::find_match(&expected, &calls).is_none());

    // Different ref allele
    let expected2 = make_expected("chr1", 12345, "G", "T", "Substitution");
    assert!(filter::reference_mode::find_match(&expected2, &calls).is_none());

    // Different alt allele
    let expected3 = make_expected("chr1", 12345, "A", "C", "Substitution");
    assert!(filter::reference_mode::find_match(&expected3, &calls).is_none());

    // Different chromosome
    let expected4 = make_expected("chr2", 12345, "A", "T", "Substitution");
    assert!(filter::reference_mode::find_match(&expected4, &calls).is_none());

    // Different type
    let expected5 = make_expected("chr1", 12345, "A", "T", "Insertion");
    assert!(filter::reference_mode::find_match(&expected5, &calls).is_none());
}

#[test]
fn test_reference_mode_chr_normalization() {
    // Call has "1", expected has "chr1"
    let calls = vec![make_call(
        "1",
        500,
        "A",
        "G",
        VariantType::Substitution,
        0.3,
        50,
        20.0,
    )];
    let expected = make_expected("chr1", 500, "A", "G", "Substitution");
    assert!(filter::reference_mode::find_match(&expected, &calls).is_some());

    // Call has "chr1", expected has "1"
    let calls2 = vec![make_call(
        "chr1",
        500,
        "A",
        "G",
        VariantType::Substitution,
        0.3,
        50,
        20.0,
    )];
    let expected2 = make_expected("1", 500, "A", "G", "Substitution");
    assert!(filter::reference_mode::find_match(&expected2, &calls2).is_some());

    // Both have "chr" prefix
    let expected3 = make_expected("chr1", 500, "A", "G", "Substitution");
    assert!(filter::reference_mode::find_match(&expected3, &calls2).is_some());

    // Neither has "chr" prefix
    let expected4 = make_expected("1", 500, "A", "G", "Substitution");
    assert!(filter::reference_mode::find_match(&expected4, &calls).is_some());

    // chrX vs X
    let calls_x = vec![make_call(
        "X",
        1000,
        "C",
        "T",
        VariantType::Substitution,
        0.2,
        30,
        10.0,
    )];
    let expected_x = make_expected("chrX", 1000, "C", "T", "Substitution");
    assert!(filter::reference_mode::find_match(&expected_x, &calls_x).is_some());
}

#[test]
fn test_reference_mode_type_case_insensitive() {
    let calls = vec![make_call(
        "chr1",
        100,
        "A",
        "T",
        VariantType::Substitution,
        0.5,
        100,
        50.0,
    )];

    // lowercase expected type
    let expected = make_expected("chr1", 100, "A", "T", "substitution");
    assert!(filter::reference_mode::find_match(&expected, &calls).is_some());

    // UPPERCASE expected type
    let expected2 = make_expected("chr1", 100, "A", "T", "SUBSTITUTION");
    assert!(filter::reference_mode::find_match(&expected2, &calls).is_some());

    // MiXeD case
    let expected3 = make_expected("chr1", 100, "A", "T", "SubStItUtIoN");
    assert!(filter::reference_mode::find_match(&expected3, &calls).is_some());
}

// ---------------------------------------------------------------------------
// Alt mode tests
// ---------------------------------------------------------------------------

#[test]
fn test_alt_mode_match() {
    let calls = vec![make_call(
        "chr1",
        100,
        "A",
        "TGCA",
        VariantType::Insertion,
        0.4,
        80,
        30.0,
    )];
    // alt_sequence is "ACGTTGCAACGT" — contains "TGCA"
    let expected = make_expected("chr1", 100, "A", "TGCA", "Insertion");
    assert!(filter::alt_mode::find_match(&expected, &calls).is_some());
}

#[test]
fn test_alt_mode_substring_match() {
    // Alt mode matches if expected alt_allele is a substring of call's alt_sequence
    let calls = vec![make_call(
        "chr1",
        100,
        "A",
        "TGCAATG",
        VariantType::Insertion,
        0.4,
        80,
        30.0,
    )];
    // alt_sequence is "ACGTTGCAATGACGT", search for just "GCAAT" substring
    let expected = make_expected("chr5", 999, "X", "GCAAT", "Insertion");
    assert!(filter::alt_mode::find_match(&expected, &calls).is_some());
}

#[test]
fn test_alt_mode_no_match() {
    let calls = vec![make_call(
        "chr1",
        100,
        "A",
        "T",
        VariantType::Substitution,
        0.5,
        100,
        50.0,
    )];
    // alt_sequence is "ACGTTACGT" — does not contain "ZZZZZ"
    let expected = make_expected("chr1", 100, "A", "ZZZZZ", "Substitution");
    assert!(filter::alt_mode::find_match(&expected, &calls).is_none());
}

// ---------------------------------------------------------------------------
// filter_results tests
// ---------------------------------------------------------------------------

#[test]
fn test_filter_passes_all_conditions() {
    let calls = vec![make_call(
        "chr1",
        12345,
        "A",
        "T",
        VariantType::Substitution,
        0.5,
        100,
        50.0,
    )];
    let expected = vec![make_expected("chr1", 12345, "A", "T", "Substitution")];
    let config = FilterConfig {
        min_coverage: 10,
        min_vaf: 0.01,
        min_expression: 1.0,
        use_alt: false,
        types: vec![],
    };

    let results = filter::filter_results(&calls, &expected, &config).unwrap();
    assert_eq!(results.len(), 1);
    assert_eq!(results[0].found, "Found");
    assert_eq!(results[0].filter_notes, "PASS");
    assert_eq!(results[0].sample, "test_sample");
    assert_eq!(results[0].chrom, "chr1");
    assert_eq!(results[0].pos, 12345);
    assert_eq!(results[0].ref_allele, "A");
    assert_eq!(results[0].alt_allele, "T");
    assert_eq!(results[0].variant_type, "Substitution");
    assert_eq!(results[0].kmer_vaf, Some(0.5));
    assert_eq!(results[0].kmer_min_coverage, Some(100));
    assert_eq!(results[0].kmer_expression, Some(50.0));
    assert!(results[0].ref_sequence.is_some());
    assert!(results[0].variant_sequence.is_some());
}

#[test]
fn test_filter_fails_coverage() {
    let calls = vec![make_call(
        "chr1",
        12345,
        "A",
        "T",
        VariantType::Substitution,
        0.5,
        5, // below threshold
        50.0,
    )];
    let expected = vec![make_expected("chr1", 12345, "A", "T", "Substitution")];
    let config = FilterConfig {
        min_coverage: 10,
        min_vaf: 0.0,
        min_expression: 0.0,
        use_alt: false,
        types: vec![],
    };

    let results = filter::filter_results(&calls, &expected, &config).unwrap();
    assert_eq!(results.len(), 1);
    assert_eq!(results[0].found, "Not Found");
    assert!(results[0].filter_notes.contains("coverage"));
    assert!(results[0].filter_notes.contains("5"));
    assert!(results[0].filter_notes.contains("10"));
}

#[test]
fn test_filter_fails_vaf() {
    let calls = vec![make_call(
        "chr1",
        12345,
        "A",
        "T",
        VariantType::Substitution,
        0.001, // below threshold
        100,
        50.0,
    )];
    let expected = vec![make_expected("chr1", 12345, "A", "T", "Substitution")];
    let config = FilterConfig {
        min_coverage: 0,
        min_vaf: 0.01,
        min_expression: 0.0,
        use_alt: false,
        types: vec![],
    };

    let results = filter::filter_results(&calls, &expected, &config).unwrap();
    assert_eq!(results.len(), 1);
    assert_eq!(results[0].found, "Not Found");
    assert!(results[0].filter_notes.contains("VAF"));
}

#[test]
fn test_filter_fails_expression() {
    let calls = vec![make_call(
        "chr1",
        12345,
        "A",
        "T",
        VariantType::Substitution,
        0.5,
        100,
        0.5, // below threshold
    )];
    let expected = vec![make_expected("chr1", 12345, "A", "T", "Substitution")];
    let config = FilterConfig {
        min_coverage: 0,
        min_vaf: 0.0,
        min_expression: 5.0,
        use_alt: false,
        types: vec![],
    };

    let results = filter::filter_results(&calls, &expected, &config).unwrap();
    assert_eq!(results.len(), 1);
    assert_eq!(results[0].found, "Not Found");
    assert!(results[0].filter_notes.contains("expression"));
}

#[test]
fn test_filter_fails_type() {
    let calls = vec![make_call(
        "chr1",
        12345,
        "A",
        "T",
        VariantType::Substitution,
        0.5,
        100,
        50.0,
    )];
    let expected = vec![make_expected("chr1", 12345, "A", "T", "Substitution")];
    let config = FilterConfig {
        min_coverage: 0,
        min_vaf: 0.0,
        min_expression: 0.0,
        use_alt: false,
        types: vec!["Insertion".into(), "Deletion".into()], // Substitution not allowed
    };

    let results = filter::filter_results(&calls, &expected, &config).unwrap();
    assert_eq!(results.len(), 1);
    assert_eq!(results[0].found, "Not Found");
    assert!(results[0].filter_notes.contains("type"));
    assert!(results[0].filter_notes.contains("not in allowed types"));
}

#[test]
fn test_filter_multiple_expected() {
    let calls = vec![
        make_call(
            "chr1",
            100,
            "A",
            "T",
            VariantType::Substitution,
            0.5,
            100,
            50.0,
        ),
        make_call(
            "chr2",
            200,
            "GG",
            "AA",
            VariantType::Substitution,
            0.3,
            80,
            30.0,
        ),
    ];
    let expected = vec![
        make_expected("chr1", 100, "A", "T", "Substitution"),     // match
        make_expected("chr3", 300, "C", "G", "Substitution"),     // no match
        make_expected("chr2", 200, "GG", "AA", "Substitution"),   // match
    ];
    let config = default_config();

    let results = filter::filter_results(&calls, &expected, &config).unwrap();
    assert_eq!(results.len(), 3);

    assert_eq!(results[0].found, "Found");
    assert_eq!(results[0].chrom, "chr1");
    assert_eq!(results[0].filter_notes, "PASS");

    assert_eq!(results[1].found, "Not Found");
    assert_eq!(results[1].chrom, "chr3");
    assert_eq!(results[1].filter_notes, "No matching variant detected");

    assert_eq!(results[2].found, "Found");
    assert_eq!(results[2].chrom, "chr2");
    assert_eq!(results[2].filter_notes, "PASS");
}

#[test]
fn test_filter_empty_calls() {
    let calls: Vec<VariantCall> = vec![];
    let expected = vec![
        make_expected("chr1", 100, "A", "T", "Substitution"),
        make_expected("chr2", 200, "C", "G", "Deletion"),
    ];
    let config = default_config();

    let results = filter::filter_results(&calls, &expected, &config).unwrap();
    assert_eq!(results.len(), 2);
    for r in &results {
        assert_eq!(r.found, "Not Found");
        assert_eq!(r.filter_notes, "No matching variant detected");
        assert!(r.kmer_vaf.is_none());
        assert!(r.kmer_min_coverage.is_none());
        assert!(r.kmer_expression.is_none());
        assert!(r.ref_sequence.is_none());
        assert!(r.variant_sequence.is_none());
        assert_eq!(r.sample, "");
    }
}

#[test]
fn test_filter_empty_expected() {
    let calls = vec![make_call(
        "chr1",
        100,
        "A",
        "T",
        VariantType::Substitution,
        0.5,
        100,
        50.0,
    )];
    let expected: Vec<ExpectedVariant> = vec![];
    let config = default_config();

    let results = filter::filter_results(&calls, &expected, &config).unwrap();
    assert!(results.is_empty());
}

#[test]
fn test_filter_use_alt_mode() {
    // Build a call where coordinate fields don't match expected,
    // but alt_sequence contains the expected alt_allele.
    let calls = vec![make_call(
        "chr5",
        9999,
        "X",
        "INSERTSEQ",
        VariantType::Insertion,
        0.4,
        60,
        25.0,
    )];
    // Expected has totally different chrom/pos — reference mode would fail
    let expected = vec![make_expected("chr1", 100, "A", "INSERTSEQ", "Insertion")];

    // Reference mode should NOT find it (different chrom/pos/ref/alt coords)
    let ref_config = FilterConfig {
        use_alt: false,
        ..default_config()
    };
    let ref_results = filter::filter_results(&calls, &expected, &ref_config).unwrap();
    assert_eq!(ref_results[0].found, "Not Found");

    // Alt mode SHOULD find it (alt_sequence contains "INSERTSEQ")
    let alt_config = FilterConfig {
        use_alt: true,
        ..default_config()
    };
    let alt_results = filter::filter_results(&calls, &expected, &alt_config).unwrap();
    assert_eq!(alt_results[0].found, "Found");
    assert_eq!(alt_results[0].filter_notes, "PASS");
    assert_eq!(alt_results[0].kmer_vaf, Some(0.4));
    assert_eq!(alt_results[0].kmer_min_coverage, Some(60));
}

// ---------------------------------------------------------------------------
// Edge cases and combined failures
// ---------------------------------------------------------------------------

#[test]
fn test_filter_multiple_failures() {
    let calls = vec![make_call(
        "chr1",
        12345,
        "A",
        "T",
        VariantType::Substitution,
        0.001, // below VAF
        3,     // below coverage
        0.1,   // below expression
    )];
    let expected = vec![make_expected("chr1", 12345, "A", "T", "Substitution")];
    let config = FilterConfig {
        min_coverage: 10,
        min_vaf: 0.01,
        min_expression: 5.0,
        use_alt: false,
        types: vec![],
    };

    let results = filter::filter_results(&calls, &expected, &config).unwrap();
    assert_eq!(results[0].found, "Not Found");
    // All three failures should be noted
    assert!(results[0].filter_notes.contains("coverage"));
    assert!(results[0].filter_notes.contains("VAF"));
    assert!(results[0].filter_notes.contains("expression"));
}
