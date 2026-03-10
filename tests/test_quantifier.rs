// Path quantification (NNLS) and variant clustering tests.

mod common;

use common::MockDb;
use kmerdet::sequence::path::KmerPath;
use kmerdet::variant::cluster::cluster_variants;
use kmerdet::variant::quantifier::quantify;
use kmerdet::variant::{VariantCall, VariantType};

// ---------------------------------------------------------------------------
// Helper to build a VariantCall with a given variant_name (positions encoded).
// ---------------------------------------------------------------------------
fn make_call(variant_name: &str, start: usize, end: usize) -> VariantCall {
    let name = format!("{}:A/T:{}", start, end);
    let name = if variant_name.is_empty() {
        name
    } else {
        variant_name.to_string()
    };
    VariantCall {
        sample: "test".to_string(),
        target: "target1".to_string(),
        variant_type: VariantType::Substitution,
        variant_name: name,
        rvaf: 0.5,
        expression: 100.0,
        min_coverage: 50,
        start_kmer_count: 100,
        ref_sequence: "ACGT".to_string(),
        alt_sequence: "TCGT".to_string(),
        info: "vs_ref".to_string(),
        chrom: None,
        pos: None,
        ref_allele: Some("A".to_string()),
        alt_allele: Some("T".to_string()),
    }
}

fn make_call_at(start: usize, end: usize) -> VariantCall {
    make_call("", start, end)
}

// ===========================================================================
// Quantifier tests
// ===========================================================================

#[test]
fn test_quantify_single_path() {
    // One path (reference only). coefficient should equal the mean k-mer count.
    // rVAF should be 1.0.
    // Path: ACGT -> CGTA -> GTAC (k=4, sequence "ACGTAC")
    let mut db = MockDb::new(4);
    db.set("ACGT", 100);
    db.set("CGTA", 100);
    db.set("GTAC", 100);

    let path = KmerPath {
        kmers: vec![
            "ACGT".to_string(),
            "CGTA".to_string(),
            "GTAC".to_string(),
        ],
        is_reference: true,
    };

    let result = quantify(&[path], &db);

    assert_eq!(result.coefficients.len(), 1);
    assert_eq!(result.rvafs.len(), 1);
    assert_eq!(result.min_coverages.len(), 1);

    // Single path: coefficient should be approximately 100.
    assert!(
        (result.coefficients[0] - 100.0).abs() < 1.0,
        "coefficient should be ~100, got {}",
        result.coefficients[0]
    );
    // rVAF of the single path must be 1.0.
    assert!(
        (result.rvafs[0] - 1.0).abs() < 1e-9,
        "rVAF should be 1.0, got {}",
        result.rvafs[0]
    );
    // min_coverage should be 100.
    assert_eq!(result.min_coverages[0], 100);
}

#[test]
fn test_quantify_two_paths_equal() {
    // Reference + variant path with equal counts.
    // Ref path:  AAAA -> AAAC -> AACG (k=4)
    // Alt path:  AAAA -> AAAT -> AATG (k=4)
    // Shared k-mer: AAAA (appears in both paths)
    // Ref-only:  AAAC, AACG
    // Alt-only:  AAAT, AATG

    let mut db = MockDb::new(4);
    // Shared k-mer at count 200 (sum of both paths contributing 100 each).
    db.set("AAAA", 200);
    // Ref-only k-mers at count 100.
    db.set("AAAC", 100);
    db.set("AACG", 100);
    // Alt-only k-mers at count 100.
    db.set("AAAT", 100);
    db.set("AATG", 100);

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

    let result = quantify(&[ref_path, alt_path], &db);

    assert_eq!(result.coefficients.len(), 2);
    assert_eq!(result.rvafs.len(), 2);

    // Both paths should have roughly equal coefficients (~100 each).
    assert!(
        (result.coefficients[0] - 100.0).abs() < 1.0,
        "ref coefficient should be ~100, got {}",
        result.coefficients[0]
    );
    assert!(
        (result.coefficients[1] - 100.0).abs() < 1.0,
        "alt coefficient should be ~100, got {}",
        result.coefficients[1]
    );

    // rVAF should be ~0.5 each.
    assert!(
        (result.rvafs[0] - 0.5).abs() < 0.01,
        "ref rVAF should be ~0.5, got {}",
        result.rvafs[0]
    );
    assert!(
        (result.rvafs[1] - 0.5).abs() < 0.01,
        "alt rVAF should be ~0.5, got {}",
        result.rvafs[1]
    );
}

#[test]
fn test_quantify_90_10_split() {
    // 90% ref, 10% variant.
    // Ref path:  CCCC -> CCCA -> CCAG (k=4)
    // Alt path:  CCCC -> CCCT -> CCTG (k=4)

    let mut db = MockDb::new(4);
    // Shared k-mer at count 100 (90 ref + 10 alt).
    db.set("CCCC", 100);
    // Ref-only k-mers at count 90.
    db.set("CCCA", 90);
    db.set("CCAG", 90);
    // Alt-only k-mers at count 10.
    db.set("CCCT", 10);
    db.set("CCTG", 10);

    let ref_path = KmerPath {
        kmers: vec![
            "CCCC".to_string(),
            "CCCA".to_string(),
            "CCAG".to_string(),
        ],
        is_reference: true,
    };
    let alt_path = KmerPath {
        kmers: vec![
            "CCCC".to_string(),
            "CCCT".to_string(),
            "CCTG".to_string(),
        ],
        is_reference: false,
    };

    let result = quantify(&[ref_path, alt_path], &db);

    // rVAFs should be approximately 0.9 and 0.1.
    assert!(
        (result.rvafs[0] - 0.9).abs() < 0.05,
        "ref rVAF should be ~0.9, got {}",
        result.rvafs[0]
    );
    assert!(
        (result.rvafs[1] - 0.1).abs() < 0.05,
        "alt rVAF should be ~0.1, got {}",
        result.rvafs[1]
    );
}

#[test]
fn test_quantify_negative_clipping() {
    // Setup where unconstrained least squares would produce a negative coefficient.
    // Three paths, but path 3 has all k-mers at zero count while overlapping with
    // path 1. This can force a negative coefficient in unconstrained solve.
    //
    // Ref path:  GGGG -> GGGA -> GGAC (k=4)
    // Alt path:  GGGG -> GGGT -> GGTC (k=4)
    // Noise path: GGGG -> GGGA -> GGAT (k=4) -- shares kmers with ref, counts = 0

    let mut db = MockDb::new(4);
    db.set("GGGG", 100);
    db.set("GGGA", 100);
    db.set("GGAC", 100);
    db.set("GGGT", 0);
    db.set("GGTC", 0);
    db.set("GGAT", 0);

    let ref_path = KmerPath {
        kmers: vec![
            "GGGG".to_string(),
            "GGGA".to_string(),
            "GGAC".to_string(),
        ],
        is_reference: true,
    };
    let alt_path = KmerPath {
        kmers: vec![
            "GGGG".to_string(),
            "GGGT".to_string(),
            "GGTC".to_string(),
        ],
        is_reference: false,
    };
    let noise_path = KmerPath {
        kmers: vec![
            "GGGG".to_string(),
            "GGGA".to_string(),
            "GGAT".to_string(),
        ],
        is_reference: false,
    };

    let result = quantify(&[ref_path, alt_path, noise_path], &db);

    // All coefficients must be >= 0 (NNLS constraint).
    for (i, &c) in result.coefficients.iter().enumerate() {
        assert!(
            c >= 0.0,
            "coefficient[{}] should be >= 0, got {}",
            i, c
        );
    }

    // All rVAFs must also be >= 0.
    for (i, &r) in result.rvafs.iter().enumerate() {
        assert!(
            r >= 0.0,
            "rVAF[{}] should be >= 0, got {}",
            i, r
        );
    }

    // The ref path should carry most of the expression.
    assert!(
        result.coefficients[0] > 50.0,
        "ref coefficient should be dominant, got {}",
        result.coefficients[0]
    );
}

#[test]
fn test_quantify_min_coverage() {
    // Verify min_coverage is the minimum k-mer count along each path.
    let mut db = MockDb::new(4);
    db.set("ACGT", 200);
    db.set("CGTA", 50); // This is the minimum for ref path
    db.set("GTAC", 300);
    db.set("CGTG", 10); // This is the minimum for alt path
    db.set("GTGA", 500);

    let ref_path = KmerPath {
        kmers: vec![
            "ACGT".to_string(),
            "CGTA".to_string(),
            "GTAC".to_string(),
        ],
        is_reference: true,
    };
    let alt_path = KmerPath {
        kmers: vec![
            "ACGT".to_string(),
            "CGTG".to_string(),
            "GTGA".to_string(),
        ],
        is_reference: false,
    };

    let result = quantify(&[ref_path, alt_path], &db);

    assert_eq!(
        result.min_coverages[0], 50,
        "ref min_coverage should be 50, got {}",
        result.min_coverages[0]
    );
    assert_eq!(
        result.min_coverages[1], 10,
        "alt min_coverage should be 10, got {}",
        result.min_coverages[1]
    );
}

#[test]
fn test_quantify_itd_contribution() {
    // Path where a k-mer appears twice (internal tandem duplication).
    // ITD path: ACGT -> CGTA -> GTAC -> ACGT -> CGTA -> GTAC
    // The k-mer ACGT appears twice, contribution should be 2.
    let mut db = MockDb::new(4);
    db.set("ACGT", 200);
    db.set("CGTA", 200);
    db.set("GTAC", 200);

    let itd_path = KmerPath {
        kmers: vec![
            "ACGT".to_string(),
            "CGTA".to_string(),
            "GTAC".to_string(),
            "ACGT".to_string(),
            "CGTA".to_string(),
            "GTAC".to_string(),
        ],
        is_reference: false,
    };

    let result = quantify(&[itd_path], &db);

    // With one path and contribution=2 for each k-mer,
    // the equation is 2*c = 200 for each k-mer, so c = 100.
    assert!(
        (result.coefficients[0] - 100.0).abs() < 1.0,
        "ITD coefficient should be ~100 (counts/2), got {}",
        result.coefficients[0]
    );
    assert!(
        (result.rvafs[0] - 1.0).abs() < 1e-9,
        "Single path rVAF should be 1.0, got {}",
        result.rvafs[0]
    );
    // min_coverage = min of query for each k-mer string in the path.
    // All are 200, so min_coverage = 200.
    assert_eq!(result.min_coverages[0], 200);
}

#[test]
fn test_quantify_empty_paths() {
    let db = MockDb::new(4);
    let result = quantify(&[], &db);

    assert!(result.coefficients.is_empty());
    assert!(result.rvafs.is_empty());
    assert!(result.min_coverages.is_empty());
}

// ===========================================================================
// Clustering tests
// ===========================================================================

#[test]
fn test_cluster_no_overlap() {
    // Three non-overlapping variants -> three separate clusters.
    let calls = vec![
        make_call_at(10, 15),
        make_call_at(20, 25),
        make_call_at(30, 35),
    ];

    let clusters = cluster_variants(&calls);

    assert_eq!(clusters.len(), 3, "expected 3 clusters, got {}", clusters.len());
    assert_eq!(clusters[0].start, 10);
    assert_eq!(clusters[0].end, 15);
    assert_eq!(clusters[0].calls.len(), 1);
    assert_eq!(clusters[1].start, 20);
    assert_eq!(clusters[1].end, 25);
    assert_eq!(clusters[1].calls.len(), 1);
    assert_eq!(clusters[2].start, 30);
    assert_eq!(clusters[2].end, 35);
    assert_eq!(clusters[2].calls.len(), 1);
}

#[test]
fn test_cluster_two_overlapping() {
    // Two variants that overlap -> one cluster with both.
    // Region [10, 20] overlaps with [15, 25].
    let calls = vec![
        make_call_at(10, 20),
        make_call_at(15, 25),
    ];

    let clusters = cluster_variants(&calls);

    assert_eq!(clusters.len(), 1, "expected 1 cluster, got {}", clusters.len());
    assert_eq!(clusters[0].calls.len(), 2);
    assert_eq!(clusters[0].start, 10);
    assert_eq!(clusters[0].end, 25);
}

#[test]
fn test_cluster_chain_overlap() {
    // A overlaps B, B overlaps C, but A doesn't directly overlap C.
    // A: [10, 20], B: [18, 28], C: [26, 36]
    // A overlaps B (18 <= 20), B overlaps C (26 <= 28), but A doesn't overlap C (26 > 20).
    // Transitive: all three should be in one cluster.
    let calls = vec![
        make_call_at(10, 20),
        make_call_at(18, 28),
        make_call_at(26, 36),
    ];

    let clusters = cluster_variants(&calls);

    assert_eq!(clusters.len(), 1, "expected 1 cluster (transitive), got {}", clusters.len());
    assert_eq!(clusters[0].calls.len(), 3);
    assert_eq!(clusters[0].start, 10);
    assert_eq!(clusters[0].end, 36);
}

#[test]
fn test_cluster_empty() {
    let clusters = cluster_variants(&[]);
    assert!(clusters.is_empty(), "empty input should produce empty output");
}
