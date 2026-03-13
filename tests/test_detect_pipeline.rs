// End-to-end detection pipeline integration tests.
//
// These tests exercise the full detection pipeline using MockDb:
//   walk -> build_graph -> find_alternative_paths -> classify -> quantify
//
// Each test sets up a MockDb with known k-mer counts, runs the pipeline,
// and verifies the final variant calls have correct types and reasonable rVAFs.

mod common;

use common::MockDb;
use kmerdet::graph::builder::build_graph;
use kmerdet::graph::pathfind::find_alternative_paths;
use kmerdet::sequence::path::KmerPath;
use kmerdet::variant::classifier::classify;
use kmerdet::variant::quantifier::quantify;
use kmerdet::variant::VariantType;
use kmerdet::walker::{walk, WalkerConfig};

/// Decompose a DNA sequence into overlapping k-mers of length k.
fn decompose(seq: &str, k: usize) -> Vec<String> {
    seq.as_bytes()
        .windows(k)
        .map(|w| std::str::from_utf8(w).unwrap().to_string())
        .collect()
}

/// Run the full detection pipeline on a MockDb and return the paths, classifications,
/// and quantification results.
struct PipelineResult {
    paths: Vec<KmerPath>,
    classifications: Vec<(VariantType, String, String)>, // (type, ref_allele, alt_allele)
    rvafs: Vec<f64>,
    min_coverages: Vec<u64>,
}

fn run_pipeline(db: &MockDb, ref_seq: &str, k: usize) -> PipelineResult {
    let ref_kmers = decompose(ref_seq, k);
    let config = WalkerConfig::default();

    // Step 1: Walk k-mers
    let walk_result = walk(db, &ref_kmers, &config);

    // Step 2: Build graph
    let graph = build_graph(&walk_result, &ref_kmers);

    // Step 3: Find paths
    let paths = find_alternative_paths(&graph);

    if paths.is_empty() {
        return PipelineResult {
            paths: vec![],
            classifications: vec![],
            rvafs: vec![],
            min_coverages: vec![],
        };
    }

    // Step 4: Classify each alt path against the reference path
    let ref_path = &paths[0];
    let mut classifications = Vec::new();
    for path in &paths {
        let class = classify(ref_path, path, k as u8);
        classifications.push((class.variant_type, class.ref_allele, class.alt_allele));
    }

    // Step 5: Quantify
    let quant = quantify(&paths, db);

    PipelineResult {
        paths: paths.clone(),
        classifications,
        rvafs: quant.rvafs,
        min_coverages: quant.min_coverages,
    }
}

// ============================================================
// Scenario A: SNV detection
// ============================================================

#[test]
fn test_pipeline_snv_detection() {
    // Reference: ATCGCAG (k=3), k-mers: ATC, TCG, CGC, GCA, CAG
    // Variant: single base change at position 3 (G->A), producing ATCACAG
    //   Variant k-mers: TCA, CAC, ACA (branching at ATC -> TCA instead of TCG)
    //   The alt path reconnects: ACA suffix=CA matches CAG prefix=CA, so ACA->CAG.
    //   Full alt path: ATC->TCA->CAC->ACA->CAG (same length as ref path).
    //
    // This is a carefully chosen reference to avoid suffix shortcuts:
    //   First ref k-mer ATC suffix=TC, Last ref k-mer CAG suffix=AG (different).
    let k = 3;
    let ref_seq = "ATCGCAG";

    let mut db = MockDb::new(k as u8);
    // Reference k-mers at 1000x coverage
    db.set_sequence(ref_seq, 1000);
    // Variant k-mers at moderate frequency
    // Threshold: 0.05 * (1000 + 100) = 55, so 100 passes easily
    db.set("TCA", 100);
    db.set("CAC", 100);
    db.set("ACA", 100);

    let result = run_pipeline(&db, ref_seq, k);

    // Should find at least 2 paths: reference + variant
    assert!(
        result.paths.len() >= 2,
        "SNV pipeline should find at least 2 paths, found {}",
        result.paths.len()
    );

    // First path is reference
    assert_eq!(result.classifications[0].0, VariantType::Reference);
    assert_eq!(result.paths[0].to_sequence(), ref_seq);

    // Find the substitution variant path
    let alt_idx = result
        .classifications
        .iter()
        .position(|(vt, _, _)| *vt == VariantType::Substitution)
        .expect("Should find a Substitution variant");

    // The variant path should produce a different sequence of the same length
    let alt_seq = result.paths[alt_idx].to_sequence();
    assert_ne!(alt_seq, ref_seq, "Alt path should differ from reference");
    assert_eq!(alt_seq.len(), ref_seq.len(), "SNV should produce same-length sequence");
    assert_eq!(alt_seq, "ATCACAG", "Alt path should be the SNV variant sequence");

    // Reference rVAF should be dominant
    assert!(
        result.rvafs[0] > 0.5,
        "Reference rVAF should be > 0.5, got {}",
        result.rvafs[0]
    );

    // rVAFs should sum to approximately 1.0
    let rvaf_sum: f64 = result.rvafs.iter().sum();
    assert!(
        (rvaf_sum - 1.0).abs() < 0.01,
        "rVAFs should sum to ~1.0, got {}",
        rvaf_sum
    );
}

// ============================================================
// Scenario B: Insertion detection
// ============================================================

#[test]
fn test_pipeline_insertion_detection() {
    // Reference: ACGTAG (k=3), k-mers: ACG, CGT, GTA, TAG
    // Insertion of T between positions 3 and 4: ACGTTAG
    //   Variant k-mers: GTT, TTA (CGT->GTT->TTA->TAG)
    let k = 3;
    let ref_seq = "ACGTAG";

    let mut db = MockDb::new(k as u8);
    db.set_sequence(ref_seq, 500);
    // Insertion variant at moderate frequency
    db.set("GTT", 100);
    db.set("TTA", 100);

    let result = run_pipeline(&db, ref_seq, k);

    assert!(
        result.paths.len() >= 2,
        "Insertion pipeline should find at least 2 paths, found {}",
        result.paths.len()
    );

    // Reference path
    assert_eq!(result.classifications[0].0, VariantType::Reference);
    assert_eq!(result.paths[0].to_sequence(), ref_seq);

    // Find the insertion path
    let alt_idx = result
        .classifications
        .iter()
        .position(|(vt, _, _)| *vt == VariantType::Insertion)
        .expect("Should find an Insertion variant");

    let alt_seq = result.paths[alt_idx].to_sequence();
    assert_eq!(
        alt_seq, "ACGTTAG",
        "Insertion path should have 1 extra base"
    );

    // Alt sequence should be longer than reference (insertion adds bases)
    assert!(
        alt_seq.len() > ref_seq.len(),
        "Insertion alt sequence ({}) should be longer than reference ({})",
        alt_seq.len(),
        ref_seq.len()
    );

    // Variant rVAF should be positive
    assert!(
        result.rvafs[alt_idx] > 0.0,
        "Insertion rVAF should be positive"
    );
}

// ============================================================
// Scenario C: Deletion detection
// ============================================================

#[test]
fn test_pipeline_deletion_detection() {
    // Reference: ACGTTAG (k=3), k-mers: ACG, CGT, GTT, TTA, TAG
    // Deletion of T at position 3: ACGTAG
    //   Variant k-mer: GTA (creating path ACG->CGT->GTA->TAG)
    let k = 3;
    let ref_seq = "ACGTTAG";

    let mut db = MockDb::new(k as u8);
    db.set_sequence(ref_seq, 500);
    // Deletion variant
    db.set("GTA", 150);

    let result = run_pipeline(&db, ref_seq, k);

    assert!(
        result.paths.len() >= 2,
        "Deletion pipeline should find at least 2 paths, found {}",
        result.paths.len()
    );

    // Reference
    assert_eq!(result.classifications[0].0, VariantType::Reference);
    assert_eq!(result.paths[0].to_sequence(), ref_seq);

    // Find the deletion path
    let alt_idx = result
        .classifications
        .iter()
        .position(|(vt, _, _)| *vt == VariantType::Deletion)
        .expect("Should find a Deletion variant");

    let alt_seq = result.paths[alt_idx].to_sequence();
    assert_eq!(
        alt_seq, "ACGTAG",
        "Deletion path should have 1 fewer base"
    );

    // Alt sequence should be shorter than reference (deletion removes bases)
    assert!(
        alt_seq.len() < ref_seq.len(),
        "Deletion alt sequence ({}) should be shorter than reference ({})",
        alt_seq.len(),
        ref_seq.len()
    );
}

// ============================================================
// Scenario D: No variant (reference only)
// ============================================================

#[test]
fn test_pipeline_reference_only() {
    // Reference: ACGTACG (k=4), no variant k-mers in database
    // Should produce exactly 1 path (the reference) with rVAF = 1.0
    let k = 4;
    let ref_seq = "ACGTACG";

    let mut db = MockDb::new(k as u8);
    db.set_sequence(ref_seq, 1000);
    // No variant k-mers -- all non-reference extensions have zero count

    let result = run_pipeline(&db, ref_seq, k);

    assert_eq!(
        result.paths.len(),
        1,
        "Reference-only should have exactly 1 path"
    );
    assert_eq!(result.classifications[0].0, VariantType::Reference);
    assert_eq!(result.paths[0].to_sequence(), ref_seq);

    // rVAF should be 1.0
    assert!(
        (result.rvafs[0] - 1.0).abs() < 1e-9,
        "Single reference path rVAF should be 1.0, got {}",
        result.rvafs[0]
    );

    // min_coverage should equal the set count
    assert_eq!(
        result.min_coverages[0], 1000,
        "Reference min coverage should match the set count"
    );
}

// ============================================================
// Scenario E: Multiple variants from same target
// ============================================================

#[test]
fn test_pipeline_multiple_variants() {
    // Reference: ACGTAG (k=3), k-mers: ACG, CGT, GTA, TAG
    // Two variant paths:
    //   Alt 1: Insertion of T -> ACGTTAG (k-mers: GTT, TTA)
    //   Alt 2: Insertion of C -> ACGCTAG (k-mers: CGC, GCT, CTA)
    let k = 3;
    let ref_seq = "ACGTAG";

    let mut db = MockDb::new(k as u8);
    db.set_sequence(ref_seq, 1000);
    // Alt 1 k-mers at 200x (20% VAF)
    db.set("GTT", 200);
    db.set("TTA", 200);
    // Alt 2 k-mers at 100x (10% VAF)
    db.set("CGC", 100);
    db.set("GCT", 100);
    db.set("CTA", 100);

    let result = run_pipeline(&db, ref_seq, k);

    // Should find reference + at least 2 alternatives
    assert!(
        result.paths.len() >= 3,
        "Multiple variant pipeline should find >= 3 paths, found {}",
        result.paths.len()
    );

    // Reference first
    assert_eq!(result.classifications[0].0, VariantType::Reference);

    // Count non-reference variants
    let variant_count = result
        .classifications
        .iter()
        .filter(|(vt, _, _)| *vt != VariantType::Reference)
        .count();
    assert!(
        variant_count >= 2,
        "Should have at least 2 non-reference variants, found {}",
        variant_count
    );

    // All rVAFs should be positive and sum to ~1.0
    for (i, &rvaf) in result.rvafs.iter().enumerate() {
        assert!(
            rvaf >= 0.0,
            "rVAF[{}] should be >= 0, got {}",
            i,
            rvaf
        );
    }
    let rvaf_sum: f64 = result.rvafs.iter().sum();
    assert!(
        (rvaf_sum - 1.0).abs() < 0.05,
        "rVAFs should sum to ~1.0, got {}",
        rvaf_sum
    );

    // Reference should have the highest rVAF (most reads support reference)
    assert!(
        result.rvafs[0] > result.rvafs[1],
        "Reference rVAF ({}) should be > first alt rVAF ({})",
        result.rvafs[0],
        result.rvafs[1]
    );
}

// ============================================================
// Scenario F: High-frequency variant (50% VAF)
// ============================================================

#[test]
fn test_pipeline_high_frequency_variant() {
    // SNV at 50% VAF: equal support for ref and alt
    let k = 4;
    let ref_seq = "ACGTACG";

    let mut db = MockDb::new(k as u8);
    // Reference k-mers at their "true" contribution (500 from ref allele)
    for kmer in &decompose(ref_seq, k) {
        db.set(kmer, 500);
    }
    // Variant k-mers at equal contribution (500 from alt allele)
    db.set("CGTT", 500);
    db.set("GTTA", 500);
    db.set("TTAC", 500);
    // Shared k-mers (ACGT, TACG) need to be at double count since they
    // are shared between ref and alt paths
    db.set("ACGT", 1000);
    db.set("TACG", 1000);

    let result = run_pipeline(&db, ref_seq, k);

    assert!(
        result.paths.len() >= 2,
        "High-freq variant should have >= 2 paths"
    );

    // Find the variant
    let alt_idx = result
        .classifications
        .iter()
        .position(|(vt, _, _)| *vt != VariantType::Reference);

    if let Some(alt_idx) = alt_idx {
        // At 50% VAF, both ref and alt should have roughly equal rVAFs
        let ref_rvaf = result.rvafs[0];
        let alt_rvaf = result.rvafs[alt_idx];
        assert!(
            (ref_rvaf - alt_rvaf).abs() < 0.2,
            "At ~50% VAF, ref ({:.3}) and alt ({:.3}) rVAFs should be roughly equal",
            ref_rvaf,
            alt_rvaf
        );
    }
}

// ============================================================
// Scenario G: Very low frequency variant (near detection limit)
// ============================================================

#[test]
fn test_pipeline_low_frequency_variant() {
    // SNV at ~1% VAF: barely above default thresholds
    let k = 4;
    let ref_seq = "ACGTACG";

    let mut db = MockDb::new(k as u8);
    db.set_sequence(ref_seq, 10000);
    // Variant at 1% frequency
    db.set("CGTT", 100);
    db.set("GTTA", 100);
    db.set("TTAC", 100);

    let config = WalkerConfig {
        count: 2,
        ratio: 0.005, // Lower ratio to detect low-frequency variants
        ..WalkerConfig::default()
    };

    let ref_kmers = decompose(ref_seq, k);
    let walk_result = walk(&db, &ref_kmers, &config);

    // The variant k-mers should still be discovered at 1% if ratio threshold is low enough
    // 100 / (10000 + 100) = ~0.99% which is above ratio=0.5%
    assert!(
        walk_result.nodes.contains_key("CGTT"),
        "1% variant should be discovered with ratio=0.005"
    );

    // Run full pipeline
    let graph = build_graph(&walk_result, &ref_kmers);
    let paths = find_alternative_paths(&graph);

    if paths.len() >= 2 {
        let quant = quantify(&paths, &db);
        // Variant rVAF should be very low
        for i in 1..quant.rvafs.len() {
            assert!(
                quant.rvafs[i] < 0.1,
                "Low-freq variant rVAF should be < 0.1, got {}",
                quant.rvafs[i]
            );
        }
    }
}

// ============================================================
// Scenario H: Pipeline with FASTA target loading
// ============================================================

#[test]
fn test_pipeline_with_fasta_target() {
    // Test that we can create a FASTA file, load it as a target,
    // decompose into k-mers, and run the pipeline
    use kmerdet::sequence::target::{RefSeq, Target};

    let target = Target {
        name: "test_target_snv".to_string(),
        sequence: "ACGTACG".to_string(),
        source: std::path::PathBuf::from("/tmp/test.fa"),
    };

    let k = 4u8;
    let ref_seq = RefSeq::from_target(target, k).unwrap();

    assert_eq!(ref_seq.kmers.len(), 4); // ACGT, CGTA, GTAC, TACG
    assert_eq!(ref_seq.kmers[0], "ACGT");
    assert_eq!(ref_seq.kmers[3], "TACG");
    assert_eq!(ref_seq.target.name, "test_target_snv");

    // Run pipeline with the decomposed k-mers
    let mut db = MockDb::new(k);
    for kmer in &ref_seq.kmers {
        db.set(kmer, 500);
    }
    db.set("CGTT", 50);
    db.set("GTTA", 50);
    db.set("TTAC", 50);

    let config = WalkerConfig::default();
    let walk_result = walk(&db, &ref_seq.kmers, &config);
    let graph = build_graph(&walk_result, &ref_seq.kmers);
    let paths = find_alternative_paths(&graph);

    // Verify we got the variant path
    assert!(
        paths.len() >= 2,
        "Pipeline with FASTA target should find variant"
    );
    let ref_path_seq = paths[0].to_sequence();
    assert_eq!(ref_path_seq, "ACGTACG");
}

// ============================================================
// Scenario I: Pipeline produces consistent min_coverage
// ============================================================

#[test]
fn test_pipeline_min_coverage_tracking() {
    // Verify that min_coverage correctly reflects the bottleneck k-mer
    let k = 3;
    let ref_seq = "ACGTAG";

    let mut db = MockDb::new(k as u8);
    // Set reference k-mers with varying coverage
    db.set("ACG", 1000);
    db.set("CGT", 500); // This is the bottleneck for reference path
    db.set("GTA", 800);
    db.set("TAG", 900);

    // Variant insertion k-mers
    db.set("GTT", 50); // This is the bottleneck for variant path
    db.set("TTA", 200);

    let result = run_pipeline(&db, ref_seq, k);

    assert!(result.paths.len() >= 2);

    // Reference path min coverage should be the bottleneck (CGT = 500)
    assert_eq!(
        result.min_coverages[0], 500,
        "Reference min_coverage should be 500 (CGT bottleneck), got {}",
        result.min_coverages[0]
    );

    // Variant path min coverage should reflect the lowest k-mer in that path
    if result.paths.len() >= 2 {
        // The alt path goes through ACG, CGT, GTT, TTA, TAG
        // Min is GTT = 50
        assert!(
            result.min_coverages[1] <= 200,
            "Variant min_coverage should reflect the variant k-mer bottleneck"
        );
    }
}

// ============================================================
// Scenario J: Empty input handling
// ============================================================

#[test]
fn test_pipeline_short_sequence() {
    // Sequence shorter than k should fail gracefully in RefSeq
    use kmerdet::sequence::target::{RefSeq, Target};

    let target = Target {
        name: "too_short".to_string(),
        sequence: "AC".to_string(),
        source: std::path::PathBuf::from("/tmp/test.fa"),
    };

    let result = RefSeq::from_target(target, 4);
    assert!(
        result.is_err(),
        "Sequence shorter than k should return an error"
    );
}

// ============================================================
// Scenario K: Walk discovers variant with default thresholds
// ============================================================

#[test]
fn test_pipeline_default_thresholds() {
    // Verify that default km thresholds (count=2, ratio=0.05) correctly filter
    let k = 4;
    let ref_seq = "ACGTACG";
    let ref_kmers = decompose(ref_seq, k);

    let mut db = MockDb::new(k as u8);
    db.set_sequence(ref_seq, 100);

    // A variant k-mer above thresholds (should be found)
    db.set("CGTT", 10); // 10/110 = ~9% > 5% ratio, and count=10 > 2

    // A variant k-mer below thresholds (should NOT be found)
    db.set("CGTG", 1); // count=1 < 2

    let config = WalkerConfig::default(); // count=2, ratio=0.05
    let result = walk(&db, &ref_kmers, &config);

    assert!(
        result.nodes.contains_key("CGTT"),
        "Variant above threshold should be discovered"
    );
    assert!(
        !result.nodes.contains_key("CGTG"),
        "Variant below threshold should NOT be discovered"
    );
}

// ============================================================
// Scenario L: Quantification rVAF accuracy
// ============================================================

#[test]
fn test_pipeline_rvaf_accuracy() {
    // Test that NNLS produces accurate rVAFs with clean data.
    // Two paths with well-separated k-mers and known counts.
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

    // Set up 75% ref, 25% alt
    let mut db = MockDb::new(4);
    db.set("AAAA", 400); // shared: 300 ref + 100 alt = 400
    db.set("AAAC", 300); // ref-only
    db.set("AACG", 300); // ref-only
    db.set("AAAT", 100); // alt-only
    db.set("AATG", 100); // alt-only

    let quant = quantify(&[ref_path, alt_path], &db);

    // rVAFs should be approximately 0.75 and 0.25
    assert!(
        (quant.rvafs[0] - 0.75).abs() < 0.05,
        "ref rVAF should be ~0.75, got {}",
        quant.rvafs[0]
    );
    assert!(
        (quant.rvafs[1] - 0.25).abs() < 0.05,
        "alt rVAF should be ~0.25, got {}",
        quant.rvafs[1]
    );
}
