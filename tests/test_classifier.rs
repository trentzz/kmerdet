// Variant classification integration tests.

use kmerdet::sequence::path::KmerPath;
use kmerdet::variant::classifier::classify;
use kmerdet::variant::VariantType;

/// Helper: build a KmerPath that produces the given DNA sequence for the given k.
///
/// For sequence "ACGT" with k=3, the k-mers are ["ACG", "CGT"].
/// to_sequence() reassembles: "ACG" + "T" = "ACGT".
fn make_path(seq: &str, k: usize, is_reference: bool) -> KmerPath {
    assert!(
        seq.len() >= k,
        "sequence length {} must be >= k={}",
        seq.len(),
        k
    );
    let kmers: Vec<String> = (0..=seq.len() - k)
        .map(|i| seq[i..i + k].to_string())
        .collect();
    let path = KmerPath {
        kmers,
        is_reference,
    };
    // Sanity check: the reconstructed sequence should match what we intended
    assert_eq!(
        path.to_sequence(),
        seq,
        "KmerPath reconstruction failed for '{}'",
        seq
    );
    path
}

#[test]
fn test_classify_reference() {
    // Identical ref and alt paths should produce Reference type
    let ref_path = make_path("ACGTACGT", 4, true);
    let alt_path = make_path("ACGTACGT", 4, false);
    let result = classify(&ref_path, &alt_path, 4);
    assert_eq!(result.variant_type, VariantType::Reference);
    assert!(result.ref_allele.is_empty());
    assert!(result.alt_allele.is_empty());
}

#[test]
fn test_classify_snv() {
    // Single base difference at position 4: A→T
    // ref: ACGTACGT
    // alt: ACGTTCGT
    //          ^
    let ref_path = make_path("ACGTACGT", 4, true);
    let alt_path = make_path("ACGTTCGT", 4, false);
    let result = classify(&ref_path, &alt_path, 4);
    assert_eq!(result.variant_type, VariantType::Substitution);
    assert_eq!(result.ref_allele, "A");
    assert_eq!(result.alt_allele, "T");
    assert_eq!(result.start, 4);
    assert_eq!(result.end, 4);
}

#[test]
fn test_classify_snv_first_base() {
    // SNV at the very first position
    let ref_path = make_path("ACGTACGT", 4, true);
    let alt_path = make_path("TCGTACGT", 4, false);
    let result = classify(&ref_path, &alt_path, 4);
    assert_eq!(result.variant_type, VariantType::Substitution);
    assert_eq!(result.ref_allele, "A");
    assert_eq!(result.alt_allele, "T");
    assert_eq!(result.start, 0);
}

#[test]
fn test_classify_snv_last_base() {
    // SNV at the very last position
    let ref_path = make_path("ACGTACGT", 4, true);
    let alt_path = make_path("ACGTACGA", 4, false);
    let result = classify(&ref_path, &alt_path, 4);
    assert_eq!(result.variant_type, VariantType::Substitution);
    assert_eq!(result.ref_allele, "T");
    assert_eq!(result.alt_allele, "A");
    assert_eq!(result.start, 7);
}

#[test]
fn test_classify_mnv() {
    // Multiple adjacent base changes: ACG→TCA at positions 2-4
    // ref: TTACGTTT
    // alt: TTTCATTT
    //        ^^^
    let ref_path = make_path("TTACGTTT", 4, true);
    let alt_path = make_path("TTTCATTT", 4, false);
    let result = classify(&ref_path, &alt_path, 4);
    assert_eq!(result.variant_type, VariantType::Substitution);
    assert_eq!(result.ref_allele, "ACG");
    assert_eq!(result.alt_allele, "TCA");
    assert_eq!(result.start, 2);
    assert_eq!(result.end, 4);
}

#[test]
fn test_classify_insertion_single_base() {
    // One base inserted: A inserted between existing bases
    // ref: ACGTACGT  (8 bases)
    // alt: ACGTAACGT (9 bases, A inserted after pos 4)
    //           ^
    let ref_path = make_path("ACGTACGT", 4, true);
    let alt_path = make_path("ACGTAACGT", 4, false);
    let result = classify(&ref_path, &alt_path, 4);
    assert!(
        result.variant_type == VariantType::Insertion,
        "expected Insertion, got {:?} (ref='{}', alt='{}')",
        result.variant_type,
        result.ref_allele,
        result.alt_allele
    );
}

#[test]
fn test_classify_insertion_multi_base() {
    // Multiple bases inserted: GGG inserted in the middle
    // ref: ACGTACGT     (8 bases)
    // alt: ACGTGGGACGT  (11 bases)
    let ref_path = make_path("ACGTACGT", 4, true);
    let alt_path = make_path("ACGTGGGACGT", 4, false);
    let result = classify(&ref_path, &alt_path, 4);
    assert!(
        result.variant_type == VariantType::Insertion,
        "expected Insertion, got {:?} (ref='{}', alt='{}')",
        result.variant_type,
        result.ref_allele,
        result.alt_allele,
    );
    assert!(
        result.alt_allele.len() > result.ref_allele.len(),
        "alt allele should be longer than ref allele in insertion"
    );
}

#[test]
fn test_classify_deletion_single_base() {
    // One base deleted: T at position 3 deleted
    // ref: ACGTACGT (8 bases)
    // alt: ACGACGT  (7 bases)
    //        ^
    let ref_path = make_path("ACGTACGT", 4, true);
    let alt_path = make_path("ACGACGT", 4, false);
    let result = classify(&ref_path, &alt_path, 4);
    assert!(
        result.variant_type == VariantType::Deletion,
        "expected Deletion, got {:?} (ref='{}', alt='{}')",
        result.variant_type,
        result.ref_allele,
        result.alt_allele,
    );
}

#[test]
fn test_classify_deletion_multi_base() {
    // Multiple bases deleted: GT at positions 2-3 deleted
    // ref: ACGTACGT (8 bases)
    // alt: ACACGT   (6 bases)
    let ref_path = make_path("ACGTACGT", 4, true);
    let alt_path = make_path("ACACGT", 4, false);
    let result = classify(&ref_path, &alt_path, 4);
    assert!(
        result.variant_type == VariantType::Deletion,
        "expected Deletion, got {:?} (ref='{}', alt='{}')",
        result.variant_type,
        result.ref_allele,
        result.alt_allele,
    );
    assert!(
        result.ref_allele.len() > result.alt_allele.len(),
        "ref allele should be longer than alt allele in deletion"
    );
}

#[test]
fn test_classify_complex_indel() {
    // Complex indel: different bases AND different lengths
    // For a true indel, the sequences must differ in length and the changed
    // bases must not be a simple prefix/suffix relationship.
    // ref: ACGTAACGT (9) → alt: ACTCCGT (7): GTA→TC (3 ref bases → 2 different alt bases)
    let ref_path = make_path("ACGTAACGT", 4, true);
    let alt_path = make_path("ACTCCGT", 4, false);
    let result = classify(&ref_path, &alt_path, 4);
    assert_eq!(
        result.variant_type,
        VariantType::Indel,
        "expected Indel, got {:?} (ref='{}', alt='{}')",
        result.variant_type,
        result.ref_allele,
        result.alt_allele,
    );
}

#[test]
fn test_classify_itd() {
    // Internal tandem duplication: ACGT is duplicated
    // ref: ACGTACGT       (8 bases)
    // alt: ACGTACGTACGT   (12 bases, ACGT repeated)
    let ref_path = make_path("ACGTACGT", 4, true);
    let alt_path = make_path("ACGTACGTACGT", 4, false);
    let result = classify(&ref_path, &alt_path, 4);
    assert_eq!(
        result.variant_type,
        VariantType::Itd,
        "expected ITD, got {:?} (ref='{}', alt='{}')",
        result.variant_type,
        result.ref_allele,
        result.alt_allele,
    );
}

#[test]
fn test_classify_itd_internal() {
    // ITD in the middle of the sequence
    // ref: GGACGTCC         (8 bases)
    // alt: GGACGTACGTCC     (12 bases, ACGT duplicated internally)
    let ref_path = make_path("GGACGTCC", 4, true);
    let alt_path = make_path("GGACGTACGTCC", 4, false);
    let result = classify(&ref_path, &alt_path, 4);
    assert_eq!(
        result.variant_type,
        VariantType::Itd,
        "expected ITD for internal duplication, got {:?} (ref='{}', alt='{}')",
        result.variant_type,
        result.ref_allele,
        result.alt_allele,
    );
}

#[test]
fn test_classify_variant_name_format() {
    // Verify variant_name follows "start:ref/alt:end" format
    let ref_path = make_path("ACGTACGT", 4, true);
    let alt_path = make_path("ACGTTCGT", 4, false);
    let result = classify(&ref_path, &alt_path, 4);

    let name = &result.variant_name;

    // Should contain '/' separating ref and alt alleles
    assert!(
        name.contains('/'),
        "variant_name '{}' should contain '/' separator",
        name
    );

    // Parse the expected format: "start:ref_allele/alt_allele:end"
    let expected = format!(
        "{}:{}/{}:{}",
        result.start, result.ref_allele, result.alt_allele, result.end
    );
    assert_eq!(
        name, &expected,
        "variant_name should follow 'start:ref/alt:end' format"
    );
}

#[test]
fn test_classify_alleles_correct() {
    // SNV: verify alleles are exactly the changed bases
    let ref_path = make_path("ACGTACGT", 4, true);
    let alt_path = make_path("ACGTTCGT", 4, false);
    let result = classify(&ref_path, &alt_path, 4);
    assert_eq!(result.ref_allele, "A", "ref_allele should be exactly 'A'");
    assert_eq!(result.alt_allele, "T", "alt_allele should be exactly 'T'");
    assert_eq!(result.start, 4, "SNV should be at position 4");
    assert_eq!(result.end, 4, "SNV end should equal start for single base");

    // MNV: verify alleles span exactly the changed region
    let ref_path2 = make_path("GGAACTTGG", 4, true);
    let alt_path2 = make_path("GGTCGATGG", 4, false);
    let result2 = classify(&ref_path2, &alt_path2, 4);
    assert_eq!(result2.variant_type, VariantType::Substitution);
    assert_eq!(result2.ref_allele.len(), result2.alt_allele.len());
    // Alleles should not include the flanking G's
    assert!(
        !result2.ref_allele.starts_with('G'),
        "ref_allele should not include flanking bases"
    );
}

#[test]
fn test_classify_reference_returns_empty_alleles() {
    let ref_path = make_path("AAGGTTCC", 4, true);
    let alt_path = make_path("AAGGTTCC", 4, false);
    let result = classify(&ref_path, &alt_path, 4);
    assert_eq!(result.variant_type, VariantType::Reference);
    assert!(
        result.ref_allele.is_empty(),
        "Reference type should have empty ref_allele"
    );
    assert!(
        result.alt_allele.is_empty(),
        "Reference type should have empty alt_allele"
    );
    assert!(
        result.variant_name.is_empty(),
        "Reference type should have empty variant_name"
    );
}

#[test]
fn test_classify_with_different_k_values() {
    // The classification result should be the same regardless of k
    // (k only affects path structure, not the sequence comparison)
    let ref_seq = "ACGTACGT";
    let alt_seq = "ACGTTCGT"; // SNV at pos 4

    for k in 3..=6 {
        let ref_path = make_path(ref_seq, k, true);
        let alt_path = make_path(alt_seq, k, false);
        let result = classify(&ref_path, &alt_path, k as u8);
        assert_eq!(
            result.variant_type,
            VariantType::Substitution,
            "k={}: expected Substitution, got {:?}",
            k,
            result.variant_type
        );
        assert_eq!(result.ref_allele, "A", "k={}: ref_allele mismatch", k);
        assert_eq!(result.alt_allele, "T", "k={}: alt_allele mismatch", k);
    }
}

#[test]
fn test_classify_longer_sequences() {
    // Test with a longer, more realistic sequence where alt is actually longer
    // ref: ATCGATCGATCG     (12 bases)
    // alt: ATCGATCCCGATCG   (14 bases, CCC inserted after pos 6)
    let ref_seq = "ATCGATCGATCG";
    let alt_seq = "ATCGATCCCGATCG";
    let ref_path = make_path(ref_seq, 4, true);
    let alt_path = make_path(alt_seq, 4, false);
    let result = classify(&ref_path, &alt_path, 4);
    assert!(
        result.variant_type == VariantType::Insertion
            || result.variant_type == VariantType::Itd
            || result.variant_type == VariantType::Indel,
        "expected insertion-like variant, got {:?} (ref='{}', alt='{}')",
        result.variant_type,
        result.ref_allele,
        result.alt_allele,
    );
}
