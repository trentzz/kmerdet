use kmerdet::kmer::Kmer;

#[test]
fn test_encode_decode_roundtrip() {
    let sequences = ["ACGT", "AAAA", "TTTT", "GCTA", "ACGTACGTACGTACGT"];
    for seq in &sequences {
        let kmer = Kmer::from_str(seq).expect("valid sequence");
        assert_eq!(kmer.to_string(), *seq);
    }
}

#[test]
fn test_reverse_complement() {
    let kmer = Kmer::from_str("ACGT").unwrap();
    let rc = kmer.reverse_complement();
    assert_eq!(rc.to_string(), "ACGT"); // ACGT is its own reverse complement
}

#[test]
fn test_reverse_complement_asymmetric() {
    let kmer = Kmer::from_str("AAAC").unwrap();
    let rc = kmer.reverse_complement();
    assert_eq!(rc.to_string(), "GTTT");
}

#[test]
fn test_canonical_form() {
    let kmer = Kmer::from_str("TTTT").unwrap();
    let canonical = kmer.canonical();
    // canonical = min(TTTT, AAAA) = AAAA
    assert_eq!(canonical.to_string(), "AAAA");
}

#[test]
fn test_extend_right() {
    let kmer = Kmer::from_str("ACGT").unwrap();
    // Extend right with A: ACGT -> CGTA (drop first, append A)
    let extended = kmer.extend_right(0b00); // A = 0b00
    assert_eq!(extended.to_string(), "CGTA");
}

#[test]
fn test_extend_left() {
    let kmer = Kmer::from_str("ACGT").unwrap();
    // Extend left with T: ACGT -> TACG (prepend T, drop last)
    let extended = kmer.extend_left(0b11); // T = 0b11
    assert_eq!(extended.to_string(), "TACG");
}

#[test]
fn test_invalid_base_returns_none() {
    assert!(Kmer::from_str("ACGN").is_none());
    assert!(Kmer::from_str("").is_none());
}
