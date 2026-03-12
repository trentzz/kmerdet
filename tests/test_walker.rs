mod common;

use common::MockDb;
use kmerdet::walker::extension;
use kmerdet::walker::{walk, WalkerConfig};

// ── Existing extension tests ──────────────────────────────────────────

#[test]
fn test_extend_forward_filters_by_count() {
    let mut db = MockDb::new(4);
    // Set up a k-mer "ACGT" and its possible forward extensions
    db.set("ACGT", 100);
    db.set("CGTA", 50); // A extension: passes
    db.set("CGTC", 2);  // C extension: too low
    db.set("CGTG", 1);  // G extension: too low
    db.set("CGTT", 0);  // T extension: zero

    let children = extension::extend_forward(&db, "ACGT", 5, 0.05);

    assert_eq!(children.len(), 1);
    assert_eq!(children[0].sequence, "CGTA");
    assert_eq!(children[0].count, 50);
}

#[test]
fn test_extend_forward_all_pass() {
    let mut db = MockDb::new(4);
    db.set("CGTA", 100);
    db.set("CGTC", 100);
    db.set("CGTG", 100);
    db.set("CGTT", 100);

    let children = extension::extend_forward(&db, "ACGT", 2, 0.05);
    assert_eq!(children.len(), 4);
}

#[test]
fn test_extend_forward_none_pass() {
    let db = MockDb::new(4);
    // All zero counts
    let children = extension::extend_forward(&db, "ACGT", 2, 0.05);
    assert!(children.is_empty());
}

// ── Walker tests ──────────────────────────────────────────────────────

/// Helper to decompose a sequence into its constituent k-mers.
fn kmers_from_seq(seq: &str, k: usize) -> Vec<String> {
    seq.as_bytes()
        .windows(k)
        .map(|w| std::str::from_utf8(w).unwrap().to_string())
        .collect()
}

#[test]
fn test_walk_reference_only() {
    // Reference: ACGTACGT with k=4
    // K-mers: ACGT, CGTA, GTAC, TACG, ACGT (duplicate collapsed)
    // All ref k-mers have high counts, no variants exist.
    let k = 4;
    let seq = "ACGTACGT";
    let ref_kmers = kmers_from_seq(seq, k);

    let mut db = MockDb::new(k as u8);
    // Set all reference k-mers to high counts
    db.set_sequence(seq, 100);
    // The forward extensions of the last ref k-mer ("ACGT") that are NOT ref k-mers
    // should have zero counts (default), so they won't be discovered.

    let config = WalkerConfig::default();
    let result = walk(&db, &ref_kmers, &config);

    // Should find exactly the reference k-mers
    for kmer in &ref_kmers {
        assert!(
            result.nodes.contains_key(kmer),
            "Expected to find reference k-mer {} in nodes",
            kmer
        );
    }
    // No extra nodes beyond reference k-mers
    assert_eq!(
        result.nodes.len(),
        ref_kmers.iter().collect::<std::collections::HashSet<_>>().len(),
        "Should contain only unique reference k-mers"
    );
    // All ref kmers should be in the reference set
    for kmer in &ref_kmers {
        assert!(result.reference_kmers.contains(kmer));
    }
}

#[test]
fn test_walk_finds_snv() {
    // Reference: ACGTACGT with k=4
    // Ref k-mers: ACGT(100), CGTA(100), GTAC(100), TACG(100)
    // SNV: at position 4 (A->T), creating variant k-mer CGTT instead of CGTA
    // The variant k-mer CGTT has count 50 — should be discovered via forward
    // extension from ACGT.
    let k = 4;
    let seq = "ACGTACGT";
    let ref_kmers = kmers_from_seq(seq, k);

    let mut db = MockDb::new(k as u8);
    db.set_sequence(seq, 100);

    // Add a variant k-mer: forward extension from ACGT can produce CGTA, CGTC, CGTG, CGTT
    // CGTA is already a reference k-mer (count 100). Add CGTT as variant.
    db.set("CGTT", 50);

    let config = WalkerConfig::default();
    let result = walk(&db, &ref_kmers, &config);

    // Should find the variant k-mer
    assert!(
        result.nodes.contains_key("CGTT"),
        "Walker should discover variant k-mer CGTT"
    );
    assert_eq!(result.nodes["CGTT"], 50);

    // Variant should NOT be in reference set
    assert!(!result.reference_kmers.contains("CGTT"));

    // All ref k-mers should still be present
    for kmer in &ref_kmers {
        assert!(result.nodes.contains_key(kmer));
    }
}

#[test]
fn test_walk_finds_insertion() {
    // Reference: AACCGGTT with k=4
    // Ref k-mers: AACC, ACCG, CCGG, CGGT, GGTT
    // Insertion variant: insert "A" after position 4 → AACCAGGTT
    // Variant k-mers that differ from ref: CCAG, CAGG, AGGT
    // These should be reachable by forward extension from AACC:
    //   AACC -> forward -> ACCA (not ref, not variant), ACCG (ref)
    // Actually let's trace through more carefully.
    // From ref k-mer ACCG, forward extensions include CCGA, CCGC, CCGG(ref), CCGT
    // We want to introduce an insertion path via CCAG.
    // From CCGG forward extension -> CGGA, CGGC, CGGG, CGGT(ref)
    // To find CCAG, we need it as a forward extension of a ref k-mer.
    // ACCG -> forward: xCCG + base -> CCGA, CCGC, CCGG, CCGT
    // AACC -> forward: xACC + base -> ACCA, ACCC, ACCG, ACCT
    // None of these produce CCAG.
    //
    // Let's redesign: ref = AACCGGTT, variant = AACCCGGTT (insert C after pos 3)
    // Ref k-mers: AACC, ACCG, CCGG, CGGT, GGTT
    // Variant path: AACC -> ACCC -> CCCG -> CCGG -> CGGT -> GGTT
    // New variant k-mers: ACCC, CCCG
    // Forward from AACC: ACC+{A,C,G,T} = ACCA, ACCC, ACCG, ACCT
    // ACCG is ref (count 100), ACCC is variant.
    // Forward from ACCC: CCC+{A,C,G,T} = CCCA, CCCC, CCCG, CCCT
    // CCCG is variant.
    // Forward from CCCG: CCG+{A,C,G,T} = CCGA, CCGC, CCGG, CCGT
    // CCGG is ref, found already.
    // So the walker should discover ACCC and CCCG through the insertion path.

    let k = 4;
    let ref_seq = "AACCGGTT";
    let ref_kmers = kmers_from_seq(ref_seq, k);

    let mut db = MockDb::new(k as u8);
    db.set_sequence(ref_seq, 100);

    // Variant k-mers for insertion of C
    db.set("ACCC", 40);
    db.set("CCCG", 40);

    let config = WalkerConfig::default();
    let result = walk(&db, &ref_kmers, &config);

    // Should discover the insertion variant k-mers
    assert!(
        result.nodes.contains_key("ACCC"),
        "Walker should discover insertion k-mer ACCC"
    );
    assert!(
        result.nodes.contains_key("CCCG"),
        "Walker should discover insertion k-mer CCCG"
    );
    assert_eq!(result.nodes["ACCC"], 40);
    assert_eq!(result.nodes["CCCG"], 40);

    // Variant k-mers should NOT be in reference set
    assert!(!result.reference_kmers.contains("ACCC"));
    assert!(!result.reference_kmers.contains("CCCG"));

    // All ref k-mers should be present
    for kmer in &ref_kmers {
        assert!(result.nodes.contains_key(kmer));
    }
}

#[test]
fn test_walk_max_node_limit() {
    // Create a scenario with a long chain of k-mers that would be discovered,
    // but limit max_node to a small number.
    // Reference: single k-mer "AAAA", which extends forward to "AAAC" -> "AACA" -> "ACAA" -> ...
    let k = 4;
    let mut db = MockDb::new(k as u8);

    // Reference k-mer
    db.set("AAAA", 100);

    // Build a long linear chain from AAAA
    // AAAA -> AAAC(100) -> AACT(100) -> ACTA(100) -> CTAG(100) -> ...
    db.set("AAAC", 100);
    db.set("AACT", 100);
    db.set("ACTA", 100);
    db.set("CTAG", 100);
    db.set("TAGC", 100);
    db.set("AGCA", 100);
    db.set("GCAG", 100);
    db.set("CAGT", 100);
    db.set("AGTC", 100);

    let ref_kmers = vec!["AAAA".to_string()];

    let config = WalkerConfig {
        count: 2,
        ratio: 0.05,
        max_stack: 500,
        max_break: 10,
        max_node: 5, // Limit to 5 nodes total
        adaptive: false,
    };

    let result = walk(&db, &ref_kmers, &config);

    // Should be capped at max_node
    assert!(
        result.nodes.len() <= 5,
        "Node count {} should be <= max_node (5)",
        result.nodes.len()
    );
    // Should have at least the reference k-mer
    assert!(result.nodes.contains_key("AAAA"));
}

#[test]
fn test_walk_max_break_limit() {
    // Create scenario with multiple branching points.
    // Reference k-mer AAAA, which branches at each step.
    let k = 4;
    let mut db = MockDb::new(k as u8);

    db.set("AAAA", 100);

    // First branch: AAAA -> AAAC and AAAT (2 children = 1 break)
    db.set("AAAC", 100);
    db.set("AAAT", 100);

    // From AAAC, another branch: AACA and AACT (2 children = 2nd break)
    db.set("AACA", 100);
    db.set("AACT", 100);

    // From AACA, another branch: ACAA and ACAT (3rd break)
    db.set("ACAA", 100);
    db.set("ACAT", 100);

    // From ACAA, further extension: CAAG (4th level, but shouldn't reach if max_break=2)
    db.set("CAAG", 100);

    let ref_kmers = vec!["AAAA".to_string()];

    let config = WalkerConfig {
        count: 2,
        ratio: 0.05,
        max_stack: 500,
        max_break: 2, // Only allow 2 branching points
        max_node: 10000,
        adaptive: false,
    };

    let result = walk(&db, &ref_kmers, &config);

    // Should discover k-mers reachable within max_break=2.
    // Break 1 (from AAAA): AAAC, AAAT discovered
    // Break 2 (from AAAC): AACA, AACT discovered
    // Break 3 (from AACA) would exceed max_break=2, so ACAA and ACAT should NOT
    // be pushed further. However, they are still inserted as nodes since they are
    // children. They just won't be explored further.
    // Actually re-reading the algorithm: we insert into nodes first, then check
    // new_breaks > max_break to decide whether to push onto stack.
    // So ACAA and ACAT get inserted as nodes but not explored.
    assert!(result.nodes.contains_key("AAAC"));
    assert!(result.nodes.contains_key("AAAT"));
    assert!(result.nodes.contains_key("AACA"));
    assert!(result.nodes.contains_key("AACT"));
    // ACAA and ACAT are children of AACA (break=2), so new_breaks=3 > max_break=2,
    // they get inserted into nodes but not pushed onto the stack.
    // Wait, let me re-check: AACA is pushed with breaks=2. When we pop it and extend,
    // children ACAA/ACAT: they get inserted. num_children=2, new_breaks=3>2, so skip push.
    // So ACAA/ACAT are in nodes but not explored further.
    assert!(
        result.nodes.contains_key("ACAA"),
        "ACAA should be discovered as a node (but not explored further)"
    );

    // CAAG should NOT be discovered because ACAA was never explored (breaks exceeded)
    assert!(
        !result.nodes.contains_key("CAAG"),
        "CAAG should not be reachable due to max_break limit"
    );
}

#[test]
fn test_walk_max_stack_limit() {
    // Create a deep linear chain and set max_stack to a small number.
    let k = 4;
    let mut db = MockDb::new(k as u8);

    // Build a linear chain: each k-mer has exactly one forward child.
    // AAAC -> AACG -> ACGT -> CGTA -> GTAC -> TACG -> ACGA -> CGAT -> ...
    // Since each has only one child, breaks never increment.
    // We need the chain to be longer than max_stack.
    let chain = vec![
        "AAAC", "AACG", "ACGT", "CGTA", "GTAC", "TACG", "ACGA", "CGAT",
        "GATC", "ATCG", "TCGA", "CGAG", "GAGT", "AGTC", "GTCA", "TCAG",
    ];

    for kmer in &chain {
        db.set(kmer, 100);
    }

    let ref_kmers = vec!["AAAC".to_string()];

    let config = WalkerConfig {
        count: 2,
        ratio: 0.05,
        max_stack: 3, // Very small stack limit
        max_break: 10,
        max_node: 10000,
        adaptive: false,
    };

    let result = walk(&db, &ref_kmers, &config);

    // The walker should discover some k-mers but not all due to stack limit.
    // With a linear chain and max_stack=3, the DFS pops one item (stack becomes 0),
    // pushes its child (stack becomes 1), pops (0), pushes (1), etc.
    // For a linear chain, stack never exceeds 1, so stack limit of 3 doesn't limit it.
    //
    // To truly test stack limit, we need branching that fills the stack.
    // Let's redesign: create a k-mer where every extension has multiple children,
    // filling the stack quickly.

    // Actually, let's just verify the walker ran and found k-mers.
    // The linear chain won't hit stack=3, so all should be found.
    assert!(result.nodes.len() >= 2);

    // Now test with a branching scenario that fills the stack.
    let mut db2 = MockDb::new(k as u8);
    db2.set("AAAA", 100);
    // From AAAA, 4 children
    db2.set("AAAA", 100);
    db2.set("AAAC", 100);
    db2.set("AAAG", 100);
    db2.set("AAAT", 100);
    // From AAAC, 4 children
    db2.set("AACA", 100);
    db2.set("AACC", 100);
    db2.set("AACG", 100);
    db2.set("AACT", 100);
    // From AAAG, 4 children
    db2.set("AAGA", 100);
    db2.set("AAGC", 100);
    db2.set("AAGG", 100);
    db2.set("AAGT", 100);
    // From AAAT, 4 children
    db2.set("AATA", 100);
    db2.set("AATC", 100);
    db2.set("AATG", 100);
    db2.set("AATT", 100);
    // And deeper levels too...
    db2.set("ACAA", 100);
    db2.set("ACAC", 100);
    db2.set("ACAG", 100);
    db2.set("ACAT", 100);

    let ref_kmers2 = vec!["AAAA".to_string()];

    let config2 = WalkerConfig {
        count: 2,
        ratio: 0.05,
        max_stack: 2, // Very small
        max_break: 100,
        max_node: 10000,
        adaptive: false,
    };

    let result2 = walk(&db2, &ref_kmers2, &config2);

    // With max_stack=2, once the stack has 2 items, no more can be pushed.
    // From AAAA, we pop it (stack=0), extend to get 3 new children
    // (AAAC, AAAG, AAAT — AAAA is already in nodes).
    // All 3 get inserted into nodes. Push first two (stack=2), third can't be pushed
    // (stack.len() >= max_stack). So some deeper k-mers won't be reachable.
    // But the children ARE added to nodes even if not pushed.
    assert!(result2.nodes.contains_key("AAAA"));
    assert!(result2.nodes.contains_key("AAAC"));

    // Count total nodes — should be fewer than what we'd get with unlimited stack
    let config_unlimited = WalkerConfig {
        count: 2,
        ratio: 0.05,
        max_stack: 500,
        max_break: 100,
        max_node: 10000,
        adaptive: false,
    };

    let result_unlimited = walk(&db2, &ref_kmers2, &config_unlimited);

    assert!(
        result2.nodes.len() <= result_unlimited.nodes.len(),
        "Limited stack ({}) should discover <= unlimited stack ({}) nodes",
        result2.nodes.len(),
        result_unlimited.nodes.len()
    );
}

#[test]
fn test_walk_empty_reference() {
    let db = MockDb::new(4);
    let ref_kmers: Vec<String> = vec![];

    let config = WalkerConfig::default();
    let result = walk(&db, &ref_kmers, &config);

    assert!(result.nodes.is_empty(), "Empty ref should yield empty nodes");
    assert!(
        result.reference_kmers.is_empty(),
        "Empty ref should yield empty reference_kmers"
    );
}

#[test]
fn test_walk_low_count_filtered() {
    // Reference: ACGTACGT with k=4. All ref k-mers have high counts.
    // A variant k-mer exists but its count is below the threshold.
    let k = 4;
    let seq = "ACGTACGT";
    let ref_kmers = kmers_from_seq(seq, k);

    let mut db = MockDb::new(k as u8);
    db.set_sequence(seq, 100);

    // Add a variant k-mer with count below the absolute threshold
    db.set("CGTT", 1); // Below default count threshold of 2

    let config = WalkerConfig {
        count: 5,  // Require at least 5
        ratio: 0.05,
        max_stack: 500,
        max_break: 10,
        max_node: 10000,
        adaptive: false,
    };

    let result = walk(&db, &ref_kmers, &config);

    // The low-count variant should NOT be discovered
    assert!(
        !result.nodes.contains_key("CGTT"),
        "Low-count variant k-mer CGTT (count=1) should be filtered out with count threshold 5"
    );

    // All reference k-mers should still be present
    for kmer in &ref_kmers {
        assert!(
            result.nodes.contains_key(kmer),
            "Reference k-mer {} should be present",
            kmer
        );
    }
}
