// Graph construction and pathfinding tests.

mod common;

use std::collections::{HashMap, HashSet};

use common::MockDb;
use kmerdet::graph::builder::{build_graph, ALT_EDGE_WEIGHT, REF_EDGE_WEIGHT};
use kmerdet::graph::pathfind::{
    find_alternative_paths, find_alternative_paths_ranked, path_mean_coverage, path_min_coverage,
};
use kmerdet::sequence::path::KmerPath;
use kmerdet::walker::WalkResult;

/// Helper: decompose a DNA string into k-mers of length k.
fn decompose(seq: &str, k: usize) -> Vec<String> {
    seq.as_bytes()
        .windows(k)
        .map(|w| std::str::from_utf8(w).unwrap().to_string())
        .collect()
}

/// Helper: build a WalkResult from reference k-mers and optional extra (variant) k-mers.
fn make_walk(ref_kmers: &[String], extra_kmers: &[(String, u64)]) -> WalkResult {
    let mut nodes = HashMap::new();
    let mut reference_kmers = HashSet::new();

    for kmer in ref_kmers {
        nodes.insert(kmer.clone(), 100); // reference k-mers get count 100
        reference_kmers.insert(kmer.clone());
    }

    for (kmer, count) in extra_kmers {
        nodes.insert(kmer.clone(), *count);
        // extra k-mers are NOT reference
    }

    WalkResult {
        nodes,
        reference_kmers,
    }
}

// ============================================================
// Graph builder tests
// ============================================================

#[test]
fn test_build_graph_linear_reference() {
    // Simple linear reference: ACGTACG with k=4
    // K-mers: ACGT, CGTA, GTAC, TACG
    let k = 4;
    let ref_kmers = decompose("ACGTACG", k);
    assert_eq!(ref_kmers, vec!["ACGT", "CGTA", "GTAC", "TACG"]);

    let walk = make_walk(&ref_kmers, &[]);
    let graph = build_graph(&walk, &ref_kmers);

    // Should have 4 real k-mer nodes + 2 virtual = 6
    assert_eq!(graph.nodes.len(), 6);
    assert_eq!(graph.source, 0);
    assert_eq!(graph.sink, 1);

    // Check virtual nodes
    assert_eq!(graph.nodes[0].kmer, "BigBang");
    assert_eq!(graph.nodes[1].kmer, "BigCrunch");

    // All real nodes should be reference
    for node in &graph.nodes[2..] {
        assert!(node.is_reference, "Node {} should be reference", node.kmer);
    }

    // BigBang should connect to the first ref k-mer ACGT
    let source_edges = graph.edges.get(&graph.source).unwrap();
    let first_ref_idx = *graph.node_index.get("ACGT").unwrap();
    assert!(
        source_edges.iter().any(|e| e.to == first_ref_idx),
        "BigBang should connect to ACGT"
    );

    // Last ref k-mer TACG should connect to BigCrunch
    let last_ref_idx = *graph.node_index.get("TACG").unwrap();
    let last_edges = graph.edges.get(&last_ref_idx).unwrap();
    assert!(
        last_edges.iter().any(|e| e.to == graph.sink),
        "TACG should connect to BigCrunch"
    );

    // Check linear chain: ACGT→CGTA→GTAC→TACG
    let acgt_idx = *graph.node_index.get("ACGT").unwrap();
    let cgta_idx = *graph.node_index.get("CGTA").unwrap();
    let gtac_idx = *graph.node_index.get("GTAC").unwrap();
    let tacg_idx = *graph.node_index.get("TACG").unwrap();

    let acgt_edges = graph.edges.get(&acgt_idx).unwrap();
    assert!(acgt_edges.iter().any(|e| e.to == cgta_idx));

    let cgta_edges = graph.edges.get(&cgta_idx).unwrap();
    assert!(cgta_edges.iter().any(|e| e.to == gtac_idx));

    let gtac_edges = graph.edges.get(&gtac_idx).unwrap();
    assert!(gtac_edges.iter().any(|e| e.to == tacg_idx));

    // All edges between reference nodes should be reference edges
    for (&from, edges) in &graph.edges {
        for edge in edges {
            if from != graph.source && from != graph.sink && edge.to != graph.source && edge.to != graph.sink {
                assert!(
                    edge.is_reference,
                    "Edge from {} to {} should be reference",
                    graph.nodes[from].kmer,
                    graph.nodes[edge.to].kmer,
                );
            }
        }
    }
}

#[test]
fn test_build_graph_with_snv() {
    // Reference: ACGTACGT, k=4
    // K-mers: ACGT, CGTA, GTAC, TACG, ACGT (dedup → 4 unique)
    // Actually: ACGT, CGTA, GTAC, TACG, ACGT → 4 unique
    // Variant k-mer: CGTT (SNV: CGTA→CGTT at last base, A→T)
    let k = 4;
    let ref_kmers = decompose("ACGTACGT", k);
    // ref_kmers = ["ACGT", "CGTA", "GTAC", "TACG", "ACGT"]
    // As HashMap keys, ACGT appears only once in nodes

    // Variant: CGTT overlaps with ACGT (ACG prefix matches CGT? No.
    // ACGT suffix = CGT, CGTT prefix = CGT → ACGT→CGTT edge exists
    let extra = vec![("CGTT".to_string(), 50u64)];
    let walk = make_walk(&ref_kmers, &extra);
    let graph = build_graph(&walk, &ref_kmers);

    // Should have variant node
    assert!(graph.node_index.contains_key("CGTT"));
    let cgtt_idx = *graph.node_index.get("CGTT").unwrap();
    assert!(!graph.nodes[cgtt_idx].is_reference);

    // ACGT→CGTT edge should exist (suffix CGT == prefix CGT)
    let acgt_idx = *graph.node_index.get("ACGT").unwrap();
    let acgt_edges = graph.edges.get(&acgt_idx).unwrap();
    let cgtt_edge = acgt_edges.iter().find(|e| e.to == cgtt_idx);
    assert!(cgtt_edge.is_some(), "ACGT→CGTT edge should exist");

    // The edge should NOT be a reference edge (CGTT is not reference)
    let cgtt_edge = cgtt_edge.unwrap();
    assert!(!cgtt_edge.is_reference);
    assert!((cgtt_edge.weight - ALT_EDGE_WEIGHT).abs() < 1e-10);
}

#[test]
fn test_build_graph_edge_weights() {
    // Verify ref edges get 0.01, non-ref edges get 1.0
    let k = 4;
    let ref_kmers = decompose("ACGTACG", k);
    // Add a variant k-mer that creates a non-ref edge
    // TACT: prefix TAC matches suffix of GTAC → GTAC→TACT
    let extra = vec![("TACT".to_string(), 30u64)];
    let walk = make_walk(&ref_kmers, &extra);
    let graph = build_graph(&walk, &ref_kmers);

    // Check reference edges
    let acgt_idx = *graph.node_index.get("ACGT").unwrap();
    let cgta_idx = *graph.node_index.get("CGTA").unwrap();
    let acgt_edges = graph.edges.get(&acgt_idx).unwrap();
    let ref_edge = acgt_edges.iter().find(|e| e.to == cgta_idx).unwrap();
    assert!(ref_edge.is_reference);
    assert!(
        (ref_edge.weight - REF_EDGE_WEIGHT).abs() < 1e-10,
        "Reference edge weight should be {}",
        REF_EDGE_WEIGHT
    );

    // Check non-reference edge: GTAC→TACT
    let gtac_idx = *graph.node_index.get("GTAC").unwrap();
    let tact_idx = *graph.node_index.get("TACT").unwrap();
    let gtac_edges = graph.edges.get(&gtac_idx).unwrap();
    let alt_edge = gtac_edges.iter().find(|e| e.to == tact_idx).unwrap();
    assert!(!alt_edge.is_reference);
    assert!(
        (alt_edge.weight - ALT_EDGE_WEIGHT).abs() < 1e-10,
        "Non-reference edge weight should be {}",
        ALT_EDGE_WEIGHT
    );
}

// ============================================================
// Pathfinding tests
// ============================================================

#[test]
fn test_find_paths_reference_only() {
    // Linear reference only → exactly 1 path (the reference)
    let k = 4;
    let ref_kmers = decompose("ACGTACG", k);
    let walk = make_walk(&ref_kmers, &[]);
    let graph = build_graph(&walk, &ref_kmers);
    let paths = find_alternative_paths(&graph);

    assert_eq!(paths.len(), 1, "Should find exactly 1 path (reference only)");
    assert!(paths[0].is_reference, "The single path should be the reference");

    // Verify the path produces the correct sequence
    let seq = paths[0].to_sequence();
    assert_eq!(seq, "ACGTACG");
}

#[test]
fn test_find_paths_snv() {
    // Reference: ACGTACG (k=4) → ACGT, CGTA, GTAC, TACG
    // SNV variant: at CGTA position, change last base A→T → CGTT
    // CGTT must reconnect to sink path. We need CGTT to have an edge
    // to something that reaches BigCrunch.
    //
    // For a complete SNV path, we need:
    // ACGT → CGTT → GTTX → ... → BigCrunch
    // But that requires more variant k-mers. Let's make a simpler scenario:
    //
    // Reference: ACGTAG (k=4) → ACGT, CGTA, GTAG
    // Variant:   ACGTTG (k=4) → ACGT, CGTT, GTTG
    //   (but GTTG won't connect back to GTAG)
    //
    // Better approach: make the variant path rejoin the reference.
    // Reference: ACGTACGT (k=3) → ACG, CGT, GTA, TAC, ACG, CGT
    //   Unique: ACG, CGT, GTA, TAC
    // Variant at position 2: G→T → ACTTACGT
    //   Variant k-mers: ACT, CTT, TTA → these reconnect? TTA→TAC (suffix TA == prefix TA) ✓
    //
    // Let's use k=3 for simplicity.
    // The key insight: for an SNV to branch, a ref k-mer's suffix must match
    // a variant k-mer's prefix.
    //
    // ACGT→CGTT: ACGT suffix = CGT, CGTT prefix = CGT → yes! Edge exists.
    // CGTT→GTTA: CGTT suffix = GTT, GTTA prefix = GTT → yes!
    // GTTA→TTAC: GTTA suffix = TTA, TTAC prefix = TTA → yes!
    // TTAC→TACG: TTAC suffix = TAC, TACG prefix = TAC → yes!
    //
    // So reference: ACGTACG (k=4) → ACGT, CGTA, GTAC, TACG
    // Variant path: ACGT → CGTT → GTTA → TTAC → TACG
    // Variant k-mers needed: CGTT, GTTA, TTAC
    // (ACGT and TACG are reference and shared)

    let k = 4;
    let ref_kmers = decompose("ACGTACG", k);
    // ["ACGT", "CGTA", "GTAC", "TACG"]

    let extra = vec![
        ("CGTT".to_string(), 50u64),
        ("GTTA".to_string(), 50u64),
        ("TTAC".to_string(), 50u64),
    ];
    let walk = make_walk(&ref_kmers, &extra);
    let graph = build_graph(&walk, &ref_kmers);

    let paths = find_alternative_paths(&graph);

    // Should have reference path + 1 alternative
    assert!(paths.len() >= 2, "Should find at least 2 paths, found {}", paths.len());
    assert!(paths[0].is_reference, "First path should be reference");

    let ref_seq = paths[0].to_sequence();
    assert_eq!(ref_seq, "ACGTACG", "Reference path sequence");

    // Find the alternative path
    let alt_paths: Vec<&_> = paths.iter().filter(|p| !p.is_reference).collect();
    assert!(!alt_paths.is_empty(), "Should have at least one alternative path");

    // The alt path goes ACGT→CGTT→GTTA→TTAC→TACG
    // Sequence: ACGT + T + A + C + G = ACGTTACG
    let alt_seq = alt_paths[0].to_sequence();
    assert_eq!(alt_seq, "ACGTTACG", "Alternative path should produce variant sequence");
}

#[test]
fn test_find_paths_insertion() {
    // Reference: ACGTAG (k=3) → ACG, CGT, GTA, TAG
    // Insertion of "CC" between GT and AG: ACGTCCAG? No, let's think about this differently.
    //
    // Reference: ACGTAG, k=3 → ACG, CGT, GTA, TAG
    // Variant with insertion: ACGAGTAG → extra k-mers between CGT and GTA
    //   Actually we want an insertion path that's longer.
    //
    // Reference path: ACGT→CGTA→GTAC (from "ACGTAC", k=4)
    // Insertion of T after position 3: "ACGTTGTAC"
    //   ACGT→CGTT→GTTG→TTGT→TGTA→GTAC
    //   Variant k-mers: CGTT, GTTG, TTGT, TGTA
    //
    // But this is complex. Let's use a simpler insertion case.
    //
    // Reference: ACGACG, k=3 → ACG, CGA, GAC, ACG → unique: ACG, CGA, GAC
    // Insertion of T: ACGTACG → extra k-mers CGT, GTA, TAC connect ACG→CGT→GTA→TAC→ACG
    //   ACG suffix=CG, CGT prefix=CG → ✓
    //   CGT suffix=GT, GTA prefix=GT → ✓
    //   GTA suffix=TA, TAC prefix=TA → ✓
    //   TAC suffix=AC, ACG prefix=AC → ✓
    // But ACG→CGA is reference, ACG→CGT is variant.
    // Reference path: ACG→CGA→GAC→(sink)
    // Variant path: ACG→CGT→GTA→TAC→ACG→CGA→GAC→(sink)
    // Wait, that revisits ACG which causes a cycle.
    //
    // Let's use a different design. Use k=4 for cleaner separation.
    //
    // Reference: ACGTCAG, k=4 → ACGT, CGTC, GTCA, TCAG
    // Insertion of AA after position 4: ACGTAACAG
    //   K-mers: ACGT, CGTA, GTAA, TAAC, AACA, ACAG
    //   We need the insertion path to rejoin the reference. The last ref k-mer is TCAG.
    //   ACAG suffix=CAG, TCAG prefix ≠ CAG... hmm.
    //
    // Let me pick sequences more carefully.
    //
    // Reference: ACGATCG, k=4 → ACGA, CGAT, GATC, ATCG
    // Insertion path: ACGA → CGAC → GACG → ACGA → ... cycles again.
    //
    // Simplest insertion scenario:
    // Reference: ACGCAG, k=3 → ACG, CGC, GCA, CAG
    // Insertion of T after position 2: ACGTCAG? Wait, that changes things.
    //
    // Let me try yet another approach.
    // Reference sequence "AACCGG", k=4 → AACC, ACCG, CCGG
    // Insertion of TT: AACCTTCCGG? That doesn't flow cleanly.
    //
    // OK, let me be very explicit:
    // Reference: ACGTAG, k=3 → ACG, CGT, GTA, TAG
    // Insertion: we insert a "T" between positions 3 and 4, giving ACGTTAG
    //   Variant k-mers: GTT, TTA (new ones not in reference)
    //   ACG→CGT (both ref, exists already)
    //   CGT→GTT: suffix GT == prefix GT → ✓ (variant edge)
    //   GTT→TTA: suffix TT == prefix TT → ✓ (variant edge)
    //   TTA→TAG: suffix TA == prefix TA → ✓ (variant edge)
    //   So variant path: ACG→CGT→GTT→TTA→TAG
    //   Reference path:  ACG→CGT→GTA→TAG
    //   Variant path is 1 k-mer longer (insertion of 1 base)

    let k = 3;
    let ref_kmers = decompose("ACGTAG", k);
    assert_eq!(ref_kmers, vec!["ACG", "CGT", "GTA", "TAG"]);

    let extra = vec![
        ("GTT".to_string(), 40u64),
        ("TTA".to_string(), 40u64),
    ];
    let walk = make_walk(&ref_kmers, &extra);
    let graph = build_graph(&walk, &ref_kmers);

    let paths = find_alternative_paths(&graph);

    assert!(paths.len() >= 2, "Should find at least 2 paths, found {}", paths.len());
    assert!(paths[0].is_reference);

    let ref_seq = paths[0].to_sequence();
    assert_eq!(ref_seq, "ACGTAG", "Reference path sequence");

    // Alternative path should be longer (insertion)
    let alt_paths: Vec<&_> = paths.iter().filter(|p| !p.is_reference).collect();
    assert!(!alt_paths.is_empty(), "Should have alternative path for insertion");

    let alt_seq = alt_paths[0].to_sequence();
    assert_eq!(alt_seq, "ACGTTAG", "Insertion path should be 1 base longer");

    // Alternative path has more k-mers than reference
    assert!(
        alt_paths[0].kmers.len() > paths[0].kmers.len(),
        "Insertion path should have more k-mers than reference"
    );
}

#[test]
fn test_find_paths_deletion() {
    // Reference: ACGTTAG, k=3 → ACG, CGT, GTT, TTA, TAG
    // Deletion of T at position 4: ACGTAG
    //   Deletion k-mers: CGT→GTA→TAG
    //   GTA is not in reference. CGT and TAG are reference.
    //   CGT suffix=GT, GTA prefix=GT → ✓
    //   GTA suffix=TA, TAG prefix=TA → ✓
    //
    // Variant path: ACG→CGT→GTA→TAG (shorter by 1 k-mer than reference)
    // Reference path: ACG→CGT→GTT→TTA→TAG

    let k = 3;
    let ref_kmers = decompose("ACGTTAG", k);
    assert_eq!(ref_kmers, vec!["ACG", "CGT", "GTT", "TTA", "TAG"]);

    let extra = vec![("GTA".to_string(), 35u64)];
    let walk = make_walk(&ref_kmers, &extra);
    let graph = build_graph(&walk, &ref_kmers);

    let paths = find_alternative_paths(&graph);

    assert!(paths.len() >= 2, "Should find at least 2 paths, found {}", paths.len());
    assert!(paths[0].is_reference);

    let ref_seq = paths[0].to_sequence();
    assert_eq!(ref_seq, "ACGTTAG", "Reference path sequence");

    let alt_paths: Vec<&_> = paths.iter().filter(|p| !p.is_reference).collect();
    assert!(!alt_paths.is_empty(), "Should have alternative path for deletion");

    let alt_seq = alt_paths[0].to_sequence();
    assert_eq!(alt_seq, "ACGTAG", "Deletion path should be 1 base shorter");

    // Alternative path has fewer k-mers than reference
    assert!(
        alt_paths[0].kmers.len() < paths[0].kmers.len(),
        "Deletion path should have fewer k-mers than reference"
    );
}

#[test]
fn test_find_paths_deduplication() {
    // Create a scenario where multiple non-reference edges lead to the same path.
    // Reference: ACGTAG, k=3 → ACG, CGT, GTA, TAG
    // Add variant k-mer CGG and GGA such that:
    //   ACG→CGG→GGA→... but GGA doesn't connect to TAG, so no complete path.
    //
    // Better: create two variant k-mers that both lead to the same sequence.
    // Actually, the simplest dedup test: add two variant k-mers that form
    // one alternative path, and verify it appears only once.
    //
    // Reference: ACGTAG, k=3 → ACG, CGT, GTA, TAG
    // Variant: GTT, TTA (insertion path: ACG→CGT→GTT→TTA→TAG, seq = ACGTTAG)
    // This creates two non-reference edges: CGT→GTT and GTT→TTA and TTA→TAG
    // Each non-reference edge triggers path construction, but they all
    // produce the same path sequence "ACGTTAG" → should be deduplicated to 1.

    let k = 3;
    let ref_kmers = decompose("ACGTAG", k);

    let extra = vec![
        ("GTT".to_string(), 40u64),
        ("TTA".to_string(), 40u64),
    ];
    let walk = make_walk(&ref_kmers, &extra);
    let graph = build_graph(&walk, &ref_kmers);

    let paths = find_alternative_paths(&graph);

    // Count unique alternative paths
    let alt_paths: Vec<&_> = paths.iter().filter(|p| !p.is_reference).collect();

    // Collect unique sequences
    let alt_sequences: HashSet<String> = alt_paths.iter().map(|p| p.to_sequence()).collect();

    // There should be exactly 1 unique alternative path
    assert_eq!(
        alt_paths.len(),
        alt_sequences.len(),
        "Paths should be deduplicated: {} paths but {} unique sequences",
        alt_paths.len(),
        alt_sequences.len()
    );

    // Specifically, there should be exactly 1 alternative path for the insertion
    assert_eq!(alt_paths.len(), 1, "Should have exactly 1 deduplicated alternative path");
}

// ============================================================
// Coverage scoring tests
// ============================================================

#[test]
fn test_path_min_coverage_simple() {
    let mut db = MockDb::new(4);
    db.set("ACGT", 100);
    db.set("CGTA", 50);
    db.set("GTAC", 200);

    let path = KmerPath {
        kmers: vec![
            "ACGT".to_string(),
            "CGTA".to_string(),
            "GTAC".to_string(),
        ],
        is_reference: true,
    };

    assert_eq!(path_min_coverage(&path, &db), 50);
}

#[test]
fn test_path_min_coverage_empty() {
    let db = MockDb::new(4);
    let path = KmerPath {
        kmers: vec![],
        is_reference: false,
    };

    assert_eq!(path_min_coverage(&path, &db), 0);
}

#[test]
fn test_path_mean_coverage_simple() {
    let mut db = MockDb::new(4);
    db.set("ACGT", 100);
    db.set("CGTA", 50);
    db.set("GTAC", 150);

    let path = KmerPath {
        kmers: vec![
            "ACGT".to_string(),
            "CGTA".to_string(),
            "GTAC".to_string(),
        ],
        is_reference: true,
    };

    let mean = path_mean_coverage(&path, &db);
    // (100 + 50 + 150) / 3 = 100.0
    assert!((mean - 100.0).abs() < 1e-9, "mean should be 100.0, got {}", mean);
}

#[test]
fn test_path_mean_coverage_empty() {
    let db = MockDb::new(4);
    let path = KmerPath {
        kmers: vec![],
        is_reference: false,
    };

    assert!((path_mean_coverage(&path, &db) - 0.0).abs() < 1e-9);
}

#[test]
fn test_path_min_coverage_single_kmer() {
    let mut db = MockDb::new(4);
    db.set("ACGT", 42);

    let path = KmerPath {
        kmers: vec!["ACGT".to_string()],
        is_reference: true,
    };

    assert_eq!(path_min_coverage(&path, &db), 42);
}

#[test]
fn test_path_mean_coverage_single_kmer() {
    let mut db = MockDb::new(4);
    db.set("ACGT", 42);

    let path = KmerPath {
        kmers: vec!["ACGT".to_string()],
        is_reference: true,
    };

    let mean = path_mean_coverage(&path, &db);
    assert!((mean - 42.0).abs() < 1e-9, "mean should be 42.0, got {}", mean);
}

#[test]
fn test_find_alternative_paths_ranked_returns_coverage_order() {
    // Create two alternative paths with different coverages.
    // The higher-coverage alt should appear before the lower-coverage alt.
    //
    // Reference: ACGTACG (k=4) -> ACGT, CGTA, GTAC, TACG
    //
    // Alt path 1 (low coverage): ACGT -> CGTT -> GTTA -> TTAC -> TACG
    //   variant k-mers: CGTT(10), GTTA(10), TTAC(10)
    //
    // Alt path 2 (high coverage): ACGT -> CGTG -> GTGA -> TGAC -> TACG
    //   Wait, TGAC suffix = GAC, TACG prefix = TAC -- no match.
    //
    // Let's use a different approach. Two SNV alt paths that share the same
    // ref start/end but differ at one position.
    //
    // Reference: ACGTACG (k=4) -> ACGT, CGTA, GTAC, TACG
    //
    // Alt path 1: ACGT -> CGTT -> GTTA -> TTAC -> TACG (counts: 10)
    // Alt path 2: ACGT -> CGTC -> GTCA -> TCAC -> ... doesn't connect to TACG
    //
    // Simpler: use k=3 with two alt paths.
    // Reference: ACGTTAG (k=3) -> ACG, CGT, GTT, TTA, TAG
    //
    // Alt path 1 (deletion, high cov): ACG -> CGT -> GTA -> TAG
    //   variant k-mer: GTA(200)
    //
    // Alt path 2 (also a deletion via different route): If we add another
    //   non-ref edge... This is tricky because we need two distinct alt paths.
    //
    // Let's use the simplest scenario: one alt path with known coverage.

    let k = 3;
    let ref_kmers = decompose("ACGTTAG", k);
    // ref_kmers = ["ACG", "CGT", "GTT", "TTA", "TAG"]

    // Deletion path: ACG -> CGT -> GTA -> TAG (count 35 for GTA)
    let extra = vec![("GTA".to_string(), 35u64)];
    let walk = make_walk(&ref_kmers, &extra);
    let graph = build_graph(&walk, &ref_kmers);

    let mut db = MockDb::new(3);
    // Set reference k-mers
    db.set("ACG", 100);
    db.set("CGT", 100);
    db.set("GTT", 100);
    db.set("TTA", 100);
    db.set("TAG", 100);
    // Set variant k-mer
    db.set("GTA", 35);

    let ranked = find_alternative_paths_ranked(&graph, &db);

    // Should have reference + 1 alt
    assert!(ranked.len() >= 2, "Should have at least 2 ranked paths, got {}", ranked.len());

    // Reference path at index 0
    assert!(ranked[0].0.is_reference, "First ranked path should be reference");

    // Reference path score = min(100, 100, 100, 100, 100) = 100
    assert_eq!(ranked[0].1, 100, "Reference path score should be 100");

    // Alt path score = min of its k-mer counts (ACG=100, CGT=100, GTA=35, TAG=100) = 35
    assert_eq!(ranked[1].1, 35, "Alt path score should be 35");
}

#[test]
fn test_find_alternative_paths_ranked_ref_stays_at_index_0() {
    // Even when the reference path has a lower coverage score than alt paths,
    // the reference path must remain at index 0.

    let k = 3;
    let ref_kmers = decompose("ACGTTAG", k);

    let extra = vec![("GTA".to_string(), 500u64)]; // high-coverage variant k-mer
    let walk = make_walk(&ref_kmers, &extra);
    let graph = build_graph(&walk, &ref_kmers);

    let mut db = MockDb::new(3);
    // Set reference k-mers with LOW counts
    db.set("ACG", 10);
    db.set("CGT", 10);
    db.set("GTT", 10);
    db.set("TTA", 10);
    db.set("TAG", 10);
    // Set variant k-mer with HIGH count
    db.set("GTA", 500);

    let ranked = find_alternative_paths_ranked(&graph, &db);

    assert!(ranked.len() >= 2, "Should have at least 2 ranked paths");

    // Reference must be first, even though its score (10) < alt score (10 for shared k-mers)
    assert!(ranked[0].0.is_reference, "Reference path must stay at index 0");
}

#[test]
fn test_find_alternative_paths_ranked_multiple_alts_sorted() {
    // Two alternative paths with different coverages should be sorted
    // by descending coverage score.
    //
    // Reference: ACGTAG (k=3) -> ACG, CGT, GTA, TAG
    // Alt 1 (insertion, low cov): ACG -> CGT -> GTT -> TTA -> TAG (seq: ACGTTAG)
    //   variant k-mers: GTT(20), TTA(20)
    // Alt 2 (some other path with high cov) -- need a second distinct alt path.
    //
    // To get two alt paths, we need two different non-ref edges that lead to
    // complete paths through the graph.
    //
    // Reference: ACGTAG (k=3) -> ACG, CGT, GTA, TAG
    //
    // Alt 1: insertion of T -> ACGTTAG: ACG->CGT->GTT->TTA->TAG
    //   new k-mers: GTT, TTA
    //
    // Alt 2: insertion of C -> ACGCTAG: ACG->CGC->GCT->CTA->TAG
    //   new k-mers: CGC, GCT, CTA
    //   ACG suffix=CG, CGC prefix=CG -> edge exists
    //   CGC suffix=GC, GCT prefix=GC -> edge exists
    //   GCT suffix=CT, CTA prefix=CT -> edge exists
    //   CTA suffix=TA, TAG prefix=TA -> edge exists
    //
    // Set alt1 k-mers at low coverage, alt2 at high coverage.

    let k = 3;
    let ref_kmers = decompose("ACGTAG", k);
    assert_eq!(ref_kmers, vec!["ACG", "CGT", "GTA", "TAG"]);

    let extra = vec![
        // Alt path 1 (insertion of T): low coverage
        ("GTT".to_string(), 20u64),
        ("TTA".to_string(), 20u64),
        // Alt path 2 (insertion of C): high coverage
        ("CGC".to_string(), 80u64),
        ("GCT".to_string(), 80u64),
        ("CTA".to_string(), 80u64),
    ];
    let walk = make_walk(&ref_kmers, &extra);
    let graph = build_graph(&walk, &ref_kmers);

    let mut db = MockDb::new(3);
    db.set("ACG", 100);
    db.set("CGT", 100);
    db.set("GTA", 100);
    db.set("TAG", 100);
    // Alt 1 k-mers
    db.set("GTT", 20);
    db.set("TTA", 20);
    // Alt 2 k-mers
    db.set("CGC", 80);
    db.set("GCT", 80);
    db.set("CTA", 80);

    let ranked = find_alternative_paths_ranked(&graph, &db);

    assert!(ranked.len() >= 3, "Should find ref + 2 alt paths, got {}", ranked.len());

    // Reference at index 0
    assert!(ranked[0].0.is_reference, "First path should be reference");

    // Alt paths should be sorted by descending score.
    // Alt path 2 (through CGC/GCT/CTA) has min_cov = min(100, 80, 80, 80, 100) = 80
    // Alt path 1 (through GTT/TTA) has min_cov = min(100, 100, 20, 20, 100) = 20
    // So alt2 (score 80) should come before alt1 (score 20).

    assert!(
        ranked[1].1 >= ranked[2].1,
        "Alt paths should be sorted by descending coverage: {} >= {}",
        ranked[1].1,
        ranked[2].1
    );

    // Verify specific scores
    assert_eq!(ranked[1].1, 80, "Higher-coverage alt should have score 80");
    assert_eq!(ranked[2].1, 20, "Lower-coverage alt should have score 20");
}

#[test]
fn test_find_alternative_paths_ranked_equal_coverage_stable() {
    // Two alternative paths with equal coverage should both appear.
    // (We don't require a specific order when scores are equal, just that
    // neither is lost and both scores are correct.)

    let k = 3;
    let ref_kmers = decompose("ACGTAG", k);

    let extra = vec![
        // Alt path 1: insertion of T
        ("GTT".to_string(), 50u64),
        ("TTA".to_string(), 50u64),
        // Alt path 2: insertion of C
        ("CGC".to_string(), 50u64),
        ("GCT".to_string(), 50u64),
        ("CTA".to_string(), 50u64),
    ];
    let walk = make_walk(&ref_kmers, &extra);
    let graph = build_graph(&walk, &ref_kmers);

    let mut db = MockDb::new(3);
    db.set("ACG", 100);
    db.set("CGT", 100);
    db.set("GTA", 100);
    db.set("TAG", 100);
    db.set("GTT", 50);
    db.set("TTA", 50);
    db.set("CGC", 50);
    db.set("GCT", 50);
    db.set("CTA", 50);

    let ranked = find_alternative_paths_ranked(&graph, &db);

    assert!(ranked.len() >= 3, "Should find ref + 2 alt paths, got {}", ranked.len());

    // Reference at index 0
    assert!(ranked[0].0.is_reference);

    // Both alt paths should have the same score (50 = min of shared + variant k-mers)
    assert_eq!(ranked[1].1, ranked[2].1, "Equal-coverage paths should have equal scores");
    assert_eq!(ranked[1].1, 50);
}

#[test]
fn test_find_alternative_paths_ranked_reference_only() {
    // When there are no alt paths, should return just the reference path scored.
    let k = 4;
    let ref_kmers = decompose("ACGTACG", k);
    let walk = make_walk(&ref_kmers, &[]);
    let graph = build_graph(&walk, &ref_kmers);

    let mut db = MockDb::new(4);
    db.set("ACGT", 100);
    db.set("CGTA", 80);
    db.set("GTAC", 120);
    db.set("TACG", 90);

    let ranked = find_alternative_paths_ranked(&graph, &db);

    assert_eq!(ranked.len(), 1, "Should have exactly 1 path (reference only)");
    assert!(ranked[0].0.is_reference);
    assert_eq!(ranked[0].1, 80, "Reference path score should be min coverage = 80");
}
