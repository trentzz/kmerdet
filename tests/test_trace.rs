/// Tests for the detection trace module.
///
/// Verifies that trace JSON output is valid, that write_trace_target and
/// write_trace_summary produce correct files, and that trace structures
/// serialize properly.

use kmerdet::trace::*;

/// Build a sample DetectionTrace for testing.
fn sample_trace(target_name: &str, outcome: &str) -> DetectionTrace {
    DetectionTrace {
        target: target_name.to_string(),
        walking: WalkingTrace {
            reference_kmers: 5,
            reference_kmer_counts: vec![100, 95, 110, 80, 90],
            total_nodes: 8,
            alt_nodes: 3,
            branching_points: 1,
            limits_hit: None,
        },
        graph: GraphTrace {
            total_nodes: 10,
            reference_nodes: 5,
            alt_nodes: 3,
            virtual_nodes: 2,
            total_edges: 12,
            reference_edges: 6,
            alt_edges: 6,
        },
        pathfinding: PathfindingTrace {
            paths_found: 2,
            reference_path_length: Some(5),
            reference_path_sequence: Some("ACGTACGT".to_string()),
            alternative_paths: vec![AltPathTrace {
                index: 0,
                length: 5,
                sequence: "ACGTTCGT".to_string(),
                shared_prefix_length: 4,
                shared_suffix_length: 3,
            }],
        },
        classifications: vec![ClassificationTrace {
            path_index: 1,
            variant_type: "Substitution".to_string(),
            variant_name: "4:A/T:4".to_string(),
            ref_allele: "A".to_string(),
            alt_allele: "T".to_string(),
            start: 4,
            end: 4,
        }],
        quantification: Some(QuantificationTrace {
            num_paths: 2,
            num_unique_kmers: 7,
            coefficients: vec![100.0, 25.0],
            rvafs: vec![0.8, 0.2],
            min_coverages: vec![80, 20],
        }),
        outcome: outcome.to_string(),
    }
}

#[test]
fn test_trace_serializes_to_valid_json() {
    let trace = sample_trace("test_target", "variant_detected");

    let json = serde_json::to_string_pretty(&trace).unwrap();
    let parsed: serde_json::Value = serde_json::from_str(&json).unwrap();

    assert_eq!(parsed["target"], "test_target");
    assert_eq!(parsed["outcome"], "variant_detected");
    assert_eq!(parsed["walking"]["reference_kmers"], 5);
    assert_eq!(parsed["walking"]["total_nodes"], 8);
    assert_eq!(parsed["walking"]["alt_nodes"], 3);
    assert_eq!(parsed["graph"]["total_nodes"], 10);
    assert_eq!(parsed["graph"]["virtual_nodes"], 2);
    assert_eq!(parsed["pathfinding"]["paths_found"], 2);
    assert_eq!(parsed["classifications"].as_array().unwrap().len(), 1);
    assert_eq!(parsed["classifications"][0]["variant_type"], "Substitution");
    assert_eq!(parsed["classifications"][0]["variant_name"], "4:A/T:4");
    assert!(parsed["quantification"].is_object());
    assert_eq!(parsed["quantification"]["num_paths"], 2);
}

#[test]
fn test_write_trace_target_produces_valid_json_file() {
    let trace = sample_trace("test_target_json", "variant_detected");
    let dir = tempfile::tempdir().unwrap();

    write_trace_target(&trace, dir.path()).unwrap();

    let target_file = dir.path().join("targets/test_target_json.json");
    assert!(target_file.exists(), "target JSON file should exist");

    let content = std::fs::read_to_string(&target_file).unwrap();
    let parsed: serde_json::Value = serde_json::from_str(&content).unwrap();

    assert_eq!(parsed["target"], "test_target_json");
    assert_eq!(parsed["outcome"], "variant_detected");
    assert_eq!(parsed["walking"]["reference_kmers"], 5);
}

#[test]
fn test_write_trace_summary_produces_valid_json_file() {
    let traces = vec![
        sample_trace("target_A", "variant_detected"),
        sample_trace("target_B", "reference_only"),
        sample_trace("target_C", "no_paths"),
    ];
    let dir = tempfile::tempdir().unwrap();

    write_trace_summary(&traces, dir.path()).unwrap();

    let summary_file = dir.path().join("summary.json");
    assert!(summary_file.exists(), "summary.json should exist");

    let content = std::fs::read_to_string(&summary_file).unwrap();
    let parsed: serde_json::Value = serde_json::from_str(&content).unwrap();

    assert_eq!(parsed["total_targets"], 3);
    assert_eq!(parsed["variants_detected"], 1);
    assert_eq!(parsed["reference_only"], 1);
    assert_eq!(parsed["no_paths"], 1);

    let targets_arr = parsed["targets"].as_array().unwrap();
    assert_eq!(targets_arr.len(), 3);
    assert_eq!(targets_arr[0]["target"], "target_A");
    assert_eq!(targets_arr[0]["outcome"], "variant_detected");
}

#[test]
fn test_write_trace_target_sanitizes_special_characters() {
    let trace = sample_trace("chr1:1000-2000/SNV", "variant_detected");
    let dir = tempfile::tempdir().unwrap();

    write_trace_target(&trace, dir.path()).unwrap();

    // Colons, slashes should be replaced with underscores.
    let expected_file = dir.path().join("targets/chr1_1000-2000_SNV.json");
    assert!(expected_file.exists(), "sanitized filename should exist");
}

#[test]
fn test_trace_without_quantification() {
    let mut trace = sample_trace("ref_only", "reference_only");
    trace.quantification = None;
    trace.classifications = vec![];

    let json = serde_json::to_string_pretty(&trace).unwrap();
    let parsed: serde_json::Value = serde_json::from_str(&json).unwrap();

    assert_eq!(parsed["outcome"], "reference_only");
    assert!(parsed["quantification"].is_null());
    assert!(parsed["classifications"].as_array().unwrap().is_empty());
}

#[test]
fn test_write_multiple_target_traces() {
    let traces = vec![
        sample_trace("target_1", "variant_detected"),
        sample_trace("target_2", "reference_only"),
    ];
    let dir = tempfile::tempdir().unwrap();

    for trace in &traces {
        write_trace_target(trace, dir.path()).unwrap();
    }
    write_trace_summary(&traces, dir.path()).unwrap();

    // All files should exist.
    assert!(dir.path().join("targets/target_1.json").exists());
    assert!(dir.path().join("targets/target_2.json").exists());
    assert!(dir.path().join("summary.json").exists());

    // Each target file should be valid JSON with the correct target name.
    for name in &["target_1", "target_2"] {
        let content =
            std::fs::read_to_string(dir.path().join(format!("targets/{}.json", name))).unwrap();
        let parsed: serde_json::Value = serde_json::from_str(&content).unwrap();
        assert_eq!(parsed["target"].as_str().unwrap(), *name);
    }
}

#[test]
fn test_trace_empty_summary() {
    let dir = tempfile::tempdir().unwrap();
    write_trace_summary(&[], dir.path()).unwrap();

    let content = std::fs::read_to_string(dir.path().join("summary.json")).unwrap();
    let parsed: serde_json::Value = serde_json::from_str(&content).unwrap();

    assert_eq!(parsed["total_targets"], 0);
    assert_eq!(parsed["variants_detected"], 0);
}
