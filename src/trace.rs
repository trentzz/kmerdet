/// Detection trace module — captures intermediate results from the variant detection pipeline.
///
/// Each target processed by the detection pipeline produces a `DetectionTrace`
/// recording the walking, graph construction, pathfinding, classification,
/// and quantification steps. Traces can be serialized to JSON for debugging
/// and audit purposes.

use std::collections::HashSet;
use std::path::Path;

use anyhow::Result;
use serde::Serialize;

use crate::graph::KmerGraph;
use crate::jellyfish::KmerDatabase;
use crate::sequence::path::KmerPath;
use crate::variant::classifier::Classification;
use crate::variant::quantifier::Quantification;
use crate::walker::{WalkResult, WalkerConfig};

/// Full trace of the detection pipeline for a single target.
#[derive(Debug, Clone, Serialize)]
pub struct DetectionTrace {
    /// Target name (from FASTA header).
    pub target: String,
    /// Trace of the k-mer walking step.
    pub walking: WalkingTrace,
    /// Trace of the graph construction step.
    pub graph: GraphTrace,
    /// Trace of the pathfinding step.
    pub pathfinding: PathfindingTrace,
    /// Trace of each variant classification.
    pub classifications: Vec<ClassificationTrace>,
    /// Trace of the quantification step (None if no paths to quantify).
    pub quantification: Option<QuantificationTrace>,
    /// Pipeline outcome: "variant_detected", "reference_only", or "no_paths".
    pub outcome: String,
}

/// Trace of the k-mer walking (DFS exploration) step.
#[derive(Debug, Clone, Serialize)]
pub struct WalkingTrace {
    /// Number of reference k-mers used as seeds.
    pub reference_kmers: usize,
    /// Database counts for each reference k-mer (in order).
    pub reference_kmer_counts: Vec<u64>,
    /// Total nodes discovered by the walk (reference + alt).
    pub total_nodes: usize,
    /// Number of alt (non-reference) nodes discovered.
    pub alt_nodes: usize,
    /// Estimated number of branching points (nodes with >1 extension).
    pub branching_points: usize,
    /// If a walk limit was hit, describes which one.
    pub limits_hit: Option<String>,
}

/// Trace of the graph construction step.
#[derive(Debug, Clone, Serialize)]
pub struct GraphTrace {
    /// Total number of nodes in the graph.
    pub total_nodes: usize,
    /// Number of reference k-mer nodes.
    pub reference_nodes: usize,
    /// Number of alt (non-reference) k-mer nodes.
    pub alt_nodes: usize,
    /// Number of virtual nodes (BigBang, BigCrunch).
    pub virtual_nodes: usize,
    /// Total number of edges in the graph.
    pub total_edges: usize,
    /// Number of reference edges (connecting reference nodes).
    pub reference_edges: usize,
    /// Number of alt edges (at least one endpoint is non-reference).
    pub alt_edges: usize,
}

/// Trace of the pathfinding step.
#[derive(Debug, Clone, Serialize)]
pub struct PathfindingTrace {
    /// Total number of paths found (including reference).
    pub paths_found: usize,
    /// Length (in k-mers) of the reference path, if found.
    pub reference_path_length: Option<usize>,
    /// DNA sequence of the reference path, if found.
    pub reference_path_sequence: Option<String>,
    /// Traces for each alternative (non-reference) path.
    pub alternative_paths: Vec<AltPathTrace>,
}

/// Trace of a single alternative path.
#[derive(Debug, Clone, Serialize)]
pub struct AltPathTrace {
    /// Index of this alternative path (0-based among alt paths).
    pub index: usize,
    /// Length in k-mers.
    pub length: usize,
    /// DNA sequence of the path.
    pub sequence: String,
    /// Number of leading bases shared with the reference path.
    pub shared_prefix_length: usize,
    /// Number of trailing bases shared with the reference path.
    pub shared_suffix_length: usize,
}

/// Trace of a single variant classification result.
#[derive(Debug, Clone, Serialize)]
pub struct ClassificationTrace {
    /// Which alternative path this classification corresponds to.
    pub path_index: usize,
    /// Variant type as a string (e.g., "Substitution", "Insertion").
    pub variant_type: String,
    /// Human-readable variant name (e.g., "41:A/T:41").
    pub variant_name: String,
    /// Reference allele at the variant position.
    pub ref_allele: String,
    /// Alternative allele at the variant position.
    pub alt_allele: String,
    /// Start position of the variant within the target sequence.
    pub start: usize,
    /// End position of the variant within the target sequence.
    pub end: usize,
}

/// Trace of the NNLS quantification step.
#[derive(Debug, Clone, Serialize)]
pub struct QuantificationTrace {
    /// Number of paths quantified.
    pub num_paths: usize,
    /// Number of unique k-mers across all paths.
    pub num_unique_kmers: usize,
    /// NNLS expression coefficients per path.
    pub coefficients: Vec<f64>,
    /// Relative VAFs per path.
    pub rvafs: Vec<f64>,
    /// Minimum k-mer coverage per path.
    pub min_coverages: Vec<u64>,
}

/// Build a walking trace from the walk result and reference k-mers.
pub fn build_walking_trace(
    walk: &WalkResult,
    ref_kmers: &[String],
    db: &dyn KmerDatabase,
    config: &WalkerConfig,
) -> WalkingTrace {
    let ref_set: HashSet<&str> = ref_kmers.iter().map(|s| s.as_str()).collect();

    let alt_nodes = walk
        .nodes
        .keys()
        .filter(|k| !ref_set.contains(k.as_str()))
        .count();

    let reference_kmer_counts: Vec<u64> = ref_kmers.iter().map(|k| db.query(k)).collect();

    // Estimate branching points: reference k-mers that have alt neighbors in the walk.
    // A simple heuristic: count walk nodes whose k-1 prefix matches the k-1 suffix of
    // a reference k-mer, but are not themselves reference k-mers.
    let branching_points = estimate_branching_points(walk, ref_kmers);

    let limits_hit = if walk.nodes.len() >= config.max_node {
        Some(format!(
            "max_node limit reached ({} >= {})",
            walk.nodes.len(),
            config.max_node
        ))
    } else {
        None
    };

    WalkingTrace {
        reference_kmers: ref_kmers.len(),
        reference_kmer_counts,
        total_nodes: walk.nodes.len(),
        alt_nodes,
        branching_points,
        limits_hit,
    }
}

/// Estimate branching points from the walk result.
///
/// A branching point is a reference k-mer whose (k-1) suffix is also the (k-1)
/// prefix of at least one non-reference k-mer. This indicates the walk diverged
/// from the reference at that position.
fn estimate_branching_points(walk: &WalkResult, ref_kmers: &[String]) -> usize {
    if ref_kmers.is_empty() {
        return 0;
    }

    let k = ref_kmers[0].len();
    if k < 2 {
        return 0;
    }

    let alt_kmers: Vec<&String> = walk
        .nodes
        .keys()
        .filter(|km| !walk.reference_kmers.contains(*km))
        .collect();

    let alt_prefixes: HashSet<&str> = alt_kmers
        .iter()
        .filter(|km| km.len() >= k)
        .map(|km| &km[..k - 1])
        .collect();

    ref_kmers
        .iter()
        .filter(|km| km.len() >= k)
        .filter(|km| alt_prefixes.contains(&km[1..]))
        .count()
}

/// Build a graph trace from a constructed KmerGraph.
pub fn build_graph_trace(graph: &KmerGraph) -> GraphTrace {
    let mut reference_nodes = 0usize;
    let mut alt_nodes = 0usize;
    let mut virtual_nodes = 0usize;

    for node in &graph.nodes {
        if node.kmer == "BigBang" || node.kmer == "BigCrunch" {
            virtual_nodes += 1;
        } else if node.is_reference {
            reference_nodes += 1;
        } else {
            alt_nodes += 1;
        }
    }

    let mut total_edges = 0usize;
    let mut reference_edges = 0usize;
    let mut alt_edges = 0usize;

    for edges in graph.edges.values() {
        for edge in edges {
            total_edges += 1;
            if edge.is_reference {
                reference_edges += 1;
            } else {
                alt_edges += 1;
            }
        }
    }

    GraphTrace {
        total_nodes: graph.nodes.len(),
        reference_nodes,
        alt_nodes,
        virtual_nodes,
        total_edges,
        reference_edges,
        alt_edges,
    }
}

/// Build a pathfinding trace from the paths found by the pathfinder.
pub fn build_pathfinding_trace(paths: &[KmerPath]) -> PathfindingTrace {
    if paths.is_empty() {
        return PathfindingTrace {
            paths_found: 0,
            reference_path_length: None,
            reference_path_sequence: None,
            alternative_paths: Vec::new(),
        };
    }

    // By convention, the first path is the reference path.
    let ref_path = &paths[0];
    let ref_sequence = ref_path.to_sequence();

    let alternative_paths: Vec<AltPathTrace> = paths[1..]
        .iter()
        .enumerate()
        .map(|(i, alt_path)| {
            let alt_sequence = alt_path.to_sequence();
            let (prefix_len, suffix_len) =
                compute_shared_affixes(&ref_sequence, &alt_sequence);

            AltPathTrace {
                index: i,
                length: alt_path.len(),
                sequence: alt_sequence,
                shared_prefix_length: prefix_len,
                shared_suffix_length: suffix_len,
            }
        })
        .collect();

    PathfindingTrace {
        paths_found: paths.len(),
        reference_path_length: Some(ref_path.len()),
        reference_path_sequence: Some(ref_sequence),
        alternative_paths,
    }
}

/// Compute the length of the shared prefix and suffix between two sequences.
fn compute_shared_affixes(ref_seq: &str, alt_seq: &str) -> (usize, usize) {
    let ref_bytes = ref_seq.as_bytes();
    let alt_bytes = alt_seq.as_bytes();

    let prefix_len = ref_bytes
        .iter()
        .zip(alt_bytes.iter())
        .take_while(|(r, a)| r == a)
        .count();

    // For the suffix, iterate from the end, but don't overlap with the prefix.
    let max_suffix = std::cmp::min(ref_bytes.len(), alt_bytes.len()) - prefix_len;
    let suffix_len = ref_bytes
        .iter()
        .rev()
        .zip(alt_bytes.iter().rev())
        .take(max_suffix)
        .take_while(|(r, a)| r == a)
        .count();

    (prefix_len, suffix_len)
}

/// Build a classification trace from a classification result.
pub fn build_classification_trace(
    classification: &Classification,
    path_index: usize,
) -> ClassificationTrace {
    ClassificationTrace {
        path_index,
        variant_type: classification.variant_type.to_string(),
        variant_name: classification.variant_name.clone(),
        ref_allele: classification.ref_allele.clone(),
        alt_allele: classification.alt_allele.clone(),
        start: classification.start,
        end: classification.end,
    }
}

/// Build a quantification trace from a quantification result.
pub fn build_quantification_trace(
    quant: &Quantification,
    paths: &[KmerPath],
) -> QuantificationTrace {
    // Count unique k-mers across all paths.
    let mut unique_kmers: HashSet<String> = HashSet::new();
    for path in paths {
        for kmer in &path.kmers {
            unique_kmers.insert(kmer.to_uppercase());
        }
    }

    QuantificationTrace {
        num_paths: quant.coefficients.len(),
        num_unique_kmers: unique_kmers.len(),
        coefficients: quant.coefficients.clone(),
        rvafs: quant.rvafs.clone(),
        min_coverages: quant.min_coverages.clone(),
    }
}

/// Write a JSON summary of all detection traces to `dir/summary.json`.
pub fn write_trace_summary(traces: &[DetectionTrace], dir: &Path) -> Result<()> {
    std::fs::create_dir_all(dir)?;

    let summary = TraceSummary {
        total_targets: traces.len(),
        variants_detected: traces
            .iter()
            .filter(|t| t.outcome == "variant_detected")
            .count(),
        reference_only: traces
            .iter()
            .filter(|t| t.outcome == "reference_only")
            .count(),
        no_paths: traces.iter().filter(|t| t.outcome == "no_paths").count(),
        targets: traces
            .iter()
            .map(|t| TargetSummaryEntry {
                target: t.target.clone(),
                outcome: t.outcome.clone(),
                total_nodes: t.walking.total_nodes,
                alt_nodes: t.walking.alt_nodes,
                paths_found: t.pathfinding.paths_found,
                classifications: t.classifications.len(),
            })
            .collect(),
    };

    let path = dir.join("summary.json");
    let file = std::fs::File::create(&path)?;
    serde_json::to_writer_pretty(file, &summary)?;

    Ok(())
}

/// Write a single target's detection trace to `dir/targets/<target>.json`.
pub fn write_trace_target(trace: &DetectionTrace, dir: &Path) -> Result<()> {
    let targets_dir = dir.join("targets");
    std::fs::create_dir_all(&targets_dir)?;

    // Sanitize target name for use as a filename.
    let safe_name: String = trace
        .target
        .chars()
        .map(|c| if c.is_alphanumeric() || c == '-' || c == '_' { c } else { '_' })
        .collect();

    let path = targets_dir.join(format!("{}.json", safe_name));
    let file = std::fs::File::create(&path)?;
    serde_json::to_writer_pretty(file, trace)?;

    Ok(())
}

/// Summary of all detection traces (written to summary.json).
#[derive(Debug, Clone, Serialize)]
struct TraceSummary {
    total_targets: usize,
    variants_detected: usize,
    reference_only: usize,
    no_paths: usize,
    targets: Vec<TargetSummaryEntry>,
}

/// Per-target entry in the summary.
#[derive(Debug, Clone, Serialize)]
struct TargetSummaryEntry {
    target: String,
    outcome: String,
    total_nodes: usize,
    alt_nodes: usize,
    paths_found: usize,
    classifications: usize,
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;

    /// Mock KmerDatabase for testing.
    struct MockDb {
        counts: HashMap<String, u64>,
        k: u8,
    }

    impl MockDb {
        fn new(counts: HashMap<String, u64>, k: u8) -> Self {
            Self { counts, k }
        }
    }

    impl KmerDatabase for MockDb {
        fn query(&self, kmer: &str) -> u64 {
            *self.counts.get(kmer).unwrap_or(&0)
        }

        fn kmer_length(&self) -> u8 {
            self.k
        }
    }

    fn make_walk_result(
        ref_kmers: &[&str],
        alt_kmers: &[(&str, u64)],
    ) -> (WalkResult, Vec<String>) {
        let mut nodes = HashMap::new();
        let mut reference_kmers = HashSet::new();

        for &kmer in ref_kmers {
            nodes.insert(kmer.to_string(), 100);
            reference_kmers.insert(kmer.to_string());
        }
        for &(kmer, count) in alt_kmers {
            nodes.insert(kmer.to_string(), count);
        }

        let ref_kmer_vec: Vec<String> = ref_kmers.iter().map(|s| s.to_string()).collect();

        (WalkResult { nodes, reference_kmers }, ref_kmer_vec)
    }

    #[test]
    fn test_build_walking_trace_basic() {
        let (walk, ref_kmers) = make_walk_result(
            &["ACGT", "CGTA", "GTAC"],
            &[("CGTT", 50), ("GTTC", 30)],
        );

        let mut counts = HashMap::new();
        counts.insert("ACGT".to_string(), 100u64);
        counts.insert("CGTA".to_string(), 120);
        counts.insert("GTAC".to_string(), 80);
        let db = MockDb::new(counts, 4);

        let config = WalkerConfig::default();
        let trace = build_walking_trace(&walk, &ref_kmers, &db, &config);

        assert_eq!(trace.reference_kmers, 3);
        assert_eq!(trace.reference_kmer_counts, vec![100, 120, 80]);
        assert_eq!(trace.total_nodes, 5);
        assert_eq!(trace.alt_nodes, 2);
        assert!(trace.limits_hit.is_none());
    }

    #[test]
    fn test_build_walking_trace_limit_hit() {
        let (walk, ref_kmers) = make_walk_result(&["ACGT"], &[]);
        let db = MockDb::new(HashMap::new(), 4);

        // Set max_node to 1 so the single node triggers the limit.
        let config = WalkerConfig {
            max_node: 1,
            ..WalkerConfig::default()
        };

        let trace = build_walking_trace(&walk, &ref_kmers, &db, &config);
        assert!(trace.limits_hit.is_some());
        assert!(trace.limits_hit.unwrap().contains("max_node"));
    }

    #[test]
    fn test_build_graph_trace() {
        use crate::graph::{Edge, GraphNode, KmerGraph};

        let nodes = vec![
            GraphNode { kmer: "BigBang".to_string(), count: 0, is_reference: true },
            GraphNode { kmer: "BigCrunch".to_string(), count: 0, is_reference: true },
            GraphNode { kmer: "ACGT".to_string(), count: 100, is_reference: true },
            GraphNode { kmer: "CGTA".to_string(), count: 80, is_reference: true },
            GraphNode { kmer: "CGTT".to_string(), count: 50, is_reference: false },
        ];

        let mut edges = HashMap::new();
        edges.insert(0, vec![Edge { to: 2, weight: 0.01, is_reference: true }]);
        edges.insert(2, vec![
            Edge { to: 3, weight: 0.01, is_reference: true },
            Edge { to: 4, weight: 1.0, is_reference: false },
        ]);
        edges.insert(3, vec![Edge { to: 1, weight: 0.01, is_reference: true }]);
        edges.insert(4, vec![Edge { to: 1, weight: 1.0, is_reference: false }]);

        let mut node_index = HashMap::new();
        for (i, n) in nodes.iter().enumerate() {
            node_index.insert(n.kmer.clone(), i);
        }

        let graph = KmerGraph {
            nodes,
            edges,
            node_index,
            source: 0,
            sink: 1,
        };

        let trace = build_graph_trace(&graph);
        assert_eq!(trace.total_nodes, 5);
        assert_eq!(trace.virtual_nodes, 2);
        assert_eq!(trace.reference_nodes, 2);
        assert_eq!(trace.alt_nodes, 1);
        assert_eq!(trace.total_edges, 5);
        assert_eq!(trace.reference_edges, 3);
        assert_eq!(trace.alt_edges, 2);
    }

    #[test]
    fn test_build_pathfinding_trace_empty() {
        let trace = build_pathfinding_trace(&[]);
        assert_eq!(trace.paths_found, 0);
        assert!(trace.reference_path_length.is_none());
        assert!(trace.alternative_paths.is_empty());
    }

    #[test]
    fn test_build_pathfinding_trace_with_paths() {
        let ref_path = KmerPath {
            kmers: vec!["ACGT".to_string(), "CGTA".to_string(), "GTAC".to_string()],
            is_reference: true,
        };
        // Alt path diverges at position 1: CGTT instead of CGTA.
        let alt_path = KmerPath {
            kmers: vec!["ACGT".to_string(), "CGTT".to_string(), "GTTC".to_string()],
            is_reference: false,
        };

        let paths = vec![ref_path, alt_path];
        let trace = build_pathfinding_trace(&paths);

        assert_eq!(trace.paths_found, 2);
        assert_eq!(trace.reference_path_length, Some(3));
        assert_eq!(trace.reference_path_sequence, Some("ACGTAC".to_string()));
        assert_eq!(trace.alternative_paths.len(), 1);

        let alt_trace = &trace.alternative_paths[0];
        assert_eq!(alt_trace.index, 0);
        assert_eq!(alt_trace.length, 3);
        assert_eq!(alt_trace.sequence, "ACGTTC");
        // "ACGTAC" vs "ACGTTC" — shared prefix is "ACGT" (4), shared suffix is "C" (1)
        assert_eq!(alt_trace.shared_prefix_length, 4);
        assert_eq!(alt_trace.shared_suffix_length, 1);
    }

    #[test]
    fn test_build_classification_trace() {
        let classification = Classification {
            variant_type: crate::variant::VariantType::Substitution,
            variant_name: "4:A/T:4".to_string(),
            ref_allele: "A".to_string(),
            alt_allele: "T".to_string(),
            start: 4,
            end: 4,
        };

        let trace = build_classification_trace(&classification, 0);
        assert_eq!(trace.path_index, 0);
        assert_eq!(trace.variant_type, "Substitution");
        assert_eq!(trace.variant_name, "4:A/T:4");
        assert_eq!(trace.ref_allele, "A");
        assert_eq!(trace.alt_allele, "T");
        assert_eq!(trace.start, 4);
        assert_eq!(trace.end, 4);
    }

    #[test]
    fn test_build_quantification_trace() {
        let quant = Quantification {
            coefficients: vec![100.0, 25.0],
            rvafs: vec![0.8, 0.2],
            min_coverages: vec![90, 20],
        };

        let paths = vec![
            KmerPath {
                kmers: vec!["ACGT".to_string(), "CGTA".to_string()],
                is_reference: true,
            },
            KmerPath {
                kmers: vec!["ACGT".to_string(), "CGTT".to_string()],
                is_reference: false,
            },
        ];

        let trace = build_quantification_trace(&quant, &paths);
        assert_eq!(trace.num_paths, 2);
        // ACGT, CGTA, CGTT = 3 unique k-mers
        assert_eq!(trace.num_unique_kmers, 3);
        assert_eq!(trace.coefficients, vec![100.0, 25.0]);
        assert_eq!(trace.rvafs, vec![0.8, 0.2]);
        assert_eq!(trace.min_coverages, vec![90, 20]);
    }

    #[test]
    fn test_compute_shared_affixes() {
        // Identical strings.
        let (p, s) = compute_shared_affixes("ACGTAC", "ACGTAC");
        assert_eq!(p, 6);
        assert_eq!(s, 0);

        // Completely different.
        let (p, s) = compute_shared_affixes("AAAA", "TTTT");
        assert_eq!(p, 0);
        assert_eq!(s, 0);

        // Shared prefix and suffix.
        let (p, s) = compute_shared_affixes("ACGTAC", "ACTTAC");
        assert_eq!(p, 2); // "AC"
        assert_eq!(s, 3); // "TAC"

        // Different lengths.
        let (p, s) = compute_shared_affixes("ACGT", "ACGGGGT");
        assert_eq!(p, 3); // "ACG"
        assert_eq!(s, 1); // "T"
    }

    #[test]
    fn test_write_trace_files() {
        let trace = DetectionTrace {
            target: "test_target".to_string(),
            walking: WalkingTrace {
                reference_kmers: 3,
                reference_kmer_counts: vec![100, 120, 80],
                total_nodes: 5,
                alt_nodes: 2,
                branching_points: 1,
                limits_hit: None,
            },
            graph: GraphTrace {
                total_nodes: 5,
                reference_nodes: 2,
                alt_nodes: 1,
                virtual_nodes: 2,
                total_edges: 5,
                reference_edges: 3,
                alt_edges: 2,
            },
            pathfinding: PathfindingTrace {
                paths_found: 2,
                reference_path_length: Some(3),
                reference_path_sequence: Some("ACGTAC".to_string()),
                alternative_paths: vec![AltPathTrace {
                    index: 0,
                    length: 3,
                    sequence: "ACGTTC".to_string(),
                    shared_prefix_length: 4,
                    shared_suffix_length: 1,
                }],
            },
            classifications: vec![ClassificationTrace {
                path_index: 0,
                variant_type: "Substitution".to_string(),
                variant_name: "4:A/T:4".to_string(),
                ref_allele: "A".to_string(),
                alt_allele: "T".to_string(),
                start: 4,
                end: 4,
            }],
            quantification: Some(QuantificationTrace {
                num_paths: 2,
                num_unique_kmers: 3,
                coefficients: vec![100.0, 25.0],
                rvafs: vec![0.8, 0.2],
                min_coverages: vec![90, 20],
            }),
            outcome: "variant_detected".to_string(),
        };

        let dir = tempfile::tempdir().unwrap();
        let dir_path = dir.path();

        // Write summary.
        write_trace_summary(&[trace.clone()], dir_path).unwrap();
        let summary_path = dir_path.join("summary.json");
        assert!(summary_path.exists());

        let summary_text = std::fs::read_to_string(&summary_path).unwrap();
        let summary: serde_json::Value = serde_json::from_str(&summary_text).unwrap();
        assert_eq!(summary["total_targets"], 1);
        assert_eq!(summary["variants_detected"], 1);

        // Write target trace.
        write_trace_target(&trace, dir_path).unwrap();
        let target_path = dir_path.join("targets/test_target.json");
        assert!(target_path.exists());

        let target_text = std::fs::read_to_string(&target_path).unwrap();
        let target: serde_json::Value = serde_json::from_str(&target_text).unwrap();
        assert_eq!(target["target"], "test_target");
        assert_eq!(target["outcome"], "variant_detected");
    }

    #[test]
    fn test_write_trace_target_sanitizes_name() {
        let trace = DetectionTrace {
            target: "chr1:1234-5678/variant".to_string(),
            walking: WalkingTrace {
                reference_kmers: 0,
                reference_kmer_counts: vec![],
                total_nodes: 0,
                alt_nodes: 0,
                branching_points: 0,
                limits_hit: None,
            },
            graph: GraphTrace {
                total_nodes: 0,
                reference_nodes: 0,
                alt_nodes: 0,
                virtual_nodes: 0,
                total_edges: 0,
                reference_edges: 0,
                alt_edges: 0,
            },
            pathfinding: PathfindingTrace {
                paths_found: 0,
                reference_path_length: None,
                reference_path_sequence: None,
                alternative_paths: vec![],
            },
            classifications: vec![],
            quantification: None,
            outcome: "no_paths".to_string(),
        };

        let dir = tempfile::tempdir().unwrap();
        write_trace_target(&trace, dir.path()).unwrap();

        // The file should exist with sanitized name.
        let expected = dir.path().join("targets/chr1_1234-5678_variant.json");
        assert!(expected.exists());
    }
}
