/// DOT format (Graphviz) exporter for KmerGraph.
///
/// Generates `.dot` files that can be rendered with `dot -Tpng graph.dot -o graph.png`
/// or similar Graphviz commands. Supports full graph export with node/edge styling,
/// simplified export for large graphs (collapsing linear chains), and a filter
/// flow diagram for visualizing the filtering cascade.

use std::collections::{HashMap, HashSet};
use std::fmt::Write as FmtWrite;
use std::path::Path;

use anyhow::Result;

use super::KmerGraph;
use crate::filter::{FilterConfig, FilterResult};
use crate::sequence::path::KmerPath;

/// Truncate a k-mer string for display: first 5 + "..." + last 5 if longer than 15 chars.
fn truncate_kmer(kmer: &str) -> String {
    if kmer.len() > 15 {
        format!("{}...{}", &kmer[..5], &kmer[kmer.len() - 5..])
    } else {
        kmer.to_string()
    }
}

/// Export a KmerGraph to DOT format with full styling.
///
/// Node styling:
/// - BigBang/BigCrunch: shape=point, width=0.2, color=gray
/// - Reference nodes: style=filled, fillcolor light green, label with truncated k-mer and count
/// - Alt nodes: style=filled, fillcolor light orange, same label format
///
/// Edge styling:
/// - Reference edges: color=green, penwidth=1
/// - Alt edges: color=red, penwidth=2, style=dashed
///
/// If `highlight_paths` is provided, edges on those paths are highlighted in blue.
pub fn export_dot(
    graph: &KmerGraph,
    target_name: &str,
    highlight_paths: Option<&[&KmerPath]>,
) -> String {
    let mut dot = String::new();

    writeln!(
        dot,
        "digraph \"{}\" {{",
        escape_dot_string(target_name)
    )
    .unwrap();
    writeln!(dot, "    rankdir=LR;").unwrap();
    writeln!(
        dot,
        "    label=\"{}\";",
        escape_dot_string(target_name)
    )
    .unwrap();
    writeln!(dot, "    fontname=\"Courier\";").unwrap();
    writeln!(dot, "    node [fontname=\"Courier\"];").unwrap();
    writeln!(dot, "    edge [fontname=\"Courier\"];").unwrap();
    writeln!(dot).unwrap();

    // Build set of highlighted edges (from_idx, to_idx) if paths are provided.
    let highlighted_edges = build_highlighted_edges(graph, highlight_paths);

    // Emit nodes.
    for (idx, node) in graph.nodes.iter().enumerate() {
        write!(dot, "    n{} [", idx).unwrap();
        if node.kmer == "BigBang" || node.kmer == "BigCrunch" {
            write!(
                dot,
                "shape=point, width=0.2, color=gray, label=\"\"",
            )
            .unwrap();
        } else {
            let label = format!("{}\\n{}", truncate_kmer(&node.kmer), node.count);
            if node.is_reference {
                write!(
                    dot,
                    "style=filled, fillcolor=\"#90EE90\", label=\"{}\"",
                    label
                )
                .unwrap();
            } else {
                write!(
                    dot,
                    "style=filled, fillcolor=\"#FFB347\", label=\"{}\"",
                    label
                )
                .unwrap();
            }
        }
        writeln!(dot, "];").unwrap();
    }

    writeln!(dot).unwrap();

    // Emit edges.
    for (&from_idx, edges) in &graph.edges {
        for edge in edges {
            let is_highlighted = highlighted_edges.contains(&(from_idx, edge.to));

            write!(dot, "    n{} -> n{} [", from_idx, edge.to).unwrap();

            if is_highlighted {
                write!(dot, "color=blue, penwidth=3, style=bold").unwrap();
            } else if edge.is_reference {
                write!(dot, "color=green, penwidth=1").unwrap();
            } else {
                write!(dot, "color=red, penwidth=2, style=dashed").unwrap();
            }

            writeln!(dot, "];").unwrap();
        }
    }

    writeln!(dot, "}}").unwrap();

    dot
}

/// Build the set of (from, to) edge pairs that should be highlighted.
fn build_highlighted_edges(
    graph: &KmerGraph,
    highlight_paths: Option<&[&KmerPath]>,
) -> HashSet<(usize, usize)> {
    let mut edges = HashSet::new();

    let paths = match highlight_paths {
        Some(p) => p,
        None => return edges,
    };

    for path in paths {
        // Map k-mer strings back to node indices and record consecutive pairs.
        let indices: Vec<Option<usize>> = path
            .kmers
            .iter()
            .map(|kmer| graph.node_index.get(kmer).copied())
            .collect();

        // Add source -> first k-mer edge.
        if let Some(Some(first_idx)) = indices.first() {
            edges.insert((graph.source, *first_idx));
        }

        // Add consecutive k-mer edges.
        for window in indices.windows(2) {
            if let (Some(from), Some(to)) = (window[0], window[1]) {
                edges.insert((from, to));
            }
        }

        // Add last k-mer -> sink edge.
        if let Some(Some(last_idx)) = indices.last() {
            edges.insert((*last_idx, graph.sink));
        }
    }

    edges
}

/// Export a simplified DOT representation for large graphs (>100 nodes).
///
/// Linear chains (sequences of nodes each with exactly 1 incoming and 1 outgoing
/// edge) are collapsed into a single node labeled with the chain length and
/// the first/last k-mers.
pub fn export_dot_simplified(graph: &KmerGraph, target_name: &str) -> String {
    // Build in-degree and out-degree maps (excluding virtual nodes from chain detection).
    let mut in_degree: HashMap<usize, usize> = HashMap::new();
    let mut out_degree: HashMap<usize, usize> = HashMap::new();

    for (&from_idx, edges) in &graph.edges {
        let edge_count = edges.len();
        *out_degree.entry(from_idx).or_default() += edge_count;
        for edge in edges {
            *in_degree.entry(edge.to).or_default() += 1;
        }
    }

    // Identify chain-eligible nodes: non-virtual nodes with exactly 1 in-edge and 1 out-edge.
    let is_chain_node = |idx: usize| -> bool {
        if idx == graph.source || idx == graph.sink {
            return false;
        }
        in_degree.get(&idx).copied().unwrap_or(0) == 1
            && out_degree.get(&idx).copied().unwrap_or(0) == 1
    };

    // Find chains by walking from non-chain nodes through chain nodes.
    let mut visited: HashSet<usize> = HashSet::new();
    let mut chains: Vec<Vec<usize>> = Vec::new();

    // For each edge, if the target is a chain node and hasn't been visited,
    // walk forward to collect the full chain.
    for (_from_idx, edges) in &graph.edges {
        for edge in edges {
            let start = edge.to;
            if visited.contains(&start) || !is_chain_node(start) {
                continue;
            }

            let mut chain = vec![start];
            visited.insert(start);

            let mut current = start;
            loop {
                // Get the single outgoing edge.
                let next = match graph.edges.get(&current) {
                    Some(out_edges) if out_edges.len() == 1 => out_edges[0].to,
                    _ => break,
                };

                if !is_chain_node(next) || visited.contains(&next) {
                    break;
                }

                visited.insert(next);
                chain.push(next);
                current = next;
            }

            if chain.len() >= 2 {
                chains.push(chain);
            } else {
                // Single-node "chains" are not collapsed.
                visited.remove(&start);
            }
        }
    }

    // Build the collapsed representation.
    // Nodes in chains are replaced by a single "collapsed" node.
    // Map from chain member index -> chain id.
    let mut chain_membership: HashMap<usize, usize> = HashMap::new();
    for (chain_id, chain) in chains.iter().enumerate() {
        for &node_idx in chain {
            chain_membership.insert(node_idx, chain_id);
        }
    }

    let mut dot = String::new();

    writeln!(
        dot,
        "digraph \"{}\" {{",
        escape_dot_string(target_name)
    )
    .unwrap();
    writeln!(dot, "    rankdir=LR;").unwrap();
    writeln!(
        dot,
        "    label=\"{} (simplified)\";",
        escape_dot_string(target_name)
    )
    .unwrap();
    writeln!(dot, "    fontname=\"Courier\";").unwrap();
    writeln!(dot, "    node [fontname=\"Courier\"];").unwrap();
    writeln!(dot, "    edge [fontname=\"Courier\"];").unwrap();
    writeln!(dot).unwrap();

    // Emit non-chain nodes.
    for (idx, node) in graph.nodes.iter().enumerate() {
        if chain_membership.contains_key(&idx) {
            continue;
        }

        write!(dot, "    n{} [", idx).unwrap();
        if node.kmer == "BigBang" || node.kmer == "BigCrunch" {
            write!(dot, "shape=point, width=0.2, color=gray, label=\"\"").unwrap();
        } else {
            let label = format!("{}\\n{}", truncate_kmer(&node.kmer), node.count);
            if node.is_reference {
                write!(
                    dot,
                    "style=filled, fillcolor=\"#90EE90\", label=\"{}\"",
                    label
                )
                .unwrap();
            } else {
                write!(
                    dot,
                    "style=filled, fillcolor=\"#FFB347\", label=\"{}\"",
                    label
                )
                .unwrap();
            }
        }
        writeln!(dot, "];").unwrap();
    }

    // Emit collapsed chain nodes.
    for (chain_id, chain) in chains.iter().enumerate() {
        let first_kmer = truncate_kmer(&graph.nodes[chain[0]].kmer);
        let last_kmer = truncate_kmer(&graph.nodes[*chain.last().unwrap()].kmer);
        let is_ref = chain.iter().all(|&idx| graph.nodes[idx].is_reference);
        let fillcolor = if is_ref { "#90EE90" } else { "#FFB347" };

        let label = format!(
            "{} nodes\\n{}\\n...\\n{}",
            chain.len(),
            first_kmer,
            last_kmer
        );

        writeln!(
            dot,
            "    chain{} [shape=box, style=filled, fillcolor=\"{}\", label=\"{}\"];",
            chain_id, fillcolor, label
        )
        .unwrap();
    }

    writeln!(dot).unwrap();

    // Emit edges, remapping chain nodes.
    let mut emitted_edges: HashSet<(String, String)> = HashSet::new();

    for (&from_idx, edges) in &graph.edges {
        for edge in edges {
            let from_label = if let Some(&chain_id) = chain_membership.get(&from_idx) {
                // Only emit edge from chain if from_idx is the last node in the chain.
                let chain = &chains[chain_id];
                if from_idx != *chain.last().unwrap() {
                    continue;
                }
                format!("chain{}", chain_id)
            } else {
                format!("n{}", from_idx)
            };

            let to_label = if let Some(&chain_id) = chain_membership.get(&edge.to) {
                // Only emit edge to chain if edge.to is the first node in the chain.
                let chain = &chains[chain_id];
                if edge.to != chain[0] {
                    continue;
                }
                format!("chain{}", chain_id)
            } else {
                format!("n{}", edge.to)
            };

            let edge_key = (from_label.clone(), to_label.clone());
            if emitted_edges.contains(&edge_key) {
                continue;
            }
            emitted_edges.insert(edge_key);

            write!(dot, "    {} -> {} [", from_label, to_label).unwrap();
            if edge.is_reference {
                write!(dot, "color=green, penwidth=1").unwrap();
            } else {
                write!(dot, "color=red, penwidth=2, style=dashed").unwrap();
            }
            writeln!(dot, "];").unwrap();
        }
    }

    writeln!(dot, "}}").unwrap();

    dot
}

/// Export a top-to-bottom flow diagram showing the filtering cascade.
///
/// The diagram shows:
/// - Input node: "N expected variants"
/// - One filter node per criterion, showing passed/rejected counts
/// - Rejected variants as red side nodes
/// - Output node: "N PASS / M filtered"
pub fn export_filter_flow_dot(
    filter_results: &[FilterResult],
    config: &FilterConfig,
) -> String {
    let mut dot = String::new();

    writeln!(dot, "digraph filter_flow {{").unwrap();
    writeln!(dot, "    rankdir=TB;").unwrap();
    writeln!(dot, "    label=\"Filter Cascade\";").unwrap();
    writeln!(dot, "    fontname=\"Courier\";").unwrap();
    writeln!(dot, "    node [fontname=\"Courier\", shape=box];").unwrap();
    writeln!(dot, "    edge [fontname=\"Courier\"];").unwrap();
    writeln!(dot).unwrap();

    let total = filter_results.len();

    // Input node.
    writeln!(
        dot,
        "    input [label=\"{} expected variants\", style=filled, fillcolor=\"#87CEEB\"];",
        total
    )
    .unwrap();

    // Categorize results by filter failure type.
    let mut not_detected = 0usize;
    let mut failed_coverage = 0usize;
    let mut failed_vaf = 0usize;
    let mut failed_expression = 0usize;
    let mut failed_type = 0usize;
    let mut passed = 0usize;

    for result in filter_results {
        if result.found == "Found" {
            passed += 1;
        } else {
            let notes = &result.filter_notes;
            if notes.contains("No matching variant detected") {
                not_detected += 1;
            } else {
                // A result can fail multiple criteria; count each.
                if notes.contains("coverage") {
                    failed_coverage += 1;
                }
                if notes.contains("VAF") {
                    failed_vaf += 1;
                }
                if notes.contains("expression") {
                    failed_expression += 1;
                }
                if notes.contains("type") {
                    failed_type += 1;
                }
            }
        }
    }

    // Running count of variants still passing after each filter.
    let after_detection = total - not_detected;

    // Detection filter.
    writeln!(
        dot,
        "    detect [label=\"Detection\\n{} passed\", style=filled, fillcolor=\"#90EE90\"];",
        after_detection
    )
    .unwrap();
    writeln!(dot, "    input -> detect;").unwrap();

    if not_detected > 0 {
        writeln!(
            dot,
            "    reject_detect [label=\"{} not detected\", style=filled, fillcolor=\"#FFB347\"];",
            not_detected
        )
        .unwrap();
        writeln!(dot, "    detect -> reject_detect [color=red, style=dashed];").unwrap();
    }

    // Coverage filter.
    let after_coverage = after_detection - failed_coverage;
    writeln!(
        dot,
        "    coverage [label=\"Coverage >= {}\\n{} passed\", style=filled, fillcolor=\"#90EE90\"];",
        config.min_coverage, after_coverage
    )
    .unwrap();
    writeln!(dot, "    detect -> coverage;").unwrap();

    if failed_coverage > 0 {
        writeln!(
            dot,
            "    reject_cov [label=\"{} failed coverage\", style=filled, fillcolor=\"#FFB347\"];",
            failed_coverage
        )
        .unwrap();
        writeln!(dot, "    coverage -> reject_cov [color=red, style=dashed];").unwrap();
    }

    // VAF filter.
    let after_vaf = after_coverage - failed_vaf;
    writeln!(
        dot,
        "    vaf [label=\"VAF >= {:.4}\\n{} passed\", style=filled, fillcolor=\"#90EE90\"];",
        config.min_vaf, after_vaf
    )
    .unwrap();
    writeln!(dot, "    coverage -> vaf;").unwrap();

    if failed_vaf > 0 {
        writeln!(
            dot,
            "    reject_vaf [label=\"{} failed VAF\", style=filled, fillcolor=\"#FFB347\"];",
            failed_vaf
        )
        .unwrap();
        writeln!(dot, "    vaf -> reject_vaf [color=red, style=dashed];").unwrap();
    }

    // Expression filter.
    let after_expr = after_vaf - failed_expression;
    writeln!(
        dot,
        "    expression [label=\"Expression >= {:.2}\\n{} passed\", style=filled, fillcolor=\"#90EE90\"];",
        config.min_expression, after_expr
    )
    .unwrap();
    writeln!(dot, "    vaf -> expression;").unwrap();

    if failed_expression > 0 {
        writeln!(
            dot,
            "    reject_expr [label=\"{} failed expression\", style=filled, fillcolor=\"#FFB347\"];",
            failed_expression
        )
        .unwrap();
        writeln!(dot, "    expression -> reject_expr [color=red, style=dashed];").unwrap();
    }

    // Type filter.
    if !config.types.is_empty() {
        let after_type = after_expr - failed_type;
        writeln!(
            dot,
            "    type_filter [label=\"Type in [{}]\\n{} passed\", style=filled, fillcolor=\"#90EE90\"];",
            config.types.join(", "),
            after_type
        )
        .unwrap();
        writeln!(dot, "    expression -> type_filter;").unwrap();

        if failed_type > 0 {
            writeln!(
                dot,
                "    reject_type [label=\"{} failed type\", style=filled, fillcolor=\"#FFB347\"];",
                failed_type
            )
            .unwrap();
            writeln!(dot, "    type_filter -> reject_type [color=red, style=dashed];").unwrap();
        }

        // Output node.
        let filtered = total - passed;
        writeln!(
            dot,
            "    output [label=\"{} PASS / {} filtered\", style=filled, fillcolor=\"#87CEEB\"];",
            passed, filtered
        )
        .unwrap();
        writeln!(dot, "    type_filter -> output;").unwrap();
    } else {
        // No type filter; connect expression directly to output.
        let filtered = total - passed;
        writeln!(
            dot,
            "    output [label=\"{} PASS / {} filtered\", style=filled, fillcolor=\"#87CEEB\"];",
            passed, filtered
        )
        .unwrap();
        writeln!(dot, "    expression -> output;").unwrap();
    }

    writeln!(dot, "}}").unwrap();

    dot
}

/// Write a DOT string to a file.
pub fn write_dot_file(dot: &str, path: &Path) -> Result<()> {
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent)?;
    }
    std::fs::write(path, dot)?;
    Ok(())
}

/// Escape a string for use inside DOT labels (double quotes).
fn escape_dot_string(s: &str) -> String {
    s.replace('\\', "\\\\")
        .replace('"', "\\\"")
        .replace('\n', "\\n")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::{Edge, GraphNode, KmerGraph};

    /// Build a small test graph with 2 reference nodes and 1 alt node.
    fn make_test_graph() -> KmerGraph {
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

        KmerGraph {
            nodes,
            edges,
            node_index,
            source: 0,
            sink: 1,
        }
    }

    #[test]
    fn test_truncate_kmer_short() {
        assert_eq!(truncate_kmer("ACGT"), "ACGT");
        assert_eq!(truncate_kmer("ACGTACGTACGTACG"), "ACGTACGTACGTACG"); // 15 chars, not truncated
    }

    #[test]
    fn test_truncate_kmer_long() {
        let long_kmer = "ACGTACGTACGTACGT"; // 16 chars
        assert_eq!(truncate_kmer(long_kmer), "ACGTA...TACGT");
    }

    #[test]
    fn test_export_dot_basic() {
        let graph = make_test_graph();
        let dot = export_dot(&graph, "test_target", None);

        // Verify DOT structure.
        assert!(dot.starts_with("digraph"));
        assert!(dot.contains("rankdir=LR"));
        assert!(dot.contains("test_target"));
        assert!(dot.ends_with("}\n"));

        // Verify node styling.
        assert!(dot.contains("shape=point")); // BigBang/BigCrunch
        assert!(dot.contains("#90EE90")); // Reference node color
        assert!(dot.contains("#FFB347")); // Alt node color

        // Verify edge styling.
        assert!(dot.contains("color=green")); // Reference edges
        assert!(dot.contains("color=red")); // Alt edges
        assert!(dot.contains("style=dashed")); // Alt edge style
    }

    #[test]
    fn test_export_dot_with_highlight() {
        let graph = make_test_graph();
        let path = KmerPath {
            kmers: vec!["ACGT".to_string(), "CGTT".to_string()],
            is_reference: false,
        };
        let paths = vec![&path];
        let dot = export_dot(&graph, "test_target", Some(&paths));

        // Highlighted path should produce blue edges.
        assert!(dot.contains("color=blue"));
    }

    #[test]
    fn test_export_dot_simplified() {
        // Build a graph with a linear chain: BigBang -> A -> B -> C -> D -> BigCrunch
        // where B and C form a chain (1 in, 1 out each).
        let nodes = vec![
            GraphNode { kmer: "BigBang".to_string(), count: 0, is_reference: true },
            GraphNode { kmer: "BigCrunch".to_string(), count: 0, is_reference: true },
            GraphNode { kmer: "AAAA".to_string(), count: 100, is_reference: true },
            GraphNode { kmer: "AAAB".to_string(), count: 90, is_reference: true },
            GraphNode { kmer: "AABC".to_string(), count: 80, is_reference: true },
            GraphNode { kmer: "ABCD".to_string(), count: 70, is_reference: true },
        ];

        let mut edges = HashMap::new();
        edges.insert(0, vec![Edge { to: 2, weight: 0.01, is_reference: true }]);
        edges.insert(2, vec![Edge { to: 3, weight: 0.01, is_reference: true }]);
        edges.insert(3, vec![Edge { to: 4, weight: 0.01, is_reference: true }]);
        edges.insert(4, vec![Edge { to: 5, weight: 0.01, is_reference: true }]);
        edges.insert(5, vec![Edge { to: 1, weight: 0.01, is_reference: true }]);

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

        let dot = export_dot_simplified(&graph, "chain_test");

        assert!(dot.contains("simplified"));
        assert!(dot.starts_with("digraph"));
        assert!(dot.ends_with("}\n"));
    }

    #[test]
    fn test_export_filter_flow_dot() {
        let results = vec![
            FilterResult {
                sample: "s1".to_string(),
                chrom: "chr1".to_string(),
                pos: 100,
                ref_allele: "A".to_string(),
                alt_allele: "T".to_string(),
                variant_type: "Substitution".to_string(),
                found: "Found".to_string(),
                filter_notes: "PASS".to_string(),
                kmer_vaf: Some(0.1),
                kmer_min_coverage: Some(50),
                kmer_expression: Some(100.0),
                ref_sequence: None,
                variant_sequence: None,
            },
            FilterResult {
                sample: "s1".to_string(),
                chrom: "chr2".to_string(),
                pos: 200,
                ref_allele: "G".to_string(),
                alt_allele: "C".to_string(),
                variant_type: "Substitution".to_string(),
                found: "Not Found".to_string(),
                filter_notes: "No matching variant detected".to_string(),
                kmer_vaf: None,
                kmer_min_coverage: None,
                kmer_expression: None,
                ref_sequence: None,
                variant_sequence: None,
            },
            FilterResult {
                sample: "s1".to_string(),
                chrom: "chr3".to_string(),
                pos: 300,
                ref_allele: "T".to_string(),
                alt_allele: "A".to_string(),
                variant_type: "Substitution".to_string(),
                found: "Not Found".to_string(),
                filter_notes: "coverage 2 < 5".to_string(),
                kmer_vaf: Some(0.01),
                kmer_min_coverage: Some(2),
                kmer_expression: Some(5.0),
                ref_sequence: None,
                variant_sequence: None,
            },
        ];

        let config = FilterConfig {
            min_coverage: 5,
            min_vaf: 0.001,
            min_expression: 1.0,
            use_alt: false,
            types: vec![],
        };

        let dot = export_filter_flow_dot(&results, &config);

        assert!(dot.contains("filter_flow"));
        assert!(dot.contains("3 expected variants"));
        assert!(dot.contains("1 not detected"));
        assert!(dot.contains("1 failed coverage"));
        assert!(dot.contains("1 PASS / 2 filtered"));
        assert!(dot.starts_with("digraph"));
        assert!(dot.ends_with("}\n"));
    }

    #[test]
    fn test_export_filter_flow_with_type_filter() {
        let results = vec![FilterResult {
            sample: "s1".to_string(),
            chrom: "chr1".to_string(),
            pos: 100,
            ref_allele: "A".to_string(),
            alt_allele: "T".to_string(),
            variant_type: "Substitution".to_string(),
            found: "Found".to_string(),
            filter_notes: "PASS".to_string(),
            kmer_vaf: Some(0.1),
            kmer_min_coverage: Some(50),
            kmer_expression: Some(100.0),
            ref_sequence: None,
            variant_sequence: None,
        }];

        let config = FilterConfig {
            min_coverage: 5,
            min_vaf: 0.001,
            min_expression: 1.0,
            use_alt: false,
            types: vec!["Substitution".to_string(), "Insertion".to_string()],
        };

        let dot = export_filter_flow_dot(&results, &config);

        assert!(dot.contains("Type in [Substitution, Insertion]"));
        assert!(dot.contains("type_filter -> output"));
    }

    #[test]
    fn test_write_dot_file() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("test.dot");

        let dot = "digraph test { a -> b; }\n";
        write_dot_file(dot, &path).unwrap();

        let content = std::fs::read_to_string(&path).unwrap();
        assert_eq!(content, dot);
    }

    #[test]
    fn test_write_dot_file_creates_parent_dirs() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("sub/dir/test.dot");

        let dot = "digraph test { a -> b; }\n";
        write_dot_file(dot, &path).unwrap();

        assert!(path.exists());
    }

    #[test]
    fn test_escape_dot_string() {
        assert_eq!(escape_dot_string("hello"), "hello");
        assert_eq!(escape_dot_string("say \"hi\""), "say \\\"hi\\\"");
        assert_eq!(escape_dot_string("line1\nline2"), "line1\\nline2");
    }

    #[test]
    fn test_export_dot_empty_graph() {
        let nodes = vec![
            GraphNode { kmer: "BigBang".to_string(), count: 0, is_reference: true },
            GraphNode { kmer: "BigCrunch".to_string(), count: 0, is_reference: true },
        ];

        let mut node_index = HashMap::new();
        node_index.insert("BigBang".to_string(), 0);
        node_index.insert("BigCrunch".to_string(), 1);

        let graph = KmerGraph {
            nodes,
            edges: HashMap::new(),
            node_index,
            source: 0,
            sink: 1,
        };

        let dot = export_dot(&graph, "empty", None);
        assert!(dot.contains("digraph"));
        assert!(dot.contains("n0"));
        assert!(dot.contains("n1"));
    }

    #[test]
    fn test_export_dot_long_kmer_truncation() {
        let long_kmer = "ACGTACGTACGTACGTACGTACGTACGTACGT"; // 32 chars
        let nodes = vec![
            GraphNode { kmer: "BigBang".to_string(), count: 0, is_reference: true },
            GraphNode { kmer: "BigCrunch".to_string(), count: 0, is_reference: true },
            GraphNode { kmer: long_kmer.to_string(), count: 42, is_reference: true },
        ];

        let mut node_index = HashMap::new();
        for (i, n) in nodes.iter().enumerate() {
            node_index.insert(n.kmer.clone(), i);
        }

        let graph = KmerGraph {
            nodes,
            edges: HashMap::new(),
            node_index,
            source: 0,
            sink: 1,
        };

        let dot = export_dot(&graph, "truncation_test", None);

        // Should contain the truncated form, not the full 32-char k-mer.
        assert!(dot.contains("ACGTA...TACGT"));
        assert!(!dot.contains(long_kmer));
    }
}
