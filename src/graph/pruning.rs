/// Graph pruning pipeline for cleaning up k-mer graphs before pathfinding.
///
/// The k-mer graph often contains noise from sequencing errors: dead-end nodes,
/// low-coverage spurious paths, and short tip branches. Pruning these before
/// pathfinding reduces false positive variant calls and speeds up path enumeration.
///
/// The pipeline runs in stages:
/// 1. Remove dead-end nodes (unreachable from source or sink)
/// 2. Remove low-coverage non-reference nodes
/// 3. Remove dead-ends created by step 2
/// 4. Clip short tips (dead-end branches shorter than a threshold)
/// 5. Collapse bubbles (remove low-coverage parallel paths)

use std::collections::{HashSet, VecDeque};

use super::KmerGraph;

/// Configuration for graph pruning.
#[derive(Debug, Clone)]
pub struct PruneConfig {
    /// Minimum node count to survive (0 = use adaptive threshold).
    pub min_count: u64,
    /// Error rate for adaptive threshold (default 0.001).
    pub error_rate: f64,
    /// Enable tip clipping.
    pub clip_tips: bool,
    /// Maximum tip length as fraction of k (default 0.5).
    pub max_tip_fraction: f64,
    /// Enable bubble collapsing.
    pub collapse_bubbles: bool,
    /// Min coverage ratio for bubble retention (default 0.10).
    pub bubble_coverage_ratio: f64,
}

impl Default for PruneConfig {
    fn default() -> Self {
        Self {
            min_count: 0,
            error_rate: 0.001,
            clip_tips: true,
            max_tip_fraction: 0.5,
            collapse_bubbles: true,
            bubble_coverage_ratio: 0.10,
        }
    }
}

/// Run the full pruning pipeline on a graph.
/// Returns the number of nodes removed.
pub fn prune_graph(graph: &mut KmerGraph, config: &PruneConfig, k: u8) -> usize {
    let initial = count_real_nodes(graph);

    // Step 1: Remove dead-end nodes (unreachable from source or sink).
    remove_dead_ends(graph);

    // Step 2: Remove low-coverage nodes.
    let min_count = if config.min_count > 0 {
        config.min_count
    } else {
        compute_adaptive_min_count(graph, config.error_rate)
    };
    remove_low_coverage(graph, min_count);

    // Step 3: Remove dead-ends again (may have been created by step 2).
    remove_dead_ends(graph);

    // Step 4: Clip tips (short dead-end branches).
    if config.clip_tips {
        let max_tip_len = ((k as f64) * config.max_tip_fraction).ceil() as usize;
        clip_tips(graph, max_tip_len);
    }

    // Step 5: Collapse bubbles.
    if config.collapse_bubbles {
        collapse_bubbles(graph, config.bubble_coverage_ratio);
    }

    initial - count_real_nodes(graph)
}

/// Count real (non-virtual) nodes that are still present in the node_index.
fn count_real_nodes(graph: &KmerGraph) -> usize {
    graph
        .node_index
        .keys()
        .filter(|k| *k != "BigBang" && *k != "BigCrunch")
        .count()
}

/// Remove a node by index: clear it from node_index and remove all edges to/from it.
///
/// Never removes the source or sink nodes.
fn remove_node(graph: &mut KmerGraph, idx: usize) {
    if idx == graph.source || idx == graph.sink {
        return;
    }

    // Remove from node_index lookup.
    let kmer = graph.nodes[idx].kmer.clone();
    graph.node_index.remove(&kmer);

    // Remove all outgoing edges from this node.
    graph.edges.remove(&idx);

    // Remove all incoming edges pointing to this node.
    for edges in graph.edges.values_mut() {
        edges.retain(|e| e.to != idx);
    }
}

/// BFS forward from source, returning the set of reachable node indices.
fn bfs_forward(graph: &KmerGraph) -> HashSet<usize> {
    let mut visited = HashSet::new();
    let mut queue = VecDeque::new();

    visited.insert(graph.source);
    queue.push_back(graph.source);

    while let Some(node) = queue.pop_front() {
        if let Some(edges) = graph.edges.get(&node) {
            for edge in edges {
                if visited.insert(edge.to) {
                    queue.push_back(edge.to);
                }
            }
        }
    }

    visited
}

/// BFS backward from sink (following reverse edges), returning reachable node indices.
fn bfs_backward(graph: &KmerGraph) -> HashSet<usize> {
    // Build reverse adjacency: for each edge a->b, record b->a.
    let mut reverse: std::collections::HashMap<usize, Vec<usize>> =
        std::collections::HashMap::new();
    for (&from, edges) in &graph.edges {
        for edge in edges {
            reverse.entry(edge.to).or_default().push(from);
        }
    }

    let mut visited = HashSet::new();
    let mut queue = VecDeque::new();

    visited.insert(graph.sink);
    queue.push_back(graph.sink);

    while let Some(node) = queue.pop_front() {
        if let Some(predecessors) = reverse.get(&node) {
            for &pred in predecessors {
                if visited.insert(pred) {
                    queue.push_back(pred);
                }
            }
        }
    }

    visited
}

/// Remove nodes that are not reachable from both source (forward) and sink (backward).
///
/// A node must be on some source-to-sink path to be useful for pathfinding.
/// Dead-end nodes are removed along with their edges.
pub fn remove_dead_ends(graph: &mut KmerGraph) {
    let forward_reachable = bfs_forward(graph);
    let backward_reachable = bfs_backward(graph);

    // Nodes to remove: present in node_index but not in both reachable sets.
    let to_remove: Vec<usize> = graph
        .node_index
        .values()
        .copied()
        .filter(|&idx| {
            idx != graph.source
                && idx != graph.sink
                && (!forward_reachable.contains(&idx) || !backward_reachable.contains(&idx))
        })
        .collect();

    for idx in to_remove {
        remove_node(graph, idx);
    }
}

/// Compute an adaptive minimum count threshold based on median reference coverage.
///
/// The idea: sequencing errors produce k-mers at roughly `error_rate * coverage`.
/// We set the threshold above this expected error level.
///
/// Returns max(2, ceil(median_ref_coverage * error_rate)).
pub fn compute_adaptive_min_count(graph: &KmerGraph, error_rate: f64) -> u64 {
    let mut ref_counts: Vec<u64> = graph
        .nodes
        .iter()
        .enumerate()
        .filter(|(idx, node)| {
            *idx != graph.source
                && *idx != graph.sink
                && node.is_reference
                && graph.node_index.contains_key(&node.kmer)
        })
        .map(|(_, node)| node.count)
        .collect();

    if ref_counts.is_empty() {
        return 2;
    }

    ref_counts.sort_unstable();
    let median = ref_counts[ref_counts.len() / 2];

    let threshold = (median as f64 * error_rate).ceil() as u64;
    std::cmp::max(2, threshold)
}

/// Remove non-reference nodes whose count is below `min_count`.
///
/// Reference nodes, source, and sink are never removed.
pub fn remove_low_coverage(graph: &mut KmerGraph, min_count: u64) {
    let to_remove: Vec<usize> = graph
        .node_index
        .values()
        .copied()
        .filter(|&idx| {
            idx != graph.source
                && idx != graph.sink
                && !graph.nodes[idx].is_reference
                && graph.nodes[idx].count < min_count
        })
        .collect();

    for idx in to_remove {
        remove_node(graph, idx);
    }
}

/// Clip short tip branches from the graph.
///
/// A "tip" is a dead-end branch: a chain of nodes ending at a leaf (no outgoing
/// edges, excluding the sink). If the chain is shorter than `max_tip_len` and
/// branches off a node with >= 2 successors, the entire chain is removed.
///
/// Reference nodes are never clipped.
pub fn clip_tips(graph: &mut KmerGraph, max_tip_len: usize) {
    // We may need multiple rounds as clipping one tip can expose another.
    loop {
        let mut clipped_any = false;

        // Find leaf nodes: nodes with no outgoing edges (excluding sink).
        let active_indices: Vec<usize> = graph
            .node_index
            .values()
            .copied()
            .filter(|&idx| idx != graph.source && idx != graph.sink)
            .collect();

        let leaves: Vec<usize> = active_indices
            .iter()
            .copied()
            .filter(|&idx| {
                // A leaf has no outgoing edges (or only edges to removed nodes).
                let out_count = graph
                    .edges
                    .get(&idx)
                    .map(|edges| {
                        edges
                            .iter()
                            .filter(|e| graph.node_index.values().any(|&v| v == e.to))
                            .count()
                    })
                    .unwrap_or(0);
                out_count == 0
            })
            .collect();

        // Build a reverse adjacency for tracing back from leaves.
        let mut predecessors_map: std::collections::HashMap<usize, Vec<usize>> =
            std::collections::HashMap::new();
        for (&from, edges) in &graph.edges {
            if !graph.node_index.values().any(|&v| v == from) {
                continue;
            }
            for edge in edges {
                if graph.node_index.values().any(|&v| v == edge.to) {
                    predecessors_map.entry(edge.to).or_default().push(from);
                }
            }
        }

        for leaf in &leaves {
            // Skip reference nodes -- never clip reference path.
            if graph.nodes[*leaf].is_reference {
                continue;
            }

            // Trace back along single-predecessor chain.
            let mut chain = vec![*leaf];
            let mut current = *leaf;

            loop {
                let preds = predecessors_map.get(&current);
                let preds = match preds {
                    Some(p) => p,
                    None => break,
                };

                // Filter to active predecessors.
                let active_preds: Vec<usize> = preds
                    .iter()
                    .copied()
                    .filter(|p| graph.node_index.values().any(|&v| v == *p))
                    .collect();

                if active_preds.len() != 1 {
                    // Either zero or multiple predecessors -- stop tracing.
                    break;
                }

                let pred = active_preds[0];

                // Don't include the branching node itself in the chain.
                if pred == graph.source || pred == graph.sink {
                    break;
                }

                // Don't clip through reference nodes.
                if graph.nodes[pred].is_reference {
                    break;
                }

                chain.push(pred);
                current = pred;

                if chain.len() > max_tip_len {
                    break;
                }
            }

            if chain.len() > max_tip_len {
                continue;
            }

            // The chain must branch off a node with >= 2 successors.
            // The "branch point" is the predecessor of the last node in our chain.
            let tail = *chain.last().unwrap();
            let tail_preds = predecessors_map.get(&tail);
            let branch_point = match tail_preds {
                Some(preds) => {
                    let active: Vec<usize> = preds
                        .iter()
                        .copied()
                        .filter(|p| graph.node_index.values().any(|&v| v == *p))
                        .collect();
                    if active.len() == 1 {
                        active[0]
                    } else {
                        continue; // No single branch point -- skip.
                    }
                }
                None => continue,
            };

            // Check that the branch point has >= 2 active successors.
            let branch_successors = graph
                .edges
                .get(&branch_point)
                .map(|edges| {
                    edges
                        .iter()
                        .filter(|e| graph.node_index.values().any(|&v| v == e.to))
                        .count()
                })
                .unwrap_or(0);

            if branch_successors < 2 {
                continue;
            }

            // Remove the entire tip chain.
            for &idx in &chain {
                remove_node(graph, idx);
            }
            clipped_any = true;
        }

        if !clipped_any {
            break;
        }
    }
}

/// Collapse bubbles by removing low-coverage non-reference parallel paths.
///
/// A simplified implementation: for each non-reference node, if its count is
/// below `min_ratio * median_reference_coverage`, remove it. This catches
/// the most common bubble pattern where a sequencing error creates a short
/// parallel path with much lower coverage than the main path.
///
/// A more sophisticated implementation would detect true bubbles (diverge and
/// reconverge within k nodes), but this simplified version is effective for
/// the common case and avoids the complexity of full bubble detection.
pub fn collapse_bubbles(graph: &mut KmerGraph, min_ratio: f64) {
    // Compute median reference coverage.
    let mut ref_counts: Vec<u64> = graph
        .nodes
        .iter()
        .enumerate()
        .filter(|(idx, node)| {
            *idx != graph.source
                && *idx != graph.sink
                && node.is_reference
                && graph.node_index.contains_key(&node.kmer)
        })
        .map(|(_, node)| node.count)
        .collect();

    if ref_counts.is_empty() {
        return;
    }

    ref_counts.sort_unstable();
    let median_ref = ref_counts[ref_counts.len() / 2];

    let threshold = (median_ref as f64 * min_ratio).ceil() as u64;

    // Remove non-reference nodes below the threshold.
    let to_remove: Vec<usize> = graph
        .node_index
        .values()
        .copied()
        .filter(|&idx| {
            idx != graph.source
                && idx != graph.sink
                && !graph.nodes[idx].is_reference
                && graph.nodes[idx].count < threshold
        })
        .collect();

    for idx in to_remove {
        remove_node(graph, idx);
    }
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;

    use super::*;
    use crate::graph::{Edge, GraphNode, KmerGraph};

    /// Helper: build a graph from a description of nodes and edges.
    ///
    /// Nodes are specified as (kmer, count, is_reference).
    /// Edges are specified as (from_kmer, to_kmer, weight, is_reference).
    fn build_test_graph(
        node_specs: &[(&str, u64, bool)],
        edge_specs: &[(&str, &str, f64, bool)],
    ) -> KmerGraph {
        let mut nodes = Vec::new();
        let mut node_index = HashMap::new();

        // Always add BigBang (0) and BigCrunch (1).
        nodes.push(GraphNode {
            kmer: "BigBang".to_string(),
            count: 0,
            is_reference: true,
        });
        node_index.insert("BigBang".to_string(), 0);

        nodes.push(GraphNode {
            kmer: "BigCrunch".to_string(),
            count: 0,
            is_reference: true,
        });
        node_index.insert("BigCrunch".to_string(), 1);

        for &(kmer, count, is_ref) in node_specs {
            let idx = nodes.len();
            nodes.push(GraphNode {
                kmer: kmer.to_string(),
                count,
                is_reference: is_ref,
            });
            node_index.insert(kmer.to_string(), idx);
        }

        let mut edges: HashMap<usize, Vec<Edge>> = HashMap::new();
        for &(from, to, weight, is_ref) in edge_specs {
            let from_idx = *node_index.get(from).expect("from node not found");
            let to_idx = *node_index.get(to).expect("to node not found");
            edges.entry(from_idx).or_default().push(Edge {
                to: to_idx,
                weight,
                is_reference: is_ref,
            });
        }

        KmerGraph {
            nodes,
            edges,
            node_index,
            source: 0,
            sink: 1,
        }
    }

    // -----------------------------------------------------------------------
    // remove_dead_ends tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_remove_dead_ends_disconnected_nodes() {
        // Graph: BigBang -> A -> B -> BigCrunch, plus disconnected node C.
        let mut graph = build_test_graph(
            &[
                ("AAAA", 100, true),
                ("AAAB", 90, true),
                ("CCCC", 50, false), // disconnected
            ],
            &[
                ("BigBang", "AAAA", 0.01, true),
                ("AAAA", "AAAB", 0.01, true),
                ("AAAB", "BigCrunch", 0.01, true),
            ],
        );

        assert!(graph.node_index.contains_key("CCCC"));

        remove_dead_ends(&mut graph);

        // CCCC should be removed (not reachable from source AND to sink).
        assert!(!graph.node_index.contains_key("CCCC"));
        // Connected nodes should remain.
        assert!(graph.node_index.contains_key("AAAA"));
        assert!(graph.node_index.contains_key("AAAB"));
        assert!(graph.node_index.contains_key("BigBang"));
        assert!(graph.node_index.contains_key("BigCrunch"));
    }

    #[test]
    fn test_remove_dead_ends_forward_only() {
        // Node D is reachable from source but has no path to sink.
        let mut graph = build_test_graph(
            &[
                ("AAAA", 100, true),
                ("AAAB", 90, true),
                ("DDDD", 30, false), // forward-reachable only
            ],
            &[
                ("BigBang", "AAAA", 0.01, true),
                ("AAAA", "AAAB", 0.01, true),
                ("AAAB", "BigCrunch", 0.01, true),
                ("AAAA", "DDDD", 1.0, false), // DDDD reachable from source
                // but no edge from DDDD to anything that reaches sink
            ],
        );

        remove_dead_ends(&mut graph);

        assert!(!graph.node_index.contains_key("DDDD"));
        assert!(graph.node_index.contains_key("AAAA"));
        assert!(graph.node_index.contains_key("AAAB"));
    }

    #[test]
    fn test_remove_dead_ends_backward_only() {
        // Node D has an edge to sink but is not reachable from source.
        let mut graph = build_test_graph(
            &[
                ("AAAA", 100, true),
                ("AAAB", 90, true),
                ("DDDD", 30, false),
            ],
            &[
                ("BigBang", "AAAA", 0.01, true),
                ("AAAA", "AAAB", 0.01, true),
                ("AAAB", "BigCrunch", 0.01, true),
                ("DDDD", "BigCrunch", 1.0, false), // can reach sink
                // but not reachable from source
            ],
        );

        remove_dead_ends(&mut graph);

        assert!(!graph.node_index.contains_key("DDDD"));
    }

    #[test]
    fn test_remove_dead_ends_preserves_valid_alt_path() {
        // Alt path: BigBang -> A -> C -> BigCrunch (should survive).
        let mut graph = build_test_graph(
            &[
                ("AAAA", 100, true),
                ("AAAB", 90, true),
                ("AAAC", 50, false), // valid alt node on a source-to-sink path
            ],
            &[
                ("BigBang", "AAAA", 0.01, true),
                ("AAAA", "AAAB", 0.01, true),
                ("AAAB", "BigCrunch", 0.01, true),
                ("AAAA", "AAAC", 1.0, false),
                ("AAAC", "BigCrunch", 1.0, false),
            ],
        );

        remove_dead_ends(&mut graph);

        // AAAC is on a valid source-to-sink path -- must survive.
        assert!(graph.node_index.contains_key("AAAC"));
        assert!(graph.node_index.contains_key("AAAA"));
        assert!(graph.node_index.contains_key("AAAB"));
    }

    // -----------------------------------------------------------------------
    // remove_low_coverage tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_remove_low_coverage_preserves_reference() {
        let mut graph = build_test_graph(
            &[
                ("AAAA", 1, true),  // low count, but reference -- must survive
                ("AAAB", 90, true),
                ("CCCC", 1, false), // low count, non-reference -- should be removed
            ],
            &[
                ("BigBang", "AAAA", 0.01, true),
                ("AAAA", "AAAB", 0.01, true),
                ("AAAB", "BigCrunch", 0.01, true),
                ("AAAA", "CCCC", 1.0, false),
                ("CCCC", "BigCrunch", 1.0, false),
            ],
        );

        remove_low_coverage(&mut graph, 5);

        // Reference node with count=1 must survive.
        assert!(graph.node_index.contains_key("AAAA"));
        // Non-reference node with count=1 (< 5) should be removed.
        assert!(!graph.node_index.contains_key("CCCC"));
    }

    #[test]
    fn test_remove_low_coverage_keeps_above_threshold() {
        let mut graph = build_test_graph(
            &[
                ("AAAA", 100, true),
                ("AAAB", 90, true),
                ("CCCC", 50, false), // above threshold
                ("DDDD", 3, false),  // below threshold
            ],
            &[
                ("BigBang", "AAAA", 0.01, true),
                ("AAAA", "AAAB", 0.01, true),
                ("AAAB", "BigCrunch", 0.01, true),
                ("AAAA", "CCCC", 1.0, false),
                ("CCCC", "BigCrunch", 1.0, false),
                ("AAAA", "DDDD", 1.0, false),
                ("DDDD", "BigCrunch", 1.0, false),
            ],
        );

        remove_low_coverage(&mut graph, 10);

        assert!(graph.node_index.contains_key("CCCC")); // 50 >= 10
        assert!(!graph.node_index.contains_key("DDDD")); // 3 < 10
    }

    // -----------------------------------------------------------------------
    // compute_adaptive_min_count tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_adaptive_min_count_basic() {
        let graph = build_test_graph(
            &[
                ("AAAA", 1000, true),
                ("AAAB", 1000, true),
                ("AAAC", 1000, true),
            ],
            &[
                ("BigBang", "AAAA", 0.01, true),
                ("AAAA", "AAAB", 0.01, true),
                ("AAAB", "AAAC", 0.01, true),
                ("AAAC", "BigCrunch", 0.01, true),
            ],
        );

        // median = 1000, error_rate = 0.001 -> threshold = ceil(1.0) = 1
        // max(2, 1) = 2
        let min_count = compute_adaptive_min_count(&graph, 0.001);
        assert_eq!(min_count, 2);
    }

    #[test]
    fn test_adaptive_min_count_high_coverage() {
        let graph = build_test_graph(
            &[
                ("AAAA", 10000, true),
                ("AAAB", 10000, true),
                ("AAAC", 10000, true),
            ],
            &[
                ("BigBang", "AAAA", 0.01, true),
                ("AAAA", "AAAB", 0.01, true),
                ("AAAB", "AAAC", 0.01, true),
                ("AAAC", "BigCrunch", 0.01, true),
            ],
        );

        // median = 10000, error_rate = 0.001 -> threshold = ceil(10.0) = 10
        // max(2, 10) = 10
        let min_count = compute_adaptive_min_count(&graph, 0.001);
        assert_eq!(min_count, 10);
    }

    #[test]
    fn test_adaptive_min_count_no_ref_nodes() {
        let graph = build_test_graph(
            &[("CCCC", 50, false)],
            &[
                ("BigBang", "CCCC", 1.0, false),
                ("CCCC", "BigCrunch", 1.0, false),
            ],
        );

        // No reference nodes -> default to 2.
        let min_count = compute_adaptive_min_count(&graph, 0.001);
        assert_eq!(min_count, 2);
    }

    // -----------------------------------------------------------------------
    // clip_tips tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_clip_tips_removes_short_branch() {
        // Graph: BigBang -> A -> B -> BigCrunch
        //                   A -> C (dead end, tip of length 1)
        let mut graph = build_test_graph(
            &[
                ("AAAA", 100, true),
                ("AAAB", 90, true),
                ("CCCC", 5, false), // tip: dead-end branch from AAAA
            ],
            &[
                ("BigBang", "AAAA", 0.01, true),
                ("AAAA", "AAAB", 0.01, true),
                ("AAAB", "BigCrunch", 0.01, true),
                ("AAAA", "CCCC", 1.0, false),
            ],
        );

        // max_tip_len = 3 (should clip CCCC, which is length 1).
        clip_tips(&mut graph, 3);

        assert!(!graph.node_index.contains_key("CCCC"));
        // Main path should be intact.
        assert!(graph.node_index.contains_key("AAAA"));
        assert!(graph.node_index.contains_key("AAAB"));
    }

    #[test]
    fn test_clip_tips_preserves_long_branch() {
        // Graph: BigBang -> A -> B -> BigCrunch
        //                   A -> C -> D -> E (dead end, tip of length 3)
        let mut graph = build_test_graph(
            &[
                ("AAAA", 100, true),
                ("AAAB", 90, true),
                ("CCCC", 20, false),
                ("DDDD", 15, false),
                ("EEEE", 10, false),
            ],
            &[
                ("BigBang", "AAAA", 0.01, true),
                ("AAAA", "AAAB", 0.01, true),
                ("AAAB", "BigCrunch", 0.01, true),
                ("AAAA", "CCCC", 1.0, false),
                ("CCCC", "DDDD", 1.0, false),
                ("DDDD", "EEEE", 1.0, false),
            ],
        );

        // max_tip_len = 2 -- the tip is 3 nodes long, so should NOT be clipped.
        clip_tips(&mut graph, 2);

        assert!(graph.node_index.contains_key("CCCC"));
        assert!(graph.node_index.contains_key("DDDD"));
        assert!(graph.node_index.contains_key("EEEE"));
    }

    #[test]
    fn test_clip_tips_never_clips_reference() {
        // Graph: BigBang -> A -> B (dead end, but both are reference)
        let mut graph = build_test_graph(
            &[
                ("AAAA", 100, true),
                ("AAAB", 90, true), // reference dead-end
            ],
            &[
                ("BigBang", "AAAA", 0.01, true),
                ("AAAA", "AAAB", 0.01, true),
                // No edge from AAAB to BigCrunch -- it's a dead end.
            ],
        );

        clip_tips(&mut graph, 5);

        // Reference nodes should never be clipped.
        assert!(graph.node_index.contains_key("AAAA"));
        assert!(graph.node_index.contains_key("AAAB"));
    }

    // -----------------------------------------------------------------------
    // collapse_bubbles tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_collapse_bubbles_removes_low_coverage_alt() {
        // Reference path: BigBang -> A -> B -> BigCrunch (counts 1000, 1000).
        // Alt node C with count 5 (well below 0.10 * 1000 = 100).
        let mut graph = build_test_graph(
            &[
                ("AAAA", 1000, true),
                ("AAAB", 1000, true),
                ("CCCC", 5, false), // very low coverage alt
            ],
            &[
                ("BigBang", "AAAA", 0.01, true),
                ("AAAA", "AAAB", 0.01, true),
                ("AAAB", "BigCrunch", 0.01, true),
                ("AAAA", "CCCC", 1.0, false),
                ("CCCC", "BigCrunch", 1.0, false),
            ],
        );

        collapse_bubbles(&mut graph, 0.10);

        // CCCC (count=5) < ceil(1000 * 0.10) = 100 -> removed.
        assert!(!graph.node_index.contains_key("CCCC"));
    }

    #[test]
    fn test_collapse_bubbles_preserves_high_coverage_alt() {
        // Alt node C with count 200 (above 0.10 * 1000 = 100).
        let mut graph = build_test_graph(
            &[
                ("AAAA", 1000, true),
                ("AAAB", 1000, true),
                ("CCCC", 200, false), // high enough coverage
            ],
            &[
                ("BigBang", "AAAA", 0.01, true),
                ("AAAA", "AAAB", 0.01, true),
                ("AAAB", "BigCrunch", 0.01, true),
                ("AAAA", "CCCC", 1.0, false),
                ("CCCC", "BigCrunch", 1.0, false),
            ],
        );

        collapse_bubbles(&mut graph, 0.10);

        // CCCC (count=200) >= ceil(1000 * 0.10) = 100 -> preserved.
        assert!(graph.node_index.contains_key("CCCC"));
    }

    // -----------------------------------------------------------------------
    // Full prune_graph pipeline tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_prune_graph_full_pipeline() {
        // Build a messy graph with multiple kinds of noise:
        // - Reference path: A -> B (counts 1000)
        // - Disconnected node D
        // - Low-coverage alt node E (count=1)
        // - Valid alt node F (count=500)
        let mut graph = build_test_graph(
            &[
                ("AAAA", 1000, true),
                ("AAAB", 1000, true),
                ("DDDD", 50, false), // disconnected
                ("EEEE", 1, false),  // low coverage alt
                ("FFFF", 500, false), // valid alt
            ],
            &[
                ("BigBang", "AAAA", 0.01, true),
                ("AAAA", "AAAB", 0.01, true),
                ("AAAB", "BigCrunch", 0.01, true),
                ("AAAA", "EEEE", 1.0, false),
                ("EEEE", "BigCrunch", 1.0, false),
                ("AAAA", "FFFF", 1.0, false),
                ("FFFF", "BigCrunch", 1.0, false),
            ],
        );

        let config = PruneConfig {
            min_count: 0,     // use adaptive
            error_rate: 0.001, // -> threshold = max(2, ceil(1000*0.001)) = 2
            clip_tips: true,
            max_tip_fraction: 0.5,
            collapse_bubbles: true,
            bubble_coverage_ratio: 0.10, // -> threshold = ceil(1000*0.10) = 100
            ..Default::default()
        };

        let removed = prune_graph(&mut graph, &config, 4);

        // DDDD: disconnected -> removed by dead_ends.
        assert!(!graph.node_index.contains_key("DDDD"));
        // EEEE: count=1 < bubble threshold 100 -> removed.
        assert!(!graph.node_index.contains_key("EEEE"));
        // FFFF: count=500 >= 100 -> survives.
        assert!(graph.node_index.contains_key("FFFF"));
        // Reference path survives.
        assert!(graph.node_index.contains_key("AAAA"));
        assert!(graph.node_index.contains_key("AAAB"));

        assert!(removed >= 2); // At least DDDD and EEEE.
    }

    #[test]
    fn test_prune_graph_clean_graph_unchanged() {
        // A clean graph with only the reference path -- pruning should change nothing.
        let mut graph = build_test_graph(
            &[
                ("AAAA", 100, true),
                ("AAAB", 90, true),
                ("AAAC", 80, true),
            ],
            &[
                ("BigBang", "AAAA", 0.01, true),
                ("AAAA", "AAAB", 0.01, true),
                ("AAAB", "AAAC", 0.01, true),
                ("AAAC", "BigCrunch", 0.01, true),
            ],
        );

        let config = PruneConfig::default();
        let removed = prune_graph(&mut graph, &config, 4);

        assert_eq!(removed, 0);
        assert!(graph.node_index.contains_key("AAAA"));
        assert!(graph.node_index.contains_key("AAAB"));
        assert!(graph.node_index.contains_key("AAAC"));
    }

    #[test]
    fn test_prune_graph_snv_alt_path_survives() {
        // Simulate a real SNV: ref path A->B, alt path A->C, where C has
        // meaningful coverage. The alt path should survive pruning.
        let mut graph = build_test_graph(
            &[
                ("AAAA", 1000, true),
                ("AAAB", 1000, true),
                ("AAAC", 100, false), // SNV alt node -- meaningful coverage
            ],
            &[
                ("BigBang", "AAAA", 0.01, true),
                ("AAAA", "AAAB", 0.01, true),
                ("AAAB", "BigCrunch", 0.01, true),
                ("AAAA", "AAAC", 1.0, false),
                ("AAAC", "BigCrunch", 1.0, false),
            ],
        );

        let config = PruneConfig::default();
        let removed = prune_graph(&mut graph, &config, 31);

        assert_eq!(removed, 0);
        assert!(graph.node_index.contains_key("AAAC"));
    }

    #[test]
    fn test_prune_graph_with_explicit_min_count() {
        // Use an explicit min_count instead of adaptive.
        let mut graph = build_test_graph(
            &[
                ("AAAA", 100, true),
                ("AAAB", 90, true),
                ("CCCC", 8, false),  // below min_count=10
                ("DDDD", 15, false), // above min_count=10
            ],
            &[
                ("BigBang", "AAAA", 0.01, true),
                ("AAAA", "AAAB", 0.01, true),
                ("AAAB", "BigCrunch", 0.01, true),
                ("AAAA", "CCCC", 1.0, false),
                ("CCCC", "BigCrunch", 1.0, false),
                ("AAAA", "DDDD", 1.0, false),
                ("DDDD", "BigCrunch", 1.0, false),
            ],
        );

        let config = PruneConfig {
            min_count: 10,
            clip_tips: false,
            collapse_bubbles: false,
            ..Default::default()
        };

        let removed = prune_graph(&mut graph, &config, 4);

        assert!(!graph.node_index.contains_key("CCCC")); // 8 < 10
        assert!(graph.node_index.contains_key("DDDD"));  // 15 >= 10
        assert_eq!(removed, 1);
    }

    #[test]
    fn test_prune_empty_graph() {
        // Graph with only source and sink -- nothing to prune.
        let mut graph = build_test_graph(&[], &[]);

        let config = PruneConfig::default();
        let removed = prune_graph(&mut graph, &config, 31);

        assert_eq!(removed, 0);
    }

    #[test]
    fn test_remove_dead_ends_keeps_source_and_sink() {
        // Even if source/sink have no edges, they should never be removed.
        let mut graph = build_test_graph(&[], &[]);

        remove_dead_ends(&mut graph);

        assert!(graph.node_index.contains_key("BigBang"));
        assert!(graph.node_index.contains_key("BigCrunch"));
    }
}
