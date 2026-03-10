/// Graph construction from walked k-mers.
///
/// Builds a directed weighted graph where:
/// - Each node is a discovered k-mer
/// - Edges connect k-mers that overlap by (k-1) bases
/// - Reference edges get weight 0.01, non-reference edges get weight 1.0
/// - Virtual source ("BigBang") and sink ("BigCrunch") nodes are added

use std::collections::HashMap;

use super::{Edge, GraphNode, KmerGraph};
use crate::walker::WalkResult;

/// Reference edge weight (low to prefer reference path in Dijkstra).
pub const REF_EDGE_WEIGHT: f64 = 0.01;
/// Non-reference edge weight.
pub const ALT_EDGE_WEIGHT: f64 = 1.0;

/// Build a KmerGraph from walk results and reference k-mer list.
///
/// The `ref_kmers` slice defines the ordered reference k-mers. The first
/// element connects to the virtual source and the last connects to the
/// virtual sink. All k-mers in `walk.nodes` are added as graph nodes,
/// and edges are created between any pair where the last (k-1) bases of
/// one k-mer equal the first (k-1) bases of another.
pub fn build_graph(walk: &WalkResult, ref_kmers: &[String]) -> KmerGraph {
    let mut nodes: Vec<GraphNode> = Vec::new();
    let mut edges: HashMap<usize, Vec<Edge>> = HashMap::new();
    let mut node_index: HashMap<String, usize> = HashMap::new();

    // 1. Add virtual source node "BigBang" at index 0
    nodes.push(GraphNode {
        kmer: "BigBang".to_string(),
        count: 0,
        is_reference: true,
    });
    node_index.insert("BigBang".to_string(), 0);
    let source = 0;

    // 2. Add virtual sink node "BigCrunch" at index 1
    nodes.push(GraphNode {
        kmer: "BigCrunch".to_string(),
        count: 0,
        is_reference: true,
    });
    node_index.insert("BigCrunch".to_string(), 1);
    let sink = 1;

    // 3. Add all walked k-mers as graph nodes
    for (kmer, &count) in &walk.nodes {
        let is_ref = walk.reference_kmers.contains(kmer);
        let idx = nodes.len();
        nodes.push(GraphNode {
            kmer: kmer.clone(),
            count,
            is_reference: is_ref,
        });
        node_index.insert(kmer.clone(), idx);
    }

    // 4. Build edges using suffix→index map for O(n) lookup
    //    For each k-mer, its suffix is kmer[1..] (last k-1 bases).
    //    For each k-mer, its prefix is kmer[..k-1] (first k-1 bases).
    //    An edge a→b exists when suffix(a) == prefix(b).
    //
    //    Build a map: suffix → list of node indices that have that suffix.
    //    Then for each node, look up its prefix in the suffix map.

    // Determine k from the first real k-mer (skip virtual nodes)
    let k = if let Some(kmer) = walk.nodes.keys().next() {
        kmer.len()
    } else {
        // No k-mers to process, return empty graph
        return KmerGraph {
            nodes,
            edges,
            node_index,
            source,
            sink,
        };
    };

    // Build suffix map: suffix string → Vec<node index>
    // A suffix is the last (k-1) characters of a k-mer.
    let mut suffix_map: HashMap<&str, Vec<usize>> = HashMap::new();
    for (kmer, &idx) in &node_index {
        // Skip virtual nodes
        if kmer == "BigBang" || kmer == "BigCrunch" {
            continue;
        }
        if kmer.len() >= k {
            let suffix = &kmer[1..]; // last k-1 bases
            suffix_map.entry(suffix).or_default().push(idx);
        }
    }

    // For each real k-mer node, find edges by looking up its prefix in the suffix map
    for (kmer, &from_idx) in &node_index {
        if kmer == "BigBang" || kmer == "BigCrunch" {
            continue;
        }
        if kmer.len() < k {
            continue;
        }

        let prefix = &kmer[..k - 1]; // first k-1 bases

        // Look up which nodes have this prefix as their suffix
        if let Some(sources) = suffix_map.get(prefix) {
            for &src_idx in sources {
                if src_idx == from_idx {
                    // Skip self-loops (shouldn't happen with distinct k-mers, but safe)
                    continue;
                }
                // Edge goes from src_idx → from_idx
                // (src's suffix matches from's prefix, so src precedes from)
                let src_is_ref = nodes[src_idx].is_reference;
                let dst_is_ref = nodes[from_idx].is_reference;
                let is_ref_edge = src_is_ref && dst_is_ref;
                let weight = if is_ref_edge {
                    REF_EDGE_WEIGHT
                } else {
                    ALT_EDGE_WEIGHT
                };

                edges.entry(src_idx).or_default().push(Edge {
                    to: from_idx,
                    weight,
                    is_reference: is_ref_edge,
                });
            }
        }
    }

    // 5. Connect BigBang to the first reference k-mer
    if let Some(first_ref) = ref_kmers.first() {
        if let Some(&first_idx) = node_index.get(first_ref) {
            edges.entry(source).or_default().push(Edge {
                to: first_idx,
                weight: REF_EDGE_WEIGHT,
                is_reference: true,
            });
        }
    }

    // 6. Connect the last reference k-mer to BigCrunch
    if let Some(last_ref) = ref_kmers.last() {
        if let Some(&last_idx) = node_index.get(last_ref) {
            edges.entry(last_idx).or_default().push(Edge {
                to: sink,
                weight: REF_EDGE_WEIGHT,
                is_reference: true,
            });
        }
    }

    // 7. Connect BigBang to any non-reference k-mer whose prefix matches
    //    the first reference k-mer's prefix (alternate starts)
    if let Some(first_ref) = ref_kmers.first() {
        if first_ref.len() >= k {
            let ref_prefix = &first_ref[..k - 1];
            for (kmer, &idx) in &node_index {
                if kmer == "BigBang" || kmer == "BigCrunch" {
                    continue;
                }
                if kmer.len() < k {
                    continue;
                }
                // Skip the first reference k-mer itself (already connected)
                if kmer == first_ref {
                    continue;
                }
                // Check if this k-mer's prefix matches the first ref k-mer's prefix
                if &kmer[..k - 1] == ref_prefix {
                    let is_ref_node = nodes[idx].is_reference;
                    // This is a non-standard connection, so it's alt weight
                    // unless the node is also reference
                    let is_ref_edge = is_ref_node;
                    let weight = if is_ref_edge {
                        REF_EDGE_WEIGHT
                    } else {
                        ALT_EDGE_WEIGHT
                    };
                    edges.entry(source).or_default().push(Edge {
                        to: idx,
                        weight,
                        is_reference: is_ref_edge,
                    });
                }
            }
        }
    }

    // 8. Also connect any non-reference k-mer that could reach the sink:
    //    k-mers whose suffix matches the last reference k-mer's suffix.
    if let Some(last_ref) = ref_kmers.last() {
        if last_ref.len() >= k {
            let ref_suffix = &last_ref[1..]; // last k-1 bases of last ref k-mer
            for (kmer, &idx) in &node_index {
                if kmer == "BigBang" || kmer == "BigCrunch" {
                    continue;
                }
                if kmer.len() < k {
                    continue;
                }
                // Skip the last reference k-mer itself (already connected)
                if kmer == last_ref {
                    continue;
                }
                // Check if this k-mer's suffix matches the last ref k-mer's suffix
                if &kmer[1..] == ref_suffix {
                    let is_ref_node = nodes[idx].is_reference;
                    let is_ref_edge = is_ref_node;
                    let weight = if is_ref_edge {
                        REF_EDGE_WEIGHT
                    } else {
                        ALT_EDGE_WEIGHT
                    };
                    edges.entry(idx).or_default().push(Edge {
                        to: sink,
                        weight,
                        is_reference: is_ref_edge,
                    });
                }
            }
        }
    }

    KmerGraph {
        nodes,
        edges,
        node_index,
        source,
        sink,
    }
}
