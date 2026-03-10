/// Pathfinding algorithms on the k-mer graph.
///
/// 1. Dijkstra from source and sink to find reachable nodes
/// 2. Enumerate alternative paths through non-reference edges
/// 3. Deduplicate by DNA sequence

use std::cmp::Ordering;
use std::collections::{BinaryHeap, HashMap, HashSet};

use super::KmerGraph;
use crate::sequence::path::KmerPath;

/// A state for Dijkstra's algorithm, ordered by distance (min-heap).
#[derive(Debug)]
struct State {
    cost: f64,
    node: usize,
}

impl PartialEq for State {
    fn eq(&self, other: &Self) -> bool {
        self.cost.to_bits() == other.cost.to_bits() && self.node == other.node
    }
}

impl Eq for State {}

impl Ord for State {
    fn cmp(&self, other: &Self) -> Ordering {
        // Reverse ordering for min-heap (BinaryHeap is a max-heap by default)
        other
            .cost
            .partial_cmp(&self.cost)
            .unwrap_or(Ordering::Equal)
            .then_with(|| self.node.cmp(&other.node))
    }
}

impl PartialOrd for State {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// Result of a Dijkstra run: distances and predecessors.
struct DijkstraResult {
    dist: HashMap<usize, f64>,
    pred: HashMap<usize, usize>,
}

/// Run Dijkstra's algorithm from a source node in the forward direction.
fn dijkstra_forward(graph: &KmerGraph, start: usize) -> DijkstraResult {
    let mut dist: HashMap<usize, f64> = HashMap::new();
    let mut pred: HashMap<usize, usize> = HashMap::new();
    let mut heap = BinaryHeap::new();

    dist.insert(start, 0.0);
    heap.push(State {
        cost: 0.0,
        node: start,
    });

    while let Some(State { cost, node }) = heap.pop() {
        // If we already found a shorter path, skip
        if let Some(&d) = dist.get(&node) {
            if cost > d {
                continue;
            }
        }

        // Explore neighbors
        if let Some(edges) = graph.edges.get(&node) {
            for edge in edges {
                let next_cost = cost + edge.weight;
                let is_shorter = match dist.get(&edge.to) {
                    Some(&d) => next_cost < d,
                    None => true,
                };
                if is_shorter {
                    dist.insert(edge.to, next_cost);
                    pred.insert(edge.to, node);
                    heap.push(State {
                        cost: next_cost,
                        node: edge.to,
                    });
                }
            }
        }
    }

    DijkstraResult { dist, pred }
}

/// Run Dijkstra's algorithm backward from a target node.
///
/// This reverses all edges: for each edge a→b in the graph, we consider
/// b→a. This finds shortest paths from any node TO the target.
fn dijkstra_backward(graph: &KmerGraph, target: usize) -> DijkstraResult {
    // Build reverse adjacency list
    let mut reverse_edges: HashMap<usize, Vec<(usize, f64)>> = HashMap::new();
    for (&from, edges) in &graph.edges {
        for edge in edges {
            reverse_edges
                .entry(edge.to)
                .or_default()
                .push((from, edge.weight));
        }
    }

    let mut dist: HashMap<usize, f64> = HashMap::new();
    let mut pred: HashMap<usize, usize> = HashMap::new();
    let mut heap = BinaryHeap::new();

    dist.insert(target, 0.0);
    heap.push(State {
        cost: 0.0,
        node: target,
    });

    while let Some(State { cost, node }) = heap.pop() {
        if let Some(&d) = dist.get(&node) {
            if cost > d {
                continue;
            }
        }

        if let Some(neighbors) = reverse_edges.get(&node) {
            for &(neighbor, weight) in neighbors {
                let next_cost = cost + weight;
                let is_shorter = match dist.get(&neighbor) {
                    Some(&d) => next_cost < d,
                    None => true,
                };
                if is_shorter {
                    dist.insert(neighbor, next_cost);
                    pred.insert(neighbor, node);
                    heap.push(State {
                        cost: next_cost,
                        node: neighbor,
                    });
                }
            }
        }
    }

    DijkstraResult { dist, pred }
}

/// Trace the path from source to a given node using the forward predecessor map.
fn trace_forward_path(pred: &HashMap<usize, usize>, source: usize, target: usize) -> Vec<usize> {
    let mut path = Vec::new();
    let mut current = target;
    while current != source {
        path.push(current);
        match pred.get(&current) {
            Some(&p) => current = p,
            None => return Vec::new(), // unreachable from source
        }
    }
    path.push(source);
    path.reverse();
    path
}

/// Trace the path from a given node to the sink using the backward predecessor map.
///
/// The backward predecessor map tells us: for each node, which node comes AFTER it
/// on the path to the sink (since Dijkstra was run on reversed edges).
fn trace_backward_path(
    pred: &HashMap<usize, usize>,
    sink: usize,
    start: usize,
) -> Vec<usize> {
    let mut path = Vec::new();
    let mut current = start;
    while current != sink {
        path.push(current);
        match pred.get(&current) {
            Some(&p) => current = p,
            None => return Vec::new(), // can't reach sink
        }
    }
    path.push(sink);
    path
}

/// Convert a path of node indices to a KmerPath, filtering out virtual nodes.
fn indices_to_kmer_path(graph: &KmerGraph, indices: &[usize], is_reference: bool) -> KmerPath {
    let kmers: Vec<String> = indices
        .iter()
        .filter(|&&idx| idx != graph.source && idx != graph.sink)
        .map(|&idx| graph.nodes[idx].kmer.clone())
        .collect();

    KmerPath {
        kmers,
        is_reference,
    }
}

/// Find all alternative paths (non-reference shortest paths) through the graph.
///
/// Algorithm:
/// 1. Run Dijkstra forward from source, backward from sink
/// 2. Find the reference path (lowest weight, all-reference edges)
/// 3. For each non-reference edge (u, v), if u is reachable from source
///    and v can reach the sink, construct the path source→...→u→v→...→sink
/// 4. Deduplicate paths by their DNA sequence
/// 5. Return reference path first, then alternative paths
pub fn find_alternative_paths(graph: &KmerGraph) -> Vec<KmerPath> {
    let source = graph.source;
    let sink = graph.sink;

    // 1. Run Dijkstra forward from source and backward from sink
    let forward = dijkstra_forward(graph, source);
    let backward = dijkstra_backward(graph, sink);

    // 2. Find the reference path by tracing forward from source to sink
    let ref_indices = trace_forward_path(&forward.pred, source, sink);
    let ref_path = if ref_indices.len() >= 2
        && *ref_indices.first().unwrap() == source
        && *ref_indices.last().unwrap() == sink
    {
        indices_to_kmer_path(graph, &ref_indices, true)
    } else {
        // No path from source to sink; return empty
        return Vec::new();
    };

    let ref_sequence = ref_path.to_sequence();

    // 3. Find alternative paths through non-reference edges
    let mut seen_sequences: HashSet<String> = HashSet::new();
    seen_sequences.insert(ref_sequence);

    let mut alt_paths: Vec<KmerPath> = Vec::new();

    // Iterate over all edges, looking for non-reference edges
    for (&from_idx, edges) in &graph.edges {
        for edge in edges {
            if edge.is_reference {
                continue; // Skip reference edges
            }

            let u = from_idx;
            let v = edge.to;

            // Check u is reachable from source and v can reach sink
            if !forward.dist.contains_key(&u) || !backward.dist.contains_key(&v) {
                continue;
            }

            // Build path: source→...→u then u→v then v→...→sink
            let path_to_u = trace_forward_path(&forward.pred, source, u);
            if path_to_u.is_empty() || *path_to_u.first().unwrap() != source {
                continue;
            }

            let path_from_v = trace_backward_path(&backward.pred, sink, v);
            if path_from_v.is_empty() || *path_from_v.last().unwrap() != sink {
                continue;
            }

            // Combine: path_to_u includes u, path_from_v includes v
            // So the full path is path_to_u + path_from_v (path_to_u ends at u,
            // and path_from_v starts at v)
            let mut full_path = path_to_u;
            full_path.extend_from_slice(&path_from_v);

            let kmer_path = indices_to_kmer_path(graph, &full_path, false);

            if kmer_path.kmers.is_empty() {
                continue;
            }

            // Deduplicate by sequence
            let sequence = kmer_path.to_sequence();
            if seen_sequences.contains(&sequence) {
                continue;
            }
            seen_sequences.insert(sequence);
            alt_paths.push(kmer_path);
        }
    }

    // 4. Return reference path first, then alternatives
    let mut result = vec![ref_path];
    result.extend(alt_paths);
    result
}
