/// Pathfinding algorithms on the k-mer graph.
///
/// 1. Dijkstra from source and sink to find reachable nodes
/// 2. Remove reference edges
/// 3. Enumerate all shortest paths through non-reference edges

use super::KmerGraph;
use crate::sequence::path::KmerPath;

/// Find all alternative paths (non-reference shortest paths) through the graph.
///
/// Algorithm:
/// 1. Run Dijkstra forward from source, backward from sink
/// 2. Identify reachable subgraph
/// 3. Remove reference edges
/// 4. Find all shortest paths from source to sink in the remaining graph
pub fn find_alternative_paths(_graph: &KmerGraph) -> Vec<KmerPath> {
    // TODO Phase 3: Implement Dijkstra + ref-edge removal + all-shortest-paths
    todo!("pathfinding not yet implemented")
}
