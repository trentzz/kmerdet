/// Graph construction from walked k-mers.
///
/// Builds a directed weighted graph where:
/// - Each node is a discovered k-mer
/// - Edges connect k-mers that overlap by (k-1) bases
/// - Reference edges get weight 0.01, non-reference edges get weight 1.0
/// - Virtual source ("BigBang") and sink ("BigCrunch") nodes are added

use super::KmerGraph;
use crate::walker::WalkResult;

/// Reference edge weight (low to prefer reference path in Dijkstra).
pub const REF_EDGE_WEIGHT: f64 = 0.01;
/// Non-reference edge weight.
pub const ALT_EDGE_WEIGHT: f64 = 1.0;

/// Build a KmerGraph from walk results and reference k-mer list.
pub fn build_graph(_walk: &WalkResult, _ref_kmers: &[String]) -> KmerGraph {
    // TODO Phase 3: Build graph from walked k-mers
    todo!("graph builder not yet implemented")
}
