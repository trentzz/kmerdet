pub mod builder;
pub mod pathfind;

use std::collections::HashMap;

/// A directed weighted k-mer graph.
///
/// Custom implementation (not petgraph) for specialized operations:
/// reference-edge removal, weighted Dijkstra, all-shortest-paths.
#[derive(Debug)]
pub struct KmerGraph {
    /// Node index → node data.
    pub nodes: Vec<GraphNode>,
    /// Adjacency list: node index → [(neighbor index, edge weight)].
    pub edges: HashMap<usize, Vec<Edge>>,
    /// K-mer string → node index lookup.
    pub node_index: HashMap<String, usize>,
    /// Index of the virtual source node ("BigBang").
    pub source: usize,
    /// Index of the virtual sink node ("BigCrunch").
    pub sink: usize,
}

#[derive(Debug, Clone)]
pub struct GraphNode {
    /// K-mer sequence (or virtual node name).
    pub kmer: String,
    /// K-mer count from the database (0 for virtual nodes).
    pub count: u64,
    /// Whether this k-mer is part of the reference sequence.
    pub is_reference: bool,
}

#[derive(Debug, Clone)]
pub struct Edge {
    /// Target node index.
    pub to: usize,
    /// Edge weight (0.01 for reference edges, 1.0 otherwise).
    pub weight: f64,
    /// Whether this is a reference edge.
    pub is_reference: bool,
}
