pub mod extension;

use std::collections::{HashMap, HashSet};

use crate::jellyfish::KmerDatabase;

/// Configuration for the DFS k-mer walker.
#[derive(Debug, Clone)]
pub struct WalkerConfig {
    /// Min absolute k-mer count for extension.
    pub count: u32,
    /// Min ratio of child count to total sibling counts.
    pub ratio: f64,
    /// Max DFS stack depth.
    pub max_stack: usize,
    /// Max branching points.
    pub max_break: usize,
    /// Max total discovered nodes.
    pub max_node: usize,
}

impl Default for WalkerConfig {
    fn default() -> Self {
        Self {
            count: 2,
            ratio: 0.05,
            max_stack: 500,
            max_break: 10,
            max_node: 10000,
        }
    }
}

/// Result of walking k-mers from a reference sequence.
#[derive(Debug)]
pub struct WalkResult {
    /// All discovered k-mer strings with their counts.
    pub nodes: HashMap<String, u64>,
    /// Set of k-mers that are part of the reference path.
    pub reference_kmers: HashSet<String>,
}

/// Walk k-mers starting from a reference sequence using iterative DFS.
pub fn walk(
    _db: &dyn KmerDatabase,
    _ref_kmers: &[String],
    _config: &WalkerConfig,
) -> WalkResult {
    // TODO Phase 2: Implement iterative DFS walking
    todo!("walker not yet implemented")
}
