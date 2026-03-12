pub mod adaptive;
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
    /// Use adaptive thresholds based on sample coverage and error rate.
    pub adaptive: bool,
}

impl Default for WalkerConfig {
    fn default() -> Self {
        Self {
            count: 2,
            ratio: 0.05,
            max_stack: 500,
            max_break: 10,
            max_node: 10000,
            adaptive: false,
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
///
/// For each reference k-mer, performs forward DFS extension to discover
/// all reachable k-mers above the count/ratio thresholds. Branching points
/// increment a break counter, and exploration is bounded by max_stack,
/// max_break, and max_node limits.
pub fn walk(
    db: &dyn KmerDatabase,
    ref_kmers: &[String],
    config: &WalkerConfig,
) -> WalkResult {
    let mut nodes: HashMap<String, u64> = HashMap::new();
    let reference_kmers: HashSet<String> = ref_kmers.iter().cloned().collect();

    // Step 1: Add all reference k-mers to nodes with their counts from db
    for kmer in ref_kmers {
        let count = db.query(kmer);
        nodes.insert(kmer.clone(), count);
    }

    // Step 2: For each reference k-mer, perform iterative DFS forward extension
    for kmer in ref_kmers {
        // Stack holds (kmer_sequence, break_count)
        let mut stack: Vec<(String, usize)> = Vec::new();
        stack.push((kmer.clone(), 0));

        while let Some((current_kmer, breaks)) = stack.pop() {
            if nodes.len() >= config.max_node {
                break;
            }

            let children = extension::extend_forward(
                db,
                &current_kmer,
                config.count,
                config.ratio,
            );

            let num_children = children.len();

            for child in children {
                // Skip already-discovered k-mers
                if nodes.contains_key(&child.sequence) {
                    continue;
                }

                // Insert the new k-mer into discovered nodes
                nodes.insert(child.sequence.clone(), child.count);

                if nodes.len() >= config.max_node {
                    break;
                }

                // A break is counted when we're at a branching point (>1 child)
                let new_breaks = if num_children > 1 {
                    breaks + 1
                } else {
                    breaks
                };

                // Respect break limit
                if new_breaks > config.max_break {
                    continue;
                }

                // Respect stack depth limit
                if stack.len() >= config.max_stack {
                    continue;
                }

                stack.push((child.sequence.clone(), new_breaks));
            }
        }
    }

    WalkResult {
        nodes,
        reference_kmers,
    }
}
