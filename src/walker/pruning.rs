/// Branch pruning for walk results.
///
/// After walking, the discovered k-mer set may contain short dead-end branches
/// (tips) caused by sequencing errors, and low-coverage bubbles where an error
/// path diverges from the true path and reconverges within a few k-mers.
/// Pruning removes these artifacts to produce cleaner graphs.

use std::collections::{HashMap, HashSet};

use super::WalkResult;

/// Remove tip branches: dead-end paths shorter than `k/2` extensions.
///
/// A "tip" is a chain of non-reference k-mers where each has exactly one
/// neighbor within the walk result (i.e., it is a dead-end chain). If the
/// chain is shorter than half the k-mer length and branches off a node with
/// multiple neighbors, it is likely a sequencing error artifact and is removed.
///
/// This iterates until no more tips are removed (fixed-point).
pub fn prune_tips(result: &mut WalkResult, k: u8) {
    let max_tip_len = (k as usize) / 2;
    if max_tip_len == 0 {
        return;
    }

    loop {
        let tips = find_tips(result, max_tip_len);
        if tips.is_empty() {
            break;
        }
        for kmer in &tips {
            result.nodes.remove(kmer);
        }
    }
}

/// Find all tip k-mers (dead-end chains of non-reference k-mers shorter than
/// `max_tip_len`).
fn find_tips(result: &WalkResult, max_tip_len: usize) -> HashSet<String> {
    // Build adjacency: for each k-mer in the walk result, find which other
    // k-mers in the result overlap by (k-1) bases (forward neighbor).
    let k = if let Some(kmer) = result.nodes.keys().next() {
        kmer.len()
    } else {
        return HashSet::new();
    };

    if k < 2 {
        return HashSet::new();
    }

    // Build suffix-to-kmer and prefix-to-kmer maps for fast neighbor lookup
    let mut suffix_map: HashMap<&str, Vec<&str>> = HashMap::new();
    let mut prefix_map: HashMap<&str, Vec<&str>> = HashMap::new();

    for kmer in result.nodes.keys() {
        let suffix = &kmer[1..];     // last k-1 bases
        let prefix = &kmer[..k - 1]; // first k-1 bases
        suffix_map.entry(suffix).or_default().push(kmer);
        prefix_map.entry(prefix).or_default().push(kmer);
    }

    // A k-mer's forward neighbors are those whose prefix matches this k-mer's suffix
    // A k-mer's backward neighbors are those whose suffix matches this k-mer's prefix
    let forward_neighbors = |kmer: &str| -> Vec<&str> {
        let suffix = &kmer[1..];
        prefix_map.get(suffix).cloned().unwrap_or_default()
    };

    let backward_neighbors = |kmer: &str| -> Vec<&str> {
        let prefix = &kmer[..k - 1];
        suffix_map.get(prefix).cloned().unwrap_or_default()
    };

    // Find leaf nodes: non-reference k-mers with no forward neighbors in the
    // walk result (dead ends).
    let mut leaves: Vec<String> = Vec::new();
    for kmer in result.nodes.keys() {
        if result.reference_kmers.contains(kmer) {
            continue;
        }
        let fwd = forward_neighbors(kmer);
        // A leaf has no forward neighbors OTHER than itself
        let has_forward = fwd.iter().any(|&n| n != kmer);
        if !has_forward {
            leaves.push(kmer.clone());
        }
    }

    // Also check backward leaves (no backward neighbors except self)
    for kmer in result.nodes.keys() {
        if result.reference_kmers.contains(kmer) {
            continue;
        }
        let bwd = backward_neighbors(kmer);
        let has_backward = bwd.iter().any(|&n| n != kmer);
        if !has_backward && !leaves.contains(kmer) {
            leaves.push(kmer.clone());
        }
    }

    let mut tips_to_remove: HashSet<String> = HashSet::new();

    for leaf in &leaves {
        // Trace back from the leaf, following the single-neighbor chain
        let mut chain: Vec<String> = vec![leaf.clone()];
        let mut current = leaf.clone();
        let is_tip = true;

        for _ in 0..max_tip_len {
            // Find the backward neighbor(s) of current that are in the walk result
            let back = backward_neighbors(&current);
            let back: Vec<&&str> = back
                .iter()
                .filter(|&&n| n != current && !tips_to_remove.contains(n))
                .collect();

            if back.is_empty() {
                // Reached a dead end backward too -- this is a disconnected fragment
                break;
            }

            if back.len() == 1 {
                let prev = back[0].to_string();
                if result.reference_kmers.contains(&prev) {
                    // Chain ends at a reference k-mer -- this is a tip
                    break;
                }
                // Check if prev has multiple forward neighbors (i.e., it's a branch point)
                let fwd_of_prev = forward_neighbors(&prev);
                let fwd_count = fwd_of_prev
                    .iter()
                    .filter(|&&n| n != prev)
                    .count();
                if fwd_count > 1 {
                    // prev is a branch point; chain from prev to leaf is a tip
                    break;
                }
                chain.push(prev.clone());
                current = prev;
            } else {
                // Multiple backward neighbors: current is at a junction
                break;
            }
        }

        // Only prune if the chain is shorter than max_tip_len
        // and none of the chain k-mers are reference k-mers
        if chain.len() <= max_tip_len
            && is_tip
            && chain
                .iter()
                .all(|km| !result.reference_kmers.contains(km))
        {
            for km in chain {
                tips_to_remove.insert(km);
            }
        }
    }

    tips_to_remove
}

/// Remove bubbles: parallel short paths with very low coverage relative to the
/// main path.
///
/// A simplified bubble pruning strategy: remove any non-reference k-mer whose
/// count is below a fraction (`min_bubble_ratio`) of the median reference
/// coverage. This catches low-coverage error paths that run parallel to the
/// reference path and reconverge.
///
/// The `min_bubble_ratio` default is 0.01 (1% of median reference coverage).
pub fn prune_bubbles(result: &mut WalkResult, min_bubble_ratio: f64) {
    // Compute median reference k-mer coverage
    let mut ref_counts: Vec<u64> = result
        .reference_kmers
        .iter()
        .filter_map(|k| result.nodes.get(k).copied())
        .collect();

    if ref_counts.is_empty() {
        return;
    }

    ref_counts.sort_unstable();
    let median = if ref_counts.len() % 2 == 0 {
        let mid = ref_counts.len() / 2;
        (ref_counts[mid - 1] + ref_counts[mid]) as f64 / 2.0
    } else {
        ref_counts[ref_counts.len() / 2] as f64
    };

    let threshold = (median * min_bubble_ratio) as u64;

    // Remove non-reference k-mers below the threshold
    let to_remove: Vec<String> = result
        .nodes
        .iter()
        .filter(|(kmer, &count)| {
            !result.reference_kmers.contains(*kmer) && count < threshold
        })
        .map(|(kmer, _)| kmer.clone())
        .collect();

    for kmer in to_remove {
        result.nodes.remove(&kmer);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::{HashMap, HashSet};

    fn make_walk_result(
        nodes: Vec<(&str, u64)>,
        ref_kmers: Vec<&str>,
    ) -> WalkResult {
        WalkResult {
            nodes: nodes
                .into_iter()
                .map(|(k, c)| (k.to_string(), c))
                .collect(),
            reference_kmers: ref_kmers.into_iter().map(|s| s.to_string()).collect(),
        }
    }

    #[test]
    fn test_prune_tips_removes_short_dead_end() {
        // Reference path: AACC -> ACCG -> CCGG -> CGGT
        // Tip branch off ACCG: CCGA (dead end, length 1)
        // With k=4, max_tip_len = 2, so this 1-kmer tip should be removed.
        let mut result = make_walk_result(
            vec![
                ("AACC", 100),
                ("ACCG", 100),
                ("CCGG", 100),
                ("CGGT", 100),
                ("CCGA", 5), // tip: forward-extends from ACCG but dead-ends
            ],
            vec!["AACC", "ACCG", "CCGG", "CGGT"],
        );

        prune_tips(&mut result, 4);

        assert!(
            !result.nodes.contains_key("CCGA"),
            "Tip k-mer CCGA should be pruned"
        );
        // Reference k-mers should remain
        assert!(result.nodes.contains_key("AACC"));
        assert!(result.nodes.contains_key("ACCG"));
        assert!(result.nodes.contains_key("CCGG"));
        assert!(result.nodes.contains_key("CGGT"));
    }

    #[test]
    fn test_prune_tips_keeps_long_branches() {
        // Reference path: AACC -> ACCG -> CCGG
        // Long branch: CCGA -> CGAT -> GATC (length 3, with k=4 max_tip_len=2)
        // This branch is longer than max_tip_len, so should NOT be pruned.
        let mut result = make_walk_result(
            vec![
                ("AACC", 100),
                ("ACCG", 100),
                ("CCGG", 100),
                ("CCGA", 50),
                ("CGAT", 50),
                ("GATC", 50),
            ],
            vec!["AACC", "ACCG", "CCGG"],
        );

        prune_tips(&mut result, 4);

        // The 3-kmer branch should survive (3 > max_tip_len=2)
        assert!(
            result.nodes.contains_key("CCGA"),
            "Long branch k-mer CCGA should survive"
        );
        assert!(
            result.nodes.contains_key("CGAT"),
            "Long branch k-mer CGAT should survive"
        );
        assert!(
            result.nodes.contains_key("GATC"),
            "Long branch k-mer GATC should survive"
        );
    }

    #[test]
    fn test_prune_tips_preserves_reference_kmers() {
        // Even if a reference k-mer looks like a dead end, it should not be pruned.
        let mut result = make_walk_result(
            vec![
                ("AACC", 100),
                ("ACCG", 100),
            ],
            vec!["AACC", "ACCG"],
        );

        prune_tips(&mut result, 4);

        assert!(result.nodes.contains_key("AACC"));
        assert!(result.nodes.contains_key("ACCG"));
    }

    #[test]
    fn test_prune_tips_empty_result() {
        let mut result = WalkResult {
            nodes: HashMap::new(),
            reference_kmers: HashSet::new(),
        };

        // Should not panic
        prune_tips(&mut result, 4);
        assert!(result.nodes.is_empty());
    }

    #[test]
    fn test_prune_bubbles_removes_low_coverage() {
        // Reference k-mers at 1000x, variant at 50x (above 1% threshold),
        // error bubble at 5x (below 1% of 1000 = 10).
        let mut result = make_walk_result(
            vec![
                ("AACC", 1000),
                ("ACCG", 1000),
                ("CCGG", 1000),
                ("CCGT", 50),  // real variant: 50 >= 10 (1% of 1000)
                ("CCGA", 5),   // error bubble: 5 < 10
            ],
            vec!["AACC", "ACCG", "CCGG"],
        );

        prune_bubbles(&mut result, 0.01);

        // Error bubble should be removed
        assert!(
            !result.nodes.contains_key("CCGA"),
            "Low-coverage bubble CCGA should be pruned"
        );
        // Real variant should remain
        assert!(
            result.nodes.contains_key("CCGT"),
            "Variant CCGT (50x) should survive bubble pruning"
        );
        // Reference k-mers should remain
        assert!(result.nodes.contains_key("AACC"));
        assert!(result.nodes.contains_key("ACCG"));
        assert!(result.nodes.contains_key("CCGG"));
    }

    #[test]
    fn test_prune_bubbles_keeps_all_above_threshold() {
        let mut result = make_walk_result(
            vec![
                ("AACC", 100),
                ("ACCG", 100),
                ("CCGT", 50), // 50 >= 1 (1% of 100)
            ],
            vec!["AACC", "ACCG"],
        );

        prune_bubbles(&mut result, 0.01);

        assert_eq!(result.nodes.len(), 3, "All k-mers should survive");
    }

    #[test]
    fn test_prune_bubbles_no_reference_kmers() {
        // Edge case: no reference k-mers
        let mut result = make_walk_result(
            vec![("AACC", 100)],
            vec![],
        );

        // Should not panic
        prune_bubbles(&mut result, 0.01);
        assert_eq!(result.nodes.len(), 1);
    }

    #[test]
    fn test_prune_bubbles_preserves_reference_even_if_low() {
        // Reference k-mer with count 0 should never be removed
        let mut result = make_walk_result(
            vec![
                ("AACC", 1000),
                ("ACCG", 0), // reference with zero count
                ("CCGT", 5), // non-ref low count
            ],
            vec!["AACC", "ACCG"],
        );

        prune_bubbles(&mut result, 0.01);

        assert!(
            result.nodes.contains_key("ACCG"),
            "Reference k-mer ACCG should survive even with count 0"
        );
    }
}
