/// Overlapping mutation clustering.
///
/// When multiple mutations overlap in k-mer space, they must be quantified
/// together as a cluster. This handles compound heterozygous mutations where
/// individual path quantification would produce zero min-coverage.

use super::VariantCall;

/// A cluster of overlapping variant calls.
#[derive(Debug)]
pub struct VariantCluster {
    /// The overlapping calls in this cluster.
    pub calls: Vec<VariantCall>,
    /// Start position of the cluster region.
    pub start: usize,
    /// End position of the cluster region.
    pub end: usize,
}

/// Parse start and end positions from a variant_name string.
///
/// The variant_name format is "start:ref/alt:end".
/// Returns (start, end) as usize values.
fn parse_positions(variant_name: &str) -> Option<(usize, usize)> {
    let parts: Vec<&str> = variant_name.split(':').collect();
    if parts.len() < 3 {
        return None;
    }
    let start = parts[0].parse::<usize>().ok()?;
    let end = parts[parts.len() - 1].parse::<usize>().ok()?;
    Some((start, end))
}

/// Group overlapping variant calls into clusters for joint quantification.
///
/// Two calls overlap if their variant regions [start, end] share any position.
/// Overlapping is transitive: if A overlaps B and B overlaps C, then A, B, and C
/// are all in the same cluster even if A and C do not directly overlap.
pub fn cluster_variants(calls: &[VariantCall]) -> Vec<VariantCluster> {
    if calls.is_empty() {
        return vec![];
    }

    // Parse positions and pair with indices.
    let mut positioned: Vec<(usize, usize, usize)> = calls
        .iter()
        .enumerate()
        .filter_map(|(idx, call)| {
            parse_positions(&call.variant_name).map(|(start, end)| (idx, start, end))
        })
        .collect();

    if positioned.is_empty() {
        return vec![];
    }

    // Sort by start position, then by end position.
    positioned.sort_by_key(|&(_, start, end)| (start, end));

    // Walk through sorted calls, grouping overlapping ones.
    let mut clusters: Vec<VariantCluster> = Vec::new();

    let (first_idx, first_start, first_end) = positioned[0];
    let mut current_calls = vec![calls[first_idx].clone()];
    let mut cluster_start = first_start;
    let mut cluster_end = first_end;

    for &(idx, start, end) in &positioned[1..] {
        if start <= cluster_end {
            // Overlapping — extend the current cluster.
            current_calls.push(calls[idx].clone());
            if end > cluster_end {
                cluster_end = end;
            }
            if start < cluster_start {
                cluster_start = start;
            }
        } else {
            // No overlap — finalize current cluster, start a new one.
            clusters.push(VariantCluster {
                calls: current_calls,
                start: cluster_start,
                end: cluster_end,
            });
            current_calls = vec![calls[idx].clone()];
            cluster_start = start;
            cluster_end = end;
        }
    }

    // Finalize the last cluster.
    clusters.push(VariantCluster {
        calls: current_calls,
        start: cluster_start,
        end: cluster_end,
    });

    clusters
}
