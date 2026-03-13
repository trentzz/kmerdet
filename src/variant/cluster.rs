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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::variant::VariantType;

    /// Helper to create a minimal VariantCall for testing.
    fn make_call(variant_name: &str, variant_type: VariantType) -> VariantCall {
        VariantCall {
            sample: "test".to_string(),
            target: "target1".to_string(),
            variant_type,
            variant_name: variant_name.to_string(),
            rvaf: 0.1,
            expression: 100.0,
            min_coverage: 10,
            path_score: 10,
            start_kmer_count: 0,
            ref_sequence: String::new(),
            alt_sequence: String::new(),
            info: "vs_ref".to_string(),
            chrom: None,
            pos: None,
            ref_allele: None,
            alt_allele: None,
            pvalue: None,
            qual: None,
            ci_lower: None,
            ci_upper: None,
        }
    }

    #[test]
    fn test_cluster_empty_calls() {
        let clusters = cluster_variants(&[]);
        assert!(clusters.is_empty());
    }

    #[test]
    fn test_cluster_single_call() {
        let calls = vec![make_call("10:A/T:10", VariantType::Substitution)];
        let clusters = cluster_variants(&calls);
        assert_eq!(clusters.len(), 1);
        assert_eq!(clusters[0].calls.len(), 1);
        assert_eq!(clusters[0].start, 10);
        assert_eq!(clusters[0].end, 10);
    }

    #[test]
    fn test_cluster_non_overlapping_calls() {
        let calls = vec![
            make_call("10:A/T:10", VariantType::Substitution),
            make_call("20:G/C:20", VariantType::Substitution),
        ];
        let clusters = cluster_variants(&calls);
        assert_eq!(clusters.len(), 2);
        assert_eq!(clusters[0].calls.len(), 1);
        assert_eq!(clusters[1].calls.len(), 1);
    }

    #[test]
    fn test_cluster_overlapping_calls() {
        // Two variants that overlap: [10,15] and [13,18]
        let calls = vec![
            make_call("10:ACGTAC/A:15", VariantType::Deletion),
            make_call("13:G/GTTTT:18", VariantType::Insertion),
        ];
        let clusters = cluster_variants(&calls);
        assert_eq!(clusters.len(), 1);
        assert_eq!(clusters[0].calls.len(), 2);
        assert_eq!(clusters[0].start, 10);
        assert_eq!(clusters[0].end, 18);
    }

    #[test]
    fn test_cluster_transitive_overlap() {
        // A overlaps B, B overlaps C, but A does not directly overlap C.
        // A=[10,14], B=[13,17], C=[16,20]
        let calls = vec![
            make_call("10:ACGT/A:14", VariantType::Deletion),
            make_call("13:G/GTTTT:17", VariantType::Insertion),
            make_call("16:A/T:20", VariantType::Substitution),
        ];
        let clusters = cluster_variants(&calls);
        assert_eq!(clusters.len(), 1);
        assert_eq!(clusters[0].calls.len(), 3);
        assert_eq!(clusters[0].start, 10);
        assert_eq!(clusters[0].end, 20);
    }

    #[test]
    fn test_cluster_mixed_overlap_and_separate() {
        // Two overlapping [10,15]+[13,18] and one separate [30,30]
        let calls = vec![
            make_call("10:ACGTAC/A:15", VariantType::Deletion),
            make_call("30:A/T:30", VariantType::Substitution),
            make_call("13:G/GTTTT:18", VariantType::Insertion),
        ];
        let clusters = cluster_variants(&calls);
        assert_eq!(clusters.len(), 2);
        // First cluster should be the overlapping pair (sorted by start)
        assert_eq!(clusters[0].calls.len(), 2);
        assert_eq!(clusters[0].start, 10);
        assert_eq!(clusters[0].end, 18);
        // Second cluster is the standalone SNV
        assert_eq!(clusters[1].calls.len(), 1);
        assert_eq!(clusters[1].start, 30);
        assert_eq!(clusters[1].end, 30);
    }

    #[test]
    fn test_cluster_unparseable_variant_names() {
        // Calls with unparseable variant names should be skipped.
        let calls = vec![
            make_call("invalid_name", VariantType::Reference),
            make_call("10:A/T:10", VariantType::Substitution),
        ];
        let clusters = cluster_variants(&calls);
        assert_eq!(clusters.len(), 1);
        assert_eq!(clusters[0].calls.len(), 1);
        assert_eq!(clusters[0].start, 10);
    }

    #[test]
    fn test_cluster_exact_boundary_overlap() {
        // Calls touching at exactly the same position: [10,15] and [15,20]
        let calls = vec![
            make_call("10:ACGT/A:15", VariantType::Deletion),
            make_call("15:G/GTTTT:20", VariantType::Insertion),
        ];
        let clusters = cluster_variants(&calls);
        assert_eq!(clusters.len(), 1);
        assert_eq!(clusters[0].calls.len(), 2);
    }

    #[test]
    fn test_cluster_adjacent_non_overlapping() {
        // Calls adjacent but not overlapping: [10,14] and [15,20]
        let calls = vec![
            make_call("10:ACGT/A:14", VariantType::Deletion),
            make_call("15:G/GTTTT:20", VariantType::Insertion),
        ];
        let clusters = cluster_variants(&calls);
        assert_eq!(clusters.len(), 2);
    }
}
