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

/// Group overlapping variant calls into clusters for joint quantification.
pub fn cluster_variants(_calls: &[VariantCall]) -> Vec<VariantCluster> {
    // TODO Phase 4: Implement overlapping mutation clustering
    todo!("clustering not yet implemented")
}
