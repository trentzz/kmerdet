/// Path quantification via NNLS linear regression.
///
/// Estimates the expression level of each path by solving:
///   contrib @ coef ≈ counts
///
/// Where:
/// - contrib[i][j] = number of times k-mer i appears in path j
/// - counts[i] = observed k-mer count from jellyfish
/// - coef[j] = estimated expression for path j
///
/// Uses non-negative least squares (NNLS) via SVD decomposition
/// with gradient descent refinement to eliminate negative coefficients.

use crate::jellyfish::KmerDatabase;
use crate::sequence::path::KmerPath;

/// Result of path quantification.
#[derive(Debug, Clone)]
pub struct Quantification {
    /// Expression coefficient per path.
    pub coefficients: Vec<f64>,
    /// Relative VAF per path (coefficient / sum).
    pub rvafs: Vec<f64>,
    /// Minimum k-mer count along each path.
    pub min_coverages: Vec<u64>,
}

/// Quantify the expression levels of multiple paths through a k-mer region.
pub fn quantify(
    _paths: &[KmerPath],
    _db: &dyn KmerDatabase,
) -> Quantification {
    // TODO Phase 4: Build contribution matrix + NNLS solve
    todo!("quantifier not yet implemented")
}
