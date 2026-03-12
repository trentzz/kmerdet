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

use indexmap::IndexSet;
use nalgebra::{DMatrix, DVector};

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
    /// Lower bound of bootstrap CI per path (empty if bootstrap not run).
    pub ci_lower: Vec<f64>,
    /// Upper bound of bootstrap CI per path (empty if bootstrap not run).
    pub ci_upper: Vec<f64>,
}

/// Quantify the expression levels of multiple paths through a k-mer region.
pub fn quantify(
    paths: &[KmerPath],
    db: &dyn KmerDatabase,
) -> Quantification {
    if paths.is_empty() {
        return Quantification {
            coefficients: vec![],
            rvafs: vec![],
            min_coverages: vec![],
            ci_lower: vec![],
            ci_upper: vec![],
        };
    }

    let n_paths = paths.len();

    // 1. Collect all unique k-mers across all paths, preserving insertion order.
    let mut all_kmers: IndexSet<String> = IndexSet::new();
    for path in paths {
        for kmer in &path.kmers {
            all_kmers.insert(kmer.to_uppercase());
        }
    }
    let n_kmers = all_kmers.len();

    if n_kmers == 0 {
        return Quantification {
            coefficients: vec![0.0; n_paths],
            rvafs: vec![0.0; n_paths],
            min_coverages: vec![0; n_paths],
            ci_lower: vec![],
            ci_upper: vec![],
        };
    }

    // 2. Build contribution matrix: contrib[i][j] = count of k-mer i in path j.
    //    For ITDs, a k-mer can appear multiple times in a single path.
    let mut contrib = vec![vec![0u32; n_paths]; n_kmers];
    for (j, path) in paths.iter().enumerate() {
        for kmer in &path.kmers {
            let upper = kmer.to_uppercase();
            if let Some(i) = all_kmers.get_index_of(&upper) {
                contrib[i][j] += 1;
            }
        }
    }

    // 3. Build observed count vector.
    let counts: Vec<f64> = all_kmers
        .iter()
        .map(|kmer| db.query(kmer) as f64)
        .collect();

    // 4. Solve NNLS via SVD with iterative negative coefficient clipping.
    let a = DMatrix::from_fn(n_kmers, n_paths, |i, j| contrib[i][j] as f64);
    let b = DVector::from_fn(n_kmers, |i, _| counts[i]);

    let coefficients = solve_nnls(&a, &b);

    // 5. Compute rVAFs.
    let total: f64 = coefficients.iter().sum();
    let rvafs: Vec<f64> = if total > 0.0 {
        coefficients.iter().map(|&c| c / total).collect()
    } else {
        vec![0.0; n_paths]
    };

    // 6. Compute min coverages.
    let min_coverages: Vec<u64> = paths
        .iter()
        .map(|path| {
            if path.kmers.is_empty() {
                return 0;
            }
            path.kmers
                .iter()
                .map(|kmer| db.query(&kmer.to_uppercase()))
                .min()
                .unwrap_or(0)
        })
        .collect();

    Quantification {
        coefficients,
        rvafs,
        min_coverages,
        ci_lower: vec![],
        ci_upper: vec![],
    }
}

/// Solve non-negative least squares: find c >= 0 minimizing ||A*c - b||^2.
///
/// Uses SVD-based pseudoinverse with iterative clipping of negative coefficients.
/// Active columns are re-solved until all remaining coefficients are non-negative.
fn solve_nnls(a: &DMatrix<f64>, b: &DVector<f64>) -> Vec<f64> {
    let n_cols = a.ncols();

    // Track which columns (paths) are still active.
    let mut active: Vec<bool> = vec![true; n_cols];
    let mut result = vec![0.0f64; n_cols];

    loop {
        // Gather active column indices.
        let active_indices: Vec<usize> = (0..n_cols).filter(|&j| active[j]).collect();

        if active_indices.is_empty() {
            break;
        }

        // Build sub-matrix with only active columns.
        let n_active = active_indices.len();
        let a_sub = DMatrix::from_fn(a.nrows(), n_active, |i, j| {
            a[(i, active_indices[j])]
        });

        // Solve via SVD pseudoinverse.
        let svd = a_sub.svd(true, true);
        let coeffs = match svd.solve(b, 1e-10) {
            Ok(c) => c,
            Err(_) => DVector::zeros(n_active),
        };

        // Check for negative coefficients.
        let mut any_negative = false;
        for (idx, &ai) in active_indices.iter().enumerate() {
            if coeffs[idx] < 0.0 {
                active[ai] = false;
                result[ai] = 0.0;
                any_negative = true;
            } else {
                result[ai] = coeffs[idx];
            }
        }

        if !any_negative {
            break;
        }
    }

    result
}
