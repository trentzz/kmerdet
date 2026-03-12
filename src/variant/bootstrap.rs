/// Bootstrap confidence intervals for rVAF estimates.
///
/// After NNLS quantification produces a point estimate of rVAF for each path,
/// bootstrap resampling provides confidence intervals by:
///
/// 1. Resampling each observed k-mer count from Poisson(observed_count)
/// 2. Re-running NNLS with the resampled counts
/// 3. Computing rVAFs for each replicate
/// 4. Taking percentiles of the replicate rVAFs as CI bounds

use indexmap::IndexSet;
use nalgebra::{DMatrix, DVector};

use crate::jellyfish::KmerDatabase;
use crate::sequence::path::KmerPath;

/// Configuration for bootstrap resampling.
#[derive(Debug, Clone)]
pub struct BootstrapConfig {
    /// Number of bootstrap replicates (default 1000).
    pub n_replicates: usize,
    /// Confidence level (default 0.95 for 95% CI).
    pub confidence_level: f64,
    /// Random seed for reproducibility (optional).
    pub seed: Option<u64>,
}

impl Default for BootstrapConfig {
    fn default() -> Self {
        Self {
            n_replicates: 1000,
            confidence_level: 0.95,
            seed: Some(42),
        }
    }
}

/// Bootstrap confidence interval result for a single path.
#[derive(Debug, Clone)]
pub struct BootstrapCI {
    /// Lower bound of CI.
    pub ci_lower: f64,
    /// Upper bound of CI.
    pub ci_upper: f64,
    /// Standard error of rVAF estimates across replicates.
    pub std_error: f64,
}

/// A simple xorshift64 pseudo-random number generator.
struct Xorshift64 {
    state: u64,
}

impl Xorshift64 {
    fn new(seed: u64) -> Self {
        Self {
            state: if seed == 0 { 1 } else { seed },
        }
    }

    fn next_u64(&mut self) -> u64 {
        let mut x = self.state;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.state = x;
        x
    }

    fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }

    /// Sample from Poisson(lambda).
    ///
    /// For lambda < 30: Knuth's inverse CDF algorithm.
    /// For lambda >= 30: normal approximation N(lambda, lambda).
    fn poisson(&mut self, lambda: f64) -> u64 {
        if lambda <= 0.0 {
            return 0;
        }
        if lambda < 30.0 {
            // Knuth's algorithm
            let l = (-lambda).exp();
            let mut k = 0u64;
            let mut p = 1.0;
            loop {
                k += 1;
                p *= self.next_f64();
                if p < l {
                    break;
                }
            }
            k - 1
        } else {
            // Normal approximation via Box-Muller
            let u1 = self.next_f64().max(f64::MIN_POSITIVE); // avoid ln(0)
            let u2 = self.next_f64();
            let z = (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos();
            (lambda + lambda.sqrt() * z).round().max(0.0) as u64
        }
    }
}

/// Run bootstrap resampling on NNLS quantification.
///
/// For each replicate:
/// 1. Resample each k-mer count from Poisson(observed_count)
/// 2. Run NNLS with resampled counts
/// 3. Compute rVAFs for this replicate
///
/// Return CIs for each path.
pub fn bootstrap_confidence_intervals(
    paths: &[KmerPath],
    db: &dyn KmerDatabase,
    config: &BootstrapConfig,
) -> Vec<BootstrapCI> {
    let n_paths = paths.len();

    if n_paths == 0 || config.n_replicates == 0 {
        return vec![];
    }

    // 1. Collect all unique k-mers across all paths.
    let mut all_kmers: IndexSet<String> = IndexSet::new();
    for path in paths {
        for kmer in &path.kmers {
            all_kmers.insert(kmer.to_uppercase());
        }
    }
    let n_kmers = all_kmers.len();

    if n_kmers == 0 {
        return paths
            .iter()
            .map(|_| BootstrapCI {
                ci_lower: 0.0,
                ci_upper: 0.0,
                std_error: 0.0,
            })
            .collect();
    }

    // 2. Build contribution matrix.
    let mut contrib = vec![vec![0u32; n_paths]; n_kmers];
    for (j, path) in paths.iter().enumerate() {
        for kmer in &path.kmers {
            let upper = kmer.to_uppercase();
            if let Some(i) = all_kmers.get_index_of(&upper) {
                contrib[i][j] += 1;
            }
        }
    }

    let a = DMatrix::from_fn(n_kmers, n_paths, |i, j| contrib[i][j] as f64);

    // 3. Get observed counts.
    let observed: Vec<f64> = all_kmers.iter().map(|kmer| db.query(kmer) as f64).collect();

    // 4. Initialize RNG.
    let mut rng = Xorshift64::new(config.seed.unwrap_or(42));

    // 5. Run bootstrap replicates.
    // rvaf_samples[j] collects all rVAF samples for path j.
    let mut rvaf_samples: Vec<Vec<f64>> = vec![Vec::with_capacity(config.n_replicates); n_paths];

    for _ in 0..config.n_replicates {
        // Resample counts from Poisson(observed_count).
        let resampled: Vec<f64> = observed.iter().map(|&c| rng.poisson(c) as f64).collect();

        let b = DVector::from_fn(n_kmers, |i, _| resampled[i]);

        // Solve NNLS with resampled counts.
        let coefficients = solve_nnls(&a, &b);

        // Compute rVAFs.
        let total: f64 = coefficients.iter().sum();
        let rvafs: Vec<f64> = if total > 0.0 {
            coefficients.iter().map(|&c| c / total).collect()
        } else {
            vec![0.0; n_paths]
        };

        for (j, &rvaf) in rvafs.iter().enumerate() {
            rvaf_samples[j].push(rvaf);
        }
    }

    // 6. Compute CIs from the collected samples.
    let alpha = 1.0 - config.confidence_level;
    let lower_pct = alpha / 2.0;
    let upper_pct = 1.0 - alpha / 2.0;

    rvaf_samples
        .iter()
        .map(|samples| {
            if samples.is_empty() {
                return BootstrapCI {
                    ci_lower: 0.0,
                    ci_upper: 0.0,
                    std_error: 0.0,
                };
            }

            let mut sorted = samples.clone();
            sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

            let ci_lower = percentile(&sorted, lower_pct);
            let ci_upper = percentile(&sorted, upper_pct);

            // Standard error = standard deviation of the bootstrap distribution.
            let mean: f64 = samples.iter().sum::<f64>() / samples.len() as f64;
            let variance: f64 =
                samples.iter().map(|&x| (x - mean).powi(2)).sum::<f64>() / samples.len() as f64;
            let std_error = variance.sqrt();

            BootstrapCI {
                ci_lower,
                ci_upper,
                std_error,
            }
        })
        .collect()
}

/// Compute a percentile from sorted data using linear interpolation.
fn percentile(sorted: &[f64], p: f64) -> f64 {
    if sorted.is_empty() {
        return 0.0;
    }
    if sorted.len() == 1 {
        return sorted[0];
    }

    let idx = p * (sorted.len() - 1) as f64;
    let lo = idx.floor() as usize;
    let hi = idx.ceil() as usize;
    let frac = idx - lo as f64;

    let lo = lo.min(sorted.len() - 1);
    let hi = hi.min(sorted.len() - 1);

    sorted[lo] * (1.0 - frac) + sorted[hi] * frac
}

/// Solve non-negative least squares: find c >= 0 minimizing ||A*c - b||^2.
///
/// Uses SVD-based pseudoinverse with iterative clipping of negative coefficients.
/// This is a copy of the NNLS solver from quantifier.rs to avoid circular
/// dependencies in the bootstrap module.
fn solve_nnls(a: &DMatrix<f64>, b: &DVector<f64>) -> Vec<f64> {
    let n_cols = a.ncols();

    let mut active: Vec<bool> = vec![true; n_cols];
    let mut result = vec![0.0f64; n_cols];

    loop {
        let active_indices: Vec<usize> = (0..n_cols).filter(|&j| active[j]).collect();

        if active_indices.is_empty() {
            break;
        }

        let n_active = active_indices.len();
        let a_sub = DMatrix::from_fn(a.nrows(), n_active, |i, j| a[(i, active_indices[j])]);

        let svd = a_sub.svd(true, true);
        let coeffs = match svd.solve(b, 1e-10) {
            Ok(c) => c,
            Err(_) => DVector::zeros(n_active),
        };

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

#[cfg(test)]
mod tests {
    use super::*;

    // -----------------------------------------------------------------------
    // RNG / Poisson sampler tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_poisson_mean_low_lambda() {
        // For Poisson(lambda=10), the mean of many samples should be close to 10.
        let mut rng = Xorshift64::new(12345);
        let n = 10_000;
        let lambda = 10.0;
        let sum: u64 = (0..n).map(|_| rng.poisson(lambda)).sum();
        let mean = sum as f64 / n as f64;
        assert!(
            (mean - lambda).abs() < 1.0,
            "Poisson mean should be ~{}, got {}",
            lambda,
            mean
        );
    }

    #[test]
    fn test_poisson_mean_high_lambda() {
        // For Poisson(lambda=100), use the normal approximation path.
        let mut rng = Xorshift64::new(54321);
        let n = 10_000;
        let lambda = 100.0;
        let sum: u64 = (0..n).map(|_| rng.poisson(lambda)).sum();
        let mean = sum as f64 / n as f64;
        assert!(
            (mean - lambda).abs() < 3.0,
            "Poisson mean should be ~{}, got {}",
            lambda,
            mean
        );
    }

    #[test]
    fn test_poisson_zero_lambda() {
        let mut rng = Xorshift64::new(42);
        for _ in 0..100 {
            assert_eq!(rng.poisson(0.0), 0);
        }
    }

    // -----------------------------------------------------------------------
    // A mock KmerDatabase for bootstrap tests
    // -----------------------------------------------------------------------

    struct MockDb {
        counts: std::collections::HashMap<String, u64>,
        k: u8,
    }

    impl MockDb {
        fn new(k: u8) -> Self {
            Self {
                counts: std::collections::HashMap::new(),
                k,
            }
        }

        fn set(&mut self, kmer: &str, count: u64) {
            self.counts.insert(kmer.to_uppercase(), count);
        }
    }

    impl KmerDatabase for MockDb {
        fn query(&self, kmer: &str) -> u64 {
            *self.counts.get(&kmer.to_uppercase()).unwrap_or(&0)
        }

        fn kmer_length(&self) -> u8 {
            self.k
        }
    }

    // -----------------------------------------------------------------------
    // Bootstrap CI tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_bootstrap_single_path_tight_ci() {
        // A single path should have rVAF ~1.0 with a very tight CI.
        let mut db = MockDb::new(4);
        db.set("ACGT", 1000);
        db.set("CGTA", 1000);
        db.set("GTAC", 1000);

        let path = KmerPath {
            kmers: vec![
                "ACGT".to_string(),
                "CGTA".to_string(),
                "GTAC".to_string(),
            ],
            is_reference: true,
        };

        let config = BootstrapConfig {
            n_replicates: 500,
            confidence_level: 0.95,
            seed: Some(42),
        };

        let cis = bootstrap_confidence_intervals(&[path], &db, &config);

        assert_eq!(cis.len(), 1);
        // Single path: rVAF is always 1.0 regardless of resampled counts.
        assert!(
            (cis[0].ci_lower - 1.0).abs() < 1e-9,
            "single path CI lower should be 1.0, got {}",
            cis[0].ci_lower
        );
        assert!(
            (cis[0].ci_upper - 1.0).abs() < 1e-9,
            "single path CI upper should be 1.0, got {}",
            cis[0].ci_upper
        );
        assert!(
            cis[0].std_error < 1e-9,
            "single path std_error should be ~0, got {}",
            cis[0].std_error
        );
    }

    #[test]
    fn test_bootstrap_two_paths_high_coverage_narrow_ci() {
        // Two paths, high coverage -> narrow CI.
        let mut db = MockDb::new(4);
        // Shared k-mer at 2000 (1000 ref + 1000 alt).
        db.set("AAAA", 2000);
        // Ref-only at 1000.
        db.set("AAAC", 1000);
        db.set("AACG", 1000);
        // Alt-only at 1000.
        db.set("AAAT", 1000);
        db.set("AATG", 1000);

        let ref_path = KmerPath {
            kmers: vec![
                "AAAA".to_string(),
                "AAAC".to_string(),
                "AACG".to_string(),
            ],
            is_reference: true,
        };
        let alt_path = KmerPath {
            kmers: vec![
                "AAAA".to_string(),
                "AAAT".to_string(),
                "AATG".to_string(),
            ],
            is_reference: false,
        };

        let config = BootstrapConfig {
            n_replicates: 1000,
            confidence_level: 0.95,
            seed: Some(42),
        };

        let cis = bootstrap_confidence_intervals(&[ref_path, alt_path], &db, &config);

        assert_eq!(cis.len(), 2);

        // Both paths should have rVAF ~0.5, CI should be narrow for high coverage.
        let ci_width_0 = cis[0].ci_upper - cis[0].ci_lower;
        let ci_width_1 = cis[1].ci_upper - cis[1].ci_lower;

        assert!(
            ci_width_0 < 0.10,
            "high coverage CI width should be <0.10, got {}",
            ci_width_0
        );
        assert!(
            ci_width_1 < 0.10,
            "high coverage CI width should be <0.10, got {}",
            ci_width_1
        );

        // CI should contain the true value (~0.5).
        assert!(
            cis[0].ci_lower <= 0.5 && cis[0].ci_upper >= 0.5,
            "CI for path 0 should contain 0.5: [{}, {}]",
            cis[0].ci_lower,
            cis[0].ci_upper
        );
        assert!(
            cis[1].ci_lower <= 0.5 && cis[1].ci_upper >= 0.5,
            "CI for path 1 should contain 0.5: [{}, {}]",
            cis[1].ci_lower,
            cis[1].ci_upper
        );
    }

    #[test]
    fn test_bootstrap_two_paths_low_coverage_wider_ci() {
        // Two paths, low coverage -> wider CI than high coverage.
        let mut db = MockDb::new(4);
        // Shared k-mer at 20 (10 ref + 10 alt).
        db.set("AAAA", 20);
        // Ref-only at 10.
        db.set("AAAC", 10);
        db.set("AACG", 10);
        // Alt-only at 10.
        db.set("AAAT", 10);
        db.set("AATG", 10);

        let ref_path = KmerPath {
            kmers: vec![
                "AAAA".to_string(),
                "AAAC".to_string(),
                "AACG".to_string(),
            ],
            is_reference: true,
        };
        let alt_path = KmerPath {
            kmers: vec![
                "AAAA".to_string(),
                "AAAT".to_string(),
                "AATG".to_string(),
            ],
            is_reference: false,
        };

        let config = BootstrapConfig {
            n_replicates: 1000,
            confidence_level: 0.95,
            seed: Some(42),
        };

        let cis = bootstrap_confidence_intervals(&[ref_path, alt_path], &db, &config);

        assert_eq!(cis.len(), 2);

        // CI should be wider than with high coverage.
        let ci_width_0 = cis[0].ci_upper - cis[0].ci_lower;
        assert!(
            ci_width_0 > 0.02,
            "low coverage CI width should be >0.02, got {}",
            ci_width_0
        );

        // But CI should still contain ~0.5.
        assert!(
            cis[0].ci_lower <= 0.55 && cis[0].ci_upper >= 0.45,
            "CI for path 0 should roughly contain 0.5: [{}, {}]",
            cis[0].ci_lower,
            cis[0].ci_upper
        );
    }

    #[test]
    fn test_bootstrap_reproducibility() {
        // Same seed should produce identical results.
        let mut db = MockDb::new(4);
        db.set("AAAA", 200);
        db.set("AAAC", 100);
        db.set("AACG", 100);
        db.set("AAAT", 100);
        db.set("AATG", 100);

        let ref_path = KmerPath {
            kmers: vec![
                "AAAA".to_string(),
                "AAAC".to_string(),
                "AACG".to_string(),
            ],
            is_reference: true,
        };
        let alt_path = KmerPath {
            kmers: vec![
                "AAAA".to_string(),
                "AAAT".to_string(),
                "AATG".to_string(),
            ],
            is_reference: false,
        };

        let config = BootstrapConfig {
            n_replicates: 200,
            confidence_level: 0.95,
            seed: Some(999),
        };

        let cis1 =
            bootstrap_confidence_intervals(&[ref_path.clone(), alt_path.clone()], &db, &config);
        let cis2 = bootstrap_confidence_intervals(&[ref_path, alt_path], &db, &config);

        assert_eq!(cis1.len(), cis2.len());
        for (a, b) in cis1.iter().zip(cis2.iter()) {
            assert!(
                (a.ci_lower - b.ci_lower).abs() < 1e-12,
                "reproducibility: ci_lower mismatch"
            );
            assert!(
                (a.ci_upper - b.ci_upper).abs() < 1e-12,
                "reproducibility: ci_upper mismatch"
            );
            assert!(
                (a.std_error - b.std_error).abs() < 1e-12,
                "reproducibility: std_error mismatch"
            );
        }
    }

    #[test]
    fn test_bootstrap_zero_replicates() {
        let mut db = MockDb::new(4);
        db.set("ACGT", 100);

        let path = KmerPath {
            kmers: vec!["ACGT".to_string()],
            is_reference: true,
        };

        let config = BootstrapConfig {
            n_replicates: 0,
            confidence_level: 0.95,
            seed: Some(42),
        };

        let cis = bootstrap_confidence_intervals(&[path], &db, &config);
        assert!(cis.is_empty(), "zero replicates should return empty CIs");
    }

    #[test]
    fn test_bootstrap_empty_paths() {
        let db = MockDb::new(4);
        let config = BootstrapConfig::default();
        let cis = bootstrap_confidence_intervals(&[], &db, &config);
        assert!(cis.is_empty(), "empty paths should return empty CIs");
    }

    #[test]
    fn test_percentile_basic() {
        let data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        assert!((percentile(&data, 0.0) - 1.0).abs() < 1e-10);
        assert!((percentile(&data, 0.5) - 3.0).abs() < 1e-10);
        assert!((percentile(&data, 1.0) - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_percentile_interpolation() {
        let data = vec![0.0, 10.0];
        assert!((percentile(&data, 0.25) - 2.5).abs() < 1e-10);
        assert!((percentile(&data, 0.75) - 7.5).abs() < 1e-10);
    }

    #[test]
    fn test_percentile_empty() {
        assert_eq!(percentile(&[], 0.5), 0.0);
    }

    #[test]
    fn test_percentile_single() {
        assert_eq!(percentile(&[42.0], 0.5), 42.0);
    }
}
