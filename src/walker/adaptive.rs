/// Adaptive walking thresholds that compute per-sample error rates and coverage
/// to set statistically-motivated count/ratio thresholds.
///
/// The fixed `--count` / `--ratio` defaults (2 / 0.05) work well at ~1000x coverage
/// but perform poorly at extremes: at 100x, real low-VAF variants are missed; at
/// 50 000x, sequencing errors pass the ratio filter. Adaptive thresholds solve this
/// by sampling reference k-mers from the database to estimate the sequencing error
/// rate, then computing a Poisson-based count threshold that controls the false
/// extension rate.

use crate::jellyfish::KmerDatabase;

/// Coverage tier for automatic parameter selection.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CoverageTier {
    /// < 500x
    UltraLow,
    /// 500-2000x
    Standard,
    /// 2000-10000x
    HighDepth,
    /// > 10000x
    UltraDeep,
}

/// Configuration for adaptive thresholds.
#[derive(Debug, Clone)]
pub struct AdaptiveConfig {
    /// Target false extension rate (default 0.01).
    pub false_extension_rate: f64,
    /// Fallback error rate when estimation fails (default 0.001).
    pub default_error_rate: f64,
    /// Number of reference k-mers to sample for estimation (default 100).
    pub sample_size: usize,
}

/// Computed adaptive thresholds for a specific sample/target.
#[derive(Debug, Clone)]
pub struct AdaptiveThresholds {
    /// Estimated per-base error rate.
    pub error_rate: f64,
    /// Median reference coverage.
    pub median_coverage: f64,
    /// Detected coverage tier.
    pub tier: CoverageTier,
    /// Computed absolute count threshold.
    pub count_threshold: u32,
    /// Computed ratio threshold.
    pub ratio_threshold: f64,
}

impl Default for AdaptiveConfig {
    fn default() -> Self {
        Self {
            false_extension_rate: 0.01,
            default_error_rate: 0.001,
            sample_size: 100,
        }
    }
}

impl std::fmt::Display for CoverageTier {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            CoverageTier::UltraLow => write!(f, "UltraLow (<500x)"),
            CoverageTier::Standard => write!(f, "Standard (500-2000x)"),
            CoverageTier::HighDepth => write!(f, "HighDepth (2000-10000x)"),
            CoverageTier::UltraDeep => write!(f, "UltraDeep (>10000x)"),
        }
    }
}

/// Classify a median coverage value into its [`CoverageTier`].
pub fn classify_tier(median_coverage: f64) -> CoverageTier {
    match median_coverage as u64 {
        0..=499 => CoverageTier::UltraLow,
        500..=1999 => CoverageTier::Standard,
        2000..=9999 => CoverageTier::HighDepth,
        _ => CoverageTier::UltraDeep,
    }
}

/// Estimate adaptive thresholds from reference k-mer coverage.
///
/// Queries reference k-mers to estimate:
/// 1. Median coverage
/// 2. Per-base error rate (from non-reference extensions)
/// 3. Coverage tier
///
/// Then computes appropriate count and ratio thresholds using a Poisson model
/// for the expected number of erroneous extensions at the observed coverage.
pub fn estimate_thresholds(
    db: &dyn KmerDatabase,
    ref_kmers: &[String],
    config: &AdaptiveConfig,
) -> AdaptiveThresholds {
    // Sample ref k-mers (use all if fewer than sample_size)
    let sample: Vec<&String> = if ref_kmers.len() <= config.sample_size {
        ref_kmers.iter().collect()
    } else {
        // Evenly spaced sampling
        let step = ref_kmers.len() / config.sample_size;
        ref_kmers.iter().step_by(step).take(config.sample_size).collect()
    };

    // Query coverage for sampled k-mers
    let mut coverages: Vec<u64> = sample.iter().map(|k| db.query(k)).collect();
    coverages.sort_unstable();

    let median_coverage = compute_median(&coverages);

    // Estimate error rate from non-reference extensions
    let error_rate = estimate_error_rate(db, &sample);
    let error_rate = if error_rate > 0.0 {
        error_rate
    } else {
        config.default_error_rate
    };

    // Determine coverage tier
    let tier = classify_tier(median_coverage);

    // Compute adaptive thresholds
    // Count threshold: Poisson 99th percentile of error distribution
    // For Poisson(lambda), 99th percentile ~ lambda + 2.33 * sqrt(lambda)
    let lambda = median_coverage * error_rate;
    let count_threshold = if lambda > 0.0 {
        let threshold = lambda + 2.33 * lambda.sqrt();
        (threshold.ceil() as u32).max(2) // Never below 2
    } else {
        2
    };

    // Ratio threshold: based on coverage tier
    let ratio_threshold = match tier {
        CoverageTier::UltraLow => 0.10,  // More permissive for low coverage
        CoverageTier::Standard => 0.05,   // Default (matches km)
        CoverageTier::HighDepth => 0.03,  // Tighter for high coverage
        CoverageTier::UltraDeep => 0.01,  // Tight for ultra-deep
    };

    AdaptiveThresholds {
        error_rate,
        median_coverage,
        tier,
        count_threshold,
        ratio_threshold,
    }
}

/// Compute the median of a sorted slice of u64 values.
fn compute_median(sorted: &[u64]) -> f64 {
    if sorted.is_empty() {
        return 0.0;
    }
    let mid = sorted.len() / 2;
    if sorted.len() % 2 == 0 {
        (sorted[mid - 1] + sorted[mid]) as f64 / 2.0
    } else {
        sorted[mid] as f64
    }
}

/// Estimate per-base error rate from non-reference k-mer extensions.
///
/// For each sampled reference k-mer, queries all four possible right-extensions.
/// The dominant extension (highest count) is assumed to be the true continuation;
/// the remaining three are treated as sequencing errors. The overall error rate
/// is `total_error_counts / (3 * total_ref_coverage)`.
fn estimate_error_rate(db: &dyn KmerDatabase, ref_kmers: &[&String]) -> f64 {
    let bases = [b'A', b'C', b'G', b'T'];
    let mut total_ref_coverage: u64 = 0;
    let mut total_error_counts: u64 = 0;

    for kmer_str in ref_kmers {
        let ref_count = db.query(kmer_str);
        if ref_count == 0 {
            continue;
        }
        total_ref_coverage += ref_count;

        let kmer_bytes = kmer_str.as_bytes();

        // Query all 4 right-extensions, find the max (likely the true extension),
        // sum the rest as error counts.
        let mut ext_counts = [0u64; 4];
        for (i, &base) in bases.iter().enumerate() {
            let mut extended = String::with_capacity(kmer_bytes.len());
            for &b in &kmer_bytes[1..] {
                extended.push(b as char);
            }
            extended.push(base as char);
            ext_counts[i] = db.query(&extended);
        }
        let max_ext = ext_counts.iter().copied().max().unwrap_or(0);
        let error_sum: u64 = ext_counts.iter().sum::<u64>() - max_ext;
        total_error_counts += error_sum;
    }

    if total_ref_coverage == 0 {
        return 0.0;
    }

    // Error rate = total error extensions / (3 * total reference coverage)
    // Factor of 3 because there are 3 non-reference extension bases
    total_error_counts as f64 / (3.0 * total_ref_coverage as f64)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;

    /// Inline mock k-mer database for unit testing.
    struct MockDb {
        counts: HashMap<String, u64>,
        k: u8,
    }

    impl MockDb {
        fn new(k: u8) -> Self {
            Self {
                counts: HashMap::new(),
                k,
            }
        }

        fn set(&mut self, kmer: &str, count: u64) {
            self.counts.insert(kmer.to_uppercase(), count);
        }

        /// Set counts for all k-mers in a sequence (sliding window of size k).
        fn set_sequence(&mut self, sequence: &str, count: u64) {
            let seq = sequence.to_uppercase();
            let k = self.k as usize;
            for window in seq.as_bytes().windows(k) {
                let kmer = std::str::from_utf8(window).unwrap().to_string();
                self.counts.insert(kmer, count);
            }
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

    // ---- Coverage tier classification ----

    #[test]
    fn test_classify_tier_ultra_low() {
        assert_eq!(classify_tier(0.0), CoverageTier::UltraLow);
        assert_eq!(classify_tier(100.0), CoverageTier::UltraLow);
        assert_eq!(classify_tier(499.0), CoverageTier::UltraLow);
    }

    #[test]
    fn test_classify_tier_standard() {
        assert_eq!(classify_tier(500.0), CoverageTier::Standard);
        assert_eq!(classify_tier(1000.0), CoverageTier::Standard);
        assert_eq!(classify_tier(1999.0), CoverageTier::Standard);
    }

    #[test]
    fn test_classify_tier_high_depth() {
        assert_eq!(classify_tier(2000.0), CoverageTier::HighDepth);
        assert_eq!(classify_tier(5000.0), CoverageTier::HighDepth);
        assert_eq!(classify_tier(9999.0), CoverageTier::HighDepth);
    }

    #[test]
    fn test_classify_tier_ultra_deep() {
        assert_eq!(classify_tier(10000.0), CoverageTier::UltraDeep);
        assert_eq!(classify_tier(50000.0), CoverageTier::UltraDeep);
    }

    // ---- Boundary tests at tier transitions ----

    #[test]
    fn test_classify_tier_boundaries() {
        // Exact boundary values
        assert_eq!(classify_tier(499.0), CoverageTier::UltraLow);
        assert_eq!(classify_tier(500.0), CoverageTier::Standard);
        assert_eq!(classify_tier(1999.0), CoverageTier::Standard);
        assert_eq!(classify_tier(2000.0), CoverageTier::HighDepth);
        assert_eq!(classify_tier(9999.0), CoverageTier::HighDepth);
        assert_eq!(classify_tier(10000.0), CoverageTier::UltraDeep);
    }

    // ---- estimate_thresholds: high coverage ----

    #[test]
    fn test_estimate_thresholds_high_coverage() {
        let k = 5;
        let mut db = MockDb::new(k);

        // Reference: ACGTACGTAC (6 5-mers at 5000x coverage)
        let ref_seq = "ACGTACGTAC";
        db.set_sequence(ref_seq, 5000);

        // Add small error counts for non-reference extensions
        // For each ref k-mer, the dominant extension is ~5000, errors are ~5 each
        // This simulates ~0.1% error rate
        // E.g., from ACGTA, forward exts: CGTAC(5000, ref), CGTAT(5), CGTAG(5), CGTAA(5)
        db.set("CGTAT", 5);
        db.set("CGTAG", 5);
        db.set("CGTAA", 5);

        let ref_kmers: Vec<String> = ref_seq
            .as_bytes()
            .windows(k as usize)
            .map(|w| std::str::from_utf8(w).unwrap().to_string())
            .collect();

        let config = AdaptiveConfig::default();
        let thresholds = estimate_thresholds(&db, &ref_kmers, &config);

        assert_eq!(thresholds.tier, CoverageTier::HighDepth);
        assert!((thresholds.median_coverage - 5000.0).abs() < 1.0);
        assert!(thresholds.count_threshold >= 2);
        assert!((thresholds.ratio_threshold - 0.03).abs() < f64::EPSILON);
    }

    // ---- estimate_thresholds: low coverage ----

    #[test]
    fn test_estimate_thresholds_low_coverage() {
        let k = 5;
        let mut db = MockDb::new(k);

        let ref_seq = "ACGTACGTAC";
        db.set_sequence(ref_seq, 100); // 100x = UltraLow

        let ref_kmers: Vec<String> = ref_seq
            .as_bytes()
            .windows(k as usize)
            .map(|w| std::str::from_utf8(w).unwrap().to_string())
            .collect();

        let config = AdaptiveConfig::default();
        let thresholds = estimate_thresholds(&db, &ref_kmers, &config);

        assert_eq!(thresholds.tier, CoverageTier::UltraLow);
        assert!((thresholds.median_coverage - 100.0).abs() < 1.0);
        assert!(thresholds.count_threshold >= 2);
        // UltraLow tier gets more permissive ratio
        assert!((thresholds.ratio_threshold - 0.10).abs() < f64::EPSILON);
    }

    // ---- Error rate: clean reference (low error) ----

    #[test]
    fn test_error_rate_clean_reference() {
        let k = 5;
        let mut db = MockDb::new(k);

        // Reference with high counts, no error extensions
        let ref_seq = "ACGTACGTAC";
        db.set_sequence(ref_seq, 1000);
        // No error k-mers set => all non-ref extensions return 0

        let ref_kmers: Vec<String> = ref_seq
            .as_bytes()
            .windows(k as usize)
            .map(|w| std::str::from_utf8(w).unwrap().to_string())
            .collect();

        let sample: Vec<&String> = ref_kmers.iter().collect();
        let error_rate = estimate_error_rate(&db, &sample);

        // With only reference extensions having counts, error rate should be very low
        // (some ref k-mers' forward extensions overlap with other ref k-mers,
        // so the "max" extension is nonzero and error sum may be tiny or zero)
        assert!(
            error_rate < 0.01,
            "Clean reference should have low error rate, got {}",
            error_rate
        );
    }

    // ---- Error rate: noisy reference (high error) ----

    #[test]
    fn test_error_rate_noisy_reference() {
        let k = 5;
        let mut db = MockDb::new(k);

        // Single reference k-mer at 1000x
        db.set("ACGTA", 1000);

        // Forward extensions: dominant = CGTAC(900), errors = CGTAT(50), CGTAG(30), CGTAA(20)
        db.set("CGTAC", 900);
        db.set("CGTAT", 50);
        db.set("CGTAG", 30);
        db.set("CGTAA", 20);

        let ref_kmers = vec!["ACGTA".to_string()];
        let sample: Vec<&String> = ref_kmers.iter().collect();
        let error_rate = estimate_error_rate(&db, &sample);

        // Error sum = 50 + 30 + 20 = 100 (the 3 non-dominant extensions)
        // Error rate = 100 / (3 * 1000) = 0.0333...
        let expected = 100.0 / 3000.0;
        assert!(
            (error_rate - expected).abs() < 1e-6,
            "Expected error rate ~{:.4}, got {:.4}",
            expected,
            error_rate
        );
    }

    // ---- Count threshold never below 2 ----

    #[test]
    fn test_count_threshold_never_below_2() {
        let k = 5;
        let mut db = MockDb::new(k);

        // Very low coverage: single k-mer at count 1
        db.set("ACGTA", 1);

        let ref_kmers = vec!["ACGTA".to_string()];
        let config = AdaptiveConfig::default();
        let thresholds = estimate_thresholds(&db, &ref_kmers, &config);

        assert!(
            thresholds.count_threshold >= 2,
            "Count threshold should never be below 2, got {}",
            thresholds.count_threshold
        );
    }

    #[test]
    fn test_count_threshold_never_below_2_zero_coverage() {
        let k = 5;
        let db = MockDb::new(k);

        // Zero coverage (no k-mers in DB)
        let ref_kmers = vec!["ACGTA".to_string()];
        let config = AdaptiveConfig::default();
        let thresholds = estimate_thresholds(&db, &ref_kmers, &config);

        assert!(
            thresholds.count_threshold >= 2,
            "Count threshold should never be below 2 even with zero coverage, got {}",
            thresholds.count_threshold
        );
    }

    // ---- Empty reference k-mers ----

    #[test]
    fn test_estimate_thresholds_empty_ref_kmers() {
        let k = 5;
        let db = MockDb::new(k);

        let ref_kmers: Vec<String> = vec![];
        let config = AdaptiveConfig::default();
        let thresholds = estimate_thresholds(&db, &ref_kmers, &config);

        assert_eq!(thresholds.median_coverage, 0.0);
        assert_eq!(thresholds.tier, CoverageTier::UltraLow);
        assert!(thresholds.count_threshold >= 2);
    }

    // ---- Sampling behavior ----

    #[test]
    fn test_estimate_thresholds_sampling_many_kmers() {
        let k = 5;
        let mut db = MockDb::new(k);

        // Create a large set of reference k-mers (more than default sample_size=100)
        let mut ref_kmers = Vec::new();
        let bases = ['A', 'C', 'G', 'T'];
        for &b1 in &bases {
            for &b2 in &bases {
                for &b3 in &bases {
                    for &b4 in &bases {
                        for &b5 in &bases {
                            let kmer = format!("{}{}{}{}{}", b1, b2, b3, b4, b5);
                            db.set(&kmer, 1500);
                            ref_kmers.push(kmer);
                        }
                    }
                }
            }
        }

        // 1024 k-mers total, all at 1500x
        let config = AdaptiveConfig {
            sample_size: 50,
            ..AdaptiveConfig::default()
        };
        let thresholds = estimate_thresholds(&db, &ref_kmers, &config);

        // Should correctly identify Standard tier
        assert_eq!(thresholds.tier, CoverageTier::Standard);
        assert!((thresholds.median_coverage - 1500.0).abs() < 1.0);
    }

    // ---- AdaptiveConfig default values ----

    #[test]
    fn test_adaptive_config_defaults() {
        let config = AdaptiveConfig::default();
        assert!((config.false_extension_rate - 0.01).abs() < f64::EPSILON);
        assert!((config.default_error_rate - 0.001).abs() < f64::EPSILON);
        assert_eq!(config.sample_size, 100);
    }

    // ---- Ratio thresholds per tier ----

    #[test]
    fn test_ratio_threshold_per_tier() {
        let k = 5;
        let ref_seq = "ACGTACGTAC";

        // Test each tier produces the correct ratio threshold
        let test_cases: &[(u64, CoverageTier, f64)] = &[
            (100, CoverageTier::UltraLow, 0.10),
            (1000, CoverageTier::Standard, 0.05),
            (5000, CoverageTier::HighDepth, 0.03),
            (20000, CoverageTier::UltraDeep, 0.01),
        ];

        for &(coverage, expected_tier, expected_ratio) in test_cases {
            let mut db = MockDb::new(k);
            db.set_sequence(ref_seq, coverage);

            let ref_kmers: Vec<String> = ref_seq
                .as_bytes()
                .windows(k as usize)
                .map(|w| std::str::from_utf8(w).unwrap().to_string())
                .collect();

            let config = AdaptiveConfig::default();
            let thresholds = estimate_thresholds(&db, &ref_kmers, &config);

            assert_eq!(
                thresholds.tier, expected_tier,
                "Coverage {} should yield tier {:?}",
                coverage, expected_tier
            );
            assert!(
                (thresholds.ratio_threshold - expected_ratio).abs() < f64::EPSILON,
                "Coverage {} (tier {:?}) should yield ratio {}, got {}",
                coverage,
                expected_tier,
                expected_ratio,
                thresholds.ratio_threshold
            );
        }
    }

    // ---- Poisson threshold computation ----

    #[test]
    fn test_poisson_threshold_increases_with_coverage() {
        let k = 5;
        let ref_seq = "ACGTACGTAC";
        let config = AdaptiveConfig::default();

        let mut prev_count = 0u32;
        for &coverage in &[100u64, 1000, 5000, 20000] {
            let mut db = MockDb::new(k);
            db.set_sequence(ref_seq, coverage);

            // Add small error extensions to ensure nonzero error rate
            db.set("CGTAT", (coverage / 200).max(1));
            db.set("CGTAG", (coverage / 200).max(1));
            db.set("CGTAA", (coverage / 200).max(1));

            let ref_kmers: Vec<String> = ref_seq
                .as_bytes()
                .windows(k as usize)
                .map(|w| std::str::from_utf8(w).unwrap().to_string())
                .collect();

            let thresholds = estimate_thresholds(&db, &ref_kmers, &config);

            assert!(
                thresholds.count_threshold >= prev_count,
                "Count threshold should increase with coverage: {} < {} at coverage {}",
                thresholds.count_threshold,
                prev_count,
                coverage
            );
            prev_count = thresholds.count_threshold;
        }
    }

    // ---- Median computation helper ----

    #[test]
    fn test_compute_median_odd_length() {
        assert!((compute_median(&[1, 2, 3]) - 2.0).abs() < f64::EPSILON);
        assert!((compute_median(&[10, 20, 30, 40, 50]) - 30.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_compute_median_even_length() {
        assert!((compute_median(&[1, 2, 3, 4]) - 2.5).abs() < f64::EPSILON);
        assert!((compute_median(&[10, 20]) - 15.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_compute_median_single_element() {
        assert!((compute_median(&[42]) - 42.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_compute_median_empty() {
        assert!((compute_median(&[]) - 0.0).abs() < f64::EPSILON);
    }
}
