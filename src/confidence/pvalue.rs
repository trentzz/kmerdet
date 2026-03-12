/// Binomial p-value computation and Fisher's method for combining p-values
/// across multiple k-mers in a variant path.
///
/// The key idea: under the null hypothesis (no variant, only sequencing errors),
/// the count of a k-mer at any position follows a binomial distribution with
/// parameters n = total coverage and p = per-base error rate.  A variant k-mer
/// with count significantly above the error expectation yields a small p-value.
///
/// We compute per-k-mer p-values and combine them with Fisher's method,
/// applying an autocorrelation correction because overlapping k-mers from the
/// same reads are not independent.

use std::collections::HashSet;

use crate::jellyfish::KmerDatabase;
use crate::sequence::path::KmerPath;

// ────────────────────────────── Binomial survival function ──────────────────

/// Compute 1 - BinomCDF(k-1, n, p) = P(X >= k | n, p).
///
/// Uses the relationship to the regularized incomplete beta function:
///   P(X >= k) = I(p, k, n - k + 1)
///
/// For k == 0 the result is 1.0 (trivially true).
/// Falls back to direct summation when n is small (n <= 1000).
pub fn binomial_sf(k: u64, n: u64, p: f64) -> f64 {
    if k == 0 {
        return 1.0;
    }
    if k > n {
        return 0.0;
    }
    if p <= 0.0 {
        return 0.0;
    }
    if p >= 1.0 {
        return 1.0;
    }

    // For small n, direct summation of the upper tail is fast and stable.
    if n <= 1000 {
        return binomial_sf_direct(k, n, p);
    }

    // For large n use the incomplete beta relationship.
    regularized_incomplete_beta(p, k as f64, (n - k + 1) as f64)
}

/// Direct summation: P(X >= k) = 1 - sum_{i=0}^{k-1} C(n,i) p^i (1-p)^{n-i}
///
/// We compute in log-space to avoid overflow then accumulate.
fn binomial_sf_direct(k: u64, n: u64, p: f64) -> f64 {
    let ln_p = p.ln();
    let ln_q = (1.0 - p).ln();

    // Sum P(X = i) for i = 0 .. k-1  (the CDF below k).
    let mut cdf = 0.0f64;
    let mut ln_binom = 0.0f64; // ln(C(n, 0)) = 0

    for i in 0..k {
        let ln_pmf = ln_binom + (i as f64) * ln_p + ((n - i) as f64) * ln_q;
        cdf += ln_pmf.exp();
        if cdf >= 1.0 {
            return 0.0;
        }
        // Update ln(C(n, i+1)) = ln(C(n, i)) + ln(n - i) - ln(i + 1)
        ln_binom += ((n - i) as f64).ln() - ((i + 1) as f64).ln();
    }

    (1.0 - cdf).max(0.0)
}

// ────────────── Regularized incomplete beta via continued fraction ──────────

/// Evaluate the regularized incomplete beta function I(x, a, b)
/// using the Lentz continued-fraction algorithm (DLMF 8.17.22).
///
/// I(x, a, b) = x^a (1-x)^b / (a B(a,b)) * CF(x, a, b)
///
/// where B(a,b) is the beta function and CF is a continued fraction.
pub fn regularized_incomplete_beta(x: f64, a: f64, b: f64) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    if x >= 1.0 {
        return 1.0;
    }

    // Use the symmetry relation when x > (a+1)/(a+b+2) for better convergence.
    if x > (a + 1.0) / (a + b + 2.0) {
        return 1.0 - regularized_incomplete_beta(1.0 - x, b, a);
    }

    // Front factor: x^a * (1-x)^b / (a * B(a,b))
    // ln(front) = a*ln(x) + b*ln(1-x) - ln(a) - ln_beta(a, b)
    let ln_front = a * x.ln() + b * (1.0 - x).ln() - a.ln() - ln_beta(a, b);
    let front = ln_front.exp();

    // Evaluate the continued fraction using the modified Lentz method.
    let cf = beta_cf(x, a, b);

    (front * cf).clamp(0.0, 1.0)
}

/// Evaluate the continued fraction for the incomplete beta function.
///
/// Uses the modified Lentz algorithm with the standard even/odd coefficient
/// recurrence (Numerical Recipes, Section 6.4).
fn beta_cf(x: f64, a: f64, b: f64) -> f64 {
    const MAX_ITER: usize = 200;
    const EPS: f64 = 1e-14;
    const TINY: f64 = 1e-30;

    let mut c = 1.0f64;
    let mut d = 1.0 - (a + b) * x / (a + 1.0);
    if d.abs() < TINY {
        d = TINY;
    }
    d = 1.0 / d;
    let mut f = d;

    for m in 1..=MAX_ITER {
        let m_f64 = m as f64;

        // Even step: d_{2m}
        let numerator_even =
            m_f64 * (b - m_f64) * x / ((a + 2.0 * m_f64 - 1.0) * (a + 2.0 * m_f64));

        d = 1.0 + numerator_even * d;
        if d.abs() < TINY {
            d = TINY;
        }
        c = 1.0 + numerator_even / c;
        if c.abs() < TINY {
            c = TINY;
        }
        d = 1.0 / d;
        f *= d * c;

        // Odd step: d_{2m+1}
        let numerator_odd = -((a + m_f64) * (a + b + m_f64) * x)
            / ((a + 2.0 * m_f64) * (a + 2.0 * m_f64 + 1.0));

        d = 1.0 + numerator_odd * d;
        if d.abs() < TINY {
            d = TINY;
        }
        c = 1.0 + numerator_odd / c;
        if c.abs() < TINY {
            c = TINY;
        }
        d = 1.0 / d;
        let delta = d * c;
        f *= delta;

        if (delta - 1.0).abs() < EPS {
            break;
        }
    }

    f
}

/// Natural log of the beta function: ln B(a, b) = ln Gamma(a) + ln Gamma(b) - ln Gamma(a+b).
fn ln_beta(a: f64, b: f64) -> f64 {
    ln_gamma(a) + ln_gamma(b) - ln_gamma(a + b)
}

/// Stirling-series approximation of ln(Gamma(x)) for x > 0.
///
/// Uses the Lanczos approximation (g = 7, n = 9) which gives ~15 digits of
/// accuracy for all x > 0.5.  For x < 0.5 the reflection formula is applied.
fn ln_gamma(x: f64) -> f64 {
    // Lanczos coefficients (g = 7, n = 9)
    const G: f64 = 7.0;
    const COEFFS: [f64; 9] = [
        0.999_999_999_999_809_93,
        676.520_368_121_885_1,
        -1259.139_216_722_403,
        771.323_428_777_653_1,
        -176.615_029_162_140_6,
        12.507_343_278_686_905,
        -0.138_571_095_265_720_12,
        9.984_369_578_019_572e-6,
        1.505_632_735_149_311_6e-7,
    ];

    if x < 0.5 {
        // Reflection: Gamma(x) Gamma(1-x) = pi / sin(pi x)
        let reflect = (std::f64::consts::PI * x).sin().abs().ln();
        return std::f64::consts::PI.ln() - reflect - ln_gamma(1.0 - x);
    }

    let x = x - 1.0;
    let mut sum = COEFFS[0];
    for (i, &c) in COEFFS.iter().enumerate().skip(1) {
        sum += c / (x + i as f64);
    }
    let t = x + G + 0.5;
    0.5 * (2.0 * std::f64::consts::PI).ln() + (t.ln() * (x + 0.5)) - t + sum.ln()
}

// ──────────────────── Per-k-mer binomial p-value ───────────────────────────

/// Compute the p-value for a single k-mer observation.
///
/// Uses the binomial survival function: P(X >= c | n, p)
/// where c = variant_count, n = total_coverage, p = error_rate.
///
/// A small p-value means the observed variant count is unlikely under the
/// null hypothesis of random sequencing errors alone.
pub fn binomial_pvalue(variant_count: u64, total_coverage: u64, error_rate: f64) -> f64 {
    if variant_count == 0 {
        return 1.0;
    }
    if total_coverage == 0 {
        return 1.0;
    }
    binomial_sf(variant_count, total_coverage, error_rate)
}

// ──────────────────── Fisher's combined p-value ────────────────────────────

/// Combine p-values across multiple k-mers using Fisher's method.
///
/// chi2 = -2 * sum(ln(p_i))
///
/// Under the null, chi2 ~ chi-squared(2k).  We compute the survival function
/// of that chi-squared distribution.
///
/// Because overlapping k-mers from the same reads are correlated, we apply an
/// autocorrelation correction: effective_k = k * kmer_length / read_length
/// (clamped to at least 1).
pub fn fisher_combined_pvalue(
    pvalues: &[f64],
    kmer_length: usize,
    read_length: usize,
) -> f64 {
    if pvalues.is_empty() {
        return 1.0;
    }

    // Clamp each p-value away from exactly 0 (which would give -inf in log).
    let clamped: Vec<f64> = pvalues.iter().map(|&p| p.max(1e-300)).collect();

    let chi2: f64 = -2.0 * clamped.iter().map(|p| p.ln()).sum::<f64>();

    // Autocorrelation correction: reduce effective degrees of freedom.
    let raw_k = pvalues.len() as f64;
    let correction = if read_length > 0 {
        (kmer_length as f64 / read_length as f64).min(1.0)
    } else {
        1.0
    };
    let effective_k = (raw_k * correction).max(1.0);

    // Degrees of freedom for Fisher's statistic.
    let df = 2.0 * effective_k;

    // P(chi2 > observed) = 1 - CDF_chisq(chi2, df)
    // = regularized upper incomplete gamma: Q(df/2, chi2/2)
    chi_squared_sf(chi2, df)
}

/// Survival function of the chi-squared distribution: P(X > x | df).
///
/// Equal to the regularized upper incomplete gamma function Q(df/2, x/2).
fn chi_squared_sf(x: f64, df: f64) -> f64 {
    if x <= 0.0 {
        return 1.0;
    }
    regularized_upper_gamma(df / 2.0, x / 2.0)
}

/// Regularized upper incomplete gamma function Q(a, x) = 1 - P(a, x).
///
/// Uses the continued fraction representation for x >= a + 1
/// and the series representation otherwise.
fn regularized_upper_gamma(a: f64, x: f64) -> f64 {
    if x < 0.0 {
        return 1.0;
    }
    if x == 0.0 {
        return 1.0;
    }

    if x < a + 1.0 {
        // Use the series for P(a, x) and return 1 - P.
        1.0 - gamma_series_p(a, x)
    } else {
        // Use the continued fraction for Q(a, x) directly.
        gamma_cf_q(a, x)
    }
}

/// Series representation of the regularized lower incomplete gamma P(a, x).
fn gamma_series_p(a: f64, x: f64) -> f64 {
    const MAX_ITER: usize = 200;
    const EPS: f64 = 1e-14;

    let ln_prefix = a * x.ln() - x - ln_gamma(a + 1.0);
    let prefix = ln_prefix.exp();

    let mut sum = 1.0;
    let mut term = 1.0;
    for n in 1..=MAX_ITER {
        term *= x / (a + n as f64);
        sum += term;
        if term.abs() < EPS * sum.abs() {
            break;
        }
    }

    (prefix * sum).clamp(0.0, 1.0)
}

/// Continued-fraction representation of Q(a, x) (Numerical Recipes 6.2.7).
fn gamma_cf_q(a: f64, x: f64) -> f64 {
    const MAX_ITER: usize = 200;
    const EPS: f64 = 1e-14;
    const TINY: f64 = 1e-30;

    // Front factor: x^a e^{-x} / Gamma(a)
    let ln_prefix = a * x.ln() - x - ln_gamma(a);

    let mut b = x + 1.0 - a;
    let mut c = 1.0 / TINY;
    let mut d = 1.0 / b;
    let mut h = d;

    for i in 1..=MAX_ITER {
        let an = -(i as f64) * (i as f64 - a);
        b += 2.0;
        d = an * d + b;
        if d.abs() < TINY {
            d = TINY;
        }
        c = b + an / c;
        if c.abs() < TINY {
            c = TINY;
        }
        d = 1.0 / d;
        let delta = d * c;
        h *= delta;
        if (delta - 1.0).abs() < EPS {
            break;
        }
    }

    (ln_prefix.exp() * h).clamp(0.0, 1.0)
}

// ──────────────────── Error rate estimation ─────────────────────────────────

const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];

/// Estimate per-sample error rate from reference k-mer non-reference extensions.
///
/// For each reference k-mer, query all 4 right extensions.  The extensions that
/// do NOT match the next reference base represent sequencing errors.
///
/// error_rate = sum(non_ref_extension_counts) / (3 * sum(ref_kmer_coverage))
///
/// Returns a default of 0.001 if insufficient data (< 10 reference k-mers with
/// non-zero coverage).
pub fn estimate_error_rate(
    db: &dyn KmerDatabase,
    ref_kmers: &[String],
) -> f64 {
    const DEFAULT_ERROR_RATE: f64 = 0.001;
    const MIN_KMERS: usize = 10;

    if ref_kmers.len() < 2 {
        return DEFAULT_ERROR_RATE;
    }

    let mut total_error_counts: u64 = 0;
    let mut total_coverage: u64 = 0;
    let mut usable_kmers: usize = 0;

    // For each ref k-mer except the last (since we need to know the next base),
    // look at all 4 right extensions.
    for i in 0..ref_kmers.len() - 1 {
        let kmer = &ref_kmers[i];
        let k = kmer.len();
        let suffix = &kmer[1..]; // last k-1 bases

        // Determine what the reference next base is.
        let next_kmer = &ref_kmers[i + 1];
        let ref_next_base = next_kmer.as_bytes().last().copied();

        let mut kmer_total: u64 = 0;
        let mut kmer_error: u64 = 0;

        for &base in &BASES {
            let mut child = String::with_capacity(k);
            child.push_str(suffix);
            child.push(base as char);
            let count = db.query(&child);
            kmer_total += count;

            if Some(base) != ref_next_base.map(|b| b.to_ascii_uppercase()) {
                kmer_error += count;
            }
        }

        if kmer_total > 0 {
            total_error_counts += kmer_error;
            total_coverage += kmer_total;
            usable_kmers += 1;
        }
    }

    if usable_kmers < MIN_KMERS || total_coverage == 0 {
        return DEFAULT_ERROR_RATE;
    }

    // The denominator accounts for the fact that we expect 3 out of 4 extensions
    // to be errors: error_rate per base = error_counts / (3 * ref_coverage).
    // But total_coverage already includes all 4 bases.  A cleaner formulation:
    //   error_rate = error_counts / total_coverage  (per-extension probability of error)
    // Then the per-base error rate would be error_counts / total_coverage.
    // We use the simpler: non_ref / total since that directly estimates the
    // probability that a random extension is erroneous.
    let rate = total_error_counts as f64 / total_coverage as f64;

    // Clamp to a reasonable range.
    rate.clamp(1e-6, 0.1)
}

// ──────────────────── Variant-level p-value ─────────────────────────────────

/// Compute the combined p-value for a variant path.
///
/// For each k-mer in the alt path that is NOT in the ref path, we compute a
/// per-k-mer binomial p-value (is this count too high to be explained by
/// sequencing error alone?).  Then we combine with Fisher's method.
///
/// Returns 1.0 if no variant-specific k-mers are found (e.g., reference call).
pub fn compute_variant_pvalue(
    db: &dyn KmerDatabase,
    alt_path: &KmerPath,
    ref_path: &KmerPath,
    error_rate: f64,
    read_length: usize,
) -> f64 {
    let ref_set: HashSet<String> = ref_path
        .kmers
        .iter()
        .map(|k| k.to_uppercase())
        .collect();

    let kmer_length = if alt_path.kmers.is_empty() {
        return 1.0;
    } else {
        alt_path.kmers[0].len()
    };

    // Estimate coverage from reference k-mers (median of ref path counts).
    let ref_counts: Vec<u64> = ref_path
        .kmers
        .iter()
        .map(|k| db.query(&k.to_uppercase()))
        .filter(|&c| c > 0)
        .collect();

    let median_coverage = if ref_counts.is_empty() {
        return 1.0;
    } else {
        let mut sorted = ref_counts.clone();
        sorted.sort_unstable();
        sorted[sorted.len() / 2]
    };

    // Collect per-k-mer p-values for variant-specific k-mers.
    let mut pvalues = Vec::new();

    for kmer in &alt_path.kmers {
        let upper = kmer.to_uppercase();
        if ref_set.contains(&upper) {
            continue; // Skip k-mers shared with the reference
        }

        let count = db.query(&upper);
        if count == 0 {
            // Zero-count variant k-mers are suspicious but don't contribute
            // a meaningful p-value (P(X >= 0) = 1).
            continue;
        }

        let p = binomial_pvalue(count, median_coverage, error_rate);
        pvalues.push(p);
    }

    if pvalues.is_empty() {
        return 1.0;
    }

    fisher_combined_pvalue(&pvalues, kmer_length, read_length)
}

// ──────────────────── Tests ────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── binomial_sf ──

    #[test]
    fn test_binomial_sf_zero_k() {
        // P(X >= 0) = 1.0 for any n, p
        assert_eq!(binomial_sf(0, 100, 0.5), 1.0);
    }

    #[test]
    fn test_binomial_sf_k_greater_than_n() {
        assert_eq!(binomial_sf(101, 100, 0.5), 0.0);
    }

    #[test]
    fn test_binomial_sf_p_zero() {
        assert_eq!(binomial_sf(5, 100, 0.0), 0.0);
    }

    #[test]
    fn test_binomial_sf_p_one() {
        assert_eq!(binomial_sf(5, 100, 1.0), 1.0);
    }

    #[test]
    fn test_binomial_sf_known_value_fair_coin() {
        // P(X >= 50 | n=100, p=0.5) should be ~0.5398
        let p = binomial_sf(50, 100, 0.5);
        assert!((p - 0.5398).abs() < 0.01, "got {}", p);
    }

    #[test]
    fn test_binomial_sf_rare_event() {
        // 5 successes out of 1000 with error rate 0.001 should be a small p-value
        // but not astronomically small.  E[X] = 1, so 5 is unusual.
        let p = binomial_sf(5, 1000, 0.001);
        assert!(p < 0.01, "expected p < 0.01, got {}", p);
        assert!(p > 1e-6, "expected p > 1e-6, got {}", p);
    }

    #[test]
    fn test_binomial_sf_large_n() {
        // Test with n > 1000 (triggers incomplete beta path).
        let p = binomial_sf(50, 5000, 0.001);
        assert!(p < 1e-10, "expected very small p-value, got {}", p);
    }

    #[test]
    fn test_binomial_sf_n_equals_k() {
        // P(X >= n | n, p) = p^n
        let p = binomial_sf(10, 10, 0.5);
        let expected = 0.5f64.powi(10);
        assert!((p - expected).abs() < 1e-10, "got {}", p);
    }

    // ── binomial_pvalue ──

    #[test]
    fn test_binomial_pvalue_zero_count() {
        assert_eq!(binomial_pvalue(0, 1000, 0.001), 1.0);
    }

    #[test]
    fn test_binomial_pvalue_zero_coverage() {
        assert_eq!(binomial_pvalue(5, 0, 0.001), 1.0);
    }

    #[test]
    fn test_binomial_pvalue_typical_variant() {
        // 10 counts at 1000x coverage with 0.001 error rate.
        // E[errors] = 1, so 10 is very significant.
        let p = binomial_pvalue(10, 1000, 0.001);
        assert!(p < 1e-5, "expected very significant, got {}", p);
    }

    // ── fisher_combined_pvalue ──

    #[test]
    fn test_fisher_single_pvalue() {
        // With one p-value, Fisher's method should give a result close to the input.
        let p = fisher_combined_pvalue(&[0.01], 31, 150);
        // With correction factor 31/150 ~= 0.207, effective_k ~= 0.207 (clamped to 1)
        // So df = 2, chi2 = -2*ln(0.01) = 9.21
        // P(chi2(2) > 9.21) ~= 0.01
        assert!(p < 0.02, "got {}", p);
        assert!(p > 0.005, "got {}", p);
    }

    #[test]
    fn test_fisher_empty() {
        assert_eq!(fisher_combined_pvalue(&[], 31, 150), 1.0);
    }

    #[test]
    fn test_fisher_all_significant() {
        // Multiple very significant p-values should combine to even more significant.
        let pvals = vec![0.001, 0.001, 0.001, 0.001, 0.001];
        let combined = fisher_combined_pvalue(&pvals, 31, 150);
        assert!(combined < 0.001, "got {}", combined);
    }

    #[test]
    fn test_fisher_all_nonsignificant() {
        let pvals = vec![0.5, 0.5, 0.5];
        let combined = fisher_combined_pvalue(&pvals, 31, 150);
        // Should not be significant.
        assert!(combined > 0.05, "got {}", combined);
    }

    // ── estimate_error_rate ──

    #[test]
    fn test_estimate_error_rate_insufficient_data() {
        // With only 1 k-mer, should return default.
        let db = MockDbForTest::new();
        let rate = estimate_error_rate(&db, &["ACGTA".to_string()]);
        assert!((rate - 0.001).abs() < 1e-10);
    }

    #[test]
    fn test_estimate_error_rate_clean_reference() {
        // Simulate a clean reference with a real overlapping k-mer chain.
        // Use a non-repeating sequence to avoid hash collisions in the mock.
        let mut db = MockDbForTest::new();

        let sequence = "ACGTACGATCAGTACGATCA";
        let k = 5usize;
        let ref_kmers: Vec<String> = sequence
            .as_bytes()
            .windows(k)
            .map(|w| std::str::from_utf8(w).unwrap().to_string())
            .collect();

        // For each consecutive pair of ref k-mers, the "correct" right extension
        // is the last base of ref_kmers[i+1].  All other extensions are errors.
        for i in 0..ref_kmers.len() - 1 {
            let kmer = &ref_kmers[i];
            let suffix = &kmer[1..]; // last k-1 bases
            let next_base = ref_kmers[i + 1].as_bytes()[k - 1];

            for &base in &[b'A', b'C', b'G', b'T'] {
                let mut child = String::from(suffix);
                child.push(base as char);
                if base == next_base {
                    db.set(&child, 1000);
                } else {
                    db.set(&child, 1); // Error count
                }
            }
        }

        let rate = estimate_error_rate(&db, &ref_kmers);
        // Expected: ~3 errors per ~1003 total per k-mer => rate ~0.003
        assert!(
            rate < 0.01,
            "expected low error rate, got {} (ref_kmers={})",
            rate,
            ref_kmers.len()
        );
        assert!(rate > 0.0, "expected positive error rate, got {}", rate);
    }

    // ── compute_variant_pvalue ──

    #[test]
    fn test_compute_variant_pvalue_reference_call() {
        // If alt path == ref path, no variant-specific k-mers => p = 1.0
        let db = MockDbForTest::new();
        let ref_path = KmerPath {
            kmers: vec!["ACGTA".to_string(), "CGTAC".to_string()],
            is_reference: true,
        };
        let alt_path = KmerPath {
            kmers: vec!["ACGTA".to_string(), "CGTAC".to_string()],
            is_reference: false,
        };
        let p = compute_variant_pvalue(&db, &alt_path, &ref_path, 0.001, 150);
        assert_eq!(p, 1.0);
    }

    #[test]
    fn test_compute_variant_pvalue_with_variant_kmers() {
        let mut db = MockDbForTest::new();

        // Reference k-mers have high counts.
        db.set("ACGTA", 1000);
        db.set("CGTAC", 1000);
        db.set("GTACG", 1000);

        // Variant k-mers have moderate counts (real variant).
        db.set("CGTTT", 50);
        db.set("GTTTC", 50);

        let ref_path = KmerPath {
            kmers: vec![
                "ACGTA".to_string(),
                "CGTAC".to_string(),
                "GTACG".to_string(),
            ],
            is_reference: true,
        };
        let alt_path = KmerPath {
            kmers: vec![
                "ACGTA".to_string(),
                "CGTTT".to_string(),
                "GTTTC".to_string(),
            ],
            is_reference: false,
        };

        let p = compute_variant_pvalue(&db, &alt_path, &ref_path, 0.001, 150);
        // 50 variant counts at 1000x coverage with 0.001 error rate should be
        // extremely significant.
        assert!(p < 1e-10, "expected very small p-value, got {}", p);
    }

    #[test]
    fn test_compute_variant_pvalue_empty_alt_path() {
        let db = MockDbForTest::new();
        let ref_path = KmerPath {
            kmers: vec!["ACGTA".to_string()],
            is_reference: true,
        };
        let alt_path = KmerPath {
            kmers: vec![],
            is_reference: false,
        };
        let p = compute_variant_pvalue(&db, &alt_path, &ref_path, 0.001, 150);
        assert_eq!(p, 1.0);
    }

    // ── Regularized incomplete beta ──

    #[test]
    fn test_incomplete_beta_zero() {
        assert_eq!(regularized_incomplete_beta(0.0, 1.0, 1.0), 0.0);
    }

    #[test]
    fn test_incomplete_beta_one() {
        assert_eq!(regularized_incomplete_beta(1.0, 1.0, 1.0), 1.0);
    }

    #[test]
    fn test_incomplete_beta_uniform() {
        // For a=1, b=1 (uniform distribution), I(x, 1, 1) = x.
        let val = regularized_incomplete_beta(0.3, 1.0, 1.0);
        assert!((val - 0.3).abs() < 1e-10, "got {}", val);
    }

    #[test]
    fn test_incomplete_beta_half() {
        // I(0.5, 1, 1) = 0.5
        let val = regularized_incomplete_beta(0.5, 1.0, 1.0);
        assert!((val - 0.5).abs() < 1e-10, "got {}", val);
    }

    // ── Helper mock DB for tests ──

    struct MockDbForTest {
        counts: std::collections::HashMap<String, u64>,
    }

    impl MockDbForTest {
        fn new() -> Self {
            Self {
                counts: std::collections::HashMap::new(),
            }
        }

        fn set(&mut self, kmer: &str, count: u64) {
            self.counts.insert(kmer.to_uppercase(), count);
        }
    }

    impl KmerDatabase for MockDbForTest {
        fn query(&self, kmer: &str) -> u64 {
            *self.counts.get(&kmer.to_uppercase()).unwrap_or(&0)
        }

        fn kmer_length(&self) -> u8 {
            5
        }
    }
}
