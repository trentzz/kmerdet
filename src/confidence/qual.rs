/// Phred-scaled quality score conversion from p-values.
///
/// QUAL = -10 * log10(p_value)
///
/// This is the standard quality encoding used in VCF files and throughout
/// bioinformatics.  A QUAL of 20 corresponds to p = 0.01, 30 to p = 0.001, etc.

/// Convert a p-value to a Phred-scaled quality score.
///
/// QUAL = -10 * log10(p_value), capped at 999.
///
/// Special cases:
/// - p_value <= 0 => 999 (maximum quality)
/// - p_value >= 1 => 0 (no evidence)
pub fn phred_qual(p_value: f64) -> f64 {
    if p_value <= 0.0 {
        return 999.0;
    }
    if p_value >= 1.0 {
        return 0.0;
    }
    let qual = -10.0 * p_value.log10();
    qual.min(999.0)
}

/// Convert a Phred quality score back to a p-value.
///
/// p_value = 10^(-QUAL / 10)
pub fn pvalue_from_qual(qual: f64) -> f64 {
    if qual <= 0.0 {
        return 1.0;
    }
    10.0f64.powf(-qual / 10.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_phred_qual_standard_values() {
        // p = 0.1 -> QUAL = 10
        let q = phred_qual(0.1);
        assert!((q - 10.0).abs() < 1e-10, "got {}", q);

        // p = 0.01 -> QUAL = 20
        let q = phred_qual(0.01);
        assert!((q - 20.0).abs() < 1e-10, "got {}", q);

        // p = 0.001 -> QUAL = 30
        let q = phred_qual(0.001);
        assert!((q - 30.0).abs() < 1e-10, "got {}", q);

        // p = 0.0001 -> QUAL = 40
        let q = phred_qual(0.0001);
        assert!((q - 40.0).abs() < 1e-10, "got {}", q);
    }

    #[test]
    fn test_phred_qual_zero_pvalue() {
        assert_eq!(phred_qual(0.0), 999.0);
    }

    #[test]
    fn test_phred_qual_negative_pvalue() {
        assert_eq!(phred_qual(-0.1), 999.0);
    }

    #[test]
    fn test_phred_qual_one() {
        assert_eq!(phred_qual(1.0), 0.0);
    }

    #[test]
    fn test_phred_qual_greater_than_one() {
        assert_eq!(phred_qual(1.5), 0.0);
    }

    #[test]
    fn test_phred_qual_cap_at_999() {
        // Very small p-value should be capped.
        let q = phred_qual(1e-200);
        assert_eq!(q, 999.0);
    }

    #[test]
    fn test_roundtrip() {
        let original = 0.005;
        let qual = phred_qual(original);
        let recovered = pvalue_from_qual(qual);
        assert!(
            (original - recovered).abs() < 1e-12,
            "original={}, recovered={}",
            original,
            recovered
        );
    }

    #[test]
    fn test_pvalue_from_qual_zero() {
        assert_eq!(pvalue_from_qual(0.0), 1.0);
    }

    #[test]
    fn test_pvalue_from_qual_standard() {
        let p = pvalue_from_qual(30.0);
        assert!((p - 0.001).abs() < 1e-12, "got {}", p);
    }
}
