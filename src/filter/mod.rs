pub mod alt_mode;
pub mod conditions;
pub mod reference_mode;

use anyhow::Result;

use crate::variant::VariantCall;

/// Result of filtering a variant call against expected variants.
#[derive(Debug, Clone, serde::Serialize)]
pub struct FilterResult {
    pub sample: String,
    pub chrom: String,
    pub pos: u64,
    pub ref_allele: String,
    pub alt_allele: String,
    pub variant_type: String,
    /// "Found" or "Not Found".
    pub found: String,
    /// Explanation of why the variant was found/filtered.
    pub filter_notes: String,
    /// rVAF from detection (if found).
    pub kmer_vaf: Option<f64>,
    /// Min coverage from detection (if found).
    pub kmer_min_coverage: Option<u64>,
    /// Expression from detection (if found).
    pub kmer_expression: Option<f64>,
    /// Full reference sequence (if found).
    pub ref_sequence: Option<String>,
    /// Full variant sequence (if found).
    pub variant_sequence: Option<String>,
}

/// Configuration for the filter engine.
#[derive(Debug, Clone)]
pub struct FilterConfig {
    pub min_coverage: u32,
    pub min_vaf: f64,
    pub min_expression: f64,
    pub use_alt: bool,
    pub types: Vec<String>,
}

/// Run tumor-informed filtering on detection results against expected variants.
pub fn filter_results(
    _calls: &[VariantCall],
    _expected: &[ExpectedVariant],
    _config: &FilterConfig,
) -> Result<Vec<FilterResult>> {
    // TODO Phase 7: Implement filtering engine
    todo!("filter engine not yet implemented")
}

/// An expected variant from the reference/targets file.
#[derive(Debug, Clone)]
pub struct ExpectedVariant {
    pub chrom: String,
    pub pos: u64,
    pub ref_allele: String,
    pub alt_allele: String,
    pub variant_type: String,
}
