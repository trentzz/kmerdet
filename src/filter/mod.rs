pub mod alt_mode;
pub mod conditions;
pub mod reference_mode;
pub mod report;

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
///
/// For each expected variant, attempts to find a matching detected call
/// (using either reference mode or alt-sequence mode depending on config),
/// then checks that the match passes all filter conditions (coverage, VAF,
/// expression, variant type).
pub fn filter_results(
    calls: &[VariantCall],
    expected: &[ExpectedVariant],
    config: &FilterConfig,
) -> Result<Vec<FilterResult>> {
    let mut results = Vec::new();

    for ev in expected {
        // Step 1: Find matching call using the configured mode
        let matched_call = if config.use_alt {
            alt_mode::find_match(ev, calls)
        } else {
            reference_mode::find_match(ev, calls)
        };

        // Step 2: Check if the match passes all filter conditions
        let (found, call, notes) = match matched_call {
            Some(call) => {
                let mut fails = Vec::new();
                if !conditions::passes_coverage(call, config.min_coverage) {
                    fails.push(format!(
                        "coverage {} < {}",
                        call.min_coverage, config.min_coverage
                    ));
                }
                if !conditions::passes_vaf(call, config.min_vaf) {
                    fails.push(format!("VAF {:.6} < {}", call.rvaf, config.min_vaf));
                }
                if !conditions::passes_expression(call, config.min_expression) {
                    fails.push(format!(
                        "expression {:.2} < {}",
                        call.expression, config.min_expression
                    ));
                }
                if !conditions::passes_type(call, &config.types) {
                    fails.push(format!(
                        "type {} not in allowed types",
                        call.variant_type
                    ));
                }

                if fails.is_empty() {
                    ("Found".to_string(), Some(call), "PASS".to_string())
                } else {
                    (
                        "Not Found".to_string(),
                        Some(call),
                        fails.join("; "),
                    )
                }
            }
            None => (
                "Not Found".to_string(),
                None,
                "No matching variant detected".to_string(),
            ),
        };

        // Step 3: Build FilterResult
        results.push(FilterResult {
            sample: call
                .map(|c| c.sample.clone())
                .unwrap_or_default(),
            chrom: ev.chrom.clone(),
            pos: ev.pos,
            ref_allele: ev.ref_allele.clone(),
            alt_allele: ev.alt_allele.clone(),
            variant_type: ev.variant_type.clone(),
            found,
            filter_notes: notes,
            kmer_vaf: call.map(|c| c.rvaf),
            kmer_min_coverage: call.map(|c| c.min_coverage),
            kmer_expression: call.map(|c| c.expression),
            ref_sequence: call.map(|c| c.ref_sequence.clone()),
            variant_sequence: call.map(|c| c.alt_sequence.clone()),
        });
    }

    Ok(results)
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
