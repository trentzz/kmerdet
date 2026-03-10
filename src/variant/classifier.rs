/// Variant classification from reference/alternative path comparison.
///
/// Implements `diff_path_without_overlap` from km:
/// 1. Left-to-right scan for first divergence point (i)
/// 2. Right-to-left scan for last divergence point (j_ref, j_var)
/// 3. Classify based on relative positions and lengths

use super::VariantType;
use crate::sequence::path::KmerPath;

/// Result of comparing a reference path to an alternative path.
#[derive(Debug)]
pub struct Classification {
    pub variant_type: VariantType,
    /// Human-readable variant name (e.g., "41:A/T:41" for SNV).
    pub variant_name: String,
    /// Reference allele at the variant position.
    pub ref_allele: String,
    /// Alternative allele at the variant position.
    pub alt_allele: String,
    /// Start position of the variant within the target.
    pub start: usize,
    /// End position of the variant within the target.
    pub end: usize,
}

/// Classify the variant by comparing reference and alternative paths.
pub fn classify(_ref_path: &KmerPath, _alt_path: &KmerPath, _k: u8) -> Classification {
    // TODO Phase 4: Implement diff_path_without_overlap algorithm
    todo!("variant classifier not yet implemented")
}
