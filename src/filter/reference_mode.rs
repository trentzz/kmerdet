/// Reference mode filtering: match by CHROM/POS/REF/ALT/TYPE.
///
/// This is the default filtering mode. Each detected variant is matched
/// against expected variants by genomic coordinates and alleles.

use crate::filter::ExpectedVariant;
use crate::variant::VariantCall;

/// Try to find a matching detected call for an expected variant
/// using reference-based matching (CHROM + POS + REF + ALT + TYPE).
pub fn find_match<'a>(
    _expected: &ExpectedVariant,
    _calls: &'a [VariantCall],
) -> Option<&'a VariantCall> {
    // TODO Phase 7: Implement reference-mode matching
    todo!("reference mode matching not yet implemented")
}
