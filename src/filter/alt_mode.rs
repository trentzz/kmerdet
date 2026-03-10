/// Alt-sequence mode filtering: match by full ALT_SEQUENCE.
///
/// From kmtools dev branch. Instead of matching by genomic coordinates,
/// match by checking if the expected variant's alternative sequence appears
/// in any detected call's alt_sequence field.

use crate::filter::ExpectedVariant;
use crate::variant::VariantCall;

/// Try to find a matching detected call for an expected variant
/// using alt-sequence matching (full ALT_SEQUENCE comparison).
pub fn find_match<'a>(
    _expected: &ExpectedVariant,
    _calls: &'a [VariantCall],
) -> Option<&'a VariantCall> {
    // TODO Phase 7: Implement alt-mode matching
    todo!("alt mode matching not yet implemented")
}
