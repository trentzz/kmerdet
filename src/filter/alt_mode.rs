/// Alt-sequence mode filtering: match by full ALT_SEQUENCE.
///
/// From kmtools dev branch. Instead of matching by genomic coordinates,
/// match by checking if the expected variant's alternative sequence appears
/// in any detected call's alt_sequence field.

use crate::filter::ExpectedVariant;
use crate::variant::VariantCall;

/// Try to find a matching detected call for an expected variant
/// using alt-sequence matching (full ALT_SEQUENCE comparison).
///
/// This is more permissive than reference mode — it checks whether the
/// expected variant's alt_allele appears anywhere within the call's
/// alt_sequence. This handles cases where coordinate-based representation
/// differs but the actual sequence is the same (common with INDELs).
pub fn find_match<'a>(
    expected: &ExpectedVariant,
    calls: &'a [VariantCall],
) -> Option<&'a VariantCall> {
    calls
        .iter()
        .find(|call| call.alt_sequence.contains(&expected.alt_allele))
}
