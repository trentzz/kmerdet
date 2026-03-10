/// Reference mode filtering: match by CHROM/POS/REF/ALT/TYPE.
///
/// This is the default filtering mode. Each detected variant is matched
/// against expected variants by genomic coordinates and alleles.

use crate::filter::ExpectedVariant;
use crate::variant::VariantCall;

/// Normalize a chromosome name by stripping the "chr" prefix if present.
fn normalize_chrom(chrom: &str) -> &str {
    chrom.strip_prefix("chr").unwrap_or(chrom)
}

/// Try to find a matching detected call for an expected variant
/// using reference-based matching (CHROM + POS + REF + ALT + TYPE).
pub fn find_match<'a>(
    expected: &ExpectedVariant,
    calls: &'a [VariantCall],
) -> Option<&'a VariantCall> {
    let expected_chrom_norm = normalize_chrom(&expected.chrom);

    calls.iter().find(|call| {
        // Match chromosome with normalization (chr1 == 1, chrX == X)
        let chrom_match = call
            .chrom
            .as_deref()
            .map(|c| normalize_chrom(c) == expected_chrom_norm)
            .unwrap_or(false);

        let pos_match = call.pos == Some(expected.pos);
        let ref_match = call.ref_allele.as_deref() == Some(expected.ref_allele.as_str());
        let alt_match = call.alt_allele.as_deref() == Some(expected.alt_allele.as_str());

        // Case-insensitive type comparison
        let type_match = call.variant_type.to_string().to_lowercase()
            == expected.variant_type.to_lowercase();

        chrom_match && pos_match && ref_match && alt_match && type_match
    })
}
