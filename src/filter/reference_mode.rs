/// Reference mode filtering: match by CHROM/POS/REF/ALT/TYPE.
///
/// This is the default filtering mode. Each detected variant is matched
/// against expected variants by genomic coordinates and alleles. When exact
/// matching fails for INDELs, INDEL normalization (left-alignment + trimming)
/// is used as a fallback to handle equivalent representations.

use crate::filter::ExpectedVariant;
use crate::variant::normalize::normalize_indel;
use crate::variant::VariantCall;

/// Normalize a chromosome name by stripping the "chr" prefix if present.
fn normalize_chrom(chrom: &str) -> &str {
    chrom.strip_prefix("chr").unwrap_or(chrom)
}

/// Check if a variant type string indicates an INDEL (insertion, deletion,
/// or complex indel) where coordinate normalization may be needed.
fn is_indel_type(variant_type: &str) -> bool {
    let lower = variant_type.to_lowercase();
    lower == "insertion" || lower == "deletion" || lower == "indel"
}

/// Try to find a matching detected call for an expected variant
/// using reference-based matching (CHROM + POS + REF + ALT + TYPE).
///
/// First attempts exact matching on all fields. If that fails and the
/// expected variant is an INDEL type, falls back to normalization-based
/// matching where both variants are left-aligned and trimmed before
/// comparison.
pub fn find_match<'a>(
    expected: &ExpectedVariant,
    calls: &'a [VariantCall],
) -> Option<&'a VariantCall> {
    find_match_with_context(expected, calls, None)
}

/// Like [`find_match`] but accepts an optional reference context for
/// INDEL left-alignment during normalization-based fallback matching.
pub fn find_match_with_context<'a>(
    expected: &ExpectedVariant,
    calls: &'a [VariantCall],
    ref_context: Option<&str>,
) -> Option<&'a VariantCall> {
    let expected_chrom_norm = normalize_chrom(&expected.chrom);

    // First pass: exact matching on all fields.
    let exact = calls.iter().find(|call| {
        let chrom_match = call
            .chrom
            .as_deref()
            .map(|c| normalize_chrom(c) == expected_chrom_norm)
            .unwrap_or(false);

        let pos_match = call.pos == Some(expected.pos);
        let ref_match = call.ref_allele.as_deref() == Some(expected.ref_allele.as_str());
        let alt_match = call.alt_allele.as_deref() == Some(expected.alt_allele.as_str());

        let type_match = call.variant_type.to_string().to_lowercase()
            == expected.variant_type.to_lowercase();

        chrom_match && pos_match && ref_match && alt_match && type_match
    });

    if exact.is_some() {
        return exact;
    }

    // Second pass: normalization-based fallback for INDELs.
    // Only attempt this if the expected variant is an INDEL type.
    if !is_indel_type(&expected.variant_type) {
        return None;
    }

    let expected_norm = normalize_indel(
        &expected.chrom,
        expected.pos,
        &expected.ref_allele,
        &expected.alt_allele,
        ref_context,
    );
    let expected_norm_chrom = normalize_chrom(&expected_norm.chrom);

    calls.iter().find(|call| {
        // Type must still match (case-insensitive).
        let type_match = call.variant_type.to_string().to_lowercase()
            == expected.variant_type.to_lowercase();
        if !type_match {
            return false;
        }

        // Chromosome must match.
        let call_chrom = match call.chrom.as_deref() {
            Some(c) => c,
            None => return false,
        };
        if normalize_chrom(call_chrom) != expected_norm_chrom {
            return false;
        }

        // Normalize the call's alleles and compare.
        let call_pos = match call.pos {
            Some(p) => p,
            None => return false,
        };
        let call_ref = match call.ref_allele.as_deref() {
            Some(r) => r,
            None => return false,
        };
        let call_alt = match call.alt_allele.as_deref() {
            Some(a) => a,
            None => return false,
        };

        let call_norm = normalize_indel(call_chrom, call_pos, call_ref, call_alt, ref_context);

        call_norm.pos == expected_norm.pos
            && call_norm.ref_allele == expected_norm.ref_allele
            && call_norm.alt_allele == expected_norm.alt_allele
    })
}
