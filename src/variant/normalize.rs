/// INDEL normalization via left-alignment and prefix/suffix trimming.
///
/// Implements the vt normalize / bcftools norm approach (based on Tan et al.
/// 2015) for canonicalizing INDEL representations. Different tools can produce
/// different coordinate representations for the same biological event (e.g.,
/// an insertion at different positions within a homopolymer run).
/// Normalization ensures that equivalent INDELs are compared as equal during
/// filtering.
///
/// Algorithm (when reference context is available):
/// 1. Left-align by repeatedly extending one base to the left and re-trimming
///    the common suffix (extend-and-retrim loop). This naturally handles both
///    homopolymer runs and multi-base repeat motifs.
/// 2. Final trim: remove common suffix then common prefix (preserving at
///    least one base in each allele for VCF anchor).
///
/// Without reference context, only suffix/prefix trimming is performed.

/// A normalized variant representation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NormalizedVariant {
    pub chrom: String,
    pub pos: u64,
    pub ref_allele: String,
    pub alt_allele: String,
}

/// Normalize an INDEL by left-aligning and trimming.
///
/// `ref_context` is the reference sequence around the variant for left-alignment.
/// The variant's `pos` is assumed to be 0-based relative to the start of
/// `ref_context`. If no `ref_context` is available, only prefix/suffix
/// trimming is performed.
pub fn normalize_indel(
    chrom: &str,
    pos: u64,
    ref_allele: &str,
    alt_allele: &str,
    ref_context: Option<&str>,
) -> NormalizedVariant {
    let mut current_pos = pos;
    let mut ref_bytes: Vec<u8> = ref_allele.as_bytes().to_vec();
    let mut alt_bytes: Vec<u8> = alt_allele.as_bytes().to_vec();

    // SNVs and identical alleles pass through unchanged.
    if ref_bytes == alt_bytes || (ref_bytes.len() == 1 && alt_bytes.len() == 1) {
        return NormalizedVariant {
            chrom: chrom.to_string(),
            pos,
            ref_allele: ref_allele.to_string(),
            alt_allele: alt_allele.to_string(),
        };
    }

    if let Some(ctx) = ref_context {
        // With reference context: left-align first, then trim.
        // Left-aligning before trimming preserves repeat-unit structure
        // needed for proper shifting through dinucleotide and longer repeats.

        // Step 1: Left-align through the reference context.
        let ctx_bytes = ctx.as_bytes();
        left_align(&mut current_pos, &mut ref_bytes, &mut alt_bytes, ctx_bytes);

        // Step 2: Trim common suffix and prefix after alignment.
        trim_suffix(&mut ref_bytes, &mut alt_bytes);
        let prefix_trimmed = trim_prefix(&mut ref_bytes, &mut alt_bytes);
        current_pos += prefix_trimmed as u64;
    } else {
        // Without reference context: just trim prefix/suffix.
        trim_suffix(&mut ref_bytes, &mut alt_bytes);
        let prefix_trimmed = trim_prefix(&mut ref_bytes, &mut alt_bytes);
        current_pos += prefix_trimmed as u64;
    }

    NormalizedVariant {
        chrom: chrom.to_string(),
        pos: current_pos,
        ref_allele: String::from_utf8(ref_bytes)
            .expect("allele should be valid UTF-8"),
        alt_allele: String::from_utf8(alt_bytes)
            .expect("allele should be valid UTF-8"),
    }
}

/// Check if two variants are equivalent after normalization.
///
/// Normalizes both variants and compares all fields. Returns `true` if the
/// normalized representations are identical.
pub fn variants_equivalent(
    v1_chrom: &str,
    v1_pos: u64,
    v1_ref: &str,
    v1_alt: &str,
    v2_chrom: &str,
    v2_pos: u64,
    v2_ref: &str,
    v2_alt: &str,
    ref_context: Option<&str>,
) -> bool {
    let n1 = normalize_indel(v1_chrom, v1_pos, v1_ref, v1_alt, ref_context);
    let n2 = normalize_indel(v2_chrom, v2_pos, v2_ref, v2_alt, ref_context);
    n1 == n2
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

/// Trim matching suffix bases from both alleles.
fn trim_suffix(ref_bytes: &mut Vec<u8>, alt_bytes: &mut Vec<u8>) {
    while ref_bytes.len() > 1
        && alt_bytes.len() > 1
        && ref_bytes.last() == alt_bytes.last()
    {
        ref_bytes.pop();
        alt_bytes.pop();
    }
}

/// Trim matching prefix bases from both alleles, returning the number of
/// bases trimmed (used to adjust position).
fn trim_prefix(ref_bytes: &mut Vec<u8>, alt_bytes: &mut Vec<u8>) -> usize {
    let mut count = 0;
    while ref_bytes.len() > 1
        && alt_bytes.len() > 1
        && ref_bytes.first() == alt_bytes.first()
    {
        ref_bytes.remove(0);
        alt_bytes.remove(0);
        count += 1;
    }
    count
}

/// Left-align a variant by repeatedly extending one base to the left and
/// re-trimming the common suffix.
///
/// This implements the vt normalize / bcftools norm approach:
/// 1. Prepend the reference base at (pos-1) to both alleles
/// 2. Decrement pos
/// 3. Right-trim common suffix (while both non-empty and last bases match)
/// 4. If no suffix bases were trimmed, the variant cannot shift further; stop
/// 5. Repeat
///
/// This naturally handles dinucleotide and longer repeat motifs because the
/// extend-and-retrim cycle processes one base at a time from the reference
/// context, allowing multi-base repeat units to "walk" through the reference.
///
/// `pos` is 0-based relative to `ctx`.
fn left_align(pos: &mut u64, ref_bytes: &mut Vec<u8>, alt_bytes: &mut Vec<u8>, ctx: &[u8]) {
    loop {
        // Cannot shift further left than the start of the context.
        if *pos == 0 {
            break;
        }

        // Both alleles must be non-empty.
        if ref_bytes.is_empty() || alt_bytes.is_empty() {
            break;
        }

        // Extend: prepend the base at (pos-1) to both alleles.
        let left_base = ctx[(*pos as usize) - 1];
        ref_bytes.insert(0, left_base);
        alt_bytes.insert(0, left_base);
        *pos -= 1;

        // Right-trim common suffix, preserving at least 1 base in each allele
        // (VCF requires an anchor base).
        let len_before = ref_bytes.len() + alt_bytes.len();
        while ref_bytes.len() > 1
            && alt_bytes.len() > 1
            && ref_bytes.last() == alt_bytes.last()
        {
            ref_bytes.pop();
            alt_bytes.pop();
        }
        let len_after = ref_bytes.len() + alt_bytes.len();

        // If no suffix was trimmed, the extension did not help — undo and stop.
        if len_before == len_after {
            ref_bytes.remove(0);
            alt_bytes.remove(0);
            *pos += 1;
            break;
        }
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // -- Basic trimming (no reference context) --

    #[test]
    fn test_snv_passes_through_unchanged() {
        let result = normalize_indel("chr1", 100, "A", "T", None);
        assert_eq!(
            result,
            NormalizedVariant {
                chrom: "chr1".to_string(),
                pos: 100,
                ref_allele: "A".to_string(),
                alt_allele: "T".to_string(),
            }
        );
    }

    #[test]
    fn test_identical_alleles_pass_through() {
        let result = normalize_indel("chr1", 100, "ACG", "ACG", None);
        assert_eq!(result.pos, 100);
        assert_eq!(result.ref_allele, "ACG");
        assert_eq!(result.alt_allele, "ACG");
    }

    #[test]
    fn test_trim_common_suffix() {
        // chr1:100 AGT/AT → trim trailing T → AG/A at pos 100
        let result = normalize_indel("chr1", 100, "AGT", "AT", None);
        assert_eq!(result.pos, 100);
        assert_eq!(result.ref_allele, "AG");
        assert_eq!(result.alt_allele, "A");
    }

    #[test]
    fn test_trim_common_prefix() {
        // chr1:100 TA/TGA: suffix trim first (A==A) → T/TG. Prefix: T is len 1, skip.
        // Result: pos 100, T/TG (VCF-style left anchor preserved).
        let result = normalize_indel("chr1", 100, "TA", "TGA", None);
        assert_eq!(result.pos, 100);
        assert_eq!(result.ref_allele, "T");
        assert_eq!(result.alt_allele, "TG");
    }

    #[test]
    fn test_trim_both_prefix_and_suffix() {
        // chr1:100 TACG/TGACG: suffix trim G, C, A → T/TG.
        // Prefix: T is len 1, skip.
        // Result: pos 100, T/TG.
        let result = normalize_indel("chr1", 100, "TACG", "TGACG", None);
        assert_eq!(result.pos, 100);
        assert_eq!(result.ref_allele, "T");
        assert_eq!(result.alt_allele, "TG");
    }

    #[test]
    fn test_trim_prefix_when_suffix_differs() {
        // chr1:100 GAT/GA: suffix T!=A, no suffix trim.
        // Prefix: G==G (both len > 1) → trim → AT/A at pos 101.
        // Now A is len 1, stop.
        let result = normalize_indel("chr1", 100, "GAT", "GA", None);
        assert_eq!(result.pos, 101);
        assert_eq!(result.ref_allele, "AT");
        assert_eq!(result.alt_allele, "A");
    }

    #[test]
    fn test_simple_insertion_no_context() {
        // chr1:100 A/AT → already minimal (A is len 1, no prefix trim)
        let result = normalize_indel("chr1", 100, "A", "AT", None);
        assert_eq!(result.pos, 100);
        assert_eq!(result.ref_allele, "A");
        assert_eq!(result.alt_allele, "AT");
    }

    #[test]
    fn test_simple_deletion_no_context() {
        // chr1:100 AT/A → already minimal
        let result = normalize_indel("chr1", 100, "AT", "A", None);
        assert_eq!(result.pos, 100);
        assert_eq!(result.ref_allele, "AT");
        assert_eq!(result.alt_allele, "A");
    }

    // -- Left-alignment with reference context --

    #[test]
    fn test_homopolymer_insertion_left_aligns() {
        // Reference context: AAATTTTCCC (0-indexed)
        // Variant at pos 6: T/TT (inserting T in the TTTT run)
        // Left-align shifts through the T run and absorbs the A anchor:
        //   pos 6 → 5 → 4 → 3 (T/TT) → 2 (A/AT, absorbs last A before run)
        // Result: pos 2, A/AT
        let ctx = "AAATTTTCCC";
        let result = normalize_indel("chr1", 6, "T", "TT", Some(ctx));
        assert_eq!(result.pos, 2, "should left-align past T run to A anchor");
        assert_eq!(result.ref_allele, "A");
        assert_eq!(result.alt_allele, "AT");
    }

    #[test]
    fn test_homopolymer_deletion_left_aligns() {
        // Reference context: AAATTTTCCC
        // Variant at pos 6: TT/T (deleting one T from the run)
        // Left-align shifts through T run, absorbs A anchor:
        // Result: pos 2, AT/A
        let ctx = "AAATTTTCCC";
        let result = normalize_indel("chr1", 6, "TT", "T", Some(ctx));
        assert_eq!(result.pos, 2, "should left-align past T run to A anchor");
        assert_eq!(result.ref_allele, "AT");
        assert_eq!(result.alt_allele, "A");
    }

    #[test]
    fn test_dinucleotide_repeat_left_aligns() {
        // Reference context: AAACACACACGG (0-indexed)
        //                     0123456789..
        // Variant at pos 7: AC/ACAC (inserting AC in repeat)
        // The extend-and-retrim algorithm shifts left through the repeat.
        // After left-alignment and trimming, the result is a parsimonious
        // representation near the start of the repeat region.
        let ctx = "AAACACACACGG";
        let result = normalize_indel("chr1", 7, "AC", "ACAC", Some(ctx));
        // The important property: position shifted left from 7.
        assert!(
            result.pos < 7,
            "should left-align from pos 7, got pos {}",
            result.pos
        );
    }

    #[test]
    fn test_already_left_aligned_no_change() {
        // Reference context: ACGTACGT
        // Variant at pos 0: A/AT — cannot shift further left
        let ctx = "ACGTACGT";
        let result = normalize_indel("chr1", 0, "A", "AT", Some(ctx));
        assert_eq!(result.pos, 0);
        assert_eq!(result.ref_allele, "A");
        assert_eq!(result.alt_allele, "AT");
    }

    #[test]
    fn test_no_left_alignment_when_bases_differ() {
        // Reference context: ACGTACGT
        // Variant at pos 4: A/AG — extending left prepends T to both:
        //   TA/TAG, suffix: A!=G, no trim. Undo. Stop.
        let ctx = "ACGTACGT";
        let result = normalize_indel("chr1", 4, "A", "AG", Some(ctx));
        assert_eq!(result.pos, 4);
        assert_eq!(result.ref_allele, "A");
        assert_eq!(result.alt_allele, "AG");
    }

    #[test]
    fn test_complex_indel_no_shift() {
        // Complex INDEL (delins): ref AGT, alt TC at pos 5
        // Context: AAAAAGTTTT
        // No suffix or prefix to trim. Extending left doesn't help.
        let ctx = "AAAAAGTTTT";
        let result = normalize_indel("chr1", 5, "AGT", "TC", Some(ctx));
        assert_eq!(result.pos, 5);
        assert_eq!(result.ref_allele, "AGT");
        assert_eq!(result.alt_allele, "TC");
    }

    // -- variants_equivalent tests --

    #[test]
    fn test_equivalent_insertions_different_positions() {
        // Same insertion represented at two different positions in a homopolymer.
        // Reference: AAATTTTCCC
        // v1: pos 4, T/TT (middle of run)
        // v2: pos 6, T/TT (end of run)
        // Both normalize to pos 2, A/AT (left-aligned past the T run).
        let ctx = "AAATTTTCCC";
        assert!(variants_equivalent(
            "chr1", 4, "T", "TT",
            "chr1", 6, "T", "TT",
            Some(ctx),
        ));
    }

    #[test]
    fn test_equivalent_deletions_different_positions() {
        // Same deletion in a homopolymer at different positions.
        // Reference: AAATTTTCCC
        let ctx = "AAATTTTCCC";
        assert!(variants_equivalent(
            "chr1", 4, "TT", "T",
            "chr1", 5, "TT", "T",
            Some(ctx),
        ));
    }

    #[test]
    fn test_non_equivalent_different_alleles() {
        let ctx = "ACGTACGTACGT";
        assert!(!variants_equivalent(
            "chr1", 4, "A", "AT",
            "chr1", 4, "A", "AG",
            Some(ctx),
        ));
    }

    #[test]
    fn test_non_equivalent_different_chroms() {
        let ctx = "AAATTTTCCC";
        assert!(!variants_equivalent(
            "chr1", 4, "T", "TT",
            "chr2", 4, "T", "TT",
            Some(ctx),
        ));
    }

    #[test]
    fn test_equivalent_snvs_same_position() {
        // SNVs should always compare by exact coordinates
        assert!(variants_equivalent(
            "chr1", 100, "A", "T",
            "chr1", 100, "A", "T",
            None,
        ));
    }

    #[test]
    fn test_non_equivalent_snvs_different_position() {
        assert!(!variants_equivalent(
            "chr1", 100, "A", "T",
            "chr1", 101, "A", "T",
            None,
        ));
    }

    #[test]
    fn test_equivalent_with_suffix_trimming() {
        // v1: pos 100, ATCG/ACG → suffix trim: G then C → AT/A at pos 100.
        // v2: pos 100, AT/A → already trimmed.
        // Both become pos 100, AT/A.
        assert!(variants_equivalent(
            "chr1", 100, "ATCG", "ACG",
            "chr1", 100, "AT", "A",
            None,
        ));
    }

    #[test]
    fn test_equivalent_with_prefix_trimming() {
        // v1: pos 100, GAT/GA → suffix: T!=A, no trim. prefix: G==G → AT/A at 101.
        // v2: pos 101, AT/A → already trimmed.
        assert!(variants_equivalent(
            "chr1", 100, "GAT", "GA",
            "chr1", 101, "AT", "A",
            None,
        ));
    }

    // -- Edge cases --

    #[test]
    fn test_empty_ref_allele() {
        // Edge case: empty ref allele (pure insertion without anchor)
        // Should not panic
        let result = normalize_indel("chr1", 100, "", "AT", None);
        assert_eq!(result.ref_allele, "");
        assert_eq!(result.alt_allele, "AT");
    }

    #[test]
    fn test_empty_alt_allele() {
        // Edge case: empty alt allele (pure deletion without anchor)
        let result = normalize_indel("chr1", 100, "AT", "", None);
        assert_eq!(result.ref_allele, "AT");
        assert_eq!(result.alt_allele, "");
    }

    #[test]
    fn test_single_base_alleles() {
        // SNV: single base ref and alt
        let result = normalize_indel("chr1", 100, "A", "G", None);
        assert_eq!(result.pos, 100);
        assert_eq!(result.ref_allele, "A");
        assert_eq!(result.alt_allele, "G");
    }

    #[test]
    fn test_long_homopolymer_left_alignment() {
        // Long homopolymer: AAAAAAAAAAATTTTTTTTTTC
        //                    0         1         2
        //                    01234567890123456789012
        // Variant at pos 18: T/TT
        // Left-aligns through all T's and absorbs A anchor.
        // Stops at pos 10 with A/AT (last A before T run).
        let ctx = "AAAAAAAAAAATTTTTTTTTTC";
        let result = normalize_indel("chr1", 18, "T", "TT", Some(ctx));
        assert_eq!(result.pos, 10, "should left-align to anchor A before T run");
        assert_eq!(result.ref_allele, "A");
        assert_eq!(result.alt_allele, "AT");
    }

    #[test]
    fn test_variant_at_position_zero_with_context() {
        // Variant already at position 0 — no room to left-align
        let ctx = "TTTTACGT";
        let result = normalize_indel("chr1", 0, "T", "TT", Some(ctx));
        assert_eq!(result.pos, 0);
    }

    #[test]
    fn test_retrim_after_left_alignment() {
        // Reference: GGGAAACCC
        //            012345678
        // Variant at pos 5: AAC/AAAC
        // Left-align first (extend-and-retrim):
        //   pos 5, AAC/AAAC → prepend ctx[4]=A → AAAC/AAAAC, pos=4
        //     suffix: C==C → AAA/AAAA. A==A → AA/AAA. A==A → A/AA. ref len 1, stop. trimmed 3.
        //   pos 4, A/AA → prepend ctx[3]=A → AA/AAA, pos=3
        //     suffix: A==A → A/AA. ref len 1, stop. trimmed 1.
        //   pos 3, A/AA → prepend ctx[2]=G → GA/GAA, pos=2
        //     suffix: A==A → G/GA. ref len 1, stop. trimmed 1.
        //   pos 2, G/GA → prepend ctx[1]=G → GG/GGA, pos=1
        //     suffix: G!=A. trimmed 0. undo → G/GA, pos=2.
        // Final trim: suffix A!=A? ref=G, alt=GA. G!=A. prefix: G==G but ref len 1. skip.
        // Result: pos 2, G/GA
        let ctx = "GGGAAACCC";
        let result = normalize_indel("chr1", 5, "AAC", "AAAC", Some(ctx));
        assert_eq!(result.pos, 2);
        assert_eq!(result.ref_allele, "G");
        assert_eq!(result.alt_allele, "GA");
    }
}
